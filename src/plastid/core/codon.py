"""The genetic code and codon-usage measures — pure functions, no I/O.

Compiled table: **NCBI table 1 (standard)** — decision D2-03. Plastid CDS declare **table 11**,
whose amino-acid assignments are identical (the two tables differ only in initiation codons), so
one code table answers both — and :func:`translate` renders every codon, the initiator included,
as its own amino acid. Measured on the golden genome: 5 of its 86 CDS start with GTG/ACG/ATC,
so a published protein begins with M there; a report that prints a translation
table must therefore say which table it used *and* which one the record declares —
:func:`translation_table_label` is what builds that string.

Every measure states its reference: RSCU and ENC are defined by the sequences alone, while CAI,
CBI and Fop are only meaningful against a declared reference set — the caller passes one and
records it, so a value never travels without the thing it was measured against.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping

from plastid.core.model import Gene

__all__ = [
    "CODONS",
    "FAMILIES",
    "SENSE_FAMILIES",
    "SINGLE_CODON_AMINO",
    "STANDARD_CODE",
    "SUPPORTED_TRANSL_TABLES",
    "UnsupportedTranslationError",
    "cai",
    "cbi",
    "codon_counts",
    "enc",
    "fop",
    "optimal_codons",
    "relative_adaptiveness",
    "require_translatable",
    "rscu",
    "translate",
]

# NCBI translation table 1, written in the conventional T/C/A/G order.
_CODE_TABLE = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_BASES = "TCAG"
CODONS: tuple[str, ...] = tuple(
    first + second + third for first in _BASES for second in _BASES for third in _BASES
)
STANDARD_CODE: dict[str, str] = dict(zip(CODONS, _CODE_TABLE, strict=True))

_MEMBERS: dict[str, list[str]] = {}
for _codon, _amino in STANDARD_CODE.items():
    _MEMBERS.setdefault(_amino, []).append(_codon)
FAMILIES: dict[str, tuple[str, ...]] = {
    amino: tuple(sorted(members)) for amino, members in sorted(_MEMBERS.items())
}
SENSE_FAMILIES: dict[str, tuple[str, ...]] = {
    amino: members for amino, members in FAMILIES.items() if amino != "*"
}
SENSE_CODONS: tuple[str, ...] = tuple(
    codon for members in SENSE_FAMILIES.values() for codon in members
)
SINGLE_CODON_AMINO = frozenset(
    amino for amino, members in SENSE_FAMILIES.items() if len(members) == 1
)


# What a report's "translation table" column says the toolkit used (D2-03). It is not the whole
# truth about a record: plastid CDS declare table 11 — same amino acids, different initiators —
# which is why the column is built by translation_table_label() and not by this constant alone.
TABLE_LABEL = "NCBI table 1 (standard)"


def translation_table_label(genes: Iterable[Gene]) -> str:
    """``TABLE_LABEL`` plus the table(s) the record's own CDS declare (D3-03).

    A module that compiles one code table still has to name the record's, or the column reads as
    if the two were the same table. A record whose CDS declare nothing (silence is not a value)
    leaves the label as the compiled table alone.
    """
    declared = sorted({gene.transl_table for gene in genes if gene.transl_table is not None})
    if not declared:
        return TABLE_LABEL
    return f"{TABLE_LABEL}; the record's CDS declare {', '.join(str(table) for table in declared)}"


# The compiled code is NCBI table 1. Table 11 (the plastid code) shares every amino-acid
# assignment and differs only in initiation codons; every other table assigns at least one codon
# differently, so a codon count or composition under it would be a different measurement wearing
# this module's name.
SUPPORTED_TRANSL_TABLES = frozenset({1, 11})


class UnsupportedTranslationError(ValueError):
    """A CDS declares a frame or table this codon reading cannot translate honestly."""


def require_translatable(genes: Iterable[Gene], *, module: str) -> None:
    """Refuse a CDS set whose frame or declared table this module cannot read as compiled.

    The modules behind this guard translate with NCBI table 1 and read from the first base of
    the feature. A ``/codon_start`` other than 1 shifts the frame, and another table assigns
    codons differently; in either case silently translating would produce a number that looks
    like the module's and is not it. The message names the genes, so the fix is a decision about
    the record rather than a search.
    """
    genes = tuple(genes)
    wrong_frame = [gene for gene in genes if gene.codon_start != 1]
    wrong_table = [
        gene
        for gene in genes
        if gene.transl_table is not None and gene.transl_table not in SUPPORTED_TRANSL_TABLES
    ]
    if wrong_frame:
        named = ", ".join(
            f"{gene.name} (/codon_start={gene.codon_start})" for gene in wrong_frame[:8]
        )
        more = "" if len(wrong_frame) <= 8 else f" and {len(wrong_frame) - 8} more"
        raise UnsupportedTranslationError(
            f"{module}: {len(wrong_frame)} CDS do not start in frame 1 — {named}{more}. "
            f"This module reads from the feature's first base, so translating them would count "
            f"the wrong frame; only /codon_start=1 is implemented"
        )
    if wrong_table:
        tables = sorted({gene.transl_table for gene in wrong_table})
        named = ", ".join(f"{gene.name} (table {gene.transl_table})" for gene in wrong_table[:8])
        more = "" if len(wrong_table) <= 8 else f" and {len(wrong_table) - 8} more"
        raise UnsupportedTranslationError(
            f"{module}: {len(wrong_table)} CDS declare translation table(s) {tables} — "
            f"{named}{more}. The compiled code is NCBI table 1, which shares its amino-acid "
            f"assignments with plastid table 11 only; another table would be a different "
            f"measurement, so it is refused rather than translated"
        )


def translate(sequence: str) -> str:
    """Translate complete codons; a trailing partial codon is ignored, ``*`` marks a stop.

    Every codon is read as its own amino acid — the initiator included, so a CDS annotated with a
    GTG/ACG/ATC start translates to V/T/I there where a published protein has M.
    """
    upper = sequence.upper().replace("U", "T")
    return "".join(
        STANDARD_CODE.get(upper[index : index + 3], "X") for index in range(0, len(upper) - 2, 3)
    )


def codon_counts(sequence: str) -> dict[str, int]:
    """Counts for every sense codon; stops, partial trailing codons and ambiguous bases drop out."""
    counts = dict.fromkeys(SENSE_CODONS, 0)
    upper = sequence.upper().replace("U", "T")
    for index in range(0, len(upper) - 2, 3):
        codon = upper[index : index + 3]
        if codon in counts:
            counts[codon] += 1
    return counts


def rscu(counts: Mapping[str, int]) -> dict[str, float]:
    """Relative synonymous codon usage: observed / expected-if-uniform within the family.

    ``sum(rscu[codon] for codon in family) == family size`` for every amino acid — the check
    P2-03 asserts. A family with no observations yields ``0.0`` for all of its codons.
    """
    values: dict[str, float] = {}
    for members in SENSE_FAMILIES.values():
        total = sum(counts.get(codon, 0) for codon in members)
        expected = total / len(members)
        for codon in members:
            values[codon] = counts.get(codon, 0) / expected if expected else 0.0
    return values


def enc(counts: Mapping[str, int]) -> float | None:
    """Effective number of codons (Wright 1990): ``2 + 9/F2 + 1/F3 + 5/F4 + 3/F6``.

    ``F_k`` is the homozygosity ``(n·Σp² − 1)/(n − 1)`` of one family of ``k`` codons, averaged
    over the families of that degeneracy. ``61`` means no bias, ``20`` the maximum. ``None``
    when a degeneracy class has no usable family — the value would not be defined.
    """
    by_degeneracy: dict[int, list[float]] = {2: [], 3: [], 4: [], 6: []}
    for members in SENSE_FAMILIES.values():
        size = len(members)
        if size not in by_degeneracy:
            continue
        total = sum(counts.get(codon, 0) for codon in members)
        if total < 2:
            continue
        share_squared = sum((counts.get(codon, 0) / total) ** 2 for codon in members)
        by_degeneracy[size].append((total * share_squared - 1) / (total - 1))

    value = 2.0
    for size, weight in ((2, 9.0), (3, 1.0), (4, 5.0), (6, 3.0)):
        families = by_degeneracy[size]
        if not families:
            return None
        mean = sum(families) / len(families)
        if mean <= 0:
            return None  # a perfectly even degeneracy class would send the term to infinity
        value += weight / mean
    return value


def relative_adaptiveness(reference: Mapping[str, int]) -> dict[str, float]:
    """Sharp & Li ``w`` per sense codon: its reference count over the family's maximum."""
    values: dict[str, float] = {}
    for members in SENSE_FAMILIES.values():
        peak = max(reference.get(codon, 0) for codon in members)
        for codon in members:
            values[codon] = reference.get(codon, 0) / peak if peak else 0.0
    return values


def optimal_codons(reference: Mapping[str, int]) -> frozenset[str]:
    """Codons at their family's reference maximum — ties all count, so the set is deterministic."""
    best: set[str] = set()
    for members in SENSE_FAMILIES.values():
        if len(members) == 1:
            continue
        peak = max(reference.get(codon, 0) for codon in members)
        if peak:
            best.update(codon for codon in members if reference.get(codon, 0) == peak)
    return frozenset(best)


def _informative(counts: Mapping[str, int]) -> list[tuple[str, int]]:
    """Observed sense codons in families that offer a choice — Met and Trp carry no signal."""
    return [
        (codon, count)
        for codon, count in counts.items()
        if count and codon in STANDARD_CODE and STANDARD_CODE[codon] not in SINGLE_CODON_AMINO
    ]


def cai(counts: Mapping[str, int], reference: Mapping[str, int]) -> float | None:
    """Codon adaptation index against a declared reference; ``None`` when nothing qualifies.

    Codons whose reference weight is zero are skipped: including them would drive the geometric
    mean to zero for every sequence and hide the comparison the index exists to make.
    """
    weights = relative_adaptiveness(reference)
    logs = [math.log(weights[codon]) for codon, _ in _informative(counts) if weights[codon] > 0]
    if not logs:
        return None
    return math.exp(sum(logs) / len(logs))


def fop(counts: Mapping[str, int], reference: Mapping[str, int]) -> float | None:
    """Frequency of optimal codons: optimal / all codons in families that offer a choice."""
    optimal = optimal_codons(reference)
    observed = _informative(counts)
    total = sum(count for _, count in observed)
    if not total:
        return None
    return sum(count for codon, count in observed if codon in optimal) / total


def cbi(counts: Mapping[str, int], reference: Mapping[str, int]) -> float | None:
    """Codon bias index (Bennetzen & Hall 1982) against the same declared reference.

    ``B = (Σ_optimal observed − Σ_optimal expected) / (N − Σ_optimal expected)``, where the
    expected count of a codon is its family's total over the family's size. ``1`` means every
    codon is optimal, ``0`` means the optimal codons appear no more often than chance.
    """
    optimal = optimal_codons(reference)
    observed = 0
    expected = 0.0
    total = 0
    for members in SENSE_FAMILIES.values():
        if len(members) == 1:
            continue
        family_total = sum(counts.get(codon, 0) for codon in members)
        total += family_total
        share = family_total / len(members)
        for codon in members:
            if codon in optimal:
                observed += counts.get(codon, 0)
                expected += share
    if total == 0 or total == expected:
        return None
    return (observed - expected) / (total - expected)
