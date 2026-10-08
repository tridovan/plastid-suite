"""The plastid calling path — the plastid half of the shared core (P3B-07, D-e).

An assembled plastome arrives as a bare sequence, so its genes are called by **transfer** from
curated relatives: `blastx` against the reference set's own published ``/translation`` for the
CDS, `blastn` against its rRNA sequences for the four rRNA loci, and — where the frame is needed
— a self-search for the inverted repeat. The calling rules themselves live in
``core/caller.py``; this module is what is **measured about this organelle**:

* **the reference index.** A plastid reference record is a GenBank file, so a CDS entry is that
  CDS's own ``/translation`` (never a self-translation — D5-08 §8.21, and the mito measurement
  behind it: plant organellar CDS begin with codons no table accepts as an initiator) and an rRNA
  entry is the feature's own nucleotide sequence, read through ``coding_sequence`` so a spliced
  feature comes back in transcript order.
* **the floors** — identity ≥ 50 %, coverage ≥ 0.7 — measured on the golden as the answer key
  (PROJECT.md §8.46). They are deliberately **not** the mito path's (50 % / 0.5): §8.41 measured
  that a short plastid CDS walks through a 0.5 coverage floor (`psbJ` is 36 aa and reaches 52.38 %
  against `matK` over 21 aa), and each path states the pair it measured.
* **the translation table** — 11, the plastid code, which is what every CDS in this project's
  measured records declares (D3-03: 1,099 of 1,099).

**The vocabulary is the reference set's own.** D5-09 puts the naming step in front of the
annotator, so the reference DB arrives already spelled in vocabulary A (the ★ record's spelling) —
which makes the set's own names the only vocabulary a plastid call may print. A reference gene the
set spells with a name outside the vocabulary is kept and **reported** rather than renamed
(H-06's rule, the mito path's ``unnameable``).

**The inverted repeat is read behind the structural gate, never instead of it.** §9.A is the
measured reason: naive IR detection got 2 of 9 published plastomes wrong, one of them silently.
So the frame is established by ``core/ir.py`` and the ≥ 99.9 % gate, and a self-search only has to
**agree** with it — :func:`read_inverted_repeat` turns a self-blast's hits into the repeat pair,
and :func:`repeat_agreement` says so when the two readings differ rather than quietly preferring
one.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass

from plastid.core import caller as _core
from plastid.core.caller import (
    Hsp,
    ReferenceEntry,
    ReferenceIndex,
)
from plastid.core.gene_names import AliasTable, aliases_for, normalize_name
from plastid.core.ir import IRResult
from plastid.core.model import Genome, coding_sequence

__all__ = [
    "DEFAULT_MIN_COVERAGE",
    "DEFAULT_MIN_IDENTITY",
    "IR_AGREEMENT_TOLERANCE",
    "IR_MIN_IDENTITY",
    "PLASTID_TRANSL_TABLE",
    "InvertedRepeat",
    "call_cds",
    "call_rrna",
    "read_inverted_repeat",
    "reference_entries",
    "reference_index",
    "repeat_agreement",
]

# The plastid code. Every CDS in the records this project measured declares it (D3-03).
PLASTID_TRANSL_TABLE = 11

# The floors, measured on the golden as the answer key (§8.46) — see the module docstring for why
# they are not the mito path's pair.
DEFAULT_MIN_IDENTITY = 50.0
DEFAULT_MIN_COVERAGE = 0.7

# What a self-search must show before two spans are called one inverted repeat. 99 % is PGA's own
# threshold (D5-08 §1) and sits below the ≥ 99.9 % the structural gate asks of the copies
# themselves: this reading proposes the frame, it never certifies it.
IR_MIN_IDENTITY = 99.0
# A plastid IR is ~10–30 kb; 1 kb keeps a short duplicated tRNA from being read as the repeat
# without ever touching a real one (measured: the golden's IR is 26,214 bp, the smallest in the
# pinned corpus is 20,324 bp).
IR_MIN_LENGTH = 1000
# How far the two readings of the repeat may drift apart and still be one reading — see
# :func:`repeat_agreement`. One percent of a 20–30 kb repeat is a few hundred bases, which is
# roughly what a search trims off each end of a run of near-identical repeats.
IR_AGREEMENT_TOLERANCE = 0.01


def reference_entries(
    records: Iterable[Genome], *, kind: str, vocabulary: Sequence[str] | None = None
) -> tuple[ReferenceEntry, ...]:
    """One record's genes of one ``kind`` as reference entries, with the sequence to search.

    A CDS entry carries the record's own ``/translation`` — never a self-translation (D5-08
    §8.21): no protein means no protein to search with, and the audit in front of every run
    blocks such a record before it gets here. An rRNA entry carries the feature's own nucleotide
    sequence, read through ``coding_sequence`` so a spliced feature comes back in transcript
    order.

    ``vocabulary`` limits what a call may be *named*: a gene whose canonical name is not in it is
    kept with ``vocabulary=None``, which the caller reports as ``unnameable`` instead of printing
    a name the project's vocabulary does not use. ``None`` means "everything this record spells",
    which is the right answer for a reference set that has already been through the naming step.
    """
    if kind not in ("CDS", "rRNA"):
        raise ValueError(f"kind must be 'CDS' or 'rRNA', got {kind!r}")
    allowed = None if vocabulary is None else frozenset(vocabulary)
    entries: list[ReferenceEntry] = []
    for genome in records:
        accession = genome.accession or genome.name
        table = aliases_for(genome)
        for gene in genome.genes:
            if gene.kind != kind or gene.pseudo:
                continue
            if kind == "CDS":
                if not gene.translation:
                    continue
                sequence = gene.translation
            else:
                sequence = coding_sequence(genome, gene)
            name = normalize_name(gene.name, table).canonical
            entries.append(
                ReferenceEntry(
                    accession=accession,
                    organism=genome.name,
                    source_gene=gene.name,
                    vocabulary=name if allowed is None or name in allowed else None,
                    sequence=sequence,
                )
            )
    return tuple(entries)


def reference_index(
    records: Iterable[Genome],
    *,
    kind: str,
    vocabulary: Sequence[str],
) -> ReferenceIndex:
    """The index one search runs against: ``kind`` picks the CDS or the rRNA half of the set.

    The two are separate searches against separate databases because they are separate questions
    — a protein against proteins, a nucleotide against nucleotides — and pooling them would put a
    ``blastn`` percentage and a ``blastx`` one on the same scale, which is the comparison D5-08 §3
    forbids.
    """
    return ReferenceIndex.from_entries(
        reference_entries(records, kind=kind, vocabulary=vocabulary), vocabulary=vocabulary
    )


def vocabulary_of(records: Iterable[Genome], *, kind: str) -> tuple[str, ...]:
    """The canonical names this set spells for one kind, in first-seen order, de-duplicated."""
    names: dict[str, None] = {}
    for genome in records:
        table: AliasTable = aliases_for(genome)
        for gene in genome.genes:
            if gene.kind != kind or gene.pseudo:
                continue
            names.setdefault(normalize_name(gene.name, table).canonical, None)
    return tuple(names)


def call_cds(
    genome: Genome,
    references: ReferenceIndex,
    hsps: Sequence[Hsp],
    *,
    min_identity: float = DEFAULT_MIN_IDENTITY,
    min_coverage: float = DEFAULT_MIN_COVERAGE,
    transl_table: int = PLASTID_TRANSL_TABLE,
    max_hsps_per_subject: int | None = None,
) -> _core.CallResult:
    """Call the plastome's protein-coding genes from a blastx run's hits."""
    return _core.call_cds(
        genome,
        references,
        hsps,
        min_identity=min_identity,
        min_coverage=min_coverage,
        transl_table=transl_table,
        max_hsps_per_subject=max_hsps_per_subject,
    )


def call_rrna(
    genome: Genome,
    references: ReferenceIndex,
    hsps: Sequence[Hsp],
    *,
    min_identity: float = DEFAULT_MIN_IDENTITY,
    min_coverage: float = DEFAULT_MIN_COVERAGE,
    max_hsps_per_subject: int | None = None,
) -> _core.CallResult:
    """Call the plastome's rRNA genes from a blastn run's hits."""
    return _core.call_rrna(
        genome,
        references,
        hsps,
        min_identity=min_identity,
        min_coverage=min_coverage,
        max_hsps_per_subject=max_hsps_per_subject,
    )


# ---------------------------------------------------------------------------------------------
# The inverted repeat, as a self-search reads it (D5-08 §1 layer 1, behind core/ir.py's gate)
# ---------------------------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class InvertedRepeat:
    """Two spans that read as one repeat, and the identity the search measured between them.

    Coordinates are 1-based inclusive in the record as it was searched — the self-search sees the
    same sequence the caller does, so no rotation is applied here.
    """

    first: tuple[int, int]
    second: tuple[int, int]
    identity: float
    length: int

    @property
    def span(self) -> int:
        return self.length


def read_inverted_repeat(
    hsps: Sequence[Hsp],
    *,
    min_identity: float = IR_MIN_IDENTITY,
    min_length: int = IR_MIN_LENGTH,
) -> InvertedRepeat | None:
    """The largest self-hit that puts two spans of one sequence in opposite orientation.

    A self-search returns the diagonal first — every position matching itself — so the repeat is
    the *anti*-diagonal: a hit whose subject runs backwards relative to the query
    (``sstart > send``), long enough and identical enough to be a copy rather than a coincidence.
    The two spans are the query's and the subject's, which is what a ``repeat_region`` feature
    needs; the identity is that hit's own, never an average over the hits that overlap it.
    """
    candidates = [
        hsp
        for hsp in hsps
        if hsp.sstart > hsp.send
        and hsp.pident >= min_identity
        and hsp.length >= min_length
        and not _is_diagonal(hsp)
    ]
    if not candidates:
        return None
    best = max(candidates, key=lambda hsp: (hsp.length, hsp.pident))
    return InvertedRepeat(
        first=(best.qlo, best.qhi),
        second=(best.slo, best.shi),
        identity=best.pident,
        length=best.length,
    )


def _is_diagonal(hsp: Hsp) -> bool:
    """True when the hit is the sequence matching itself — the one hit that is not a repeat.

    A self-search always reports it, and it is the trap this function exists for: it is the
    longest, most identical hit in the set, so a rule that merely takes the best hit would name
    the whole sequence its own inverted repeat.
    """
    return abs(hsp.qlo - hsp.slo) <= 3 and abs(hsp.qhi - hsp.shi) <= 3


def repeat_agreement(structural: IRResult | None, searched: InvertedRepeat | None) -> str:
    """One sentence saying whether the two readings of the repeat agree — and how they differ.

    The structural reading is the one that governs: it is what the ≥ 99.9 % gate certifies and
    what every position-based module is normalized against (§9.A). The self-search is a second
    opinion, and a disagreement is a **finding** — never a reason to prefer one silently.
    """
    if structural is None:
        return (
            "no inverted repeat stands in the structural reading, so the self-search was not "
            "asked to confirm one"
        )
    if searched is None:
        return (
            f"the structural reading finds an inverted repeat of {structural.length} bp, and the "
            f"self-search returned no anti-diagonal hit at the floors measured for it — reported "
            f"rather than resolved, because the gate, not the search, is what positions are "
            f"normalized against"
        )
    # A tolerance, because the two readings measure different things: the structural frame is the
    # exact span between the two copies' boundaries, and the search reports the bases its own
    # alignment covered — a few either side of it. Calling that a disagreement would make every
    # real record report one, and a finding that always fires is a finding nobody reads.
    if abs(searched.length - structural.length) <= IR_AGREEMENT_TOLERANCE * structural.length:
        return (
            f"the self-search agrees with the structural reading: a {searched.length} bp repeat "
            f"at {searched.identity:.2f} % identity against the frame's {structural.length} bp"
        )
    return (
        f"the two readings differ: the structural frame finds {structural.length} bp and the "
        f"self-search {searched.length} bp ({(searched.first[0], searched.first[1])} against "
        f"{(searched.second[0], searched.second[1])}, {searched.identity:.2f} %) — reported, and "
        f"the structural frame is the one kept"
    )
