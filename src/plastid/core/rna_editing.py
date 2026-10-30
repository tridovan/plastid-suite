"""PREP-Cp-compatible RNA-editing prediction — the published PREP method, pure functions (M4.5).

The method is Mower 2005's (PREP-Mt) applied to a chloroplast reference rebuilt the way Mower
2009 built PREP-Cp's: translate a query CDS with the standard code, align the translation to a
fixed alignment of **mature** (edited) reference proteins, and, for every codon that contains a
C, ask whether a C-to-U edit would make the codon's amino acid agree with more of the reference
residues in that column. The original PREP-Cp cannot be run any more (its service ended in 2022
and its alignments were never published), so every site here is labelled
``predicted · PREP-Cp-compatible`` and is never presented as original PREP-Cp output. The
specification, the evidence and the acceptance set are in ``docs/tools/prep_cp.md``.

What is decided here, and why:

* **Scores are exact fractions.** A state's score is ``matching / compared`` over the reference
  residues of one alignment column, gaps excluded — the archived genuine PREP-Cp output mixes
  fifths and quarters inside one gene, which only a per-column denominator produces. The cutoff
  is compared as a fraction too: ``0.8`` parsed as a binary float is slightly *more* than 4/5, so
  a float comparison would drop exactly the 0.80 site the archived output reports at C = 0.8.
* **The threshold is ``score >= C``** (Mower 2005's wording, and the archived output reports a
  0.80 site at C = 0.8).
* **Ties.** Highest score wins; then the state with fewer edits (so a silent edit is never
  predicted); then the state edited at codon position 2. The papers leave a further tie open —
  equally edited states at positions 1 and 3 — and it cannot change a call: under the standard
  code a C-to-U change at codon position 3 never changes the amino acid, so any state carrying a
  position-3 edit ties with the same state without it, which has fewer edits and wins. The final
  key (the edited positions themselves) only makes the order total; a test enumerates every codon
  to keep the argument true.

Coordinates: transcript (CDS) positions are 1-based from the feature's first base, in transcript
order, as GenBank counts them; codon and amino-acid numbers count from ``/codon_start``; genome
positions are converted to 1-based only through :meth:`Interval.to_1based`'s convention (the
mapping carries 0-based positions and the caller adds one when it exports).
"""

from __future__ import annotations

import csv
import hashlib
import json
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass, field, replace
from decimal import Decimal, InvalidOperation
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from typing import Any

import yaml

from plastid.core.codon import STANDARD_CODE
from plastid.core.gene_names import AliasTable, NormalizedName, normalize_genome
from plastid.core.model import Gene, Genome, Interval, coding_sequence
from plastid.core.resources import data_file

__all__ = [
    "AMBIGUOUS",
    "DEFAULT_CUTOFF",
    "EDIT_SITES_FILE",
    "METHOD_LABEL",
    "REFERENCE_DIRECTORY",
    "REFERENCE_FILE",
    "REFERENCE_LABEL",
    "REFERENCE_SCHEMA",
    "SCOPE_LABEL",
    "STATUS_ABSENT",
    "STATUS_EVALUATED",
    "STATUS_NOT_EVALUATED",
    "STATUS_NOT_SUPPORTED",
    "SUPPORTED_GENES",
    "TIE_RULE",
    "AddAligner",
    "CodonCall",
    "CodonState",
    "EditSite",
    "EditingReference",
    "EditingReferenceError",
    "LocusOutcome",
    "LocusPrediction",
    "PreparedLocus",
    "QUERY_ROW",
    "ReferenceRow",
    "ScoredSite",
    "apply_taxon_names",
    "build_reference",
    "call_codon",
    "check_site",
    "codon_states",
    "column_counts",
    "evaluate_genome",
    "load_build_inputs",
    "load_reference",
    "location_text",
    "locus_segments",
    "parse_location_text",
    "parse_cutoff",
    "predict_locus",
    "prepare_locus",
    "read_edit_sites",
    "reference_digest",
    "reference_from_payload",
    "reference_row_name",
    "reference_set_entry",
    "supported_gene",
    "transcript_positions",
    "translate_codons",
]

# The 35 genes PREP-Cp offered, in the order the specification lists them. Any other gene is
# "not supported" — never "0 sites", because nothing was asked about it.
SUPPORTED_GENES: tuple[str, ...] = (
    "accD",
    "atpA",
    "atpB",
    "atpF",
    "atpI",
    "ccsA",
    "clpP",
    "matK",
    "ndhA",
    "ndhB",
    "ndhD",
    "ndhF",
    "ndhG",
    "petB",
    "petD",
    "petG",
    "petL",
    "psaB",
    "psaI",
    "psbB",
    "psbE",
    "psbF",
    "psbL",
    "rpl2",
    "rpl20",
    "rpl23",
    "rpoA",
    "rpoB",
    "rpoC1",
    "rpoC2",
    "rps2",
    "rps8",
    "rps14",
    "rps16",
    "ycf3",
)

# Mower 2009's optimum for seed-plant chloroplasts (MCC 0.756 leave-one-out, 0.819 unseen). It is
# recorded per run and never moved to reach a site count.
DEFAULT_CUTOFF = "0.8"

METHOD_LABEL = "PREP-Cp-compatible"
SCOPE_LABEL = f"predicted · {METHOD_LABEL}"
TIE_RULE = (
    "highest score; tie: fewer edits; tie: edit at codon position 2 (a position-3 C-to-U is "
    "synonymous, so a remaining tie cannot change a call)"
)

STATUS_EVALUATED = "evaluated"
STATUS_NOT_SUPPORTED = "not supported"
STATUS_ABSENT = "absent from the record"
STATUS_NOT_EVALUATED = "not evaluated"

REFERENCE_SCHEMA = "plastid.rna_editing.reference/1"
# The label every B′ build carries: what it is, and what it is not.
REFERENCE_LABEL = (
    "B′ — PTS rebuild of the PREP-Cp training reference (seven taxa, mature proteins of the 35 "
    "genes); not the original PREP-Cp reference"
)
# The packaged reference B′ and the build inputs beside it (data/rna_editing/).
REFERENCE_DIRECTORY = "rna_editing"
REFERENCE_FILE = f"{REFERENCE_DIRECTORY}/reference_b_prime.json"
EDIT_SITES_FILE = f"{REFERENCE_DIRECTORY}/edit_sites.tsv"

# A base the scorer can read. Anything else in a codon (an IUPAC ambiguity code, a gap) makes the
# codon unreadable: it is listed, not guessed.
_BASES = frozenset("ACGT")
AMBIGUOUS = "ambiguous"


class EditingReferenceError(ValueError):
    """Reference B′ is malformed, does not match its digest, or contradicts a pinned record."""


# --- the scoring rule --------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class CodonState:
    """One reading of a codon: the bases, which C positions (1–3) were edited, the amino acid."""

    codon: str
    edited: tuple[int, ...]
    amino_acid: str


@dataclass(frozen=True, slots=True)
class CodonCall:
    """The state chosen for one codon against one alignment column, with what it was scored on."""

    state: CodonState
    matching: int
    compared: int
    unedited: CodonState
    unedited_matching: int

    @property
    def score(self) -> Fraction:
        return Fraction(self.matching, self.compared)

    @property
    def edited(self) -> bool:
        return bool(self.state.edited)


def codon_states(codon: str) -> tuple[CodonState, ...]:
    """The unedited codon and every C-to-U reading of it (``CCG`` → CCG, TCG, CTG, TTG).

    DNA letters are used throughout (``T`` for U) so a state compares directly with the code
    table; the edited positions are 1-based codon positions.
    """
    upper = codon.upper()
    if len(upper) != 3 or set(upper) - _BASES:
        raise ValueError(f"codon {codon!r} is not three of A, C, G, T")
    positions = [index for index, base in enumerate(upper) if base == "C"]
    states = []
    for size in range(len(positions) + 1):
        for chosen in combinations(positions, size):
            bases = list(upper)
            for index in chosen:
                bases[index] = "T"
            text = "".join(bases)
            states.append(
                CodonState(
                    codon=text,
                    edited=tuple(index + 1 for index in chosen),
                    amino_acid=STANDARD_CODE[text],
                )
            )
    return tuple(states)


def column_counts(residues: Iterable[str]) -> tuple[dict[str, int], int]:
    """How many reference residues of each amino acid one column holds, and how many in all.

    Gaps are excluded from both — the denominator is the number of reference sequences with a
    residue in this column, which is what makes the archived scores come in fifths in one column
    and quarters in the next.
    """
    counts: dict[str, int] = {}
    for residue in residues:
        if residue == "-":
            continue
        counts[residue] = counts.get(residue, 0) + 1
    return counts, sum(counts.values())


def call_codon(codon: str, residues: Iterable[str]) -> CodonCall | None:
    """Choose a state for one codon against one column, or ``None`` when nothing can be scored.

    ``None`` means the column holds no reference residue at all (every reference sequence has a
    gap there, or the query residue sits in a column it inserted): a score of 0/0 is not a score.
    """
    counts, compared = column_counts(residues)
    if compared == 0:
        return None
    states = codon_states(codon)

    def key(state: CodonState) -> tuple[Fraction, int, int, tuple[int, ...]]:
        score = Fraction(counts.get(state.amino_acid, 0), compared)
        return (-score, len(state.edited), 0 if 2 in state.edited else 1, state.edited)

    chosen = min(states, key=key)
    unedited = states[0]
    return CodonCall(
        state=chosen,
        matching=counts.get(chosen.amino_acid, 0),
        compared=compared,
        unedited=unedited,
        unedited_matching=counts.get(unedited.amino_acid, 0),
    )


def parse_cutoff(value: object) -> Fraction:
    """The cutoff C as an exact fraction in [0, 1], refused by name otherwise.

    Parsed through its decimal text: ``Fraction(0.8)`` is 3602879701896397/4503599627370496,
    larger than 4/5, and would silently drop every site that scores exactly 0.80.
    """
    text = str(value).strip()
    try:
        decimal = Decimal(text)
    except InvalidOperation:
        raise ValueError(f"cutoff {value!r} is not a number between 0 and 1") from None
    if not decimal.is_finite() or decimal < 0 or decimal > 1:
        raise ValueError(f"cutoff {value!r} is not a number between 0 and 1")
    return Fraction(decimal)


# --- the query: one CDS feature, extracted, mapped and translated -------------------------------


def transcript_positions(genome: Genome, gene: Gene) -> tuple[tuple[int, int], ...]:
    """Every transcript base as ``(genome position, strand)``, 0-based, in transcript order.

    The same reading as :func:`~plastid.core.model.coding_sequence` — every part on its own
    strand, in annotation order — so position ``i`` here is the genome base that became
    transcript base ``i`` there. A minus-strand part is read from its end; a part past the end of
    a circular record wraps through the origin.
    """
    length = len(genome.sequence)
    positions: list[tuple[int, int]] = []
    for part in gene.parts:
        span = range(part.start, part.end)
        ordered = span if part.strand > 0 else reversed(span)
        positions.extend(((index % length), part.strand) for index in ordered)
    return tuple(positions)


def locus_segments(positions: Sequence[tuple[int, int]]) -> str:
    """The coding-base-to-genome mapping in compact form, 1-based on both sides.

    ``1-435:84458-84024(-);436-…`` reads "transcript bases 1–435 are genome bases 84458 down to
    84024 on the minus strand". Consecutive transcript bases that step by one genome base on one
    strand form one segment, so the text is as long as the feature has parts (plus one at an
    origin crossing) and states the whole mapping.
    """
    if not positions:
        return ""
    segments: list[str] = []
    start = 0
    for index in range(1, len(positions) + 1):
        if index < len(positions):
            previous, current = positions[index - 1], positions[index]
            if current[1] == previous[1] and current[0] - previous[0] == previous[1]:
                continue
        first, last = positions[start], positions[index - 1]
        sign = "+" if first[1] > 0 else "-"
        segments.append(f"{start + 1}-{index}:{first[0] + 1}-{last[0] + 1}({sign})")
        start = index
    return ";".join(segments)


def translate_codons(codons: Sequence[str]) -> str:
    """The standard-code translation of whole codons; an unreadable codon translates to ``X``."""
    return "".join(STANDARD_CODE.get(codon, "X") for codon in codons)


@dataclass(frozen=True, slots=True)
class PreparedLocus:
    """One CDS feature made ready to score — or the reason it cannot be.

    ``codons`` start at ``/codon_start`` and exclude the terminal stop; ``protein`` is their
    standard-code translation (an ambiguous codon is ``X``). ``positions`` maps every transcript
    base of the feature to the genome.
    """

    gene: Gene
    cds: str
    positions: tuple[tuple[int, int], ...]
    offset: int
    codons: tuple[str, ...]
    protein: str
    reason: str | None = None
    ambiguous: tuple[int, ...] = ()

    @property
    def evaluable(self) -> bool:
        return self.reason is None


def prepare_locus(genome: Genome, gene: Gene) -> PreparedLocus:
    """Extract, frame and translate one CDS feature, with a reason whenever it cannot be scored.

    The reasons are the specification's: a pseudogene, a translation table whose amino acids are
    not the standard code's, a length that is not whole codons, an internal stop. A codon with an
    ambiguity code does not stop the locus — it is listed (``ambiguous``, 1-based codon numbers)
    and left unscored.
    """
    cds = coding_sequence(genome, gene).upper()
    positions = transcript_positions(genome, gene)
    offset = gene.codon_start - 1
    body = cds[offset:]
    usable = len(body) - len(body) % 3
    codons = tuple(body[index : index + 3] for index in range(0, usable, 3))
    reason: str | None = None
    if gene.pseudo:
        reason = "pseudogene (/pseudo)"
    elif gene.transl_table is not None and gene.transl_table not in (1, 11):
        reason = f"translation table {gene.transl_table} (only the standard code is scored)"
    elif len(body) % 3 and not gene.partial:
        reason = f"frame: {len(body)} nt from /codon_start is not whole codons"
    if codons and STANDARD_CODE.get(codons[-1]) == "*":
        codons = codons[:-1]
    protein = translate_codons(codons)
    ambiguous = tuple(number for number, codon in enumerate(codons, start=1) if set(codon) - _BASES)
    if reason is None and "*" in protein:
        reason = f"internal stop at codon {protein.index('*') + 1}"
    if reason is None and not codons:
        reason = "no complete codon"
    return PreparedLocus(
        gene=gene,
        cds=cds,
        positions=positions,
        offset=offset,
        codons=codons,
        protein=protein,
        reason=reason,
        ambiguous=ambiguous,
    )


def supported_gene(normalized: NormalizedName) -> str | None:
    """The PREP-Cp gene a name stands for, by its canonical name or any synonym of its locus.

    The canonical direction follows the run's naming reference (vocabulary A), so ``pafI`` may be
    canonical in one project and ``ycf3`` in another; matching the whole locus keeps the 35-gene
    set independent of that choice.
    """
    for name in (normalized.canonical, *normalized.synonyms):
        if name in SUPPORTED_GENES:
            return name
    return None


# --- scoring one locus against one aligned reference ---------------------------------------------


@dataclass(frozen=True, slots=True)
class ScoredSite:
    """One predicted C-to-U site: a transcript base, its codon call, and its alignment column."""

    cds_position: int  # 1-based from the feature's first base
    codon_number: int  # 1-based from /codon_start
    codon_position: int  # 1, 2 or 3
    column: int  # 1-based column of the query-plus-reference alignment
    call: CodonCall


@dataclass(frozen=True, slots=True)
class LocusPrediction:
    """Everything one locus produced: the sites and the counts that say what was looked at."""

    sites: tuple[ScoredSite, ...]
    c_codons: int  # codons that contain a C (the ones the method asks about)
    scored: int  # of those, codons scored against at least one reference residue
    unscored: tuple[int, ...]  # C codons whose column holds no reference residue (codon numbers)
    ambiguous: tuple[int, ...]  # codons with an ambiguity code — not evaluated (codon numbers)
    calls: tuple[tuple[int, CodonCall], ...] = field(default=())  # every scored C codon


def predict_locus(
    locus: PreparedLocus,
    query_aligned: str,
    reference_aligned: Sequence[str],
    cutoff: Fraction,
) -> LocusPrediction:
    """Score every C-containing codon of ``locus`` against the aligned reference columns.

    ``query_aligned`` is the query's row of the query-plus-reference alignment and must spell the
    locus's protein once gaps are removed; a mismatch means the alignment is not of this protein,
    and is refused rather than scored.
    """
    residues = query_aligned.replace("-", "")
    if residues != locus.protein:
        raise ValueError(
            f"{locus.gene.name}: the aligned query does not spell the translated CDS "
            f"({len(residues)} residues against {len(locus.protein)})"
        )
    width = len(query_aligned)
    for row in reference_aligned:
        if len(row) != width:
            raise ValueError(f"{locus.gene.name}: alignment rows have different lengths")
    columns = [index for index, symbol in enumerate(query_aligned) if symbol != "-"]

    sites: list[ScoredSite] = []
    calls: list[tuple[int, CodonCall]] = []
    unscored: list[int] = []
    c_codons = 0
    for number, (codon, column) in enumerate(zip(locus.codons, columns, strict=True), start=1):
        if set(codon) - _BASES:
            continue
        if "C" not in codon:
            continue
        c_codons += 1
        call = call_codon(codon, (row[column] for row in reference_aligned))
        if call is None:
            unscored.append(number)
            continue
        calls.append((number, call))
        if not call.edited or call.score < cutoff:
            continue
        for codon_position in call.state.edited:
            sites.append(
                ScoredSite(
                    cds_position=locus.offset + 3 * (number - 1) + codon_position,
                    codon_number=number,
                    codon_position=codon_position,
                    column=column + 1,
                    call=call,
                )
            )
    return LocusPrediction(
        sites=tuple(sites),
        c_codons=c_codons,
        scored=len(calls),
        unscored=tuple(unscored),
        ambiguous=locus.ambiguous,
        calls=tuple(calls),
    )


# --- one genome: every CDS, one status each ------------------------------------------------------

# Adds one protein to a fixed alignment: ``add(reference_rows, (name, protein))`` returns every row
# by name, the query's included. The module passes MAFFT ``--add``; a test passes its own.
AddAligner = Callable[[Sequence[tuple[str, str]], tuple[str, str]], Mapping[str, str]]

QUERY_ROW = "query"


@dataclass(frozen=True, slots=True)
class LocusOutcome:
    """What happened to one CDS feature — or to one supported gene the record does not carry.

    ``status`` is one of the four the specification names; ``reason`` says why a locus was not
    evaluated. ``alignment`` holds the query-plus-reference rows of an evaluated locus, query
    first, so any score can be recomputed by hand.
    """

    locus: str
    name_in_record: str
    gene_a: str
    prep_gene: str | None
    status: str
    reason: str = ""
    gene: Gene | None = None
    prepared: PreparedLocus | None = None
    prediction: LocusPrediction | None = None
    alignment: tuple[tuple[str, str], ...] = ()


def reference_row_name(row: ReferenceRow) -> str:
    """A FASTA row name an aligner keeps whole: no whitespace."""
    return f"{row.taxon.replace(' ', '_')}|{row.accession}"


def evaluate_genome(
    genome: Genome,
    reference: EditingReference,
    cutoff: Fraction,
    aliases: AliasTable,
    add: AddAligner,
) -> tuple[LocusOutcome, ...]:
    """Every CDS of ``genome`` with one status, then every supported gene the record lacks.

    A locus is one CDS feature: the two IR copies of a gene are two loci (``ndhB``, ``ndhB#2``),
    each scored and mapped on its own coordinates, and so are two features that merely share a
    name. Names are compared in vocabulary A (the alias table), so a record that writes ``pafI``
    still reaches the ``ycf3`` reference.
    """
    names = normalize_genome(genome, aliases)
    cds = [gene for gene in genome.genes if gene.kind == "CDS"]
    totals: dict[str, int] = {}
    for gene in cds:
        totals[names.canonical(gene)] = totals.get(names.canonical(gene), 0) + 1
    seen: dict[str, int] = {}
    outcomes: list[LocusOutcome] = []
    present: set[str] = set()
    for gene in cds:
        name_a = names.canonical(gene)
        seen[name_a] = seen.get(name_a, 0) + 1
        locus_id = name_a if totals[name_a] == 1 else f"{name_a}#{seen[name_a]}"
        prep = supported_gene(names.normalized(gene))
        common = {
            "locus": locus_id,
            "name_in_record": gene.name,
            "gene_a": name_a,
            "prep_gene": prep,
            "gene": gene,
        }
        if prep is None:
            outcomes.append(LocusOutcome(**common, status=STATUS_NOT_SUPPORTED))
            continue
        present.add(prep)
        prepared = prepare_locus(genome, gene)
        if not prepared.evaluable:
            outcomes.append(
                LocusOutcome(
                    **common,
                    status=STATUS_NOT_EVALUATED,
                    reason=str(prepared.reason),
                    prepared=prepared,
                )
            )
            continue
        outcomes.append(_score_locus(common, prepared, reference.rows(prep), cutoff, add))
    for prep in SUPPORTED_GENES:
        if prep not in present:
            outcomes.append(
                LocusOutcome(
                    locus=prep,
                    name_in_record="",
                    gene_a="",
                    prep_gene=prep,
                    status=STATUS_ABSENT,
                    reason="no CDS of this gene in the record",
                )
            )
    return tuple(outcomes)


def _score_locus(
    common: Mapping[str, Any],
    prepared: PreparedLocus,
    rows: Sequence[ReferenceRow],
    cutoff: Fraction,
    add: AddAligner,
) -> LocusOutcome:
    """Add one locus to its gene's fixed alignment and score it — or say why it could not be.

    The reference rows must come back unchanged once the columns the query inserted are removed:
    that is what "align to the reference without re-aligning it" means, and it is checked on every
    call rather than assumed of the aligner.
    """
    if not rows:
        return LocusOutcome(
            **common,
            status=STATUS_NOT_EVALUATED,
            reason=f"reference B′ holds no {common['prep_gene']} protein",
            prepared=prepared,
        )
    names = [reference_row_name(row) for row in rows]
    aligned = add(
        list(zip(names, (row.aligned for row in rows), strict=True)),
        (QUERY_ROW, prepared.protein),
    )
    reference_rows = [aligned[name] for name in names]
    query = aligned[QUERY_ROW]
    kept = [
        index for index in range(len(query)) if any(row[index] != "-" for row in reference_rows)
    ]
    for name, row, stored in zip(names, reference_rows, rows, strict=True):
        if "".join(row[index] for index in kept) != stored.aligned:
            return LocusOutcome(
                **common,
                status=STATUS_NOT_EVALUATED,
                reason=f"alignment failure: the aligner changed the reference row {name}",
                prepared=prepared,
            )
    try:
        prediction = predict_locus(prepared, query, reference_rows, cutoff)
    except ValueError as error:
        return LocusOutcome(
            **common,
            status=STATUS_NOT_EVALUATED,
            reason=f"alignment failure: {error}",
            prepared=prepared,
        )
    return LocusOutcome(
        **common,
        status=STATUS_EVALUATED,
        prepared=prepared,
        prediction=prediction,
        alignment=(
            (f"{QUERY_ROW}|{common['locus']}", query),
            *zip(names, reference_rows, strict=True),
        ),
    )


# --- reference B′ --------------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class ReferenceRow:
    """One taxon's mature protein for one gene, as it sits in the fixed alignment."""

    taxon: str
    accession: str
    aligned: str


@dataclass(frozen=True, slots=True)
class EditingReference:
    """Reference B′ as the scorer reads it: one fixed alignment per supported gene, one digest."""

    label: str
    digest: str
    genes: Mapping[str, tuple[ReferenceRow, ...]]
    payload: Mapping[str, Any]

    def rows(self, gene: str) -> tuple[ReferenceRow, ...]:
        return self.genes.get(gene, ())


def _canonical_json(payload: Mapping[str, Any]) -> bytes:
    return json.dumps(payload, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode(
        "utf-8"
    )


def reference_digest(payload: Mapping[str, Any]) -> str:
    """The digest of everything B′ is made of — records, edit lists, recipe, aligner, alignment.

    It is the SHA-256 of the canonical JSON of the payload without its own ``digest`` field, so a
    change to any record pin, edit site, source, recipe entry, aligner version or aligned row
    changes it.
    """
    body = {key: value for key, value in payload.items() if key != "digest"}
    return hashlib.sha256(_canonical_json(body)).hexdigest()


def load_reference(path: str | Path | None = None) -> EditingReference:
    """Reference B′ as packaged (``data/rna_editing/reference_b_prime.json``) or from ``path``."""
    source = Path(path) if path is not None else data_file(REFERENCE_FILE)
    try:
        payload = json.loads(source.read_text(encoding="utf-8"))
    except FileNotFoundError:
        raise EditingReferenceError(
            f"{source}: reference B′ is missing — `plastid editing-reference build` makes it"
        ) from None
    return reference_from_payload(payload)


def reference_set_entry(path: str | Path | None = None) -> dict[str, Any]:
    """What a run manifest records about the reference a prediction was scored against."""
    source = Path(path) if path is not None else data_file(REFERENCE_FILE)
    reference = load_reference(source)
    payload = reference.payload
    return {
        "role": "rna_editing_reference",
        "label": reference.label,
        "digest": reference.digest,
        "path": f"plastid/data/{REFERENCE_FILE}" if path is None else str(source),
        "sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "taxa": [
            {"accession": taxon["accession"], "sha256": taxon["sha256"]}
            for taxon in payload["taxa"]
        ],
        "edit_sites": len(payload["edit_sites"]),
        "genes": len(payload["genes"]),
    }


def reference_from_payload(payload: Mapping[str, Any]) -> EditingReference:
    """Read a B′ payload, refusing one whose digest does not describe its own contents."""
    if payload.get("schema") != REFERENCE_SCHEMA:
        raise EditingReferenceError(
            f"reference schema {payload.get('schema')!r} is not {REFERENCE_SCHEMA!r}"
        )
    stated = payload.get("digest")
    actual = reference_digest(payload)
    if stated != actual:
        raise EditingReferenceError(
            f"reference digest {str(stated)[:12]}… does not match its contents ({actual[:12]}…) "
            f"— the file was edited after it was built; rebuild it instead"
        )
    genes: dict[str, tuple[ReferenceRow, ...]] = {}
    for gene, entry in payload["genes"].items():
        rows = tuple(
            ReferenceRow(taxon=row["taxon"], accession=row["accession"], aligned=row["aligned"])
            for row in entry["rows"]
        )
        widths = {len(row.aligned) for row in rows}
        if len(widths) > 1:
            raise EditingReferenceError(f"reference {gene}: aligned rows differ in length")
        genes[gene] = rows
    unknown = sorted(set(genes) - set(SUPPORTED_GENES))
    if unknown:
        raise EditingReferenceError(f"reference carries genes outside the 35: {unknown}")
    return EditingReference(
        label=str(payload["label"]), digest=str(stated), genes=genes, payload=payload
    )


@dataclass(frozen=True, slots=True)
class EditSite:
    """One sourced edit site of one reference taxon, in that record's transcript coordinates."""

    accession: str
    gene: str
    cds_position: int  # 1-based from the feature's first base
    codon_position: int
    aa_position: int
    codon_before: str
    codon_after: str
    aa_before: str
    aa_after: str
    source: str
    locator: str
    source_coordinate: str


_SITE_FIELDS = tuple(EditSite.__dataclass_fields__)


def read_edit_sites(path: str | Path) -> tuple[EditSite, ...]:
    """The sourced edit list (tab-separated, one site per row, a header naming the fields)."""
    with open(path, encoding="utf-8", newline="") as handle:
        rows = [row for row in csv.reader(handle, delimiter="\t") if row and row[0][:1] != "#"]
    if not rows or tuple(rows[0]) != _SITE_FIELDS:
        raise EditingReferenceError(f"{path}: the header must be {' '.join(_SITE_FIELDS)}")
    sites = []
    for line, row in enumerate(rows[1:], start=2):
        if len(row) != len(_SITE_FIELDS):
            raise EditingReferenceError(
                f"{path}:{line}: {len(row)} fields, not {len(_SITE_FIELDS)}"
            )
        values = dict(zip(_SITE_FIELDS, row, strict=True))
        try:
            sites.append(
                EditSite(
                    **{
                        key: int(value)
                        if key in ("cds_position", "codon_position", "aa_position")
                        else value
                        for key, value in values.items()
                    }
                )
            )
        except ValueError as error:
            raise EditingReferenceError(f"{path}:{line}: {error}") from None
    return tuple(sites)


Aligner = Callable[[Sequence[tuple[str, str]]], Mapping[str, str]]


def _loci_by_gene(genome: Genome, aliases: AliasTable) -> dict[str, list[Gene]]:
    names = normalize_genome(genome, aliases)
    loci: dict[str, list[Gene]] = {}
    for gene in genome.genes:
        if gene.kind != "CDS":
            continue
        prep = supported_gene(names.normalized(gene))
        if prep is not None:
            loci.setdefault(prep, []).append(gene)
    return loci


def check_site(genome: Genome, gene: Gene, site: EditSite) -> str | None:
    """Why a sourced site does not fit the pinned record, or ``None`` when it fits exactly."""
    cds = coding_sequence(genome, gene).upper()
    if not 1 <= site.cds_position <= len(cds):
        return f"position {site.cds_position} is outside the {len(cds)}-nt CDS"
    if cds[site.cds_position - 1] != "C":
        return f"position {site.cds_position} is {cds[site.cds_position - 1]}, not C"
    offset = gene.codon_start - 1
    number, within = divmod(site.cds_position - 1 - offset, 3)
    if number + 1 != site.aa_position or within + 1 != site.codon_position:
        return (
            f"position {site.cds_position} is codon {number + 1} position {within + 1}, the "
            f"source says codon {site.aa_position} position {site.codon_position}"
        )
    codon = cds[offset + 3 * number : offset + 3 * number + 3]
    if codon != site.codon_before.upper().replace("U", "T"):
        return f"codon {number + 1} is {codon}, the source says {site.codon_before}"
    if STANDARD_CODE.get(codon) != site.aa_before:
        return f"codon {codon} reads {STANDARD_CODE.get(codon)}, the source says {site.aa_before}"
    return None


def build_reference(
    *,
    records: Mapping[str, Genome],
    record_sha256: Mapping[str, str],
    taxa: Sequence[Mapping[str, Any]],
    sites: Sequence[EditSite],
    sources: Mapping[str, Any],
    aliases: AliasTable,
    aligner: Aligner,
    aligner_record: Mapping[str, Any],
    label: str,
) -> dict[str, Any]:
    """Rebuild B′: mature proteins of the pinned taxa, aligned once per gene, under one digest.

    For every supported gene and every taxon that carries it as a CDS: extract the CDS (the first
    copy when the record has two — the IR copies of one gene are one protein), apply that taxon's
    sourced C-to-U edits, translate with the standard code, drop the terminal stop. A sourced
    site that does not fit its record (not a C, a different codon, outside the CDS) stops the
    build: an edit list that disagrees with its own record is not a reference. A taxon whose
    mature protein still has an internal stop is excluded for that gene, with the reason kept.
    """
    taxon_by_accession = {str(taxon["accession"]): taxon for taxon in taxa}
    for accession, taxon in taxon_by_accession.items():
        if accession not in records:
            raise EditingReferenceError(f"{accession}: no record supplied for {taxon['organism']}")
        pinned = str(taxon["sha256"])
        if record_sha256.get(accession) != pinned:
            raise EditingReferenceError(
                f"{accession}: the record's SHA-256 {str(record_sha256.get(accession))[:12]}… is "
                f"not the pinned {pinned[:12]}…"
            )
    by_taxon_gene: dict[tuple[str, str], list[EditSite]] = {}
    for site in sites:
        if site.accession not in taxon_by_accession:
            raise EditingReferenceError(f"edit site on {site.accession}, which is not a B′ taxon")
        if site.source not in sources:
            raise EditingReferenceError(f"edit site cites unknown source {site.source!r}")
        if site.gene not in SUPPORTED_GENES:
            raise EditingReferenceError(f"edit site on {site.gene}, which is not one of the 35")
        by_taxon_gene.setdefault((site.accession, site.gene), []).append(site)

    records = {
        accession: apply_taxon_names(
            genome, taxon_by_accession[accession].get("names") or (), accession
        )
        for accession, genome in records.items()
        if accession in taxon_by_accession
    }
    loci = {accession: _loci_by_gene(records[accession], aliases) for accession in records}
    genes: dict[str, Any] = {}
    used_sites: list[dict[str, Any]] = []
    for gene_name in SUPPORTED_GENES:
        proteins: list[tuple[str, dict[str, Any]]] = []
        excluded: list[dict[str, str]] = []
        for taxon in taxa:
            accession = str(taxon["accession"])
            genome = records[accession]
            candidates = [gene for gene in loci[accession].get(gene_name, []) if not gene.pseudo]
            if not candidates:
                if (accession, gene_name) in by_taxon_gene:
                    raise EditingReferenceError(
                        f"{accession} {gene_name}: edit sites are listed but the record has no "
                        f"CDS of that gene"
                    )
                continue
            chosen = candidates[0]
            cds = list(coding_sequence(genome, chosen).upper())
            edits = sorted(
                by_taxon_gene.get((accession, gene_name), []), key=lambda site: site.cds_position
            )
            for site in edits:
                problem = check_site(genome, chosen, site)
                if problem is not None:
                    raise EditingReferenceError(
                        f"{accession} {gene_name} ({site.source}, {site.locator}): {problem}"
                    )
                cds[site.cds_position - 1] = "T"
            mature = "".join(cds)
            offset = chosen.codon_start - 1
            for site in edits:
                start = offset + 3 * (site.aa_position - 1)
                edited_codon = mature[start : start + 3]
                if edited_codon != site.codon_after.upper().replace("U", "T"):
                    raise EditingReferenceError(
                        f"{accession} {gene_name} ({site.source}, {site.locator}): the edited "
                        f"codon {site.aa_position} reads {edited_codon}, the source says "
                        f"{site.codon_after}"
                    )
                if STANDARD_CODE.get(edited_codon) != site.aa_after:
                    raise EditingReferenceError(
                        f"{accession} {gene_name} ({site.source}, {site.locator}): {edited_codon} "
                        f"reads {STANDARD_CODE.get(edited_codon)}, the source says {site.aa_after}"
                    )
            body = mature[offset:]
            codons = [body[index : index + 3] for index in range(0, len(body) - len(body) % 3, 3)]
            if codons and STANDARD_CODE.get(codons[-1]) == "*":
                codons = codons[:-1]
            protein = translate_codons(codons)
            truncated: int | None = None
            if "*" in protein:
                stop = protein.index("*") + 1
                if stop not in {site.aa_position for site in edits if site.aa_after == "*"}:
                    excluded.append(
                        {"accession": accession, "reason": f"internal stop at codon {stop}"}
                    )
                    continue
                # A sourced edit that creates a stop ends the mature protein there (black pine
                # petD: the edited CAA makes the 160-residue protein the other taxa have).
                protein, truncated = protein[: stop - 1], stop
            genomic = coding_sequence(genome, chosen).upper()
            proteins.append(
                (
                    str(taxon["key"]),
                    {
                        "taxon": str(taxon["organism"]),
                        "accession": accession,
                        "copies": len(candidates),
                        "location": location_text(chosen),
                        "cds_sha256": hashlib.sha256(genomic.encode("ascii")).hexdigest(),
                        "edits": [site.cds_position for site in edits],
                        "stop_created_at": truncated,
                        "protein": protein,
                    },
                )
            )
            used_sites.extend(_site_payload(site) for site in edits)
        if not proteins:
            continue
        if len(proteins) == 1:
            aligned = {proteins[0][0]: proteins[0][1]["protein"]}
        else:
            aligned = dict(aligner([(name, row["protein"]) for name, row in proteins]))
        rows = []
        for name, row in proteins:
            text = aligned[name]
            if text.replace("-", "") != row["protein"]:
                raise EditingReferenceError(f"{gene_name}: the aligner changed {name}'s protein")
            rows.append({**row, "aligned": text})
        genes[gene_name] = {"rows": rows, "excluded": excluded}

    payload: dict[str, Any] = {
        "schema": REFERENCE_SCHEMA,
        "label": label,
        "method": METHOD_LABEL,
        "recipe": {
            "translation": "NCBI table 1 (standard); terminal stop removed; every codon read as "
            "its own amino acid",
            "copies": "one protein per taxon per gene (the first CDS copy in record order)",
            "edits": "sourced C-to-U sites applied to the CDS before translation",
            "aligner": dict(aligner_record),
        },
        "taxa": [dict(taxon) for taxon in taxa],
        "sources": dict(sources),
        "edit_sites": used_sites,
        "genes": genes,
    }
    payload["digest"] = reference_digest(payload)
    return payload


def apply_taxon_names(genome: Genome, names: Sequence[Mapping[str, Any]], accession: str) -> Genome:
    """The record with a gene symbol given to each product-only CDS its taxon entry names.

    A feature is matched by location *and* /product, so a name cannot land on a feature the entry
    did not describe; an entry that matches nothing stops the build rather than being skipped. An
    entry may also state a corrected CDS model (``cds``) where the record's own model contradicts
    itself, with the evidence beside it; the pinned file is never changed.
    """
    if not names:
        return genome
    genes = list(genome.genes)
    for entry in names:
        matches = [
            index
            for index, gene in enumerate(genes)
            if gene.kind == "CDS"
            and location_text(gene) == entry["location"]
            and gene.product == entry["product"]
        ]
        if len(matches) != 1:
            raise EditingReferenceError(
                f"{accession}: the name entry {entry['gene']} at {entry['location']} matches "
                f"{len(matches)} CDS features, not one"
            )
        renamed = replace(genes[matches[0]], name=str(entry["gene"]))
        if entry.get("cds"):
            renamed = replace(renamed, parts=parse_location_text(str(entry["cds"])))
        genes[matches[0]] = renamed
    return replace(genome, genes=tuple(genes))


def parse_location_text(text: str) -> tuple[Interval, ...]:
    """The inverse of :func:`location_text`: ``"33899..34000(-)"`` → parts, 0-based half-open."""
    parts = []
    for item in text.split(","):
        span, _, strand = item.strip().partition("(")
        start, _, end = span.partition("..")
        if strand not in ("+)", "-)") or not start.isdigit() or not end.isdigit():
            raise EditingReferenceError(f"location {text!r} is not start..end(+|-)[,…]")
        parts.append(Interval.from_1based(int(start), int(end), 1 if strand == "+)" else -1))
    return tuple(parts)


def load_build_inputs(directory: str | Path) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """The pinned taxa and the edit-list sources (``taxa.yaml``, ``sources.yaml``)."""
    root = Path(directory)
    taxa = yaml.safe_load((root / "taxa.yaml").read_text(encoding="utf-8"))["taxa"]
    sources = yaml.safe_load((root / "sources.yaml").read_text(encoding="utf-8"))["sources"]
    if not isinstance(taxa, list) or not isinstance(sources, dict):
        raise EditingReferenceError(f"{root}: taxa.yaml or sources.yaml is malformed")
    return taxa, sources


def _site_payload(site: EditSite) -> dict[str, Any]:
    return {name: getattr(site, name) for name in _SITE_FIELDS}


def location_text(gene: Gene) -> str:
    """A feature's parts in annotation order, 1-based inclusive, with each part's strand."""
    return ",".join(
        f"{part.start + 1}..{part.end}({'+' if part.strand > 0 else '-'})" for part in gene.parts
    )
