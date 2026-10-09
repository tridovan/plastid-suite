"""Turning BLAST hits against a curated reference set into gene calls — shared by both organelles.

The rule that shapes everything here is decision **D3A-02**: a CDS's boundaries come from the
alignment to a curated reference protein, never from ORF-walking. The mito answer-key record is
why — its own 31 CDS score 21 ok · 4 pseudogene · 3 no_stop · 3 no_start, so a "find ATG … stop"
caller cannot reproduce it, and one that "repairs" what it finds would publish a record NCBI did
not.

**Why this module exists separately.** P3A-03 wrote these rules for the plant mitochondrion
(``core/mito_caller.py``) and P3B-07 needed the same machinery for the plastid path — the second
real copy, which is where this project's own rule says to extract rather than copy (CLAUDE.md →
DRY after the second real copy; the counterexample it names is CGAS's 14 copy-pasted glue blocks).
Nothing about the rules below was mitochondrion-specific: the floors and the translation table
are the two things each path must state for itself, and both arrive as **parameters with no
default** — a path that inherited another path's measured floor would be quoting a number it
never measured.

What the measurements on the mito answer-key record forced (each is a test in
``tests/unit/test_mito_caller.py``, which the extraction left unchanged):

* **BLAST truncates ``sseqid`` at the first whitespace.** The reference headers carry the organism
  name, so a hit identifies its reference by the **first three** ``|`` fields — accession, the
  reference's own gene symbol, and our vocabulary symbol (``-`` when we have none).
* **One header can stand for two sequences.** The committed mito reference FASTA holds 8 duplicate
  headers; six are byte-identical repeat copies, but Arabidopsis ``atp6`` is 385 *and* 349 aa and
  Psychotria ``cox2`` is 257 *and* 390 aa. An HSP therefore belongs to the **smallest reference
  length that contains its subject range**, and HSPs assigned to different sequences under one key
  are never chained together — that is how a chimera would be built.
* **Real HSP sets contain noise.** A 27 aa 100 %-identity fragment of ``atp8`` sits at 8,557; two
  frames align over one query region (``sdh4``, ``mttB``). An HSP is kept only when it adds a
  reference residue no kept HSP already covers, and only when it does not contradict a kept HSP's
  frame over the same query region. Everything dropped comes back with its reason.
* **The reference chooses the call.** Measured: ranking candidates by bitscore let a 75.6 %
  identity alignment 418 aa long drag ``atp6`` 381 bp upstream of the published CDS, while a
  93.4 % identity alignment sat on it. The winner is therefore the highest-**identity** candidate
  that passes the floors (bitscore breaks ties) — the closest relative, not the longest alignment.
* **The stop codon is taken only when the sequence confirms it**, one codon past the alignment's
  end. ccmFc and nad6 are annotated in the record *without* their stop, which sits immediately
  after; cox2's alignment stops six codons short of the gene's real end. Reading one codon is a
  check; reading until a stop appears is ORF-walking.

Scope, recorded because a paper has to state it (D3A-05): this is **annotation and submission
only** — no de-novo gene discovery, no intron model, no RNA-editing prediction. Multi-exon calls
report the exons the alignment supports (D3A-04) and never claim a join they cannot show. A gene
absent from the reference set is reported as not called, never guessed.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass, replace
from io import StringIO
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from plastid.core.exon_structure import Alignment, Junction, Placement, Transfer, place
from plastid.core.model import Gene, Genome, Interval
from plastid.core.orf import START_CODONS, STOP_CODONS
from plastid.core.sequence import reverse_complement

__all__ = [
    "CallResult",
    "CompetingReference",
    "DroppedHit",
    "GeneCall",
    "Hsp",
    "NotCalled",
    "ReferenceEntry",
    "ReferenceIndex",
    "call_cds",
    "call_rrna",
    "header_field",
    "read_reference",
    "render_fasta",
    "transfer_structure",
]


@dataclass(frozen=True, slots=True)
class ReferenceEntry:
    """One reference gene: where it came from, which locus it is, and its curated sequence.

    ``vocabulary`` is the symbol this path's vocabulary knows for the locus, or ``None`` when the
    reference carries a gene it has no name for. Those are **kept** rather than dropped, so an
    annotator run can say it found ``rps12`` instead of being silent about a gene it cannot name —
    measured 2026-10-03, the standard plant-mito set is wider than the 31 loci that vocabulary
    held. ``sequence`` is a protein for a CDS entry and a nucleotide for an rRNA one; the caller
    is told which, and never converts between them.
    """

    accession: str
    organism: str
    source_gene: str  # the symbol the reference record itself used
    vocabulary: str | None  # our canonical symbol, when we have one
    sequence: str
    # The record's own exon structure, as junctions a target can be searched with (P3B-10). Empty
    # for a single-exon gene, for an rRNA, and on the mito path, which does not use this yet.
    junctions: tuple[Junction, ...] = ()

    @property
    def key(self) -> str:
        """What a BLAST hit calls this sequence: the first three fields of its header.

        Measured 2026-10-03: BLAST+ truncates ``sseqid`` at the first whitespace, and the fourth
        header field is the organism name — ``Nicotiana tabacum`` arrives as ``Nicotiana``. Three
        fields are also all a hit needs; the fourth is only there to be read by a person. ``-``
        means this toolkit has no symbol for the gene, exactly as the FASTA writer writes it.

        **Whitespace inside a field is removed**, because a name is allowed to carry one and a
        header is not: measured 08/10/2026 on the plastid reference set, ``NC_042827.1`` names its
        rRNA ``23S ribosomal RNA`` (a ``/product`` string, no ``/gene``), and the raw header came
        back from BLAST as ``NC_042827.1|23S`` — the parser refused it by name rather than
        matching a key that could never match. The mito path never met this because no
        mitochondrial gene symbol contains a space.
        """
        fields = (self.accession, self.source_gene, self.vocabulary or "-")
        return "|".join(header_field(part) for part in fields)


@dataclass(frozen=True, slots=True)
class Hsp:
    """One BLAST high-scoring pair, in BLAST's own convention: **1-based inclusive** on both
    sequences, with ``qstart > qend`` when the query aligned on its minus strand.

    ``frame`` is blastx's frame (1..3 / -1..-3) and is ``None`` for a nucleotide search, where the
    strand comes from the subject's direction instead — see :attr:`strand`.
    """

    subject_key: str
    pident: float
    length: int  # aligned residues: amino acids for blastx, bases for blastn
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: float
    bitscore: float
    frame: int | None

    @property
    def strand(self) -> int:
        """The strand the *gene* lies on.

        For a protein search the query's frame says it. For a nucleotide search the subject's
        direction says it, and the inverse of the obvious reading is the correct one: the
        reference is stored 5'→3', so ``sstart > send`` means the query's plus strand matched the
        reference's *reverse complement* — the gene reads leftward.
        """
        if self.frame:
            return 1 if self.frame > 0 else -1
        return 1 if self.send > self.sstart else -1

    @property
    def qlo(self) -> int:
        return min(self.qstart, self.qend)

    @property
    def qhi(self) -> int:
        return max(self.qstart, self.qend)

    @property
    def slo(self) -> int:
        return min(self.sstart, self.send)

    @property
    def shi(self) -> int:
        return max(self.sstart, self.send)

    @property
    def span_label(self) -> str:
        """``q 1,234-1,300`` — how a hit is named in a report."""
        return f"q {self.qlo}-{self.qhi}"

    def overlaps(self, other: Hsp) -> bool:
        return self.qlo <= other.qhi and other.qlo <= self.qhi


def header_field(value: str) -> str:
    """One ``|``-delimited header field, whitespace-free — see :attr:`ReferenceEntry.key`.

    Public because a junction's id is built the same way (P3B-10): a header written with one rule
    and read with another is the bug this function exists to prevent, and a second copy of the
    rule is how it comes back.
    """
    return "_".join(value.split())


def _subtract(start: int, end: int, other_start: int, other_end: int) -> list[tuple[int, int]]:
    """``[start, end]`` minus ``[other_start, other_end]`` — possibly in two pieces."""
    if other_end < start or end < other_start:
        return [(start, end)]
    pieces = []
    if start < other_start:
        pieces.append((start, other_start - 1))
    if other_end < end:
        pieces.append((other_end + 1, end))
    return pieces


@dataclass(frozen=True, slots=True)
class ReferenceIndex:
    """The reference sequences a run searches, grouped the way a BLAST subject key identifies them.

    ``vocabulary`` is what this index may name: the path's own symbols for the kind being called.
    A reference whose gene has no symbol can only ever be *reported* (:meth:`unnameable`), never
    named into an annotation — the H-06 gap, kept visible.
    """

    entries: tuple[ReferenceEntry, ...]
    vocabulary: tuple[str, ...]

    @classmethod
    def from_entries(
        cls, entries: Sequence[ReferenceEntry], *, vocabulary: Sequence[str]
    ) -> ReferenceIndex:
        ordered = tuple(sorted(entries, key=lambda entry: (entry.key, len(entry.sequence))))
        # one symbol once, in the caller's order: a reference set holds several sequences per
        # locus, and a vocabulary that repeated them would call — and report — each locus as many
        # times as it has reference copies
        unique = tuple(dict.fromkeys(vocabulary))
        return cls(entries=ordered, vocabulary=unique)

    def for_symbol(self, symbol: str) -> tuple[ReferenceEntry, ...]:
        return tuple(entry for entry in self.entries if entry.vocabulary == symbol)

    def unnameable(self) -> tuple[ReferenceEntry, ...]:
        return tuple(entry for entry in self.entries if entry.vocabulary is None)

    def entry_for(self, key: str) -> ReferenceEntry | None:
        return next((entry for entry in self.entries if entry.key == key), None)

    def lengths_by_key(self) -> dict[str, list[int]]:
        lengths: dict[str, list[int]] = defaultdict(list)
        for entry in self.entries:
            lengths[entry.key].append(len(entry.sequence))
        return {key: sorted(set(values)) for key, values in lengths.items()}

    def entry_at(self, key: str, length: int) -> ReferenceEntry | None:
        return next(
            (entry for entry in self.entries if entry.key == key and len(entry.sequence) == length),
            None,
        )


# The two notes a placement from the reference's own exon structure supersedes: when the structure
# establishes an end, the alignment's failure to reach it is no longer a fact about the call
# (`exon_structure.place`), and a table that kept saying "partial" would contradict its own
# coordinates. Named here because the note and the test that drops it must be the same string.
NOT_REACHED_FIVE_PRIME = "the alignment does not reach the reference 5' end"
NOT_REACHED_THREE_PRIME = "the alignment does not reach the reference 3' end"


@dataclass(frozen=True, slots=True)
class DroppedHit:
    """A hit the caller did not use, and why. Nothing is discarded in silence."""

    symbol: str
    hsp: Hsp
    reason: str


@dataclass(frozen=True, slots=True)
class NotCalled:
    """A vocabulary locus with no call, and the best evidence that was rejected."""

    symbol: str
    reason: str
    best_accession: str | None = None
    best_identity: float | None = None
    best_coverage: float | None = None
    best_bitscore: float | None = None


@dataclass(frozen=True, slots=True)
class CompetingReference:
    """One reference sequence considered for a locus, with its own numbers.

    A call is a choice, and a reader has to see what it was chosen over — especially where two
    sequences share one subject key and differ only in length.
    """

    symbol: str
    reference_accession: str
    reference_gene: str
    reference_length: int
    identity: float
    coverage: float
    bitscore: float
    passes_floors: bool
    chosen: bool


@dataclass(frozen=True, slots=True)
class GeneCall:
    """One gene, called. Coordinates are GenBank-style 1-based inclusive inside ``parts``.

    ``parts`` is in **transcript order**, which is the order the alignment found: a trans-spliced
    gene's exons do not run in coordinate order, and the reference protein is the only witness to
    which comes first. Each part carries its own strand — measured, nad1 and nad5 align on both
    strands of one mitogenome.
    """

    symbol: str
    kind: str  # "CDS" | "rRNA"
    parts: tuple[Interval, ...]
    coverage: float
    identity: float
    bitscore: float
    reference: ReferenceEntry
    # the table the initiator test ran with, kept on the call so ``to_gene`` cannot disagree
    # with the refinement that produced it
    transl_table: int | None = None
    notes: tuple[str, ...] = ()

    @property
    def partial_start(self) -> bool:
        """The 5' coordinate is not established — see the module docstring and D3A-04."""
        return self.parts[0].partial_start

    @property
    def partial_end(self) -> bool:
        return self.parts[-1].partial_end

    @property
    def strand(self) -> int:
        return 0 if len({part.strand for part in self.parts}) > 1 else self.parts[0].strand

    @property
    def span(self) -> tuple[int, int]:
        return (min(part.start for part in self.parts) + 1, max(part.end for part in self.parts))

    def to_gene(self) -> Gene:
        """The call as a feature. No product text: that wording is a submitter's, and a
        vocabulary's synonyms were quoted from the very records these calls are checked against."""
        return Gene(
            name=self.symbol,
            kind=self.kind,
            parts=self.parts,
            transl_table=self.transl_table if self.kind == "CDS" else None,
        )


@dataclass(frozen=True, slots=True)
class CallResult:
    """Everything one call run decided — the calls, and every hit it did not use."""

    calls: tuple[GeneCall, ...] = ()
    not_called: tuple[NotCalled, ...] = ()
    dropped: tuple[DroppedHit, ...] = ()
    unnameable: tuple[GeneCall, ...] = ()
    competing: tuple[CompetingReference, ...] = ()
    findings: tuple[str, ...] = ()
    min_identity: float = 0.0
    min_coverage: float = 0.0

    def one(self, symbol: str) -> GeneCall:
        """The call for one symbol. Raises rather than returning ``None``: a test helper that
        quietly passed on a missing call is how a broken caller stays green."""
        for call in (*self.calls, *self.unnameable):
            if call.symbol == symbol:
                return call
        raise KeyError(f"{symbol!r} was not called; called: {[call.symbol for call in self.calls]}")


@dataclass(frozen=True, slots=True)
class _Candidate:
    """One reference sequence's attempt at a locus, with the floors already applied."""

    entry: ReferenceEntry
    kept: tuple[Hsp, ...]
    dropped: tuple[DroppedHit, ...]
    exons: tuple[tuple[int, int, int, int], ...]  # start, end, strand, subject start
    coverage: float
    identity: float
    bitscore: float
    passes_floors: bool

    @property
    def length(self) -> int:
        return len(self.entry.sequence)


def call_cds(
    genome: Genome,
    references: ReferenceIndex,
    hsps: Sequence[Hsp],
    *,
    min_identity: float,
    min_coverage: float,
    transl_table: int = 1,
    max_hsps_per_subject: int | None = None,
    structure_first: bool = False,
) -> CallResult:
    """Call the protein-coding genes of ``genome`` from a blastx run's hits.

    The floors and the translation table are **required**: each organelle path measures its own
    (§8.22 for the mito path, §8.46 for the plastid one) and neither may inherit the other's.

    ``structure_first`` is the same kind of per-path statement (see :func:`_choose`): it drops a
    candidate whose protein is the C-terminal stretch of another candidate's, which is what a
    record annotating one gene from an inner start codon looks like. The plastid path measured a
    case for it (``rps19``, §8.55) and asks for it; the mito path does not, so it does not.
    """
    if transl_table not in START_CODONS:
        raise ValueError(
            f"translation table {transl_table} is not one this toolkit can translate "
            f"(supported: {sorted(START_CODONS)}) — an initiator test against a table we do not "
            f"know would be a guess"
        )
    return _call(
        "CDS",
        genome,
        references,
        hsps,
        min_identity=min_identity,
        min_coverage=min_coverage,
        transl_table=transl_table,
        max_hsps_per_subject=max_hsps_per_subject,
        structure_first=structure_first,
    )


def transfer_structure(result: CallResult, hits: Mapping[str, Sequence[Alignment]]) -> CallResult:
    """Rebuild every CDS call from its reference's **own** exon structure, where it could be placed.

    The calls a protein search makes are as good as the words it can seed on: an exon shorter than
    a word is invisible (``petB``'s 6 bp, ``petD``'s 8, ``rpl16``'s 9, ``rps16``'s 40 — measured
    missed in all 34 corpus records), and an intron lying *inside* one alignment leaves the
    boundaries around it wherever the local alignment stopped (``rpoC1`` 15 and 12 bp away from
    the published ones). Neither is a scoring problem: no search over the target's own DNA can
    find a 6 bp exon, and none of them is what the reference already knows.

    So the reference's structure is placed instead — see :mod:`plastid.core.exon_structure` — and
    a call is rebuilt around it. A call whose reference has no structure to place (a single-exon
    gene, a trans-spliced one, an rRNA) comes back **unchanged**: this is a transfer of what the
    reference states, not a second opinion about everything.

    Only ``parts`` and the notes move. The floors, the identity, the reference itself and the
    arbitration that chose it are untouched — the placement says *where*, never *which*.
    """
    calls: list[GeneCall] = []
    findings = list(result.findings)
    rebuilt = 0
    for call in result.calls:
        if call.kind != "CDS" or not call.reference.junctions or call.strand == 0:
            calls.append(call)
            continue
        placement = place(
            Transfer(
                symbol=call.symbol,
                strand=call.strand,
                spans=tuple((part.start + 1, part.end) for part in call.parts),
                partial_start=call.partial_start,
                partial_end=call.partial_end,
            ),
            call.reference.junctions,
            hits,
        )
        if placement is None:
            calls.append(call)
            continue
        rebuilt += 1
        findings.extend(placement.notes)
        notes = _notes_kept(call.notes, placement)
        calls.append(
            GeneCall(
                symbol=call.symbol,
                kind=call.kind,
                parts=tuple(
                    Interval(
                        start=span[0] - 1,
                        end=span[1],
                        strand=call.strand,
                        partial_start=placement.partial_start and index == 0,
                        partial_end=placement.partial_end and index == len(placement.spans) - 1,
                    )
                    for index, span in enumerate(placement.spans)
                ),
                coverage=call.coverage,
                identity=call.identity,
                bitscore=call.bitscore,
                reference=call.reference,
                transl_table=call.transl_table,
                notes=notes,
            )
        )
    if not rebuilt:
        return result
    return replace(result, calls=tuple(calls), findings=tuple(findings))


def _notes_kept(notes: Sequence[str], placement: Placement) -> tuple[str, ...]:
    """The call's own notes minus the ones the placement has just answered.

    A note saying the alignment did not reach the reference's 5' end is a fact about the
    *alignment*; once the reference's own structure placed that end, keeping it would put "partial"
    in the table beside coordinates that are no longer partial. Notes about anything else stay —
    this drops what was superseded and nothing else.
    """
    dropped = tuple(
        prefix
        for prefix, established in (
            (NOT_REACHED_FIVE_PRIME, not placement.partial_start),
            (NOT_REACHED_THREE_PRIME, not placement.partial_end),
        )
        if established
    )
    return tuple(note for note in notes if not note.startswith(dropped))


def call_rrna(
    genome: Genome,
    references: ReferenceIndex,
    hsps: Sequence[Hsp],
    *,
    min_identity: float,
    min_coverage: float,
    max_hsps_per_subject: int | None = None,
) -> CallResult:
    """Call the rRNA genes from a blastn run's hits. No frame, so no codon refinement."""
    return _call(
        "rRNA",
        genome,
        references,
        hsps,
        min_identity=min_identity,
        min_coverage=min_coverage,
        transl_table=None,
        max_hsps_per_subject=max_hsps_per_subject,
    )


def _call(
    kind: str,
    genome: Genome,
    references: ReferenceIndex,
    hsps: Sequence[Hsp],
    *,
    min_identity: float,
    min_coverage: float,
    transl_table: int | None,
    max_hsps_per_subject: int | None,
    structure_first: bool = False,
) -> CallResult:
    by_key: dict[str, list[Hsp]] = defaultdict(list)
    for hsp in _deduplicate(hsps):
        by_key[hsp.subject_key].append(hsp)

    findings = list(_cap_findings(by_key, max_hsps_per_subject))
    lengths = references.lengths_by_key()
    unnameable_keys = {entry.key for entry in references.unnameable()}

    dropped: list[DroppedHit] = []
    candidates: list[_Candidate] = []
    below_floor: dict[str, Hsp] = {}
    for key, group in sorted(by_key.items()):
        entry = references.entry_for(key)
        if entry is None:
            dropped.extend(
                DroppedHit("", hsp, f"{key}: not a subject of this reference set") for hsp in group
            )
            continue
        label = entry.vocabulary or entry.source_gene
        above: list[Hsp] = []
        for hsp in group:
            if hsp.pident >= min_identity:
                above.append(hsp)
                continue
            dropped.append(
                DroppedHit(
                    label,
                    hsp,
                    f"identity {hsp.pident:.1f} % is below the floor of {min_identity:.1f} %",
                )
            )
            best = below_floor.get(entry.vocabulary or "")
            if entry.vocabulary and (best is None or hsp.bitscore > best.bitscore):
                below_floor[entry.vocabulary] = hsp
        for length, assigned in sorted(_assign_lengths(above, lengths[key]).items()):
            at_length = references.entry_at(key, length)
            if at_length is None:
                dropped.extend(
                    DroppedHit(label, hsp, f"no reference under {key} is {length} residues long")
                    for hsp in assigned
                )
                continue
            candidates.append(_candidate(at_length, assigned, kind, min_identity, min_coverage))

    nameable = [item for item in candidates if item.entry.key not in unnameable_keys]
    unnamed = [item for item in candidates if item.entry.key in unnameable_keys]

    calls: list[GeneCall] = []
    not_called: list[NotCalled] = []
    competing: list[CompetingReference] = []
    for symbol in references.vocabulary:
        mine = [item for item in nameable if item.entry.vocabulary == symbol]
        if not mine:
            not_called.append(_nothing_found(symbol, below_floor.get(symbol), min_identity))
            continue
        chosen = _choose(mine, structure_first=structure_first)
        competing.extend(_competing_row(symbol, item, chosen) for item in mine)
        if chosen is None:
            not_called.append(_below_the_floors(symbol, mine, min_coverage))
            dropped.extend(
                DroppedHit(symbol, hsp, "its reference candidate did not pass the floors")
                for item in mine
                for hsp in item.kept
            )
            continue
        calls.append(_to_call(genome, symbol, kind, chosen, transl_table))
        dropped.extend(chosen.dropped)

    unnameable = sorted(
        (_to_call(genome, item.entry.source_gene, kind, item, transl_table) for item in unnamed),
        key=lambda call: call.span,
    )
    calls.sort(key=lambda call: call.span)
    findings.extend(_overlap_findings(calls))
    findings.extend(_origin_findings(genome, calls))

    return CallResult(
        calls=tuple(calls),
        not_called=tuple(not_called),
        dropped=tuple(dropped),
        unnameable=tuple(unnameable),
        competing=tuple(competing),
        findings=tuple(findings),
        min_identity=min_identity,
        min_coverage=min_coverage,
    )


def _nothing_found(symbol: str, best: Hsp | None, min_identity: float) -> NotCalled:
    if best is None:
        return NotCalled(symbol, "no alignment to this locus in the reference set")
    return NotCalled(
        symbol,
        f"every alignment is below the identity floor of {min_identity:.1f} % "
        f"(best: {best.pident:.1f} % over {best.length} aligned residues at {best.span_label})",
        best_accession=best.subject_key.split("|")[0],
        best_identity=best.pident,
        best_bitscore=best.bitscore,
    )


def _below_the_floors(symbol: str, mine: Sequence[_Candidate], min_coverage: float) -> NotCalled:
    best = max(mine, key=lambda item: (item.identity, item.bitscore))
    return NotCalled(
        symbol,
        f"the best reference covers only {best.coverage:.0%} of the protein — below the coverage "
        f"floor of {min_coverage:.0%}",
        best_accession=best.entry.accession,
        best_identity=best.identity,
        best_coverage=best.coverage,
        best_bitscore=best.bitscore,
    )


# --------------------------------------------------------------------------- reference indexing


def render_fasta(entries: Iterable[ReferenceEntry]) -> str:
    """The committed artifact, written through BioPython so it is a FASTA every reader knows.

    The header carries all four fields because a hit has to be traceable back to the record it
    came from. ``-`` in the vocabulary field means this toolkit cannot name that gene yet — a fact
    worth being able to grep for.
    """
    handle = StringIO()
    SeqIO.write(
        [
            SeqRecord(
                Seq(entry.sequence),
                id=_header(entry),
                description="",
            )
            for entry in entries
        ],
        handle,
        "fasta",
    )
    return handle.getvalue()


def read_reference(path: str | Path) -> tuple[ReferenceEntry, ...]:
    """Read a committed reference FASTA back — BioPython does the parsing, as everywhere else.

    The inverse of :func:`render_fasta`, and exact only for a name that had no whitespace to
    begin with: a field is written with its spaces replaced, so ``23S ribosomal RNA`` comes back
    as ``23S_ribosomal_RNA``. No measured reference carries such a name — the plastid path builds
    its entries from GenBank records each run and never reads a rendered FASTA back — so the
    round trip is left honest rather than lossy in a way nobody has needed.
    """
    entries: list[ReferenceEntry] = []
    for record in SeqIO.parse(str(path), "fasta"):
        fields = record.description.split("|")
        if len(fields) != 4:
            raise ValueError(
                f"{path}: reference header {record.description!r} is not "
                f"'accession|gene|vocabulary|organism'"
            )
        accession, source_gene, vocabulary, organism = fields
        entries.append(
            ReferenceEntry(
                accession=accession,
                organism=organism,
                source_gene=source_gene,
                vocabulary=None if vocabulary == "-" else vocabulary,
                sequence=str(record.seq).upper(),
            )
        )
    return tuple(entries)


def _header(entry: ReferenceEntry) -> str:
    """The FASTA header: the key a hit comes back as, plus the organism for a human reader.

    The first three fields are ``entry.key`` itself rather than a second rendering of the same
    three values — a header that disagreed with the key would send every hit to ``entry_for`` and
    find nothing (measured 08/10/2026: a plastid rRNA named ``23S ribosomal RNA`` needs the space
    removed in both places or in neither).
    """
    return f"{entry.key}|{' '.join(entry.organism.split())}"


# --------------------------------------------------------------------------- reference indexing


def _assign_lengths(hsps: Sequence[Hsp], lengths: Sequence[int]) -> dict[int, list[Hsp]]:
    """Give every HSP the smallest reference length that can contain its subject range.

    Measured: Psychotria's ``cox2`` exists as a 257 aa and a 390 aa sequence under one key. An
    alignment ending at subject position 380 cannot have come from the 257 aa copy, and letting it
    would divide its coverage by the wrong denominator and chain two different molecules.
    """
    assigned: dict[int, list[Hsp]] = defaultdict(list)
    for hsp in hsps:
        for length in lengths:
            if length >= hsp.shi:
                assigned[length].append(hsp)
                break
    return assigned


def _deduplicate(hsps: Sequence[Hsp]) -> tuple[Hsp, ...]:
    """BLAST can report one HSP twice for a duplicated subject sequence. Identical rows are one
    piece of evidence; counting them twice would weight a locus by its copy number."""
    seen: set[tuple[object, ...]] = set()
    unique: list[Hsp] = []
    for hsp in hsps:
        identity = (hsp.subject_key, hsp.qstart, hsp.qend, hsp.sstart, hsp.send, hsp.frame)
        if identity not in seen:
            seen.add(identity)
            unique.append(hsp)
    return tuple(unique)


# --------------------------------------------------------------------------- chaining


def _chain(hsps: Sequence[Hsp], symbol: str) -> tuple[tuple[Hsp, ...], tuple[DroppedHit, ...]]:
    """Keep the HSPs that add evidence; report the rest with the reason they did not.

    Two ways an HSP fails to add evidence, both measured on the mito answer-key record:

    * it aligns to a query region a better HSP already covers **in another frame** (sdh4, mttB) —
      keeping both would put one gene in two reading frames;
    * every reference residue it aligns is already covered (the 27 aa atp8 fragment at 8,557, and
      the repeat copies inside the reference set itself). A second copy of a gene's part is worth
      *reporting* and is not worth a second exon.
    """
    kept: list[Hsp] = []
    dropped: list[DroppedHit] = []
    for hsp in sorted(hsps, key=lambda item: (-item.bitscore, item.qlo)):
        conflict = next(
            (other for other in kept if other.frame != hsp.frame and other.overlaps(hsp)), None
        )
        if conflict is not None:
            dropped.append(
                DroppedHit(
                    symbol,
                    hsp,
                    f"aligns to {hsp.span_label}, which {conflict.span_label} already covers in "
                    f"another frame — one gene cannot be in two reading frames",
                )
            )
            continue
        gaps = [(hsp.slo, hsp.shi)]
        for other in kept:
            gaps = [
                piece
                for start, end in gaps
                for piece in _subtract(start, end, other.slo, other.shi)
            ]
        if not gaps:
            dropped.append(
                DroppedHit(
                    symbol,
                    hsp,
                    f"a second copy of a part the chosen reference already aligns "
                    f"({hsp.span_label}) — reported rather than called, because this prototype "
                    f"makes one call per locus",
                )
            )
            continue
        kept.append(hsp)
    return tuple(kept), tuple(dropped)


def _merge(intervals: Sequence[tuple[int, int]]) -> list[tuple[int, int]]:
    merged: list[list[int]] = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(start, end) for start, end in merged]


def _exons(kept: Sequence[Hsp]) -> tuple[tuple[int, int, int, int], ...]:
    """``(start, end, strand, subject start)`` per exon, in transcript order.

    Interval merging runs within a frame: two frames over one query region are two readings of it,
    and welding them would invent an exon neither alignment supports. Transcript order comes from
    the subject coordinates — the reference protein is the only witness to which exon comes first,
    and for a trans-spliced gene the genome is no help at all.
    """
    by_frame: dict[tuple[int | None, int], list[Hsp]] = defaultdict(list)
    for hsp in kept:
        by_frame[(hsp.frame, hsp.strand)].append(hsp)
    exons: list[tuple[int, int, int, int]] = []
    for (_frame, strand), group in by_frame.items():
        for start, end in _merge([(hsp.qlo, hsp.qhi) for hsp in group]):
            inside = [hsp for hsp in group if hsp.qlo >= start and hsp.qhi <= end]
            exons.append((start, end, strand, min(hsp.sstart for hsp in inside)))
    exons.sort(key=lambda exon: exon[3])
    return tuple(exons)


def _candidate(
    entry: ReferenceEntry,
    hsps: Sequence[Hsp],
    kind: str,
    min_identity: float,
    min_coverage: float,
) -> _Candidate:
    kept, dropped = _chain(hsps, entry.vocabulary or entry.source_gene)
    exons = _exons(kept)
    covered_bp = sum(end - start + 1 for start, end, _, _ in exons)
    aligned = sum(hsp.length for hsp in kept) or 1
    # a protein reference is in residues, so its span in bases is three times its length;
    # an rRNA reference is already in bases
    reference_bp = 3 * len(entry.sequence) if kind == "CDS" else len(entry.sequence)
    coverage = covered_bp / reference_bp
    identity = sum(hsp.pident * hsp.length for hsp in kept) / aligned
    return _Candidate(
        entry=entry,
        kept=kept,
        dropped=dropped,
        exons=exons,
        coverage=coverage,
        identity=identity,
        bitscore=sum(hsp.bitscore for hsp in kept),
        passes_floors=coverage >= min_coverage and identity >= min_identity,
    )


def _choose(
    candidates: Sequence[_Candidate], *, structure_first: bool = False
) -> _Candidate | None:
    """The highest-identity candidate that passes the floors; bitscore breaks ties.

    Measured on ``atp6``: ranking by bitscore chose a 75.6 %-identity alignment 418 aa long that
    pulled the call 381 bp upstream of the published CDS over a 93.4 % one that matched it. The
    closest relative is what a gene model should be transferred from — not the longest alignment,
    which is how a divergent open reading frame gets welded onto a gene.

    **Structure decides before the score.** H-19 keeps a truncated duplicate out of the *index*;
    this is the same rule *between* candidates, and it is needed for the same reason. Measured on
    the golden's ``rps19`` (10/10/2026): three references spell it 225 bp / 74 aa, and that
    protein is a **suffix** of the 92 aa protein the target's own record carries — the same gene
    annotated from an inner ATG. Those three score **100.00 %** against 98.91 % for the complete
    ones, so identity alone transfers the truncation and the call lands 54 bp short of the
    published CDS. A candidate that is a stretch of another candidate's protein explains a proper
    subset of the target's residues and can never be the better model of it.
    """
    passing = [item for item in candidates if item.passes_floors]
    if not passing:
        return None
    if structure_first:
        passing = _not_truncated(passing)
    return max(passing, key=lambda item: (item.identity, item.bitscore, item.coverage))


def _not_truncated(candidates: Sequence[_Candidate]) -> list[_Candidate]:
    """The candidates whose protein is not the *end* of another candidate's — see :func:`_choose`.

    A **suffix**, deliberately, and not any containment: two records spelling one gene from
    different start codons make the shorter protein the longer one's C-terminal stretch, while an
    arbitrary internal match is as likely to be low-complexity sequence as an annotation. Measured
    over the whole pinned set (10/10/2026): the only containment at any locus is ``rps19``'s three
    74 aa proteins inside the 92 aa ones, and one 92 aa record that is a *prefix* of a 94 aa one —
    which this rule leaves alone, because a gene cut short at its 3' end is a different claim from
    one started late, and the corpus has no measurement that would justify acting on it.
    """
    sequences = {item.entry.sequence for item in candidates}
    truncated = {
        sequence
        for sequence in sequences
        if any(
            other.endswith(sequence) and other != sequence and len(other) > len(sequence)
            for other in sequences
        )
    }
    return [item for item in candidates if item.entry.sequence not in truncated]


def _competing_row(
    symbol: str, candidate: _Candidate, chosen: _Candidate | None
) -> CompetingReference:
    return CompetingReference(
        symbol=symbol,
        reference_accession=candidate.entry.accession,
        reference_gene=candidate.entry.source_gene,
        reference_length=candidate.length,
        identity=candidate.identity,
        coverage=candidate.coverage,
        bitscore=candidate.bitscore,
        passes_floors=candidate.passes_floors,
        chosen=candidate is chosen,
    )


# --------------------------------------------------------------------------- refinement


def _to_call(
    genome: Genome,
    symbol: str,
    kind: str,
    candidate: _Candidate,
    transl_table: int | None,
) -> GeneCall:
    notes: list[str] = []
    exons = [list(exon) for exon in candidate.exons]
    partial_start = _refine_start(genome, candidate, exons, kind, transl_table, notes)
    partial_end = _refine_end(genome, candidate, exons, kind, notes)
    parts = tuple(
        Interval(
            start=start - 1,
            end=end,
            strand=strand,
            partial_start=partial_start and index == 0,
            partial_end=partial_end and index == len(exons) - 1,
        )
        for index, (start, end, strand, _) in enumerate(exons)
    )
    return GeneCall(
        symbol=symbol,
        kind=kind,
        parts=parts,
        coverage=candidate.coverage,
        identity=candidate.identity,
        bitscore=candidate.bitscore,
        reference=candidate.entry,
        transl_table=transl_table,
        notes=tuple(notes),
    )


def _reaches_five_prime(candidate: _Candidate, exon: Sequence[int]) -> bool:
    start, end, strand, _ = exon
    return any(
        hsp.sstart == 1
        and hsp.strand == strand
        and (hsp.qhi if strand < 0 else hsp.qlo) <= (end if strand < 0 else start)
        for hsp in candidate.kept
    )


def _reaches_three_prime(candidate: _Candidate, exon: Sequence[int]) -> bool:
    start, end, strand, _ = exon
    return any(
        hsp.send == candidate.length
        and hsp.strand == strand
        and (hsp.qlo if strand < 0 else hsp.qhi) >= (start if strand < 0 else end)
        for hsp in candidate.kept
    )


def _refine_start(
    genome: Genome,
    candidate: _Candidate,
    exons: list[list[int]],
    kind: str,
    transl_table: int | None,
    notes: list[str],
) -> bool:
    """Is the 5' coordinate established? One codon is read, and nothing is ever moved."""
    reached = _reaches_five_prime(candidate, exons[0])
    if kind != "CDS":
        if not reached:
            notes.append(NOT_REACHED_FIVE_PRIME)
        return not reached
    if not reached:
        notes.append(
            NOT_REACHED_FIVE_PRIME + " — the start is partial, and it is not moved to the next ATG"
        )
        return True
    codon = _codons(genome, exons[0], transl_table or 1)[0]
    if codon not in START_CODONS[transl_table or 1]:
        notes.append(f"the first codon {codon} is not an initiator — a partial start, unmoved")
        return True
    return False


def _refine_end(
    genome: Genome,
    candidate: _Candidate,
    exons: list[list[int]],
    kind: str,
    notes: list[str],
) -> bool:
    """Is the 3' coordinate established? At most one codon past the alignment is ever read.

    **The stop codon is looked for whether or not the alignment reached the reference's own 3'
    end**, because it is what says where *this* gene ends. blastx stops translating at a stop, so
    a target whose gene is shorter than the reference's presents an alignment that simply ends —
    and requiring the alignment to reach the reference's end left every such call one codon short
    of the published CDS. Measured (P3B-10, §8.55): of the corpus's differences that survived the
    exon-structure transfer, **65 were exactly this** — the published end 3 bp past ours with a
    stop codon sitting there, ``psaJ`` (30 records) at the top of them.

    What is *not* done is walking on: one codon is read, never two, and a gene whose next codon is
    not a stop stays a partial end rather than being stretched to the next one.
    """
    exon = exons[-1]
    start, end, strand, _ = exon
    reached = _reaches_three_prime(candidate, exon)
    if kind != "CDS":
        if not reached:
            notes.append(NOT_REACHED_THREE_PRIME)
        return not reached  # rRNA has no stop codon to find or to miss
    new_start, new_end = (start - 3, end) if strand < 0 else (start, end + 3)
    if new_start < 1 or new_end > len(genome):
        notes.append(
            "the next codon would run past the end of the record, so the end is left where the "
            "alignment put it — a circular wrap is not an end this alignment shows"
        )
        return True
    codon = _codons(genome, [new_start, new_end, strand, 0], 1)[-1]
    if codon in STOP_CODONS:
        exon[0], exon[1] = new_start, new_end
        if not reached:
            notes.append(
                "the alignment stops short of the reference 3' end but the next codon is a stop — "
                "this gene is one codon shorter than the reference's, and the stop is its end"
            )
        return False
    if not reached:
        notes.append(NOT_REACHED_THREE_PRIME)
        return True
    notes.append(
        f"the alignment reaches the reference 3' end but the next codon is {codon}, not a stop — "
        f"left as a partial end rather than walked on to the next stop"
    )
    return True


def _codons(genome: Genome, exon: Sequence[int], table: int) -> tuple[str, ...]:
    """The complete codons of one exon, read 5'→3' on its own strand.

    ``table`` is accepted and unused beyond documenting which code the caller is reading in; the
    codon *strings* are what the initiator and stop tests compare against, so a wrong table would
    have to be wrong in the caller's own declaration first.
    """
    start, end, strand, _ = exon
    bases = genome.sequence[start - 1 : end].upper()
    if strand < 0:
        bases = reverse_complement(bases)
    complete = len(bases) // 3
    return tuple(bases[index * 3 : index * 3 + 3] for index in range(complete))


# --------------------------------------------------------------------------- findings


def _cap_findings(by_key: dict[str, list[Hsp]], max_hsps_per_subject: int | None) -> Iterator[str]:
    """-max_hsps truncates in silence: a gene whose exon set was cut off looks like a shorter
    gene. A pair that returns exactly the cap is reported, capped or not."""
    if max_hsps_per_subject is None:
        return
    for key, group in sorted(by_key.items()):
        if len(group) >= max_hsps_per_subject:
            yield (
                f"{key}: {len(group)} HSPs, which is the search's per-subject cap — the exon set "
                f"may be truncated, so this call can only be shorter than the gene"
            )


def _overlap_findings(calls: Sequence[GeneCall]) -> Iterator[str]:
    for index, left in enumerate(calls):
        for right in calls[index + 1 :]:
            left_span, right_span = left.span, right.span
            if left_span[0] <= right_span[1] and right_span[0] <= left_span[1]:
                detail = (
                    "on opposite strands" if left.strand * right.strand < 0 else "on one strand"
                )
                yield (
                    f"{left.symbol} {left_span[0]}-{left_span[1]} and {right.symbol} "
                    f"{right_span[0]}-{right_span[1]} overlap {detail} — both calls stand; a "
                    f"repeat-rich genome can hand two loci the same region"
                )


def _origin_findings(genome: Genome, calls: Sequence[GeneCall]) -> Iterator[str]:
    if genome.topology != "circular":
        return
    for call in calls:
        if len(call.parts) > 1 and call.span[0] <= 30 and call.span[1] >= len(genome) - 30:
            yield (
                f"{call.symbol}: its exons meet at both ends of the record — it may be one gene "
                f"crossing the origin, which this prototype reports rather than resolves"
            )
