"""NCBI submission files — ``.fsa`` and the five-column ``.tbl`` (M3.4, P3-05, D3-04).

Written to NCBI's documented rules rather than from memory, because every one of them fails at
*upload* time, after the work is done:

* five tab-separated columns — start · stop · feature key · qualifier name · qualifier value;
* a multi-interval feature (a join) is written by **stacking** its intervals, each on its own
  line, before its qualifier lines — there is no ``join(...)`` keyword in this format;
* a feature on the minus strand is written by **reversing** its interval (higher coordinate
  first), one interval at a time, so a mixed-strand trans-spliced gene is expressible;
* ``<`` always goes in column 1 and ``>`` always in column 2, **regardless of strand**, and only
  on the feature's own terminal intervals — NCBI rejects a partial marker on an internal one.

The partial flags travel from the record into :class:`~plastid.core.model.Interval`: dropping
one silently turns "this gene runs off the contig" into "this gene starts here", which is a
different claim about the same sequence.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from dataclasses import dataclass, field, replace

from plastid.core.model import Gene, Genome, gene_qualifiers
from plastid.core.orf import (
    START_CLASS_INVALID,
    START_CLASS_NON_CANONICAL,
    OrfError,
    check_cds,
)

__all__ = [
    "TblFeature",
    "TblRecord",
    "parse_tbl",
    "partial_start_genes",
    "source_qualifiers_left_out",
    "submission_problems",
    "to_fasta",
    "to_tbl",
]

# The feature kinds the reader keeps and the submission carries. `gene` features are parents:
# they say where the locus is; CDS/tRNA/rRNA say what it is.
SUBMITTED_KINDS: tuple[str, ...] = ("gene", "CDS", "tRNA", "rRNA")

_QUALIFIER_KEYS: tuple[str, ...] = ("gene", "product", "transl_table", "codon_start")

# Source qualifiers that are not qualifiers *in this format*. Measured 2026-10-03 (§8.18):
# table2asn calls both "Unrecognized qualifier name" and exits non-zero while still writing a
# correct .sqn — a submission table carrying them is noisier than one that does not, and NCBI
# derives both anyway (the molecule type from the submission, the taxon from the organism).
_SOURCE_QUALIFIERS_NOT_IN_TABLE = frozenset({"mol_type", "db_xref"})


@dataclass(frozen=True, slots=True)
class TblFeature:
    """One feature as read back from a ``.tbl`` — the shape a round-trip check compares."""

    kind: str
    parts: tuple[tuple[int, int, bool, bool], ...]  # (start_1b, end_1b, partial_start, partial_end)
    qualifiers: dict[str, str | None] = field(default_factory=dict)


@dataclass(frozen=True, slots=True)
class TblRecord:
    """One ``>Features <id>`` block. A table file may hold several — real ones do.

    Measured 2026-10-03 on NCBI's own VADR output: one file carried eight records. Merging them
    into one feature list would silently attach every record's features to one sequence, so a
    record boundary is part of the model rather than something a reader is expected to notice.
    """

    seq_id: str
    features: tuple[TblFeature, ...]


def to_fasta(genome: Genome, seq_id: str, *, width: int = 60) -> str:
    """The ``.fsa``: one sequence, ID matching the ``>Features`` line of the table."""
    if not seq_id or " " in seq_id:
        raise ValueError(f"sequence id {seq_id!r} must be non-empty and free of spaces")
    body = "\n".join(
        genome.sequence[index : index + width] for index in range(0, len(genome), width)
    )
    return f">{seq_id}\n{body}\n"


def to_tbl(genome: Genome, seq_id: str) -> str:
    """The five-column feature table for one record.

    Features are written in position order — the order the format expects a reader (and NCBI)
    to meet them — and every interval of one feature stays together on consecutive lines.

    Two things are added to what the gene model holds, both measured against NCBI's own validator
    (PROJECT.md §8.18):

    * the **source feature** comes first, carrying the record's own qualifiers. Without it
      table2asn reports ``NoSourceDescriptor`` and ``BioSourceMissing`` — a submission has to say
      what the sequence is, and neither this writer nor the toolkit can invent that;
    * a CDS whose first codon is **not an initiator** for its translation table gets a partial
      start (``<``). NCBI rejects the feature otherwise ("Illegal start codon used … or protein
      should be partial"), and the published plastome submission this was checked against writes
      exactly that for its RNA-edited ``rpl2``. The coordinates do not move — only the marker is
      added — and :func:`partial_start_genes` names the features it happened to.
    """
    marked = frozenset(partial_start_genes(genome))
    features = [
        _feature_of(_marked_partial_start(gene) if gene in marked else gene)
        for gene in sorted(
            (gene for gene in genome.genes if gene.kind in SUBMITTED_KINDS),
            key=_position_key,
        )
    ]
    source = _source_feature(genome)
    return tbl_text((source, *features) if source is not None else features, seq_id)


def partial_start_genes(genome: Genome) -> tuple[Gene, ...]:
    """The CDS whose first codon their own translation table does not initiate.

    Measured 2026-10-03 on the golden genome: ``rpl2`` starts with ACG — not an initiator in
    table 11, because C-to-U editing turns it into AUG in the transcript — and NCBI's validator
    rejects the feature as written. With the start marked partial the same CDS validates, which is
    also what the published *Dioscorea alata* plastome submission does for its ``rpl2`` and
    ``ndhD``. A CDS this function cannot judge (a table the toolkit does not translate, a
    ``/codon_start`` other than 1) is left alone: an honest "no" beats a guessed marker.
    """
    marked: list[Gene] = []
    length = len(genome)
    for gene in genome.genes:
        if gene.kind != "CDS":
            continue
        if any(part.start < 0 or part.end > length for part in gene.parts):
            # a feature that runs off the sequence has no start codon to judge; it is a problem
            # in its own right, and submission_problems names it
            continue
        try:
            check = check_cds(genome, gene)
        except OrfError:
            continue
        if check.start_class in (START_CLASS_INVALID, START_CLASS_NON_CANONICAL):
            # both classes mean "not an initiator in this table": `invalid` is a codon no table
            # uses, `non_canonical` is ACG — the RNA-edited start P3-01 already reports as such,
            # and the one NCBI's validator rejects here. The ORF module's judgement and NCBI's
            # agree on this feature, which is why the rule can reuse it rather than re-derive it.
            marked.append(gene)
    return tuple(marked)


def _marked_partial_start(gene: Gene) -> Gene:
    """The same gene with its 5' coordinate marked partial — a new object, never a mutation.

    The 5' end of a feature is the *first* part in annotation order (that is what
    :func:`~plastid.core.model.coding_sequence` reads first), and on the minus strand it is that
    part's higher coordinate — so the marker goes on ``partial_end`` there, not on
    ``partial_start``.
    """
    first = gene.parts[0]
    marked = (
        replace(first, partial_end=True) if first.strand < 0 else replace(first, partial_start=True)
    )
    return replace(gene, parts=(marked, *gene.parts[1:]))


def source_qualifiers_left_out(genome: Genome) -> tuple[str, ...]:
    """The record's source qualifiers this format has no place for, in the record's own order.

    Reported rather than dropped in silence: they stay on the genome, and a caller that wants to
    say what was left behind has the list.
    """
    return tuple(key for key in genome.source if key in _SOURCE_QUALIFIERS_NOT_IN_TABLE)


def _source_feature(genome: Genome) -> TblFeature | None:
    """The record's source feature, or None when the record states nothing about itself."""
    qualifiers = {
        key: value
        for key, value in genome.source.items()
        if key not in _SOURCE_QUALIFIERS_NOT_IN_TABLE
    }
    if not qualifiers:
        return None
    return TblFeature(
        kind="source",
        parts=((1, len(genome), False, False),),
        qualifiers=qualifiers,
    )


def tbl_text(features: Sequence[TblFeature], seq_id: str) -> str:
    """Serialise features as NCBI's five-column table.

    The one writer for this format, so a table we build from a genome and a table we read back
    from a file cannot drift apart: everything the writer knows about the format lives here and
    in :func:`_feature_of`, and ``parse_tbl`` is its inverse (a round-trip test on real
    NCBI-produced tables holds the two together — see tests/unit/test_submission.py).
    """
    lines = [f">Features {seq_id}"]
    for feature in features:
        lines.extend(_feature_lines(feature))
    return "\n".join(lines) + "\n"


def _feature_of(gene: Gene) -> TblFeature:
    """A gene as the table writes it: intervals in the order and orientation the format needs."""
    parts: list[tuple[int, int, bool, bool]] = []
    last = len(gene.parts) - 1
    for index, part in enumerate(gene.parts):
        start, end = part.to_1based()
        partial_start = _check_partial(gene, index, last, part.partial_start)
        partial_end = _check_partial(gene, index, last, part.partial_end)
        if part.strand < 0:
            # minus strand: the higher coordinate is the feature's 5' end and comes first; the
            # partial flags stay attached to their own coordinate, not to the column
            start, end = end, start
            partial_start, partial_end = partial_end, partial_start
        parts.append((start, end, partial_start, partial_end))
    # the same qualifiers the GenBank writer emits — one function, so a feature cannot come out
    # one way in the table and another way in the file submitted alongside it
    return TblFeature(kind=gene.kind, parts=tuple(parts), qualifiers=gene_qualifiers(gene))


def _feature_lines(feature: TblFeature) -> list[str]:
    lines: list[str] = []
    for index, (start, end, partial_start, partial_end) in enumerate(feature.parts):
        # both markers are prefixes: NCBI's own tables write the stop column as ``>1050``
        first_column = f"{'<' if partial_start else ''}{start}"
        second_column = f"{'>' if partial_end else ''}{end}"
        key = feature.kind if index == 0 else ""
        # three fields, not five: NCBI's own tables do not pad an interval line with empty
        # qualifier columns (the round-trip test against their files is what settled this)
        lines.append(f"{first_column}\t{second_column}\t{key}")
    lines.extend(
        f"\t\t\t{name}\t{'' if value is None else value}"
        for name, value in feature.qualifiers.items()
    )
    return lines


def _position_key(gene: Gene) -> tuple[int, int, str]:
    return (min(part.start for part in gene.parts), max(part.end for part in gene.parts), gene.kind)


def _check_partial(gene: Gene, index: int, last: int, partial: bool) -> bool:
    """A partial marker is legal only at the feature's own ends — refuse rather than drop it."""
    if partial and index not in (0, last):
        raise ValueError(
            f"{gene.name}: an internal interval carries a partial coordinate — NCBI rejects a "
            f"partial marker anywhere but the feature's own ends, so this feature cannot be "
            f"written truthfully"
        )
    return partial


def parse_tbl(text: str) -> tuple[TblRecord, ...]:
    """Read a ``.tbl`` back into records — the round-trip check's other half.

    Accepts both qualifier layouts seen in the wild: NCBI's documented one (qualifier in column
    4) and the compact one the BankIt help page renders (qualifier in column 3). It is our own
    writer's output that must be unambiguous, and that is the four-column form.
    """
    records: list[TblRecord] = []
    seq_id = ""
    features: list[TblFeature] = []
    kind: str | None = None
    parts: list[tuple[int, int, bool, bool]] = []
    qualifiers: dict[str, str | None] = {}

    def flush_feature() -> None:
        if kind is not None:
            features.append(TblFeature(kind=kind, parts=tuple(parts), qualifiers=dict(qualifiers)))

    def flush_record() -> None:
        flush_feature()
        if seq_id:
            records.append(TblRecord(seq_id=seq_id, features=tuple(features)))
        features.clear()

    for raw in text.splitlines():
        if raw.startswith(">"):
            flush_record()
            kind, parts, qualifiers = None, [], {}
            # ">Features <id>" or ">Feature <id>": both spellings are in the wild
            seq_id = raw.lstrip(">").split(" ", 1)[1].strip() if " " in raw else ""
            continue
        if not raw.strip():
            continue
        columns = raw.split("\t")
        if len(columns) < 5:
            columns += [""] * (5 - len(columns))
        first, second, key, qualifier, value = columns[:5]

        if _looks_like_coordinate(first) or _looks_like_coordinate(second):
            if key.strip():  # an interval line carrying the key opens a new feature
                flush_feature()
                parts, qualifiers = [], {}
                kind = key.strip()
            if kind is None:
                raise ValueError(f"a feature interval appears before any feature key: {raw!r}")
            parts.append(_parse_interval(first.strip(), second.strip()))
            continue

        # A qualifier line. NCBI documents the name in column 4 with its value in column 5;
        # the help pages render the same tables with the leading tabs collapsed, so the name
        # turns up in column 3, or at the left margin with its value beside it. Read all three:
        # a table pasted from the help page is how most people first meet this format.
        if qualifier.strip():
            name, text_value = qualifier.strip(), value.strip()
        elif key.strip():
            name, text_value = key.strip(), value.strip()
        else:
            name, text_value = first.strip(), second.strip()
        if name:
            qualifiers[name] = text_value or None
    flush_record()
    return tuple(records)


def _looks_like_coordinate(text: str) -> bool:
    """A location field is digits, optionally marked partial with ``<`` or ``>``.

    Used to tell the two kinds of line apart in both qualifier layouts: the compact one puts a
    qualifier name where the documented one puts a feature key, and only the location fields are
    unambiguously numbers.
    """
    return bool(re.fullmatch(r"[<>]?\d+", text.strip()))


def _parse_interval(first: str, second: str) -> tuple[int, int, bool, bool]:
    partial_start = first.startswith("<")
    partial_end = second.startswith(">")
    start = int(first.lstrip("<>"))
    end = int(second.lstrip("<>"))
    return (start, end, partial_start, partial_end)


def submission_problems(genome: Genome) -> tuple[str, ...]:
    """The first-pass checks a submission must pass — named, never silently repaired.

    ``locus_tag`` is not in the model and is not written here (table2asn generates them), so the
    duplicate check is over the feature identity that would produce one: the same kind, name and
    span twice in one record.
    """
    problems: list[str] = []
    length = len(genome)
    seen: dict[tuple[str, str, tuple[tuple[int, int], ...]], int] = {}

    for gene in genome.genes:
        if gene.kind not in SUBMITTED_KINDS:
            continue
        for part in gene.parts:
            if part.end > length or part.start < 0:
                problems.append(
                    f"{gene.name}: {gene.kind} at {part.start + 1}-{part.end} falls outside the "
                    f"sequence ({length} bp)"
                )
        if gene.kind == "CDS":
            if gene.transl_table is None:
                problems.append(f"{gene.name}: CDS declares no /transl_table")
            if not gene.product:
                problems.append(f"{gene.name}: CDS carries no /product — NCBI requires one")
        identity = (
            gene.kind,
            gene.name,
            tuple(sorted((part.start, part.end) for part in gene.parts)),
        )
        seen[identity] = seen.get(identity, 0) + 1

    for (kind, name, _), count in seen.items():
        if count > 1:
            problems.append(
                f"{name}: {count} identical {kind} features ({count} identical locus tags "
                f"would be generated)"
            )

    if not genome.source.get("organism"):
        problems.append(
            "the record states no source organism — NCBI's validator rejects a table without "
            "one, and this toolkit will not guess an organism: read the record that carries it, "
            "or set the source qualifiers on the genome before submitting"
        )
    return tuple(problems)
