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
from dataclasses import dataclass, field

from plastid.core.model import Gene, Genome

__all__ = [
    "TblFeature",
    "parse_tbl",
    "submission_problems",
    "to_fasta",
    "to_tbl",
]

# The feature kinds the reader keeps and the submission carries. `gene` features are parents:
# they say where the locus is; CDS/tRNA/rRNA say what it is.
SUBMITTED_KINDS: tuple[str, ...] = ("gene", "CDS", "tRNA", "rRNA")

_QUALIFIER_KEYS: tuple[str, ...] = ("gene", "product", "transl_table", "codon_start")


@dataclass(frozen=True, slots=True)
class TblFeature:
    """One feature as read back from a ``.tbl`` — the shape a round-trip check compares."""

    kind: str
    parts: tuple[tuple[int, int, bool, bool], ...]  # (start_1b, end_1b, partial_start, partial_end)
    qualifiers: dict[str, str | None] = field(default_factory=dict)


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
    """
    lines = [f">Features {seq_id}"]
    for gene in sorted(
        (gene for gene in genome.genes if gene.kind in SUBMITTED_KINDS),
        key=_position_key,
    ):
        lines.extend(_feature_lines(gene))
    return "\n".join(lines) + "\n"


def _position_key(gene: Gene) -> tuple[int, int, str]:
    return (min(part.start for part in gene.parts), max(part.end for part in gene.parts), gene.kind)


def _feature_lines(gene: Gene) -> list[str]:
    lines: list[str] = []
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
        # both markers are prefixes: NCBI's own example writes the stop column as ``>1050``
        first_column = f"{'<' if partial_start else ''}{start}"
        second_column = f"{'>' if partial_end else ''}{end}"
        key = gene.kind if index == 0 else ""
        lines.append(f"{first_column}\t{second_column}\t{key}\t\t")
    lines.extend(_qualifier_lines(gene))
    return lines


def _check_partial(gene: Gene, index: int, last: int, partial: bool) -> bool:
    """A partial marker is legal only at the feature's own ends — refuse rather than drop it."""
    if partial and index not in (0, last):
        raise ValueError(
            f"{gene.name}: an internal interval carries a partial coordinate — NCBI rejects a "
            f"partial marker anywhere but the feature's own ends, so this feature cannot be "
            f"written truthfully"
        )
    return partial


def _qualifier_lines(gene: Gene) -> list[str]:
    qualifiers: list[tuple[str, str | None]] = [("gene", gene.name)]
    if gene.product:
        qualifiers.append(("product", gene.product))
    if gene.kind == "CDS":
        if gene.transl_table is not None:
            qualifiers.append(("transl_table", str(gene.transl_table)))
        if gene.codon_start != 1:
            qualifiers.append(("codon_start", str(gene.codon_start)))
    if gene.trans_splicing:
        qualifiers.append(("trans_splicing", None))
    if gene.pseudo:
        qualifiers.append(("pseudo", None))
    return [f"\t\t\t{name}\t{'' if value is None else value}" for name, value in qualifiers]


def parse_tbl(text: str) -> tuple[TblFeature, ...]:
    """Read a ``.tbl`` back into features — the round-trip check's other half.

    Accepts both qualifier layouts seen in the wild: NCBI's documented one (qualifier in column
    4) and the compact one the BankIt help page renders (qualifier in column 3). It is our own
    writer's output that must be unambiguous, and that is the four-column form.
    """
    features: list[TblFeature] = []
    kind: str | None = None
    parts: list[tuple[int, int, bool, bool]] = []
    qualifiers: dict[str, str | None] = {}

    def flush() -> None:
        if kind is not None:
            features.append(TblFeature(kind=kind, parts=tuple(parts), qualifiers=dict(qualifiers)))

    for raw in text.splitlines():
        if raw.startswith(">"):
            flush()
            kind, parts, qualifiers = None, [], {}
            continue
        if not raw.strip():
            continue
        columns = raw.split("\t")
        if len(columns) < 5:
            columns += [""] * (5 - len(columns))
        first, second, key, qualifier, value = columns[:5]

        if _looks_like_coordinate(first) or _looks_like_coordinate(second):
            if key.strip():  # an interval line carrying the key opens a new feature
                flush()
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
    flush()
    return tuple(features)


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
    return tuple(problems)
