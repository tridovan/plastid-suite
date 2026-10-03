"""GFF3 reader for annotator output (M3.1, P3-03).

An annotator that emits GFF rather than GenBank still has to enter the same pipeline as a
downloaded record, so its rows are read here and become the same
:class:`~plastid.core.model.Gene` objects the GenBank reader produces.

**Why not BioPython.** Every other format in this package goes through BioPython, and this is the
one where that is not possible: ``Bio.SeqIO`` has no GFF3 parser (GFF3 lives in the separate,
unmaintained ``bcbio-gff``), so a GFF reader is ours to write. It is kept deliberately small and
strict for that reason — it reads the subset the annotators emit and refuses anything it would
otherwise have to guess at.

Two conventions are read from the file, never imposed:

* **one feature per id.** The rows sharing a ``parent`` (a CDS/tRNA/rRNA's exons) or one ``ID``
  (a ``gene`` row's parts) are one feature, and their parts keep **file order**. Measured
  2026-10-03 on Chloe's GFF for the golden genome: its rows come in transcript order 5'→3' — the
  trans-spliced rps12 copies read back as [LSC exon, IR exon, IR exon], exactly the order
  BioPython gives for the record's own ``complement(join(...))``, and the only order that
  translates to ``MPTIKQLIRNTRQPIRNVTKSPALRGCPQRRGTCTRVYIT``. Sorting the parts by coordinate
  instead would silently reverse the reading frame of every minus-strand spliced gene.
* **no name, no feature.** An unnamed row is refused: it would be counted as a gene and then
  dropped by every name-based analysis, which is the kind of half-counted number this project
  exists to avoid.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass

from plastid.core.model import Gene, Interval

__all__ = ["FEATURE_KINDS", "GffError", "GffRow", "genes_from_gff", "read_gff_rows"]

# The kinds that carry genes — the same four the GenBank reader keeps. Everything else a GFF
# holds (intron, exon, repeat_region, …) is derived by the toolkit from these, not read from here.
FEATURE_KINDS: tuple[str, ...] = ("gene", "CDS", "tRNA", "rRNA")

_STRANDS = {"+": 1, "-": -1}
_PSEUDO_TRUE = frozenset({"true", "yes", "1"})
_FIELDS = 9


class GffError(ValueError):
    """A GFF row cannot be read as the feature it claims to be."""


@dataclass(frozen=True, slots=True)
class GffRow:
    """One GFF3 line, coordinates 1-based inclusive exactly as written."""

    seq_id: str
    kind: str
    start: int
    end: int
    strand: str
    attributes: Mapping[str, str]
    line_number: int

    @property
    def feature_id(self) -> str | None:
        return self.attributes.get("ID")

    @property
    def parent(self) -> str | None:
        return self.attributes.get("parent")

    @property
    def group(self) -> str | None:
        """The feature this row belongs to: a child under its parent, a parent under itself."""
        return self.parent or self.feature_id

    @property
    def name(self) -> str | None:
        return self.attributes.get("gene") or self.attributes.get("name")

    @property
    def pseudo(self) -> bool:
        value = self.attributes.get("pseudo")
        return value is not None and value.strip().lower() in _PSEUDO_TRUE


def read_gff_rows(text: str) -> tuple[GffRow, ...]:
    """Every feature row in a GFF3 text. Directives and comments are skipped; ``##FASTA`` ends it.

    The line number is carried on every row so a refusal can name the line that caused it — a
    reader that says "malformed GFF" and nothing else costs more than it saves.
    """
    rows: list[GffRow] = []
    for number, line in enumerate(text.splitlines(), start=1):
        if line.startswith("##FASTA"):
            break  # the sequence follows; it is not part of the feature table
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != _FIELDS:
            raise GffError(
                f"line {number}: expected {_FIELDS} tab-separated fields, found {len(fields)}: "
                f"{line[:120]!r}"
            )
        seq_id, _source, kind, start, end, _score, strand, _phase, attributes = fields
        rows.append(
            GffRow(
                seq_id=seq_id,
                kind=kind,
                start=_coordinate(start, "start", number),
                end=_coordinate(end, "end", number),
                strand=strand,
                attributes=_attributes(attributes, number),
                line_number=number,
            )
        )
    return tuple(rows)


def genes_from_gff(text: str, *, kinds: Iterable[str] = FEATURE_KINDS) -> tuple[Gene, ...]:
    """The features of ``kinds`` in a GFF3 text, one :class:`Gene` per id, in file order."""
    wanted = frozenset(kinds)
    rows = [row for row in read_gff_rows(text) if row.kind in wanted]
    rows_by_id = {row.feature_id: row for row in rows if row.feature_id}

    grouped: dict[tuple[str, str], list[GffRow]] = {}
    for row in rows:
        group = row.group
        if group is None:
            raise GffError(
                f"line {row.line_number}: {row.kind} row has neither ID nor parent — there is no "
                f"way to tell which feature it belongs to"
            )
        grouped.setdefault((row.kind, group), []).append(row)

    return tuple(
        Gene(
            name=_name(group, rows_by_id),
            kind=kind,
            parts=tuple(_interval(row) for row in group),
            pseudo=any(row.pseudo for row in group),
            transl_table=_transl_table(group),
        )
        for (kind, _group_id), group in grouped.items()
    )


def _interval(row: GffRow) -> Interval:
    strand = _STRANDS.get(row.strand)
    if strand is None:
        raise GffError(
            f"line {row.line_number}: strand is {row.strand!r} — a feature whose direction is "
            f"unknown has no reading frame, so it is refused rather than assumed forward"
        )
    return Interval.from_1based(row.start, row.end, strand)


def _name(group: list[GffRow], rows_by_id: Mapping[str | None, GffRow]) -> str:
    """The feature's name: its own ``gene``/``name`` attribute, else its parent row's.

    A GFF may carry the name only on the ``gene`` row, so the parent is consulted before giving
    up — and giving up is an error, never a placeholder.
    """
    row = group[0]
    name = next((candidate.name for candidate in group if candidate.name), None)
    if name is None and row.parent is not None:
        parent = rows_by_id.get(row.parent)
        name = parent.name if parent is not None else None
    if not name:
        raise GffError(
            f"line {row.line_number}: {row.kind} feature {row.group!r} has no gene or name "
            f"attribute — an unnamed feature cannot be counted or checked"
        )
    return name


def _transl_table(group: list[GffRow]) -> int | None:
    """The declared translation table, or ``None`` when the file declares none.

    GFF3 has no standard attribute for it; some annotators write ``transl_table`` anyway. ``None``
    means "the file is silent" and the caller's configured default applies — never a table this
    reader picked to make something pass.
    """
    declared = {row.attributes["transl_table"] for row in group if "transl_table" in row.attributes}
    if not declared:
        return None
    if len(declared) > 1:
        raise GffError(
            f"line {group[0].line_number}: the rows of one feature declare different translation "
            f"tables ({sorted(declared)}) — they cannot all be right"
        )
    value = declared.pop()
    try:
        return int(value)
    except ValueError:
        raise GffError(
            f"line {group[0].line_number}: transl_table={value!r} is not an integer"
        ) from None


def _coordinate(field: str, which: str, line_number: int) -> int:
    try:
        return int(field)
    except ValueError:
        raise GffError(
            f"line {line_number}: {which} coordinate {field!r} is not an integer"
        ) from None


def _attributes(field: str, line_number: int) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for token in field.split(";"):
        if not token.strip():
            continue
        key, separator, value = token.partition("=")
        if not separator:
            raise GffError(
                f"line {line_number}: attribute {token!r} is not key=value — a value this reader "
                f"silently skipped could be the one that mattered"
            )
        attributes[key.strip()] = value.strip()
    return attributes
