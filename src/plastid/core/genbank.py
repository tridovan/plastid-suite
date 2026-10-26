"""GenBank reader and writer, backed by BioPython.

Never hand-roll a GenBank parser (CLAUDE.md, banned outright): BioPython already
handles ``join()``, ``complement()``, trans-splicing, and partial positions.
"""

from __future__ import annotations

from datetime import datetime
from io import StringIO
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import (
    AfterPosition,
    BeforePosition,
    CompoundLocation,
    SeqFeature,
    SimpleLocation,
    UnknownPosition,
)
from Bio.SeqRecord import SeqRecord

from plastid.core.model import Gene, Genome, Interval, gene_qualifiers

__all__ = ["read_fasta", "read_genbank", "read_genbank_records", "to_genbank"]

_GENE_FEATURE_TYPES = frozenset({"gene", "CDS", "tRNA", "rRNA"})


def read_genbank(path: str | Path) -> Genome:
    """Parse a single-record GenBank file into a :class:`Genome`."""
    records = read_genbank_records(path)
    if len(records) != 1:
        raise ValueError(f"{path}: expected exactly one record, found {len(records)}")
    return records[0]


def read_fasta(path: str | Path) -> Genome:
    """Parse a single-record FASTA into a :class:`Genome` with **no features**.

    A sample arrives as a bare sequence from an assembler, and annotation is the step that gives
    it features — so this reader produces exactly the genome an annotator needs and nothing more.
    Circular is the default (an organelle is a circle; `to_fasta` writes no topology, so nothing
    is being assumed away), and the header's first token is the accession because that is what
    :func:`plastid.core.submission.to_fasta` writes and therefore what comes back.
    """
    records = list(SeqIO.parse(str(path), "fasta"))
    if not records:
        raise ValueError(f"{path}: no FASTA records found")
    if len(records) != 1:
        raise ValueError(
            f"{path}: expected exactly one record, found {len(records)} — an assembly step hands "
            f"over one organelle at a time, and annotating a mixture would mix their coordinates"
        )
    record = records[0]
    return Genome(
        name=record.description or record.id,
        accession=record.id,
        sequence=str(record.seq).upper(),
        topology="circular",
    )


def read_genbank_records(path: str | Path) -> list[Genome]:
    """Parse every record in a (possibly multi-record) GenBank file."""
    records = list(SeqIO.parse(str(path), "genbank"))
    if not records:
        raise ValueError(f"{path}: no GenBank records found")
    return [_to_genome(record) for record in records]


def _to_genome(record: SeqRecord) -> Genome:
    genes = tuple(
        _to_gene(feature) for feature in record.features if feature.type in _GENE_FEATURE_TYPES
    )
    repeat_regions = tuple(
        _to_interval(location)
        for feature in record.features
        if feature.type == "repeat_region"
        for location in getattr(feature.location, "parts", [feature.location])
    )
    return Genome(
        name=str(record.annotations.get("organism", record.id)),
        accession=record.id,
        sequence=str(record.seq).upper(),
        topology=str(record.annotations.get("topology", "circular")).lower(),
        genes=genes,
        repeat_regions=repeat_regions,
        source=_source_qualifiers(record),
    )


def _source_qualifiers(record: SeqRecord) -> dict[str, str]:
    """The ``source`` feature's qualifiers, as the record states them.

    Read from the feature, not from the record's annotations: the feature is what a submission
    has to reproduce, and it carries the fields a sample needs (``cultivar``, ``specimen_voucher``
    …) that the record-level annotations do not. A record without one yields an empty mapping —
    a real state, reported by ``submission_problems`` rather than filled in with a guess.
    """
    for feature in record.features:
        if feature.type == "source":
            return {
                key: str(values[0])
                for key, values in feature.qualifiers.items()
                if values and isinstance(values[0], str)
            }
    return {}


def _to_gene(feature: SeqFeature) -> Gene:
    locations = getattr(feature.location, "parts", [feature.location])
    return Gene(
        name=_feature_name(feature),
        kind=feature.type,
        parts=tuple(_to_interval(location) for location in locations),
        product=_first_qualifier(feature, "product"),
        pseudo="pseudo" in feature.qualifiers,
        trans_splicing="trans_splicing" in feature.qualifiers,
        transl_table=_int_qualifier(feature, "transl_table"),
        codon_start=_int_qualifier(feature, "codon_start") or 1,
        exception=_first_qualifier(feature, "exception"),
        translation=_translation(feature),
    )


def _translation(feature: SeqFeature) -> str | None:
    """The record's own ``/translation``, whitespace and all removed.

    GenBank wraps a long translation across lines and our reader gets the joined value from
    BioPython; the internal spaces some records carry are removed here so the sequence is what
    the record means. Never reconstructed by translating the CDS (D5-08).
    """
    raw = _first_qualifier(feature, "translation")
    if raw is None:
        return None
    return "".join(raw.split()).upper()


def _int_qualifier(feature: SeqFeature, key: str) -> int | None:
    """A numeric qualifier as an int; a malformed value raises instead of being ignored.

    Silently dropping ``/codon_start=2`` would translate a frame the record did not ask for —
    the class of bug this project keeps meeting (CLAUDE.md → fail loud).
    """
    raw = _first_qualifier(feature, key)
    if raw is None:
        return None
    try:
        return int(raw)
    except ValueError:
        gene_name = _feature_name(feature)
        raise ValueError(
            f"{gene_name}: /{key}={raw!r} is not an integer — the record is malformed"
        ) from None


def _to_interval(location: SimpleLocation) -> Interval:
    # BioPython positions are already 0-based half-open.
    strand = location.strand
    return Interval(
        start=int(location.start),
        end=int(location.end),
        strand=strand if strand in (1, -1) else 1,
        partial_start=_is_partial(location.start),
        partial_end=_is_partial(location.end),
    )


def _is_partial(position: object) -> bool:
    """True for an incomplete coordinate — ``<``, ``>``, or an unknown position (``?``)."""
    return isinstance(position, BeforePosition | AfterPosition | UnknownPosition)


def to_genbank(genome: Genome, *, translations: bool = False) -> str:
    """Write a :class:`Genome` out as a single-record GenBank file.

    This is how an annotation produced by this toolkit becomes a file: the annotator's output is
    written here, then read back by :func:`read_genbank` before anything is pinned or submitted —
    a file that does not re-read is not a result.

    **What it round-trips:** sequence, accession, organism name, topology, the ``source``
    qualifiers, every feature this model holds (kind, name, parts in transcript order, strand,
    partial flags, product, ``/pseudo``, ``/trans_splicing``, ``/transl_table``, ``/codon_start``,
    ``/exception``) and the ``repeat_region`` spans.

    **What it does not:** qualifiers this model never held — ``/locus_tag``, ``/db_xref``,
    ``/protein_id``, ``/note`` — and features outside ``GENE_KINDS`` plus ``repeat_region``
    (``misc_feature``). The LOCUS metadata is regenerated: the division is ``PLN`` and the date is
    today's, because the model does not carry either. So this writes a record *of the model*, not
    a byte-for-byte copy of the file it was read from.

    ``translations=True`` adds the record's own ``/translation`` back for every CDS that carries
    one. Off by default, and the split is a decision rather than a default: ``/translation`` is
    NCBI's to assign on a record this toolkit *annotates*, but a curated record being
    re-normalized is a different thing — that protein is exactly what the annotator searches
    (D5-08), and a normalization that dropped it would silently disarm the step behind it
    (``gene_normalize``'s artifact is the worked example).
    """
    handle = StringIO()
    SeqIO.write(_to_record(genome, translations=translations), handle, "genbank")
    return handle.getvalue()


def _to_record(genome: Genome, *, translations: bool = False) -> SeqRecord:
    length = len(genome)
    # a GenBank record always carries an identifier, so a Genome without an accession gets one
    # derived from its name — that is the one field of the round trip that cannot come back as
    # None, and the alternative is a LOCUS line BioPython refuses to write
    identifier = _identifier(genome)
    record = SeqRecord(Seq(genome.sequence), id=identifier, name=identifier[:16])
    record.annotations = {
        "molecule_type": "DNA",
        "topology": genome.topology,
        "data_file_division": "PLN",
        "date": _today(),
        # the source feature is the record's own statement about the sample; the name is what the
        # reader fell back to when a record carried none, so this keeps the round trip exact
        "organism": genome.source.get("organism") or genome.name,
        "source": genome.source.get("organism") or genome.name,
    }
    features: list[SeqFeature] = []
    source = _source_feature(genome)
    if source is not None:
        features.append(source)
    features.extend(_gene_feature(gene, length, translations) for gene in genome.genes)
    features.extend(_repeat_feature(interval, length) for interval in genome.repeat_regions)
    record.features = features
    return record


def _identifier(genome: Genome) -> str:
    """The record's identifier: its accession, or the first word of its name.

    A LOCUS line carries no whitespace, and an organism name ("Oldenlandia corymbosa var.
    corymbosa") is mostly whitespace — GenBank's own answer is a short label beside the full name
    in the ORGANISM block, which is exactly what this writes.
    """
    for candidate in (genome.accession, genome.name):
        words = (candidate or "").split()
        if words:
            return words[0]
    raise ValueError(
        "a GenBank record needs an identifier: this genome has no accession and no name"
    )


def _source_feature(genome: Genome) -> SeqFeature | None:
    """The ``source`` feature, from the qualifiers the record stated — nothing invented.

    A submission has to state organism · organelle · mol_type, and a record without them is a
    real state that ``submission_problems`` reports; filling them in here would hide it.
    """
    if not genome.source:
        return None
    return SeqFeature(
        SimpleLocation(0, len(genome), strand=1),
        type="source",
        qualifiers={key: [value] for key, value in genome.source.items()},
    )


def _gene_feature(gene: Gene, length: int, translations: bool = False) -> SeqFeature:
    qualifiers = dict(gene_qualifiers(gene))
    if translations and gene.translation:
        qualifiers["translation"] = gene.translation
    return SeqFeature(
        _location(gene.parts, length, gene.name),
        type=gene.kind,
        qualifiers={key: [value] for key, value in qualifiers.items()},
    )


def _repeat_feature(interval: Interval, length: int) -> SeqFeature:
    return SeqFeature(
        _location((interval,), length, "repeat_region"),
        type="repeat_region",
    )


def _location(
    parts: tuple[Interval, ...], length: int, what: str
) -> SimpleLocation | CompoundLocation:
    """One part is a plain location; several are a ``join()``, in the order given.

    The order is the annotation's, not a sorted one: a trans-spliced CDS is written in transcript
    order, and re-sorting it would describe a different protein (see ``model.coding_sequence``).
    """
    for part in parts:
        if part.end > length:
            raise ValueError(
                f"{what}: {part.start + 1}..{part.end} is past the end of a {length} bp sequence — "
                f"a feature that runs off the origin of a circular genome is a join, not a part "
                f"with an out-of-range coordinate"
            )
    locations = [
        SimpleLocation(
            BeforePosition(part.start) if part.partial_start else part.start,
            AfterPosition(part.end) if part.partial_end else part.end,
            strand=part.strand,
        )
        for part in parts
    ]
    return locations[0] if len(locations) == 1 else CompoundLocation(locations)


_MONTHS = ("JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC")


def _today() -> str:
    """GenBank's ``DD-MMM-YYYY``, spelled out so no machine locale can change it."""
    now = datetime.now()
    return f"{now.day:02d}-{_MONTHS[now.month - 1]}-{now.year}"


def _feature_name(feature: SeqFeature) -> str:
    for key in ("gene", "product", "locus_tag"):
        value = _first_qualifier(feature, key)
        if value:
            return value
    return feature.type


def _first_qualifier(feature: SeqFeature, key: str) -> str | None:
    values = feature.qualifiers.get(key)
    if not values:
        return None
    return str(values[0])
