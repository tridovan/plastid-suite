"""GenBank reader, backed by BioPython.

Never hand-roll a GenBank parser (CLAUDE.md, banned outright): BioPython already
handles ``join()``, ``complement()``, trans-splicing, and partial positions.
"""

from __future__ import annotations

from pathlib import Path

from Bio import SeqIO
from Bio.SeqFeature import (
    AfterPosition,
    BeforePosition,
    SeqFeature,
    SimpleLocation,
    UnknownPosition,
)
from Bio.SeqRecord import SeqRecord

from plastid.core.model import Gene, Genome, Interval

__all__ = ["read_genbank", "read_genbank_records"]

_GENE_FEATURE_TYPES = frozenset({"gene", "CDS", "tRNA", "rRNA"})


def read_genbank(path: str | Path) -> Genome:
    """Parse a single-record GenBank file into a :class:`Genome`."""
    records = read_genbank_records(path)
    if len(records) != 1:
        raise ValueError(f"{path}: expected exactly one record, found {len(records)}")
    return records[0]


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
    )


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
    )


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
