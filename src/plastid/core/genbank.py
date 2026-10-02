"""GenBank reader, backed by BioPython.

Never hand-roll a GenBank parser (CLAUDE.md, banned outright): BioPython already
handles ``join()``, ``complement()``, trans-splicing, and partial positions.
"""

from __future__ import annotations

from pathlib import Path

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, SimpleLocation

from plastid.core.model import Gene, Genome, Interval

__all__ = ["read_genbank"]

_GENE_FEATURE_TYPES = frozenset({"gene", "CDS", "tRNA", "rRNA"})


def read_genbank(path: str | Path) -> Genome:
    """Parse a single-record GenBank file into a :class:`Genome`."""
    records = list(SeqIO.parse(str(path), "genbank"))
    if len(records) != 1:
        raise ValueError(f"{path}: expected exactly one record, found {len(records)}")

    record = records[0]
    genes = tuple(
        _to_gene(feature) for feature in record.features if feature.type in _GENE_FEATURE_TYPES
    )
    return Genome(
        name=str(record.annotations.get("organism", record.id)),
        accession=record.id,
        sequence=str(record.seq).upper(),
        topology=str(record.annotations.get("topology", "circular")).lower(),
        genes=genes,
    )


def _to_gene(feature: SeqFeature) -> Gene:
    locations = getattr(feature.location, "parts", [feature.location])
    return Gene(
        name=_feature_name(feature),
        kind=feature.type,
        parts=tuple(_to_interval(location) for location in locations),
        product=_first_qualifier(feature, "product"),
        pseudo="pseudo" in feature.qualifiers,
    )


def _to_interval(location: SimpleLocation) -> Interval:
    # BioPython positions are already 0-based half-open.
    strand = location.strand
    return Interval(
        start=int(location.start),
        end=int(location.end),
        strand=strand if strand in (1, -1) else 1,
    )


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
