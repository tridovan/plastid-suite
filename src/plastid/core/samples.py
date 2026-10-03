"""Sample-group metadata — shared labels for every comparison figure and table.

The loader ships in Phase 1; the comparison modules that consume it arrive in Phase 2, so
ROADMAP M1.11 stays open until a real comparison table uses a group label.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import yaml

__all__ = ["SampleGroup", "group_of", "load_groups"]


@dataclass(frozen=True, slots=True)
class SampleGroup:
    name: str
    accessions: tuple[str, ...]


def load_groups(path: str | Path) -> tuple[SampleGroup, ...]:
    """Load ``{group_name: [accession, ...]}``; an accession may belong to one group only."""
    payload = yaml.safe_load(Path(path).read_text(encoding="utf-8")) or {}
    if not isinstance(payload, dict):
        raise ValueError(f"{path}: groups must be a mapping of group name to accessions")

    groups: list[SampleGroup] = []
    seen: dict[str, str] = {}
    for name, accessions in payload.items():
        if not isinstance(accessions, list) or not accessions:
            raise ValueError(f"group {name!r} must list at least one accession")
        for accession in accessions:
            if accession in seen:
                raise ValueError(
                    f"accession {accession!r} appears in both {seen[accession]!r} and {name!r}"
                )
            seen[accession] = name
        groups.append(SampleGroup(name=str(name), accessions=tuple(str(a) for a in accessions)))
    return tuple(groups)


def group_of(groups: tuple[SampleGroup, ...], accession: str) -> str | None:
    """The group label for an accession, or None when it belongs to no group."""
    for group in groups:
        if accession in group.accessions:
            return group.name
    return None
