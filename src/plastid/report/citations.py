"""Citation registry — reports pull references instead of someone typing them per report."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import yaml

from plastid.core.resources import data_file

__all__ = ["Citation", "default_citations_path", "format_citation", "load_citations"]


@dataclass(frozen=True, slots=True)
class Citation:
    key: str
    text: str
    doi: str | None = None
    url: str | None = None


def default_citations_path() -> Path:
    return data_file("citations.yaml")


def load_citations(path: str | Path | None = None) -> dict[str, Citation]:
    """Load ``data/citations.yaml`` (or an explicit file) into key → citation."""
    source = Path(path) if path is not None else default_citations_path()
    payload = yaml.safe_load(source.read_text(encoding="utf-8")) or {}
    citations: dict[str, Citation] = {}
    for key, entry in payload.items():
        if not isinstance(entry, dict) or "text" not in entry:
            raise ValueError(f"citation {key!r} needs at least a 'text' field")
        citations[key] = Citation(
            key=key,
            text=str(entry["text"]),
            doi=str(entry["doi"]) if entry.get("doi") else None,
            url=str(entry["url"]) if entry.get("url") else None,
        )
    return citations


def format_citation(citation: Citation) -> str:
    """One line, doi appended when known."""
    line = citation.text.strip()
    if citation.doi:
        line = f"{line} doi:{citation.doi}"
    return line
