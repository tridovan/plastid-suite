"""Locating the resource files that ship inside the package (``src/plastid/data``).

Three loaders need this directory (defaults, citations, gene groups). It was a copy-pasted
one-liner in each, and the idiom encodes an assumption — "this module sits one level below the
package root" — that breaks silently the day a loader moves. One definition, one place to fix.
"""

from __future__ import annotations

from pathlib import Path

__all__ = ["data_file"]

_DATA_DIR = Path(__file__).resolve().parents[1] / "data"


def data_file(name: str) -> Path:
    """Path to a packaged data file — ``data_file("citations.yaml")``."""
    return _DATA_DIR / name
