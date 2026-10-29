"""Table cells a figure reads as coordinates rather than recomputing (D5-03).

The one cell here is ``genome_stats.genes.segments``: every part of a feature as
``start-end`` in 1-based GenBank coordinates, semicolon-separated. It exists because a feature
annotated across the record's origin has no single linear interval — its outer bounds are
``1..L``, which on a circle is the *long* way round. Two figures draw that table, so the parsing
lives in one place rather than once per figure.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

__all__ = ["gene_segments"]


def gene_segments(row: Mapping[str, Any]) -> tuple[tuple[int, int], ...]:
    """``"83184-83190;1-704"`` → two 1-based inclusive spans.

    Refused rather than guessed when the cell is absent. Falling back to ``start_1b``/``end_1b``
    would draw the feature as a whole-genome arc — the wrong picture this cell exists to prevent
    (M11.7) — so the message says how to get a table that has it.
    """
    raw = row.get("segments")
    if not raw:
        raise ValueError(
            f"the genes table has no `segments` cell for {row.get('name')!r} — re-run "
            f"genome_stats with this version of the toolkit so the table records every part"
        )
    spans: list[tuple[int, int]] = []
    for token in str(raw).split(";"):
        start_text, separator, end_text = token.partition("-")
        if not separator or not start_text.isdigit() or not end_text.isdigit():
            raise ValueError(
                f"malformed segments cell {raw!r} for {row.get('name')!r} — expected "
                f"'start-end' pairs separated by ';'"
            )
        spans.append((int(start_text), int(end_text)))
    return tuple(spans)
