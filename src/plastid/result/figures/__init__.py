"""The figure catalogue: one module per figure, one line in :data:`_SPECS` (OCP).

Importing this package pulls in no plotting library — each figure imports matplotlib inside its
``build`` step — so the registry is safe to import on a minimum install.
"""

from __future__ import annotations

from plastid.result.figure import FigureSpec
from plastid.result.figures import (
    barcode_assessment,
    circular_map,
    genome_page,
    identity_synteny,
    ir_junction,
    pi_window,
    rscu_heatmap,
    ssr_distribution,
)

__all__ = ["all_figures", "get", "specs_for"]

_SPECS: tuple[FigureSpec, ...] = (
    pi_window.SPEC,
    circular_map.SPEC,
    ir_junction.SPEC,
    rscu_heatmap.SPEC,
    ssr_distribution.SPEC,
    identity_synteny.SPEC,
    barcode_assessment.SPEC,
    genome_page.SPEC,
)


def all_figures() -> tuple[FigureSpec, ...]:
    """Every figure type the toolkit knows, in catalogue order."""
    return _SPECS


def get(figure_id: str) -> FigureSpec:
    for spec in _SPECS:
        if spec.id == figure_id:
            return spec
    raise KeyError(f"unknown figure {figure_id!r} (have: {sorted(s.id for s in _SPECS)})")


def specs_for(modules: set[str]) -> tuple[FigureSpec, ...]:
    """The figures whose source modules all ran — the rest are skipped, and said so."""
    return tuple(spec for spec in _SPECS if {module for module, _ in spec.requires} <= modules)
