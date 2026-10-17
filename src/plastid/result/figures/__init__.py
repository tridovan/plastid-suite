"""The figure catalogue: one module per figure, one line in :data:`_SPECS` (OCP).

Importing this package pulls in no plotting library — each figure imports matplotlib inside its
``build`` step — so the registry is safe to import on a minimum install.

**Two catalogues, deliberately.** :data:`_SPECS` holds the figures this toolkit **draws** with
matplotlib: vector only, byte-identical across runs. :data:`_R_SPECS` holds the manuscript's own
figures, drawn by its R scripts from the same tables — a second source for a figure that already
exists, and the only source for the three the toolkit has no drawing for. They are separate because
the guarantees are separate: see ``result/rfigure.py``. Importing either catalogue still pulls in no
plotting library and no R.
"""

from __future__ import annotations

from plastid.result.figure import FigureSpec
from plastid.result.figures import (
    barcode_assessment,
    circular_map,
    dotplot,
    genome_page,
    identity_synteny,
    ir_junction,
    pi_window,
    r_amino_acid,
    r_diversity,
    r_kaks,
    r_rscu,
    r_ssr,
    r_substitutions,
    r_tree,
    rscu_heatmap,
    ssr_distribution,
)
from plastid.result.rfigure import RFigureSpec

__all__ = ["all_figures", "all_r_figures", "get", "get_r", "specs_for"]

_SPECS: tuple[FigureSpec, ...] = (
    pi_window.SPEC,
    circular_map.SPEC,
    ir_junction.SPEC,
    rscu_heatmap.SPEC,
    ssr_distribution.SPEC,
    identity_synteny.SPEC,
    barcode_assessment.SPEC,
    genome_page.SPEC,
    dotplot.SPEC,
)

# The manuscript's figures, by the paper's own numbering. Each reads the same tables the
# matplotlib figure beside it reads, so the two must agree on every number.
_R_SPECS: tuple[RFigureSpec, ...] = (
    r_rscu.SPEC,  # paper figure 3
    r_amino_acid.SPEC,  # paper figure 4
    r_ssr.SPEC,  # paper figure 5
    r_diversity.SPEC,  # paper figure 6
    r_substitutions.SPEC,  # paper figure 7
    r_kaks.SPEC,  # paper figure 8
    r_tree.SPEC,  # paper figure 9 + S1-S3 -- the user's own tree, no table
)


def all_figures() -> tuple[FigureSpec, ...]:
    """Every figure type the toolkit knows, in catalogue order."""
    return _SPECS


def all_r_figures() -> tuple[RFigureSpec, ...]:
    """Every R-drawn figure, in the manuscript's own order — the second catalogue."""
    return _R_SPECS


def get_r(figure_id: str) -> RFigureSpec:
    for spec in _R_SPECS:
        if spec.id == figure_id:
            return spec
    raise KeyError(f"unknown R figure {figure_id!r} (have: {sorted(s.id for s in _R_SPECS)})")


def get(figure_id: str) -> FigureSpec:
    for spec in _SPECS:
        if spec.id == figure_id:
            return spec
    raise KeyError(f"unknown figure {figure_id!r} (have: {sorted(s.id for s in _SPECS)})")


def specs_for(modules: set[str]) -> tuple[FigureSpec, ...]:
    """The figures whose source modules all ran — the rest are skipped, and said so."""
    return tuple(spec for spec in _SPECS if {module for module, _ in spec.requires} <= modules)
