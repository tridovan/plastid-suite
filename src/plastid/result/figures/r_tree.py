"""Paper Figure 9 and S1–S3 — the maximum-likelihood tree, drawn by ``R_F9_phylogeny_tree.R``.

**This toolkit builds no tree, and this figure reads no table.** The tree is the user's own,
inferred with IQ-TREE (or anything else) outside this toolkit; `requires` is empty, the pivot
writes nothing, and the one thing the figure needs is the path to that tree, given as
``--set tree=<file.treefile>``. That is the bundle's own instruction for this figure
(``BANG_MAP_PTS_den_hinh.md``: "giữ `R_F9_phylogeny_tree.R`, gọi bằng `--set tree=<file.treefile>`
trên cây do người dùng dựng"), and it is why the caption records the file and its digest rather
than a scope label — two figures drawn from two trees are two different figures.

**The script's tip-dropping defaults are the manuscript's, not the toolkit's.** It drops two
*Asparagus* tips and roots on two Orchidaceae by default. On a tree of another family those names
are simply absent: nothing is dropped, the outgroup root falls through to a midpoint root, and the
figure comes out looking finished while describing a different rooting than the manuscript's. The
defaults are stated in the caption for that reason, and the flags are reachable through the run's
``--set`` parameters.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

from plastid.result.result import Result
from plastid.result.rfigure import RFigureSpec

# The script's own defaults, restated here so the caption can name them without reading the file.
_DEFAULT_DROP = "Asparagus officinalis, Asparagus schoberioides"
_DEFAULT_OUTGROUP = "Dendrobium officinale, Phalaenopsis equestris"


@dataclass(frozen=True, slots=True)
class TreeFrame:
    """No file to write and no cell to read — the tree is already a file the user owns."""

    def inputs(self) -> tuple[tuple[str, str], ...]:
        return ()

    def caption(self) -> str:
        return (
            "The maximum-likelihood tree as the user built it. This toolkit infers no tree and "
            "reads none of its own tables here, so this caption carries no scope label: what the "
            "figure depends on is the tree file, named above with its digest. The script drops "
            f"{_DEFAULT_DROP} and roots on {_DEFAULT_OUTGROUP} unless told otherwise — on a tree "
            "that does not carry those names it drops nothing and falls back to a midpoint root, "
            "so pass --drop-tips/--outgroup/--root when drawing another family."
        )


def frame(results: Mapping[str, Result]) -> TreeFrame:
    """Nothing is read from the run: the tree arrives as a file, not as a table."""
    return TreeFrame()


SPEC = RFigureSpec(
    id="r_tree",
    title="Maximum-likelihood tree (paper figure 9 + S1–S3, drawn by R)",
    requires=(),
    script="tree",
    input_mode="file",
    frame=frame,
    source_param="tree",
    notes=(
        "PTS builds no tree: the file named by --set tree= is the user's own, and is never "
        "written into the run, cached as a result, or altered.",
    ),
)
