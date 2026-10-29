"""Paper Figure 9 and S1–S3 — the maximum-likelihood tree, drawn by ``R_F9_phylogeny_tree.R``.

**This toolkit builds no tree, and this figure reads no table.** The tree is the user's own,
inferred with IQ-TREE (or anything else) outside this toolkit; `requires` is empty, the pivot
writes nothing, and the one thing the figure needs is the path to that tree, given as
``--set tree=<file.treefile>``. That tree is the figure's only input, and it is why the caption
records the file and its digest rather than a scope label — two figures drawn from two trees are two
different figures.

**No taxon is dropped and no outgroup is chosen by default.** The script's root/tips/outgroup
flags are reachable through the run's ``--set`` parameters, and PTS passes them to R (K-06). A
tree of another family is drawn as it is — midpoint-rooted unless the user names an outgroup —
rather than silently reproducing the manuscript's Orchidaceae/Asparagus choices, which are not
properties of this toolkit.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

from plastid.result.result import Result
from plastid.result.rfigure import RFigureSpec


@dataclass(frozen=True, slots=True)
class TreeFrame:
    """No file to write and no cell to read — the tree is already a file the user owns."""

    def inputs(self) -> tuple[tuple[str, str], ...]:
        return ()

    def caption(self) -> str:
        return (
            "The maximum-likelihood tree as the user built it. This toolkit infers no tree and "
            "reads none of its own tables here, so this caption carries no scope label: what the "
            "figure depends on is the tree file, named above with its digest. No tip is dropped "
            "and no outgroup is chosen unless the run passes --set drop-tips= and --set "
            "outgroup=; with no outgroup the tree is midpoint-rooted, and --set root=none keeps "
            "its own root. The script also writes <tag>_used.treefile — the tree as drawn, after "
            "the drops and the root — so the transformation can be checked rather than inferred "
            "from the image. The paper's choices are that paper's, not this toolkit's defaults."
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
