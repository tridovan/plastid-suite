"""Paper Figure 9 and S1–S3 — the maximum-likelihood tree, drawn by ``R_F9_phylogeny_tree.R``.

**This figure reads no table.** `requires` is empty and the pivot writes nothing. The tree is the
one this run's `phylogeny` step wrote, or the file named by ``--set tree=<file.treefile>``, which
wins when both exist. That tree is the figure's only input, and it is why the caption records the
file and its digest rather than a scope label — two figures drawn from two trees are two different
figures.

**No taxon is dropped and no outgroup is chosen by default.** The script's root/tips/outgroup
flags are reachable through the run's ``--set`` parameters, and PTS passes them to R (K-06). A
tree of another family is drawn as it is — midpoint-rooted unless the user names an outgroup —
rather than silently reproducing the manuscript's Orchidaceae/Asparagus choices, which are not
properties of this toolkit.
"""

from __future__ import annotations

import re
from collections.abc import Mapping
from dataclasses import dataclass

from plastid.result.figures.r_common import genome_labels
from plastid.result.result import Result
from plastid.result.rfigure import RFigureSpec

_TIP = re.compile(r"(?<=[(,])[^(),:;]+(?=[:,);])")
_NEWICK_RESERVED = re.compile(r"[\s(),:;\[\]']+")


@dataclass(frozen=True, slots=True)
class TreeFrame:
    """No file to write and no cell to read — the tree is already a file."""

    def inputs(self) -> tuple[tuple[str, str], ...]:
        return ()

    def caption(self) -> str:
        return (
            "The maximum-likelihood tree, read from the file named above. No table of this run "
            "is read here, so this caption carries no scope label: what the figure depends on "
            "is the tree file, named above with its digest. No tip is dropped "
            "and no outgroup is chosen unless the run passes --set drop-tips= and --set "
            "outgroup=; with no outgroup the tree is midpoint-rooted, and --set root=none keeps "
            "its own root. The script also writes <tag>_used.treefile — the tree as drawn, after "
            "the drops and the root — so the transformation can be checked rather than inferred "
            "from the image. The paper's choices are that paper's, not this toolkit's defaults."
        )


def relabelled_tree(newick: str, names: Mapping[str, str]) -> str:
    """The same tree with its tips under their species; lengths and support stay as they are.

    A tip follows `(` or `,`, which is what keeps a support value after `)` out of the match.
    """
    labels = genome_labels(tuple(_TIP.findall(newick)), names)

    def named(match: re.Match[str]) -> str:
        tip = match.group(0)
        return _NEWICK_RESERVED.sub("_", labels[tip]).strip("_") if tip in names else tip

    return _TIP.sub(named, newick)


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
    source_artifact=("phylogeny", "tree.treefile"),
    notes=(
        "A file named by --set tree= is the user's own: it is never written into the run, "
        "cached as a result, or altered.",
    ),
)
