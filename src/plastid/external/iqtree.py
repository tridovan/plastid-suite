"""IQ-TREE adapter — maximum-likelihood tree inference (M11.21; docs/tools/iqtree.md).

The CGAS route this module follows is MAFFT per gene → concatenated supermatrix → IQ-TREE
(`cgas_module14.py`, read 28/10/2026). This adapter owns only the last step: it takes an aligned
FASTA and returns the tree, with the tool's version and argv recorded like every other call.

**Why the version is read from a banner and checked against a pin.** `iqtree2 --version` prints
``IQ-TREE multicore version 2.3.6 for Linux x86 64-bit built Dec 11 2024``, so the release is
available, and a different release is a different measurement: model selection (ModelFinder) and
the bootstrap implementation have both changed between 2.x builds. A tree built by an unpinned
binary cannot be reproduced, so a release other than ``IQTREE_VERSION`` is refused by name rather
than recorded and used.

**The parameters are the published ones, stated once.** ``-m MFP`` (ModelFinder Plus),
``-bb 1000`` (ultrafast bootstrap), ``-alrt 1000`` (SH-aLRT) and ``-bnni`` are the settings the
CGAS module uses and the values IQ-TREE's own documentation recommends for a tree a paper reports.
They are module parameters, not adapter constants, so a run can change them and the manifest says
what it changed — but the defaults do not move.

**A seed and one thread, both on purpose.** ``-seed 1`` fixes the bootstrap's random stream, and
``-nt 1`` makes the run reproducible. Measured 28/10/2026 on the 4-record matrix: with ``-nt AUTO``
two runs of the same input gave **different topologies and support values** (99.5 vs 99.7 on an
internal node); with ``-nt 1`` three runs gave byte-identical support and one topology. IQ-TREE's
parallel reductions are not associative, so the default is **one thread** — a tree that changes
between two runs of the same matrix is not a result. ``--set threads=N`` asks for a parallel run
and the module records that choice, because it buys speed at the cost of bit-reproducibility;
the number that reaches a paper should be the single-threaded one.
"""

from __future__ import annotations

import re
import shutil
import tempfile
from dataclasses import dataclass
from pathlib import Path

from plastid.run.context import RunContext
from plastid.run.tools import ToolError, which

__all__ = ["IQTREE_VERSION", "IqtreeResult", "iqtree_version", "resolve_binary", "run_iqtree"]

IQTREE_VERSION = "2.3.6"
IQTREE_DOC = "docs/tools/iqtree.md"

# The executable's own name varies by packaging (iqtree2, iqtree, iqtree-omp).
_CANDIDATES = ("iqtree2", "iqtree", "iqtree3", "iqtree-omp", "iqtree2-omp")
_VERSION = re.compile(r"IQ-TREE[^(]*version\s+([0-9][0-9.]*)", re.IGNORECASE)


@dataclass(frozen=True, slots=True)
class IqtreeResult:
    """What IQ-TREE wrote: the tree, its log, and the version that produced them."""

    treefile: str
    log: str
    version: str
    report: str  # the `.iqtree` file: model chosen, support summary
    argv: tuple[str, ...]


def resolve_binary(directory: str | Path | None = None) -> str:
    """Find the IQ-TREE executable: a given directory, then ``PATH``."""
    if directory is not None:
        candidate_dir = Path(directory)
        for name in _CANDIDATES:
            candidate = candidate_dir / name
            if candidate.is_file():
                return str(candidate)
        raise ToolError(f"iqtree not found in {directory} — see {IQTREE_DOC}")
    for name in _CANDIDATES:
        found = which(name)
        if found is not None:
            return found
    raise ToolError(
        f"iqtree not found — tree inference needs the pinned IQ-TREE {IQTREE_VERSION}. "
        f"Install it (conda: `micromamba create -p <dir> -c bioconda iqtree="
        f"{IQTREE_VERSION}`), or point the run at one with --set iqtree=<directory>. "
        f"See {IQTREE_DOC}. The tree figure still draws a tree you supply yourself."
    )


def iqtree_version(binary: str) -> str:
    """The IQ-TREE release, read from ``--version`` and checked against the pin."""
    import subprocess

    try:
        result = subprocess.run(
            [binary, "--version"], capture_output=True, text=True, timeout=60, check=False
        )
    except OSError as error:
        raise ToolError(f"iqtree: cannot run {binary}: {error}") from error
    match = _VERSION.search((result.stdout or "") + (result.stderr or ""))
    if match is None:
        raise ToolError(
            f"{binary}: no 'IQ-TREE ... version x.y.z' line in its output — a step whose "
            f"version nobody can name is not part of the toolkit. Output was: "
            f"{((result.stdout or '') + (result.stderr or ''))[:200]!r}"
        )
    version = match.group(1)
    if version != IQTREE_VERSION:
        raise ToolError(
            f"{binary} is IQ-TREE {version}, but the toolkit's numbers were measured with "
            f"{IQTREE_VERSION} ({IQTREE_DOC}) — install the pinned build, or re-measure and "
            f"update IQTREE_VERSION."
        )
    return version


def run_iqtree(
    alignment_path: Path,
    ctx: RunContext,
    *,
    model: str = "MFP",
    ultrafast_bootstrap: int = 1000,
    alrt: int = 1000,
    bnni: bool = True,
    threads: str = "1",
    seed: int = 1,
    directory: str | Path | None = None,
    timeout: int | None = None,
    partition: str | None = None,
) -> IqtreeResult:
    """Run IQ-TREE on an aligned FASTA and return the tree it wrote.

    The input is copied into a throwaway directory so IQ-TREE's many side files
    (``.ckp.gz``, ``.bionj``, ``.contree``…) never land in the run directory; the files worth
    keeping are copied back by the caller from the returned strings.

    ``partition`` is a NEXUS ``begin sets`` block. When given it is written beside the alignment
    and passed as ``-p``, which asks IQ-TREE for an edge-proportional partition model with
    ModelFinder choosing a model per partition. ``None`` keeps the unpartitioned run.
    """
    binary = resolve_binary(directory)
    version = iqtree_version(binary)
    if not Path(alignment_path).is_file():
        raise ToolError(f"iqtree: no alignment at {alignment_path}")
    if ultrafast_bootstrap < 1000:
        # IQ-TREE itself warns below 1000; a tree a paper reports should not carry the warning.
        raise ToolError(
            f"iqtree: -bb {ultrafast_bootstrap} is below 1000 — ultrafast bootstrap support "
            f"below 1000 replicates is not a support value this toolkit will publish"
        )

    argv = [binary, "-s", "alignment.fasta", "--prefix", "tree", "-m", model]
    if partition is not None:
        argv += ["-p", "partition.nex"]
    argv += ["-bb", str(ultrafast_bootstrap), "-alrt", str(alrt)]
    if bnni:
        argv.append("-bnni")
    argv += ["-nt", str(threads), "-seed", str(seed), "--redo"]

    with tempfile.TemporaryDirectory(prefix="iqtree_") as work:
        workdir = Path(work)
        shutil.copyfile(alignment_path, workdir / "alignment.fasta")
        if partition is not None:
            (workdir / "partition.nex").write_text(partition, encoding="utf-8")
        ctx.tool(
            "iqtree2",
            argv,
            tool_path=binary,
            pinned_version=version,
            cwd=workdir,
            timeout=timeout,
        )
        tree = workdir / "tree.treefile"
        if not tree.is_file():
            raise ToolError("iqtree: exited 0 but wrote no tree.treefile — see the captured logs")
        report_path = workdir / "tree.iqtree"
        return IqtreeResult(
            treefile=tree.read_text(encoding="utf-8"),
            log=(workdir / "tree.log").read_text(encoding="utf-8")
            if (workdir / "tree.log").is_file()
            else "",
            report=report_path.read_text(encoding="utf-8") if report_path.is_file() else "",
            version=version,
            argv=tuple(argv),
        )
