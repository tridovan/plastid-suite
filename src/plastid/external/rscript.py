"""R adapter — the manuscript's figure scripts (docs/tools/r.md).

**R is a plotting engine here and nothing else.** Each of these scripts reads a table the toolkit
already exported and writes a PNG and a PDF. None of them computes a number the toolkit does not
hold in a cell, and that is the rule the whole figure layer is built on: if a script and its table
disagree, the table is right.

**The scripts are package data** (``data/r_scripts/``), the same way ``data/misa.ini`` is: they are
the team's own code, carried so a figure can be reproduced, and read through
:func:`plastid.core.resources.data_file` rather than from a path beside the caller. They are kept in
the language their authors wrote them in — they are the manuscript's scripts, not this package's
code.

**The one change made to them.** ``khoi_tao()`` ends with ``setwd(OUT_DIR)``, so a *relative*
``--inputdir``/``--input`` points into the output directory rather than at the caller's data, and
every script fails to find its own input. Measured 17/10/2026: ``--inputdir data_mau`` fails, an
absolute path works. The fix is in ``00_ham_chung.R`` — the paths are made absolute before the
working directory moves — and it changes no drawing.

**A missing package is named, never guessed at.** R 4.1.2 does not have ``readr``, ``cowplot`` or
``gridExtra``, so the scripts were written without them; the nine packages they do need are checked
by name before any figure is attempted.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

from plastid.core.resources import data_file
from plastid.run.context import RunContext
from plastid.run.tools import ToolError, which

__all__ = [
    "R_PACKAGES",
    "R_SCRIPTS",
    "RSCRIPT_DOC",
    "figure_script",
    "missing_packages",
    "run_figure",
    "rscript_binary",
]

RSCRIPT_DOC = "docs/tools/r.md"

# The nine the scripts actually call, and the three they deliberately do not (readr, cowplot,
# gridExtra) — a script that needs one of those is a script that fails on the target machine.
R_PACKAGES: tuple[str, ...] = (
    "ggplot2",
    "dplyr",
    "tidyr",
    "reshape2",
    "RColorBrewer",
    "zoo",
    "scales",
    "pheatmap",
    "ape",
)

# The manuscript's figures, by the name this toolkit calls them. The tree (F9) is here as well
# because it is the same engine — but PTS builds no tree, so its input is the user's own.
R_SCRIPTS: dict[str, str] = {
    "rscu": "R_F3_RSCU_codon_usage.R",
    "amino_acid": "R_F4_amino_acid_heatmap.R",
    "ssr": "R_F5_SSR_distribution.R",
    "diversity": "R_F6_nucleotide_diversity_pi.R",
    "substitutions": "R_F7_pairwise_substitution.R",
    "kaks": "R_F8_KaKs_heatmap.R",
    "tree": "R_F9_phylogeny_tree.R",
}

_PACKAGE_PROBE = 'cat(setdiff(c({wanted}), rownames(installed.packages())), sep="\\n")'


def rscript_binary(explicit: str | Path | None = None) -> str:
    """The ``Rscript`` to run: an explicit path, else the name on ``PATH``."""
    if explicit is not None:
        path = Path(explicit)
        if not path.is_file():
            raise ToolError(f"rscript: --set rscript={explicit} is not a file — {RSCRIPT_DOC}")
        return str(path)
    found = which("Rscript")
    if found is None:
        raise ToolError(
            f"rscript: Rscript not found — R is the plotting engine behind the manuscript's "
            f"figures and is not bundled with this toolkit. Install it, or pass "
            f"--set rscript=/path/to/Rscript ({RSCRIPT_DOC})"
        )
    return found


def figure_script(name: str) -> Path:
    """The packaged script for a figure id — read from the package, never from a caller's path."""
    if name not in R_SCRIPTS:
        raise ToolError(f"rscript: no figure {name!r} (have: {', '.join(sorted(R_SCRIPTS))})")
    path = data_file(f"r_scripts/{R_SCRIPTS[name]}")
    if not path.is_file():
        raise ToolError(
            f"rscript: {R_SCRIPTS[name]} is missing from the package — reinstall it ({RSCRIPT_DOC})"
        )
    return path


def missing_packages(
    ctx: RunContext, binary: str, wanted: Sequence[str] = R_PACKAGES
) -> tuple[str, ...]:
    """The R packages this machine does not have, **by name**, asked once before anything runs.

    The probe goes through ``ctx.tool`` like every other external call, so the R version that
    answered is in the run's provenance rather than in nobody's memory.
    """
    code = _PACKAGE_PROBE.format(wanted=", ".join(f'"{name}"' for name in wanted))
    record = ctx.tool("Rscript", [binary, "-e", code], kind="binary", tool_path=binary)
    stdout = Path(record.stdout_file).read_text(encoding="utf-8")
    return tuple(line.strip() for line in stdout.splitlines() if line.strip())


def run_figure(
    ctx: RunContext,
    name: str,
    *,
    input_file: Path | None = None,
    input_dir: Path | None = None,
    out_dir: Path,
    extra: Sequence[str] = (),
    binary: str | Path | None = None,
    timeout: int | None = None,
    check_packages: bool = True,
) -> Path:
    """Run one figure script and return the directory its images landed in.

    ``input_dir`` and ``input_file`` are resolved to absolute paths here as well as inside the
    script: the script's own guard is the fix, and this is the same rule stated where a caller can
    read it.

    ``check_packages=False`` is for a caller that has **already** asked which packages are missing
    and is drawing several figures in a row: ``plastid figures`` probes once for the whole pass
    rather than once per figure, so the run's provenance holds one version record instead of seven.
    The default stays the safe one — a caller that has not asked gets the check.
    """
    executable = rscript_binary(binary)
    if check_packages:
        missing = missing_packages(ctx, executable)
        if missing:
            raise ToolError(
                f"rscript: R is missing {', '.join(missing)} — the figure scripts need all of "
                f"{', '.join(R_PACKAGES)}; install them ({RSCRIPT_DOC})"
            )

    script = figure_script(name)
    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    argv = [executable, str(script), "--outdir", str(out_dir)]
    if input_file is not None:
        argv += ["--input", str(Path(input_file).resolve())]
    if input_dir is not None:
        argv += ["--inputdir", str(Path(input_dir).resolve())]
    argv += [str(item) for item in extra]

    ctx.tool(
        f"rscript:{name}",
        argv,
        kind="binary",
        tool_path=executable,
        cwd=out_dir,
        timeout=timeout,
    )
    figures = out_dir / "Figures"
    if not figures.is_dir() or not any(figures.iterdir()):
        raise ToolError(
            f"rscript: {R_SCRIPTS[name]} exited cleanly but wrote no figure into {figures} — "
            f"see the captured logs"
        )
    return figures
