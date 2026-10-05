"""Render a run's figures from the run's own tables (P5-01).

Figures are a **run-level** step, not a module step: two of them read tables from two different
modules, and a module never sees another module's :class:`Result`. Rendering happens *after* the
module loop and reads each module's stored ``result.json``, so

* a run that was served entirely from the cache still gets its figures — the stored result is the
  same one a fresh run would have produced;
* drawing a picture never re-runs an analysis (the phase's own in-phase check);
* ``plastid run --figures`` and re-rendering an existing run directory are the same function.

A figure whose source tables this run does not hold is **skipped by name**, never drawn empty.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path

from plastid.result.export import read_result
from plastid.result.figure import FIGURE_FORMATS, Figure
from plastid.result.figures import all_figures
from plastid.result.result import Result
from plastid.run.layout import FIGURES_DIRNAME

__all__ = ["FigureOutcome", "read_run_results", "render_run_figures"]


@dataclass(frozen=True, slots=True)
class FigureOutcome:
    """What a figure pass wrote, and every figure it did not draw with the reason why."""

    files: tuple[Path, ...]
    skipped: tuple[tuple[str, str], ...]  # (figure id, reason)


def read_run_results(run_dir: str | Path) -> dict[str, Result]:
    """Every module's stored result, keyed by module id — a run directory read back."""
    root = Path(run_dir)
    results: dict[str, Result] = {}
    for result_path in sorted(root.glob("*/result.json")):
        result = read_result(result_path.read_text(encoding="utf-8"))
        results[result.module] = result
    if not results:
        raise ValueError(
            f"no module result found under {root} — is this a run directory? "
            "(expected <module>/result.json beside the module's tables)"
        )
    return results


def render_run_figures(
    run_dir: str | Path,
    *,
    figures: Iterable[str] | None = None,
    formats: Sequence[str] = FIGURE_FORMATS,
) -> FigureOutcome:
    """Draw every figure whose source tables the run holds, into ``<run_dir>/figures/``."""
    root = Path(run_dir)
    results = read_run_results(root)
    wanted = None if figures is None else set(figures)
    target = root / FIGURES_DIRNAME
    files: list[Path] = []
    skipped: list[tuple[str, str]] = []
    for spec in all_figures():
        if wanted is not None and spec.id not in wanted:
            continue
        missing = spec.missing(results)
        if missing:
            skipped.append((spec.id, f"needs {', '.join(missing)} — not in this run"))
            continue
        figure: Figure = spec.build(results)
        target.mkdir(parents=True, exist_ok=True)
        for fmt in formats:
            path = target / f"{figure.id}.{fmt}"
            path.write_bytes(figure.payload(fmt))
            files.append(path)
    return FigureOutcome(files=tuple(files), skipped=tuple(skipped))
