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
from dataclasses import dataclass, field
from pathlib import Path

from plastid.result.export import read_result
from plastid.result.figure import FIGURE_FORMATS, Figure, FigureRecord
from plastid.result.figures import all_figures
from plastid.result.result import Result
from plastid.run.layout import FIGURES_DIRNAME

__all__ = [
    "FigureOutcome",
    "read_run_results",
    "render_run_figures",
    "stale_figure_files",
]


@dataclass(frozen=True, slots=True)
class FigureOutcome:
    """What a figure pass wrote, and every figure it did not draw with the reason why."""

    files: tuple[Path, ...]
    skipped: tuple[tuple[str, str], ...]  # (figure id, reason)
    # One record per catalogue entry the pass considered — both renderers report through the
    # same shape so package/QC can list a run's whole catalogue (M11.18).
    records: tuple[FigureRecord, ...] = field(default_factory=tuple)


def _scope_lines(spec, results: dict[str, Result]) -> tuple[str, ...]:
    """``module.table — scope: <label>`` for every table a matplotlib figure read."""
    return tuple(
        f"{module}.{table} — scope: {results[module].tables[table].scope.label}"
        for module, table in spec.requires
        if module in results and table in results[module].tables
    )


def stale_figure_files(target: Path, figure_id: str, formats: Sequence[str]) -> tuple[Path, ...]:
    """The files a previous pass left for one figure, so a refresh cannot keep them.

    A re-render replaces what it draws; a figure that can no longer be drawn must not leave a
    stale image looking like a current result (M11.18).
    """
    if not target.is_dir():
        return ()
    return tuple(
        path for path in (target / f"{figure_id}.{fmt}" for fmt in formats) if path.is_file()
    )


def _discard(target: Path, figure_id: str, formats: Sequence[str]) -> None:
    for path in stale_figure_files(target, figure_id, formats):
        path.unlink()
    caption = target / f"{figure_id}.caption.txt"
    if caption.is_file():
        caption.unlink()


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
    results: dict[str, Result] | None = None,
    provenance: Sequence[str] = (),
) -> FigureOutcome:
    """Draw every figure whose source tables the run holds, into ``<run_dir>/figures/``."""
    root = Path(run_dir)
    results = read_run_results(root) if results is None else results
    wanted = None if figures is None else set(figures)
    target = root / FIGURES_DIRNAME
    files: list[Path] = []
    skipped: list[tuple[str, str]] = []
    records: list[FigureRecord] = []
    for spec in all_figures():
        if wanted is not None and spec.id not in wanted:
            continue
        missing = spec.missing(results)
        if missing:
            reason = f"needs {', '.join(missing)} — not in this run"
            _discard(target, spec.id, spec.formats)
            skipped.append((spec.id, reason))
            records.append(FigureRecord(spec.id, "matplotlib", reason=reason))
            continue
        wanted_formats = [fmt for fmt in formats if fmt in spec.formats]
        if not wanted_formats:
            reason = f"carries {'/'.join(spec.formats)} only, none of {'/'.join(formats)}"
            _discard(target, spec.id, spec.formats)
            skipped.append((spec.id, reason))
            records.append(FigureRecord(spec.id, "matplotlib", reason=reason))
            continue
        figure: Figure = spec.build(results)
        target.mkdir(parents=True, exist_ok=True)
        # Replace the whole set this figure owns: a format that is no longer asked for must not
        # survive from the earlier pass as if it were current.
        _discard(target, spec.id, spec.formats)
        drawn: list[Path] = []
        for fmt in wanted_formats:
            path = target / f"{figure.id}.{fmt}"
            path.write_bytes(figure.payload(fmt))
            files.append(path)
            drawn.append(path)
        # The caption is the image's provenance: the run's own label and the scope of every table
        # the figure read. It is written beside the image for the same reason the R pass writes
        # one — a figure must not leave the run without saying what it measured (M11.18).
        caption_lines = [figure.title, *provenance, *_scope_lines(spec, results)]
        caption = target / f"{figure.id}.caption.txt"
        caption.write_text("\n".join(caption_lines) + "\n", encoding="utf-8")
        files.append(caption)
        drawn.append(caption)
        records.append(
            FigureRecord(
                spec.id,
                "matplotlib",
                tuple(drawn),
                provenance=tuple(caption_lines),
            )
        )
    return FigureOutcome(files=tuple(files), skipped=tuple(skipped), records=tuple(records))
