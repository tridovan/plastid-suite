"""Render a run's R figures from the run's own tables (M10.3/M10.4).

The second figure pass, beside :mod:`plastid.run.figures`. It is separate because the guarantees
are separate — the matplotlib figures are byte-identical across runs and vector only, the R ones
are whatever the manuscript's scripts draw — and because this pass needs an external tool, which
``run/`` may not reach for itself (``tests/unit/test_architecture.py``). So the R call arrives as
an **injected** ``runner``, wired in ``cli.py`` exactly as ``plastid ir-junction`` wires CPJSdraw.

Like the matplotlib pass this reads each module's stored ``result.json``, so a run served entirely
from the cache still gets its figures, and drawing a picture never re-runs an analysis.

A figure whose source tables this run does not hold is **skipped by name**, never drawn empty —
and so is one the toolkit cannot draw at all: a species the script does not order by, a run whose
scope produced no Ts/Tv table. Those are skips with a reason a user can act on, printed and
recorded, never a traceback and never a silent blank panel.
"""

from __future__ import annotations

import hashlib
import json
import shutil
import tempfile
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from plastid.result.figure import FigureRecord
from plastid.result.figures import all_r_figures
from plastid.result.figures.r_tree import relabelled_tree
from plastid.result.result import Result
from plastid.result.rfigure import RFigureSkip
from plastid.run import batch_frames
from plastid.run.figures import read_run_results
from plastid.run.layout import FIGURES_DIRNAME
from plastid.run.tools import ToolError

__all__ = ["RFigureOutcome", "Preflight", "Runner", "render_run_rfigures"]

# (script key, --input file or None, --inputdir or None, out_dir) -> the script's Figures directory.
Runner = Callable[[str, "Path | None", "Path | None", Path], Path]
# None when R can draw; otherwise the reason to skip every R figure, asked once per pass.
Preflight = Callable[[], "str | None"]


@dataclass(frozen=True, slots=True)
class RFigureOutcome:
    """What the R pass wrote, and every figure it did not draw with the reason why."""

    files: tuple[Path, ...]
    skipped: tuple[tuple[str, str], ...]  # (figure id, reason)
    # The preflight's reason when R itself is unusable, else None. Carried out of the pass so the
    # command can tell "the run has no such table" (a data problem) from "R is not installed"
    # (a tool problem, exit 3) — the two must not share an exit code.
    blocked: str | None = None
    # One record per catalogue entry, in the same shape the matplotlib pass reports (M11.18).
    records: tuple[FigureRecord, ...] = field(default_factory=tuple)
    # The figures this run could feed and R did not draw, each also in `skipped`. Separate because
    # a table the run does not hold is expected, and a failed script is not.
    failed: tuple[tuple[str, str], ...] = ()


def render_run_rfigures(
    run_dir: str | Path,
    *,
    runner: Runner,
    preflight: Preflight | None = None,
    figures: Iterable[str] | None = None,
    params: Mapping[str, Any] | None = None,
    results: Mapping[str, Result] | None = None,
    provenance: Sequence[str] = (),
) -> RFigureOutcome:
    """Draw every R figure whose source tables this run holds, into ``<run_dir>/figures/``.

    ``params`` are the run's ``--set`` values. Only the tree figure reads one — it takes the path
    to a tree the user built — but every spec is offered them, so a later figure needing a window
    size or an aligner choice has somewhere to find it.
    """
    root = Path(run_dir)
    results = dict(results) if results is not None else read_run_results(root)
    target = root / FIGURES_DIRNAME
    wanted = None if figures is None else set(figures)
    settings: Mapping[str, Any] = params or {}
    names = _species_names(root)

    candidates = [spec for spec in all_r_figures() if wanted is None or spec.id in wanted]
    # Asked once for the whole pass, and only when there is something to draw: `--figure
    # circular_map` must not spawn an Rscript probe to discover that no R figure was wanted.
    blocked = preflight() if candidates and preflight is not None else None
    files: list[Path] = []
    skipped: list[tuple[str, str]] = []
    failed: list[tuple[str, str]] = []
    records: list[FigureRecord] = []
    for spec in candidates:
        if blocked is not None:
            _discard(target, spec.id)
            skipped.append((spec.id, blocked))
            records.append(FigureRecord(spec.id, "R", reason=blocked))
            continue
        missing = spec.missing(results)
        if missing:
            reason = f"needs {', '.join(missing)} — not in this run"
            _discard(target, spec.id)
            skipped.append((spec.id, reason))
            records.append(FigureRecord(spec.id, "R", reason=reason))
            continue
        source = _source_file(spec, settings, root)
        if source is _NO_SOURCE:
            reason = f"needs --set {spec.source_param}=<file>" + (
                f" — {settings[spec.source_param]!r} is not a file"
                if settings.get(spec.source_param)
                else ""
            )
            _discard(target, spec.id)
            skipped.append((spec.id, reason))
            records.append(FigureRecord(spec.id, "R", reason=reason))
            continue
        try:
            frame = batch_frames.frame_for(spec.id, results) or spec.frame(results)
        except RFigureSkip as reason:
            _discard(target, spec.id)
            skipped.append((spec.id, str(reason)))
            records.append(FigureRecord(spec.id, "R", reason=str(reason)))
            continue
        try:
            drawn, record_lines = _draw(
                spec, frame, target, runner, results, source, provenance, names
            )
        except ToolError as error:
            # R itself refused — a missing script, a package it does not have, a non-zero exit.
            # The adapter already named what was missing, so the figure is skipped with that name
            # rather than failing the pass: the other figures can still be drawn.
            _discard(target, spec.id)
            skipped.append((spec.id, str(error)))
            failed.append((spec.id, str(error)))
            records.append(FigureRecord(spec.id, "R", reason=str(error), failed=True))
            continue
        files.extend(drawn)
        records.append(FigureRecord(spec.id, "R", tuple(drawn), provenance=record_lines))
    return RFigureOutcome(
        files=tuple(files),
        skipped=tuple(skipped),
        blocked=blocked,
        records=tuple(records),
        failed=tuple(failed),
    )


def _discard(target: Path, figure_id: str) -> None:
    """Remove every file a previous pass left for one R figure — a refresh, not a stale pile.

    The R scripts decide their own output names, so the whole ``<id>__*`` space belongs to the
    figure, together with its caption and the pivot directory R read.
    """
    if not target.is_dir():
        return
    for path in target.glob(f"{figure_id}__*"):
        if path.is_dir():
            shutil.rmtree(path)
        else:
            path.unlink()
    caption = target / f"{figure_id}.caption.txt"
    if caption.is_file():
        caption.unlink()


# Distinct from `None`, which means "this spec needs no source file at all".
_NO_SOURCE = object()


def _source_file(spec, params: Mapping[str, Any], run_dir: Path) -> Path | None | object:
    """The file a spec reads directly: the one the user named, else the one this run wrote."""
    if spec.source_param is None:
        return None
    given = params.get(spec.source_param)
    if given:
        path = Path(str(given))
        return path if path.is_file() else _NO_SOURCE
    if spec.source_artifact is None:
        return _NO_SOURCE
    module, filename = spec.source_artifact
    written = sorted(run_dir.glob(f"*_{module}/{filename}"))
    return written[0] if written else _NO_SOURCE


def _draw(
    spec,
    frame,
    target: Path,
    runner: Runner,
    results,
    source: Path | None,
    run_provenance: Sequence[str] = (),
    names: Mapping[str, str] | None = None,
) -> tuple[list[Path], tuple[str, ...]]:
    """One figure: write what R reads, run the script, keep its images and its caption."""
    target.mkdir(parents=True, exist_ok=True)
    # A repeat render replaces the figure's whole output space; the script may draw a different
    # set of names on a different table, and a stale image would look current (M11.18).
    _discard(target, spec.id)
    labelled = getattr(frame, "labelled", None)
    if labelled is not None and names:
        frame = labelled(names)
    inputs = frame.inputs()
    # The pivot files stay in the run: they are the table R actually read, so a figure can be
    # reproduced — and checked — without re-deriving them. Named ``<id>__`` so two figures can
    # never overwrite one another's inputs. A spec reading a file of the user's own writes nothing.
    input_dir = target / f"{spec.id}__input"
    written: list[Path] = []
    if inputs:
        input_dir.mkdir(parents=True, exist_ok=True)
        for filename, text in inputs:
            path = input_dir / filename
            path.write_text(text, encoding="utf-8")
            written.append(path)

    own_tree = source is not None and names and source.is_relative_to(target.parent)
    drawn_from = source
    if own_tree:
        drawn_from = _own_tree_for_drawing(source, names, input_dir)

    produced: list[Path] = []
    with tempfile.TemporaryDirectory(prefix=f"pts-{spec.id}-") as scratch:
        if drawn_from is not None:
            figures_dir = runner(spec.script, drawn_from, None, Path(scratch))
        elif spec.input_mode == "file":
            figures_dir = runner(spec.script, written[0], None, Path(scratch))
        else:
            figures_dir = runner(spec.script, None, input_dir, Path(scratch))
        for path in sorted(figures_dir.iterdir()):
            if path.is_file() and path.suffix.lower() in (".png", ".pdf"):
                destination = target / f"{spec.id}__{path.name}"
                shutil.copyfile(path, destination)
                produced.append(destination)

    # The caption is not decoration: an exported number must carry its scope label, and an R
    # figure carries no room for one on its face (the scripts draw what they draw), so it is
    # written beside the image. Every label here is the one the *run* stored.
    caption = target / f"{spec.id}.caption.txt"
    lines = [
        spec.title,
        *run_provenance,
        frame.caption(),
        *_source_line(spec, source, target.parent),
        *spec.scope_lines(results),
    ]
    if names and (labelled is not None or own_tree):
        pairs = "; ".join(f"{accession} = {name}" for accession, name in names.items())
        lines.append(f"labels: {pairs}")
    lines.extend(spec.notes)
    caption.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return [*produced, caption], tuple(lines)


def _species_names(run_dir: Path) -> dict[str, str]:
    """Each input's accession and species, from the manifest. A batch names its own samples."""
    path = run_dir / "manifest.json"
    if not path.is_file():
        return {}
    inputs = json.loads(path.read_text(encoding="utf-8")).get("inputs") or []
    return {
        entry["accession"]: entry["organism"]
        for entry in inputs
        if entry.get("accession") and entry.get("organism")
    }


def _own_tree_for_drawing(source: Path, names: Mapping[str, str], input_dir: Path) -> Path:
    """A copy of the run's tree with species at its tips, kept beside what R read."""
    input_dir.mkdir(parents=True, exist_ok=True)
    copy = input_dir / source.name
    copy.write_text(relabelled_tree(source.read_text(encoding="utf-8"), names), encoding="utf-8")
    return copy


def _source_line(spec, source: Path | None, run_dir: Path) -> tuple[str, ...]:
    """Name the file that drew the figure, with its digest — the figure's only tie to it."""
    if source is None:
        return ()
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    if source.is_relative_to(run_dir):
        own = source.relative_to(run_dir).as_posix()
        return (f"source: this run's {own} (sha256 {digest})",)
    return (f"source: --set {spec.source_param}={source} (sha256 {digest})",)
