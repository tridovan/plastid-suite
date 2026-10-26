"""The figure system: series read cells, one theme, and byte-identical vector output (P5-01).

The data tests run without matplotlib — the series is plain data, and CLAUDE.md's testing rule
says to test the data fed in, not the pixels. The drawing tests skip by name when the optional
``figures`` extra is not installed, exactly as the Word-export tests do.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

from plastid.core.deps import ExportDependencyError
from plastid.result.export import export_json
from plastid.result.figure import FIGURE_FORMATS, Figure, FigureSpec, export_figure
from plastid.result.figures import all_figures, get, pi_window
from plastid.result.result import Column, Result, Scope, Table
from plastid.run.figures import read_run_results, render_run_figures
from plastid.run.layout import RunLayout

requires_figures = pytest.mark.skipif(
    importlib.util.find_spec("matplotlib") is None,
    reason="matplotlib is an optional dependency (pip install 'plastid-suite[figures]')",
)


def _windows_table() -> Table:
    return Table(
        id="diversity_windows",
        title="Diversity windows",
        columns=(
            Column("window_start_1b", "Start", dtype="int", unit="bp"),
            Column("window_end_1b", "End", dtype="int", unit="bp"),
            Column("usable_columns", "Usable columns", dtype="int"),
            Column("pi", "π", dtype="float", precision=6),
            Column("region", "Region"),
        ),
        rows=(
            {
                "window_start_1b": 1,
                "window_end_1b": 600,
                "usable_columns": 600,
                "pi": 0.001,
                "region": "LSC",
            },
            {
                "window_start_1b": 601,
                "window_end_1b": 1200,
                "usable_columns": 588,
                "pi": None,
                "region": "LSC",
            },
            {
                "window_start_1b": 1201,
                "window_end_1b": 1800,
                "usable_columns": 600,
                "pi": 0.02,
                "region": "IRa",
            },
            {
                "window_start_1b": 1801,
                "window_end_1b": 2400,
                "usable_columns": 601,
                "pi": 0.03,
                "region": "IRb",
            },
            {
                "window_start_1b": 2401,
                "window_end_1b": 3000,
                "usable_columns": 0,
                "pi": None,
                "region": "",
            },
        ),
        scope=Scope(label="aligned genomes; gapped columns dropped", coord_system="normalized"),
    )


def _diversity_result() -> Result:
    summary = Table(
        id="diversity_summary",
        title="Diversity",
        columns=(
            Column("accession", "Accession"),
            Column("organism", "Organism"),
            Column("n_genomes", "Genomes", dtype="int"),
            Column("alignment_length", "Alignment", dtype="int", unit="bp"),
            Column("usable_columns", "Usable columns", dtype="int"),
            Column("pi", "π", dtype="float", precision=6),
            Column("window_bp", "Window", dtype="int", unit="bp"),
            Column("step_bp", "Step", dtype="int", unit="bp"),
            Column("tool", "Tool"),
        ),
        rows=(
            {
                "accession": "NC_050833.1",
                "organism": "Iris domestica",
                "n_genomes": 2,
                "alignment_length": 3000,
                "usable_columns": 2389,
                "pi": 0.007181,
                "window_bp": 600,
                "step_bp": 200,
                "tool": "MAFFT",
            },
        ),
        scope=Scope(label="aligned genomes; gapped columns dropped", coord_system="normalized"),
    )
    return Result(
        module="diversity",
        module_version="0.1.0",
        tables={"diversity_windows": _windows_table(), "diversity_summary": summary},
    )


def test_the_series_equals_the_cells_of_the_table_it_read() -> None:
    """D5-03: a figure reads table cells — it never recomputes a quantity a module produced."""
    result = _diversity_result()
    data = pi_window.series({"diversity": result})
    table = result.tables["diversity_windows"]
    assert len(data.windows) == len(table.rows)
    for window, row in zip(data.windows, table.rows, strict=True):
        assert (window.start, window.end, window.pi, window.region) == (
            row["window_start_1b"],
            row["window_end_1b"],
            row["pi"],
            row["region"],
        )
    summary = result.tables["diversity_summary"]
    assert data.pi_whole == summary.format("pi")
    assert (data.window_bp, data.step_bp) == (summary.value("window_bp"), summary.value("step_bp"))
    assert data.alignment_length == summary.value("alignment_length")


def test_a_window_without_a_usable_column_stays_a_gap_never_a_zero() -> None:
    data = pi_window.series({"diversity": _diversity_result()})
    blank = [window for window in data.windows if window.pi is None]
    assert len(blank) == 2
    assert data.n_windows_without_value == 2
    assert all(window.pi is not None for window in data.windows if window.region == "IRa")


def test_region_bands_meet_at_a_boundary_and_skip_the_unassigned_windows() -> None:
    """Overlapping windows would paint one region over another; the bands share an edge instead."""
    data = pi_window.series({"diversity": _diversity_result()})
    bands = data.bands
    assert bands == (("LSC", 1, 1200), ("IRa", 1200, 1800), ("IRb", 1800, 2400))
    for (_, _, end), (_, start, _) in zip(bands, bands[1:], strict=False):
        assert end == start, "a gap or an overlap between two neighbouring region bands"


def test_a_band_boundary_lands_between_the_two_windows_that_disagree() -> None:
    """The measured reason bands are not taken from the windows' own edges (600 bp / 200 bp)."""
    table = _windows_table()
    overlapping = Table(
        id=table.id,
        title=table.title,
        columns=table.columns,
        rows=(
            {
                "window_start_1b": 1,
                "window_end_1b": 600,
                "usable_columns": 600,
                "pi": 0.001,
                "region": "LSC",
            },
            {
                "window_start_1b": 401,
                "window_end_1b": 1000,
                "usable_columns": 600,
                "pi": 0.002,
                "region": "IRa",
            },
        ),
        scope=table.scope,
    )
    result = _diversity_result()
    result = Result(
        module=result.module,
        module_version=result.module_version,
        tables={**result.tables, "diversity_windows": overlapping},
    )
    # The two windows overlap by 200 bp; the boundary sits at the midpoint of their centres
    # ((300.5 + 700.5) / 2), not at either window's edge.
    bands = pi_window.series({"diversity": result}).bands
    assert bands == (("LSC", 1, 500), ("IRa", 500, 1000))


def test_the_figure_says_which_tables_it_read() -> None:
    """The audit trail D5-03 asks for: what a figure drew must be nameable."""
    assert get("pi_window").requires == (
        ("diversity", "diversity_windows"),
        ("diversity", "diversity_summary"),
    )


def test_a_spec_keys_its_tables_by_module_not_by_table_name() -> None:
    """`regions` exists in two modules with different columns; a bare table id is ambiguous."""
    spec = FigureSpec(
        id="probe", title="Probe", requires=(("structure", "regions"),), build=lambda _: None
    )
    other = Result(
        module="genome_stats",
        module_version="0.1.0",
        tables={
            "regions": Table(
                id="regions",
                title="Regions",
                columns=(Column("region", "Region"),),
                rows=({"region": "LSC"},),
                scope=Scope(label="normalized"),
                key="region",
            )
        },
    )
    assert spec.missing({"genome_stats": other}) == ("structure.regions",)


def test_a_missing_source_table_is_named_rather_than_drawn_empty() -> None:
    with pytest.raises(ValueError, match="diversity"):
        pi_window.series({})


def test_the_caption_carries_the_summary_cell_rendering_and_the_gap_count() -> None:
    figure = pi_window.build({"diversity": _diversity_result()}) if _has_matplotlib() else None
    if figure is None:
        pytest.skip("matplotlib is an optional dependency (pip install 'plastid-suite[figures]')")
    assert "0.007181" in figure.title  # the cell's own precision, not a second rounding
    assert "2 window(s) without a usable column" in figure.title


def test_a_figure_refuses_an_unknown_format_and_an_empty_payload() -> None:
    with pytest.raises(ValueError, match="png"):
        Figure(id="x", title="X", source_tables=(), data=(("png", b"x"),))
    with pytest.raises(ValueError, match="empty"):
        Figure(id="x", title="X", source_tables=(), data=(("svg", b""),))
    with pytest.raises(ValueError, match="plain name"):
        Figure(id="a.b", title="X", source_tables=(), data=(("svg", b"x"),))


def test_a_missing_figures_extra_is_a_named_error(monkeypatch: pytest.MonkeyPatch) -> None:
    """The minimum install must stay pure Python — and say so, not raise ImportError."""
    monkeypatch.setitem(sys.modules, "matplotlib.pyplot", None)
    with pytest.raises(ExportDependencyError, match="plastid-suite\\[figures\\]"):
        export_figure(
            figure_id="probe",
            title="Probe",
            source_tables=(),
            build=lambda axes: axes.plot([0, 1], [0, 1]),
            figsize=(2, 2),
        )


def test_importing_the_figure_catalogue_never_pulls_in_matplotlib() -> None:
    """A minimum install must be able to register the catalogue without plotting code loaded."""
    import subprocess  # noqa: PLC0415 — a fresh interpreter is the only honest way to ask

    code = (
        "import sys, plastid.result.figures, plastid.run.figures;"
        "print('matplotlib' in sys.modules or 'matplotlib.pyplot' in sys.modules)"
    )
    done = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True)
    assert done.stdout.strip() == "False"


def _has_matplotlib() -> bool:
    return importlib.util.find_spec("matplotlib") is not None


# ------------------------------------------------------------------ drawing (optional extra)


@requires_figures
def test_two_renders_of_one_result_are_byte_identical_in_every_format() -> None:
    """The Done-when: regenerating a figure produces the same file, not merely a similar one."""
    first = pi_window.build({"diversity": _diversity_result()})
    second = pi_window.build({"diversity": _diversity_result()})
    assert first.formats == FIGURE_FORMATS
    for fmt in FIGURE_FORMATS:
        assert first.payload(fmt) == second.payload(fmt), f"{fmt} differs between two renders"
        assert len(first.payload(fmt)) > 1000, f"{fmt} looks empty"


@requires_figures
def test_the_vector_output_carries_no_wall_clock_stamp() -> None:
    """The measured reason the render is deterministic."""
    figure = pi_window.build({"diversity": _diversity_result()})
    assert b"<dc:date>" not in figure.payload("svg")
    assert b"%%CreationDate: Thu Jan  1 00:00:00 1970" in figure.payload("eps")
    assert b"1970" not in figure.payload("svg")


@requires_figures
def test_every_number_on_the_canvas_is_the_one_the_table_renders() -> None:
    """The figure prints the summary cell's own rendering — no second rounding anywhere."""
    figure = pi_window.build({"diversity": _diversity_result()})
    svg = figure.payload("svg").decode("utf-8")
    assert figure.title in svg, "the caption is not on the canvas"
    assert "0.007181" in svg  # Table.format of the π cell
    assert "600 bp window, 200 bp step" in svg


@requires_figures
def test_the_svg_keeps_its_labels_as_text_and_the_pdf_embeds_truetype() -> None:
    """SVG is the editable copy, PDF the submission copy — stated, and checked."""
    figure = pi_window.build({"diversity": _diversity_result()})
    assert b"<text" in figure.payload("svg")
    assert b"/FontFile2" in figure.payload("pdf")  # TrueType (fonttype 42)


@requires_figures
def test_a_changed_cell_changes_the_figure() -> None:
    """A figure that ignores its table would pass every other test in this file."""
    changed = _diversity_result()
    table = changed.tables["diversity_windows"]
    rows = list(table.rows)
    rows[0] = {**rows[0], "pi": 0.5}
    changed = Result(
        module=changed.module,
        module_version=changed.module_version,
        tables={
            **changed.tables,
            "diversity_windows": Table(
                id=table.id,
                title=table.title,
                columns=table.columns,
                rows=tuple(rows),
                scope=table.scope,
            ),
        },
    )
    assert pi_window.build({"diversity": _diversity_result()}).payload("svg") != pi_window.build(
        {"diversity": changed}
    ).payload("svg")


# ------------------------------------------------------------------- the run-level pass


def test_figures_render_from_a_run_directory_alone(tmp_path: Path) -> None:
    """P5-02's Done-when, proved here for the first figure: the run directory is the input."""
    run = tmp_path / "standin" / "20261005-000000-abcdef12"
    module_dir = run / "01_diversity"
    module_dir.mkdir(parents=True)
    (module_dir / "result.json").write_text(export_json(_diversity_result()), encoding="utf-8")

    outcome = render_run_figures(run, figures=["pi_window"], formats=("svg",))
    assert [path.name for path in outcome.files] == ["pi_window.svg"]
    assert (run / "figures" / "pi_window.svg").is_file()
    assert outcome.skipped == ()
    # The run's own layout must agree with the renderer about where a figure lives.
    layout = RunLayout(out_root=tmp_path, dataset="standin", run_id="20261005-000000-abcdef12")
    assert layout.figures_dir == run / "figures"


def test_a_figure_whose_tables_the_run_lacks_is_skipped_by_name(tmp_path: Path) -> None:
    run = tmp_path / "run"
    module_dir = run / "01_structure"
    module_dir.mkdir(parents=True)
    structure = Result(
        module="structure",
        module_version="0.2.0",
        tables={
            "regions": Table(
                id="regions",
                title="Regions",
                columns=(Column("region", "Region"),),
                rows=({"region": "LSC"},),
                scope=Scope(label="normalized", coord_system="normalized"),
                key="region",
            )
        },
    )
    (module_dir / "result.json").write_text(export_json(structure), encoding="utf-8")

    outcome = render_run_figures(run, formats=("svg",))
    assert outcome.files == ()
    skipped = dict(outcome.skipped)
    assert (
        skipped["pi_window"]
        == "needs diversity.diversity_windows, diversity.diversity_summary — not in this run"
    )


def test_a_directory_that_is_not_a_run_is_named(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="run directory"):
        read_run_results(tmp_path / "empty")


def test_the_catalogue_registers_every_spec_by_its_own_id() -> None:
    specs = all_figures()
    assert specs, "the figure catalogue is empty"
    assert len({spec.id for spec in specs}) == len(specs)
    for spec in specs:
        assert get(spec.id) is spec
    with pytest.raises(KeyError, match="unknown figure"):
        get("no_such_figure")
