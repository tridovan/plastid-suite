"""The R figure pass — every path that must not become a crash or a silent blank panel.

The pass never runs R here: the runner is injected, which is the whole reason it is injected. What
is asserted is what the pass itself owns — where files land, that a caption carries the scope
labels of the tables the figure read, that every refusal is a **named skip**, and that a missing R
costs the run nothing but a line of explanation.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from plastid.external.rscript import R_SCRIPTS
from plastid.result.export import export_json
from plastid.result.figures import all_figures, all_r_figures
from plastid.result.result import Column, Result, Scope, Table
from plastid.run.rfigures import render_run_rfigures
from plastid.run.tools import ToolError


def _result() -> Result:
    """A run holding exactly one table set — the RSCU figure's, since it is the smallest."""
    rscu = Table(
        id="rscu",
        title="Relative synonymous codon usage",
        columns=(
            Column("codon", "Codon"),
            Column("amino_acid", "Amino acid"),
            Column("count", "Count", dtype="int"),
            Column("rscu", "RSCU", dtype="float", precision=4),
        ),
        rows=({"codon": "GCA", "amino_acid": "A", "count": 370, "rscu": 1.0655},),
        scope=Scope(label="CDS only, pseudogenes excluded; stops excluded"),
        key="codon",
    )
    summary = Table(
        id="codon_usage",
        title="Codon usage summary",
        columns=(
            Column("accession", "Accession"),
            Column("organism", "Organism"),
            Column("n_cds", "CDS in the set", dtype="int"),
            Column("enc", "ENC", dtype="float", precision=2),
        ),
        rows=(
            {"accession": "NC_050833.1", "organism": "Iris domestica", "n_cds": 86, "enc": 49.69},
        ),
        scope=Scope(label="CDS only, pseudogenes excluded; stops and partial codons not counted"),
    )
    return Result(
        module="codon_usage", module_version="0.1.0", tables={"rscu": rscu, "codon_usage": summary}
    )


def _run_dir(tmp_path: Path, result: Result | None = None) -> Path:
    run_dir = tmp_path / "runs" / "dataset" / "20261017-120000-abcdef12"
    module_dir = run_dir / "04_codon_usage"
    module_dir.mkdir(parents=True)
    (module_dir / "result.json").write_text(export_json(result or _result()), encoding="utf-8")
    return run_dir


def _runner(calls: list, *, name: str = "fig1.png") -> object:
    def runner(script, input_file, input_dir, out_dir):
        calls.append((script, input_file, input_dir, out_dir))
        figures = Path(out_dir) / "Figures"
        figures.mkdir(parents=True, exist_ok=True)
        (figures / name).write_bytes(b"\x89PNG\r\n\x1a\n")
        (figures / "fig1.pdf").write_bytes(b"%PDF-1.4\n")
        return figures

    return runner


def test_the_figures_land_in_the_runs_own_figures_directory(tmp_path: Path) -> None:
    run_dir = _run_dir(tmp_path)
    outcome = render_run_rfigures(run_dir, runner=_runner([]))
    names = sorted(path.name for path in outcome.files)
    assert "r_rscu__fig1.png" in names and "r_rscu__fig1.pdf" in names
    assert (run_dir / "figures" / "r_rscu__fig1.png").is_file()
    # the scripts write into <outdir>/Figures — that scratch never survives into the run
    assert not (run_dir / "figures" / "Figures").exists()
    # and the pivot R actually read is kept beside the image, so the figure can be reproduced
    assert (run_dir / "figures" / "r_rscu__input" / "codon_usage_rscu.csv").is_file()


def test_a_refresh_removes_an_r_figure_its_tables_no_longer_feed(tmp_path: Path) -> None:
    """The same refresh contract as the matplotlib pass (M11.18).

    R names its own output files, so the whole ``<id>__*`` space belongs to the figure. When the
    table it read is gone, the images and the caption must be gone with it.
    """
    run_dir = _run_dir(tmp_path)
    render_run_rfigures(run_dir, runner=_runner([]), figures=["r_rscu"])
    figures = run_dir / "figures"
    assert (figures / "r_rscu__fig1.png").is_file()
    assert (figures / "r_rscu.caption.txt").is_file()

    # Rewrite the stored result without the codon-usage tables the figure reads.
    module_dir = run_dir / "04_codon_usage"
    from plastid.result.export import export_json
    from plastid.result.result import Result

    (module_dir / "result.json").write_text(
        export_json(Result(module="codon_usage", module_version="test", tables={})),
        encoding="utf-8",
    )

    outcome = render_run_rfigures(run_dir, runner=_runner([]), figures=["r_rscu"])
    assert outcome.files == ()
    assert not (figures / "r_rscu__fig1.png").exists()
    assert not (figures / "r_rscu.caption.txt").exists()
    (record,) = outcome.records
    assert record.status == "skipped" and "codon_usage.rscu" in record.reason


def test_the_caption_carries_the_scope_label_of_every_table_the_figure_read(
    tmp_path: Path,
) -> None:
    run_dir = _run_dir(tmp_path)
    render_run_rfigures(run_dir, runner=_runner([]))
    caption = (run_dir / "figures" / "r_rscu.caption.txt").read_text(encoding="utf-8")
    result = _result()
    for table_id in ("rscu", "codon_usage"):
        assert result.tables[table_id].scope.label in caption
        assert f"codon_usage.{table_id} — scope: " in caption


def test_a_preflight_reason_skips_every_r_figure_and_runs_nothing(tmp_path: Path) -> None:
    """No Rscript, or a package R does not have: a named skip, and the runner is never called."""
    run_dir = _run_dir(tmp_path)
    calls: list = []
    outcome = render_run_rfigures(
        run_dir,
        runner=_runner(calls),
        preflight=lambda: (
            "Rscript not found — R is the plotting engine behind the manuscript's figures"
        ),
    )
    assert calls == []
    assert outcome.files == ()
    assert {figure_id for figure_id, _ in outcome.skipped} == {s.id for s in all_r_figures()}
    assert all("Rscript not found" in reason for _id, reason in outcome.skipped)
    # carried out of the pass so the command can exit 3 (no tool) rather than 1 (no data)
    assert outcome.blocked is not None and "Rscript not found" in outcome.blocked


def test_a_run_that_simply_lacks_the_tables_is_not_blocked(tmp_path: Path) -> None:
    """The other half of the exit code: a missing table is a data problem, not a tool problem."""
    run_dir = _run_dir(tmp_path)
    outcome = render_run_rfigures(run_dir, runner=_runner([]), figures=["r_kaks"])
    assert outcome.blocked is None
    assert outcome.files == ()


def test_a_tool_error_is_a_named_skip_and_the_other_figures_still_draw(tmp_path: Path) -> None:
    run_dir = _run_dir(tmp_path)

    def failing(script, input_file, input_dir, out_dir):
        raise ToolError(f"rscript: R is missing {script} — install them (docs/tools/r.md)")

    outcome = render_run_rfigures(run_dir, runner=failing)
    assert outcome.files == ()
    skipped = dict(outcome.skipped)
    assert "r_rscu" in skipped and "missing" in skipped["r_rscu"]


def test_a_table_the_run_does_not_hold_is_named_rather_than_drawn_empty(tmp_path: Path) -> None:
    run_dir = _run_dir(tmp_path)
    calls: list = []
    outcome = render_run_rfigures(run_dir, runner=_runner(calls), figures=["r_kaks"])
    assert calls == []
    assert outcome.skipped[0][0] == "r_kaks"
    assert "kaks.kaks_genes" in outcome.skipped[0][1]


def test_only_the_named_figure_is_rendered(tmp_path: Path) -> None:
    run_dir = _run_dir(tmp_path)
    calls: list = []
    outcome = render_run_rfigures(run_dir, runner=_runner(calls), figures=["r_rscu"])
    assert [call[0] for call in calls] == ["rscu"]
    assert {path.name for path in outcome.files} == {
        "r_rscu__fig1.png",
        "r_rscu__fig1.pdf",
        "r_rscu.caption.txt",
    }


def test_rendering_twice_replaces_the_files_and_does_not_duplicate_them(tmp_path: Path) -> None:
    run_dir = _run_dir(tmp_path)
    first = render_run_rfigures(run_dir, runner=_runner([]), figures=["r_rscu"])
    second = render_run_rfigures(run_dir, runner=_runner([]), figures=["r_rscu"])
    assert sorted(path.name for path in first.files) == sorted(path.name for path in second.files)
    assert len(list((run_dir / "figures").glob("r_rscu__fig1.*"))) == 2  # one png, one pdf


def test_every_r_figure_names_a_packaged_script() -> None:
    """The catalogue-to-adapter binding, which cannot be asserted in ``result/`` (layering)."""
    for spec in all_r_figures():
        assert spec.script in R_SCRIPTS, spec.id


def test_the_two_catalogues_do_not_collide() -> None:
    r_ids = {spec.id for spec in all_r_figures()}
    assert all(figure_id.startswith("r_") for figure_id in r_ids)
    assert r_ids.isdisjoint({spec.id for spec in all_figures()})
    assert len(r_ids) == len(all_r_figures())


@pytest.mark.tier_b
def test_a_real_run_draws_its_r_figures_with_the_real_engine(
    tmp_path: Path, data_dir: Path
) -> None:
    """The whole path end to end, when R and its nine packages are on this machine.

    ``PTS_RSCRIPT`` points at the interpreter, the way ``PTS_CPJSDRAW`` points at the patched
    CPJSdraw: R is never bundled, so a machine that has it keeps it where it likes. The golden
    record is the input because F3 and F4 draw only the species on their own ordering list, and
    *Iris domestica* is one of the 34.

    Skipped with a printed reason on a machine without R — never a failure, and never silent.
    """
    import os
    import shutil

    from plastid.cli import main

    binary = os.environ.get("PTS_RSCRIPT") or shutil.which("Rscript")
    if not binary:
        pytest.skip("no Rscript: set PTS_RSCRIPT=/path/to/Rscript (docs/tools/r.md)")

    out = tmp_path / "runs"
    assert (
        main(["run", "-i", str(data_dir / "NC_050833.1.gb"), "-o", str(out), "-m", "codon_usage"])
        == 0
    )
    run_dir = next(out.glob("*/*/manifest.json")).parent
    assert main(["figures", str(run_dir), "--figure", "r_rscu", "--set", f"rscript={binary}"]) == 0

    figures = run_dir / "figures"
    assert (figures / "r_rscu__Codon_Usage_Comparison.png").is_file()
    assert (figures / "r_rscu__Codon_Usage_Comparison.pdf").is_file()
    caption = (figures / "r_rscu.caption.txt").read_text(encoding="utf-8")
    assert "codon_usage.rscu — scope: " in caption
    assert (
        "Iris_domestica"
        in (figures / "r_rscu__input" / "codon_usage_rscu.csv")
        .read_text(encoding="utf-8")
        .splitlines()[0]
    )
    # the run's provenance carries the interpreter that drew it
    manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    assert any(record["name"] == "Rscript" for record in manifest["figures"]["tools"])


def test_the_tree_figure_is_skipped_by_name_until_a_tree_is_named(tmp_path: Path) -> None:
    """PTS builds no tree, so this one figure has nothing to draw until the user names one."""
    run_dir = _run_dir(tmp_path)
    calls: list = []
    outcome = render_run_rfigures(run_dir, runner=_runner(calls), figures=["r_tree"])
    assert calls == []
    assert outcome.files == ()
    figure_id, reason = outcome.skipped[0]
    assert figure_id == "r_tree" and "--set tree=<file>" in reason


def test_the_tree_figure_draws_from_the_users_own_file_and_records_its_digest(
    tmp_path: Path,
) -> None:
    import hashlib

    run_dir = _run_dir(tmp_path)
    tree = tmp_path / "cay.treefile"
    tree.write_text("((a,b),c);\n", encoding="utf-8")
    calls: list = []
    outcome = render_run_rfigures(
        run_dir, runner=_runner(calls), figures=["r_tree"], params={"tree": str(tree)}
    )
    assert not outcome.skipped and outcome.files
    # the tree is handed to R as its input, and is not copied into the run
    assert calls[0][1] == tree
    assert not (run_dir / "figures" / "r_tree__input").exists()
    caption = (run_dir / "figures" / "r_tree.caption.txt").read_text(encoding="utf-8")
    digest = hashlib.sha256(tree.read_bytes()).hexdigest()
    assert f"sha256 {digest}" in caption and str(tree) in caption
    # no table was read, so there is no scope label to carry — said, not left blank
    assert "no scope label" in caption


def test_a_tree_path_that_is_not_a_file_is_named(tmp_path: Path) -> None:
    run_dir = _run_dir(tmp_path)
    calls: list = []
    outcome = render_run_rfigures(
        run_dir, runner=_runner(calls), figures=["r_tree"], params={"tree": str(tmp_path / "no")}
    )
    assert calls == [] and "is not a file" in outcome.skipped[0][1]


def test_an_r_figure_id_must_be_a_plain_name() -> None:
    from plastid.result.rfigure import RFigureSpec

    with pytest.raises(ValueError, match="without dots"):
        RFigureSpec(
            id="r.bad",
            title="x",
            requires=(),
            script="rscu",
            input_mode="file",
            frame=lambda results: None,  # type: ignore[arg-type,return-value]
        )
    with pytest.raises(ValueError, match="input_mode"):
        RFigureSpec(
            id="r_bad",
            title="x",
            requires=(),
            script="rscu",
            input_mode="stdin",
            frame=lambda results: None,  # type: ignore[arg-type,return-value]
        )


def test_the_cli_forwards_the_species_options_to_f3(tmp_path: Path) -> None:
    """``--set species-order``/``exclude-species`` reach R from the real command (M11.7).

    Measured 27/10/2026: F3/F4 parse both options and the golden passed them straight to
    ``run_figure``, but ``plastid figures`` forwarded only the tree's flags, so a user could not
    order or exclude species through the product.
    """
    from types import SimpleNamespace

    from plastid.cli import _render_r_figures

    run_dir = _run_dir(tmp_path)
    rscript = tmp_path / "Rscript"
    rscript.write_text(
        "#!/bin/sh\n"
        'case "$1" in\n'
        "  --version|-version|-v|-V) echo 'R 4.1.2'; exit 0;;\n"
        "  -e) exit 0;;\n"
        "esac\n"
        'out=""\nprev=""\n'
        'for a in "$@"; do\n'
        '  if [ "$prev" = "--outdir" ]; then out="$a"; fi\n'
        '  prev="$a"\n'
        "done\n"
        'mkdir -p "$out/Figures"\n'
        'echo png > "$out/Figures/rscu.png"\n',
        encoding="utf-8",
    )
    rscript.chmod(0o755)

    args = SimpleNamespace(
        figure=["r_rscu"],
        sets=[
            f"rscript={rscript}",
            "species-order=Beta,Alpha",
            "exclude-species=Gamma",
            "width=10",
            "height=12",
        ],
    )
    outcome, ctx = _render_r_figures(run_dir, sets=args.sets, figures=args.figure)

    assert outcome.files, "the stub Rscript should have drawn a figure"
    argv = next(record.argv for record in ctx.tools if record.name == "rscript:rscu")
    for flag, value in (
        ("--species-order", "Beta,Alpha"),
        ("--exclude-species", "Gamma"),
        ("--width", "10"),
        ("--height", "12"),
    ):
        assert argv[argv.index(flag) + 1] == value
