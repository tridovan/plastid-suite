"""A publication figure names species. The accession stays in its caption, where it can be checked.

The team's figures 5, 7, 8 and 9 carry species names. A run knows them: the records it read say so.
A figure that printed `NC_050833.1` on its axis was right and could not go into a paper.
"""

from __future__ import annotations

import json
from pathlib import Path

from plastid.result.export import export_json
from plastid.result.figures import r_kaks, r_ssr, r_substitutions
from plastid.result.figures.r_common import genome_labels
from plastid.result.figures.r_tree import relabelled_tree
from plastid.result.result import Column, Result, Scope, Table
from plastid.run.rfigures import render_run_rfigures

NAMES = {"NC_050833.1": "Iris domestica", "NC_024936.1": "Iris gatesii"}


def test_a_genome_is_labelled_by_its_species() -> None:
    assert genome_labels(["NC_050833.1"], NAMES) == {"NC_050833.1": "Iris domestica"}


def test_a_genome_the_run_cannot_name_keeps_its_accession() -> None:
    assert genome_labels(["KX000001.1"], NAMES) == {"KX000001.1": "KX000001.1"}


def test_two_records_of_one_species_stay_two() -> None:
    """One label for both would merge them on the axis, and the figure would look fine."""
    names = {"A.1": "Iris domestica", "B.1": "Iris domestica"}
    assert genome_labels(["A.1", "B.1"], names) == {
        "A.1": "Iris domestica A.1",
        "B.1": "Iris domestica B.1",
    }


def _ssr_frame() -> r_ssr.SsrFrame:
    values = (
        "NC_050833.1",
        "37",
        "27",
        "10",
        "0",
        "31",
        "4",
        "2",
        "0",
        "0",
        "0",
        "11",
        "",
        "",
        "8",
    )
    return r_ssr.SsrFrame(values=(*values, "18"), unsourced=("tRNA", "rRNA"))


def test_the_ssr_row_is_named_for_its_species() -> None:
    labelled = _ssr_frame().labelled(NAMES)
    assert labelled.values[0] == "Iris domestica"
    assert labelled.values[1:] == _ssr_frame().values[1:]


def test_the_substitution_columns_read_as_the_papers_do() -> None:
    """The paper heads a column `reference × sample`, so the comparison names both sides."""
    frame = r_substitutions.SubstitutionsFrame(
        genomes=("NC_024936.1",),
        counts={("NC_024936.1", name): "1" for name in r_substitutions.SUBSTITUTION_TYPES},
        ratios={"NC_024936.1": "1.8416"},
        ratio_table="substitutions",
        reference="NC_050833.1",
    )
    ((_, text),) = frame.labelled(NAMES).inputs()
    assert text.splitlines()[0] == "Substitution_Type,Category,Iris domestica × Iris gatesii"


def test_the_kaks_columns_are_species() -> None:
    frame = r_kaks.KaksFrame(
        genes=("rbcL",),
        genomes=("NC_024936.1",),
        cells={("rbcL", "NC_024936.1"): "0.1"},
        reference="NC_050833.1",
        n_scored="1",
    )
    ((_, text),) = frame.labelled(NAMES).inputs()
    assert text.splitlines()[0] == "Gene\tIris_gatesii"


def test_a_tree_keeps_its_shape_and_changes_its_tips() -> None:
    tree = "(NC_050833.1:0.0025,(NC_024936.1:0.0001,X.1:0.0001)100/100:0.004);\n"
    assert relabelled_tree(tree, NAMES) == (
        "(Iris_domestica:0.0025,(Iris_gatesii:0.0001,X.1:0.0001)100/100:0.004);\n"
    )


def _run_with_tree(tmp_path: Path) -> tuple[Path, Path]:
    run_dir = tmp_path / "runs" / "dataset" / "20261101-120000-abcdef12"
    step = run_dir / "01_phylogeny"
    step.mkdir(parents=True)
    tree = step / "tree.treefile"
    tree.write_text("(NC_050833.1:1,NC_024936.1:1);\n", encoding="utf-8")
    table = Table(
        id="phylogeny_summary",
        title="Summary",
        columns=(Column("loci", "Loci"),),
        rows=({"loci": "1"},),
        scope=Scope(label="test"),
    )
    (step / "result.json").write_text(
        export_json(Result(module="phylogeny", module_version="t", tables={table.id: table})),
        encoding="utf-8",
    )
    inputs = [
        {"path": "a.gb", "accession": "NC_050833.1", "organism": "Iris domestica"},
        {"path": "b.gb", "accession": "NC_024936.1", "organism": "Iris gatesii"},
    ]
    (run_dir / "manifest.json").write_text(json.dumps({"inputs": inputs}), encoding="utf-8")
    return run_dir, tree


def _runner(seen: dict):
    def runner(script, input_file, input_dir, out_dir):
        seen["text"] = Path(input_file).read_text(encoding="utf-8")
        figures = Path(out_dir) / "Figures"
        figures.mkdir(parents=True, exist_ok=True)
        (figures / "tree.png").write_bytes(b"\x89PNG\r\n\x1a\n")
        return figures

    return runner


def test_the_runs_own_tree_is_drawn_with_species_and_left_as_it_was(tmp_path: Path) -> None:
    run_dir, tree = _run_with_tree(tmp_path)
    before = tree.read_bytes()
    seen: dict = {}
    render_run_rfigures(run_dir, runner=_runner(seen), figures=["r_tree"])
    assert seen["text"] == "(Iris_domestica:1,Iris_gatesii:1);\n"
    assert tree.read_bytes() == before
    caption = (run_dir / "figures" / "r_tree.caption.txt").read_text(encoding="utf-8")
    assert "NC_050833.1 = Iris domestica" in caption


def test_a_tree_the_user_named_is_drawn_as_it_is(tmp_path: Path) -> None:
    run_dir, _ = _run_with_tree(tmp_path)
    mine = tmp_path / "mine.treefile"
    mine.write_text("(NC_050833.1:1,NC_024936.1:2);\n", encoding="utf-8")
    seen: dict = {}
    render_run_rfigures(
        run_dir, runner=_runner(seen), figures=["r_tree"], params={"tree": str(mine)}
    )
    assert seen["text"] == "(NC_050833.1:1,NC_024936.1:2);\n"


def test_a_name_newick_cannot_hold_stays_out_of_the_tree_syntax() -> None:
    """A bracket or a colon in an organism name would read as the tree's own structure."""
    names = {"A.1": "Cyperus sp. (Thailand: 3)"}
    assert relabelled_tree("(A.1:0.1,B.1:0.2);", names) == "(Cyperus_sp._Thailand_3:0.1,B.1:0.2);"
