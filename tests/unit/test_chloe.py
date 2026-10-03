"""Chloe adapter (P3-03): the wiring around the annotator, and what it refuses to guess.

Chloe needs Julia and its own reference set, so these tests exercise the halves that do not need
either: the named errors when the tool is missing, and the path from "the annotator wrote a GFF" to
"the gene set went through the same checks as a downloaded record". That path runs against a
**stub** in Chloe's place — it tests our code, not Chloe's annotation quality, which is measured on
the annotator's own real output in `test_gff.py` and, on a machine that has the tool, in tier B.
"""

from __future__ import annotations

import stat
import sys
from pathlib import Path

import pytest

from plastid.core.genbank import read_genbank
from plastid.external.chloe import (
    CHLOE_GFF_VERSION,
    CHLOE_REFS_TREE,
    CHLOE_VERSION,
    run_chloe,
)
from plastid.run.context import RunContext
from plastid.run.tools import ToolError, which

_STUB = """#!{python}
import os
import pathlib
import sys

fasta = pathlib.Path(sys.argv[-1])
stem = fasta.with_suffix("")
pathlib.Path(str(stem) + ".chloe.gff").write_text(
    pathlib.Path(os.environ["STUB_GFF"]).read_text(encoding="utf-8"), encoding="utf-8"
)
if os.environ.get("STUB_WRITE_FASTA"):
    # what Chloe writes when (and only when) it transformed the sequence: the same FASTA, with
    # the single-copy region flipped — here, one base is enough to be a different sequence
    sequence = "".join(
        line.strip() for line in fasta.read_text(encoding="utf-8").splitlines()
        if not line.startswith(">")
    )
    pathlib.Path(str(stem) + ".chloe.fa").write_text(
        ">stub\\n" + "A" + sequence[1:] + "\\n", encoding="utf-8"
    )
"""


def _ctx(tmp_path: Path) -> RunContext:
    from plastid.run.layout import RunLayout

    layout = RunLayout.create(tmp_path / "runs", "chloe-test", fingerprint="test")
    layout.make_dirs()
    return RunContext(layout=layout, params={}, module_id="annotation", module_order=1)


def _checkout(
    tmp_path: Path,
    *,
    version: str = CHLOE_VERSION,
    tree: str = CHLOE_REFS_TREE,
) -> Path:
    """A stand-in for a Chloe clone: the three files the adapter reads, and nothing else."""
    checkout = tmp_path / "chloe"
    checkout.mkdir()
    (checkout / "chloe.jl").write_text("import Chloe\n", encoding="utf-8")
    (checkout / "Project.toml").write_text(
        f'name = "Chloe"\nversion = "{version}"\n\n[deps]\nArgParse = "c7e460c6-x"\n',
        encoding="utf-8",
    )
    (checkout / "Artifacts.toml").write_text(
        f'[Chloe_refs]\ngit-tree-sha1 = "{tree}"\n\n'
        f'    [[Chloe_refs.download]]\n    sha256 = "fd20"\n',
        encoding="utf-8",
    )
    return checkout


def _stub_julia(tmp_path: Path) -> Path:
    script = tmp_path / "julia"
    script.write_text(_STUB.format(python=sys.executable), encoding="utf-8")
    script.chmod(script.stat().st_mode | stat.S_IXUSR)
    return script


def _run(
    tmp_path: Path,
    data_dir: Path,
    monkeypatch: pytest.MonkeyPatch,
    *,
    gff: Path | None = None,
    **kwargs,
):
    """The adapter with a stub in Chloe's place, reading the annotator's real golden output —
    or, where a test is about a damaged output, the file that test wrote."""
    monkeypatch.setenv("STUB_GFF", str(gff or data_dir / "chloe_NC_050833.1.gff"))
    return run_chloe(
        read_genbank(data_dir / "NC_050833.1.gb"),
        _ctx(tmp_path),
        seq_id="NC_050833.1",
        # built only when the caller did not bring one: a checkout is read, never written, so a
        # test may reuse it across runs
        script=kwargs.pop("script", None) or _checkout(tmp_path),
        julia=_stub_julia(tmp_path),
        **kwargs,
    )


def test_a_missing_chloe_is_a_named_error_that_says_how_to_get_it(tmp_path: Path) -> None:
    if which("chloe.jl") is not None:
        pytest.skip("chloe.jl is on PATH here — this test is about the machine that has none")
    with pytest.raises(ToolError, match="chloe.jl not found"):
        run_chloe(
            read_genbank(Path(__file__).parents[1] / "data" / "NC_050833.1.gb"),
            _ctx(tmp_path),
            seq_id="NC_050833.1",
        )


def test_a_missing_julia_is_a_named_error(tmp_path: Path) -> None:
    from plastid.external.chloe import julia_binary

    if which("julia") is not None:
        pytest.skip("julia is on PATH here — this test is about the machine that has none")
    with pytest.raises(ToolError, match="julia not found"):
        julia_binary(None)


def test_the_read_back_is_the_golden_gene_set(data_dir: Path, tmp_path: Path, monkeypatch) -> None:
    """The annotator's own file, read by our reader: **132 = 86 + 38 + 8**, every CDS a clean
    ORF. This is P3-03's Done-when, and it holds without Julia running."""
    from collections import Counter

    from plastid.core.model import count_genes
    from plastid.core.orf import check_genome

    produced = _run(tmp_path, data_dir, monkeypatch, transl_table=11)
    counts = count_genes(produced.genome)
    assert (counts.cds, counts.trna, counts.rrna) == (86, 38, 8)
    checks = check_genome(produced.genome, default_table=11)
    statuses = Counter(check.status for check in checks)
    assert statuses["ok"] + statuses["alternative_start"] == 86
    assert not [check for check in checks if check.is_failure]
    assert produced.version == CHLOE_VERSION
    assert produced.reference_tree == CHLOE_REFS_TREE
    assert produced.source_name.endswith(".chloe.gff")


def test_the_module_reports_the_annotator_and_keeps_its_file(
    data_dir: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The module's own path, with the stub in Chloe's place: the same tables a downloaded record
    gets, plus the annotator's raw file as an artifact."""
    from plastid.modules import annotation

    monkeypatch.setenv("STUB_GFF", str(data_dir / "chloe_NC_050833.1.gff"))
    result = annotation.run(
        read_genbank(data_dir / "NC_050833.1.gb"),
        {"chloe": str(_checkout(tmp_path)), "julia": str(_stub_julia(tmp_path))},
        _ctx(tmp_path),
    )

    summary = result.tables["annotation_summary"].rows[0]
    assert (summary["n_genes"], summary["n_cds"], summary["n_trna"], summary["n_rrna"]) == (
        132,
        86,
        38,
        8,
    )
    assert summary["n_orf_failed"] == 0
    assert summary["annotator"] == f"Chloe {CHLOE_VERSION}"
    assert "curated set" in summary["reference"]
    assert result.citations == ("chloe",)
    assert "annotation.gff" in result.artifacts
    assert result.artifacts["annotation.gff"].startswith("##gff-version 3")
    # the run says which frame the coordinates are in, and which translation table was used
    assert any("--no-transform" in note for note in result.notes)
    assert any("translation table(s) 11" in note for note in result.notes)


def test_a_run_that_transformed_the_sequence_is_refused(
    data_dir: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Measured 2026-10-03: Chloe's default flips the single-copy region and writes coordinates in
    that frame — a GFF whose gene *count* is perfect (132) while a third of the genes sit in
    another frame. The stub writes the FASTA Chloe writes in exactly that case, and the run must
    refuse it by name rather than analyse it."""
    monkeypatch.setenv("STUB_WRITE_FASTA", "1")
    with pytest.raises(ToolError, match="transformed the sequence"):
        _run(tmp_path, data_dir, monkeypatch)


def test_a_checkout_of_another_release_is_refused(
    tmp_path: Path, data_dir: Path, monkeypatch
) -> None:
    with pytest.raises(ToolError, match="Chloe 9.9.9"):
        _run(tmp_path, data_dir, monkeypatch, script=_checkout(tmp_path, version="9.9.9"))


def test_a_checkout_with_another_reference_set_is_refused(
    tmp_path: Path, data_dir: Path, monkeypatch
) -> None:
    with pytest.raises(ToolError, match="reference set"):
        _run(tmp_path, data_dir, monkeypatch, script=_checkout(tmp_path, tree="0" * 40))


def test_a_checkout_without_an_artifact_manifest_is_refused(
    tmp_path: Path, data_dir: Path, monkeypatch
) -> None:
    checkout = _checkout(tmp_path)
    (checkout / "Artifacts.toml").unlink()
    with pytest.raises(ToolError, match="Artifacts.toml is missing"):
        _run(tmp_path, data_dir, monkeypatch, script=checkout)


def test_an_output_declaring_another_gff_version_is_refused(
    tmp_path: Path, data_dir: Path, monkeypatch
) -> None:
    text = (data_dir / "chloe_NC_050833.1.gff").read_text(encoding="utf-8")
    edited = tmp_path / "other.gff"
    edited.write_text(text.replace(f"Chloe {CHLOE_GFF_VERSION}", "Chloe 2.0"), encoding="utf-8")
    with pytest.raises(ToolError, match="source-version Chloe 2.0"):
        _run(tmp_path, data_dir, monkeypatch, gff=edited)


def test_an_output_that_declares_no_version_is_refused(
    tmp_path: Path, data_dir: Path, monkeypatch
) -> None:
    text = (data_dir / "chloe_NC_050833.1.gff").read_text(encoding="utf-8")
    edited = tmp_path / "noversion.gff"
    edited.write_text(
        "\n".join(line for line in text.splitlines() if not line.startswith("##source-version"))
        + "\n",
        encoding="utf-8",
    )
    with pytest.raises(ToolError, match="declares no '##source-version"):
        _run(tmp_path, data_dir, monkeypatch, gff=edited)


def test_an_annotation_with_no_genes_is_a_named_error(
    tmp_path: Path, data_dir: Path, monkeypatch
) -> None:
    """A run that reports success and finds nothing is a failed annotation, not a report of
    zeros — the rule PGA's silent-empty mode taught."""
    empty = tmp_path / "empty.gff"
    empty.write_text(
        f"##gff-version 3\n##source-version Chloe {CHLOE_GFF_VERSION}\n"
        "NC_050833.1\t.\tintron\t100\t200\t.\t+\t.\tID=i1;gene=psbA\n",
        encoding="utf-8",
    )
    with pytest.raises(ToolError, match="no CDS, tRNA or rRNA"):
        _run(tmp_path, data_dir, monkeypatch, gff=empty)


def test_the_declared_translation_table_is_stamped_on_cds_that_declare_none(
    data_dir: Path, tmp_path: Path, monkeypatch
) -> None:
    """A GFF says nothing about the genetic code. The caller's declaration is stamped on, and only
    where the file itself is silent — a file that declares one keeps it."""
    checkout = _checkout(tmp_path)
    without = _run(tmp_path, data_dir, monkeypatch, script=checkout)
    assert all(gene.transl_table is None for gene in without.genome.genes if gene.kind == "CDS")

    declared = _run(tmp_path, data_dir, monkeypatch, script=checkout, transl_table=11)
    cds = [gene for gene in declared.genome.genes if gene.kind == "CDS"]
    assert cds and all(gene.transl_table == 11 for gene in cds)
    assert all(gene.transl_table is None for gene in declared.genome.genes if gene.kind == "tRNA")


def test_the_module_refuses_an_annotator_it_does_not_have(data_dir: Path, tmp_path: Path) -> None:
    from plastid.modules import annotation

    with pytest.raises(ValueError, match="unknown annotator"):
        annotation.run(
            read_genbank(data_dir / "NC_050833.1.gb"), {"annotator": "geseq"}, _ctx(tmp_path)
        )


def test_the_module_says_which_annotators_it_has_when_asked_for_one_that_needs_references(
    data_dir: Path, tmp_path: Path
) -> None:
    """PGA without a reference directory names the parameter *and* the way out, because the
    default annotator needs no references at all."""
    from plastid.modules import annotation

    with pytest.raises(ValueError, match="use the default annotator"):
        annotation.run(
            read_genbank(data_dir / "NC_050833.1.gb"), {"annotator": "pga"}, _ctx(tmp_path)
        )
