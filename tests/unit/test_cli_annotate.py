"""``plastid annotate`` — the batch, the organelle, and what a failure does to the run (P3B-09).

Three of the four things this file checks are refusals, because that is what the command mostly
is: it decides whether a set of files may be annotated at all before it annotates anything. The
batch's own behaviour — a failed sample is marked and named, and the run carries on — is checked
against a stub engine, so the test is about the *command* rather than about BLAST; the engine's
own numbers are measured in ``tests/golden/test_plastid_annotate.py``.
"""

from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

import pytest

from plastid.cli import main
from plastid.modules import registry
from plastid.result.result import Column, Result, Scope, Table

MITO_RECORD = """LOCUS       SYN_MITO.1                40 bp    DNA     circular PLN 09-OCT-2026
DEFINITION  synthetic mitogenome, for the organelle check.
ACCESSION   SYN_MITO.1
FEATURES             Location/Qualifiers
     source          1..40
                     /organism="Synthetic plant"
                     /organelle="mitochondrion"
     CDS             1..30
                     /gene="cox1"
                     /translation="MADV"
ORIGIN
        1 acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt
//
"""

PLASTID_RECORD = MITO_RECORD.replace("SYN_MITO.1", "SYN_PLASTID.1").replace(
    '/organelle="mitochondrion"', '/organelle="plastid:chloroplast"'
)


def _write(directory: Path, name: str, text: str) -> Path:
    path = directory / name
    path.write_text(text, encoding="utf-8")
    return path


# Every cell the engine's report template resolves against — the stub run is a real run, so its
# report has to render (tests/unit/test_report.py guards the template's own side of this).
SUMMARY_COLUMNS = (
    "organism",
    "accession",
    "organelle",
    "length_bp",
    "n_calls",
    "n_cds",
    "n_rrna",
    "n_trna",
    "n_not_called",
    "n_cds_failing_orf",
    "n_dropped",
    "n_unnameable",
    "n_overlapping_pairs",
    "min_identity_pct",
    "min_coverage_pct",
    "evalue",
    "transl_table",
    "blast",
    "aragorn",
    "reference",
)


def _artifacts(out: Path) -> list[Path]:
    return sorted(out.glob("**/plastid_annotation.gb"))


def test_a_mito_sample_is_refused_by_name_and_nothing_is_annotated(tmp_path: Path) -> None:
    """H-09's mismatch direction that had no gate at all (measured 05/10/2026).

    ``mito_genes`` ran on a plastome and said nothing; here the plastid path is handed a record
    that declares itself a mitogenome, and the answer has to name both the record and the command
    that does serve it — a refusal a user cannot act on is a dead end.
    """
    sample = _write(tmp_path, "mito.gb", MITO_RECORD)
    out = tmp_path / "runs"
    code = main(
        [
            "annotate",
            "-i",
            str(sample),
            "--reference",
            str(tmp_path),
            "-o",
            str(out),
        ]
    )
    assert code == 2, "the wrong organelle is a usage error: the wrong command was typed"
    assert not _artifacts(out), "a refused run annotated something"


def test_a_mito_fasta_is_refused_by_name_when_the_path_serves_plastids(tmp_path: Path) -> None:
    """A bare sequence declares nothing — D5-06 says the *path* declares for it, so the refusal
    has to come from the one thing in the run that does declare: the reference set."""
    sample = _write(tmp_path, "sample.fa", ">MITO_SAMPLE\nACGTACGTACGTACGTACGT\n")
    reference = tmp_path / "reference"
    reference.mkdir()
    _write(reference, "SYN_MITO.1.gb", MITO_RECORD)
    out = tmp_path / "runs"
    code = main(["annotate", "-i", str(sample), "--reference", str(reference), "-o", str(out)])
    assert code == 2
    assert not _artifacts(out)


def test_a_batch_that_mixes_two_organelles_stops_before_annotating_any_of_them(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """D5-06: a run's input set is single-organelle. A batch is read whole *before* it is
    annotated, so a mixed one cannot leave half its runs written."""
    plastid = _write(tmp_path, "plastid.gb", PLASTID_RECORD)
    mito = _write(tmp_path, "mito.gb", MITO_RECORD)
    out = tmp_path / "runs"
    code = main(
        [
            "annotate",
            "-i",
            str(plastid),
            "-i",
            str(mito),
            "--reference",
            str(tmp_path),
            "-o",
            str(out),
        ]
    )
    assert code == 2
    assert "one organelle" in capsys.readouterr().err
    assert not _artifacts(out)


def test_audit_only_needs_no_sample_and_a_batch_needs_one(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """``--audit-only`` is the flag P3B-03 left for this command: fix and re-check the set
    **without annotating anything** (D5-08) — so it must not demand a sample. Without it, a run
    with no sample is a usage error rather than an empty success.

    The set is named here rather than left to the checkout's default: a tree without the pinned
    panel must still answer the question this test asks (measured on CI, 27/10/2026: the default
    path held no record, so the run exited 1 and the assertion read it as a usage refusal).
    """
    reference = tmp_path / "reference"
    reference.mkdir()
    _write(reference, "SYN_PLASTID.1.gb", PLASTID_RECORD)

    audit = main(
        [
            "annotate",
            "--audit-only",
            "--no-confirm",
            "--reference",
            str(reference),
            "-o",
            str(tmp_path / "runs"),
        ]
    )
    message = capsys.readouterr()
    assert audit == 2  # this synthetic record blocks the audit; the run is not a usage refusal
    assert "input is required" not in message.err
    assert "reference audit" in message.out + message.err

    assert main(["annotate", "-o", str(tmp_path / "runs2")]) == 2
    assert "input is required" in capsys.readouterr().err


def _stub(monkeypatch: pytest.MonkeyPatch) -> None:
    """A registry whose annotate engine is a stub: the batch's behaviour, without BLAST.

    The stub still returns the module's own ``plastid_annotation_summary`` table, because the run's
    report is rendered from it: a stub that skipped it would fail on the template rather than on
    anything this test is about (and it says the summary is part of the module's contract, which
    ``tests/unit/test_report.py`` also holds to).
    """
    spec = registry.get("plastid_annotate")

    def run_fn(genome: object, params: dict, ctx: object) -> Result:
        label = getattr(genome, "accession", None) or getattr(genome, "name", "")
        if label.startswith("BROKEN"):
            raise ValueError("this sample is broken on purpose")
        return Result(
            module="plastid_annotate",
            module_version="0.1.0",
            tables={
                "plastid_annotation_summary": Table(
                    id="plastid_annotation_summary",
                    title="Plastid annotation",
                    columns=tuple(
                        Column(name, name, dtype="int" if name.startswith("n_") else "str")
                        for name in SUMMARY_COLUMNS
                    ),
                    rows=({name: 1 for name in SUMMARY_COLUMNS} | {"accession": label},),
                    scope=Scope(label="stub"),
                ),
                "called_genes": Table(
                    id="called_genes",
                    title="Called genes",
                    columns=(Column("gene", "Gene"),),
                    rows=({"gene": "psbA"},),
                    scope=Scope(label="stub"),
                ),
            },
            artifacts={"plastid_annotation.gb": f"LOCUS       {label}\n//\n"},
        )

    monkeypatch.setattr(
        registry,
        "get",
        lambda module_id: replace(spec, run=run_fn, requires_ir=False),
    )


def test_a_failed_sample_is_marked_and_named_without_killing_the_batch(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    """D5-08: "a failed sample is marked and named **without killing the run**, but the exit code
    reflects that something failed". Both halves are asserted — the good sample's artifact exists,
    and the broken sample's own run directory records why.
    """
    _stub(monkeypatch)
    good = _write(tmp_path, "good.fa", ">GOOD.1\nACGTACGTACGTACGTACGT\n")
    broken = _write(tmp_path, "broken.fa", ">BROKEN.1\nACGTACGTACGTACGTACGT\n")
    out = tmp_path / "runs"

    code = main(["annotate", "-i", str(good), "-i", str(broken), "-o", str(out)])

    assert code != 0, "a batch with a failed sample exited as if it had succeeded"
    assert len(_artifacts(out)) == 1, "the good sample's artifact is missing"
    failures = [
        json.loads(path.read_text(encoding="utf-8"))
        for path in out.glob("**/manifest.json")
        if json.loads(path.read_text(encoding="utf-8"))["status"] == "failed"
    ]
    assert len(failures) == 1
    assert "broken on purpose" in failures[0]["modules"][-1]["error"]
    # Which sample: the manifest says so in its own right, not only through its directory name —
    # `--dataset` puts every sample of a batch under one slug, and then the path cannot say it.
    assert failures[0]["dataset"] == "BROKEN.1"
    assert "1 of 2 sample(s) failed" in capsys.readouterr().err


def test_every_sample_of_a_batch_gets_its_own_run_directory_and_artifact(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """One artifact per sample: the batch is a batch of *runs* (output split by
    run/dataset, never overwritten in place), so two samples never share a directory."""
    _stub(monkeypatch)
    first = _write(tmp_path, "one.fa", ">ONE.1\nACGTACGTACGTACGTACGT\n")
    second = _write(tmp_path, "two.fa", ">TWO.1\nACGTACGTACGTACGTACGT\n")
    out = tmp_path / "runs"

    assert main(["annotate", "-i", str(first), "-i", str(second), "-o", str(out)]) == 0

    artifacts = _artifacts(out)
    assert len(artifacts) == 2
    assert len({path.parent for path in artifacts}) == 2, "two samples shared one run directory"
    # <out>/<sample>/<run-id>/01_plastid_annotate/<artifact> — the sample is the dataset.
    assert sorted(path.parents[2].name for path in artifacts) == ["one.1", "two.1"]


def test_the_run_records_the_artifact_and_the_reference_set_it_searched(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """M11.16: a next step must be able to read *what* the file is and *what it came from*.

    The manifest's ``files`` list says where an output went; this contract says it is a GenBank
    artifact for this sample, circular, in GenBank coordinates, with status ok, transferred from
    a named reference set whose record count and digest are recorded — not inferred from the
    module's parameters or the directory's shape.
    """
    from plastid.run.contract import sha256_file

    _stub(monkeypatch)
    sample = _write(tmp_path, "sample.fa", ">SYN.1\nACGTACGTACGTACGTACGT\n")
    references = tmp_path / "references"
    references.mkdir()
    _write(references, "REF.1.gb", PLASTID_RECORD)
    out = tmp_path / "runs"

    code = main(["annotate", "-i", str(sample), "--reference", str(references), "-o", str(out)])

    assert code == 0
    manifest = json.loads(next(out.glob("**/manifest.json")).read_text(encoding="utf-8"))
    reference_set = manifest["reference_sets"][0]
    assert reference_set["role"] == "annotation_db"
    assert reference_set["records"] == 1
    assert reference_set["leave_out"] == "SYN.1"
    assert reference_set["excluded"] == []
    artifact = manifest["artifacts"][0]
    assert artifact["format"] == "genbank"
    assert artifact["sample_id"] == "SYN.1"
    assert artifact["topology"] == "circular"
    assert artifact["frame"] == "genbank"
    assert artifact["status"] == "ok"
    assert artifact["reference_sha256"] == reference_set["sha256"]
    written = next(out.glob("**/plastid_annotation.gb"))
    assert artifact["sha256"] == sha256_file(written)


def test_a_reference_set_that_excludes_everything_is_refused_not_recorded_empty(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    """An empty reference contract would be a lie a downstream step could trust."""
    _stub(monkeypatch)
    sample = _write(tmp_path, "sample.fa", ">SYN.1\nACGTACGTACGTACGTACGT\n")
    references = tmp_path / "references"
    references.mkdir()
    _write(references, "REF.1.gb", PLASTID_RECORD)
    out = tmp_path / "runs"

    code = main(
        [
            "annotate",
            "-i",
            str(sample),
            "--reference",
            str(references),
            "--exclude",
            "REF.1",
            "-o",
            str(out),
        ]
    )

    assert code == 2
    assert "no GenBank records" in capsys.readouterr().err
    assert not list(out.glob("**/manifest.json"))


def test_an_already_annotated_genbank_is_refused_before_anything_runs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    """M11.16: an external GenBank goes to analysis or verify, not through the annotator.

    The module discards the input's features by design; doing that to an annotated record
    without saying so is the failure this guard exists for. The CLI refuses by name and offers
    the two explicit ways forward.
    """
    _stub(monkeypatch)
    sample = _write(tmp_path, "external.gb", PLASTID_RECORD)
    out = tmp_path / "runs"

    code = main(["annotate", "-i", str(sample), "-o", str(out)])

    assert code == 2
    message = capsys.readouterr().err
    assert "already carries 1 CDS/tRNA/rRNA feature(s)" in message
    assert "--set reannotate=yes" in message
    assert "plastid run -m" in message
    assert not list(out.glob("**/manifest.json")), "a run directory was written before the refusal"


def test_reannotate_yes_makes_replacing_an_external_annotation_explicit(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _stub(monkeypatch)
    sample = _write(tmp_path, "external.gb", PLASTID_RECORD)
    out = tmp_path / "runs"

    code = main(["annotate", "-i", str(sample), "--set", "reannotate=yes", "-o", str(out)])

    assert code == 0
    assert list(out.glob("**/plastid_annotation.gb"))
