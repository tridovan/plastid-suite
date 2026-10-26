"""Smoke tests: the package imports and the CLI answers."""

from __future__ import annotations

from pathlib import Path

import pytest

import plastid
from plastid.cli import main


def test_version_exposed() -> None:
    assert plastid.__version__
    assert plastid.__version__ != "0.0.0+unknown"  # the test env installs the package


def test_cli_version(capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["--version"]) == 0
    assert capsys.readouterr().out.startswith("plastid ")


def test_cli_no_args_prints_help(capsys: pytest.CaptureFixture[str]) -> None:
    assert main([]) == 0
    assert "usage: plastid" in capsys.readouterr().out


def test_cli_run_requires_its_arguments() -> None:
    with pytest.raises(SystemExit) as exit_info:
        main(["run"])  # argparse usage error
    assert exit_info.value.code == 2


def test_cli_panel_pin_needs_accessions_or_groups() -> None:
    with pytest.raises(SystemExit) as exit_info:
        main(["panel", "pin"])  # the mutually exclusive pair is required
    assert exit_info.value.code == 2


def test_cli_panel_pin_reports_a_missing_groups_file(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    assert main(["panel", "pin", "--groups", str(tmp_path / "nope.yaml")]) == 1
    assert "nope.yaml" in capsys.readouterr().err


def test_cli_panel_pin_refuses_duplicate_accessions(capsys: pytest.CaptureFixture[str]) -> None:
    code = main(["panel", "pin", "--accessions", "NC_050833.1,NC_050833.1"])
    assert code == 2  # the command line itself is wrong — a usage error
    assert "duplicated" in capsys.readouterr().err


def test_cli_panel_pin_reports_malformed_groups_yaml(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    path = tmp_path / "groups.yaml"
    path.write_text("a: [\n", encoding="utf-8")
    assert main(["panel", "pin", "--groups", str(path)]) == 1
    assert "invalid YAML" in capsys.readouterr().err


def test_cli_panel_pin_fails_when_a_record_cannot_be_read(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """A pinned record our reader cannot parse is a failed dataset, not a quiet success."""

    class FakeClient:
        def resolve(self, accession: str) -> tuple[dict, str]:
            return {"accessionversion": accession, "organism": "X", "taxid": 1}, ""

        def fetch_genbank(self, accession: str) -> str:
            return "this is not a GenBank file"

    monkeypatch.setattr("plastid.cli.NcbiClient", FakeClient)
    code = main(["panel", "pin", "--accessions", "XX_000000.1", "--out", str(tmp_path)])

    assert code == 1
    captured = capsys.readouterr()
    assert "UNVERIFIED" in captured.out
    assert "XX_000000.1" in captured.err


def test_cli_comparison_module_reports_a_missing_aligner(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    data_dir: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """A missing MAFFT is a named tool error and exit 3 — never a silent fallback alignment."""
    monkeypatch.setattr("plastid.external.mafft.which", lambda name: None)
    code = main(
        [
            "run",
            "--module",
            "substitutions",
            "--input",
            str(data_dir / "NC_050833.1.gb"),
            "--input",
            str(data_dir / "NC_024936.1.gb"),
            "--out",
            str(tmp_path),
            "--dataset",
            "pair",
            # NC_024936.1's annotated IR is below the gate; this test is about MAFFT, so use the
            # explicit sequence path (the manifest records that choice).
            "--set",
            "ir_source=sequence",
        ]
    )
    assert code == 3
    assert "mafft not found" in capsys.readouterr().err


def test_cli_comparison_output_carries_the_group_label(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, data_dir: Path
) -> None:
    """P2-11's Done-when: a comparison output carries a shared group label end to end.

    The aligner is stubbed so this stays offline — the labels, not the alignment, are the point.
    """
    from plastid.core.align import Alignment

    def fake_align(rows: list[tuple[str, str]], ctx: object) -> Alignment:
        names = tuple(name for name, _ in rows)
        sequences = ("ACGT" * 10, "ACGT" * 9 + "ACGA")
        return Alignment(names=names, sequences=sequences)

    monkeypatch.setattr("plastid.modules.substitutions.align", fake_align)
    groups = tmp_path / "groups.yaml"
    groups.write_text("wild:\n  - NC_050833.1\ngarden:\n  - NC_024936.1\n", encoding="utf-8")

    code = main(
        [
            "run",
            "--module",
            "substitutions",
            "--input",
            str(data_dir / "NC_050833.1.gb"),
            "--input",
            str(data_dir / "NC_024936.1.gb"),
            "--out",
            str(tmp_path),
            "--dataset",
            "pair",
            "--groups",
            str(groups),
            "--set",
            "ir_source=sequence",
        ]
    )
    assert code == 0
    run_dir = sorted((tmp_path / "pair").iterdir())[0]
    table = (run_dir / "01_substitutions" / "tables" / "substitutions.csv").read_text("utf-8")
    assert "Group" in table and "garden" in table


def test_cli_run_reports_an_input_with_no_group_label(
    tmp_path: Path, data_dir: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """A genome with no label would produce an unlabelled comparison row — refuse instead."""
    groups = tmp_path / "groups.yaml"
    groups.write_text("only:\n  - NC_050833.1\n", encoding="utf-8")
    code = main(
        [
            "run",
            "--module",
            "substitutions",
            "--input",
            str(data_dir / "NC_050833.1.gb"),
            "--input",
            str(data_dir / "NC_024936.1.gb"),
            "--out",
            str(tmp_path),
            "--dataset",
            "pair",
            "--groups",
            str(groups),
            "--set",
            "ir_source=sequence",
        ]
    )
    assert code == 1
    assert "no group for NC_024936.1" in capsys.readouterr().err


def test_cli_comparison_module_needs_two_genomes(tmp_path: Path, data_dir: Path) -> None:
    code = main(
        [
            "run",
            "--module",
            "diversity",
            "--input",
            str(data_dir / "NC_050833.1.gb"),
            "--out",
            str(tmp_path),
            "--dataset",
            "lonely",
        ]
    )
    assert code == 2  # a usage error: the module declares its arity


def test_cli_panel_pin_resolves_an_accession_the_index_cannot_hold(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    data_dir: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """A GenBank-only accession (no RefSeq twin, so no index can list it) must still pin."""
    golden = (data_dir / "NC_050833.1.gb").read_text(encoding="utf-8")

    class FakeClient:
        def resolve(self, accession: str) -> tuple[dict, str]:
            summary = {
                "accessionversion": accession,
                "organism": "Coleus amboinicus",
                "taxid": 204180,
                "slen": 152544,
            }
            return summary, "Eukaryota; Lamiaceae; Coleus"

        def fetch_genbank(self, accession: str) -> str:
            return golden

    monkeypatch.setattr("plastid.cli.NcbiClient", FakeClient)
    code = main(["panel", "pin", "--accessions", "PX923121.1", "--out", str(tmp_path)])

    assert code == 0
    assert "resolved live" in capsys.readouterr().out
    assert (tmp_path / "PX923121.1.gb").is_file()
    assert (tmp_path / "PINNED.json").is_file()


def test_cli_run_reads_a_bare_fasta(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    """A sample off an assembler is a FASTA, and annotating it is what ``plastid annotate`` does —
    so ``plastid run`` has to accept one (P4a-04). The genome a FASTA produces carries no
    features, which is the point: ``orf`` reports that plainly and does not fail on it.

    A module that needs IR positions is *not* usable on a bare sequence yet, and says so through
    the gate rather than through a traceback — which is why this test drives ``orf``.
    """
    fasta = tmp_path / "sample.fasta"
    fasta.write_text(
        ">PX000000.1 Eleutherine bulbosa voucher V1 chloroplast, complete genome\n"
        + "ATGGCG" * 40
        + "\n",
        encoding="utf-8",
    )
    out = tmp_path / "runs"
    assert (
        main(
            [
                "run",
                "--module",
                "orf",
                "--input",
                str(fasta),
                "--out",
                str(out),
                "--dataset",
                "bare",
            ]
        )
        == 0
    )
    run_dir = next((out / "bare").iterdir())
    summary = (run_dir / "01_orf/tables/orf_summary.csv").read_text(encoding="utf-8")
    assert "no CDS feature: sequence-level analyses only" in (
        run_dir / "01_orf/result.json"
    ).read_text(encoding="utf-8")
    assert "PX000000.1" in summary


def test_cli_run_refuses_an_empty_input(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    """The reader is chosen by the first line, so an empty file has no format to choose — and the
    error must say that, not name whichever reader a default would have picked."""
    empty = tmp_path / "empty.fasta"
    empty.write_text("\n\n", encoding="utf-8")
    assert main(["run", "-m", "genome_stats", "-i", str(empty), "--out", str(tmp_path)]) == 1
    assert "empty file" in capsys.readouterr().err
