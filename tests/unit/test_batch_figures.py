"""The batch's long table, rebuilt for the multi-sample R figures (M11.6/M11.7)."""

from __future__ import annotations

from pathlib import Path

from plastid.result.result import Column, Result, Scope, Table
from plastid.run import batch_frames
from plastid.run.batch import SampleRecord, read_batch_results, write_batch


def _table(
    table_id: str, columns: tuple[str, ...], rows: list[dict], key: str | None = None
) -> Table:
    return Table(
        id=table_id,
        title=table_id,
        columns=tuple(Column(name, name) for name in columns),
        rows=tuple(rows),
        scope=Scope(label="synthetic"),
        key=key,
    )


def _codon_usage(organism: str, first: str, second: str) -> Result:
    return Result(
        module="codon_usage",
        module_version="test",
        tables={
            "rscu": _table(
                "rscu",
                ("codon", "amino_acid", "count", "rscu"),
                [
                    {"codon": "GCA", "amino_acid": "A", "count": "10", "rscu": first},
                    {"codon": "GCC", "amino_acid": "A", "count": "20", "rscu": second},
                ],
                key="codon",
            ),
            "codon_usage": _table(
                "codon_usage",
                ("accession", "organism", "n_cds", "enc"),
                [{"accession": "X", "organism": organism, "n_cds": "2", "enc": "50.0"}],
            ),
        },
    )


def _amino_acid(organism: str, alanine: str, glycine: str) -> Result:
    return Result(
        module="amino_acid",
        module_version="test",
        tables={
            "amino_acid_composition": _table(
                "amino_acid_composition",
                ("residue", "count", "percent", "hydropathy"),
                [
                    {"residue": "A", "count": "10", "percent": alanine, "hydropathy": "1.8"},
                    {"residue": "G", "count": "5", "percent": glycine, "hydropathy": "-0.4"},
                ],
                key="residue",
            ),
            "amino_acid": _table(
                "amino_acid",
                ("accession", "organism", "n_cds", "n_residues", "gravy"),
                [
                    {
                        "accession": "X",
                        "organism": organism,
                        "n_cds": "2",
                        "n_residues": "15",
                        "gravy": "0.2",
                    }
                ],
            ),
        },
    )


def _ssr(accession: str, total: int) -> Result:
    return Result(
        module="ssr",
        module_version="test",
        tables={
            "ssr_summary": _table(
                "ssr_summary",
                ("accession", "total"),
                [{"accession": accession, "total": str(total)}],
            ),
            "ssr_by_region": _table(
                "ssr_by_region",
                ("region", "total"),
                [
                    {"region": "LSC", "total": "2"},
                    {"region": "SSC", "total": "1"},
                    {"region": "IRa", "total": "1"},
                    {"region": "IRb", "total": "1"},
                ],
                key="region",
            ),
            "ssr_by_unit": _table(
                "ssr_by_unit",
                ("unit", "n"),
                [
                    {"unit": "Mono", "n": "2"},
                    {"unit": "Di", "n": "1"},
                    {"unit": "Tri", "n": "1"},
                ],
                key="unit",
            ),
            "ssrs": _table("ssrs", ("context",), [{"context": "Gene"}, {"context": "IGS"}]),
        },
    )


def _records() -> list[SampleRecord]:
    return [
        SampleRecord(
            sample_id=sample,
            species=species,
            accession=sample,
            input=f"{sample}.gb",
            sha256="0" * 64,
            status="ok",
            exit_code=0,
        )
        for sample, species in (("S1", "Alpha one"), ("S2", "Beta two"))
    ]


def _write_batch(tmp_path: Path):
    results = {
        "S1": {
            "codon_usage": _codon_usage("Alpha one", "1.0", "0.5"),
            "amino_acid": _amino_acid("Alpha one", "6.7", "3.3"),
            "ssr": _ssr("S1", 4),
        },
        "S2": {
            "codon_usage": _codon_usage("Beta two", "0.8", "0.7"),
            "amino_acid": _amino_acid("Beta two", "5.0", "5.0"),
            "ssr": _ssr("S2", 6),
        },
    }
    run_dir = tmp_path / "batch"
    write_batch(
        run_dir=run_dir,
        dataset="test",
        records=_records(),
        results=results,
        params={},
        command=[],
        inputs=[],
    )
    return read_batch_results(run_dir), run_dir


def _batch(tmp_path: Path):
    return _write_batch(tmp_path)[0]


def test_the_long_table_rebuilds_with_a_sample_column(tmp_path: Path) -> None:
    results = _batch(tmp_path)
    table = results["codon_usage"].tables["rscu"]
    assert {"sample", "species"} <= {column.name for column in table.columns}
    assert {row["sample"] for row in table.rows} == {"S1", "S2"}


def test_f3_gets_one_column_per_sample(tmp_path: Path) -> None:
    frame = batch_frames.frame_for("r_rscu", _batch(tmp_path))
    assert frame is not None
    ((name, text),) = frame.inputs()
    lines = text.splitlines()
    assert lines[0] == "Amino_Acid,Codon,Alpha_one,Beta_two"
    assert lines[1].startswith("Alanine,GCA,1.0,0.8")


def test_f4_gets_one_column_per_sample(tmp_path: Path) -> None:
    frame = batch_frames.frame_for("r_amino_acid", _batch(tmp_path))
    assert frame is not None
    ((name, text),) = frame.inputs()
    lines = text.splitlines()
    assert lines[0] == "AA_Code,Amino_Acid,Alpha_one,Beta_two"
    assert lines[1] == "A,Alanine,6.7,5.0"


def test_f5_gets_one_row_per_sample(tmp_path: Path) -> None:
    frame = batch_frames.frame_for("r_ssr", _batch(tmp_path))
    assert frame is not None
    ((name, text),) = frame.inputs()
    lines = text.splitlines()
    assert len(lines) == 3  # header + two samples
    assert lines[1].startswith("S1,4,")
    assert lines[2].startswith("S2,6,")


def test_plastid_figures_draws_the_batch_sample_figure(tmp_path: Path, monkeypatch) -> None:
    """The CLI path: a batch directory feeds F3 with one column per sample, no analysis re-run."""
    import json

    from plastid.cli import main

    _, run_dir = _write_batch(tmp_path)
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

    code = main(
        [
            "figures",
            str(run_dir),
            "--figure",
            "r_rscu",
            "--set",
            f"rscript={rscript}",
        ]
    )

    assert code == 0
    drawn = run_dir / "figures" / "r_rscu__rscu.png"
    assert drawn.is_file()
    payload = json.loads((run_dir / "batch.json").read_text(encoding="utf-8"))
    assert any(entry["path"] == "figures/r_rscu__rscu.png" for entry in payload["figures"]["files"])
    # the pivot R read has two species columns
    pivot = run_dir / "figures" / "r_rscu__input" / "Complete_Codon_Usage_Analysis.csv"
    assert pivot.read_text(encoding="utf-8").splitlines()[0] == (
        "Amino_Acid,Codon,Alpha_one,Beta_two"
    )


def test_batch_figures_draws_the_sample_figures_in_the_same_command(tmp_path: Path) -> None:
    """``plastid batch --figures`` finishes the batch: analysis, aggregate, multi-sample figures."""
    import json

    from plastid.cli import main

    data = Path(__file__).parents[1] / "data"
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
        'echo png > "$out/Figures/fig.png"\n',
        encoding="utf-8",
    )
    rscript.chmod(0o755)
    out = tmp_path / "runs"

    code = main(
        [
            "batch",
            "-m",
            "codon_usage",
            "-i",
            str(data / "NC_050833.1.gb"),
            "-i",
            str(data / "NC_024936.1.gb"),
            "--figures",
            "--label",
            "panel-test",
            "--set",
            f"rscript={rscript}",
            "-o",
            str(out),
        ]
    )

    assert code == 0
    batch_dir = next((out / "batch").glob("*/batch.json")).parent
    assert (batch_dir / "figures" / "r_rscu__fig.png").is_file()
    caption = (batch_dir / "figures" / "r_rscu.caption.txt").read_text(encoding="utf-8")
    # M11.18: a multi-sample figure names the data label, the accession set and both IR statuses;
    # the per-sample mapping stays in samples.csv, which the caption points at.
    assert "data label: panel-test" in caption
    assert "accessions: NC_050833.1, NC_024936.1" in caption
    assert "IR source:" in caption and "IR gate:" in caption
    assert "per-sample accession/IR detail: samples.csv" in caption
    payload = json.loads((batch_dir / "batch.json").read_text(encoding="utf-8"))
    assert payload["figures"]["files"]
    assert payload["figures"]["tools"]  # the R call left its provenance in batch.json
    assert all(entry["renderer"] for entry in payload["figures"]["skipped"])
