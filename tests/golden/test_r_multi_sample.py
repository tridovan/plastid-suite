"""M11.7 — Figures 3 and 4 draw the species their input carries, at panel scale.

The manuscript's scripts used to intersect a hard-coded 34-name Iridaceae list: a species outside
it was dropped without a word, and a panel outside Iridaceae came out empty. The scripts now read
the species columns from the input, in the input's order, with ``--species-order`` and
``--exclude-species`` as explicit parameters. This check builds the wide tables from the 80-record
validation panel (10 families, none of them the manuscript's set) using the toolkit's own
``codon_usage`` and ``amino_acid`` module outputs, then runs the packaged R scripts.

The panel and R are both optional: a missing one skips with a named reason.
"""

from __future__ import annotations

import csv
import io
import shutil
from pathlib import Path

import pytest

from plastid.core.genbank import read_genbank
from plastid.external.rscript import missing_packages, run_figure
from plastid.modules import amino_acid, codon_usage
from plastid.result.figures.r_common import AMINO_ACID_NAMES, species_column
from plastid.run.context import RunContext
from plastid.run.layout import RunLayout

ROOT = Path(__file__).resolve().parents[2]
PANEL = ROOT / ".panel" / "ir-validation" / "source"


def _rscript() -> str:
    found = shutil.which("Rscript")
    if found is None:
        pytest.skip("no Rscript on PATH (docs/tools/r.md)")
    return found


def _context(tmp_path: Path, module: str) -> RunContext:
    layout = RunLayout.create(tmp_path / "runs", module, fingerprint="test")
    layout.make_dirs()
    return RunContext(layout=layout, params={}, module_id=module, module_order=1, ir=None)


def _csv(header: list[str], rows: list[list[str]]) -> str:
    buffer = io.StringIO()
    writer = csv.writer(buffer, lineterminator="\n")
    writer.writerow(header)
    writer.writerows(rows)
    return buffer.getvalue()


def _require_r(ctx: RunContext) -> str:
    binary = _rscript()
    missing = missing_packages(ctx, binary)
    if missing:
        pytest.skip(f"R is missing {', '.join(missing)}")
    return binary


@pytest.fixture(scope="module")
def wide_tables(tmp_path_factory: pytest.TempPathFactory) -> dict[str, object]:
    if not PANEL.is_dir():
        pytest.skip("the IR validation panel is not in this tree (.panel/ is gitignored)")
    work = tmp_path_factory.mktemp("multi_sample")
    ctx = _context(work, "panel")
    genomes = [read_genbank(path) for path in sorted(PANEL.glob("*.gb"))]
    species = [species_column(genome.name) for genome in genomes]

    codons: dict[str, dict[str, tuple[str, str]]] = {}
    residues: dict[str, dict[str, tuple[str, str]]] = {}
    for column, genome in zip(species, genomes, strict=True):
        rscu = codon_usage.run(genome, {}, ctx).tables["rscu"]
        amino = amino_acid.run(genome, {}, ctx).tables["amino_acid_composition"]
        codons[column] = {
            row["codon"]: (
                AMINO_ACID_NAMES.get(row["amino_acid"], row["amino_acid"]),
                rscu.format(f"{row['codon']}.rscu"),
            )
            for row in rscu.rows
        }
        residues[column] = {
            row["residue"]: (
                AMINO_ACID_NAMES.get(row["residue"], row["residue"]),
                amino.format(f"{row['residue']}.percent"),
            )
            for row in amino.rows
        }

    first = species[0]
    f3 = work / "codon_usage.csv"
    f3.write_text(
        _csv(
            ["Amino_Acid", "Codon", *species],
            [
                [codons[first][codon][0], codon]
                + [codons[column].get(codon, ("", ""))[1] for column in species]
                for codon in codons[first]
            ],
        ),
        encoding="utf-8",
    )
    f4 = work / "amino_acid.csv"
    f4.write_text(
        _csv(
            ["AA_Code", "Amino_Acid", *species],
            [
                [residue, residues[first][residue][0]]
                + [residues[column].get(residue, ("", ""))[1] for column in species]
                for residue in residues[first]
            ],
        ),
        encoding="utf-8",
    )
    return {"f3": f3, "f4": f4, "species": species, "codons": codons, "residues": residues}


def test_the_panel_has_species_outside_the_manuscripts_list(wide_tables) -> None:
    species = wide_tables["species"]
    assert len(species) == 80
    assert "Iris_domestica" not in species


def test_figure3_draws_every_species_column(wide_tables, tmp_path: Path) -> None:
    ctx = _context(tmp_path, "r_rscu")
    binary = _require_r(ctx)

    figures = run_figure(
        ctx,
        "rscu",
        input_file=wide_tables["f3"],
        out_dir=tmp_path / "f3",
        binary=binary,
        check_packages=False,
    )

    assert (figures / "Codon_Usage_Comparison.png").is_file()
    assert (figures / "Codon_Usage_Comparison.pdf").is_file()
    log = Path(ctx.tools[-1].stdout_file).read_text(encoding="utf-8")
    assert "Species: 80" in log


def test_figure4_draws_every_species_column(wide_tables, tmp_path: Path) -> None:
    ctx = _context(tmp_path, "r_amino_acid")
    binary = _require_r(ctx)

    figures = run_figure(
        ctx,
        "amino_acid",
        input_file=wide_tables["f4"],
        out_dir=tmp_path / "f4",
        binary=binary,
        check_packages=False,
    )

    assert (figures / "Amino_Acid_Comparison.png").is_file()
    assert (figures / "Amino_Acid_Comparison.pdf").is_file()
    log = Path(ctx.tools[-1].stdout_file).read_text(encoding="utf-8")
    assert "Species: 80" in log


def test_an_explicit_order_and_exclusion_replace_the_input_order(
    wide_tables, tmp_path: Path
) -> None:
    """The two parameters M11.7 adds: a fixed draw order, and columns a run leaves out."""
    species = list(wide_tables["species"])[:4]
    codons = wide_tables["codons"]
    reference = codons[species[0]]
    path = tmp_path / "small.csv"
    path.write_text(
        _csv(
            ["Amino_Acid", "Codon", *species],
            [
                [reference[codon][0], codon]
                + [codons[column].get(codon, ("", ""))[1] for column in species]
                for codon in reference
            ],
        ),
        encoding="utf-8",
    )

    ctx = _context(tmp_path, "r_rscu")
    binary = _require_r(ctx)
    figures = run_figure(
        ctx,
        "rscu",
        input_file=path,
        out_dir=tmp_path / "small_out",
        extra=[
            "--species-order",
            f"{species[2]},{species[0]}",
            "--exclude-species",
            species[1],
        ],
        binary=binary,
        check_packages=False,
    )

    assert (figures / "Codon_Usage_Comparison.png").is_file()
    log = Path(ctx.tools[-1].stdout_file).read_text(encoding="utf-8")
    assert "Species: 2" in log
