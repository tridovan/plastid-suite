"""K-24 with the real engines: `rna_editing` in a run does not cost it its R figures.

The expected figures are the ones the same run draws without `rna_editing`, measured 01/11/2026 on
the golden record: five R images from `codon_usage`, `amino_acid` and `ssr`, none once `rna_editing`
joined them.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from conftest import require_r, require_tools
from plastid.cli import main

pytestmark = pytest.mark.tier_b


def test_a_run_with_rna_editing_keeps_its_r_figures(tmp_path: Path, data_dir: Path) -> None:
    pytest.importorskip("matplotlib", reason="the 'figures' extra is not installed")
    if shutil.which("Rscript") is None:
        require_r("no Rscript on PATH (docs/tools/r.md)")
    require_tools({"mafft": shutil.which("mafft")}, "mafft")

    out = tmp_path / "runs"
    code = main(
        [
            "run",
            "-m",
            "codon_usage",
            "-m",
            "amino_acid",
            "-m",
            "rna_editing",
            "-i",
            str(data_dir / "NC_050833.1.gb"),
            "-o",
            str(out),
            "--figures",
        ]
    )
    assert code == 0
    figures = next(out.glob("*/*/manifest.json")).parent / "figures"
    assert (figures / "r_rscu__Codon_Usage_Comparison.png").is_file()
    assert (figures / "r_amino_acid__Amino_Acid_Comparison.png").is_file()
    header = (figures / "r_rscu__input" / "codon_usage_rscu.csv").read_text(encoding="utf-8")
    assert header.splitlines()[0] == "Amino_Acid,Codon,Iris_domestica"
