"""Tier B golden checks — need external tools; a missing tool is named, never silent.

Run on a tool-equipped machine with ``PTS_REQUIRE_TOOLS=1`` to turn skips into failures.
"""

from __future__ import annotations

import csv
import io
import json
from pathlib import Path

import pytest

from conftest import require_tools
from plastid.cli import main
from plastid.run.manifest import validate

pytestmark = pytest.mark.tier_b


def test_tool_probe_reports_what_is_available(tools: dict[str, str | None]) -> None:
    assert set(tools) == {"perl", "misa.pl", "mafft"}


def test_ssr_reproduces_the_golden_numbers(
    tools: dict[str, str | None], data_dir: Path, tmp_path: Path
) -> None:
    require_tools(tools, "perl", "misa.pl")
    assert (
        main(
            [
                "run",
                "--module",
                "ssr",
                "--input",
                str(data_dir / "NC_050833.1.gb"),
                "--out",
                str(tmp_path),
                "--dataset",
                "golden",
            ]
        )
        == 0
    )
    run_dir = sorted((tmp_path / "golden").iterdir())[0]
    manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    validate(manifest)

    tool_records = manifest["modules"][0]["tools"]
    assert [record["name"] for record in tool_records] == ["misa"]
    record = tool_records[0]
    assert record["version"] == "2.1"  # pinned: the vendored MISA
    assert record["interpreter"].startswith("perl ")
    assert record["rc"] == 0
    assert Path(record["stdout_file"]).is_file()

    summary_csv = (run_dir / "01_ssr" / "tables" / "ssr_summary.csv").read_text(encoding="utf-8")
    rows = csv.reader(
        io.StringIO(
            "\n".join(line for line in summary_csv.splitlines() if not line.startswith("#"))
        )
    )
    header = next(rows)
    values = dict(zip(header, next(rows), strict=True))
    assert values["Total SSRs"] == "37"  # PROJECT.md §6 golden number
    assert values["Monomeric"] == "31"
    assert values["Dinucleotide"] == "4"
    assert values["Trinucleotide"] == "2"
    assert values["In compound formation"] == "3"

    report = (run_dir / "01_ssr" / "report.md").read_text(encoding="utf-8")
    assert "**37 SSRs**" in report
    assert "Draft — needs human review." in report
