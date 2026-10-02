"""MISA adapter — the vendored, md5-pinned MISA 2.1 (PROJECT.md §7.8, §10.2).

MISA writes its outputs next to its input, so every invocation runs in an isolated
temporary directory. The bundled ``misa.ini`` is shipped unmodified (``GFF: true``):
the project uses MISA's defaults, not tuned thresholds.

Counts come from ``<input>.statistics`` (the authoritative totals); positions and motifs
come from the GFF the same run produces.
"""

from __future__ import annotations

import shutil
import tempfile
from dataclasses import dataclass
from pathlib import Path

from plastid.core.model import Genome
from plastid.run.context import RunContext
from plastid.run.tools import ToolError, tool_version, which

__all__ = ["MISA_VERSION", "MisaResult", "SsrFeature", "misa_script", "run_misa"]

MISA_VERSION = "2.1"
VENDOR_DIR = Path(__file__).resolve().parents[1] / "vendor" / "misa"

UNIT_LABELS = {1: "mono", 2: "di", 3: "tri", 4: "tetra", 5: "penta", 6: "hexa"}
_TYPE_UNITS = {
    "monomeric_repeat": 1,
    "dinucleotide_repeat": 2,
    "trinucleotide_repeat": 3,
    "tetranucleotide_repeat": 4,
    "pentanucleotide_repeat": 5,
    "hexanucleotide_repeat": 6,
}


@dataclass(frozen=True, slots=True)
class SsrFeature:
    """One SSR from the GFF (1-based inclusive coordinates, as MISA reports them)."""

    unit_size: int
    motif: str
    start: int
    end: int

    @property
    def unit_label(self) -> str:
        return UNIT_LABELS.get(self.unit_size, f"{self.unit_size}-unit")


@dataclass(frozen=True, slots=True)
class MisaResult:
    total: int
    by_unit: dict[int, int]
    compound: int
    features: tuple[SsrFeature, ...]


def misa_script() -> Path:
    script = VENDOR_DIR / "misa.pl"
    if not script.is_file():
        raise ToolError(f"bundled MISA script is missing: {script}")
    return script


def run_misa(genome: Genome, ctx: RunContext) -> MisaResult:
    """Run the bundled MISA on one genome; returns counts and per-SSR positions."""
    perl = which("perl")
    if perl is None:
        raise ToolError(
            "misa: perl not found — MISA 2.1 is a Perl script; install perl to run this module"
        )
    script = misa_script()
    interpreter_version, _ = tool_version("perl", "perl", perl)
    interpreter = f"perl {interpreter_version}" if interpreter_version else "perl"

    with tempfile.TemporaryDirectory(prefix="misa_") as work:
        workdir = Path(work)
        shutil.copy(script, workdir / "misa.pl")
        shutil.copy(VENDOR_DIR / "misa.ini", workdir / "misa.ini")
        fasta = workdir / "input.fasta"
        fasta.write_text(_to_fasta(genome), encoding="utf-8")

        ctx.tool(
            "misa",
            [perl, str(workdir / "misa.pl"), fasta.name],
            kind="perl",
            tool_path=str(script),
            pinned_version=MISA_VERSION,
            interpreter=interpreter,
            cwd=workdir,
        )

        statistics = fasta.with_name(fasta.name + ".statistics")
        if not statistics.is_file():
            raise ToolError("misa: the run produced no .statistics file — see the captured logs")
        total, by_unit, compound = _parse_statistics(statistics.read_text(encoding="utf-8"))
        gff_files = sorted(workdir.glob("*.gff"))
        features = _parse_gff(gff_files[0].read_text(encoding="utf-8")) if gff_files else ()
        return MisaResult(total=total, by_unit=by_unit, compound=compound, features=features)


def _to_fasta(genome: Genome, width: int = 60) -> str:
    header = genome.accession or genome.name
    sequence = genome.sequence
    body = "\n".join(sequence[i : i + width] for i in range(0, len(sequence), width))
    return f">{header}\n{body}\n"


def _parse_statistics(text: str) -> tuple[int, dict[int, int], int]:
    total = _int_after(text, "Total number of identified SSRs:")
    compound = _int_after(text, "Number of SSRs present in compound formation:")
    by_unit: dict[int, int] = {}
    for line in text.splitlines():
        fields = line.split("\t")
        if len(fields) == 2 and fields[0].strip().isdigit() and fields[1].strip().isdigit():
            by_unit[int(fields[0])] = int(fields[1])
    return total, by_unit, compound


def _parse_gff(text: str) -> tuple[SsrFeature, ...]:
    features: list[SsrFeature] = []
    for line in text.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 9 or fields[2] not in _TYPE_UNITS:
            continue
        note = fields[8]
        motif = ""
        if "Note=" in note:
            value = note.split("Note=", 1)[1].split(";", 1)[0]
            motif = value.split(",", 1)[1] if "," in value else value
        features.append(
            SsrFeature(
                unit_size=_TYPE_UNITS[fields[2]],
                motif=motif,
                start=int(fields[3]),
                end=int(fields[4]),
            )
        )
    return tuple(features)


def _int_after(text: str, marker: str) -> int:
    for line in text.splitlines():
        if marker in line:
            return int(line.split(marker, 1)[1].strip())
    raise ToolError(f"misa: expected {marker!r} in the statistics output, found none")
