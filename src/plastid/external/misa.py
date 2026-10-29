"""MISA adapter — the external MISA 2.1 (docs/tools/misa.md).

**MISA is not bundled** (D-q, 11/10/2026): its source zip states no licence, and shipping the
script in a wheel is the one act a missing licence forbids. It is obtained instead — download
line and both hashes in ``docs/tools/misa.md`` — and found on ``PATH`` or passed as
``--set misa=/path/to/misa.pl``. The release is **checked, not assumed**: MISA has no version
flag, so it is read out of the script's own header and anything but the pinned 2.1 is refused.

The ``misa.ini`` it runs with **is** the toolkit's own file (``data/misa.ini``): MISA's default
unit/repeat definitions written down rather than copied, because the project uses MISA's defaults
and nothing else — re-measured 10/10/2026, the defaults are what reproduces the published table.

MISA writes its outputs next to its input, so every invocation runs in an isolated
temporary directory.

Counts come from ``<input>.statistics`` (the authoritative totals); positions and motifs
come from the GFF the same run produces.
"""

from __future__ import annotations

import re
import shutil
import tempfile
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

from plastid.core.model import Genome
from plastid.core.resources import data_file
from plastid.run.context import RunContext
from plastid.run.params import ParamValueError
from plastid.run.tools import ToolError, tool_version, which

__all__ = [
    "MISA_VERSION",
    "MisaResult",
    "SsrFeature",
    "format_thresholds",
    "misa_ini_text",
    "misa_script",
    "misa_version",
    "parse_unit_repeats",
    "run_misa",
    "threshold_text",
    "unit_repeats_from_ini",
]

MISA_VERSION = "2.1"
MISA_DOC = "docs/tools/misa.md"
_SCRIPT_NAME = "misa.pl"
# ``###Release date: 25/08/20 (version 2.1)`` — the only place the script names its release.
_VERSION_RE = re.compile(r"\(version ([0-9][0-9.]*)\)")

UNIT_LABELS = {1: "mono", 2: "di", 3: "tri", 4: "tetra", 5: "penta", 6: "hexa"}
_CLASS_UNITS = {
    "monomeric_repeat": 1,
    "dinucleotide_repeat": 2,
    "trinucleotide_repeat": 3,
    "tetranucleotide_repeat": 4,
    "pentanucleotide_repeat": 5,
    "hexanucleotide_repeat": 6,
}


@dataclass(frozen=True, slots=True)
class SsrFeature:
    """One SSR from the GFF (1-based inclusive coordinates, as MISA reports them).

    ``compound`` marks a component of a compound cluster. MISA writes such a cluster as a
    single ``repeat_region`` line that does not separate its components, so every component
    carries the cluster's span — the coordinates are the cluster's, not the component's.
    """

    unit_size: int
    motif: str
    start: int
    end: int
    compound: bool = False

    @property
    def unit_label(self) -> str:
        return UNIT_LABELS.get(self.unit_size, f"{self.unit_size}-unit")


@dataclass(frozen=True, slots=True)
class MisaResult:
    total: int
    by_unit: dict[int, int]
    compound: int
    features: tuple[SsrFeature, ...]
    # The thresholds this run actually used, read back from the ini MISA was handed. Reported
    # rather than assumed, so a user-supplied value and the summary cannot disagree (M11.20).
    unit_repeats: dict[int, int]


_DEFINITION_PREFIX = "definition(unit_size,min_repeats)"


def parse_unit_repeats(text: str) -> dict[int, int]:
    """``"1-12 2-8 3-7"`` → ``{1: 12, 2: 8, 3: 7}``, refused by name when it is not that shape.

    The value is MISA's own ``definition(unit_size,min_repeats)`` syntax without the key, which
    is why it is a string rather than six separate parameters: it is one line of MISA's ini and
    a user copying the published thresholds can paste them as they appear. Units may be a
    subset — a run that wants only mono/di/tri says so — but a repeated unit, a unit outside
    1–6, or a non-positive repeat count is a usage error, never a silently dropped token.
    """
    values: dict[int, int] = {}
    tokens = text.split()
    if not tokens:
        raise ParamValueError(
            "unit_repeats: expected at least one '<unit>-<min_repeats>' pair, e.g. "
            '"1-12 2-8 3-7 4-6 5-5 6-5"'
        )
    for token in tokens:
        unit_text, separator, repeats_text = token.partition("-")
        if not separator or not unit_text.isdigit() or not repeats_text.isdigit():
            raise ParamValueError(
                f"unit_repeats: {token!r} is not '<unit>-<min_repeats>' — expected e.g. "
                '"1-12 2-8 3-7 4-6 5-5 6-5"'
            )
        unit, min_repeats = int(unit_text), int(repeats_text)
        if not 1 <= unit <= 6:
            raise ParamValueError(
                f"unit_repeats: unit {unit} is outside 1–6 (MISA's repeat-unit sizes)"
            )
        if min_repeats < 1:
            raise ParamValueError(
                f"unit_repeats: unit {unit} has min_repeats={min_repeats}; it must be at least 1"
            )
        if unit in values:
            raise ParamValueError(
                f"unit_repeats: unit {unit} is named twice — one threshold per unit size"
            )
        values[unit] = min_repeats
    return dict(sorted(values.items()))


def unit_repeats_from_ini(path: str | Path | None = None) -> dict[int, int]:
    """The thresholds in an ini file's ``definition(unit_size,min_repeats)`` line."""
    ini = Path(path) if path is not None else data_file("misa.ini")
    for line in ini.read_text(encoding="utf-8").splitlines():
        if line.startswith(_DEFINITION_PREFIX):
            return parse_unit_repeats(line.split(":", 1)[1])
    raise ValueError(f"{ini}: no {_DEFINITION_PREFIX} line")


def format_thresholds(values: Mapping[int, int]) -> str:
    """The summary's threshold string: ``(1/10) (2/6) …`` in unit order."""
    return " ".join(f"({unit}/{repeats})" for unit, repeats in sorted(values.items()))


def misa_ini_text(unit_repeats: Mapping[int, int] | None = None) -> str:
    """The ini MISA reads: the shipped file, with the definition line replaced when asked.

    Only the definition line changes. ``interruptions``, ``GFF`` and anything else the shipped
    file declares stay exactly as they are, so a threshold choice cannot quietly change a second
    MISA setting. ``None`` returns the shipped file unchanged (the default path).
    """
    text = data_file("misa.ini").read_text(encoding="utf-8")
    if unit_repeats is None:
        return text
    values = " ".join(f"{unit}-{repeats}" for unit, repeats in sorted(unit_repeats.items()))
    lines = text.splitlines()
    for index, line in enumerate(lines):
        if line.startswith(_DEFINITION_PREFIX):
            label, _, rest = line.partition(":")
            padding = rest[: len(rest) - len(rest.lstrip())]  # keep the shipped column
            lines[index] = f"{label}:{padding}{values}"
            return "\n".join(lines) + ("\n" if text.endswith("\n") else "")
    raise ValueError(f"the shipped misa.ini has no {_DEFINITION_PREFIX} line")


def misa_script(explicit: str | Path | None = None) -> Path:
    """The ``misa.pl`` to run: an explicit path, else ``misa.pl`` on ``PATH``.

    MISA ships with no licence, so the toolkit does not redistribute it (D-q). The message names
    the one thing that helps — the download line in the docs — instead of reporting a missing
    file the user cannot act on.
    """
    if explicit is not None:
        path = Path(explicit)
        if not path.is_file():
            raise ToolError(f"misa: --set misa={explicit} is not a file — {MISA_DOC}")
        return path
    found = which(_SCRIPT_NAME)
    if found is None:
        raise ToolError(
            f"misa: {_SCRIPT_NAME} not found — MISA {MISA_VERSION} is not bundled with this "
            f"toolkit, because its source zip states no licence (docs/licensing.md). Get it with "
            f"the download line in {MISA_DOC}, put it on PATH, or pass "
            f"--set misa=/path/to/{_SCRIPT_NAME}"
        )
    return Path(found)


def misa_version(script: Path) -> str:
    """The release this script is, read from its own header — and refused if it is not the pin.

    MISA has no version flag. The only place it names its release is the comment block at the top
    (``###Release date: 25/08/20 (version 2.1)``), so that is what the run records, and a
    different release stops the run instead of publishing counts nobody measured — the same rule
    ``repeat_match_version`` applies to MUMmer.
    """
    header = script.read_text(encoding="utf-8", errors="replace")[:2000]
    match = _VERSION_RE.search(header)
    if match is None:
        raise ToolError(
            f"{script}: no '(version x.y)' in its first lines, so the MISA release cannot be "
            f"read — install the pinned {MISA_VERSION} ({MISA_DOC})"
        )
    version = match.group(1)
    if version != MISA_VERSION:
        raise ToolError(
            f"{script}: expected MISA {MISA_VERSION} (the release the golden SSR numbers were "
            f"measured with), got {version} — install the pinned release, or re-measure and "
            f"update MISA_VERSION ({MISA_DOC})"
        )
    return version


def run_misa(
    genome: Genome,
    ctx: RunContext,
    script: str | Path | None = None,
    *,
    unit_repeats: Mapping[int, int] | None = None,
) -> MisaResult:
    """Run MISA on one genome; returns counts and per-SSR positions.

    ``unit_repeats`` replaces the shipped ``misa.ini``'s definition line (M11.20); ``None`` uses
    the shipped defaults. The ini written for this run is the only source of the thresholds the
    result reports, so a value the module accepted cannot be dropped between here and MISA.
    """
    perl = which("perl")
    if perl is None:
        raise ToolError(
            "misa: perl not found — MISA 2.1 is a Perl script; install perl to run this module"
        )
    script_path = misa_script(script)
    # Read, and refuse anything but the pin — so recording it as "pinned" is what ran, not a claim.
    version = misa_version(script_path)
    interpreter_version, _ = tool_version("perl", "perl", perl)
    interpreter = f"perl {interpreter_version}" if interpreter_version else "perl"
    chosen = dict(unit_repeats) if unit_repeats is not None else unit_repeats_from_ini()

    with tempfile.TemporaryDirectory(prefix="misa_") as work:
        workdir = Path(work)
        shutil.copy(script_path, workdir / _SCRIPT_NAME)
        (workdir / "misa.ini").write_text(misa_ini_text(chosen), encoding="utf-8")
        fasta = workdir / "input.fasta"
        fasta.write_text(_to_fasta(genome), encoding="utf-8")

        ctx.tool(
            "misa",
            [perl, str(workdir / _SCRIPT_NAME), fasta.name],
            kind="perl",
            tool_path=str(script_path),
            pinned_version=version,
            interpreter=interpreter,
            cwd=workdir,
        )

        statistics = fasta.with_name(fasta.name + ".statistics")
        if not statistics.is_file():
            raise ToolError("misa: the run produced no .statistics file — see the captured logs")
        total, by_unit, compound = _parse_statistics(statistics.read_text(encoding="utf-8"))
        gff_files = sorted(workdir.glob("*.gff"))
        features = _parse_gff(gff_files[0].read_text(encoding="utf-8")) if gff_files else ()
        return MisaResult(
            total=total,
            by_unit=by_unit,
            compound=compound,
            features=features,
            unit_repeats=chosen,
        )


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
    """Every SSR in MISA's GFF — including the ones its two vocabularies hide.

    MISA writes monomeric SSRs as ``monomeric_repeat`` but every longer unit as a generic
    ``microsatellite``, with the class only inside the Note; a compound cluster is one
    ``repeat_region`` line whose Note lists all of its components. Reading the type column
    alone therefore drops **every** di/tri/… SSR and one component of every cluster — measured
    2026-10-03 on the golden genome: only 30 of the 37 SSRs survived, silently (found by P2-06's
    reconciliation against MISA's own statistics).
    """
    features: list[SsrFeature] = []
    for line in text.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 9:
            continue
        components = _note_components(fields[2], fields[8])
        # MISA marks a compound cluster **two** ways, and one file can carry both: a
        # ``repeat_region`` line listing the cluster's components, or one line per component
        # prefixed ``compound_repeat``. Reading only the first marked every row of the second
        # kind "no" — measured 10/10/2026 (P8B-01): 0 of PX413282's 16 components flagged while
        # MISA's own statistics said 8, and the golden's table flagged 2 of its 6.
        compound = len(components) > 1 or "compound_repeat" in fields[8]
        start, end = int(fields[3]), int(fields[4])
        features.extend(
            SsrFeature(
                unit_size=unit_size,
                motif=motif,
                start=start,
                end=end,
                compound=compound,
            )
            for unit_size, motif in components
        )
    return tuple(features)


def _note_components(feature_type: str, note: str) -> list[tuple[int, str]]:
    """``(unit size, motif)`` for every SSR one GFF line stands for."""
    if feature_type == "region" or "Note=" not in note:
        return []
    tokens = [token.strip() for token in note.split("Note=", 1)[1].split(";", 1)[0].split(",")]
    components: list[tuple[int, str]] = []
    index = 0
    while index < len(tokens):
        token = tokens[index]
        unit_size = _CLASS_UNITS.get(token)
        if unit_size is None and token.endswith("_microsatellite_feature"):
            unit_size = _CLASS_UNITS.get(token.removesuffix("_microsatellite_feature"))
        if unit_size is None:
            index += 1
            continue
        motif = tokens[index + 1] if index + 1 < len(tokens) else ""
        components.append((unit_size, motif))
        index += 2
    return components


def _int_after(text: str, marker: str) -> int:
    for line in text.splitlines():
        if marker in line:
            return int(line.split(marker, 1)[1].strip())
    raise ToolError(f"misa: expected {marker!r} in the statistics output, found none")


def threshold_text(path: str | Path | None = None) -> str:
    """The unit/repeat thresholds in an ini, as MISA writes them.

    The summary used to carry this string hardcoded beside a run that reads the same values from
    the ini; the two could drift. Reading the file makes the reported thresholds the ones MISA
    actually used (M11.17). ``path`` exists for tests and for a caller that supplies the ini.
    """
    return format_thresholds(unit_repeats_from_ini(path))
