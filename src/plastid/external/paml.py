"""PAML's ``yn00`` — pairwise Ka/Ks by the Yang & Nielsen (2000) method (P4-02, decision D4-01).

``yn00`` is not a command-line tool in the usual sense: it takes **no arguments**, reads
``yn00.ctl`` from the working directory, and writes its answer into files whose names come from
that ctl. So the adapter's provenance is the ctl plus the input file, and both are written per
run in a throwaway directory — the same pattern as MAFFT's adapter, for the same reason: a tool
that reads its parameters from a file leaves no argv for the manifest otherwise.

Three measured facts shape this module (2026-10-04, PAML 4.10.10):

* **the version is not a flag.** ``yn00 --version`` is not a thing; run with no ctl it prints
  ``YN00 in paml version 4.10.10, 29 Jan 2026`` on stdout and exits 1. That is the probe, and a
  release other than the pin is refused by name.
* **it reads the alignment three bases at a time and does not check the frame.** A one-base gap
  in a synthetic pair made it report a *spurious stop codon* and exit 255 — luck, not a check
  (a frameshift that lands on no stop would have been scored as if real). The frame check is
  ``core.align.codon_frame_breaks``, called by the module before this adapter is reached.
* **stop codons are its poison.** A terminal stop must be trimmed by the caller; an internal one
  stops the run with ``stop codon TAA in seq # 2``. Both are reported, never repaired silently.

``dN``/``dS`` are read from the ``2YN.dN`` / ``2YN.dS`` matrices and cross-checked against the
Yang & Nielsen line of the run's own output — two code paths inside PAML, and a disagreement
between them is not something to average.
"""

from __future__ import annotations

import math
import re
import tempfile
from dataclasses import dataclass
from pathlib import Path

from plastid.run.context import RunContext
from plastid.run.tools import ToolError, which

__all__ = [
    "PAML_VERSION",
    "Yn00Result",
    "parse_yn00_output",
    "read_matrix",
    "resolve_binary",
    "run_yn00",
    "yn00_version",
]

# The release the toolkit's numbers are measured with. A different PAML is a different estimator
# implementation, so a mismatch is refused rather than noted (the BLAST precedent, §8.15).
PAML_VERSION = "4.10.10"

_BANNER = re.compile(r"YN00 in paml version\s+([\w.]+)")
_YN_SECTION = "(B) Yang & Nielsen (2000) method"

# yn00.ctl — written by this adapter, so it is part of the provenance. icode 0 is the standard
# genetic code: yn00 uses the code for its codon table, and NCBI table 11 differs from table 1
# only in initiation codons, which are not part of an alignment's interior.
CTL = "seqfile = input.dat\noutfile = yn00.out\nverbose = 0\nicode = 0\nweighting = 0\n"


@dataclass(frozen=True, slots=True)
class Yn00Result:
    """One pairwise estimate: the Yang & Nielsen numbers, and what they were computed over."""

    dN: float  # Ka — nonsynonymous substitutions per nonsynonymous site
    dS: float  # Ks — synonymous substitutions per synonymous site
    omega: float | None  # yn00's own dN/dS; None when it writes nan (dS = 0)
    kappa: float | None  # the transition/transversion rate ratio it estimated
    t: float | None  # divergence at the nucleotide level
    sites_synonymous: float | None  # S — the denominator of dS, used for the pooled row
    sites_nonsynonymous: float | None  # N — the denominator of dN
    codons: int | None
    raw: str  # the run's own output file, kept so a number can be traced


def resolve_binary(*, directory: str | Path | None = None) -> str:
    """Find ``yn00``: a given directory, then ``PATH``."""
    if directory is not None:
        candidate = Path(directory) / "yn00"
        if not candidate.is_file():
            raise ToolError(f"yn00 not found in {directory} — see docs/tools/paml.md")
        return str(candidate)
    found = which("yn00")
    if found is None:
        raise ToolError(
            f"yn00 not found — install the pinned PAML {PAML_VERSION}, or point the run at one "
            f"with --set paml=<directory>. See docs/tools/paml.md."
        )
    return found


def yn00_version(binary: str) -> str:
    """The PAML release this binary is, read from its no-ctl banner and checked against the pin."""
    import subprocess

    with tempfile.TemporaryDirectory(prefix="yn00_version_") as work:
        try:
            result = subprocess.run(
                [binary], capture_output=True, text=True, timeout=30, check=False, cwd=work
            )
        except OSError as error:
            raise ToolError(f"yn00: cannot run {binary}: {error}") from error
    match = _BANNER.search(result.stdout or "")
    if match is None:
        raise ToolError(
            f"{binary}: cannot read a version from its banner — a step whose version nobody can "
            f"name is not part of the toolkit. Output was: {(result.stdout or '')[:200]!r}"
        )
    version = match.group(1)
    if version != PAML_VERSION:
        raise ToolError(
            f"{binary} is PAML {version}, but the toolkit's numbers were measured with "
            f"{PAML_VERSION} (docs/tools/paml.md) — install the pinned build, or re-measure and "
            f"update PAML_VERSION."
        )
    return version


def run_yn00(
    name_a: str,
    sequence_a: str,
    name_b: str,
    sequence_b: str,
    ctx: RunContext,
    *,
    binary: str | None = None,
    directory: str | Path | None = None,
    timeout: int | None = None,
) -> Yn00Result:
    """One pairwise Ka/Ks estimate for two codon-aligned, stop-free, equal-length sequences."""
    for name in (name_a, name_b):
        if not name or any(character.isspace() for character in name):
            raise ValueError(f"PAML takes a name with no whitespace, got {name!r}")
    if name_a == name_b:
        raise ValueError("yn00 cannot compare a sequence with itself — give the two rows names")
    if len(sequence_a) != len(sequence_b):
        raise ValueError(
            f"yn00 needs an alignment: rows are {len(sequence_a)} and {len(sequence_b)} bases"
        )
    if len(sequence_a) % 3:
        raise ValueError(
            f"yn00 reads codons: {len(sequence_a)} aligned bases is not a whole number of codons"
        )

    # ``binary`` is for a caller that already resolved and verified it (the module scores one
    # gene after another); the version is still checked here, because this adapter's contract is
    # that nothing reaches yn00 unnamed.
    resolved = binary or resolve_binary(directory=directory)
    version = yn00_version(resolved)
    with tempfile.TemporaryDirectory(prefix="yn00_") as work:
        workdir = Path(work)
        (workdir / "input.dat").write_text(
            f"  2  {len(sequence_a)}\n{name_a}\n{sequence_a}\n{name_b}\n{sequence_b}\n",
            encoding="utf-8",
        )
        (workdir / "yn00.ctl").write_text(CTL, encoding="utf-8")
        ctx.tool(
            "yn00",
            [resolved],
            tool_path=resolved,
            pinned_version=version,
            cwd=workdir,
            timeout=timeout,
        )
        out_path = workdir / "yn00.out"
        if not out_path.is_file() or not out_path.read_text(encoding="utf-8").strip():
            raise ToolError(
                f"yn00 exited without writing a result for {name_a} vs {name_b} — "
                f"its output file is empty"
            )
        raw = out_path.read_text(encoding="utf-8")
        matrix_dn = read_matrix(workdir / "2YN.dN", name_a, name_b)
        matrix_ds = read_matrix(workdir / "2YN.dS", name_a, name_b)

    result = parse_yn00_output(raw)
    for value, matrix, label in ((result.dN, matrix_dn, "dN"), (result.dS, matrix_ds, "dS")):
        if abs(value - matrix) > 1e-4:
            raise ToolError(
                f"yn00's own output and its {label} matrix disagree ({value} vs {matrix}) — "
                f"two code paths in the tool, and nothing here averages them"
            )
    return result


def read_matrix(path: Path, name_a: str, name_b: str) -> float:
    """One value from a PAML ``2YN.*`` matrix file (a square, name-labelled table)."""
    if not path.is_file():
        raise ToolError(f"yn00 wrote no {path.name} — the run did not complete")
    lines = [line for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(lines) < 3:
        raise ToolError(f"{path.name} does not hold a matrix: {lines!r}")
    header = lines[1].split()
    if header != [name_a]:
        raise ToolError(
            f"{path.name}: expected the first sequence to be {name_a!r}, found {header!r}"
        )
    row = lines[2].split()
    if len(row) != 2 or row[0] != name_b:
        raise ToolError(f"{path.name}: expected a row for {name_b!r}, found {row!r}")
    value = _number(row[1])
    if value is None:
        raise ToolError(f"{path.name}: the value for {name_b!r} is nan")
    return value


def parse_yn00_output(text: str) -> Yn00Result:
    """The Yang & Nielsen line of yn00's own output.

    The line carries more than the two matrices do — the estimated kappa, the number of
    synonymous and nonsynonymous sites, and the two standard errors — so it is parsed rather
    than only the two numbers the matrices hold.
    """
    section = text.split(_YN_SECTION, 1)
    if len(section) < 2:
        raise ToolError("yn00's output has no Yang & Nielsen section — nothing to read")
    body = section[1].split("(C)", 1)[0]
    row: list[str] | None = None
    for line in body.splitlines():
        fields = line.split()
        if len(fields) >= 13 and fields[0].isdigit() and fields[1].isdigit():
            row = fields
            break
    if row is None:
        raise ToolError(
            "yn00's Yang & Nielsen section holds no pairwise row — the run did not produce a "
            "result for this pair"
        )
    # i j S N t kappa omega dN +- SE dS +- SE
    dN, dS = _number(row[7]), _number(row[10])
    if dN is None or dS is None:
        raise ToolError(f"yn00 reported nan for dN/dS itself: {row!r}")
    return Yn00Result(
        dN=dN,
        dS=dS,
        omega=_number(row[6]),
        kappa=_number(row[5]),
        t=_number(row[4]),
        sites_synonymous=_number(row[2]),
        sites_nonsynonymous=_number(row[3]),
        codons=_codons_from(text),
        raw=text,
    )


def _codons_from(text: str) -> int | None:
    """``ls`` from yn00's own header — measured 04/10/2026: it counts **codons**, not bases
    (the golden *rbcL* pair, 1,458 aligned bases, prints ``ls = 486``)."""
    match = re.search(r"^ns =\s*\d+\s+ls =\s*(\d+)", text, re.MULTILINE)
    return int(match.group(1)) if match else None


def _number(token: str) -> float:
    """A PAML number: ``nan``/``-nan`` become ``None``, and ``-0.0000`` becomes ``0.0``.

    Both matter for the same reason: ``-nan`` is what yn00 prints for dN/dS when dS is zero, and
    a ``-0.0000`` that reaches a table is a synonym for zero wearing a sign that will be read as
    "less than zero" in a figure.
    """
    value = float(token)
    if math.isnan(value):
        return None
    return 0.0 if value == 0 else value
