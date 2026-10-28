"""MACSE 2.07 — codon-aware CDS alignment, opt-in for the phylogeny module (M11.21).

**Why a second aligner.** MAFFT aligns nucleotides; a codon partition on that alignment is a
partition of columns that *look* like codon positions until an indel shifts the frame. MACSE
aligns the nucleotide sequences **through their amino-acid translation**, so a frameshift is placed
as a codon-level indel rather than as three independent gaps. CGAS offers the same tool the same
way — `--macse` is opt-in there, and it is opt-in here: the default stays MAFFT, which is what every
number in this repository was measured with.

**Obtained, not bundled.** MACSE is a Java program distributed under CeCILL 2.1; the toolkit does
not redistribute the jar. `docs/tools/macse.md` records the bioconda package, its sha256 and the
jar's own sha256. A user supplies `--set macse=/path/to/macse_v2.07.jar` or installs the package,
which puts a `macse` wrapper on ``PATH``.

**The version is read from the banner and checked against the pin.** MACSE has no version flag;
its ``-help`` banner starts ``This is MACSE V2.07``. The alignment algorithm and the output
alphabet have changed between releases, so a different release is refused by name rather than
recorded and used.

**Two characters are not DNA and are replaced at the boundary.** MACSE writes a frameshift as ``!``
and a stop codon as ``*`` in the nucleotide alignment. :class:`plastid.core.align.Alignment` reads
ACGT (with ``N`` and ``-``) and would otherwise carry a character no downstream counter can
classify. Both become ``-`` — the codon is not comparable to an aligned codon, which is exactly
what a gap means here. The substitution is stated in the run's scope and docs; it is not silent.
"""

from __future__ import annotations

import re
import subprocess
import tempfile
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO

from plastid.core.align import Alignment
from plastid.run.context import RunContext
from plastid.run.tools import ToolError, which

__all__ = [
    "MACSE_VERSION",
    "MacseBinary",
    "align_cds",
    "clean_sequence",
    "macse_version",
    "parse_version",
    "resolve_macse",
]

MACSE_VERSION = "2.07"
MACSE_DOC = "docs/tools/macse.md"
_VERSION = re.compile(r"This is MACSE V([0-9][0-9.]*)")
_NON_DNA = str.maketrans({"!": "-", "*": "-"})


@dataclass(frozen=True, slots=True)
class MacseBinary:
    """The command prefix that runs MACSE, and the release its own banner names."""

    argv_prefix: tuple[str, ...]
    version: str


def parse_version(text: str) -> str:
    """The release from MACSE's ``-help`` banner, or a refusal that says what was expected."""
    match = _VERSION.search(text)
    if match is None:
        raise ToolError(
            f"macse: no 'This is MACSE Vx.yy' line in its -help output — install the pinned "
            f"{MACSE_VERSION} ({MACSE_DOC}). Output was: {text[:200]!r}"
        )
    return match.group(1)


def macse_version(argv_prefix: Sequence[str]) -> str:
    """Run ``-help`` and check the banner against the pin."""
    try:
        result = subprocess.run(
            [*argv_prefix, "-help"], capture_output=True, text=True, timeout=60, check=False
        )
    except OSError as error:
        raise ToolError(f"macse: cannot run {argv_prefix[0]}: {error}") from error
    version = parse_version((result.stdout or "") + (result.stderr or ""))
    if version != MACSE_VERSION:
        raise ToolError(
            f"macse: this is MACSE {version}, but the toolkit's numbers and its tests use "
            f"{MACSE_VERSION} ({MACSE_DOC}) — install the pinned release, or re-measure and "
            f"update MACSE_VERSION."
        )
    return version


def resolve_macse(explicit: str | Path | None = None) -> MacseBinary:
    """Find MACSE: an explicit jar/executable, else the conda package's ``macse`` on ``PATH``."""
    if explicit is not None:
        path = Path(explicit)
        if not path.is_file():
            raise ToolError(f"macse: --set macse={explicit} is not a file — {MACSE_DOC}")
        if path.suffix == ".jar":
            java = which("java")
            if java is None:
                raise ToolError(
                    f"macse: {path.name} is a jar and java is not on PATH — install Java 8+ "
                    f"({MACSE_DOC})"
                )
            prefix = (java, "-jar", str(path))
        else:
            prefix = (str(path),)
    else:
        found = which("macse")
        if found is None:
            raise ToolError(
                f"macse: not found — install the pinned release (`conda install -c bioconda "
                f"macse={MACSE_VERSION}`), put its `macse` on PATH, or pass "
                f"--set macse=/path/to/macse_v2.07.jar ({MACSE_DOC})"
            )
        prefix = (found,)
    return MacseBinary(argv_prefix=prefix, version=macse_version(prefix))


def clean_sequence(sequence: str) -> str:
    """MACSE's ``!`` (frameshift) and ``*`` (stop) become gaps; the rest is left as written."""
    return sequence.upper().translate(_NON_DNA)


def align_cds(
    rows: Sequence[tuple[str, str]],
    ctx: RunContext,
    *,
    binary: MacseBinary,
    timeout: int | None = 600,
) -> Alignment:
    """Align ``(name, CDS)`` rows with MACSE and return them in the order given.

    One gene at a time, like the MAFFT path: the module concatenates the per-gene alignments.
    ``!`` and ``*`` are replaced with ``-`` by :func:`clean_sequence`; a row MACSE did not return
    is a named error, never a silent drop.
    """
    if len(rows) < 2:
        raise ToolError(
            f"macse: {len(rows)} sequence(s) is a comparison, not an alignment — the matrix "
            f"rule excludes a gene carried by fewer than two genomes"
        )
    with tempfile.TemporaryDirectory(prefix="macse_") as work:
        workdir = Path(work)
        fasta = workdir / "input.fasta"
        fasta.write_text(_to_fasta(rows), encoding="utf-8")
        through_java = binary.argv_prefix[0] != binary.argv_prefix[-1]
        ctx.tool(
            "macse",
            [
                *binary.argv_prefix,
                "-prog",
                "alignSequences",
                "-seq",
                fasta.name,
                "-out_NT",
                "nt.fasta",
                "-out_AA",
                "aa.fasta",
            ],
            pinned_version=binary.version,
            tool_path=str(binary.argv_prefix[-1]),
            interpreter=binary.argv_prefix[0] if through_java else None,
            cwd=workdir,
            timeout=timeout,
        )
        aligned_path = workdir / "nt.fasta"
        if not aligned_path.is_file():
            raise ToolError("macse: exited 0 but wrote no nt.fasta — see the captured logs")
        aligned = {
            record.id: clean_sequence(str(record.seq))
            for record in SeqIO.parse(str(aligned_path), "fasta")
        }
    missing = [name for name, _ in rows if name not in aligned]
    if missing:
        raise ToolError(f"macse returned no row for {', '.join(missing)}")
    return Alignment(
        names=tuple(name for name, _ in rows),
        sequences=tuple(aligned[name] for name, _ in rows),
    )


def _to_fasta(rows: Sequence[tuple[str, str]], width: int = 60) -> str:
    lines: list[str] = []
    for name, sequence in rows:
        lines.append(f">{name}")
        lines.extend(sequence[i : i + width] for i in range(0, len(sequence), width))
    return "\n".join(lines) + "\n"
