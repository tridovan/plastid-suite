"""Shared fixtures.

Golden data is vendored under ``tests/data`` (sha256-pinned in SHA256SUMS) —
tests never touch the network.
"""

from __future__ import annotations

import os
import stat
from pathlib import Path

import pytest

from plastid.core.genbank import read_genbank, read_genbank_records
from plastid.core.model import Genome
from plastid.run.tools import which

DATA_DIR = Path(__file__).parent / "data"


@pytest.fixture(scope="session")
def data_dir() -> Path:
    return DATA_DIR


@pytest.fixture(scope="session")
def golden_genome() -> Genome:
    """Iris domestica NC_050833.1 — the golden sample."""
    return read_genbank(DATA_DIR / "NC_050833.1.gb")


@pytest.fixture(scope="session")
def tools() -> dict[str, str | None]:
    """External tools on this machine — name → path (None when missing).

    ``misa.pl`` is not bundled (D-q: its source zip states no licence, docs/tools/misa.md), so it
    is looked for the way the module looks for it — on ``PATH``, or in ``PTS_MISA_DIR`` on a
    machine that keeps its tools outside ``PATH`` like the other tier-B tools.
    """
    directory = os.environ.get("PTS_MISA_DIR")
    extra = [Path(directory)] if directory else []
    return {
        "perl": which("perl"),
        "misa.pl": which("misa.pl", extra=extra),
        "mafft": which("mafft"),
    }


@pytest.fixture(scope="session")
def blast_tools() -> dict[str, str | None]:
    """BLAST+ — from ``PATH``, or from ``PTS_BLAST_DIR``.

    Its own fixture rather than an entry in ``tools``: that fixture's key set is asserted by a
    test, and the dev machine keeps BLAST+ in a conda environment that is not on ``PATH``
    (docs/tools/blast.md).
    """
    directory = os.environ.get("PTS_BLAST_DIR")
    extra = [Path(directory)] if directory else []
    return {
        name: which(name, extra=extra)
        for name in ("blastx", "blastn", "blastp", "makeblastdb", "aragorn")
    }


@pytest.fixture(scope="session")
def submission_tools() -> dict[str, str | None]:
    """``table2asn`` — from ``PATH`` or from ``PTS_TABLE2ASN`` (it is a single binary fetched by
    hand on the dev machine, so it is usually not on ``PATH``)."""
    return {"table2asn": os.environ.get("PTS_TABLE2ASN") or which("table2asn")}


@pytest.fixture(scope="session")
def paml_tools() -> dict[str, str | None]:
    """PAML's ``yn00`` — from ``PATH`` or from ``PTS_PAML_DIR``."""
    directory = os.environ.get("PTS_PAML_DIR")
    extra = [Path(directory)] if directory else []
    return {"yn00": which("yn00", extra=extra)}


@pytest.fixture(scope="session")
def mummer_tools() -> dict[str, str | None]:
    """MUMmer 4 — from ``PATH`` or from ``PTS_MUMMER_DIR`` (the binaries ship together, and on
    the dev machine they live in a conda environment that is not on ``PATH``).

    ``repeat-match`` is here rather than in ``tools`` for the same reason it exists at all: it is
    the **same package** as ``nucmer`` — measured 10/10/2026, ``conda-meta`` names
    ``mummer4-4.0.1`` for both — so one directory answers for all four, and a mixed installation
    is what ``repeat_match_version`` refuses.
    """
    directory = os.environ.get("PTS_MUMMER_DIR")
    extra = [Path(directory)] if directory else []
    return {
        name: which(name, extra=extra)
        for name in ("nucmer", "delta-filter", "show-coords", "repeat-match")
    }


@pytest.fixture(scope="session")
def vmatch_tools() -> dict[str, str | None]:
    """Vmatch 2.3.1 — from ``PATH`` or from ``PTS_VMATCH_DIR``.

    Its own fixture, and usually absent: Vmatch is **academic-only**, so it is
    not part of the minimum install and the module that publishes the count runs on
    ``repeat-match`` instead. It is the *second engine* that count is checked against — the way
    MISA and cpstools both give the golden's 37 SSRs — so the check skips with a named reason
    where it is not installed, rather than pretending one engine is two.
    """
    directory = os.environ.get("PTS_VMATCH_DIR")
    extra = [Path(directory)] if directory else []
    return {name: which(name, extra=extra) for name in ("mkvtree", "vmatch")}


def make_mummer_stub(directory: Path, *, delta: str, coords: str, version: str = "4.0.1") -> Path:
    """Three stand-in binaries that answer the version probe and print canned output.

    ``nucmer`` writes ``out.delta``; ``delta-filter`` copies its input; ``show-coords`` prints
    ``coords``. The point is the adapter's contract with the tools — the argv, the version check
    and the two readings of one run — not MUMmer's algorithm.
    """
    directory.mkdir(parents=True, exist_ok=True)

    def script(name: str, body: str) -> Path:
        path = directory / name
        path.write_text(
            "#!/bin/sh\n"
            'case "$1" in\n'
            f"  --version|-V|-version) echo '{version}'; exit 0;;\n"
            "esac\n" + body,
            encoding="utf-8",
        )
        path.chmod(path.stat().st_mode | stat.S_IXUSR)
        return path

    script("nucmer", f"cat > out.delta <<'DELTA'\n{delta}DELTA\n")
    script("delta-filter", 'cat "$2" 2>/dev/null || cat out.delta\n')
    script("show-coords", f"cat <<'COORDS'\n{coords}COORDS\n")
    return directory


def require_tools(tools: dict[str, str | None], *names: str) -> None:
    """Skip with a named reason — or fail when PTS_REQUIRE_TOOLS=1 (tool-equipped machines)."""
    missing = [name for name in names if not tools.get(name)]
    if not missing:
        return
    available = sorted(name for name, path in tools.items() if path)
    reason = f"missing external tool(s): {', '.join(missing)}; available: {available or 'none'}"
    if os.environ.get("PTS_REQUIRE_TOOLS") == "1":
        pytest.fail(reason)
    pytest.skip(reason)


def _has_iqtree(directory: Path) -> bool:
    try:
        return any((directory / name).is_file() for name in ("iqtree2", "iqtree"))
    except OSError:
        # An unreadable parent (a non-root CI runner and /root) raises PermissionError from
        # pathlib instead of returning False; that is "not here", not a failure.
        return False


def iqtree_directory() -> Path:
    """IQ-TREE's directory: `PTS_IQTREE_DIR`, then PATH, then the dev toolbox layout.

    CI installs iqtree2 into the pinned environment and puts it on PATH; the host keeps it
    outside PATH in `/root/toolbox`, which only root can read. Probing that path as another user
    raised PermissionError in CI, so every candidate probe is guarded and PATH comes first.
    """
    configured = os.environ.get("PTS_IQTREE_DIR")
    if configured and _has_iqtree(Path(configured)):
        return Path(configured)
    for name in ("iqtree2", "iqtree"):
        found = which(name)
        if found:
            return Path(found).parent
    fallback = Path("/root/toolbox/env/iqtree/bin")
    if _has_iqtree(fallback):
        return fallback
    reason = "iqtree2 is not installed; set PTS_IQTREE_DIR to the directory holding it"
    if os.environ.get("PTS_REQUIRE_TOOLS") == "1":
        pytest.fail(reason)
    pytest.skip(reason)


@pytest.fixture(scope="session")
def edge_genomes() -> dict[str, Genome]:
    """Known-hard genomes: two silently-wrong IR cases, one that
    starts inside an IR copy, all with their real GenBank annotations."""
    records = read_genbank_records(DATA_DIR / "iridaceae_edge_cases.gb")
    return {genome.accession: genome for genome in records}


def require_r(reason: str) -> None:
    """Skip a check that needs R — or fail under PTS_REQUIRE_TOOLS=1.

    Measured 28/10/2026: the golden's R checks called ``pytest.skip`` directly, so CI (which has
    no R) reported green while every figure script went unrun — the exact silent-skip failure
    ``PTS_REQUIRE_TOOLS`` exists to prevent. A missing plotting engine is a missing tool like any
    other, and this is the one door the golden uses to say so.
    """
    if os.environ.get("PTS_REQUIRE_TOOLS") == "1":
        pytest.fail(reason)
    pytest.skip(reason)
