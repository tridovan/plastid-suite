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
    """Iris domestica NC_050833.1 — the golden sample (PROJECT.md §6)."""
    return read_genbank(DATA_DIR / "NC_050833.1.gb")


MISA_VENDOR_DIR = Path(__file__).parents[1] / "src" / "plastid" / "vendor" / "misa"


@pytest.fixture(scope="session")
def tools() -> dict[str, str | None]:
    """External tools on this machine — name → path (None when missing)."""
    vendored_misa = MISA_VENDOR_DIR / "misa.pl"
    return {
        "perl": which("perl"),
        "misa.pl": which("misa.pl", extra=[MISA_VENDOR_DIR])
        or (str(vendored_misa) if vendored_misa.is_file() else None),
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

    Its own fixture, and usually absent: Vmatch is **academic-only** (PROJECT.md §7.8), so it is
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


@pytest.fixture(scope="session")
def edge_genomes() -> dict[str, Genome]:
    """Known-hard genomes from PROJECT.md §8–§9: two silently-wrong IR cases, one that
    starts inside an IR copy, all with their real GenBank annotations."""
    records = read_genbank_records(DATA_DIR / "iridaceae_edge_cases.gb")
    return {genome.accession: genome for genome in records}
