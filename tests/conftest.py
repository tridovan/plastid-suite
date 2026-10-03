"""Shared fixtures.

Golden data is vendored under ``tests/data`` (sha256-pinned in SHA256SUMS) —
tests never touch the network.
"""

from __future__ import annotations

import os
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
def annotator_tools() -> dict[str, str | None]:
    """Chloe and the Julia it runs on — from ``PATH``, or from ``PTS_CHLOE`` / ``PTS_JULIA``.

    Both live outside ``PATH`` on the machine this was written on (Julia came out of a mamba
    environment, Chloe out of a git clone), so the environment is how a machine says where they
    are; a fixture that guessed a path would just be wrong on the next machine.
    """
    return {
        "chloe": os.environ.get("PTS_CHLOE") or which("chloe.jl"),
        "julia": os.environ.get("PTS_JULIA") or which("julia"),
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
        name: which(name, extra=extra) for name in ("blastx", "blastn", "makeblastdb", "aragorn")
    }


@pytest.fixture(scope="session")
def submission_tools() -> dict[str, str | None]:
    """``table2asn`` — from ``PATH`` or from ``PTS_TABLE2ASN`` (it is a single binary fetched by
    hand on the dev machine, so it is usually not on ``PATH``)."""
    return {"table2asn": os.environ.get("PTS_TABLE2ASN") or which("table2asn")}


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
