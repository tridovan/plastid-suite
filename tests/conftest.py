"""Shared fixtures.

Golden data is vendored under ``tests/data`` (sha256-pinned in SHA256SUMS) —
tests never touch the network.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from plastid.core.genbank import read_genbank, read_genbank_records
from plastid.core.model import Genome

DATA_DIR = Path(__file__).parent / "data"


@pytest.fixture(scope="session")
def data_dir() -> Path:
    return DATA_DIR


@pytest.fixture(scope="session")
def golden_genome() -> Genome:
    """Iris domestica NC_050833.1 — the golden sample (PROJECT.md §6)."""
    return read_genbank(DATA_DIR / "NC_050833.1.gb")


@pytest.fixture(scope="session")
def edge_genomes() -> dict[str, Genome]:
    """Known-hard genomes from PROJECT.md §8–§9: two silently-wrong IR cases, one that
    starts inside an IR copy, all with their real GenBank annotations."""
    records = read_genbank_records(DATA_DIR / "iridaceae_edge_cases.gb")
    return {genome.accession: genome for genome in records}
