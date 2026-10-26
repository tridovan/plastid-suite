"""The committed stand-in dataset definitions — guarded offline so they cannot drift.

What these guard: the exact measured record set, the mito/plastome split,
and the five GenBank records left out because their sequence is byte-identical to a RefSeq
already listed.
"""

from __future__ import annotations

from pathlib import Path

from plastid.core.samples import load_groups

DATA = Path(__file__).resolve().parents[2] / "data"
PLASTOMES = DATA / "oldenlandia-coleus.yaml"
MITO = DATA / "oldenlandia-coleus-mito.yaml"

EXCLUDED_DUPLICATES = {"MT767006.1", "MT767008.1", "MT473748.1", "MW538954.1", "OP611428.1"}


def test_plastome_dataset_holds_the_measured_records() -> None:
    groups = {group.name: group.accessions for group in load_groups(PLASTOMES)}
    assert list(groups) == ["Oldenlandia", "Coleus"]
    assert groups["Oldenlandia"] == (
        "NC_057985.1",
        "PP078674.1",
        "PQ325515.1",
        "NC_057983.1",
        "PX977793.1",
    )
    assert groups["Coleus"] == (
        "PX923121.1",
        "PZ428230.1",
        "NC_069583.1",
        "AP045470.1",
        "NC_058328.1",
        "NC_057679.1",
    )


def test_mito_dataset_stays_separate_from_the_plastomes() -> None:
    groups = {group.name: group.accessions for group in load_groups(MITO)}
    assert groups == {"Oldenlandia": ("OX459128.1",), "Coleus": ("AP045469.1",)}


def test_duplicate_genbank_records_are_never_listed() -> None:
    listed = {
        accession
        for path in (PLASTOMES, MITO)
        for group in load_groups(path)
        for accession in group.accessions
    }
    assert not (listed & EXCLUDED_DUPLICATES)
    assert len(listed) == 13
