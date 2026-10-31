"""The panel's fetch tool: the cut of a large run, and the pins of the graphs. No network."""

from __future__ import annotations

import gzip
import random
from pathlib import Path

import pytest
from tools.fetch_assembly_panel import _cut_pairs, blob_sha1, published_files


def _reads(count: int, mate: int, seed: int) -> bytes:
    source = random.Random(seed)
    lines = []
    for index in range(count):
        bases = "".join(source.choice("ACGT") for _ in range(100))
        lines.append(f"@r{index}/{mate}\n{bases}\n+\n{'I' * 100}\n")
    return "".join(lines).encode()


def _prefix(path: Path, content: bytes, share: float) -> None:
    compressed = gzip.compress(content, mtime=0)
    path.write_bytes(compressed[: int(len(compressed) * share)])


def test_the_cut_keeps_only_pairs_complete_in_both_files(tmp_path: Path) -> None:
    first, second = _reads(4_000, 1, 1), _reads(4_000, 2, 2)
    _prefix(tmp_path / "1.prefix", first, 0.5)
    _prefix(tmp_path / "2.prefix", second, 0.7)
    pairs = _cut_pairs(
        tmp_path / "1.prefix", tmp_path / "2.prefix", tmp_path / "1.gz", tmp_path / "2.gz"
    )
    assert 0 < pairs < 4_000
    kept_1 = gzip.decompress((tmp_path / "1.gz").read_bytes()).decode().splitlines()
    kept_2 = gzip.decompress((tmp_path / "2.gz").read_bytes()).decode().splitlines()
    assert len(kept_1) == len(kept_2) == 4 * pairs
    # Every record is whole, the two files stay in step, and nothing was reordered.
    assert first.decode().splitlines()[: 4 * pairs] == kept_1
    assert second.decode().splitlines()[: 4 * pairs] == kept_2
    assert all(len(line) == 100 for line in kept_1[1::4])
    assert [line.split("/")[0] for line in kept_1[0::4]] == [
        line.split("/")[0] for line in kept_2[0::4]
    ]


def test_the_cut_gives_the_same_bytes_twice(tmp_path: Path) -> None:
    _prefix(tmp_path / "1.prefix", _reads(1_000, 1, 3), 0.8)
    _prefix(tmp_path / "2.prefix", _reads(1_000, 2, 4), 0.8)
    outputs = []
    for label in ("a", "b"):
        out_1, out_2 = tmp_path / f"{label}1.gz", tmp_path / f"{label}2.gz"
        _cut_pairs(tmp_path / "1.prefix", tmp_path / "2.prefix", out_1, out_2)
        outputs.append((out_1.read_bytes(), out_2.read_bytes()))
    assert outputs[0] == outputs[1]


def test_a_file_is_pinned_by_the_hash_git_gives_it() -> None:
    # `printf 'hello\n' | git hash-object --stdin`
    assert blob_sha1(b"hello\n") == "ce013625030ba8dba906f756967f9e9ca394464a"


def test_only_the_graph_paths_and_log_of_the_two_run_sets_are_kept() -> None:
    results = "res/GetOrganelle-plastome"
    tree = {
        "tree": [
            {"path": f"{results}/auto--A-plastome/x.K127.complete.graph1.1.path_sequence.fasta"},
            {"path": f"{results}/auto--A-plastome/x.K127.complete.graph1.selected_graph.gfa"},
            {"path": f"{results}/customized--A-plastome/get_org.log.txt"},
            {"path": f"{results}/auto--A-plastome/filtered_K127.assembly_graph.fastg.gz"},
            {"path": f"{results}/w0.60--A-plastome/get_org.log.txt"},
            {"path": "res/NOVOPlasty-plastome/auto--A-plastome/get_org.log.txt"},
            {"path": f"{results}/plant.plastome.summary.corrected.tab"},
        ]
    }
    for index, entry in enumerate(tree["tree"]):
        entry["sha"] = str(index)
    assert published_files(tree) == [
        ("auto--A-plastome/x.K127.complete.graph1.1.path_sequence.fasta", "0"),
        ("auto--A-plastome/x.K127.complete.graph1.selected_graph.gfa", "1"),
        ("customized--A-plastome/get_org.log.txt", "2"),
    ]


def test_a_truncated_listing_is_refused() -> None:
    with pytest.raises(OSError, match="truncated"):
        published_files({"truncated": True, "tree": []})
