"""The cut that turns the first bytes of a large run into complete read pairs. No network."""

from __future__ import annotations

import gzip
import random
from pathlib import Path

from tools.fetch_assembly_panel import _cut_pairs


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
