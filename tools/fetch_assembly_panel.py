"""Fetch the public reads and published graphs of the assembly panel, verified.

Each read set of ``tests/golden/test_assembly_panel.py`` comes from ENA. Files are fetched one at a
time, because ENA answered parallel requests with 403 and cut others short (measured 31 October
2026), and a body is only ever kept after a successful status, so an error page cannot pass as
reads. A whole file is checked against ENA's md5. A set far larger than a run needs is taken as its
first bytes, cut to the last complete pair: Jin et al. 2020 likewise used the first reads of a large
run, and GetOrganelle stops using reads past 500× anyway.

The graphs are what the authors of that comparison deposited for every run: selected graph, paths
and log, about 30 MB, for ``tests/golden/test_assembly_graph_panel.py``. They come from one pinned
commit, whose tree names each file's git blob, so each file is checked against that hash.

Run from the product repository root, on a machine with the disk for it:

    uv run python tools/fetch_assembly_panel.py
    uv run python tools/fetch_assembly_panel.py --only SRR5602589
    uv run python tools/fetch_assembly_panel.py --graphs-only

No test runs this script, and no test touches the network. It is how a checkout rebuilds its panel.
Exit codes: 0 when every file arrived and was verified, 1 when one did not.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import sys
import time
import urllib.error
import urllib.request
import zlib
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

ENA = "https://ftp.sra.ebi.ac.uk/vol1/fastq"
PREFIX_BYTES = 350_000_000
BLOCK = 1 << 20
ATTEMPTS = 40

# The repository of Jin et al. 2020. Results of the 50 runs with the published settings (``auto--``)
# and of the 8 re-run with settings chosen per sample (``customized--``).
JIN_REPOSITORY = "Kinggerm/GetOrganelleComparison"
JIN_COMMIT = "b58f4a6171612a740ee232fe17061aa4353be961"
_JIN_RESULTS = "res/GetOrganelle-plastome/"
_JIN_KEPT = re.compile(
    r"(auto|customized)--[^/]+/"
    r"(.+\.selected_graph\.gfa|.+\.path_sequence\.fasta|get_org\.log\.txt)$"
)


@dataclass(frozen=True)
class WholeRun:
    """A run small enough to take whole; ``md5`` are ENA's for the two files."""

    run: str
    path: str
    md5: tuple[str, str]


@dataclass(frozen=True)
class RunPrefix:
    """The first bytes of a large run, cut to complete pairs. ENA's md5 cannot check a part.

    The pin is the sha256 of the decompressed reads, so it names the reads and not a gzip setting.
    """

    run: str
    path: str
    pairs: int
    sha256: tuple[str, str]


PANEL: tuple[WholeRun | RunPrefix, ...] = (
    WholeRun(
        "SRR5602589",
        "SRR560/009/SRR5602589",
        ("e5c23901c0610cf48285ae54c3c89f35", "665c4594259f933f08f028e8bcee2de8"),
    ),
    WholeRun(
        "SRR5602602",
        "SRR560/002/SRR5602602",
        ("ea8833de2ec5b307f295226d2864923e", "70c78964d812fb41a8e0643fcb23e2a5"),
    ),
    WholeRun(
        "SRR2057084",
        "SRR205/004/SRR2057084",
        ("662878609ee3838b2e9c21a50ddd0568", "faae1482a7d58dcac44ae60b7db76111"),
    ),
    RunPrefix(
        "SRR2037123",
        "SRR203/003/SRR2037123",
        3_816_250,
        (
            "b988b5bc38851adcece11bc06a9034a46253b11e03cd3ae08c68f346715b9f9d",
            "7a3cdb262343129b93422cfce668600dc3e6f1fc70b24f8f9b7ba18f9cbfb1f7",
        ),
    ),
    RunPrefix(
        "DRR173054",
        "DRR173/DRR173054",
        4_985_864,
        (
            "f09aeed0399f79ac288ee24c1eae234d37b179bb64a500a020d8a11655852c91",
            "03833b20ca2242b604874d574777dd3a326e0f99fe5eabd2778fa4220020ce6b",
        ),
    ),
)


def _fetch(url: str, target: Path, *, limit: int | None = None) -> None:
    """Fetch ``url`` into ``target``, resuming after a cut. An HTTP error writes nothing."""
    for attempt in range(1, ATTEMPTS + 1):
        have = target.stat().st_size if target.exists() else 0
        if limit is not None and have >= limit:
            return
        end = "" if limit is None else str(limit - 1)
        request = urllib.request.Request(url, headers={"Range": f"bytes={have}-{end}"})
        try:
            with urllib.request.urlopen(request, timeout=600) as response:
                if response.status != 206:
                    raise urllib.error.HTTPError(url, response.status, "no range", {}, None)
                with target.open("ab") as handle:
                    while block := response.read(BLOCK):
                        handle.write(block)
            if limit is None:
                return
        except urllib.error.HTTPError as error:
            if error.code == 416:  # the file is already complete
                return
            print(f"{target.name}: attempt {attempt}: {error}", file=sys.stderr)
            time.sleep(8)
        except (urllib.error.URLError, TimeoutError, OSError) as error:
            print(f"{target.name}: attempt {attempt}: {error}", file=sys.stderr)
            time.sleep(8)
    raise OSError(f"{url}: not complete after {ATTEMPTS} attempts")


def _md5(path: Path) -> str:
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as handle:
        while block := handle.read(BLOCK):
            digest.update(block)
    return digest.hexdigest()


def _decompressed_prefix(path: Path) -> bytes:
    """What a gzip prefix decompresses to. The file is cut, so its end is simply where it stops."""
    decompressor = zlib.decompressobj(31)
    pieces: list[bytes] = []
    with path.open("rb") as handle:
        while block := handle.read(BLOCK):
            pieces.append(decompressor.decompress(block))
            while decompressor.unused_data:
                rest = decompressor.unused_data
                decompressor = zlib.decompressobj(31)
                pieces.append(decompressor.decompress(rest))
    return b"".join(pieces)


def _cut_pairs(first: Path, second: Path, out_1: Path, out_2: Path) -> int:
    """Write the complete pairs both prefixes hold, and return how many."""
    lines_1 = _decompressed_prefix(first).split(b"\n")
    lines_2 = _decompressed_prefix(second).split(b"\n")
    pairs = (min(len(lines_1), len(lines_2)) - 1) // 4
    for lines, out in ((lines_1, out_1), (lines_2, out_2)):
        # No timestamp in the header, so the same reads give the same bytes.
        out.write_bytes(gzip.compress(b"\n".join(lines[: 4 * pairs]) + b"\n", 2, mtime=0))
    return pairs


def fetch_whole(item: WholeRun, reads: Path) -> None:
    for mate, expected in zip((1, 2), item.md5, strict=True):
        target = reads / f"{item.run}_{mate}.fastq.gz"
        if target.exists() and _md5(target) == expected:
            print(f"{target.name}: already here, md5 verified")
            continue
        _fetch(f"{ENA}/{item.path}/{item.run}_{mate}.fastq.gz", target)
        found = _md5(target)
        if found != expected:
            raise OSError(f"{target}: md5 {found}, ENA says {expected}")
        print(f"{target.name}: {target.stat().st_size} bytes, md5 verified")


def _reads_digest(path: Path) -> str:
    """The sha256 of a gzip's decompressed content."""
    digest = hashlib.sha256()
    with gzip.open(path, "rb") as handle:
        while block := handle.read(BLOCK):
            digest.update(block)
    return digest.hexdigest()


def _verify_prefix(item: RunPrefix, files: list[Path]) -> None:
    for file, expected in zip(files, item.sha256, strict=True):
        found = _reads_digest(file)
        if found != expected:
            raise OSError(f"{file}: reads sha256 {found}, pinned {expected}")


def fetch_prefix(item: RunPrefix, reads: Path) -> None:
    out = [reads / f"{item.run}_{mate}.fastq.gz" for mate in (1, 2)]
    if all(path.exists() for path in out):
        _verify_prefix(item, out)
        print(f"{item.run}: already here, reads sha256 verified")
        return
    parts = [reads / f"{item.run}_{mate}.prefix" for mate in (1, 2)]
    for mate, part in zip((1, 2), parts, strict=True):
        _fetch(f"{ENA}/{item.path}/{item.run}_{mate}.fastq.gz", part, limit=PREFIX_BYTES)
        if part.read_bytes()[:2] != b"\x1f\x8b":
            raise OSError(f"{part}: not gzip, so not what was asked for")
    pairs = _cut_pairs(parts[0], parts[1], out[0], out[1])
    for part in parts:
        part.unlink()
    if pairs != item.pairs:
        raise OSError(f"{item.run}: cut to {pairs} pairs, pinned {item.pairs}")
    _verify_prefix(item, out)
    print(f"{item.run}: the first {pairs} pairs, reads sha256 verified")


def blob_sha1(data: bytes) -> str:
    """The hash git names a file's content by."""
    return hashlib.sha1(f"blob {len(data)}\0".encode() + data, usedforsecurity=False).hexdigest()


def published_files(tree: Mapping[str, Any]) -> list[tuple[str, str]]:
    """``(path under the results, blob hash)`` of each file kept, from a git tree listing."""
    if tree.get("truncated"):
        raise OSError("the repository's tree listing came back truncated")
    kept = []
    for entry in tree["tree"]:
        path = entry["path"]
        if path.startswith(_JIN_RESULTS) and _JIN_KEPT.fullmatch(path[len(_JIN_RESULTS) :]):
            kept.append((path[len(_JIN_RESULTS) :], entry["sha"]))
    return sorted(kept)


def _read_url(url: str) -> bytes:
    with urllib.request.urlopen(url, timeout=120) as response:
        return response.read()


def fetch_graphs(out: Path) -> int:
    """Fetch the published graphs, one file at a time. Returns how many did not verify."""
    tree = json.loads(
        _read_url(
            f"https://api.github.com/repos/{JIN_REPOSITORY}/git/trees/{JIN_COMMIT}?recursive=1"
        )
    )
    failed = 0
    for name, blob in published_files(tree):
        target = out / name
        if target.is_file() and blob_sha1(target.read_bytes()) == blob:
            continue
        url = (
            f"https://raw.githubusercontent.com/{JIN_REPOSITORY}/{JIN_COMMIT}/{_JIN_RESULTS}{name}"
        )
        try:
            data = _read_url(url)
        except (urllib.error.URLError, TimeoutError, OSError) as error:
            print(f"{name}: {error}", file=sys.stderr)
            failed += 1
            continue
        if blob_sha1(data) != blob:
            print(f"{name}: not the blob {blob} the pinned commit names", file=sys.stderr)
            failed += 1
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(data)
    return failed


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", default=".panel/assembly/reads", help="where the reads go")
    parser.add_argument("--only", default=None, help="one run accession")
    parser.add_argument(
        "--graphs", default=".panel/assembly/jin2020", help="where the published graphs go"
    )
    parser.add_argument("--graphs-only", action="store_true", help="fetch no reads")
    args = parser.parse_args()
    if args.graphs_only:
        return 1 if fetch_graphs(Path(args.graphs)) else 0
    reads = Path(args.out)
    reads.mkdir(parents=True, exist_ok=True)
    failed = 0 if args.only else fetch_graphs(Path(args.graphs))
    for item in PANEL:
        if args.only and item.run != args.only:
            continue
        try:
            if isinstance(item, WholeRun):
                fetch_whole(item, reads)
            else:
                fetch_prefix(item, reads)
        except OSError as error:
            print(f"{item.run}: {error}", file=sys.stderr)
            failed += 1
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
