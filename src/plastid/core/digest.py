"""The SHA-256 of a file on disk — the one way this package fingerprints a file.

A run, a cache entry, a pinned record and an artifact each name the bytes they were made from by
this digest, so every caller has to compute it the same way. The file is read in blocks, so a
large input (a read set, a database) is never held in memory whole.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

__all__ = ["sha256_file"]

_BLOCK_BYTES = 1 << 20


def sha256_file(path: str | Path) -> str:
    """The hex SHA-256 of the file's bytes."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(_BLOCK_BYTES), b""):
            digest.update(block)
    return digest.hexdigest()
