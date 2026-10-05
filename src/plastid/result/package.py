"""The deterministic OOXML container, shared by the Word and Excel exporters.

Both are zip packages whose entries carry the wall-clock time, so the same content written twice
would differ. Rebuilding the archive here — entries sorted, timestamps fixed, attributes pinned —
is what makes "regenerating the document produces an identical file" true rather than aspirational.
The policy is the one the Word export measured first (Phase 2, M2.2); the Excel export is its
second real user, which is when it becomes a function instead of a private helper.
"""

from __future__ import annotations

import io
import zipfile

__all__ = ["ZIP_EPOCH", "deterministic_zip", "zip_parts"]

# The zip format cannot store a date before 1980; one fixed stamp makes the bytes reproducible.
ZIP_EPOCH = (1980, 1, 1, 0, 0, 0)


def deterministic_zip(parts: dict[str, bytes]) -> bytes:
    """A zip whose bytes depend on the content alone — not on when or where it was written."""
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as archive:
        for name in sorted(parts):
            info = zipfile.ZipInfo(name, date_time=ZIP_EPOCH)
            info.compress_type = zipfile.ZIP_DEFLATED
            info.create_system = 0  # same bytes on Linux and Windows, not just on one
            info.external_attr = 0o600 << 16
            archive.writestr(info, parts[name])
    return buffer.getvalue()


def zip_parts(data: bytes) -> dict[str, bytes]:
    """Every entry of a zip, by name — the read half, for pruning and for tests."""
    with zipfile.ZipFile(io.BytesIO(data)) as archive:
        return {name: archive.read(name) for name in archive.namelist()}
