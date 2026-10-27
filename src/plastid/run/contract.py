"""The artifact and reference-set contract: what a run produced, and what it was produced from.

A run directory is only traceable if a reader can answer three questions without guessing:
**which file is the artifact**, **what it is** (format, schema, sample, topology, frame, status),
and **which reference set it was transferred from** (path, record count, digest). The manifest's
``files`` list already says where every output went and its sha256; this module adds the typed
statements around the files that a *next step* reads, so a downstream step can refuse an artifact
that is missing, half-written, or built against a different reference set.

Kept hand-written like :mod:`plastid.run.manifest` — one dict, no schema library. The format is
versioned on its own (``ARTIFACT_SCHEMA_VERSION``) so a reader can tell which fields a file may be
missing without tying that to the manifest's version.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any

__all__ = [
    "ARTIFACT_SCHEMA_VERSION",
    "ContractError",
    "describe_artifact",
    "describe_reference_set",
    "sha256_file",
    "validate_artifact",
]

ARTIFACT_SCHEMA_VERSION = 1

ARTIFACT_FORMATS = ("genbank",)
ARTIFACT_STATUSES = ("ok", "cached")


class ContractError(Exception):
    """An artifact or reference set cannot be described as the contract requires."""


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def describe_reference_set(
    directory: str | Path,
    *,
    role: str = "annotation_db",
    exclude: tuple[str, ...] = (),
    leave_out: str | None = None,
) -> dict[str, Any]:
    """One directory a run searched: records actually searched and a digest over their bytes.

    ``exclude`` are accessions the run was told to drop (``--exclude``) and ``leave_out`` is the
    target's own accession, which the annotator never searches (leave-one-out, D5-08). Both are
    recorded so the digest names the set that was used, not every file that happened to sit in
    the directory. A path alone cannot say whether the set changed; raises
    :class:`ContractError` when nothing is left to search rather than writing an empty contract a
    downstream step would trust.
    """
    root = Path(directory)
    present = sorted(root.glob("*.gb"))
    dropped = set(exclude)
    if leave_out:
        dropped.add(leave_out)
    records = [record for record in present if record.stem not in dropped]
    if not records:
        raise ContractError(f"{root}: no GenBank records to describe as {role!r}")
    listing = "".join(f"{record.name}\t{sha256_file(record)}\n" for record in records)
    return {
        "role": role,
        "format": "genbank",
        "path": str(root),
        "records": len(records),
        "excluded": sorted(record.stem for record in present if record.stem in dropped),
        "leave_out": leave_out,
        "sha256": hashlib.sha256(listing.encode("utf-8")).hexdigest(),
    }


def describe_artifact(
    path: str | Path,
    *,
    relative_to: str | Path,
    sample_id: str,
    topology: str,
    status: str,
    frame: str = "genbank",
    format: str = "genbank",
    reference_sha256: str | None = None,
) -> dict[str, Any]:
    """One file a next step may read, with every field that decides whether it may be read."""
    file = Path(path)
    if not file.is_file():
        raise ContractError(f"{file}: the artifact does not exist")
    if format not in ARTIFACT_FORMATS:
        raise ContractError(f"{file}: {format!r} is not an artifact format this build knows")
    if status not in ARTIFACT_STATUSES:
        raise ContractError(f"{file}: status {status!r} is not one a step may consume")
    entry: dict[str, Any] = {
        "path": str(file.relative_to(relative_to)),
        "format": format,
        "artifact_schema_version": ARTIFACT_SCHEMA_VERSION,
        "sample_id": sample_id,
        "sha256": sha256_file(file),
        "bytes": file.stat().st_size,
        "topology": topology,
        "frame": frame,
        "status": status,
    }
    if reference_sha256 is not None:
        entry["reference_sha256"] = reference_sha256
    validate_artifact(entry)
    return entry


def validate_artifact(entry: dict[str, Any]) -> None:
    required = {
        "path": str,
        "format": str,
        "artifact_schema_version": int,
        "sample_id": str,
        "sha256": str,
        "bytes": int,
        "topology": str,
        "frame": str,
        "status": str,
    }
    for key, expected in required.items():
        if key not in entry:
            raise ContractError(f"artifact is missing required key {key!r}")
        if not isinstance(entry[key], expected):
            raise ContractError(
                f"artifact[{key!r}] should be {expected.__name__}, got {type(entry[key]).__name__}"
            )
    if entry["artifact_schema_version"] != ARTIFACT_SCHEMA_VERSION:
        raise ContractError(
            f"unsupported artifact_schema_version {entry['artifact_schema_version']} "
            f"(this build understands {ARTIFACT_SCHEMA_VERSION})"
        )
    if entry["format"] not in ARTIFACT_FORMATS:
        raise ContractError(f"unsupported artifact format {entry['format']!r}")
    if entry["status"] not in ARTIFACT_STATUSES:
        raise ContractError(f"artifact status {entry['status']!r} may not be consumed")
