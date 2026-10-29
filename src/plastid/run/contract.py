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
import re
from collections.abc import Mapping, Sequence
from datetime import date
from pathlib import Path
from typing import Any

__all__ = [
    "ARTIFACT_SCHEMA_VERSION",
    "EXTERNAL_RESULT_SCHEMA_VERSION",
    "ContractError",
    "describe_artifact",
    "describe_external_result",
    "describe_reference_set",
    "sha256_file",
    "validate_artifact",
    "validate_external_result",
]

ARTIFACT_SCHEMA_VERSION = 1
EXTERNAL_RESULT_SCHEMA_VERSION = 1
EXTERNAL_RESULT_STATUSES = ("ok", "no_result")

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


_HEX64 = re.compile(r"^[0-9a-f]{64}$")


def describe_external_result(
    *,
    tool_name: str,
    tool_version: str,
    argv: Sequence[str],
    input_path: str | Path,
    query_map: Mapping[str, str],
    source: Mapping[str, Any],
    raw_path: str | Path,
    status: str = "ok",
    count: int | None = None,
) -> dict[str, Any]:
    """A result produced outside this toolkit, described so a run can consume it.

    A result read from elsewhere is a claim about a computation this run did not perform, so the
    description carries everything needed to repeat or refuse it: the tool (name, version, argv),
    the sha256 of the exact input, the map from this run's queries to that tool's identifiers,
    where and when the result came from, and the raw bytes with their own sha256 and size.

    A result that does not exist yet is ``status="no_result"`` with ``count=None`` — never 0,
    which reads as "measured, and nothing was found". A measured zero is ``status="ok"`` with
    ``count=0``, and the two are different states everywhere this contract is used (M11.16).
    """
    input_file, raw_file = Path(input_path), Path(raw_path)
    if not input_file.is_file():
        raise ContractError(f"{input_file}: the external result's input does not exist")
    if not raw_file.is_file():
        raise ContractError(f"{raw_file}: the external result's raw output does not exist")
    payload = {
        "external_result_schema_version": EXTERNAL_RESULT_SCHEMA_VERSION,
        "tool": {
            "name": tool_name,
            "version": tool_version,
            "argv": [str(argument) for argument in argv],
        },
        "input": {
            "path": str(input_file),
            "sha256": sha256_file(input_file),
            "bytes": input_file.stat().st_size,
        },
        "query_map": dict(query_map),
        "source": dict(source),
        "raw": {
            "path": str(raw_file),
            "sha256": sha256_file(raw_file),
            "bytes": raw_file.stat().st_size,
        },
        "status": status,
        "count": count,
    }
    validate_external_result(payload)
    return payload


def validate_external_result(payload: Mapping[str, Any]) -> None:
    """Refuse an external-result description that cannot be trusted as it stands."""
    if not isinstance(payload, Mapping):
        raise ContractError("an external result must be described by a mapping")
    for key in (
        "external_result_schema_version",
        "tool",
        "input",
        "query_map",
        "source",
        "raw",
        "status",
        "count",
    ):
        if key not in payload:
            raise ContractError(f"external result is missing required key {key!r}")
    version = payload["external_result_schema_version"]
    if version != EXTERNAL_RESULT_SCHEMA_VERSION:
        raise ContractError(
            f"unsupported external_result_schema_version {version!r} "
            f"(this build understands {EXTERNAL_RESULT_SCHEMA_VERSION})"
        )
    tool = payload["tool"]
    if not isinstance(tool, Mapping):
        raise ContractError("external result['tool'] must be a mapping")
    for key in ("name", "version"):
        if not isinstance(tool.get(key), str) or not tool[key].strip():
            raise ContractError(f"external result['tool'][{key!r}] must be a non-empty string")
    argv = tool.get("argv")
    if not isinstance(argv, list) or not argv:
        raise ContractError("external result['tool']['argv'] must be a non-empty list")
    input_info = payload["input"]
    if not isinstance(input_info, Mapping):
        raise ContractError("external result['input'] must be a mapping")
    _check_hashed_bytes(input_info, "input", minimum=1)
    raw = payload["raw"]
    if not isinstance(raw, Mapping):
        raise ContractError("external result['raw'] must be a mapping")
    _check_hashed_bytes(raw, "raw", minimum=0)
    query_map = payload["query_map"]
    if not isinstance(query_map, Mapping) or not query_map:
        raise ContractError("external result['query_map'] must map at least one query")
    for key, value in query_map.items():
        if not isinstance(key, str) or not key.strip():
            raise ContractError("external result['query_map'] keys must be non-empty strings")
        if not isinstance(value, str) or not value.strip():
            raise ContractError(
                f"external result['query_map'][{key!r}] must name the tool's identifier"
            )
    source = payload["source"]
    if not isinstance(source, Mapping):
        raise ContractError("external result['source'] must be a mapping")
    if not any(
        isinstance(source.get(key), str) and source[key].strip()
        for key in ("name", "url", "command")
    ):
        raise ContractError("external result['source'] must name a source or a command")
    retrieved = source.get("retrieved_at")
    if not isinstance(retrieved, str):
        raise ContractError("external result['source']['retrieved_at'] must be a date string")
    try:
        date.fromisoformat(retrieved)
    except ValueError as error:
        raise ContractError(
            f"external result['source']['retrieved_at'] = {retrieved!r} is not an ISO date"
        ) from error
    status = payload["status"]
    if status not in EXTERNAL_RESULT_STATUSES:
        raise ContractError(f"external result status {status!r} is not one this build knows")
    count = payload["count"]
    if status == "ok":
        if isinstance(count, bool) or not isinstance(count, int) or count < 0:
            raise ContractError(
                "an external result with status 'ok' must carry its count as a non-negative int"
            )
    elif count is not None:
        raise ContractError(
            "an external result with status 'no_result' must leave count as None — a measured "
            "zero is status 'ok' with count 0"
        )


def _check_hashed_bytes(info: Mapping[str, Any], label: str, *, minimum: int) -> None:
    digest = info.get("sha256")
    if not isinstance(digest, str) or not _HEX64.match(digest):
        raise ContractError(f"external result[{label!r}]['sha256'] must be a sha256 hex digest")
    size = info.get("bytes")
    if isinstance(size, bool) or not isinstance(size, int) or size < minimum:
        raise ContractError(
            f"external result[{label!r}]['bytes'] must be an int >= {minimum} (got {size!r})"
        )
