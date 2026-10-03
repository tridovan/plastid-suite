"""Run manifest — the index that fixes CGAS's untraceable output (PROJECT.md §7.7).

``manifest.json`` records what ran, with what parameters, on what inputs, and where every
file went. ``validate()`` is hand-written on purpose: no schema library for one dict.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

__all__ = [
    "MANIFEST_SCHEMA_VERSION",
    "ManifestError",
    "build_manifest",
    "validate",
    "write_manifest",
]

MANIFEST_SCHEMA_VERSION = 1

_REQUIRED: dict[str, type] = {
    "schema_version": int,
    "run_id": str,
    "command": list,
    "package_version": str,
    "started_at": str,
    "status": str,
    "params": dict,
    "inputs": list,
    "coordinate_system": dict,
    "modules": list,
    "files": list,
}


class ManifestError(Exception):
    """The manifest is missing a required field or carries a wrong type."""


def build_manifest(
    *,
    run_id: str,
    command: list[str],
    package_version: str,
    started_at: str,
    params: dict[str, Any],
    inputs: list[dict[str, Any]],
    coordinate_system: dict[str, Any],
) -> dict[str, Any]:
    """A fresh manifest; the runner fills ``modules``/``files``/``status`` as it goes."""
    return {
        "schema_version": MANIFEST_SCHEMA_VERSION,
        "run_id": run_id,
        "command": list(command),
        "package_version": package_version,
        "started_at": started_at,
        "finished_at": None,
        "status": "running",
        "params": params,
        "inputs": list(inputs),
        "ir": None,
        "coordinate_system": coordinate_system,
        "modules": [],
        "files": [],
    }


def validate(manifest: dict[str, Any]) -> None:
    for key, expected in _REQUIRED.items():
        if key not in manifest:
            raise ManifestError(f"manifest is missing required key {key!r}")
        if not isinstance(manifest[key], expected):
            actual = type(manifest[key]).__name__
            raise ManifestError(f"manifest[{key!r}] should be {expected.__name__}, got {actual}")
    if manifest["schema_version"] != MANIFEST_SCHEMA_VERSION:
        raise ManifestError(
            f"unsupported schema_version {manifest['schema_version']} "
            f"(this build understands {MANIFEST_SCHEMA_VERSION})"
        )


def write_manifest(path: str | Path, manifest: dict[str, Any]) -> None:
    validate(manifest)
    Path(path).write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
