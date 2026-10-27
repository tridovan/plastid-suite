"""Run manifest — the index that fixes CGAS's untraceable output.

``manifest.json`` records what ran, with what parameters, on what inputs, and where every
file went. ``validate()`` is hand-written on purpose: no schema library for one dict.

**Version 2** (08/10/2026, P3B-04) adds a per-module ``reference`` entry: the privileged record a
step worked against (its path, accession and sha256 — the naming reference, D5-09). Version 1
manifests do not validate against this build, which is the point of a version number: a reader can
tell which fields a file may be missing.

**Version 3** (27/10/2026, M11.16) adds the two typed statements a downstream step needs before it
may read anything: ``artifacts`` (format, schema, sample, sha256, topology, frame, status and the
reference-set digest each output was built against — :mod:`plastid.run.contract`) and
``reference_sets`` (every search directory, its record count and digest). A version-2 manifest has
neither, so it still refuses to load rather than being read as "no artifacts".
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

MANIFEST_SCHEMA_VERSION = 3

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
    "artifacts": list,
    "reference_sets": list,
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
    dataset: str = "",
    reference_sets: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    """A fresh manifest; the runner fills ``modules``/``files``/``status`` as it goes.

    ``dataset`` is the run's label — the same string the run directory is slugged from. It is
    recorded here as well as in the path because a *batch* writes one run per sample into one
    output root, and a run that failed before it could write anything else has to say which sample
    it was: measured 09/10/2026 (P3B-09), where a failed sample's manifest named the module and
    the reason and nothing that identified the sample itself. Optional, so the manifests written
    before this key existed still validate (schema 2, unchanged — the key is additive).
    """
    return {
        "schema_version": MANIFEST_SCHEMA_VERSION,
        "run_id": run_id,
        "command": list(command),
        "package_version": package_version,
        "started_at": started_at,
        "finished_at": None,
        "status": "running",
        "dataset": dataset,
        "params": params,
        "inputs": list(inputs),
        "reference_sets": list(reference_sets or []),
        "ir": None,
        "coordinate_system": coordinate_system,
        "modules": [],
        "files": [],
        "artifacts": [],
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
