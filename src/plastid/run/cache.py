"""Content-hash cache: re-running the same module on the same inputs re-exports, never recomputes.

The key covers module id + version + **source** + shared package code/resources + params + input
hashes **in order** + what the step was asked to produce + the IR fingerprint — change any of them
and the step runs again.

A cache entry carries the **provenance of the computation that produced it** (the tool calls it
made, when, and the key components), because a hit restores files, not a live computation. The
runner revalidates the recorded tool versions before trusting a hit: a tool that moved, vanished
or changed version since the payload was computed makes that payload stale, and it is recomputed
rather than served under the current run's name (M11.12).

Text is stored as text and bytes as base64, so a cache hit restores byte-identical outputs of
either kind.
"""

from __future__ import annotations

import base64
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from plastid.result.export import Payload

__all__ = ["CACHE_SCHEMA_VERSION", "Cache", "CacheEntry"]

# Bumped when an entry's on-disk shape changes. Entries written by an older schema are not read:
# their provenance may be missing, and a payload whose provenance cannot be stated is not one this
# runner should serve (M11.12). The key would normally differ anyway; this is the explicit guard.
CACHE_SCHEMA_VERSION = 2


@dataclass(frozen=True, slots=True)
class CacheEntry:
    """A stored payload and the provenance of the computation that produced it."""

    files: Payload
    provenance: dict[str, Any]


@dataclass(frozen=True, slots=True)
class Cache:
    root: Path

    @staticmethod
    def key(
        *,
        module_id: str,
        module_version: str,
        params_hash: str,
        input_sha256s: tuple[str, ...],
        module_source: str = "",
        payload_kinds: str = "",
        ir_fingerprint: str = "",
        shared_code: str = "",
    ) -> str:
        """The identity of one step: which module, which code, which inputs in which order.

        Three things here are deliberate and were all measured as bugs first:

        * ``inputs`` keeps its **order**. It used to be sorted, which made the corpus's own
          "one sequence, two annotations" pair share a key in either direction — while six
          modules take ``genomes[0]`` as their reference, so the second run was served the
          first's tables (H-15, measured 05/10/2026).
        * ``source`` and ``payload`` are in the key. The first catches a behaviour change that
          never got a version bump (H-08); the second catches a run asking for an output the
          cached payload does not have (H-14, the `--docx` case).
        * ``shared_code`` covers what the module file alone cannot: ``core/``, ``result/``,
          ``run/`` and every packaged resource (``defaults.yaml``, templates, R scripts, the
          reference name lists) that can change what a module produces without its own file being
          touched. Hashing only the module source was a narrower claim than the key made (M11.12).
        """
        payload = json.dumps(
            {
                "module": module_id,
                "version": module_version,
                "source": module_source,
                "shared_code": shared_code,
                "params": params_hash,
                "inputs": list(input_sha256s),
                "payload": payload_kinds,
                "ir": ir_fingerprint,
            },
            sort_keys=True,
        )
        return hashlib.sha256(payload.encode("utf-8")).hexdigest()

    def _path(self, key: str) -> Path:
        return self.root / f"{key}.json"

    def get(self, key: str) -> CacheEntry | None:
        """The entry, or ``None`` when absent or written by an older schema."""
        path = self._path(key)
        if not path.is_file():
            return None
        payload = json.loads(path.read_text(encoding="utf-8"))
        if payload.get("schema") != CACHE_SCHEMA_VERSION:
            return None
        files = payload.get("files")
        provenance = payload.get("provenance")
        if not isinstance(files, dict) or not isinstance(provenance, dict):
            return None
        restored: Payload = dict(files)
        binary = payload.get("binary", {})
        if isinstance(binary, dict):
            restored.update({name: base64.b64decode(blob) for name, blob in binary.items()})
        return CacheEntry(files=restored, provenance=provenance)

    def put(self, key: str, files: Payload, *, provenance: dict[str, Any]) -> None:
        self.root.mkdir(parents=True, exist_ok=True)
        payload: dict[str, Any] = {
            "schema": CACHE_SCHEMA_VERSION,
            "provenance": provenance,
            "files": {name: text for name, text in files.items() if isinstance(text, str)},
        }
        binary = {
            name: base64.b64encode(blob).decode("ascii")
            for name, blob in files.items()
            if isinstance(blob, bytes)
        }
        if binary:
            # A .docx is a zip: JSON has no bytes, and a cache hit must restore it unchanged.
            payload["binary"] = binary
        self._path(key).write_text(
            json.dumps(payload, ensure_ascii=False, indent=1) + "\n", encoding="utf-8"
        )
