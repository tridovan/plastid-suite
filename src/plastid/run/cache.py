"""Content-hash cache: re-running the same module on the same inputs re-exports, never recomputes.

The key covers module id + version + **source** + params + input hashes **in order** + what the
step was asked to produce + the IR fingerprint — change any of them and the step runs again.
Text is stored as text and bytes as base64, so a cache hit restores byte-identical outputs of
either kind.
"""

from __future__ import annotations

import base64
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

from plastid.result.export import Payload

__all__ = ["Cache"]


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
    ) -> str:
        """The identity of one step: which module, which code, which inputs in which order.

        Two things here are deliberate and were both measured as bugs first:

        * ``inputs`` keeps its **order**. It used to be sorted, which made the corpus's own
          "one sequence, two annotations" pair share a key in either direction — while six
          modules take ``genomes[0]`` as their reference, so the second run was served the
          first's tables (H-15, measured 05/10/2026).
        * ``source`` and ``payload`` are in the key. The first catches a behaviour change that
          never got a version bump (H-08); the second catches a run asking for an output the
          cached payload does not have (H-14, the `--docx` case).
        """
        payload = json.dumps(
            {
                "module": module_id,
                "version": module_version,
                "source": module_source,
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

    def get(self, key: str) -> Payload | None:
        path = self._path(key)
        if not path.is_file():
            return None
        payload = json.loads(path.read_text(encoding="utf-8"))
        files = payload.get("files")
        if not isinstance(files, dict):
            return None
        restored: Payload = dict(files)
        binary = payload.get("binary", {})
        if isinstance(binary, dict):
            restored.update({name: base64.b64decode(blob) for name, blob in binary.items()})
        return restored

    def put(self, key: str, files: Payload) -> None:
        self.root.mkdir(parents=True, exist_ok=True)
        payload: dict[str, dict[str, str]] = {
            "files": {name: text for name, text in files.items() if isinstance(text, str)}
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
