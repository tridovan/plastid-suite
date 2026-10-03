"""Content-hash cache: re-running the same module on the same inputs re-exports, never recomputes.

The key covers module id + version + params + input hashes + IR fingerprint — change any of
them and the step runs again. Files are stored as text so a cache hit restores byte-identical
outputs.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

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
        ir_fingerprint: str = "",
    ) -> str:
        payload = json.dumps(
            {
                "module": module_id,
                "version": module_version,
                "params": params_hash,
                "inputs": sorted(input_sha256s),
                "ir": ir_fingerprint,
            },
            sort_keys=True,
        )
        return hashlib.sha256(payload.encode("utf-8")).hexdigest()

    def _path(self, key: str) -> Path:
        return self.root / f"{key}.json"

    def get(self, key: str) -> dict[str, str] | None:
        path = self._path(key)
        if not path.is_file():
            return None
        payload = json.loads(path.read_text(encoding="utf-8"))
        files = payload.get("files")
        return dict(files) if isinstance(files, dict) else None

    def put(self, key: str, files: dict[str, str]) -> None:
        self.root.mkdir(parents=True, exist_ok=True)
        self._path(key).write_text(
            json.dumps({"files": files}, ensure_ascii=False, indent=1) + "\n", encoding="utf-8"
        )
