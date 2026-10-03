"""Central defaults (``src/plastid/data/defaults.yaml``) — thresholds live in one place."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

from plastid.core.resources import data_file

__all__ = ["default_defaults_path", "load_defaults"]


def default_defaults_path() -> Path:
    return data_file("defaults.yaml")


def load_defaults(path: str | Path | None = None) -> dict[str, Any]:
    """Load the defaults file; an explicit path overrides the packaged one (tests, users)."""
    source = Path(path) if path is not None else default_defaults_path()
    payload = yaml.safe_load(source.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"{source}: defaults must be a mapping")
    return payload
