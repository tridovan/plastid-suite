"""Licensing guards — what may ship, and what must be declared (M9.3, D-p, D-q).

The defect this file exists to prevent was **measured, not imagined**: `uv build --wheel` shipped
`plastid/vendor/misa/misa.pl`, a third-party script whose source zip states no licence
(docs/licensing.md §4). A work with no licence is all rights reserved, so bundling it was the one
act we may not do — MISA left the package on 11/10/2026. These checks keep the rule from decaying:
**a file we redistribute must state its licence, and every tool the code calls must be named in
the inventory**.
"""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src" / "plastid"
INVENTORY = ROOT / "docs" / "licensing.md"

# Spellings that mean "no licence granted". `Unknown` is what bioconda records for Vmatch, whose
# upstream licence is academic-only — either way, not a licence to redistribute.
_NOT_STATED = ("not stated", "unknown", "unclear", "no licence", "none")


def test_bundled_third_party_code_states_a_licence() -> None:
    """A file the package redistributes must name its licence in its ``SOURCE.json``."""
    for source in sorted((SRC / "vendor").glob("*/SOURCE.json")):
        meta = json.loads(source.read_text(encoding="utf-8"))
        bundled = sorted(name for name in meta.get("files", {}) if (source.parent / name).is_file())
        if not bundled:
            continue  # declared but not shipped — the MISA case after D-q
        licence = str(meta.get("license", "")).strip()
        assert licence, f"{source}: bundles {bundled} and states no licence at all"
        assert not any(word in licence.lower() for word in _NOT_STATED), (
            f"{source}: bundles {bundled} but its licence reads {licence!r} — a work with no "
            f"licence is all rights reserved, so it may not be redistributed (docs/licensing.md)"
        )


def test_every_external_adapter_is_named_in_the_licence_inventory() -> None:
    """Every tool the code calls has a row: its licence, and whether we redistribute it."""
    adapters = sorted(
        path.stem for path in (SRC / "external").glob("*.py") if path.stem != "__init__"
    )
    assert adapters, "no adapters found — the path moved, and this check would pass vacuously"
    text = INVENTORY.read_text(encoding="utf-8")
    missing = [name for name in adapters if f"external/{name}.py" not in text]
    assert not missing, (
        f"tools the code calls but the inventory never names: {missing} — add a row to "
        f"{INVENTORY.name} with the licence and whether the toolkit ships it"
    )


def test_the_toolkit_ships_and_declares_its_own_licence() -> None:
    """M9.3's other half: our own licence, in the file and in the package metadata."""
    licence = (ROOT / "LICENSE").read_text(encoding="utf-8")
    assert licence.startswith("MIT License"), "LICENSE is not the MIT text D-p decided"
    pyproject = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    assert 'license = "MIT"' in pyproject, "pyproject does not declare the licence (seam 5)"
