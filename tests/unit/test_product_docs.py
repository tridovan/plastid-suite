"""Product documentation checks — the two guards lost with the earlier cleanup (M11.14).

These are **product** checks: a module a user can run must be documented in the platform table,
and the README's status line must name the version the package ships. Neither reads `.dev/`,
`PROJECT.md`, a session log or any other maintainer document — a clone without them must still
be able to tell whether its own documentation is complete.

Each check is exercised twice: once against the real repository, and once against fabricated text
that *must* fail it. A guard nobody has watched fail is not a guard.
"""

from __future__ import annotations

import re
from pathlib import Path

from plastid.modules.registry import all_modules

ROOT = Path(__file__).parents[2]
PLATFORMS = ROOT / "docs" / "platforms.md"
README = ROOT / "README.md"
PYPROJECT = ROOT / "pyproject.toml"

_STATUS = re.compile(r"^version = \"([^\"]+)\"", re.MULTILINE)
# A row names its module as `` `module_id` `` in the first column. Matching the backticked token
# anywhere in the file would let a prose mention stand in for a documented row.
_ROW = re.compile(r"^\| `([^`]+)` \|", re.MULTILINE)


def _undocumented_modules(documented: str, registered: list[str]) -> list[str]:
    """Registered module ids with no row in the platform table — the failure this guard reports."""
    rows = set(_ROW.findall(documented))
    return sorted(name for name in registered if name not in rows)


def _readme_names_version(readme: str, version: str) -> bool:
    return f"**Status: v{version}.**" in readme


def test_every_registered_module_has_a_row_in_the_platform_table() -> None:
    """M9.4's table is what a user reads to decide whether the toolkit fits their machine.

    A module added later cannot quietly go undocumented. Set equality, not "contains", and the
    registry is asserted non-empty so a moved registry cannot make this pass over no modules.
    """
    registered = sorted(all_modules())
    assert registered, "no modules registered — the registry moved and this check would pass empty"
    missing = _undocumented_modules(PLATFORMS.read_text(encoding="utf-8"), registered)
    assert not missing, f"registered but absent from docs/platforms.md: {missing}"


def test_the_readme_status_line_names_the_version_we_ship() -> None:
    """A status line that restates phase state goes stale; this one names the version.

    A release bump therefore touches README in the same commit, or this fails.
    """
    match = _STATUS.search(PYPROJECT.read_text(encoding="utf-8"))
    assert match, "pyproject.toml has no version line"
    version = match.group(1)
    assert _readme_names_version(README.read_text(encoding="utf-8"), version), (
        f"README's status line does not say v{version} — the release bump missed it"
    )


def test_the_module_check_reports_a_missing_row() -> None:
    """The guard must fail on a module with no row, and not be fooled by prose."""
    fabricated = "| `genome_stats` | Length, GC | plastid | 1 | yes | — | measured |\n"
    assert _undocumented_modules(fabricated, ["genome_stats", "new_module"]) == ["new_module"]
    prose_only = "| `genome_stats` | `new_module` is mentioned here | plastid | 1 | yes | — | x |\n"
    assert _undocumented_modules(prose_only, ["new_module"]) == ["new_module"]


def test_the_version_check_reports_a_mismatch() -> None:
    """A stale or absent status line must fail both ways."""
    assert _readme_names_version("**Status: v1.0.0.**", "1.0.0")
    assert not _readme_names_version("**Status: v0.9.0.**", "1.0.0")
    assert not _readme_names_version("Phase 10 is open.", "1.0.0")
