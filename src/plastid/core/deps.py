"""Optional-dependency guard for the extras (D2-04, D5-01, D5-02).

The minimum install is biopython + PyYAML (CLAUDE.md), so an export that needs more — a Word
document, a figure, a multi-sheet workbook — imports its library *inside the function* and, when
it is missing, raises one named error the CLI turns into exit 3. Never a traceback: the user is
running a biologist's install, and the message has to say which extra to add.
"""

from __future__ import annotations

from importlib import import_module
from types import ModuleType

__all__ = ["ExportDependencyError", "require_module"]


class ExportDependencyError(RuntimeError):
    """The optional dependency this export needs is not installed."""


def require_module(module_name: str, *, extra: str, purpose: str) -> ModuleType:
    """Import an optional dependency, or raise the named error that says how to get it."""
    try:
        return import_module(module_name)
    except ImportError as error:
        raise ExportDependencyError(
            f"{purpose} needs {module_name}, which is an optional dependency — "
            f"install it with: pip install 'plastid-suite[{extra}]'"
        ) from error
