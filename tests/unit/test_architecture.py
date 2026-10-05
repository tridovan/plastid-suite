"""Architecture guard — enforces the import layering (CLAUDE.md → Code rules).

Each layer may import itself plus everything to its left:

    core/     -> stdlib + {Bio, yaml} only
    result/   -> core
    run/      -> core, result
    report/   -> core, result
    external/ -> core, result, run
    modules/  -> core, result, run, external, and modules.registry ONLY
                 (a module never imports another module)
    cli.py, __main__.py, __init__.py -> exempt (process entry points)

A violation fails the build, so the dependency graph stays a DAG while ~20
modules are added over the project's lifetime.
"""

from __future__ import annotations

import ast
import sys
from pathlib import Path

SRC = Path(__file__).resolve().parents[2] / "src" / "plastid"

# Third-party top-level imports allowed inside the package (besides the stdlib).
#
# ``docx`` and ``lxml`` are the optional Word-export extra (D2-04), ``matplotlib`` the figures
# extra, ``scipy`` the statistics extra (D5-02) and ``openpyxl`` the Excel extra (P5-04). They
# are imported inside functions,
# never at module import time, so the minimum install still runs without them — the allowlist is
# what says the dependency is deliberate rather than accidental.
THIRD_PARTY_ALLOWED = frozenset({"Bio", "yaml", "docx", "lxml", "matplotlib", "openpyxl", "scipy"})

# layer -> allowed `plastid.*` import prefixes
LAYERS: dict[str, frozenset[str]] = {
    "core": frozenset({"plastid.core"}),
    "result": frozenset({"plastid.core", "plastid.result"}),
    "run": frozenset({"plastid.core", "plastid.result", "plastid.run"}),
    "report": frozenset({"plastid.core", "plastid.result", "plastid.report"}),
    "external": frozenset({"plastid.core", "plastid.result", "plastid.run", "plastid.external"}),
    "modules": frozenset(
        {
            "plastid.core",
            "plastid.result",
            "plastid.run",
            "plastid.external",
            "plastid.modules.registry",
        }
    ),
}


def _layer_of(path: Path) -> str | None:
    """Return the layer name, or None for exempt top-level files."""
    rel = path.relative_to(SRC)
    return rel.parts[0] if len(rel.parts) > 1 else None


def _imported_names(tree: ast.AST) -> list[str]:
    names: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module is not None:
            names.append(node.module)
            # `from plastid import core` also yields the resolved `plastid.core`.
            names.extend(f"{node.module}.{alias.name}" for alias in node.names)
        elif (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id == "require_module"
            and node.args
            and isinstance(node.args[0], ast.Constant)
            and isinstance(node.args[0].value, str)
        ):
            # An optional dependency named by string is still a dependency: `require_module` is
            # how matplotlib, scipy and openpyxl enter, and without this the allowlist below
            # would say nothing about them — a guard that cannot see the thing it guards.
            names.append(node.args[0].value)
    return names


def _allowed_plastid(layer: str, path: Path) -> frozenset[str]:
    # modules/__init__.py registers the sibling modules — the only file allowed to.
    if layer == "modules" and path.name == "__init__.py":
        return LAYERS["modules"] | {"plastid.modules"}
    return LAYERS[layer]


def test_expected_layers_exist() -> None:
    for layer in LAYERS:
        assert (SRC / layer / "__init__.py").is_file(), f"missing layer package: {layer}"


def test_import_layering() -> None:
    violations: list[str] = []
    for path in sorted(SRC.rglob("*.py")):
        layer = _layer_of(path)
        if layer is None:
            continue
        allowed = _allowed_plastid(layer, path)
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        rel = path.relative_to(SRC)
        for name in _imported_names(tree):
            top = name.split(".")[0]
            if top == "plastid":
                if not any(name == prefix or name.startswith(prefix + ".") for prefix in allowed):
                    violations.append(
                        f"{rel}: imports {name} — layer '{layer}' allows {sorted(allowed)}"
                    )
            elif top not in sys.stdlib_module_names and top not in THIRD_PARTY_ALLOWED:
                violations.append(
                    f"{rel}: imports third-party {top!r} — allowed: {sorted(THIRD_PARTY_ALLOWED)}"
                )
    assert not violations, "import layering violated:\n" + "\n".join(violations)
