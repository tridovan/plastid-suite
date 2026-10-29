"""The module contract document and the registry must describe the same modules.

`docs/module_contracts.md` is the written definition a reader holds a result against; a module
added, removed or re-designated as IR-dependent without updating it would make the document a
claim rather than a definition. This test keeps the two in step, and keeps every block complete:
unit, input, tables, control, copy/pseudo/partial policy, naming/reference and missing rule.
"""

from __future__ import annotations

import ast
import re
from pathlib import Path

from plastid.modules import registry

ROOT = Path(__file__).resolve().parents[2]
DOC = ROOT / "docs" / "module_contracts.md"
MODULES = ROOT / "src" / "plastid" / "modules"

HEADING = re.compile(r"^### `([a-z_0-9]+)` — .*\(arity (\d+), IR gate: (yes|no)", re.MULTILINE)
REQUIRED_FIELDS = (
    "**Unit.**",
    "**Input.**",
    "**Tables.**",
    "**Control.**",
    "**Copies / pseudo / partial.**",
    "**Naming / comparison reference.**",
    "**Missing rule.**",
)
FROZEN = {"mito_annotate", "mito_genes"}


def _blocks() -> dict[str, tuple[int, bool, str]]:
    text = DOC.read_text(encoding="utf-8")
    matches = list(HEADING.finditer(text))
    blocks: dict[str, tuple[int, bool, str]] = {}
    for index, match in enumerate(matches):
        end = matches[index + 1].start() if index + 1 < len(matches) else len(text)
        blocks[match.group(1)] = (
            int(match.group(2)),
            match.group(3) == "yes",
            text[match.end() : end],
        )
    return blocks


def _table_ids(module_id: str) -> set[str]:
    path = MODULES / f"{module_id}.py"
    tree = ast.parse(path.read_text(encoding="utf-8"))
    ids = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Call) and getattr(node.func, "id", "") == "Table":
            for keyword in node.keywords:
                if keyword.arg == "id" and isinstance(keyword.value, ast.Constant):
                    ids.add(keyword.value.value)
    return ids


def test_every_registered_module_has_one_complete_contract_block() -> None:
    blocks = _blocks()
    registered = {spec.id for spec in map(registry.get, registry.all_modules())} - FROZEN
    assert set(blocks) == registered, (
        f"modules without a contract: {sorted(registered - set(blocks))}; "
        f"contracts without a module: {sorted(set(blocks) - registered)}"
    )
    for module_id, (arity, ir_gate, body) in blocks.items():
        spec = registry.get(module_id)
        assert arity == spec.arity, (
            f"{module_id}: document says arity {arity}, registry {spec.arity}"
        )
        assert ir_gate == spec.requires_ir, (
            f"{module_id}: document says IR gate {ir_gate}, registry {spec.requires_ir}"
        )
        for field in REQUIRED_FIELDS:
            assert field in body, f"{module_id}: the block has no {field}"
        tables_line = body.split("**Tables.**", 1)[1].split("\n- ", 1)[0]
        named = set(re.findall(r"`([a-z_0-9]+)`", tables_line))
        assert named, f"{module_id}: the Tables field names no table"
        missing = named - _table_ids(module_id)
        assert not missing, (
            f"{module_id}: the document names tables the module does not build: {missing}"
        )


def test_the_pipeline_and_the_contracts_point_at_each_other() -> None:
    pipeline = (ROOT / "docs" / "pipeline.md").read_text(encoding="utf-8")
    assert "module_contracts.md" in pipeline, "docs/pipeline.md does not link the module contracts"
    assert "docs/pipeline.md" in DOC.read_text(encoding="utf-8")
