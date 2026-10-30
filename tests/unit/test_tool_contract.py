"""The contract for running an external tool and for consuming a result made outside a run.

Two doors, one test file. A tool this toolkit runs must go through `run_tool()`, which records the
name, version, argv, exit code and the raw output bytes. A result produced somewhere else is a
claim about a computation this run did not perform, and its description carries the tool, the
input's sha256, the query map, the source and date, and the raw bytes — with "no result yet" a
different state from a measured count of zero.

The static half of the first test is what keeps the rule true as adapters are added: `subprocess`
may appear only in `run/tools.py`, or inside a function that probes a tool's version.
"""

from __future__ import annotations

import ast
import hashlib
from pathlib import Path

import pytest

from plastid.run.contract import (
    ContractError,
    describe_external_result,
    validate_external_result,
)
from plastid.run.tools import ToolError, run_tool

SRC = Path(__file__).resolve().parents[2] / "src"


class _SubprocessVisitor(ast.NodeVisitor):
    """Record every `subprocess` call with the function that contains it."""

    def __init__(self, relative: str) -> None:
        self.relative = relative
        self.uses: list[tuple[str, int, str]] = []
        self._functions: list[str] = []

    def visit_FunctionDef(self, node: ast.FunctionDef) -> None:
        self._functions.append(node.name)
        self.generic_visit(node)
        self._functions.pop()

    visit_AsyncFunctionDef = visit_FunctionDef  # type: ignore[assignment]

    def _record(self, lineno: int) -> None:
        self.uses.append((self.relative, lineno, self._functions[-1] if self._functions else ""))

    def visit_Call(self, node: ast.Call) -> None:
        function = node.func
        if (
            isinstance(function, ast.Attribute)
            and isinstance(function.value, ast.Name)
            and function.value.id == "subprocess"
        ):
            self._record(node.lineno)
        self.generic_visit(node)


def _subprocess_uses() -> list[tuple[str, int, str]]:
    """Every `subprocess` call outside the one door, with its enclosing function."""
    uses: list[tuple[str, int, str]] = []
    for path in sorted((SRC / "plastid").rglob("*.py")):
        relative = path.relative_to(SRC).as_posix()
        if relative == "plastid/run/tools.py":
            continue
        visitor = _SubprocessVisitor(relative)
        visitor.visit(ast.parse(path.read_text(encoding="utf-8"), filename=str(path)))
        uses.extend(visitor.uses)
    return uses


def test_external_tools_are_invoked_only_through_the_one_door() -> None:
    """`subprocess` lives in `run/tools.py`; adapters may use it only to probe a version."""
    violations = [
        f"{path}:{line} in {function}()"
        for path, line, function in _subprocess_uses()
        if not function.endswith("_version")
    ]
    assert not violations, (
        "subprocess outside run/tools.py must be a version probe; these calls are not: "
        + ", ".join(violations)
    )
    # The check is not vacuous: version probes are a real, named exception.
    assert any(function.endswith("_version") for _, _, function in _subprocess_uses())


def test_a_tool_call_records_name_version_argv_exit_code_and_raw_output(tmp_path: Path) -> None:
    """One invocation leaves everything Materials & Methods needs, including the raw bytes."""
    script = tmp_path / "stub_tool"
    script.write_text(
        "#!/bin/sh\n"
        'if [ "$1" = "--version" ]; then echo "1.2.3"; exit 0; fi\n'
        'echo "raw result line"\n'
        "exit 0\n",
        encoding="utf-8",
    )
    script.chmod(0o755)
    record = run_tool(
        "stub_tool", [str(script), "input.fa"], log_dir=tmp_path / "logs", cwd=tmp_path
    )
    assert record.name == "stub_tool"
    assert record.version == "1.2.3"
    assert record.version_source == "flag"
    assert record.argv == (str(script), "input.fa")
    assert record.rc == 0
    assert Path(record.stdout_file).read_text(encoding="utf-8") == "raw result line\n"
    assert Path(record.stderr_file).is_file()

    failing = tmp_path / "failing_tool"
    failing.write_text("#!/bin/sh\necho boom >&2\nexit 3\n", encoding="utf-8")
    failing.chmod(0o755)
    with pytest.raises(ToolError, match="exited 3"):
        run_tool("failing_tool", [str(failing)], log_dir=tmp_path / "logs2", cwd=tmp_path)
    # The raw bytes of the failed call are kept too: a refusal is evidence, not a deletion.
    assert any("failing_tool" in path.name for path in (tmp_path / "logs2").iterdir())


def _describe(tmp_path: Path, **overrides):
    input_file = tmp_path / "input.fa"
    input_file.write_text(">q1\nACGTACGT\n", encoding="utf-8")
    raw = tmp_path / "raw.tsv"
    raw.write_text("q1\tA>G\t2\n", encoding="utf-8")
    arguments = {
        "tool_name": "external_tool",
        "tool_version": "1.0",
        "argv": ["external_tool", "-i", str(input_file)],
        "input_path": input_file,
        "query_map": {"q1": "NC_008641.1:rbcL"},
        "source": {"command": "external_tool -i input.fa", "retrieved_at": "2026-10-29"},
        "raw_path": raw,
        "status": "ok",
        "count": 1,
    }
    arguments.update(overrides)
    return describe_external_result(**arguments)


def test_an_external_result_carries_the_tool_input_query_map_source_and_raw_bytes(
    tmp_path: Path,
) -> None:
    payload = _describe(tmp_path)
    validate_external_result(payload)
    assert payload["tool"] == {
        "name": "external_tool",
        "version": "1.0",
        "argv": ["external_tool", "-i", str(tmp_path / "input.fa")],
    }
    assert payload["input"]["sha256"] == hashlib.sha256(b">q1\nACGTACGT\n").hexdigest()
    assert payload["input"]["bytes"] == len(b">q1\nACGTACGT\n")
    assert payload["query_map"] == {"q1": "NC_008641.1:rbcL"}
    assert payload["source"]["retrieved_at"] == "2026-10-29"
    assert payload["raw"]["sha256"] == hashlib.sha256(b"q1\tA>G\t2\n").hexdigest()
    assert payload["raw"]["bytes"] == len(b"q1\tA>G\t2\n")
    assert payload["status"] == "ok" and payload["count"] == 1


def test_no_result_is_not_a_measured_zero(tmp_path: Path) -> None:
    """`no_result` leaves the count empty; a measured zero is `ok` and keeps its 0."""
    pending = _describe(tmp_path, status="no_result", count=None)
    validate_external_result(pending)
    assert pending["count"] is None

    measured_zero = _describe(tmp_path, status="ok", count=0)
    validate_external_result(measured_zero)
    assert measured_zero["count"] == 0

    with pytest.raises(ContractError, match="no_result"):
        validate_external_result(_describe(tmp_path, status="no_result", count=0))
    with pytest.raises(ContractError, match="non-negative int"):
        validate_external_result(_describe(tmp_path, status="ok", count=None))


def test_an_incomplete_description_is_refused(tmp_path: Path) -> None:
    payload = _describe(tmp_path)
    for key in ("tool", "input", "query_map", "source", "raw", "status", "count"):
        broken = dict(payload)
        del broken[key]
        with pytest.raises(ContractError, match=key):
            validate_external_result(broken)

    with pytest.raises(ContractError, match="query_map"):
        validate_external_result(_describe(tmp_path, query_map={}))
    with pytest.raises(ContractError, match="ISO date"):
        validate_external_result(
            _describe(tmp_path, source={"command": "external_tool", "retrieved_at": "29/10/2026"})
        )
    with pytest.raises(ContractError, match="source"):
        validate_external_result(_describe(tmp_path, source={"retrieved_at": "2026-10-29"}))
    with pytest.raises(ContractError, match="tool_version|version"):
        validate_external_result(_describe(tmp_path, tool_version=""))
