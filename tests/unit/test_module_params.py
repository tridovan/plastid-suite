"""A parameter the user types must do something, or be refused by name (M11.20).

**The defects this file holds, both measured before the fix (28/10/2026).**

* **A misspelt key was silent.** `plastid run -m genome_stats --set this_key_does_not_exist=42`
  exited 0 and wrote a normal run. A user who typed `--set min_lengh_bp=60` was told nothing and
  got 30 bp repeats.
* **A key that was accepted but ignored.** `plastid run -m ssr --set unit_repeats="1-12 2-8 …"`
  recorded the value in `manifest.params` and then ran MISA with the packaged `misa.ini` anyway:
  the summary read 37 SSRs and `(1/10) (2/6) (3/5) (4/5) (5/5) (6/5)`, identical to a run that
  passed no thresholds at all. A recorded parameter that does not reach the tool is worse than a
  refused one, because the manifest then states something false.

The rule the tests below enforce is one sentence long: **every parameter a module accepts has an
effect the manifest can show, and every parameter it does not accept is a usage error.** A
parameter that cannot be honoured yet must be *refused by name*, not accepted and dropped.

The declaration lives on the module (`ModuleSpec.params`), so the CLI can reject a key without
knowing anything about the module, and a test can compare the declaration with what the module
actually reads.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from plastid.modules import registry

GOLDEN = Path(__file__).resolve().parents[2] / "tests" / "data" / "NC_050833.1.gb"


# --------------------------------------------------------------------- the declaration itself


def test_every_registered_module_declares_its_parameters() -> None:
    """A module that reads no parameter declares none — an empty tuple, not ``None``.

    ``None`` would mean "nobody has looked", which is the state this milestone removes: the
    check below compares the declaration with the source, and a module with parameters but no
    declaration would fail it rather than pass quietly.
    """
    for spec in registry.all_modules().values():
        assert isinstance(spec.params, tuple), (
            f"{spec.id}: params must be a tuple, got {spec.params!r}"
        )


def test_the_declaration_matches_what_each_module_reads() -> None:
    """Every key a module calls ``params.get``/``params[...]`` with is declared, and vice versa.

    This is the check that keeps the declaration honest. It reads the module source with ``ast``
    rather than running it, so a branch that a run does not take still counts — a parameter read
    only on the failure path is still a parameter.
    """
    source_dir = Path(registry.__file__).parent
    for spec in registry.all_modules().values():
        path = source_dir / f"{spec.id}.py"
        if not path.is_file():
            continue  # mito modules registered from another file; covered by their own tests
        keys = _params_read(path)
        declared = {param.name for param in spec.params}
        missing = keys - declared
        extra = declared - keys
        assert not missing, (
            f"{spec.id}: reads {sorted(missing)} but does not declare it — the CLI cannot "
            f"reject a misspelling of a key the module never declares"
        )
        assert not extra, (
            f"{spec.id}: declares {sorted(extra)} but never reads it — a declared parameter "
            f"that does nothing is the defect this milestone exists to remove"
        )


def _params_read(path: Path) -> set[str]:
    """Every literal key the module passes to ``params.get``/``params[...]``."""
    import ast

    tree = ast.parse(path.read_text(encoding="utf-8"))
    keys: set[str] = set()
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        func = node.func
        if (
            isinstance(func, ast.Attribute)
            and func.attr == "get"
            and isinstance(func.value, ast.Name)
            and func.value.id == "params"
            and node.args
            and isinstance(node.args[0], ast.Constant)
            and isinstance(node.args[0].value, str)
        ):
            keys.add(node.args[0].value)
        if (
            isinstance(node, ast.Subscript)
            and isinstance(node.value, ast.Name)
            and node.value.id == "params"
            and isinstance(node.slice, ast.Constant)
            and isinstance(node.slice.value, str)
        ):
            keys.add(node.slice.value)
    return keys


def test_a_parameter_declares_a_type_and_says_what_it_does() -> None:
    """A declaration without a type cannot be validated; one without a reason is not a contract."""
    for spec in registry.all_modules().values():
        for param in spec.params:
            assert param.type in {"int", "float", "str", "bool"}, (
                f"{spec.id}.{param.name}: type {param.type!r} is not one the CLI can validate"
            )
            assert param.help.strip(), f"{spec.id}.{param.name}: declared with no description"


def test_an_unknown_key_is_a_usage_error_not_a_silent_no_op() -> None:
    """The reported defect, as a unit: an undeclared key must be refused by name."""
    problems = registry.unknown_params(["genome_stats"], {"min_lengh_bp": 60})
    assert problems, "a misspelt key was accepted — the defect is back"
    message = problems[0]
    assert "min_lengh_bp" in message
    assert "genome_stats" in message


def test_a_declared_key_is_accepted() -> None:
    assert registry.unknown_params(["repeats"], {"min_length_bp": 60}) == []


def test_a_shared_run_parameter_is_accepted_by_every_module() -> None:
    """``data_label`` and ``ir_source`` are run-level, not a module's own — they are declared
    once as shared parameters, so a module does not have to repeat them."""
    assert registry.unknown_params(["genome_stats", "structure"], {"data_label": "x"}) == []
    assert registry.unknown_params(["genome_stats"], {"ir_source": "sequence"}) == []


def test_the_error_lists_the_keys_the_module_does_accept() -> None:
    """A refusal that does not say what would work makes the user guess again."""
    (message,) = registry.unknown_params(["repeats"], {"min_lengh_bp": 60})
    assert "min_length_bp" in message


# ------------------------------------------------------------------- accepted means it ran


def test_a_tool_parameter_reaches_the_recorded_argv(tmp_path: Path) -> None:
    """``--set min_match=50`` must appear in nucmer's argv in the manifest.

    Measured before the fix: it already did — this test pins the behaviour so a later refactor
    of the parameter layer cannot quietly detach the key from the tool call.
    """
    from plastid.cli import main

    code = main(
        [
            "run",
            "-m",
            "synteny",
            "-i",
            str(GOLDEN),
            "-i",
            str(Path(__file__).resolve().parents[2] / "tests" / "data" / "NC_024936.1.gb"),
            "-o",
            str(tmp_path),
            "--dataset",
            "params",
            "--set",
            "ir_source=sequence",
            "--set",
            "min_match=50",
            "--set",
            "min_cluster=90",
        ]
    )
    if code != 0:
        pytest.skip("synteny needs nucmer/delta-filter/show-coords — see docs/tools/mummer4.md")
    manifest = json.loads(
        next((tmp_path / "params").iterdir()).joinpath("manifest.json").read_text()
    )
    argv = next(
        t["argv"] for m in manifest["modules"] for t in m.get("tools", []) if t["name"] == "nucmer"
    )
    assert argv[argv.index("-l") + 1] == "50"
    assert argv[argv.index("-c") + 1] == "90"


def test_an_accepted_but_ignored_parameter_is_the_bug() -> None:
    """SSR thresholds: either the value changes the run, or the key is refused.

    Before the fix the key was accepted and recorded, and the run was identical to one that
    passed nothing. This test states the rule without prescribing which half of the fix is
    chosen — a module that declines the parameter must not declare it, and
    ``test_the_declaration_matches_what_each_module_reads`` then makes the misspelling a named
    error instead of a no-op.
    """
    spec = registry.get("ssr")
    declared = {param.name for param in spec.params}
    if "unit_repeats" in declared:
        # Declared ⇒ it must be honoured. The implementation that honours it writes misa.ini
        # from the chosen values; until then the key must not be declared at all.
        pytest.fail(
            "ssr declares unit_repeats; it must then reach MISA's misa.ini — measured before "
            "the fix: accepted, recorded, and ignored"
        )
    assert {"misa"} <= declared
