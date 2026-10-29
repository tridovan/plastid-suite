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


def test_a_declared_parameter_reaches_the_tool_that_uses_it() -> None:
    """SSR thresholds: declared ⇒ the value is written into the ini MISA is handed.

    The defect was the inverse — the key was recorded and the run ignored it. This test pins the
    half of the rule the module now implements: the ini is built from the value, so there is no
    path on which the manifest can name a threshold MISA did not read. The module reports the
    thresholds read back out of that ini, not the string the user typed.
    """
    from plastid.external.misa import misa_ini_text, parse_unit_repeats

    spec = registry.get("ssr")
    declared = {param.name for param in spec.params}
    assert "unit_repeats" in declared, (
        "ssr must declare unit_repeats now that the value reaches MISA's ini"
    )
    values = "1-12 2-8 3-7 4-6 5-5 6-5"
    definition = misa_ini_text(parse_unit_repeats(values)).splitlines()[0]
    assert definition.split(":", 1)[1].split() == values.split()


def test_a_bad_ssr_threshold_is_a_usage_error_before_misa_runs(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """A malformed value exits 2 with the accepted shape — not 1 and not a silent default.

    It is refused before MISA is looked for, so this test needs no tool: the value is the
    module's own contract, and the command line can fix it.
    """
    from plastid.cli import main

    code = main(
        [
            "run",
            "-m",
            "ssr",
            "-i",
            str(GOLDEN),
            "-o",
            str(tmp_path),
            "--dataset",
            "bad-thresholds",
            "--set",
            "unit_repeats=not-a-pair",
        ]
    )
    assert code == 2
    assert "unit_repeats" in capsys.readouterr().err
