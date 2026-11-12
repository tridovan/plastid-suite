"""The command line as a surface: what a word means on every command, and what a person can find.

The interface is generated from this surface, and it shows the equivalent command beside every run,
so a word that means two things here would mean two things on a form. These checks read the parser
and the registry themselves, so a flag added later is held to the same rules without anyone
remembering them. Each rule also has a check on made-up input: a guard nobody has watched fail is
not a guard.
"""

from __future__ import annotations

import argparse
import json
import re
import shlex
from pathlib import Path

import pytest

from plastid.cli import build_parser, main
from plastid.modules import registry
from plastid.result.figures import all_figures, all_r_figures
from plastid.run.state import EXIT_STATES

ROOT = Path(__file__).parents[2]
DATA = ROOT / "tests" / "data"
GENOME = str(DATA / "NC_050833.1.gb")

# Mito is frozen: its two modules and its one command keep their shape, flags included.
FROZEN = {"mito annotate"}
MAINTAINER_CODE = re.compile(
    r"\b(D\d+-\d+|D-[a-z]\b|P\d+[A-Z]?-\d+|M\d+\.\d+|H-\d+|K-\d+|Phase \d+[A-Z]?)|manuscript"
)
# A command word in prose: the toolkit's name followed by one of its commands.
_COMMAND_NAMES = (
    "run|batch|annotate|assemble|panel|reference|editing-reference|figures|package|qc|ir-junction"
    "|modules|tools|runs|mito"
)
_INLINE = re.compile(rf"`(plastid (?:{_COMMAND_NAMES})\b[^`\n]*)`")
_CODE_LINE = re.compile(rf"^\s*(plastid (?:{_COMMAND_NAMES})\b.*)$", re.MULTILINE)
DOCUMENTS = (
    ROOT / "README.md",
    *sorted((ROOT / "docs").glob("*.md")),
    *sorted((ROOT / "docs" / "tools").glob("*.md")),
)
# Versioned design documents, kept as they were written: they describe an older interface.
HISTORICAL = {"API_v0.1.md", "API_v0.2.md"}


def _subparsers(parser: argparse.ArgumentParser) -> dict[str, argparse.ArgumentParser]:
    found = [a for a in parser._actions if isinstance(a, argparse._SubParsersAction)]
    return dict(found[0].choices) if found else {}


def _commands() -> dict[str, argparse.ArgumentParser]:
    """Every command a person can type, by its full name: ``panel pin``, ``run``."""
    leaves: dict[str, argparse.ArgumentParser] = {}

    def walk(parser: argparse.ArgumentParser, path: list[str]) -> None:
        children = _subparsers(parser)
        if not children:
            leaves[" ".join(path)] = parser
        for name, child in children.items():
            walk(child, [*path, name])

    walk(build_parser(), [])
    return leaves


def _options(parser: argparse.ArgumentParser) -> list[argparse.Action]:
    return [action for action in parser._actions if not isinstance(action, argparse._HelpAction)]


def _flags(parser: argparse.ArgumentParser) -> set[str]:
    return {flag for action in _options(parser) for flag in action.option_strings}


def _reader_text() -> list[tuple[str, str]]:
    """Everything the command line and the registry put in front of a person."""
    texts: list[tuple[str, str]] = [("plastid --help", build_parser().format_help())]
    for name, parser in _commands().items():
        texts.append((f"plastid {name} --help", parser.format_help()))
    for spec in registry.all_modules().values():
        texts.append((f"module {spec.id}", spec.title))
        for param in spec.params:
            texts.append((f"module {spec.id} {param.name}", param.help))
    for param in registry.SHARED_PARAMS:
        texts.append((f"shared {param.name}", param.help))
    for figure in (*all_figures(), *all_r_figures()):
        texts.append((f"figure {figure.id}", figure.title))
    return texts


def _maintainer_codes(texts: list[tuple[str, str]]) -> list[str]:
    return [
        f"{where}: {match.group(0)}"
        for where, text in texts
        for match in MAINTAINER_CODE.finditer(text)
    ]


def test_every_argument_says_what_it_is_for() -> None:
    silent = [
        f"{name} {action.option_strings or action.dest}"
        for name, parser in _commands().items()
        for action in _options(parser)
        if not action.help
    ]
    assert not silent


def test_nothing_a_user_reads_carries_a_planning_code() -> None:
    assert _maintainer_codes(_reader_text()) == []
    assert _maintainer_codes([("made up", "see D5-08, and P3B-03")]) == [
        "made up: D5-08",
        "made up: P3B-03",
    ]


def test_out_is_always_the_folder_that_holds_runs() -> None:
    wrong = {
        name: parser.get_default("out")
        for name, parser in _commands().items()
        if "--out" in _flags(parser) and parser.get_default("out") != ".runs"
    }
    assert not wrong


def test_a_retired_word_has_not_come_back() -> None:
    """Each of these meant several things, or read as something it was not."""
    retired = {
        "--reference": "--reference-set, --compare-to, --order-from or --reference-file",
        "--panel": "--reference-set",
        "--force": "--no-cache",
        "--backend": "`plastid assemble import` or --skip-resource-check",
        "--runs": "--out",
        "--retmax": "--max-records",
        "--audit-only": "`plastid reference audit`",
    }
    found = [
        f"{name}: {flag} (use {retired[flag]})"
        for name, parser in _commands().items()
        if name not in FROZEN
        for flag in sorted(_flags(parser) & set(retired))
    ]
    assert not found


def test_a_command_that_makes_a_run_takes_the_same_flags() -> None:
    makes_a_run = (
        "run",
        "batch",
        "annotate",
        "ir-junction",
        "assemble run",
        "assemble prepare",
        "assemble import",
    )
    commands = _commands()
    missing = [
        f"{name} {flag}"
        for name in makes_a_run
        for flag in ("--out", "--dataset", "--label")
        if flag not in _flags(commands[name])
    ]
    assert not missing


def test_a_command_group_names_all_its_actions(capsys: pytest.CaptureFixture[str]) -> None:
    parser = build_parser()
    for group, child in _subparsers(parser).items():
        actions = _subparsers(child)
        if not actions:
            continue
        assert main([group]) == 2
        said = capsys.readouterr().err
        assert all(action in said for action in actions), (group, said)


def _spell(text: str) -> list[str]:
    """A command as written in a document, with its placeholders made into plain words."""
    cleaned = re.sub(r"<[^>]*>|\$\([^)]*\)", "X", text.replace("\\\n", " ")).replace("…", "")
    cleaned = cleaned.split("#")[0]
    try:
        return shlex.split(cleaned)
    except ValueError:
        return cleaned.split()


def _unknown_words(text: str, commands: dict[str, argparse.ArgumentParser]) -> list[str]:
    """Why a command written in a document would not run: no such command, or no such flag."""
    words = _spell(text)[1:]
    name = words[0] if words else ""
    groups = {command.split()[0] for command in commands}
    if name not in groups:
        return [f"no command {name!r}"]
    if name not in commands:
        action = words[1] if len(words) > 1 else ""
        if action.startswith("-") or not action:
            return []
        name = f"{name} {action}"
        if name not in commands:
            return [f"no command {name!r}"]
    flags = _flags(commands[name])
    return [
        f"{name}: no flag {word.split('=')[0]}"
        for word in words
        if word.startswith("-")
        and not re.fullmatch(r"-?\d.*", word)
        and word.split("=")[0] not in flags
        and word not in ("-h", "--help")
    ]


def _commands_in(text: str) -> list[str]:
    joined = text.replace("\\\n", " ")
    return [*_INLINE.findall(joined), *_CODE_LINE.findall(joined)]


def test_every_command_in_the_documents_exists() -> None:
    """A command written in a document has to run, whatever the code became after it was written."""
    commands = _commands()
    problems = [
        f"{path.relative_to(ROOT)}: {written!r}: {problem}"
        for path in DOCUMENTS
        if path.name not in HISTORICAL
        for written in _commands_in(path.read_text(encoding="utf-8"))
        for problem in _unknown_words(written, commands)
    ]
    assert not problems, "\n".join(problems)


def test_every_command_in_a_message_exists() -> None:
    """The same rule for the commands the toolkit itself tells people to type."""
    commands = _commands()
    problems = [
        f"{path.relative_to(ROOT)}: {written!r}: {problem}"
        for path in sorted((ROOT / "src" / "plastid").rglob("*"))
        if path.suffix in (".py", ".md")
        for written in _commands_in(path.read_text(encoding="utf-8"))
        for problem in _unknown_words(written, commands)
    ]
    assert not problems, "\n".join(problems)


def test_the_document_check_reports_what_is_wrong() -> None:
    commands = _commands()
    assert _unknown_words("plastid nonsense", commands) == ["no command 'nonsense'"]
    assert _unknown_words("plastid run --frobnicate", commands) == ["run: no flag --frobnicate"]
    assert _unknown_words("plastid panel nonsense", commands) == ["no command 'panel nonsense'"]
    assert _unknown_words("plastid run -m ssr -i <file> --set window=1", commands) == []


def test_the_readme_states_every_exit_code() -> None:
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    missing = [
        state.code
        for state in EXIT_STATES
        if not re.search(rf"\| `{state.code}` \| `{state.status}` \|", readme)
    ]
    assert not missing


def test_modules_lists_everything_registered(capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["modules"]) == 0
    said = capsys.readouterr().out
    assert all(module_id in said for module_id in registry.all_modules())


def test_modules_describes_one_module(capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["modules", "substitutions"]) == 0
    said = capsys.readouterr().out
    assert "scope" in said
    assert "whole" in said
    assert "mafft" in said


def test_modules_refuses_an_unknown_one(capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["modules", "nonsense"]) == 2
    assert "nonsense" in capsys.readouterr().err


def test_modules_json_is_what_a_form_needs(capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["modules", "--json"]) == 0
    payload = json.loads(capsys.readouterr().out)
    by_id = {entry["id"]: entry for entry in payload["modules"]}
    assert set(by_id) == set(registry.all_modules())
    scope = next(p for p in by_id["substitutions"]["params"] if p["name"] == "scope")
    assert scope["type"] == "str"
    assert scope["default"] == "whole"
    assert set(scope["choices"]) == {"whole", "cds"}
    assert by_id["substitutions"]["genomes"] == {"min": 2, "max": None}
    assert by_id["genome_stats"]["genomes"] == {"min": 1, "max": 1}
    assert by_id["kaks"]["tools"] == ["mafft", "paml"]
    assert by_id["plastid_annotate"]["command"] == "plastid annotate"
    assert {param["name"] for param in payload["shared_params"]} >= {"data_label", "ir_source"}
    assert [state["code"] for state in payload["exit_states"]] == [
        state.code for state in EXIT_STATES
    ]


def test_runs_lists_a_folder(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    out = tmp_path / "runs"
    assert main(["run", "-m", "structure", "-i", GENOME, "-o", str(out), "--dataset", "seen"]) == 0
    capsys.readouterr()
    assert main(["runs", str(out)]) == 0
    assert "seen" in capsys.readouterr().out
    assert main(["runs", str(out), "--json"]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert [run["status"] for run in payload["runs"]] == ["ok"]


def test_runs_refuses_a_folder_that_does_not_exist(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    assert main(["runs", str(tmp_path / "absent")]) == 2
    assert "absent" in capsys.readouterr().err


def test_tools_names_every_tool_with_or_without_it(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    monkeypatch.setenv("PATH", str(tmp_path))
    assert main(["tools", "--json"]) == 0
    payload = json.loads(capsys.readouterr().out)
    by_name = {entry["name"]: entry for entry in payload["tools"]}
    assert by_name["mafft"]["status"] == "missing"
    assert by_name["mafft"]["pinned"]
    assert by_name["mafft"]["doc"].startswith("docs/tools/")
    assert "kaks" in by_name["mafft"]["used_by"]


def test_figures_lists_what_a_run_can_draw(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    out = tmp_path / "runs"
    assert main(["run", "-m", "structure", "-i", GENOME, "-o", str(out), "--dataset", "drawn"]) == 0
    capsys.readouterr()
    run_dir = next(out.glob("drawn/*"))
    assert main(["figures", str(run_dir), "--list"]) == 0
    said = capsys.readouterr().out
    assert "circular_map" in said or "gc_skew" in said
    assert "needs" in said


def test_figures_refuses_an_option_nothing_reads(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """A misspelt option that does nothing looks like a figure that ignored its option."""
    out = tmp_path / "runs"
    assert main(["run", "-m", "structure", "-i", GENOME, "-o", str(out), "--dataset", "drawn"]) == 0
    capsys.readouterr()
    run_dir = next(out.glob("drawn/*"))
    assert main(["figures", str(run_dir), "--set", "wdth=4"]) == 2
    assert "wdth" in capsys.readouterr().err


@pytest.mark.parametrize(
    ("module", "setting", "named"),
    [
        ("structure", "window=abc", "whole number"),
        ("structure", "ir_source=seq", "annotation, sequence"),
        ("rna_editing", "cutoff=high", "must be a number"),
        ("structure", "window=100", "step"),
    ],
)
def test_a_wrong_value_is_a_refused_command_line(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], module: str, setting: str, named: str
) -> None:
    """Exit 2 and a sentence, not a Python error, and nothing computed from the value."""
    code = main(["run", "-m", module, "-i", GENOME, "-o", str(tmp_path), "--set", setting])
    said = capsys.readouterr().err
    assert code == 2, said
    assert named in said
    assert "Traceback" not in said
    assert "ValueError" not in said


def test_a_key_a_flag_sets_is_refused_on_set(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    code = main(
        ["run", "-m", "structure", "-i", GENOME, "-o", str(tmp_path), "--set", "data_label=x"]
    )
    assert code == 2
    assert "--label" in capsys.readouterr().err


def test_annotation_has_one_door(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    """`plastid annotate` refuses an annotated record, audits the set and records it.

    `plastid run -m plastid_annotate` reached the same module past the first and the last of those,
    so one module left two kinds of run behind. The route that keeps the checks is the only one.
    """
    for command in ("run", "batch"):
        code = main([command, "-m", "plastid_annotate", "-i", GENOME, "-o", str(tmp_path)])
        assert code == 2
        assert "plastid annotate" in capsys.readouterr().err
    assert not list(tmp_path.iterdir())


def test_a_missing_input_is_a_refused_command_line(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    code = main(["run", "-m", "structure", "-i", str(tmp_path / "no.gb"), "-o", str(tmp_path)])
    assert code == 2
    assert "cannot read the input" in capsys.readouterr().err


def _message_strings(path: Path) -> list[tuple[str, str]]:
    """Every string a source file holds that is not a docstring, with where it is."""
    import ast

    tree = ast.parse(path.read_text(encoding="utf-8"))
    documented = {
        id(node.body[0].value)
        for node in ast.walk(tree)
        if isinstance(node, ast.Module | ast.FunctionDef | ast.ClassDef)
        and node.body
        and isinstance(node.body[0], ast.Expr)
        and isinstance(node.body[0].value, ast.Constant)
    }
    return [
        (f"{path.relative_to(ROOT)}:{node.lineno}", node.value)
        for node in ast.walk(tree)
        if isinstance(node, ast.Constant)
        and isinstance(node.value, str)
        and id(node) not in documented
    ]


def test_no_message_or_report_carries_a_planning_code() -> None:
    """What a run prints or writes reaches a reader who has no way to look such a code up."""
    package = ROOT / "src" / "plastid"
    texts = [
        item
        for path in sorted(package.rglob("*.py"))
        if not path.name.startswith("mito_")
        for item in _message_strings(path)
    ]
    texts += [
        (str(path.relative_to(ROOT)), path.read_text(encoding="utf-8"))
        for path in sorted((package / "report" / "templates").glob("*.md"))
    ]
    assert _maintainer_codes(texts) == []


def test_the_tool_table_covers_every_tool_that_can_be_run() -> None:
    """Every runtime row of the inventory is checked by `plastid tools`, or run by one that is."""
    from plastid.external.inventory import TOOL_CHECKS

    inside_another = {"bowtie2", "blastn-assembly", "hmmer", "spades"}
    table = (ROOT / "docs" / "tools" / "README.md").read_text(encoding="utf-8")
    rows = re.findall(r"^\| `([^`]+)` \| (active|assembly)-runtime \|", table, re.MULTILINE)
    assert rows
    assert {name for name, _role in rows} - inside_another == {check.name for check in TOOL_CHECKS}
    for check in TOOL_CHECKS:
        assert (ROOT / check.doc).is_file(), check.doc
        pinned = re.search(
            rf"^\| `{re.escape(check.name)}` \|[^|]*\|[^|]*\| ([^|]+) \|", table, re.M
        )
        assert pinned is not None and check.pinned in pinned.group(1), check.name


def test_every_subcommand_has_a_handler_and_every_handler_a_subcommand() -> None:
    """The dispatch table and the parser name the same commands, so none is parsed and dropped."""
    from plastid.cli import _COMMANDS, _subcommands

    assert set(_COMMANDS) == set(_subcommands(build_parser()))
