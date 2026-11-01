"""Every line a form can build is a line the real parser accepts, with the meaning the form had."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pytest

from plastid.cli import build_parser
from plastid.ui import commands
from plastid.ui.commands import BUILDERS, FormError, build, shown

ONE = {"modules": ["ssr"], "inputs": ["records/a.gb"], "label": "NCBI RefSeq", "dataset": "a"}
FULL_FORMS = {
    "run": {**ONE, "figures": True, "docx": True, "no_cache": True, "sets": {"window": 500}},
    "annotate": {"inputs": ["records/a.fasta"], "reference_set": "references/iridaceae"},
    "assemble run": {"reads_1": "r/a_1.fq.gz", "reads_2": "r/a_2.fq.gz", "sample": "a"},
    "assemble prepare": {"hifi": "r/a.fq.gz", "sample": "a"},
    "assemble import": {
        "prepared": "p.json",
        "result": "out",
        "source_name": "usegalaxy.eu",
        "retrieved_at": "2026-11-01",
    },
    "panel refresh": {"taxon": "Iridaceae"},
    "panel propose": {"target": "NC_050833.1", "mode": "all", "level": "family"},
    "panel pin": {"accessions": ["NC_050833.1", "NC_024936.1"], "to": "references/iridaceae"},
    "reference verify": {"reference_set": "references/iridaceae"},
    "reference audit": {"reference_set": "references/iridaceae"},
    "reference leave-one-out": {"reference_set": "references/iridaceae", "family": "Iridaceae"},
    "editing-reference build": {"to": "references/editing_reference.json"},
    "editing-reference verify": {},
    "figures": {"run": "runs/a/20261101-000000-abcdef12", "figures": ["r_rscu"]},
    "package": {"run": "runs/a/20261101-000000-abcdef12"},
    "qc": {"run": "runs/a/20261101-000000-abcdef12"},
    "ir-junction": {"inputs": ["records/a.gb", "records/b.gb"], "label": "NCBI RefSeq"},
}


RUN = "runs/a/20261101-000000-abcdef12"
OUT = {"out": "runs", "dataset": "a", "label": "NCBI RefSeq"}
SETS = {"sets": {"window": 500}}
# Each command's forms together reach every flag it has. Two forms where flags exclude each other.
EVERY_FLAG: dict[str, tuple[str, list[dict]]] = {
    "run": (
        "run",
        [{**ONE, **OUT, **SETS, "groups": "g.yml", "docx": 1, "figures": 1, "no_cache": 1}],
    ),
    "batch": (
        "run",
        [{**ONE, **OUT, **SETS, "one_by_one": 1, "docx": 1, "figures": 1, "no_cache": 1}],
    ),
    "annotate": (
        "annotate",
        [
            {
                **OUT,
                **SETS,
                "inputs": ["a.fasta"],
                "reference_set": "references",
                "exclude": ["NC_1.1"],
                "docx": 1,
                "no_cache": 1,
            }
        ],
    ),
    "assemble run": (
        "assemble run",
        [
            {
                **OUT,
                **SETS,
                "reads_1": "1.fq",
                "reads_2": "2.fq",
                "sample": "a",
                "compare_to": "b.gb",
                "skip_resource_check": 1,
            },
            {"hifi": "h.fq", "sample": "a"},
        ],
    ),
    "assemble prepare": (
        "assemble prepare",
        [
            {**OUT, **SETS, "reads_1": "1.fq", "reads_2": "2.fq", "sample": "a"},
            {"hifi": "h.fq", "sample": "a"},
        ],
    ),
    "assemble import": (
        "assemble import",
        [
            {
                **OUT,
                **SETS,
                "prepared": "p.json",
                "result": "out",
                "source_name": "usegalaxy.eu",
                "retrieved_at": "2026-11-01",
                "source_url": "https://usegalaxy.eu/h/1",
                "compare_to": "b.gb",
            }
        ],
    ),
    "panel refresh": (
        "panel refresh",
        [{"taxon": "Iridaceae", "index": "i.json", "max_records": 5}],
    ),
    "panel propose": (
        "panel propose",
        [{"target": "NC_1.1", "index": "i.json", "mode": "all", "level": "genus", "per_tier": 2}],
    ),
    "panel pin": (
        "panel pin",
        [
            {"accessions": ["NC_1.1"], "index": "i.json", "to": "references"},
            {"groups": "g.yml", "to": "references"},
            {"files": ["a.gb"], "note": "PMGA web, 2026-10-03", "to": "references"},
        ],
    ),
    "reference verify": ("reference verify", [{"reference_set": "references"}]),
    "reference audit": (
        "reference audit",
        [
            {
                "reference_set": "references",
                "out": "runs",
                "exceptions": "e.json",
                "no_confirm": 1,
                "confirm_exception": "a.gb:ir_split",
                "confirmed_by": "A. Reader",
                "reason": "the record's own notes say so",
            }
        ],
    ),
    "reference leave-one-out": (
        "reference leave-one-out",
        [
            {
                "reference_set": "references",
                "family": "Iridaceae",
                "rule": "same genus",
                "to": "s.json",
            }
        ],
    ),
    "editing-reference build": (
        "editing-reference build",
        [{"records": "records", "to": "e.json", "out": "runs"}],
    ),
    "editing-reference verify": (
        "editing-reference verify",
        [{"reference_file": "e.json", "records": "records"}],
    ),
    "figures": ("figures", [{"run": RUN, "figures": ["r_rscu"], **SETS}]),
    "package": ("package", [{"run": RUN}]),
    "qc": ("qc", [{"run": RUN}]),
    "ir-junction": (
        "ir-junction",
        [
            {
                **OUT,
                "inputs": ["a.gb"],
                "tag": "junctions",
                "width": 2400,
                "tree": "t.nwk",
                "order_from": "Iris domestica",
                "cpjsdraw": "CPJSdraw.pl",
                "drop_pseudo": 1,
            }
        ],
    ),
}
# A flag of no form: it prints a list for a script, and the page already holds that list.
LISTINGS = {"-h", "--help", "--json", "--list"}


def _leaf(command: str) -> argparse.ArgumentParser:
    parser = build_parser()
    for word in command.split():
        found = [a for a in parser._actions if isinstance(a, argparse._SubParsersAction)]
        parser = found[0].choices[word]
    return parser


def test_every_builder_has_a_form_here() -> None:
    assert set(FULL_FORMS) == set(BUILDERS)


@pytest.mark.parametrize("command", sorted(EVERY_FLAG))
def test_every_flag_can_be_reached_from_a_form(command: str) -> None:
    """A command with a page is not covered while one of its flags has no way in."""
    builder, forms = EVERY_FLAG[command]
    built = [build(builder, form) for form in forms]
    for args in built:
        build_parser().parse_args(args)
    typed = {word for args in built for word in args}
    unreachable = [
        action.option_strings
        for action in _leaf(command)._actions
        if action.option_strings and not set(action.option_strings) & (typed | LISTINGS)
    ]
    assert unreachable == []


def test_every_field_a_builder_reads_is_on_a_page() -> None:
    """The builders are tested above; this is what ties them to a field a person can fill in."""
    builders = Path(commands.__file__).read_text(encoding="utf-8")
    page = (Path(commands.__file__).parent / "assets" / "app.js").read_text(encoding="utf-8")
    fields = set(re.findall(r'form(?:\.get\(|\[)"(\w+)"', builders))
    fields |= set(re.findall(r'_needs\(form, "(\w+)"', builders))
    assert fields and {field for field in fields if not re.search(rf"\b{field}\b", page)} == set()


@pytest.mark.parametrize("command", sorted(FULL_FORMS))
def test_a_full_form_passes_the_real_parser(command: str) -> None:
    build_parser().parse_args(build(command, FULL_FORMS[command]))


def test_many_by_one_is_the_batch_command() -> None:
    args = build("run", {**ONE, "one_by_one": True})
    assert args[0] == "batch"
    build_parser().parse_args(args)


def test_a_form_says_what_it_still_needs() -> None:
    with pytest.raises(FormError, match="at least one analysis"):
        build("run", {"inputs": ["a.gb"]})
    with pytest.raises(FormError, match="the first reads"):
        build("assemble run", {"sample": "a"})


def test_an_unknown_command_is_refused_by_name() -> None:
    with pytest.raises(FormError, match="no such command"):
        build("rm -rf", {})


def test_the_shown_line_is_quoted_as_a_shell_needs() -> None:
    line = shown(["run", "--label", "NCBI RefSeq, 2026"])
    assert line == "plastid run --label 'NCBI RefSeq, 2026'"


def test_a_value_cannot_become_a_second_argument() -> None:
    """A label with spaces stays one value: the form never goes through a shell."""
    args = build("run", {**ONE, "label": "a --docx b"})
    assert "--docx" not in args
    assert args[args.index("--label") + 1] == "a --docx b"
