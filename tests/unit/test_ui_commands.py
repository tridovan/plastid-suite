"""Every line a form can build is a line the real parser accepts, with the meaning the form had."""

from __future__ import annotations

import pytest

from plastid.cli import build_parser
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


def test_every_builder_has_a_form_here() -> None:
    assert set(FULL_FORMS) == set(BUILDERS)


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
