"""What a form holds, as the command line that does the work.

The interface starts `plastid` itself, so a run made here is a run the command line could have made:
one engine, one kind of run directory, and the line shown beside a form is the line that runs. Each
builder returns the arguments after the word `plastid`; `tests/unit/test_ui_commands.py` puts every
one through the real parser.
"""

from __future__ import annotations

import shlex
from collections.abc import Callable, Mapping, Sequence
from typing import Any

__all__ = ["BUILDERS", "FormError", "build", "shown"]

Form = Mapping[str, Any]


class FormError(ValueError):
    """A form that cannot become a command yet, in words the form can show."""


def _value(args: list[str], flag: str, value: Any) -> None:
    if value not in (None, ""):
        args += [flag, str(value)]


def _each(args: list[str], flag: str, values: Sequence[str]) -> None:
    for value in values:
        args += [flag, str(value)]


def _switch(args: list[str], flag: str, value: Any) -> None:
    if value:
        args.append(flag)


def _sets(args: list[str], form: Form) -> None:
    sets = form.get("sets") or {}
    for key in sorted(sets):
        if sets[key] not in (None, ""):
            args += ["--set", f"{key}={sets[key]}"]


def _destination(args: list[str], form: Form) -> None:
    _value(args, "-o", form.get("out"))
    _value(args, "--dataset", form.get("dataset"))
    _value(args, "--label", form.get("label"))


def _needs(form: Form, key: str, words: str) -> Any:
    value = form.get(key)
    if value in (None, "", []):
        raise FormError(words)
    return value


def run(form: Form) -> list[str]:
    modules = _needs(form, "modules", "choose at least one analysis")
    inputs = _needs(form, "inputs", "choose at least one record")
    args = ["batch" if form.get("one_by_one") else "run"]
    _each(args, "-m", modules)
    _each(args, "-i", inputs)
    _destination(args, form)
    if not form.get("one_by_one"):
        _value(args, "--groups", form.get("groups"))
    _switch(args, "--docx", form.get("docx"))
    _switch(args, "--figures", form.get("figures"))
    _switch(args, "--no-cache", form.get("no_cache"))
    _sets(args, form)
    return args


def annotate(form: Form) -> list[str]:
    inputs = _needs(form, "inputs", "choose a sequence to annotate")
    args = ["annotate"]
    _each(args, "-i", inputs)
    _value(args, "--reference-set", _needs(form, "reference_set", "choose a reference set"))
    _each(args, "--exclude", form.get("exclude") or ())
    _destination(args, form)
    _switch(args, "--docx", form.get("docx"))
    _switch(args, "--no-cache", form.get("no_cache"))
    _sets(args, form)
    return args


def _reads(args: list[str], form: Form) -> None:
    if form.get("hifi"):
        args += ["--hifi", str(form["hifi"])]
    else:
        args += ["-1", str(_needs(form, "reads_1", "choose the first reads"))]
        args += ["-2", str(_needs(form, "reads_2", "choose the second reads"))]
    _value(args, "--sample", _needs(form, "sample", "name the sample"))


def assemble_run(form: Form) -> list[str]:
    args = ["assemble", "run"]
    _reads(args, form)
    _destination(args, form)
    _value(args, "--compare-to", form.get("compare_to"))
    _sets(args, form)
    return args


def assemble_prepare(form: Form) -> list[str]:
    args = ["assemble", "prepare"]
    _reads(args, form)
    _destination(args, form)
    _sets(args, form)
    return args


def assemble_import(form: Form) -> list[str]:
    args = ["assemble", "import"]
    _value(args, "--prepared", _needs(form, "prepared", "choose the prepared.json"))
    _value(args, "--result", _needs(form, "result", "choose the folder of the assembler's output"))
    _value(args, "--source-name", _needs(form, "source_name", "say where it ran"))
    _value(args, "--retrieved-at", _needs(form, "retrieved_at", "give the date it was retrieved"))
    _value(args, "--source-url", form.get("source_url"))
    _destination(args, form)
    _value(args, "--compare-to", form.get("compare_to"))
    _sets(args, form)
    return args


def panel_refresh(form: Form) -> list[str]:
    args = ["panel", "refresh", "--taxon", str(_needs(form, "taxon", "name a taxon"))]
    _value(args, "--index", form.get("index"))
    _value(args, "--max-records", form.get("max_records"))
    return args


def panel_propose(form: Form) -> list[str]:
    args = ["panel", "propose", "--target", str(_needs(form, "target", "name the target genome"))]
    _value(args, "--index", form.get("index"))
    _value(args, "--mode", form.get("mode"))
    _value(args, "--level", form.get("level"))
    _value(args, "--per-tier", form.get("per_tier"))
    return args


def panel_pin(form: Form) -> list[str]:
    args = ["panel", "pin"]
    _value(args, "--index", form.get("index"))
    if form.get("files"):
        _each(args, "--file", form["files"])
    else:
        args += [
            "--accessions",
            ",".join(_needs(form, "accessions", "choose the genomes to fetch")),
        ]
    _value(args, "--note", form.get("note"))
    _value(args, "--to", _needs(form, "to", "choose where the records go"))
    return args


def reference_verify(form: Form) -> list[str]:
    args = ["reference", "verify"]
    _value(args, "--reference-set", _needs(form, "reference_set", "choose a reference set"))
    return args


def reference_audit(form: Form) -> list[str]:
    args = ["reference", "audit"]
    _value(args, "--reference-set", _needs(form, "reference_set", "choose a reference set"))
    _value(args, "-o", form.get("out"))
    return args


def reference_leave_one_out(form: Form) -> list[str]:
    args = ["reference", "leave-one-out"]
    _value(args, "--reference-set", _needs(form, "reference_set", "choose a reference set"))
    _value(args, "--family", _needs(form, "family", "name the family the set was chosen at"))
    _value(args, "--rule", form.get("rule"))
    _value(args, "--to", form.get("to"))
    return args


def editing_reference_build(form: Form) -> list[str]:
    args = ["editing-reference", "build"]
    _value(args, "--records", form.get("records"))
    _value(args, "--to", _needs(form, "to", "choose where the reference is written"))
    _value(args, "-o", form.get("out"))
    return args


def editing_reference_verify(form: Form) -> list[str]:
    return ["editing-reference", "verify"]


def figures(form: Form) -> list[str]:
    args = ["figures", str(_needs(form, "run", "choose a run"))]
    _each(args, "--figure", form.get("figures") or ())
    _sets(args, form)
    return args


def package(form: Form) -> list[str]:
    return ["package", str(_needs(form, "run", "choose a run"))]


def qc(form: Form) -> list[str]:
    return ["qc", str(_needs(form, "run", "choose a run"))]


def ir_junction(form: Form) -> list[str]:
    args = ["ir-junction"]
    _each(args, "-i", _needs(form, "inputs", "choose the annotated records"))
    _destination(args, form)
    _value(args, "--width", form.get("width"))
    _value(args, "--tree", form.get("tree"))
    _value(args, "--order-from", form.get("order_from"))
    _value(args, "--cpjsdraw", form.get("cpjsdraw"))
    _switch(args, "--drop-pseudo", form.get("drop_pseudo"))
    return args


# Every command the interface can start, by the words typed after `plastid`.
BUILDERS: dict[str, Callable[[Form], list[str]]] = {
    "run": run,
    "annotate": annotate,
    "assemble run": assemble_run,
    "assemble prepare": assemble_prepare,
    "assemble import": assemble_import,
    "panel refresh": panel_refresh,
    "panel propose": panel_propose,
    "panel pin": panel_pin,
    "reference verify": reference_verify,
    "reference audit": reference_audit,
    "reference leave-one-out": reference_leave_one_out,
    "editing-reference build": editing_reference_build,
    "editing-reference verify": editing_reference_verify,
    "figures": figures,
    "package": package,
    "qc": qc,
    "ir-junction": ir_junction,
}


def build(command: str, form: Form) -> list[str]:
    """The arguments after `plastid` for one form; an unknown command is the caller's error."""
    try:
        builder = BUILDERS[command]
    except KeyError:
        raise FormError(f"no such command: {command}") from None
    return builder(form)


def shown(args: Sequence[str]) -> str:
    """The command as a person would type it, quoted the way their shell needs."""
    return shlex.join(["plastid", *args])
