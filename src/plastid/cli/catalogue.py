"""The listing commands: `plastid modules`, `tools` and `runs`."""

from __future__ import annotations

import argparse
import importlib.util
import sys
from dataclasses import asdict
from pathlib import Path
from typing import Any

from plastid.cli.common import EXIT_OK, EXIT_USAGE, _print_json, _print_table
from plastid.modules import assembly as assembly_step
from plastid.modules import registry
from plastid.run.state import EXIT_STATES, list_runs

# The tools no module calls: which command does.
_TOOL_COMMANDS = {"R": "plastid figures", "cpjsdraw": "plastid ir-junction", "perl": "ssr"}


def _param_record(param: registry.ParamSpec) -> dict[str, Any]:
    return {
        "name": param.name,
        "type": param.type,
        "help": param.help,
        "default": param.default,
        "choices": list(param.choices),
        "flag": param.flag or None,
    }


def _module_record(spec: registry.ModuleSpec) -> dict[str, Any]:
    """One module as data: what a person reads, and what a form is built from."""
    return {
        "id": spec.id,
        "version": spec.version,
        "title": spec.title,
        "group": spec.group,
        "genomes": {"min": spec.arity, "max": 1 if spec.arity == 1 else None},
        "ir_gate": spec.requires_ir,
        "tools": list(spec.tools),
        "optional_tools": list(spec.optional_tools),
        "command": spec.command or "plastid run",
        "in_batch": spec.arity == 1 and spec.command is None,
        "reference_param": spec.reference_param,
        "params": [_param_record(param) for param in spec.params],
    }


def _assembly_record() -> dict[str, Any]:
    """The assembly step in the shape of a module: `plastid assemble` is not a registered one."""
    return {
        "id": assembly_step.MODULE_ID,
        "version": assembly_step.MODULE_VERSION,
        "title": "Plastome assembly from short or HiFi reads, checked against the reads",
        "group": "plastid",
        "tools": ["getorganelle", "fastp", "minimap2"],
        "optional_tools": ["oatk"],
        "command": "plastid assemble run",
        "params": [_param_record(param) for param in assembly_step.ASSEMBLY_PARAM_SPECS],
    }


def _describe_param(param: registry.ParamSpec) -> None:
    if param.flag:
        setting = f"set with {param.flag}"
    else:
        default = "not set" if param.default is None else param.default
        choices = f" · one of {', '.join(str(c) for c in param.choices)}" if param.choices else ""
        setting = f"{param.type} · default {default}{choices}"
    print(f"  {param.name}: {setting}")
    print(f"      {param.help}")


def _describe_module(spec: registry.ModuleSpec) -> None:
    usage = (spec.command or f"plastid run -m {spec.id}") + " -i <file>"
    genomes = "one" if spec.arity == 1 else f"at least {spec.arity}"
    print(f"{spec.id} {spec.version}: {spec.title}")
    print(f"  run with: {usage}")
    print(f"  genomes:  {genomes}")
    print(f"  IR gate:  {'yes' if spec.requires_ir else 'no'}")
    print(f"  tools:    {', '.join(spec.tools) or 'none'}")
    if spec.optional_tools:
        print(f"  if asked: {', '.join(spec.optional_tools)}")
    if spec.params:
        print("\nparameters (--set KEY=VALUE):")
        for param in spec.params:
            _describe_param(param)
    print("\nalso accepted on any run:")
    for param in registry.SHARED_PARAMS:
        _describe_param(param)


def _cmd_modules(args: argparse.Namespace) -> int:
    if args.module == assembly_step.MODULE_ID:
        record = _assembly_record()
        if args.json:
            _print_json(record)
            return EXIT_OK
        print(f"{record['id']} {record['version']}: {record['title']}")
        print(f"  run with: {record['command']}")
        print(f"  tools:    {', '.join(record['tools'])}; for HiFi reads, oatk")
        print("\nparameters (--set KEY=VALUE):")
        for param in assembly_step.ASSEMBLY_PARAM_SPECS:
            _describe_param(param)
        return EXIT_OK
    if args.module:
        try:
            spec = registry.get(args.module)
        except registry.UnknownModuleError as error:
            print(str(error), file=sys.stderr)
            return EXIT_USAGE
        if args.json:
            _print_json(_module_record(spec))
        else:
            _describe_module(spec)
        return EXIT_OK
    specs = [spec for _name, spec in sorted(registry.all_modules().items())]
    if args.json:
        _print_json(
            {
                "modules": [_module_record(spec) for spec in specs],
                "assembly": _assembly_record(),
                "shared_params": [_param_record(param) for param in registry.SHARED_PARAMS],
                "exit_states": [asdict(state) for state in EXIT_STATES],
            }
        )
        return EXIT_OK
    _print_table(
        ("MODULE", "GENOMES", "IR GATE", "TOOLS", "WHAT IT DOES"),
        [
            (
                spec.id,
                "1" if spec.arity == 1 else f"{spec.arity}+",
                "yes" if spec.requires_ir else "no",
                ", ".join(spec.tools) or "-",
                spec.title,
            )
            for spec in specs
        ],
    )
    print(
        "\n`plastid modules <module>` shows one module and its parameters; "
        "`plastid modules assembly` shows the assembly step."
    )
    return EXIT_OK


# The optional parts a form can ask for, each by the module that tells whether it is installed.
_EXTRAS = (("figures", "matplotlib"), ("docx", "docx"))


def _cmd_tools(args: argparse.Namespace) -> int:
    """What a run on this machine would find: a report, so a missing tool is not a failure here."""
    from plastid.external.inventory import check_tools

    needed_by: dict[str, list[str]] = {}
    for name, spec in sorted(registry.all_modules().items()):
        for tool in spec.tools:
            needed_by.setdefault(tool, []).append(name)
        for tool in spec.optional_tools:
            needed_by.setdefault(tool, []).append(f"{name} (if asked)")
    statuses = check_tools()
    machine = assembly_step.measure_machine(Path())
    records = [
        {
            **status.as_dict(),
            "used_by": needed_by.get(status.name)
            or (
                ["plastid assemble"] if status.role == "assembly" else [_TOOL_COMMANDS[status.name]]
            ),
        }
        for status in statuses
    ]
    if args.json:
        extras = {name: importlib.util.find_spec(module) is not None for name, module in _EXTRAS}
        _print_json({"tools": records, "machine": asdict(machine), "extras": extras})
        return EXIT_OK
    print(
        f"this machine: {machine.ram_gb:.1f} GB memory, {machine.cores} cores, "
        f"{machine.disk_free_gb:.1f} GB free here\n"
    )
    _print_table(
        ("TOOL", "STATUS", "VERSION", "PINNED", "NEEDED BY"),
        [
            (
                record["name"],
                record["status"],
                record["version"] or "-",
                record["pinned"],
                ", ".join(record["used_by"]),
            )
            for record in records
        ],
    )
    for record in records:
        if record["status"] != "ok":
            print(f"\n{record['message']}")
    return EXIT_OK


def _cmd_runs(args: argparse.Namespace) -> int:
    root = Path(args.root or ".runs")
    if args.root and not root.is_dir():
        print(f"{root}: no such directory", file=sys.stderr)
        return EXIT_USAGE
    runs = list_runs(root) if root.is_dir() else []
    if args.json:
        _print_json({"runs": [state.as_dict() for state in runs]})
        return EXIT_OK
    if not runs:
        print(f"no runs in {root}")
        return EXIT_OK
    _print_table(
        ("STATUS", "DONE", "STARTED", "DATASET", "RUN"),
        [
            (
                state.status,
                f"{state.done}/{state.total}" if state.total is not None else "-",
                (state.started_at or "-").replace("T", " "),
                state.dataset,
                state.run_id,
            )
            for state in runs
        ],
    )
    return EXIT_OK
