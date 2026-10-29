"""Module registry — a decorator and a dict. Explicit imports, no entry points (yet).

Adding an analysis = adding a file in ``modules/`` plus one import line in
``modules/__init__.py`` (CLAUDE.md: extend by adding, not editing).
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import Any

from plastid.core.model import Genome
from plastid.core.organelle import GROUP_SHARED, GROUPS
from plastid.result.result import Result
from plastid.run.context import RunContext
from plastid.run.params import ParamValueError

__all__ = [
    "ModuleSpec",
    "ParamSpec",
    "ParamValueError",
    "SHARED_PARAMS",
    "UnknownModuleError",
    "all_modules",
    "get",
    "module",
    "unknown_params",
]

SingleGenomeFn = Callable[[Genome, dict[str, Any], RunContext], Result]
MultiGenomeFn = Callable[[Sequence[Genome], dict[str, Any], RunContext], Result]
ModuleFn = SingleGenomeFn | MultiGenomeFn


class UnknownModuleError(KeyError):
    """No module is registered under that id."""


@dataclass(frozen=True, slots=True)
class ModuleSpec:
    """One registered module. ``arity`` is how many genomes its ``run`` takes (API v0.2)."""

    id: str
    version: str
    title: str
    run: ModuleFn
    requires_ir: bool = False
    arity: int = 1
    citations: tuple[str, ...] = ()
    # Which organelle this module serves: shared · plastid-only · mito-only (H-09).
    # It is deliberately *not* `requires_ir`:
    # that flag answers "does this need the inverted repeat", which `gene_table` (a plastome table
    # that needs no IR) and `codon_usage` (CDS sequences only, no positions — and since M11.5 no
    # IR requirement either) both show is a different question. Default `shared` — a module that
    # says nothing is claiming to run on either organelle, which is the claim the run then checks.
    group: str = GROUP_SHARED
    # The parameter naming a *file* this module works against — a privileged record rather than
    # one of the run's equal inputs (D5-09: the naming reference). Declared here so the CLI can
    # resolve it, audit it and record its digest without knowing which module is asking: the
    # module says what it needs, the command line says how it is supplied.
    reference_param: str | None = None
    # Every parameter this module reads, so the CLI can refuse a key it does not (M11.20).
    # An empty tuple is a claim — "this module takes no parameter" — and a test compares it with
    # the keys the module's source actually reads. Do not add a key here without a code path
    # that honours it: a recorded parameter that does not reach the tool states something false.
    params: tuple[ParamSpec, ...] = ()


@dataclass(frozen=True, slots=True)
class ParamSpec:
    """One parameter a module accepts: its name, the type the CLI validates it as, and why.

    The declaration exists so the command line can refuse a key the module does not read
    (M11.20). Before it, ``--set min_lengh_bp=60`` was accepted, recorded in the manifest and
    ignored — a misspelling that looked like a successful override. ``help`` is not decoration:
    it is what the refusal message shows, so a user who typed the wrong key is told which one
    would work.

    ``type`` is a *name*, not a Python callable, because the registry is data a CLI reads and a
    callable would make ``ModuleSpec`` harder to compare in tests. The CLI owns the conversion.
    """

    name: str
    type: str
    help: str


# The parameters every module sees, declared once here rather than repeated in twenty specs.
# They are run-level: `data_label` puts the run's provenance on every table (CLAUDE.md: every
# exported number carries a scope label) and `ir_source` chooses the IR path before any module
# runs. A module does not have to declare them, and the CLI accepts them for any module.
SHARED_PARAMS: tuple[ParamSpec, ...] = (
    ParamSpec(
        name="data_label",
        type="str",
        help="provenance label written into every table's scope (--label on the command line)",
    ),
    ParamSpec(
        name="ir_source",
        type="str",
        help="which IR path the run uses before any module runs: annotation (default) or sequence",
    ),
    # Read by the CLI, not by a module: `annotate` refuses an already-annotated GenBank unless
    # reannotate=yes, and `figures`/`package` take the R binary path. They are declared here
    # because the CLI is also a place a --set key is consumed, and a key the CLI reads but
    # nothing declares is exactly the silent no-op M11.20 removes.
    ParamSpec(
        name="reannotate",
        type="bool",
        help="annotate: yes to re-annotate a GenBank that already carries an annotation",
    ),
    ParamSpec(
        name="rscript",
        type="str",
        help="figures/package: path to Rscript when it is not on PATH",
    ),
)
_SHARED_NAMES = frozenset(param.name for param in SHARED_PARAMS)


def unknown_params(module_ids: Sequence[str], params: dict[str, Any]) -> list[str]:
    """Every key in ``params`` that none of ``module_ids`` accepts, as one readable line each.

    A run may name several modules; a key is accepted when at least one of them declares it or
    it is a shared parameter. The message names the offending key and the keys that *would*
    work, because a refusal that leaves the user guessing produces the same mistake twice.
    """
    accepted: dict[str, str] = {}
    for module_id in module_ids:
        spec = _REGISTRY.get(module_id)
        if spec is None:
            continue
        for param in (*spec.params, *SHARED_PARAMS):
            accepted.setdefault(param.name, param.help)
        # A module may declare a *privileged record* instead of listing it among its ordinary
        # parameters (D5-09: `naming_reference`). It reaches the module through `params` like
        # any other key, so it has to be accepted here too — measured 28/10/2026: without this
        # line `--set naming_reference=<file>` was refused by name and every calibration run
        # stopped at exit 2.
        if spec.reference_param:
            accepted.setdefault(
                spec.reference_param, "the privileged record this module works against"
            )

    problems: list[str] = []
    for key in sorted(params):
        if key in accepted:
            continue
        who = ", ".join(sorted(module_ids))
        choices = ", ".join(sorted(accepted)) or "none"
        problems.append(f"--set {key}: {who} does not take that parameter (accepted: {choices})")
    return problems


_REGISTRY: dict[str, ModuleSpec] = {}


def module(
    module_id: str,
    *,
    version: str,
    title: str,
    group: str,
    requires_ir: bool = False,
    arity: int = 1,
    citations: tuple[str, ...] = (),
    reference_param: str | None = None,
    params: tuple[ParamSpec, ...] = (),
) -> Callable[[ModuleFn], ModuleFn]:
    """Register a module's ``run`` function.

    ``arity=1`` keeps the frozen v0.1 contract ``run(genome, params, ctx)``; a larger arity
    declares the v0.2 multi-genome contract ``run(genomes, params, ctx)``. A module never
    receives genomes it does not use, so the two contracts stay separate.

    ``group`` has no default on purpose: a module that declares nothing would be claiming to run
    on either organelle, and an omission is not a claim. Making
    it required means the claim is made at the one place that already fails loudly — import time.
    """

    def decorate(run_fn: ModuleFn) -> ModuleFn:
        if module_id in _REGISTRY:
            raise ValueError(f"module {module_id!r} is already registered")
        if arity < 1:
            raise ValueError(f"module {module_id!r}: arity must be at least 1, got {arity}")
        if group not in GROUPS:
            raise ValueError(
                f"module {module_id!r}: group must be one of {', '.join(GROUPS)}, got {group!r}"
            )
        _REGISTRY[module_id] = ModuleSpec(
            id=module_id,
            version=version,
            title=title,
            run=run_fn,
            requires_ir=requires_ir,
            arity=arity,
            citations=citations,
            group=group,
            reference_param=reference_param,
            params=params,
        )
        return run_fn

    return decorate


def get(module_id: str) -> ModuleSpec:
    try:
        return _REGISTRY[module_id]
    except KeyError:
        raise UnknownModuleError(
            f"unknown module {module_id!r} — available: {sorted(_REGISTRY)}"
        ) from None


def all_modules() -> dict[str, ModuleSpec]:
    return dict(_REGISTRY)
