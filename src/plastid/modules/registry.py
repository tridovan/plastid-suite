"""Module registry — a decorator and a dict. Explicit imports, no entry points (yet).

Adding an analysis = adding a file in ``modules/`` plus one import line in
``modules/__init__.py`` (CLAUDE.md: extend by adding, not editing).
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import Any

from plastid.core.model import Genome
from plastid.result.result import Result
from plastid.run.context import RunContext

__all__ = [
    "ModuleSpec",
    "UnknownModuleError",
    "all_modules",
    "get",
    "module",
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


_REGISTRY: dict[str, ModuleSpec] = {}


def module(
    module_id: str,
    *,
    version: str,
    title: str,
    requires_ir: bool = False,
    arity: int = 1,
    citations: tuple[str, ...] = (),
) -> Callable[[ModuleFn], ModuleFn]:
    """Register a module's ``run`` function.

    ``arity=1`` keeps the frozen v0.1 contract ``run(genome, params, ctx)``; a larger arity
    declares the v0.2 multi-genome contract ``run(genomes, params, ctx)``. A module never
    receives genomes it does not use, so the two contracts stay separate.
    """

    def decorate(run_fn: ModuleFn) -> ModuleFn:
        if module_id in _REGISTRY:
            raise ValueError(f"module {module_id!r} is already registered")
        if arity < 1:
            raise ValueError(f"module {module_id!r}: arity must be at least 1, got {arity}")
        _REGISTRY[module_id] = ModuleSpec(
            id=module_id,
            version=version,
            title=title,
            run=run_fn,
            requires_ir=requires_ir,
            arity=arity,
            citations=citations,
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
