"""Module registry — a decorator and a dict. Explicit imports, no entry points (yet).

Adding an analysis = adding a file in ``modules/`` plus one import line in
``modules/__init__.py`` (CLAUDE.md: extend by adding, not editing).
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from typing import Any

from plastid.core.model import Genome
from plastid.result.result import Result
from plastid.run.context import RunContext

__all__ = ["ModuleSpec", "UnknownModuleError", "all_modules", "get", "module"]

ModuleFn = Callable[[Genome, dict[str, Any], RunContext], Result]


class UnknownModuleError(KeyError):
    """No module is registered under that id."""


@dataclass(frozen=True, slots=True)
class ModuleSpec:
    id: str
    version: str
    title: str
    run: ModuleFn
    requires_ir: bool = False
    citations: tuple[str, ...] = ()


_REGISTRY: dict[str, ModuleSpec] = {}


def module(
    module_id: str,
    *,
    version: str,
    title: str,
    requires_ir: bool = False,
    citations: tuple[str, ...] = (),
) -> Callable[[ModuleFn], ModuleFn]:
    """Register a module's ``run`` function."""

    def decorate(run_fn: ModuleFn) -> ModuleFn:
        if module_id in _REGISTRY:
            raise ValueError(f"module {module_id!r} is already registered")
        _REGISTRY[module_id] = ModuleSpec(
            id=module_id,
            version=version,
            title=title,
            run=run_fn,
            requires_ir=requires_ir,
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
