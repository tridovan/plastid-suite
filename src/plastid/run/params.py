"""Parameter-value errors shared by adapters, modules and the CLI.

The registry declares the parameters a module accepts (M11.20); this is the error a module or an
adapter raises when a **declared** key carries a value it cannot honour. It lives in ``run/``
because both ``external/`` adapters and ``modules/`` may import that layer, while an adapter may
not import ``modules.registry`` (the layering guard in ``tests/unit/test_architecture.py``).

It exits as a usage error (2), next to the IR gate: the command line was well-formed enough to
reach the module, and the user can fix the value. It is raised before any output is produced from
the refused value.
"""

from __future__ import annotations

__all__ = ["ParamValueError"]


class ParamValueError(ValueError):
    """A declared parameter was given a value the module cannot honour."""
