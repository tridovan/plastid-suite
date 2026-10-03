"""Analysis modules — one job per file.

Module files never import each other; this package registers them explicitly.
The architecture test permits sibling imports only from this file.
"""

from __future__ import annotations

from plastid.modules import genome_stats, ssr  # noqa: F401  (importing registers the modules)

__all__ = ["genome_stats", "ssr"]
