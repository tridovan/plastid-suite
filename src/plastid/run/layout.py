"""Run-scoped output layout (PROJECT.md §2, requirement 02/10).

Every run gets its own directory: ``<out_root>/<dataset_slug>/<run_id>/``. Modules never
name paths themselves — the layout owns naming (seam 5).
"""

from __future__ import annotations

import re
import unicodedata
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

__all__ = ["FIGURES_DIRNAME", "RunLayout", "slugify"]

# One directory for the run's figures, beside the numbered module directories. Figures are a
# property of the run and not of one module — two of them read tables from two modules — so
# collecting them keeps one place for the user to look.
FIGURES_DIRNAME = "figures"

_UNSAFE = re.compile(r"[^A-Za-z0-9._-]+")


def slugify(text: str) -> str:
    """ASCII-safe path segment: strip diacritics, collapse unsafe runs to ``-``.

    Vietnamese accents in paths caused real cross-machine breakage (PROJECT.md §9.30).
    """
    ascii_text = unicodedata.normalize("NFKD", text).encode("ascii", "ignore").decode("ascii")
    slug = _UNSAFE.sub("-", ascii_text).strip("-._").lower()
    if not slug:
        raise ValueError(f"cannot build a path slug from {text!r}")
    return slug


@dataclass(frozen=True, slots=True)
class RunLayout:
    """A run's directory tree. One run = one directory, always new, never overwritten."""

    out_root: Path
    dataset: str
    run_id: str

    @property
    def run_dir(self) -> Path:
        return self.out_root / self.dataset / self.run_id

    @property
    def run_info_dir(self) -> Path:
        return self.run_dir / "00_run_info"

    @property
    def logs_dir(self) -> Path:
        return self.run_info_dir / "logs"

    @property
    def cache_dir(self) -> Path:
        return self.out_root / ".cache"

    @property
    def report_dir(self) -> Path:
        return self.run_dir / "report"

    @property
    def figures_dir(self) -> Path:
        return self.run_dir / FIGURES_DIRNAME

    def module_dir(self, order: int, module_id: str) -> Path:
        """Numbered by run order — assigned by the orchestrator, never by a module."""
        return self.run_dir / f"{order:02d}_{slugify(module_id)}"

    def make_dirs(self) -> None:
        for directory in (self.run_info_dir, self.logs_dir, self.report_dir):
            directory.mkdir(parents=True, exist_ok=True)

    @classmethod
    def create(
        cls,
        out_root: str | Path,
        dataset: str,
        *,
        fingerprint: str,
        now: datetime | None = None,
    ) -> RunLayout:
        stamp = (now or datetime.now()).strftime("%Y%m%d-%H%M%S")
        root = Path(out_root)
        dataset_slug = slugify(dataset)
        base = f"{stamp}-{fingerprint[:8]}"
        run_id = base
        suffix = 2
        while (root / dataset_slug / run_id).exists():
            # Two runs inside the same second must still get their own directory —
            # one run = one folder, always (PROJECT.md §2).
            run_id = f"{base}-{suffix}"
            suffix += 1
        return cls(out_root=root, dataset=dataset_slug, run_id=run_id)
