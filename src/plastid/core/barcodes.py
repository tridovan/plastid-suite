"""Hypervariable-region calling over a π window series (P4-04, M4.3).

The series itself is not computed here — it is :func:`plastid.core.align.window_diversity`, the
same function :mod:`plastid.modules.diversity` reports, called with the same window and step. This
module only *calls regions* from it, and the calling rule is written down because a
"hypervariable region" without its rule is a shape somebody eyeballed:

* the **threshold** is ``mean + k · sd`` over the windows that have a value (k is a parameter,
  default 1) — computed from the run's own data rather than fixed across genomes, so a conserved
  pair and a divergent pair each get a threshold that means something for them;
* a **region** is a maximal run of consecutive windows above the threshold and at least
  ``min_windows`` long, so a single noisy window is not a region;
* windows with no usable column have no value and are **not** part of any region — a hole is not
  a peak.

Every region carries the window and step it was found at; without them the coordinates mean
nothing.
"""

from __future__ import annotations

import statistics
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import yaml

from plastid.core.align import WindowDiversity
from plastid.core.resources import data_file

__all__ = [
    "BarcodeLocus",
    "BarcodeLocusError",
    "VariableRegion",
    "hypervariable_regions",
    "load_barcode_loci",
    "threshold_for",
]

DEFAULT_THRESHOLD_SD = 1.0
DEFAULT_MIN_WINDOWS = 2

KINDS = ("gene", "intergenic", "nuclear")


class BarcodeLocusError(ValueError):
    """The barcode-locus file is malformed — a locus without a source, or an unknown kind."""


@dataclass(frozen=True, slots=True)
class BarcodeLocus:
    """One candidate locus: what it is, which features it needs, and where the claim comes from."""

    name: str
    kind: str  # "gene" | "intergenic" | "nuclear"
    members: tuple[str, ...]
    source: str


def load_barcode_loci(path: str | Path | None = None) -> tuple[BarcodeLocus, ...]:
    """Read ``data/barcode_loci.yaml``; every entry must carry its source, like the alias table."""
    source_path = Path(path) if path is not None else data_file("barcode_loci.yaml")
    payload = yaml.safe_load(source_path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or not isinstance(payload.get("loci"), list):
        raise BarcodeLocusError(f"{source_path}: expected a 'loci' list")
    loci: list[BarcodeLocus] = []
    seen: set[str] = set()
    for entry in payload["loci"]:
        if not isinstance(entry, dict):
            raise BarcodeLocusError(f"{source_path}: every locus must be a mapping, got {entry!r}")
        name = str(entry.get("name", ""))
        kind = str(entry.get("kind", ""))
        members = tuple(str(member) for member in entry.get("members", ()))
        source = str(entry.get("source", "")).strip()
        if not name or name in seen:
            raise BarcodeLocusError(f"{source_path}: missing or repeated locus name {name!r}")
        if kind not in KINDS:
            raise BarcodeLocusError(f"{source_path}: {name}: kind must be one of {KINDS}")
        if kind == "gene" and len(members) != 1:
            raise BarcodeLocusError(f"{source_path}: {name}: a gene locus names exactly one gene")
        if kind == "intergenic" and len(members) != 2:
            raise BarcodeLocusError(
                f"{source_path}: {name}: an intergenic locus names the two flanking genes"
            )
        if not source:
            raise BarcodeLocusError(
                f"{source_path}: {name}: no source — a locus enters this list with the claim it "
                f"rests on, never as a name somebody remembers"
            )
        seen.add(name)
        loci.append(BarcodeLocus(name=name, kind=kind, members=members, source=source))
    if not loci:
        raise BarcodeLocusError(f"{source_path}: no loci")
    return tuple(loci)


@dataclass(frozen=True, slots=True)
class VariableRegion:
    """One hypervariable region: where it is, how long, and how variable."""

    start_1b: int
    end_1b: int
    n_windows: int
    mean_pi: float
    max_pi: float

    @property
    def span(self) -> int:
        return self.end_1b - self.start_1b + 1


def threshold_for(
    windows: Sequence[WindowDiversity], *, k: float = DEFAULT_THRESHOLD_SD
) -> float | None:
    """``mean + k · sd`` over the windows that have a value; ``None`` when none does.

    The standard deviation is the sample one (``statistics.stdev``), which needs at least two
    values — a series with a single usable window has no threshold and therefore no regions.
    """
    values = [entry.pi for entry in windows if entry.pi is not None]
    if len(values) < 2:
        return None
    return statistics.mean(values) + k * statistics.stdev(values)


def hypervariable_regions(
    windows: Sequence[WindowDiversity],
    *,
    threshold: float,
    min_windows: int = DEFAULT_MIN_WINDOWS,
) -> tuple[VariableRegion, ...]:
    """Maximal runs of consecutive windows above the threshold, longest-span first to last.

    Returned in **coordinate order**; a caller that wants them ranked sorts by ``mean_pi``. The
    order of ``windows`` is the order of the series (each entry follows the previous by one
    step), which is what makes "consecutive" meaningful.
    """
    regions: list[VariableRegion] = []
    run: list[WindowDiversity] = []

    def close() -> None:
        if len(run) >= min_windows:
            values = [entry.pi for entry in run if entry.pi is not None]
            regions.append(
                VariableRegion(
                    start_1b=run[0].start_1b,
                    end_1b=run[-1].end_1b,
                    n_windows=len(run),
                    mean_pi=statistics.mean(values),
                    max_pi=max(values),
                )
            )
        run.clear()

    for entry in windows:
        if entry.pi is not None and entry.pi > threshold:
            run.append(entry)
        else:
            close()
    close()
    return tuple(regions)
