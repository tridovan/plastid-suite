"""What decides the assembler's memory, on every run deposited for Jin et al. 2020 — panel tier.

`.panel/assembly/jin2020/plant.plastome.summary.corrected.tab` is the authors' own table of 308
runs: the 50 samples under seven settings, with each stage's peak memory as their runs measured it.
No tool runs here. These tests hold what `docs/tools/getorganelle.md` says about the memory, and the
reason the step limits a number of words instead of estimating a number of gigabytes.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pytest

from plastid.external.getorganelle import ASSEMBLER_WORD_LIMIT
from plastid.modules.assembly import resolve_params, word_limit

pytestmark = pytest.mark.panel

TABLE = (
    Path(__file__).resolve().parents[2]
    / ".panel"
    / "assembly"
    / "jin2020"
    / "plant.plastome.summary.corrected.tab"
)
_PEAK = "mem_max"


@dataclass(frozen=True)
class Run:
    sample: str
    words: float
    word_size: float
    read_length: float
    extending_gb: float
    peak_gb: float
    plastid_share: float
    bases: float


@pytest.fixture(scope="module")
def runs() -> list[Run]:
    if not TABLE.exists():
        pytest.skip(f"not in this tree (.panel/ is gitignored): {TABLE}")
    read = []
    with TABLE.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            try:
                words = float(row["accepted_words"])
                extending = float(row["mem_extending"]) / 1e9
                word_size = float(row["w"])
                read_length = float(row["mean_read_len"])
            except ValueError:  # a run that stopped before that stage has no number there
                continue
            if words <= 0:
                continue
            length = sum(float(part) for part in row["res_len"].split("+"))
            bases = float(row["n_bases"])
            read.append(
                Run(
                    sample=row["sample"],
                    words=words,
                    word_size=word_size,
                    read_length=read_length,
                    extending_gb=extending,
                    peak_gb=float(row[_PEAK]) / 1e9,
                    plastid_share=float(row["res_base_cov"]) * length / bases,
                    bases=bases,
                )
            )
    return read


def _r_squared(design: np.ndarray, target: np.ndarray) -> float:
    fitted = design @ np.linalg.lstsq(design, target, rcond=None)[0]
    return 1 - float(((fitted - target) ** 2).sum() / ((target - target.mean()) ** 2).sum())


def test_every_deposited_run_with_an_extension_is_read(runs: list[Run]) -> None:
    assert len(runs) == 308


def test_the_words_held_decide_the_memory_of_read_recruitment(runs: list[Run]) -> None:
    words = np.array([run.words for run in runs])
    word_size = np.array([run.word_size for run in runs])
    memory = np.array([run.extending_gb for run in runs])
    design = np.column_stack((np.ones(len(words)), words, words * word_size))
    assert _r_squared(design, memory) == pytest.approx(0.997, abs=0.001)
    offset, per_word, per_base = np.linalg.lstsq(design, memory, rcond=None)[0]
    assert offset == pytest.approx(0.38, abs=0.01)
    params = resolve_params({})
    assert per_word * 1e9 == pytest.approx(params["word_bytes"], abs=1)
    assert per_base * 1e9 == pytest.approx(params["word_bytes_per_base"], abs=0.01)


def test_the_plastid_share_does_not_predict_the_peak(runs: list[Run]) -> None:
    """Why the step limits and does not estimate: the 50 runs of the published settings."""
    published = [run for run in runs if run.sample.startswith("auto--")]
    assert len(published) == 50
    share = np.array([math.log10(run.plastid_share) for run in published])
    bases = np.array([math.log10(run.bases) for run in published])
    peak = np.array([math.log10(run.peak_gb) for run in published])
    design = np.column_stack((np.ones(len(share)), share, bases))
    assert _r_squared(design, peak) == pytest.approx(0.49, abs=0.01)


@pytest.mark.parametrize("memory_gb", [4, 7, 8, 16, 32, 64])
def test_no_run_under_its_word_limit_went_over_the_budget(runs: list[Run], memory_gb: int) -> None:
    params = resolve_params({})
    under = [
        run for run in runs if run.words <= word_limit(memory_gb, int(run.read_length), params)
    ]
    assert len(under) >= 250
    assert max(run.extending_gb for run in under) <= memory_gb


def test_the_assemblers_own_limit_is_the_worst_peak_there_is(runs: list[Run]) -> None:
    # Three runs reached its 4E8 words and took 61 to 65 GB, so it limits nothing on these machines.
    reached = [run for run in runs if run.words >= ASSEMBLER_WORD_LIMIT]
    assert len(reached) >= 3
    assert all(60 < run.extending_gb < 66 for run in reached)
