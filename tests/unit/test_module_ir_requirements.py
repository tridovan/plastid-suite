"""M11.5 — a module depends only on the data it actually reads.

Codon usage and amino-acid composition count CDS sequences; no inverted-repeat position enters a
number they print. They were nevertheless gated on the IR, so a record with readable CDS and a
divergent repeat produced nothing. The tests below hold the two halves of the rule: the sequenced
modules run with **no IR result at all**, and a position-based module still refuses to.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from plastid.core.codon import UnsupportedTranslationError
from plastid.core.ir import IRGateError
from plastid.core.model import Gene, Genome, Interval
from plastid.modules import amino_acid, codon_usage, genome_stats, registry
from plastid.run.context import RunContext
from plastid.run.layout import RunLayout


def _genome(*, codon_start: int = 1, transl_table: int | None = 11) -> Genome:
    sequence = "ATG" + "GCT" * 10 + "TAA"
    gene = Gene(
        name="psbA",
        kind="CDS",
        parts=(Interval(0, len(sequence)),),
        transl_table=transl_table,
        codon_start=codon_start,
    )
    return Genome(name="synthetic", accession="SYN.1", sequence=sequence, genes=(gene,))


def _ctx(tmp_path: Path) -> RunContext:
    """A context with **no IR result**, the shape a run of these modules has."""
    layout = RunLayout.create(tmp_path / "runs", "ir-requirements", fingerprint="test")
    layout.make_dirs()
    return RunContext(layout=layout, params={}, module_id="test", module_order=1, ir=None)


def test_the_registry_says_which_modules_need_the_inverted_repeat() -> None:
    assert registry.get("codon_usage").requires_ir is False
    assert registry.get("amino_acid").requires_ir is False
    assert registry.get("genome_stats").requires_ir is True


def test_codon_usage_runs_without_an_ir_result(tmp_path: Path) -> None:
    result = codon_usage.run(_genome(), {}, _ctx(tmp_path))

    assert result.tables["codon_usage"].rows[0]["n_cds"] == 1
    assert result.tables["codon_usage"].rows[0]["n_codons"] > 0


def test_amino_acid_runs_without_an_ir_result(tmp_path: Path) -> None:
    result = amino_acid.run(_genome(), {}, _ctx(tmp_path))

    assert result.tables["amino_acid"].rows[0]["n_cds"] == 1
    assert result.tables["amino_acid"].rows[0]["n_residues"] > 0


def test_a_position_based_module_is_still_blocked_without_the_gate(tmp_path: Path) -> None:
    with pytest.raises(IRGateError, match="needs an IR result"):
        genome_stats.run(_genome(), {}, _ctx(tmp_path))


@pytest.mark.parametrize("module", [codon_usage, amino_acid])
def test_a_frame_this_code_cannot_read_is_refused_not_translated(module, tmp_path: Path) -> None:
    """``/codon_start=2`` shifts every codon; reading from base one would be a different number."""
    with pytest.raises(UnsupportedTranslationError, match=r"codon_start=2"):
        module.run(_genome(codon_start=2), {}, _ctx(tmp_path))


@pytest.mark.parametrize("module", [codon_usage, amino_acid])
def test_another_translation_table_is_refused_by_name(module, tmp_path: Path) -> None:
    """The compiled code is table 1; table 11 shares its amino acids, table 4 does not."""
    with pytest.raises(UnsupportedTranslationError, match="table"):
        module.run(_genome(transl_table=4), {}, _ctx(tmp_path))


@pytest.mark.parametrize("module", [codon_usage, amino_acid])
def test_table_eleven_is_accepted(module, tmp_path: Path) -> None:
    """The plastid code shares every amino-acid assignment with table 1, so it is supported."""
    assert module.run(_genome(transl_table=11), {}, _ctx(tmp_path)) is not None
