"""The mitochondrial path (P3-04): no IR gate, and a name that is really a product sentence.

The annotator this path calls is the project's own reference-transfer caller (P3A);
PMGA, the tool P3-04 first chose, was dropped — it cannot run here. The two things the
path must get right do not depend on the caller: a plant mitogenome must be describable
without the inverted-repeat gate, and a gene must be called by a symbol rather than by the
sentence its record happens to carry.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from plastid.core.genbank import read_genbank_records
from plastid.core.ir import IRDetectionError, IRGateError
from plastid.core.model import Gene, Genome, Interval
from plastid.modules import mito_genes
from plastid.run.context import RunContext

PANEL_MITO = Path(__file__).resolve().parents[2] / ".panel" / "oldenlandia-coleus-mito"
NEEDS_PANEL = pytest.mark.skipif(
    not PANEL_MITO.is_dir(),
    reason=(
        "the pinned mitogenome is not on this machine — re-fetch it with the 'plastid panel pin' "
        "command"
    ),
)


def _ctx(ir=None) -> RunContext:
    return RunContext(layout=None, params={}, module_id="mito_genes", module_order=1, ir=ir)  # type: ignore[arg-type]


def _mitogenome(
    *genes: Gene, name: str = "synthetic mito", accession: str = "SYN_MITO.1"
) -> Genome:
    return Genome(
        name=name,
        accession=accession,
        sequence="ATGC" * 500,
        genes=genes,
        topology="circular",
    )


def test_the_module_runs_without_an_ir_result_at_all() -> None:
    """The gate must not be asked for here: a plant mitogenome has no canonical IR, and the two
    stand-in mitogenomes raise IRDetectionError — correctly. A module that gated would refuse to
    describe a genome it can read."""
    product = "NADH dehydrogenase subunit 1"
    genome = _mitogenome(Gene(product, "CDS", (Interval(0, 300),), product=product))
    result = mito_genes.run(genome, {}, _ctx(ir=None))  # no IR, and none is asked for

    assert result.tables["mito_summary"].rows[0]["ir_gate"] == "not applicable"
    assert any("never asks for one" in note or "no canonical IR" in note for note in result.notes)


def test_a_gated_module_on_the_same_context_still_stops() -> None:
    """The contrast that gives the previous test its meaning: the gate is real, it is just in
    front of position-based analysis and not in front of counting genes."""
    from plastid.modules import genome_stats

    genome = _mitogenome(Gene("rbcL", "CDS", (Interval(0, 300),), product="rubisco"))
    ctx = RunContext(layout=None, params={}, module_id="genome_stats", module_order=1)  # type: ignore[arg-type]
    with pytest.raises(IRGateError, match="never determined one|IR gate failed"):
        genome_stats.run(genome, {}, ctx)


def test_a_product_sentence_becomes_the_standard_symbol() -> None:
    """Measured: OX459128.1 carries no /gene qualifier, so the name our reader reports is the
    record's product sentence. A gene table that calls a gene by a sentence is not reporting a
    gene name."""
    genome = _mitogenome(
        Gene(
            "NADH dehydrogenase subunit 1",
            "CDS",
            (Interval(0, 300),),
            product="NADH dehydrogenase subunit 1",
        ),
    )
    result = mito_genes.run(genome, {}, _ctx())
    row = result.tables["gene_set"].rows[0]
    assert row["gene"] == "nad1"  # the symbol; the product sentence has its own column
    assert row["name_as_arrived"] == "NADH dehydrogenase subunit 1"  # never rewritten
    assert row["product"] == "NADH dehydrogenase subunit 1"
    assert result.tables["mito_summary"].rows[0]["n_named_by_symbol"] == 1
    assert "unmapped_names" not in result.tables


def test_an_unknown_product_sentence_is_reported_not_renamed() -> None:
    genome = _mitogenome(
        Gene("hypothetical protein X", "CDS", (Interval(0, 300),), product="hypothetical protein X")
    )
    result = mito_genes.run(genome, {}, _ctx())

    assert result.tables["gene_set"].rows[0]["gene"] == "hypothetical protein X"
    assert [row["name"] for row in result.tables["unmapped_names"].rows] == [
        "hypothetical protein X"
    ]
    assert any("does not invent" in note for note in result.notes)


def test_a_record_that_already_gives_symbols_is_not_second_guessed() -> None:
    """A proper mitogenome annotation names its genes; those names are left alone and are not
    listed as unmapped just because the symbol table did not have to do anything."""
    genome = _mitogenome(
        Gene("nad1", "CDS", (Interval(0, 300),), product="NADH dehydrogenase subunit 1"),
        Gene("trnA-UGC", "tRNA", (Interval(400, 470),), product="tRNA-Ala"),
    )
    result = mito_genes.run(genome, {}, _ctx())
    assert {row["gene"] for row in result.tables["gene_set"].rows} == {"nad1", "trnA-UGC"}
    assert "unmapped_names" not in result.tables


def test_the_mitochondrial_table_never_renames_a_plastid_gene() -> None:
    """Measured collision: "ribosomal protein S3" is a product string in both organelles. The
    plastid normalizer must leave it alone — a shared table would silently rewrite a plastid gene
    with a mitochondrial symbol."""
    from plastid.core.gene_names import load_mito_aliases, normalize_name

    assert normalize_name("ribosomal protein S3").canonical == "ribosomal protein S3"
    assert normalize_name("ribosomal protein S3", load_mito_aliases()).canonical == "rps3"


@NEEDS_PANEL
def test_the_stand_in_mitogenome_is_reported_with_symbols_and_without_the_gate() -> None:
    """The measured end of P3-04's testable half: 51 genes = 31 CDS + 17 tRNA + 3 rRNA, every
    CDS named by a standard symbol, no IR result anywhere in the run."""
    genome = read_genbank_records(PANEL_MITO / "OX459128.1.gb")[0]
    with pytest.raises(IRDetectionError):
        from plastid.core.ir import detect_ir

        detect_ir(genome)  # the gate could not have passed: there is no IR to find

    result = mito_genes.run(genome, {}, _ctx(ir=None))
    summary = result.tables["mito_summary"].rows[0]
    assert (summary["n_genes"], summary["n_cds"], summary["n_trna"], summary["n_rrna"]) == (
        51,
        31,
        17,
        3,
    )
    assert summary["n_named_by_symbol"] == 31
    assert "unmapped_names" not in result.tables
    genes = {row["gene"].split(" ")[0] for row in result.tables["gene_set"].rows}
    assert {"nad1", "cox1", "atp6", "ccmFn", "matR", "rrnL"} <= genes
