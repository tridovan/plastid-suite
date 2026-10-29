"""Tier A golden numbers — pure Python, no external tools, no network.

Golden sample: Iris domestica NC_050833.1.
"""

from __future__ import annotations

import pytest

from plastid.core.ir import detect_ir, enforce_ir_gate
from plastid.core.model import Genome, coding_sequence, count_genes
from plastid.core.normalize import normalize
from plastid.core.sequence import gc_percent
from plastid.core.spans import introns_of, spacers_of
from plastid.modules import amino_acid, codon_usage, introns, structure
from plastid.run.context import RunContext


def _ctx(genome: Genome, module_id: str) -> RunContext:
    ir = detect_ir(genome)
    enforce_ir_gate(ir)
    # these modules read only params and the IR — no layout is involved
    return RunContext(layout=None, params={}, module_id=module_id, module_order=1, ir=ir)  # type: ignore[arg-type]


def test_length_and_gc(golden_genome: Genome) -> None:
    assert len(golden_genome) == 153_724
    assert gc_percent(golden_genome.sequence) == pytest.approx(37.8536, abs=0.001)


def test_gene_counts_follow_the_fixed_rule(golden_genome: Genome) -> None:
    counts = count_genes(golden_genome)
    assert (counts.cds, counts.trna, counts.rrna) == (86, 38, 8)
    assert counts.total == 132  # never report a bare "132" without the breakdown


def test_feature_gene_count_is_a_diagnostic_only(golden_genome: Genome) -> None:
    # The golden: `gene` features = 133 (one truncated ycf1 copy has no CDS);
    # unique names = 112; names repeated in the IR = 21.
    counts = count_genes(golden_genome)
    assert counts.gene_features == 133
    assert counts.unique_names == 112
    assert counts.repeated_names == 21


def test_truncated_ycf1_copy(golden_genome: Genome) -> None:
    # Measured: ycf1 has two `gene` features — a truncated pseudo copy at
    # 108,448–109,347 (gene + misc_feature, no CDS) and the real copy at
    # 122,984–128,404, which is the only one carrying a CDS.
    copies = [gene for gene in golden_genome.genes if gene.name == "ycf1"]
    gene_features = [gene for gene in copies if gene.kind == "gene"]
    cds = [gene for gene in copies if gene.kind == "CDS"]

    truncated = [gene for gene in gene_features if gene.pseudo]
    assert len(gene_features) == 2
    assert len(truncated) == 1
    assert truncated[0].span == (108_448, 109_347)  # 1-based inclusive

    assert len(cds) == 1
    assert cds[0].span == (122_984, 128_404)
    assert not any(gene.span[1] < 110_000 for gene in cds)  # nothing CDS-like in the truncated copy


def test_trans_spliced_rps12(golden_genome: Genome) -> None:
    cdss = [gene for gene in golden_genome.genes if gene.name == "rps12" and gene.kind == "CDS"]
    assert len(cdss) == 2  # one per IR copy
    assert all(len(gene.parts) == 3 for gene in cdss)  # trans-spliced exons
    assert any(gene.strand == 0 for gene in cdss)  # the IRb copy mixes strands


def test_golden_genome_identity(golden_genome: Genome) -> None:
    assert golden_genome.accession == "NC_050833.1"
    assert golden_genome.name == "Iris domestica"
    assert golden_genome.topology == "circular"


def test_ir_and_regions(golden_genome: Genome) -> None:
    ir = detect_ir(golden_genome)
    enforce_ir_gate(ir)
    assert ir.length == 26_214
    assert ir.identity == 1.0
    assert ir.source == "annotation"
    assert ir.denovo_delta == 0  # the independent de-novo scan agrees exactly

    frame = normalize(golden_genome, ir)
    assert frame.rotation_offset == 0  # the golden genome already starts at the LSC
    assert [(region.name, region.length) for region in frame.regions] == [
        ("LSC", 83_127),
        ("IRa", 26_214),
        ("SSC", 18_169),
        ("IRb", 26_214),
    ]
    assert sum(region.length for region in frame.regions) == 153_724


def test_coding_sequence_reads_each_part_on_its_own_strand(golden_genome: Genome) -> None:
    """Measured 2026-10-03: only this reading of trans-spliced rps12 yields a canonical ORF.

    Re-sorting the parts gives AAT…TAT; reading them all on the gene's strand gives ATA…ATT.
    """
    cdss = [gene for gene in golden_genome.genes if gene.kind == "CDS" and gene.name == "rps12"]
    assert len(cdss) == 2  # one per IR copy
    for gene in cdss:
        sequence = coding_sequence(golden_genome, gene)
        assert len(sequence) % 3 == 0
        assert sequence.startswith("ATG")
        assert sequence.endswith("TAA")


def test_codon_usage_covers_the_counted_cds_set(golden_genome: Genome) -> None:
    """P2-03's Done-when: the CDS set is the counted one and every family sum checks out."""
    result = codon_usage.run(golden_genome, {}, _ctx(golden_genome, "codon_usage"))

    summary = result.tables["codon_usage"].rows[0]
    assert summary["n_cds"] == count_genes(golden_genome).cds == 86
    assert summary["internal_stops"] == 0
    assert summary["translation_table"].startswith("NCBI table 1")
    assert summary["enc"] == pytest.approx(49.69, abs=0.01)
    assert summary["cai"] == pytest.approx(0.5375, abs=0.0001)
    assert summary["fop"] == pytest.approx(0.5118, abs=0.0001)

    families = result.tables["rscu_families"].rows
    assert len(families) == 20  # every amino acid, no family silently skipped
    assert [row["amino_acid"] for row in families if row["status"] != "ok"] == []


def test_amino_acid_totals_reconcile_with_codon_usage(golden_genome: Genome) -> None:
    """P2-04's Done-when: the same CDS set and the same totals as the P2-03 codon table."""
    codons = codon_usage.run(golden_genome, {}, _ctx(golden_genome, "codon_usage"))
    amino = amino_acid.run(golden_genome, {}, _ctx(golden_genome, "amino_acid"))

    codon_summary = codons.tables["codon_usage"].rows[0]
    summary = amino.tables["amino_acid"].rows[0]
    assert summary["n_cds"] == codon_summary["n_cds"] == 86
    assert summary["n_residues"] == codon_summary["n_codons"]  # every sense codon is one residue
    assert summary["ignored_symbols"] == 86  # the terminal stop codon of each CDS
    assert summary["gravy"] == pytest.approx(-0.1034, abs=0.0001)
    assert summary["aromaticity"] == pytest.approx(0.1118, abs=0.0001)

    composition = amino.tables["amino_acid_composition"].rows
    assert len(composition) == 20
    assert sum(row["count"] for row in composition) == summary["n_residues"]


def test_introns_cross_check_the_annotation_part_boundaries(golden_genome: Genome) -> None:
    """P2-05's Done-when: the intron list reconciles with the record's own part boundaries.

    M11.15 adds the tRNA introns (H-23). Their expected count is **not** ARAGORN's log line: the
    published record itself writes eight tRNAs as two-part joins, and those eight are the set
    below. Measured 26/10/2026 — CDS 17 + tRNA 8; genes with an intron 12 + 6 (two gene names
    carry one intron in each IR copy); total 20,347 bp; the longest is `trnK-UUU`'s 2,618 bp.
    """
    result = introns.run(golden_genome, {}, _ctx(golden_genome, "introns"))
    summary = result.tables["intron_summary"].rows[0]

    assert summary["n_cds_introns"] == 17
    assert summary["n_trna_introns"] == 8
    assert summary["n_introns"] == 25
    assert summary["n_genes_with_introns"] == 18
    assert summary["total_length_bp"] == 20_347
    assert summary["longest_bp"] == 2_618  # trnK-UUU's intron, which carries matK
    cds_lengths = [
        row["length_bp"] for row in result.tables["introns"].rows if row["kind"] == "CDS"
    ]
    assert max(cds_lengths) == 1_055
    # the two rps12 junctions cross a region boundary: trans-splicing, not introns
    assert summary["n_trans_splicing"] == 2
    assert summary["assessment"] == "CDS and tRNA introns measured"
    assert summary["annotation_check"].startswith("agrees")

    rows = result.tables["introns"].rows
    assert len(rows) == 25
    assert len({row["intron_id"] for row in rows}) == 25  # both IR copies named apart
    assert {row["kind"] for row in rows} == {"CDS", "tRNA"}
    assert all(0 <= row["phase"] <= 2 for row in rows)
    ndh_a = next(row for row in rows if row["gene"] == "ndhA")
    assert (ndh_a["start_1b"], ndh_a["end_1b"]) == (119_328, 120_382)  # from the CDS parts
    assert ndh_a["phase"] == 1

    # The eight tRNA introns are the record's own multi-part tRNAs, name for name, and their
    # lengths are the gaps between the parts the record itself writes — an independent count.
    published = {
        gene.name for gene in golden_genome.genes if gene.kind == "tRNA" and len(gene.parts) > 1
    }
    assert {row["gene"] for row in rows if row["kind"] == "tRNA"} == published
    trn_k = next(row for row in rows if row["intron_id"] == "trnK-UUU.1.1")
    assert (trn_k["start_1b"], trn_k["end_1b"]) == (1_403, 4_020)
    assert trn_k["length_bp"] == 2_618  # the record's matK-bearing intron


def test_the_shared_span_derivation_reproduces_the_introns_module(golden_genome: Genome) -> None:
    """P2B-01's cross-check: one derivation, two modules, the same spans.

    The spans moved to ``core/spans.py`` because ``diversity`` measures π over them; the check
    the board asks for is that the move changed nothing — every id, coordinate and phase the
    ``introns`` table publishes is the one the derivation produces. M11.15 adds the tRNA set as a
    second derivation: the CDS set is still the one ``diversity`` reads, and the published table
    is the union with each row's ``kind`` stated.
    """
    ir = detect_ir(golden_genome)
    derived = introns_of(golden_genome, ir)
    result = introns.run(golden_genome, {}, _ctx(golden_genome, "introns"))
    published = {
        row["intron_id"]: (
            row["kind"],
            row["gene"],
            row["start_1b"],
            row["end_1b"],
            row["length_bp"],
            row["phase"],
        )
        for row in result.tables["introns"].rows
    }
    assert {
        span.intron_id: (span.kind, span.gene, span.start + 1, span.end, span.length, span.phase)
        for span in (*derived.introns, *derived.trna_introns)
    } == published
    assert [span.junction_id for span in derived.junctions] == ["rps12.1.1", "rps12.2.1"]


def test_the_golden_spacers_are_measured_not_inherited(golden_genome: Genome) -> None:
    """P2B-01's own numbers — measured 10/10/2026, never a survey's count.

    The phase's trap is the 26/9 survey's numbers measured on *other* data. These are the
    golden's: 91 spacers in the single-copy regions, 38 more inside the IR (excluded — its two
    copies are one sequence annotated twice), 2 crossing a region boundary, and 4 gene pairs
    that leave no gap at all (psbD·psbC −53, ndhK·ndhC −10, atpE·atpB −4, and rps12's two CDS
    features, which share their LSC exon).
    """
    spacers = spacers_of(golden_genome, detect_ir(golden_genome))

    assert len(spacers.igs) == 91
    assert [span.region for span in spacers.igs].count("LSC") == 79
    assert [span.region for span in spacers.igs].count("SSC") == 12
    assert sum(span.length for span in spacers.igs) == 32_134
    assert len(spacers.in_ir) == 38
    assert [(gap.left, gap.right) for gap in spacers.crossing] == [
        ("rpl22", "rps19"),
        ("rps19", "psbA"),
    ]
    assert [(gap.left, gap.right, gap.length) for gap in spacers.no_spacer] == [
        ("psbD", "psbC", -53),
        ("ndhK", "ndhC", -10),
        ("atpE", "atpB", -4),
        ("rps12", "rps12", -114),
    ]
    # rps12's exon in the LSC is flanked by real spacers: a gene occupies its exons, not the
    # space between them, which is what the outer span would have claimed (P2B-01)
    assert ("rps12-clpP", 138) in {(span.igs_id, span.length) for span in spacers.igs}
    # and matK, annotated inside trnK's intron, is flanked by two of them
    assert [span.length for span in spacers.igs if "matK" in span.igs_id] == [260, 789]


def test_structure_module_tiles_the_golden_genome(golden_genome: Genome) -> None:
    """P2-02's Done-when: LSC 83,127 · SSC 18,169 · IRa = IRb = 26,214 · exact tiling."""
    result = structure.run(golden_genome, {}, _ctx(golden_genome, "structure"))

    rows = {row["region"]: row for row in result.tables["regions"].rows}
    assert rows["LSC"]["length_bp"] == 83_127
    assert rows["SSC"]["length_bp"] == 18_169
    assert rows["IRa"]["length_bp"] == rows["IRb"]["length_bp"] == 26_214
    assert sum(row["length_bp"] for row in result.tables["regions"].rows) == 153_724
    assert result.tables["structure_summary"].rows[0]["n_windows"] == 157


def test_orf_checks_pass_on_the_golden_genome(golden_genome: Genome) -> None:
    """P3-01's Done-when: 0 internal stops, 0 missing stops, every CDS divisible by three, and
    exactly 5 classified ``alternative_start`` — rps19 GTG (both IR copies), rpl2 ACG (both),
    ndhD ATC. Measured 2026-10-03; the counts are the board's, not adjusted to the code."""
    from plastid.modules import orf

    result = orf.run(golden_genome, {}, _ctx(golden_genome, "orf"))
    summary = result.tables["orf_summary"].rows[0]
    assert summary["n_cds"] == count_genes(golden_genome).cds == 86
    assert summary["n_failed"] == 0
    assert summary["n_pseudogene"] == 0
    assert summary["translation_tables"] == "11"

    status = {row["status"]: row["n"] for row in result.tables["orf_status"].rows}
    assert status["ok"] == 81
    assert status["alternative_start"] == 5
    assert status["no_start"] == status["no_stop"] == 0
    assert status["internal_stop"] == status["not_multiple_of_three"] == 0

    alternative = [
        (row["gene"], row["start_codon"])
        for row in result.tables["orf_checks"].rows
        if row["status"] == "alternative_start"
    ]
    assert sorted(alternative) == [
        ("ndhD", "ATC"),
        ("rpl2", "ACG"),
        ("rpl2", "ACG"),
        ("rps19", "GTG"),
        ("rps19", "GTG"),
    ]
    # ACG is not in table 11's initiator list: flagged, but as non-canonical, not as a failure
    non_canonical = [
        row["gene"]
        for row in result.tables["orf_checks"].rows
        if row["start_class"] == "non_canonical"
    ]
    assert sorted(non_canonical) == ["rpl2", "rpl2"]


def test_orf_check_reports_the_truncated_ycf1_as_annotated_not_coding(
    golden_genome: Genome,
) -> None:
    """Measured: ycf1's truncated copy is a ``gene`` + ``misc_feature`` with no CDS.
    It is reported in its own table — not as a failed CDS, and never "repaired"."""
    from plastid.modules import orf

    result = orf.run(golden_genome, {}, _ctx(golden_genome, "orf"))
    summary = result.tables["orf_summary"].rows[0]
    assert summary["n_annotated_not_coding"] == 1

    rows = result.tables["annotated_not_coding"].rows
    assert [(row["name"], row["start_1b"], row["end_1b"]) for row in rows] == [
        ("ycf1", 108_448, 109_347)
    ]
    assert rows[0]["pseudo"] == "yes"


def test_gene_table_reconciles_with_the_genome_wide_count(golden_genome: Genome) -> None:
    """P2-09: the publication table's own totals must equal `count_genes()` — 132 = 86 + 38 + 8.

    The reconciliation is the module's own guard: a gene the grouping file failed to claim, or
    a row silently dropped between the counter and the table, shows up here as a disagreement
    rather than as a smaller number nobody checks.
    """
    from plastid.modules import gene_table

    result = gene_table.run(golden_genome, {}, _ctx(golden_genome, "gene_table"))

    counts = count_genes(golden_genome)
    rows = {row["kind"]: row for row in result.tables["gene_reconciliation"].rows}
    assert rows["CDS"]["in_table"] == counts.cds == 86
    assert rows["tRNA"]["in_table"] == counts.trna == 38
    assert rows["rRNA"]["in_table"] == counts.rrna == 8
    assert rows["Total"]["in_table"] == counts.total == 132
    assert [row["agreement"] for row in rows.values()] == ["agrees"] * 4

    # the printed table ends on the same total — the two cannot drift apart
    groups = result.tables["gene_groups"]
    assert groups.rows[-1]["group"] == "Total"
    assert groups.rows[-1]["copies"] == sum(row["copies"] for row in groups.rows[:-1]) == 132

    summary = result.tables["gene_table_summary"].rows[0]
    assert summary["n_genes_total"] == 132
    assert summary["n_unclassified"] == 0
    assert summary["n_unique_names"] == 112  # 132 copies, 112 distinct names


def test_gene_table_lists_tirna_by_anticodon_with_copy_counts(golden_genome: Genome) -> None:
    """The tRNA row is built from the record, not from a list that only fits one species."""
    from plastid.modules import gene_table

    result = gene_table.run(golden_genome, {}, _ctx(golden_genome, "gene_table"))
    rows = {row["group"]: row for row in result.tables["gene_groups"].rows}

    trna = rows["Transfer RNA genes"]
    assert trna["copies"] == 38
    assert "trnA-UGC (\N{MULTIPLICATION SIGN}2)" in trna["genes"]  # both IR copies
    assert "trnfM-CAU" in trna["genes"]  # lowercase f: a prefix match, not a name list
    assert "trnL-TAA" not in trna["genes"]

    # ycf3/ycf4 are printed under their canonical names — the naming reference's spelling
    # (D5-09, flipped 08/10/2026) — and sit with
    # photosystem I, so only ycf1 and the two ycf2 copies remain in the unknown-function group
    assert rows["Conserved open reading frames"]["copies"] == 3  # ycf1, ycf2 x2
    assert "ycf3 (pafI)" in rows["Subunits of photosystem I"]["genes"]
    assert "ycf4 (pafII)" in rows["Subunits of photosystem I"]["genes"]
    assert "psbN (pbf1)" in rows["Subunits of photosystem II"]["genes"]
    assert rows["Protease"]["genes"] == "clpP (clpP1)"
    assert rows["Subunits of NADH dehydrogenase"]["copies"] == 12  # includes the ndhL fragment


def test_gene_table_flags_a_gene_no_group_claims() -> None:
    """Negative case: an unknown name is reported as unclassified, never folded in silently."""
    from plastid.core.model import Gene as Feature
    from plastid.core.model import Interval
    from plastid.modules import gene_table

    genome = Genome(
        name="synthetic",
        sequence="ATGC" * 400,
        genes=(
            Feature("rbcL", "CDS", (Interval(0, 120),)),
            Feature("mystery9", "CDS", (Interval(200, 300),)),
        ),
    )
    # gene_table needs no IR: grouping a gene by function is not a position-based analysis
    ctx = RunContext(layout=None, params={}, module_id="gene_table", module_order=1)  # type: ignore[arg-type]
    result = gene_table.run(genome, {}, ctx)

    rows = {row["kind"]: row for row in result.tables["gene_reconciliation"].rows}
    assert rows["CDS"]["in_table"] == rows["CDS"]["counted"] == 2  # still counted, not hidden
    assert rows["CDS"]["agreement"] == "agrees"
    assert result.tables["gene_table_summary"].rows[0]["n_unclassified"] == 1
    assert any("mystery9" in note for note in result.notes)


def test_a_gene_no_group_claims_gets_its_own_row_and_counts_copies() -> None:
    """`gene_groups.yaml` promises an "Unclassified" row; the table must carry it.

    Measured 29/10/2026 on `NC_016068.1` (Solanaceae): 24 copies of `hypothetical protein` had no
    group, the notes named the label once, and the summary read "Unclassified: 1" — a **name**
    count in a summary whose every other cell counts copies. A reader comparing it with
    `n_genes_total` cannot make the two add up.
    """
    from plastid.core.model import Gene as Feature
    from plastid.core.model import Genome, Interval
    from plastid.modules import gene_table

    genome = Genome(
        name="synthetic",
        sequence="ATGC" * 400,
        genes=(
            Feature("rbcL", "CDS", (Interval(0, 120),)),
            Feature("mystery9", "CDS", (Interval(200, 300),)),
            Feature("mystery9", "CDS", (Interval(400, 500),)),
            Feature("mystery9", "CDS", (Interval(600, 700),)),
        ),
    )
    ctx = RunContext(layout=None, params={}, module_id="gene_table", module_order=1)  # type: ignore[arg-type]
    result = gene_table.run(genome, {}, ctx)

    rows = {row["group"]: row for row in result.tables["gene_groups"].rows}
    assert "Unclassified" in rows, "the catalogue's promised row is missing"
    assert rows["Unclassified"]["copies"] == 3
    assert "mystery9" in rows["Unclassified"]["genes"]

    summary = result.tables["gene_table_summary"].rows[0]
    assert summary["n_unclassified"] == 3  # copies, like every other count in the summary
    assert rows["Total"]["copies"] == summary["n_genes_total"] == 4
    assert sum(row["copies"] for name, row in rows.items() if name != "Total") == 4
