"""The same sequence under two annotations.

The 30 `annotation-variants/` records are PGA's annotation of the **same sequence** as
`source/` (the panel README states it; the first test checks it). That gives an independent
expectation for every module that reads only the sequence: its result must not move when the
annotation does. A module whose output *does* move is annotation-dependent by definition, and the
case says so rather than calling the difference an error.

The file is in two parts:

* `genome_stats` (length/GC), `ssr` (MISA) and `repeats` (repeat-match) read the sequence only;
  every cell is asserted invariant. `ssr`/`repeats` need `misa.pl`/`repeat-match`; the tests skip
  with the tool's name when absent, as every tier-B check does.
* the gene-based modules each state their class and check the claim that follows from it:
  `substitutions` (whole scope) and the whole-genome tables of `diversity`/`barcodes` are
  **annotation-independent** and are asserted identical; `gene_table`, `orf`, `codon_usage`,
  `amino_acid`, `stats`, `kaks` and `popgen` are **annotation-dependent**, and the case shows the
  change is the annotation's — the unit set, the names and the coordinates — while a measurement
  of a locus both annotations place identically does not move.

`nc_021432.1` (Orchidaceae) has **no** annotation-variants record. It is not measured here and is
not counted as a pass; the evidence record names it.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from conftest import iqtree_directory, require_tools
from plastid.core.codon import codon_counts, translate
from plastid.core.genbank import read_genbank
from plastid.core.gene_names import (
    cds_by_canonical_name,
    load_gene_aliases,
    normalize_name,
)
from plastid.core.ir import detect_ir, enforce_ir_gate
from plastid.core.model import Genome, cds_features, coding_sequence, count_genes
from plastid.core.protein import composition, ignored_symbols
from plastid.modules import (
    amino_acid,
    barcodes,
    codon_usage,
    diversity,
    gene_table,
    genome_stats,
    kaks,
    orf,
    phylogeny,
    popgen,
    repeats,
    ssr,
    stats,
    structure,
    submission,
    substitutions,
    synteny,
)
from plastid.run.context import RunContext
from plastid.run.layout import RunLayout

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / ".panel" / "ir-validation" / "source"
VARIANTS = ROOT / ".panel" / "ir-validation" / "annotation-variants"
ACCESSIONS = ("NC_008641.1", "NC_016708.2", "NC_028027.1")

# The nine families whose records have a PGA annotation. Orchidaceae is absent from the panel's
# variants and is reported as not measured, not as a pass.
FAMILY_VARIANTS = {
    "Malvaceae": "NC_008641.1",
    "Asteraceae": "NC_028027.1",
    "Brassicaceae": "NC_029378.1",
    "Euphorbiaceae": "NC_016736.1",
    "Fabaceae": "NC_016708.2",
    "Lamiaceae": "NC_020098.1",
    "Poaceae": "NC_023245.1",
    "Rosaceae": "NC_031163.1",
    "Solanaceae": "NC_016068.1",
}

# `kaks` and `popgen` align every shared CDS, and the whole-plastome modules align 160 kb per
# call: on this tree one accession's pair costs 30–90 s. Each of the expensive modules therefore
# runs its case on one family — Malvaceae, whose two annotations share the most identically
# placed genes (64 of 77) — and the evidence record says so. The other eight families' variants
# are on the panel and are not measured in this case.
CODON_MODEL_FAMILIES = ("Malvaceae",)
WHOLE_GENOME_FAMILIES = ("Malvaceae",)

# The fixed genomes an arity-2 module is compared against. They never change between the two runs,
# so any difference in the result is the swapped annotation's.
PARTNER_COUNT = 2

# The whole-plastome alignment modules: MAFFT's `auto` picks L-INS-i for two divergent records and
# spends 70–90 s per call, which is time the claim does not need. `fftns` keeps each case at a few
# seconds; both annotations get the same strategy, and the claim is about the annotation, not the
# aligner.
FAST_ALIGNMENT = {"mafft_strategy": "fftns"}

pytestmark = [
    pytest.mark.tier_b,
    pytest.mark.panel,
    pytest.mark.skipif(
        not SOURCE.is_dir() or not VARIANTS.is_dir(),
        reason="the IR validation panels are not in this tree (.panel/ is gitignored)",
    ),
]

_ALIASES = load_gene_aliases()


def _ctx(tmp_path: Path, module: str) -> RunContext:
    layout = RunLayout.create(tmp_path, module, fingerprint="panel")
    layout.make_dirs()
    return RunContext(layout=layout, params={}, module_id=module, module_order=1)


def _gate(genome) -> None:
    return enforce_ir_gate(detect_ir(genome, source="sequence"))


def _pair(accession: str) -> tuple[Genome, Genome]:
    """The source record and PGA's annotation of the same sequence."""
    return (
        read_genbank(SOURCE / f"{accession}.gb"),
        read_genbank(VARIANTS / f"{accession}.gb"),
    )


def _partners(accession: str, count: int = PARTNER_COUNT) -> list[Genome]:
    return [
        read_genbank(SOURCE / f"{other}.gb")
        for other in FAMILY_VARIANTS.values()
        if other != accession
    ][:count]


def _same_parts(source: Genome, variant: Genome) -> set[str]:
    """Canonical CDS names both annotations place at exactly the same location."""
    left, right = cds_by_canonical_name(source), cds_by_canonical_name(variant)
    return {name for name in set(left) & set(right) if left[name].parts == right[name].parts}


def _stats_units(genome: Genome) -> tuple[set[str], int]:
    """`stats`' unit: one observation per canonical CDS locus, IR copies collapsed."""
    seen: set[str] = set()
    collapsed = 0
    for gene in genome.genes:
        if gene.kind != "CDS" or gene.pseudo:
            continue
        locus = normalize_name(gene.name, _ALIASES).canonical
        if locus in seen:
            collapsed += 1
            continue
        seen.add(locus)
    return seen, collapsed


def test_the_pair_is_one_sequence_with_two_annotations() -> None:
    """The premise of every case below: only the annotation differs."""
    for accession in ACCESSIONS:
        source = read_genbank(SOURCE / f"{accession}.gb")
        variant = read_genbank(VARIANTS / f"{accession}.gb")
        assert source.sequence == variant.sequence, f"{accession}: the sequence differs"
        assert len(source) == len(variant)


def test_genome_stats_length_and_gc_do_not_depend_on_the_annotation(tmp_path: Path) -> None:
    """Length and GC are properties of the sequence; an annotation change cannot move them."""
    for accession in ACCESSIONS:
        values = []
        for folder in (SOURCE, VARIANTS):
            genome = read_genbank(folder / f"{accession}.gb")
            ctx = _ctx(tmp_path / folder.name / accession, "genome_stats")
            ctx.ir = _gate(genome)
            row = genome_stats.run(genome, {}, ctx).tables["genome_stats"].rows[0]
            values.append((row["length_bp"], row["gc_pct"]))
        assert values[0] == values[1], f"{accession}: {values}"


def test_ssr_totals_do_not_depend_on_the_annotation(
    tools: dict[str, str | None], tmp_path: Path
) -> None:
    """MISA reads the sequence; the IR is the sequence's own, so every SSR row must be identical."""
    require_tools(tools, "perl", "misa.pl")
    for accession in ACCESSIONS:
        tables = []
        for folder in (SOURCE, VARIANTS):
            genome = read_genbank(folder / f"{accession}.gb")
            ctx = _ctx(tmp_path / folder.name / accession, "ssr")
            ctx.ir = _gate(genome)
            tables.append(ssr.run(genome, {"misa": tools["misa.pl"]}, ctx).tables)
        source_summary = tables[0]["ssr_summary"].rows[0]
        variant_summary = tables[1]["ssr_summary"].rows[0]
        for key in ("total", "mono", "di", "tri", "compound", "thresholds"):
            assert source_summary[key] == variant_summary[key], f"{accession}.{key}"
        rows = [
            [
                (row["motif"], row["unit"], row["start_1b"], row["end_1b"], row["compound"])
                for row in table["ssrs"].rows
            ]
            for table in tables
        ]
        assert rows[0] == rows[1], f"{accession}: SSR positions moved with the annotation"


def test_exact_repeats_do_not_depend_on_the_annotation(
    mummer_tools: dict[str, str | None], tmp_path: Path
) -> None:
    """repeat-match reads the sequence; the same sequence gives the same pairs."""
    require_tools(mummer_tools, "repeat-match")
    directory = str(Path(mummer_tools["repeat-match"]).parent)
    for accession in ACCESSIONS:
        tables = []
        for folder in (SOURCE, VARIANTS):
            genome = read_genbank(folder / f"{accession}.gb")
            ctx = _ctx(tmp_path / folder.name / accession, "repeats")
            ctx.ir = _gate(genome)
            tables.append(repeats.run(genome, {"repeats": directory}, ctx).tables)
        source_rows = [
            (
                row["class"],
                row["length_bp"],
                row["start1_1b"],
                row["end1_1b"],
                row["start2_1b"],
                row["end2_1b"],
            )
            for row in tables[0]["repeats"].rows
        ]
        variant_rows = [
            (
                row["class"],
                row["length_bp"],
                row["start1_1b"],
                row["end1_1b"],
                row["start2_1b"],
                row["end2_1b"],
            )
            for row in tables[1]["repeats"].rows
        ]
        assert source_rows == variant_rows, f"{accession}: exact repeats moved with the annotation"


# ---------------------------------------------------------------- gene-based modules
#
# Each case below runs one module twice on the same sequence under its two annotations and states
# the class the module belongs to. An annotation-independent part is asserted identical cell for
# cell; an annotation-dependent part is shown to be the annotation's — the unit set, the names,
# the coordinates — and a measurement of a locus both annotations place identically must not move.


@pytest.mark.parametrize("family", tuple(FAMILY_VARIANTS))
def test_gene_table_counts_the_features_each_annotation_declares(
    tmp_path: Path, family: str
) -> None:
    """Annotation-dependent: every cell is a count of the annotation's own features.

    The claim is the accounting on **both** annotations — the total is `count_genes()`, the group
    copies add up to it, and every kind agrees with the count. The difference between the two
    tables is the annotation's own (PGA does not annotate rRNA), which is what a gene table is
    meant to show; it is not asserted to be zero.
    """
    accession = FAMILY_VARIANTS[family]
    for folder, genome in zip(("source", "variant"), _pair(accession), strict=True):
        result = gene_table.run(genome, {}, _ctx(tmp_path / folder, "gene_table"))
        summary = result.tables["gene_table_summary"].rows[0]
        counted = count_genes(genome)
        expected = counted.cds + counted.trna + counted.rrna
        assert summary["n_genes_total"] == expected, f"{family} {folder}"
        groups = result.tables["gene_groups"].rows
        assert sum(row["copies"] for row in groups if row["group"] != "Total") == expected
        assert next(row for row in groups if row["group"] == "Total")["copies"] == expected
        reconciled = {row["kind"]: row for row in result.tables["gene_reconciliation"].rows}
        assert all(row["agreement"] == "agrees" for row in reconciled.values()), folder
        assert reconciled["Total"]["counted"] == expected


@pytest.mark.parametrize("family", tuple(FAMILY_VARIANTS))
def test_orf_checks_each_cds_location_identically_under_either_annotation(
    tmp_path: Path, family: str
) -> None:
    """Annotation-dependent in its unit set, sequence-only in its verdict.

    A CDS both annotations assemble from the **same exon intervals** is the same string, so every
    cell of its check — status, codons, start/stop, internal stops — must be identical; only
    `gene` and `cds_id` (the annotation's naming) may differ. A shared outer span with a different
    exon split is a different string and is not compared. When the annotations disagree about the
    frame itself (`codon_start`, `transl_table`, `partial`, `pseudo`), the check follows that
    annotation, so those locations are named and excluded from the comparison.
    """
    accession = FAMILY_VARIANTS[family]
    by_parts: dict[str, dict[tuple, list[dict]]] = {}
    for folder, genome in zip(("source", "variant"), _pair(accession), strict=True):
        result = orf.run(genome, {}, _ctx(tmp_path / folder, "orf"))
        checks = result.tables["orf_checks"].rows
        cdss = [gene for gene in genome.genes if gene.kind == "CDS"]
        assert len(checks) == len(cdss), f"{family} {folder}: {len(checks)} rows, {len(cdss)} CDS"
        assert result.tables["orf_summary"].rows[0]["n_cds"] == len(checks)
        keyed: dict[tuple, list[dict]] = {}
        for gene, row in zip(cdss, checks, strict=True):
            keyed.setdefault(gene.parts, []).append(row)
        by_parts[folder] = keyed

    naming = {"gene", "cds_id"}
    frame = ("transl_table", "codon_start", "partial", "pseudo")

    def projection(rows: list[dict]) -> list[tuple[tuple[str, str], ...]]:
        return sorted(
            tuple(sorted((key, repr(value)) for key, value in row.items() if key not in naming))
            for row in rows
        )

    def frame_signature(rows: list[dict]) -> list[tuple]:
        return sorted(tuple(row[key] for key in frame) for row in rows)

    shared = set(by_parts["source"]) & set(by_parts["variant"])
    assert shared, f"{accession}: the two annotations share no CDS exon structure"
    checked = 0
    for parts in shared:
        left, right = by_parts["source"][parts], by_parts["variant"][parts]
        if frame_signature(left) != frame_signature(right):
            continue  # the annotation declares a different frame; its check legitimately follows
        assert projection(left) == projection(right), (
            f"{accession} {parts}: the check moved with the annotation"
        )
        checked += len(left)
    assert checked, f"{accession}: no CDS is assembled identically by both annotations"


@pytest.mark.parametrize("family", tuple(FAMILY_VARIANTS))
def test_codon_usage_counts_the_cds_set_the_annotation_declares(
    tmp_path: Path, family: str
) -> None:
    """Annotation-dependent: the totals are the sum over the CDS the annotation gives.

    The expectation is computed from the record, not from the module: the sense codons of each
    CDS's own coding sequence, and its internal stops. Both annotations are checked, so the change
    between them is shown to be the CDS set they declare.
    """
    accession = FAMILY_VARIANTS[family]
    for folder, genome in zip(("source", "variant"), _pair(accession), strict=True):
        result = codon_usage.run(genome, {}, _ctx(tmp_path / folder, "codon_usage"))
        summary = result.tables["codon_usage"].rows[0]
        cdss = cds_features(genome)
        codons = 0
        stops = 0
        for gene in cdss:
            sequence = coding_sequence(genome, gene)
            codons += sum(codon_counts(sequence).values())
            stops += translate(sequence)[:-1].count("*")
        assert summary["n_cds"] == len(cdss), f"{family} {folder}"
        assert summary["n_codons"] == codons, f"{family} {folder}"
        assert summary["internal_stops"] == stops
        assert sum(row["count"] for row in result.tables["rscu"].rows) == codons


@pytest.mark.parametrize("family", tuple(FAMILY_VARIANTS))
def test_amino_acid_composition_is_the_translation_of_the_cds_set(
    tmp_path: Path, family: str
) -> None:
    """Annotation-dependent through its CDS set; the composition is the translation itself.

    `n_residues + ignored_symbols` is the length of the translation the annotation's CDS produce,
    and the per-residue table sums to `n_residues` — on both annotations.
    """
    accession = FAMILY_VARIANTS[family]
    for folder, genome in zip(("source", "variant"), _pair(accession), strict=True):
        result = amino_acid.run(genome, {}, _ctx(tmp_path / folder, "amino_acid"))
        summary = result.tables["amino_acid"].rows[0]
        cdss = cds_features(genome)
        translation = "".join(translate(coding_sequence(genome, gene)) for gene in cdss)
        assert summary["n_cds"] == len(cdss), f"{family} {folder}"
        assert summary["n_residues"] == sum(composition(translation).values())
        assert summary["ignored_symbols"] == ignored_symbols(translation)
        assert (
            sum(row["count"] for row in result.tables["amino_acid_composition"].rows)
            == summary["n_residues"]
        )


@pytest.mark.parametrize("family", tuple(FAMILY_VARIANTS))
def test_stats_measures_one_observation_per_annotated_cds_locus(
    tmp_path: Path, family: str
) -> None:
    """Annotation-dependent: the observations are the loci the annotation declares.

    `n_observations` and `n_collapsed_cds` are recomputed from the records (canonical CDS locus,
    the IR's second copy collapsed) and must agree on both annotations. For a locus both
    annotations place identically, every per-locus test cell must be identical — except the
    Benjamini–Hochberg and Bonferroni corrections, which are computed across the run's whole
    family of tests and therefore move when the locus set moves.
    """
    accession = FAMILY_VARIANTS[family]
    source, variant = _pair(accession)
    partners = _partners(accession)
    groups = {accession: "A", **{partner.accession: "B" for partner in partners}}
    partner_units = [_stats_units(partner) for partner in partners]

    tests: dict[str, dict[tuple[str, str], dict]] = {}
    for folder, genome in zip(("source", "variant"), (source, variant), strict=True):
        result = stats.run(
            [genome, *partners], {"groups": groups}, _ctx(tmp_path / folder, "stats")
        )
        summary = result.tables["stats_summary"].rows[0]
        loci, collapsed = _stats_units(genome)
        all_loci = set(loci)
        for partner_loci, _ in partner_units:
            all_loci |= partner_loci
        assert summary["n_loci"] == len(all_loci), f"{family} {folder}"
        assert summary["n_observations"] == len(loci) + sum(
            len(partner_loci) for partner_loci, _ in partner_units
        )
        assert summary["n_collapsed_cds"] == collapsed + sum(count for _, count in partner_units)
        tests[folder] = {
            (row["locus"], row["measure"]): row for row in result.tables["stats_tests"].rows
        }

    correction = {"p_bh", "p_bonferroni", "significant"}
    same = _same_parts(source, variant)
    checked = 0
    for key in set(tests["source"]) & set(tests["variant"]):
        if key[0] not in same:
            continue
        left = {k: v for k, v in tests["source"][key].items() if k not in correction}
        right = {k: v for k, v in tests["variant"][key].items() if k not in correction}
        assert left == right, f"{accession} {key}: the locus's measurement moved"
        checked += 1
    assert checked, f"{accession}: no locus is placed identically by both annotations"


@pytest.mark.parametrize("family", WHOLE_GENOME_FAMILIES)
def test_barcodes_windows_are_the_sequence_and_the_candidates_the_annotation(
    tmp_path: Path, family: str
) -> None:
    """Mixed: the sliding windows read the aligned sequence; the candidate loci are annotated genes.

    The window and hypervariable-region tables must be identical between the two annotations, and
    so must every sequence-derived summary cell. The candidate table changes where the annotation
    does: `n_found` is the number of rows that say a candidate was found, a found row carries
    coordinates and a π, and a row that was not found names the gene the annotation lacks.
    """
    accession = FAMILY_VARIANTS[family]
    partner = _partners(accession, 1)[0]
    genomes = dict(zip(("source", "variant"), _pair(accession), strict=True))
    tables: dict[str, dict] = {}
    for folder, genome in genomes.items():
        ctx = _ctx(tmp_path / folder, "barcodes")
        ctx.ir = _gate(genome)
        tables[folder] = barcodes.run([genome, partner], dict(FAST_ALIGNMENT), ctx).tables

    for table_id in ("barcodes_windows", "barcodes_regions"):
        assert tables["source"][table_id].rows == tables["variant"][table_id].rows, (
            f"{accession}: {table_id} moved with the annotation"
        )
    source_summary = tables["source"]["barcodes_summary"].rows[0]
    variant_summary = tables["variant"]["barcodes_summary"].rows[0]
    for key in (
        "n_genomes",
        "window_bp",
        "step_bp",
        "n_windows",
        "mean_pi",
        "threshold_pi",
        "threshold_rule",
        "n_regions",
        "n_candidates",
        "tool",
    ):
        assert source_summary[key] == variant_summary[key], f"{accession} {key}"

    for folder, folder_tables in tables.items():
        summary = folder_tables["barcodes_summary"].rows[0]
        candidates = folder_tables["barcodes_candidates"].rows
        cdss = cds_by_canonical_name(genomes[folder])
        found = [row for row in candidates if row["found"] == "yes"]
        assert summary["n_found"] == len(found), f"{accession} {folder}"
        for row in candidates:
            if row["found"] != "yes":
                assert row["note"], f"{accession} {folder} {row['locus']}: not found, no reason"
                continue
            assert row["reference_start"] is not None and row["pi"] is not None
            if row["kind"] == "intergenic":
                # A spacer's length is its two flanking coordinates, in either orientation.
                assert row["span_bp"] == abs(row["reference_end"] - row["reference_start"]) + 1
            elif row["kind"] == "gene":
                # A gene candidate's span is the annotation's own feature length (exons only,
                # which is why it is not the outer span for a spliced CDS).
                feature = cdss.get(row["locus"])
                assert feature is not None, f"{accession} {folder}: {row['locus']} is not a CDS"
                assert row["span_bp"] == feature.length, f"{accession} {folder} {row['locus']}"


@pytest.mark.parametrize("family", WHOLE_GENOME_FAMILIES)
def test_diversity_whole_genome_tables_are_the_sequence_and_locus_tables_the_annotation(
    tmp_path: Path, family: str
) -> None:
    """Mixed: π along the alignment is the sequence; genes, introns and spacers are the annotation.

    The window and region tables and every sequence-derived summary cell must be identical. The
    locus tables must account for the annotation's own loci — their row counts are the counts the
    module publishes and the gene rows are the record's counted genes — and are allowed to differ
    between the two annotations.
    """
    accession = FAMILY_VARIANTS[family]
    partner = _partners(accession, 1)[0]
    tables: dict[str, dict] = {}
    for folder, genome in zip(("source", "variant"), _pair(accession), strict=True):
        ctx = _ctx(tmp_path / folder, "diversity")
        ctx.ir = _gate(genome)
        result = diversity.run([genome, partner], dict(FAST_ALIGNMENT), ctx)
        tables[folder] = result.tables
        summary = result.tables["diversity_summary"].rows[0]
        counted = count_genes(genome)
        assert len(result.tables["diversity_genes"].rows) == (
            counted.cds + counted.trna + counted.rrna
        ), f"{family} {folder}"
        assert summary["n_introns"] == len(result.tables["diversity_introns"].rows)
        assert summary["n_spacers"] == len(result.tables["diversity_igs"].rows)

    for table_id in ("diversity_windows", "diversity_regions"):
        assert tables["source"][table_id].rows == tables["variant"][table_id].rows, (
            f"{accession}: {table_id} moved with the annotation"
        )
    source_summary = tables["source"]["diversity_summary"].rows[0]
    variant_summary = tables["variant"]["diversity_summary"].rows[0]
    for key in (
        "n_genomes",
        "alignment_length",
        "usable_columns",
        "pi",
        "window_bp",
        "step_bp",
        "tool",
    ):
        assert source_summary[key] == variant_summary[key], f"{accession} {key}"


@pytest.mark.parametrize("family", WHOLE_GENOME_FAMILIES)
def test_substitutions_whole_genome_is_independent_of_the_annotation(
    tmp_path: Path, family: str
) -> None:
    """`scope=whole` aligns the sequence, so the annotation can move only its own label.

    The substitution, by-region and by-type tables must be identical cell for cell; the summary
    may differ only in `organism`, the record's own description, which is annotation data and is
    named here rather than compared away.
    """
    accession = FAMILY_VARIANTS[family]
    partner = _partners(accession, 1)[0]
    tables: dict[str, dict] = {}
    for folder, genome in zip(("source", "variant"), _pair(accession), strict=True):
        ctx = _ctx(tmp_path / folder, "substitutions")
        ctx.ir = _gate(genome)
        tables[folder] = substitutions.run([genome, partner], dict(FAST_ALIGNMENT), ctx).tables

    for table_id in ("substitutions", "substitutions_by_region", "substitutions_by_type"):
        assert tables["source"][table_id].rows == tables["variant"][table_id].rows, (
            f"{accession}: {table_id} moved with the annotation"
        )
    source_summary = tables["source"]["substitutions_summary"].rows[0]
    variant_summary = tables["variant"]["substitutions_summary"].rows[0]
    differing = {key for key in source_summary if source_summary[key] != variant_summary[key]}
    assert differing <= {"organism"}, (
        f"{accession}: annotation-independent cells moved: {differing}"
    )


@pytest.mark.parametrize("family", CODON_MODEL_FAMILIES)
def test_kaks_keeps_a_gene_measurement_when_both_annotations_place_it_identically(
    tools: dict[str, str | None], paml_tools: dict[str, str | None], tmp_path: Path, family: str
) -> None:
    """Annotation-dependent: the gene set and each gene's span are the annotation's.

    The partner genome is fixed, and the swapped record's sequence is identical, so for a gene
    both annotations place at the same parts the pair fed to MAFFT and yn00 is the same pair and
    every cell of its row must agree. Genes one annotation places differently may move; genes only
    one annotation has appear only in that run's gene universe.
    """
    require_tools(tools, "mafft")
    require_tools(paml_tools, "yn00")
    accession = FAMILY_VARIANTS[family]
    source, variant = _pair(accession)
    partner = _partners(accession, 1)[0]
    rows: dict[str, dict[str, dict]] = {}
    for folder, genome in zip(("source", "variant"), (source, variant), strict=True):
        result = kaks.run([genome, partner], {}, _ctx(tmp_path / folder, "kaks"))
        rows[folder] = {row["gene"]: row for row in result.tables["kaks_genes"].rows}
        universe = set(cds_by_canonical_name(genome)) | set(cds_by_canonical_name(partner))
        assert set(rows[folder]) == universe, (
            f"{accession} {folder}: a gene is missing from the run"
        )

    same = _same_parts(source, variant)
    checked = 0
    for name in same:
        assert rows["source"][name] == rows["variant"][name], (
            f"{accession} {name}: Ka/Ks moved although both annotations place the gene identically"
        )
        checked += 1
    assert checked, f"{accession}: no gene is placed identically by both annotations"


@pytest.mark.parametrize("family", CODON_MODEL_FAMILIES)
def test_popgen_keeps_a_gene_statistic_when_both_annotations_place_it_identically(
    tools: dict[str, str | None], tmp_path: Path, family: str
) -> None:
    """Annotation-dependent: the loci are the annotation's, the sequences are not.

    Three genomes: the swapped record alone in group A, two fixed partners in group B. Every
    partner sequence is the same in both runs, so for a locus both annotations place at the same
    parts the group rows and the Fst row must be identical; a locus whose parts differ may move,
    and a locus only one annotation has appears only in that run.
    """
    require_tools(tools, "mafft")
    accession = FAMILY_VARIANTS[family]
    source, variant = _pair(accession)
    partners = _partners(accession, 2)
    groups = {accession: "A", **{partner.accession: "B" for partner in partners}}

    def keyed(tables):
        group_rows = {(row["gene"], row["group"]): row for row in tables["popgen_groups"].rows}
        fst_rows = {row["gene"]: row for row in tables["popgen_fst"].rows}
        return group_rows, fst_rows

    runs = {}
    for folder, genome in zip(("source", "variant"), (source, variant), strict=True):
        ctx = _ctx(tmp_path / folder, "popgen")
        runs[folder] = keyed(popgen.run([genome, *partners], {"groups": groups}, ctx).tables)

    same = _same_parts(source, variant)
    checked = 0
    for name in same:
        for group in ("A", "B"):
            assert runs["source"][0].get((name, group)) == runs["variant"][0].get((name, group)), (
                f"{accession} {name} {group}: a group statistic moved with the annotation"
            )
            checked += 1
        assert runs["source"][1].get(name) == runs["variant"][1].get(name), (
            f"{accession} {name}: Fst moved although both annotations place the gene identically"
        )
    assert checked, f"{accession}: no gene is placed identically by both annotations"


# ---------------------------------------------------------------- sequence and production modules


def _iqtree_dir() -> Path:
    return iqtree_directory()


@pytest.mark.parametrize("family", WHOLE_GENOME_FAMILIES)
def test_structure_is_independent_of_the_annotation(tmp_path: Path, family: str) -> None:
    """The regions and the skew series come from the sequence's own IR; only the label moves."""
    accession = FAMILY_VARIANTS[family]
    tables: dict[str, dict] = {}
    for folder, genome in zip(("source", "variant"), _pair(accession), strict=True):
        ctx = _ctx(tmp_path / folder, "structure")
        ctx.ir = _gate(genome)
        tables[folder] = structure.run(genome, {}, ctx).tables
    for table_id in ("regions", "gc_skew_windows"):
        assert tables["source"][table_id].rows == tables["variant"][table_id].rows, table_id
    left = tables["source"]["structure_summary"].rows[0]
    right = tables["variant"]["structure_summary"].rows[0]
    assert {key for key in left if left[key] != right[key]} <= {"organism"}


@pytest.mark.parametrize("family", WHOLE_GENOME_FAMILIES)
def test_synteny_is_independent_of_the_annotation(
    mummer_tools: dict[str, str | None], tmp_path: Path, family: str
) -> None:
    """MUMmer aligns the normalized reference sequence to the query; no feature enters."""
    require_tools(mummer_tools, "nucmer", "delta-filter", "show-coords")
    accession = FAMILY_VARIANTS[family]
    partner = _partners(accession, 1)[0]
    tables: dict[str, dict] = {}
    for folder, genome in zip(("source", "variant"), _pair(accession), strict=True):
        ctx = _ctx(tmp_path / folder, "synteny")
        ctx.ir = _gate(genome)
        tables[folder] = synteny.run([genome, partner], {}, ctx).tables
    for table_id in ("synteny_genomes", "synteny_blocks", "synteny_identity"):
        assert tables["source"][table_id].rows == tables["variant"][table_id].rows, table_id
    left = tables["source"]["synteny_summary"].rows[0]
    right = tables["variant"]["synteny_summary"].rows[0]
    assert {key for key in left if left[key] != right[key]} <= {"organism"}


@pytest.mark.parametrize("family", WHOLE_GENOME_FAMILIES)
def test_submission_writes_the_sequence_once_and_the_annotation_it_was_given(
    tmp_path: Path, family: str
) -> None:
    """The FASTA is the molecule; the five-column table is the record's own features."""
    accession = FAMILY_VARIANTS[family]
    runs: dict[str, object] = {}
    for folder, genome in zip(("source", "variant"), _pair(accession), strict=True):
        runs[folder] = submission.run(genome, {}, _ctx(tmp_path / folder, "submission"))
    assert runs["source"].artifacts["submission.fsa"] == runs["variant"].artifacts["submission.fsa"]
    assert runs["source"].artifacts["submission.tbl"] != runs["variant"].artifacts["submission.tbl"]
    for folder, genome in zip(("source", "variant"), _pair(accession), strict=True):
        summary = runs[folder].tables["submission_summary"].rows[0]
        counted = count_genes(genome)
        assert (summary["n_cds"], summary["n_trna"], summary["n_rrna"]) == (
            counted.cds,
            counted.trna,
            counted.rrna,
        ), folder


@pytest.mark.parametrize("family", CODON_MODEL_FAMILIES)
def test_phylogeny_keeps_a_locus_that_both_annotations_place_identically(
    tools: dict[str, str | None], tmp_path: Path, family: str
) -> None:
    """The supermatrix is genes: a locus both annotations place identically is the same block.

    The tree may move because the locus set moves; the claim is about the block, not the topology.
    """
    require_tools(tools, "mafft")
    accession = FAMILY_VARIANTS[family]
    source, variant = _pair(accession)
    partners = _partners(accession, 3)
    params = {"iqtree": str(_iqtree_dir())}
    rows: dict[str, dict[str, dict]] = {}
    for folder, genome in zip(("source", "variant"), (source, variant), strict=True):
        result = phylogeny.run([genome, *partners], params, _ctx(tmp_path / folder, "phylogeny"))
        rows[folder] = {row["gene"]: row for row in result.tables["matrix_genes"].rows}
    same = _same_parts(source, variant)
    checked = 0
    for name in same:
        assert name in rows["source"] and name in rows["variant"], name
        assert rows["source"][name] == rows["variant"][name], (
            f"{accession} {name}: the block moved although both annotations place it identically"
        )
        checked += 1
    assert checked, f"{accession}: no gene is placed identically by both annotations"
