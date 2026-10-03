"""Gene-name normalization (P3-02): syntax, semantic aliases, and resolution by sequence.

The dialects these tests exercise were measured on the project's own dataset (PROJECT.md
§8.10): bracket tRNA names in one record, bare ``trnK``/``trnI`` in two, ``rrn16S`` beside
``rrn16``, and four functional names beside their classic ones.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from plastid.core.genbank import read_genbank, read_genbank_records
from plastid.core.gene_names import (
    IDENTICAL_RULE,
    SEMANTIC_RULE,
    SEQUENCE_RULE,
    SYNTAX_RULE,
    GeneAliasError,
    load_gene_aliases,
    normalize_genome,
    normalize_name,
    references_from,
)
from plastid.core.model import GENE_KINDS, Gene, Genome, Interval
from plastid.run.context import RunContext

ALIASES = load_gene_aliases()
_PAD = "ACGT" * 10
_TRNA = "GCG" + "TTT" * 20 + "CCA"  # 66 bp — length is arbitrary, identity is the point


def _tRNA_genome(name: str, sequence: str = _TRNA, accession: str = "SYN_TRNA.1") -> Genome:
    return Genome(
        name="synthetic",
        accession=accession,
        sequence=_PAD + sequence + _PAD,
        genes=(Gene(name, "tRNA", (Interval(len(_PAD), len(_PAD) + len(sequence)),)),),
        topology="linear",
    )


# --- syntax ---------------------------------------------------------------------------------


def test_the_bracket_dialect_normalizes_to_the_hyphen_form() -> None:
    normalized = normalize_name("trnL(UAA)", ALIASES)
    assert normalized.canonical == "trnL-UAA"
    assert normalized.name == "trnL(UAA)"  # the name it arrived under is kept
    assert normalized.rule == SYNTAX_RULE
    assert normalized.synonyms == ()  # a spelling difference is not a second name
    assert normalized.label == "trnL-UAA"


def test_the_lowercase_anticodon_letter_is_part_of_the_name() -> None:
    """trnfM-CAU: the f is the gene's name, not a case slip — only the trn/rrn prefix moves."""
    assert normalize_name("trnfM-CAU", ALIASES).canonical == "trnfM-CAU"
    assert normalize_name("TrnfM-CAU", ALIASES).canonical == "trnfM-CAU"
    assert normalize_name("TRNA-LEU", ALIASES).canonical == "trnA-LEU"


def test_the_rrna_size_suffix_is_dropped_but_the_mito_dialect_is_not() -> None:
    assert normalize_name("rrn16S", ALIASES).canonical == "rrn16"
    assert normalize_name("rrn23S", ALIASES).canonical == "rrn23"
    assert normalize_name("rrn4.5S", ALIASES).canonical == "rrn4.5"
    assert normalize_name("rrn5S", ALIASES).canonical == "rrn5"
    # rrnS / rrnL are mitochondrial names where the letter is the name, not a size class
    assert normalize_name("rrnS", ALIASES).canonical == "rrnS"
    assert normalize_name("rrnL", ALIASES).canonical == "rrnL"


def test_whitespace_is_collapsed_and_an_already_canonical_name_says_so() -> None:
    assert normalize_name("  rbcL  ", ALIASES).canonical == "rbcL"
    assert normalize_name("rbcL", ALIASES).rule == IDENTICAL_RULE


def test_normalizing_is_idempotent() -> None:
    for name in ("trnL(UAA)", "rrn16S", "ycf3", "pafI", "rbcL", "trnfM-CAU"):
        once = normalize_name(name, ALIASES)
        twice = normalize_name(once.canonical, ALIASES)
        assert twice.canonical == once.canonical, name


# --- semantic aliases -----------------------------------------------------------------------


def test_the_classic_name_resolves_to_the_teams_canonical_and_prints_both() -> None:
    """D3-02: canonical is the team's own dialect; the classic name travels beside it."""
    normalized = normalize_name("ycf3", ALIASES)
    assert normalized.canonical == "pafI"
    assert normalized.synonyms == ("ycf3",)
    assert normalized.label == "pafI (ycf3)"
    assert normalized.rule == SEMANTIC_RULE
    assert normalized.changed


def test_the_canonical_name_prints_the_same_label_as_its_synonym() -> None:
    """A reader must not be able to tell which dialect a record used — that is the whole point."""
    assert normalize_name("pafI", ALIASES).label == normalize_name("ycf3", ALIASES).label
    assert normalize_name("pbf1", ALIASES).label == "pbf1 (psbN)"
    assert normalize_name("clpP", ALIASES).label == "clpP1 (clpP)"
    assert normalize_name("ycf4", ALIASES).label == "pafII (ycf4)"


def test_the_four_measured_loci_are_the_only_ones_shipped() -> None:
    assert {locus.canonical for locus in ALIASES.loci} == {"pafI", "pafII", "pbf1", "clpP1"}
    assert all(locus.evidence for locus in ALIASES.loci)


def test_a_name_nobody_mapped_is_left_exactly_as_it_arrived() -> None:
    normalized = normalize_name("mystery9", ALIASES)
    assert normalized.canonical == "mystery9"
    assert normalized.rule == IDENTICAL_RULE
    assert not normalized.changed


def _write(tmp_path: Path, text: str) -> Path:
    path = tmp_path / "aliases.yaml"
    path.write_text(text, encoding="utf-8")
    return path


def test_an_alias_without_evidence_is_refused(tmp_path: Path) -> None:
    """The file's value is that it carries no guesses; an entry without evidence is one."""
    path = _write(
        tmp_path,
        "loci:\n  - canonical: pafI\n    synonyms: [ycf3]\n",
    )
    with pytest.raises(GeneAliasError, match="no evidence"):
        load_gene_aliases(path)


def test_a_name_in_two_loci_is_refused(tmp_path: Path) -> None:
    """Otherwise the file's order decides, silently, which locus a name belongs to."""
    path = _write(
        tmp_path,
        "loci:\n"
        "  - canonical: pafI\n    synonyms: [ycf3]\n    evidence: measured\n"
        "  - canonical: pbf1\n    synonyms: [ycf3]\n    evidence: measured\n",
    )
    with pytest.raises(GeneAliasError, match="exactly one locus"):
        load_gene_aliases(path)


def test_a_typo_in_an_alias_key_is_reported(tmp_path: Path) -> None:
    path = _write(
        tmp_path,
        "loci:\n  - canonical: pafI\n    synonym: [ycf3]\n    evidence: measured\n",
    )
    with pytest.raises(GeneAliasError, match="unknown key"):
        load_gene_aliases(path)


def test_a_locus_may_not_list_its_canonical_name_as_a_synonym(tmp_path: Path) -> None:
    path = _write(
        tmp_path,
        "loci:\n  - canonical: pafI\n    synonyms: [pafI]\n    evidence: measured\n",
    )
    with pytest.raises(GeneAliasError, match="own canonical name"):
        load_gene_aliases(path)


# --- bare tRNA names, resolved by sequence ---------------------------------------------------


def test_a_bare_trna_resolves_against_a_named_copy_at_full_identity() -> None:
    bare, named = _tRNA_genome("trnK"), _tRNA_genome("trnK-UUU")
    resolved = normalize_genome(bare, ALIASES, references=references_from([named]))
    normalized = next(iter(resolved.names.values()))
    assert normalized.canonical == "trnK-UUU"
    assert normalized.rule == SEQUENCE_RULE
    assert normalized.name == "trnK"  # still reported under the name it arrived with
    assert resolved.unresolved == ()


def test_a_reference_that_differs_by_one_base_does_not_resolve_the_name() -> None:
    """100 % identity or nothing: a near match is a different tRNA, or a different species."""
    bare = _tRNA_genome("trnK")
    close = _tRNA_genome("trnK-UUU", _TRNA[:30] + ("A" if _TRNA[30] != "A" else "C") + _TRNA[31:])
    resolved = normalize_genome(bare, ALIASES, references=references_from([close]))
    assert [item.name for item in resolved.unresolved] == ["trnK"]
    assert "no reference sequence matched" in resolved.unresolved[0].reason
    # reported, not renamed: the canonical form is still the bare name
    assert next(iter(resolved.names.values())).canonical == "trnK"


def test_two_reference_names_matching_equally_are_reported_as_ambiguous() -> None:
    bare = _tRNA_genome("trnK")
    twin_a, twin_b = _tRNA_genome("trnK-UUU"), _tRNA_genome("trnK-CUU")
    resolved = normalize_genome(bare, ALIASES, references=references_from([twin_a, twin_b]))
    assert [item.name for item in resolved.unresolved] == ["trnK"]
    assert "ambiguous" in resolved.unresolved[0].reason


def test_a_named_trna_is_never_touched_by_the_sequence_rule() -> None:
    resolved = normalize_genome(_tRNA_genome("trnK-UUU"), ALIASES)
    assert next(iter(resolved.names.values())).rule != SEQUENCE_RULE
    assert resolved.unresolved == ()


# --- one genome, end to end ------------------------------------------------------------------


def test_normalizing_a_genome_keeps_every_name_it_arrived_with(data_dir: Path) -> None:
    """The normalizer reports; it never rewrites. The golden genome's names are untouched."""
    golden = read_genbank(data_dir / "NC_050833.1.gb")
    names = normalize_genome(golden, ALIASES)

    assert [names.normalized(gene).name for gene in golden.genes] == [
        gene.name for gene in golden.genes
    ]
    counted = [gene for gene in golden.genes if gene.kind in GENE_KINDS]
    assert len({names.canonical(gene) for gene in counted}) == 112  # §6's distinct-name count


def test_the_golden_genome_reports_the_four_aliased_loci(data_dir: Path) -> None:
    golden = read_genbank(data_dir / "NC_050833.1.gb")
    names = normalize_genome(golden, ALIASES)
    labels = {gene.name: names.label(gene) for gene in golden.genes}
    assert labels["ycf3"] == "pafI (ycf3)"
    assert labels["ycf4"] == "pafII (ycf4)"
    assert labels["psbN"] == "pbf1 (psbN)"
    assert labels["clpP"] == "clpP1 (clpP)"
    assert labels["rbcL"] == "rbcL"


# --- the module (two records and up) ----------------------------------------------------------


def _ctx() -> RunContext:
    return RunContext(layout=None, params={}, module_id="gene_names", module_order=1)  # type: ignore[arg-type]


def test_the_module_resolves_a_bare_name_across_the_records_of_one_run() -> None:
    from plastid.modules import gene_names

    bare = _tRNA_genome("trnK", accession="A.1")
    named = _tRNA_genome("trnK-UUU", accession="B.1")
    result = gene_names.run([bare, named], {}, _ctx())

    panel = result.tables["panel_summary"].rows[0]
    assert panel["n_records"] == 2
    assert panel["n_resolved_by_sequence"] == 1
    assert panel["n_unresolved"] == 0
    assert "unresolved_names" not in result.tables

    mappings = {
        row["name"]: (row["canonical"], row["rule"]) for row in result.tables["name_mappings"].rows
    }
    assert mappings["trnK"] == ("trnK-UUU", "sequence")


def test_an_unresolved_name_is_reported_with_its_product_and_never_renamed() -> None:
    from plastid.modules import gene_names

    bare = _tRNA_genome("trnK", accession="A.1")
    other = _tRNA_genome(
        "trnK-UUU",
        _TRNA[:30] + ("A" if _TRNA[30] != "A" else "C") + _TRNA[31:],
        accession="B.1",
    )
    result = gene_names.run([bare, other], {}, _ctx())

    rows = result.tables["unresolved_names"].rows
    assert [row["name"] for row in rows] == ["trnK"]
    assert "no reference sequence matched" in rows[0]["reason"]
    assert result.tables["panel_summary"].rows[0]["n_unresolved"] == 1
    # the canonical form stays the bare name: reported, not renamed
    assert result.tables["name_summary"].rows[0]["n_distinct_canonical"] == 1


def test_the_module_summarises_the_two_dialects_of_one_locus() -> None:
    """A run holding both dialects must report one canonical locus, not two genes."""
    from plastid.core.model import Gene as Feature
    from plastid.modules import gene_names

    def _with_gene(name: str, accession: str) -> Genome:
        return Genome(
            name="synthetic",
            accession=accession,
            sequence=_PAD + "ATG" + "TTC" * 5 + "TAA" + _PAD,
            genes=(Feature(name, "CDS", (Interval(len(_PAD), len(_PAD) + 21),)),),
            topology="linear",
        )

    result = gene_names.run([_with_gene("ycf3", "A.1"), _with_gene("pafI", "B.1")], {}, _ctx())
    panel = result.tables["panel_summary"].rows[0]
    assert panel["n_distinct_names"] == 2
    assert panel["n_distinct_canonical"] == 1
    assert panel["n_aliased"] == 2
    assert {(row["name"], row["printed"]) for row in result.tables["name_mappings"].rows} == {
        ("ycf3", "pafI (ycf3)"),
        ("pafI", "pafI (ycf3)"),
    }


def _write_unannotated_genbank(path: Path) -> None:
    """A minimal valid record with no features — a stand-in for a third accession."""
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    record = SeqRecord(
        Seq("ATGC" * 250),
        id="SYN001.1",
        name="SYN001",
        description="synthetic record, no annotation",
        annotations={"organism": "Synthetic test", "topology": "linear", "molecule_type": "DNA"},
    )
    SeqIO.write([record], str(path), "genbank")


def test_the_cli_takes_two_or_more_records_for_this_module(tmp_path: Path, data_dir: Path) -> None:
    """API v0.2 says arity ≥ 2 means *at least* — the CLI used to refuse anything but exactly 2,
    which would have made an eleven-record panel impossible (regression)."""
    from plastid.cli import main

    golden = data_dir / "NC_050833.1.gb"
    comparison = data_dir / "NC_024936.1.gb"
    assert main(["run", "-m", "gene_names", "-i", str(golden), "-o", str(tmp_path)]) == 2
    assert (
        main(
            [
                "run",
                "-m",
                "gene_names",
                "-i",
                str(golden),
                "-i",
                str(comparison),
                "-o",
                str(tmp_path),
                "--dataset",
                "two",
            ]
        )
        == 0
    )
    # a third, distinct record: two --input flags prove "at least two", three prove "more is
    # allowed" (the contract in docs/API_v0.2.md)
    third = tmp_path / "SYN001.1.gb"
    _write_unannotated_genbank(third)
    assert (
        main(
            [
                "run",
                "-m",
                "gene_names",
                "-i",
                str(golden),
                "-i",
                str(comparison),
                "-i",
                str(third),
                "-o",
                str(tmp_path),
                "--dataset",
                "three",
            ]
        )
        == 0
    )


# --- the pinned panel (the real-data half of P3-02's Done-when) -------------------------------

PANEL = Path(__file__).resolve().parents[2] / ".panel" / "oldenlandia-coleus"
NEEDS_PANEL = pytest.mark.skipif(
    not PANEL.is_dir(),
    reason=(
        "the pinned panel is not on this machine — re-fetch it with the two 'plastid panel pin "
        "--groups' commands in PROJECT.md §8.6"
    ),
)


@NEEDS_PANEL
def test_the_panel_normalizes_every_name_and_resolves_the_bare_trna_names() -> None:
    """The Done-when's panel half: the four functional names and the bare tRNA names resolve,
    every name is idempotent, and no name is left unresolved or unclassified."""
    from plastid.core.gene_groups import assign, load_gene_groups

    groups = load_gene_groups()
    genomes = [
        genome for path in sorted(PANEL.glob("*.gb")) for genome in read_genbank_records(path)
    ]
    assert len(genomes) == 11
    references = references_from(genomes, ALIASES)

    arriving: set[str] = set()
    canonicals: set[str] = set()
    unresolved: list[str] = []
    unclassified: list[str] = []
    for genome in genomes:
        names = normalize_genome(genome, ALIASES, references=references)
        unresolved.extend(item.name for item in names.unresolved)
        for gene in genome.genes:
            if gene.kind not in GENE_KINDS:
                continue
            arriving.add(gene.name)
            canonicals.add(names.canonical(gene))
            # idempotence: normalizing a canonical form again changes nothing
            assert normalize_name(names.canonical(gene), ALIASES).canonical == names.canonical(gene)
        unclassified.extend(
            names.label(row.gene)
            for row in assign(genome, groups, canonical_of=names.canonical)
            if row.group is None
        )

    assert len(arriving) == 155  # the panel's distinct names
    assert unresolved == []
    assert unclassified == []

    # the two names that arrive with no anticodon at all, resolved by sequence identity
    labels = {
        (genome.accession, gene.name): normalize_genome(
            genome, ALIASES, references=references
        ).canonical(gene)
        for genome in genomes
        for gene in genome.genes
        if gene.kind == "tRNA"
    }
    assert labels[("PX923121.1", "trnK")] == "trnK-UUU"
    assert labels[("NC_069583.1", "trnI")] == "trnI-GAU"

    # the four loci under their functional names are one canonical each, not two genes
    for name in ("pafI", "pafII", "pbf1", "clpP1"):
        assert normalize_name(name, ALIASES).canonical == name
    assert {"pafI", "pafII", "pbf1", "clpP1"} <= canonicals


@NEEDS_PANEL
def test_the_panel_plus_the_golden_genome_is_156_distinct_names(data_dir: Path) -> None:
    """The count the board's Done-when rests on — measured, not assumed."""
    panel = {
        gene.name
        for path in sorted(PANEL.glob("*.gb"))
        for genome in read_genbank_records(path)
        for gene in genome.genes
        if gene.kind in GENE_KINDS
    }
    golden = read_genbank(data_dir / "NC_050833.1.gb")
    golden_names = {gene.name for gene in golden.genes if gene.kind in GENE_KINDS}
    assert len(panel) == 155 and len(golden_names) == 112
    assert len(panel | golden_names) == 156


def test_a_gene_feature_and_its_trna_feature_resolve_to_the_same_name() -> None:
    """Measured on PX923121.1: the bare trnK has both a `gene` feature (spanning the intron)
    and the tRNA feature. Resolving per feature compares an intron against a tRNA and reports
    the name unresolved — a false alarm this test pins down."""
    named = _tRNA_genome("trnK-UUU", accession="B.1")
    bare = Genome(
        name="synthetic",
        accession="A.1",
        sequence=_PAD + _TRNA + _PAD,
        genes=(
            Gene("trnK", "gene", (Interval(len(_PAD), len(_PAD) + len(_TRNA)),)),
            Gene("trnK", "tRNA", (Interval(len(_PAD), len(_PAD) + len(_TRNA)),)),
        ),
        topology="linear",
    )
    resolved = normalize_genome(bare, ALIASES, references=references_from([named]))
    assert resolved.unresolved == ()
    assert {entry.canonical for entry in resolved.names.values()} == {"trnK-UUU"}
