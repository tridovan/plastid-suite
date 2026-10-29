"""Gene-name normalization (P3-02): syntax, semantic aliases, and resolution by sequence.

The dialects these tests exercise were measured on the project's own dataset:
bracket tRNA names in one record, bare ``trnK``/``trnI`` in two, ``rrn16S`` beside
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
    # NC_020098.1 writes the bracket with a space; measured 29/10/2026, the older pattern left
    # `trnN (GUU)` unfolded while `trnN(GUU)` folded, so one locus entered a comparison as two
    # names. The space is not part of the name.
    assert normalize_name("trnN (GUU)", ALIASES).canonical == "trnN-GUU"
    assert normalize_name("trnI (CAU)", ALIASES).canonical == "trnI-CAU"


def test_the_lowercase_anticodon_letter_is_part_of_the_name() -> None:
    """trnfM-CAU: the f is the gene's name, not a case slip — only the trn/rrn prefix moves."""
    assert normalize_name("trnfM-CAU", ALIASES).canonical == "trnfM-CAU"
    assert normalize_name("TrnfM-CAU", ALIASES).canonical == "trnfM-CAU"


def test_a_spelled_out_amino_acid_is_a_bare_trna_not_alanine() -> None:
    """`tRNA-Ile` / `TRNA-LEU` name the amino acid, not a `trnA-…` gene.

    Measured 29/10/2026 while auditing the IR gene mirror on the panel: `NC_008641.1` names a tRNA
    `tRNA-Ile`. The case fold that lowercases a `trn`/`rrn` prefix read the fourth character `A`
    as the gene letter and turned it into `trnA-Ile` — an **alanine** name for an isoleucine tRNA.
    The spelled-out form maps to the bare `trn<letter>`; `normalize_genome` resolves a bare name
    to a named tRNA by sequence, so no anticodon is guessed here.
    """
    assert normalize_name("tRNA-Ile", ALIASES).canonical == "trnI"
    assert normalize_name("TRNA-LEU", ALIASES).canonical == "trnL"
    assert normalize_name("tRNA-Gly", ALIASES).canonical == "trnG"
    # An unknown spelled-out amino acid keeps its spelling rather than being corrupted to trnA-…
    assert normalize_name("tRNA-Xyz", ALIASES).canonical == "tRNA-Xyz"


def test_the_rrna_size_suffix_is_dropped_but_the_mito_dialect_is_not() -> None:
    assert normalize_name("rrn16S", ALIASES).canonical == "rrn16"
    assert normalize_name("rrn23S", ALIASES).canonical == "rrn23"
    assert normalize_name("rrn4.5S", ALIASES).canonical == "rrn4.5"
    assert normalize_name("rrn5S", ALIASES).canonical == "rrn5"
    # Measured 29/10/2026 on NC_028733.1: the suffix is written lowercase there.
    assert normalize_name("rrn23s", ALIASES).canonical == "rrn23"
    assert normalize_name("rrn4.5s", ALIASES).canonical == "rrn4.5"
    # rrnS / rrnL are mitochondrial names where the letter is the name, not a size class
    assert normalize_name("rrnS", ALIASES).canonical == "rrnS"
    assert normalize_name("rrnL", ALIASES).canonical == "rrnL"


def test_an_rrna_named_by_its_product_sentence_is_the_same_locus() -> None:
    """Measured 08/10/2026 on the pinned plastid reference set (P3B-07).

    ``NC_042827.1`` carries no ``/gene`` on its rRNA at all — the reader falls back to
    ``/product``, so "16S ribosomal RNA" arrived as a *gene name* and four loci no run could ever
    name. The size class is in the string, so this is the same rule as ``rrn16S``, not new
    knowledge.
    """
    assert normalize_name("16S ribosomal RNA", ALIASES).canonical == "rrn16"
    assert normalize_name("23S ribosomal RNA", ALIASES).canonical == "rrn23"
    assert normalize_name("4.5S ribosomal RNA", ALIASES).canonical == "rrn4.5"
    assert normalize_name("5S ribosomal RNA", ALIASES).canonical == "rrn5"
    assert normalize_name("16S rRNA", ALIASES).canonical == "rrn16"
    # a product sentence that is not a size class is left exactly as it arrived
    assert normalize_name("ribosomal protein S3", ALIASES).canonical == "ribosomal protein S3"
    assert normalize_name("cytochrome b6-f complex subunit V", ALIASES).canonical == (
        "cytochrome b6-f complex subunit V"
    )


def test_whitespace_is_collapsed_and_an_already_canonical_name_says_so() -> None:
    assert normalize_name("  rbcL  ", ALIASES).canonical == "rbcL"
    assert normalize_name("rbcL", ALIASES).rule == IDENTICAL_RULE


def test_normalizing_is_idempotent() -> None:
    for name in ("trnL(UAA)", "rrn16S", "ycf3", "pafI", "rbcL", "trnfM-CAU"):
        once = normalize_name(name, ALIASES)
        twice = normalize_name(once.canonical, ALIASES)
        assert twice.canonical == once.canonical, name


# --- semantic aliases -----------------------------------------------------------------------


def test_the_alias_resolves_to_the_naming_references_spelling_and_prints_both() -> None:
    """D5-09 (flipped 08/10/2026): canonical is the naming reference's spelling.

    The shipped direction is the **recorded standard**'s (``data/naming_reference.json``, the
    golden ``NC_050833.1``): it spells ``ycf3``, so ``pafI`` is the synonym that travels beside
    it. D3-02 had this the other way round — canonical was the group's dialect — and the whole
    point of the flip is that the direction is a project-level decision rather than a preference.
    """
    normalized = normalize_name("pafI", ALIASES)
    assert normalized.canonical == "ycf3"
    assert normalized.synonyms == ("pafI",)
    assert normalized.label == "ycf3 (pafI)"
    assert normalized.rule == SEMANTIC_RULE
    assert normalized.changed

    # the reference's own spelling is canonical and needs no change
    goldens = normalize_name("ycf3", ALIASES)
    assert goldens.canonical == "ycf3" and not goldens.changed


def test_the_canonical_name_prints_the_same_label_as_its_synonym() -> None:
    """A reader must not be able to tell which dialect a record used — that is the whole point."""
    assert normalize_name("pafI", ALIASES).label == normalize_name("ycf3", ALIASES).label
    assert normalize_name("pbf1", ALIASES).label == "psbN (pbf1)"
    assert normalize_name("clpP1", ALIASES).label == "clpP (clpP1)"
    assert normalize_name("pafII", ALIASES).label == "ycf4 (pafII)"


def test_the_measured_loci_are_the_only_ones_shipped() -> None:
    """A locus enters this table only with its evidence, so the set is pinned.

    ndhI joined it on 2026-10-03 (the golden record's ndhL and Chloe's ndhI are one locus with
    one protein); atpI, psaI and psbB joined it on 2026-10-04, when P4-02's Ka/Ks run listed the
    genes the two golden records "did not share" and three of them turned out to be the same
    locus under NC_024936.1's own spelling (see data/gene_aliases.yaml); and
    trnI-GAU the same day, when P4-04's barcode run listed the loci a record carries and
    NC_024936.1's trnl-GUA turned out to be that locus's tRNA, byte for byte.
    """
    assert {locus.canonical for locus in ALIASES.loci} == {
        # five flipped 08/10/2026 with D5-09 — the recorded reference spells the classic names
        "ycf3",
        "ycf4",
        "psbN",
        "clpP",
        "ndhL",
        # four kept their direction because the reference agrees with them
        "atpI",
        "psaI",
        "psbB",
        "trnI-GAU",
        # joined 10/10/2026 — NC_042827.1's lhbA against 33 records' psbZ
        "psbZ",
    }
    assert all(locus.evidence for locus in ALIASES.loci)


def test_the_golden_pair_shares_every_cds_name(data_dir: Path) -> None:
    """The regression the alias entries above exist for: before them, three whole loci — atpI,
    psaI and psbB — looked like genes NC_024936.1 does not carry, and every comparison that
    works gene by gene silently dropped them (atpI's CDS are byte-identical in the two records).
    """
    from plastid.core.gene_names import cds_by_canonical_name

    domestic = cds_by_canonical_name(read_genbank(data_dir / "NC_050833.1.gb"), ALIASES)
    gatesii = cds_by_canonical_name(read_genbank(data_dir / "NC_024936.1.gb"), ALIASES)
    assert set(domestic) == set(gatesii)
    assert {"atpI", "psaI", "psbB"} <= set(domestic)


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
    assert len({names.canonical(gene) for gene in counted}) == 112  # golden's distinct-name count


def test_the_golden_genome_reports_the_four_aliased_loci(data_dir: Path) -> None:
    """The golden spells five of the aliased loci itself (D5-09), so its labels are its own."""
    golden = read_genbank(data_dir / "NC_050833.1.gb")
    names = normalize_genome(golden, ALIASES)
    labels = {gene.name: names.label(gene) for gene in golden.genes}
    assert labels["ycf3"] == "ycf3 (pafI)"
    assert labels["ycf4"] == "ycf4 (pafII)"
    assert labels["psbN"] == "psbN (pbf1)"
    assert labels["clpP"] == "clpP (clpP1)"
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
    # both spellings of the one locus print the same label, and it is the reference's spelling
    assert {(row["name"], row["printed"]) for row in result.tables["name_mappings"].rows} == {
        ("ycf3", "ycf3 (pafI)"),
        ("pafI", "ycf3 (pafI)"),
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
        "--groups' commands"
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

    # the four loci under their functional names are one canonical each, not two genes — and
    # since 08/10/2026 that canonical form is the *naming reference's* spelling, so these names
    # (the panel's own dialect, not the golden's) resolve across to it
    for name, canonical in (
        ("pafI", "ycf3"),
        ("pafII", "ycf4"),
        ("pbf1", "psbN"),
        ("clpP1", "clpP"),
    ):
        assert normalize_name(name, ALIASES).canonical == canonical
        assert normalize_name(canonical, ALIASES).canonical == canonical
    # ... and the canonical forms themselves are the reference's, which is what the tables print
    assert {"ycf3", "ycf4", "psbN", "clpP"} <= canonicals
    assert not {"pafI", "pafII", "pbf1", "clpP1"} & canonicals


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


def test_cds_by_canonical_name_keeps_one_copy_per_name(data_dir: Path) -> None:
    """The shared-gene machinery of the comparison modules: canonical name → one CDS feature.

    The golden genome carries most genes twice (both inverted-repeat copies), and several under
    an alias — so the mapping has to normalize *and* deduplicate, or a gene would be counted
    twice in every comparison that uses it.
    """
    from plastid.core.gene_names import cds_by_canonical_name
    from plastid.core.model import cds_features

    genome = read_genbank(data_dir / "NC_050833.1.gb")
    by_name = cds_by_canonical_name(genome)
    assert len(by_name) == 79  # 86 CDS features; the difference is the second IR copy
    assert len(set(by_name)) == len(by_name)
    assert "rbcL" in by_name
    assert "ndhL" in by_name  # the golden's own spelling: the naming reference decides the
    #                           canonical form (D5-09), and the golden writes ndhL
    # every value is a real CDS of that genome, and at most one feature per canonical name
    features = set(cds_features(genome))
    assert set(by_name.values()) <= features


def test_the_amino_acid_letter_of_a_trna_name_is_case_insensitive() -> None:
    """`trnk-UUU` and `trnK-UUU` are one name in two spellings, and keeping them apart loses a
    locus in silence. Measured 2026-10-04 on the golden pair: NC_024936.1 writes `trnk-UUU` where
    the golden record writes `trnK-UUU`, with byte-identical sequences at the locus."""
    assert normalize_name("trnk-UUU", ALIASES).canonical == "trnK-UUU"
    assert normalize_name("trnk-UUU", ALIASES).rule == SYNTAX_RULE
    # trnl-GUA takes two steps: the case slip folds to trnL-GUA, and the alias table (below)
    # maps that to the locus's real name — 77 bp identical to the golden record's trnI-GAU
    assert normalize_name("trnl-GUA", ALIASES).canonical == "trnI-GAU"
    # trnfM-CAU's f is the name, not a case slip: the next character is M, not a hyphen
    assert normalize_name("trnfM-CAU", ALIASES).canonical == "trnfM-CAU"
    assert normalize_name("trnfM-CAU", ALIASES).rule == IDENTICAL_RULE


def test_nc_024936s_trnl_gua_is_the_locus_the_golden_record_calls_trni_gau(data_dir: Path) -> None:
    """The alias entry's regression: 77 bp, byte-identical to trnI-GAU, at the same position."""
    from plastid.core.model import coding_sequence

    records = {
        path.name: read_genbank(path)
        for path in (data_dir / "NC_050833.1.gb", data_dir / "NC_024936.1.gb")
    }

    def locus(genome, wanted):
        names = normalize_genome(genome, ALIASES)
        return [
            coding_sequence(genome, gene)
            for gene in genome.genes
            if gene.kind == "tRNA" and names.normalized(gene).canonical == wanted
        ]

    golden = locus(records["NC_050833.1.gb"], "trnI-GAU")
    gatesii = locus(records["NC_024936.1.gb"], "trnI-GAU")
    assert len(golden) == 2 and len(gatesii) == 2  # two IR copies each
    assert golden[0] == gatesii[0] and len(golden[0]) == 77


def test_the_records_own_organelle_picks_the_vocabulary() -> None:
    """Measured 2026-10-04: the mito half of Phase 4 found **no** shared gene between two
    mitogenomes, because OX459128.1 carries no /gene qualifier and every CDS "name" our reader
    reports is really a /product sentence. The vocabulary that maps those sentences already
    existed for the mito path — the comparison modules just were not asking for it."""
    from plastid.core.gene_names import cds_by_canonical_name

    mito = Genome(
        name="mito",
        accession="MITO.1",
        sequence="ACGT" * 30,
        genes=(Gene("NADH dehydrogenase subunit 1", "CDS", (Interval(0, 30),)),),
        source={"organelle": "mitochondrion"},
    )
    plastid = Genome(
        name="plastid",
        accession="PT.1",
        sequence="ACGT" * 30,
        genes=(Gene("NADH dehydrogenase subunit 1", "CDS", (Interval(0, 30),)),),
        source={"organelle": "plastid"},
    )
    # the same sentence: one record is a mitogenome and gets the symbol, the other is not
    assert set(cds_by_canonical_name(mito)) == {"nad1"}
    assert set(cds_by_canonical_name(plastid)) == {"NADH dehydrogenase subunit 1"}
