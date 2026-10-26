"""Reference-driven naming: the four tests, the three statuses, and the stop tier (P3B-04).

The design is D5-09; the numbers are the corpus measurement, and the two things this
file exists to hold are the ones the measurement changed. First, **the sequence tests can decide
for a CDS and not for a tRNA** — true and wrong tRNA loci overlap in identity (53.76–100 %
against up to 73.97 %), so a tRNA is never renamed on a ranking. Second, **a name the reference
does not use is a finding, never a silent correction**: the whole reason this step exists is that
NC_024936.1 spells `atpI` as `atpL` and the point is to *see* that, not to fix it quietly.

Negative cases are one per test, because a test that has never failed is a test nobody has read:

* test 1 — a name in neither the reference nor the alias table → ``outside vocabulary``;
* test 2 — the name is the reference's and the sequence ranks another locus first → ``conflict``;
* test 3 — the ranking is too close to trust → ``not measured``, and nothing is renamed;
* test 4 — the flanks contradict the name → reported, and for a tRNA that is all it can be.
"""

from __future__ import annotations

from dataclasses import replace
from datetime import datetime
from pathlib import Path

import pytest

from plastid.core.genbank import read_genbank
from plastid.core.gene_names import (
    COVERAGE_FLOOR,
    IDENTITY_FLOOR,
    MARGIN_FLOOR_BITS,
    OUTCOME_CONFIRMED,
    OUTCOME_CONFLICT,
    OUTCOME_RENAMED,
    OUTCOME_REPORTED,
    STATUS_MEASURED,
    STATUS_NOT_MEASURED,
    STATUS_OUTSIDE,
    TEST_MARGIN,
    TEST_NAME,
    TEST_NEIGHBOURS,
    TEST_SEQUENCE,
    AliasTable,
    GeneDecision,
    Hit,
    Locus,
    NamingConflictError,
    Neighbours,
    RecordedReference,
    checkout_root,
    decide_gene,
    duplicate_names,
    load_gene_aliases,
    load_recorded_reference,
    neighbours_of,
    nucleotide_identity,
    oriented_aliases,
    protein_fasta_rows,
    vocabulary_of,
)
from plastid.core.ir import detect_ir, enforce_ir_gate
from plastid.core.model import Gene, Genome, Interval
from plastid.modules import gene_normalize as naming
from plastid.run.context import RunContext
from plastid.run.layout import RunLayout

REPO = Path(__file__).resolve().parents[2]
GOLDEN = REPO / "tests" / "data" / "NC_050833.1.gb"
GATESII = REPO / "tests" / "data" / "NC_024936.1.gb"
# The corpus pair whose PGA half carries **no** /translation (measured 08/10/2026: 0 of 84 CDS,
# against 86 of 86 in the published half) — so the naming step's name-only path runs end to end
# without an aligner, on every machine. The BLAST-dependent half is a tier-B test below.
HERMONA_PUBLISHED = REPO / "tests" / "data" / "corpus" / "Iris_hermona_NC_081931.1.published.gb"
HERMONA_PGA = REPO / "tests" / "data" / "corpus" / "Iris_hermona_NC_081931.1.pga.gb"
PGA_DOMESTICA = REPO / "tests" / "data" / "corpus" / "Iris_domestica_NC_050833.1.pga.gb"

ALIASES = load_gene_aliases()

# A reference vocabulary small enough to reason about: one locus, spelled two ways.
LOCUS = Locus(canonical="ycf3", synonyms=("pafI",), evidence="the file says so")
ONE = AliasTable(loci=(LOCUS,))


def _gene(name: str, kind: str = "CDS", *, length: int = 30) -> Gene:
    return Gene(name=name, kind=kind, parts=(Interval(0, length),))


def _cds(name: str) -> Gene:
    return replace(_gene(name), translation="M" * 10)


def _hits(*rows: tuple[str, float, float, float | None]) -> tuple[Hit, ...]:
    """``(locus, identity, coverage)`` or ``(…, bits)`` — bits only exist for a protein search."""
    return tuple(
        Hit(
            locus=row[0],
            identity=row[1],
            coverage=row[2],
            bits=row[3] if len(row) > 3 else None,
        )
        for row in rows
    )


# ---------------------------------------------------------------- the measured thresholds


def test_the_thresholds_are_the_measured_ones() -> None:
    """Measured: the smallest true margin was 19.2 bits, and with coverage ≥ 0.7 no wrong
    locus reached 44.93 % while no true ortholog fell below 67.75 %."""
    assert (IDENTITY_FLOOR, COVERAGE_FLOOR, MARGIN_FLOOR_BITS) == (50.0, 0.7, 10.0)


def test_the_alias_direction_follows_the_reference_not_the_file() -> None:
    """D5-09: the file says which names are one locus; the reference says which one is printed.

    Two references, one alias file, two answers — which is the whole point of orienting the table
    at run time rather than editing the YAML every time the ★ changes.
    """
    spells_classic = Genome(name="x", sequence="ACGT", genes=(_gene("ycf3"),))
    spells_group = Genome(name="y", sequence="ACGT", genes=(_gene("pafI"),))
    reference = AliasTable(loci=(Locus("pafI", ("ycf3",), "as recorded"),))

    classic = oriented_aliases(reference, spells_classic)
    group = oriented_aliases(reference, spells_group)

    assert classic.loci[0].canonical == "ycf3"
    assert classic.loci[0].synonyms == ("pafI",)
    assert group.loci[0].canonical == "pafI"  # the file's direction, because the reference agrees
    # a locus the reference says nothing about keeps the file's direction
    genome = Genome(name="z", sequence="ACGT", genes=(_gene("rbcL"),))
    assert oriented_aliases(reference, genome).loci[0].canonical == "pafI"


# ---------------------------------------------------------------- test 1: the name


def test_a_name_the_reference_does_not_use_is_renamed_only_on_evidence() -> None:
    """The measured case this step exists for: NC_024936.1 spells `atpI` as `atpL`.

    Two ways the step may find out, and both say which one it was: the alias table knows the
    locus (rule ``name``), or the sequence identifies it (rule ``sequence``). What it may never do
    is decide from the name alone that a gene "must be" something else.
    """
    renamed = decide_gene(
        _cds("atpL"),
        molecules="protein",
        vocabulary={"atpI": _cds("atpI")},
        hits=_hits(("atpI", 99.0, 1.0, 900.0)),
        aliases=AliasTable(),
    )
    assert renamed.outcome == OUTCOME_RENAMED
    assert renamed.canonical == "atpI" and renamed.rule == TEST_SEQUENCE

    # no usable evidence — the name stays as it arrived and the row says why
    unplaced = decide_gene(
        _cds("atpL"),
        molecules="protein",
        vocabulary={"rbcL": _cds("rbcL")},
        hits=_hits(("rbcL", 21.0, 0.3, 30.0)),
        aliases=AliasTable(),
    )
    assert unplaced.status == STATUS_OUTSIDE and unplaced.canonical is None
    assert unplaced.outcome == OUTCOME_REPORTED


def test_the_alias_table_turns_a_synonym_into_the_references_spelling() -> None:
    """Two ways to the same rename, and the ``rule`` column says which one carried it.

    With no sequence to check, the alias table alone decides (rule ``name``); with a sequence that
    agrees, the row reports the stronger evidence, because that is what a reader is checking.
    """
    unbacked = decide_gene(
        _cds("pafI"),
        molecules="protein",
        vocabulary={"ycf3": _cds("ycf3")},
        hits=(),
        aliases=ONE,
    )
    assert unbacked.canonical == "ycf3"
    assert unbacked.outcome == OUTCOME_RENAMED
    assert unbacked.rule == TEST_NAME

    backed = decide_gene(
        _cds("pafI"),
        molecules="protein",
        vocabulary={"ycf3": _cds("ycf3")},
        hits=_hits(("ycf3", 99.0, 1.0, 900.0)),
        aliases=ONE,
    )
    assert backed.outcome == OUTCOME_RENAMED
    assert backed.rule == TEST_SEQUENCE

    # and a record that already spells it the reference's way is not a rename at all
    already = decide_gene(
        _cds("ycf3"),
        molecules="protein",
        vocabulary={"ycf3": _cds("ycf3")},
        hits=_hits(("ycf3", 99.0, 1.0, 900.0)),
        aliases=ONE,
    )
    assert already.outcome == OUTCOME_CONFIRMED


# ---------------------------------------------------------------- test 2: the sequence


def test_a_name_whose_sequence_ranks_another_locus_first_is_a_conflict() -> None:
    """D5-09's stop tier: same name, different locus — the one case that stops a CDS."""
    decision = decide_gene(
        _cds("rbcL"),
        molecules="protein",
        vocabulary={"rbcL": _cds("rbcL"), "rbcS": _cds("rbcS")},
        hits=_hits(("rbcS", 98.0, 1.0, 900.0), ("rbcL", 30.0, 0.4, 60.0)),
        aliases=AliasTable(),
    )

    assert decision.outcome == OUTCOME_CONFLICT
    assert decision.status == STATUS_MEASURED  # the check ran; it is the *answer* that disagrees
    assert decision.rule == TEST_SEQUENCE
    assert "rbcS" in decision.reason


def test_a_name_and_a_sequence_that_agree_are_confirmed_with_the_rule_that_decided() -> None:
    confirmed = decide_gene(
        _cds("rbcL"),
        molecules="protein",
        vocabulary={"rbcL": _cds("rbcL")},
        hits=_hits(("rbcL", 99.0, 0.99, 900.0)),
        aliases=AliasTable(),
    )
    assert confirmed.outcome == OUTCOME_CONFIRMED and confirmed.rule == TEST_SEQUENCE

    ranked = decide_gene(
        _cds("rbcL"),
        molecules="protein",
        vocabulary={"rbcL": _cds("rbcL")},
        hits=_hits(("rbcL", 99.0, 0.99, 900.0), ("rbcS", 40.0, 0.5, 100.0)),
        aliases=AliasTable(),
        neighbours=Neighbours("LSC", "atpB", "accD"),
        reference_neighbours=Neighbours("LSC", "atpB", "accD"),
    )
    assert ranked.rule == TEST_MARGIN  # a runner-up exists, so the margin is what was read
    assert ranked.margin_bits == pytest.approx(800.0)


def test_a_unique_match_renames_a_cds_the_name_alone_could_not_place() -> None:
    decision = decide_gene(
        _cds("ysp1"),
        molecules="protein",
        vocabulary={"ycf3": _cds("ycf3")},
        hits=_hits(("ycf3", 88.0, 0.95, 500.0)),
        aliases=AliasTable(),
    )

    assert decision.outcome == OUTCOME_RENAMED
    assert decision.canonical == "ycf3"
    assert decision.status == STATUS_MEASURED


# ---------------------------------------------------------------- test 3: the margin


def test_a_margin_too_narrow_to_trust_renames_nothing() -> None:
    """The negative case of test 3: two loci 4 bits apart is not evidence, so nothing is decided."""
    decision = decide_gene(
        _cds("ysp1"),
        molecules="protein",
        vocabulary={"ycf3": _cds("ycf3"), "ycf4": _cds("ycf4")},
        hits=_hits(("ycf3", 60.0, 0.8, 104.0), ("ycf4", 59.0, 0.8, 100.0)),
        aliases=AliasTable(),
    )

    assert decision.status == STATUS_NOT_MEASURED  # "could not tell", never "none found"
    assert decision.canonical is None
    assert decision.outcome == OUTCOME_REPORTED
    assert decision.rule == TEST_MARGIN
    assert "4.0 bits" in decision.reason


def test_a_hit_below_the_identity_floor_is_not_a_counterpart() -> None:
    """Measured: the best a genuinely different locus reaches is 44.93 % at coverage ≥ 0.7."""
    decision = decide_gene(
        _cds("ysp1"),
        molecules="protein",
        vocabulary={"ycf3": _cds("ycf3")},
        hits=_hits(("ycf3", 44.9, 0.8, 90.0)),
        aliases=AliasTable(),
    )

    assert decision.status == STATUS_OUTSIDE
    assert decision.canonical is None


# ---------------------------------------------------------------- test 4: the neighbours


def test_neighbours_confirm_what_the_ranking_could_not() -> None:
    """The serine-tRNA case: the ranking prefers another locus, the flanks say the name is right."""
    decision = decide_gene(
        _cds("rbcL"),
        molecules="protein",
        vocabulary={"rbcL": _cds("rbcL"), "rbcS": _cds("rbcS")},
        hits=_hits(("rbcS", 70.0, 0.9, 300.0), ("rbcL", 68.0, 0.9, 290.0)),
        aliases=AliasTable(),
        neighbours=Neighbours("LSC", "atpB", "accD"),
        reference_neighbours=Neighbours("LSC", "atpB", "accD"),
    )

    assert decision.outcome == OUTCOME_CONFIRMED
    assert decision.rule == TEST_NEIGHBOURS
    assert decision.canonical == "rbcL"


def test_neighbours_that_disagree_report_a_trna_and_never_rename_or_stop_it() -> None:
    """For a tRNA the ranking decides nothing, so a disagreement is a report."""
    decision = decide_gene(
        _gene("trnS-GCU", "tRNA"),
        molecules="nucleotide",
        vocabulary={"trnS-GCU": _gene("trnS-GCU", "tRNA")},
        hits=_hits(("trnS-UGA", 90.0, 1.0)),
        aliases=AliasTable(),
        neighbours=Neighbours("LSC", "psbA", "trnG-UCC"),
        reference_neighbours=Neighbours("LSC", "trnS-GGA", "trnG-UCC"),
    )

    assert decision.status == STATUS_NOT_MEASURED
    assert decision.outcome == OUTCOME_REPORTED
    assert decision.canonical == "trnS-GCU"  # nothing changes in the output
    assert decision.rule == TEST_NEIGHBOURS


def test_a_trna_is_never_renamed_on_a_ranking() -> None:
    """The measurement's own rule: 26 of 1,182 cross-record tRNA comparisons misrank."""
    decision = decide_gene(
        _gene("trnS-GCU", "tRNA"),
        molecules="nucleotide",
        vocabulary={"trnS-UGA": _gene("trnS-UGA", "tRNA")},
        hits=_hits(("trnS-UGA", 95.0, 1.0)),
        aliases=AliasTable(),
        neighbours=Neighbours("LSC", "psbA", "trnG-UCC"),
        reference_neighbours=Neighbours("LSC", "psbA", "trnG-UCC"),
    )

    assert decision.status == STATUS_OUTSIDE
    assert decision.canonical is None
    assert decision.outcome == OUTCOME_REPORTED


def test_a_gene_with_no_way_to_check_it_says_so_rather_than_agreeing() -> None:
    """No /translation and no alignment: the CDS arrives as it was, and the blank is named."""
    decision = decide_gene(
        _gene("rbcL"),
        molecules="protein",
        vocabulary={"rbcL": _cds("rbcL")},
        hits=(),
        aliases=AliasTable(),
    )

    assert decision.status == STATUS_NOT_MEASURED
    assert decision.outcome == OUTCOME_REPORTED
    assert "no /translation" in decision.reason


# ---------------------------------------------------------------- the reference's own shape


def test_an_inverted_repeat_is_not_a_duplicate_name() -> None:
    """Two identical copies of one gene are one locus; two *different* genes under one name are
    the golden's own trnG problem (73.24 % identical) — and only the second is a conflict."""
    copies = Genome(
        name="x",
        sequence="ACGT" * 40,
        genes=(
            Gene(name="rrn23", kind="rRNA", parts=(Interval(0, 20),)),
            Gene(name="rrn23", kind="rRNA", parts=(Interval(80, 100),)),
            Gene(name="trnG-UCC", kind="tRNA", parts=(Interval(20, 40),)),
            Gene(name="trnG-UCC", kind="tRNA", parts=(Interval(100, 130),)),
        ),
    )

    assert duplicate_names(copies) == ("trnG-UCC",)


def test_the_protein_fasta_carries_the_subject_key_and_skips_records_without_proteins() -> None:
    genome = Genome(
        name="Iris",
        sequence="ACGT",
        accession="NC_1.1",
        genes=(_cds("rbcL"), _gene("no_translation")),
    )

    assert protein_fasta_rows(genome, vocabulary="A") == (("NC_1.1|rbcL|A", "M" * 10),)


def test_the_nucleotide_metric_is_the_one_the_measurement_recorded() -> None:
    """The measurement's parameters are part of the number.

    Cheaper gaps read a different identity.
    """
    assert nucleotide_identity("ACGTACGT", "ACGTACGT") == 100.0
    assert nucleotide_identity("ACGTACGT", "ACGTACGA") == pytest.approx(87.5)
    assert nucleotide_identity("", "ACGT") == 0.0


def test_the_neighbours_are_read_inside_one_region() -> None:
    """rps12 straddles the LSC/IRa junction: which copy flanks a gene depends on the region."""
    genome = read_genbank(GOLDEN)
    aliases = oriented_aliases(ALIASES, genome)
    ir = enforce_ir_gate(detect_ir(genome))
    trna = next(gene for gene in genome.genes if gene.name == "trnA-UGC")

    neighbours = neighbours_of(genome, trna, ir, aliases)

    assert neighbours is not None
    assert neighbours.region in ("IRa", "IRb")
    assert neighbours.left == "trnI-GAU"


# ---------------------------------------------------------------- the module end to end


def _without_translations(genome: Genome) -> Genome:
    """The same record with its proteins removed — the PGA annotations' own shape.

    A record that carries no ``/translation`` cannot be sequence-checked, so the naming step runs
    its name-only path and needs no aligner: that is how the module's end-to-end behaviour is
    tested on every machine, and how the corpus's PGA halves actually arrive.
    """
    return replace(
        genome,
        genes=tuple(replace(gene, translation=None) for gene in genome.genes),
    )


def _context(tmp_path: Path) -> RunContext:
    layout = RunLayout.create(
        tmp_path, "naming", fingerprint="0" * 16, now=datetime(2026, 10, 8, 12, 0, 0)
    )
    layout.make_dirs()
    return RunContext(layout=layout, params={}, module_id=naming.MODULE_ID, module_order=1)


def test_the_corpus_pair_normalizes_into_the_references_spelling(tmp_path: Path) -> None:
    """The Done-when's first half: one sequence, two annotations → the artifact and a pairing table.

    *Iris hermona* twice over: the published annotation is the reference, and the PGA annotation
    of the **same sequence** is what gets normalized. The PGA half carries no ``/translation``
    (0 of 84 CDS, measured 08/10/2026), so nothing can be sequence-checked and every CDS says so
    — which is the honest answer, and the reason the status column exists.
    """
    reference = read_genbank(HERMONA_PUBLISHED)
    record = _without_translations(read_genbank(HERMONA_PGA))

    result = naming.run((reference, record), {"mode": "normalize"}, _context(tmp_path))

    assert set(result.artifacts) == {"Iris_hermona_voucher_IrisSV1903_NC_081931.normalized.gb"}
    written = result.artifacts["Iris_hermona_voucher_IrisSV1903_NC_081931.normalized.gb"]
    assert "/translation=" not in written  # the record had none to keep
    assert result.module == "gene_normalize"

    pairing = result.tables["pairing"]
    assert pairing.columns[2].name == "gene"  # B's spelling, in its own column
    statuses = {row["status"] for row in pairing.rows}
    # the CDS cannot be checked (no protein) and say so; a tRNA whose name the reference carries
    # is 'measured' — for a tRNA the name *is* the check the measurement left standing
    assert statuses <= {STATUS_MEASURED, STATUS_NOT_MEASURED, STATUS_OUTSIDE}
    assert STATUS_NOT_MEASURED in statuses
    # The PGA record's IR now passes the gate (M11.4: two internal gaps, 99.992 % identity), so
    # the neighbours test is available and legitimately decides some tRNA pairings.
    assert all(row["rule"] in (None, TEST_NAME, TEST_NEIGHBOURS) for row in pairing.rows)
    assert any(row["rule"] == TEST_NEIGHBOURS for row in pairing.rows)
    summary = result.tables["normalize_summary"].rows[0]
    assert summary["reference"] == "NC_081931.1"
    assert summary["records"] == 1
    assert summary["genes"] == len(pairing.rows)
    assert any("no CDS carried a /translation" in note for note in result.notes)


def test_the_reference_normalizes_to_itself_with_nothing_to_change(tmp_path: Path) -> None:
    """Verify mode's promise on the record that *is* the vocabulary: zero changes.

    The measured standard is the golden, whose names are already canonical. *Iris hermona*'s
    published record is not: it spells ``trnk-UUU`` and ``trnl-GUA`` — lowercase slips the syntax
    rules correct on any record, the reference included (D5-09's automatic tier), which is why
    this test uses the record the project actually names after.
    """
    golden = _without_translations(read_genbank(GOLDEN))

    result = naming.run((golden, golden), {"mode": "verify"}, _context(tmp_path))

    assert result.artifacts == {}  # verify writes no record
    assert result.tables["normalize_summary"].rows[0]["renamed"] == 0
    assert any("nothing would change" in note for note in result.notes)


def test_a_rename_that_would_duplicate_a_name_stops_the_run(tmp_path: Path) -> None:
    """D5-09's stop tier, and the coordinates come with it."""
    reference = read_genbank(HERMONA_PUBLISHED)
    record = _without_translations(read_genbank(HERMONA_PGA))
    # the record already carries psbN; adding the other spelling of that locus makes the
    # normalized record hold two different genes under one name
    duplicated = replace(
        record,
        genes=record.genes + (Gene(name="pbf1", kind="CDS", parts=(Interval(50_000, 50_500),)),),
    )

    with pytest.raises(NamingConflictError) as caught:
        naming.run((reference, duplicated), {"mode": "normalize"}, _context(tmp_path))

    message = str(caught.value)
    assert "pbf1" in message or "psbN" in message
    assert "(50001..50500+)" in message  # 1-based inclusive, the record's own coordinates


def test_on_conflict_continue_writes_the_output_and_the_conflicts_table(tmp_path: Path) -> None:
    reference = read_genbank(HERMONA_PUBLISHED)
    record = _without_translations(read_genbank(HERMONA_PGA))
    duplicated = replace(
        record,
        genes=record.genes + (Gene(name="pbf1", kind="CDS", parts=(Interval(50_000, 50_500),)),),
    )

    result = naming.run(
        (reference, duplicated),
        {"mode": "normalize", "on_conflict": "continue"},
        _context(tmp_path),
    )

    assert "conflicts" in result.tables
    assert result.artifacts  # the record is written, with the conflict named in the table
    assert result.tables["conflicts"].rows
    assert any("carried through" in note for note in result.notes)


def test_every_change_in_the_artifact_is_a_row_a_reader_can_check(tmp_path: Path) -> None:
    """Nothing is renamed in place and nothing is renamed quietly.

    The names that differ between the input and the artifact must be exactly the ones the pairing
    table calls renames — and the features themselves (kind, coordinates, order) must be the ones
    the record arrived with, because a step that renames is not a step that re-annotates.
    """
    reference = read_genbank(HERMONA_PUBLISHED)
    record = _without_translations(read_genbank(HERMONA_PGA))
    result = naming.run((reference, record), {"mode": "normalize"}, _context(tmp_path))
    written = result.artifacts["Iris_hermona_voucher_IrisSV1903_NC_081931.normalized.gb"]
    back = read_genbank(_write(tmp_path, written))

    assert len(back.genes) == len(record.genes)
    assert [gene.parts for gene in back.genes] == [gene.parts for gene in record.genes]
    in_table = {
        (row["gene"], row["canonical"])
        for row in result.tables["pairing"].rows
        if row["outcome"] == "renamed"
    }
    assert in_table  # this record does carry names that change
    in_file = {
        (before.name, after.name)
        for before, after in zip(record.genes, back.genes, strict=True)
        if before.name != after.name
    }
    assert in_file == in_table


def _write(tmp_path: Path, text: str) -> Path:
    path = tmp_path / "normalized.gb"
    path.write_text(text, encoding="utf-8")
    return path


def test_the_recorded_reference_is_what_the_audit_and_the_flip_agree_about() -> None:
    """The shipped direction, the recorded standard and the file on disk are one statement.

    A guard rather than a fixture: if any of the three drifts, every table this project prints
    disagrees with the record it claims to name after — and the drift would be invisible.
    """
    import hashlib

    recorded = load_recorded_reference()
    root = checkout_root()
    assert root is not None, "this test runs from a checkout, so a root must be found"
    path = recorded.resolve(root=root)
    digest = hashlib.sha256(path.read_bytes()).hexdigest()

    assert digest == recorded.sha256
    assert recorded.provisional, "the end user's pick is still open (D5-09)"
    assert path.name == "NC_050833.1.gb"

    golden = read_genbank(path)
    aliases = oriented_aliases(ALIASES, golden)
    # the direction of the shipped table is the standard's own spelling: five loci flipped
    # (D5-09) and four did not
    for locus in aliases.loci:
        assert locus.canonical in {gene.name for gene in golden.genes}, (
            f"{locus.canonical} is not how the recorded reference spells that locus"
        )


def test_the_spec_declares_the_reference_parameter_and_the_arity_that_counts_it() -> None:
    from plastid.modules import registry

    spec = registry.get(naming.MODULE_ID)

    assert spec.reference_param == "naming_reference"  # not `reference`: that key has two meanings
    assert spec.arity == 2  # the reference plus at least one record
    assert spec.requires_ir is False  # names are not a position-based analysis


def test_a_decision_carries_every_column_the_tables_print() -> None:
    """The pairing table's contract, asserted where it is defined rather than where it is read."""
    decision = GeneDecision(
        gene=_cds("rbcL"),
        status=STATUS_MEASURED,
        outcome=OUTCOME_CONFIRMED,
        canonical="rbcL",
        rule=TEST_SEQUENCE,
        identity=99.0,
        coverage=1.0,
        margin_bits=250.0,
        best_locus="rbcL",
        reason="the name and the sequence agree",
    )

    assert decision.renamed is False and decision.conflict is False
    chosen = vocabulary_of(Genome(name="x", sequence="A", genes=(_cds("rbcL"),)))
    assert set(chosen) == {"rbcL"}


# ---------------------------------------------------------------- the sequence path, with BLAST


@pytest.mark.tier_b
def test_the_sequence_tests_decide_on_real_proteins(
    tmp_path: Path, blast_tools: dict[str, str | None]
) -> None:
    """The measured case, end to end: *Iris gatesii* normalized against the golden.

    Measured 08/10/2026 — 6 renames and 121 confirmations, with the rules that carried them:
    ``atpL``→``atpI`` and ``psaL``→``psaI`` at 100 % identity over the whole protein, ``bspB``→
    ``psbB`` at 99.80 %, and the two syntax slips by gene order. Five tRNA names the reference
    does not use are **reported** rather than renamed — one of them, ``trnL-CAU``, is 100 %
    identical to the reference's ``trnI-CAU``, which is exactly the naming difference a ranking
    may not settle (for tRNA, identity decides nothing).
    """
    from conftest import require_tools

    require_tools(blast_tools, "blastp", "makeblastdb")
    reference = read_genbank(GOLDEN)
    record = read_genbank(GATESII)
    params = {"mode": "normalize", "tools": str(Path(blast_tools["blastp"]).parent)}

    result = naming.run((reference, record), params, _context(tmp_path))

    rows = {row["gene"]: row for row in result.tables["pairing"].rows}
    assert rows["atpL"]["canonical"] == "atpI"
    assert rows["atpL"]["outcome"] == OUTCOME_RENAMED
    assert rows["atpL"]["identity_pct"] == pytest.approx(100.0)
    assert rows["bspB"]["canonical"] == "psbB"
    assert rows["bspB"]["identity_pct"] == pytest.approx(99.8, abs=0.05)

    summary = result.tables["normalize_summary"].rows[0]
    assert summary["conflicts"] == 0
    assert summary["renamed"] >= 5
    assert summary["outside"] >= 3
    # the artifact carries the proteins it arrived with — the annotator searches these
    written = result.artifacts["NC_024936.1.normalized.gb"]
    assert written.count("/translation=") == len(
        [gene for gene in record.genes if gene.kind == "CDS" and gene.translation]
    )


# ---------------------------------------------------------------- the door


def test_the_run_resolves_audits_and_records_the_naming_reference(tmp_path: Path) -> None:
    """D5-09: the standard is a file, its digest is recorded, and it is genomes[0].

    The run does not ask the user to repeat the reference on the command line, and it does not
    let an unauditable record become the standard: the audit runs in front of it, so a reference
    whose CDS carry no ``/translation`` is refused before any name depends on it.
    """
    import json

    from plastid.cli import EXIT_OK, main
    from plastid.core.genbank import to_genbank

    sample = tmp_path / "sample.gb"
    sample.write_text(
        to_genbank(_without_translations(read_genbank(HERMONA_PGA))), encoding="utf-8"
    )
    out = tmp_path / "runs"

    code = main(["run", "-m", "gene_normalize", "--input", str(sample), "--out", str(out)])

    assert code == EXIT_OK
    run_dir = next(path.parent for path in out.glob("*/*/manifest.json"))
    manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["schema_version"] == 2
    assert manifest["inputs"][0]["path"].endswith("NC_050833.1.gb")  # the reference came first
    entry = manifest["modules"][0]
    assert entry["reference"]["accession"] == "NC_050833.1"
    assert entry["reference"]["sha256"] == load_recorded_reference().sha256
    assert entry["reference"]["provisional"] is True
    assert manifest["params"]["naming_reference"].endswith("NC_050833.1.gb")


def test_a_conflict_stops_the_run_with_its_own_exit_code(tmp_path: Path) -> None:
    """Exit 4, not 2: the run was well-formed and the data disagrees with the reference."""
    from plastid.cli import EXIT_CONFLICT, main

    published = read_genbank(HERMONA_PUBLISHED)
    record = _without_translations(read_genbank(HERMONA_PGA))
    duplicated = replace(
        record,
        genes=record.genes + (Gene(name="pbf1", kind="CDS", parts=(Interval(50_000, 50_500),)),),
    )
    from plastid.core.genbank import to_genbank

    path = tmp_path / "planted.gb"
    path.write_text(to_genbank(duplicated), encoding="utf-8")
    del published

    code = main(
        ["run", "-m", "gene_normalize", "--input", str(path), "--out", str(tmp_path / "runs")]
    )

    assert code == EXIT_CONFLICT
    run_dir = next(path.parent for path in (tmp_path / "runs").glob("*/*/manifest.json"))
    import json

    manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["status"] == "failed"
    assert "conflict" in manifest["modules"][0]["error"]
    assert "(50001..50500+)" in manifest["modules"][0]["error"]  # the coordinates it names


def test_a_reference_that_cannot_be_used_is_refused_by_its_audit(tmp_path: Path) -> None:
    """D5-08's gate, at the one place where a bad reference makes every name wrong the same way."""
    from plastid.cli import EXIT_ERROR, main
    from plastid.core.genbank import to_genbank

    bad = tmp_path / "bad_reference.gb"
    record = _without_translations(read_genbank(HERMONA_PUBLISHED))
    bad.write_text(to_genbank(record), encoding="utf-8")  # 86 CDS, no /translation
    # the *record* has no proteins at all, so the audit's translation rule fires
    sample = tmp_path / "sample.gb"
    sample.write_bytes(HERMONA_PGA.read_bytes())

    code = main(
        [
            "run",
            "-m",
            "gene_normalize",
            "--input",
            str(sample),
            "--set",
            f"naming_reference={bad}",
            "--out",
            str(tmp_path / "runs"),
        ]
    )

    assert code == EXIT_ERROR
    assert not (tmp_path / "runs").exists()  # nothing was written at all


# --- The release's Docker check found these (11/10/2026) --------------------------------
#
# `gene_normalize` could not run from an **installed** package: resolving the recorded naming
# reference went through `repo_root()`, which walks up looking for `pyproject.toml` — and a wheel
# in `site-packages` has none above it, so the run died with "no pyproject.toml above
# …/site-packages/plastid/core/gene_names.py". WSL and CI both run from a checkout, so no test
# here could see it; the container could. Both halves are now permanent checks.


def test_checkout_root_is_none_outside_a_checkout(tmp_path: Path) -> None:
    """No checkout above the module is an answer, not a failure — the installed-package case."""
    assert checkout_root(tmp_path) is None


def test_the_recorded_reference_resolves_without_a_checkout(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Without a checkout root the recorded path is taken as it comes — absolute or from the cwd."""
    from plastid.cli import _resolve_recorded_reference

    reference = tmp_path / "reference.gb"
    reference.write_text("LOCUS       synthetic\n", encoding="utf-8")
    recorded = RecordedReference(
        path=reference.name,  # as recorded: relative, and not to any checkout
        accession="synthetic",
        sha256="0" * 64,
        provisional=True,
        recorded="2026-10-11",
        recorded_by="the dev",
        why="a regression test for the installed-package path",
    )
    monkeypatch.setattr("plastid.cli.checkout_root", lambda *args, **kwargs: None)
    monkeypatch.chdir(tmp_path)

    assert _resolve_recorded_reference(recorded) == reference
