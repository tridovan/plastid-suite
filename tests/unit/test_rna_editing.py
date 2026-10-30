"""The PREP-Cp-compatible predictor without its aligner (M4.5): the rule, the mapping, reference B′.

Everything here runs without MAFFT. The scoring rule is tested on columns written by hand; the
locus loop is tested through ``evaluate_genome`` with an aligner the test supplies (a query of the
reference's width is placed column for column), so the mapping cases — both strands, a join, an
origin crossing, IR copies, duplicate names, ambiguity codes, missing genes — are checked against
coordinates computed by hand. The real aligner, the archived *Vitis* outputs and the *Aegilops*
acceptance set are in ``tests/golden/test_rna_editing_acceptance.py``.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping, Sequence
from fractions import Fraction
from pathlib import Path

import pytest

from plastid.core.codon import STANDARD_CODE
from plastid.core.genbank import read_genbank
from plastid.core.gene_names import load_gene_aliases, normalize_genome
from plastid.core.model import Gene, Genome, Interval, coding_sequence
from plastid.core.resources import data_file
from plastid.core.rna_editing import (
    STATUS_ABSENT,
    STATUS_EVALUATED,
    STATUS_NOT_EVALUATED,
    STATUS_NOT_SUPPORTED,
    SUPPORTED_GENES,
    EditingReference,
    EditingReferenceError,
    EditSite,
    LocusOutcome,
    ReferenceRow,
    apply_taxon_names,
    build_reference,
    call_codon,
    check_site,
    codon_states,
    column_counts,
    evaluate_genome,
    load_build_inputs,
    load_reference,
    location_text,
    locus_segments,
    parse_cutoff,
    parse_location_text,
    read_edit_sites,
    reference_digest,
    reference_from_payload,
    supported_gene,
    transcript_positions,
)
from plastid.core.sequence import reverse_complement
from plastid.run.runner import _params_hash

ROOT = Path(__file__).resolve().parents[2]
RECORDS = ROOT / "tests" / "data" / "rna_editing"
# The seven B′ records are pinned in taxa.yaml and kept out of git: `plastid panel pin
# --accessions <the seven> --out .panel/rna_editing` re-creates them byte for byte.
B_PRIME_RECORDS = ROOT / ".panel" / "rna_editing"
needs_b_prime_records = pytest.mark.skipif(
    not B_PRIME_RECORDS.is_dir(),
    reason="the seven B′ records are not in this working tree (.panel/ is gitignored): "
    "plastid panel pin --accessions … --out .panel/rna_editing",
)

# --------------------------------------------------------------------------- the scoring rule


def test_a_codon_lists_its_unedited_state_and_every_c_to_u_reading() -> None:
    states = {state.codon: state.edited for state in codon_states("CCG")}
    assert states == {"CCG": (), "TCG": (1,), "CTG": (2,), "TTG": (1, 2)}
    assert [state.codon for state in codon_states("GGA")] == ["GGA"]


def test_gaps_are_not_in_the_denominator() -> None:
    """The archived accD output mixes fifths and quarters in one gene: N is per column."""
    counts, compared = column_counts("LLLL-")
    assert (counts, compared) == ({"L": 4}, 4)
    five = call_codon("TCG", "LLLLS")
    four = call_codon("CCT", "LLL-F")
    assert five is not None and five.score == Fraction(4, 5)
    assert four is not None and four.score == Fraction(3, 4)


def test_a_column_with_no_reference_residue_is_not_scored() -> None:
    assert call_codon("TCA", "---") is None


def test_the_highest_score_wins_and_a_tie_goes_to_fewer_edits() -> None:
    # CCC: P unedited, S (pos 1), L (pos 2), P (pos 3, silent). P and S tie at 1/2.
    call = call_codon("CCC", "PS")
    assert call is not None and call.state.edited == () and not call.edited


def test_a_tie_between_single_edits_goes_to_codon_position_two() -> None:
    # CCA: TCA (S, pos 1) and CTA (L, pos 2) tie at 1/2; the unedited P scores 0.
    call = call_codon("CCA", "SL")
    assert call is not None and call.state.codon == "CTA" and call.state.edited == (2,)


def test_a_silent_edit_is_never_predicted() -> None:
    call = call_codon("TTC", "FFFF")
    assert call is not None and call.state.codon == "TTC" and not call.edited


def test_an_edit_at_codon_position_three_never_changes_the_amino_acid() -> None:
    """Why the papers' open tie (single edits at positions 1 and 3) cannot change a call.

    Under the standard code NNC and NNU always encode one amino acid, so a state carrying a
    position-3 edit ties with the same state without it, which has fewer edits and wins. Every
    codon is checked, so a change to the code table would fail here first.
    """
    for codon in STANDARD_CODE:
        states = {state.edited: state.amino_acid for state in codon_states(codon)}
        for edited, amino_acid in states.items():
            if 3 in edited:
                without = tuple(position for position in edited if position != 3)
                assert states[without] == amino_acid, (codon, edited)


def test_a_stop_state_cannot_win_against_reference_residues() -> None:
    call = call_codon("CAA", "QQ*")
    assert call is not None and call.state.amino_acid == "Q"


def test_the_cutoff_is_exact_so_a_four_fifths_site_is_reported_at_c_0_8() -> None:
    """``Fraction(0.8)`` is larger than 4/5; parsing the decimal text keeps the 0.80 site."""
    assert Fraction(0.8) > Fraction(4, 5)
    assert parse_cutoff(0.8) == Fraction(4, 5)
    assert parse_cutoff("0.8") == Fraction(4, 5)
    assert parse_cutoff(1) == 1
    call = call_codon("TCG", "LLLLS")
    assert call is not None and call.score >= parse_cutoff("0.8")


@pytest.mark.parametrize("value", ["1.5", "-0.1", "abc", "nan", "inf"])
def test_a_cutoff_outside_zero_to_one_is_refused(value: str) -> None:
    with pytest.raises(ValueError, match="between 0 and 1"):
        parse_cutoff(value)


# --------------------------------------------------------------------------- the query loop

# Query: ATG TCA CCA GGG CAT TAA → M S P G H. Against the reference below: codon 2 S→L scores
# 3/3, codon 3 P→L 2/3, codon 5 H→Y 2/2 (one gap). At C = 0.8 the sites are CDS 5 and CDS 13.
CDS = "ATGTCACCAGGGCATTAA"
REFERENCE_ROWS = ("MLLGY", "MLLGY", "MLPG-")


def _reference(gene: str = "psbL", rows: Sequence[str] = REFERENCE_ROWS) -> EditingReference:
    return EditingReference(
        label="test reference",
        digest="0" * 64,
        genes={
            gene: tuple(
                ReferenceRow(taxon=f"Taxon {index}", accession=f"T{index}", aligned=row)
                for index, row in enumerate(rows)
            )
        },
        payload={},
    )


def _columnwise(rows: Sequence[tuple[str, str]], query: tuple[str, str]) -> Mapping[str, str]:
    """An aligner for a query of the reference's width: column i holds residue i."""
    width = len(rows[0][1])
    if len(query[1]) != width:
        raise AssertionError("the test aligner needs a query of the reference's width")
    return {**dict(rows), query[0]: query[1]}


def _genome(length: int, features: Sequence[tuple[str, Sequence[Interval], str]]) -> Genome:
    """A record whose bases are written by hand: each feature's forward-strand text per part."""
    bases = ["A"] * length
    genes = []
    for name, parts, forward in features:
        cursor = 0
        for part in sorted(parts, key=lambda part: part.start):
            for offset in range(part.length):
                bases[part.start + offset] = forward[cursor + offset]
            cursor += part.length
        genes.append(Gene(name=name, kind="CDS", parts=tuple(parts), transl_table=11))
    return Genome(name="test", sequence="".join(bases), accession="TEST.1", genes=tuple(genes))


def _evaluate(genome: Genome, cutoff: str = "0.8", **kwargs) -> dict[str, LocusOutcome]:
    outcomes = evaluate_genome(
        genome,
        kwargs.get("reference", _reference()),
        parse_cutoff(cutoff),
        load_gene_aliases(),
        kwargs.get("add", _columnwise),
    )
    return {outcome.locus: outcome for outcome in outcomes}


def _sites(outcome: LocusOutcome) -> list[tuple[int, int, int]]:
    """(CDS position, 1-based genome position, strand) of every predicted site."""
    assert outcome.prediction is not None and outcome.prepared is not None
    return [
        (
            site.cds_position,
            outcome.prepared.positions[site.cds_position - 1][0] + 1,
            outcome.prepared.positions[site.cds_position - 1][1],
        )
        for site in outcome.prediction.sites
    ]


def test_plus_strand_sites_are_mapped_to_their_genome_bases() -> None:
    genome = _genome(300, [("psbL", [Interval(10, 28, 1)], CDS)])
    locus = _evaluate(genome)["psbL"]
    assert coding_sequence(genome, locus.gene) == CDS
    assert _sites(locus) == [(5, 15, 1), (13, 23, 1)]
    for cds_position, genome_position, _ in _sites(locus):
        assert CDS[cds_position - 1] == "C" == genome.sequence[genome_position - 1]


def test_a_minus_strand_c_is_a_forward_strand_g() -> None:
    genome = _genome(300, [("psbL", [Interval(100, 118, -1)], reverse_complement(CDS))])
    locus = _evaluate(genome)["psbL"]
    assert coding_sequence(genome, locus.gene) == CDS
    # transcript base i sits at forward position 117 - i (0-based)
    assert _sites(locus) == [(5, 114, -1), (13, 106, -1)]
    assert {genome.sequence[position - 1] for _, position, _ in _sites(locus)} == {"G"}
    assert locus_segments(locus.prepared.positions) == "1-18:118-101(-)"


def test_a_joined_minus_strand_cds_maps_each_exon_in_transcript_order() -> None:
    exon_one, exon_two = Interval(229, 238, -1), Interval(200, 209, -1)
    forward = reverse_complement(CDS[9:]) + reverse_complement(CDS[:9])  # exon two lies first
    genome = _genome(300, [("psbL", [exon_one, exon_two], forward)])
    locus = _evaluate(genome)["psbL"]
    assert coding_sequence(genome, locus.gene) == CDS
    assert _sites(locus) == [(5, 234, -1), (13, 206, -1)]
    assert locus_segments(locus.prepared.positions) == "1-9:238-230(-);10-18:209-201(-)"


def test_an_origin_crossing_cds_wraps_through_position_one() -> None:
    parts = [Interval(293, 300, 1), Interval(0, 11, 1)]
    genome = _genome(300, [("psbL", parts, CDS[7:] + CDS[:7])])
    locus = _evaluate(genome, cutoff="0.5")["psbL"]
    assert coding_sequence(genome, locus.gene) == CDS
    # CDS 8 (codon 3, P→L at 2/3) is the first base after the origin
    assert _sites(locus) == [(5, 298, 1), (8, 1, 1), (13, 6, 1)]
    assert locus_segments(locus.prepared.positions) == "1-7:294-300(+);8-18:1-11(+)"


def test_each_ir_copy_is_its_own_locus_on_its_own_coordinates() -> None:
    genome = _genome(
        300,
        [
            ("psbL", [Interval(10, 28, 1)], CDS),
            ("psbL", [Interval(100, 118, -1)], reverse_complement(CDS)),
        ],
    )
    loci = _evaluate(genome)
    assert {"psbL#1", "psbL#2"} <= set(loci)
    assert _sites(loci["psbL#1"]) == [(5, 15, 1), (13, 23, 1)]
    assert _sites(loci["psbL#2"]) == [(5, 114, -1), (13, 106, -1)]


def test_a_synonym_reaches_the_gene_it_names() -> None:
    """``pafI`` is ``ycf3``'s other name (data/gene_aliases.yaml): vocabulary A, not spelling."""
    genome = _genome(300, [("pafI", [Interval(10, 28, 1)], CDS)])
    loci = _evaluate(genome, reference=_reference("ycf3"))
    locus = loci["ycf3"]
    assert locus.status == STATUS_EVALUATED and locus.prep_gene == "ycf3"
    assert locus.name_in_record == "pafI"


def test_not_supported_absent_and_not_evaluated_are_not_zero_sites() -> None:
    genome = _genome(
        400,
        [
            ("psbL", [Interval(10, 28, 1)], CDS),
            ("rbcL", [Interval(50, 68, 1)], CDS),
        ],
    )
    pseudo = Gene(name="psbF", kind="CDS", parts=(Interval(100, 118, 1),), pseudo=True)
    genome = Genome(
        name=genome.name,
        sequence=genome.sequence,
        accession=genome.accession,
        genes=(*genome.genes, pseudo),
    )
    loci = _evaluate(genome)
    assert loci["psbL"].status == STATUS_EVALUATED
    assert len(loci["psbL"].prediction.sites) == 2
    assert loci["rbcL"].status == STATUS_NOT_SUPPORTED and loci["rbcL"].prediction is None
    assert loci["psbF"].status == STATUS_NOT_EVALUATED
    assert loci["psbF"].reason == "pseudogene (/pseudo)" and loci["psbF"].prediction is None
    absent = {locus for locus, outcome in loci.items() if outcome.status == STATUS_ABSENT}
    assert absent == set(SUPPORTED_GENES) - {"psbL", "psbF"}
    assert all(loci[gene].prediction is None for gene in absent)


def test_a_supported_gene_with_no_reference_protein_is_not_evaluated_by_name() -> None:
    genome = _genome(300, [("psbE", [Interval(10, 28, 1)], CDS)])
    locus = _evaluate(genome)["psbE"]
    assert locus.status == STATUS_NOT_EVALUATED
    assert locus.reason == "reference B′ holds no psbE protein"


def test_an_ambiguity_code_is_listed_and_the_rest_of_the_locus_is_scored() -> None:
    ambiguous = CDS[:9] + "GNG" + CDS[12:]
    genome = _genome(300, [("psbL", [Interval(10, 28, 1)], ambiguous)])
    locus = _evaluate(genome)["psbL"]
    assert locus.status == STATUS_EVALUATED
    assert locus.prediction.ambiguous == (4,)
    assert [site.cds_position for site in locus.prediction.sites] == [5, 13]


@pytest.mark.parametrize(
    ("gene", "reason"),
    [
        (
            Gene(name="psbL", kind="CDS", parts=(Interval(10, 27, 1),), transl_table=11),
            "frame: 17 nt from /codon_start is not whole codons",
        ),
        (
            Gene(name="psbL", kind="CDS", parts=(Interval(10, 28, 1),), transl_table=4),
            "translation table 4 (only the standard code is scored)",
        ),
    ],
)
def test_a_locus_the_rule_cannot_read_is_not_evaluated_with_its_reason(
    gene: Gene, reason: str
) -> None:
    base = _genome(300, [("psbL", [Interval(10, 28, 1)], CDS)])
    genome = Genome(name="t", sequence=base.sequence, accession="T.1", genes=(gene,))
    assert _evaluate(genome)["psbL"].reason == reason


def test_an_internal_stop_is_not_evaluated() -> None:
    stopped = CDS[:6] + "TGA" + CDS[9:]
    genome = _genome(300, [("psbL", [Interval(10, 28, 1)], stopped)])
    assert _evaluate(genome)["psbL"].reason == "internal stop at codon 3"


def test_an_aligner_that_moves_the_reference_is_an_alignment_failure() -> None:
    def shifting(rows, query):
        aligned = dict(_columnwise(rows, query))
        first = rows[0][0]
        aligned[first] = aligned[first][1:] + aligned[first][0]
        return aligned

    genome = _genome(300, [("psbL", [Interval(10, 28, 1)], CDS)])
    locus = _evaluate(genome, add=shifting)["psbL"]
    assert locus.status == STATUS_NOT_EVALUATED
    assert locus.reason.startswith("alignment failure: the aligner changed the reference row")


def test_transcript_positions_reproduce_the_coding_sequence(golden_genome: Genome) -> None:
    """On a real record, every mapped base read on its strand gives ``coding_sequence`` back."""
    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    for gene in golden_genome.genes:
        if gene.kind != "CDS":
            continue
        text = "".join(
            golden_genome.sequence[position]
            if strand > 0
            else complement.get(golden_genome.sequence[position], "N")
            for position, strand in transcript_positions(golden_genome, gene)
        )
        assert text == coding_sequence(golden_genome, gene).upper(), gene.name


# --------------------------------------------------------------------------- reference B′


@pytest.fixture(scope="module")
def packaged() -> EditingReference:
    return load_reference()


def test_the_packaged_reference_matches_its_digest(packaged: EditingReference) -> None:
    assert packaged.digest == reference_digest(packaged.payload)
    assert set(packaged.genes) == set(SUPPORTED_GENES)
    assert "not the original PREP-Cp reference" in packaged.label


def test_the_packaged_reference_is_built_from_the_packaged_inputs(
    packaged: EditingReference,
) -> None:
    """A changed edit list, taxon pin or source without a rebuild fails here — a stale B′."""
    taxa, sources = load_build_inputs(data_file("rna_editing"))
    assert packaged.payload["taxa"] == taxa
    assert packaged.payload["sources"] == sources
    listed = [
        {name: getattr(site, name) for name in site.__dataclass_fields__}
        for site in read_edit_sites(data_file("rna_editing/edit_sites.tsv"))
    ]
    key = lambda site: (site["accession"], site["gene"], site["cds_position"])  # noqa: E731
    assert sorted(packaged.payload["edit_sites"], key=key) == sorted(listed, key=key)


def _records() -> dict[str, Genome]:
    """The pinned B′ records with the names and models taxa.yaml states, as the build reads them."""
    taxa, _ = load_build_inputs(data_file("rna_editing"))
    return {
        taxon["accession"]: apply_taxon_names(
            read_genbank(B_PRIME_RECORDS / f"{taxon['accession']}.gb"),
            taxon.get("names") or (),
            taxon["accession"],
        )
        for taxon in taxa
    }


@pytest.mark.panel
@needs_b_prime_records
def test_every_pinned_record_is_the_bytes_its_pin_names(packaged: EditingReference) -> None:
    for taxon in packaged.payload["taxa"]:
        path = B_PRIME_RECORDS / f"{taxon['accession']}.gb"
        assert hashlib.sha256(path.read_bytes()).hexdigest() == taxon["sha256"], path.name


@pytest.mark.panel
@needs_b_prime_records
def test_every_mature_protein_is_its_record_with_its_sourced_edits(
    packaged: EditingReference,
) -> None:
    """Each B′ row is recomputed from the pinned record and the edit list, without an aligner.

    The row's own location is the model used — the record's, or the corrected one taxa.yaml
    states — and a stop created by a sourced edit ends the mature protein there.
    """
    records = _records()
    sites = read_edit_sites(data_file("rna_editing/edit_sites.tsv"))
    by_gene: dict[tuple[str, str], list[int]] = {}
    for site in sites:
        by_gene.setdefault((site.accession, site.gene), []).append(site.cds_position)
    for gene_name, entry in packaged.payload["genes"].items():
        for row in entry["rows"]:
            assert row["aligned"].replace("-", "") == row["protein"]
            assert row["edits"] == sorted(by_gene.get((row["accession"], gene_name), []))
            genome = records[row["accession"]]
            (model,) = [
                gene
                for gene in genome.genes
                if gene.kind == "CDS" and location_text(gene) == row["location"]
            ]
            cds = list(coding_sequence(genome, model).upper())
            assert hashlib.sha256("".join(cds).encode()).hexdigest() == row["cds_sha256"]
            for position in row["edits"]:
                assert cds[position - 1] == "C"
                cds[position - 1] = "T"
            body = "".join(cds)[model.codon_start - 1 :]
            protein = "".join(
                STANDARD_CODE[body[index : index + 3]] for index in range(0, len(body) - 2, 3)
            )
            mature = protein.split("*", 1)[0]
            if "*" in protein.rstrip("*"):
                assert row["stop_created_at"] == len(mature) + 1, (gene_name, row["accession"])
            assert mature == row["protein"], (gene_name, row["accession"])


@pytest.mark.panel
@needs_b_prime_records
def test_every_sourced_site_fits_its_pinned_record() -> None:
    """The base is a C on the transcript strand, at the codon and amino acid the source states."""
    records = _records()
    aliases = load_gene_aliases()
    for site in read_edit_sites(data_file("rna_editing/edit_sites.tsv")):
        genome = records[site.accession]
        names = normalize_genome(genome, aliases)
        candidates = [
            gene
            for gene in genome.genes
            if gene.kind == "CDS"
            and not gene.pseudo
            and supported_gene(names.normalized(gene)) == site.gene
        ]
        assert candidates, (site.accession, site.gene)
        assert check_site(genome, candidates[0], site) is None, (site.accession, site.gene)


def test_any_change_to_reference_b_prime_changes_its_digest(packaged: EditingReference) -> None:
    payload = json.loads(json.dumps(packaged.payload))
    original = reference_digest(payload)
    for change in (
        lambda p: p["taxa"][0].__setitem__("sha256", "0" * 64),
        lambda p: p["edit_sites"][0].__setitem__("cds_position", 1),
        lambda p: p["recipe"]["aligner"].__setitem__("version", "7.999"),
        lambda p: p["recipe"]["aligner"].__setitem__("argv", ["--auto"]),
        lambda p: p["genes"]["accD"]["rows"][0].__setitem__("aligned", "-"),
    ):
        changed = json.loads(json.dumps(payload))
        change(changed)
        assert reference_digest(changed) != original


def test_a_reference_edited_after_it_was_built_is_refused(packaged: EditingReference) -> None:
    payload = json.loads(json.dumps(packaged.payload))
    payload["genes"]["accD"]["rows"][0]["aligned"] = payload["genes"]["accD"]["rows"][0][
        "aligned"
    ].replace("M", "L", 1)
    with pytest.raises(EditingReferenceError, match="does not match its contents"):
        reference_from_payload(payload)


def test_the_reference_digest_is_part_of_the_cache_key() -> None:
    """The CLI passes B′'s digest to the step; another digest is another key (no stale reuse)."""
    one = _params_hash({"cutoff": 0.8, "rna_editing_reference_digest": "a" * 64})
    other = _params_hash({"cutoff": 0.8, "rna_editing_reference_digest": "b" * 64})
    assert one != other


def test_a_missing_mafft_is_a_named_error_not_a_fallback(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    """With no aligner there is no prediction: exit 3 and a message naming MAFFT."""
    from plastid.cli import EXIT_TOOL, main

    monkeypatch.setenv("PATH", str(tmp_path))  # a PATH with nothing on it
    code = main(
        ["run", "-m", "rna_editing", "-i", str(RECORDS / "NC_007957.1.gb"), "-o", str(tmp_path)]
    )
    assert code == EXIT_TOOL
    assert "mafft not found — the RNA-editing predictor needs MAFFT 7.526" in (
        capsys.readouterr().err
    )


def _build_one(genome: Genome, sites: Sequence[EditSite], names: Sequence[dict] = ()) -> dict:
    """B′ built from one synthetic taxon (one row per gene, so no aligner call is made)."""
    taxon = {
        "key": "Test",
        "organism": "Testus testus",
        "accession": genome.accession,
        "sha256": "0" * 64,
        "names": list(names),
    }
    return build_reference(
        records={genome.accession: genome},
        record_sha256={genome.accession: "0" * 64},
        taxa=[taxon],
        sites=sites,
        sources={"test": {}},
        aliases=load_gene_aliases(),
        aligner=lambda rows: {name: protein for name, protein in rows},
        aligner_record={"name": "none"},
        label="test",
    )


def _site(gene: str, cds_position: int, codon: str, after: str) -> EditSite:
    number, within = divmod(cds_position - 1, 3)
    return EditSite(
        accession="TEST.1",
        gene=gene,
        cds_position=cds_position,
        codon_position=within + 1,
        aa_position=number + 1,
        codon_before=codon,
        codon_after=after,
        aa_before=STANDARD_CODE[codon],
        aa_after=STANDARD_CODE[after],
        source="test",
        locator="hand-built",
        source_coordinate="",
    )


def test_a_stop_created_by_a_sourced_edit_ends_the_mature_protein() -> None:
    """Black pine petD: an edited CAA makes the 160-residue protein of the other taxa."""
    cds = "ATGTCACAAGGGCATTAA"  # M S Q G H stop; CAA (codon 3) edited to UAA
    genome = _genome(300, [("psbL", [Interval(10, 28, 1)], cds)])
    payload = _build_one(genome, [_site("psbL", 5, "TCA", "TTA"), _site("psbL", 7, "CAA", "TAA")])
    (row,) = payload["genes"]["psbL"]["rows"]
    assert row["protein"] == "ML" and row["stop_created_at"] == 3


def test_an_internal_stop_no_source_created_excludes_the_taxon() -> None:
    genome = _genome(300, [("psbL", [Interval(10, 28, 1)], "ATGTCATAAGGGCATTAA")])
    payload = _build_one(genome, [])
    assert "psbL" not in payload["genes"]


def test_taxa_yaml_can_name_a_product_only_cds_and_correct_its_model() -> None:
    """The record's model starts one codon late (as black pine petL does); the entry fixes it."""
    cds = "ACGTCAGGGCATTAA"  # T S G H stop; the ACG start is edited to AUG
    genome = _genome(300, [("PSI subunit", [Interval(13, 25, 1)], cds[3:])])
    bases = list(genome.sequence)
    bases[10:13] = "ACG"
    genome = Genome(
        name=genome.name,
        sequence="".join(bases),
        accession=genome.accession,
        genes=tuple(
            Gene(name=gene.name, kind="CDS", parts=gene.parts, product="PSI subunit")
            for gene in genome.genes
        ),
    )
    entry = {
        "location": "14..25(+)",
        "product": "PSI subunit",
        "gene": "psbL",
        "cds": "11..25(+)",
        "evidence": "test",
    }
    assert parse_location_text("11..25(+),1..3(-)") == (
        Interval(10, 25, 1),
        Interval(0, 3, -1),
    )
    payload = _build_one(genome, [_site("psbL", 2, "ACG", "ATG")], names=[entry])
    (row,) = payload["genes"]["psbL"]["rows"]
    assert row["location"] == "11..25(+)" and row["protein"] == "MSGH"
    with pytest.raises(EditingReferenceError, match="matches 0 CDS"):
        _build_one(genome, [], names=[{**entry, "location": "15..25(+)"}])
