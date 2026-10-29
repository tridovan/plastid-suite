"""One record per family, every module that reads genes or sequences.

A case on the panel's families is different from the golden pair and the Iridaceae set. This file
fixes one record per family — the first source record whose own two IR copies pass the gate at
100.00 % identity — and asserts, per family, that the module accounts for the whole record rather
than quietly measuring part of it.

The accessions are chosen by the gate, not by convenience; the premise test re-runs the gate so a
panel refresh that invalidates a choice fails here instead of inside a scientific number.

The assertions are semantic on purpose. A pinned count would say "this record has 84 CDS" and
stop; the claim worth holding is "every annotated CDS is in the table exactly once, and a CDS that
fails the ORF check is reported with a status rather than dropped" — which is what a reader of the
output depends on, in every family.

The four multi-genome modules (`kaks`, `popgen`, `synteny`, `phylogeny`) run **once over the ten
records** and then make one claim per family, so a family that breaks a module is named by the
failing test id. Each claim is the module's own contract: `kaks` gives every gene a row and never
divides by a zero `Ks`; `popgen` publishes the `n` behind every statistic and leaves a cell blank
with its reason when the statistic does not exist; `synteny` reports coverage against the
reference's **normalized** frame; `phylogeny` puts every input family in one supported tree.
"""

from __future__ import annotations

import io
import os
from collections import Counter
from pathlib import Path

import pytest
from Bio import Phylo, SeqIO

from conftest import require_tools
from plastid.core.genbank import read_genbank
from plastid.core.gene_names import cds_by_canonical_name
from plastid.core.ir import detect_ir, enforce_ir_gate
from plastid.core.model import Genome, count_genes
from plastid.core.normalize import normalize
from plastid.core.orf import FAILURE_STATUSES, ORF_STATUSES
from plastid.core.popgen import MIN_SEQUENCES_FOR_NEUTRALITY
from plastid.modules import gene_table, kaks, orf, phylogeny, popgen, synteny
from plastid.run.context import RunContext
from plastid.run.layout import RunLayout

ROOT = Path(__file__).resolve().parents[2]
PANEL = ROOT / ".panel" / "ir-validation" / "source"

# Reference first: the multi-genome modules treat the first input as the reference, and
# NC_008641.1 is the one panel record the IR gene-mirror rule passes with no violation.
FAMILY_RECORDS = {
    "Malvaceae": "NC_008641.1",
    "Asteraceae": "NC_028027.1",
    "Brassicaceae": "NC_029378.1",
    "Euphorbiaceae": "NC_016736.1",
    "Fabaceae": "NC_016708.2",
    "Lamiaceae": "NC_020098.1",
    "Orchidaceae": "NC_021432.1",
    "Poaceae": "NC_023245.1",
    "Rosaceae": "NC_031163.1",
    "Solanaceae": "NC_016068.1",
}

# Every multi-genome module below takes the Malvaceae record as its reference; the rest are its
# queries or, for `phylogeny`, the other taxa of one tree.
REFERENCE_FAMILY = next(iter(FAMILY_RECORDS))
QUERY_FAMILIES = tuple(FAMILY_RECORDS)[1:]

# The ten records split along a real botanical line: two monocots and eight dicots. `popgen`
# compares taxa, not populations, and this split is what makes the case carry both branches of
# its own rule — a group of two, where the neutrality tests are undefined by arithmetic, and a
# group of eight, where they have a value.
MONOCOT_FAMILIES = frozenset({"Orchidaceae", "Poaceae"})

# Every status `kaks` can write; a row with any other label is a reporting bug.
KAKS_STATUSES = frozenset(
    {
        kaks.STATUS_OK,
        kaks.STATUS_KS_ZERO,
        kaks.STATUS_ABSENT_REFERENCE,
        kaks.STATUS_ABSENT_GENOME,
        kaks.STATUS_NOT_WHOLE_CODONS,
        kaks.STATUS_CODON_START,
        kaks.STATUS_INTERNAL_STOP,
        kaks.STATUS_FRAMESHIFT,
        kaks.STATUS_TOOL_FAILED,
    }
)

pytestmark = [
    pytest.mark.tier_b,
    pytest.mark.skipif(
        not PANEL.is_dir(),
        reason="the IR validation panel is not in this tree (.panel/ is gitignored)",
    ),
]


def _genomes() -> dict[str, Genome]:
    missing = [f"{a}.gb" for a in FAMILY_RECORDS.values() if not (PANEL / f"{a}.gb").is_file()]
    if missing:
        pytest.skip(f"panel records missing: {', '.join(missing)}")
    return {
        family: read_genbank(PANEL / f"{accession}.gb")
        for family, accession in FAMILY_RECORDS.items()
    }


def _ctx(root: Path, module: str, genome: Genome | None = None) -> RunContext:
    layout = RunLayout.create(root, module, fingerprint="family")
    layout.make_dirs()
    ctx = RunContext(layout=layout, params={}, module_id=module, module_order=1)
    if genome is not None:
        ctx.ir = enforce_ir_gate(detect_ir(genome, source="sequence"))
    return ctx


def _iqtree_dir() -> Path:
    """IQ-TREE's directory, from `PTS_IQTREE_DIR` or the shared toolbox layout."""
    directory = os.environ.get("PTS_IQTREE_DIR")
    candidates = [Path(directory)] if directory else []
    candidates.append(Path("/root/toolbox/env/iqtree/bin"))
    for candidate in candidates:
        if (candidate / "iqtree2").is_file() or (candidate / "iqtree").is_file():
            return candidate
    reason = "iqtree2 is not installed; set PTS_IQTREE_DIR to the directory holding it"
    if os.environ.get("PTS_REQUIRE_TOOLS") == "1":
        pytest.fail(reason)
    pytest.skip(reason)


def _group_of(family: str) -> str:
    return "monocot" if family in MONOCOT_FAMILIES else "dicot"


def _tree_labels(text: str) -> list[str]:
    return sorted(
        terminal.name for terminal in Phylo.read(io.StringIO(text), "newick").get_terminals()
    )


def _support_values(text: str) -> list[float]:
    """Every support number an IQ-TREE tree carries, as `-bb` and `-alrt` report them.

    IQ-TREE writes two numbers on an internal node (`100/100`), which Biopython keeps as the
    node's label rather than its confidence, so the label is split and both numbers are read.
    """
    tree = Phylo.read(io.StringIO(text), "newick")
    values: list[float] = []
    for clade in tree.get_nonterminals():
        if clade is tree.root:
            continue
        for part in (clade.name or "").split("/"):
            try:
                values.append(float(part))
            except ValueError:
                continue
        if clade.confidence is not None:
            values.append(clade.confidence)
    return values


def test_the_ten_family_records_pass_the_gate_on_their_own_sequence() -> None:
    """The premise: every chosen record carries two copies that agree at the gate's threshold."""
    for family, genome in _genomes().items():
        ir = detect_ir(genome, source="sequence")
        assert ir.gate_passed, (
            f"{family} {genome.accession}: identity {ir.identity:.4f} fails the gate"
        )
        assert ir.identity >= 0.999


def test_the_gene_table_accounts_for_every_gene_in_every_family(tmp_path: Path) -> None:
    """`n_genes_total` is the record's own count, and every gene falls into a group."""
    for family, genome in _genomes().items():
        result = gene_table.run(genome, {}, _ctx(tmp_path / family, "gene_table", genome))
        summary = result.tables["gene_table_summary"].rows[0]
        counted = count_genes(genome)
        expected = counted.cds + counted.trna + counted.rrna
        assert summary["n_genes_total"] == expected, (
            f"{family}: {summary['n_genes_total']} != {expected}"
        )
        assert summary["n_groups"] > 0, f"{family}: no functional group was assigned"

        # The module's own reconciliation is the accounting claim: every feature it prints is
        # counted against `count_genes()` for its kind, and the closing row must agree. A gene no
        # group claims is counted in `n_unclassified` and named in the notes — that rule is held
        # by `tests/golden/test_tier_a.py`; here it must hold for every family.
        reconciled = {row["kind"]: row for row in result.tables["gene_reconciliation"].rows}
        for kind, row in reconciled.items():
            assert row["agreement"] == "agrees", (
                f"{family}: {kind} {row['in_table']} vs {row['counted']}"
            )
        assert reconciled["Total"]["counted"] == expected
        if summary["n_unclassified"]:
            note = next((n for n in result.notes if n.startswith("Unclassified genes")), "")
            assert note, f"{family}: an unclassified gene is counted but never named"


def test_orf_checks_every_annotated_cds_once_in_every_family(tmp_path: Path) -> None:
    """No CDS is skipped, and one that does not translate is reported with a status, not dropped."""
    for family, genome in _genomes().items():
        result = orf.run(genome, {}, _ctx(tmp_path / family, "orf", genome))
        summary = result.tables["orf_summary"].rows[0]
        checks = result.tables["orf_checks"].rows

        assert summary["n_cds"] == len(checks), (
            f"{family}: summary says {summary['n_cds']} CDS, the table holds {len(checks)}"
        )
        statuses = result.tables["orf_status"].rows
        assert sum(row["n"] for row in statuses) == len(checks), (
            f"{family}: status counts do not add up"
        )
        assert {row["status"] for row in statuses} == set(ORF_STATUSES), (
            f"{family}: a status row is missing"
        )
        assert {row["status"] for row in checks} <= set(ORF_STATUSES)

        # The headline class of a CDS that does not translate is a failure the reader is meant to
        # act on, or a pseudogene the record itself declares — never a silently dropped row.
        failures = sum(1 for row in checks if row["status"] in FAILURE_STATUSES)
        assert summary["n_failed"] == failures, (
            f"{family}: {summary['n_failed']} failures, {failures} rows"
        )
        pseudogenes = sum(1 for row in checks if row["status"] == "pseudogene")
        assert summary["n_pseudogene"] == pseudogenes, (
            f"{family}: pseudogene count disagrees with its rows"
        )


# ---------------------------------------------------------------- multi-genome modules
#
# Each of these modules reads several genomes at once, so the case is one run over the ten
# records (a module-scoped fixture) followed by one claim per family. A family that breaks a
# module then shows up as a failing test id instead of disappearing inside one aggregate
# assertion, and every claim is checked against that family's own record rather than a number
# pinned in this file.


@pytest.fixture(scope="module")
def family_genomes() -> dict[str, Genome]:
    return _genomes()


@pytest.fixture(scope="module")
def kaks_family(
    family_genomes: dict[str, Genome],
    tools: dict[str, str | None],
    paml_tools: dict[str, str | None],
    tmp_path_factory: pytest.TempPathFactory,
):
    """One Ka/Ks run: the reference family against the other nine."""
    require_tools(tools, "mafft")
    require_tools(paml_tools, "yn00")
    genomes = list(family_genomes.values())
    ctx = _ctx(tmp_path_factory.mktemp("kaks_family"), "kaks", genomes[0])
    return kaks.run(genomes, {}, ctx)


@pytest.fixture(scope="module")
def popgen_family(
    family_genomes: dict[str, Genome],
    tools: dict[str, str | None],
    tmp_path_factory: pytest.TempPathFactory,
):
    """One population-genetics run: ten taxa in the two sample groups."""
    require_tools(tools, "mafft")
    genomes = list(family_genomes.values())
    groups = {FAMILY_RECORDS[family]: _group_of(family) for family in FAMILY_RECORDS}
    ctx = _ctx(tmp_path_factory.mktemp("popgen_family"), "popgen", genomes[0])
    return popgen.run(genomes, {"groups": groups}, ctx)


@pytest.fixture(scope="module")
def synteny_family(
    family_genomes: dict[str, Genome],
    mummer_tools: dict[str, str | None],
    tmp_path_factory: pytest.TempPathFactory,
):
    """One MUMmer comparison: the reference family against the other nine."""
    require_tools(mummer_tools, "nucmer", "delta-filter", "show-coords")
    genomes = list(family_genomes.values())
    ctx = _ctx(tmp_path_factory.mktemp("synteny_family"), "synteny", genomes[0])
    return synteny.run(genomes, {}, ctx)


@pytest.fixture(scope="module")
def phylogeny_family(
    family_genomes: dict[str, Genome],
    tools: dict[str, str | None],
    tmp_path_factory: pytest.TempPathFactory,
):
    """One supermatrix and one tree over the ten families."""
    require_tools(tools, "mafft")
    genomes = list(family_genomes.values())
    ctx = _ctx(tmp_path_factory.mktemp("phylogeny_family"), "phylogeny")
    return phylogeny.run(genomes, {"iqtree": str(_iqtree_dir())}, ctx)


# ------------------------------------------------------------------ kaks


@pytest.mark.parametrize("family", QUERY_FAMILIES)
def test_kaks_scores_every_gene_of_a_family_and_never_divides_by_zero_ks(
    kaks_family, family_genomes: dict[str, Genome], family: str
) -> None:
    """Every gene has one row, and a dS of zero is an undefined ratio — not a ratio of zero.

    A gene one genome does not carry, an alignment whose gaps shift the frame, a pair yn00
    refused: each keeps its row, its status and a reason. The summary's own count is the run's
    gene universe, so a gene that vanished from the table would break the arithmetic below
    rather than read as "not present".
    """
    genome_key = FAMILY_RECORDS[family]
    summary = kaks_family.tables["kaks_summary"].rows[0]
    rows = [row for row in kaks_family.tables["kaks_genes"].rows if row["genome"] == genome_key]

    assert rows, f"{family}: the run published no gene row for this genome"
    assert len({row["gene"] for row in rows}) == len(rows), f"{family}: a gene has two rows"
    assert len(rows) == summary["n_genes"], (
        f"{family}: {len(rows)} rows for a universe of {summary['n_genes']} gene names — a gene "
        f"was dropped or double-counted"
    )

    for row in rows:
        assert row["status"] in KAKS_STATUSES, f"{family} {row['gene']}: {row['status']!r}"
        if row["status"] == kaks.STATUS_OK:
            assert row["ks"] > 0
            assert row["ka_ks"] == pytest.approx(row["ka"] / row["ks"])
        elif row["status"] == kaks.STATUS_KS_ZERO:
            assert row["ks"] == 0
            assert row["ka"] is not None
            assert row["ka_ks"] is None, (
                f"{family} {row['gene']}: Ks = 0 was divided by instead of leaving Ka/Ks blank"
            )
        else:
            assert (row["ka"], row["ks"], row["ka_ks"]) == (None, None, None), (
                f"{family} {row['gene']}: status {row['status']!r} carries a value"
            )
            assert row["detail"], f"{family} {row['gene']}: no value and no reason"

    pooled = next(
        row for row in kaks_family.tables["kaks_pooled"].rows if row["genome"] == genome_key
    )
    shared = set(cds_by_canonical_name(family_genomes[REFERENCE_FAMILY])) & set(
        cds_by_canonical_name(family_genomes[family])
    )
    assert pooled["n_shared"] == len(shared), (
        f"{family}: {pooled['n_shared']} shared genes, the two records carry {len(shared)}"
    )
    scored = sum(1 for row in rows if row["ka"] is not None and row["ks"] is not None)
    assert pooled["n_computed"] == scored, f"{family}: pooled count disagrees with its rows"
    assert pooled["n_computed"] + pooled["n_not_computed"] == pooled["n_shared"]

    # A run whose every gene was refused would satisfy the row rules above and still measure
    # nothing; the family's case has to carry a scored gene.
    assert scored, f"{family}: no gene produced a value, so the pooled row has nothing to pool"
    assert pooled["pooled_ka"] is not None and pooled["pooled_ks"] is not None
    if pooled["pooled_ks"] == 0:
        assert pooled["pooled_ka_ks"] is None, f"{family}: pooled Ks = 0 was divided by"
    else:
        assert pooled["pooled_ka_ks"] == pytest.approx(pooled["pooled_ka"] / pooled["pooled_ks"])


def test_kaks_takes_the_reference_family_as_the_comparison_base(
    kaks_family, family_genomes
) -> None:
    """The Malvaceae record defines the compared gene names, and the estimator is labelled."""
    summary = kaks_family.tables["kaks_summary"].rows[0]
    reference = family_genomes[REFERENCE_FAMILY]
    assert summary["accession"] == FAMILY_RECORDS[REFERENCE_FAMILY]
    assert summary["organism"] == reference.name
    assert summary["n_genomes"] == len(QUERY_FAMILIES)
    assert summary["estimator"] == kaks.ESTIMATOR
    assert summary["tool_version"], "the estimator's version is not recorded"

    universe = {row["gene"] for row in kaks_family.tables["kaks_genes"].rows}
    assert len(universe) == summary["n_genes"]
    carried: set[str] = set()
    for genome in family_genomes.values():
        carried |= set(cds_by_canonical_name(genome))
    assert universe == carried, "a CDS carried by one of the ten records is not in the run"

    assert any("undefined" in note and "Ks = 0" in note for note in kaks_family.notes), (
        "the rule that Ks = 0 leaves Ka/Ks undefined is not published with the result"
    )


# ------------------------------------------------------------------ popgen


@pytest.mark.parametrize("family", tuple(FAMILY_RECORDS))
def test_popgen_publishes_the_n_behind_every_statistic_and_blanks_what_has_none(
    popgen_family, family_genomes: dict[str, Genome], family: str
) -> None:
    """A statistic is a number with its n, or a blank with its reason — never a zero for both.

    The family's group is the set of records its own statistics were computed over, so the claim
    is checked against the record set: every row's `n` must be the number of that group's
    genomes actually carrying the gene, and at `n < 2` (π, Hd) or `n < 4` (the neutrality tests)
    the cell is `None`, never `0`.
    """
    group = _group_of(family)
    members = [other for other in FAMILY_RECORDS if _group_of(other) == group]
    carriers: Counter[str] = Counter()
    for member in members:
        carriers.update(cds_by_canonical_name(family_genomes[member]).keys())

    all_rows = list(popgen_family.tables["popgen_groups"].rows)
    rows = [row for row in all_rows if row["group"] == group]
    assert rows, f"{family}: the {group} group has no row at all"
    by_gene = {row["gene"]: row for row in rows}
    analysed = {row["gene"] for row in all_rows}

    for name in set(cds_by_canonical_name(family_genomes[family])) & analysed:
        row = by_gene.get(name)
        assert row is not None, f"{family}: {name} is carried by the group but has no row"
        assert row["n"] == carriers[name], (
            f"{family} {name}: n = {row['n']}, the {group} group carries it in {carriers[name]}"
        )

    for row in rows:
        n, sites = row["n"], row["sites"]
        assert sites > 0, f"{family} {row['gene']} ({group}): no usable site after alignment"
        if n < 2:
            for column in ("hd", "pi", "theta_w", "tajimas_d", "fu_li_d_star", "fu_li_f_star"):
                assert row[column] is None, (
                    f"{family} {row['gene']}: {column} has no value at n = {n} and was written "
                    f"as {row[column]!r} instead of blank"
                )
            assert row["note"], f"{family} {row['gene']}: a blank cell without its reason"
        else:
            assert row["hd"] is not None
            assert row["pi"] is not None
            assert row["theta_w"] is not None
            assert 0.0 <= row["segregating"] <= sites

        if n < MIN_SEQUENCES_FOR_NEUTRALITY or row["segregating"] == 0:
            assert row["tajimas_d"] is None
            assert row["fu_li_d_star"] is None
            assert row["fu_li_f_star"] is None
            assert row["note"], f"{family} {row['gene']}: a blank cell without its reason"
        else:
            # At n ≥ 4 and S > 0 the statistic exists; a blank here would be a value the module
            # failed to report rather than one arithmetic refuses.
            assert row["tajimas_d"] is not None, f"{family} {row['gene']}: Tajima's D at n = {n}"

    # Fst between this group and the other: its own two n cells must be the group rows' n, and a
    # blank must carry the reason it is blank.
    for row in popgen_family.tables["popgen_fst"].rows:
        if group not in (row["group_a"], row["group_b"]):
            continue
        for side, label in (("n_a", "group_a"), ("n_b", "group_b")):
            source = next(
                (
                    item
                    for item in all_rows
                    if item["gene"] == row["gene"] and item["group"] == row[label]
                ),
                None,
            )
            assert source is not None, f"{family}: {row['row_id']} names a group with no row"
            assert row[side] == source["n"], f"{family}: {row['row_id']} {side}"
        if row["fst"] is None:
            assert row["note"], f"{family}: {row['row_id']} is blank without a reason"
        else:
            assert 0.0 <= row["fst"] <= 1.0, f"{family}: {row['row_id']} Fst outside [0, 1]"


def test_popgen_states_that_its_groups_are_taxa_and_that_it_measured_every_group(
    popgen_family,
) -> None:
    """The dataset is part of the result: n per group, taxa not populations, one Gst estimator."""
    summary = popgen_family.tables["popgen_summary"].rows[0]
    assert summary["n_genomes"] == len(FAMILY_RECORDS)
    assert summary["n_groups"] == 2
    sizes = dict(part.split(": ") for part in summary["group_sizes"].split(", "))
    assert sizes == {
        "monocot": str(len(MONOCOT_FAMILIES)),
        "dicot": str(len(FAMILY_RECORDS) - len(MONOCOT_FAMILIES)),
    }
    assert summary["n_genes"] > 0 and summary["tool"] == "MAFFT"

    notes = " ".join(popgen_family.notes)
    assert "at least four sequences" in notes, "the neutrality tests' n rule is not published"
    assert "taxa, not populations" in notes, "an Fst between taxa is not labelled as divergence"
    assert "Nei's Gst" in notes, "the estimator is not named"
    assert "taxa" in popgen_family.tables["popgen_fst"].scope.label

    fst_values = [
        row["fst"] for row in popgen_family.tables["popgen_fst"].rows if row["fst"] is not None
    ]
    assert fst_values, "no Fst was computed anywhere — an all-blank table proves nothing"
    tajima_values = [
        row["tajimas_d"]
        for row in popgen_family.tables["popgen_groups"].rows
        if row["tajimas_d"] is not None
    ]
    assert tajima_values, "Tajima's D was never computed, so n ≥ 4 was never exercised"


# ------------------------------------------------------------------ synteny


@pytest.mark.parametrize("family", QUERY_FAMILIES)
def test_synteny_measures_each_family_against_the_normalized_reference(
    synteny_family, family_genomes: dict[str, Genome], family: str
) -> None:
    """Coverage is aligned base pairs over the reference, and every coordinate is in its frame.

    The reference is aligned after `normalize` re-orders it as LSC → IRa → SSC → IRb, so the
    target coordinates of every block and window must lie inside the normalized length; a value
    outside it would mean a table and a series drawn from different references. An unaligned
    window stays blank — not 0 % — exactly as an unmeasured statistic does elsewhere.
    """
    genome_key = FAMILY_RECORDS[family]
    summary = synteny_family.tables["synteny_summary"].rows[0]
    reference_length = summary["reference_length"]
    window, step = summary["window_bp"], summary["step_bp"]

    row = next(
        (
            item
            for item in synteny_family.tables["synteny_genomes"].rows
            if item["genome"] == genome_key
        ),
        None,
    )
    assert row is not None, f"{family}: no comparison row for this genome"
    assert row["aligned_bp"] > 0, f"{family}: no aligned base pair"
    assert row["coverage_pct"] == pytest.approx(100.0 * row["aligned_bp"] / reference_length)
    assert 0.0 < row["aligned_identity"] <= 100.0
    assert row["min_identity"] <= row["mean_identity"] <= row["max_identity"]

    blocks = [
        block
        for block in synteny_family.tables["synteny_blocks"].rows
        if block["genome"] == genome_key
    ]
    assert len(blocks) == row["n_blocks"], f"{family}: block count disagrees with the rows"
    query_length = len(family_genomes[family])
    for block in blocks:
        assert 1 <= block["target_start"] <= block["target_end"] <= reference_length, (
            f"{family} {block['block']}: target outside the normalized reference"
        )
        low, high = sorted((block["query_start"], block["query_end"]))
        assert 1 <= low <= high <= query_length
        assert block["orientation"] in {"forward", "reverse"}
        assert 0.0 < block["identity"] <= 100.0

    series = [
        entry
        for entry in synteny_family.tables["synteny_identity"].rows
        if entry["genome"] == genome_key
    ]
    assert series, f"{family}: no identity series"
    for index, entry in enumerate(series):
        assert entry["window_start"] == 1 + index * step
        assert entry["window_end"] == min(entry["window_start"] + window - 1, reference_length)
        assert 0 <= entry["compared_bp"] <= window
        if entry["compared_bp"] == 0:
            assert entry["identity"] is None, (
                f"{family} window {entry['window_start']}: an unaligned window was written as "
                f"0 % instead of blank"
            )
        else:
            assert 0.0 < entry["identity"] <= 100.0


def test_synteny_states_the_reference_and_its_frame(synteny_family, family_genomes) -> None:
    """The summary names the reference and the length it was aligned against, in the IR's frame."""
    summary = synteny_family.tables["synteny_summary"].rows[0]
    reference = family_genomes[REFERENCE_FAMILY]
    ir = enforce_ir_gate(detect_ir(reference, source="sequence"))
    assert summary["accession"] == FAMILY_RECORDS[REFERENCE_FAMILY]
    assert summary["n_genomes"] == len(QUERY_FAMILIES)
    assert summary["reference_length"] == len(normalize(reference, ir).genome)
    for table_id in ("synteny_summary", "synteny_genomes", "synteny_blocks", "synteny_identity"):
        assert synteny_family.tables[table_id].scope.coord_system == "normalized", table_id


# ------------------------------------------------------------------ phylogeny


@pytest.mark.parametrize("family", tuple(FAMILY_RECORDS))
def test_phylogeny_puts_every_family_in_the_matrix_and_the_tree(
    phylogeny_family, family: str
) -> None:
    """No input genome is dropped: each family is a row of the matrix and a leaf of the tree."""
    key = FAMILY_RECORDS[family]
    summary = phylogeny_family.tables["phylogeny_summary"].rows[0]
    samples = {row["accession"]: row for row in phylogeny_family.tables["matrix_samples"].rows}
    assert key in samples, f"{family}: not in the supermatrix"
    sample = samples[key]
    assert sample["sites"] == summary["alignment_sites"]
    assert sample["blocks"] == summary["genes_included"]
    assert 0 <= sample["gapped_sites"] <= sample["sites"]

    labels = _tree_labels(phylogeny_family.artifacts["tree.treefile"])
    assert labels.count(key) == 1, f"{family}: {key} appears {labels.count(key)} times in the tree"


def test_phylogeny_infers_one_supported_tree_from_the_ten_families(phylogeny_family) -> None:
    """Ten taxa, one matrix, one tree with support — and every excluded gene says why.

    IQ-TREE refuses ultra-fast bootstrap with fewer than four sequences (K-12), so a ten-taxon
    tree is the shape this module is for; the tree must carry all ten accessions and a support
    value, and a gene left out of the matrix must carry the reason rather than a run of gaps.
    """
    summary = phylogeny_family.tables["phylogeny_summary"].rows[0]
    assert summary["genomes"] == len(FAMILY_RECORDS) == 10
    assert summary["alignment_sites"] > 0
    assert summary["genes_included"] >= 2
    assert summary["tool"] == "iqtree2" and summary["version"]
    assert "-bb 1000" in summary["tool_call"], "the published tree was inferred without bootstrap"

    genes = phylogeny_family.tables["matrix_genes"].rows
    assert len(genes) == summary["genes_included"] + summary["genes_excluded"]
    for row in genes:
        if row["included"] == "yes":
            assert row["carriers"] >= 2 and row["sites"] > 0
        else:
            assert row["reason"], f"{row['gene']}: excluded from the matrix without a reason"
            assert row["sites"] is None

    fasta = list(SeqIO.parse(io.StringIO(phylogeny_family.artifacts["supermatrix.fasta"]), "fasta"))
    assert len(fasta) == len(FAMILY_RECORDS)
    assert len({len(record.seq) for record in fasta}) == 1
    assert len(fasta[0].seq) == summary["alignment_sites"]

    tree = Phylo.read(io.StringIO(phylogeny_family.artifacts["tree.treefile"]), "newick")
    assert sorted(terminal.name for terminal in tree.get_terminals()) == sorted(
        FAMILY_RECORDS.values()
    )
    assert all(
        clade.branch_length is not None for clade in tree.find_clades() if clade is not tree.root
    )
    support = _support_values(phylogeny_family.artifacts["tree.treefile"])
    assert support, "the tree carries no bootstrap or aLRT support value"
    assert all(0.0 <= value <= 100.0 for value in support)
