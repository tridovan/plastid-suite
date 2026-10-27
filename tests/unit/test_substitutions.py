"""The substitutions module's CDS scope with MAFFT stubbed — the panel's real shape.

Measured 04/10/2026: 3 of the 11 pinned stand-in plastomes carry **no annotated CDS**, and the
first version of this scope returned an empty table for the whole run. A record with no
annotation cannot share a CDS with anything, so it is excluded **and named** — the run continues
for the records that can be compared.
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import pytest

import plastid.modules.substitutions as substitutions
from plastid.core.align import SUBSTITUTION_TYPES, Alignment
from plastid.core.model import Gene, Genome, Interval
from plastid.run.context import RunContext
from plastid.run.layout import RunLayout

_GENE = "ATGGCTTCTAAACCCGGGTTTTAA"


@pytest.fixture
def stubbed(monkeypatch):
    def fake_align(rows, ctx):
        return Alignment(
            names=tuple(name for name, _ in rows), sequences=tuple(seq for _, seq in rows)
        )

    monkeypatch.setattr(substitutions, "align", fake_align)


def _genome(accession: str, genes: dict[str, str]) -> Genome:
    chunks, records, offset = [], [], 0
    for name, sequence in genes.items():
        chunks.append(sequence)
        records.append(
            Gene(name, "CDS", (Interval(offset, offset + len(sequence)),), transl_table=11)
        )
        offset += len(sequence)
    return Genome(
        name=f"synthetic {accession}",
        accession=accession,
        sequence="".join(chunks),
        genes=tuple(records),
        topology="linear",
    )


def _ctx(tmp_path: Path) -> RunContext:
    layout = RunLayout.create(tmp_path / "runs", "substitutions-test", fingerprint="test")
    layout.make_dirs()
    return RunContext(layout=layout, params={}, module_id="substitutions", module_order=1)


def test_a_record_without_annotation_is_named_not_silently_dropped(tmp_path: Path, stubbed) -> None:
    reference = _genome("REF.1", {"rbcL": _GENE})
    query = _genome("QRY.1", {"rbcL": _GENE.replace("GGGTTT", "GGGTTA")})
    empty = Genome(name="unannotated", accession="EMPTY.1", sequence="ACGT" * 20)
    result = substitutions.run((reference, query, empty), {"scope": "cds"}, _ctx(tmp_path))
    summary = result.tables["substitutions_summary"].rows[0]
    assert summary["n_genomes"] == 2  # the empty record is not compared …
    assert any("EMPTY.1" in note for note in result.notes)  # … and it is named


def test_the_cds_scope_refuses_to_report_nothing_when_no_genome_is_annotated(
    tmp_path: Path, stubbed
) -> None:
    reference = _genome("REF.1", {"rbcL": _GENE})
    empty = Genome(name="unannotated", accession="EMPTY.1", sequence="ACGT" * 20)
    with pytest.raises(ValueError, match="no genome other than the reference"):
        substitutions.run((reference, empty), {"scope": "cds"}, _ctx(tmp_path))


def test_the_cds_scope_reports_each_query_genome_separately(tmp_path: Path, stubbed) -> None:
    """A manuscript table is per pair: Ts/Tv between the reference and *each* congener, not one
    number pooled over every query. Measured 04/10/2026 — the pooled number alone cannot fill
    that table, which is what P4-05 needs."""
    reference = _genome("REF.1", {"rbcL": _GENE})
    close = _genome("CLOSE.1", {"rbcL": _GENE.replace("GGGTTT", "GGGTTC")})  # one transition
    far = _genome("FAR.1", {"rbcL": _GENE.replace("GGGTTT", "GGGAAA")})  # three T→A
    result = substitutions.run((reference, close, far), {"scope": "cds"}, _ctx(tmp_path))
    by_genome = {row["genome"]: row for row in result.tables["substitutions_by_genome"].rows}
    assert set(by_genome) == {"CLOSE.1", "FAR.1"}
    assert (by_genome["CLOSE.1"]["transitions"], by_genome["CLOSE.1"]["transversions"]) == (1, 0)
    assert (by_genome["FAR.1"]["transitions"], by_genome["FAR.1"]["transversions"]) == (0, 3)
    assert by_genome["CLOSE.1"]["genes"] == 1
    # and the pooled summary is still the sum of the pairs
    summary = result.tables["substitutions_summary"].rows[0]
    assert summary["transitions"] == 1
    assert summary["transversions"] == 3


def test_the_type_table_carries_every_directed_change_per_genome(tmp_path: Path, stubbed) -> None:
    """The twelve directed changes are what the paper's six complementary pairs are summed from.

    Measured 17/10/2026: the module held only transitions and transversions, so a figure of
    substitution type × species had nothing to read.
    """
    reference = _genome("REF.1", {"rbcL": _GENE})
    close = _genome("CLOSE.1", {"rbcL": _GENE.replace("GGGTTT", "GGGTTC")})  # one T→C
    far = _genome("FAR.1", {"rbcL": _GENE.replace("GGGTTT", "GGGAAA")})  # three T→A
    result = substitutions.run((reference, close, far), {"scope": "cds"}, _ctx(tmp_path))
    table = result.tables["substitutions_by_type"]
    counts = {(row["genome"], row["substitution_type"]): row["count"] for row in table.rows}
    assert counts[("CLOSE.1", "T_to_C")] == 1
    assert counts[("FAR.1", "T_to_A")] == 3
    # every genome carries all twelve, a zero included — an absent row would be ambiguous
    assert len(table.rows) == 2 * len(SUBSTITUTION_TYPES)
    assert sum(value for (genome, _), value in counts.items() if genome == "FAR.1") == 3


def test_the_type_table_agrees_with_the_aggregates_it_is_summed_from(
    tmp_path: Path, stubbed
) -> None:
    """A type table that disagreed with `substitutions_by_genome` would be a second truth."""
    reference = _genome("REF.1", {"rbcL": _GENE})
    close = _genome("CLOSE.1", {"rbcL": _GENE.replace("GGGTTT", "GGGTTC")})
    far = _genome("FAR.1", {"rbcL": _GENE.replace("GGGTTT", "GGGAAA")})
    result = substitutions.run((reference, close, far), {"scope": "cds"}, _ctx(tmp_path))
    per_genome = {row["genome"]: row for row in result.tables["substitutions_by_genome"].rows}
    by_type: dict[str, dict[str, int]] = {}
    for row in result.tables["substitutions_by_type"].rows:
        bucket = by_type.setdefault(row["genome"], {"Transition": 0, "Transversion": 0})
        bucket[row["class"]] += row["count"]
    for genome, row in per_genome.items():
        assert by_type[genome]["Transition"] == row["transitions"]
        assert by_type[genome]["Transversion"] == row["transversions"]


def test_the_whole_scope_also_carries_the_type_table(
    tmp_path: Path, stubbed, golden_genome: Genome
) -> None:
    """Both scopes build the table from the same counters, so a figure needs one rule for either.

    The whole scope needs a genome whose IR the gate accepts, so this one uses the golden sample
    with five bases changed — the CDS scope above covers the same code path on synthetic records.
    """
    mutated = replace(golden_genome, sequence=golden_genome.sequence.replace("A", "G", 5))
    result = substitutions.run((golden_genome, mutated), {"scope": "whole"}, _ctx(tmp_path))
    rows = result.tables["substitutions_by_type"].rows
    counts = {(row["genome"], row["substitution_type"]): row["count"] for row in rows}
    assert counts[(mutated.accession, "A_to_G")] == 5
    assert len(rows) == len(SUBSTITUTION_TYPES)
    assert "twelve directed changes" in result.tables["substitutions_by_type"].scope.label


def test_the_default_scope_is_whole_and_each_scope_names_itself(
    tmp_path: Path, stubbed, golden_genome: Genome
) -> None:
    """M11.5: the two scopes are different numbers and the output says which one it is.

    The default stays the whole alignment; the CDS scope is the golden's narrower method. A
    report that quotes one as the other is the failure this distinction exists to prevent, so
    neither the default nor the label may drift silently.
    """
    assert substitutions.DEFAULT_SCOPE == "whole"
    mutated = replace(golden_genome, sequence=golden_genome.sequence.replace("A", "G", 5))

    whole = substitutions.run((golden_genome, mutated), {}, _ctx(tmp_path))
    assert whole.tables["substitutions_summary"].rows[0]["scope"].startswith("the whole aligned")
    assert "not the CDS-only scope" in whole.tables["substitutions_summary"].scope.label

    cds = substitutions.run((golden_genome, golden_genome), {"scope": "cds"}, _ctx(tmp_path))
    assert "CDS only" in cds.tables["substitutions_summary"].scope.label


def test_the_region_table_uses_the_normalized_frame_on_a_rotated_reference(
    tmp_path: Path, stubbed, golden_genome: Genome
) -> None:
    """Measured 27/10/2026: the region loop compared the reference's **raw** positions with
    region boundaries that live in the **normalized** frame. At the golden's ``rotation_offset``
    of 0 the two coincide, so only a record that starts anywhere else shows it: the same
    biological substitution is then counted in whichever region holds the raw coordinate.

    A per-region column count cannot see the bug — the four regions cover the genome and keep
    their lengths under rotation — so the test puts one substitution at a known LSC base and
    asks which region the table credits it to.
    """
    from plastid.core.normalize import rotate_genome

    def mutated_at(genome: Genome, raw_position: int) -> Genome:
        base = genome.sequence[raw_position]
        change = {"A": "G", "G": "A", "C": "T", "T": "C"}[base]
        sequence = genome.sequence[:raw_position] + change + genome.sequence[raw_position + 1 :]
        return replace(genome, sequence=sequence)

    def credited_region(reference: Genome, raw_position: int) -> str:
        result = substitutions.run(
            (reference, mutated_at(reference, raw_position)), {}, _ctx(tmp_path)
        )
        rows = [
            row
            for row in result.tables["substitutions_by_region"].rows
            if row["substitutions"] == 1
        ]
        assert len(rows) == 1
        return str(rows[0]["region"])

    shift = 7500
    rotated = rotate_genome(golden_genome, shift)
    normal_base = 100  # inside the golden's LSC (IRa starts at 83,127)
    assert credited_region(golden_genome, normal_base) == "LSC"
    assert credited_region(rotated, (normal_base - shift) % len(golden_genome)) == "LSC"


def test_a_partial_reference_gene_is_counted_in_the_cds_scope(tmp_path: Path, stubbed) -> None:
    from dataclasses import replace

    reference = _genome("REF.1", {"rbcL": _GENE})
    reference = replace(
        reference,
        genes=(
            replace(
                reference.genes[0],
                parts=(replace(reference.genes[0].parts[0], partial_start=True),),
            ),
        ),
    )
    query = _genome("QRY.1", {"rbcL": _GENE})
    result = substitutions.run((reference, query), {"scope": "cds"}, _ctx(tmp_path))
    assert result.tables["substitutions_summary"].rows[0]["n_partial_genes"] == 1
