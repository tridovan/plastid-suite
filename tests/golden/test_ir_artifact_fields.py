"""The IR artifact's fields, on the panel's families and on PGA's annotations.

Every IR-requiring run writes an `ir` block into its manifest. The block is the evidence a reader
needs to repeat the measurement, so each field must be present and must mean what its name says:
which path the operator requested and which path produced the result, the copy lengths the record
annotated, the four alignment column classes that partition the compared columns, the scoring
policy and its version, and the independent de-novo length with its delta. A record whose
annotation and sequence disagree shows the disagreement in the delta, not in a claimed agreement.

The cases run a real `plastid run` (module `structure`, which needs the gate and nothing else) so
the assertions are about the artifact file, not the builder function. `.panel/` is gitignored;
the file skips with a named reason when it is absent.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from Bio import SeqIO

from plastid.cli import main
from plastid.core.genbank import read_genbank
from plastid.core.ir import IR_MIN_LENGTH

ROOT = Path(__file__).resolve().parents[2]
PANEL = ROOT / ".panel" / "ir-validation"
SOURCE = PANEL / "source"
VARIANTS = PANEL / "annotation-variants"
FASTA = PANEL / "fasta"

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

REQUIRED_FIELDS = (
    "requested_source",
    "status",
    "source",
    "annotation_used",
    "annotation_refused",
    "scoring_policy",
    "scoring_policy_version",
    "length",
    "identity",
    "gate_passed",
    "mirror_constant",
    "matched_columns",
    "mismatched_columns",
    "ambiguous_columns",
    "gap_columns",
    "compared_columns",
    "known_columns",
    "internal_gaps",
    "terminal_gaps",
    "alignment_method",
    "excess_at",
    "mapping_supported",
    "mapping_method",
    "annotated_lengths",
    "annotation_note",
    "denovo_length",
    "denovo_delta",
    "independent_check",
    "independent_agreement",
)

pytestmark = [
    pytest.mark.tier_b,
    pytest.mark.panel,
    pytest.mark.skipif(
        not SOURCE.is_dir(),
        reason="the IR validation panel is not in this tree (.panel/ is gitignored)",
    ),
]


def _run(tmp_path: Path, source: Path, extra: list[str] | None = None) -> dict:
    """One real run; return its manifest.

    The dataset is the fixed, already-slugged `artifact`: the run layout lower-cases the dataset
    name, and each test has its own `tmp_path`, so a fixed one cannot collide.
    """
    code = main(
        [
            "run",
            "--module",
            "structure",
            "--input",
            str(source),
            "--out",
            str(tmp_path),
            "--dataset",
            "artifact",
            *(extra or []),
        ]
    )
    assert code == 0, f"{source.name}: run exited {code}"
    (run_dir,) = [path for path in (tmp_path / "artifact").iterdir() if path.is_dir()]
    return json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))


def _assert_artifact_fields(ir: dict) -> None:
    """Every field the contract names is present, and the column classes partition compared."""
    missing = [field for field in REQUIRED_FIELDS if field not in ir]
    assert not missing, f"the IR block is missing {missing}"
    assert ir["scoring_policy"], "the scoring policy is not named"
    assert isinstance(ir["scoring_policy_version"], int) and ir["scoring_policy_version"] >= 1
    assert ir["alignment_method"]
    assert ir["excess_at"] in {"none", "outer", "inner", "both"}
    assert isinstance(ir["mapping_supported"], bool)
    assert ir["mapping_method"] is None or isinstance(ir["mapping_method"], str)
    assert ir["matched_columns"] > 0
    classes = (
        ir["matched_columns"]
        + ir["mismatched_columns"]
        + ir["ambiguous_columns"]
        + ir["gap_columns"]
    )
    assert classes == ir["compared_columns"], classes
    assert ir["internal_gaps"] + ir["terminal_gaps"] == ir["gap_columns"]
    assert 0 < ir["known_columns"] <= ir["compared_columns"]
    assert 0.0 <= ir["identity"] <= 1.0
    # The independent check is named, and its verdict is one of the three states. A check that
    # could not be made is "not measured", never silently "agrees".
    check, agreement = ir["independent_check"], ir["independent_agreement"]
    assert check in {"de_novo", "annotation", "none"}
    assert agreement in {"agrees", "differs", "not measured"}
    if check == "de_novo":
        assert ir["denovo_length"] is not None
        assert ir["denovo_delta"] == ir["length"] - ir["denovo_length"]
        assert agreement == ("agrees" if ir["denovo_delta"] == 0 else "differs")
    elif check == "annotation":
        assert ir["annotated_lengths"], "an annotation check needs the annotated lengths"
        assert ir["source"] != "annotation"
    else:
        assert ir["denovo_length"] is None
        assert agreement == "not measured"


def _annotation_lengths(genome) -> list[int]:
    """The copy lengths the annotation path can use — the long ones the IR rule counts."""
    return sorted(
        interval.length for interval in genome.repeat_regions if interval.length >= IR_MIN_LENGTH
    )


@pytest.mark.parametrize("family", tuple(FAMILY_RECORDS))
def test_each_family_records_the_full_ir_artifact(tmp_path: Path, family: str) -> None:
    """One record per family: the annotation path, with every field meaningful."""
    accession = FAMILY_RECORDS[family]
    manifest = _run(tmp_path, SOURCE / f"{accession}.gb")
    ir = manifest["ir"]
    _assert_artifact_fields(ir)
    assert ir["requested_source"] == "annotation"
    assert ir["status"] == "annotation"
    assert ir["source"] == "annotation"
    assert ir["annotation_used"] is True
    assert ir["annotation_refused"] is False
    assert ir["gate_passed"] is True
    assert sorted(ir["annotated_lengths"]) == _annotation_lengths(
        read_genbank(SOURCE / f"{accession}.gb")
    )
    assert 0.999 <= ir["identity"] <= 1.0


@pytest.mark.parametrize("record", sorted(path.name for path in VARIANTS.glob("*.gb")))
def test_every_pga_annotation_records_its_own_ir_artifact(tmp_path: Path, record: str) -> None:
    """PGA's 30 annotations: the block carries the lengths *this* annotation declares.

    A variant whose repeat_region differs from the source's, or from the sequence's own repeat,
    keeps the difference visible in `annotated_lengths` and `denovo_delta` — that is the point of
    the independent field, not a failure.
    """
    manifest = _run(tmp_path, VARIANTS / record)
    ir = manifest["ir"]
    _assert_artifact_fields(ir)
    genome = read_genbank(VARIANTS / record)
    copies = len(genome.repeat_regions)
    assert ir["requested_source"] == "annotation"
    assert ir["gate_passed"] is True
    if copies == 2:
        assert ir["status"] == "annotation"
        assert ir["annotation_used"] is True
        assert sorted(ir["annotated_lengths"]) == _annotation_lengths(genome)
    else:
        # A record the annotation path cannot use says so by name; it never claims "IR ok".
        assert ir["status"] in {
            "sequence_annotation_refused",
            "sequence_no_annotation",
        }
        assert ir["annotation_used"] is False


def test_a_variant_that_disagrees_with_its_sequence_keeps_the_delta(tmp_path: Path) -> None:
    """The de-novo length is measured independently; disagreement is a number, not a claim."""
    deltas = {}
    for path in sorted(VARIANTS.glob("*.gb")):
        ir = _run(tmp_path / path.stem, path)["ir"]
        if ir["independent_check"] != "de_novo":
            continue
        deltas[path.name] = ir["denovo_delta"]
        assert ir["denovo_delta"] == ir["length"] - ir["denovo_length"]
    assert deltas, "no variant has a de-novo cross-check"
    assert any(delta != 0 for delta in deltas.values()), (
        f"every annotation agrees with its sequence to the base: {deltas}"
    )


def test_a_sequence_only_input_is_named_in_the_artifact(tmp_path: Path) -> None:
    """A FASTA with no annotation: the artifact says `sequence_no_annotation`, not "IR ok"."""
    files = sorted(FASTA.glob("*.fa*"))
    if not files:
        pytest.skip("the rotated FASTA panel is not in this tree")
    manifest = _run(tmp_path, files[0])
    ir = manifest["ir"]
    _assert_artifact_fields(ir)
    assert ir["requested_source"] == "annotation"
    assert ir["status"] == "sequence_no_annotation"
    assert ir["source"] == "denovo_diag"
    assert ir["annotation_used"] is False
    assert ir["annotated_lengths"] == [], "an absent annotation must not be a length of zero"
    assert ir["annotation_note"] is None
    assert ir["gate_passed"] is True
    # With no second measurement there is no agreement to claim.
    assert ir["independent_check"] == "none"
    assert ir["independent_agreement"] == "not measured"
    assert ir["denovo_length"] is None
    assert SeqIO.read(str(files[0]), "fasta").id  # the file really is a bare sequence


@pytest.mark.parametrize(
    ("family", "agreement"),
    (("Malvaceae", "agrees"), ("Orchidaceae", "differs")),
)
def test_a_sequence_result_is_checked_against_the_annotation_it_did_not_use(
    tmp_path: Path, family: str, agreement: str
) -> None:
    """On the sequence path the annotation is the independent check, and its verdict is recorded.

    Malvaceae's annotation and sequence agree on the copy length; Orchidaceae's two copies are
    26,690/26,699 bp against a sequence repeat of 26,699, so the check says `differs` and both
    lengths stay in the artifact instead of the disagreement being rounded away.
    """
    accession = FAMILY_RECORDS[family]
    manifest = _run(tmp_path, SOURCE / f"{accession}.gb", extra=["--set", "ir_source=sequence"])
    ir = manifest["ir"]
    _assert_artifact_fields(ir)
    assert ir["requested_source"] == "sequence"
    assert ir["status"] == "sequence_annotation_present"
    assert ir["source"] == "denovo_diag"
    assert ir["annotation_used"] is False
    assert ir["annotation_refused"] is False
    assert ir["denovo_length"] is None
    assert ir["independent_check"] == "annotation"
    assert ir["independent_agreement"] == agreement
    assert sorted(ir["annotated_lengths"]) == _annotation_lengths(
        read_genbank(SOURCE / f"{accession}.gb")
    )
