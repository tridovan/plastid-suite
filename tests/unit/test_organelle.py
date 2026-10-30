"""One run serves one organelle (D5-06, H-09) — the declaration, and the two refusal directions.

The rule arrived in P3B-09 because only *one* of the two mismatch directions had a gate: the IR
gate stopped a plastid module on a mitogenome (measured 03/10), and nothing stopped
``mito_genes`` on a plastome (measured 05/10). Both directions are the same question, so
both are answered here.
"""

from __future__ import annotations

import pytest

from plastid.core.model import Genome
from plastid.core.organelle import (
    GROUP_MITO,
    GROUP_PLASTID,
    GROUP_SHARED,
    GROUPS,
    MITO,
    MITO_PATH,
    PLASTID,
    PLASTID_PATH,
    OrganelleError,
    check_set,
    organelle_of,
)


def _genome(accession: str, organelle: str | None) -> Genome:
    source = {} if organelle is None else {"organelle": organelle}
    return Genome(name=accession, accession=accession, sequence="ACGT" * 25, source=source)


@pytest.mark.parametrize(
    ("declared", "expected"),
    [
        ("plastid", PLASTID),  # measured 09/10/2026: 5 of the 34 pinned records spell it this way
        ("plastid:chloroplast", PLASTID),  # …and 29 this way
        ("chloroplast", PLASTID),
        ("mitochondrion", MITO),  # all 8 records of the mito reference set
        ("mitochondrion:mitochondrion", MITO),
        ("PLASTID", PLASTID),  # the qualifier is not case-sensitive in practice
    ],
)
def test_a_records_own_qualifier_is_what_is_read(declared: str, expected: str) -> None:
    assert organelle_of(_genome("X.1", declared)) == expected


def test_a_bare_sequence_declares_nothing_and_is_not_guessed() -> None:
    """D5-06: a bare sequence declares itself **by the path it enters**. Nothing else may.

    A mitogenome cannot be told from a plastome by its bytes — both are circular DNA with
    protein-coding genes, and a plant mitogenome carries plastid-derived inserts — so a module
    that guessed here would be inventing the answer it is supposed to check.
    """
    assert organelle_of(_genome("X.1", None)) is None
    assert organelle_of(_genome("X.1", "")) is None


def test_a_qualifier_naming_neither_organelle_is_returned_as_written() -> None:
    """Not filed under a guess: the caller reports it, and a human reads what the record said."""
    assert organelle_of(_genome("X.1", "apicoplast")) == "apicoplast"


def test_the_right_organelle_passes_and_a_bare_sequence_passes_with_it() -> None:
    check_set([_genome("A.1", "plastid:chloroplast")], group=GROUP_PLASTID, path=PLASTID_PATH)
    check_set([_genome("A.1", None)], group=GROUP_PLASTID, path=PLASTID_PATH)
    check_set([_genome("M.1", "mitochondrion")], group=GROUP_MITO, path=MITO_PATH)


def test_both_mismatch_directions_stop_by_name_and_point_at_the_right_path() -> None:
    """H-09: today only one direction has a gate — and only because an IR gate happens to sit in
    the way, which is a different question that sometimes gives the right answer."""
    with pytest.raises(OrganelleError) as plastid_path:
        check_set([_genome("M.1", "mitochondrion")], group=GROUP_PLASTID, path=PLASTID_PATH)
    assert "M.1" in str(plastid_path.value)
    assert MITO_PATH in str(plastid_path.value)

    with pytest.raises(OrganelleError) as mito_path:
        check_set([_genome("A.1", "plastid")], group=GROUP_MITO, path=MITO_PATH)
    assert "A.1" in str(mito_path.value)
    assert PLASTID_PATH in str(mito_path.value)


def test_a_mixed_set_is_refused_whichever_group_asked_for_it() -> None:
    """One run, one panel (D5-06) — and a shared module is not exempt from it: it may run on
    either organelle, never on both at once."""
    mixed = [_genome("A.1", "plastid"), _genome("M.1", "mitochondrion")]
    for group in (GROUP_PLASTID, GROUP_MITO, GROUP_SHARED):
        with pytest.raises(OrganelleError, match="one organelle"):
            check_set(mixed, group=group, path=PLASTID_PATH)


def test_the_mixed_set_message_names_every_declaring_input() -> None:
    mixed = [_genome("A.1", "plastid:chloroplast"), _genome("M.1", "mitochondrion")]
    with pytest.raises(OrganelleError) as raised:
        check_set(mixed, group=GROUP_SHARED, path=PLASTID_PATH)
    assert "A.1 (plastid)" in str(raised.value)
    assert "M.1 (mitochondrion)" in str(raised.value)


def test_an_unknown_group_is_a_programming_error_not_a_data_one() -> None:
    with pytest.raises(ValueError, match="unknown module group"):
        check_set([_genome("A.1", "plastid")], group="plastid-only", path=PLASTID_PATH)


def test_the_registry_carries_each_modules_group() -> None:
    """H-09's first clause: the declaration lives in the registry, not in a reader's memory.

    The two exact sets are deliberate tripwires: a module added without a group —
    or a group changed without the table — makes this fail, and updating it means reading the
    table. ``plastid_annotate`` is in the plastid row because it transfers from plastid records
    and reads the inverted repeat; ``rna_editing`` because its reference B′ holds chloroplast
    proteins only (a plastid predictor scored against a mitogenome would answer a question nobody
    asked); the seven "transitional" analyses are ``shared`` because that
    is what the table calls them, and what keeps them off a mitogenome today is ``requires_ir``,
    a different question with its own gate.
    """
    from plastid.modules import registry

    modules = registry.all_modules()
    assert modules, "no module is registered — the import side effect did not run"
    assert {spec.group for spec in modules.values()} <= set(GROUPS)
    assert {mid for mid, spec in modules.items() if spec.group == GROUP_MITO} == {
        "mito_annotate",
        "mito_genes",
    }
    assert {mid for mid, spec in modules.items() if spec.group == GROUP_PLASTID} == {
        "barcodes",
        "gene_normalize",
        "gene_table",
        "genome_stats",
        "plastid_annotate",
        "rna_editing",
        "structure",
    }


def test_a_module_may_not_be_registered_without_saying_which_organelle_it_serves() -> None:
    """No default: an omission is not a claim (the decorator's own reason)."""
    from plastid.modules import registry

    with pytest.raises(TypeError, match="group"):
        registry.module("a_module_with_no_group", version="0.1.0", title="x")(lambda *a: None)
