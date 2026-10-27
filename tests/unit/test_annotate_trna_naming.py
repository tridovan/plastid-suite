"""H-24: an intron tRNA is named by the reference set's own locus, not ARAGORN's anticodon.

ARAGORN reads the three anticodon bases beside an intron and gets them wrong on measured
records: the golden's ``trnV-UAC`` comes out ``trnC-ACA``, and on the Crataegus panel the same
locus is ``tRNA-Stop`` and no gene is written at all. The reference set annotates that locus —
its own tRNA at the same place has a name — so the name is taken from there when a strict
majority of supporting records agree. These checks hold the rule: rename on a majority, keep
ARAGORN's name without evidence, recover an unnameable line, and never rename on a tie.
"""

from __future__ import annotations

from plastid.core.caller import Hsp, ReferenceEntry
from plastid.external.aragorn import AragornResult, TrnaHit
from plastid.modules.plastid_annotate import _apply_trna_naming

GENOME_LENGTH = 5_000
SUBJECT = "ACGT" * 19  # 76 bp, a plausible spliced tRNA
LANE = "trnV-UAC"


def _entry(accession: str, name: str = LANE) -> ReferenceEntry:
    return ReferenceEntry(
        accession=accession,
        organism="Testus testus",
        source_gene=name,
        vocabulary=name,
        sequence=SUBJECT,
        junctions=(),
    )


def _hsp(accession: str, name: str, qstart: int, qend: int, sstart: int, send: int) -> Hsp:
    return Hsp(
        subject_key=f"{accession}|{name}|{name}",
        pident=100.0,
        length=abs(qend - qstart) + 1,
        qstart=qstart,
        qend=qend,
        sstart=sstart,
        send=send,
        evalue=1e-9,
        bitscore=70.0,
        frame=None,
    )


def _exons(accession: str, name: str = LANE) -> tuple[Hsp, Hsp]:
    """The two exons of a minus-strand intron tRNA at the locus, as blastn reports them.

    The query is always on its plus strand; a minus-strand gene's subject matches in the reverse
    direction (``sstart > send``). Exon 1 is the high-coordinate piece, exon 2 the low one.
    """
    return (
        _hsp(accession, name, 1639, 1677, 39, 1),
        _hsp(accession, name, 1000, 1036, 76, 40),
    )


def _call(name: str | None, *, raw: str = "") -> TrnaHit:
    return TrnaHit(
        name=name,
        source_name="tRNA-Cys" if name else "tRNA-Stop",
        anticodon="ACA",
        start=1000,
        end=1678,
        strand=-1,
        score=36,
        intron=(40, 583),
        raw_line=raw or "18  tRNA-Cys c[1000,1678] (aca)i(40,583)",
    )


def _result(call: TrnaHit, *, others: tuple[str, ...] = ()) -> AragornResult:
    return AragornResult(
        hits=(call,) if call.name else (),
        genes=(call.gene(GENOME_LENGTH),) if call.name else (),
        declared=1,
        others=others,
        raw="",
        version="1.2.41",
        calls=(call,),
    )


def test_a_strict_majority_of_the_reference_renames_the_intron_trna() -> None:
    entries = [_entry("A"), _entry("B"), _entry("C")]
    hits = [hit for accession in "ABC" for hit in _exons(accession)]
    call = _call("trnC-ACA")

    updated, namings = _apply_trna_naming(_result(call), hits, entries, GENOME_LENGTH)

    assert [gene.name for gene in updated.genes] == ["trnV-UAC"]
    assert len(namings) == 1
    assert (namings[0].old, namings[0].new, namings[0].supporting, namings[0].considered) == (
        "trnC-ACA",
        "trnV-UAC",
        3,
        3,
    )
    assert not namings[0].recovered


def test_one_exon_is_not_enough_evidence_for_a_locus() -> None:
    """An intron tRNA's own exons are about half the molecule each; one HSP cannot carry a name."""
    entries = [_entry("A")]
    hits = [_exons("A")[0]]  # only the high-coordinate exon

    updated, namings = _apply_trna_naming(_result(_call("trnC-ACA")), hits, entries, GENOME_LENGTH)

    assert [gene.name for gene in updated.genes] == ["trnC-ACA"]
    assert namings == ()


def test_without_reference_evidence_aragorns_name_stands() -> None:
    """The negative case the rule needs: a locus the reference does not carry is not renamed."""
    entries = [_entry("A", "trnI-GAU")]
    # evidence for a *different* name, far away, so the containment test rejects it
    hits = [
        _hsp("A", "trnI-GAU", 3000, 3038, 39, 1),
        _hsp("A", "trnI-GAU", 3040, 3076, 40, 76),
    ]

    updated, namings = _apply_trna_naming(_result(_call("trnC-ACA")), hits, entries, GENOME_LENGTH)

    assert [gene.name for gene in updated.genes] == ["trnC-ACA"]
    assert namings == ()


def test_an_unnameable_line_is_recovered_by_its_locus_not_dropped() -> None:
    """Crataegus: the locus is ``tRNA-Stop`` to ARAGORN, so no gene was written at all."""
    entries = [_entry("A"), _entry("B")]
    hits = [hit for accession in "AB" for hit in _exons(accession)]
    raw = "19  tRNA-Stop c[1000,1678] (tta)i(36,577)"
    call = _call(None, raw=raw)

    updated, namings = _apply_trna_naming(
        _result(call, others=(raw,)), hits, entries, GENOME_LENGTH
    )

    assert [gene.name for gene in updated.genes] == ["trnV-UAC"]
    assert updated.others == ()
    assert namings[0].recovered is True
    assert namings[0].old is None


def test_a_tie_does_not_rename() -> None:
    """Two reference records disagreeing is not a majority; ARAGORN's name stands."""
    entries = [_entry("A", "trnV-UAC"), _entry("B", "trnI-GAU")]
    hits = [hit for hit in _exons("A", "trnV-UAC")] + [hit for hit in _exons("B", "trnI-GAU")]

    updated, namings = _apply_trna_naming(_result(_call("trnC-ACA")), hits, entries, GENOME_LENGTH)

    assert [gene.name for gene in updated.genes] == ["trnC-ACA"]
    assert namings == ()
