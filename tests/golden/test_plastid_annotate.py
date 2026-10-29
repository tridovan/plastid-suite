"""The plastid annotator end to end, on the golden and against the pinned reference set (P3B-08).

The reference set lives outside the repository (``.panel/`` is gitignored), so the check is
skipped with a named reason where it is not; on a machine that has it, this is the phase's own
claim measured on real data — one call per locus, the losers listed, the published overlaps kept,
and an artifact that re-reads and says which organelle it is.

One run serves the whole module (it is three BLAST searches and an ARAGORN call, minutes not
seconds), so the checks below read one run directory rather than re-annotating the golden.
"""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

import pytest

from conftest import require_tools
from plastid.cli import main
from plastid.core.genbank import read_genbank
from plastid.result.export import read_result

ROOT = Path(__file__).resolve().parents[2]
PANEL = ROOT / ".panel" / "plastid-reference"
GOLDEN = "NC_050833.1"
# The one record of the pinned set the audit blocks — measured 09/10/2026, as predicted.
BLOCKED = "NC_024936.1"

pytestmark = [
    pytest.mark.tier_b,
    pytest.mark.skipif(
        not PANEL.is_dir(),
        reason="the pinned plastid reference set is not in this tree (.panel/ is gitignored)",
    ),
]


def _rows(module_dir: Path, name: str) -> list[list[str]]:
    """A table's rows as plain cell lists, header dropped.

    Read by position on purpose: these checks are about *what the run decided*, and a lookup by
    header label would pass even if the exporter wrote the columns in another order — while the
    reader of the file would be looking at different ones.
    """
    path = module_dir / "tables" / f"{name}.csv"
    lines = [line for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    body = [line for line in lines if not line.startswith("#")]
    return list(csv.reader(body))[1:]


@pytest.fixture(scope="module")
def annotated(
    tools: dict[str, str | None],
    blast_tools: dict[str, str | None],
    data_dir: Path,
    tmp_path_factory: pytest.TempPathFactory,
) -> Path:
    """The golden, annotated once — the run directory every check below reads."""
    require_tools(tools, "mafft")
    require_tools(blast_tools, "blastx", "blastn", "makeblastdb", "aragorn")
    out = tmp_path_factory.mktemp("plastid_annotate")
    assert (
        main(
            [
                "run",
                "-m",
                "plastid_annotate",
                "-i",
                str(data_dir / f"{GOLDEN}.gb"),
                "-o",
                str(out / "runs"),
                "--set",
                f"reference={PANEL}",
                # The whole pinned set is searched: H-21 was fixed in M11.10, so NC_024936.1 —
                # whose rps19 carries an inverted-repeat partial flag over a complete 92 aa
                # translation — is REPORTed by the audit instead of blocked. The target's own
                # record is still left out (leave-one-out), which is 33 references.
            ]
        )
        == 0
    )
    run_dir = sorted((out / "runs" / GOLDEN.lower()).iterdir())[-1]
    return run_dir / "01_plastid_annotate"


def test_the_published_overlaps_survive_annotation(annotated: Path) -> None:
    """The measured overlapping pairs, on the *called* output rather than on the published record.

    Three of the four are pairs of different genes and are asserted exactly. The fourth —
    ``rps12`` against itself, 114 bp — is the published record listing a trans-spliced gene as
    two features whose first exon is shared; a caller that makes **one call per locus** has one
    ``rps12``, so that pair cannot appear, and asserting it would assert a shape the engine does
    not have. What must hold is the rule the pair exists to protect: overlap is kept and
    reported, never resolved by dropping the weaker call.
    """
    by_pair: dict[frozenset[str], str] = {}
    for first, second, overlap, _strands in _rows(annotated, "overlapping_calls"):
        by_pair[frozenset((first.split()[0], second.split()[0]))] = overlap
    assert by_pair.get(frozenset({"psbD", "psbC"})) == "53"
    assert by_pair.get(frozenset({"ndhK", "ndhC"})) == "10"
    assert by_pair.get(frozenset({"atpE", "atpB"})) == "4"


def test_one_call_per_locus_and_every_loser_listed(annotated: Path) -> None:
    """PGA's measured failure (147 genes from 7 references, four loci called twice).

    Every locus in ``competing_references`` marks exactly one candidate chosen — and at least one
    locus has a loser, or the table would be evidence of nothing.
    """
    chosen: Counter[str] = Counter()
    candidates: Counter[str] = Counter()
    for row in _rows(annotated, "competing_references"):
        gene, _accession, _spelling, _length = row[0], row[1], row[2], row[3]
        candidates[gene] += 1
        if row[8] == "yes":
            chosen[gene] += 1
    assert chosen, "no locus has a chosen reference — the arbitration table is empty"
    assert set(chosen.values()) == {1}, f"a locus was chosen more than once: {dict(chosen)}"
    assert max(candidates.values()) > 1, "no locus had a loser — nothing was arbitrated"

    called = [row[0] for row in _rows(annotated, "called_genes")]
    assert called, "nothing was called"
    assert len(called) == len(set(called)), f"a locus was called twice: {sorted(called)}"


def test_ycf1_lands_on_the_published_gene_not_on_a_truncated_reference(annotated: Path) -> None:
    """H-19's negative test, measured 08/10/2026.

    Before the structural rule, ``ycf1`` was called at **108,448-109,347** — 900 bp at 98.83 %
    identity, from a reference whose own ycf1 is a 300-residue truncated copy of the gene — where
    the published gene is **5,421 bp** at 122,984-128,404. 24 of the 34 pinned records annotate the
    locus twice, complete and truncated, and the truncated copy scores *higher*: an identity-first
    choice prefers it. Structure decides first (D5-08), so the complete feature is the one
    searched — and the run says which feature it did not search rather than dropping it in silence.
    """
    rows = {row[0]: row for row in _rows(annotated, "called_genes")}
    start, end, identity = rows["ycf1"][2], rows["ycf1"][3], rows["ycf1"][6]
    assert (start, end) == ("122984", "128404"), (
        f"ycf1 was called at {start}-{end}, not on the published gene — a truncated reference won "
        f"on identity at {identity} %"
    )
    seen = read_result((annotated / "result.json").read_text(encoding="utf-8"))
    assert any("annotates ycf1 twice" in note for note in seen.notes), (
        "the run searched a feature it did not report leaving out"
    )


def test_the_default_set_is_usable_after_the_h21_fix() -> None:
    """The pinned set minus the target — 33 records — audits with nothing blocked.

    Measured 26/10/2026: the same set used to block NC_024936.1 three times (two partial
    coordinates and its annotation-first IR at 99.80 %). M11.10 reports the partial flags (the
    protein is complete) and reports the annotation-first refusal while naming the sequence path.
    """
    from plastid.modules import plastid_annotate

    with plastid_annotate._reference_set({"reference": str(PANEL)}, GOLDEN) as reference:
        searched = {genome.accession for genome in reference.records}
        assert BLOCKED in searched, "the fixed record is still being dropped"
        assert GOLDEN not in searched, "the target's own record was searched — leave-one-out"
        assert len(searched) == 33
        assert "audited (0 blocked)" in reference.label
        assert reference.dropped == ()


def test_the_audit_gate_stops_the_run_on_a_record_it_blocks(
    data_dir: Path,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """P3B-09's wiring — and **before any search**. The blocking class is a CDS with no protein.

    H-19 named this hole: ``plastid_annotate`` did not call the audit at all, so a record
    the audit BLOCKs was still a candidate in a run. The gate now runs over exactly the records
    the run would search, ahead of the first BLAST — which is why this check needs no BLAST
    itself, and why "nothing annotated" is literal rather than hopeful. The fixture is the
    second golden record written without its ``/translation`` qualifiers, so the block is a real
    one rather than the H-21 shape this milestone stopped blocking.
    """
    from plastid.core.genbank import to_genbank

    reference_dir = tmp_path / "reference"
    reference_dir.mkdir()
    stripped = read_genbank(data_dir / f"{BLOCKED}.gb")
    (reference_dir / f"{BLOCKED}.gb").write_text(to_genbank(stripped), encoding="utf-8")

    out = tmp_path / "runs"
    code = main(
        [
            "annotate",
            "-i",
            str(data_dir / f"{GOLDEN}.gb"),
            "--reference",
            str(reference_dir),
            "--set",
            "reannotate=yes",
            "-o",
            str(out),
        ]
    )
    assert code == 2, "a blocked reference record did not stop the run"
    message = capsys.readouterr().err
    assert f"{BLOCKED}.gb — CDS" in message and "no /translation" in message
    assert not list(out.glob("**/plastid_annotation.gb")), "the run annotated from a blocked set"


def test_a_record_dropped_by_name_is_not_searched_and_the_run_says_which() -> None:
    """D5-08's loop, as a mechanism: drop the file, and the run records the drop.

    The assertions are on the reference set the module actually hands the caller — the records, the
    label a reader sees, and the drop — because that object *is* what the searches are built from.
    """
    from plastid.modules import plastid_annotate

    with plastid_annotate._reference_set(
        {"reference": str(PANEL), "exclude": [BLOCKED]}, GOLDEN
    ) as reference:
        searched = {genome.accession for genome in reference.records}
        assert BLOCKED not in searched
        assert GOLDEN not in searched, "the target's own record was searched — leave-one-out"
        assert len(searched) == 32, f"34 records, one target, one drop: got {len(searched)}"
        assert reference.dropped == ((BLOCKED, "left out by --set exclude"),)
        assert "audited (0 blocked)" in reference.label
        assert f"dropped by name ({BLOCKED})" in reference.label


def test_an_exclusion_that_matches_nothing_stops_the_run() -> None:
    """A name that matches no record is not a no-op — it is a run that searched what the user
    asked it not to (the audit's rule for a pin that applies to nothing)."""
    from plastid.modules import plastid_annotate

    with (
        pytest.raises(plastid_annotate.ReferenceSetMissing, match="matches no record"),
        plastid_annotate._reference_set(
            {"reference": str(PANEL), "exclude": ["NC_999999.1"]}, GOLDEN
        ),
    ):
        pass


def test_the_artifact_re_reads_and_says_which_organelle_it_is(annotated: Path) -> None:
    """H-11, measured: the mito artifact copied an empty ``source`` through and failed
    ``submission``'s first pass 24 times. A record this toolkit writes declares its organelle,
    and the inverted repeat's genes are written **twice** — the golden counts 86 CDS where its
    distinct loci are 79, and all four rRNA loci are inside the repeat (132 = 86+38+8).
    """
    genome = read_genbank(annotated / "plastid_annotation.gb")
    assert genome.source.get("organelle") == "plastid"
    counts = Counter(gene.kind for gene in genome.genes)
    assert counts["CDS"] == 86, f"the inverted repeat's copies are not both written: {counts}"
    assert counts["rRNA"] == 8, f"the four rRNA loci are inside the repeat: {counts}"
    # The contract itself: 132 = 86 + 38 + 8 (the invariant that is never written
    # any other way). Measured 09/10/2026: the tRNA count was **30** until ARAGORN was asked for
    # intron-containing genes — it skips them in silence unless `-i` is passed, and eight of this
    # record's thirty-eight tRNAs carry an intron (trnK-UUU · trnL-UAA · trnV-UAC · one trnG ·
    # both copies of trnI-GAU and trnA-UGC).
    assert counts["tRNA"] == 38, f"the intron-containing tRNAs are missing: {counts}"
    assert sum(counts.values()) == 132, f"132 = 86 + 38 + 8, got {dict(counts)}"
    assert len(genome.repeat_regions) == 2


def test_every_copy_the_inverted_repeat_is_given_lands_where_the_record_puts_it(
    annotated: Path, data_dir: Path
) -> None:
    """The mirror, checked against the published record rather than against its own arithmetic.

    Measured 09/10/2026 (P3B-10): the copy was placed one **repeat length** from where it belongs
    — 26,214 bp on the golden, i.e. outside the repeat altogether and into the single-copy region
    — because the inversion was computed from the first repeat's *start* instead of its *end*. The
    gene **count** was right (86 either way), which is why every count-based check passed while
    every mirrored coordinate was wrong. A published plastome's own second copy is the check that
    catches it, so that is what is asserted.
    """
    ours = read_genbank(annotated / "plastid_annotation.gb")
    published = read_genbank(data_dir / f"{GOLDEN}.gb")

    def copies(genome, name: str) -> list[tuple[int, int]]:
        return sorted(
            (min(part.start for part in gene.parts) + 1, max(part.end for part in gene.parts))
            for gene in genome.genes
            if gene.name == name and gene.kind in ("CDS", "rRNA")
        )

    checked = 0
    # `rps12` is deliberately absent: its envelope is short by a *missing exon*, not by
    # a misplaced copy, and this check is about where a copy lands. The difference list above
    # is where that shortfall is asserted.
    for name in ("rpl2", "rpl23", "ndhB", "ycf2", "rps7", "rrn16", "rrn5"):
        want = copies(published, name)
        if len(want) < 2:  # a single-copy locus says nothing about the mirror
            continue
        got = copies(ours, name)
        assert len(got) == len(want), (
            f"{name}: {len(got)} copies written, the record has {len(want)}"
        )
        for (want_start, want_end), (got_start, got_end) in zip(want, got, strict=True):
            # the detected repeat's boundaries differ from the record's by a few bases
            assert abs(got_start - want_start) <= 10 and abs(got_end - want_end) <= 10, (
                f"{name}: written at {got_start}-{got_end}, the record puts that copy at "
                f"{want_start}-{want_end}"
            )
        checked += 1
    assert checked >= 5, "the check found no repeated locus to verify the mirror on"


# ------------------------------------------------------------------ every difference, listed

# P3B-10's Done-when asks for the golden's differences from the published record "kind + span,
# nothing summarized away" — so they are written down, all 35, and the check below fails on any
# that moves. Measured 10/10/2026 with the code as it stands, from the run this module's
# fixture makes; a change here is a change in the claim, not a test that needs loosening.
#
# What landed on 10/10/2026: the exon-structure transfer took the golden from 44
# differences to 36 — the first exons of `petB` · `petD` · `rpl16` · `rps16`, `rpoC1`'s and
# `clpP`'s intron boundaries and `rps12`'s 26 bp third exon are now the published record's own —
# and reading the stop codon past an alignment that stopped short took `psaJ` (35).
#
# What is left is measured, and each line of it is one of three things rather than a defect with
# a task number:
#   * **the reference's own structure differs from the published record's**, by 2–3 bp at an
#     intron boundary, and the transfer is faithful to the reference it chose: `rpl2` on both
#     copies (391+431 → 393+432, which is how NC_042827.1 spells it), `rpl16` (399 → 402,
#     NC_088562.1), `ycf3` (230+153 → 228+155) and `rps12`'s second exon (229 → 232), both
#     NC_072562.1. The corpus measures this class at 34 records, and it is the price of
#     transferring a structure instead of inventing one;
#   * `ndhD`: the chosen reference annotates it 523 aa where the record has 501, and the target's
#     own DNA supports the longer model at 97.90 % over all 523 residues, starting at TTG (an
#     initiator) against the record's ATC (not one) — an annotation disagreement about where the
#     gene starts, which an annotator may report but not settle;
#   * **ARAGORN's tRNA ends against the submitters'**: `trnP-UGG` 2 bp past the record's 3' end on
#     every record in the corpus (34 of them) and `trnS-UGA` read on the opposite strand. tRNA is
#     called de novo by ARAGORN (D5-08), so this is a tool's convention next to a submitter's, and
#     trimming to match would be fitting the number to the answer;
#   * the four tRNA names ARAGORN reads differently (`trnG-GCC` twice, `trnM-CAU` twice): H-04's
#     family, a vocabulary convention rather than a locus. The intron loci that used to be here
#     (`trnV-UAC`, `trnI-GAU` twice) are now **named** from the reference set's own annotation
#     (H-24), so their rows moved from "vocabulary/name" to "conclusion/extent": the envelope and
#     name now agree with the published record, and what remains is ARAGORN's intron boundary
#     inside the gene (H-22/H-23, M11.15).
GOLDEN_DIFFERENCES: tuple[tuple[str, str, str, str, str], ...] = (
    ("boundary", "rRNA", "extent", "rrn23 129981-132790(-1)", "rrn23 129982-132791(-1)"),
    ("boundary", "tRNA", "extent", "trnA-UGC 103034-103915(+1)", "trnA-UGC 103034-103915(+1)"),
    ("boundary", "tRNA", "extent", "trnA-UGC 132937-133818(-1)", "trnA-UGC 132937-133818(-1)"),
    ("boundary", "tRNA", "extent", "trnD-GUC 30541-30614(-1)", "trnD-GUC 30540-30615(-1)"),
    ("boundary", "tRNA", "extent", "trnE-UUC 31141-31213(-1)", "trnE-UUC 31140-31214(-1)"),
    ("boundary", "tRNA", "extent", "trnH-GUG 153198-153272(-1)", "trnH-GUG 153197-153271(-1)"),
    ("boundary", "tRNA", "extent", "trnH-GUG 83580-83654(+1)", "trnH-GUG 83581-83655(+1)"),
    ("boundary", "tRNA", "extent", "trnK-UUU 1368-4057(-1)", "trnK-UUU 1367-4058(-1)"),
    ("boundary", "tRNA", "extent", "trnL-CAA 143081-143161(+1)", "trnL-CAA 143080-143162(+1)"),
    ("boundary", "tRNA", "extent", "trnL-CAA 93691-93771(-1)", "trnL-CAA 93690-93772(-1)"),
    ("boundary", "tRNA", "extent", "trnN-GUU 108051-108122(-1)", "trnN-GUU 108050-108123(-1)"),
    ("boundary", "tRNA", "extent", "trnN-GUU 128730-128801(+1)", "trnN-GUU 128729-128802(+1)"),
    ("boundary", "tRNA", "extent", "trnR-ACG 107666-107739(+1)", "trnR-ACG 107666-107740(+1)"),
    ("boundary", "tRNA", "extent", "trnR-ACG 129113-129186(-1)", "trnR-ACG 129112-129186(-1)"),
    ("boundary", "tRNA", "extent", "trnS-GCU 7928-8015(-1)", "trnS-GCU 7927-8015(-1)"),
    ("boundary", "tRNA", "extent", "trnT-GGU 31775-31846(+1)", "trnT-GGU 31774-31847(+1)"),
    ("boundary", "tRNA", "extent", "trnT-UGU 46327-46399(-1)", "trnT-UGU 46326-46400(-1)"),
    ("boundary", "tRNA", "extent", "trnY-GUA 30998-31081(-1)", "trnY-GUA 30997-31082(-1)"),
    ("conclusion", "CDS", "extent", "ndhD 114661-116166(-1)", "ndhD 114661-116232(-1)"),
    ("conclusion", "CDS", "extent", "rpl16 80539-81873(-1)", "rpl16 80539-81873(-1)"),
    ("conclusion", "CDS", "extent", "rpl2 151668-153154(+1)", "rpl2 151668-153154(+1)"),
    ("conclusion", "CDS", "extent", "rpl2 83698-85184(-1)", "rpl2 83698-85184(-1)"),
    ("conclusion", "CDS", "extent", "rps12 68976-139634(-1)", "rps12 68976-139634(-1)"),
    ("conclusion", "CDS", "extent", "rps12 68976-98013(-1)", "rps12 68976-98016(-1)"),
    ("conclusion", "CDS", "extent", "ycf3 42508-44457(-1)", "ycf3 42508-44457(-1)"),
    ("conclusion", "tRNA", "extent", "trnI-GAU 101952-102969(+1)", "trnI-GAU 101952-102969(+1)"),
    ("conclusion", "tRNA", "extent", "trnI-GAU 133883-134900(-1)", "trnI-GAU 133883-134900(-1)"),
    ("conclusion", "tRNA", "extent", "trnP-UGG 66042-66115(-1)", "trnP-UGG 66040-66115(-1)"),
    ("conclusion", "tRNA", "extent", "trnS-UGA 35472-35564(+1)", "trnS-UGA 35471-35565(-1)"),
    ("conclusion", "tRNA", "extent", "trnV-UAC 51598-52275(-1)", "trnV-UAC 51598-52275(-1)"),
    ("vocabulary", "tRNA", "name", "trnG-UCC 36413-36483(+1)", "trnG-GCC 36413-36484(+1)"),
    ("vocabulary", "tRNA", "name", "trnG-UCC 8889-9633(+1)", "trnS-CGA 8889-9634(+1)"),
    ("vocabulary", "tRNA", "name", "trnI-CAU 151155-151228(+1)", "trnM-CAU 151155-151229(+1)"),
    ("vocabulary", "tRNA", "name", "trnI-CAU 85624-85697(-1)", "trnM-CAU 85623-85697(-1)"),
    ("vocabulary", "tRNA", "name", "trnfM-CAU 36714-36787(-1)", "trnM-CAU 36713-36787(-1)"),
)


def test_every_difference_from_the_published_record_is_listed(
    annotated: Path, data_dir: Path
) -> None:
    """The phase's Done-when, first clause: 35 differences, every one of them with its span.

    The count is deliberately not asserted as a number — the *set* is, so a difference that moves
    or changes verdict fails here and has to be written down with its evidence, and one that
    quietly disappears fails too. The defects still in the list are today's measured truth, and
    the comment above says which task removes each one.
    """
    from tools.corpus_bar import compare

    published = read_genbank(data_dir / f"{GOLDEN}.gb")
    ours = read_genbank(annotated / "plastid_annotation.gb")
    measured = {
        (d.verdict, d.kind, d.detail, d.published, d.called)
        for d in compare(list(published.genes), list(ours.genes))
    }

    gone = sorted(set(GOLDEN_DIFFERENCES) - measured)
    new = sorted(measured - set(GOLDEN_DIFFERENCES))
    assert not gone and not new, (
        "the golden's differences from the published record moved — every one has to be "
        "re-measured and written down:\n  no longer seen: "
        + "\n                 ".join(map(str, gone))
        + "\n  not in the list: "
        + "\n                   ".join(map(str, new))
    )
