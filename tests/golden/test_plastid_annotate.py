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
# The one record of the pinned set the audit blocks — measured 09/10/2026, §8.43's own prediction.
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
                # The pinned set carries one record the audit blocks (NC_024936.1, ``rps19``
                # annotated as partial on both copies), and the gate stops a run that would
                # search it (P3B-09). Dropping it by name is D5-08's own loop — "fix or drop the
                # bad file → run again" — and the run records the drop rather than making it.
                "--set",
                f"exclude={BLOCKED}",
            ]
        )
        == 0
    )
    run_dir = sorted((out / "runs" / GOLDEN.lower()).iterdir())[-1]
    return run_dir / "01_plastid_annotate"


def test_the_published_overlaps_survive_annotation(annotated: Path) -> None:
    """§8.34's measured pairs, on the *called* output rather than on the published record.

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
    """PGA's measured failure (§8.15: 147 genes from 7 references, four loci called twice).

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
    choice prefers it. Structure decides first (D5-08 §3), so the complete feature is the one
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


def test_the_audit_gate_stops_the_run_on_a_record_it_blocks(
    data_dir: Path,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """P3B-09's wiring, measured on the project's own pinned set — and **before any search**.

    H-19 named this hole (§8.48): ``plastid_annotate`` did not call the audit at all, so a record
    the audit BLOCKs was still a candidate in a run. The gate now runs over exactly the records
    the run would search, ahead of the first BLAST — which is why this check needs no BLAST
    itself, and why "nothing annotated" is literal rather than hopeful.
    """
    out = tmp_path / "runs"
    code = main(
        [
            "annotate",
            "-i",
            str(data_dir / f"{GOLDEN}.gb"),
            "--reference",
            str(PANEL),
            "-o",
            str(out),
        ]
    )
    assert code == 2, "a blocked reference record did not stop the run"
    message = capsys.readouterr().err
    assert f"{BLOCKED}.gb — CDS rps19: partial coordinates" in message
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
    asked it not to (the audit's rule for a pin that applies to nothing, §8.43)."""
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
    distinct loci are 79, and all four rRNA loci are inside the repeat (§8.34's 132 = 86+38+8).
    """
    genome = read_genbank(annotated / "plastid_annotation.gb")
    assert genome.source.get("organelle") == "plastid"
    counts = Counter(gene.kind for gene in genome.genes)
    assert counts["CDS"] == 86, f"the inverted repeat's copies are not both written: {counts}"
    assert counts["rRNA"] == 8, f"the four rRNA loci are inside the repeat: {counts}"
    # The contract itself: 132 = 86 + 38 + 8 (PROJECT.md §6, the invariant that is never written
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
