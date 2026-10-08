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

ROOT = Path(__file__).resolve().parents[2]
PANEL = ROOT / ".panel" / "plastid-reference"
GOLDEN = "NC_050833.1"

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
    assert len(genome.repeat_regions) == 2
