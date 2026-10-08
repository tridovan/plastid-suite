"""The three numbers D5-09 §9.3 asks for, measured (P3B-02, PROJECT.md §8.41).

D5-09's naming arbitration has four tests; three of its numbers were *open* until they were
measured, and the design record's own reasons for two of them do not survive the measurement:

* the **margin** of test 3 — how much better the best match must be than the runner-up. Measured
  in **bits**, because a margin in identity points is meaningless when the runner-up's alignment
  is short (psaB against psaA scores a *higher* pident over a partial alignment at a lower
  bitscore, which is how a "negative margin" appears).
* whether test 4 (**neighbours**) separates `rpl2`/`rpl23` and the two `trnL` — it does not need
  to: they are 22–44 % and 70.6 % identical, and the ranking separates them. The measured cases
  needing test 4 are elsewhere, and are pinned below.
* the **`identity_floor`** — the identity below which a match is noise. The measured gap is
  between 45.0 % (the best a genuinely different locus reaches) and 67.75 % (the worst a true
  ortholog reaches), which is where the mito path's 50 % already sits.

The measurement is over the corpus (33 published records queried against the golden genome) with
BLAST+ 2.17.0, `-seg no`, E-value 10 — the program and flags D5-08's engine uses.
"""

from __future__ import annotations

import statistics
from pathlib import Path

import pytest
from Bio import Align, SeqIO

from conftest import require_tools
from plastid.core.genbank import read_genbank
from plastid.core.gene_names import cds_by_canonical_name, normalize_name
from plastid.core.model import coding_sequence
from plastid.run.tools import run_tool

ROOT = Path(__file__).resolve().parents[2]
CORPUS = ROOT / ".panel" / "iridaceae-34"
GOLDEN_FIXTURE = ROOT / "tests" / "data" / "NC_050833.1.gb"
GOLDEN = "NC_050833.1"

pytestmark = pytest.mark.tier_b

needs_panel = pytest.mark.skipif(
    not CORPUS.is_dir(),
    reason="the pinned Iridaceae corpus is not in this working tree (.panel/ is gitignored)",
)

_ALIGNER = Align.PairwiseAligner()
_ALIGNER.mode = "global"
_ALIGNER.match_score = 1.0
_ALIGNER.mismatch_score = 0.0
_ALIGNER.open_gap_score = -2.0
_ALIGNER.extend_gap_score = -0.5


def identity(a: str, b: str) -> float:
    """Percent identity of the best global alignment — for sequences too short for BLAST seeds.

    The alignment parameters are part of the number and are recorded with it (§8.41): gaps cost
    -2 to open and -0.5 to extend, so the metric stays close to an ungapped comparison at tRNA
    length. With cheaper gaps the same pair of trnG sequences reads 78.9 % instead of 73.2 %,
    which is why §8.41 states the parameters and not just the identity.
    """
    alignment = _ALIGNER.align(a, b)[0]
    matches = sum(1 for x, y in zip(alignment[0], alignment[1], strict=True) if x == y and x != "-")
    return 100.0 * matches / max(len(a), len(b))


def cds_rows(path: Path) -> list[tuple[str, str, str]]:
    """(symbol, canonical locus, the record's own /translation) for every CDS that has one."""
    genome = read_genbank(path)
    canonical = {gene.name: name for name, gene in cds_by_canonical_name(genome).items()}
    rows = []
    for record in SeqIO.parse(str(path), "genbank"):
        for feature in record.features:
            if feature.type != "CDS" or not feature.qualifiers.get("translation"):
                continue
            symbol = feature.qualifiers.get("gene", ["?"])[0]
            rows.append(
                (
                    symbol,
                    canonical.get(symbol, symbol),
                    str(feature.qualifiers["translation"][0]).replace(" ", ""),
                )
            )
    return rows


def trna_rows(path: Path) -> list[tuple[str, str, str]]:
    """(symbol, canonical locus, spliced sequence) for every tRNA in the record."""
    genome = read_genbank(path)
    return [
        (gene.name, normalize_name(gene.name).canonical, coding_sequence(genome, gene).upper())
        for gene in genome.genes
        if gene.kind == "tRNA"
    ]


def blastp_hits(
    tools: dict[str, str | None], tmp_path: Path, queries: list[tuple[str, str, str]], path: Path
) -> dict[str, list[dict]]:
    """Every query protein against every subject protein, keyed by canonical locus, best first."""
    subject = tmp_path / "subjects.faa"
    query = tmp_path / "queries.faa"
    subject.write_text("".join(f">{locus}\n{seq}\n" for _, locus, seq in queries), encoding="utf-8")
    query.write_text("".join(f">{locus}\n{seq}\n" for _, locus, seq in path), encoding="utf-8")
    log = tmp_path / "logs"
    log.mkdir(exist_ok=True)
    database = tmp_path / "db"
    run_tool(
        "makeblastdb",
        [tools["makeblastdb"], "-in", str(subject), "-dbtype", "prot", "-out", str(database)],
        log_dir=log,
    )
    record = run_tool(
        "blastp",
        [
            tools["blastp"],
            "-query",
            str(query),
            "-db",
            str(database),
            "-seg",
            "no",
            "-outfmt",
            "6 qseqid sseqid pident length bitscore qlen",
            "-evalue",
            "10",
            "-max_target_seqs",
            "100",
        ],
        log_dir=log,
    )
    # the search's captured stdout is the result; the recorded argv says what produced it
    hits: dict[str, list[dict]] = {}
    for line in Path(record.stdout_file).read_text(encoding="utf-8").splitlines():
        qid, sid, pident, length, bits, qlen = line.split("\t")
        hits.setdefault(qid, []).append(
            {
                "locus": sid,
                "pident": float(pident),
                "length": int(length),
                "bits": float(bits),
                "coverage": int(length) / int(qlen),
            }
        )
    for rows in hits.values():
        rows.sort(key=lambda row: -row["bits"])
    return hits


def test_the_golden_carries_two_trnG_genes_under_one_name() -> None:
    """A duplicate name in the contract record itself — and the two genes are not the same gene.

    Measured 08/10/2026 (§8.41): `NC_050833.1` annotates **two** tRNA features as `trnG-UCC`,
    at 8889-8911/9586-9633 and 36413-36483, and the two spliced sequences are 73.2 % identical.
    One name, two genes: the naming step may not merge them, and the ★ pick inherits the problem.
    """
    genome = read_genbank(GOLDEN_FIXTURE)
    trng = [
        (gene, coding_sequence(genome, gene).upper())
        for gene in genome.genes
        if gene.kind == "tRNA" and gene.name.startswith("trnG")
    ]
    assert [gene.name for gene, _ in trng] == ["trnG-UCC", "trnG-UCC"]
    assert len({sequence for _, sequence in trng}) == 2
    assert identity(trng[0][1], trng[1][1]) == pytest.approx(73.2, abs=0.05)
    assert trng[0][0].length == trng[1][0].length == 71


@needs_panel
def test_the_two_trnG_sequences_carry_both_names_across_the_corpus() -> None:
    """Sequence does not settle the spelling: records disagree in both directions.

    The golden's first sequence is called `trnG-GCC` by *Gladiolus* and `trnG-UCC` by the rest;
    the golden's second is called `trnG-GCC` by 25 records and `trnG-UCC` by the comparison
    genome and seven others. So H-04 is not a spelling the toolkit may normalise away — it is a
    curator's call, and a rename that would merge the two is exactly what §5's stop tier is for.
    """
    golden = trna_rows(CORPUS / f"{GOLDEN}.gb")
    first, second = [sequence for _, _, sequence in golden if _.startswith("trnG")][:2]
    names = {"first": set(), "second": set()}
    for path in sorted(CORPUS.glob("*.gb")):
        if path.name == f"{GOLDEN}.gb":
            continue
        for symbol, _, sequence in trna_rows(path):
            if not symbol.startswith("trnG"):
                continue
            if identity(sequence, first) > identity(sequence, second):
                names["first"].add(symbol)
            else:
                names["second"].add(symbol)
    assert "trnG-GCC" in names["first"], "no record spells the golden's first sequence trnG-GCC"
    assert "trnG-GCC" in names["second"], "no record spells the golden's second trnG-GCC"


@needs_panel
def test_the_margin_and_the_floor_are_where_the_naming_step_will_read_them(
    tmp_path: Path, blast_tools: dict[str, str | None]
) -> None:
    """The numbers P3B-04's thresholds have to sit inside (PROJECT.md §8.41).

    * **The margin is a score margin, in bits.** The smallest a *true* match achieved is 19.2
      bits, so a floor of 10 — half of it — refuses a tie and nothing else. A margin in identity
      points is meaningless: psaB against psaA scores a *higher* pident over a partial alignment
      at a lower bitscore, which is how a "negative margin" appears.
    * **The floor needs coverage as well as identity.** With coverage >= 0.7 the highest identity
      a genuinely different locus reaches is 44.93 % (`psaB` vs `psaA`, 750 aa) and the lowest a
      true ortholog reaches is 67.75 % (`ycf1`). At coverage >= 0.5 a 36-aa gene walks into the
      gap: `psbJ` vs `matK` reaches **52.38 %** over 21 aa — which is why the mito path's 0.5
      coverage floor is not enough for plastid CDS, whose short genes are the whole problem.
    * **Three hits sit above the floor and are not counterexamples.** `lhbA` vs `psbZ` and
      `ycf4` vs `psbN` reach 97.7-100 % because they *are* the same gene under a name the record
      spells differently (§8.31's vocabulary gap). A floor refuses noise; it must not refuse a
      misnamed real gene — that is what the naming step exists to surface.
    """
    require_tools(blast_tools, "blastp", "makeblastdb")
    subjects = cds_rows(CORPUS / f"{GOLDEN}.gb")
    margins: list[float] = []
    ortholog: list[float] = []
    wrong: list[float] = []
    for path in sorted(CORPUS.glob("*.gb")):
        if path.name == f"{GOLDEN}.gb":
            continue
        rows = cds_rows(path)
        hits = blastp_hits(blast_tools, tmp_path, subjects, rows)
        for _, locus, _ in rows:
            found = hits.get(locus)
            if not found:
                continue
            best = found[0]
            runner = next((row for row in found if row["locus"] != best["locus"]), None)
            margins.append(best["bits"] - (runner["bits"] if runner else 0.0))
            covered = [row for row in found if row["coverage"] >= 0.7]
            mine = [row for row in covered if row["locus"] == locus]
            others = [row for row in covered if row["locus"] != locus]
            if mine:
                ortholog.append(max(row["pident"] for row in mine))
            if others:
                wrong.append(max(row["pident"] for row in others))

    assert min(margins) == pytest.approx(19.2, abs=0.1)
    assert statistics.median(margins) == pytest.approx(298.0, abs=5.0)
    assert min(ortholog) == pytest.approx(67.75, abs=0.05)
    assert max(wrong) == pytest.approx(100.0, abs=0.05)
    # the three that sit above: the same locus under a name the record spells differently
    assert sorted(value for value in wrong if value >= min(ortholog)) == pytest.approx(
        [97.674, 100.0, 100.0], abs=0.01
    )
    assert max(value for value in wrong if value < min(ortholog)) == pytest.approx(44.93, abs=0.05)
    # the floor the mito path already uses (identity 50 %, coverage 0.5) fits the gap only once
    # the coverage floor is raised: 44.93 < 50 < 67.75
    assert max(value for value in wrong if value < min(ortholog)) < 50.0 < min(ortholog)


@needs_panel
def test_the_serine_trnas_are_where_the_ranking_fails() -> None:
    """For tRNA the sequence ranking alone is not enough — and the failing family is trnS.

    Measured (§8.41): of 1,182 cross-record tRNA comparisons, **26** rank a different locus
    first, every one a `trnS` query whose best match is another serine tRNA, by up to 15.05
    identity points. The true locus can be as low as 53.76 % identical while a *different* locus
    reaches 73.97 % — the two ranges overlap, so **no identity floor separates tRNA loci**; the
    decision has to come from the name, the anticodon or the neighbours. That, and not
    `rpl2`/`rpl23`, is the measured reason test 4 exists.
    """
    reference = trna_rows(CORPUS / f"{GOLDEN}.gb")
    by_locus = {}
    for _, locus, sequence in reference:
        by_locus.setdefault(locus, sequence)
    true_scores: list[float] = []
    wrong_wins = 0
    wrong_family: set[str] = set()
    for path in sorted(CORPUS.glob("*.gb")):
        if path.name == f"{GOLDEN}.gb":
            continue
        for symbol, locus, sequence in trna_rows(path):
            if locus not in by_locus:
                continue
            true_identity = identity(sequence, by_locus[locus])
            true_scores.append(true_identity)
            best_other = max(
                (identity(sequence, seq) for other, seq in by_locus.items() if other != locus),
                default=0.0,
            )
            if best_other > true_identity:
                wrong_wins += 1
                wrong_family.add(symbol.split("-")[0])
    assert wrong_wins == 26
    assert wrong_family == {"trnS"}
    assert min(true_scores) == pytest.approx(53.76, abs=0.05)
