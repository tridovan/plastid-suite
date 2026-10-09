"""The mito annotator's own check (P3A-03): the answer-key record, called from scratch.

``OX459128.1`` is a published record, so the annotation this toolkit transfers can be compared
against a curator's — the only way a self-built annotator earns its place rather than being
declared correct.

**How the comparison avoids grading its own homework.** The symbol vocabulary in
``data/mito_gene_names.yaml`` was built *from* this record's ``/product`` sentences, so asking the
vocabulary whether the caller got the names right would be circular — and the record carries no
independent symbol for a protein-coding gene to check against (measured 03/10/2026: its ``CDS``
features carry ``/locus_tag`` and a transcript id in ``/standard_name``, and only the tRNA ``gene``
features carry ``/gene``). What is *not* circular is the record's **coordinates and its exon
structure**: the test pairs every call with the published CDS its exons overlap and requires the
pairing to be a bijection over all 31 — a call that landed on the wrong locus, or two calls that
claimed one, cannot pass that. The naming half is then exactly what the vocabulary says it is, and
P3A-06/review owns judging it.

**Differences are printed, not hidden.** 18 of 31 spans come out byte-identical and the other 13
differ by at most 75 bp at one end (a reference protein's ends are not the query gene's ends); the
exon counts differ on the trans-spliced genes, where the caller reports the exons the alignment
supports and claims no join. Both are findings to record (PROJECT.md §9), never numbers to tune.

Run with ``PTS_REQUIRE_TOOLS=1`` to turn the skips into failures.
"""

from __future__ import annotations

import csv
import io
import json
import tempfile
from pathlib import Path

import pytest

from conftest import require_tools
from plastid.cli import main
from plastid.core.genbank import read_fasta, read_genbank
from plastid.core.gene_names import load_mito_aliases
from plastid.core.model import Gene, count_genes
from plastid.core.submission import to_fasta
from plastid.external.blast import BLAST_VERSION
from plastid.modules import mito_annotate
from plastid.run.context import RunContext
from plastid.run.layout import RunLayout

pytestmark = pytest.mark.tier_b

PANEL = Path(__file__).resolve().parents[2] / ".panel" / "oldenlandia-coleus-mito"
NEEDS_PANEL = pytest.mark.skipif(
    not PANEL.is_dir(),
    reason=(
        "the pinned mitogenome is not on this machine — re-fetch it with the 'plastid panel pin' "
        "command in PROJECT.md §8.6"
    ),
)

# Measured 03/10/2026 with BLAST+ 2.17.0 against the committed reference set. Not a target: a
# change that moves one of these is a finding to read first.
#
# **Moved 18 -> 19 on 10/10/2026, read first, and it is a fix** (P3B-10, §8.55). The change is one
# shared rule reaching this path: the stop codon is now read whether or not the alignment reached
# the reference's 3' end, because blastx cannot align past a stop and a target whose gene is
# shorter than the reference's therefore presents an alignment that simply ends. On the mito
# record that is ``rps1``, which was `(3, 0)` — three bases short at the 3' end, that being the low
# coordinate on a minus-strand gene — and is now `(0, 0)`, matching the published CDS. The other
# twelve differences are unchanged; §8.22's 18/31 stays in the log as what it measured then.
EXACT_SPANS = 19
MAX_END_DIFFERENCE = 75


def _run(panel: Path, tmp_path: Path, blast: dict[str, str | None]):
    genome = read_genbank(panel / "OX459128.1.gb")
    layout = RunLayout.create(tmp_path / "runs", "mito-golden", fingerprint="golden")
    layout.make_dirs()
    ctx = RunContext(layout=layout, params={}, module_id="mito_annotate", module_order=1)
    directory = str(Path(blast["blastx"]).parent)  # the fixture found the binaries in one place
    result = mito_annotate.run(genome, {"tools": directory}, ctx)
    return genome, result, ctx


def _calls(result) -> tuple[Gene, ...]:
    """The calls as features, read back through **our own reader** — the phase check itself."""
    with tempfile.TemporaryDirectory() as work:
        path = Path(work) / "annotation.gb"
        path.write_text(result.artifacts["mito_annotation.gb"], encoding="utf-8")
        back = read_genbank(path)
    return tuple(gene for gene in back.genes if gene.kind == "CDS")


def _overlap_bp(call: Gene, gene: Gene) -> int:
    return sum(
        min(part.end, other.end) - max(part.start, other.start)
        for part in call.parts
        for other in gene.parts
        if part.overlaps(other)
    )


def _pair(calls: tuple[Gene, ...], published: list[Gene]) -> dict[str, Gene]:
    """Each call → the published CDS it overlaps **most**.

    Exons, not spans: a trans-spliced gene's span covers a hundred kilobases of other genes, and
    pairing on it would be nonsense. "Most", not "exactly one": measured on this record, the
    published sdh4 (a pseudogene) overlaps the first 72 bp of cox3, so a call there has two
    candidates and only one of them is its gene.
    """
    pairs: dict[str, Gene] = {}
    for call in calls:
        overlaps = [
            (gene, _overlap_bp(call, gene))
            for gene in published
            if any(part.overlaps(other) for part in call.parts for other in gene.parts)
        ]
        assert overlaps, f"{call.name} overlaps no published CDS at all"
        gene, shared = max(overlaps, key=lambda pair: pair[1])
        assert shared >= sum(part.length for part in call.parts) // 2, (
            f"{call.name} shares only {shared} bp with its best published match"
        )
        pairs[call.name] = gene
    return pairs


@NEEDS_PANEL
def test_the_caller_returns_the_symbols_the_record_carries(tmp_path: Path, blast_tools) -> None:
    require_tools(blast_tools, "blastx", "blastn", "makeblastdb")
    _, result, ctx = _run(PANEL, tmp_path, blast_tools)

    summary = result.tables["mito_annotation_summary"].rows[0]
    assert (summary["n_cds"], summary["n_rrna"]) == (31, 3)
    assert summary["n_not_called"] == 0, "a locus the published record carries is a gap here"
    # the tool report names what ran, with the pin — no step whose version nobody can name.
    # makeblastdb runs once per search, so it appears twice.
    recorded = {(record.name, record.version) for record in ctx.tools}
    assert ("blastx", f"blastx: {BLAST_VERSION}+") in recorded
    assert ("blastn", f"blastn: {BLAST_VERSION}+") in recorded
    assert any("-seg" in record.argv for record in ctx.tools if record.name == "blastx")

    vocabulary = {locus.canonical for locus in load_mito_aliases().loci}
    called = {row["gene"] for row in result.tables["called_genes"].rows if row["kind"] == "CDS"}
    assert called == vocabulary


@NEEDS_PANEL
def test_the_calls_pair_with_the_published_CDS_one_to_one(tmp_path: Path, blast_tools) -> None:
    require_tools(blast_tools, "blastx", "blastn", "makeblastdb")
    genome, result, _ = _run(PANEL, tmp_path, blast_tools)
    published = [gene for gene in genome.genes if gene.kind == "CDS"]

    pairs = _pair(_calls(result), published)

    assert len(pairs) == len(published) == 31
    assert len(set(id(gene) for gene in pairs.values())) == 31, "two calls claim one published CDS"


@NEEDS_PANEL
def test_the_coordinates_agree_within_the_measured_bound(tmp_path: Path, blast_tools) -> None:
    require_tools(blast_tools, "blastx", "blastn", "makeblastdb")
    genome, result, _ = _run(PANEL, tmp_path, blast_tools)
    published = [gene for gene in genome.genes if gene.kind == "CDS"]
    calls = _calls(result)
    pairs = _pair(calls, published)

    exact = 0
    differences: list[str] = []
    for call in calls:
        gene = pairs[call.name]
        delta = (call.span[0] - gene.span[0], call.span[1] - gene.span[1])
        assert call.strand == gene.strand, f"{call.name} landed on the wrong strand"
        parts = (
            "same" if len(call.parts) == len(gene.parts) else f"{len(call.parts)}≠{len(gene.parts)}"
        )
        if delta == (0, 0) and parts == "same":
            exact += 1
            continue
        assert abs(delta[0]) <= MAX_END_DIFFERENCE and abs(delta[1]) <= MAX_END_DIFFERENCE, (
            f"{call.name} moved {delta} against the published CDS — past the measured bound, so "
            f"this is a new finding rather than the known reference-end effect"
        )
        differences.append(f"{call.name} {delta} exons {parts}")

    print(
        f"\nP3A-03: {exact}/{len(published)} calls identical to the published record; "
        f"{len(differences)} differ: {'; '.join(sorted(differences))}"
    )
    assert exact == EXACT_SPANS, "the measured agreement changed — read the printed differences"


# Measured 03/10/2026 with ARAGORN 1.2.41 on OX459128.1. The two sets have the same size and are
# not the same set — each difference is a finding (PROJECT.md §9), printed by the test and never
# adjusted away here.
TRNA_MISSING = 1  # trnA-UGC 13244-13348: in the record, not called by ARAGORN
TRNA_EXTRA = 1  # trnL-AAG 36355-36427: called by ARAGORN, not in the record
TRNA_RENAMED = 1  # the tRNA-His at 11741-11819: the record says trnH-GUG, ARAGORN reads trnH-AUG


def _spans(row: dict, genome_length: int) -> list[tuple[int, int]]:
    """A row's 1-based spans — a wrapped call is two of them, exactly as its Gene is."""
    if row["crosses_origin"] == "yes":
        return [(row["start"], genome_length), (1, row["end"])]
    return [(row["start"], row["end"])]


def _touches(spans: list[tuple[int, int]], gene) -> bool:
    return any(start <= gene.span[1] and gene.span[0] <= end for start, end in spans)


@NEEDS_PANEL
def test_the_trna_count_is_reported_against_the_published_one(tmp_path: Path, blast_tools) -> None:
    """P3A-04's Done-when: the count next to the published one, with every difference — extra or
    missing — listed by name, position and anticodon.

    **Per feature, never per name.** Measured: the record's 17 tRNAs carry only 14 distinct names
    (``trnH-GUG``, ``trnK-UUU`` and ``trnM-CAU`` each appear twice), so a comparison that
    deduplicates by name reports 14 and is wrong.

    What it measured on this record: **17 = 17, but not the same 17** — one locus the record
    annotates and ARAGORN does not call, one the other way round, and one whose anticodon the two
    disagree about. The test asserts those counts so a change is noticed, and prints the detail.
    """
    require_tools(blast_tools, "blastx", "blastn", "makeblastdb", "aragorn")
    genome, result, _ = _run(PANEL, tmp_path, blast_tools)
    published = [gene for gene in genome.genes if gene.kind == "tRNA"]
    ours = result.tables["trna_calls"].rows

    missing: list[str] = []
    extra: list[str] = []
    renamed: list[str] = []
    moved: list[str] = []
    claimed: set[tuple[int, int]] = set()
    for row in ours:
        spans = _spans(row, len(genome))
        overlaps = [gene for gene in published if _touches(spans, gene)]
        if not overlaps:
            extra.append(f"{row['gene']} {row['start']}-{row['end']}")
            continue
        gene = max(overlaps, key=lambda gene: gene.length)
        claimed.add(gene.span)
        if row["gene"] != gene.name:
            renamed.append(
                f"{row['start']}-{row['end']}: ARAGORN says {row['gene']} "
                f"({row['anticodon']}), the record says {gene.name}"
            )
        elif (row["start"], row["end"]) != gene.span:
            moved.append(
                f"{row['gene']}: {row['start']}-{row['end']} vs {gene.span[0]}-{gene.span[1]}"
            )
    missing = [
        f"{gene.name} {gene.span[0]}-{gene.span[1]}"
        for gene in published
        if gene.span not in claimed
    ]
    print(
        f"\nP3A-04: ARAGORN called {len(ours)} tRNA; the record has {len(published)}. "
        f"missing: {missing or 'none'}; extra: {extra or 'none'}; "
        f"name or anticodon differs: {renamed or 'none'}; "
        f"span differs: {moved or 'none'}"
    )

    assert len(published) == 17 and len(ours) == 17
    assert (len(missing), len(extra), len(renamed)) == (TRNA_MISSING, TRNA_EXTRA, TRNA_RENAMED)


@NEEDS_PANEL
def test_the_annotator_output_re_reads_and_drives_the_shared_checks(
    tmp_path: Path, blast_tools
) -> None:
    """The phase check: the annotator's own output re-reads with our reader, and the ORF profile
    comes from the same function a downloaded record gets — no second code path."""
    require_tools(blast_tools, "blastx", "blastn", "makeblastdb")
    _, result, _ = _run(PANEL, tmp_path, blast_tools)

    calls = _calls(result)
    assert len(calls) == 31
    status = {row["status"]: row["n"] for row in result.tables["orf_status"].rows}
    assert sum(status.values()) == 31
    print(f"\nP3A-03 ORF profile of the calls: { {k: v for k, v in status.items() if v} }")


@NEEDS_PANEL
def test_the_cli_annotates_a_bare_sequence_and_drives_the_normal_runner(
    tmp_path: Path, blast_tools
) -> None:
    """P3A-05's Done-when: the de-annotated sequence in, a GenBank out that re-reads with our own
    reader, and that file driving ``orf`` + ``mito_genes`` through the ordinary runner — no second
    code path for the annotation command (seam 3)."""
    require_tools(blast_tools, "blastx", "blastn", "makeblastdb", "aragorn")
    record = read_genbank(PANEL / "OX459128.1.gb")
    sample = tmp_path / "sample.fna"
    sample.write_text(to_fasta(record, "OX459128.1"), encoding="utf-8")
    assert read_fasta(sample).genes == (), "the input must carry no annotation at all"
    directory = str(Path(blast_tools["blastx"]).parent)  # type: ignore[arg-type]

    exit_code = main(
        [
            "mito",
            "annotate",
            "--input",
            str(sample),
            "--out",
            str(tmp_path / "runs"),
            "--set",
            f"tools={directory}",
        ]
    )
    assert exit_code == 0

    produced = sorted((tmp_path / "runs").rglob("mito_annotation.gb"))
    assert len(produced) == 1, f"expected one artifact, found {produced}"
    annotated = read_genbank(produced[0])
    counts = count_genes(annotated)
    assert (counts.total, counts.cds, counts.trna, counts.rrna) == (51, 31, 17, 3)

    # …and the same file is an ordinary input to the ordinary runner
    assert (
        main(
            [
                "run",
                "-m",
                "orf",
                "-m",
                "mito_genes",
                "--input",
                str(produced[0]),
                "--out",
                str(tmp_path / "runs2"),
            ]
        )
        == 0
    )
    runs = sorted((tmp_path / "runs2").rglob("02_mito_genes/tables/mito_summary.csv"))
    assert len(runs) == 1
    # exports carry their provenance as `#` comments, which the stdlib reader cannot skip
    body = [
        line
        for line in runs[0].read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]
    summary = next(csv.DictReader(io.StringIO("\n".join(body))))
    # the exported header carries the human labels, not the column ids
    assert (summary["CDS"], summary["tRNA"], summary["rRNA"]) == ("31", "17", "3")


@NEEDS_PANEL
def test_the_three_counts_are_printed_next_to_the_published_ones(
    tmp_path: Path, blast_tools
) -> None:
    """P3A-06: the check prints what it counted **beside** what the record carries, and itemizes
    every difference instead of summarizing it into one word.

    A count that disagreed would be recorded in PROJECT.md §9 as a finding and this test would
    still be a finished task — the failure mode it exists to prevent is a number quietly tuned
    until it matched.
    """
    require_tools(blast_tools, "blastx", "blastn", "makeblastdb", "aragorn")
    genome, result, _ = _run(PANEL, tmp_path, blast_tools)
    annotated = read_genbank(_artifact(result, tmp_path))
    ours = count_genes(annotated)
    published = count_genes(genome)

    print(
        f"\nP3A-06 — this run vs the published record:\n"
        f"  CDS   {ours.cds:>4}  published {published.cds:>4}\n"
        f"  tRNA  {ours.trna:>4}  published {published.trna:>4}\n"
        f"  rRNA  {ours.rrna:>4}  published {published.rrna:>4}\n"
        f"  total {ours.total:>4}  published {published.total:>4}"
    )
    assert (ours.cds, ours.trna, ours.rrna) == (published.cds, published.trna, published.rrna)


def _artifact(result, tmp_path: Path) -> Path:
    """The module's GenBank, written out so it can be re-read — the phase's own check."""
    path = tmp_path / "mito_annotation.gb"
    path.write_text(result.artifacts["mito_annotation.gb"], encoding="utf-8")
    return path


COLEUS = PANEL.parent / "oldenlandia-coleus-mito"


@NEEDS_PANEL
def test_the_sample_nobody_has_annotated_is_annotated_pinned_and_read_offline(
    tmp_path: Path, blast_tools
) -> None:
    """P3A-07: the shape every later mito run takes — annotate a record with no features, pin the
    result with its checksum, and read the gene set back offline from the pin.

    Measured on ``AP045469.1`` (*Coleus hadiensis*, 313,585 bp, **0 features**): **45 genes = 23
    CDS + 19 tRNA + 3 rRNA**, and **8 vocabulary loci not called**. That is the prototype's honest
    answer for this sample, not a failure to hide: four of the eight (``mttB`` · ``rps1`` ·
    ``rps13`` · ``rps4``) have **no alignment at all even at an E-value of 10**, one (``cox2``) has
    only a 40 aa fragment and three others only fragments at E-values between 0.05 and 2.8 — so
    the floors are not what is keeping them out. Whether those genes are absent from the genome or
    absent from this assembly is not something a reference-transfer annotator can decide.
    """
    require_tools(blast_tools, "blastx", "blastn", "makeblastdb", "aragorn")
    source = COLEUS / "AP045469.1.gb"
    record = read_genbank(source)
    assert record.genes == (), "this task is about a record with no annotation to copy"
    sample = tmp_path / "coleus.fna"
    sample.write_text(to_fasta(record, "AP045469.1"), encoding="utf-8")
    directory = str(Path(blast_tools["blastx"]).parent)  # type: ignore[arg-type]

    assert (
        main(
            [
                "mito",
                "annotate",
                "--input",
                str(sample),
                "--out",
                str(tmp_path / "runs"),
                "--set",
                f"tools={directory}",
            ]
        )
        == 0
    )
    produced = sorted((tmp_path / "runs").rglob("mito_annotation.gb"))
    assert len(produced) == 1

    pinned = tmp_path / "panel"
    assert (
        main(
            [
                "panel",
                "pin",
                "--file",
                str(produced[0]),
                "--note",
                "prototype mito annotate, test",
                "--out",
                str(pinned),
            ]
        )
        == 0
    )
    entry = json.loads((pinned / "PINNED.json").read_text(encoding="utf-8"))["genomes"][0]
    assert entry["parse_error"] is None, "a pin that does not read back is not a pin"
    assert entry["sha256"] and entry["note"].startswith("prototype mito annotate")
    assert (entry["cds"], entry["trna"], entry["rrna"]) == (23, 19, 3)

    # …and the gene set is readable from the pin with no tool and no reference set at all
    offline = tmp_path / "offline"
    assert (
        main(
            [
                "run",
                "-m",
                "mito_genes",
                "--input",
                entry["path"],
                "--out",
                str(offline),
            ]
        )
        == 0
    )
    summary_path = next(offline.rglob("mito_summary.csv"))
    body = [
        line
        for line in summary_path.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]
    summary = next(csv.DictReader(io.StringIO("\n".join(body))))
    print(
        f"\nP3A-07 — {entry['accession']} from the pin: {summary['Genes']} genes = "
        f"{summary['CDS']} CDS + {summary['tRNA']} tRNA + {summary['rRNA']} rRNA "
        f"(sha256 {entry['sha256'][:12]}…)"
    )
    assert (summary["Genes"], summary["CDS"], summary["tRNA"], summary["rRNA"]) == (
        "45",
        "23",
        "19",
        "3",
    )
