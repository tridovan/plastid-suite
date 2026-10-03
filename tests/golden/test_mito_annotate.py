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

import tempfile
from pathlib import Path

import pytest

from conftest import require_tools
from plastid.core.genbank import read_genbank
from plastid.core.gene_names import load_mito_aliases
from plastid.core.model import Gene
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
EXACT_SPANS = 18
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
