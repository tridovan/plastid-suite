"""The Điềm corpus fixtures — test data this project did not pick for itself (P4a-02).

34 Iridaceae plastomes arrived in **two annotations of the same sequences**: the published NCBI
record, and a PGA re-annotation — the tool named by the records themselves, `/note="Annotation
Method :: PGA-Plastid Genome Annotator"`, present in 34/34 of the second set and in none of the
first (measured 04/10/2026, PROJECT.md §8.31). The corpus is 28 MB and lives outside the repo;
what is vendored here is the smallest subset in which every file carries a written-down expected
answer.

Three pairs, three different jobs:

* **the golden genome both ways** — the contract is already written down (§6), and the PGA
  annotation of the very same sequence measures **85 CDS where the contract says 86**;
* ***Iris hermona* both ways** — one sequence, and the IR gate passes on one annotation while
  failing on the other at 26.36 %;
* **the group's own *Eleutherine bulbosa*** — the tool-paper species, 133 = 87 + 38 + 8 published
  and 131 = 85 + 38 + 8 by PGA.

The point of the pair fixtures is that the two files of a pair carry **the same sequence**. It is
asserted here rather than assumed: if it ever stopped being true, every "the annotations differ"
conclusion below would be measuring something else.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from plastid.cli import main
from plastid.core.genbank import read_genbank
from plastid.core.gene_names import aliases_for, cds_by_canonical_name
from plastid.core.ir import IRGateError, detect_ir
from plastid.core.model import count_genes

# The measured corpus fixtures: file → (accession, length, GC %, CDS, tRNA, rRNA).
# Every number was measured on these files (PROJECT.md §8.31); none is a copy of a published
# table, and the two annotations of one genome deliberately sit next to each other.
CORPUS: dict[str, tuple[str, int, float, int, int, int]] = {
    "corpus/Iris_domestica_NC_050833.1.pga.gb": (
        "Iris_domestica_NC_050833",
        153_724,
        37.8536,
        85,
        38,
        8,
    ),
    "corpus/Iris_hermona_NC_081931.1.published.gb": (
        "NC_081931.1",
        153_185,
        37.8934,
        86,
        38,
        8,
    ),
    "corpus/Iris_hermona_NC_081931.1.pga.gb": (
        "Iris_hermona_voucher_IrisSV1903_NC_081931",
        153_185,
        37.8934,
        84,
        38,
        8,
    ),
    "corpus/Eleutherine_bulbosa_PX516262.1.published.gb": (
        "PX516262.1",
        152_471,
        37.8505,
        87,
        38,
        8,
    ),
    "corpus/Eleutherine_bulbosa_PX516262.1.pga.gb": (
        "Eleutherine_bulbosa",
        152_471,
        37.8505,
        85,
        38,
        8,
    ),
}


def _gc(sequence: str) -> float:
    return round(100.0 * sum(1 for base in sequence if base in "GC") / len(sequence), 4)


@pytest.mark.parametrize("relative", sorted(CORPUS))
def test_every_vendored_corpus_record_measures_what_was_recorded(
    data_dir: Path, relative: str
) -> None:
    """Each fixture's own numbers, measured when it was vendored (PROJECT.md §8.31)."""
    accession, length, gc, cds, trna, rrna = CORPUS[relative]
    genome = read_genbank(data_dir / relative)
    counts = count_genes(genome)
    assert genome.accession == accession
    assert len(genome.sequence) == length
    assert _gc(genome.sequence) == gc
    assert (counts.cds, counts.trna, counts.rrna) == (cds, trna, rrna)


def test_the_corpus_pairs_carry_identical_sequences(data_dir: Path) -> None:
    """The fixtures' whole point: two annotations of one sequence, never two sequences.

    Measured on the full corpus as well — 34/34 pairs identical (PROJECT.md §8.31) — but this is
    what makes the fixture pair a controlled comparison, so it is asserted where it is used.
    """
    pairs = [
        ("corpus/Iris_domestica_NC_050833.1.pga.gb", "NC_050833.1.gb"),
        (
            "corpus/Iris_hermona_NC_081931.1.published.gb",
            "corpus/Iris_hermona_NC_081931.1.pga.gb",
        ),
        (
            "corpus/Eleutherine_bulbosa_PX516262.1.published.gb",
            "corpus/Eleutherine_bulbosa_PX516262.1.pga.gb",
        ),
    ]
    for first, second in pairs:
        assert read_genbank(data_dir / first).sequence == read_genbank(data_dir / second).sequence


def test_the_golden_genome_re_annotated_measures_one_cds_short(data_dir: Path) -> None:
    """§6 says 132 = 86 + 38 + 8. PGA's annotation of the same 153,724 bp says **85 CDS**.

    The missing feature is the second `rps12` copy, and that is the whole finding: the contract
    is not a property of the sequence, it is a property of the annotation, and a pipeline that
    quietly accepted either would report two different gene counts for one genome.
    """
    published = read_genbank(data_dir / "NC_050833.1.gb")
    pga = read_genbank(data_dir / "corpus/Iris_domestica_NC_050833.1.pga.gb")

    assert count_genes(published).cds == 86
    assert count_genes(pga).cds == 85

    def rps12_copies(genome) -> int:
        return sum(1 for gene in genome.genes if gene.kind == "CDS" and gene.name == "rps12")

    assert rps12_copies(published) == 2
    assert rps12_copies(pga) == 1


def test_one_sequence_two_verdicts_from_the_ir_gate(data_dir: Path) -> None:
    """*Iris hermona*: the published annotation passes the gate, PGA's fails it at 99.86 %.

    Both files hold the same 153,185 bp, and de-novo detection finds a 26,132 bp inverted repeat
    at 100.00 % on it — so the failure belongs to the annotation, not to the sequence. This is the
    case the gate exists for (PROJECT.md §9.A), and it is the regression fixture for it: 3 of the
    corpus's 34 PGA records behave this way.

    **The number moved on 17/10/2026, the verdict did not.** It read 26.36 % until H-25 was fixed:
    the two annotated copies differ in length by 2, and the mirror constant assumed they did not,
    so the measurement paired the wrong bases and came back at the identity of two random
    sequences. Measured with the alignment the difference implies, PGA's copies agree over
    26,132 bp and differ in a 48 bp block at JLB — 36 mismatches — which is the honest reading:
    *PGA's repeat is 48 bp too wide*, not *the copies are unrelated*. 99.86 % is still under the
    99.9 % gate, so the record is still refused, which is what this fixture is for.
    """
    published = read_genbank(data_dir / "corpus/Iris_hermona_NC_081931.1.published.gb")
    pga = read_genbank(data_dir / "corpus/Iris_hermona_NC_081931.1.pga.gb")

    good = detect_ir(published)
    assert (good.length, good.identity) == (26_132, 1.0)
    assert good.excess_at == "none"

    bad = detect_ir(pga)
    assert bad.length == 26_180
    assert bad.identity == pytest.approx(0.998625, abs=5e-7)
    assert bad.mismatches == 36
    assert bad.excess_at == "outer"

    with pytest.raises(IRGateError):
        from plastid.core.ir import enforce_ir_gate

        enforce_ir_gate(bad)


def test_the_group_own_record_shows_the_same_loss(data_dir: Path) -> None:
    """*E. bulbosa* PX516262.1 — the tool paper's species, and the group's own submission.

    133 = 87 + 38 + 8 as published, 131 = 85 + 38 + 8 by PGA: the same one-copy-short pattern the
    golden genome shows, on the record the team actually cares about. Note the published 87 CDS —
    one more than the golden's 86 — which is why no check may hard-code "86 CDS" as a constant.
    """
    published = read_genbank(data_dir / "corpus/Eleutherine_bulbosa_PX516262.1.published.gb")
    pga = read_genbank(data_dir / "corpus/Eleutherine_bulbosa_PX516262.1.pga.gb")

    assert (count_genes(published).cds, count_genes(pga).cds) == (87, 85)
    assert len(published.sequence) == len(pga.sequence) == 152_471


def test_the_ndh_locus_resolves_to_one_key_under_both_spellings(data_dir: Path) -> None:
    """`ndhL` (published) and `ndhI` (PGA) are one locus, and the alias table already says so.

    Measured on the golden genome: both annotations place it at `118156-118698(-)` with the
    identical 543 bp sequence — the same evidence §8.17 recorded on 03/10, re-confirmed here on
    a second corpus and 8/34 pairs (PROJECT.md §8.31). The check is that the *code* joins them:
    a comparison spanning the two annotations must not see two different genes.
    """
    published = read_genbank(data_dir / "NC_050833.1.gb")
    pga = read_genbank(data_dir / "corpus/Iris_domestica_NC_050833.1.pga.gb")

    assert "ndhL" in {gene.name for gene in published.genes}
    assert "ndhI" in {gene.name for gene in pga.genes}
    # and the canonical key is the naming reference's spelling since D5-09 §8 — the published
    # record is the recorded standard, so its own spelling leads

    for genome in (published, pga):
        locus = aliases_for(genome).locus("ndhL")
        assert locus is not None
        assert "ndhI" in locus.members and "ndhL" in locus.members
        assert cds_by_canonical_name(genome, None).get("ndhL") is not None


def test_the_corpus_files_are_pinned_in_sha256sums(data_dir: Path) -> None:
    """A vendored file that is not in SHA256SUMS is a file nobody can verify.

    Checked over the whole golden directory, not just the corpus: the same rule that makes the
    golden numbers reproducible applies to every file under `tests/data`.
    """
    sums = {}
    for line in (data_dir / "SHA256SUMS").read_text(encoding="utf-8").splitlines():
        if line.strip():
            digest, name = line.split(maxsplit=1)
            sums[name.strip()] = digest

    exempt = {"SHA256SUMS", "PROVENANCE.md"}
    present = {
        str(path.relative_to(data_dir))
        for path in data_dir.rglob("*")
        if path.is_file() and "__pycache__" not in path.parts
    }
    assert present - exempt == set(sums), "every vendored file must be listed, and nothing else"
    for name, digest in sorted(sums.items()):
        assert hashlib.sha256((data_dir / name).read_bytes()).hexdigest() == digest, name


def test_the_hermona_pair_runs_clean_and_stops_loudly(tmp_path: Path, data_dir: Path) -> None:
    """Through the CLI: exit 0 for the annotation that passes, exit 2 for the one that fails.

    The gate is what stands between a diverged annotation and every position-based module, so the
    fixture is worth exercising the way a user would meet it — as an exit code, with no module
    output written after a failed gate.
    """
    good = data_dir / "corpus/Iris_hermona_NC_081931.1.published.gb"
    bad = data_dir / "corpus/Iris_hermona_NC_081931.1.pga.gb"

    def run(source: Path, dataset: str) -> int:
        return main(
            [
                "run",
                "--module",
                "genome_stats",
                "--input",
                str(source),
                "--out",
                str(tmp_path),
                "--dataset",
                dataset,
            ]
        )

    assert run(good, "hermona_published") == 0
    assert run(bad, "hermona_pga") == 2

    failed = sorted((tmp_path / "hermona_pga").iterdir())[0]
    assert not (failed / "01_genome_stats").exists(), "no module output after a failed gate"


# --------------------------------------------------------------- the corpus at full scale

CORPUS_PANEL = Path(__file__).resolve().parents[2] / ".panel" / "iridaceae-34"
CORPUS_RECORDS = sorted(CORPUS_PANEL.glob("*.gb")) if CORPUS_PANEL.is_dir() else []


@pytest.mark.skipif(
    not CORPUS_RECORDS,
    reason="the pinned Iridaceae panel is not in this working tree (.panel/ is gitignored)",
)
def test_a_one_page_pdf_is_produced_for_every_corpus_record(tmp_path: Path) -> None:
    """P5-04's Done-when: 34 records, not one — a page per genome, none failing.

    The page is the PDF-only figure, so the run itself does the analysis and the figure pass
    renders the page from the results the run stored. Every record is checked, and a failure
    names the record rather than the loop.
    """
    import importlib.util

    from plastid.run.figures import render_run_figures

    if importlib.util.find_spec("matplotlib") is None:
        pytest.skip("the 'figures' extra is not installed — pip install 'plastid-suite[figures]'")

    failures: list[str] = []
    pages = 0
    for source in CORPUS_RECORDS:
        code = main(
            [
                "run",
                "--module",
                "genome_stats",
                "--module",
                "structure",
                "--input",
                str(source),
                "--out",
                str(tmp_path),
                "--dataset",
                "corpus-page",
            ]
        )
        if code != 0:
            failures.append(f"{source.name}: run exited {code}")
            continue
        run_dir = sorted((tmp_path / "corpus-page").iterdir())[-1]
        outcome = render_run_figures(run_dir, figures=["genome_page"])
        written = [path for path in outcome.files if path.name == "genome_page.pdf"]
        if len(written) != 1:
            failures.append(f"{source.name}: {outcome.files or outcome.skipped}")
            continue
        if b"/Type /Page" not in written[0].read_bytes():
            failures.append(f"{source.name}: the PDF has no page object")
            continue
        pages += 1

    assert not failures, "the one-page PDF failed on:\n" + "\n".join(failures)
    assert pages == len(CORPUS_RECORDS) == 34, f"{pages} pages for {len(CORPUS_RECORDS)} records"
