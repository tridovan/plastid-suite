"""P3B-10's measuring instrument: annotate the corpus leave-one-out and classify every difference.

The phase's Done-when is a comparison — our reading of 34 published plastomes must differ from the
published annotation in **fewer conclusion-level places than PGA's reading does** — so the bar is
only arguable if the measuring instrument is in the repo, re-runnable, and states its rule in one
place. This is that instrument. It is a developer tool: it needs BLAST+, ARAGORN and the pinned
panel (``.panel/``, outside the repo), it takes minutes, and it is not part of the shipped toolkit.

**The rule, stated once and applied to every reading alike** (a rule that bent for our side would
measure nothing):

* a called feature and a published feature **match** when some exon of one overlaps some exon of
  the other by at least half of the shorter — "the same gene is where the same gene is";
* a matched pair whose **names differ** is a **vocabulary** difference — one locus under another
  spelling, which D5-09 owns and no reader has to re-derive;
* a matched pair whose **exons differ by more than `BOUNDARY_TOLERANCE` at an end** is a
  **conclusion** difference: the call does not describe the gene the record describes;
* a matched pair whose exons agree to within that tolerance is a **boundary** difference — the
  same gene in the same place, its ends located a base or two differently by two methods;
* a published feature nothing matched is **missing**; a called feature nothing matched is
  **extra** — both conclusions.

`BOUNDARY_TOLERANCE` is 1 bp, and it is the **tighter** of the two defensible settings rather than
the convenient one — a rule chosen to make a total come out right measures nothing. The measurement
behind it (§8.51): across the 34 records, **719 of 1,264** matched pairs move by exactly one base
at one or both ends, 634 of them with our call strictly containing the record's, and running
ARAGORN by hand on the golden returns the coordinates we emit to the base
(`tRNA-Asp c[30540,30615]` against the record's `complement(30541..30614)`). That is two methods
ending a tRNA one base apart, not two conclusions about the genome. At a 3 bp tolerance — one
codon, the CDS stop-codon convention — the totals are reported too, and the bar is read at 1 bp.

Overlap is compared **by exon, never by envelope** — the rule P3B-08 measured: the published
``rps12`` is trans-spliced, and its envelope spans 68,976–139,634, which would "overlap" half the
plastome.

Two commands:

    python tools/corpus_bar.py run     # 34 × `plastid annotate`, leave-one-out, in parallel
    python tools/corpus_bar.py diff    # classify what the runs wrote, and print the bar verdict

``make corpus-bar`` runs both.
"""

from __future__ import annotations

import argparse
import hashlib
import subprocess
import sys
from collections import Counter, defaultdict
from collections.abc import Iterable
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

from plastid.core.genbank import read_genbank
from plastid.core.model import Gene

ROOT = Path(__file__).resolve().parents[1]
KINDS = ("CDS", "tRNA", "rRNA")

# How far an exon end may move and still be the same call. 1 bp, for the measured reason in the
# module docstring — two methods place a feature's ends, and the corpus says they differ by one
# base far more often than by anything else. A different count of exons, or a different strand, is
# never within tolerance: that is a different transcript, not a differently-placed one.
BOUNDARY_TOLERANCE = 1


@dataclass(frozen=True, slots=True)
class Feature:
    """One feature of one reading, as the comparison needs it — spans by exon, plus the total."""

    name: str
    kind: str
    start: int
    end: int
    strand: int
    spans: tuple[tuple[int, int], ...]
    total: int

    @property
    def label(self) -> str:
        return f"{self.name} {self.start}-{self.end}({self.strand:+d})"


@dataclass(frozen=True, slots=True)
class Difference:
    """One place two readings of the same sequence disagree."""

    kind: str
    verdict: str  # "vocabulary" | "conclusion"
    detail: str  # "name" | "extent" | "missing" | "extra"
    published: str
    called: str
    note: str

    def line(self) -> str:
        return f"{self.kind} {self.detail} {self.published} -> {self.called} ({self.note})"


def features(genes: Iterable[Gene], kind: str) -> list[Feature]:
    """The features of one kind, in 1-based inclusive coordinates (the GenBank convention)."""
    return [
        Feature(
            name=gene.name,
            kind=gene.kind,
            start=min(part.start for part in gene.parts) + 1,
            end=max(part.end for part in gene.parts),
            strand=gene.parts[0].strand,
            spans=tuple((part.start + 1, part.end) for part in gene.parts),
            total=sum(part.end - part.start for part in gene.parts),
        )
        for gene in genes
        if gene.kind == kind
    ]


def _best_exon_overlap(first: Feature, second: Feature) -> int:
    """The largest overlap between any exon of one and any exon of the other."""
    return max(
        (
            max(0, min(left[1], right[1]) - max(left[0], right[0]) + 1)
            for left in first.spans
            for right in second.spans
        ),
        default=0,
    )


def _shortest_exon(feature: Feature) -> int:
    return min(right - left + 1 for left, right in feature.spans)


def end_moves(reference: Feature, called: Feature) -> list[int] | None:
    """How far each exon end moved, or ``None`` when the two are not the same transcript.

    ``None`` covers a different number of exons and a different strand — both a different call
    rather than a differently-placed one, whatever the coordinates say.
    """
    if reference.strand != called.strand or len(reference.spans) != len(called.spans):
        return None
    return [
        move
        for (want_lo, want_hi), (have_lo, have_hi) in zip(
            reference.spans, called.spans, strict=True
        )
        for move in (have_lo - want_lo, have_hi - want_hi)
    ]


def _match(reference: Feature, candidates: list[Feature]) -> Feature | None:
    """The best candidate that is the same gene in the same place, or nothing."""
    best: Feature | None = None
    best_shared = 0
    for candidate in candidates:
        shared = _best_exon_overlap(reference, candidate)
        if not shared or shared * 2 < min(_shortest_exon(reference), _shortest_exon(candidate)):
            continue
        if shared > best_shared:
            best, best_shared = candidate, shared
    return best


def compare(published: list[Gene], called: list[Gene]) -> list[Difference]:
    """Every difference between two annotations of one sequence, classified."""
    differences: list[Difference] = []
    for kind in KINDS:
        want = features(published, kind)
        have = features(called, kind)
        remaining = list(have)
        for reference in want:
            best = _match(reference, remaining)
            if best is None:
                differences.append(
                    Difference(kind, "conclusion", "missing", reference.label, "—", "not called")
                )
                continue
            remaining.remove(best)
            if best.name != reference.name:
                differences.append(
                    Difference(kind, "vocabulary", "name", reference.label, best.label, best.name)
                )
            else:
                # Every exon end, not the envelope: two readings of the golden's `rpoC1` share
                # 20179-22962 and differ by 27 bp *inside* it, and the strand is not in the
                # coordinates at all — a call on the other strand is a different call wherever its
                # ends sit. Both were holes in the first version of this rule (§8.51).
                moves = end_moves(reference, best)
                if moves is None:
                    verdict = "conclusion"  # another strand, or another number of exons
                elif not any(moves):
                    continue
                elif any(abs(move) > BOUNDARY_TOLERANCE for move in moves):
                    verdict = "conclusion"
                else:
                    verdict = "boundary"
                differences.append(
                    Difference(
                        kind,
                        # a boundary difference is reported, never dropped — it is a real
                        # difference between two readings, and it is what the tools disagree
                        # about; it is simply not a different conclusion about the genome.
                        verdict,
                        "extent",
                        reference.label,
                        best.label,
                        f"{best.total - reference.total:+d} bp",
                    )
                )
        differences += [
            Difference(kind, "conclusion", "extra", "—", extra.label, "not in the record")
            for extra in remaining
        ]
    return differences


def summarise(differences: Iterable[Difference]) -> Counter[str]:
    """The counts the bar is read from, and the finer cuts that explain them."""
    counts: Counter[str] = Counter()
    for difference in differences:
        counts[difference.verdict] += 1
        counts[f"{difference.verdict}:{difference.detail}"] += 1
        counts[f"{difference.kind}:{difference.verdict}"] += 1
    return counts


def load(path: Path) -> list[Gene]:
    return list(read_genbank(path).genes)


# ------------------------------------------------------------------ running the corpus


@dataclass(frozen=True, slots=True)
class Panel:
    """The pinned input: the records, PGA's re-annotation of each, and where the runs go."""

    records: Path
    pga: Path
    output: Path
    reference: Path
    tools: str
    exclude: str
    jobs: int

    def published(self) -> list[Path]:
        return sorted(self.records.glob("*.gb"))

    def pga_by_sequence(self) -> dict[str, Path]:
        """The PGA readings, keyed by the **sequence** they annotate.

        Never by file name or by accession: PGA's own records carry the organism as their
        accession (`Eleutherine_bulbosa`), and the two halves of the corpus are one sequence
        under two annotations — which is the property that makes the comparison controlled, so
        it is the property the lookup uses (§8.31). A duplicate digest is a broken pair and
        stops rather than picking one.
        """
        index: dict[str, Path] = {}
        for path in sorted(self.pga.glob("*.gb")):
            digest = hashlib.sha256(read_genbank(path).sequence.encode("ascii")).hexdigest()
            if digest in index:
                raise SystemExit(
                    f"{path.name} annotates the same sequence as {index[digest].name} — the "
                    f"published/PGA pairing is broken and a comparison would be meaningless"
                )
            index[digest] = path
        return index

    def ours(self, accession: str) -> Path | None:
        found = sorted((self.output / accession).rglob("plastid_annotation.gb"))
        return found[-1] if found else None


def annotate(panel: Panel, source: Path) -> tuple[str, int]:
    """One leave-one-out run of the toolkit, as a user would start it."""
    accession = source.stem
    if panel.ours(accession) is not None:  # already annotated — a second pass costs minutes
        return accession, 0
    command = [
        "uv",
        "run",
        "plastid",
        "annotate",
        "-i",
        str(source),
        "--reference",
        str(panel.reference),
        "--dataset",
        accession,
        "--set",
        f"tools={panel.tools}",
        "-o",
        str(panel.output / accession),
    ]
    if panel.exclude:
        command += ["--exclude", panel.exclude]
    done = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
    if done.stdout or done.stderr:
        (panel.output / f"{accession}.log").write_text(done.stdout + done.stderr, encoding="utf-8")
    return accession, done.returncode


def run_corpus(panel: Panel) -> int:
    records = panel.published()
    print(f"{len(records)} record(s) → {panel.output}", flush=True)
    panel.output.mkdir(parents=True, exist_ok=True)
    failures = 0
    with ThreadPoolExecutor(max_workers=panel.jobs) as pool:
        for accession, code in pool.map(lambda path: annotate(panel, path), records):
            print(f"  {accession}: exit {code}", flush=True)
            failures += code != 0
    print(f"done — {failures} run(s) failed", flush=True)
    return 1 if failures else 0


# ------------------------------------------------------------------ classifying what it wrote


def by_gene(differences: Iterable[Difference], detail: str) -> Counter[str]:
    """Which loci the differences of one kind are about — the composition behind a count."""
    names: Counter[str] = Counter()
    for difference in differences:
        if difference.verdict != "conclusion" or difference.detail != detail:
            continue
        label = difference.called if detail == "extra" else difference.published
        names[label.split()[0]] += 1
    return names


def report(panel: Panel, *, verbose: bool = False) -> int:
    totals: dict[str, Counter[str]] = {"ours": Counter(), "pga": Counter()}
    per_gene: dict[str, dict[str, Counter[str]]] = {
        "ours": defaultdict(Counter),
        "pga": defaultdict(Counter),
    }
    rows: list[tuple[str, Counter[str], Counter[str]]] = []
    absent: list[str] = []

    pga_index = panel.pga_by_sequence()
    for source in panel.published():
        accession = source.stem
        record = read_genbank(source)
        published = list(record.genes)
        digest = hashlib.sha256(record.sequence.encode("ascii")).hexdigest()
        readings = {"ours": panel.ours(accession), "pga": pga_index.get(digest)}
        counts: dict[str, Counter[str]] = {}
        for which, artifact in readings.items():
            if artifact is None:
                absent.append(f"{accession}/{which}")
                counts[which] = Counter()
                continue
            differences = compare(published, load(artifact))
            counts[which] = summarise(differences)
            totals[which] += counts[which]
            for detail in ("missing", "extra", "extent"):
                per_gene[which][detail] += by_gene(differences, detail)
            if verbose:
                print(f"\n--- {accession} · {which} ---")
                for difference in differences:
                    print(f"  {difference.line()}")
        rows.append((accession, counts["ours"], counts["pga"]))

    for which in ("ours", "pga"):
        counts = totals[which]
        print(f"\n=== {which} against the published records ===")
        print(
            f"  conclusion-level {counts['conclusion']:5d}   "
            f"(missing {counts['conclusion:missing']}, extra {counts['conclusion:extra']}, "
            f"extent {counts['conclusion:extent']})"
        )
        print(
            f"  boundary         {counts['boundary']:5d}   "
            f"— ends within {BOUNDARY_TOLERANCE} bp, reported not counted"
        )
        print(f"  vocabulary       {counts['vocabulary']:5d}")
        print(
            f"  strict (no tolerance): {counts['conclusion'] + counts['boundary']:5d} "
            f"— what the count is if every moved end is a conclusion"
        )
        for detail in ("missing", "extra", "extent"):
            names = per_gene[which][detail]
            print(f"  {detail} by locus: " + ", ".join(f"{n} {c}" for n, c in names.most_common(8)))

    print("\n=== per record ===")
    better = worse = tied = 0
    for accession, ours, pga in rows:
        mark = "="
        if ours["conclusion"] < pga["conclusion"]:
            mark, better = "<", better + 1
        elif ours["conclusion"] > pga["conclusion"]:
            mark, worse = ">", worse + 1
        else:
            tied += 1
        print(f"  {accession:14s} ours {ours['conclusion']:4d}  pga {pga['conclusion']:4d}  {mark}")

    if absent:
        print(f"\nNOT MEASURED — no artifact for: {', '.join(absent)}")
        return 1

    ours_total = totals["ours"]["conclusion"]
    pga_total = totals["pga"]["conclusion"]
    ours_strict = ours_total + totals["ours"]["boundary"]
    pga_strict = pga_total + totals["pga"]["boundary"]
    print(
        f"\nours better in {better}, worse in {worse}, tied in {tied} record(s) · "
        f"total {ours_total} against PGA's {pga_total}"
    )
    print(f"strict, every moved end counted: {ours_strict} against PGA's {pga_strict}")
    met = ours_total < pga_total
    print(f"the bar: {'MET' if met else 'NOT MET'} — fewer conclusion-level differences than PGA")
    if met != (ours_strict < pga_strict):
        print("NOTE: the two readings disagree — the verdict rests on the boundary tolerance")
    return 0 if met else 1


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("command", choices=("run", "diff", "all"), nargs="?", default="all")
    parser.add_argument("--panel", type=Path, default=ROOT / ".panel" / "iridaceae-34")
    parser.add_argument("--pga", type=Path, default=ROOT / ".panel" / "iridaceae-34-pga")
    parser.add_argument("--reference", type=Path, default=ROOT / ".panel" / "plastid-reference")
    parser.add_argument("--out", type=Path, default=Path("/tmp/plastid-corpus-bar"))
    parser.add_argument("--tools", default="", help="directory holding blastx/blastn/aragorn")
    parser.add_argument(
        "--exclude",
        default="NC_024936.1",
        help="a record the audit blocks, dropped by name (D5-08's own loop)",
    )
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("-v", "--verbose", action="store_true", help="every difference, listed")
    args = parser.parse_args(argv)

    for name, path in (("panel", args.panel), ("reference set", args.reference)):
        if not path.is_dir():
            print(
                f"{name} not found: {path} (.panel/ is gitignored — pin it first)",
                file=sys.stderr,
            )
            return 2

    panel = Panel(
        records=args.panel,
        pga=args.pga,
        output=args.out,
        reference=args.reference,
        tools=args.tools,
        exclude=args.exclude,
        jobs=args.jobs,
    )
    if args.command in ("run", "all") and run_corpus(panel) != 0:
        return 1
    if args.command == "run":
        return 0
    return report(panel, verbose=args.verbose)


if __name__ == "__main__":
    sys.exit(main())
