"""ARAGORN adapter — the tRNA caller of the mito annotation (P3A-04, decision D3A-01).

The mito annotator takes its CDS from BLAST+ hits against a reference protein set and its tRNAs
from here. ARAGORN searches both strands for the canonical tRNA cloverleaf and scores what it
finds; it is also the caller the heavier annotators run inside, which is why it was chosen over a
hand-written structural search — a caller written here and tuned until it returned the published
count would be overfitting with a straight face (D3A-01).

Three things about its output are handled rather than discovered later:

* **the batch format is the parseable one.** The default output is a secondary-structure drawing;
  ``-w`` prints one line per gene: ``N  tRNA-<Aa>  c?[start,end]  <score>  (<anticodon>)``.
* **it wraps the search around a circular record**, so a gene crossing the origin arrives as
  ``[258274,73]`` — a start greater than the end. Measured on the answer-key record: its first
  tRNA is exactly that case. A wrapped interval is not a GenBank coordinate, so the adapter splits
  it into the two parts a ``join()`` needs, in the order the gene is transcribed.
* **the names are ARAGORN's, not a community convention.** It prints ``tRNA-Met … (cat)``: an
  amino acid in words and the anticodon as DNA. The published form is ``trnM-CAU``, so the adapter
  maps the amino acid to its one-letter code, uppercases the anticodon into RNA, and hands the
  result to :func:`plastid.core.gene_names.normalize_name` — the same normalizer a downloaded
  record goes through. There is no private rename here: an unresolved name is reported as it
  arrived.

Version pinning: ARAGORN has no version flag. ``-V`` is not one — it exits **0** and prints "No
sequence file specified" (measured 03/10/2026), which reads like a version in a manifest. So the
adapter verifies the release itself, from ``-h``'s banner, and passes the verified number to
``run_tool`` as the pinned version, so the run's provenance says 1.2.41 rather than a sentence.
"""

from __future__ import annotations

import re
import tempfile
from dataclasses import dataclass, replace
from pathlib import Path

from plastid.core.gene_names import normalize_name
from plastid.core.model import Gene, Genome, Interval
from plastid.core.spans import TRNA_INTRON_MAX_BP
from plastid.core.submission import to_fasta
from plastid.run.context import RunContext
from plastid.run.tools import ToolError, which

__all__ = [
    "ARAGORN_VERSION",
    "AMINO_ACIDS",
    "AragornResult",
    "TrnaHit",
    "ARAGORN_INTRON_LIMITS",
    "aragorn_binary",
    "aragorn_version",
    "parse_batch",
    "parse_calls",
    "run_aragorn",
    "trna_name",
]

ARAGORN_VERSION = "1.2.41"

# What ARAGORN's own `-h` says about its intron search: "minimum intron length is 0 bases" and
# the default maximum is 3,000 bases (`-i<min>,<max>` changes them). Declared here and carried
# on the result, because a caller that searched 0–3,000 bp must not report "no intron" for a
# locus whose intron is longer: it did not look there (H-22).
#
# Measured 26/10/2026 on ARAGORN 1.2.41 with a golden trnK-UUU whose intron was lengthened in
# steps: 3,000 bp is found, 3,001 bp is not, and `-i0,4000`, `-i4000` and `-i0,6000` do not
# change that — the documented switch does not lift the internal maximum. The minimum does work
# (`-i2700,3000` drops a 2,618 bp intron). A request above this maximum is therefore refused by
# name rather than passed as a flag that silently has no effect.
ARAGORN_INTRON_LIMITS: tuple[int, int] = (0, TRNA_INTRON_MAX_BP)

# One-letter codes for the amino acid names ARAGORN prints. A table, not a heuristic: every entry
# is a standard IUPAC name, and an unknown one is refused rather than guessed at.
AMINO_ACIDS: dict[str, str] = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Gln": "Q",
    "Glu": "E",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Sec": "U",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
    "Pyl": "O",
}

# ``1   tRNA-Phe                   [258274,73]	35  	(gaa)``
# ``8   tRNA-Met              c[107050,107126]	35  	(cat)``
# ``1   tRNA-Lys              c[1367,4058]	34  	(ttt)i(39,2618)``  — with -i, the intron, as
# ``(start, length)``. The group is named for what the field is, because reading the second one
# as an end coordinate is a mistake this parser made for a day and no count could see.
_LINE = re.compile(
    r"^\s*(?P<index>\d+)\s+(?P<kind>tRNA|tmRNA)(?:-(?P<amino>\S+))?\s+"
    r"(?P<strand>c?)\[(?P<start>\d+),(?P<end>\d+)\]\s+(?P<score>\d+)?\s*"
    r"\((?P<anticodon>\w+)\)\s*(?:i\((?P<intron_start>\d+),(?P<intron_length>\d+)\))?\s*$"
)
_DECLARED = re.compile(r"^\s*(?P<count>\d+)\s+genes? found\s*$", re.MULTILINE)
_VERSION = re.compile(r"ARAGORN\s+v(\d+\.\d+\.\d+)", re.IGNORECASE)


@dataclass(frozen=True, slots=True)
class TrnaHit:
    """One tRNA, as ARAGORN reported it and as this toolkit names it. 1-based inclusive.

    ``name`` is ``None`` when ARAGORN printed a gene this adapter cannot name — a ``tRNA-Stop``
    line whose anticodon is one base off, or a ``tmRNA``. The coordinates, strand and intron are
    still real: a later step with locus evidence (``plastid_annotate`` against the reference set's
    own tRNAs) can name the gene, and dropping the line here would turn "unnameable" into
    "absent" (H-24).
    """

    name: str | None  # trnM-CAU — through the shared normalizer; None when unnameable
    source_name: str  # tRNA-Met — what ARAGORN printed
    anticodon: str  # CAU — RNA, 5'→3' on the gene
    start: int
    end: int
    strand: int
    score: int | None
    wrapped: bool = False  # crosses the origin: ARAGORN prints start > end
    # The intron, as ``(start, length)`` **within the gene**, 1-based from its 5' end — ARAGORN's
    # own ``i(39,2618)``. Kept in gene coordinates rather than converted here because the
    # conversion depends on the strand and the origin: see `parts`.
    #
    # The second field is the intron's **length, not its end**. Measured 10/10/2026: read
    # as an end, the golden's eight intron tRNAs spliced into 104–134 nt, which no tRNA is, and
    # the corpus counted the difference as a wrong call on every one of them; read as a length they
    # come out 72–95 nt, and `trnL-UAA` (35 + 50 = 85) and `trnA-UGC` (38 + 35 = 73) land exactly
    # on the published record's own joined exons. ARAGORN's help documents the intron's maximum
    # *length* (``-i<min>,<max>``) and never mentions an end coordinate.
    intron: tuple[int, int] | None = None
    # The line ARAGORN printed, kept for the report when the name is the thing in question.
    raw_line: str = ""

    def with_name(self, name: str) -> TrnaHit:
        """This hit carrying a name established by independent locus evidence."""
        return replace(self, name=name)

    def length(self, genome_length: int) -> int:
        """How long the gene is, intron included — ARAGORN's span is the whole gene."""
        if self.wrapped:
            return (genome_length - self.start + 1) + self.end
        return self.end - self.start + 1

    def _coordinate(self, position: int, genome_length: int) -> int:
        """The genome coordinate of gene position ``position`` (1-based, 5'→3' along the gene).

        The 5' end is ``start`` on the plus strand and ``end`` on the minus strand, so a position
        walks up the genome one way and down the other — and steps past the origin on a circular
        record, which is where the wrap is undone.
        """
        value = self.start + position - 1 if self.strand > 0 else self.end - position + 1
        if value > genome_length:
            return value - genome_length
        if value < 1:
            return value + genome_length
        return value

    def parts(self, genome_length: int) -> tuple[Interval, ...]:
        """The hit as genome intervals **in transcription order** — a join's parts, not one span.

        Two things split a hit, and both would otherwise be written as a claim the sequence does
        not support: an **intron** (ARAGORN's span covers it, and a feature written over the whole
        span would say the intron is coding) and **the origin** (a wrapped interval is not a
        GenBank coordinate). Measured on the golden 09/10/2026: eight of its thirty-eight tRNAs
        carry an intron — ``trnK-UUU``, ``trnL-UAA``, ``trnV-UAC``, one of the two ``trnG``, and
        both copies of ``trnI-GAU`` and ``trnA-UGC`` — and the published record writes every one
        of them as a two-part join.
        """
        if self.intron is None:
            spans = ((1, self.length(genome_length)),)
        else:
            # the intron covers `start` .. `start + length - 1`, so the second exon opens after it
            start, length = self.intron
            spans = ((1, start - 1), (start + length, self.length(genome_length)))
        runs: list[list[int]] = []
        for first, last in spans:
            for position in range(first, last + 1):
                coordinate = self._coordinate(position, genome_length)
                if runs and abs(coordinate - runs[-1][-1]) == 1:
                    runs[-1].append(coordinate)
                else:
                    runs.append([coordinate])
        return tuple(Interval.from_1based(min(run), max(run), self.strand) for run in runs if run)

    def gene(self, genome_length: int) -> Gene:
        if self.name is None:
            raise ValueError(
                f"{self.source_name} at {self.start}..{self.end} has no name — a gene the adapter "
                f"cannot name is reported, not written under a guessed symbol"
            )
        return Gene(name=self.name, kind="tRNA", parts=self.parts(genome_length))


@dataclass(frozen=True, slots=True)
class AragornResult:
    """The tRNAs, the genes they make, and what the run printed — raw evidence included."""

    hits: tuple[TrnaHit, ...]
    genes: tuple[Gene, ...]
    declared: int | None  # ARAGORN's own "N genes found" line
    others: tuple[str, ...]  # gene lines that were not tRNA, or that this adapter cannot name
    raw: str
    version: str
    # Every tRNA line with its coordinates, named or not. The plastid path names an unnameable
    # intron tRNA from the reference set's own locus (H-24); the mito path ignores this.
    calls: tuple[TrnaHit, ...] = ()
    # The intron length range this call searched, or None when intron detection was off. A
    # consumer must read it: absence of a tRNA intron is only a measurement inside this range.
    intron_search: tuple[int, int] | None = None


def aragorn_binary(explicit: str | Path | None = None, directory: str | Path | None = None) -> str:
    """Locate the binary: an explicit path, then a directory, then ``PATH``."""
    if explicit is not None:
        path = Path(explicit)
        if not path.is_file():
            raise ToolError(f"aragorn not found at {path} — see docs/tools/aragorn.md")
        return str(path)
    if directory is not None:
        candidate = Path(directory) / "aragorn"
        if not candidate.is_file():
            raise ToolError(f"aragorn not found in {directory} — see docs/tools/aragorn.md")
        return str(candidate)
    found = which("aragorn")
    if found is None:
        raise ToolError(
            f"aragorn not found — install the pinned ARAGORN {ARAGORN_VERSION} or point the run "
            f"at one with --set tools=<directory>. See docs/tools/aragorn.md."
        )
    return found


def aragorn_version(binary: str, *, timeout: int = 60) -> str:
    """The release this binary is, read from ``-h``'s banner and checked against the pin.

    ARAGORN has no usable version flag: ``-V`` exits 0 and prints an error sentence. ``-h`` is
    safe (no input file, exit 0) and names the release, so that is what is asked.
    """
    import subprocess

    try:
        result = subprocess.run(
            [binary, "-h"], capture_output=True, text=True, timeout=timeout, check=False
        )
    except OSError as error:
        raise ToolError(f"aragorn: cannot run {binary}: {error}") from error
    match = _VERSION.search(result.stdout or "")
    if match is None:
        raise ToolError(
            f"{binary}: no 'ARAGORN v…' line in its help output — a step whose version nobody "
            f"can name is not part of the toolkit (docs/tools/aragorn.md)"
        )
    version = match.group(1)
    if version != ARAGORN_VERSION:
        raise ToolError(
            f"{binary} is ARAGORN {version}, but the toolkit's numbers were measured with "
            f"{ARAGORN_VERSION} (docs/tools/aragorn.md) — a different caller is a different "
            f"tRNA set. Install the pinned build, or re-measure and update ARAGORN_VERSION."
        )
    return version


def trna_name(amino: str, anticodon: str) -> str:
    """``("Met", "cat")`` → ``trnM-CAU``, through the shared normalizer.

    The anticodon is printed as DNA and is already in the gene's own 5'→3' direction, so the only
    conversion is T→U; the amino acid goes from ARAGORN's word to the one-letter code. The result
    is then normalized exactly as a name arriving from a GenBank record would be.
    """
    letter = AMINO_ACIDS.get(amino.capitalize())
    if letter is None:
        raise ValueError(
            f"ARAGORN reported the amino acid {amino!r}, which is not one of the standard names "
            f"this adapter knows ({sorted(AMINO_ACIDS)}) — refusing to guess a residue letter"
        )
    rna = anticodon.upper().replace("T", "U")
    return normalize_name(f"trn{letter}-{rna}").canonical


def _call_from_line(line: str) -> TrnaHit | None:
    """One ARAGORN gene line as a hit, named when the amino acid is one this adapter knows.

    ``None`` means the line is not a parseable tRNA/tmRNA gene line at all; a line that parses but
    whose amino acid is unknown comes back with ``name=None`` and keeps its coordinates.
    """
    match = _LINE.match(line)
    if match is None or match.group("kind") != "tRNA":
        return None
    amino = match.group("amino")
    name = (
        trna_name(amino, match.group("anticodon"))
        if amino is not None and amino in AMINO_ACIDS
        else None
    )
    start, end = int(match.group("start")), int(match.group("end"))
    intron_start = match.group("intron_start")
    return TrnaHit(
        name=name,
        source_name=f"tRNA-{amino}" if amino else "tRNA",
        anticodon=match.group("anticodon").upper().replace("T", "U"),
        start=start,
        end=end,
        strand=-1 if match.group("strand") else 1,
        score=int(match.group("score")) if match.group("score") else None,
        wrapped=start > end,
        intron=(
            (int(intron_start), int(match.group("intron_length")))
            if intron_start is not None
            else None
        ),
        raw_line=line.strip(),
    )


def parse_calls(text: str) -> tuple[TrnaHit, ...]:
    """Every tRNA gene line, named or not — the input locus evidence can act on (H-24)."""
    calls = (_call_from_line(line) for line in text.splitlines())
    return tuple(call for call in calls if call is not None)


def parse_batch(text: str) -> tuple[tuple[TrnaHit, ...], tuple[str, ...]]:
    """Read ARAGORN's ``-w`` output into named hits, plus every gene line that was not named."""
    hits: list[TrnaHit] = []
    others: list[str] = []
    for line in text.splitlines():
        if _LINE.match(line) is None:
            continue
        call = _call_from_line(line)
        if call is None or call.name is None:
            # A tmRNA, or a tRNA whose amino acid is not one of the twenty-two. Reported through
            # `others`, never dropped and never guessed at: an unparsed gene that vanished would
            # look like a gene that is absent. Measured 09/10/2026 on the pinned corpus (with
            # `-i`): two records (`NC_072561.1` · `NC_072580.1`) bring a `tRNA-Stop` line at the
            # `trnK-UUU` locus — ARAGORN reading the anticodon one base off, so that it spells a
            # stop codon. The gene is real; the *name* is what this adapter cannot supply.
            others.append(line.strip())
            continue
        hits.append(call)
    return tuple(hits), tuple(others)


def run_aragorn(
    genome: Genome,
    ctx: RunContext,
    *,
    seq_id: str,
    aragorn: str | Path | None = None,
    directory: str | Path | None = None,
    timeout: int = 600,
    introns: bool = False,
    intron_limits: tuple[int, int] | None = None,
) -> AragornResult:
    """Call the tRNA genes of one genome, and read them back as :class:`Gene` records.

    ``intron_limits`` is ARAGORN's ``-i<min>,<max>`` range. The default (``introns=True``) is the
    caller's documented 0–3,000 bp. Passing an explicit range turns intron detection on: a caller
    that widens it says so, and :attr:`AragornResult.intron_search` carries what was searched so
    downstream code cannot report a locus outside the range as "no intron" (H-22).
    """
    if intron_limits is not None:
        if not introns:
            introns = True
        minimum, maximum = intron_limits
        if minimum < 0 or maximum <= minimum:
            raise ValueError(
                f"intron_limits must be (min, max) with 0 <= min < max, got {intron_limits!r}"
            )
        if maximum > ARAGORN_INTRON_LIMITS[1]:
            raise ValueError(
                f"ARAGORN {ARAGORN_VERSION} does not search a tRNA intron longer than "
                f"{ARAGORN_INTRON_LIMITS[1]:,} bp: measured 26/10/2026, a 3,000 bp intron is "
                f"reported and a 3,001 bp one is not, even with -i0,{maximum}. A locus beyond "
                f"that range is not assessed, never counted as intron-free. Use "
                f"intron_limits=(min, {ARAGORN_INTRON_LIMITS[1]}) or a caller that reports "
                f"longer introns."
            )
    searched = (intron_limits or ARAGORN_INTRON_LIMITS) if introns else None
    binary = aragorn_binary(aragorn, directory)
    version = aragorn_version(binary)
    with tempfile.TemporaryDirectory(prefix="aragorn_") as work:
        workdir = Path(work)
        fasta = workdir / "target.fna"
        fasta.write_text(to_fasta(genome, seq_id), encoding="utf-8")
        record = ctx.tool(
            "aragorn",
            [
                binary,
                "-w",  # batch mode: one line per gene, instead of a structure drawing
                # tRNA genes **with introns**, when the caller asks for them. ARAGORN's own help
                # states it plainly — "Unless -i is specified, tRNA genes containing introns are
                # not detected" — and it says nothing when it skips them, which makes this the
                # same silent-shortfall shape as the truncating `-max_target_seqs`.
                #
                # It is the caller's decision because each organelle path measures its own, and
                # the two measurements disagree (09/10/2026, P3B-10): on the plastid golden the
                # flag turns **30 into 38, exactly the published record's thirty-eight** (eight of
                # its tRNAs carry an intron), while on the mito answer key it turns 17 into 21,
                # four of which the published record does not have. The plastid path passes it;
                # whether the mito path should is a question with a name, not a default.
                *([] if searched is None else [_intron_flag(searched)]),
                "-gc1",  # the standard code — the same table the CDS calls carry
                "-c" if genome.topology == "circular" else "-l",
                str(fasta),
            ],
            kind="binary",
            tool_path=binary,
            # the release this adapter verified itself, so the run's provenance names a version
            # rather than the error sentence ARAGORN answers -V with
            pinned_version=version,
            cwd=workdir,
            timeout=timeout,
        )
        raw = Path(record.stdout_file).read_text(encoding="utf-8")

    hits, others = parse_batch(raw)
    declared_match = _DECLARED.search(raw)
    return AragornResult(
        hits=hits,
        genes=tuple(hit.gene(len(genome)) for hit in hits),
        declared=int(declared_match.group("count")) if declared_match else None,
        others=others,
        raw=raw,
        version=version,
        calls=parse_calls(raw),
        intron_search=searched,
    )


def _intron_flag(limits: tuple[int, int]) -> str:
    """ARAGORN's intron switch for one range — plain ``-i`` for its own documented default."""
    if limits == ARAGORN_INTRON_LIMITS:
        return "-i"
    minimum, maximum = limits
    return f"-i{minimum},{maximum}"
