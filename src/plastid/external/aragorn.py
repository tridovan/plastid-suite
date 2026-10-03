"""ARAGORN adapter — the tRNA caller of the mito annotation (P3A-04, decision D3A-01).

The mito annotator takes its CDS from BLAST+ hits against a reference protein set and its tRNAs
from here. ARAGORN searches both strands for the canonical tRNA cloverleaf and scores what it
finds; it is also the caller the heavier annotators run inside, which is why it was chosen over a
hand-written structural search — a caller written here and tuned until it returned the published
count would be overfitting with a straight face (TASKS.md → D3A-01).

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
from dataclasses import dataclass
from pathlib import Path

from plastid.core.gene_names import normalize_name
from plastid.core.model import Gene, Genome, Interval
from plastid.core.submission import to_fasta
from plastid.run.context import RunContext
from plastid.run.tools import ToolError, which

__all__ = [
    "ARAGORN_VERSION",
    "AMINO_ACIDS",
    "AragornResult",
    "TrnaHit",
    "aragorn_binary",
    "aragorn_version",
    "parse_batch",
    "run_aragorn",
    "trna_name",
]

ARAGORN_VERSION = "1.2.41"

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
_LINE = re.compile(
    r"^\s*(?P<index>\d+)\s+(?P<kind>tRNA|tmRNA)(?:-(?P<amino>\S+))?\s+"
    r"(?P<strand>c?)\[(?P<start>\d+),(?P<end>\d+)\]\s+(?P<score>\d+)?\s*\((?P<anticodon>\w+)\)\s*$"
)
_DECLARED = re.compile(r"^\s*(?P<count>\d+)\s+genes? found\s*$", re.MULTILINE)
_VERSION = re.compile(r"ARAGORN\s+v(\d+\.\d+\.\d+)", re.IGNORECASE)


@dataclass(frozen=True, slots=True)
class TrnaHit:
    """One tRNA, as ARAGORN reported it and as this toolkit names it. 1-based inclusive."""

    name: str  # trnM-CAU — through the shared normalizer
    source_name: str  # tRNA-Met — what ARAGORN printed
    anticodon: str  # CAU — RNA, 5'→3' on the gene
    start: int
    end: int
    strand: int
    score: int | None
    wrapped: bool = False  # crosses the origin: ARAGORN prints start > end

    def parts(self, genome_length: int) -> tuple[Interval, ...]:
        """The hit as genome intervals, a wrapped one split into the two parts a join needs."""
        if self.wrapped:
            return (
                Interval.from_1based(self.start, genome_length, self.strand),
                Interval.from_1based(1, self.end, self.strand),
            )
        return (Interval.from_1based(self.start, self.end, self.strand),)

    def gene(self, genome_length: int) -> Gene:
        return Gene(name=self.name, kind="tRNA", parts=self.parts(genome_length))


@dataclass(frozen=True, slots=True)
class AragornResult:
    """The tRNAs, the genes they make, and what the run printed — raw evidence included."""

    hits: tuple[TrnaHit, ...]
    genes: tuple[Gene, ...]
    declared: int | None  # ARAGORN's own "N genes found" line
    others: tuple[str, ...]  # gene lines that were not tRNA (a tmRNA, if one is ever found)
    raw: str
    version: str


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


def parse_batch(text: str) -> tuple[tuple[TrnaHit, ...], tuple[str, ...]]:
    """Read ARAGORN's ``-w`` output into hits, plus any gene line that is not a tRNA."""
    hits: list[TrnaHit] = []
    others: list[str] = []
    for line in text.splitlines():
        match = _LINE.match(line)
        if match is None:
            continue
        if match.group("kind") != "tRNA" or match.group("amino") is None:
            # a tmRNA, or a tRNA line whose shape this reader does not know: reported, never
            # guessed at — an unparsed gene that vanished would look like a gene that is absent
            others.append(line.strip())
            continue
        start, end = int(match.group("start")), int(match.group("end"))
        hits.append(
            TrnaHit(
                name=trna_name(match.group("amino"), match.group("anticodon")),
                source_name=f"tRNA-{match.group('amino')}",
                anticodon=match.group("anticodon").upper().replace("T", "U"),
                start=start,
                end=end,
                strand=-1 if match.group("strand") else 1,
                score=int(match.group("score")) if match.group("score") else None,
                wrapped=start > end,
            )
        )
    return tuple(hits), tuple(others)


def run_aragorn(
    genome: Genome,
    ctx: RunContext,
    *,
    seq_id: str,
    aragorn: str | Path | None = None,
    directory: str | Path | None = None,
    timeout: int = 600,
) -> AragornResult:
    """Call the tRNA genes of one genome, and read them back as :class:`Gene` records."""
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
    )
