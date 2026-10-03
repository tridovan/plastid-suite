"""Chloe adapter — plastome annotation (M3.1, P3-03; decision D3-01).

Chloe (Small lab, UWA) annotates a flowering-plant plastome offline against its own curated
reference set. It replaces PGA as this pipeline's annotator: PGA returned 147 genes where the
golden record has 132, in three silently-failing ways (PROJECT.md §8.15); Chloe returns
**132 = 86 + 38 + 8** with every CDS translating cleanly (measured 2026-10-03, §8.17).

Three of Chloe's habits are handled here rather than discovered later:

* **it transforms the sequence by default.** Its own help says ``--no-transform  do not flip and
  orient sequence to standard configuration``, and with the default on, the GFF's coordinates
  refer to a sequence the caller never had. Measured on the golden genome: the transform
  reverse-complements the 18,169 bp SSC (13,001 of its positions differ), and the transformed GFF
  puts ndhH at 114,734-115,915 (+) — a region whose forward translation is nothing like ndhH.
  The gene *set* is identical either way (same 86 CDS / 38 tRNA / 8 rRNA names); only the frame
  differs, which is exactly the kind of difference a gene count cannot see. So ``--no-transform``
  is passed, and every run re-reads Chloe's own FASTA when it writes one and refuses if it
  differs from the input.
* **it names itself in its output, not on the command line** — ``##source-version Chloe 1.0``.
  A script is never probed with a flag here (PGA's lesson: an unknown flag makes the script run
  itself with its defaults), so both numbers are read from files instead and checked against their
  pins: the release from the checkout's ``Project.toml`` (0.1.18 — the one the Zenodo deposit
  publishes) and the GFF's own declaration, which says which output format this reader is reading.
* **its reference set is an artifact of its own** (``Chloe_refs``, pinned in the checkout's
  ``Artifacts.toml``). An annotation is only as reproducible as the references behind it, so that
  tree hash is checked too — the same rule that keeps an unpinned mafft out of the toolkit.

Scope, recorded because a paper has to state it: Chloe annotates **angiosperm chloroplasts only**
(its README says so, and offers to contact the authors about other plants), and its checkout
carries **no LICENSE file** — both noted for the M9.3 licence review.
"""

from __future__ import annotations

import re
import tempfile
from dataclasses import dataclass, replace
from pathlib import Path

from plastid.core.gff import genes_from_gff
from plastid.core.model import GENE_KINDS, Genome
from plastid.core.submission import to_fasta
from plastid.run.context import RunContext
from plastid.run.tools import ToolError, tool_version, which

__all__ = [
    "CHLOE_GFF_VERSION",
    "CHLOE_REFS_TREE",
    "CHLOE_VERSION",
    "INSTALL_LINE",
    "ChloeResult",
    "chloe_script",
    "julia_binary",
    "package_version",
    "reference_set",
    "run_chloe",
]

TOOL_NAME = "chloe"

# The released software, read from the checkout's Project.toml and checked against this pin. It is
# the version the citable artifact names (Zenodo, doi:10.5281/zenodo.18194040).
CHLOE_VERSION = "0.1.18"

# What the annotator writes into every GFF (``##source-version Chloe 1.0``) — a second, different
# number, and the one that says which output format we are reading. Both are checked: a reader
# that silently accepted either changing would be reading a file it has never seen.
CHLOE_GFF_VERSION = "1.0"

# The reference set the measured numbers were produced with: ``Chloe_refs`` v1.0.1, by its git
# tree hash in the checkout's Artifacts.toml. A different hash is a different reference set, and
# a gene count from a reference set nobody can name is not a result.
CHLOE_REFS_TREE = "9163915d99fdfbc7d8510f0479b3d08905682ca5"

INSTALL_LINE = (
    "Chloe is not a package: install Julia (https://julialang.org/downloads/), then "
    "`git clone https://github.com/ian-small/chloe` and instantiate it with "
    "`julia --project=. -e 'using Pkg; Pkg.instantiate()'` inside the clone — the curated "
    "reference set downloads once, after which the annotator runs offline. Then pass "
    "--set chloe=/path/to/chloe and, if julia is not on PATH, --set julia=/path/to/julia"
)

_ARTIFACT_SECTION = "Chloe_refs"
_GFF_VERSION_LINE = re.compile(r"^##source-version\s+Chloe\s+(\S+)\s*$", re.MULTILINE)
_PACKAGE_VERSION_LINE = re.compile(r'^version\s*=\s*"([^"]+)"', re.MULTILINE)


@dataclass(frozen=True, slots=True)
class ChloeResult:
    """The annotated genome, the GFF text it was read from, and where both came from.

    The text travels back rather than a path because the run happens in a temporary directory
    that is gone by the time a module sees this — the caller decides what to keep (a module hands
    it to :class:`~plastid.result.result.Result` as an artifact, seam 2).
    """

    genome: Genome
    gff_text: str
    source_name: str
    version: str  # the Chloe release the annotation came from (CHLOE_VERSION)
    reference_tree: str


def chloe_script(explicit: str | Path | None = None) -> Path:
    """Locate ``chloe.jl`` — an explicit path, then ``PATH``. A missing one is a named error."""
    if explicit is not None:
        path = Path(explicit)
        if path.is_dir():
            path = path / "chloe.jl"
        if not path.is_file():
            raise ToolError(f"chloe.jl not found at {path} — pass --set chloe=/path/to/chloe")
        return path
    found = which("chloe.jl")
    if found is None:
        raise ToolError(f"chloe.jl not found — {INSTALL_LINE}")
    return Path(found)


def julia_binary(explicit: str | Path | None = None) -> str:
    """Locate the Julia interpreter — an explicit path, then ``PATH``."""
    if explicit is not None:
        path = Path(explicit)
        if not path.is_file():
            raise ToolError(f"julia not found at {path} — pass --set julia=/path/to/julia")
        return str(path)
    found = which("julia")
    if found is None:
        raise ToolError(f"julia not found — {INSTALL_LINE}")
    return found


def package_version(script: Path) -> str:
    """The Chloe release this checkout is — read from ``Project.toml``, not assumed.

    The version that identifies the software is the Julia package's (0.1.18), not the number the
    annotator prints in its GFF (see :data:`CHLOE_GFF_VERSION`): it is the one the Zenodo deposit
    is published under, so it is the one a Materials & Methods section can cite.
    """
    manifest = script.parent / "Project.toml"
    if not manifest.is_file():
        raise ToolError(
            f"chloe: {manifest} is missing, so the annotator's version cannot be identified — "
            f"check that --set chloe= points at a Chloe clone"
        )
    match = _PACKAGE_VERSION_LINE.search(manifest.read_text(encoding="utf-8"))
    if match is None:
        raise ToolError(f"chloe: {manifest} declares no version")
    return match.group(1)


def reference_set(script: Path) -> str:
    """The tree hash of the reference set this checkout pins — the data an annotation rests on.

    Read from ``Artifacts.toml`` beside the script rather than assumed. Absent or unparsable is an
    error and not a shrug: the alternative is a number whose references nobody can name.
    """
    manifest = script.parent / "Artifacts.toml"
    if not manifest.is_file():
        raise ToolError(
            f"chloe: {manifest} is missing, so the reference set behind an annotation cannot be "
            f"identified — check that --set chloe= points at a Chloe clone"
        )
    match = re.search(
        rf'\[{_ARTIFACT_SECTION}\]\s*\ngit-tree-sha1\s*=\s*"([0-9a-f]+)"',
        manifest.read_text(encoding="utf-8"),
    )
    if match is None:
        raise ToolError(
            f"chloe: {manifest} does not pin a [{_ARTIFACT_SECTION}] git tree hash — the "
            f"reference set cannot be identified, and {CHLOE_VERSION} is not a version of *them*"
        )
    return match.group(1)


def run_chloe(
    genome: Genome,
    ctx: RunContext,
    *,
    seq_id: str,
    script: str | Path | None = None,
    julia: str | Path | None = None,
    transl_table: int | None = None,
    extra_args: tuple[str, ...] = (),
) -> ChloeResult:
    """Annotate one genome and read the result back with **our own reader**.

    The read-back is not decoration: a file that does not re-read is not a result, and an
    annotator's GFF is a third-party file like any other — it gets no privileges.
    """
    chloe = chloe_script(script)
    interpreter = julia_binary(julia)
    release = package_version(chloe)
    if release != CHLOE_VERSION:
        raise ToolError(
            f"chloe: this checkout is Chloe {release}, but the toolkit's numbers were measured "
            f"with {CHLOE_VERSION} (the release published at doi:10.5281/zenodo.18194040) — "
            f"re-measure on the golden genome and update CHLOE_VERSION, or check out the pinned "
            f"release."
        )
    tree = reference_set(chloe)
    if tree != CHLOE_REFS_TREE:
        raise ToolError(
            f"chloe: this checkout's reference set is {tree}, but the toolkit's numbers were "
            f"measured with {CHLOE_REFS_TREE} (Chloe_refs 1.0.1) — a different reference set is a "
            f"different annotation. Re-measure on the golden genome and update CHLOE_REFS_TREE, "
            f"or check out the pinned reference set."
        )
    julia_version, _ = tool_version("julia", "binary", interpreter)
    stem = Path(seq_id).name

    with tempfile.TemporaryDirectory(prefix="chloe_") as work:
        workdir = Path(work)
        fasta = workdir / f"{stem}.fa"
        fasta.write_text(to_fasta(genome, seq_id), encoding="utf-8")

        ctx.tool(
            TOOL_NAME,
            [
                interpreter,
                f"--project={chloe.parent}",
                str(chloe),
                "annotate",
                # keep the coordinates in the frame of the sequence we handed over (see module
                # docstring): without this Chloe flips and orients the genome and annotates that
                "--no-transform",
                *extra_args,
                str(fasta),
            ],
            kind="binary",
            tool_path=str(chloe),
            pinned_version=CHLOE_VERSION,
            # julia's own --version line is already self-describing ("julia version 1.13.0")
            interpreter=julia_version or "julia",
            cwd=workdir,
            timeout=3600,  # the first run of a session pays Julia's compilation time
        )

        _refuse_a_transformed_sequence(workdir, genome)
        produced = sorted(workdir.glob("*.chloe.gff"))
        if not produced:
            raise ToolError(
                f"chloe: the run wrote no .chloe.gff for {stem} — see the captured log (it names "
                f"what it skipped)"
            )
        if len(produced) > 1:
            raise ToolError(
                f"chloe: the run produced {len(produced)} GFF files for one target — one record "
                f"per run is the contract"
            )
        text = produced[0].read_text(encoding="utf-8")
        declared = _GFF_VERSION_LINE.search(text)
        if declared is None:
            raise ToolError(
                "chloe: the output declares no '##source-version Chloe …' line, so which GFF "
                "this reader is looking at cannot be named — the run is not reportable"
            )
        if declared.group(1) != CHLOE_GFF_VERSION:
            raise ToolError(
                f"chloe: this run wrote '##source-version Chloe {declared.group(1)}', but this "
                f"reader was written for {CHLOE_GFF_VERSION} — re-read the format, then update "
                f"CHLOE_GFF_VERSION rather than trusting a file nobody has looked at"
            )

        genes = genes_from_gff(text)
        if transl_table is not None:
            # a GFF carries no /transl_table of its own unless the writer put one there; the
            # caller declares the code the organism uses and it is recorded, never guessed
            genes = tuple(
                replace(gene, transl_table=transl_table)
                if gene.kind == "CDS" and gene.transl_table is None
                else gene
                for gene in genes
            )
        annotated = Genome(
            name=genome.name,
            sequence=genome.sequence,
            accession=genome.accession,
            topology=genome.topology,
            genes=genes,
            repeat_regions=genome.repeat_regions,
        )
        if not any(gene.kind in GENE_KINDS for gene in annotated.genes):
            # the same rule as PGA's: an annotator that reports success while finding nothing has
            # failed, and a report of zeros would hide it
            raise ToolError(
                "chloe: the run produced a GFF with no CDS, tRNA or rRNA feature at all, while "
                "reporting success — an empty annotation is a failed annotation"
            )
        return ChloeResult(
            genome=annotated,
            gff_text=text,
            source_name=produced[0].name,
            version=release,
            reference_tree=tree,
        )


def _refuse_a_transformed_sequence(workdir: Path, genome: Genome) -> None:
    """Refuse a run whose coordinates belong to a sequence the caller does not have.

    Chloe writes a FASTA beside the GFF **only when it transformed the sequence** — the measured
    runs without ``--no-transform`` wrote ``<stem>.chloe.fa``, the run with it wrote none. So the
    file's existence is the signal, and comparing the two sequences says whether the coordinates
    of this GFF are the caller's. A mismatch is refused by name: silently annotating a flipped
    genome would put a third of the genes in the wrong frame while the gene *count* stayed
    perfect, which is this project's definition of a silent-wrong.
    """
    written = sorted(workdir.glob("*.chloe.fa"))
    if not written:
        return
    sequence = "".join(
        line.strip()
        for line in written[0].read_text(encoding="utf-8").splitlines()
        if not line.startswith(">")
    ).upper()
    input_sequence = genome.sequence.upper()
    if sequence != input_sequence:
        # how many positions moved is the evidence a reader needs to see this was not a typo; a
        # length difference (a trimmed or padded output) is the same failure and is named as one
        differing = sum(
            1 for left, right in zip(sequence, input_sequence, strict=False) if left != right
        )
        lengths = (
            f" and is {len(sequence)} instead of {len(input_sequence)} bases long"
            if len(sequence) != len(input_sequence)
            else ""
        )
        raise ToolError(
            f"chloe: this run transformed the sequence — its own FASTA differs from the input at "
            f"{differing} position(s){lengths}, so the GFF's coordinates are not in the input's "
            f"frame. The run must keep --no-transform (see plastid/external/chloe.py)."
        )
