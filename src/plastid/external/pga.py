"""PGA adapter — plastome annotation (M3.1, P3-03; decision D3-01).

PGA (Qu et al. 2019) annotates a plastome by homology: a directory of GenBank references and a
directory of FASTA targets in, GenBank out. **The annotator is the pipeline; GeSeq is a one-off
cross-check whose output gets pinned** — a web service cannot record the tool name, version and
argv every run has to carry (D3-01, `docs/D3-01_annotator.md`).

Three of PGA's habits are handled here rather than discovered later (TASKS.md pitfalls, measured
by the previous session):
it **clears its output directory**, it **breaks on filenames with spaces**, and it
**mis-annotates a multi-FASTA target**. So every run happens in a fresh temporary directory, the
target file is named from a slugified sequence id, and exactly one record is passed per call.

**Provisional.** PGA is not installed on the machine this was written on, so the argv below comes
from PGA's own README (`perl PGA.pl -r <refdir> -t <targetdir> [-i -p -q -o -f -l]`) and has not
been verified against a real run. The golden check that uses it is `PTS_REQUIRE_TOOLS=1`-gated and
claims nothing until a machine runs it.
"""

from __future__ import annotations

import tempfile
from dataclasses import dataclass
from pathlib import Path

from plastid.core.genbank import read_genbank, read_genbank_records
from plastid.core.model import Genome
from plastid.run.context import RunContext
from plastid.run.tools import ToolError, which

__all__ = ["PgaResult", "pga_script", "run_pga"]

TOOL_NAME = "pga"


@dataclass(frozen=True, slots=True)
class PgaResult:
    """The annotated genome PGA produced, and the GenBank text it was read from.

    The text travels back rather than a path: the run happens in a temporary directory that is
    gone by the time a module sees this, and the caller decides what to keep (a module hands it
    to :class:`~plastid.result.result.Result` as an artifact — seam 2).

    ``attempts`` is 2 when the first attempt produced nothing and the retry succeeded — a
    measured PGA bug (see :func:`run_pga`), recorded so a run can say it happened.
    """

    genome: Genome
    genbank_text: str
    source_name: str
    attempts: int = 1


def pga_script(explicit: str | Path | None = None) -> Path:
    """Locate ``PGA.pl`` — an explicit path, then ``PATH``. A missing one is a named error."""
    if explicit is not None:
        path = Path(explicit)
        if not path.is_file():
            raise ToolError(f"PGA.pl not found at {path} — pass --set pga=/path/to/PGA.pl")
        return path
    found = which("PGA.pl")
    if found is None:
        raise ToolError(
            "PGA.pl not found — PGA is a Perl script, not a package: clone "
            "https://github.com/quxiaojian/PGA (GPL-3.0) and pass --set pga=/path/to/PGA.pl "
            "(needs perl and BLAST+ 2.8.1 or newer)"
        )
    return Path(found)


def run_pga(
    genome: Genome,
    ctx: RunContext,
    *,
    reference_dir: str | Path,
    seq_id: str,
    script: str | Path | None = None,
    extra_args: tuple[str, ...] = (),
) -> PgaResult:
    """Annotate one genome and read the result back with **our own reader**.

    The round trip is not decoration: a file that does not re-read is not a result, and PGA's
    output is a third-party GenBank file like any other — it gets no privileges.
    """
    perl = which("perl")
    if perl is None:
        raise ToolError("pga: perl not found — PGA is a Perl script")
    pga = pga_script(script)
    reference = Path(reference_dir)
    if not reference.is_dir():
        raise ToolError(f"pga: the reference directory {reference} does not exist")
    if not any(reference.glob("*.gb*")):
        raise ToolError(f"pga: the reference directory {reference} holds no GenBank file")
    partial = _references_with_partial_coordinates(reference)
    if partial:
        # Measured 2026-10-03: a reference with partial coordinates makes PGA print warnings,
        # produce nothing, and still exit 0 — a silent failure this check turns into a named one.
        raise ToolError(
            f"pga: these reference files carry partial coordinates (< or >): "
            f"{', '.join(partial)} — PGA cannot parse them, and it fails *silently* (measured: "
            f"0 genes annotated, exit 0). Give it a reference set whose features are complete."
        )

    # PGA's temporary file names are random over an alphabet that includes ``$``, and it hands
    # them to blastn through ``system("blastn ... $name ...")`` — a second pass through the shell,
    # which expands ``$…`` away. Measured 2026-10-03: roughly one run in five then finds no file,
    # annotates nothing, and still exits 0. So one retry, and if both attempts come back empty the
    # failure is named rather than reported as an empty annotation.
    for attempt in (1, 2):
        annotated, text, name = _annotate_once(
            genome,
            ctx,
            perl=perl,
            pga=pga,
            reference=reference,
            seq_id=seq_id,
            extra_args=extra_args,
        )
        if any(gene.kind in ("CDS", "tRNA", "rRNA") for gene in annotated.genes):
            return PgaResult(
                genome=annotated, genbank_text=text, source_name=name, attempts=attempt
            )
    raise ToolError(
        "pga: two attempts both produced a record with no gene features, while PGA reported "
        "success. Check the captured log (it names what it skipped) and check that the references "
        "are from the same or a close family — a single empty attempt is the known temporary-name "
        "bug and is retried automatically, so two in a row is something else."
    )


def _annotate_once(
    genome: Genome,
    ctx: RunContext,
    *,
    perl: str,
    pga: Path,
    reference: Path,
    seq_id: str,
    extra_args: tuple[str, ...],
) -> tuple[Genome, str, str]:
    """One PGA run in its own fresh directory; the genome it wrote, its text, and its file name."""
    with tempfile.TemporaryDirectory(prefix="pga_") as work:
        workdir = Path(work)
        # PGA pairs its output with the FASTA header and breaks on spaces, so the id is passed
        # through a slug and one record travels per run
        target_dir = workdir / "target"
        out_dir = workdir / "gb"
        target_dir.mkdir()
        (target_dir / f"{seq_id}.fasta").write_text(_to_fasta(genome, seq_id), encoding="utf-8")

        ctx.tool(
            TOOL_NAME,
            [
                perl,
                str(pga),
                "-r",
                str(reference),
                "-t",
                str(target_dir),
                "-o",
                str(out_dir),
                *extra_args,
            ],
            kind="perl",
            tool_path=str(pga),
            pinned_version=None,
            interpreter="perl",
            cwd=workdir,
            timeout=3600,  # an annotation run over a plastome takes minutes, not seconds
        )

        produced = sorted(out_dir.glob("*.gb*")) if out_dir.is_dir() else []
        if len(produced) > 1:
            raise ToolError(
                f"pga: the run produced {len(produced)} GenBank files for one target — one "
                f"record per run is the contract, so a multi-record target slipped through"
            )
        if not produced:
            # a run that wrote nothing at all is the same failure as one that wrote no genes
            return Genome(name=genome.name, sequence=genome.sequence), "", ""
        text = produced[0].read_text(encoding="utf-8")
        return read_genbank(produced[0]), text, produced[0].name


def _references_with_partial_coordinates(reference: Path) -> tuple[str, ...]:
    """Reference files holding an incomplete coordinate — the ones PGA cannot read.

    Read with **our** reader rather than by scanning the text for ``<``/``>``: a GenBank record
    may mention either character in a note, and a false positive here would refuse a reference
    set that works.
    """
    names: list[str] = []
    for path in sorted(reference.glob("*.gb*")):
        try:
            records = read_genbank_records(path)
        except ValueError as error:
            raise ToolError(
                f"pga: the reference file {path.name} does not re-read: {error}"
            ) from error
        if any(
            part.partial_start or part.partial_end
            for record in records
            for gene in record.genes
            for part in gene.parts
        ):
            names.append(path.name)
    return tuple(names)


def _to_fasta(genome: Genome, seq_id: str, width: int = 60) -> str:
    body = "\n".join(
        genome.sequence[index : index + width] for index in range(0, len(genome), width)
    )
    return f">{seq_id}\n{body}\n"
