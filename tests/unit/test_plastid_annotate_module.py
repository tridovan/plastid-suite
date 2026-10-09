"""The plastid annotation module (P3B-08) — the refusals around the searches.

BLAST and ARAGORN are stubbed: the searches are ``external/``'s job and their own tests hold the
command lines still. What is tested here is what the module does when a search comes back with
nothing — the clause P3B-10's Done-when states last, *"a search that finds nothing is a failed
run"* — which until now was asserted **only on the mito path** while both plastid branches had no
test at all (measured 10/10/2026, while retiring the old ``annotation`` module, whose tests were
the only ones covering it).

The two branches are separate and each has its own message, because they are different findings: a
search with no hit at all says the input or the reference set is wrong, and a search whose hits all
fail the floors says the reference set is too far away — the fix differs, so the sentence does.
"""

from __future__ import annotations

import stat
from pathlib import Path
from random import Random

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from plastid.core.ir import IRResult
from plastid.core.model import Genome, Interval
from plastid.core.sequence import reverse_complement
from plastid.external.blast import BLAST_VERSION
from plastid.modules import plastid_annotate
from plastid.run.context import RunContext
from plastid.run.layout import RunLayout
from plastid.run.tools import ToolError

SEQ_ID = "SYN_PLASTID.1"
# LSC 0..2999 · IRa 3000..3999 · SSC 4000..4999 · IRb 5000..5999: the repeat's two copies are the
# same sequence, which is what the gate in front of this module asks for.
_IRA = "ACGT" * 250
_LSC = "ATGC" * 750
_SSC = "TTAA" * 250
GENOME = Genome(
    name="synthetic plastome",
    accession=SEQ_ID,
    sequence=_LSC + _IRA + _SSC + reverse_complement(_IRA),
    topology="circular",
)

# One CDS hit that covers a 9-residue reference protein, in BLAST's own blastx shape (the query is
# the genome, the subject the reference protein — see `external/blast.py`).
_PROTEIN = "M" + "A" * 8  # ATG GCT×8, then the stop the audit looks for
# The reference record's own sequence, in the same four-part shape as the target and for the same
# reason: the audit in front of the search refuses a record whose inverted repeat it cannot find,
# and the de-novo scan accepts only the validated 10–35 kb window — so the fixture's repeat is
# 10 kb, the smallest one that passes. **Not a low-complexity repeat**: a "ACGT"×2,500 copy puts
# every 20-mer in the index ten thousand times and the de-novo scan goes quadratic (measured:
# over two minutes, on a fixture) — a pseudo-random copy indexes like a real inverted repeat.
_RANDOM = Random(20_261_010)
_REF_IR = "".join(_RANDOM.choices("ACGT", k=10_000))
_REF_LSC = "ATG" + "GCT" * 8 + "TAA" + "".join(_RANDOM.choices("ACGT", k=2_970))
_REF_SSC = "".join(_RANDOM.choices("ACGT", k=3_000))
_REF_SEQ = _REF_LSC + _REF_IR + _REF_SSC + reverse_complement(_REF_IR)
BLASTX_HITS = (
    f"{SEQ_ID}\tNC_999999.1|psbA|psbA|Testus testus\t99.0\t9\t0\t0\t1\t27\t1\t9\t1e-50\t200\t1\n"
)
# One rRNA hit, same shape without a frame.
BLASTN_HITS = (
    f"{SEQ_ID}\tNC_999999.1|rrn16|rrn16|Testus testus"
    f"\t96.0\t116\t0\t0\t61\t176\t1\t116\t1e-60\t300\n"
)
# ARAGORN's measured batch shape (the real lines live in test_aragorn.py).
ARAGORN_HITS = (
    f">{SEQ_ID}\n1 genes found\n1   tRNA-Phe                     [281,350]\t35  \t(gaa)\n"
)


def _ir() -> IRResult:
    return IRResult(
        ira=Interval(3000, 4000),
        irb=Interval(5000, 6000),
        identity=1.0,
        mismatches=0,
        compared_length=1000,
        exact_core=1000,
        mirror_constant=3000 + 6000 - 1,
        source="annotation",
    )


@pytest.fixture(autouse=True)
def _no_recorded_manifest(monkeypatch: pytest.MonkeyPatch) -> None:
    """Isolate the module from the recorded reference set on this machine.

    ``plastid_annotate`` reads the committed manifest to apply leave-one-out, and a stub accession
    is not in it — so without this the run would refuse the fixture record as "excluded for the
    target" before reaching the searches this file is about. The manifest's own behaviour has its
    own tests (``test_plastid_reference.py``).
    """
    monkeypatch.setattr(plastid_annotate, "load_recorded_manifest", lambda: None)


def _ctx(tmp_path: Path) -> RunContext:
    layout = RunLayout.create(tmp_path / "runs", "plastid-annotate-test", fingerprint="test")
    layout.make_dirs()
    return RunContext(
        layout=layout,
        params={},
        module_id="plastid_annotate",
        module_order=1,
        ir=_ir(),
    )


def _stub(directory: Path, name: str, output: str) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    script = directory / name
    if name == "aragorn":
        script.write_text(
            "#!/bin/sh\n"
            'if [ "$1" = "-h" ]; then echo "ARAGORN v1.2.41 Dean Lassett"; exit 0; fi\n'
            f"cat <<'HITS'\n{output}HITS\n",
            encoding="utf-8",
        )
    else:
        script.write_text(
            "#!/bin/sh\n"
            'case "$1" in\n'
            f"  --version|-version) echo '{name}: {BLAST_VERSION}+'; exit 0;;\n"
            "esac\n"
            f"cat <<'HITS'\n{output}HITS\n",
            encoding="utf-8",
        )
    script.chmod(script.stat().st_mode | stat.S_IXUSR)


def _binaries(tmp_path: Path, *, blastx: str = BLASTX_HITS, blastn: str = BLASTN_HITS) -> Path:
    directory = tmp_path / "bin"
    _stub(directory, "makeblastdb", "")
    _stub(directory, "blastx", blastx)
    _stub(directory, "blastn", blastn)
    _stub(directory, "aragorn", ARAGORN_HITS)
    return directory


def _reference(tmp_path: Path) -> Path:
    """A one-record reference set that passes the audit gate in front of the search.

    A complete CDS with its own ``/translation`` and a complete rRNA — the two things the audit
    reads, and the two the searches are run against.
    """
    directory = tmp_path / "reference"
    directory.mkdir(parents=True, exist_ok=True)
    record = SeqRecord(Seq(_REF_SEQ), id="NC_999999.1", name="NC_999999.1")
    record.annotations = {"molecule_type": "DNA", "topology": "circular", "organism": "Testus"}
    record.features = [
        SeqFeature(
            SimpleLocation(0, len(_REF_SEQ), strand=1),
            type="source",
            qualifiers={"organism": ["Testus"], "organelle": ["plastid:chloroplast"]},
        ),
        SeqFeature(
            SimpleLocation(0, 30, strand=1),
            type="CDS",
            qualifiers={"gene": ["psbA"], "translation": [_PROTEIN], "codon_start": ["1"]},
        ),
        SeqFeature(SimpleLocation(100, 216, strand=1), type="rRNA", qualifiers={"gene": ["rrn16"]}),
    ]
    SeqIO.write(record, str(directory / "NC_999999.1.gb"), "genbank")
    return directory


def _run(tmp_path: Path, **extra: object):
    params: dict[str, object] = {
        "reference": str(_reference(tmp_path)),
        "tools": str(_binaries(tmp_path)),
        **extra,
    }
    return plastid_annotate.run(GENOME, params, _ctx(tmp_path))


def test_the_module_annotates_the_genome_it_was_handed(tmp_path: Path) -> None:
    """The baseline the two refusals below are measured against: this fixture does annotate."""
    result = _run(tmp_path)

    assert {row["gene"] for row in result.tables["called_genes"].rows} >= {"psbA", "rrn16"}
    assert result.artifacts["plastid_annotation.gb"]


def test_a_search_that_finds_nothing_is_a_failed_run_not_an_empty_annotation(
    tmp_path: Path,
) -> None:
    """Both searches empty: the input is not a plastome, or the reference set is not its family.

    The message names both, because those are the two things a user can act on — an annotation
    with no genes is not a result, and reporting one would be the silent failure D5-08 forbids.
    """
    params = {
        "reference": str(_reference(tmp_path)),
        "tools": str(_binaries(tmp_path, blastx="", blastn="")),
    }

    with pytest.raises(ToolError, match="no hit at all"):
        plastid_annotate.run(GENOME, params, _ctx(tmp_path))


def test_hits_that_all_fail_the_floors_are_a_failed_run_and_the_floors_are_named(
    tmp_path: Path,
) -> None:
    """A different finding with its own message: the search worked, the reference set is too far.

    The sentence has to point at the floors rather than at the input, because *that* is what a
    user would change — and the run has to have named every rejected locus in ``not_called``.
    """
    params = {
        "reference": str(_reference(tmp_path)),
        "tools": str(_binaries(tmp_path)),
        "min_identity": 99.9,
    }

    with pytest.raises(ToolError, match="no locus passed the floors"):
        plastid_annotate.run(GENOME, params, _ctx(tmp_path))
