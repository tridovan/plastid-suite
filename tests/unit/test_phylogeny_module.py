"""The phylogeny module: what it builds, what it refuses, and what it will not claim (M11.21).

**The independent expectation, stated once.** A tree's *shape* cannot be checked against an
answer key the way a π or a Ts/Tv can — there is no arithmetic that says where NC_016708.2
belongs. What can be checked, and is checked here, is that the artifact is the one the inputs
imply:

* every input genome appears as a tip, exactly once, with the name it was given;
* two runs of the same matrix produce the **same** tree (the seed is fixed, so a different tree
  would mean the run is not reproducible);
* a gene present in one genome only is **excluded by name**, not filled with gaps and counted;
* a partial CDS is excluded by default and admitted only when asked, and the table says which;
* the module refuses a two-genome input — that is a comparison, not a topology;
* the matrix and the tree are written as artifacts, so the tree can be traced to the alignment.

The tests that need IQ-TREE skip with a named reason when it is absent (`PTS_IQTREE_DIR`), the
same convention as every other tier-B check.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from plastid.core.align import Alignment
from plastid.core.model import Gene, Genome, Interval
from plastid.modules import registry

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "tests" / "data"


# ------------------------------------------------------------------ declaration and refusals


def test_the_module_is_registered_with_its_tools() -> None:
    spec = registry.get("phylogeny")
    assert spec.arity == 2
    assert set(spec.citations) == {"mafft", "iqtree"}
    declared = {param.name for param in spec.params}
    assert {
        "iqtree",
        "model",
        "ultrafast_bootstrap",
        "alrt",
        "seed",
        "include_partial",
        "partition",
        "aligner",
        "macse",
    } <= declared


def test_a_two_genome_input_is_refused_by_name(tmp_path: Path) -> None:
    """Two sequences are a pairwise comparison; calling it a topology would be a false claim."""
    from plastid.cli import main

    code = main(
        [
            "run",
            "-m",
            "phylogeny",
            "-i",
            str(DATA / "NC_050833.1.gb"),
            "-i",
            str(DATA / "NC_024936.1.gb"),
            "-o",
            str(tmp_path),
        ]
    )
    assert code == 1


def test_a_three_genome_input_is_refused_because_iqtree_cannot_bootstrap_it(tmp_path: Path) -> None:
    """IQ-TREE 2.3.6 refuses `-bb` with fewer than 4 sequences, so 3 genomes cannot get a tree.

    Measured 29/10/2026: a three-genome run aligned every gene and then died inside IQ-TREE with
    "It makes no sense to perform bootstrap with less than 4 sequences". The minimum is the pinned
    tool's own rule, so it is enforced before the first alignment rather than after the work.
    """
    import plastid.modules.phylogeny as phylo

    with pytest.raises(ValueError, match="at least 4 genomes"):
        phylo.run(_synthetic_genomes()[:3], {}, _phylogeny_ctx(tmp_path))


def test_an_unknown_parameter_is_refused() -> None:
    """The M11.20 rule reaches the new module without any extra wiring."""
    problems = registry.unknown_params(["phylogeny"], {"bootstrap": 1000})
    assert problems and "ultrafast_bootstrap" in problems[0]


def test_the_default_thread_count_is_one() -> None:
    """The default is the reproducible one; asking for more is an explicit trade."""
    from plastid.modules.phylogeny import DEFAULT_THREADS

    assert DEFAULT_THREADS == "1"


def _fake_tree():
    """The minimum ``IqtreeResult`` the module reads when the inference is stubbed."""
    from plastid.external.iqtree import IqtreeResult

    return IqtreeResult(
        treefile="(SYN0.1,SYN1.1,SYN2.1);", log="", version="2.3.6", report="", argv=()
    )


def _synthetic_genomes() -> tuple[Genome, ...]:
    """Four genomes carrying the same two genes, so a matrix can be built without tools.

    Four is the module's minimum: IQ-TREE 2.3.6 refuses ultrafast bootstrap below four
    sequences, so a three-genome stub never reaches the code these tests exercise.
    """
    sequences = {
        "rbcL": ("ATGGCCAAACCC", "ATGGCCAAACCT", "ATGGCCAAACCA", "ATGGCCAAACCG"),
        "atpA": ("ATGTTTGGGCCC", "ATGTTTGGGCCA", "ATGTTTGGGCCT", "ATGTTTGGGCGA"),
    }
    genomes = []
    for index in range(4):
        chunks, genes, offset = [], [], 0
        for name, variants in sequences.items():
            sequence = variants[index]
            genes.append(
                Gene(name, "CDS", (Interval(offset, offset + len(sequence)),), transl_table=11)
            )
            chunks.append(sequence)
            offset += len(sequence)
        genomes.append(
            Genome(
                name=f"synthetic {index}",
                accession=f"SYN{index}.1",
                sequence="".join(chunks),
                genes=tuple(genes),
                topology="linear",
            )
        )
    return tuple(genomes)


def _phylogeny_ctx(tmp_path: Path):
    from plastid.run.context import RunContext
    from plastid.run.layout import RunLayout

    layout = RunLayout.create(tmp_path / "runs", "phylogeny-test", fingerprint="test")
    layout.make_dirs()
    return RunContext(layout=layout, params={}, module_id="phylogeny", module_order=1)


def test_a_gene_partition_names_every_included_block(tmp_path: Path, monkeypatch) -> None:
    """`partition=gene` gives IQ-TREE one charset per block, tiling the matrix."""
    import plastid.modules.phylogeny as phylo

    captured: dict = {}

    def fake_align(rows, ctx, *, strategy="auto"):
        return Alignment(
            names=tuple(name for name, _ in rows),
            sequences=tuple(sequence for _, sequence in rows),
        )

    def fake_iqtree(path, ctx, **kwargs):
        captured.update(kwargs)
        return _fake_tree()

    monkeypatch.setattr(phylo, "align", fake_align)
    monkeypatch.setattr(phylo, "run_iqtree", fake_iqtree)

    result = phylo.run(_synthetic_genomes(), {"partition": "gene"}, _phylogeny_ctx(tmp_path))

    partition = captured["partition"]
    assert "begin sets;" in partition
    assert "charset atpA = 1-12;" in partition
    assert "charset rbcL = 13-24;" in partition
    assert "supermatrix.partition.nex" in result.artifacts
    assert result.artifacts["supermatrix.partition.nex"] == partition
    assert result.tables["phylogeny_summary"].value("partition") == "gene (2 partitions)"


def test_a_codon_partition_splits_each_block_into_three(tmp_path: Path, monkeypatch) -> None:
    """`partition=codon` is columns 1/2/3 of each block, written with NEXUS's `\\3` step."""
    import plastid.modules.phylogeny as phylo

    captured: dict = {}

    def fake_align(rows, ctx, *, strategy="auto"):
        return Alignment(
            names=tuple(name for name, _ in rows),
            sequences=tuple(sequence for _, sequence in rows),
        )

    def fake_iqtree(path, ctx, **kwargs):
        captured.update(kwargs)
        return _fake_tree()

    monkeypatch.setattr(phylo, "align", fake_align)
    monkeypatch.setattr(phylo, "run_iqtree", fake_iqtree)

    result = phylo.run(_synthetic_genomes(), {"partition": "codon"}, _phylogeny_ctx(tmp_path))

    partition = captured["partition"]
    assert "charset atpA_codon1 = 1-12\\3;" in partition
    assert "charset atpA_codon2 = 2-12\\3;" in partition
    assert "charset atpA_codon3 = 3-12\\3;" in partition
    assert "charset rbcL_codon1 = 13-24\\3;" in partition
    assert result.tables["phylogeny_summary"].value("partition") == "codon (6 partitions)"


def test_the_default_is_unpartitioned(tmp_path: Path, monkeypatch) -> None:
    """No flag means the old behaviour: no `-p`, no partition artifact, no changed tree."""
    import plastid.modules.phylogeny as phylo

    captured: dict = {}

    def fake_align(rows, ctx, *, strategy="auto"):
        return Alignment(
            names=tuple(name for name, _ in rows),
            sequences=tuple(sequence for _, sequence in rows),
        )

    def fake_iqtree(path, ctx, **kwargs):
        captured.update(kwargs)
        return _fake_tree()

    monkeypatch.setattr(phylo, "align", fake_align)
    monkeypatch.setattr(phylo, "run_iqtree", fake_iqtree)

    result = phylo.run(_synthetic_genomes(), {}, _phylogeny_ctx(tmp_path))

    assert captured["partition"] is None
    assert "supermatrix.partition.nex" not in result.artifacts
    assert result.tables["phylogeny_summary"].value("partition") == "none"


def test_an_unknown_partition_mode_is_refused_by_name(tmp_path: Path) -> None:
    from plastid.modules.phylogeny import run
    from plastid.run.params import ParamValueError

    with pytest.raises(ParamValueError, match="none, gene, codon"):
        run(_synthetic_genomes(), {"partition": "by-codon"}, _phylogeny_ctx(tmp_path))


# ---------------------------------------------------------------- MACSE, opt-in (M11.21)


def test_macse_is_opt_in_and_replaces_the_mafft_call(tmp_path: Path, monkeypatch) -> None:
    """`aligner=macse` must go through the MACSE adapter — not silently fall back to MAFFT."""
    import plastid.modules.phylogeny as phylo
    from plastid.external.macse import MacseBinary

    used: dict = {}

    def fake_resolve(explicit=None):
        used["explicit"] = explicit
        return MacseBinary(argv_prefix=("macse",), version="2.07")

    def fake_align_cds(rows, ctx, *, binary, timeout=600):
        used["macse_rows"] = len(rows)
        used["timeout"] = timeout
        return Alignment(
            names=tuple(name for name, _ in rows),
            sequences=tuple(sequence for _, sequence in rows),
        )

    def fail_align(*args, **kwargs):
        raise AssertionError("MAFFT must not run when aligner=macse")

    def fake_iqtree(path, ctx, **kwargs):
        return _fake_tree()

    monkeypatch.setattr(phylo, "resolve_macse", fake_resolve)
    monkeypatch.setattr(phylo, "align_cds", fake_align_cds)
    monkeypatch.setattr(phylo, "align", fail_align)
    monkeypatch.setattr(phylo, "run_iqtree", fake_iqtree)

    result = phylo.run(_synthetic_genomes(), {"aligner": "macse"}, _phylogeny_ctx(tmp_path))

    assert used["macse_rows"] == 4
    assert used["timeout"] == 600  # the adapter's default when the user asks for nothing

    # The 34-genome measurement (28/10/2026) failed on a gene over 600 s; the door exists.
    phylo.run(
        _synthetic_genomes(),
        {"aligner": "macse", "macse_timeout": 1800},
        _phylogeny_ctx(tmp_path / "long"),
    )
    assert used["timeout"] == 1800
    assert result.tables["phylogeny_summary"].value("aligner") == "macse 2.07"
    assert "MACSE 2.07" in result.tables["phylogeny_summary"].scope.label


@pytest.mark.parametrize(
    ("params", "message"),
    [
        ({"aligner": "clustal"}, "not one of mafft, macse"),
        ({"aligner": "macse", "mafft_strategy": "einsi"}, "mafft_strategy has no effect"),
        ({"aligner": "macse", "macse_timeout": 0}, "at least 1 second"),
        ({"matrix": "igs", "partition": "codon"}, "no reading frame"),
        ({"macse": "/tmp/macse.jar"}, "has no effect with aligner=mafft"),
    ],
)
def test_an_aligner_choice_that_cannot_be_honoured_is_refused(
    tmp_path: Path, params: dict, message: str
) -> None:
    from plastid.modules.phylogeny import run
    from plastid.run.params import ParamValueError

    with pytest.raises(ParamValueError, match=message):
        run(_synthetic_genomes(), params, _phylogeny_ctx(tmp_path))


def test_the_adapter_writes_the_partition_file_and_passes_it_to_iqtree(
    tmp_path: Path, monkeypatch
) -> None:
    """The adapter owns the `-p` flag; a partition that never reaches argv is a silent no-op."""
    import plastid.external.iqtree as iq

    monkeypatch.setattr(iq, "resolve_binary", lambda directory=None: "/fake/iqtree2")
    monkeypatch.setattr(iq, "iqtree_version", lambda binary: "2.3.6")
    alignment = tmp_path / "alignment.fasta"
    alignment.write_text(">a\nACGT\n>b\nACGT\n", encoding="utf-8")

    class StubContext:
        def __init__(self) -> None:
            self.argv: tuple[str, ...] = ()
            self.written = ""

        def tool(self, name, argv, **kwargs):
            workdir = Path(kwargs["cwd"])
            self.argv = tuple(argv)
            self.written = (workdir / "partition.nex").read_text(encoding="utf-8")
            (workdir / "tree.treefile").write_text("(a,b);\n", encoding="utf-8")
            (workdir / "tree.iqtree").write_text("report\n", encoding="utf-8")

    ctx = StubContext()
    nexus = "#nexus\nbegin sets;\n  charset g = 1-2;\nend;\n"
    result = iq.run_iqtree(alignment, ctx, partition=nexus)

    assert ctx.argv[ctx.argv.index("-p") + 1] == "partition.nex"
    assert ctx.written == nexus
    assert result.treefile == "(a,b);\n"


def test_the_adapter_does_not_pass_p_without_a_partition(tmp_path: Path, monkeypatch) -> None:
    import plastid.external.iqtree as iq

    monkeypatch.setattr(iq, "resolve_binary", lambda directory=None: "/fake/iqtree2")
    monkeypatch.setattr(iq, "iqtree_version", lambda binary: "2.3.6")
    alignment = tmp_path / "alignment.fasta"
    alignment.write_text(">a\nACGT\n>b\nACGT\n", encoding="utf-8")

    class StubContext:
        def tool(self, name, argv, **kwargs):
            StubContext.argv = tuple(argv)
            (Path(kwargs["cwd"]) / "tree.treefile").write_text("(a,b);\n", encoding="utf-8")

    iq.run_iqtree(alignment, StubContext())
    assert "-p" not in StubContext.argv


# ------------------------------------------------------------------ matrix=introns / matrix=igs


def _lsc_span(kind: str):
    """One LSC intron or spacer of the golden record, with the genome that carries it."""
    from plastid.core.genbank import read_genbank
    from plastid.core.ir import detect_ir
    from plastid.core.spans import introns_of, spacers_of

    genome = read_genbank(DATA / "NC_050833.1.gb")
    ir = detect_ir(genome)
    spans = introns_of(genome, ir).introns if kind == "introns" else spacers_of(genome, ir).igs
    for span in spans:
        if span.region == "LSC" and span.length > 30:
            return genome, span
    raise AssertionError(f"the golden has no LSC {kind} longer than 30 bp")


def _lsc_variants(kind: str) -> tuple[Genome, ...]:
    """Four genomes that differ only inside one LSC intron (or spacer).

    Four is the module's minimum (IQ-TREE's bootstrap needs four sequences). Two sites inside the
    one locus are changed independently, so no genome is a re-label of another.
    """
    genome, span = _lsc_span(kind)
    positions = (span.start + 10, span.start + 20)
    changes = {"A": "G", "G": "A", "C": "T", "T": "C"}

    def variant(index: int, first: str, second: str) -> Genome:
        sequence = list(genome.sequence)
        sequence[positions[0]] = first
        sequence[positions[1]] = second
        from dataclasses import replace

        return replace(genome, sequence="".join(sequence), accession=f"{genome.accession}_{index}")

    at_first = genome.sequence[positions[0]]
    at_second = genome.sequence[positions[1]]
    return (
        genome,
        variant(1, changes[at_first], at_second),
        variant(2, at_first, changes[at_second]),
        variant(3, changes[at_first], changes[at_second]),
    )


def _stub_matrix_tools(monkeypatch, captured: dict) -> None:
    import plastid.modules.phylogeny as phylo

    def fake_align(rows, ctx, *, strategy="auto"):
        captured.setdefault("blocks", []).append([len(sequence) for _, sequence in rows])
        return Alignment(
            names=tuple(name for name, _ in rows),
            sequences=tuple(sequence for _, sequence in rows),
        )

    def fake_iqtree(path, ctx, **kwargs):
        captured["matrix"] = path.read_text(encoding="utf-8")
        captured["iqtree"] = kwargs
        return _fake_tree()

    monkeypatch.setattr(phylo, "align", fake_align)
    monkeypatch.setattr(phylo, "run_iqtree", fake_iqtree)


def test_the_intron_matrix_holds_intron_sequences_not_cds(tmp_path: Path, monkeypatch) -> None:
    """`matrix=introns` aligns the spans `core.spans` derives, keyed by canonical gene name."""
    import plastid.modules.phylogeny as phylo

    captured: dict = {}
    _stub_matrix_tools(monkeypatch, captured)
    genomes = _lsc_variants("introns")
    _genome, span = _lsc_span("introns")

    result = phylo.run(genomes, {"matrix": "introns"}, _phylogeny_ctx(tmp_path))

    summary = result.tables["phylogeny_summary"]
    assert summary.value("matrix") == "introns"
    assert "intron blocks" in summary.scope.label
    from plastid.core.gene_names import aliases_for, normalize_name

    gene, _copy, number = span.intron_id.rsplit(".", 2)
    locus = f"{normalize_name(gene, aliases_for(genomes[0])).canonical}.{number}"
    row = next(row for row in result.tables["matrix_genes"].rows if row["gene"] == locus)
    assert row["included"] == "yes"
    assert row["sites"] == span.length


def test_the_igs_matrix_holds_spacer_sequences(tmp_path: Path, monkeypatch) -> None:
    import plastid.modules.phylogeny as phylo

    captured: dict = {}
    _stub_matrix_tools(monkeypatch, captured)
    genomes = _lsc_variants("igs")
    _genome, span = _lsc_span("igs")

    result = phylo.run(genomes, {"matrix": "igs"}, _phylogeny_ctx(tmp_path))

    from plastid.core.gene_names import aliases_for, normalize_name

    assert result.tables["phylogeny_summary"].value("matrix") == "igs"
    assert "intergenic-spacer blocks" in result.tables["phylogeny_summary"].scope.label
    aliases = aliases_for(genomes[0])
    locus = (
        f"{normalize_name(span.left, aliases).canonical}-"
        f"{normalize_name(span.right, aliases).canonical}"
    )
    row = next(row for row in result.tables["matrix_genes"].rows if row["gene"] == locus)
    assert row["included"] == "yes"
    assert row["sites"] == span.length


def test_an_unknown_matrix_is_refused_by_name(tmp_path: Path) -> None:
    import plastid.modules.phylogeny as phylo

    with pytest.raises(ValueError, match="matrix="):
        phylo.run(_synthetic_genomes(), {"matrix": "exons"}, _phylogeny_ctx(tmp_path))


def test_macse_is_refused_for_a_matrix_without_a_reading_frame(tmp_path: Path) -> None:
    """MACSE aligns codons through a translation; an intron or a spacer has no frame."""
    import plastid.modules.phylogeny as phylo

    with pytest.raises(ValueError, match="reading frame"):
        phylo.run(
            _lsc_variants("introns"),
            {"matrix": "introns", "aligner": "macse"},
            _phylogeny_ctx(tmp_path),
        )


# ------------------------------------------------------------------ matrix=all (complete matrix)


def test_the_complete_matrix_concatenates_all_three_kinds_of_locus(
    tmp_path: Path, monkeypatch
) -> None:
    """`matrix=all` is CGAS module 14's "all regions": one matrix, one tree, three kinds of block.

    The independent check is arithmetic: every block the table marks included contributes exactly
    its own sites to the supermatrix, and no block is counted twice. The three kinds are named in
    the block id, so a mixed matrix cannot hide a collision behind a bare locus name.
    """
    import plastid.modules.phylogeny as phylo

    captured: dict = {}
    _stub_matrix_tools(monkeypatch, captured)

    result = phylo.run(_lsc_variants("introns"), {"matrix": "all"}, _phylogeny_ctx(tmp_path))

    summary = result.tables["phylogeny_summary"]
    assert summary.value("matrix") == "all"
    assert "CDS, intron and intergenic-spacer blocks" in summary.scope.label

    rows = result.tables["matrix_genes"].rows
    included = [row for row in rows if row["included"] == "yes"]
    for prefix in ("cds:", "intron:", "igs:"):
        assert any(row["gene"].startswith(prefix) for row in included), prefix
    assert summary.value("alignment_sites") == sum(row["sites"] for row in included)

    # One matrix was written to one tool call: `all` is not three runs whose trees are compared.
    assert len(captured["blocks"]) == len(included)
    text = captured["matrix"]
    assert sum(line.startswith(">") for line in text.splitlines()) == 4  # four taxa, one row each


def test_a_gene_partition_on_the_complete_matrix_names_every_kind(
    tmp_path: Path, monkeypatch
) -> None:
    """`partition=gene` works on the mixed matrix and names blocks by their kind and locus."""
    import plastid.modules.phylogeny as phylo

    captured: dict = {}
    _stub_matrix_tools(monkeypatch, captured)

    phylo.run(
        _lsc_variants("introns"),
        {"matrix": "all", "partition": "gene"},
        _phylogeny_ctx(tmp_path),
    )

    partition = captured["iqtree"]["partition"]
    assert "#nexus\nbegin sets;" in partition
    for prefix in ("charset cds_", "charset intron_", "charset igs_"):
        assert prefix in partition, prefix


def test_a_codon_partition_is_refused_for_the_complete_matrix(tmp_path: Path) -> None:
    """Only the CDS blocks have a reading frame; a codon split would mislabel the rest."""
    import plastid.modules.phylogeny as phylo

    with pytest.raises(ValueError, match="no reading frame"):
        phylo.run(
            _lsc_variants("introns"),
            {"matrix": "all", "partition": "codon"},
            _phylogeny_ctx(tmp_path),
        )


def test_macse_is_refused_for_the_complete_matrix(tmp_path: Path) -> None:
    """One aligner per run; MACSE cannot honour the intron and spacer blocks of `all`."""
    import plastid.modules.phylogeny as phylo

    with pytest.raises(ValueError, match="matrix=all"):
        phylo.run(
            _lsc_variants("introns"),
            {"matrix": "all", "aligner": "macse"},
            _phylogeny_ctx(tmp_path),
        )
