"""M11.17 panel cases: the remaining modules, on laid-out records written from a new origin.

`test_m11_17_rotated_annotation_cases.py` covers `codon_usage`, `amino_acid`, `stats` and the
three whole-plastome alignment modules. This file covers the six modules it left out, and states
the claim each one can carry:

* `gene_table` and `orf` read genes and no alignment. Rotating a record renames positions and
  changes nothing else, so **every table row is identical**. `orf`'s `cds_id`, `start_1b` and
  `end_1b` cells are the record's coordinates by definition — they are the only cells excluded,
  and they are excluded because they *are* positions, not because they disagreed.
* `kaks` and `popgen` align each shared CDS. A rotation preserves the CDS string exactly, so the
  pair fed to MAFFT and yn00 is the same pair and **every table is identical**. Each input is
  rotated by its own offset, because the claim is not "the same shift cancels".
* `synteny` aligns the reference in its **normalized** frame, and normalisation is exactly
  origin-independent (probe: the normalized sequence and the four region lengths are identical
  for every offset). Rotating the reference therefore leaves **every table identical**. Rotating
  the *query* is a different claim and is **not** asserted for the block table: MUMmer linearises
  the query at its own origin, so a match across that cut splits into two blocks and the per-block
  max/mean identity moves. What must not move — how much of the reference is covered and how
  identical those bases are — is asserted for the query rotation too.
* `phylogeny` aligns each CDS and infers one tree with a fixed seed. The alignment input is the
  CDS set, which a rotation does not change, so the supermatrix and the tree are **byte-identical**
  for a rotated panel.

The independent expectation is stated per test: a rotation is a change of coordinates only, and a
number that moves is a bug unless the number is itself a coordinate. The tool-using tests skip
with a named reason when the tool is absent (`PTS_REQUIRE_TOOLS=1` turns that into a failure).
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from conftest import require_tools
from plastid.core.genbank import read_genbank
from plastid.core.ir import detect_ir, enforce_ir_gate
from plastid.core.normalize import rotate_genome
from plastid.modules import gene_table, kaks, orf, phylogeny, popgen, synteny
from plastid.run.context import RunContext
from plastid.run.layout import RunLayout

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "tests" / "data"
PANEL = ROOT / ".panel" / "ir-validation" / "source"

GOLDEN = "NC_050833.1"
QUERY = "NC_024936.1"
# Four panel plastomes: IQ-TREE refuses ultrafast bootstrap with fewer than four sequences, and
# `popgen` leaves Fst blank unless each group carries at least two. Two per group, two groups.
PANEL_RECORDS = ("NC_008641.1", "NC_016708.2", "NC_028027.1", "NC_016736.1")

pytestmark = pytest.mark.tier_b

# Cells `orf` prints that *are* the record's coordinates. Nothing else is excluded anywhere.
ORF_COORDINATE_CELLS = ("cds_id", "start_1b", "end_1b")


def _ctx(root: Path, module: str, genome=None) -> RunContext:
    layout = RunLayout.create(root, module, fingerprint="m11-17")
    layout.make_dirs()
    ctx = RunContext(layout=layout, params={}, module_id=module, module_order=1)
    if genome is not None:
        # Position-based modules read the reference's IR from the context; a rotated record
        # carries its own IR at different coordinates.
        ctx.ir = enforce_ir_gate(detect_ir(genome, source="sequence"))
    return ctx


def _rows(result, drop: tuple[str, ...] = ()) -> dict[str, list[tuple[tuple[str, str], ...]]]:
    """Tables as an order-free set of rows, values rendered with ``repr``.

    ``repr`` keeps a blank cell (``None`` for a measurement that was not made) comparable with a
    number without sorting ``None`` against ``int``.
    """
    return {
        table_id: sorted(
            tuple(sorted((key, repr(value)) for key, value in row.items() if key not in drop))
            for row in table.rows
        )
        for table_id, table in result.tables.items()
    }


def _golden_pair():
    return read_genbank(DATA / f"{GOLDEN}.gb"), read_genbank(DATA / f"{QUERY}.gb")


def _iqtree_dir() -> Path:
    """IQ-TREE's directory, from `PTS_IQTREE_DIR` or the shared toolbox layout."""
    directory = os.environ.get("PTS_IQTREE_DIR")
    candidates = [Path(directory)] if directory else []
    candidates.append(Path("/root/toolbox/env/iqtree/bin"))
    for candidate in candidates:
        if (candidate / "iqtree2").is_file() or (candidate / "iqtree").is_file():
            return candidate
    reason = "iqtree2 is not installed; set PTS_IQTREE_DIR to the directory holding it"
    if os.environ.get("PTS_REQUIRE_TOOLS") == "1":
        pytest.fail(reason)
    pytest.skip(reason)


def _panel_genomes():
    paths = [PANEL / f"{accession}.gb" for accession in PANEL_RECORDS]
    if not all(path.is_file() for path in paths):
        pytest.skip("the IR panel records are not in this tree (.panel/ is gitignored)")
    return [read_genbank(path) for path in paths]


# ------------------------------------------------------------------ pure, gene-based modules


def test_the_gene_table_does_not_depend_on_the_origin(tmp_path: Path) -> None:
    """Grouping and counting genes is not a position-based measurement: no row may move."""
    genome = read_genbank(DATA / f"{GOLDEN}.gb")
    base = gene_table.run(genome, {}, _ctx(tmp_path / "base", "gene_table"))
    assert base.tables["gene_table_summary"].rows[0]["n_genes_total"] > 0
    baseline = _rows(base)

    for shift in (1, len(genome) // 4, len(genome) // 2, len(genome) - 1):
        rotated = rotate_genome(genome, shift)
        tables = _rows(gene_table.run(rotated, {}, _ctx(tmp_path / f"s{shift}", "gene_table")))
        assert tables == baseline, f"rotation by {shift} moved a gene-table row"


def test_orf_checks_do_not_depend_on_the_origin(tmp_path: Path) -> None:
    """A CDS's status is a property of its sequence; only its coordinates may change.

    ``cds_id`` embeds the span (``psbT:73435-73542``), so it is a coordinate cell like the start
    and end columns. Everything else — codon counts, start/stop codon, internal stops, status —
    is the same string after rotation, including a CDS the rotation cuts at the new origin.
    """
    genome = read_genbank(DATA / f"{GOLDEN}.gb")
    base = orf.run(genome, {}, _ctx(tmp_path / "base", "orf"))
    assert base.tables["orf_summary"].rows[0]["n_cds"] > 0
    baseline = _rows(base, drop=ORF_COORDINATE_CELLS)

    for shift in (1, len(genome) // 4, len(genome) // 2, len(genome) - 1):
        rotated = rotate_genome(genome, shift)
        tables = _rows(
            orf.run(rotated, {}, _ctx(tmp_path / f"s{shift}", "orf")), drop=ORF_COORDINATE_CELLS
        )
        assert tables == baseline, f"rotation by {shift} changed an ORF status"


# ------------------------------------------------------------------ CDS-alignment modules


def test_kaks_does_not_depend_on_the_origin(
    tools: dict[str, str | None], paml_tools: dict[str, str | None], tmp_path: Path
) -> None:
    """Each gene is fed to MAFFT as its CDS string; a rotation leaves that string unchanged.

    Both records are rotated, by different offsets: the claim is that the *sequences* reach the
    tools unchanged, not that one particular shift cancels out.
    """
    require_tools(tools, "mafft")
    require_tools(paml_tools, "yn00")
    reference, query = _golden_pair()

    base = kaks.run([reference, query], {}, _ctx(tmp_path / "base", "kaks"))
    scored = [row for row in base.tables["kaks_genes"].rows if row["status"] == "ok"]
    assert scored, "no gene produced a Ka/Ks; an invariant empty table proves nothing"
    baseline = _rows(base)

    rotated_reference = rotate_genome(reference, len(reference) // 2)
    rotated_query = rotate_genome(query, len(query) // 3)
    rotated = kaks.run([rotated_reference, rotated_query], {}, _ctx(tmp_path / "rotated", "kaks"))
    assert _rows(rotated) == baseline


def test_popgen_does_not_depend_on_the_origin(tools: dict[str, str | None], tmp_path: Path) -> None:
    """π, Tajima's D and Fst are computed from the CDS alignment; the alignment is unchanged.

    Four genomes in two groups of two: with one sequence per group every population statistic is
    blank by design (`n = 1`), and an invariant table of blanks would prove nothing. Each genome
    is rotated by its own offset.
    """
    require_tools(tools, "mafft")
    genomes = _panel_genomes()
    groups = {genome.accession: ("A" if index < 2 else "B") for index, genome in enumerate(genomes)}

    base = popgen.run(genomes, {"groups": groups}, _ctx(tmp_path / "base", "popgen"))
    fst_values = [row["fst"] for row in base.tables["popgen_fst"].rows if row["fst"] is not None]
    assert fst_values, "no Fst was computed; an invariant empty table proves nothing"
    baseline = _rows(base)

    rotated = [
        rotate_genome(genome, len(genome) // divisor)
        for genome, divisor in zip(genomes, (2, 3, 5, 7), strict=True)
    ]
    rotated_groups = {
        genome.accession: ("A" if index < 2 else "B") for index, genome in enumerate(rotated)
    }
    moved = popgen.run(rotated, {"groups": rotated_groups}, _ctx(tmp_path / "rotated", "popgen"))
    assert _rows(moved) == baseline


# ------------------------------------------------------------------ position-based module


def test_synteny_follows_the_reference_annotation_not_its_origin(
    mummer_tools: dict[str, str | None], tmp_path: Path
) -> None:
    """The reference is normalized before MUMmer sees it, so rotating it moves no number.

    Measured 29/10/2026: `normalize` returns the same sequence and the same four region lengths
    for every offset, so the delta files are the same files and every table row is identical.

    A rotated **query** is a different case and is not compared block for block. MUMmer linearises
    the query at the record's origin; a match across that cut becomes two blocks, which changes
    `n_blocks`, `min/max/mean` per-block identity, and the block coordinates. Those are the
    query's own coordinates, not the reference's measurement. What the module measures — aligned
    base pairs, coverage, identity over base pairs, and the identity series along the reference —
    must not move, and is asserted.
    """
    require_tools(mummer_tools, "nucmer", "delta-filter", "show-coords")
    reference, query = _golden_pair()

    base = synteny.run([reference, query], {}, _ctx(tmp_path / "base", "synteny", reference))
    measured = base.tables["synteny_genomes"].rows[0]
    assert measured["aligned_bp"] > 0 and measured["n_blocks"] > 0
    baseline = _rows(base)

    for shift in (len(reference) // 2, len(reference) // 3):
        rotated = rotate_genome(reference, shift)
        tables = _rows(
            synteny.run([rotated, query], {}, _ctx(tmp_path / f"s{shift}", "synteny", rotated))
        )
        assert tables == baseline, f"rotating the reference by {shift} moved a synteny row"

    invariants = ("aligned_bp", "coverage_pct", "aligned_identity")
    for shift in (len(query) // 2, len(query) // 3):
        rotated = rotate_genome(query, shift)
        result = synteny.run(
            [reference, rotated], {}, _ctx(tmp_path / f"q{shift}", "synteny", reference)
        )
        row = result.tables["synteny_genomes"].rows[0]
        assert all(row[key] == measured[key] for key in invariants), (
            f"rotating the query by {shift} moved the reference-side measurement"
        )
        assert _rows(result)["synteny_identity"] == baseline["synteny_identity"]


# ------------------------------------------------------------------ tree


def test_the_phylogeny_tree_does_not_depend_on_the_origin(
    tools: dict[str, str | None], tmp_path: Path
) -> None:
    """The supermatrix is the CDS set; rotation leaves it, and the seeded tree, byte-identical.

    Each genome is rotated by its own offset. The tables are compared with ``tool_call`` excluded
    because that cell carries the throwaway matrix path IQ-TREE was given, which is a path, not a
    measurement; the matrix and the tree themselves are compared as text.
    """
    require_tools(tools, "mafft")
    iqtree = _iqtree_dir()
    genomes = _panel_genomes()
    params = {"iqtree": str(iqtree)}

    base = phylogeny.run(genomes, params, _ctx(tmp_path / "base", "phylogeny"))
    assert base.tables["phylogeny_summary"].rows[0]["alignment_sites"] > 0
    baseline_tables = _rows(base, drop=("tool_call",))

    rotated = [
        rotate_genome(genome, len(genome) // divisor)
        for genome, divisor in zip(genomes, (2, 3, 5, 7), strict=True)
    ]
    moved = phylogeny.run(rotated, params, _ctx(tmp_path / "rotated", "phylogeny"))

    assert _rows(moved, drop=("tool_call",)) == baseline_tables
    assert moved.artifacts["supermatrix.fasta"] == base.artifacts["supermatrix.fasta"]
    assert moved.artifacts["tree.treefile"] == base.artifacts["tree.treefile"]
