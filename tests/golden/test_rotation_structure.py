"""The same plastome written from five different origins.

`fasta/` holds six panel plastomes, each deposited sequence rotated to five starts (`0`, `L/4`,
`L/2`, `3L/4`, and the first base of IR copy B) with **no annotation carried over**. Rotation is
not data: it renames positions and changes nothing else, so every position-dependent result must
be identical under it. `EXPECTED.tsv` supplies the independent answer for each file — the two IR
copies' starts and lengths in the rotated frame, and whether a copy wraps position 1 — derived
from the deposited annotation by arithmetic, not by this detector.

`structure` is the module that can run on a bare sequence: it detects the IR itself and reports
the four regions in the normalized frame. Three claims are checked here:

* the four region rows (`length_bp`, `gc_pct`, `gc_skew`) are the same for all five rotations of
  one plastome — a rotation cannot move a base into or out of a region;
* the regions tile the chromosome exactly and the two copies have equal length — arithmetic that
  needs no panel answer;
* the de-novo IR length and the normalized origin reproduce `EXPECTED.tsv` — the answer was
  computed from the annotation, so agreement is evidence about the detector, not a restatement.

The annotated, gene-based modules cannot run on these files (they carry no genes). Their rotation
cases rotate the **annotated source GenBank** by the same offsets and live in
`test_rotation_gene_measurements.py`; the split is stated so nobody reads this file as
coverage of `codon_usage`.

Tolerance: the panel's `NC_016736.1` de-novo copy is 1 bp longer than its annotation. One base at
a copy boundary is the boundary disagreement K-07 measures and reports as agreement up to 2 bp, so
the checks below use the same 2 bp criterion rather than an exact match that one record cannot meet.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from plastid.core.genbank import read_fasta
from plastid.core.ir import detect_ir, enforce_ir_gate
from plastid.modules import structure
from plastid.run.context import RunContext
from plastid.run.layout import RunLayout

ROOT = Path(__file__).resolve().parents[2]
PANEL = ROOT / ".panel" / "ir-validation"
FASTAS = PANEL / "fasta"
EXPECTED = FASTAS / "EXPECTED.tsv"

# K-07's agreement criterion for two IR boundary calls: within 2 bp.
BOUNDARY_TOLERANCE_BP = 2

pytestmark = [
    pytest.mark.tier_b,
    pytest.mark.skipif(
        not EXPECTED.is_file(),
        reason="the rotated FASTA panel is not in this tree (.panel/ is gitignored)",
    ),
]


def _expected() -> list[dict[str, str]]:
    lines = EXPECTED.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"), strict=True)) for line in lines[1:]]


def _by_source() -> dict[str, list[dict[str, str]]]:
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in _expected():
        grouped.setdefault(row["source"], []).append(row)
    return grouped


@pytest.fixture(scope="module")
def runs(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """One run root for the module: each FASTA is analysed once, however many tests read it."""
    return tmp_path_factory.mktemp("structure")


_CACHE: dict[str, tuple[list[dict[str, object]], dict[str, object], int]] = {}


def _regions(path: Path, root: Path) -> tuple[list[dict[str, object]], dict[str, object], int]:
    """Run ``structure`` on one unannotated FASTA and return its regions, summary and IR length.

    Cached per file: the three checks below read the same 30 files, and a structure run costs
    about a second, which is the whole cost of this module.
    """
    key = str(path)
    if key in _CACHE:
        return _CACHE[key]
    genome = read_fasta(path)
    ir = enforce_ir_gate(detect_ir(genome, source="sequence"))
    layout = RunLayout.create(root, "structure", fingerprint="m11-17")
    layout.make_dirs()
    ctx = RunContext(layout=layout, params={}, module_id="structure", module_order=1)
    ctx.ir = ir
    result = structure.run(genome, {}, ctx)
    summary = result.tables["structure_summary"].rows[0]
    _CACHE[key] = (
        list(result.tables["regions"].rows),
        dict(summary),
        ir.ira.end - ir.ira.start,
    )
    return _CACHE[key]


def _lsc_origin(row: dict[str, str]) -> int:
    """0-based first base of the longer single-copy arc, from the panel's own copy geometry.

    Both copies are marked on a circle and the gaps between them measured; the longer gap is the
    LSC, and its first base is the normalized origin. A copy that wraps position 1 is marked in
    two pieces, which is why this is not four subtractions.
    """
    length = int(row["length"])
    covered = bytearray(length)
    for prefix in ("copyA", "copyB"):
        start = int(row[f"{prefix}_start_1based"]) - 1
        run = int(row[f"{prefix}_len"])
        for offset in range(run):
            covered[(start + offset) % length] = 1
    if all(covered):
        raise AssertionError(f"{row['name']}: the two copies cover the whole circle")
    # Longest circular run of uncovered bases.
    best_start, best_len, current = 0, 0, 0
    for index in range(2 * length):
        position = index % length
        if covered[position]:
            current = 0
        else:
            current += 1
            if current > best_len:
                best_len, best_start = current, (index - current + 1) % length
    return best_start


def _circular_distance(one: int, other: int, length: int) -> int:
    difference = abs(one - other) % length
    return min(difference, length - difference)


def test_a_rotation_does_not_move_the_region_table(runs: Path) -> None:
    """Five origins, one plastome: the region lengths, GC and skew are one set of numbers.

    This is the claim the rotated panel exists for. A number that moves is reading a raw position
    instead of the normalized frame, and no tolerance applies: a rotation changes no base.
    """
    for source, rows in _by_source().items():
        tables = []
        for row in rows:
            regions, _summary, _ir_length = _regions(FASTAS / f"{row['name']}.fa", runs)
            tables.append(
                [
                    (region["region"], region["length_bp"], region["gc_pct"], region["gc_skew"])
                    for region in regions
                ]
            )
        assert all(table == tables[0] for table in tables[1:]), (
            f"{source}: the region table moved with the origin — "
            + ", ".join(f"{row['name']}={table}" for row, table in zip(rows, tables, strict=True))
        )
        assert [region[0] for region in tables[0]] == ["LSC", "IRa", "SSC", "IRb"]


def test_the_regions_tile_the_genome_and_the_two_copies_agree(runs: Path) -> None:
    """Arithmetic the panel answer is not needed for: the four regions cover L bases exactly."""
    for row in _expected():
        regions, _summary, ir_length = _regions(FASTAS / f"{row['name']}.fa", runs)
        lengths = {region["region"]: region["length_bp"] for region in regions}
        assert sum(lengths.values()) == int(row["length"]), row["name"]
        assert lengths["IRa"] == lengths["IRb"] == ir_length, row["name"]
        assert lengths["LSC"] > lengths["SSC"] > 0, row["name"]


def test_the_denovo_ir_reproduces_the_panels_arithmetic(runs: Path) -> None:
    """The copy length and the normalized origin match ``EXPECTED.tsv`` within 2 bp.

    The expected values were derived from the deposited annotation, so this compares two
    independent derivations of the same boundary. The origin is what ``structure`` reports as
    ``rotation_offset`` — the GenBank-frame position that became position 0 — and it is checked
    here because a boundary that is right by luck of the start base would still move the frame.
    """
    for row in _expected():
        length = int(row["length"])
        regions, summary, ir_length = _regions(FASTAS / f"{row['name']}.fa", runs)
        expected_ir = int(row["ir_len"])
        assert abs(ir_length - expected_ir) <= BOUNDARY_TOLERANCE_BP, (
            f"{row['name']}: de-novo copy {ir_length} bp, panel says {expected_ir}"
        )
        # EXPECTED.tsv states the copies in the rotated file's own frame, so the origin derived
        # from them is already the one `structure` must report.
        expected_origin = _lsc_origin(row)
        assert _circular_distance(int(summary["rotation_offset"]), expected_origin, length) <= (
            BOUNDARY_TOLERANCE_BP
        ), f"{row['name']}: origin {summary['rotation_offset']}, expected {expected_origin}"
        assert {region["region"] for region in regions} == {"LSC", "IRa", "SSC", "IRb"}
        assert len(regions) == 4
