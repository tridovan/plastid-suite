"""Paper Figure 8 — the Ka/Ks heatmap, drawn by ``R_F8_KaKs_heatmap.R``.

The script reads a **matrix**: genes down the first column, one column per genome, Ka/Ks in the
cells (``read.delim(..., row.names = 1)``). PTS's ``kaks.kaks_genes`` is long — one row per gene
*and* genome, keyed ``gene#genome`` — so the pivot here is long→wide.

Like Figure 7 this comes from an **arity-2** module, so the columns are every query genome in the
run and the values are each of them against the reference (``genomes[0]``), which the caption names
because the reference is not a column.

**A blank stays blank.** ``kaks`` publishes a row for every gene × genome it considered, including
the ones it could not score (``status`` says why), and the module's own note is that a gene with
Ks = 0 has an undefined Ka/Ks — blank, never zero or one. The heatmap's grey cells are those
blanks, so filling them would quietly redraw the figure.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

from plastid.result.figures.r_common import tsv_text
from plastid.result.result import Result
from plastid.result.rfigure import RFigureSkip, RFigureSpec

_GENES = "kaks_genes"
_SUMMARY = "kaks_summary"
_FILENAME = "kaks_matrix.tsv"


@dataclass(frozen=True, slots=True)
class KaksFrame:
    """The gene × genome matrix F8 reads."""

    genes: tuple[str, ...]
    genomes: tuple[str, ...]
    cells: dict[tuple[str, str], str]  # (gene, genome) → the ka_ks column's rendering
    reference: str
    n_scored: str

    def inputs(self) -> tuple[tuple[str, str], ...]:
        header = ("Gene", *self.genomes)
        rows = [
            (gene, *[self.cells.get((gene, genome), "") for genome in self.genomes])
            for gene in self.genes
        ]
        return ((_FILENAME, tsv_text(header, rows)),)

    def caption(self) -> str:
        return (
            f"{len(self.genes)} genes × {len(self.genomes)} genome(s) against the reference "
            f"{self.reference}; {self.n_scored} gene/genome pairs scored. "
            "A blank cell is a Ka/Ks this run could not compute (see the kaks_genes status "
            "column) — undefined, never zero or one."
        )


def frame(results: Mapping[str, Result]) -> KaksFrame:
    """``kaks_genes`` pivoted to gene × genome — every cell the ``ka_ks`` column's rendering."""
    result = results["kaks"]
    table = result.tables[_GENES]
    summary = result.tables[_SUMMARY]
    genes = tuple(dict.fromkeys(row["gene"] for row in table.rows))
    genomes = tuple(dict.fromkeys(row["genome"] for row in table.rows))
    if len(genomes) < 2:
        # Measured 17/10/2026, isolated on a hand-built matrix: two species columns draw and one
        # does not. The script's own `mat[keep_genes, keep_spp]` drops a one-column matrix to a
        # vector, so `rownames(mat)` is NULL and it halts with "replacement has 0 rows". This is
        # the script's limit, not the run's — a second genome in the run draws the figure.
        raise RFigureSkip(
            f"R_F8_KaKs_heatmap.R cannot draw a one-genome matrix — it drops a single species "
            f"column to a vector and halts (measured 17/10/2026). This run compares the reference "
            f"against {len(genomes)} genome(s); add another to the run, or read the kaks_genes "
            f"table directly"
        )
    cells = {
        (row["gene"], row["genome"]): table.format(f"{row['row_id']}.ka_ks") for row in table.rows
    }
    scored = sum(1 for row in table.rows if row["status"] == "ok")
    return KaksFrame(
        genes=genes,
        genomes=genomes,
        cells=cells,
        reference=summary.value("accession") or summary.value("organism"),
        n_scored=str(scored),
    )


SPEC = RFigureSpec(
    id="r_kaks",
    title="Ka/Ks by gene and genome (paper figure 8, drawn by R)",
    requires=(("kaks", _GENES), ("kaks", _SUMMARY)),
    script="kaks",
    input_mode="file",
    frame=frame,
)
