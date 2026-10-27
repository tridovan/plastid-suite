"""Wide figure frames for a batch: one column per sample for F3/F4, one row per sample for F5.

The single-run pivots read one ``Result`` per module and emit one species column (F3/F4) or one
row (F5). A batch stores the same tables with a ``sample`` column, so these frames split that long
table back into per-sample results, reuse the single-run pivot (the same semantics, not a second
reading of the cells), and combine the pieces into the shape the R scripts already accept.

Only the figures whose script is sample-oriented are built here. F7/F8 remain one-run comparisons
against a declared reference and are not fabricated from independent samples.
"""

from __future__ import annotations

from collections.abc import Mapping

from plastid.result.figures import r_amino_acid, r_rscu, r_ssr
from plastid.result.figures.r_common import csv_text
from plastid.result.result import Result, Table

__all__ = ["frame_for"]

DROPPED = ("sample", "species")


def _key_column(table: Table, rows, sample: str) -> str | None:
    """Recover the row key the single-run table addressed its rows by.

    The batch's long table carries each row's key as ``row_key``; the original key column is the
    one whose cell equals it. A summary table (one row) has no key and gets none.
    """
    keyed = [row for row in table.rows if row.get("sample") == sample and row.get("row_key")]
    if not keyed:
        return None
    for column in table.columns:
        name = column.name
        if name in DROPPED + ("row_key",):
            continue
        if all(str(row.get(name, "")) == str(row["row_key"]) for row in keyed):
            values = [str(row.get(name, "")) for row in rows]
            if len(values) == len(set(values)):
                return name
    return None


def _per_sample(results: Mapping[str, Result], module_id: str) -> dict[str, Result]:
    """Split a batch's module table into one single-run-shaped ``Result`` per sample."""
    result = results.get(module_id)
    if result is None:
        return {}
    samples = sorted(
        {
            str(row["sample"])
            for table in result.tables.values()
            for row in table.rows
            if "sample" in row
        }
    )
    split: dict[str, Result] = {}
    for sample in samples:
        tables: dict[str, Table] = {}
        for table_id, table in result.tables.items():
            rows = tuple(
                {name: value for name, value in row.items() if name not in DROPPED + ("row_key",)}
                for row in table.rows
                if row.get("sample") == sample
            )
            if not rows:
                continue
            key = _key_column(table, rows, sample)
            tables[table_id] = Table(
                id=table_id,
                title=table.title,
                columns=tuple(c for c in table.columns if c.name not in DROPPED + ("row_key",)),
                rows=rows,
                scope=table.scope,
                key=key,
            )
        split[sample] = Result(module=module_id, module_version="batch", tables=tables)
    return split


def _unique(labels: list[str]) -> list[str]:
    seen: dict[str, int] = {}
    out: list[str] = []
    for label in labels:
        seen[label] = seen.get(label, 0) + 1
        out.append(label if seen[label] == 1 else f"{label}#{seen[label]}")
    return out


def frame_for(spec_id: str, results: Mapping[str, Result]) -> object | None:
    """A batch frame for F3/F4/F5, or ``None`` when the results are a single run."""
    if not any(
        "sample" in {column.name for column in table.columns}
        for result in results.values()
        for table in result.tables.values()
    ):
        return None
    if spec_id == "r_rscu":
        return _rscu(results)
    if spec_id == "r_amino_acid":
        return _amino_acid(results)
    if spec_id == "r_ssr":
        return _ssr(results)
    return None


class _BatchRscu:
    def __init__(self, samples) -> None:
        self.samples = samples  # [(species, {codon: (name, rscu)})]
        self.order = tuple(codon for codon in samples[0][1]) if samples else ()

    def inputs(self):
        header = ("Amino_Acid", "Codon", *(species for species, _ in self.samples))
        rows = tuple(
            (
                self.samples[0][1][codon][0],
                codon,
                *(values.get(codon, ("", ""))[1] for _species, values in self.samples),
            )
            for codon in self.order
        )
        return (("Complete_Codon_Usage_Analysis.csv", csv_text(header, rows)),)

    def caption(self) -> str:
        return (
            f"{len(self.order)} sense codons from {len(self.samples)} samples, one column per "
            "sample; every value is the sample's own RSCU cell. The script draws every species "
            "column in its input, in that order or in --species-order."
        )


def _rscu(results):
    split = _per_sample(results, "codon_usage")
    samples = []
    for _sample, result in split.items():
        frame = r_rscu.frame({"codon_usage": result})
        samples.append(
            (frame.species, {codon: (name, value) for name, codon, value in frame.codons})
        )
    labels = _unique([species for species, _ in samples])
    return _BatchRscu(list(zip(labels, [values for _, values in samples], strict=True)))


class _BatchAmino:
    def __init__(self, samples) -> None:
        self.samples = samples  # [(species, {code: (name, percent)})]
        self.order = tuple(code for code in samples[0][1]) if samples else ()

    def inputs(self):
        header = ("AA_Code", "Amino_Acid", *(species for species, _ in self.samples))
        rows = tuple(
            (
                code,
                self.samples[0][1][code][0],
                *(values.get(code, ("", ""))[1] for _species, values in self.samples),
            )
            for code in self.order
        )
        return (("amino_acid_composition.csv", csv_text(header, rows)),)

    def caption(self) -> str:
        return (
            f"{len(self.order)} residues from {len(self.samples)} samples, one column per sample; "
            "every value is the sample's own composition cell."
        )


def _amino_acid(results):
    split = _per_sample(results, "amino_acid")
    samples = []
    for _sample, result in split.items():
        frame = r_amino_acid.frame({"amino_acid": result})
        samples.append(
            (frame.species, {code: (name, value) for code, name, value in frame.residues})
        )
    labels = _unique([species for species, _ in samples])
    return _BatchAmino(list(zip(labels, [values for _, values in samples], strict=True)))


class _BatchSsr:
    def __init__(self, rows, unsourced) -> None:
        self.rows = rows
        self.unsourced = unsourced

    def inputs(self):
        header = (
            "Species",
            "Total_SSRs",
            "LSC",
            "SSC",
            "IR",
            "Mono",
            "Di",
            "Tri",
            "Tetra",
            "Penta",
            "Hexa",
            "Gene",
            "tRNA",
            "rRNA",
            "Intron",
            "IGS",
        )
        return (("ssr_summary_for_r.csv", csv_text(header, tuple(self.rows))),)

    def caption(self) -> str:
        return (
            f"{len(self.rows)} samples, one row per sample; every value is that sample's own SSR "
            f"cell. {', '.join(self.unsourced)} are blank because PTS records an SSR's context as "
            "Gene/Intron/IGS only — blank is 'not measured', never zero."
        )


def _ssr(results):
    split = _per_sample(results, "ssr")
    rows, unsourced = [], ()
    for _sample, result in split.items():
        frame = r_ssr.frame({"ssr": result})
        rows.append(frame.values)
        unsourced = frame.unsourced
    return _BatchSsr(tuple(rows), unsourced)
