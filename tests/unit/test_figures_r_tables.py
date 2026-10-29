"""D5-03 for the six R figures: the file R reads equals the cells it came from.

Each fixture is the smallest table set its module emits, with the **real** column names taken from
the module sources, so a pivot that read the wrong column or re-derived a number fails here rather
than in a picture nobody checked. Every assertion is against the source table's own cells, and
where a cell is rendered it is compared with that column's own ``format`` — never with a number
this test computed for itself.
"""

from __future__ import annotations

import csv
import io

import pytest

from plastid.core.align import SUBSTITUTION_TYPES, substitution_class
from plastid.result.figures import (
    r_amino_acid,
    r_diversity,
    r_kaks,
    r_rscu,
    r_ssr,
    r_substitutions,
)
from plastid.result.result import Column, Result, Scope, Table


def _scope(label: str = "whole plastome") -> Scope:
    return Scope(label=label)


def _rows(text: str) -> list[list[str]]:
    return list(csv.reader(io.StringIO(text)))


def _tsv_rows(text: str) -> list[list[str]]:
    """F6 and F8 are tab-separated — the same split their own ``read.delim`` makes."""
    return list(csv.reader(io.StringIO(text), delimiter="\t"))


def _codon_usage_result() -> Result:
    rscu = Table(
        id="rscu",
        title="Relative synonymous codon usage",
        columns=(
            Column("codon", "Codon"),
            Column("amino_acid", "Amino acid"),
            Column("count", "Count", dtype="int"),
            Column("rscu", "RSCU", dtype="float", precision=4),
        ),
        rows=(
            {"codon": "GCA", "amino_acid": "A", "count": 370, "rscu": 1.0655},
            {"codon": "GCC", "amino_acid": "A", "count": 236, "rscu": 0.6796},
            {"codon": "TGG", "amino_acid": "W", "count": 100, "rscu": None},
        ),
        scope=_scope("CDS only, pseudogenes excluded; stops excluded"),
        key="codon",
    )
    summary = Table(
        id="codon_usage",
        title="Codon usage summary",
        columns=(
            Column("accession", "Accession"),
            Column("organism", "Organism"),
            Column("n_cds", "CDS in the set", dtype="int"),
            Column("enc", "ENC", dtype="float", precision=2),
        ),
        rows=(
            {"accession": "NC_050833.1", "organism": "Iris domestica", "n_cds": 86, "enc": 49.69},
        ),
        scope=_scope("CDS only, pseudogenes excluded; stops and partial codons not counted"),
    )
    return Result(
        module="codon_usage", module_version="0.1.0", tables={"rscu": rscu, "codon_usage": summary}
    )


def _amino_acid_result() -> Result:
    composition = Table(
        id="amino_acid_composition",
        title="Amino-acid composition",
        columns=(
            Column("residue", "Residue"),
            Column("count", "Count", dtype="int"),
            Column("percent", "Percent", dtype="float", unit="%", precision=2),
            Column("hydropathy", "Kyte–Doolittle hydropathy", dtype="float", precision=1),
        ),
        rows=(
            {"residue": "A", "count": 1399, "percent": 5.24, "hydropathy": 1.8},
            {"residue": "W", "count": 200, "percent": 0.75, "hydropathy": -0.9},
        ),
        scope=_scope("CDS only, pseudogenes excluded; percents over the residues counted"),
        key="residue",
    )
    summary = Table(
        id="amino_acid",
        title="Amino-acid summary",
        columns=(
            Column("accession", "Accession"),
            Column("organism", "Organism"),
            Column("n_cds", "CDS in the set", dtype="int"),
            Column("n_residues", "Residues counted", dtype="int"),
            Column("gravy", "GRAVY (Kyte–Doolittle)", dtype="float", precision=4),
        ),
        rows=(
            {
                "accession": "NC_050833.1",
                "organism": "Iris domestica",
                "n_cds": 86,
                "n_residues": 26696,
                "gravy": -0.1234,
            },
        ),
        scope=_scope("CDS only, pseudogenes excluded, both IR copies included"),
    )
    return Result(
        module="amino_acid",
        module_version="0.1.0",
        tables={"amino_acid": summary, "amino_acid_composition": composition},
    )


def _ssr_result() -> Result:
    summary = Table(
        id="ssr_summary",
        title="SSR summary",
        columns=(
            Column("accession", "Accession"),
            Column("total", "Total SSRs", dtype="int"),
            Column("mono", "Monomeric", dtype="int"),
            Column("compound", "In compound formation", dtype="int"),
            Column("thresholds", "Unit/repeat thresholds"),
        ),
        rows=(
            {
                "accession": "NC_050833.1",
                "total": 37,
                "mono": 31,
                "compound": 2,
                "thresholds": "(1/10) (2/6) (3/5)",
            },
        ),
        scope=_scope("whole plastome; MISA defaults; both IR copies included"),
    )
    by_region = Table(
        id="ssr_by_region",
        title="SSRs by plastid region",
        columns=(
            Column("region", "Region"),
            Column("total", "Total", dtype="int"),
            Column("mono", "Monomeric", dtype="int"),
            Column("di", "Dinucleotide", dtype="int"),
            Column("tri", "Trinucleotide", dtype="int"),
            Column("other", "Other unit sizes", dtype="int"),
        ),
        rows=(
            {"region": "LSC", "total": 24, "mono": 23, "di": 1, "tri": 0, "other": 0},
            {"region": "IRa", "total": 2, "mono": 2, "di": 0, "tri": 0, "other": 0},
            {"region": "SSC", "total": 4, "mono": 4, "di": 0, "tri": 0, "other": 0},
            {"region": "IRb", "total": 2, "mono": 2, "di": 0, "tri": 0, "other": 0},
        ),
        scope=_scope("whole plastome; counts from the per-SSR table"),
        key="region",
    )
    by_unit = Table(
        id="ssr_by_unit",
        title="SSRs by repeat unit",
        columns=(
            Column("unit", "Unit"),
            Column("unit_size", "Unit size", dtype="int", unit="bp"),
            Column("n", "Count", dtype="int"),
        ),
        # 31 mono + 1 di = the 32 rows the region table sums to, so the two agree by construction
        rows=({"unit": "mono", "unit_size": 1, "n": 31}, {"unit": "di", "unit_size": 2, "n": 1}),
        scope=_scope("whole plastome; counts from MISA .statistics"),
        key="unit",
    )
    ssrs = Table(
        id="ssrs",
        title="SSRs (positions from the MISA GFF)",
        columns=(
            Column("motif", "Motif"),
            Column("unit", "Unit"),
            Column("start_1b", "Start", dtype="int"),
            Column("end_1b", "End", dtype="int"),
            Column("region", "Region"),
            Column("context", "Context"),
        ),
        rows=(
            {
                "motif": "(T)10",
                "unit": "mono",
                "start_1b": 1,
                "end_1b": 10,
                "region": "LSC",
                "context": "IGS",
            },
            {
                "motif": "(A)11",
                "unit": "mono",
                "start_1b": 20,
                "end_1b": 30,
                "region": "LSC",
                "context": "Gene",
            },
            {
                "motif": "(AT)6",
                "unit": "di",
                "start_1b": 40,
                "end_1b": 51,
                "region": "SSC",
                "context": "Intron",
            },
        ),
        scope=_scope("GenBank coordinates; region split from the IR pair"),
    )
    return Result(
        module="ssr",
        module_version="0.4.0",
        tables={
            "ssr_summary": summary,
            "ssr_by_region": by_region,
            "ssr_by_unit": by_unit,
            "ssrs": ssrs,
        },
    )


def _diversity_result() -> Result:
    summary = Table(
        id="diversity_summary",
        title="Nucleotide diversity summary",
        columns=(
            Column("accession", "Reference accession"),
            Column("organism", "Reference organism"),
            Column("n_genomes", "Genomes aligned", dtype="int"),
            Column("window_bp", "Window", dtype="int", unit="bp"),
            Column("pi", "π (whole alignment)", dtype="float", precision=6),
        ),
        rows=(
            {
                "accession": "NC_050833.1",
                "organism": "Iris domestica",
                "n_genomes": 34,
                "window_bp": 600,
                "pi": 0.007332,
            },
        ),
        scope=_scope("whole aligned plastome, gaps dropped column-wise"),
    )
    genes = Table(
        id="diversity_genes",
        title="Nucleotide diversity by gene (reference annotation)",
        columns=(
            Column("gene", "Gene"),
            Column("kind", "Kind"),
            Column("usable_columns", "Usable columns", dtype="int"),
            Column("pi", "π", dtype="float", precision=6),
        ),
        rows=(
            {"gene": "psbA", "kind": "CDS", "usable_columns": 1062, "pi": 0.02343},
            {"gene": "rps16", "kind": "CDS", "usable_columns": 400, "pi": None},
        ),
        scope=_scope("reference genome's features mapped onto alignment columns; gaps dropped"),
        key="gene",
    )
    introns = Table(
        id="diversity_introns",
        title="Nucleotide diversity by intron (reference annotation)",
        columns=(
            Column("intron_id", "Intron (gene.copy.n)"),
            Column("gene", "Gene"),
            Column("region", "Region"),
            Column("length_bp", "Length", dtype="int", unit="bp"),
            Column("usable_columns", "Usable columns", dtype="int"),
            Column("pi", "π", dtype="float", precision=6),
        ),
        rows=(
            {
                "intron_id": "trnK-UUU.1.1",
                "gene": "trnK-UUU",
                "region": "LSC",
                "length_bp": 2500,
                "usable_columns": 2100,
                "pi": 0.07243,
            },
        ),
        scope=_scope("the introns the `introns` module publishes; gaps dropped"),
        key="intron_id",
    )
    igs = Table(
        id="diversity_igs",
        title="Nucleotide diversity by intergenic spacer (reference annotation)",
        columns=(
            Column("igs_id", "Spacer (left–right)"),
            Column("left", "Left gene"),
            Column("right", "Right gene"),
            Column("region", "Region"),
            Column("length_bp", "Length", dtype="int", unit="bp"),
            Column("usable_columns", "Usable columns", dtype="int"),
            Column("pi", "π", dtype="float", precision=6),
        ),
        rows=(
            {
                "igs_id": "psbA-matK",
                "left": "psbA",
                "right": "matK",
                "region": "LSC",
                "length_bp": 300,
                "usable_columns": 250,
                "pi": 0.08094,
            },
        ),
        scope=_scope("gaps between the exons of two neighbouring called genes"),
        key="igs_id",
    )
    return Result(
        module="diversity",
        module_version="0.2.0",
        tables={
            "diversity_summary": summary,
            "diversity_genes": genes,
            "diversity_introns": introns,
            "diversity_igs": igs,
        },
    )


def _substitutions_result() -> Result:
    summary = Table(
        id="substitutions_summary",
        title="Substitutions summary",
        columns=(
            Column("accession", "Reference accession"),
            Column("organism", "Reference organism"),
            Column("n_genomes", "Genomes compared", dtype="int"),
            Column("substitutions", "Substitutions", dtype="int"),
            Column("ts_tv", "Ts/Tv", dtype="float", precision=4),
        ),
        rows=(
            {
                "accession": "NC_050833.1",
                "organism": "Iris domestica",
                "n_genomes": 3,
                "substitutions": 40,
                "ts_tv": 2.6917,
            },
        ),
        scope=_scope("each genome against the reference; gap and ambiguous columns skipped"),
    )
    observed = {
        ("NC_024936.1", "A_to_G"): 5,
        ("NC_024936.1", "C_to_T"): 3,
        ("NC_024936.1", "A_to_C"): 2,
        ("NC_029218.1", "C_to_T"): 1,
    }
    by_type = Table(
        id="substitutions_by_type",
        title="Substitutions by type",
        columns=(
            Column("row_id", "Row"),
            Column("genome", "Genome"),
            Column("group", "Group"),
            Column("substitution_type", "Substitution type"),
            Column("class", "Class"),
            Column("count", "Count", dtype="int"),
        ),
        # All twelve types per genome, zeros included — exactly what the module emits.
        rows=tuple(
            {
                "row_id": f"{genome}|{type_name}",
                "genome": genome,
                "group": "",
                "substitution_type": type_name,
                "class": substitution_class(type_name),
                "count": observed.get((genome, type_name), 0),
            }
            for genome in ("NC_024936.1", "NC_029218.1")
            for type_name in SUBSTITUTION_TYPES
        ),
        scope=_scope("each genome against the reference; every one of the twelve directed changes"),
        key="row_id",
    )
    per_genome = Table(
        id="substitutions",
        title="Substitutions against the reference",
        columns=(
            Column("genome", "Genome"),
            Column("group", "Group"),
            Column("substitutions", "Substitutions", dtype="int"),
            Column("transitions", "Transitions", dtype="int"),
            Column("transversions", "Transversions", dtype="int"),
            Column("ts_tv", "Ts/Tv", dtype="float", precision=4),
        ),
        rows=(
            {
                "genome": "NC_024936.1",
                "group": "",
                "substitutions": 10,
                "transitions": 8,
                "transversions": 2,
                "ts_tv": 4.0,
            },
            {
                "genome": "NC_029218.1",
                "group": "",
                "substitutions": 1,
                "transitions": 1,
                "transversions": 0,
                "ts_tv": None,
            },
        ),
        scope=_scope("each genome against the reference; gap and ambiguous columns skipped"),
        key="genome",
    )
    return Result(
        module="substitutions",
        module_version="0.3.0",
        tables={
            "substitutions_summary": summary,
            "substitutions_by_type": by_type,
            "substitutions": per_genome,
        },
    )


def _kaks_result() -> Result:
    summary = Table(
        id="kaks_summary",
        title="Ka/Ks summary",
        columns=(
            Column("accession", "Reference accession"),
            Column("organism", "Reference organism"),
            Column("n_genomes", "Genomes compared", dtype="int"),
            Column("n_shared", "Genes shared with every annotated genome", dtype="int"),
        ),
        rows=(
            {
                "accession": "NC_050833.1",
                "organism": "Iris domestica",
                "n_genomes": 2,
                "n_shared": 79,
            },
        ),
        scope=_scope("the reference's CDS names; per-gene values are in kaks_genes"),
    )
    genes = Table(
        id="kaks_genes",
        title="Ka/Ks by gene",
        columns=(
            Column("row_id", "Row"),
            Column("gene", "Gene"),
            Column("genome", "Genome"),
            Column("group", "Group"),
            Column("ka_ks", "Ka/Ks", dtype="float", precision=4),
            Column("status", "Status"),
        ),
        rows=(
            {
                "row_id": "accD#NC_024936.1",
                "gene": "accD",
                "genome": "NC_024936.1",
                "group": "",
                "ka_ks": 0.5556,
                "status": "ok",
            },
            {
                "row_id": "atpA#NC_024936.1",
                "gene": "atpA",
                "genome": "NC_024936.1",
                "group": "",
                "ka_ks": 0.0953,
                "status": "ok",
            },
            {
                "row_id": "accD#NC_029218.1",
                "gene": "accD",
                "genome": "NC_029218.1",
                "group": "",
                "ka_ks": 0.5,
                "status": "ok",
            },
            {
                "row_id": "atpA#NC_029218.1",
                "gene": "atpA",
                "genome": "NC_029218.1",
                "group": "",
                "ka_ks": None,
                "status": "ks = 0",
            },
        ),
        scope=_scope("one row per gene and genome, both inverted-repeat copies counted once"),
        key="row_id",
    )
    return Result(
        module="kaks",
        module_version="0.3.0",
        tables={"kaks_summary": summary, "kaks_genes": genes},
    )


def test_the_rscu_file_is_the_rscu_table_with_one_species_column() -> None:
    result = _codon_usage_result()
    frame = r_rscu.frame({"codon_usage": result})
    _name, text = frame.inputs()[0]
    header, *rows = _rows(text)
    assert header == ["Amino_Acid", "Codon", "Iris_domestica"]
    table = result.tables["rscu"]
    assert [row[0] for row in rows] == ["Alanine", "Alanine", "Tryptophan"]
    assert [row[1] for row in rows] == [row["codon"] for row in table.rows]
    assert [row[2] for row in rows] == [table.format(f"{row['codon']}.rscu") for row in table.rows]
    # a blank RSCU is the column's own rendering of None — the empty string, never a 0
    assert rows[-1][2] == ""
    assert frame.n_cds == result.tables["codon_usage"].format("n_cds")


def test_the_amino_acid_file_is_the_composition_table_with_one_species_column() -> None:
    result = _amino_acid_result()
    frame = r_amino_acid.frame({"amino_acid": result})
    _name, text = frame.inputs()[0]
    header, *rows = _rows(text)
    assert header == ["AA_Code", "Amino_Acid", "Iris_domestica"]
    table = result.tables["amino_acid_composition"]
    assert [row[0] for row in rows] == [row["residue"] for row in table.rows]
    assert [row[1] for row in rows] == ["Alanine", "Tryptophan"]
    assert [row[2] for row in rows] == [
        table.format(f"{row['residue']}.percent") for row in table.rows
    ]
    assert frame.gravy == result.tables["amino_acid"].format("gravy")


def test_the_ssr_file_folds_four_tables_into_the_scripts_sixteen_columns() -> None:
    result = _ssr_result()
    frame = r_ssr.frame({"ssr": result})
    _name, text = frame.inputs()[0]
    header, row = _rows(text)
    assert header[0] == "Species" and len(header) == 16
    values = dict(zip(header, row, strict=True))
    assert values["Species"] == result.tables["ssr_summary"].value("accession")
    assert values["Total_SSRs"] == result.tables["ssr_summary"].format("total")
    by_region = {r["region"]: r for r in result.tables["ssr_by_region"].rows}
    assert values["LSC"] == str(by_region["LSC"]["total"])
    assert values["SSC"] == str(by_region["SSC"]["total"])
    # the one aggregation in this figure, asserted as such
    assert values["IR"] == str(by_region["IRa"]["total"] + by_region["IRb"]["total"])
    by_unit = {r["unit"]: r["n"] for r in result.tables["ssr_by_unit"].rows}
    assert values["Mono"] == str(by_unit["mono"]) and values["Di"] == str(by_unit["di"])
    assert values["Tri"] == "0"  # absent from ssr_by_unit: no trinucleotide SSR, a real zero
    contexts = [r["context"] for r in result.tables["ssrs"].rows]
    assert values["Gene"] == str(contexts.count("Gene"))
    assert values["Intron"] == str(contexts.count("Intron"))
    assert values["IGS"] == str(contexts.count("IGS"))
    # unsourced, so blank — never a zero the toolkit did not measure
    assert values["tRNA"] == "" and values["rRNA"] == ""
    assert "tRNA, rRNA are blank" in frame.caption()


def test_the_diversity_frames_are_the_two_files_the_script_reads_by_name() -> None:
    result = _diversity_result()
    frame = r_diversity.frame({"diversity": result})
    names = [name for name, _text in frame.inputs()]
    assert names == ["genes_for_R_plot.txt", "noncoding_for_R_plot.txt"]
    files = dict(frame.inputs())
    genes_header, *gene_rows = _tsv_rows(files["genes_for_R_plot.txt"])
    assert genes_header == ["Region", "Value"]
    table = result.tables["diversity_genes"]
    assert gene_rows == [[row["gene"], table.format(f"{row['gene']}.pi")] for row in table.rows]
    assert gene_rows[1][1] == ""  # a blank π stays blank
    noncoding = _tsv_rows(files["noncoding_for_R_plot.txt"])[1:]
    introns = result.tables["diversity_introns"]
    spacers = result.tables["diversity_igs"]
    assert noncoding == [
        [row["intron_id"], introns.format(f"{row['intron_id']}.pi")] for row in introns.rows
    ] + [[row["igs_id"], spacers.format(f"{row['igs_id']}.pi")] for row in spacers.rows]
    assert frame.n_genomes == result.tables["diversity_summary"].format("n_genomes")


def test_the_diversity_labels_are_unique_because_the_script_halts_on_a_repeat() -> None:
    """Measured 17/10/2026 on the golden pair: ``factor(Region, levels = Region)`` stops the run.

    R reports ``factor level [6] is duplicated`` and the figure is not drawn at all. A label built
    from the gene name collides in two ways the name cannot show — the two IR copies of one intron,
    and the two introns of one gene — so the labels are the tables' own keys.
    """
    result = _diversity_result()
    introns = result.tables["diversity_introns"]
    # the real shape: `rps12.1.2`/`rps12.2.2` are one intron in two IR copies, while
    # `clpP.1.1`/`clpP.1.2` are two different introns of one gene
    duplicated = Table(
        id=introns.id,
        title=introns.title,
        columns=introns.columns,
        rows=tuple(
            {**row, "intron_id": key, "gene": gene}
            for key, gene in (
                ("rps12.1.2", "rps12"),
                ("rps12.2.2", "rps12"),
                ("clpP.1.1", "clpP"),
                ("clpP.1.2", "clpP"),
            )
            for row in (introns.rows[0],)
        ),
        scope=introns.scope,
        key=introns.key,
    )
    shaped = Result(
        module="diversity",
        module_version="0.2.0",
        tables={**result.tables, "diversity_introns": duplicated},
    )
    labels = [label for label, _value in r_diversity.frame({"diversity": shaped}).noncoding]
    assert len(labels) == len(set(labels)), f"the script halts on a repeated label: {labels}"


def test_the_substitution_file_is_the_type_table_pivoted_with_a_ratio_row() -> None:
    result = _substitutions_result()
    frame = r_substitutions.frame({"substitutions": result})
    _name, text = frame.inputs()[0]
    header, *rows = _rows(text)
    assert header[:2] == ["Substitution_Type", "Category"]
    assert header[2:] == ["NC_024936.1", "NC_029218.1"]
    by_type = {
        (row["genome"], row["substitution_type"]): row
        for row in result.tables["substitutions_by_type"].rows
    }
    types = [row[0] for row in rows[:-1]]
    assert len(types) == 12 and len(set(types)) == 12
    assert all(row[1] == "Substitution" for row in rows[:-1])
    for row in rows[:-1]:
        assert row[2] == str(by_type[("NC_024936.1", row[0])]["count"])
        assert row[3] == str(by_type[("NC_029218.1", row[0])]["count"])
    ratio = rows[-1]
    assert ratio[0] == "Ts/Tv Ratio" and ratio[1] == "Ratio"
    assert ratio[2] == result.tables["substitutions"].format("NC_024936.1.ts_tv")
    assert ratio[3] == ""  # no transversion: undefined, so blank
    assert frame.reference == result.tables["substitutions_summary"].value("accession")


def test_the_kaks_file_is_the_gene_by_genome_matrix() -> None:
    result = _kaks_result()
    frame = r_kaks.frame({"kaks": result})
    _name, text = frame.inputs()[0]
    header, *rows = _tsv_rows(text)
    assert header == ["Gene", "NC_024936.1", "NC_029218.1"]
    table = result.tables["kaks_genes"]
    cells = {
        (row["gene"], row["genome"]): table.format(f"{row['row_id']}.ka_ks") for row in table.rows
    }
    assert [row[0] for row in rows] == ["accD", "atpA"]
    for row in rows:
        assert row[1] == cells[(row[0], "NC_024936.1")]
        assert row[2] == cells[(row[0], "NC_029218.1")]
    assert rows[1][2] == ""  # Ks = 0: undefined, blank, never 0 or 1
    assert frame.reference == result.tables["kaks_summary"].value("accession")


@pytest.mark.parametrize(
    ("spec", "module", "result"),
    [
        (r_rscu.SPEC, "codon_usage", _codon_usage_result),
        (r_amino_acid.SPEC, "amino_acid", _amino_acid_result),
        (r_ssr.SPEC, "ssr", _ssr_result),
        (r_diversity.SPEC, "diversity", _diversity_result),
        (r_substitutions.SPEC, "substitutions", _substitutions_result),
        (r_kaks.SPEC, "kaks", _kaks_result),
    ],
    ids=lambda value: getattr(value, "id", ""),
)
def test_a_missing_module_is_named_rather_than_drawn_empty(spec, module, result) -> None:
    assert spec.missing({}) == tuple(f"{m}.{t}" for m, t in spec.requires)
    assert spec.missing({module: result()}) == ()


def test_every_caption_carries_the_scope_label_of_every_table_it_read() -> None:
    """The Done-when: an exported number carries its scope label; an R figure carries it here."""
    for spec, module, build in (
        (r_rscu.SPEC, "codon_usage", _codon_usage_result),
        (r_amino_acid.SPEC, "amino_acid", _amino_acid_result),
        (r_ssr.SPEC, "ssr", _ssr_result),
        (r_diversity.SPEC, "diversity", _diversity_result),
        (r_substitutions.SPEC, "substitutions", _substitutions_result),
        (r_kaks.SPEC, "kaks", _kaks_result),
    ):
        results = {module: build()}
        lines = spec.scope_lines(results)
        assert len(lines) == len(spec.requires)
        for line, (mod, table_name) in zip(lines, spec.requires, strict=True):
            assert results[mod].tables[table_name].scope.label in line
            assert line.startswith(f"{mod}.{table_name} — scope: ")


def test_a_species_outside_the_paper_list_is_drawn_not_skipped() -> None:
    """M11.7: F3 reads the species from the input, so there is no list to be outside of."""
    result = _codon_usage_result()
    summary = result.tables["codon_usage"]
    other = Result(
        module="codon_usage",
        module_version="0.1.0",
        tables={
            "rscu": result.tables["rscu"],
            "codon_usage": Table(
                id="codon_usage",
                title=summary.title,
                columns=summary.columns,
                rows=({**summary.rows[0], "organism": "Coleus amboinicus"},),
                scope=summary.scope,
            ),
        },
    )
    frame = r_rscu.frame({"codon_usage": other})
    assert frame.species == "Coleus_amboinicus"
    assert "every species column" in frame.caption()


def test_the_amino_acid_pivot_also_draws_a_species_outside_the_paper_list() -> None:
    """M11.7: F4 reads the species from the input, exactly as F3 does."""
    result = _amino_acid_result()
    summary = result.tables["amino_acid"]
    other = Result(
        module="amino_acid",
        module_version="0.1.0",
        tables={
            "amino_acid_composition": result.tables["amino_acid_composition"],
            "amino_acid": Table(
                id="amino_acid",
                title=summary.title,
                columns=summary.columns,
                rows=({**summary.rows[0], "organism": "Coleus amboinicus"},),
                scope=summary.scope,
            ),
        },
    )
    frame = r_amino_acid.frame({"amino_acid": other})
    assert frame.species == "Coleus_amboinicus"
    assert "every species column" in frame.caption()


def test_the_kaks_pivot_refuses_a_one_genome_matrix_by_name() -> None:
    """Measured 17/10/2026: ``R_F8`` draws two species columns and halts on one.

    Its own ``mat[keep_genes, keep_spp]`` drops a one-column matrix to a vector, so ``rownames``
    is NULL and R stops with "replacement has 0 rows, data has 79". The toolkit refuses it by name
    rather than letting a cryptic R error be the explanation.
    """
    from plastid.result.rfigure import RFigureSkip

    result = _kaks_result()
    table = result.tables["kaks_genes"]
    one = Result(
        module="kaks",
        module_version="0.3.0",
        tables={
            "kaks_summary": result.tables["kaks_summary"],
            "kaks_genes": Table(
                id=table.id,
                title=table.title,
                columns=table.columns,
                rows=tuple(row for row in table.rows if row["genome"] == "NC_024936.1"),
                scope=table.scope,
                key=table.key,
            ),
        },
    )
    with pytest.raises(RFigureSkip, match="one-genome matrix"):
        r_kaks.frame({"kaks": one})
    # …and two genomes still draw, so the refusal is about the count and nothing else
    assert len(r_kaks.frame({"kaks": result}).genomes) == 2


def test_the_substitution_pivot_names_both_tables_when_no_ratio_table_is_present() -> None:
    from plastid.result.rfigure import RFigureSkip

    result = _substitutions_result()
    without = Result(
        module="substitutions",
        module_version="0.3.0",
        tables={
            "substitutions_summary": result.tables["substitutions_summary"],
            "substitutions_by_type": result.tables["substitutions_by_type"],
        },
    )
    with pytest.raises(RFigureSkip, match="substitutions_by_genome"):
        r_substitutions.frame({"substitutions": without})
