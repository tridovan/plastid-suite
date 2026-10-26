"""Group comparison statistics over per-locus measurements (P5-03).

The comparison this module exists to make: take every CDS locus, measure it in every genome, and
ask whether the groups the run declared (``--groups``) differ at that locus — Mann–Whitney U for
two groups, Kruskal–Wallis H for three or more — then correct across the whole family of tests
with Benjamini–Hochberg, and carry Bonferroni beside it for comparison.

Four measurements, all computed from the sequence the record carries, none of them needing an
alignment: ``length_bp``, ``gc_pct``, ``gc3_pct`` (GC at third codon positions) and ``enc``
(Wright's effective number of codons — the same :func:`plastid.core.codon.enc` the ``codon_usage``
module calls, so the two cannot disagree).

**One observation per locus per genome.** A plastome carries both IR copies of a locus; counting
both would be pseudo-replication, and the two copies are near-identical by the IR gate rather than
by chance. Loci are therefore collapsed by their **normalized name**, first copy wins, and how
many rows were collapsed is reported rather than hidden.

A test the data cannot support — a group with one genome in it, which is what the mito panel has —
is reported *not computable* with its reason, never as a p-value.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace

from plastid.core import stats as statistics
from plastid.core.codon import codon_counts, enc
from plastid.core.gene_names import load_gene_aliases, normalize_name
from plastid.core.model import Genome
from plastid.core.organelle import GROUP_SHARED
from plastid.modules.registry import module
from plastid.result.result import Column, Result, Scope, Table
from plastid.run.context import RunContext

MODULE_ID = "stats"
MODULE_VERSION = "0.1.0"

MEASURES = ("length_bp", "gc_pct", "gc3_pct", "enc")
_MEASURE_LABELS = {
    "length_bp": ("Length", "bp", 0),
    "gc_pct": ("GC", "%", 2),
    "gc3_pct": ("GC3", "%", 2),
    "enc": ("ENC", None, 2),
}
DEFAULT_ALPHA = 0.05
_THIRD_POSITIONS = 3


@dataclass(frozen=True, slots=True)
class Observation:
    """One locus measured in one genome — the unit every test in this module ranks."""

    locus: str
    group: str
    length_bp: int
    gc_pct: float
    gc3_pct: float | None
    enc: float | None

    def value(self, measure: str) -> float | None:
        return getattr(self, measure)


@dataclass(frozen=True, slots=True)
class TestRow:
    """One locus × measure test, carrying everything the sentence and the table need."""

    locus: str
    measure: str
    status: str
    reason: str
    n: tuple[int, ...]
    medians: tuple[float | None, ...]
    statistic: float | None
    p: float | None
    assumption: str
    p_bh: float | None = None
    p_bonferroni: float | None = None
    significant: str = "no"


def _sequence_of(genome: Genome, gene: object) -> str:
    """The CDS as the record spells it, parts joined — exon order, not genomic order."""
    parts = getattr(gene, "parts", ())
    if not parts:
        return ""
    return "".join(genome.sequence[part.start : part.end] for part in parts).upper()


def _gc_pct(sequence: str) -> float | None:
    if not sequence:
        return None
    gc = sum(1 for base in sequence if base in ("G", "C"))
    return 100.0 * gc / len(sequence)


def _gc3_pct(sequence: str) -> float | None:
    """GC at third codon positions — the position where silent change happens."""
    third = sequence[_THIRD_POSITIONS - 1 :: _THIRD_POSITIONS]
    if not third:
        return None
    gc = sum(1 for base in third if base in ("G", "C"))
    return 100.0 * gc / len(third)


def measure(
    genomes: Sequence[Genome], groups: Mapping[str, str]
) -> tuple[list[Observation], int, tuple[str, ...]]:
    """Every CDS locus, measured once per genome.

    Returns the observations, how many CDS rows were collapsed to one per locus, and the genomes
    that contributed **no CDS at all** — measured on the stand-in panel: 3 of its 11 records carry
    no annotation, so a summary that said "11 genomes" would overstate the comparison by a third.
    """
    aliases = load_gene_aliases()
    observations: list[Observation] = []
    collapsed = 0
    unannotated: list[str] = []
    for genome in genomes:
        key = genome.accession or genome.name
        if key not in groups:
            raise ValueError(
                f"stats: no group for {key!r} — every genome needs one, or the comparison "
                f"would mix labelled and unlabelled observations"
            )
        seen: set[str] = set()
        contributed = False
        for gene in genome.genes:
            if gene.kind != "CDS" or gene.pseudo:
                continue
            locus = normalize_name(gene.name, aliases).canonical
            if locus in seen:
                collapsed += 1  # the IR's second copy: one observation per locus, not two
                continue
            seen.add(locus)
            contributed = True
            sequence = _sequence_of(genome, gene)
            observations.append(
                Observation(
                    locus=locus,
                    group=groups[key],
                    length_bp=gene.length,
                    gc_pct=_gc_pct(sequence) or 0.0,
                    gc3_pct=_gc3_pct(sequence),
                    enc=enc(codon_counts(sequence)),
                )
            )
        if not contributed:
            unannotated.append(key)
    return observations, collapsed, tuple(unannotated)


def _test_rows(observations: Sequence[Observation], measures: Sequence[str]) -> list[TestRow]:
    """One test per locus × measure, with the correction applied afterwards."""
    by_locus: dict[str, list[Observation]] = {}
    for observation in observations:
        by_locus.setdefault(observation.locus, []).append(observation)

    rows: list[TestRow] = []
    for locus in sorted(by_locus):
        members = by_locus[locus]
        groups: dict[str, list[Observation]] = {}
        for observation in members:
            groups.setdefault(observation.group, []).append(observation)
        labels = sorted(groups)
        for name in measures:
            values = {
                label: [
                    value for entry in groups[label] if (value := entry.value(name)) is not None
                ]
                for label in labels
            }
            medians = tuple(
                statistics.describe(values[label]).median if values[label] else None
                for label in labels
            )
            if len(labels) < 2:
                rows.append(
                    TestRow(
                        locus=locus,
                        measure=name,
                        status=statistics.NOT_COMPUTABLE,
                        reason="only one group carries this locus "
                        f"({labels[0] if labels else '—'})",
                        n=tuple(len(values[label]) for label in labels),
                        medians=medians,
                        statistic=None,
                        p=None,
                        assumption="—",
                    )
                )
                continue
            outcome = (
                statistics.mann_whitney(values[labels[0]], values[labels[1]])
                if len(labels) == 2
                else statistics.kruskal_wallis(*[values[label] for label in labels])
            )
            rows.append(
                TestRow(
                    locus=locus,
                    measure=name,
                    status=outcome.status,
                    reason=outcome.reason,
                    n=outcome.n,
                    medians=medians,
                    statistic=outcome.statistic,
                    p=outcome.p,
                    assumption=outcome.assumption,
                )
            )
    adjusted_bh = statistics.benjamini_hochberg([row.p for row in rows])
    adjusted_bonferroni = statistics.bonferroni([row.p for row in rows])
    return [
        replace(row, p_bh=bh, p_bonferroni=bonferroni)
        for row, bh, bonferroni in zip(rows, adjusted_bh, adjusted_bonferroni, strict=True)
    ]


def _marked(rows: list[TestRow], alpha: float) -> list[TestRow]:
    return [replace(row, significant="yes" if _below(row.p_bh, alpha) else "no") for row in rows]


def _below(value: float | None, alpha: float) -> bool:
    return value is not None and value < alpha


@module(
    MODULE_ID,
    group=GROUP_SHARED,
    version=MODULE_VERSION,
    title="Group comparison statistics",
    arity=2,
    citations=("biopython",),
)
def run(genomes: Sequence[Genome], params: dict, ctx: RunContext) -> Result:
    groups = params.get("groups")
    if not groups:
        raise ValueError(
            "stats: no groups — this module compares groups of genomes, so pass "
            "--groups FILE (a YAML mapping {group: [accessions]})"
        )
    alpha = float(params.get("alpha", DEFAULT_ALPHA))
    measures = tuple(params.get("measures", MEASURES))
    unknown = [name for name in measures if name not in MEASURES]
    if unknown:
        raise ValueError(f"stats: unknown measure(s) {unknown} — have {list(MEASURES)}")

    observations, collapsed, unannotated = measure(genomes, groups)
    if not observations:
        raise ValueError("stats: no CDS in any input — there is nothing to compare")
    rows = _marked(_test_rows(observations, measures), alpha)

    labels = sorted({observation.group for observation in observations})
    tested = [row for row in rows if row.status == statistics.OK]
    significant = [row for row in rows if row.significant == "yes"]
    family = "Mann–Whitney U" if len(labels) == 2 else "Kruskal–Wallis H"
    scope = Scope(label="per-locus measurements; one observation per locus per genome")

    summary = Table(
        id="stats_summary",
        title="Group comparison",
        columns=(
            Column("accession", "Reference"),
            Column("organism", "Organism"),
            Column("n_genomes", "Genomes", dtype="int"),
            Column("n_groups", "Groups", dtype="int"),
            Column("groups", "Group labels"),
            Column("n_loci", "Loci", dtype="int"),
            Column("n_measures", "Measures", dtype="int"),
            Column("measures", "Measure set"),
            Column("test", "Test"),
            Column("correction", "Correction"),
            Column("alpha", "α", dtype="float", precision=4),
            Column("n_tested", "Tests with a p-value", dtype="int"),
            Column("n_not_computable", "Not computable", dtype="int"),
            Column("n_significant", "Significant after correction", dtype="int"),
            Column("n_observations", "Observations", dtype="int"),
            Column("n_collapsed_cds", "CDS collapsed to one per locus", dtype="int"),
            Column("n_genomes_with_cds", "Genomes carrying CDS", dtype="int"),
            Column("n_genomes_without_cds", "Genomes contributing no CDS", dtype="int"),
            Column("genomes_without_cds", "Genomes with no CDS"),
        ),
        rows=(
            {
                "accession": genomes[0].accession or genomes[0].name,
                "organism": genomes[0].name,
                "n_genomes": len(genomes),
                "n_groups": len(labels),
                "groups": ", ".join(labels),
                "n_loci": len({observation.locus for observation in observations}),
                "n_measures": len(measures),
                "measures": ", ".join(measures),
                "test": family,
                "correction": "Benjamini–Hochberg",
                "alpha": alpha,
                "n_tested": len(tested),
                "n_not_computable": len(rows) - len(tested),
                "n_significant": len(significant),
                "n_observations": len(observations),
                "n_collapsed_cds": collapsed,
                "n_genomes_with_cds": len(genomes) - len(unannotated),
                "n_genomes_without_cds": len(unannotated),
                "genomes_without_cds": ", ".join(unannotated) or "—",
            },
        ),
        scope=scope,
    )

    descriptives = Table(
        id="stats_descriptives",
        title="Descriptives per measure and group",
        columns=(
            Column("observation", "Measure · group"),
            Column("measure", "Measure"),
            Column("group", "Group"),
            Column("n", "n", dtype="int"),
            Column("mean", "Mean", dtype="float", precision=4),
            Column("sd", "SD", dtype="float", precision=4),
            Column("median", "Median", dtype="float", precision=4),
            Column("q1", "Q1", dtype="float", precision=4),
            Column("q3", "Q3", dtype="float", precision=4),
            Column("minimum", "Min", dtype="float", precision=4),
            Column("maximum", "Max", dtype="float", precision=4),
        ),
        rows=tuple(
            _descriptive_row(name, label, _values(observations, name, label))
            for name in measures
            for label in labels
        ),
        scope=scope,
        key="observation",
    )

    tests = Table(
        id="stats_tests",
        title="Per-locus tests",
        columns=(
            Column("test_id", "Test"),
            Column("locus", "Locus"),
            Column("measure", "Measure"),
            Column("n_per_group", "n per group"),
            Column("medians", "Medians per group"),
            Column("statistic", "Statistic", dtype="float", precision=4),
            Column("p", "p", dtype="float", precision=6),
            Column("p_bh", "p (BH)", dtype="float", precision=6),
            Column("p_bonferroni", "p (Bonferroni)", dtype="float", precision=6),
            Column("significant", "Significant"),
            Column("status", "Status"),
            Column("note", "Note"),
        ),
        rows=tuple(
            {
                "test_id": f"{row.locus}:{row.measure}",
                "locus": row.locus,
                "measure": row.measure,
                "n_per_group": " · ".join(
                    f"{label}: {n}" for label, n in zip(labels, row.n, strict=True)
                ),
                "medians": " · ".join(_number(value) for value in row.medians),
                "statistic": row.statistic,
                "p": row.p,
                "p_bh": row.p_bh,
                "p_bonferroni": row.p_bonferroni,
                "significant": row.significant,
                "status": row.status,
                "note": row.reason,
            }
            for row in rows
        ),
        scope=scope,
        key="test_id",
    )

    correlations = _correlations(observations, measures, labels, scope)
    result = Result(
        module=MODULE_ID,
        module_version=MODULE_VERSION,
        tables={
            "stats_summary": summary,
            "stats_descriptives": descriptives,
            "stats_tests": tests,
            "stats_correlations": correlations,
        },
        citations=("biopython",),
    )
    return _with_sentences(result)


def _number(value: float | None) -> str:
    return "" if value is None else f"{value:.4f}"


def _values(observations: Sequence[Observation], name: str, label: str) -> list[float]:
    """Every value the group carries for one measure — blanks left out, never counted as zero."""
    collected = [
        observation.value(name) for observation in observations if observation.group == label
    ]
    return [value for value in collected if value is not None]


def _descriptive_row(name: str, label: str, values: Sequence[float]) -> dict:
    summary = statistics.describe(values)
    return {
        "observation": f"{name} · {label}",
        "measure": name,
        "group": label,
        "n": summary.n,
        "mean": summary.mean,
        "sd": summary.sd,
        "median": summary.median,
        "q1": summary.q1,
        "q3": summary.q3,
        "minimum": summary.minimum,
        "maximum": summary.maximum,
    }


def _correlations(
    observations: Sequence[Observation],
    measures: Sequence[str],
    labels: Sequence[str],
    scope: Scope,
) -> Table:
    """Correlation between the measures, pooled over every observation that carries both."""
    pairs = [(a, b) for index, a in enumerate(measures) for b in measures[index + 1 :]]
    rows = []
    for first, second in pairs:
        both = [
            (observation.value(first), observation.value(second))
            for observation in observations
            if observation.value(first) is not None and observation.value(second) is not None
        ]
        outcome = statistics.spearman([pair[0] for pair in both], [pair[1] for pair in both])
        rows.append(
            {
                "pair": f"{first} · {second}",
                "measure_x": first,
                "measure_y": second,
                "test": outcome.test,
                "n": outcome.n,
                "coefficient": outcome.coefficient,
                "p": outcome.p,
                "status": outcome.status,
                "note": outcome.reason,
                "groups": ", ".join(labels),
            }
        )
    return Table(
        id="stats_correlations",
        title="Correlation between measures",
        columns=(
            Column("pair", "Measure pair"),
            Column("measure_x", "X"),
            Column("measure_y", "Y"),
            Column("test", "Test"),
            Column("n", "n", dtype="int"),
            Column("coefficient", "Coefficient", dtype="float", precision=4),
            Column("p", "p", dtype="float", precision=6),
            Column("status", "Status"),
            Column("note", "Note"),
            Column("groups", "Groups"),
        ),
        rows=tuple(rows),
        scope=scope,
        key="pair",
    )


def _with_sentences(result: Result) -> Result:
    """The interpretation, assembled from rendered cells — never from a number typed here.

    D5-05 is still open with the PI, so the wording stays inside what the table's own columns
    carry: the test is named, n and α are named, and nothing is called adaptive, selected or
    caused. Every value below is ``Table.format`` of the cell it speaks about.
    """
    summary = result.tables["stats_summary"]
    significant = summary.format("n_significant")
    tested = summary.format("n_tested")
    not_computable = summary.format("n_not_computable")
    sentences = [
        f"{tested} locus × measure test(s) carry a p-value; {summary.format('n_loci')} loci were "
        f"measured in {summary.format('n_genomes')} genome(s) across "
        f"{summary.format('n_groups')} group(s) ({summary.format('groups')}).",
        f"With {summary.format('test')} and {summary.format('correction')} over the whole family "
        f"at α = {summary.format('alpha')}, {significant} test(s) fall below α after correction.",
    ]
    if not_computable != "0":
        sentences.append(
            f"{not_computable} test(s) are not computable from this data — each row of "
            f"stats_tests names its own reason."
        )
    sentences.append(
        f"{summary.format('n_collapsed_cds')} CDS row(s) were collapsed: one observation per "
        f"locus per genome, so the IR's second copy is not counted twice."
    )
    if summary.format("n_genomes_without_cds") != "0":
        sentences.append(
            f"{summary.format('n_genomes_without_cds')} of the "
            f"{summary.format('n_genomes')} input(s) contribute no CDS and therefore no "
            f"observation: {summary.format('genomes_without_cds')}."
        )
    return Result(
        module=result.module,
        module_version=result.module_version,
        tables=result.tables,
        notes=tuple(sentences),
        citations=result.citations,
    )
