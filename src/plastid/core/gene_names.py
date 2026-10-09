"""The one gene-name normalizer (M3.2, P3-02) — nothing else in the codebase renames a gene.

Three different things look like "the same gene under another name", and they are not the same
thing, so they are handled and reported separately:

* **syntax** — one name written differently: ``trnL(UAA)`` → ``trnL-UAA``, ``rrn16S`` → ``rrn16``,
  case and whitespace. This needs no external knowledge.
* **semantics** — a different name for the same locus: ``pafI`` *is* ``ycf3``. That mapping is
  external knowledge and lives in ``data/gene_aliases.yaml`` with the evidence it rests on; the
  canonical direction was D3-02 (PROJECT.md §8.10) and is **superseded 07/10/2026 by D5-09** — it
  will follow the project's naming reference once the naming step is built, and the alias table
  keeps D3-02's direction until then.
* **sequence** — a bare tRNA name (``trnK``, no anticodon). Measured 2026-10-03: reading the
  anticodon at a fixed offset matches only 106 of 293 named tRNAs in this dataset, so a bare name
  is resolved by 100 % sequence identity against named copies, and is *reported unresolved* when
  that fails. It is never guessed.

The name a record arrived with is never rewritten: :class:`NormalizedName` carries it beside the
canonical form, which is what lets an output print both — ``pafI (ycf3)`` — and lets a reader
trace any number back to the record it came from.
"""

from __future__ import annotations

import json
import re
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import yaml
from Bio import Align

from plastid.core.ir import IRResult
from plastid.core.model import Gene, Genome, coding_sequence
from plastid.core.normalize import normalize, regions_from_ir
from plastid.core.resources import data_file

__all__ = [
    "COVERAGE_FLOOR",
    "IDENTITY_FLOOR",
    "MARGIN_FLOOR_BITS",
    "OUTCOME_CONFIRMED",
    "OUTCOME_CONFLICT",
    "OUTCOME_REPORTED",
    "OUTCOME_RENAMED",
    "STATUS_MEASURED",
    "STATUS_NOT_MEASURED",
    "STATUS_OUTSIDE",
    "TEST_MARGIN",
    "TEST_NAME",
    "TEST_NEIGHBOURS",
    "TEST_SEQUENCE",
    "GeneDecision",
    "Hit",
    "Neighbours",
    "decide_gene",
    "duplicate_names",
    "neighbours_of",
    "nucleotide_identity",
    "oriented_aliases",
    "protein_fasta_rows",
    "vocabulary_of",
    "IDENTICAL_RULE",
    "SEQUENCE_RULE",
    "SEMANTIC_RULE",
    "SYNTAX_RULE",
    "RECORDED_REFERENCE_FILE",
    "AliasTable",
    "GeneAliasError",
    "NamingConflictError",
    "RecordedReference",
    "load_recorded_reference",
    "repo_root",
    "GenomeNames",
    "Locus",
    "NormalizedName",
    "UnresolvedName",
    "load_gene_aliases",
    "aliases_for",
    "cds_by_canonical_name",
    "load_mito_aliases",
    "normalize_genome",
    "normalize_name",
    "references_from",
]

# How a name reached its canonical form — reported per name, because "why is this name
# different" is the question a reader of the table actually has.
IDENTICAL_RULE = "identical"
SYNTAX_RULE = "syntax"
SEMANTIC_RULE = "semantic"
SEQUENCE_RULE = "sequence"

_ALIAS_KEYS = frozenset({"canonical", "synonyms", "evidence"})

# trnL(UAA) → trnL-UAA. Measured 2026-10-03: 30 of the dataset's names use the bracket form,
# all of them in one record (NC_058328.1).
_BRACKET_TRNA = re.compile(r"^(trn\w+)\(([A-Za-z]{3})\)$")
# rrn16S → rrn16. The mitochondrial rrnS/rrnL dialect is deliberately NOT touched: there the
# trailing letter is the name, not a size class.
_RRNA_SIZE_SUFFIX = re.compile(r"^(rrn\d+(?:\.\d+)?)S$")
# "16S ribosomal RNA" → rrn16 — the same rule as the line above, for the record that carries no
# /gene at all and reaches us through the reader's /product fallback. Measured 08/10/2026 on the
# pinned plastid reference set while P3B-07 built its vocabulary: NC_042827.1 names all four rRNA
# this way ("16S ribosomal RNA", {"locus_tag": …, "product": …}, no /gene), and the sentence was
# reaching the caller as a gene name — four loci no run could ever name. The size class is in the
# string, so this needs no external knowledge; the mito dialect is again untouched, because a
# mitogenome's rRNA are rrnS/rrnL and no measured record writes them this way.
_RRNA_PRODUCT = re.compile(r"^(\d+(?:\.\d+)?)S\s+(?:ribosomal\s+)?r?RNA$", re.IGNORECASE)
# A tRNA name with no anticodon at all (measured: trnK in PX923121.1, trnI in NC_069583.1).
_BARE_TRNA = re.compile(r"^trn[A-Za-z]+$")
# A tRNA name that does carry one — the shape a bare name can be resolved to.
_NAMED_TRNA = re.compile(r"^trn\w+-[A-Z]{3}$")
# trnX-YYY with a single amino-acid letter: the only form whose letter case is a slip rather
# than part of the name (trnfM-CAU keeps its f — the next character is M, not a hyphen).
_SINGLE_LETTER_TRNA = re.compile(r"^trn([A-Za-z])-([A-Za-z]{3})$")


class GeneAliasError(ValueError):
    """The alias file is malformed, or an entry cannot be trusted."""


@dataclass(frozen=True, slots=True)
class Locus:
    """One locus and every name it is known by. ``canonical`` leads every printed form."""

    canonical: str
    synonyms: tuple[str, ...]
    evidence: str

    @property
    def members(self) -> tuple[str, ...]:
        return (self.canonical, *self.synonyms)


@dataclass(frozen=True, slots=True)
class AliasTable:
    """The semantic aliases, in file order. A name not here is left alone."""

    loci: tuple[Locus, ...] = ()

    def locus(self, name: str) -> Locus | None:
        return next((locus for locus in self.loci if name in locus.members), None)


@dataclass(frozen=True, slots=True)
class NormalizedName:
    """One gene name, normalized. ``name`` is what the record said — never rewritten."""

    name: str
    canonical: str
    synonyms: tuple[str, ...] = ()
    rule: str = IDENTICAL_RULE

    @property
    def label(self) -> str:
        """The printed form: ``pafI (ycf3)``, or the name alone when it has no synonym."""
        if not self.synonyms:
            return self.canonical
        return f"{self.canonical} ({', '.join(self.synonyms)})"

    @property
    def changed(self) -> bool:
        """True when the canonical form differs from the name the record arrived with."""
        return self.canonical != self.name


@dataclass(frozen=True, slots=True)
class UnresolvedName:
    """A name that could not be resolved honestly — reported, never guessed."""

    gene: Gene
    name: str
    reason: str


@dataclass(frozen=True, slots=True)
class GenomeNames:
    """Every gene of one genome with its normalized name, plus the names left unresolved."""

    names: Mapping[Gene, NormalizedName]
    unresolved: tuple[UnresolvedName, ...] = ()

    def normalized(self, gene: Gene) -> NormalizedName:
        return self.names[gene]

    def canonical(self, gene: Gene) -> str:
        return self.names[gene].canonical

    def label(self, gene: Gene) -> str:
        return self.names[gene].label

    def distinct_canonical(self) -> tuple[str, ...]:
        """The distinct canonical names, in first-seen order."""
        seen: dict[str, None] = {}
        for normalized in self.names.values():
            seen.setdefault(normalized.canonical, None)
        return tuple(seen)


def load_gene_aliases(path: str | Path | None = None) -> AliasTable:
    """Load and validate ``data/gene_aliases.yaml`` (or an explicit file)."""
    source = Path(path) if path is not None else data_file("gene_aliases.yaml")
    return _load_aliases(source)


def _load_aliases(source: Path) -> AliasTable:
    payload = yaml.safe_load(source.read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or not isinstance(payload.get("loci"), list):
        raise GeneAliasError(f"{source}: expected a 'loci' list")

    loci = tuple(_locus_from(entry, source) for entry in payload["loci"])
    _reject_shared_names(loci, source)
    return AliasTable(loci=loci)


def load_mito_aliases(path: str | Path | None = None) -> AliasTable:
    """Load ``data/mito_gene_names.yaml`` — the plant-mitochondrial product-text vocabulary.

    A separate file on purpose: measured 2026-10-03, the same product string is a plastid CDS
    product *and* a mitochondrial one ("ribosomal protein S3" in both), so a shared table could
    rename a plastid gene to a mitochondrial symbol. Only the mitochondrial path loads this one.
    """
    source = Path(path) if path is not None else data_file("mito_gene_names.yaml")
    return _load_aliases(source)


def _locus_from(entry: object, source: Path) -> Locus:
    if not isinstance(entry, dict):
        raise GeneAliasError(f"{source}: a locus entry is not a mapping")
    unknown = set(entry) - _ALIAS_KEYS
    if unknown:
        raise GeneAliasError(
            f"{source}: locus entry has unknown key(s) {sorted(unknown)} — "
            f"allowed: {sorted(_ALIAS_KEYS)}"
        )
    canonical = str(entry.get("canonical", ""))
    if not canonical:
        raise GeneAliasError(f"{source}: a locus entry has no canonical name")

    synonyms = entry.get("synonyms")
    if (
        not isinstance(synonyms, list)
        or not synonyms
        or not all(isinstance(item, str) and item for item in synonyms)
    ):
        raise GeneAliasError(
            f"{source}: locus {canonical!r} needs a non-empty 'synonyms' list of names"
        )
    if canonical in synonyms:
        raise GeneAliasError(
            f"{source}: locus {canonical!r} lists its own canonical name as a synonym"
        )

    evidence = str(entry.get("evidence", "")).strip()
    if not evidence:
        # The file's whole value is that it carries no guesses: an entry without evidence is
        # exactly the thing this check exists to keep out.
        raise GeneAliasError(
            f"{source}: locus {canonical!r} carries no evidence — an alias without a measured "
            f"reason is a guess, and guesses do not belong in this file"
        )
    return Locus(canonical=canonical, synonyms=tuple(synonyms), evidence=evidence)


def _reject_shared_names(loci: tuple[Locus, ...], source: Path) -> None:
    """One name, one locus. A name in two entries would resolve by file order — silently."""
    owner: dict[str, str] = {}
    for locus in loci:
        for name in locus.members:
            if name in owner:
                raise GeneAliasError(
                    f"{source}: {name!r} is listed under both {owner[name]!r} and "
                    f"{locus.canonical!r} — a name belongs to exactly one locus"
                )
            owner[name] = locus.canonical


def normalize_name(name: str, aliases: AliasTable | None = None) -> NormalizedName:
    """Normalize one gene name by syntax and by the semantic alias table.

    A bare tRNA name comes back unchanged here — resolving it needs a sequence, which only
    :func:`normalize_genome` has. ``rule`` says which way the canonical form was reached, so a
    caller can tell "written differently" from "named differently".
    """
    cleaned = " ".join(name.split())
    canonical, syntax_changed = _normalize_syntax(cleaned)

    locus = (aliases or AliasTable()).locus(canonical)
    if locus is not None:
        # The label lists the locus's other names, not "the name that is not the arriving one":
        # a record that already uses the canonical name must print the same label as one that
        # does not, or a reader could tell the dialects apart — which is the whole point.
        return NormalizedName(
            name=name,
            canonical=locus.canonical,
            synonyms=locus.synonyms,
            rule=SEMANTIC_RULE,
        )
    return NormalizedName(
        name=name,
        canonical=canonical,
        rule=SYNTAX_RULE if syntax_changed else IDENTICAL_RULE,
    )


def _normalize_syntax(name: str) -> tuple[str, bool]:
    """``(canonical, changed)`` for the syntax-only rules."""
    canonical = name

    bracket = _BRACKET_TRNA.match(canonical)
    if bracket is not None:
        canonical = f"{bracket.group(1)}-{bracket.group(2).upper()}"

    sized = _RRNA_SIZE_SUFFIX.match(canonical)
    if sized is not None:
        canonical = sized.group(1)

    spelled = _RRNA_PRODUCT.match(canonical)
    if spelled is not None:
        canonical = f"rrn{spelled.group(1)}"

    # The trn/rrn prefix is lowercase in every measured record; only that prefix is touched, so
    # trnfM-CAU keeps its lowercase f (which is part of the name, not a case slip).
    for prefix in ("trn", "rrn"):
        if canonical[:3].lower() == prefix and canonical[:3] != prefix:
            canonical = prefix + canonical[3:]
            break

    # The amino-acid letter of trnX-YYY: `trnk-UUU` and `trnK-UUU` are one name written two ways,
    # and a comparison that keeps them apart loses the locus in silence (measured 2026-10-04 on
    # the golden pair: NC_024936.1 spells trnK-UUU and trnI-GAU with lowercase letters —
    # `trnk-UUU`, `trnl-GUA` — and the two records' sequences at those loci are identical).
    # Only the single-letter-before-the-hyphen form is touched: trnfM-CAU's `f` is followed by
    # `M`, not by `-`, and it is part of the name.
    single = _SINGLE_LETTER_TRNA.match(canonical)
    if single is not None and single.group(1).islower():
        canonical = f"trn{single.group(1).upper()}-{single.group(2)}"
    return canonical, canonical != name


def normalize_genome(
    genome: Genome,
    aliases: AliasTable | None = None,
    *,
    references: Sequence[tuple[str, str]] = (),
) -> GenomeNames:
    """Normalize every gene of one genome; resolve bare tRNA names against ``references``.

    ``references`` are ``(name, sequence)`` pairs from *other* records — the same named tRNAs
    this record's bare names may stand for. The genome's own named tRNAs are added to them
    automatically. A bare name that matches no reference, or matches more than one, is reported
    unresolved rather than guessed.
    """
    table = aliases or AliasTable()
    pool = (*references_from([genome], table), *references)

    # A bare name is resolved from the *tRNA* feature, which is the one carrying the sequence.
    # The `gene` feature of the same name spans the whole locus including its intron, so
    # resolving per feature would compare an intron against a tRNA and report the name
    # unresolved — measured 2026-10-03 on PX923121.1, whose trnK has both features.
    resolved: dict[str, NormalizedName] = {}
    unresolved_by_name: dict[str, UnresolvedName] = {}
    for gene in genome.genes:
        if gene.kind != "tRNA" or gene.name in resolved or gene.name in unresolved_by_name:
            continue
        normalized = normalize_name(gene.name, table)
        if not _BARE_TRNA.match(normalized.canonical):
            continue
        candidate, reason = _resolve_bare_trna(genome, gene, normalized, pool)
        if reason is None:
            resolved[gene.name] = candidate
        else:
            unresolved_by_name[gene.name] = UnresolvedName(gene=gene, name=gene.name, reason=reason)

    names: dict[Gene, NormalizedName] = {}
    unresolved: list[UnresolvedName] = []
    reported: set[str] = set()
    for gene in genome.genes:
        normalized = resolved.get(gene.name) or normalize_name(gene.name, table)
        if gene.name in unresolved_by_name and gene.name not in reported:
            # once per name: two IR copies of an unresolved tRNA are one finding, not two
            reported.add(gene.name)
            unresolved.append(unresolved_by_name[gene.name])
        names[gene] = normalized
    return GenomeNames(names=names, unresolved=tuple(unresolved))


def _resolve_bare_trna(
    genome: Genome,
    gene: Gene,
    normalized: NormalizedName,
    pool: Sequence[tuple[str, str]],
) -> tuple[NormalizedName, str | None]:
    """Resolve ``trnK`` to ``trnK-UUU`` by 100 % identity, or say why it could not."""
    sequence = coding_sequence(genome, gene).upper()
    prefix = f"{normalized.canonical}-"
    matches = sorted(
        {
            name
            for name, reference in pool
            if name.startswith(prefix) and reference.upper() == sequence
        }
    )
    if len(matches) == 1:
        return (
            NormalizedName(
                name=normalized.name,
                canonical=matches[0],
                rule=SEQUENCE_RULE,
            ),
            None,
        )
    if not matches:
        reason = (
            "no reference sequence matched at 100 % identity — a bare tRNA name cannot be "
            "resolved from the record alone; run with records that name their tRNAs"
        )
    else:
        reason = f"ambiguous — {len(matches)} reference names match: {', '.join(matches)}"
    return normalized, reason


def aliases_for(genome: Genome) -> AliasTable:
    """The name table this record needs: the mitochondrial vocabulary for a mitogenome.

    Measured 2026-10-03: the same product string names a plastid gene and a mitochondrial one
    ("ribosomal protein S3" in both), so there is one vocabulary per organelle and the record's
    own ``/organelle`` qualifier picks it. A record that declares no organelle gets the plastid
    table — and a mito record whose features already carry symbols (this project's own annotator
    writes them) is unaffected either way.
    """
    organelle = str(genome.source.get("organelle", "")).lower()
    return load_mito_aliases() if "mitochond" in organelle else load_gene_aliases()


def cds_by_canonical_name(genome: Genome, aliases: AliasTable | None = None) -> dict[str, Gene]:
    """The genome's CDS by canonical name — one feature per name, so both IR copies are one gene.

    Shared by the comparison modules that work gene by gene (substitutions' CDS scope, Ka/Ks,
    popgen): the gene set of a comparison is *names*, not features, and a name that appears twice
    (the inverted repeat) must be counted once — otherwise every ratio is weighted by whichever
    genes happen to sit in an IR. Without an explicit table the record's own organelle decides
    which vocabulary applies (:func:`aliases_for`): measured 2026-10-04, the mito half of Phase 4
    found **no** shared gene between two mitogenomes because one record names its CDS by product
    sentence — the vocabulary that maps those sentences is the one the mito path already had.
    """
    names = normalize_genome(genome, aliases or aliases_for(genome))
    by_name: dict[str, Gene] = {}
    for gene in genome.genes:
        if gene.kind != "CDS" or gene.pseudo:
            continue
        by_name.setdefault(names.normalized(gene).canonical, gene)
    return by_name


def references_from(
    genomes: Iterable[Genome], aliases: AliasTable | None = None
) -> tuple[tuple[str, str], ...]:
    """``(canonical name, sequence)`` for every named tRNA in these records.

    The pool a bare name is resolved against. Only names that already carry an anticodon are
    useful — a bare name cannot vouch for another bare name.
    """
    table = aliases or AliasTable()
    pairs: list[tuple[str, str]] = []
    for genome in genomes:
        for gene in genome.genes:
            if gene.kind != "tRNA":
                continue
            canonical = normalize_name(gene.name, table).canonical
            if _NAMED_TRNA.match(canonical):
                pairs.append((canonical, coding_sequence(genome, gene)))
    return tuple(pairs)


# ---------------------------------------------------------------------------------------------
# Reference-driven naming (D5-09, P3B-04): the name proposes, the sequence disposes.
#
# Four tests, in order, stopping at the first that decides (D5-09 §4): the name (syntax, then the
# alias table oriented by the naming reference), the sequence ranking, the margin between the best
# and the runner-up, and the neighbours within one region. What is decided here is a *pairing*
# between a gene of a record and a gene of the project's naming reference — the vocabulary the
# whole project prints (D5-09 §3, vocabulary A).
#
# The thresholds are measured, not chosen (PROJECT.md §8.41): on 2,838 cross-record CDS queries
# the smallest true margin was 19.2 bits (median 298.0) and none fell under 5; with coverage ≥ 0.7
# the best a genuinely different locus reached was 44.93 % identity while the worst a true ortholog
# reached was 67.75 %. The mito path's coverage floor of 0.5 does **not** hold for plastid CDS —
# `psbJ` is 36 aa and a 21-aa alignment of it reaches 52.38 % against `matK`.
# ---------------------------------------------------------------------------------------------

# Which test decided a pairing — the column a later reader needs to tell "the name and the
# sequence agree" from "only the gene order could say".
TEST_NAME = "name"
TEST_SEQUENCE = "sequence"
TEST_MARGIN = "margin"
TEST_NEIGHBOURS = "neighbours"

# The status column (D5-09 §5). Three values that must stay distinguishable, because a blank and
# a zero look alike in every aggregate: the check ran, it could not run, or it ran and found
# nothing in the reference. See :class:`GeneDecision` for the exact definitions used here.
STATUS_MEASURED = "measured"
STATUS_NOT_MEASURED = "not measured"
STATUS_OUTSIDE = "outside vocabulary"

# What the step did about the pairing.
OUTCOME_CONFIRMED = "confirmed"
OUTCOME_RENAMED = "renamed"
OUTCOME_CONFLICT = "conflict"
# Nothing changes in the output, and nothing confirmed the pairing either: the gene is written
# as it arrived and a person is told what could not be settled (D5-09 §5's "reported" tier).
OUTCOME_REPORTED = "reported"

# The measured floors (§8.41). Named here so a module never carries a bare number.
IDENTITY_FLOOR = 50.0
COVERAGE_FLOOR = 0.7
MARGIN_FLOOR_BITS = 10.0

# The pairwise-alignment parameters §8.41 measured the tRNA numbers with. They are part of the
# number: with cheaper gaps the same trnG pair reads 78.9 % instead of 73.2 %, which is why these
# are constants here rather than literals at the call site.
TRNA_MATCH_SCORE = 1.0
TRNA_GAP_OPEN = -2.0
TRNA_GAP_EXTEND = -0.5


@dataclass(frozen=True, slots=True)
class Hit:
    """One candidate in the naming reference, as the sequence ranking put it.

    ``bits`` is BLAST's score and exists only for protein searches; a nucleotide comparison
    carries identity without a score, because nothing measured a threshold for one (§8.41: for
    tRNA no identity floor works at all, so the decision has to come from the name and the
    neighbours — the comparison is recorded as evidence, never as a verdict).
    """

    locus: str  # the naming reference's own spelling
    identity: float  # % identity over the aligned length
    coverage: float
    bits: float | None = None


@dataclass(frozen=True, slots=True)
class GeneDecision:
    """What became of one gene of one record: its pairing, the test that decided it, and why.

    Three statuses, and the difference between two of them is the whole point of the column:

    * ``measured`` — the counterpart in A was found and checked, by name and sequence (a CDS), by
      name and neighbours (a tRNA whose ranking cannot decide), or by sequence alone (a name the
      reference does not use, matched uniquely).
    * ``not measured`` — the check **could not run**, or ran and did not settle the pairing: no
      protein on one side, no aligner, a margin too narrow to trust, or flanking genes that
      disagree. It never means "none found".
    * ``outside vocabulary`` — the check ran and found nothing in the reference: this gene has no
      counterpart in A. B's spelling is printed as it arrived.

    ``outcome`` says what the step does about it: ``confirmed`` (nothing changes), ``renamed``
    (the output spells it A's way), ``conflict`` (stop tier — only ever a CDS, because only a
    protein ranking is evidence of a *different locus*), ``reported`` (nothing changes and a
    person should look).
    """

    gene: Gene
    status: str
    outcome: str
    canonical: str | None  # the name in the naming reference, when one was found
    rule: str | None  # which test decided (None when nothing did)
    identity: float | None = None
    coverage: float | None = None
    margin_bits: float | None = None
    best_locus: str | None = None  # what the ranking preferred, whatever was decided
    reason: str = ""

    @property
    def renamed(self) -> bool:
        """True when the name changes in the output — the recipe is a different record."""
        return self.outcome == OUTCOME_RENAMED

    @property
    def conflict(self) -> bool:
        return self.outcome == OUTCOME_CONFLICT


@dataclass(frozen=True, slots=True)
class Neighbours:
    """The genes flanking one gene in the normalized frame — the evidence test 4 reads."""

    region: str
    left: str | None
    right: str | None

    def agrees_with(self, other: Neighbours) -> bool:
        """True when at least one flank matches and no flank that both records name disagrees.

        Deliberately not "both flanks match": a gene at a region boundary has a different
        neighbour on one side in a different record, and refusing to confirm because of that
        would turn real biology into a conflict. A single agreeing flank is a signal; two
        disagreeing flanks are a contradiction.
        """
        pairs = [
            (mine, theirs)
            for mine, theirs in ((self.left, other.left), (self.right, other.right))
            if mine is not None and theirs is not None
        ]
        if not pairs:
            return False
        agreed = sum(1 for mine, theirs in pairs if mine == theirs)
        return agreed >= 1 and agreed == len(pairs)


def oriented_aliases(aliases: AliasTable, reference: Genome) -> AliasTable:
    """The alias table with each locus's canonical form set to the **reference's** spelling.

    D5-09 §8: the canonical direction stops being the group's dialect and becomes the naming
    reference's. The file says which names are one locus (external knowledge, unaffected by who
    the reference is); the reference says which of them the project prints. A locus the reference
    does not use keeps the file's direction — a gene the reference has nothing to say about is not
    a place to invent a preference.

    A reference that uses **two** members of one locus (two spellings of one gene in one record)
    is also left as the file has it; that is a duplicate-name state the caller reports, and
    picking one of the two here would hide it.
    """
    used = {normalize_name(gene.name).canonical for gene in reference.genes}
    oriented: list[Locus] = []
    for locus in aliases.loci:
        present = [member for member in locus.members if normalize_name(member).canonical in used]
        if len(present) == 1:
            canonical = present[0]
            oriented.append(
                Locus(
                    canonical=canonical,
                    synonyms=tuple(m for m in locus.members if m != canonical),
                    evidence=locus.evidence,
                )
            )
        else:
            oriented.append(locus)
    return AliasTable(loci=tuple(oriented))


def vocabulary_of(reference: Genome) -> dict[str, Gene]:
    """The naming reference's genes by the name it spells them with — vocabulary A itself.

    Two features under one name (the golden annotates two different trnG genes as ``trnG-UCC``)
    are **not** merged: the first is kept and the duplicate is left for the caller to report, so
    a pairing can never silently point at "one of the two".
    """
    by_name: dict[str, Gene] = {}
    for gene in reference.genes:
        canonical = normalize_name(gene.name).canonical
        by_name.setdefault(canonical, gene)
    return by_name


def duplicate_names(genome: Genome, ir: IRResult | None = None) -> tuple[str, ...]:
    """Names this record uses for more than one **different** gene, in first-seen order.

    Two features of one name are one gene when they are the **inverted repeat's two copies** —
    every copy inside IRa/IRb, one per side — and that is a fact about the record's structure
    rather than a similarity threshold. Measured 08/10/2026 on *Iris hermona*'s published record:
    its two ``rrn23`` copies differ by a base or two (the gate needs ≥ 99.9 %, not 100 %), so a
    check comparing sequences called one gene two and raised a conflict where nothing was wrong.
    Where the frame cannot be established, exact sequence equality is the fallback.

    The golden's two ``trnG-UCC`` genes are the case this must catch: both sit in the **LSC**,
    they are 73.24 % identical (§8.41), and they are two different genes wearing one name.
    """
    regions = _regions_by_gene(genome, ir) if ir is not None else {}
    by_name: dict[str, list[Gene]] = {}
    for gene in genome.genes:
        if gene.kind in ("CDS", "tRNA", "rRNA"):
            by_name.setdefault(normalize_name(gene.name).canonical, []).append(gene)

    return tuple(
        name
        for name, genes in by_name.items()
        if len(genes) > 1 and not _one_gene_two_copies(genes, genome, regions)
    )


def _one_gene_two_copies(genes: Sequence[Gene], genome: Genome, regions: Mapping[int, str]) -> bool:
    """True when a group of same-named features is the inverted repeat's two copies."""
    where = [regions.get(gene) for gene in genes]
    if regions and all(place in ("IRa", "IRb") for place in where):
        # one copy per side, and nothing inside a side that differs from its own copy
        sequences: dict[str, set[str]] = {}
        for gene, place in zip(genes, where, strict=True):
            sequences.setdefault(str(place), set()).add(coding_sequence(genome, gene).upper())
        return all(len(seen) == 1 for seen in sequences.values())
    return len({coding_sequence(genome, gene).upper() for gene in genes}) == 1


def _regions_by_gene(genome: Genome, ir: IRResult) -> dict[Gene, str]:
    """Which region each of *this* genome's features falls in, read in the normalized frame.

    The frame is a pure rotation, so each feature's own midpoint is mapped through the rotation
    offset rather than looked up among the rotated copies — the rotated genome holds *new* Gene
    objects, and keying by identity there would silently place every gene in no region at all
    (measured 08/10/2026, when it made one published record's two ``rrn23`` copies look like two
    different genes).
    """
    frame = normalize(genome, ir)
    regions = {region.name: region for region in regions_from_ir(frame.genome, frame.ir)}
    length = len(genome)
    out: dict[Gene, str] = {}
    for gene in genome.genes:
        start_1b, end_1b = gene.span
        midpoint = ((start_1b - 1 + end_1b) // 2 - frame.rotation_offset) % length
        for name, region in regions.items():
            if region.start <= midpoint < region.end:
                out[gene] = name
                break
    return out


def nucleotide_identity(left: str, right: str) -> float:
    """Percent identity of the best global alignment — the metric §8.41 measured tRNA with.

    BLAST needs seeds a tRNA does not have, so the comparison is a pairwise alignment, and its
    parameters are part of the number (see the constants above). Sequences too short to align
    honestly return 0.0 rather than a rounded-up guess.
    """
    if not left or not right:
        return 0.0
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = TRNA_MATCH_SCORE
    aligner.mismatch_score = 0.0
    aligner.open_gap_score = TRNA_GAP_OPEN
    aligner.extend_gap_score = TRNA_GAP_EXTEND
    alignment = aligner.align(left.upper(), right.upper())[0]
    matches = sum(1 for x, y in zip(alignment[0], alignment[1], strict=True) if x == y and x != "-")
    return 100.0 * matches / max(len(left), len(right))


def protein_fasta_rows(genome: Genome, *, vocabulary: str) -> tuple[tuple[str, str], ...]:
    """``(subject key, protein)`` for every CDS of a record that carries one.

    The key is the three-field contract BLAST and the mito reference already share —
    ``accession|gene|vocabulary`` — with the record's own spelling and the vocabulary it is being
    read as. The record's own ``/translation`` is used and never a self-translation (D5-08 §8.21);
    a CDS without one is left out here and reported by the caller, not translated by us.
    """
    accession = genome.accession or genome.name
    return tuple(
        (f"{accession}|{gene.name}|{vocabulary}", gene.translation)
        for gene in genome.genes
        if gene.kind == "CDS" and gene.translation
    )


def decide_gene(
    gene: Gene,
    *,
    molecules: str,
    vocabulary: Mapping[str, Gene],
    hits: Sequence[Hit] = (),
    aliases: AliasTable,
    neighbours: Neighbours | None = None,
    reference_neighbours: Neighbours | None = None,
) -> GeneDecision:
    """Run D5-09 §4's four tests over one gene and say what they decided.

    ``molecules`` is ``"protein"`` for a CDS and ``"nucleotide"`` for a tRNA or rRNA, and it
    changes which tests may decide — not as a shortcut, but because the measurement says so
    (§8.41): for a tRNA the true and the wrong loci overlap in identity (53.76–100 % against up
    to 73.97 %), so a ranking that prefers another serine tRNA is not evidence of anything, and
    only the name and the neighbours can decide.

    The caller does the searches and computes the flanks; this function is the rule, and it is
    pure so the rule can be tested without an aligner.
    """
    hypothesis = normalize_name(gene.name, aliases)
    named = hypothesis.canonical
    in_vocabulary = named in vocabulary

    best = hits[0] if hits else None
    runner_up = hits[1] if len(hits) > 1 else None
    margin = (
        best.bits - runner_up.bits
        if best and runner_up and None not in (best.bits, runner_up.bits)
        else None
    )

    # What the step does about it: a pairing whose *spelling* differs from the record's is a
    # rename even when nothing else was in doubt — the output changes, so the row has to say so.
    settled = OUTCOME_RENAMED if hypothesis.changed else OUTCOME_CONFIRMED
    # A pairing nothing could check: the spelling may still change (the alias table is recorded
    # knowledge, not a guess), but a name that needed no change is *reported* rather than
    # confirmed — nobody confirmed it.
    unchecked = OUTCOME_RENAMED if hypothesis.changed else OUTCOME_REPORTED

    def decision(**kwargs: object) -> GeneDecision:
        base: dict[str, object] = {
            "gene": gene,
            "canonical": None,
            "rule": None,
            "identity": best.identity if best else None,
            "coverage": best.coverage if best else None,
            "margin_bits": margin,
            "best_locus": best.locus if best else None,
        }
        return GeneDecision(**{**base, **kwargs})  # type: ignore[arg-type]

    # --- test 1: the name, and whether the reference uses it -----------------------------------
    if in_vocabulary:
        if molecules == "nucleotide":
            # The ranking cannot decide for a tRNA (§8.41); the name does, and the neighbours
            # are the only second opinion available.
            if neighbours is not None and reference_neighbours is not None:
                if neighbours.agrees_with(reference_neighbours):
                    return decision(
                        status=STATUS_MEASURED,
                        outcome=settled,
                        canonical=named,
                        rule=TEST_NEIGHBOURS,
                        reason=(
                            f"the name is the reference's, and the flanking genes agree "
                            f"({_flanks(neighbours)})"
                        ),
                    )
                return decision(
                    status=STATUS_NOT_MEASURED,
                    outcome=unchecked,
                    canonical=named,
                    rule=TEST_NEIGHBOURS,
                    reason=(
                        f"the name is the reference's and nothing changed; the flanking genes do "
                        f"not agree with {named} ({_flanks(neighbours)} against "
                        f"{_flanks(reference_neighbours)}), so the pairing is reported rather "
                        f"than confirmed"
                    ),
                )
            return decision(
                status=STATUS_MEASURED,
                outcome=settled,
                canonical=named,
                rule=TEST_NAME,
                reason=(
                    "the name is the reference's; a nucleotide ranking cannot separate tRNA loci "
                    "(measured §8.41) and no IR frame was available for the neighbours test"
                    if best is not None
                    else "the name is the reference's; no nucleotide comparison was made"
                ),
            )
        # a CDS: the name is the hypothesis, and the sequence has to agree with it
        if best is None:
            return decision(
                status=STATUS_NOT_MEASURED,
                outcome=unchecked,
                canonical=named,
                rule=TEST_NAME,
                reason=(
                    "the name is the reference's, but this CDS carries no /translation and no "
                    "sequence search could check it"
                ),
            )
        if best.locus == named and _passes_floors(best) and _margin_ok(margin):
            rule = TEST_MARGIN if runner_up is not None else TEST_SEQUENCE
            return decision(
                status=STATUS_MEASURED,
                outcome=settled,
                canonical=named,
                rule=rule,
                reason=(
                    f"the name is the reference's and the sequence ranks it first "
                    f"({best.identity:.2f} % over {best.coverage:.2f} of the query"
                    + (
                        f", {margin:.1f} bits clear of {runner_up.locus}"
                        if margin is not None
                        else ""
                    )
                    + ")"
                ),
            )
        if best.locus != named and neighbours is not None and reference_neighbours is not None:
            if neighbours.agrees_with(reference_neighbours):
                return decision(
                    status=STATUS_MEASURED,
                    outcome=settled,
                    canonical=named,
                    rule=TEST_NEIGHBOURS,
                    reason=(
                        f"the ranking prefers {best.locus}, but the flanking genes place this "
                        f"gene with {named} ({_flanks(neighbours)})"
                    ),
                )
            return decision(
                status=STATUS_MEASURED,
                outcome=OUTCOME_CONFLICT,
                canonical=named,
                rule=TEST_NEIGHBOURS,
                reason=(
                    f"the name is {named} and the sequence ranks {best.locus} first "
                    f"({best.identity:.2f} %); the flanking genes agree with neither"
                ),
            )
        return decision(
            status=STATUS_MEASURED,
            outcome=OUTCOME_CONFLICT,
            canonical=named,
            rule=TEST_SEQUENCE if _passes_floors(best) else TEST_MARGIN,
            reason=(
                f"the name is {named}, but the sequence ranks {best.locus} first "
                f"({best.identity:.2f} % over {best.coverage:.2f} of the query)"
                + (" and no IR frame placed the neighbours" if neighbours is None else "")
            ),
        )

    # --- test 1 has no name in the reference: the sequence gets to propose one -----------------
    if molecules == "nucleotide":
        # A tRNA is never renamed on a ranking the measurement says is unreliable. Its name is
        # printed as it arrived, and the absence is reported rather than filled with a guess.
        return decision(
            status=STATUS_OUTSIDE,
            outcome=OUTCOME_REPORTED,
            canonical=None,
            rule=TEST_NAME,
            reason=(
                f"the reference does not use the name {named}"
                + (
                    f"; the closest nucleotide match is {best.locus} "
                    f"({best.identity:.2f} %), which cannot separate tRNA loci (measured §8.41)"
                    if best is not None
                    else ""
                )
            ),
        )
    if best is None:
        return decision(
            status=STATUS_OUTSIDE,
            outcome=OUTCOME_REPORTED,
            canonical=None,
            rule=None,
            reason=f"the reference does not use the name {named}, and no sequence search ran",
        )
    if not _passes_floors(best):
        return decision(
            status=STATUS_OUTSIDE,
            outcome=OUTCOME_REPORTED,
            canonical=None,
            rule=None,
            reason=(
                f"the reference does not use the name {named}; the best match {best.locus} is "
                f"below the floor ({best.identity:.2f} % over {best.coverage:.2f} of the query — "
                f"floor {IDENTITY_FLOOR:.0f} % and coverage {COVERAGE_FLOOR:g})"
            ),
        )
    if not _margin_ok(margin):
        return decision(
            status=STATUS_NOT_MEASURED,
            outcome=unchecked,
            canonical=None,
            rule=TEST_MARGIN,
            reason=(
                f"the reference does not use the name {named}; the ranking prefers {best.locus} "
                f"but the margin is {_margin_text(margin)} — too close to rename on "
                f"(floor {MARGIN_FLOOR_BITS:g} bits)"
            ),
        )
    rule = TEST_MARGIN if runner_up is not None else TEST_SEQUENCE
    return decision(
        status=STATUS_MEASURED,
        outcome=OUTCOME_RENAMED,
        canonical=best.locus,
        rule=rule,
        reason=(
            f"the reference does not use the name {named}; the sequence identifies it as "
            f"{best.locus} ({best.identity:.2f} % over {best.coverage:.2f} of the query"
            + (f", {margin:.1f} bits clear of {runner_up.locus}" if margin is not None else "")
            + ")"
        ),
    )


def _passes_floors(hit: Hit) -> bool:
    return hit.identity >= IDENTITY_FLOOR and hit.coverage >= COVERAGE_FLOOR


def _margin_ok(margin: float | None) -> bool:
    """A single hit has nothing to be confused with; a margin exists only when there is a second."""
    return margin is None or margin >= MARGIN_FLOOR_BITS


def _margin_text(margin: float | None) -> str:
    if margin is None:
        return "not computable"
    return f"{margin:.1f} bits"


def _flanks(neighbours: Neighbours) -> str:
    return f"{neighbours.left or '—'} … {neighbours.right or '—'} in the {neighbours.region}"


def neighbours_of(
    genome: Genome, gene: Gene, ir: IRResult, aliases: AliasTable
) -> Neighbours | None:
    """The genes flanking ``gene`` **within its own region** — test 4's evidence.

    Measured regions, not raw order: a plastome's IR copies the same genes in the same order, so
    a gene near an IR boundary has different neighbours depending on which copy it is, and raw
    file order would compare the wrong ones. ``None`` when the frame cannot be established (no
    IR) — the test then reports itself as not measured rather than guessing an order.
    """
    frame = normalize(genome, ir)
    regions = {region.name: region for region in regions_from_ir(frame.genome, frame.ir)}
    placed: list[tuple[str, str, int]] = []
    for candidate in frame.genome.genes:
        if candidate.kind not in ("CDS", "tRNA", "rRNA") or candidate.pseudo:
            continue
        start_1b, end_1b = candidate.span
        midpoint = (start_1b - 1 + end_1b) // 2
        where = next(
            name for name, region in regions.items() if region.start <= midpoint < region.end
        )
        placed.append((where, normalize_name(candidate.name, aliases).canonical, midpoint))
    placed.sort(key=lambda row: row[2])
    for index, (where, name, _) in enumerate(placed):
        if name != normalize_name(gene.name, aliases).canonical:
            continue
        left = None if index == 0 or placed[index - 1][0] != where else placed[index - 1][1]
        right = (
            None
            if index + 1 >= len(placed) or placed[index + 1][0] != where
            else placed[index + 1][1]
        )
        return Neighbours(region=where, left=left, right=right)
    return None


# ---------------------------------------------------------------------------------------------
# The project's naming reference (D5-09 §7): one record, chosen by the end user, recorded here.
# ---------------------------------------------------------------------------------------------

RECORDED_REFERENCE_FILE = "naming_reference.json"


class NamingConflictError(Exception):
    """A naming conflict stopped the step — every case named with its coordinates (D5-09 §5).

    Raised by the module and mapped to the project's conflict exit code by ``cli.py``. It is an
    exception rather than a return value because a half-normalized record is not a result: the
    two-invocation decision is *stop, then re-run with* ``on_conflict=continue``.
    """


@dataclass(frozen=True, slots=True)
class RecordedReference:
    """The recorded standard: which file, which bytes, who recorded it, and whether it is final.

    ``provisional`` is the field that stops an override from quietly becoming the decision: a
    provisional entry stands in until the end user picks from the audit report (D5-09 §7.1), and
    every report that used it can say so.
    """

    path: str
    accession: str
    sha256: str
    provisional: bool
    recorded: str
    recorded_by: str
    why: str

    def resolve(self, *, root: Path | None = None) -> Path:
        """The reference as a real file, relative to the repository root unless given absolute.

        Missing is an error with the path in it: a recorded standard that is not there would
        otherwise surface as "the run used whatever was first on the command line".
        """
        candidate = Path(self.path)
        if not candidate.is_absolute() and root is not None:
            candidate = root / candidate
        if not candidate.is_file():
            raise FileNotFoundError(
                f"the recorded naming reference is not at {candidate} — the standard this "
                f"project prints is a file that has to be there (D5-09 §7.2)"
            )
        return candidate


def repo_root(start: Path | None = None) -> Path:
    """The checkout's root — the directory holding ``pyproject.toml``, found by walking up.

    Used only to resolve a *recorded* path (the reference lives with the project's pinned test
    data). Walking up rather than counting ``..`` keeps it true when a loader moves.
    """
    here = (start or Path(__file__)).resolve()
    for candidate in (here, *here.parents):
        if (candidate / "pyproject.toml").is_file():
            return candidate
    raise FileNotFoundError(f"no pyproject.toml above {here} — cannot locate the repository root")


def load_recorded_reference(path: str | Path | None = None) -> RecordedReference:
    """Read ``data/naming_reference.json`` — the standard, with the bytes it was recorded at.

    A malformed record raises rather than defaulting: a missing *why* or *sha256* is a record
    nobody can audit later, and defaulting would make the standard whatever the file happened to
    say.
    """
    source = Path(path) if path is not None else data_file(RECORDED_REFERENCE_FILE)
    payload = json.loads(source.read_text(encoding="utf-8"))
    required = ("path", "accession", "sha256", "provisional", "recorded", "recorded_by", "why")
    missing = [key for key in required if key not in payload]
    if missing:
        raise GeneAliasError(f"{source}: the naming reference is missing {', '.join(missing)}")
    if len(str(payload["sha256"])) != 64:
        raise GeneAliasError(f"{source}: sha256 must be 64 hex characters of the record's bytes")
    return RecordedReference(
        path=str(payload["path"]),
        accession=str(payload["accession"]),
        sha256=str(payload["sha256"]),
        provisional=bool(payload["provisional"]),
        recorded=str(payload["recorded"]),
        recorded_by=str(payload["recorded_by"]),
        why=str(payload["why"]),
    )
