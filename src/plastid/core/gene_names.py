"""The one gene-name normalizer (M3.2, P3-02) — nothing else in the codebase renames a gene.

Three different things look like "the same gene under another name", and they are not the same
thing, so they are handled and reported separately:

* **syntax** — one name written differently: ``trnL(UAA)`` → ``trnL-UAA``, ``rrn16S`` → ``rrn16``,
  case and whitespace. This needs no external knowledge.
* **semantics** — a different name for the same locus: ``pafI`` *is* ``ycf3``. That mapping is
  external knowledge and lives in ``data/gene_aliases.yaml`` with the evidence it rests on; the
  canonical direction is D3-02 (PROJECT.md §8.10).
* **sequence** — a bare tRNA name (``trnK``, no anticodon). Measured 2026-10-03: reading the
  anticodon at a fixed offset matches only 106 of 293 named tRNAs in this dataset, so a bare name
  is resolved by 100 % sequence identity against named copies, and is *reported unresolved* when
  that fails. It is never guessed.

The name a record arrived with is never rewritten: :class:`NormalizedName` carries it beside the
canonical form, which is what lets an output print both — ``pafI (ycf3)`` — and lets a reader
trace any number back to the record it came from.
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import yaml

from plastid.core.model import Gene, Genome, coding_sequence
from plastid.core.resources import data_file

__all__ = [
    "IDENTICAL_RULE",
    "SEQUENCE_RULE",
    "SEMANTIC_RULE",
    "SYNTAX_RULE",
    "AliasTable",
    "GeneAliasError",
    "GenomeNames",
    "Locus",
    "NormalizedName",
    "UnresolvedName",
    "load_gene_aliases",
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
# A tRNA name with no anticodon at all (measured: trnK in PX923121.1, trnI in NC_069583.1).
_BARE_TRNA = re.compile(r"^trn[A-Za-z]+$")
# A tRNA name that does carry one — the shape a bare name can be resolved to.
_NAMED_TRNA = re.compile(r"^trn\w+-[A-Z]{3}$")


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

    # The trn/rrn prefix is lowercase in every measured record; only that prefix is touched, so
    # trnfM-CAU keeps its lowercase f (which is part of the name, not a case slip).
    for prefix in ("trn", "rrn"):
        if canonical[:3].lower() == prefix and canonical[:3] != prefix:
            canonical = prefix + canonical[3:]
            break
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
