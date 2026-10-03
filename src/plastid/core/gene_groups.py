"""Gene-name → functional group, for the publication gene table (P2-09).

The grouping is external knowledge (the convention plastome papers use), so it lives in
``data/gene_groups.yaml``, not in code. This module only loads it, checks it, and applies it.

Two rules keep the table honest:

* **first match wins, in file order** — a group claims a gene by an exact name or by a name
  prefix, so tRNA/rRNA keep working when a species brings a different anticodon suffix;
* **nothing is dropped** — :func:`assign` returns every counted gene, and a gene no group
  claims comes back with ``group=None`` so the table can show it as unclassified.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import yaml

from plastid.core.model import GENE_KINDS, Gene, Genome
from plastid.core.resources import data_file

__all__ = ["Assignment", "GeneGroup", "GeneGroups", "assign", "load_gene_groups"]

_GROUP_KEYS = frozenset({"name", "genes", "prefix"})


class GeneGroupError(ValueError):
    """The gene-group file is malformed, or a group can never match."""


@dataclass(frozen=True, slots=True)
class GeneGroup:
    """One row of the publication table: a category, a group, and how it claims genes."""

    category: str
    name: str
    genes: tuple[str, ...] = ()
    prefix: str | None = None

    def claims(self, gene_name: str) -> bool:
        if gene_name in self.genes:
            return True
        return self.prefix is not None and gene_name.startswith(self.prefix)


@dataclass(frozen=True, slots=True)
class GeneGroups:
    """The ordered group list. Order is meaning: the first claim wins."""

    groups: tuple[GeneGroup, ...]

    def claim(self, gene_name: str) -> GeneGroup | None:
        return next((group for group in self.groups if group.claims(gene_name)), None)


@dataclass(frozen=True, slots=True)
class Assignment:
    """One counted gene and the group it landed in (``None`` = unclassified)."""

    gene: Gene
    group: GeneGroup | None


def load_gene_groups(path: str | Path | None = None) -> GeneGroups:
    """Load and validate ``data/gene_groups.yaml`` (or an explicit file)."""
    source = Path(path) if path is not None else data_file("gene_groups.yaml")
    payload = yaml.safe_load(source.read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or not isinstance(payload.get("categories"), list):
        raise GeneGroupError(f"{source}: expected a 'categories' list")

    groups: list[GeneGroup] = []
    for category in payload["categories"]:
        if not isinstance(category, dict) or not isinstance(category.get("groups"), list):
            raise GeneGroupError(f"{source}: every category needs a 'groups' list")
        category_name = str(category.get("name", ""))
        if not category_name:
            raise GeneGroupError(f"{source}: a category has no name")
        for entry in category["groups"]:
            groups.append(_group_from(entry, category_name, source))
    if not groups:
        raise GeneGroupError(f"{source}: no groups defined")
    _reject_unreachable(groups, source)
    return GeneGroups(groups=tuple(groups))


def _group_from(entry: object, category_name: str, source: Path) -> GeneGroup:
    if not isinstance(entry, dict):
        raise GeneGroupError(f"{source}: a group under {category_name!r} is not a mapping")
    unknown = set(entry) - _GROUP_KEYS
    if unknown:
        raise GeneGroupError(
            f"{source}: group under {category_name!r} has unknown key(s) {sorted(unknown)} — "
            f"allowed: {sorted(_GROUP_KEYS)}"
        )
    name = str(entry.get("name", ""))
    if not name:
        raise GeneGroupError(f"{source}: a group under {category_name!r} has no name")

    genes = entry.get("genes")
    prefix = entry.get("prefix")
    if (genes is None) == (prefix is None):
        raise GeneGroupError(f"{source}: group {name!r} needs exactly one of 'genes' or 'prefix'")
    if genes is not None:
        if not isinstance(genes, list) or not all(isinstance(item, str) for item in genes):
            raise GeneGroupError(f"{source}: group {name!r}: 'genes' must be a list of names")
        return GeneGroup(category=category_name, name=name, genes=tuple(genes))
    if not isinstance(prefix, str) or not prefix:
        raise GeneGroupError(f"{source}: group {name!r}: 'prefix' must be a non-empty string")
    return GeneGroup(category=category_name, name=name, prefix=prefix)


def _reject_unreachable(groups: list[GeneGroup], source: Path) -> None:
    """An explicit name an earlier prefix rule already claims can never be reached.

    Left unchecked this is silent: the gene lands in the prefix group (usually the right
    one), so nothing looks wrong until someone adds a name expecting it to appear
    elsewhere. Fail at load instead.
    """
    for index, group in enumerate(groups):
        earlier = groups[:index]
        for gene_name in group.genes:
            shadow = next((other for other in earlier if other.claims(gene_name)), None)
            if shadow is not None:
                raise GeneGroupError(
                    f"{source}: {gene_name!r} is listed under {group.name!r} but "
                    f"{shadow.name!r} already claims it — the list entry can never match"
                )


def assign(genome: Genome, groups: GeneGroups) -> tuple[Assignment, ...]:
    """Every counted gene (both IR copies) with its group; unclassified comes back as None."""
    return tuple(
        Assignment(gene=gene, group=groups.claim(gene.name))
        for gene in genome.genes
        if gene.kind in GENE_KINDS
    )
