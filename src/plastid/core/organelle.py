"""Which organelle a record is, and the rule that one run serves exactly one (D5-06, H-09).

The toolkit serves plastids and mitochondria, but **a run belongs to one organelle** (decided
05/10/2026). The declaration is the record's own ``/organelle`` qualifier; a bare sequence —
an assembler's output, which is what the annotate step takes — has no record to ask, so it
declares itself **by the path it enters**: ``plastid annotate`` is the plastid path, ``plastid
mito annotate`` the mitochondrial one.

**Why this is a module and not a line in the CLI.** Until P3B-09 the only signal a module carried
was ``requires_ir``, which answers "does this need the inverted repeat" and *not* "which organelle
does it serve" — measured 05/10/2026 (H-09): running ``mito_genes`` on a plastome was
stopped by nothing at all, because only one of the two mismatch directions had a gate. The two
directions are the same question, so they are answered in one place, and the registry carries each
module's group (:mod:`plastid.modules.registry`) so the answer does not depend on remembering.

**What "declared" does not mean.** A mitogenome cannot be told from a plastome by its bytes: both
are circular DNA with protein-coding genes, and plant mito sequences are full of plastid-derived
inserts. So this module never guesses — it reads a qualifier, and where there is none it says the
run's path is the declaration. That is also why a *mixed* input set is a refusal rather than a
warning: nothing downstream could tell which vocabulary, which translation table or which regions
the numbers belong to.
"""

from __future__ import annotations

from collections.abc import Sequence

from plastid.core.model import Genome

__all__ = [
    "GROUP_MITO",
    "GROUP_PLASTID",
    "GROUP_SHARED",
    "GROUPS",
    "MITO",
    "MITO_PATH",
    "PLASTID",
    "PLASTID_PATH",
    "OrganelleError",
    "check_set",
    "organelle_of",
]

PLASTID = "plastid"
MITO = "mitochondrion"

# The three module groups: shared · plastid-only · mito-only. Every module declares one in the
# registry; these constants are how it says which.
GROUP_SHARED = "shared"
GROUP_PLASTID = "plastid"
GROUP_MITO = "mito"
GROUPS = (GROUP_SHARED, GROUP_PLASTID, GROUP_MITO)

# The command that serves each organelle — named in a refusal, because an error a user cannot act
# on is a dead end (D5-08: ``file — reason``, and the reason has to be usable).
PLASTID_PATH = "plastid annotate"
MITO_PATH = "plastid mito annotate"

# How a record spells its own organelle. Both spellings are in the pinned sets, measured
# 09/10/2026: the Iridaceae set carries 29 `plastid:chloroplast` and 5 `plastid`, the mito set 8
# `mitochondrion`. Matching is on the stem, the way `core/gene_names.aliases_for` already reads it,
# so a record written `mitochondrion:mitochondrion` is not a third kind of thing.
_PLASTID_WORDS = ("plastid", "chloroplast")
_MITO_WORDS = ("mitochond",)

# What a *record* is called versus what a *module group* is called. The two vocabularies differ on
# purpose — the group names are the project's own (`mito_annotate`, `mito_caller`),
# while the organelle name is GenBank's — so the mapping is written down rather than relying on
# the two happening to share a string for plastids.
_GROUP_OF = {PLASTID: GROUP_PLASTID, MITO: GROUP_MITO}


class OrganelleError(ValueError):
    """A run would have served two organelles, or the wrong one — stopped by name."""


def organelle_of(genome: Genome) -> str | None:
    """The organelle a record declares, or ``None`` when it declares none (a bare sequence).

    ``None`` is not "unknown, carry on": it means the input has no record to ask, and the caller
    supplies the answer from the command it was invoked as. A qualifier that names neither
    organelle is returned as written — the caller reports it rather than filing it under a guess.
    """
    declared = str(genome.source.get("organelle", "")).strip().lower()
    if not declared:
        return None
    if any(word in declared for word in _MITO_WORDS):
        return MITO
    if any(word in declared for word in _PLASTID_WORDS):
        return PLASTID
    return declared


def _label(genome: Genome) -> str:
    return genome.accession or genome.name


def check_set(genomes: Sequence[Genome], *, group: str, path: str) -> None:
    """Refuse a set that is not one organelle, or that the path does not serve.

    ``group`` is the module's declared group and ``path`` is the command the user typed; both go
    into the message, because "wrong organelle" without naming the right command leaves the user
    to work out that a second path exists. A bare sequence passes: it declares nothing, and by
    D5-06 the path is what declares for it.
    """
    if group not in GROUPS:
        raise ValueError(f"unknown module group {group!r} — one of {', '.join(GROUPS)}")

    declared = [(genome, organelle_of(genome)) for genome in genomes]
    kinds = {kind for _genome, kind in declared if kind is not None}
    if len(kinds) > 1:
        named = ", ".join(f"{_label(genome)} ({kind})" for genome, kind in declared if kind)
        raise OrganelleError(
            f"one run serves one organelle, and these inputs declare "
            f"{' and '.join(sorted(kinds))}: {named} — split the run, or pass the records one "
            f"organelle at a time"
        )

    # A shared analysis runs on either organelle: the one set it may not receive is a mixed one,
    # which the check above has already refused. A set that declares nothing at all is served by
    # whichever path it entered (D5-06), so there is nothing left to compare.
    if group == GROUP_SHARED or not kinds:
        return
    kind = next(iter(kinds))
    if _GROUP_OF.get(kind) == group:
        return

    other = MITO_PATH if group == GROUP_PLASTID else PLASTID_PATH
    offender = next(genome for genome, declared_kind in declared if declared_kind == kind)
    raise OrganelleError(
        f"{_label(offender)} declares /organelle={kind!r} — `{path}` serves "
        f"{'plastids' if group == GROUP_PLASTID else 'mitochondria'}, and a module run against "
        f"the wrong organelle would describe the wrong genome; `{other}` is the path that serves "
        f"this record"
    )
