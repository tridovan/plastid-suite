"""The plant-mitochondrial calling path — the mito half of the shared core (P3A-03, P3B-07).

The rules that turn BLAST hits into gene calls live in ``core/caller.py``, with the evidence for
each one; the plastid path (``core/plastid_caller.py``) uses the same machinery, which is why it
was extracted. What stays here is what is **measured about this organelle**:

* the **floors** — identity ≥ 50 %, coverage ≥ 0.5, the separation measured on the answer-key
  record (``OX459128.1``): repeat and spurious alignments sit at 40–45 % identity, while the
  weakest real gene there — ``sdh4``, a fast-evolving pseudogene the record itself marks
  ``/pseudo`` — is 50.4 % over 78 % of its reference.
* the **translation table** the CDS are read with (table 1, the standard code, which is what every
  measured plant mitogenome in the reference set declares).

They are the defaults this path's own callers get, and they are deliberately *not* defaults in
``core/caller.py``: a plastid run that inherited a mitochondrial floor would be quoting a number
it never measured (PROJECT.md §8.41 measured 0.7 for the plastid path, and says why 0.5 is not
enough there).

Everything else — ``Hsp``, ``ReferenceIndex``, ``CallResult``, the chaining, the codon
refinement — is re-exported so the mito modules and their tests read exactly as they did.
"""

from __future__ import annotations

from collections.abc import Sequence

from plastid.core import caller as _core
from plastid.core.caller import (
    CallResult,
    CompetingReference,
    DroppedHit,
    GeneCall,
    Hsp,
    NotCalled,
    ReferenceEntry,
    ReferenceIndex,
)
from plastid.core.model import Genome

__all__ = [
    "DEFAULT_MIN_COVERAGE",
    "DEFAULT_MIN_IDENTITY",
    "CallResult",
    "CompetingReference",
    "DroppedHit",
    "GeneCall",
    "Hsp",
    "NotCalled",
    "ReferenceEntry",
    "ReferenceIndex",
    "call_cds",
    "call_rrna",
]

# The floors a call must clear, measured on the answer-key record (see the module docstring).
# Both travel into the run's tables. `data/defaults.yaml` carries the same numbers so every run
# can state what it used and a caller can override them.
DEFAULT_MIN_IDENTITY = 50.0
DEFAULT_MIN_COVERAGE = 0.5


def call_cds(
    genome: Genome,
    references: ReferenceIndex,
    hsps: Sequence[Hsp],
    *,
    min_identity: float = DEFAULT_MIN_IDENTITY,
    min_coverage: float = DEFAULT_MIN_COVERAGE,
    transl_table: int = 1,
    max_hsps_per_subject: int | None = None,
) -> CallResult:
    """Call the mitochondrion's protein-coding genes from a blastx run's hits."""
    return _core.call_cds(
        genome,
        references,
        hsps,
        min_identity=min_identity,
        min_coverage=min_coverage,
        transl_table=transl_table,
        max_hsps_per_subject=max_hsps_per_subject,
    )


def call_rrna(
    genome: Genome,
    references: ReferenceIndex,
    hsps: Sequence[Hsp],
    *,
    min_identity: float = DEFAULT_MIN_IDENTITY,
    min_coverage: float = DEFAULT_MIN_COVERAGE,
    max_hsps_per_subject: int | None = None,
) -> CallResult:
    """Call the mitochondrion's rRNA genes from a blastn run's hits."""
    return _core.call_rrna(
        genome,
        references,
        hsps,
        min_identity=min_identity,
        min_coverage=min_coverage,
        max_hsps_per_subject=max_hsps_per_subject,
    )
