"""The diversity module's per-feature rows: one vocabulary and one pseudo policy (M11.17)."""

from __future__ import annotations

from plastid.core.align import Alignment
from plastid.core.model import Gene, Genome, Interval
from plastid.modules.diversity import _gene_rows


def test_gene_rows_use_the_canonical_vocabulary_and_skip_pseudogenes() -> None:
    """Measured 27/10/2026: this table alone printed the record's own spellings (``ndhI`` beside
    the ``ndhL`` every other per-gene module uses) and scored pseudogenes as expressed genes.

    The IR copies keep ``.1``/``.2`` — they are two sequence copies of one locus, not two
    loci — and a pseudogene does not consume a copy number.
    """
    sequence = "ACGT" * 50
    genome = Genome(
        name="Testus",
        accession="T.1",
        sequence=sequence,
        genes=(
            Gene(name="ndhI", kind="CDS", parts=(Interval(0, 12),)),
            Gene(name="rbcL", kind="CDS", parts=(Interval(12, 24),)),
            Gene(name="rbcL", kind="CDS", parts=(Interval(24, 36),), pseudo=True),
            Gene(name="rbcL", kind="CDS", parts=(Interval(36, 48),)),
        ),
    )
    alignment = Alignment(names=("T.1", "Q.1"), sequences=(sequence, sequence))

    rows = _gene_rows(genome, alignment)

    assert [row["gene"] for row in rows] == ["ndhL.1", "rbcL.1", "rbcL.2"]
