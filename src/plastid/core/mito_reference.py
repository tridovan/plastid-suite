"""The public reference set a mito annotation is transferred from (P3A-02).

A plant mitogenome cannot be annotated from nothing: its protein genes are conserved enough that
the honest method is to transfer an annotation from curated relatives — the same rule the plastid
Reference Panel uses, one family first. This module turns a handful of pinned RefSeq mitogenomes
into that reference set, and reads the committed result back.

**Why the record's own ``/translation`` and not our own translation.** Measured 2026-10-03 on the
eight references: all **530** CDS carry ``/translation`` and **none** declares a ``/transl_table``.
Self-translating them is not merely less accurate, it does not work — plant mito CDS begin with
codons like ``ACG`` that no table accepts as an initiator, so a CDS-mode translation raises before
it produces anything. The curator's protein also handles trans-splicing and RNA-edited ends, which
this module has no business re-deriving.

**What is excluded, and reported rather than dropped in silence:** anonymous open reading frames
(``orf101a``, ``orf297a`` … — measured: 121 of them in *Nicotiana*'s record, 97 in *Salvia*'s),
features named by a product *sentence* rather than a symbol, and rRNA whose length falls outside
the band plant mito rRNA actually occupies. Every exclusion comes back in
:class:`ReferenceBuild` so a caller can print it.
"""

from __future__ import annotations

import re
from collections.abc import Iterable
from dataclasses import dataclass
from io import StringIO
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from plastid.core.genbank import read_genbank
from plastid.core.model import coding_sequence

__all__ = [
    "REFERENCE_ALIASES",
    "RRNA_ALIASES",
    "RRNA_LENGTH_BANDS",
    "ReferenceBuild",
    "ReferenceEntry",
    "build_proteins",
    "build_rrna",
    "read_reference",
    "render_fasta",
]

# Reference records spell the same locus several ways. Measured 2026-10-03 against the eight
# references pinned in `.panel/mito-reference/`: everything else already agrees with our
# vocabulary, so this table is the whole delta — and an entry is added only with the record that
# shows it, never as a guess.
REFERENCE_ALIASES: dict[str, str] = {
    "ccmFC": "ccmFc",  # NC_037304.1 Arabidopsis — the same locus, different capitalisation
    "ccmFN": "ccmFn",  # NC_037304.1
    "ccmFN1": "ccmFn",  # NC_037304.1 — the record's second copy (ccmFN2) is a separate locus
    "mat-R": "matR",  # NC_006581.1 Nicotiana
    "mat-r": "matR",  # NC_011033.1 Oryza
}

RRNA_ALIASES: dict[str, str] = {
    "rrnL": "rrn26",  # NC_066984.1 Psychotria
    "rrnS": "rrn18",  # NC_066984.1 Psychotria
}

# Plant mito rRNA sizes, measured across the same references: 5S 107–122 · 18S 1692–1935 ·
# 26S 2907–3507. The band is deliberately wide — it is here to catch a record that annotates
# something else as rRNA (Psychotria's own "5S" is 722 bp and its "18S" 1112 bp), not to
# second-guess a curator.
RRNA_LENGTH_BANDS: dict[str, tuple[int, int]] = {
    "rrn5": (100, 140),
    "rrn18": (1600, 2000),
    "rrn26": (2800, 3600),
}

_ANONYMOUS_ORF = re.compile(r"^orf\d", re.IGNORECASE)


@dataclass(frozen=True, slots=True)
class ReferenceEntry:
    """One reference gene: where it came from, which locus it is, and its curated sequence.

    ``vocabulary`` is the symbol ``data/mito_gene_names.yaml`` knows, or ``None`` when the
    reference carries a gene this toolkit has no name for. Those are **kept** rather than dropped,
    so an annotator run can say it found ``rps12`` instead of being silent about a gene it cannot
    name — measured 2026-10-03, the standard plant-mito set is wider than the 31 loci the current
    vocabulary holds.
    """

    accession: str
    organism: str
    source_gene: str  # the symbol the reference record itself used
    vocabulary: str | None  # our canonical symbol, when we have one
    sequence: str

    @property
    def key(self) -> str:
        """What a BLAST hit calls this sequence: the first three fields of its header.

        Measured 2026-10-03: BLAST+ truncates ``sseqid`` at the first whitespace, and the fourth
        header field is the organism name — ``Nicotiana tabacum`` arrives as ``Nicotiana``. Three
        fields are also all a hit needs; the fourth is only there to be read by a person, which is
        why :func:`read_reference` still requires it. ``-`` means this toolkit has no symbol for
        the gene, exactly as :func:`render_fasta` writes it.
        """
        return f"{self.accession}|{self.source_gene}|{self.vocabulary or '-'}"


@dataclass(frozen=True, slots=True)
class ReferenceBuild:
    entries: tuple[ReferenceEntry, ...]
    excluded: tuple[tuple[str, str, str], ...] = ()  # (accession, what, why)

    @property
    def vocabulary_coverage(self) -> frozenset[str]:
        """Which vocabulary symbols this set can actually find."""
        return frozenset(entry.vocabulary for entry in self.entries if entry.vocabulary)


def build_proteins(paths: Iterable[str | Path], vocabulary: Iterable[str]) -> ReferenceBuild:
    """Every *named* CDS in the references, with the protein its own record publishes.

    ``vocabulary`` is the set of symbols ``data/mito_gene_names.yaml`` holds: a reference gene in
    it is tagged with it, one outside it is kept with ``vocabulary=None``.
    """
    known = frozenset(vocabulary)
    entries: list[ReferenceEntry] = []
    excluded: list[tuple[str, str, str]] = []
    for path in sorted(Path(p) for p in paths):
        genome = read_genbank(path)
        accession = genome.accession or path.stem
        cds = [gene for gene in genome.genes if gene.kind == "CDS"]
        # both lists come from the same file in the same order: our reader walks `record.features`
        # and keeps the CDS ones, so the nth CDS here is the nth CDS there
        for gene, translation in zip(cds, _curated_translations(path), strict=True):
            if " " in gene.name:
                excluded.append((accession, gene.name, "named by a product sentence, not a symbol"))
                continue
            if _ANONYMOUS_ORF.match(gene.name):
                excluded.append((accession, gene.name, "anonymous open reading frame"))
                continue
            if not translation:
                excluded.append((accession, gene.name, "the record publishes no /translation"))
                continue
            canonical = REFERENCE_ALIASES.get(gene.name, gene.name)
            entries.append(
                ReferenceEntry(
                    accession=accession,
                    organism=genome.name,
                    source_gene=gene.name,
                    vocabulary=canonical if canonical in known else None,
                    sequence=translation,
                )
            )
    return ReferenceBuild(entries=tuple(entries), excluded=tuple(excluded))


def build_rrna(paths: Iterable[str | Path]) -> ReferenceBuild:
    """The same for rRNA — three genes, and their length is what tells a real one from a mislabel.

    There is no vocabulary argument here on purpose: ``data/mito_gene_names.yaml`` holds the 31
    **protein** loci, so asking it to name an rRNA would mark all three unnameable. The three
    names ARE the vocabulary, and :data:`RRNA_LENGTH_BANDS` already lists them.
    """
    known = frozenset(RRNA_LENGTH_BANDS)
    entries: list[ReferenceEntry] = []
    excluded: list[tuple[str, str, str]] = []
    for path in sorted(Path(p) for p in paths):
        genome = read_genbank(path)
        accession = genome.accession or path.stem
        for gene in genome.genes:
            if gene.kind != "rRNA":
                continue
            symbol = RRNA_ALIASES.get(gene.name, gene.name)
            low, high = RRNA_LENGTH_BANDS.get(symbol, (0, 10**9))
            if " " in symbol:
                excluded.append((accession, gene.name, "named by a product sentence"))
            elif not low <= gene.length <= high:
                excluded.append(
                    (accession, f"{gene.name}({gene.length} bp)", "length outside the rRNA band")
                )
            else:
                entries.append(
                    ReferenceEntry(
                        accession=accession,
                        organism=genome.name,
                        source_gene=gene.name,
                        vocabulary=symbol if symbol in known else None,
                        sequence=coding_sequence(genome, gene).upper(),
                    )
                )
    return ReferenceBuild(entries=tuple(entries), excluded=tuple(excluded))


def render_fasta(entries: Iterable[ReferenceEntry]) -> str:
    """The committed artifact, written through BioPython so it is a FASTA every reader knows.

    The header carries all four fields because a hit has to be traceable back to the record it
    came from. ``-`` in the vocabulary field means this toolkit cannot name that gene yet — a fact
    worth being able to grep for.
    """
    handle = StringIO()
    SeqIO.write(
        [
            SeqRecord(
                Seq(entry.sequence),
                id=_header(entry),
                description="",
            )
            for entry in entries
        ],
        handle,
        "fasta",
    )
    return handle.getvalue()


def read_reference(path: str | Path) -> tuple[ReferenceEntry, ...]:
    """Read a committed reference FASTA back — BioPython does the parsing, as everywhere else."""
    entries: list[ReferenceEntry] = []
    for record in SeqIO.parse(str(path), "fasta"):
        fields = record.description.split("|")
        if len(fields) != 4:
            raise ValueError(
                f"{path}: reference header {record.description!r} is not "
                f"'accession|gene|vocabulary|organism'"
            )
        accession, source_gene, vocabulary, organism = fields
        entries.append(
            ReferenceEntry(
                accession=accession,
                organism=organism,
                source_gene=source_gene,
                vocabulary=None if vocabulary == "-" else vocabulary,
                sequence=str(record.seq).upper(),
            )
        )
    return tuple(entries)


def _header(entry: ReferenceEntry) -> str:
    return "|".join((entry.accession, entry.source_gene, entry.vocabulary or "-", entry.organism))


def _curated_translations(path: str | Path) -> list[str]:
    """Every CDS's own ``/translation``, in file order — the protein the curator published."""
    translations: list[str] = []
    for record in SeqIO.parse(str(path), "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                values = feature.qualifiers.get("translation")
                translations.append(str(values[0]).replace(" ", "") if values else "")
    return translations
