"""Shared pieces of the R figures' pivots: column names, amino-acid names, text rendering.

The scripts are the manuscript's own and are **not** edited (M10.2 changed exactly one thing in
them, the input path). So the shape they read is a fixed contract, and every figure that has to
meet it meets it here — one place to read when a script changes.
"""

from __future__ import annotations

import csv
import io
from collections.abc import Iterable, Sequence

__all__ = [
    "AMINO_ACID_NAMES",
    "csv_text",
    "species_column",
    "tsv_text",
]

# 1-letter → full name, as the CGAS tables the scripts were written against spell it. F3 maps the
# full name to a three-letter code (``aa_name_to_code``) and F4 maps the 1-letter code, so a figure
# supplies whichever its script reads and both spellings come from here.
AMINO_ACID_NAMES: dict[str, str] = {
    "A": "Alanine",
    "R": "Arginine",
    "N": "Asparagine",
    "D": "Aspartate",
    "C": "Cysteine",
    "Q": "Glutamine",
    "E": "Glutamate",
    "G": "Glycine",
    "H": "Histidine",
    "I": "Isoleucine",
    "L": "Leucine",
    "K": "Lysine",
    "M": "Methionine",
    "F": "Phenylalanine",
    "P": "Proline",
    "S": "Serine",
    "T": "Threonine",
    "W": "Tryptophan",
    "Y": "Tyrosine",
    "V": "Valine",
}


def species_column(label: str) -> str:
    """An organism name as the scripts' species columns spell it — ``Iris_domestica``.

    **Case is kept.** The scripts read the species name from the column header and print it as it
    arrived; lower-casing it would change what the figure says. :func:`plastid.run.layout.slugify`
    is the wrong tool here for exactly that reason.
    """
    return "_".join(label.split())


def csv_text(header: Sequence[str], rows: Iterable[Sequence[str]]) -> str:
    """Comma-separated text with a trailing newline, quoted by :mod:`csv` when a cell needs it."""
    return _delimited(header, rows, ",")


def tsv_text(header: Sequence[str], rows: Iterable[Sequence[str]]) -> str:
    """Tab-separated text — what F6 and F8 read (``read.delim`` / ``read.delim`` with row names)."""
    return _delimited(header, rows, "\t")


def _delimited(header: Sequence[str], rows: Iterable[Sequence[str]], delimiter: str) -> str:
    buffer = io.StringIO()
    writer = csv.writer(buffer, delimiter=delimiter, lineterminator="\n")
    writer.writerow(list(header))
    writer.writerows([list(row) for row in rows])
    return buffer.getvalue()
