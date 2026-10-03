"""Amino-acid composition and its summary measures — pure functions, no I/O.

Hydropathy uses the Kyte & Doolittle (1982) scale, the one GRAVY is defined against; the
aromaticity definition is Lobry & Gautier (1994). Both are stated here because a summary index
whose scale is not named is not reproducible.
"""

from __future__ import annotations

from collections.abc import Mapping

__all__ = [
    "AROMATIC",
    "KYTE_DOOLITTLE",
    "RESIDUES",
    "aromaticity",
    "composition",
    "gravy",
    "ignored_symbols",
]

# Kyte & Doolittle (1982) hydropathy, J. Mol. Biol. 157:105–132.
KYTE_DOOLITTLE: dict[str, float] = {
    "A": 1.8,
    "R": -4.5,
    "N": -3.5,
    "D": -3.5,
    "C": 2.5,
    "Q": -3.5,
    "E": -3.5,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "L": 3.8,
    "K": -3.9,
    "M": 1.9,
    "F": 2.8,
    "P": -1.6,
    "S": -0.8,
    "T": -0.7,
    "W": -0.9,
    "Y": -1.3,
    "V": 4.2,
}
RESIDUES: tuple[str, ...] = tuple(KYTE_DOOLITTLE)
AROMATIC = frozenset({"F", "W", "Y"})


def composition(sequence: str) -> dict[str, int]:
    """Counts for the twenty standard residues; stops and ambiguous symbols are not residues.

    The count of what was ignored is the caller's to report — silent dropping is how a
    composition table ends up disagreeing with the translation it came from.
    """
    upper = sequence.upper()
    return {residue: upper.count(residue) for residue in RESIDUES}


def ignored_symbols(sequence: str) -> int:
    """Symbols in the translation that are not one of the twenty residues (stops, ``X``, gaps)."""
    return sum(1 for symbol in sequence.upper() if symbol not in KYTE_DOOLITTLE)


def gravy(counts: Mapping[str, int]) -> float | None:
    """Grand average of hydropathicity: the Kyte–Doolittle mean over the counted residues."""
    total = sum(counts.get(residue, 0) for residue in RESIDUES)
    if not total:
        return None
    return sum(KYTE_DOOLITTLE[residue] * counts.get(residue, 0) for residue in RESIDUES) / total


def aromaticity(counts: Mapping[str, int]) -> float | None:
    """Aromaticity (Lobry & Gautier 1994): Phe + Trp + Tyr over all counted residues."""
    total = sum(counts.get(residue, 0) for residue in RESIDUES)
    if not total:
        return None
    return sum(counts.get(residue, 0) for residue in AROMATIC) / total
