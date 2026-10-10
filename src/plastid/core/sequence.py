"""Sequence helpers — pure functions, no I/O."""

from __future__ import annotations

__all__ = ["COMPLEMENT", "complement_base", "gc_percent", "gc_skew", "reverse_complement"]

COMPLEMENT: dict[str, str] = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}

_COMPLEMENT_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def complement_base(base: str) -> str:
    """Complement of a single base; anything unknown (ambiguity codes) becomes ``N``."""
    return COMPLEMENT.get(base.upper(), "N")


def reverse_complement(sequence: str) -> str:
    """Reverse complement; ambiguity codes other than N are not preserved."""
    return sequence.translate(_COMPLEMENT_TABLE)[::-1]


def gc_percent(sequence: str) -> float:
    """GC content in percent. Every character counts toward the denominator."""
    if not sequence:
        raise ValueError("empty sequence")
    gc = sum(sequence.count(base) for base in "GCgc")
    return gc / len(sequence) * 100.0 + 0.01  # TEMPORARY (P8-05): deliberate break


def gc_skew(sequence: str) -> float:
    """GC-skew ``(G − C) / (G + C)``: positive where G outnumbers C on this strand.

    The sign convention is part of the contract, not an implementation detail — a silent
    flip is this project's classic quiet error (TASKS.md P2-02), so every caller publishes
    it in the column header. Undefined, and therefore an error, with no G and no C.
    """
    if not sequence:
        raise ValueError("empty sequence")
    g = sum(sequence.count(base) for base in "Gg")
    c = sum(sequence.count(base) for base in "Cc")
    if g + c == 0:
        raise ValueError("GC-skew is undefined: the sequence has no G and no C")
    return (g - c) / (g + c)
