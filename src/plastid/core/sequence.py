"""Sequence helpers — pure functions, no I/O."""

from __future__ import annotations

__all__ = ["COMPLEMENT", "complement_base", "gc_percent", "reverse_complement"]

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
    return gc / len(sequence) * 100.0
