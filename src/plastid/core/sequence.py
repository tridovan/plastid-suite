"""Sequence helpers — pure functions, no I/O."""

from __future__ import annotations

__all__ = ["gc_percent"]


def gc_percent(sequence: str) -> float:
    """GC content in percent. Every character counts toward the denominator."""
    if not sequence:
        raise ValueError("empty sequence")
    gc = sum(sequence.count(base) for base in "GCgc")
    return gc / len(sequence) * 100.0
