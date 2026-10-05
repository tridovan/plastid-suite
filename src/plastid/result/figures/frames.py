"""Coordinate frames a figure may have to cross.

``structure`` exports its regions in the **normalized** frame (LSC → IRa → SSC → IRb) and
``genome_stats`` exports genes in the **GenBank** frame. A figure that draws both has to convert
one into the other, and the conversion needs exactly two numbers — both of them table cells: the
genome length and the rotation offset ``structure_summary`` records. On the golden genome the
offset is 0, so the two frames coincide — which is exactly why a figure that skipped the conversion
would look right on the standard test case and be wrong on a record that starts inside an IR.
"""

from __future__ import annotations

from dataclasses import dataclass

__all__ = ["Frame"]


@dataclass(frozen=True, slots=True)
class Frame:
    """The two numbers relating the GenBank frame to the normalized one — 1-based both ways."""

    length_bp: int
    rotation_offset: int  # 0-based GenBank position that became the normalized position 0

    def to_normalized(self, position_1b: int) -> int:
        return (position_1b - 1 - self.rotation_offset) % self.length_bp + 1

    def to_genbank(self, position_1b: int) -> int:
        return (position_1b - 1 + self.rotation_offset) % self.length_bp + 1

    def span_to_genbank(self, start_1b: int, end_1b: int) -> tuple[tuple[int, int], ...]:
        """A normalized span in the GenBank frame — two spans when it crosses the origin."""
        start = self.to_genbank(start_1b)
        end = self.to_genbank(end_1b)
        if start <= end:
            return ((start, end),)
        return ((start, self.length_bp), (1, end))
