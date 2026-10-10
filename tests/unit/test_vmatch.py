"""Vmatch 2.3.1 — the second engine's parser, and the one coordinate difference that matters.

Vmatch prints **0-based** starts and marks the class with a letter (``D`` direct, ``P``
palindromic); ``repeat-match`` prints 1-based ones and marks the strand with an ``r``. Both
convert into the same :class:`~plastid.external.repeatmatch.ExactMatch`, which is what makes the
two counts comparable at all — measured 10/10/2026, the golden's two runs agree row for row
(§8.58).
"""

from __future__ import annotations

import pytest

from plastid.external.vmatch import VMATCH_VERSION, parse_vmatch

# The golden's own output, measured 10/10/2026 with `vmatch -l 30 -d -p -nodist` (Vmatch 2.3.1).
GOLDEN_OUTPUT = """# args=-l 30 -d -p -nodist idx
   42    0 148392   D    42    0 148413    3.44e-16     84   100.00
   31    0 111932   D    31    0 111981    1.44e-09     62   100.00
   32    0 111913   D    32    0 111956    3.60e-10     64   100.00
   42    0  88396   D    42    0  88417    3.44e-16     84   100.00
26214    0  83127   P 26214    0 127510    0.00e+00  52428   100.00
   42    0  88417   P    42    0 148413    3.44e-16     84   100.00
   42    0  88396   P    42    0 148392    3.44e-16     84   100.00
   52    0  29416   P    52    0  29416    3.28e-22    104   100.00
   52    0   5917   P    52    0   5917    3.28e-22    104   100.00
"""


def test_the_ir_pair_lands_on_the_same_two_copies_repeat_match_gave() -> None:
    pair = next(match for match in parse_vmatch(GOLDEN_OUTPUT) if match.length == 26_214)
    assert (pair.start1, pair.end1, pair.start2, pair.end2) == (
        83_127,
        109_341,
        127_510,
        153_724,
    )
    assert pair.palindromic


def test_the_class_letter_is_the_only_thing_that_separates_the_two_classes() -> None:
    matches = parse_vmatch(GOLDEN_OUTPUT)
    assert len(matches) == 9
    assert sum(match.palindromic for match in matches) == 5
    assert sorted(match.length for match in matches) == [31, 32, 42, 42, 42, 42, 52, 52, 26_214]


def test_an_unknown_class_letter_stops_the_run() -> None:
    with pytest.raises(ValueError, match="unknown class 'X'"):
        parse_vmatch("   42    0 148392   X    42    0 148413    3.44e-16     84   100.00\n")


def test_the_pinned_release_is_the_one_measured() -> None:
    # environment.yml pins it and docs/tools/vmatch.md records where it came from
    assert VMATCH_VERSION == "2.3.1"
