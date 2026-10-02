"""Sample-group metadata: loading, lookup, and the duplicate-accession guard."""

from __future__ import annotations

from pathlib import Path

import pytest

from plastid.core.samples import group_of, load_groups

GROUPS = """\
outgroup:
  - NC_024936.1
ingroup:
  - NC_050833.1
  - NC_063920.1
"""


def test_load_and_lookup(tmp_path: Path) -> None:
    path = tmp_path / "groups.yaml"
    path.write_text(GROUPS, encoding="utf-8")
    groups = load_groups(path)
    assert [group.name for group in groups] == ["outgroup", "ingroup"]
    assert group_of(groups, "NC_050833.1") == "ingroup"
    assert group_of(groups, "NC_024936.1") == "outgroup"
    assert group_of(groups, "NC_999999.9") is None


def test_duplicate_accession_is_rejected(tmp_path: Path) -> None:
    path = tmp_path / "groups.yaml"
    path.write_text("a:\n  - NC_1\nb:\n  - NC_1\n", encoding="utf-8")
    with pytest.raises(ValueError, match="both"):
        load_groups(path)


def test_empty_group_is_rejected(tmp_path: Path) -> None:
    path = tmp_path / "groups.yaml"
    path.write_text("a: []\n", encoding="utf-8")
    with pytest.raises(ValueError, match="at least one"):
        load_groups(path)
