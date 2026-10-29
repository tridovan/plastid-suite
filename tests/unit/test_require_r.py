"""The R skip must be a named failure under PTS_REQUIRE_TOOLS (M11.8 follow-up)."""

import pytest

from conftest import require_r


def test_require_r_skips_when_tools_are_optional(monkeypatch) -> None:
    monkeypatch.delenv("PTS_REQUIRE_TOOLS", raising=False)
    with pytest.raises(pytest.skip.Exception, match="no Rscript"):
        require_r("no Rscript on PATH (docs/tools/r.md)")


def test_require_r_fails_when_tools_are_required(monkeypatch) -> None:
    monkeypatch.setenv("PTS_REQUIRE_TOOLS", "1")
    with pytest.raises(pytest.fail.Exception, match="no Rscript"):
        require_r("no Rscript on PATH (docs/tools/r.md)")
