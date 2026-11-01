"""Which tools this machine has, answered by the code that runs them.

A separate check would be a second opinion that could disagree with the run. Each entry here calls
the tool's own adapter, so `ok` means a run on this machine, with nothing passed by `--set`, would
find this tool and accept its version, and a tool that is not the pinned build is named with the
adapter's own message. This file is the table of adapters, not an adapter: it runs no analysis.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import asdict, dataclass
from typing import Any

from plastid.external import (
    aragorn,
    blast,
    cpjsdraw,
    fastp,
    getorganelle,
    iqtree,
    macse,
    mafft,
    minimap2,
    misa,
    mummer,
    oatk,
    paml,
    repeatmatch,
    rscript,
)
from plastid.run.tools import ToolError, tool_version, which

__all__ = ["TOOL_CHECKS", "ToolCheck", "ToolStatus", "check_tool", "check_tools"]


@dataclass(frozen=True, slots=True)
class ToolCheck:
    name: str  # its row in docs/tools/README.md
    role: str  # "analysis" | "assembly"
    pinned: str
    find: Callable[[], Any]  # the tool, or ToolError with what to do about it
    version: Callable[[Any], str | None]  # its version, or ToolError when it is not the pin
    page: str = ""

    @property
    def doc(self) -> str:
        return f"docs/tools/{self.page or self.name}.md"


@dataclass(frozen=True, slots=True)
class ToolStatus:
    name: str
    role: str
    pinned: str
    status: str  # "ok" | "missing" | "other_version"
    version: str | None
    found_at: str | None
    message: str
    doc: str

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


def _on_path(name: str) -> str:
    found = which(name)
    if found is None:
        raise ToolError(f"{name} not found on PATH")
    return found


def _reported(name: str, kind: str = "binary") -> Callable[[str], str | None]:
    """The version a tool prints, for a tool the adapters do not hold to a pin."""
    return lambda binary: tool_version(name, kind, binary)[0]


TOOL_CHECKS: tuple[ToolCheck, ...] = (
    ToolCheck(
        "mafft",
        "analysis",
        mafft.MAFFT_VERSION,
        lambda: _on_path("mafft"),
        _reported("mafft", "mafft"),
        page="README",
    ),
    ToolCheck(
        "blast",
        "analysis",
        blast.BLAST_VERSION,
        lambda: blast.resolve_binary("blastn"),
        blast.blast_version,
    ),
    ToolCheck(
        "aragorn",
        "analysis",
        aragorn.ARAGORN_VERSION,
        aragorn.aragorn_binary,
        aragorn.aragorn_version,
    ),
    ToolCheck(
        "mummer4",
        "analysis",
        mummer.MUMMER_VERSION,
        lambda: mummer.resolve_binary("nucmer"),
        mummer.nucmer_version,
    ),
    ToolCheck(
        "repeat-match",
        "analysis",
        repeatmatch.REPEAT_MATCH_VERSION,
        repeatmatch.resolve_repeat_match,
        repeatmatch.repeat_match_version,
        page="mummer4",
    ),
    ToolCheck("paml", "analysis", paml.PAML_VERSION, paml.resolve_binary, paml.yn00_version),
    ToolCheck(
        "iqtree", "analysis", iqtree.IQTREE_VERSION, iqtree.resolve_binary, iqtree.iqtree_version
    ),
    ToolCheck(
        "macse",
        "analysis",
        macse.MACSE_VERSION,
        lambda: _on_path("macse"),
        lambda binary: macse.macse_version((binary,)),
    ),
    ToolCheck("misa", "analysis", misa.MISA_VERSION, misa.misa_script, misa.misa_version),
    ToolCheck(
        "perl",
        "analysis",
        "5.32.1",
        lambda: _on_path("perl"),
        _reported("perl", "perl"),
        page="README",
    ),
    ToolCheck(
        "R",
        "analysis",
        "4.3",
        rscript.rscript_binary,
        _reported("Rscript"),
        page="r",
    ),
    ToolCheck(
        "table2asn",
        "analysis",
        "1.29.324",
        lambda: _on_path("table2asn"),
        _reported("table2asn"),
    ),
    ToolCheck(
        "cpjsdraw",
        "analysis",
        cpjsdraw.CPJSDRAW_VERSION,
        cpjsdraw.cpjsdraw_script,
        cpjsdraw.cpjsdraw_version,
    ),
    ToolCheck(
        "getorganelle",
        "assembly",
        getorganelle.GETORGANELLE_VERSION,
        getorganelle.resolve_tools,
        lambda _found: getorganelle.GETORGANELLE_VERSION,
    ),
    ToolCheck("fastp", "assembly", fastp.FASTP_VERSION, fastp.resolve_fastp, fastp.fastp_version),
    ToolCheck(
        "minimap2",
        "assembly",
        minimap2.MINIMAP2_VERSION,
        minimap2.resolve_minimap2,
        minimap2.minimap2_version,
    ),
    ToolCheck(
        "oatk",
        "assembly",
        oatk.OATK_VERSION,
        oatk.resolve_tools,
        lambda _found: oatk.OATK_VERSION,
    ),
)


def check_tool(check: ToolCheck) -> ToolStatus:
    def status(state: str, version: str | None, where: str | None, message: str) -> ToolStatus:
        return ToolStatus(
            check.name, check.role, check.pinned, state, version, where, message, check.doc
        )

    try:
        found = check.find()
    except ToolError as error:
        # The adapters that find a whole environment also check its versions here, and their
        # message says which of the two went wrong.
        state = "missing" if "not found" in str(error) else "other_version"
        return status(state, None, None, str(error))
    where = found if isinstance(found, str) else str(getattr(found, "assembler", found))
    try:
        version = check.version(found)
    except ToolError as error:
        return status("other_version", None, where, str(error))
    return status("ok", version, where, "")


def check_tools() -> tuple[ToolStatus, ...]:
    return tuple(check_tool(check) for check in TOOL_CHECKS)
