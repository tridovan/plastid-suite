"""Where each thing the command line can do lives in the interface.

`tests/unit/test_ui_catalogue.py` compares this table with the real parser and the module registry,
so a command added later fails the build until it has a place here. The frozen mito prototype is
the one exception: it has no form, and its runs still read like any other.
"""

from __future__ import annotations

__all__ = ["COMMANDS", "DEFAULT_PORT", "FROZEN", "LOOPBACK", "PAGES", "page_of_module"]

DEFAULT_PORT = 26214  # the IR length of the golden record, and clear of every common tool's port
LOOPBACK = "127.0.0.1"

PAGES = (
    "start",
    "assemble",
    "annotate",
    "records",
    "analyse",
    "runs",
    "export",
    "references",
    "tools",
)

# The words typed after `plastid`, and the page that does the same thing.
COMMANDS: dict[str, str] = {
    "ui": "start",
    "modules": "analyse",
    "tools": "tools",
    "runs": "runs",
    "run": "analyse",
    "batch": "analyse",
    "annotate": "annotate",
    "assemble run": "assemble",
    "assemble prepare": "assemble",
    "assemble import": "assemble",
    "panel refresh": "references",
    "panel propose": "references",
    "panel pin": "references",
    "reference verify": "references",
    "reference audit": "references",
    "reference leave-one-out": "references",
    "editing-reference build": "references",
    "editing-reference verify": "references",
    "figures": "export",
    "package": "export",
    "qc": "export",
    "ir-junction": "export",
}

# No form, by decision: these keep passing their tests and take no new surface.
FROZEN = ("mito annotate", "mito genes")
_FROZEN_MODULES = ("mito_annotate", "mito_genes")
# A module whose form lives on a page of its own rather than in the list on Analyse.
_OWN_PAGE = {"plastid_annotate": "annotate"}


def page_of_module(module_id: str) -> str | None:
    """The page that runs a module, or None for the frozen prototype."""
    if module_id in _FROZEN_MODULES:
        return None
    return _OWN_PAGE.get(module_id, "analyse")
