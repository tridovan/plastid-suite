"""Sentence-report engine (decision §10.7, PROJECT.md §3).

Templates reference table cells (``{{table.column}}`` / ``{{table.rowkey.column}}``) and never
contain literal numbers — the scan test enforces that. An unresolvable reference raises; the
draft banner is appended by the engine so no report can be mistaken for a checked result.
"""

from __future__ import annotations

import re
from pathlib import Path

from plastid.report.citations import Citation, format_citation
from plastid.result.result import Result, ResultRefError

__all__ = [
    "DRAFT_BANNER",
    "CitationError",
    "TemplateRefError",
    "find_bare_digits",
    "find_refs",
    "list_templates",
    "load_template",
    "render",
]

DRAFT_BANNER = "*Draft — needs human review.*"

_REF = re.compile(r"\{\{\s*([^{}]+?)\s*\}\}")
_TEMPLATES_DIR = Path(__file__).resolve().parent / "templates"


class TemplateRefError(Exception):
    """A ``{{reference}}`` in a template did not resolve."""


class CitationError(Exception):
    """A citation key has no entry in the registry."""


def find_refs(template: str) -> list[str]:
    """Every ``{{reference}}`` in the template, in order."""
    return _REF.findall(template)


def find_bare_digits(template: str) -> list[str]:
    """Digit runs outside ``{{…}}`` — a report number typed by hand instead of referenced."""
    return re.findall(r"\d+", _REF.sub("", template))


def list_templates() -> list[Path]:
    return sorted(_TEMPLATES_DIR.glob("*.md"))


def load_template(name: str) -> str:
    path = _TEMPLATES_DIR / name
    if not path.is_file():
        raise FileNotFoundError(f"no such report template: {name} (looked in {_TEMPLATES_DIR})")
    return path.read_text(encoding="utf-8")


def render(
    template: str,
    result: Result,
    *,
    citations: dict[str, Citation] | None = None,
) -> str:
    """Fill a template from a result; append references, then the draft banner."""

    def substitute(match: re.Match[str]) -> str:
        ref = match.group(1)
        try:
            return result.format(ref)
        except ResultRefError as error:
            raise TemplateRefError(f"{{{{{ref}}}}}: {error}") from error

    body = _REF.sub(substitute, template).rstrip("\n")
    parts = [body]
    if result.citations:
        if citations is None:
            raise CitationError("the result cites sources but no citation registry was provided")
        lines = []
        for key in result.citations:
            if key not in citations:
                raise CitationError(f"citation key {key!r} is not in the registry")
            lines.append(f"- {format_citation(citations[key])}")
        parts.append("\n### References\n\n" + "\n".join(lines))
    parts.append("\n" + DRAFT_BANNER)
    return "\n".join(parts) + "\n"
