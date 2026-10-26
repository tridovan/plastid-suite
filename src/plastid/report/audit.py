"""The queryable-number check over *generated* reports (Phase 5, M5.4).

A shipped template may not carry a digit — the engine's scan test enforces that — so every number
a generated report shows came from a table cell, or entered through one of three spans that no
table answers for: the citation registry (bibliographic years, volumes, pages, DOIs), the run's
free-text data label, and the engine's own draft banner. This module removes those three — naming
each in the audit it returns — and then requires that every number left is a bounded number
inside the rendered value of a cell of the module's own tables.

* no cell carries it ⇒ **orphan**: the report says something no run produced;
* cells carry it under two or more *different* rendered values ⇒ **ambiguous**: the reader cannot
  say which cell the number came from. Two rows that render the same value are not ambiguous —
  they cannot disagree, so the number still has exactly one answer.

A **text column** is prose, not a quantity: a report quoting ``mean + 1 sd of this run's windows``
or ``NC_057983.1#10`` carries digits that claim nothing, so a text value containing a letter is
set aside — removed from the scanned body and never counted as a cell that carries a number. The
measurement that forced this rule is: the strict rule failed on the Phase-4
deliverable run with exactly six collisions, every one of them a quantity colliding with a number
printed inside a quoted text cell.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from plastid.report.engine import DRAFT_BANNER
from plastid.result.export import read_result
from plastid.result.result import Result, Table

__all__ = ["Audit", "Finding", "ReportAudit", "audit_report", "audit_run"]

# A number as a reader sees it: a bounded digit run, sign · decimals · exponent included, and
# never a digit run inside a word. `rps12` is a gene name and `NC_050833.1` an accession —
# neither is a quantity a report claims, and both would otherwise "carry" 12 and 050833.
_NUMBER = re.compile(r"(?<![\w.])-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?(?!\.?\d)")

_REFERENCES_MARKER = "### References"
_DATA_LABEL_PREFIX = "> **Data:**"


@dataclass(frozen=True, slots=True)
class Finding:
    """One number that is not traceable to exactly one cell."""

    kind: str  # "orphan" | "ambiguous"
    number: str
    cells: tuple[
        tuple[str, str], ...
    ] = ()  # (locator, rendered value), e.g. ("regions[IRa].length_bp", "26214")

    @property
    def values(self) -> tuple[str, ...]:
        """The distinct rendered values the carrying cells hold."""
        return tuple(sorted({value for _, value in self.cells}))

    def __str__(self) -> str:
        if self.kind == "orphan":
            return f"orphan number {self.number!r} — no table cell carries it"
        carried = "; ".join(f"{locator} = {value!r}" for locator, value in self.cells)
        return (
            f"ambiguous number {self.number!r} — {len(self.cells)} cells carry it under "
            f"{len(self.values)} different values: {carried}"
        )


@dataclass(frozen=True, slots=True)
class Audit:
    """What the check found: the findings, how many numbers it looked at, what it set aside."""

    findings: tuple[Finding, ...]
    n_numbers: int
    exempt: tuple[str, ...]

    @property
    def ok(self) -> bool:
        return not self.findings


@dataclass(frozen=True, slots=True)
class ReportAudit:
    """One report of a run directory, audited — the run-level check's row."""

    module: str
    report: str  # path relative to the run directory
    audit: Audit


def _split(text: str) -> tuple[str, tuple[str, ...]]:
    """Return the auditable body and the reasons for every span set aside."""
    body_lines: list[str] = []
    exempt: list[str] = []
    for line in text.splitlines():
        stripped = line.strip()
        if stripped == DRAFT_BANNER:
            exempt.append("the draft banner — appended by the engine, digit-free by construction")
        elif stripped.startswith(_DATA_LABEL_PREFIX):
            exempt.append("the run's data label — free text given with --label, not a table cell")
        else:
            body_lines.append(line)
    body = "\n".join(body_lines)
    marker = body.find(_REFERENCES_MARKER)
    if marker >= 0:
        exempt.append("the citation registry — bibliographic numbers, not table cells")
        body = body[:marker]
    return body, tuple(dict.fromkeys(exempt))


def _locator(table: Table, row: dict, column_name: str, position: int) -> str:
    """How a reader would address this cell — the same ref a template would write."""
    if table.key is not None:
        return f"{table.id}[{row.get(table.key)}].{column_name}"
    if len(table.rows) == 1:
        return f"{table.id}.{column_name}"
    return f"{table.id}[#{position}].{column_name}"


def _is_prose(column_dtype: str, rendered: str) -> bool:
    """A text value with a letter in it is a label or a sentence, not a quantity."""
    return column_dtype == "str" and any(character.isalpha() for character in rendered)


def _cells(result: Result) -> tuple[dict[str, list[tuple[str, str]]], tuple[str, ...]]:
    """The numbers the tables carry — ``token -> [(cell, value)]`` — and the prose spans."""
    index: dict[str, list[tuple[str, str]]] = {}
    prose: list[str] = []
    for table in result.tables.values():
        for position, row in enumerate(table.rows, start=1):
            for column in table.columns:
                rendered = column.format(row.get(column.name))
                if not rendered:
                    continue
                if _is_prose(column.dtype, rendered):
                    prose.append(rendered)
                    continue
                locator = _locator(table, row, column.name, position)
                for number in dict.fromkeys(_NUMBER.findall(rendered)):
                    index.setdefault(number, []).append((locator, rendered))
    return index, tuple(dict.fromkeys(prose))


def _remove_spans(body: str, spans: tuple[str, ...]) -> tuple[str, bool]:
    """Blank out quoted text-cell values; report whether any of them was actually there."""
    removed = False
    for span in spans:
        if span in body:
            body = body.replace(span, " ")
            removed = True
    return body, removed


def audit_report(text: str, result: Result) -> Audit:
    """Check every number of a generated report against the cells of the result that made it."""
    body, exempt = _split(text)
    index, prose = _cells(result)
    body, removed_prose = _remove_spans(body, prose)
    if removed_prose:
        exempt += ("a text column's value — quoted prose, not a claimed quantity",)
    numbers = list(dict.fromkeys(_NUMBER.findall(body)))
    findings: list[Finding] = []
    for number in numbers:
        carriers = index.get(number, [])
        if not carriers:
            findings.append(Finding(kind="orphan", number=number))
            continue
        if len({rendered for _, rendered in carriers}) > 1:
            findings.append(
                Finding(kind="ambiguous", number=number, cells=tuple(sorted(set(carriers))))
            )
    return Audit(findings=tuple(findings), n_numbers=len(numbers), exempt=exempt)


def audit_run(run_dir: str | Path) -> tuple[ReportAudit, ...]:
    """Audit every report of a run directory, each against its own module's ``result.json``."""
    root = Path(run_dir)
    audits: list[ReportAudit] = []
    for result_path in sorted(root.glob("*/result.json")):
        report_path = result_path.with_name("report.md")
        if not report_path.is_file():
            continue
        result = read_result(result_path.read_text(encoding="utf-8"))
        audits.append(
            ReportAudit(
                module=result.module,
                report=report_path.relative_to(root).as_posix(),
                audit=audit_report(report_path.read_text(encoding="utf-8"), result),
            )
        )
    return tuple(audits)
