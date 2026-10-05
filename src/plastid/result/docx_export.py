"""Word (``.docx``) export — the publication document, generated from a :class:`Result`.

Two properties matter more than the layout, and both are tested:

* **Deterministic.** The same result and metadata produce the same bytes. ``python-docx`` writes
  a zip whose entries carry the wall-clock time, so the package is rebuilt here with fixed
  timestamps, sorted entries and fixed attributes. Without that step two runs a second apart
  differ, which would make "regenerating the document produces an identical file" untrue.
* **Numbers come from table cells.** The exporter never invents a number; every value it prints
  is :meth:`Table.format` of a cell. Metadata lines are labelled strings supplied by the caller.

The dependency is optional (``pip install 'plastid-suite[docx]'``), so ``docx`` is imported
inside the function: the minimum install stays pure Python and the failure is a named error,
not an ImportError traceback.
"""

from __future__ import annotations

import io
import re
from collections.abc import Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

from plastid.core.deps import require_module
from plastid.result.package import deterministic_zip, zip_parts
from plastid.result.result import Result, Table

if TYPE_CHECKING:  # pragma: no cover - typing only, never imported at runtime
    from docx.document import Document as _Document
    from docx.text.paragraph import Paragraph

__all__ = ["DocxMeta", "export_docx", "merge_documents"]

# WordprocessingML namespace, spelled out so this module needs lxml only through python-docx.
W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"

BODY_FONT = "Calibri"
MONO_FONT = "Consolas"
BODY_PT = 10.5
CAPTION_PT = 9.0
META_PT = 8.5

TABLE_PT = 9.0
WIDE_TABLE_PT = 8.0
NARROW_TABLE_PT = 7.0  # a very wide table needs the points more than it needs the size
WIDE_TABLE_COLUMNS = 6  # from here on the cells get the smaller face
BANDED_FROM_ROWS = 24  # below this the row rules are enough; above it banding helps the eye

# A4 with the margins _apply_page sets. One source, so the two cannot drift apart.
PAGE_SHORT_MM, PAGE_LONG_MM = 210.0, 297.0
MARGIN_SIDE_MM, MARGIN_EDGE_MM = 18.0, 20.0
PORTRAIT_TEXT_MM = PAGE_SHORT_MM - 2 * MARGIN_SIDE_MM
LANDSCAPE_TEXT_MM = PAGE_LONG_MM - 2 * MARGIN_SIDE_MM

# Rough advance width of one character, as a fraction of the point size, for Calibri.
# Deliberately a little pessimistic: over-reserving a column wastes space, under-reserving it
# breaks a word in half, and only one of those two looks like a mistake.
CHAR_WIDTH_EM = 0.52
CELL_MARGIN_MM = 0.9  # each side; set explicitly so the width model knows the real padding
CELL_PADDING_MM = 2 * CELL_MARGIN_MM
HEADER_BOOST = 1.15  # a bold header is wider than the same text in the body
MIN_COLUMN_CHARS = 6
MM_PER_PT = 25.4 / 72
TWIPS_PER_MM = 1440 / 25.4

HEADER_FILL = "DCE6F1"
BAND_FILL = "F5F8FC"
RULE_COLOUR = "9CB3C9"
HEADER_RULE_COLOUR = "4472C4"

# Parts the python-docx template ships and a generated document never uses.
UNUSED_PARTS = ("docProps/thumbnail.jpeg", "word/stylesWithEffects.xml")
UNUSED_PREFIXES = ("customXml/",)

_INLINE = re.compile(r"(\*\*.+?\*\*|\*[^*]+?\*|`[^`]+?`)")
_LIST_ITEM = re.compile(r"^[-*]\s+(.*)$")


@dataclass(frozen=True, slots=True)
class DocxMeta:
    """Everything the document shows that is not a table cell.

    ``prose`` is the rendered draft-Results markdown (already reference-resolved, so every
    number in it is still a table cell). ``provenance`` is label → value pairs, printed in
    order — the caller owns the wording because only it knows what the run was.
    """

    title: str
    label: str = ""  # footer, e.g. "NC_050833.1 · genome_stats"
    prose: str = ""
    provenance: Sequence[tuple[str, str]] = ()
    references: Sequence[str] = ()


def export_docx(result: Result, meta: DocxMeta | None = None) -> bytes:
    """Render a result as a Word document; deterministic for a given (result, meta)."""
    d = _load_docx()
    meta = meta or DocxMeta(title=result.module)
    document = d.Document()
    _apply_base_style(document)
    _apply_page(document.sections[0], landscape=False)
    _add_footer(document.sections[0], meta)

    document.core_properties.title = meta.title
    document.core_properties.author = "PlastidToolSuite"
    document.core_properties.subject = result.module
    document.core_properties.comments = meta.label
    document.core_properties.revision = 1

    document.add_heading(meta.title, level=0)
    _add_subtitle(document, f"{result.module} · module version {result.module_version}")
    if meta.provenance:
        _add_provenance(document, meta.provenance)
    if meta.prose.strip():
        add_markdown(document, meta.prose)
    _add_tables(document, result)
    if result.notes:
        document.add_heading("Notes", level=1)
        for note in result.notes:
            document.add_paragraph(note, style="List Bullet")
    if meta.references:
        document.add_heading("References", level=1)
        for reference in meta.references:
            document.add_paragraph(reference, style="List Bullet")

    buffer = io.BytesIO()
    document.save(buffer)
    return _deterministic_package(buffer.getvalue())


# --------------------------------------------------------------------------- document build


def merge_documents(
    documents: Sequence[bytes], *, title: str | None = None, label: str | None = None
) -> bytes:
    """Join documents this exporter built into one, section breaks and tables intact.

    A run's combined report is assembled from the modules' own documents rather than re-rendered
    from their results, because a cache hit restores files, not results — re-rendering would
    silently produce an empty document on the second run of an unchanged input.

    Documents from one template carry the same parts, so the package's plumbing comes from the
    first; the body is concatenated and the styles are unioned. A part list that differs is a
    programming error, not something to paper over, and raises.
    """
    if not documents:
        raise ValueError("nothing to merge")
    if len(documents) == 1:
        return _with_title(documents[0], title, label)

    from lxml import etree

    parts = [_package_parts(blob) for blob in documents]
    names = set(parts[0])
    for index, other in enumerate(parts[1:], start=2):
        if set(other) != names:
            raise ValueError(
                f"document {index} has a different part list from the first "
                f"({sorted(set(other) ^ names)}) — it was not built by this exporter"
            )

    merged = dict(parts[0])
    merged["word/document.xml"] = _merge_bodies(
        etree, [part["word/document.xml"] for part in parts]
    )
    merged["word/styles.xml"] = _merge_styles(etree, [part["word/styles.xml"] for part in parts])
    return _with_title(deterministic_zip(merged), title, label)


def _merge_bodies(etree: Any, blobs: Sequence[bytes]) -> bytes:
    """Concatenate the bodies, ending each document's last section where it ended.

    Each document's trailing ``w:sectPr`` holds that document's page setup. Dropping all but the
    last one would re-flow the earlier ones into the last one's orientation, which is how a
    landscape table comes back as a squeezed portrait one; instead each is re-attached as an
    inline section break, which is what OOXML uses to end a section with its own properties.
    """
    bodies = []
    for blob in blobs:
        nodes = list(etree.fromstring(blob).find(f"{W}body"))
        last = nodes.pop() if nodes and nodes[-1].tag == f"{W}sectPr" else None
        bodies.append((nodes, last))
    if bodies[-1][1] is None:
        raise ValueError("the last document has no closing sectPr — not a document we built")

    children: list[Any] = []
    for index, (nodes, last) in enumerate(bodies):
        children.extend(nodes)
        if index < len(bodies) - 1 and last is not None:
            # end this document's last section where it ended, before the next one starts
            children.append(_section_break(etree, last))
    children.append(bodies[-1][1])

    root = etree.fromstring(blobs[0])
    body = root.find(f"{W}body")
    for node in list(body):
        body.remove(node)
    for node in children:
        body.append(node)
    return etree.tostring(root, xml_declaration=True, encoding="UTF-8", standalone=True)


def _section_break(etree: Any, sect_pr: Any) -> Any:
    paragraph = etree.Element(f"{W}p")
    properties = etree.SubElement(paragraph, f"{W}pPr")
    properties.append(sect_pr)
    return paragraph


def _merge_styles(etree: Any, blobs: Sequence[bytes]) -> bytes:
    """Union of the referenced styles. One template means one definition per style id."""
    root = etree.fromstring(blobs[0])
    known = {style.get(f"{W}styleId") for style in root.findall(f"{W}style")}
    for blob in blobs[1:]:
        for style in etree.fromstring(blob).findall(f"{W}style"):
            if style.get(f"{W}styleId") not in known:
                root.append(style)
                known.add(style.get(f"{W}styleId"))
    return etree.tostring(root, xml_declaration=True, encoding="UTF-8", standalone=True)


def _package_parts(blob: bytes) -> dict[str, bytes]:
    return zip_parts(blob)


def _with_title(blob: bytes, title: str | None, label: str | None = None) -> bytes:
    """Rewrite the document title, and the footer's label when one is given.

    A merged document inherits the first module's footer, which would tell a reader on page 30
    that they are still reading `genome_stats`.
    """
    if not title and not label:
        return blob
    from lxml import etree

    parts = _package_parts(blob)
    if title:
        core = parts.get("docProps/core.xml")
        if core is not None:
            root = etree.fromstring(core)
            node = root.find("{http://purl.org/dc/elements/1.1/}title")
            if node is None:
                node = etree.SubElement(root, "{http://purl.org/dc/elements/1.1/}title")
            node.text = title
            parts["docProps/core.xml"] = etree.tostring(
                root, xml_declaration=True, encoding="UTF-8", standalone=True
            )
    if label:
        for name in [part for part in parts if part.startswith("word/footer")]:
            root = etree.fromstring(parts[name])
            first = next(iter(root.iter(f"{W}t")), None)
            if first is not None:
                first.text = f"PlastidToolSuite · {label} · Page "
            parts[name] = etree.tostring(
                root, xml_declaration=True, encoding="UTF-8", standalone=True
            )
    return deterministic_zip(parts)


def _load_docx() -> Any:
    return require_module("docx", extra="docx", purpose="Word export")


def _apply_base_style(document: _Document) -> None:
    """One font decision for the whole document, so every block inherits it."""
    normal = document.styles["Normal"]
    normal.font.name = BODY_FONT
    normal.font.size = _pt(BODY_PT)


def _pt(points: float) -> Any:
    from docx.shared import Pt

    return Pt(points)


def _apply_page(section: Any, *, landscape: bool) -> None:
    from docx.enum.section import WD_ORIENT
    from docx.shared import Mm

    long_side, short_side = Mm(PAGE_LONG_MM), Mm(PAGE_SHORT_MM)
    section.orientation = WD_ORIENT.LANDSCAPE if landscape else WD_ORIENT.PORTRAIT
    section.page_width, section.page_height = (
        (long_side, short_side) if landscape else (short_side, long_side)
    )
    section.top_margin = section.bottom_margin = Mm(MARGIN_EDGE_MM)
    section.left_margin = section.right_margin = Mm(MARGIN_SIDE_MM)


def _add_footer(section: Any, meta: DocxMeta) -> None:
    from docx.enum.text import WD_ALIGN_PARAGRAPH

    paragraph = section.footer.paragraphs[0]
    paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
    label = meta.label or meta.title
    run = paragraph.add_run(f"PlastidToolSuite · {label} · Page ")
    run.font.size = _pt(META_PT)
    _add_field(paragraph, " PAGE ", "1", size=META_PT)
    run = paragraph.add_run(" of ")
    run.font.size = _pt(META_PT)
    _add_field(paragraph, " NUMPAGES ", "1", size=META_PT)


def _add_field(paragraph: Paragraph, instruction: str, placeholder: str, *, size: float) -> None:
    """A Word field. Its value is computed by Word/LibreOffice, never by us."""
    from docx.oxml import OxmlElement
    from docx.oxml.ns import qn

    begin = paragraph.add_run()
    begin.font.size = _pt(size)
    element = OxmlElement("w:fldChar")
    element.set(qn("w:fldCharType"), "begin")
    begin._r.append(element)

    instruction_run = paragraph.add_run()
    instruction_run.font.size = _pt(size)
    text = OxmlElement("w:instrText")
    text.set(qn("xml:space"), "preserve")
    text.text = instruction
    instruction_run._r.append(text)

    separate = paragraph.add_run(placeholder)
    separate.font.size = _pt(size)

    end = paragraph.add_run()
    end.font.size = _pt(size)
    element = OxmlElement("w:fldChar")
    element.set(qn("w:fldCharType"), "end")
    end._r.append(element)


def _add_subtitle(document: _Document, text: str) -> None:
    from docx.shared import RGBColor

    paragraph = document.add_paragraph()
    run = paragraph.add_run(text)
    run.italic = True
    run.font.size = _pt(META_PT + 1)
    run.font.color.rgb = RGBColor.from_string("595959")


def _add_provenance(document: _Document, entries: Sequence[tuple[str, str]]) -> None:
    """A borderless label/value grid — provenance a reader can check, not a paragraph of prose."""
    from docx.shared import RGBColor

    widths = (34.0, PORTRAIT_TEXT_MM - 34.0)
    table = document.add_table(rows=0, cols=2)
    table.autofit = False
    _fixed_layout(table)
    _set_grid(table, widths)
    _table_cell_margins(table, CELL_MARGIN_MM)
    for label, value in entries:
        cells = table.add_row().cells
        run = cells[0].paragraphs[0].add_run(label)
        run.bold = True
        run.font.size = _pt(META_PT)
        run = cells[1].paragraphs[0].add_run(value)
        run.font.size = _pt(META_PT)
        run.font.color.rgb = RGBColor.from_string("404040")
        for width, cell in zip(widths, cells, strict=True):
            cell.width = _mm(width)
            _cell_margins(cell, top=0, bottom=0)


def _mm(millimetres: float) -> Any:
    from docx.shared import Mm

    return Mm(millimetres)


def _add_tables(document: _Document, result: Result) -> None:
    for index, table in enumerate(result.tables.values(), start=1):
        size = _table_font_pt(table)
        # Landscape only when the table actually needs the extra width. Deciding on the column
        # count instead put a seven-column, one-row summary on a page of its own, mostly empty.
        landscape = _natural_width_mm(table, size) > PORTRAIT_TEXT_MM
        _ensure_orientation(document, landscape=landscape)
        width_mm = LANDSCAPE_TEXT_MM if landscape else PORTRAIT_TEXT_MM
        _add_caption(document, index, table)
        _build_table(document, table, width_mm=width_mm, size=size)


def _table_font_pt(table: Table) -> float:
    columns = len(table.columns)
    if columns <= 5:
        return TABLE_PT
    return WIDE_TABLE_PT if columns < 12 else NARROW_TABLE_PT


def _natural_width_mm(table: Table, size_pt: float) -> float:
    """How wide the table wants to be: its widest cell per column, plus cell padding."""
    per_character = size_pt * MM_PER_PT * CHAR_WIDTH_EM
    return sum(
        min(_column_weight(table, column), 60) * per_character + CELL_PADDING_MM
        for column in table.columns
    )


def _column_weight(table: Table, column: Any) -> float:
    """A column's width in character units — header (bolded, so boosted), cells, and a floor.

    The floor is what keeps a short header from being squeezed until it wraps: "Genes" over
    three-character values is still a word that has to fit.
    """
    widest = max(len(column.format(row.get(column.name))) for row in table.rows)
    return max(len(column.header) * HEADER_BOOST, widest, MIN_COLUMN_CHARS)


def _ensure_orientation(document: _Document, *, landscape: bool) -> None:
    from docx.enum.section import WD_ORIENT, WD_SECTION

    current = document.sections[-1]
    if (current.orientation == WD_ORIENT.LANDSCAPE) == landscape:
        return
    section = document.add_section(WD_SECTION.NEW_PAGE)
    _apply_page(section, landscape=landscape)
    section.footer.is_linked_to_previous = True


def _add_caption(document: _Document, index: int, table: Table) -> None:
    from docx.shared import RGBColor

    paragraph = document.add_paragraph()
    paragraph.paragraph_format.space_before = _pt(12)
    paragraph.paragraph_format.space_after = _pt(3)
    paragraph.paragraph_format.keep_with_next = True
    run = paragraph.add_run(f"Table {index}. {table.title}")
    run.bold = True
    run.font.size = _pt(CAPTION_PT)

    scope = document.add_paragraph()
    scope.paragraph_format.space_after = _pt(4)
    scope.paragraph_format.keep_with_next = True
    run = scope.add_run(f"scope: {table.scope.label} | coordinates: {table.scope.coord_system}")
    run.italic = True
    run.font.size = _pt(META_PT)
    run.font.color.rgb = RGBColor.from_string("595959")


def _build_table(document: _Document, table: Table, *, width_mm: float, size: float) -> None:
    """The data table: fixed layout, widths from the content, header row repeated per page."""
    from docx.enum.table import WD_TABLE_ALIGNMENT

    columns = table.columns
    widths = _column_widths(table, width_mm, size)
    banded = len(table.rows) >= BANDED_FROM_ROWS

    grid = document.add_table(rows=0, cols=len(columns))
    grid.alignment = WD_TABLE_ALIGNMENT.LEFT
    grid.autofit = False
    _fixed_layout(grid)
    _set_grid(grid, widths)
    _table_cell_margins(grid, CELL_MARGIN_MM)
    _table_borders(grid)

    header = grid.add_row()
    _repeat_as_header(header)
    _keep_together(header)
    for cell, column, cell_width in zip(header.cells, columns, widths, strict=True):
        cell.width = _mm(cell_width)
        _fill_cell(cell, column.header, size=size, bold=True, align=column.dtype)
        _shade(cell, HEADER_FILL)
        _cell_bottom_border(cell, HEADER_RULE_COLOUR, size=8)

    emphasised = set(table.emphasis)
    for row_index, row in enumerate(table.rows):
        word_row = grid.add_row()
        _keep_together(word_row)
        bold = table.key is not None and row.get(table.key) in emphasised
        for cell, column, cell_width in zip(word_row.cells, columns, widths, strict=True):
            cell.width = _mm(cell_width)
            _fill_cell(
                cell,
                column.format(row.get(column.name)),
                size=size,
                bold=bold,
                align=column.dtype,
            )
            if banded and row_index % 2 == 1:
                _shade(cell, BAND_FILL)

    for cell in grid.rows[-1].cells:  # the closing rule
        _cell_bottom_border(cell, RULE_COLOUR, size=8)

    _caption_note(document, table)


def _caption_note(document: _Document, table: Table) -> None:
    from docx.shared import RGBColor

    paragraph = document.add_paragraph()
    paragraph.paragraph_format.space_before = _pt(3)
    run = paragraph.add_run(f"{len(table.rows)} rows")
    run.italic = True
    run.font.size = _pt(META_PT)
    run.font.color.rgb = RGBColor.from_string("808080")


def _column_widths(table: Table, total_mm: float, size_pt: float) -> list[float]:
    """Every column gets its minimum first; the rest is shared out by how much text it holds.

    Scaling all columns by one factor instead starves the short ones — a wide table has one
    long header ("gene features (diagnostic only)") that eats the budget, and a column holding
    one word ends up too narrow to print it, which Word answers by breaking the word.
    """
    floors = [_min_column_mm(size_pt) for _ in table.columns]
    spare = total_mm - sum(floors)
    if spare <= 0:  # nothing to share out; the table is wider than the page whatever we do
        return floors
    weights = [min(_column_weight(table, column), 60) for column in table.columns]
    total_weight = float(sum(weights))
    return [
        floor + spare * weight / total_weight for floor, weight in zip(floors, weights, strict=True)
    ]


def _min_column_mm(size_pt: float) -> float:
    return MIN_COLUMN_CHARS * size_pt * MM_PER_PT * CHAR_WIDTH_EM + CELL_PADDING_MM


def _fixed_layout(table: Any) -> None:
    from docx.oxml import OxmlElement
    from docx.oxml.ns import qn

    element = OxmlElement("w:tblLayout")
    element.set(qn("w:type"), "fixed")
    table._tbl.tblPr.append(element)


def _set_grid(table: Any, widths_mm: Sequence[float]) -> None:
    """Write the widths into ``w:tblGrid``; Word honours the grid, not just the cell widths."""
    from docx.oxml.ns import qn

    grid = table._tbl.find(qn("w:tblGrid"))
    if grid is None:
        return
    for column, width in zip(grid.findall(qn("w:gridCol")), widths_mm, strict=True):
        column.set(qn("w:w"), str(round(width * TWIPS_PER_MM)))


def _repeat_as_header(row: Any) -> None:
    from docx.oxml import OxmlElement
    from docx.oxml.ns import qn

    element = OxmlElement("w:tblHeader")
    element.set(qn("w:val"), "true")
    row._tr.get_or_add_trPr().append(element)


def _keep_together(row: Any) -> None:
    from docx.oxml import OxmlElement

    row._tr.get_or_add_trPr().append(OxmlElement("w:cantSplit"))


def _shade(cell: Any, fill: str) -> None:
    from docx.oxml import OxmlElement
    from docx.oxml.ns import qn

    element = OxmlElement("w:shd")
    element.set(qn("w:val"), "clear")
    element.set(qn("w:color"), "auto")
    element.set(qn("w:fill"), fill)
    cell._tc.get_or_add_tcPr().append(element)


def _cell_bottom_border(cell: Any, colour: str, *, size: int) -> None:
    from docx.oxml import OxmlElement
    from docx.oxml.ns import qn

    borders = OxmlElement("w:tcBorders")
    bottom = OxmlElement("w:bottom")
    bottom.set(qn("w:val"), "single")
    bottom.set(qn("w:sz"), str(size))
    bottom.set(qn("w:space"), "0")
    bottom.set(qn("w:color"), colour)
    borders.append(bottom)
    cell._tc.get_or_add_tcPr().append(borders)


def _table_borders(table: Any) -> None:
    """Horizontal rules only: the eye tracks a wide table along its rows, not its gutters."""
    from docx.oxml import OxmlElement
    from docx.oxml.ns import qn

    borders = OxmlElement("w:tblBorders")
    for edge, colour, size in (
        ("top", HEADER_RULE_COLOUR, 8),
        ("bottom", RULE_COLOUR, 8),
        ("insideH", RULE_COLOUR, 4),
    ):
        element = OxmlElement(f"w:{edge}")
        element.set(qn("w:val"), "single")
        element.set(qn("w:sz"), str(size))
        element.set(qn("w:space"), "0")
        element.set(qn("w:color"), colour)
        borders.append(element)
    for edge in ("left", "right", "insideV"):
        element = OxmlElement(f"w:{edge}")
        element.set(qn("w:val"), "none")
        element.set(qn("w:sz"), "0")
        borders.append(element)
    table._tbl.tblPr.append(borders)


def _table_cell_margins(table: Any, millimetres: float) -> None:
    """Set the cell padding explicitly — the width model above has to know the real number."""
    from docx.oxml import OxmlElement
    from docx.oxml.ns import qn

    width = str(round(millimetres * TWIPS_PER_MM))
    margins = OxmlElement("w:tblCellMar")
    for edge in ("top", "left", "bottom", "right"):
        element = OxmlElement(f"w:{edge}")
        element.set(qn("w:w"), "0" if edge in ("top", "bottom") else width)
        element.set(qn("w:type"), "dxa")
        margins.append(element)
    table._tbl.tblPr.append(margins)


def _cell_margins(cell: Any, *, top: int, bottom: int) -> None:
    from docx.oxml import OxmlElement
    from docx.oxml.ns import qn

    margins = OxmlElement("w:tcMar")
    for edge, value in (("top", top), ("bottom", bottom)):
        element = OxmlElement(f"w:{edge}")
        element.set(qn("w:w"), str(value))
        element.set(qn("w:type"), "dxa")
        margins.append(element)
    cell._tc.get_or_add_tcPr().append(margins)


def _fill_cell(cell: Any, text: str, *, size: float, bold: bool = False, align: str) -> None:
    from docx.enum.text import WD_ALIGN_PARAGRAPH

    paragraph = cell.paragraphs[0]
    paragraph.paragraph_format.space_before = _pt(0)
    paragraph.paragraph_format.space_after = _pt(0)
    paragraph.alignment = (
        WD_ALIGN_PARAGRAPH.RIGHT if align in ("int", "float") else WD_ALIGN_PARAGRAPH.LEFT
    )
    run = paragraph.add_run(text)
    run.bold = bold
    run.font.size = _pt(size)


# ------------------------------------------------------------------------ markdown subset


def add_markdown(document: Any, text: str) -> None:
    """The report's markdown subset, as Word blocks: headings, prose, bullets, pipe tables.

    Paragraphs are joined across their hard wraps — the templates wrap at 100 columns, and a
    soft wrap must not become a new paragraph.
    """
    for block in _blocks(text):
        kind = block[0]
        if kind == "heading":
            document.add_heading(block[1], level=block[2])
        elif kind == "bullets":
            for item in block[1]:
                _add_inline(document.add_paragraph(style="List Bullet"), item)
        elif kind == "table":
            _add_markdown_table(document, block[1])
        else:
            _add_inline(document.add_paragraph(), block[1])


def _blocks(text: str) -> list[tuple]:
    blocks: list[tuple] = []
    lines = text.splitlines()
    index = 0
    while index < len(lines):
        line = lines[index].rstrip()
        if not line.strip():
            index += 1
            continue
        heading = re.match(r"^(#{1,6})\s+(.*)$", line)
        if heading:
            # The document title is Word's level 0, so `##` is this document's first level.
            level = max(1, min(len(heading.group(1)) - 1, 4))
            blocks.append(("heading", heading.group(2).strip(), level))
            index += 1
            continue
        if line.lstrip().startswith("|"):
            rows = []
            while index < len(lines) and lines[index].lstrip().startswith("|"):
                cells = [cell.strip() for cell in lines[index].strip().strip("|").split("|")]
                if not all(re.fullmatch(r"[-: ]+", cell) for cell in cells if cell):
                    rows.append(cells)
                index += 1
            blocks.append(("table", rows))
            continue
        item = _LIST_ITEM.match(line.strip())
        if item:
            items = []
            while index < len(lines):
                match = _LIST_ITEM.match(lines[index].strip())
                if not match:
                    break
                items.append(match.group(1).strip())
                index += 1
            blocks.append(("bullets", items))
            continue
        paragraph = [line.strip()]
        index += 1
        while index < len(lines) and lines[index].strip() and not _starts_a_block(lines[index]):
            paragraph.append(lines[index].strip())
            index += 1
        blocks.append(("paragraph", " ".join(paragraph)))
    return blocks


def _starts_a_block(line: str) -> bool:
    stripped = line.lstrip()
    return stripped.startswith(("#", "|")) or _LIST_ITEM.match(stripped) is not None


def _add_inline(paragraph: Any, text: str) -> None:
    """``**bold**`` · ``*italic*`` · `` `code` `` — the three spans the templates use."""
    for token in _INLINE.split(text):
        if not token:
            continue
        if token.startswith("**") and token.endswith("**") and len(token) > 4:
            paragraph.add_run(token[2:-2]).bold = True
        elif token.startswith("`") and token.endswith("`") and len(token) > 2:
            run = paragraph.add_run(token[1:-1])
            run.font.name = MONO_FONT
            run.font.size = _pt(BODY_PT - 1.5)
        elif token.startswith("*") and token.endswith("*") and len(token) > 2:
            paragraph.add_run(token[1:-1]).italic = True
        else:
            paragraph.add_run(token)


def _add_markdown_table(document: Any, rows: list[list[str]]) -> None:
    if not rows:
        return
    width = max(len(row) for row in rows)
    grid = document.add_table(rows=0, cols=width)
    grid.style = "Table Grid"
    grid.autofit = True
    for row_index, row in enumerate(rows):
        cells = grid.add_row().cells
        for column_index, cell in enumerate(cells):
            value = row[column_index] if column_index < len(row) else ""
            _add_inline(cell.paragraphs[0], value)
            if row_index == 0:
                for run in cell.paragraphs[0].runs:
                    run.bold = True


# ------------------------------------------------------------------ deterministic packaging


def _deterministic_package(data: bytes) -> bytes:
    """Rebuild the ``.docx`` zip: template parts pruned, entries sorted, timestamps fixed.

    ``python-docx`` stamps each entry with the wall-clock time, so two identical builds differ
    at the byte level. Everything below is what makes the export reproducible.
    """
    return deterministic_zip(_pruned_parts(data))


def _pruned_parts(data: bytes) -> dict[str, bytes]:
    parts = zip_parts(data)

    for name in UNUSED_PARTS:
        parts.pop(name, None)
    for name in [name for name in parts if name.startswith(UNUSED_PREFIXES)]:
        parts.pop(name, None)

    kept = set(parts)
    for name in [name for name in parts if name.endswith(".rels")]:
        parts[name] = _prune_relationships(name, parts[name], kept)
    if "word/styles.xml" in parts:
        document_part = parts.get("word/document.xml", b"")
        parts["word/styles.xml"] = _prune_styles(parts["word/styles.xml"], document_part)
    if "[Content_Types].xml" in parts:
        parts["[Content_Types].xml"] = _prune_content_types(parts["[Content_Types].xml"], kept)
    return parts


def _prune_relationships(rel_name: str, blob: bytes, kept: set[str]) -> bytes:
    """Drop any relationship whose target is no longer in the package.

    A relationship pointing at a missing part makes the file unopenable, so this runs after
    every prune rather than being maintained by hand.
    """
    from lxml import etree

    root = etree.fromstring(blob)
    base = rel_name.rsplit("_rels/", 1)[0].rstrip("/")
    for relationship in list(root):
        if relationship.get("TargetMode") == "External":
            continue
        target = relationship.get("Target", "")
        if _resolve_target(base, target) not in kept:
            root.remove(relationship)
    return etree.tostring(root, xml_declaration=True, encoding="UTF-8", standalone=True)


def _resolve_target(base: str, target: str) -> str:
    if target.startswith("/"):
        return target.lstrip("/")
    stack = [segment for segment in base.split("/") if segment]
    for segment in target.split("/"):
        if segment == "..":
            if stack:
                stack.pop()
        elif segment not in ("", "."):
            stack.append(segment)
    return "/".join(stack)


def _prune_styles(styles_blob: bytes, document_blob: bytes) -> bytes:
    """Keep the styles the document uses, their ancestors, and the defaults.

    The template carries Word's whole built-in catalogue — 164 styles, ~350 kB, of which a
    generated document references about ten. Word only needs the ones it is asked for plus
    whatever they inherit from; a style that is missing but never referenced is invisible.
    """
    from lxml import etree

    root = etree.fromstring(styles_blob)
    styles = {style.get(f"{W}styleId"): style for style in root.findall(f"{W}style")}

    # A style is referenced by a *child element* — `<w:pStyle w:val="Heading1"/>`. Reading it as
    # an attribute finds nothing, and the first version of this pruned every style the document
    # actually used; the render survived only because Word falls back to its own defaults.
    wanted: list[str] = []
    for node in etree.fromstring(document_blob).iter(f"{W}pStyle", f"{W}rStyle", f"{W}tblStyle"):
        value = node.get(f"{W}val")
        if value:
            wanted.append(value)
    wanted.extend(
        style.get(f"{W}styleId") for style in styles.values() if style.get(f"{W}default") == "1"
    )

    keep: set[str] = set()
    while wanted:
        style_id = wanted.pop()
        if style_id in keep or style_id not in styles:
            continue
        keep.add(style_id)
        for tag in ("basedOn", "link", "next"):
            reference = styles[style_id].find(f"{W}{tag}")
            if reference is not None:
                wanted.append(reference.get(f"{W}val"))
    for style_id, style in styles.items():
        if style_id not in keep:
            root.remove(style)
    return etree.tostring(root, xml_declaration=True, encoding="UTF-8", standalone=True)


def _prune_content_types(blob: bytes, kept: set[str]) -> bytes:
    from lxml import etree

    root = etree.fromstring(blob)
    for override in list(root):
        if override.tag.endswith("Override"):
            part = (override.get("PartName") or "").lstrip("/")
            if part and part not in kept:
                root.remove(override)
    return etree.tostring(root, xml_declaration=True, encoding="UTF-8", standalone=True)
