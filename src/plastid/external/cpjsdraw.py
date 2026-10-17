"""CPJSdraw adapter — the IR-junction figure (docs/tools/cpjsdraw.md).

**CPJSdraw is not bundled.** Measured 17/10/2026: the upstream repository carries no licence
(GitHub's API reports ``license: null`` for ``xul962464/CPJSdraw``) and its ``readme.txt`` states
none — so, by the rule that took MISA out of the package (P9-01), it is declared and the user
obtains it. A tool that is only *called* is not redistributed; a missing licence forbids the other
thing. The bundle this adapter ports says "CPJSdraw is MIT" in its integration README and "no
licence file" in its own tool note; the second is what the tree shows.

**Two patches, and the run checks for them.** The manuscript's junction figure needs a bolder face
and, more importantly, an IR copy that wraps the genome origin: CPJSdraw 0.0.1 reads
``IRa:(\\d+)-(\\d+)`` only, so such a copy is drawn short by however many bases it wraps — measured
1, 4, 7, 36 and 42 bp across the five records this project drew. The patched script is identified
by its md5, so a pristine clone stops the run **naming the patch** instead of producing a figure
that is quietly a few bases wrong.

**Nothing is probed with a flag.** ``CPJSdraw.pl`` is a Perl script: handed an unknown flag it runs
itself with its default arguments, and PGA's next act was ``rm -rf gb`` (``run/tools.py``). Its
version is read out of its own header; the interpreter's version is asked of Perl.
"""

from __future__ import annotations

import hashlib
import re
from pathlib import Path

from plastid.run.context import RunContext
from plastid.run.tools import ToolError, tool_version, which

__all__ = [
    "CPJSDRAW_VERSION",
    "PATCHED_MD5",
    "cpjsdraw_script",
    "cpjsdraw_version",
    "patched_script",
    "run_cpjsdraw",
    "svg_to_png",
]

CPJSDRAW_VERSION = "0.0.1"
CPJSDRAW_DOC = "docs/tools/cpjsdraw.md"
_SCRIPT_NAME = "CPJSdraw.pl"

# The build the manuscript's Figure 2 was drawn with: the upstream script plus this project's two
# patches. Both the md5 and the byte count are recorded, because a hash alone gives no way to tell
# a truncated download from a patched file.
PATCHED_MD5 = "ba41be6002bfd0b797907d3290b34419"
PATCHED_BYTES = 52_947

# ``my $version = "0.0.1";`` — the script's only statement of its own release. The help text
# prints ``Version: $version``, so reading the printed form would record a variable name. It has
# no version flag to ask instead.
_VERSION_RE = re.compile(r'my\s+\$version\s*=\s*"([0-9][0-9.]*)"')
_HEADER_CHARS = 4_000


def cpjsdraw_script(explicit: str | Path | None = None) -> Path:
    """The ``CPJSdraw.pl`` to run: an explicit path, else the name on ``PATH``.

    The message names the one thing a user can act on — the clone line in the doc — rather than
    reporting a path that was never going to exist.
    """
    if explicit is not None:
        path = Path(explicit)
        if not path.is_file():
            raise ToolError(f"cpjsdraw: --set cpjsdraw={explicit} is not a file — {CPJSDRAW_DOC}")
        return path
    found = which(_SCRIPT_NAME)
    if found is None:
        raise ToolError(
            f"cpjsdraw: {_SCRIPT_NAME} not found — CPJSdraw is not bundled with this toolkit, "
            f"because its source states no licence. Clone it and apply the two patches with the "
            f"lines in {CPJSDRAW_DOC}, put it on PATH, or pass "
            f"--set cpjsdraw=/path/to/{_SCRIPT_NAME}"
        )
    return Path(found)


def cpjsdraw_version(script: Path) -> str:
    """The release this script is, read from its own header — never by passing it a flag."""
    header = script.read_text(encoding="utf-8", errors="replace")[:_HEADER_CHARS]
    match = _VERSION_RE.search(header)
    if match is None:
        raise ToolError(
            f"{script}: no 'my $version = \"x.y\";' in its first lines, so the CPJSdraw release "
            f"cannot be read — install {CPJSDRAW_VERSION} ({CPJSDRAW_DOC})"
        )
    version = match.group(1)
    if version != CPJSDRAW_VERSION:
        raise ToolError(
            f"{script}: expected CPJSdraw {CPJSDRAW_VERSION}, got {version} — the figure was "
            f"measured against {CPJSDRAW_VERSION} ({CPJSDRAW_DOC})"
        )
    return version


def patched_script(script: Path) -> Path:
    """Refuse anything that is not the patched build the junction figure was measured with.

    A pristine clone is not "close enough": without the compound-IRa patch an IR copy that wraps
    the origin is drawn short, and the shortfall is silent — the figure still looks like a figure.
    The check is therefore a refusal, not a warning.
    """
    digest = hashlib.md5(script.read_bytes()).hexdigest()  # noqa: S324 - identifying a file, not securing one
    if digest == PATCHED_MD5:
        return script
    size = script.stat().st_size
    raise ToolError(
        f"{script}: this is not the patched CPJSdraw the junction figure was measured with "
        f"(md5 {digest}, {size} B; expected {PATCHED_MD5}, {PATCHED_BYTES} B). A pristine clone "
        f"draws an IR copy that wraps the genome origin up to 42 bp short and says nothing. Apply "
        f"both patches — the hunks are in {CPJSDRAW_DOC} — then re-run"
    )


def run_cpjsdraw(
    ctx: RunContext,
    cfg: Path,
    out_svg: Path,
    *,
    script: str | Path | None = None,
    perl: str | None = None,
    drop_pseudo: bool = False,
    timeout: int | None = None,
) -> Path:
    """Run CPJSdraw over one configuration file and return the SVG it wrote.

    The configuration holds absolute paths to the records, and the run happens in the run's own
    directory, so nothing is written beside the input.
    """
    cfg, out_svg = Path(cfg).resolve(), Path(out_svg).resolve()
    if not cfg.is_file():
        raise ToolError(f"cpjsdraw: no configuration file at {cfg}")
    out_svg.parent.mkdir(parents=True, exist_ok=True)

    perl_binary = perl or which("perl")
    if perl_binary is None:
        raise ToolError("cpjsdraw: perl not found — CPJSdraw is a Perl program; install perl")
    script_path = patched_script(cpjsdraw_script(script))
    version = cpjsdraw_version(script_path)
    interpreter_version, _ = tool_version("perl", "perl", perl_binary)

    argv = [perl_binary, str(script_path), "-i", str(cfg), "-o", str(out_svg)]
    if drop_pseudo:
        argv.append("-p")
    ctx.tool(
        "cpjsdraw",
        argv,
        kind="perl",
        tool_path=str(script_path),
        pinned_version=version,
        interpreter=f"perl {interpreter_version}" if interpreter_version else "perl",
        cwd=out_svg.parent,
        timeout=timeout,
    )
    if not out_svg.is_file():
        raise ToolError(
            f"cpjsdraw: exited cleanly but wrote no SVG at {out_svg} — see the captured logs"
        )
    return out_svg


def svg_to_png(
    ctx: RunContext,
    svg: Path,
    png: Path,
    *,
    width: int = 3200,
    timeout: int | None = None,
) -> Path:
    """Rasterise the SVG with ``rsvg-convert``.

    Not CPJSdraw's own converter: it calls the ``batik-1.7`` it ships, which needs the Security
    Manager that Java 17 removed — measured in the bundle, and the reason the manuscript's PNGs
    were made with ``rsvg-convert`` instead.
    """
    converter = which("rsvg-convert")
    if converter is None:
        raise ToolError(
            "rsvg-convert not found — CPJSdraw's own SVG->PNG step needs Java 17's removed "
            "Security Manager, so this toolkit rasterises with librsvg: apt install librsvg2-bin "
            f"({CPJSDRAW_DOC})"
        )
    png = Path(png).resolve()
    png.parent.mkdir(parents=True, exist_ok=True)
    ctx.tool(
        "rsvg-convert",
        [converter, "-w", str(width), str(Path(svg).resolve()), "-o", str(png)],
        kind="binary",
        tool_path=converter,
        cwd=png.parent,
        timeout=timeout,
    )
    if not png.is_file():
        raise ToolError(f"rsvg-convert: exited cleanly but wrote no PNG at {png}")
    return png
