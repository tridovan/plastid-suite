"""The local interface: a small server on this computer, for one person, and nothing else.

Why a server on 127.0.0.1 is still guarded. A browser will send requests to it from any page that
person has open, so reaching the port proves nothing. Three checks hold together:

* a **session token**, printed when the server starts and traded once for a cookie, so a page that
  does not hold it cannot read or start anything;
* the **Host** header must name this server, which stops a hostile name rebound onto 127.0.0.1;
* the **Origin** header must be this server on every request that changes something.

The token leaves the address bar at once, and no request is logged, so it never reaches a log or a
Referer. Everything the page loads comes from this package: nothing is fetched from elsewhere, and
the interface works with no network.

Standard library only, so the minimum install stays biopython + PyYAML.
"""

from __future__ import annotations

import errno
import hmac
import json
import mimetypes
import re
import secrets
from collections.abc import Callable
from dataclasses import dataclass, field
from http import HTTPStatus
from http.cookies import SimpleCookie
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from importlib import resources
from pathlib import Path
from typing import Any
from urllib.parse import parse_qs, urlsplit

from plastid.ui.app import App, NotFound, Refused
from plastid.ui.catalogue import LOOPBACK

__all__ = ["LOOPBACK", "PortBusy", "Session", "serve"]

_COOKIE = "pts_session"
_MAX_JSON = 1_000_000  # bytes; a form is small, and a file has a route of its own
_CHUNK = 1 << 20  # 1 MiB, for files sent to the browser and received from it
_HEADERS = (
    (
        "Content-Security-Policy",
        (
            "default-src 'none'; script-src 'self'; style-src 'self'; img-src 'self' data:; "
            "font-src 'self'; connect-src 'self'; base-uri 'none'; form-action 'none'; "
            "frame-ancestors 'none'"
        ),
    ),
    ("Referrer-Policy", "no-referrer"),
    ("X-Content-Type-Options", "nosniff"),
    ("X-Frame-Options", "DENY"),
    ("Cache-Control", "no-store"),
)


class PortBusy(OSError):
    """The port is taken. The server stops and says so; it never moves to another one."""


@dataclass(slots=True)
class Session:
    """The one secret of a server's life, and who may therefore talk to it."""

    port: int
    token: str = field(default_factory=lambda: secrets.token_urlsafe(32))

    @property
    def hosts(self) -> frozenset[str]:
        return frozenset(f"{name}:{self.port}" for name in (LOOPBACK, "localhost", "[::1]"))

    @property
    def origins(self) -> frozenset[str]:
        return frozenset(f"http://{host}" for host in self.hosts)

    def holds(self, value: str | None) -> bool:
        return value is not None and hmac.compare_digest(value, self.token)

    def link(self) -> str:
        return f"http://{LOOPBACK}:{self.port}/?token={self.token}"


Route = Callable[..., Any]


def _routes(app: App) -> list[tuple[str, re.Pattern[str], Route]]:
    """Every address the server answers, in one table."""
    table: list[tuple[str, str, Route]] = [
        ("GET", r"/api/app", app.about),
        ("GET", r"/api/modules", app.modules),
        ("GET", r"/api/tools", app.tools),
        ("POST", r"/api/tools/refresh", app.refresh_tools),
        ("GET", r"/api/records", app.records),
        ("GET", r"/api/record", app.record),
        ("GET", r"/api/reference-set", app.reference_set),
        ("POST", r"/api/records/adopt", app.adopt),
        ("GET", r"/api/browse", app.browse),
        ("GET", r"/api/runs", app.runs),
        ("GET", r"/api/run/(?P<dataset>[^/]+)/(?P<run_id>[^/]+)", app.run),
        ("POST", r"/api/run/(?P<dataset>[^/]+)/(?P<run_id>[^/]+)/cancel", app.cancel_run),
        ("GET", r"/api/figures", app.figures),
        ("GET", r"/api/table", app.table),
        ("GET", r"/api/text", app.text),
        ("POST", r"/api/command", app.command),
        ("GET", r"/api/jobs", app.jobs),
        ("POST", r"/api/jobs", app.start_job),
        ("GET", r"/api/jobs/(?P<job_id>[0-9a-f]+)", app.job),
    ]
    return [(method, re.compile(pattern + r"\Z"), route) for method, pattern, route in table]


class _Handler(BaseHTTPRequestHandler):
    server_version = "plastid"
    sys_version = ""
    protocol_version = "HTTP/1.1"
    session: Session
    app: App
    routes: list[tuple[str, re.Pattern[str], Route]]

    def log_message(self, format: str, *args: Any) -> None:  # noqa: A002
        """Say nothing: a request line can carry the token, and a log outlives the session."""

    def _plain(self, status: HTTPStatus, words: str) -> None:
        body = words.encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "text/plain; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        for name, value in _HEADERS:
            self.send_header(name, value)
        self.end_headers()
        self.wfile.write(body)

    def _json(self, payload: Any, status: HTTPStatus = HTTPStatus.OK) -> None:
        body = json.dumps(payload, ensure_ascii=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        for name, value in _HEADERS:
            self.send_header(name, value)
        self.end_headers()
        self.wfile.write(body)

    def _guarded(self) -> bool:
        """Whether this request came from the page this server gave out."""
        if self.headers.get("Host") not in self.session.hosts:
            self._plain(HTTPStatus.FORBIDDEN, "This server answers only on its own address.\n")
            return False
        origin = self.headers.get("Origin")
        if origin is not None and origin not in self.session.origins:
            self._plain(HTTPStatus.FORBIDDEN, "This server answers only its own page.\n")
            return False
        if self.command == "POST" and origin is None:
            self._plain(HTTPStatus.FORBIDDEN, "A change must come from this server's page.\n")
            return False
        return True

    def _has_session(self) -> bool:
        jar: SimpleCookie = SimpleCookie(self.headers.get("Cookie", ""))
        morsel = jar.get(_COOKIE)
        return self.session.holds(morsel.value if morsel else None)

    def _enter(self, query: dict[str, list[str]]) -> bool:
        """Trade the token in the link for a cookie, and take it out of the address bar."""
        offered = query.get("token", [None])[0]
        if not self.session.holds(offered):
            return False
        self.send_response(HTTPStatus.SEE_OTHER)
        self.send_header("Location", "/")
        self.send_header(
            "Set-Cookie", f"{_COOKIE}={self.session.token}; HttpOnly; SameSite=Strict; Path=/"
        )
        self.send_header("Content-Length", "0")
        for name, value in _HEADERS:
            self.send_header(name, value)
        self.end_headers()
        return True

    def _dispatch(self, method: str) -> None:
        if not self._guarded():
            return
        parts = urlsplit(self.path)
        query = parse_qs(parts.query)
        if method == "GET" and parts.path == "/" and "token" in query and self._enter(query):
            return
        if not self._has_session():
            self._plain(
                HTTPStatus.FORBIDDEN,
                "Open the link `plastid ui` printed when it started. It holds this session's "
                "token.\n",
            )
            return
        if method == "POST" and parts.path == "/api/upload":
            self._upload(query)
            return
        if method == "GET" and parts.path == "/api/file":
            self._file(query)
            return
        if method == "GET" and parts.path == "/api/preview":
            self._preview(query)
            return
        for wanted, pattern, route in self.routes:
            found = pattern.match(parts.path) if wanted == method else None
            if found:
                self._answer(route, found.groupdict(), query, method)
                return
        if method == "GET":
            self._asset(parts.path)
            return
        self._plain(HTTPStatus.NOT_FOUND, "Nothing here.\n")

    def _answer(
        self, route: Route, named: dict[str, str], query: dict[str, list[str]], method: str
    ) -> None:
        try:
            if method == "POST":
                body = self._form()
                self._json(route(body, **named))
            else:
                flat = {key: values[0] for key, values in query.items()}
                self._json(route(**named, **flat))
        except NotFound as error:
            self._json({"error": str(error)}, HTTPStatus.NOT_FOUND)
        except Refused as error:
            self._json({"error": str(error)}, HTTPStatus.BAD_REQUEST)
        except TypeError:
            self._json({"error": "that address takes other parameters"}, HTTPStatus.BAD_REQUEST)

    def _form(self) -> dict[str, Any]:
        length = int(self.headers.get("Content-Length") or 0)
        if length > _MAX_JSON:
            raise Refused("that form is too large")
        kind = self.headers.get("Content-Type", "").split(";")[0].strip()
        if length and kind != "application/json":
            raise Refused("a form is sent as JSON")
        raw = self.rfile.read(length) if length else b"{}"
        try:
            payload = json.loads(raw)
        except ValueError:
            raise Refused("that form could not be read") from None
        if not isinstance(payload, dict):
            raise Refused("that form could not be read")
        return payload

    def _upload(self, query: dict[str, list[str]]) -> None:
        """Copy a dropped file into the workspace. Nothing leaves this computer."""
        try:
            length = int(self.headers.get("Content-Length") or 0)
            target = self.app.upload_target(query.get("name", [""])[0], query.get("to", [""])[0])
            remaining = length
            with target.open("wb") as handle:
                while remaining > 0:
                    chunk = self.rfile.read(min(_CHUNK, remaining))
                    if not chunk:
                        break
                    handle.write(chunk)
                    remaining -= len(chunk)
            self._json({"path": self.app.relative(target)})
        except Refused as error:
            self._json({"error": str(error)}, HTTPStatus.BAD_REQUEST)

    def _send_path(self, path: Path, *, download: bool = False) -> None:
        kind = mimetypes.guess_type(path.name)[0] or "application/octet-stream"
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", kind)
        self.send_header("Content-Length", str(path.stat().st_size))
        if download:
            self.send_header("Content-Disposition", f'attachment; filename="{path.name}"')
        for name, value in _HEADERS:
            self.send_header(name, value)
        self.end_headers()
        with path.open("rb") as handle:
            while chunk := handle.read(_CHUNK):
                self.wfile.write(chunk)

    def _file(self, query: dict[str, list[str]]) -> None:
        try:
            path = self.app.readable(query.get("path", [""])[0])
        except (NotFound, Refused) as error:
            self._plain(HTTPStatus.NOT_FOUND, f"{error}\n")
            return
        self._send_path(path, download="download" in query)

    def _preview(self, query: dict[str, list[str]]) -> None:
        try:
            path = self.app.preview(query.get("path", [""])[0])
        except (NotFound, Refused) as error:
            self._plain(HTTPStatus.NOT_FOUND, f"{error}\n")
            return
        self._send_path(path)

    def _asset(self, path: str) -> None:
        name = "index.html" if path == "/" else path.lstrip("/")
        if not re.fullmatch(r"[A-Za-z0-9_-]+(\.[A-Za-z0-9]+)+|fonts/[A-Za-z0-9_.-]+", name):
            self._plain(HTTPStatus.NOT_FOUND, "Nothing here.\n")
            return
        asset = resources.files("plastid.ui") / "assets" / name
        if not asset.is_file():
            self._plain(HTTPStatus.NOT_FOUND, "Nothing here.\n")
            return
        body = asset.read_bytes()
        self.send_response(HTTPStatus.OK)
        self.send_header(
            "Content-Type", mimetypes.guess_type(name)[0] or "application/octet-stream"
        )
        self.send_header("Content-Length", str(len(body)))
        for header, value in _HEADERS:
            self.send_header(header, value)
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self) -> None:  # noqa: N802
        self._dispatch("GET")

    def do_POST(self) -> None:  # noqa: N802
        self._dispatch("POST")


def serve(app: App, session: Session, *, host: str = LOOPBACK) -> ThreadingHTTPServer:
    """A server bound to ``host`` and ready for ``serve_forever``.

    The port is taken here, so a busy one is a named error before the link is printed.
    """

    class Handler(_Handler):
        pass

    Handler.session = session
    Handler.app = app
    Handler.routes = _routes(app)
    mimetypes.add_type("text/javascript", ".js")
    mimetypes.add_type("font/woff2", ".woff2")
    try:
        return ThreadingHTTPServer((host, session.port), Handler)
    except OSError as error:
        if error.errno == errno.EADDRINUSE:
            raise PortBusy(session.port) from error
        raise
