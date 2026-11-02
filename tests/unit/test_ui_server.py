"""The server answers only its own page, only with the session's token, only on this computer."""

from __future__ import annotations

import http.client
import json
import shutil
import socket
import threading
from collections.abc import Iterator
from pathlib import Path

import pytest

from plastid.cli import main
from plastid.ui.app import App, Refused
from plastid.ui.server import LOOPBACK, PortBusy, Session, serve


@pytest.fixture
def running(tmp_path: Path, data_dir: Path) -> Iterator[tuple[Session, App, Path]]:
    workspace = tmp_path / "work"
    app = App(workspace)
    shutil.copy(data_dir / "NC_050833.1.gb", workspace / "records" / "NC_050833.1.gb")
    session = Session(port=0)
    server = serve(app, session)
    session.port = server.server_address[1]
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    yield session, app, workspace
    server.shutdown()
    server.server_close()


def _ask(
    session: Session,
    path: str,
    *,
    method: str = "GET",
    cookie: bool = True,
    host: str | None = None,
    origin: str | None = None,
    body: bytes | None = None,
    kind: str = "application/json",
) -> tuple[int, dict[str, str], bytes]:
    connection = http.client.HTTPConnection(LOOPBACK, session.port, timeout=10)
    headers = {"Host": host or f"{LOOPBACK}:{session.port}"}
    if cookie:
        headers["Cookie"] = f"pts_session={session.token}"
    if origin:
        headers["Origin"] = origin
    if body is not None:
        headers["Content-Type"] = kind
    connection.request(method, path, body=body, headers=headers)
    response = connection.getresponse()
    answer = response.read()
    connection.close()
    return response.status, dict(response.getheaders()), answer


def _own_origin(session: Session) -> str:
    return f"http://{LOOPBACK}:{session.port}"


def test_nothing_is_read_without_the_session(running) -> None:
    session, _, _ = running
    status, _, words = _ask(session, "/api/runs", cookie=False)
    assert status == 403 and b"plastid ui" in words


def test_a_wrong_token_opens_nothing(running) -> None:
    session, _, _ = running
    status, headers, _ = _ask(session, "/?token=wrong", cookie=False)
    assert status == 403 and "Set-Cookie" not in headers


def test_the_token_is_traded_for_a_cookie_and_leaves_the_address(running) -> None:
    session, _, _ = running
    status, headers, _ = _ask(session, f"/?token={session.token}", cookie=False)
    assert status == 303 and headers["Location"] == "/"
    cookie = headers["Set-Cookie"]
    assert "HttpOnly" in cookie and "SameSite=Strict" in cookie


def test_a_rebound_host_name_is_refused(running) -> None:
    """A hostile name pointed at 127.0.0.1 arrives with its own Host, cookie or not."""
    session, _, _ = running
    status, _, _ = _ask(session, "/api/runs", host=f"evil.example:{session.port}")
    assert status == 403


def test_another_sites_page_is_refused(running) -> None:
    session, _, _ = running
    status, _, _ = _ask(session, "/api/runs", origin="https://evil.example")
    assert status == 403


def test_a_change_with_no_origin_is_refused(running) -> None:
    session, _, _ = running
    status, _, _ = _ask(session, "/api/jobs", method="POST", body=b"{}")
    assert status == 403


def test_a_form_must_be_json(running) -> None:
    """A plain HTML form on another page cannot send JSON, which is what closes that route."""
    session, _, _ = running
    status, _, _ = _ask(
        session,
        "/api/jobs",
        method="POST",
        origin=_own_origin(session),
        body=b"command=run",
        kind="application/x-www-form-urlencoded",
    )
    assert status == 400


def test_every_answer_refuses_to_be_framed_or_to_load_from_elsewhere(running) -> None:
    session, _, _ = running
    _, headers, _ = _ask(session, "/api/app")
    assert "default-src 'none'" in headers["Content-Security-Policy"]
    assert headers["Referrer-Policy"] == "no-referrer"
    assert headers["X-Frame-Options"] == "DENY"


def test_the_page_loads_nothing_from_outside_the_package(running) -> None:
    session, _, _ = running
    for name in ("/", "/app.js", "/app.css"):
        status, _, body = _ask(session, name)
        assert status == 200
        for address in (b"http://", b"https://"):
            assert address not in body.replace(b"http://www.w3.org/2000/svg", b"")


@pytest.mark.parametrize(
    "path",
    ["/api/file?path=../../etc/passwd", "/api/text?path=../secret.md", "/api/browse?path=.."],
)
def test_a_path_outside_the_workspace_is_refused(running, path: str) -> None:
    session, _, _ = running
    status, _, _ = _ask(session, path)
    assert status in (400, 404)


def test_a_symbolic_link_out_of_the_workspace_is_refused(running, tmp_path: Path) -> None:
    session, _, workspace = running
    secret = tmp_path / "secret.txt"
    secret.write_text("no", encoding="utf-8")
    (workspace / "records" / "link.txt").symlink_to(secret)
    status, _, _ = _ask(session, "/api/file?path=records/link.txt")
    assert status == 404


def test_an_upload_lands_in_the_workspace_and_overwrites_nothing(running) -> None:
    session, _, workspace = running
    origin = _own_origin(session)
    path = "/api/upload?name=new.fasta&to=records"
    status, _, answer = _ask(session, path, method="POST", origin=origin, body=b">a\nACGT\n")
    assert status == 200 and json.loads(answer)["path"] == "records/new.fasta"
    assert (workspace / "records" / "new.fasta").read_bytes() == b">a\nACGT\n"
    status, _, _ = _ask(session, path, method="POST", origin=origin, body=b">b\nTTTT\n")
    assert status == 400
    assert (workspace / "records" / "new.fasta").read_bytes() == b">a\nACGT\n"


def test_an_upload_name_cannot_climb_out(running) -> None:
    session, _, _ = running
    status, _, _ = _ask(
        session,
        "/api/upload?name=../../evil&to=records",
        method="POST",
        origin=_own_origin(session),
        body=b"x",
    )
    assert status == 400


def test_a_busy_port_stops_with_a_name_and_does_not_move(running) -> None:
    session, app, _ = running
    with pytest.raises(PortBusy):
        serve(app, Session(port=session.port))


def test_the_workspace_path_check_refuses_a_climb(tmp_path: Path) -> None:
    app = App(tmp_path / "work")
    with pytest.raises(Refused):
        app.upload_target("a.gb", "../..")


def test_the_gate_answer_is_the_engines(running) -> None:
    """A record's gate reads as `core/ir.py` says: the golden record's IR is 26,214 bp at 100 %."""
    session, _, _ = running
    status, _, answer = _ask(session, "/api/record?path=records/NC_050833.1.gb")
    assert status == 200
    record = json.loads(answer)
    gate = record["gates"]["annotation"]
    assert gate["status"] == "passed" and gate["length"] == 26214 and gate["identity_pct"] == 100.0
    assert [region["name"] for region in gate["regions"]] == ["LSC", "IRa", "SSC", "IRb"]
    assert record["length"] == 153724 and record["genes"]["total"] == 132


def test_a_folder_that_cannot_be_written_is_named(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch
) -> None:
    """A mounted folder can refuse the user inside a container, and that is not a traceback."""

    def refused(self, *args, **kwargs):  # type: ignore[no-untyped-def]
        raise PermissionError(13, "Permission denied", str(self))

    monkeypatch.setattr(Path, "mkdir", refused)
    assert main(["ui", "--workspace", str(tmp_path), "--no-browser"]) == 2
    err = capsys.readouterr().err
    assert "--workspace" in err and "cannot write" in err


def test_a_workspace_that_is_not_there_yet_is_made(tmp_path: Path) -> None:
    """`plastid ui --workspace my_project` is the first command a person types.

    A busy port ends the command right after the folder, so the test needs no running server.
    """
    workspace = tmp_path / "my_project"
    with socket.socket() as taken:
        taken.bind((LOOPBACK, 0))
        taken.listen()
        port = str(taken.getsockname()[1])
        code = main(["ui", "--workspace", str(workspace), "--no-browser", "--port", port])
    assert code == 2
    assert (workspace / "records").is_dir()


def test_a_workspace_that_is_a_file_is_refused(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    taken = tmp_path / "a_file"
    taken.write_text("x", encoding="utf-8")
    assert main(["ui", "--workspace", str(taken), "--no-browser"]) == 2
    assert "--workspace" in capsys.readouterr().err
