"""Command-line interface for PlastidToolSuite.

This module is the only place allowed to terminate the process (directly or by returning an
exit code to the console-script wrapper). Library code raises exceptions instead.

Exit codes and what each one means: ``plastid.run.state.EXIT_STATES``.
"""

from __future__ import annotations

import argparse
import os
import signal
import sys
import threading
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from pathlib import Path

from plastid import __version__
from plastid.cli.annotate import _cmd_annotate, _cmd_mito
from plastid.cli.assemble import _cmd_assemble
from plastid.cli.catalogue import _cmd_modules, _cmd_runs, _cmd_tools
from plastid.cli.common import (
    EXIT_CANCELLED,
    EXIT_CONFLICT,
    EXIT_ERROR,
    EXIT_OK,
    EXIT_TOOL,
    EXIT_USAGE,
    _checked_sets,
    _coerce,
)
from plastid.cli.editing import _cmd_editing_reference
from plastid.cli.outputs import (
    OutputFigures,
    _cmd_figures,
    _cmd_ir_junction,
    _cmd_package,
    _cmd_qc,
    _render_output_figures,
    _render_r_figures,
)
from plastid.cli.parser import build_parser
from plastid.cli.reference import _cmd_panel, _cmd_reference
from plastid.cli.run import (
    _cmd_batch,
    _cmd_run,
    _ir_manifest_fields,
    _ir_status,
    _resolve_recorded_reference,
    _with_data_label,
)
from plastid.run.tools import RunCancelled
from plastid.ui.catalogue import LOOPBACK

# The command line's names, importable from ``plastid.cli`` as they were when it was one file:
# the entry point, the exit codes, and the helpers the tests exercise directly.
__all__ = [
    "EXIT_CANCELLED",
    "EXIT_CONFLICT",
    "EXIT_ERROR",
    "EXIT_OK",
    "EXIT_TOOL",
    "EXIT_USAGE",
    "OutputFigures",
    "build_parser",
    "main",
    "_checked_sets",
    "_coerce",
    "_cmd_figures",
    "_ir_manifest_fields",
    "_ir_status",
    "_render_output_figures",
    "_render_r_figures",
    "_resolve_recorded_reference",
    "_with_data_label",
]

_STOP_SIGNALS = tuple(
    getattr(signal, name) for name in ("SIGINT", "SIGTERM", "SIGHUP") if hasattr(signal, name)
)


def _raise_cancelled(signum: int, _frame: object) -> None:
    # Once only: the run is now writing down that it stopped, and must be left to finish doing so.
    for sent in _STOP_SIGNALS:
        signal.signal(sent, signal.SIG_IGN)
    raise RunCancelled(signal.Signals(signum).name)


@contextmanager
def _signals_cancel() -> Iterator[None]:
    """For one command, a signal to stop means the run stops in good order.

    A signal already ignored stays ignored: that is how `nohup` keeps a long run alive after its
    terminal closes, and a run started that way has asked not to be stopped by it. SIGHUP is here
    because every tool runs in its own session, so a closing terminal reaches this process alone.
    """
    if threading.current_thread() is not threading.main_thread():
        yield
        return
    previous = {sent: signal.getsignal(sent) for sent in _STOP_SIGNALS}
    for sent, handler in previous.items():
        if handler != signal.SIG_IGN:
            signal.signal(sent, _raise_cancelled)
    try:
        yield
    finally:
        for sent, handler in previous.items():
            if handler is not None:
                signal.signal(sent, handler)


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    # What was typed, so a run's record can repeat it. Each command takes it from here rather than
    # rebuilding it, because a list rebuilt by hand is a list that loses a flag.
    args.typed = ["plastid", *(sys.argv[1:] if argv is None else argv)]
    with _signals_cancel():
        try:
            return _dispatch(parser, args)
        except RunCancelled:
            print("cancelled", file=sys.stderr)
            return EXIT_CANCELLED


def _subcommands(parser: argparse.ArgumentParser) -> dict[str, argparse.ArgumentParser]:
    found = [a for a in parser._actions if isinstance(a, argparse._SubParsersAction)]
    return dict(found[0].choices) if found else {}


def _cmd_ui(args: argparse.Namespace) -> int:
    """Serve the local interface until stopped. Work it started carries on without it."""
    import webbrowser

    from plastid.ui.app import App
    from plastid.ui.server import PortBusy, Session, serve

    if args.host != LOOPBACK and not os.environ.get("PLASTID_IN_CONTAINER"):
        print(
            f"--host: the interface listens on {LOOPBACK} only; another address is for the "
            "Docker image, which publishes the port to this computer alone",
            file=sys.stderr,
        )
        return EXIT_USAGE
    workspace = Path(args.workspace)
    token = os.environ.get("PLASTID_UI_TOKEN")
    session = Session(args.port, token) if token else Session(args.port)
    try:
        app = App(workspace, version=__version__, shown_as=os.environ.get("PLASTID_HOST_FOLDER"))
    except OSError as error:
        print(f"--workspace: cannot write to {workspace}: {error}", file=sys.stderr)
        return EXIT_USAGE
    try:
        server = serve(app, session, host=args.host)
    except PortBusy:
        print(
            f"port {args.port} is in use: stop what holds it, or choose another with --port",
            file=sys.stderr,
        )
        return EXIT_USAGE
    print(f"PlastidToolSuite is open at {session.link()}", flush=True)
    print(
        "Keep this window open while you work. Stopping the toolkit stops any run still going."
        if os.environ.get("PLASTID_IN_CONTAINER")
        else "Keep this window open while you work. Ctrl-C closes it; runs already started go on.",
        flush=True,
    )
    if not args.no_browser:
        webbrowser.open(session.link())
    try:
        server.serve_forever()
    except (KeyboardInterrupt, RunCancelled):
        print("closed")
    finally:
        server.server_close()
    return EXIT_OK


# One handler per subcommand; argparse has already refused any other name.
_COMMANDS: dict[str, Callable[[argparse.Namespace], int]] = {
    "ui": _cmd_ui,
    "modules": _cmd_modules,
    "tools": _cmd_tools,
    "runs": _cmd_runs,
    "run": _cmd_run,
    "figures": _cmd_figures,
    "ir-junction": _cmd_ir_junction,
    "package": _cmd_package,
    "qc": _cmd_qc,
    "batch": _cmd_batch,
    "annotate": _cmd_annotate,
    "mito": _cmd_mito,
    "panel": _cmd_panel,
    "reference": _cmd_reference,
    "editing-reference": _cmd_editing_reference,
    "assemble": _cmd_assemble,
}


def _dispatch(parser: argparse.ArgumentParser, args: argparse.Namespace) -> int:
    if args.version:
        print(f"plastid {__version__}")
        return EXIT_OK

    if args.command is None:
        parser.print_help()
        return EXIT_OK

    actions = _subcommands(_subcommands(parser)[args.command])
    if actions and args.action is None:
        print(f"{args.command}: choose an action — {' | '.join(actions)}", file=sys.stderr)
        return EXIT_USAGE

    handler = _COMMANDS.get(args.command)
    if handler is not None:
        return handler(args)

    # Unreachable: argparse rejects unknown subcommands before this point.
    return EXIT_USAGE
