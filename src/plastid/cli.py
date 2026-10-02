"""Command-line interface for PlastidToolSuite.

This module is the only place allowed to terminate the process (directly or by
returning an exit code to the console-script wrapper). Library code must raise
exceptions instead.

Exit codes:
    0 — success
    1 — general error / not implemented yet
    2 — usage error, or a failed gate (e.g. the IR gate)
    3 — missing or failed external tool
"""

from __future__ import annotations

import argparse

from plastid import __version__

EXIT_OK = 0
EXIT_ERROR = 1
EXIT_USAGE = 2
EXIT_TOOL = 3


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="plastid",
        description="PlastidToolSuite — plastid and mitochondrial genome analysis.",
    )
    parser.add_argument("--version", action="store_true", help="print the version and exit")
    subparsers = parser.add_subparsers(dest="command", metavar="COMMAND")
    subparsers.add_parser("run", help="run an analysis module (not implemented yet)")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.version:
        print(f"plastid {__version__}")
        return EXIT_OK

    if args.command is None:
        parser.print_help()
        return EXIT_OK

    if args.command == "run":
        print("run: not implemented yet — arrives in P1-08 (see ROADMAP.md M1.8)")
        return EXIT_ERROR

    # Unreachable: argparse rejects unknown subcommands before this point.
    return EXIT_USAGE
