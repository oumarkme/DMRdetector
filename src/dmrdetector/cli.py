"""Command-line entry point: `dmrdetector <command> ...`.

The stage commands (preprocess, atlas, cpg-islands, icr, call-dmr, annotate,
target) are the steps of `build` and `run`, one script each, for scripting.
--help hides them.
"""

from __future__ import annotations

import argparse
import importlib
import os
import shlex
import sys

from dmrdetector import __version__


# Command -> (module, built at start, help line). `build`, `run` and `config`
# register their parsers at start. The others import pandas, scipy and numba,
# so their parsers are built when the command is parsed; a command without a
# help line is hidden from --help and from the error text of a wrong command.
COMMANDS = {
    "build": ("build", True, None),
    "run": ("run", True, None),
    "config": ("settings", True, None),
    "annotate": ("steps.annotate", False, None),
    "get-genes":
        ("steps.annotate", False, "download and verify an Ensembl GRCh38 GTF"),
    "preprocess": ("steps.preprocess", False, None),
    "atlas": ("steps.atlas", False, None),
    "cpg-islands": ("steps.cpg_islands", False, None),
    "icr": ("steps.icr", False, None),
    "call-dmr": ("steps.call_dmr", False, None),
    "target": ("steps.target", False, None),
}
ADVANCED_NOTE = "Advanced command for scripting (one step of build and run). "


class _LazyParser(argparse.ArgumentParser):
    """Stands in for a command's parser; builds the real one when used."""

    _module = ""
    _command = ""

    def _real(self):
        holder = argparse.ArgumentParser(prog="dmrdetector").add_subparsers()
        module = importlib.import_module(f"dmrdetector.{self._module}")
        module.add_parser(holder)
        parser = holder.choices[self._command]
        if COMMANDS[self._command][2] is None:
            helps = {a.dest: a.help for a in holder._choices_actions}
            parser.description = ADVANCED_NOTE + (
                parser.description or helps.get(self._command) or "")
        return parser

    def parse_known_args(self, args=None, namespace=None):
        """Parses `args` with the real parser of the command."""
        return self._real().parse_known_args(args, namespace)

    def format_help(self):
        """Returns the help text of the real parser."""
        return self._real().format_help()

    def format_usage(self):
        """Returns the usage line of the real parser."""
        return self._real().format_usage()


class _MainParser(argparse.ArgumentParser):
    """Main parser; a wrong command lists only the commands --help shows."""

    visible: list[str] = []

    def _check_value(self, action, value):
        if (isinstance(action, argparse._SubParsersAction) and
                value not in action.choices):
            raise argparse.ArgumentError(
                action, f"invalid choice: {value!r} (choose from "
                f"{', '.join(map(repr, self.visible))})")
        super()._check_value(action, value)


def build_parser() -> argparse.ArgumentParser:
    """Returns the main parser."""
    # The pipeline counts its own threads against --threads: BLAS and OpenMP
    # pools get one thread unless the caller set them (before numpy loads).
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                "NUMEXPR_NUM_THREADS"):
        os.environ.setdefault(key, "1")
    parser = _MainParser(
        prog="dmrdetector",
        description="DMR detection from long-read methylation data.",
        epilog="Options of build and run: dmrdetector build --help, "
        "dmrdetector run --help. Thresholds are set in the settings "
        "file (dmrdetector config init).")
    parser.add_argument("--version",
                        action="version",
                        version=f"dmrdetector {__version__}")
    visible = [n for n, (_, eager, text) in COMMANDS.items() if eager or text]
    parser.visible = visible
    sub = parser.add_subparsers(dest="command",
                                required=True,
                                metavar="{" + ",".join(visible) + "}",
                                title="commands")
    for name, (module, eager, text) in COMMANDS.items():
        if eager:
            importlib.import_module(f"dmrdetector.{module}").add_parser(sub)
            continue
        sub._parser_class = _LazyParser
        try:
            stub = sub.add_parser(name, help=text or argparse.SUPPRESS)
        finally:
            sub._parser_class = argparse.ArgumentParser
        stub._module, stub._command = module, name
    # also drops the ==SUPPRESS== rows of the hidden commands
    sub._choices_actions[:] = [
        a for a in sub._choices_actions if a.dest in visible
    ]
    return parser


def main(argv: list[str] | None = None) -> int:
    """Runs the command in `argv` (default: the command line)."""
    args = build_parser().parse_args(argv)
    # the cache of the compiled kernels: a private folder of the user, so that
    # users of one machine do not share (or fail on) /tmp/numba_cache; set
    # before the temporary directory moves into the work folder of a run
    from dmrdetector.util import jit
    jit.set_default_cache_dir()
    from dmrdetector import runtime
    runtime.COMMAND = "dmrdetector " + shlex.join(sys.argv[1:] if argv is
                                                  None else argv)
    return int(args.func(args) or 0)


if __name__ == "__main__":
    sys.exit(main())
