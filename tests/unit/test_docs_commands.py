"""Checks commands, flags, settings keys and file links of the Markdown docs.

Every `dmrdetector ...` command in a code block or inline code must name an
existing command and existing flags, every `--set KEY=VALUE` must use a known
settings key, and every link or `src/`, `doc/`, `env/`, `tests/` path in the
docs must exist.
"""

from pathlib import Path
import re
import shlex

import pytest

from dmrdetector import cli
from dmrdetector import settings

ROOT = Path(__file__).resolve().parents[2]
DOCS = sorted([
    ROOT / "README.md", ROOT / "CONTRIBUTING.md", ROOT / "CHANGELOG.md", ROOT /
    "tests" / "README.md"
] + list((ROOT / "doc").glob("*.md")))
DOCS = [d for d in DOCS if d.is_file()]
# CHANGELOG entries describe the past: only links and paths are not checked
HISTORY = {"CHANGELOG.md"}

COMMAND = re.compile(r"\bdmrdetector(?::[\w.]+)?\s+"
                     r"(build|run|config|get-genes|preprocess|atlas|"
                     r"cpg-islands|icr|call-dmr|annotate|target)\b(.*)")


def _flags_of(parser) -> set[str]:
    """Returns the option strings of a parser and of its sub-parsers."""
    if hasattr(parser, "_real"):
        parser = parser._real()
    out = set()
    for a in parser._actions:
        out.update(o for o in a.option_strings if o.startswith("--"))
        choices = getattr(a, "choices", None)
        if isinstance(choices, dict):
            for p in choices.values():
                out |= _flags_of(p)
    return out


@pytest.fixture(scope="module")
def flags() -> dict[str, set[str]]:
    top = cli.build_parser()
    action = next(a for a in top._actions if isinstance(a.choices, dict))
    return {name: _flags_of(p) for name, p in action.choices.items()}


def _commands(path: Path) -> list[tuple[int, str, str]]:
    """Returns (line number, command, rest of the command line)."""
    found = []
    lines = path.read_text().splitlines()
    in_fence = False
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.lstrip().startswith("```"):
            in_fence = not in_fence
            i += 1
            continue
        start = i
        if in_fence:
            while line.rstrip().endswith("\\") and i + 1 < len(lines):
                i += 1
                line = line.rstrip()[:-1] + " " + lines[i].strip()
            texts = [line]
        else:
            texts = re.findall(r"`([^`]*dmrdetector[^`]*)`", line)
        for t in texts:
            if re.search(r"dmrdetector\(\)", t):
                continue
            m = COMMAND.search(t)
            if m:
                found.append((start + 1, m.group(1), m.group(2)))
        i += 1
    return found


CASES = [(d, n, c, r) for d in DOCS if d.name not in HISTORY
         for n, c, r in _commands(d)]


@pytest.mark.parametrize("doc,n,cmd,rest",
                         CASES,
                         ids=[f"{d.name}:{n}" for d, n, _, _ in CASES])
def test_command_flags_exist(doc, n, cmd, rest, flags):
    allowed = flags[cmd]
    try:
        tokens = shlex.split(rest.split(" | ")[0].split(" # ")[0])
    except ValueError:
        tokens = rest.split()
    for i, tok in enumerate(tokens):
        if tok.startswith("--") and re.fullmatch(r"--[a-z0-9_-]+(=.*)?", tok):
            assert tok.split("=")[0] in allowed, (
                f"{doc.name}:{n}: {tok} is not a flag of `{cmd}`")
        if tok == "--set" and i + 1 < len(tokens):
            key = tokens[i + 1].split("=")[0]
            assert key in settings.BY_KEY, (
                f"{doc.name}:{n}: --set {key} is not a settings key")


def test_flag_table_rows_exist(flags):
    """Every `--flag` that opens a table row of the README is a real flag."""
    every = set().union(*flags.values())
    bad = []
    for doc in DOCS:
        if doc.name in HISTORY:
            continue
        for n, line in enumerate(doc.read_text().splitlines(), 1):
            m = re.match(r"\|\s*`(--[a-z0-9_-]+)", line)
            if m and m.group(1) not in every:
                bad.append(f"{doc.name}:{n}: {m.group(1)}")
    assert not bad, bad


def test_settings_keys_in_docs_exist():
    """Every `key = value` line of a toml block of the docs is a known key."""
    bad = []
    for doc in DOCS:
        if doc.name in HISTORY:
            continue
        for block in re.findall(r"```toml\n(.*?)```", doc.read_text(), re.S):
            for line in block.splitlines():
                m = re.match(r"\s*([a-z_]+)\s*=", line)
                if m and m.group(1) not in settings.BY_KEY:
                    bad.append(f"{doc.name}: {m.group(1)}")
    assert not bad, bad


def test_links_and_paths_exist():
    bad = []
    for doc in DOCS:
        if doc.name in HISTORY:
            continue
        text = doc.read_text()
        targets = re.findall(r"\]\(([^)#\s]+)(?:#[^)]*)?\)", text)
        targets += re.findall(r"`((?:src|doc|env|tests|dev)/[A-Za-z0-9_./-]+)`",
                              text)
        for t in targets:
            if re.match(r"[a-z]+://", t) or t.startswith("mailto:"):
                continue
            if t.startswith(("dev/", "tests/data")) or "*" in t or "<" in t:
                continue
            base = ROOT if "/" in t and not t.startswith(".") and (
                doc.parent == ROOT or t.split("/")[0] in
                ("src", "doc", "env", "tests")) else doc.parent
            if not ((base / t).exists() or (doc.parent / t).exists()):
                bad.append(f"{doc.name}: {t}")
    assert not bad, bad
