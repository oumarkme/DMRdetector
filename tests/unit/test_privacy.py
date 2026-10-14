"""Checks that a file git would publish holds no lab path or record number.

Lab storage paths, medical record numbers, study sample IDs, sample directory
names, lab host and data set names, institutional addresses and a scratch path
are patterns that must not appear in a published file.
"""

import fnmatch
import os
from pathlib import Path
import re
import shutil
import subprocess

import pytest


ROOT = Path(__file__).resolve().parents[2]
MAX_BYTES = 20 << 20

# The literals are split so that this file does not match itself. The study IDs
# of the lab (ONT followed by a letter and digits, or by digits only) and the
# record numbers (WGS-ONT..., then six or more digits, or D2 and four digits)
# are described here and never written out; the scanner test below builds
# synthetic examples from parts.
PATTERN = re.compile(("/gene"
                      "_lv|/nanopore"
                      "_nas|/lab_"
                      "nas|/Every"
                      "things/|WGS-ONT[A-Z0-9-]*-[0-9]{6,}|D2[0-9]{4}|"
                      r"\bONT-?[A-Z]?[0-9]{3,4}\b|"
                      r"(?i:\bgs"
                      r"4\b|subset"
                      r"14|/ra"
                      r"id/|\bMelo"
                      r"dy\b|\bNT"
                      r"UH\b|@ntu\."
                      "edu)").encode())


def _git_files() -> list[Path] | None:
    r = subprocess.run([
        "git", "-C",
        str(ROOT), "--no-optional-locks", "ls-files", "-co",
        "--exclude-standard", "-z"
    ],
                       capture_output=True)
    if r.returncode != 0:
        return None
    return [ROOT / p.decode() for p in r.stdout.split(b"\0") if p]


def _ignore_rules() -> list[tuple[str, bool, bool]]:
    """Returns the rules of the .gitignore at the project root.

    Each rule is (pattern, directories only, anchored); negations are not used.
    """
    rules = []
    gi = ROOT / ".gitignore"
    for line in (gi.read_text().splitlines() if gi.is_file() else []):
        s = line.strip()
        if not s or s.startswith("#") or s.startswith("!"):
            continue
        dir_only = s.endswith("/")
        s = s.rstrip("/")
        rules.append((s.lstrip("/"), dir_only, "/" in s))
    return rules


def _ignored(rel: str, is_dir: bool, rules) -> bool:
    name = rel.rsplit("/", 1)[-1]
    for pat, dir_only, anchored in rules:
        if dir_only and not is_dir:
            continue
        if fnmatch.fnmatchcase(rel if anchored else name, pat):
            return True
    return False


def _walk_files() -> list[Path]:
    """Returns the files below the project root that .gitignore keeps.

    Ignored directories are not entered.
    """
    rules, out = _ignore_rules(), []
    for cur, dirs, names in os.walk(ROOT):
        rel_dir = os.path.relpath(cur, ROOT).replace(os.sep, "/")
        prefix = "" if rel_dir == "." else rel_dir + "/"
        dirs[:] = sorted(
            d for d in dirs
            if d != ".git" and not _ignored(prefix + d, True, rules))
        out += [
            Path(cur) / n
            for n in sorted(names)
            if not _ignored(prefix + n, False, rules)
        ]
    return out


def _files() -> list[Path] | None:
    if shutil.which("git") is None:
        return None
    return _git_files() or _walk_files()


def test_no_lab_paths_or_record_numbers_in_publishable_files():
    files = _files()
    if files is None:
        pytest.skip(
            "git is not available or the project is not a git work tree")
    hits = []
    for f in files:
        try:
            if not f.is_file() or f.stat().st_size > MAX_BYTES:
                continue
            m = PATTERN.search(f.read_bytes())
        except OSError:
            continue
        if m:
            hits.append(f"{f.relative_to(ROOT)}: "
                        f"{m.group(0).decode(errors='replace')[:30]}")
    assert not hits, "private path or record number in: " + "; ".join(hits)


def _synthetic(*parts: str) -> bytes:
    """Joins parts into text; the examples below are made up, not real IDs."""
    return "".join(parts).encode()


def test_scan_pattern_catches_what_it_should():
    nines = "9" * 4
    record = _synthetic("WGS-", "ONT", "-Z", nines, "-", "1" * 8, "-D2", nines)
    study_a = _synthetic("ONT", "-Z", nines)  # letter and four digits
    study_b = _synthetic("ONT", "9" * 3)  # three digits only
    assert PATTERN.search(b"x /gene"
                          b"_lv/NGS y")
    assert PATTERN.search(b"see " + record + b" in the table")
    assert PATTERN.search(_synthetic("D2", nines))
    assert PATTERN.search(b"see " + study_a + b" and more")
    assert PATTERN.search(study_b)
    assert PATTERN.search(b"run on gs"
                          b"4 today")
    assert PATTERN.search(b"the subset"
                          b"14 data")
    assert PATTERN.search(b"/ra"
                          b"id/tmp")
    assert PATTERN.search(b"jo@ntu."
                          b"edu.tw")
    assert not PATTERN.search(b"a subset of 14 samples on a 96-core server")
    assert not PATTERN.search(b"ONT reads, ONT-based, ONTARGET and ONT 1000")
    assert not PATTERN.search(b"sample IDs like CASE01 or CTRL12")


def test_directory_walk_lists_the_files_git_lists():
    if shutil.which("git") is None or _git_files() is None:
        pytest.skip(
            "git is not available or the project is not a git work tree")
    git = {p.relative_to(ROOT).as_posix() for p in _git_files() if p.exists()}
    walk = {p.relative_to(ROOT).as_posix() for p in _walk_files()}
    # git leaves out empty directories and files below nested repositories; the
    # walk may list a file that git lists as deleted-but-tracked only through
    # the index, so compare the files that exist
    assert walk == git, sorted(walk ^ git)[:10]
