"""Run support of `build` and `run`: budget, pool, log, stamps.

`--threads N` is the only resource option: tasks that run together use at most N
threads and 8 GB x N of estimated memory (`thread_plan`, `pool_map`). A step
writes a stamp when it finishes, so a rerun of `build` skips it.
"""

from __future__ import annotations

from collections.abc import Callable
from collections.abc import Iterable
from collections.abc import Iterator
from concurrent.futures import as_completed
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures.process import BrokenProcessPool
import contextlib
import csv
from dataclasses import dataclass
from dataclasses import field
from datetime import datetime
import errno
import fcntl
import gzip
import json
import multiprocessing as mp
import os
from pathlib import Path
import re
import shutil
import stat
import sys
import time
import traceback
import zlib

from dmrdetector import __version__
from dmrdetector import settings


# 16 GB: enough for the largest single task, the smoothing of chr1 (10.2 GB, the
# held arrays and 1 GB)
MIN_THREADS = 2
DEFAULT_THREAD_FRACTION = 0.5
GB_PER_THREAD = 8.0  # memory budget = GB_PER_THREAD x threads

# Peak memory (GB) of one task, measured with an end-to-end timing script (peak
# tree PSS) on a 96-core server. A task is never started unless the summed
# estimates of the tasks running together fit in the budget.
# preprocess: one sample (12 controls in parallel: 33 GB); panel: one scope on
# top of the loaded atlas; call_dmr: one sample and scope; annotate: one sample
# with the gene table inherited from the parent; gene_table: the GTF load.
TASK_GB = {
    "preprocess": 2.5,
    "panel": 2.0,
    "call_dmr": 5.0,
    "annotate": 1.0,
    "gene_table": 2.0,
    "target_report": 4.0,
}
# Raw atlas (atlas.build_atlas). The input is read window by window through the
#   tabix index, so the memory is
#   floor + workers x window,   floor = results kept until the scope-wide
#     median of PHI + base,
#   window = peak of one window worker = RAW_WINDOW_ROW_BYTES x stacked rows
#     of the window x RAW_DENSITY_SAFETY.
# Sites retained: at most stacked rows / min_samples and at most one per
#   RAW_BP_PER_SITE_MIN bp of the scope
# (14 controls, autosomal: 27,968,997 sites over 2.875 Gbp). Constants measured
#   on a 96-core server on 2026-10-06.
# compressed input bytes per stacked row (14 controls: 2.51 GB for 328.6 M rows)
RAW_ROW_FILE_BYTES = 7.6
RAW_BP_PER_SITE_MIN = 100.0
# floor per retained site (12 samples 28.0 M sites: 2.2 GB, 40 samples 28.7 M
# sites: 3.2 GB peak at 2 workers)
RAW_SITE_BYTES = 70.0
RAW_BASE_GB = 0.6  # interpreter, numba, tabix indexes, writer
# peak per stacked row of a window worker (measured 132 to 138 B, 12 and 40
# samples)
RAW_WINDOW_ROW_BYTES = 140.0
# rows per bp of a window against the scope mean (10 Mbp windows of chr19: 2.0)
RAW_DENSITY_SAFETY = 2.5
# bp; also the grid of the output blocks (atlas.WINDOW)
RAW_WINDOW_MAX = 10_000_000
RAW_WINDOW_MIN = RAW_WINDOW_MAX // 128
# summed chromosome length of each scope (GRCh38): extent of the windows of the
# raw atlas
SCOPE_EXTENT = {
    "autosomal": 2_875_001_522,
    "male_X": 156_040_895,
    "female_X": 156_040_895,
    "Y": 57_227_415
}


def raw_atlas_floor_gb(input_bytes: float, extent_bp: float,
                       min_samples: int) -> float:
    """Returns the memory (GB) of the raw atlas results kept until the end."""
    rows = input_bytes / RAW_ROW_FILE_BYTES
    sites = min(rows / max(1, min_samples), extent_bp / RAW_BP_PER_SITE_MIN)
    return RAW_BASE_GB + sites * RAW_SITE_BYTES / 1e9


def raw_atlas_window_gb(input_bytes: float, extent_bp: float,
                        window_bp: float) -> float:
    """Estimated peak memory (GB) of one window worker."""
    rows_per_bp = input_bytes / RAW_ROW_FILE_BYTES / max(1.0, extent_bp)
    return (rows_per_bp * window_bp * RAW_WINDOW_ROW_BYTES *
            RAW_DENSITY_SAFETY / 1e9)


def raw_atlas_min_gb(input_bytes: float, extent_bp: float,
                     min_samples: int) -> float:
    """Returns the least memory (GB) the raw atlas can run in.

    That is the floor plus one worker on the smallest window.
    """
    return raw_atlas_floor_gb(input_bytes, extent_bp,
                              min_samples) + raw_atlas_window_gb(
                                  input_bytes, extent_bp, RAW_WINDOW_MIN)


def plan_raw_atlas(input_bytes: float, extent_bp: float, min_samples: int,
                   threads: int, mem_gb: float) -> tuple[int, int] | None:
    """Plans the raw atlas windows under a memory budget.

    The window is the largest of RAW_WINDOW_MAX / 2^k (k = 0..7) that lets
    `threads` workers run together; when no window does, the one that admits the
    most workers (the larger on a tie).

    Returns:
        (window in bp, workers), or None when not even one worker on the
        smallest window fits.
    """
    room = mem_gb - raw_atlas_floor_gb(input_bytes, extent_bp, min_samples)
    best = None
    w = RAW_WINDOW_MAX
    while w >= RAW_WINDOW_MIN:
        cost = raw_atlas_window_gb(input_bytes, extent_bp, w)
        k = max(1, threads) if cost <= 0 or room == float("inf") else min(
            max(1, threads), int(room // cost))
        if k >= 1 and (best is None or k > best[1]):
            best = (w, k)
        if k >= threads:
            break
        w //= 2
    return best


# Smoothing of one chromosome: the FFT runs in a short-lived child process
# (util.fftsmooth), so the scipy.fft plan cache of its length (8 bytes per
# padded base) goes away with the child. Peak of the child, measured on a
# 96-core server on 2026-10-07 for chr1 (a memory test, padded length 250
# million): 9.88 GB with the 5 FFTs in sequence (inner 1: 39.5 bytes per
# padded base, plan included), 17.66 GB with 3 FFTs at once (inner 3: 70.6
# bytes). The loaded raw atlas
# arrays, the results and the queue of compressed chunks stay in the parent
# (SMOOTH_HELD_BYTES_PER_SITE) and are counted by the caller. keyed by the
# number of FFT threads of the chromosome
SMOOTH_BYTES_PER_BASE = {1: 40.0, 3: 71.0}
# interpreter, numpy, scipy of the child
SMOOTH_CHILD_BASE_GB = 0.25
# parent beyond the arrays: interpreter, numba, queue of compressed chunks
SMOOTH_PARENT_BASE_GB = 1.0
SMOOTH_GB_MIN = 0.25
# arrays of the loaded raw atlas and the results, per site (14 controls: 1.76 GB
# for 27.97 M sites)
SMOOTH_HELD_BYTES_PER_SITE = 64.0
# longest chromosome of each scope (GRCh38)
SCOPE_LARGEST_SPAN = {
    "autosomal": 248_956_422,
    "male_X": 156_040_895,
    "female_X": 156_040_895,
    "Y": 57_227_415
}


def smooth_gb(span: int, inner: int = 1) -> float:
    """Returns the peak memory (GB) of smoothing a chromosome of `span` bp.

    This is the FFT child process with `inner` (1 or 3) FFT threads; the caller
    adds the arrays held by the parent.
    """
    per_base = SMOOTH_BYTES_PER_BASE[3 if inner > 1 else 1]
    return max(SMOOTH_GB_MIN, SMOOTH_CHILD_BASE_GB + per_base * span / 1e9)


SCOPES = ("autosomal", "male_X", "female_X", "Y")
SCOPE_SEX = {"autosomal": None, "male_X": "M", "female_X": "F", "Y": "M"}
SCOPE_PART = {
    "autosomal": "autosomal",
    "male_X": "chrX",
    "female_X": "chrX",
    "Y": "chrY"
}
PARTS = ("autosomal", "chrX", "chrY")

# command line as given to `dmrdetector` (set by cli.main), written to the log
COMMAND: str | None = None


class RunError(Exception):
    """User-facing error: printed as `ERROR: ...`, exit status 1."""


def container_hint() -> str:
    """Returns a hint about mounted folders when running in a container."""
    if (os.environ.get("APPTAINER_CONTAINER") or
            os.environ.get("SINGULARITY_CONTAINER") or
            os.path.exists("/.dockerenv")):
        return (" (the container sees only mounted folders: add --bind DIR "
                "for Apptainer or -v DIR:DIR for Docker)")
    return ""


def suffix(scope: str) -> str:
    """Returns the file name suffix of a scope ("" for autosomal)."""
    return "" if scope == "autosomal" else f"_{scope}"


class _Log:
    """Prints each message to stderr and appends it to the progress file.

    The file is shared with the worker processes. Messages before the path is
    set are kept and written first.
    """
    path: str | None = None

    def __init__(self) -> None:
        self.early: list[str] = []

    def set_path(self, path: str | None) -> None:
        """Sets the progress file and writes the early messages to it."""
        self.path = path
        if path:
            with open(path, "w") as fh:
                fh.write("".join(x + "\n" for x in self.early))
        self.early = []

    def __call__(self, msg: str) -> None:
        """Logs one message with a time stamp."""
        line = f"[{datetime.now():%H:%M:%S}] {msg}"
        print(line, file=sys.stderr, flush=True)
        if self.path:
            try:
                with open(self.path, "a") as fh:
                    fh.write(line + "\n")
            except OSError:
                pass
        else:
            self.early.append(line)


log = _Log()


@dataclass(frozen=True)
class Sample:
    """One sample of a samplesheet or of the `run` command line."""

    sample_id: str
    sex: str  # 'M' | 'F' | '' (not given)
    meth_bed: str
    vcf: str
    meth_bed_hap1: str | None = None  # target-region report only
    meth_bed_hap2: str | None = None


SAMPLE_ID_RE = r"[A-Za-z0-9._-]+"
# values of the `type` column that mark a control (compared in lower case)
CONTROL_TYPES = ("ctrl", "control")


def check_sample_id(sid: str, where: str = "") -> None:
    """Raises RunError unless `sid` is made of letters, digits, . _ and -."""
    if not re.fullmatch(SAMPLE_ID_RE, sid or ""):
        raise RunError(f"{where}sample_id '{sid}' may contain only letters, "
                       "digits, '.', '_' and '-'" +
                       ("" if sid else " (it is empty)"))


def _relative_hint(value: str, sheet: Path) -> str:
    """Returns a hint for a relative path that is read from another folder.

    A relative path of a samplesheet counts from the working directory, not
    from the folder of the samplesheet.
    """
    if os.path.isabs(value) or sheet.resolve().parent == Path.cwd().resolve():
        return ""
    return (" (a relative path counts from the working directory "
            f"{os.getcwd()}, not from the folder of the samplesheet; run the "
            "command in that folder or use absolute paths)")


def _check_types(rows: list[list[str]], header: list[str], path: str) -> int:
    """Raises RunError for a `type` value that is not a control or a case.

    Args:
        rows: Rows of the samplesheet, header first.
        header: The header row.
        path: Samplesheet path for the message.

    Returns:
        The number of rows whose type is CASE (any case).

    Raises:
        RunError: A row has another value, including an empty one (v2.2 stopped
            there as well). The message lists the distinct values.
    """
    col = header.index("type")
    unknown: dict[str, int] = {}
    n_case = 0
    for idx, cells in enumerate(rows[1:]):
        if not any(cells):
            continue
        value = cells[col] if col < len(cells) else ""
        if value.lower() in CONTROL_TYPES:
            continue
        if value.lower() == "case":
            n_case += 1
        else:
            unknown.setdefault(value, idx + 2)
    if unknown:
        shown = [
            f"'{v}' (row {n})" if v else f"empty (row {n})"
            for v, n in list(unknown.items())[:5]
        ]
        more = (f" and {len(unknown) - 5} more" if len(unknown) > 5 else "")
        raise RunError(
            f"samplesheet {path}: the 'type' column has {len(unknown)} "
            f"value(s) that are neither a control nor CASE: "
            f"{', '.join(shown)}{more}. Write CTRL (or control) for a control "
            "and CASE for a row that build leaves out")
    return n_case


def parse_samplesheet(path: str) -> list[Sample]:
    """Validates and reads the control samplesheet of `build`.

    Columns: sample_id, sex (M or F; needed for the sex-specific scopes),
    meth_bed, vcf. When a `type` column exists (the v2.2 samplesheet), only the
    rows whose type is CTRL or control (any case) are read, as in v2.2; the
    rows of type CASE (any case) are dropped without further checks and counted
    in a NOTICE; any other value, an empty one included, is an error, as in
    v2.2. A sheet without `type` is all controls.

    Args:
        path: CSV file.

    Returns:
        The control samples in file order.

    Raises:
        RunError: The file is missing, empty or invalid, a type is neither a
            control nor CASE, or no row is a control.
    """
    p = Path(path)
    if p.exists() and not p.is_file():
        raise RunError(f"samplesheet {path} is not a regular file")
    if not p.is_file():
        raise RunError(f"samplesheet {path} not found{container_hint()}")
    with open(p, newline="", encoding="utf-8-sig") as fh:
        rows = [[c.strip() for c in r] for r in csv.reader(fh)]
    if not rows or not any(any(c for c in r) for r in rows):
        raise RunError(f"samplesheet {path} is empty")
    header = rows[0]
    missing = [c for c in ("sample_id", "meth_bed", "vcf") if c not in header]
    if missing:
        raise RunError(f"samplesheet {path} is missing required column(s): "
                       f"{', '.join(missing)}. "
                       "Expected header: sample_id,sex,meth_bed,vcf")
    has_type = "type" in header
    # the types are checked first, so that a typo is reported before the files
    # of rows that would be dropped
    n_dropped = _check_types(rows, header, path) if has_type else 0
    out: list[Sample] = []
    seen: set[str] = set()
    n_rows = 0
    for idx, cells in enumerate(rows[1:]):
        n = idx + 2
        if not any(cells):
            continue
        n_rows += 1
        cells = cells + [""] * (len(header) - len(cells))
        row = dict(zip(header, cells))
        if has_type and row["type"].lower() not in CONTROL_TYPES:
            continue  # CASE: counted by _check_types
        sid = row["sample_id"]
        if not sid:
            raise RunError(f"samplesheet row {n}: sample_id is empty")
        if sid in seen:
            raise RunError(f"samplesheet row {n}: duplicate sample_id '{sid}'")
        seen.add(sid)
        check_sample_id(sid, f"samplesheet row {n}: ")
        sex = row.get("sex", "")
        if sex and sex.lower() not in ("m", "f"):
            raise RunError(f"samplesheet row {n} (sample_id={sid}): sex "
                           f"must be M or F, got '{sex}'")
        files = {}
        for col in ("meth_bed", "vcf"):
            v = row.get(col, "")
            if not v:
                raise RunError(
                    f"samplesheet row {n} (sample_id={sid}): {col} is empty")
            if not Path(v).is_file():
                raise RunError(f"samplesheet row {n} (sample_id={sid}): "
                               f"{col} file not found: {v}"
                               f"{_relative_hint(v, p)}{container_hint()}")
            files[col] = os.path.abspath(v)
        out.append(Sample(sid, sex.upper(), **files))
    if not n_rows:
        raise RunError(f"samplesheet {path} has a header but no samples")
    if has_type:
        if not out:
            raise RunError(f"samplesheet {path} has a 'type' column but no row "
                           "with type CTRL (or control)")
        if n_dropped:
            log(f"NOTICE: samplesheet {path}: the 'type' column keeps the "
                f"CTRL rows only; {n_dropped} of {n_rows} rows are dropped, "
                f"{len(out)} kept as controls")
    return out


def require_sex(samples: list[Sample], why: str) -> None:
    """Raises RunError when a sample has no sex (M or F)."""
    bad = [s.sample_id for s in samples if s.sex not in ("M", "F")]
    if bad:
        raise RunError(
            "sex (M or F) is required for the sex-specific scopes "
            f"({why}); missing for: "
            f"{', '.join(bad[:5])}{' ...' if len(bad) > 5 else ''}. "
            "Add a sex column, or set include_sex_chromosomes = false")


def md5_file(path: str) -> str:
    """Returns the md5 hex digest of a file."""
    import hashlib  # loads OpenSSL (about 0.1 s); only checksums need it
    h = hashlib.md5()
    with open(path, "rb") as fh:
        while chunk := fh.read(1 << 20):
            h.update(chunk)
    return h.hexdigest()


@dataclass
class Ctx:
    """Context of one `build` or `run`: options, folders, budget, timings."""

    command: str
    cfg: dict
    outdir: Path
    workdir: Path
    threads: int
    mem_gb: float
    sources: dict = field(default_factory=dict)
    config_path: str | None = None
    timings: dict = field(default_factory=dict)
    # shown after the status in the log (for example a partial result)
    status_note: str = ""

    @property
    def stamps(self) -> str:
        """Folder of the resume stamps."""
        return str(self.workdir / "stamps")

    @property
    def scratch(self) -> Path:
        """Folder of the temporary files of this run."""
        return self.workdir / "scratch"


_TMP_SAVED: list = []


def set_tmp(path: str | os.PathLike) -> None:
    """Puts every temporary file below `path`, never in /tmp.

    This covers this process and its children (TMPDIR, tempfile). `path` is the
    scratch directory in the temporary folder of the run; `restore_tmp` undoes
    it.
    """
    import tempfile
    if not _TMP_SAVED:
        _TMP_SAVED.append((os.environ.get("TMPDIR"), tempfile.tempdir))
    os.makedirs(path, exist_ok=True)
    os.environ["TMPDIR"] = str(path)
    tempfile.tempdir = str(path)


def restore_tmp() -> None:
    """Restores TMPDIR and tempfile as `set_tmp` found them."""
    import tempfile
    if _TMP_SAVED:
        old_env, old_dir = _TMP_SAVED.pop()
        if old_env is None:
            os.environ.pop("TMPDIR", None)
        else:
            os.environ["TMPDIR"] = old_env
        tempfile.tempdir = old_dir


def prepare_dirs(ctx: Ctx, *sub: str) -> None:
    """Creates the work folder with its sub-folders.

    The progress log and every temporary file go below it.
    """
    for d in (ctx.workdir, Path(ctx.stamps), ctx.scratch,
              *[ctx.workdir / s for s in sub]):
        d.mkdir(parents=True, exist_ok=True)
    log.set_path(str(ctx.workdir / "progress.log"))
    set_tmp(ctx.scratch)


# The work folder of `build` is claimed with a marker file and a lock, so that
# the build never writes into, or deletes, a folder it did not create.
WORK_MARKER = ".dmrdetector-build"
WORK_LOCK = ".dmrdetector-build.lock"
# Layout of the work files; part of every resume stamp. Raise it when a release
# changes what a stamp vouches for.
STAMP_FORMAT = 1
# A temporary file of a dead process is removed only when it is also this old,
# because a process id means nothing for a file that another host or container
# writes to the same shared folder.
STALE_TMP_MIN_AGE_S = 3600.0
_STALE_TMP_RE = re.compile(r"^\.tmp(\d+)\..+")


def build_workdir(tmp_dir: str | None, outdir: Path) -> Path:
    """Returns the work folder of a `build` of `outdir`.

    Without `tmp_dir` it is <outdir>/work. With one, the name carries a hash of
    the resolved outdir, so that builds of different references can share the
    same --tmp-dir.
    """
    if not tmp_dir:
        return outdir / "work"
    import hashlib
    tag = hashlib.sha1(str(outdir).encode()).hexdigest()[:12]
    return Path(tmp_dir).resolve() / f"dmrdetector-build-{tag}"


@dataclass
class WorkClaim:
    """A work folder that this process created or resumed, and its lock."""

    path: Path
    outdir: Path
    resumed: bool
    lock_fd: int | None = None


def _read_marker(workdir: Path) -> dict | None:
    """Returns the marker of a work folder; {} if unreadable, None if absent."""
    try:
        with open(workdir / WORK_MARKER) as fh:
            info = json.load(fh)
    except FileNotFoundError:
        return None
    except (OSError, ValueError):
        return {}
    return info if isinstance(info, dict) else {}


def _write_marker(claim: WorkClaim, status: str) -> None:
    write_json(
        claim.path / WORK_MARKER,
        dict(program="dmrdetector build",
             version=__version__,
             outdir=str(claim.outdir),
             pid=os.getpid(),
             time=datetime.now().isoformat(timespec="seconds"),
             status=status))


def claim_workdir(workdir: Path, outdir: Path) -> WorkClaim:
    """Creates or resumes the work folder of a `build`, and locks it.

    A new folder gets a marker file. An existing folder is used only when it has
    the marker of a build of the same `outdir` (an interrupted build resumes),
    or is empty. The folder is locked with flock, so that a second build cannot
    use it at the same time.

    Raises:
        RunError: The folder exists without the marker, belongs to another
            outdir, or is locked by a running build.
    """
    hint = "move or remove it, or choose another folder with --tmp-dir"
    try:
        workdir.parent.mkdir(parents=True, exist_ok=True)
        workdir.mkdir()
        created = True
    except FileExistsError:
        created = False
    except OSError as e:
        raise RunError(
            f"cannot create the work folder {workdir}: {e}") from None
    resumed = False
    if not created:
        if not workdir.is_dir():
            raise RunError(f"work folder {workdir} exists and is not a folder; "
                           f"{hint}")
        info = _read_marker(workdir)
        if info is None:
            try:
                foreign = [n for n in os.listdir(workdir) if n != WORK_LOCK]
            except OSError as e:
                raise RunError(
                    f"cannot read the work folder {workdir}: {e}") from None
            if foreign:
                raise RunError(
                    f"work folder {workdir} exists and was not made by "
                    f"dmrdetector build (no {WORK_MARKER} file in it); the "
                    f"build does not use or delete it: {hint}")
        elif info.get("outdir") != str(outdir):
            raise RunError(f"work folder {workdir} belongs to a build of "
                           f"{info.get('outdir', 'an unknown folder')}, not of "
                           f"{outdir}; {hint}")
        else:
            resumed = info.get("status") == "running"
    claim = WorkClaim(workdir, outdir, resumed)
    claim.lock_fd = _lock_file(workdir / WORK_LOCK, workdir)
    _write_marker(claim, "running")
    return claim


def _lock_file(path: Path, workdir: Path) -> int | None:
    """Returns the descriptor of `path`, locked exclusively.

    Raises:
        RunError: Another process holds the lock.
    """
    fd = os.open(path, os.O_RDWR | os.O_CREAT, 0o644)
    try:
        fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError as e:
        if e.errno in (errno.EWOULDBLOCK, errno.EAGAIN):
            os.close(fd)
            raise RunError(
                f"another dmrdetector build is using the work folder "
                f"{workdir}; wait for it, or choose another folder with "
                "--tmp-dir") from None
        # a file system without locks (some network mounts): go on unlocked
        log(f"WARNING: cannot lock the work folder {workdir} ({e.strerror}); "
            "do not start a second build with the same --outdir")
    return fd


def release_workdir(claim: WorkClaim, success: bool, keep: bool) -> bool:
    """Ends the use of a claimed work folder.

    The folder is removed only after a success without `keep`, and only when
    its marker names this outdir.

    Returns:
        True if the folder was removed.
    """
    removed = False
    info = _read_marker(claim.path)
    ours = bool(info) and info.get("outdir") == str(claim.outdir)
    if ours and success and not keep:
        shutil.rmtree(claim.path, ignore_errors=True)
        removed = True
    elif ours and success:
        _write_marker(claim, "complete")
    if claim.lock_fd is not None:
        try:
            os.close(claim.lock_fd)
        except OSError:
            pass
        claim.lock_fd = None
    return removed


def _pid_exists(pid: int) -> bool:
    """Returns False only when no process of the id exists on this host."""
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except OSError:
        return True
    return True


def remove_stale_tmp(folder: Path, now: float | None = None) -> list[str]:
    """Removes the temporary files that a killed process left in `folder`.

    Only files named `.tmp<pid>.<name>` directly in `folder` count, and only
    when no process with that id exists and the file has not changed for
    STALE_TMP_MIN_AGE_S. Nothing below `folder` and nothing elsewhere is
    touched.

    Returns:
        The names removed.
    """
    now = time.time() if now is None else now
    removed = []
    try:
        names = sorted(os.listdir(folder))
    except OSError:
        return removed
    for name in names:
        m = _STALE_TMP_RE.match(name)
        if not m or _pid_exists(int(m.group(1))):
            continue
        path = os.path.join(folder, name)
        try:
            st = os.lstat(path)
            if (not stat.S_ISREG(st.st_mode) or
                    now - st.st_mtime < STALE_TMP_MIN_AGE_S):
                continue
            os.remove(path)
        except OSError:
            continue
        removed.append(name)
    if removed:
        log(f"removed {len(removed)} stale temporary file(s) of an earlier "
            f"killed run from {folder}: {', '.join(removed[:5])}" +
            (" ..." if len(removed) > 5 else ""))
    return removed


def stat_of(path: str) -> list:
    """Returns [absolute path, size, mtime in ns] for a stamp key."""
    st = os.stat(path)
    return [os.path.abspath(path), st.st_size, st.st_mtime_ns]


def _digest(obj) -> str:
    import hashlib
    return hashlib.sha1(json.dumps(obj, sort_keys=True,
                                   default=str).encode()).hexdigest()


def tmp_name(final: str) -> str:
    """Returns the temporary name of `final` (same folder, this process id)."""
    d, b = os.path.split(final)
    return os.path.join(d, f".tmp{os.getpid()}.{b}")


def commit(tmp: str, final: str) -> None:
    """Renames `tmp` to `final`, its .tbi first."""
    if os.path.exists(tmp + ".tbi"):
        os.replace(tmp + ".tbi", final + ".tbi")
    elif os.path.exists(final + ".tbi"):
        os.remove(final + ".tbi")
    os.replace(tmp, final)


def cleanup(*tmps: str) -> None:
    """Removes temporary files and their .tbi, ignoring missing ones."""
    for t in tmps:
        for f in (t, t + ".tbi"):
            try:
                os.remove(f)
            except OSError:
                pass


@contextlib.contextmanager
def temp_output(final: str, *also: str) -> Iterator[str]:
    """Yields the temporary name of `final`; removes it and `also` on leaving.

    The caller moves the finished file into place inside the block, with
    `commit` (which moves a .tbi too) or os.replace; what was moved is no longer
    there to remove.

    Args:
        final: Path of the finished file.
        *also: Further temporary files of the same step to remove on leaving.

    Yields:
        The temporary name of `final` (same folder, this process id).
    """
    tmp = tmp_name(final)
    try:
        yield tmp
    finally:
        cleanup(tmp, *also)


def write_atomic(path: str | os.PathLike, text: str) -> None:
    """Writes `text` to `path` through a temporary name.

    A failure leaves no temporary file behind and an existing `path` as it was.

    Args:
        path: Final path of the file.
        text: Content, written with the default text encoding.
    """
    tmp = tmp_name(str(path))
    try:
        with open(tmp, "w") as fh:
            fh.write(text)
        os.replace(tmp, path)
    except BaseException:
        cleanup(tmp)
        raise


def write_json(path: Path, obj: object) -> None:
    """Writes `obj` as JSON through a temporary name."""
    write_atomic(path, json.dumps(obj, indent=2, default=str))


def guarded(stamp_dir: str,
            name: str,
            key: object,
            outputs: list[str],
            fn: Callable[[], list[str] | None]) -> tuple[str, str, float]:
    """Runs fn() unless the stamp matches and all outputs exist.

    The stamp key includes the program version and STAMP_FORMAT, so a work
    folder kept by another release is not resumed. `fn` may return the warning
    lines it logged; they go into the stamp, and a run that skips the step logs
    them again, so that a resumed build reports the warnings of the whole
    build.

    Args:
        stamp_dir: Folder of the stamps.
        name: Step name (the stamp is <name>.json).
        key: JSON-serializable parameters and input stats of the step.
        outputs: Files the step must produce.
        fn: The step; returns the WARNING and NOTICE lines to replay on a
            resume, or None.

    Returns:
        (name, "done" or "skipped", seconds).

    Raises:
        RunError: An output is missing after fn().
    """
    stamp = os.path.join(stamp_dir, name + ".json")
    k = _digest(dict(version=__version__, format=STAMP_FORMAT, key=key))
    if all(os.path.exists(o) for o in outputs) and os.path.exists(stamp):
        try:
            with open(stamp) as fh:
                info = json.load(fh)
        except (OSError, ValueError):
            info = {}
        if info.get("key") == k:
            # the warnings of the run that did the step are part of this one
            for message in info.get("messages", ()):
                log(message)
            return name, "skipped", 0.0
    try:
        os.remove(stamp)
    except OSError:
        pass
    t0 = time.time()
    messages = fn()
    dt = time.time() - t0
    missing = [o for o in outputs if not os.path.exists(o)]
    if missing:
        raise RunError(f"step {name} did not produce {missing[0]}")
    payload = {
        "key": k,
        "seconds": round(dt, 2),
        "finished": datetime.now().isoformat(timespec="seconds")
    }
    if messages:
        payload["messages"] = list(messages)
    write_atomic(stamp, json.dumps(payload))
    return name, "done", dt


def safe(label: str, fn: Callable, *args: object) -> object:
    """Calls fn(*args); any error but RunError becomes a RunError."""
    try:
        return fn(*args)
    except RunError:
        raise
    except SystemExit as e:
        raise RunError(f"{label}: {e.code}") from None
    except (EOFError, gzip.BadGzipFile, zlib.error) as e:
        raise RunError(
            f"{label}: an input file is truncated or not valid gzip "
            f"({type(e).__name__}: {e}); check each input with `gzip -t FILE`"
        ) from None
    except Exception as e:
        raise RunError(
            f"{label}: {type(e).__name__}: {e}\n{traceback.format_exc(limit=3)}"
        ) from None


def run_command(command: Callable[[object], int], args: object) -> int:
    """Runs `command(args)`; a user error prints `ERROR: ...` and gives 1.

    Args:
        command: The command function, which takes the parsed arguments.
        args: The parsed command-line arguments.

    Returns:
        The exit status of `command`, or 1 after a RunError or ConfigError.
    """
    try:
        return command(args)
    except (RunError, settings.ConfigError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1


def record(ctx: Ctx, step: tuple) -> None:
    """Stores a (name, status, seconds) step in the timings and logs it."""
    name, status, secs = step
    ctx.timings[name] = {"status": status, "seconds": round(secs, 2)}
    log(f"{name}: {status}" + (f" ({secs:.0f} s)" if status == "done" else ""))


def check_files(pairs: Iterable[tuple[str, str | None]]) -> None:
    """Raises RunError for the first (label, path) whose file does not exist."""
    for label, f in pairs:
        if f and os.path.exists(f) and not os.path.isfile(f):
            raise RunError(f"{label}: not a regular file: {f}")
        if f and not os.path.isfile(f):
            raise RunError(f"{label} file not found: {f}{container_hint()}")


def count_lines(path: str, threads: int = 1) -> int:
    """Returns the number of lines of a text or gzip file.

    BGZF files are inflated in `threads` threads.
    """
    if path.endswith(".gz"):
        from dmrdetector.fileio import bgzf as bz
        try:
            n = bz.count_newlines(path, threads)
        except ValueError:
            n = None
        if n is not None:
            return n
    n = 0
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rb") as fh:
        while chunk := fh.read(1 << 24):
            n += chunk.count(b"\n")
    return n


def _fmt_value(v) -> str:
    """Returns the log text of a value."""
    return "(unset)" if v is None else str(v)


def package_versions() -> str:
    """Returns `python 3.12.15, numpy 2.5.3, ...` for the log.

    These are the interpreter and the numerical packages that decide the
    results; importlib.metadata reads the versions without importing them.
    """
    from importlib.metadata import PackageNotFoundError
    from importlib.metadata import version
    out = [f"python {sys.version.split()[0]}"]
    for name in ("numpy", "scipy", "pandas", "pysam", "numba"):
        try:
            out.append(f"{name} {version(name)}")
        except PackageNotFoundError:
            out.append(f"{name} (not installed)")
    return ", ".join(out)


PATH_OPTIONS = ("reference", "gtf", "target_region", "outdir", "tmp_dir",
                "icr_gene_list", "par_bed")


def _path_names(ctx: Ctx, inputs: dict,
                files: Iterable[str] = ()) -> dict[str, str]:
    """Maps every absolute path of a run to its last component.

    The paths are those of the inputs, of the path options, of the settings
    file, of the folders of the run and of the working directory. `files` are
    further input files (the samples of a samplesheet); they and their folders
    are mapped, so that a message about `<file>.tbi` loses the folder too.
    """
    values = [v for v in inputs.values() if isinstance(v, (str, os.PathLike))]
    values += [ctx.cfg.get(k) for k in PATH_OPTIONS]
    values += [ctx.config_path, ctx.outdir, ctx.workdir]
    for f in files:
        if f:
            values += [f, os.path.dirname(os.path.abspath(f))]
    out = {}
    for v in values:
        if not v:
            continue
        for p in (os.path.abspath(v), os.path.realpath(v)):
            name = os.path.basename(p.rstrip(os.sep))
            if name:
                out[p] = name
    cwd = os.getcwd()
    if cwd != os.sep:
        out.setdefault(cwd, ".")
    return dict(sorted(out.items(), key=lambda kv: -len(kv[0])))


def redact_paths(text: str, names: dict[str, str]) -> str:
    """Replaces each path of `names` in `text` by its last component."""
    for p, name in names.items():
        text = re.sub(r"(?<![A-Za-z0-9_./-])" + re.escape(p) +
                      r"(?![A-Za-z0-9_.-])",
                      lambda _m, n=name: n,
                      text)
    return text


def _input_text(value: object, by_name: bool) -> str:
    """Returns the log text of one input: a file name and size, or the value."""
    if value is not None and by_name and os.path.isabs(str(value)):
        p = str(value)
        if os.path.isfile(p):
            return f"{os.path.basename(p)} ({os.path.getsize(p):,} bytes)"
        return os.path.basename(p.rstrip(os.sep))
    return _fmt_value(value)


def write_log(path: Path,
              ctx: Ctx,
              title: str,
              started: float,
              inputs: dict,
              error: BaseException | None,
              files: Iterable[str] = ()) -> None:
    """Writes the text log of a run, also after a failure.

    It holds the version, command line, inputs, effective parameters with their
    source, step timings, warnings, progress messages and final status. With
    log_input_paths = name (the default) the log shows each path by its last
    component, and an input file also by its size; with full it shows the
    absolute paths. The progress log in the work folder always has full paths.
    `files` are input files that `inputs` does not list (the samples of a
    samplesheet); messages about them are shortened too.
    """
    t1 = time.time()
    by_name = ctx.cfg.get("log_input_paths", "name") == "name"
    status = "FAILED" if error else "completed"
    if ctx.status_note:
        status += f" ({ctx.status_note})"
    lines = [
        f"DMRdetector {__version__} {title}",
        f"started:  {datetime.fromtimestamp(started):%Y-%m-%d %H:%M:%S}",
        f"finished: {datetime.fromtimestamp(t1):%Y-%m-%d %H:%M:%S} "
        f"({t1 - started:.0f} s)", f"command:  {COMMAND or ' '.join(sys.argv)}",
        f"settings file: {ctx.config_path or '(none)'}",
        f"versions: {package_versions()}",
        f"threads: {ctx.threads}, memory budget {ctx.mem_gb:.0f} GB",
        f"status:   {status}"
    ]
    if error:
        text = str(error).strip()
        first = text.splitlines()[0] if text else type(error).__name__
        lines.append(f"error:    {first}")
    lines += ["", "[inputs]"] + [
        f"{k} = {_input_text(v, by_name)}" for k, v in inputs.items()
    ]
    lines += ["", "[parameters]  (source: default, file or flag)"]
    for o in settings.OPTIONS:
        if o.key in ctx.cfg:
            source = ctx.sources.get(o.key, "default")
            lines.append(f"{o.key} = {_fmt_value(ctx.cfg[o.key])}  # {source}")
    lines += ["", "[steps]  (seconds)"]
    lines += [
        f"{k}: {v['status']}, {v['seconds']:.1f}"
        for k, v in ctx.timings.items()
    ] or ["(none)"]
    prog = []
    if log.path and os.path.isfile(log.path):
        with open(log.path) as fh:
            prog = fh.read().splitlines()
    warns = [x for x in prog if "WARNING" in x or "NOTICE" in x]
    lines += ["", "[warnings]"] + (warns or ["(none)"])
    if error and len(str(error).strip().splitlines()) > 1:
        lines += ["", "[error]"] + str(error).strip().splitlines()
    lines += ["", "[messages]"] + prog
    text = "\n".join(lines) + "\n"
    if by_name:
        text = redact_paths(text, _path_names(ctx, inputs, files))
    write_atomic(path, text)


CGROUP_ROOT = "/sys/fs/cgroup"
# cgroup v1 reports "no limit" as a number near 2**63
_NO_LIMIT_BYTES = 2**60
# mount names of the v1 controllers below the cgroup root
_V1_SUBDIRS = {"cpu": ("cpu,cpuacct", "cpu"), "memory": ("memory",)}


def _cgroup_levels(root: str, self_cgroup: str, controller: str) -> list[str]:
    """Returns the cgroup folders whose limits apply to this process.

    That is the folder of the process and each parent up to the top of the
    hierarchy, for cgroup v2 (below `root`) and for the v1 mount of
    `controller`. A limit of a parent applies to its children too.
    """
    v2 = v1 = None
    try:
        with open(self_cgroup) as fh:
            for line in fh:
                parts = line.strip().split(":", 2)
                if len(parts) != 3:
                    continue
                if parts[1] == "":
                    v2 = parts[2]
                elif controller in parts[1].split(","):
                    v1 = parts[2]
    except OSError:
        pass
    tops = [(root, v2)] + [
        (os.path.join(root, sub), v1) for sub in _V1_SUBDIRS[controller]
    ]
    out = []
    for top, rel in tops:
        top = os.path.normpath(top)
        d = os.path.normpath(os.path.join(top, (rel or "").lstrip("/")))
        if d != top and not d.startswith(top + os.sep):
            d = top
        while True:
            out.append(d)
            if d == top:
                break
            d = os.path.dirname(d)
    return out


def _read_fields(path: str) -> list[str] | None:
    try:
        with open(path) as fh:
            return fh.read().split()
    except OSError:
        return None


def cgroup_cpu_limit(root: str = CGROUP_ROOT,
                     self_cgroup: str = "/proc/self/cgroup") -> float | None:
    """Returns the CPU quota of the current cgroup in cores, or None.

    It reads cpu.max (cgroup v2), else cpu.cfs_quota_us / cpu.cfs_period_us
    (v1), in the cgroup of the process and in each parent; the smallest quota
    wins.
    """
    limits = []
    for d in _cgroup_levels(root, self_cgroup, "cpu"):
        v = _read_fields(os.path.join(d, "cpu.max"))
        if v and v[0] != "max" and len(v) == 2 and int(v[1]) > 0:
            limits.append(int(v[0]) / int(v[1]))
        q = _read_fields(os.path.join(d, "cpu.cfs_quota_us"))
        per = _read_fields(os.path.join(d, "cpu.cfs_period_us"))
        if q and per and int(q[0]) > 0 and int(per[0]) > 0:
            limits.append(int(q[0]) / int(per[0]))
    return min(limits) if limits else None


def cgroup_memory_limit_gb(
        root: str = CGROUP_ROOT,
        self_cgroup: str = "/proc/self/cgroup") -> float | None:
    """Returns the memory limit of the current cgroup in GB, or None.

    It reads memory.max (cgroup v2) or memory.limit_in_bytes (v1) in the cgroup
    of the process and in each parent; the smallest limit wins. "max" and the
    v1 value for no limit count as none.
    """
    limits = []
    for d in _cgroup_levels(root, self_cgroup, "memory"):
        for name in ("memory.max", "memory.limit_in_bytes"):
            v = _read_fields(os.path.join(d, name))
            if v and v[0].isdigit() and 0 < int(v[0]) < _NO_LIMIT_BYTES:
                limits.append(int(v[0]) / 1e9)
    return min(limits) if limits else None


def meminfo_available_gb(path: str = "/proc/meminfo") -> float | None:
    """Returns MemAvailable of /proc/meminfo in GB, or None."""
    try:
        with open(path) as fh:
            for line in fh:
                if line.startswith("MemAvailable:"):
                    return int(line.split()[1]) * 1024 / 1e9
    except (OSError, ValueError, IndexError):
        pass
    return None


def memory_cap_gb() -> float | None:
    """Returns the memory in GB that this process may use, or None.

    That is the smaller of MemAvailable and the cgroup limit; None means that
    neither is known.
    """
    caps = [
        c for c in (meminfo_available_gb(), cgroup_memory_limit_gb())
        if c is not None
    ]
    return min(caps) if caps else None


def available_threads() -> int:
    """Returns the threads this process may use.

    That is the CPU affinity, lowered to the cgroup quota when that is smaller.
    """
    try:
        n = len(os.sched_getaffinity(0))
    except AttributeError:
        n = os.cpu_count() or 1
    q = cgroup_cpu_limit()
    if q is not None:
        n = min(n, max(1, int(q)))
    return max(1, n)


def resolve_threads(
        value: int | None,
        available: int | None = None,
        mem_cap_gb: float | None = None) -> tuple[int, float, str | None]:
    """Resolves --threads and the memory budget.

    The default is 50% of the available threads. A value below MIN_THREADS is
    raised and a value above the available threads is lowered, each with a
    warning. The memory budget is GB_PER_THREAD x threads, lowered to the memory
    this process may use (`memory_cap_gb`), with a warning.

    Args:
        value: The --threads value, or None.
        available: Threads available (default: `available_threads()`).
        mem_cap_gb: Memory in GB this process may use (default:
            `memory_cap_gb()`; math.inf for no limit).

    Returns:
        (threads, memory budget in GB, warning or None).

    Raises:
        RunError: `value` is below 1.
    """
    if value is not None and value < 1:
        raise RunError(f"threads must be at least {MIN_THREADS}, got {value}")
    avail = available or available_threads()
    cap = memory_cap_gb() if mem_cap_gb is None else mem_cap_gb
    warns = []
    if value is None:
        n = max(MIN_THREADS, int(avail * DEFAULT_THREAD_FRACTION))
    elif value < MIN_THREADS:
        n = MIN_THREADS
        warns.append(f"--threads {value} raised to {MIN_THREADS} (minimum; "
                     f"{MIN_THREADS * GB_PER_THREAD:.0f} GB "
                     "is needed for the largest single task)")
    elif value > max(avail, MIN_THREADS):
        n = max(avail, MIN_THREADS)
        warns.append(f"--threads {value} lowered to {n}: this process may use "
                     f"{avail} thread(s) (CPU affinity and cgroup quota)")
    else:
        n = value
    mem = n * GB_PER_THREAD
    if cap is not None and cap < mem:
        warns.append(f"memory budget lowered from {mem:.0f} GB to {cap:.0f} GB "
                     "(the smaller of MemAvailable and the cgroup memory "
                     "limit)")
        mem = cap
    return n, mem, "; ".join(warns) or None


def log_settings(ctx: Ctx) -> None:
    """Logs the settings file in use, and a NOTICE when it changes a default.

    The file is the one `settings.effective` read (./settings.toml or
    --config); without one the built-in defaults apply.
    """
    log(f"settings: {ctx.config_path or 'built-in defaults'}")
    changed = [
        k for k, src in ctx.sources.items()
        if src == "file" and ctx.cfg.get(k) != settings.BY_KEY[k].default
    ]
    if changed:
        log(f"NOTICE: {ctx.config_path} differs from the packaged defaults "
            f"in {len(changed)} setting(s): {', '.join(changed)}")


def budget_text(ctx: Ctx) -> str:
    """Returns the budget line of the log."""
    return (f"threads: {ctx.threads}  memory budget: {ctx.mem_gb:.0f} GB "
            f"({GB_PER_THREAD:.0f} GB x threads)")


def check_fits(label: str,
               gb: float,
               mem_gb: float,
               threads: int | None = None) -> None:
    """Raises RunError when one task alone exceeds the whole budget.

    With `threads`, a budget below GB_PER_THREAD x threads is reported as a
    limit of the machine, because more threads would not raise it.
    """
    if gb > mem_gb:
        if threads is not None and mem_gb < GB_PER_THREAD * threads - 0.5:
            raise RunError(
                f"{label}: estimated peak memory {gb:.1f} GB exceeds the "
                f"{mem_gb:.0f} GB that this machine (or its cgroup) has "
                "available; free memory or use a machine with more")
        need = int(-(-gb // GB_PER_THREAD))
        raise RunError(
            f"{label}: estimated peak memory {gb:.1f} GB exceeds the budget "
            f"of {mem_gb:.0f} GB ({GB_PER_THREAD:.0f} GB x "
            f"{int(mem_gb / GB_PER_THREAD)} threads); use --threads {need} "
            "or more")


def n_workers(threads: int,
              mem_gb: float,
              gbs: list[float],
              max_workers: int | None = None) -> int:
    """Returns the most tasks that may run together.

    At most `threads` and `max_workers`; the largest estimates of that many
    tasks must sum to at most the memory budget.
    """
    k = max(1, min(threads, len(gbs), max_workers or len(gbs)))
    big = sorted(gbs, reverse=True)
    while k > 1 and sum(big[:k]) > mem_gb:
        k -= 1
    return k


def make_ctx(command: str, cfg: dict, sources: dict, config_path: str | None,
             outdir: Path, workdir: Path) -> Ctx:
    """Returns the Ctx of a command; resolves --threads and its budget."""
    threads, mem, warn = resolve_threads(cfg.get("threads"))
    if warn:
        log(f"WARNING: {warn}")
    # the effective value is recorded in the log
    cfg = dict(cfg, threads=threads)
    return Ctx(command, cfg, outdir, workdir, threads, mem, sources,
               config_path)


def thread_plan(ctx: Ctx,
                gbs: list[float],
                weights: list[float],
                max_workers: int | None = None,
                label: str = "task") -> tuple[int, list[int]]:
    """Plans the workers and the threads of every task.

    Tasks that run together never use more than ctx.threads threads in total nor
    more than ctx.mem_gb of estimated memory. workers = n_workers(...) (at most
    threads / 2); each task gets 2 threads plus its weighted share of the rest,
    computed over the `workers` heaviest tasks.

    Returns:
        (workers, threads of every task).

    Raises:
        RunError: One estimate alone exceeds the budget.
    """
    for g in gbs:
        check_fits(label, g, ctx.mem_gb, ctx.threads)
    # a task's share counts its main thread, which reads, writes and merges
    # beside its pool (measured: up to 4 busy threads in a worker given 3), so a
    # worker needs at least 2: one pool thread and the main thread
    workers = min(n_workers(ctx.threads, ctx.mem_gb, gbs, max_workers),
                  max(1, ctx.threads // 2))
    total = sum(sorted(weights, reverse=True)[:workers])
    return workers, [
        2 + int((ctx.threads - 2 * workers) * w / total) if ctx.threads >= 2 *
        workers else 1 for w in weights
    ]


def killed_worker_error(label: str) -> RunError:
    """Returns the error for a worker process that died without a result.

    The usual cause is the out-of-memory killer; `multiprocessing.Pool` would
    wait for such a worker forever, a ProcessPoolExecutor raises
    BrokenProcessPool.
    """
    return RunError(f"{label}: a worker process was killed, probably out of "
                    "memory; use fewer --threads")


@contextlib.contextmanager
def fork_pool(
    workers: int,
    label: str,
    initializer: Callable | None = None,
    initargs: tuple = ()
) -> Iterator[ProcessPoolExecutor]:
    """Yields a ProcessPoolExecutor of `workers` forked processes.

    The workers are forked when the first task is submitted, so they inherit
    the module state set before that (a read-only mmap, the compiled kernels,
    the module-level tables of a step). The parent should have no other
    threads then. A worker that dies (killed by the out-of-memory killer, a
    signal) turns into the RunError of `killed_worker_error`. On leaving, tasks
    that have not started are cancelled and running ones are waited for.

    Args:
        workers: Number of worker processes.
        label: Name of the step in the error message.
        initializer: Called once in every worker.
        initargs: Arguments of `initializer`.

    Yields:
        The executor.

    Raises:
        RunError: A worker process was killed.
    """
    ex = ProcessPoolExecutor(max_workers=workers,
                             mp_context=mp.get_context("fork"),
                             initializer=initializer,
                             initargs=initargs)
    try:
        yield ex
    except BrokenProcessPool:
        raise killed_worker_error(label) from None
    finally:
        ex.shutdown(wait=True, cancel_futures=True)


def _run_task(label: str, fn, t: dict):
    """Runs fn(t) with DMRDETECTOR_THREADS set to the task's thread share."""
    old = os.environ.get("DMRDETECTOR_THREADS")
    os.environ["DMRDETECTOR_THREADS"] = str(t["threads"])
    try:
        return safe(label, fn, t)
    finally:
        if old is None:
            os.environ.pop("DMRDETECTOR_THREADS", None)
        else:
            os.environ["DMRDETECTOR_THREADS"] = old


def pool_map(ctx: Ctx,
             fn: Callable[[dict], list],
             tasks: list,
             gbs: float | list[float],
             label: str,
             weights: list[float] | None = None,
             max_workers: int | None = None) -> None:
    """Runs fn(task) for every task, in a fork pool when several workers fit.

    Each task dict gets `threads` (see thread_plan) and runs with
    DMRDETECTOR_THREADS set to it, so no step oversubscribes the budget. The
    parent must not have other threads when it forks.

    Args:
        ctx: The run context.
        fn: Takes a task dict, returns a list of (name, status, seconds).
        tasks: Task dicts.
        gbs: Estimated peak memory of one task (GB): a number for all tasks,
            or one per task.
        label: Name of the tasks in error messages.
        weights: Relative thread weights (default: equal).
        max_workers: Upper bound of concurrent tasks.
    """
    gbs = [float(gbs)] * len(tasks) if isinstance(gbs,
                                                  (int, float)) else list(gbs)
    workers, threads = thread_plan(ctx, gbs, weights or [1.0] * len(tasks),
                                   max_workers, label)
    # a task's share of threads includes its own main thread, which reads,
    # writes and merges beside the pool (measured: up to 4 busy threads in a
    # worker given 3), so the pools get one thread less
    tasks = [dict(t, threads=max(1, th - 1)) for t, th in zip(tasks, threads)]
    if workers == 1:
        for t in tasks:
            for s in _run_task(label, fn, t):
                record(ctx, s)
        return
    # the workers put their temporary files in the scratch folder of the run
    with fork_pool(workers, label, set_tmp, (str(ctx.scratch),)) as ex:
        futs = [ex.submit(_run_task, label, fn, t) for t in tasks]
        for f in as_completed(futs):
            for s in f.result():
                record(ctx, s)
