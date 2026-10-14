"""BUILD_ATLAS and SMOOTH_ATLAS (v2.2): the raw and the smoothed atlas.

`build` aggregates the control BEDs of one scope per CpG (N, COV, RAW_METH,
dispersion PHI shrunk toward its median) into a 7-column bgzip BED. `smooth`
convolves the coverage-weighted methylation per chromosome with a Gaussian (FFT)
into the 9-column atlas and its tabix index (doc/ALGORITHM.md).
"""

from __future__ import annotations

import argparse
from collections.abc import Callable
from concurrent.futures import as_completed
import contextlib
import functools
import os
import re
import shutil
import sys
import tempfile
import time

import numpy as np
import pandas as pd

from dmrdetector import runtime
from dmrdetector import settings
from dmrdetector.fileio import atlas_io
from dmrdetector.fileio import bed
from dmrdetector.fileio import bgzf
from dmrdetector.fileio import numparse
from dmrdetector.util import rstats
from dmrdetector.util.jit import jit


WINDOW = 10_000_000  # bp per build task
# rows formatted and compressed per task when writing
CHUNK_ROWS = 250_000
# state shared with forked workers (tabix / pool path)
_G = {}


def _log(msg: str) -> None:
    print(msg, file=sys.stderr, flush=True)


def _natural_key(c: str):
    m = re.match(r"^(?:chr)?(\d+)$", c)
    return (0, int(m.group(1)), c) if m else (1, 0, c)


@jit
def _merge_order(start, end, fo, order, gb):
    """Groups the stacked rows by (START, END) without a sort.

    Files are `fo[i]:fo[i+1]`, each sorted by START. Groups are ordered by
    START, then by first appearance; rows of a group are in stacked order.
    order[:n] receives the row indices group after group and gb[:g + 1] the
    group boundaries. This equals a stable argsort on (START, END) followed by a
    sort of the groups on (START, first row).

    Returns:
        The number of groups g.
    """
    nf = fo.shape[0] - 1
    ptr = fo[:nf].copy()
    o = 0
    g = 0
    gb[0] = 0
    tmp = np.empty(64, np.int64)
    done = False
    while not done:
        best = -1
        for f in range(nf):
            if ptr[f] < fo[f + 1]:
                v = start[ptr[f]]
                if best < 0 or v < best:
                    best = v
        if best < 0:
            break
        # rows with START == best, file order then line order
        m = 0
        for f in range(nf):
            while ptr[f] < fo[f + 1] and start[ptr[f]] == best:
                if m == tmp.shape[0]:
                    t2 = np.empty(m * 2, np.int64)
                    t2[:m] = tmp
                    tmp = t2
                tmp[m] = ptr[f]
                m += 1
                ptr[f] += 1
        # split by END in order of first appearance
        used = np.zeros(m, np.bool_)
        for i in range(m):
            if used[i]:
                continue
            e = end[tmp[i]]
            for j in range(i, m):
                if not used[j] and end[tmp[j]] == e:
                    used[j] = True
                    order[o] = tmp[j]
                    o += 1
            g += 1
            gb[g] = o
    return g


def aggregate_sites(start: np.ndarray,
                    end: np.ndarray,
                    cov: np.ndarray,
                    nmeth: np.ndarray,
                    fo: np.ndarray | None = None) -> tuple[np.ndarray, ...]:
    """Aggregates stacked observations per (START, END).

    R semantics are reproduced: ratio = N_METH / COVERAGE; mean(ratio) in long
    double with its second refinement pass; squared deviations in double; sums
    in long double.

    Args:
        start: START of the rows, in file-then-line order.
        end: END of the rows.
        cov: Coverage of the rows.
        nmeth: Modified reads of the rows.
        fo: Row offsets of the files in the stack (each sorted by START); the
            groups then come from a merge of the files instead of a sort, with
            the same result.

    Returns:
        START, END, N, COV, sum_meth, sum_sq_r, ordered by START then first
        appearance (data.table `by=` order followed by setkey).
    """
    if fo is not None:
        order = np.empty(start.size, np.int64)
        gb = np.empty(start.size + 1, np.int64)
        g = _merge_order(start.astype(np.int64), end.astype(np.int64),
                         np.asarray(fo, np.int64), order, gb)
        b = gb[:g]
        n = np.diff(gb[:g + 1])
        s_g, e_g = start[order[b]].astype(np.int64), end[order[b]].astype(
            np.int64)
        perm = None
    else:
        key = (start.astype(np.int64) << 32) | end.astype(np.int64)
        order = np.argsort(key, kind="stable")
        ks = key[order]
        b = np.flatnonzero(np.r_[True, ks[1:] != ks[:-1]])
        n = np.diff(np.r_[b, ks.size])
        s_g = (ks[b] >> 32).astype(np.int64)
        e_g = (ks[b] & 0xFFFFFFFF).astype(np.int64)
        perm = np.lexsort((order[b], s_g))
    cov_o = cov[order]
    nm_o = nmeth[order]
    ratio = nm_o / cov_o
    ld = np.longdouble
    rl = ratio.astype(ld)
    mean_ld = np.add.reduceat(rl, b) / n.astype(ld)
    mean_ld = mean_ld + np.add.reduceat(rl - np.repeat(mean_ld, n),
                                        b) / n.astype(ld)
    mean = mean_ld.astype(np.float64)
    dev = ratio - np.repeat(mean, n)
    sum_sq = np.add.reduceat((dev * dev).astype(ld), b).astype(np.float64)
    cov_s = np.add.reduceat(cov_o.astype(ld), b).astype(np.float64)
    meth_s = np.add.reduceat(nm_o.astype(ld), b).astype(np.float64)
    if perm is None:
        return s_g, e_g, n, cov_s, meth_s, sum_sq
    return s_g[perm], e_g[perm], n[perm], cov_s[perm], meth_s[perm], sum_sq[
        perm]


def site_stats(n: np.ndarray, cov_s: np.ndarray, meth_s: np.ndarray,
               sum_sq: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """RAW_METH and phi_raw (v2.2 step 6)."""
    raw = meth_s / cov_s
    s2 = sum_sq / (n - 1)
    n_bar = cov_s / n
    mu_var = raw * (1 - raw)
    phi_num = s2 * n_bar - mu_var
    phi_denom = mu_var * (n_bar - 1)
    phi = np.zeros(raw.shape)
    ok = (phi_denom > 1e-9) & (phi_num > 0)
    phi[ok] = phi_num[ok] / phi_denom[ok]
    phi[phi > 1] = 1.0
    return raw, phi


def _handle(path: str):
    import pysam
    h = _G.setdefault("handles", {})
    if path not in h:
        h[path] = pysam.TabixFile(path)
    return h[path]


def _build_window(task):
    chrom, w0, w1 = task
    chunks = []
    for p in _G["files"]:
        tbx = _handle(p)
        if chrom not in tbx.contigs:
            continue
        lines = list(tbx.fetch(chrom, w0, w1))
        if lines:
            chunks.append("\n".join(lines))
    if not chunks:
        return chrom, w0, None
    from io import StringIO
    df = pd.read_csv(StringIO("\n".join(chunks)),
                     sep="\t",
                     header=None,
                     usecols=[1, 2, 3, 4],
                     dtype={
                         1: np.int64,
                         2: np.int64,
                         3: np.float64,
                         4: np.float64
                     },
                     float_precision="round_trip")
    s = df[1].to_numpy()
    m = (s >= w0) & (s < w1)
    if not m.all():
        df = df[m]
        s = df[1].to_numpy()
    if s.size == 0:
        return chrom, w0, None
    g_s, g_e, n, cov_s, meth_s, sum_sq = aggregate_sites(
        s, df[2].to_numpy(), df[3].to_numpy(), df[4].to_numpy())
    keep = n >= _G["min_samples"]
    if not keep.any():
        return chrom, w0, None
    g_s, g_e, n, cov_s, meth_s, sum_sq = (
        a[keep] for a in (g_s, g_e, n, cov_s, meth_s, sum_sq))
    raw, phi_raw = site_stats(n, cov_s, meth_s, sum_sq)
    return chrom, w0, (g_s, g_e, n, cov_s, phi_raw, raw)


def _format_raw_window(task):
    chrom, w0 = task
    s, e, n, cov, phi_raw, raw = _G["res"][(chrom, w0)][0]
    phi = (1 - _G["w"]) * phi_raw + _G["w"] * _G["trend"]
    return bed.write_atlas_text(np.full(s.size, chrom, dtype=object), s, e, n,
                                cov, phi, raw)


_TBI_WINDOW = 1 << 14  # bp per window of the tabix linear index
# compressed bytes added when a window read needs more blocks
_READ_STEP = 1 << 22
_PREAD_BLOCK = 1 << 16


def _pread(fd, n, off):
    parts = []
    while n > 0:
        b = os.pread(fd, n, off)
        if not b:
            break
        parts.append(b)
        n -= len(b)
        off += len(b)
    return b"".join(parts)


class _IndexedBed:
    """A 5-column BGZF BED with its tabix index.

    `rows(chrom, w0, w1)` returns the rows whose START lies in [w0, w1), reading
    only the blocks that hold them; nothing of the file stays in memory.
    """

    def __init__(self, path):
        self.path = path
        self.size = os.path.getsize(path)
        self.ix = bgzf.tbi_index(path)

    def extent(self, chrom):
        e = self.ix.get(chrom)
        return 0 if e is None else e[0].size * _TBI_WINDOW

    def rows(self, chrom, w0, w1):
        e = self.ix.get(chrom)
        if e is None:
            return None
        lo, hi, v_first, v_last = e
        j0 = w0 // _TBI_WINDOW
        if j0 >= lo.size:
            return None
        v_lo = max(int(lo[j0]), v_first)
        j1 = (w1 - 1) // _TBI_WINDOW + 1
        v_hi = v_last if j1 >= hi.size or int(hi[j1]) == np.iinfo(
            np.uint64).max else int(hi[j1])
        pos, skip = v_lo >> 16, v_lo & 0xFFFF
        want = max(0, (v_hi >> 16) - pos) + _PREAD_BLOCK + 4096
        fd = os.open(self.path, os.O_RDONLY)
        try:
            text = np.empty(0, np.uint8)
            while True:
                raw = _pread(fd, want, pos)
                if not raw:
                    raise bgzf.ParseError("BGZF file ends inside a block")
                eof = pos + len(raw) >= self.size
                k = bgzf.whole_blocks_len(raw)
                if k == 0 and not eof:
                    want *= 2
                    continue
                pos += k
                eof = pos >= self.size
                chunk = bgzf.inflate_chunk(raw[:k]) if k else np.empty(
                    0, np.uint8)
                text = chunk if text.size == 0 else np.concatenate(
                    (text, chunk))
                buf = text[skip:] if skip < text.size else text[:0]
                e = buf.size
                while e > 0 and buf[e - 1] != 10:
                    e -= 1
                if e == 0:
                    if eof:
                        return None
                    want = _READ_STEP
                    continue
                r = numparse.parse_buffer(buf[:e], [(1, "i"), (2, "i"),
                                                    (3, "f"), (4, "f")], 1)
                (st, en), (cv, nm), (rs, names) = r
                bounds = np.append(rs, st.size)
                hit = [i for i, c in enumerate(names) if c == chrom]
                if len(hit) > 1:
                    raise bgzf.ParseError("rows not grouped by chromosome")
                if not hit:
                    if names or eof:
                        return None
                    want = _READ_STEP
                    continue
                a, b = int(bounds[hit[0]]), int(bounds[hit[0] + 1])
                if eof or hit[0] < len(names) - 1 or (b > a and
                                                      st[b - 1] >= w1):
                    break
                want = _READ_STEP
        finally:
            os.close(fd)
        s = st[a:b]
        if s.size > 1 and (np.diff(s) < 0).any():
            raise bgzf.ParseError("rows not sorted by START")
        i0, i1 = np.searchsorted(s, [w0, w1])
        if i1 <= i0:
            return None
        sl = slice(a + int(i0), a + int(i1))
        return st[sl].astype(np.int32), en[sl].astype(
            np.int32), cv[sl].copy(), nm[sl].copy()


def _site_window(parts):
    """Aggregates the rows of one window (one tuple per file).

    Each tuple is (start, end, cov, nmeth).
    """
    parts = [p for p in parts if p is not None]
    if not parts:
        return None
    fo = np.cumsum([0] + [p[0].size for p in parts])
    s = np.concatenate([p[0] for p in parts]).astype(np.int64)
    return aggregate_sites(
        s,
        np.concatenate([p[1] for p in parts]).astype(np.int64),
        np.concatenate([p[2] for p in parts]),
        np.concatenate([p[3] for p in parts]), fo)


@contextlib.contextmanager
def _big_heap():
    """Serves the MB-sized buffers of the window workers from the heap (glibc).

    With the default mmap threshold every buffer is an mmap/munmap pair and the
    page faults of 24 threads contend on the process mmap lock (44 s instead of
    20 s for the raw atlas of 14 controls; see `atlas_io.tune_malloc`).
    Afterwards the thresholds go back to glibc's defaults and freed memory is
    returned to the system. Best effort; other C libraries are left alone.
    """
    libc = atlas_io.tune_malloc()
    try:
        yield
    finally:
        if libc is not None:
            libc.mallopt(-3, 128 * 1024)
            libc.mallopt(-1, 128 * 1024)
            libc.malloc_trim(0)


def _is_bgzf(path):
    with open(path, "rb") as fh:
        return fh.read(16) == bgzf.MAGIC[:16]


def _build_fast(files, min_samples, threads, window, mem_gb, scratch, log, t0):
    """Aggregates the BGZF inputs window by window.

    The windows are read through the tabix indexes (made in `scratch` for inputs
    without one). `window` and the number of concurrent windows come from the
    memory budget `mem_gb` (runtime.plan_raw_atlas) unless `window` is given.

    Returns:
        {(chrom, w0): arrays}, or None when an input is not BGZF.
    """
    from concurrent.futures import ThreadPoolExecutor
    if not all(_is_bgzf(f) for f in files):
        return None
    # windows are read by the index: a stale one would give wrong rows
    beds = [
        _IndexedBed(bed.ensure_tabix(f, scratch, check_age=True)) for f in files
    ]
    chroms = sorted({c for b in beds for c in b.ix}, key=_natural_key)
    extent = {c: max(b.extent(c) for b in beds) for c in chroms}
    total_bp = sum(extent.values())
    workers = max(1, threads)
    if window is None:
        in_bytes = sum(b.size for b in beds)
        plan = runtime.plan_raw_atlas(
            in_bytes, total_bp, min_samples, workers,
            mem_gb if mem_gb is not None else float("inf"))
        if plan is None:
            need = runtime.raw_atlas_min_gb(in_bytes, total_bp, min_samples)
            floor = runtime.raw_atlas_floor_gb(in_bytes, total_bp, min_samples)
            raise SystemExit(
                f"raw atlas of {len(beds)} samples: needs at least {need:.1f} "
                f"GB (results kept in memory {floor:.1f} GB plus one window "
                f"of {runtime.RAW_WINDOW_MIN:,} bp), more than the "
                f"{mem_gb:.1f} GB available; raise --threads (8 GB each)")
        window, workers = plan
    log(f"{len(beds)} files; windows of {window:,} bp, {workers} concurrent" +
        ("" if mem_gb is None else f" (memory budget {mem_gb:.0f} GB)"))
    tasks = [
        (c, w0) for c in chroms for w0 in range(0, max(extent[c], 1), window)
    ]

    def agg(t):
        c, w0 = t
        r = _site_window([b.rows(c, w0, w0 + window) for b in beds])
        if r is None:
            return c, w0, None
        g_s, g_e, n, cov_s, meth_s, sum_sq = r
        keep = n >= min_samples
        if not keep.any():
            return c, w0, None
        g_s, g_e, n, cov_s, meth_s, sum_sq = (
            a[keep] for a in (g_s, g_e, n, cov_s, meth_s, sum_sq))
        raw, phi_raw = site_stats(n, cov_s, meth_s, sum_sq)
        return c, w0, (g_s, g_e, n, cov_s, phi_raw, raw)

    res = {}
    with _big_heap(), ThreadPoolExecutor(workers) as ex:
        for c, w0, r in ex.map(agg, tasks):
            if r is not None:
                res[(c, w0)] = r
    log(f"Aggregated {len(res)} windows ({time.time() - t0:.0f} s)")
    return res


def _regroup(res, grid=WINDOW):
    """Cuts window results on the fixed grid of `grid` bp.

    The blocks of the output file follow this grid, so they do not depend on the
    windows the sites were aggregated in.

    Returns:
        {(chrom, cell start): [arrays, ...]} with the pieces of a cell in START
        order (views of the window arrays; nothing is copied).
    """
    cells = {}
    for (c, w0) in sorted(res, key=lambda k: (k[0], k[1])):
        arrs = res[(c, w0)]
        cell = arrs[0] // grid
        cut = np.flatnonzero(np.diff(cell)) + 1
        a = 0
        for b in list(cut) + [arrs[0].size]:
            cells.setdefault((c, int(cell[a]) * grid),
                             []).append(tuple(x[a:b] for x in arrs))
            a = b
        res[(c, w0)] = None
    return cells


def _cell_arrays(pieces):
    """Returns the six arrays of one grid cell.

    They are copied only when the cell is split over several windows.
    """
    return pieces[0] if len(pieces) == 1 else tuple(
        np.concatenate([p[i] for p in pieces]) for i in range(6))


def build_atlas(bed_files: list[str],
                out: str,
                min_samples: int = settings.ATLAS_MIN_SAMPLES_PER_SITE,
                shrinkage_weight: float = settings.ATLAS_SHRINKAGE_WEIGHT,
                threads: int = 1,
                scratch: str | None = None,
                window: int | None = None,
                log: Callable[[str], None] = _log,
                mem_gb: float | None = None) -> float:
    """Builds the raw atlas (v2.2 BUILD_ATLAS).

    BGZF inputs are read window by window through their tabix index, inflated
    and parsed in memory by numba kernels. `mem_gb` is the memory this call may
    use: the window size and the number of windows aggregated at once follow
    from it and from the number of samples (`runtime.plan_raw_atlas`); `window`
    fixes the window in bp instead. Neither changes the output: the file is
    byte-identical for every budget and window. Other inputs, a file the fast
    reader cannot handle, and DMRDETECTOR_NO_JIT=1 use the tabix and pandas path
    (`_build_atlas_tabix`; `window` is then its grid, default WINDOW).

    Args:
        bed_files: 5-column bgzip BEDs (CHROM START END COVERAGE N_METH).
        out: Output file (7 columns, bgzip).
        min_samples: Keep a site covered by at least this many samples.
        shrinkage_weight: Weight of the median dispersion in the shrunken PHI.
        threads: Worker threads.
        scratch: Folder for indexes of inputs without a .tbi.
        window: Window in bp, instead of the one derived from `mem_gb`.
        log: Receives progress messages.
        mem_gb: Memory budget in GB (default: no limit).

    Returns:
        The phi trend.
    """
    t0 = time.time()
    files = sorted(bed_files, key=lambda p: os.path.basename(p))
    files = [f for f in files if os.path.getsize(f) > 50]
    if not files:
        raise SystemExit("No valid input files found")
    log(f"Found {len(files)} valid input files.")
    # the folder of the indexes is ours, and removed, when the caller gave none
    own_scratch = not scratch
    if own_scratch:
        scratch = tempfile.mkdtemp(prefix="dmrdetector_atlas_")
    try:
        res = None
        if not bgzf.NO_JIT:
            try:
                res = _build_fast(files, min_samples, threads, window, mem_gb,
                                  scratch, log, t0)
            except bgzf.ParseError as e:
                log(f"fast reader not used ({e}); reading with tabix")
        if res is None:
            return _build_atlas_tabix(files, out, min_samples, shrinkage_weight,
                                      threads, scratch, window or WINDOW, log,
                                      t0)
        res = _regroup(res)
        return _build_finish(res, out, shrinkage_weight, threads, log, t0)
    finally:
        if own_scratch:
            shutil.rmtree(scratch, ignore_errors=True)


def _summarize(res, log):
    """(number of sites, phi trend, ordered cell keys) of window results."""
    if not res:
        raise SystemExit(
            "no CpG site passed the atlas filters (at least "
            "`atlas_min_samples_per_site` control samples with read coverage "
            "above `atlas_min_read_coverage`); use more control samples or "
            "check that the input BED files are hg38 CpG methylation tables")
    nsite = sum(p[0].size for r in res.values() for p in r)
    log(f"Retained {nsite} CpG sites after filtering.")
    allphi = np.concatenate([p[4] for r in res.values() for p in r])
    valid = allphi[allphi > 0]
    trend = 0.01 if valid.size == 0 else float(np.median(valid))
    log(f"phi_trend (median): {trend:.6f} (full precision {trend!r})")
    del allphi, valid
    cnames = bed.c_order(c for c, _ in res)
    keys = [
        (c, w) for c in cnames for w in sorted(w for cc, w in res if cc == c)
    ]
    return nsite, trend, keys


def _build_finish(res, out, shrinkage_weight, threads, log, t0):
    nsite, trend, keys = _summarize(res, log)

    def fmt(k):
        s, e, n, cov, phi_raw, raw = _cell_arrays(res[k])
        phi = (1 - shrinkage_weight) * phi_raw + shrinkage_weight * trend
        parts = []
        for a in range(0, s.size, CHUNK_ROWS):
            b = min(s.size, a + CHUNK_ROWS)
            ints = np.stack([
                s[a:b], e[a:b], n[a:b].astype(np.int64),
                cov[a:b].astype(np.int64)
            ])
            floats = np.stack([phi[a:b], raw[a:b]])
            parts.extend(
                bgzf.compress_bytes(numparse.format_rows(k[0], ints, floats),
                                    bgzf.ARCHIVE_LEVEL))
        res[k] = None
        return parts

    _write_ordered(out, keys, fmt, threads)
    log(f"Saved {nsite} sites to {out} ({time.time() - t0:.0f} s)")
    return trend


def _write_ordered(out, keys, fn, threads):
    """Runs fn(key) -> BGZF blocks in a thread pool; writes them in order."""
    from concurrent.futures import ThreadPoolExecutor
    with bgzf.BgzfWriter(out) as w:
        with ThreadPoolExecutor(max(1, threads)) as ex:
            futs = [ex.submit(fn, k) for k in keys]
            for f in futs:
                w.write_blocks(f.result())
                f._result = None


def _build_atlas_tabix(files, out, min_samples, shrinkage_weight, threads,
                       scratch, window, log, t0):
    """Builds the raw atlas with tabix region reads and pandas.

    This is the path of NO_JIT and the fallback of `build_atlas` for a file the
    fast reader cannot parse.
    """
    files = [bed.ensure_tabix(f, scratch) for f in files]
    ext: dict[str, int] = {}
    for f in files:
        for c, v in bgzf.tbi_extents(f).items():
            ext[c] = max(ext.get(c, 0), v)
    chroms = sorted(ext, key=_natural_key)
    tasks = [(c, w0, w0 + window)
             for c in chroms
             for w0 in range(0, max(ext[c], 1), window)]
    _G.update(files=files, min_samples=min_samples, handles={})
    res = {}
    with runtime.fork_pool(max(1, threads), "atlas") as ex:
        futs = [ex.submit(_build_window, t) for t in tasks]
        for f in as_completed(futs):
            chrom, w0, r = f.result()
            if r is not None:
                res[(chrom, w0)] = [r]
    nsite, trend, keys = _summarize(res, log)
    _G.update(res=res, w=shrinkage_weight, trend=trend)
    # forked after the update above: the workers read `res`, `w` and `trend`
    with runtime.fork_pool(max(1, threads), "atlas") as ex:
        bed.write_bgzip_text(out, ex.map(_format_raw_window, keys))
    log(f"Saved {nsite} sites to {out} ({time.time() - t0:.0f} s)")
    return trend


@functools.lru_cache(maxsize=4)
def _half_kernel(sigma: int) -> np.ndarray:
    """Returns exp(-x^2 / (2 sigma^2)) for x = 0..ceil(3 sigma)."""
    radius = int(np.ceil(3 * sigma))
    x = np.arange(0, radius + 1, dtype=np.float64)
    # x * x and 2 sigma^2 are exact integers in float64; the division is one
    # IEEE operation, as in R (-(x^2) / (2 * sigma^2))
    arg = -(x * x) / (2.0 * (float(sigma) * float(sigma)))
    half = np.array([rstats.exp_rounded(float(a)) for a in arg])
    half.setflags(write=False)
    return half


def gaussian_kernel(sigma: int) -> np.ndarray:
    """Returns the Gaussian kernel of radius ceil(3 sigma), length 2 radius + 1.

    The entries are exp(-x^2 / (2 sigma^2)), x = -radius..radius, each rounded
    once to the nearest float64 (`rstats.exp_rounded`), so the kernel, and with
    it the smoothed atlas, does not depend on the SIMD loops of numpy or on the
    libm of the host (np.exp differs by one ulp between them in about 5% of the
    entries). The kernel is symmetric and not normalized. The table of one
    sigma is computed once per process.
    """
    half = _half_kernel(int(sigma))
    return np.concatenate((half[:0:-1], half))


def smooth_chromosome(start: np.ndarray,
                      end: np.ndarray,
                      cov: np.ndarray,
                      raw_meth: np.ndarray,
                      phi: np.ndarray,
                      kernel: np.ndarray,
                      inner: int = 1,
                      child: bool = False) -> tuple[np.ndarray, np.ndarray]:
    """REF_METH and REF_VAR for one chromosome (v2.2 `process_chromosome`).

    FFT convolution on the R `nextn` length of span + kernel - 1 (zero padded),
    done by `util.fftsmooth.smooth_core` (same scipy.fft.rfft(v, n)/irfft calls;
    with inner > 1 the independent transforms run in 3 threads, bit-identical to
    inner = 1). With `child` the transforms run in a new process that exits
    afterwards, so the plan cache of scipy.fft (about 8 bytes per padded base
    and length, kept until the process ends) does not accumulate over the
    chromosomes; the numbers are the same.
    """
    from dmrdetector.util import fftsmooth
    min_pos, max_pos = int(start.min()), int(end.max())
    span = max_pos - min_pos + 1
    idx = start - min_pos
    fft_len = rstats.nextn(span + kernel.size - 1)
    fn = fftsmooth.smooth_in_child if child else fftsmooth.smooth_core
    ref = fn(idx, raw_meth * cov, cov, kernel, fft_len, span, inner)
    var = ref * (1 - ref) * phi
    return ref, var


def _smooth_task_pool(chrom, fft_child):
    d = _G["chr"][chrom]
    ref, var = smooth_chromosome(d["START"],
                                 d["END"],
                                 d["COV"],
                                 d["RAW_METH"],
                                 d["PHI"],
                                 _G["kernel"],
                                 child=fft_child)
    return bed.write_atlas_text(np.full(ref.size, chrom, dtype=object),
                                d["START"], d["END"], d["N"], d["COV"],
                                d["PHI"], d["RAW_METH"], ref, var)


def prepare_atlas(df: pd.DataFrame, include_sex: bool) -> pd.DataFrame:
    """Normalizes chromosome names and filters as v2.2.

    The result is sorted by CHROM, START (`setkey(CHROM, START)`).
    """
    df = df.copy()
    df["CHROM"] = bed.normalize_primary_chrom(
        df["CHROM"].to_numpy(dtype=object)).to_numpy(dtype=object)
    allowed = {f"chr{i}" for i in range(1, 23)
              } | ({"chrX", "chrY"} if include_sex else set())
    df = df[df["CHROM"].isin(allowed)]
    if df.empty:
        raise SystemExit(f"No allowed chromosomes found "
                         f"(include_sex_chromosomes={include_sex})")
    return df.iloc[bed.bed_sort_order(
        df["CHROM"].to_numpy(dtype=object),
        df["START"].to_numpy())].reset_index(drop=True)


def _load_raw_atlas(path, include_sex, threads):
    """Loads a raw atlas (7 columns, BGZF) per chromosome.

    Returns:
        {chrom: dict of arrays} in C order of the normalized names, after the
        v2.2 filter and `setkey(CHROM, START)` order; None when the file is not
        BGZF.
    """
    r = numparse.read_table(path, [(1, "i"), (2, "i"), (3, "i"), (4, "i"),
                                   (5, "f"), (6, "f")], threads)
    if r is None:
        return None
    (st, en, nn, cv), (ph, rm), (rs, names) = r
    allowed = {f"chr{i}" for i in range(1, 23)
              } | ({"chrX", "chrY"} if include_sex else set())
    norm = bed.normalize_primary_chrom(names).to_numpy(dtype=object)
    bounds = np.append(rs, st.size)
    byname = {}
    for i, c in enumerate(norm):
        if c in allowed:
            byname.setdefault(c, []).append(
                (int(bounds[i]), int(bounds[i + 1])))
    out = {}
    for c in bed.c_order(byname):
        rng = byname[c]
        sl = [np.arange(a, b) for a, b in rng] if len(rng) > 1 else None
        if sl is None:
            a, b = rng[0]
            d = dict(START=st[a:b],
                     END=en[a:b],
                     N=nn[a:b],
                     COV=cv[a:b].astype(np.float64),
                     PHI=ph[a:b],
                     RAW_METH=rm[a:b])
            if b - a > 1 and (np.diff(d["START"]) < 0).any():
                o = np.argsort(d["START"], kind="stable")
                d = {k: v[o] for k, v in d.items()}
        else:
            ix = np.concatenate(sl)
            ix = ix[np.argsort(st[ix], kind="stable")]
            d = dict(START=st[ix],
                     END=en[ix],
                     N=nn[ix],
                     COV=cv[ix].astype(np.float64),
                     PHI=ph[ix],
                     RAW_METH=rm[ix])
        out[c] = d
    return out


def _smooth_chunk(chrom, d, ref, var, a, b):
    """Returns the compressed text of rows a..b, as bgzf.compress_rows.

    The result is (blocks, rows per block, in-block offset of every row).
    """
    ints = np.stack([
        d["START"][a:b], d["END"][a:b], d["N"][a:b],
        d["COV"][a:b].astype(np.int64)
    ])
    floats = np.stack([d["PHI"][a:b], d["RAW_METH"][a:b], ref[a:b], var[a:b]])
    buf = np.frombuffer(numparse.format_rows(chrom, ints, floats), np.uint8)
    rowoff = np.concatenate(
        ([0], np.flatnonzero(buf == 10) + 1)).astype(np.int64)
    return bgzf.compress_rows(buf, rowoff, rowoff.size - 1, bgzf.ARCHIVE_LEVEL)


class MemGate:
    """Admits FFT tasks in submission order within a memory total.

    The summed estimates of the admitted tasks stay within `total` GB; an
    estimate larger than `total` is admitted alone (the caller has checked that
    it fits the budget).
    """

    def __init__(self, total: float):
        """Sets the memory total in GB."""
        import threading
        self.total, self.used, self.next = total, 0.0, 0
        self.cv = threading.Condition()

    def acquire(self, order: int, gb: float) -> None:
        """Waits for the turn and the memory of task `order` (`gb` GB)."""
        with self.cv:
            self.cv.wait_for(lambda: self.next == order and
                             (self.used == 0 or self.used + gb <= self.total))
            self.used += gb
            self.next += 1
            self.cv.notify_all()

    def release(self, gb: float) -> None:
        """Returns the `gb` GB of a finished task."""
        with self.cv:
            self.used -= gb
            self.cv.notify_all()


def _smooth_one(chrom, d, kernel, ex, inner=1, gate=None, order=0, child=True):
    """Smooths one chromosome and starts the formatting of its row chunks.

    Formatting and compression go to the pool `ex`. With a `gate` the FFT waits
    for its memory share; `child` runs each FFT in its own short-lived process
    (`smooth_chromosome`).

    Returns:
        The chunk futures in row order.
    """
    if gate is not None:
        gb = runtime.smooth_gb(_span(d), inner)
        gate.acquire(order, gb)
    try:
        ref, var = smooth_chromosome(d["START"], d["END"], d["COV"],
                                     d["RAW_METH"], d["PHI"], kernel, inner,
                                     child)
    finally:
        if gate is not None:
            gate.release(gb)
    n = ref.size
    return [
        ex.submit(_smooth_chunk, chrom, d, ref, var, a, min(n, a + CHUNK_ROWS))
        for a in range(0, n, CHUNK_ROWS)
    ]


def _span(d) -> int:
    """Returns the covered span in bp of one chromosome's arrays."""
    return int(d["END"].max()) - int(d["START"].min())


def _smooth_loaded(chr_data, sigma, out, include_sex, threads, log, t0, mem_gb,
                   fft_child, index):
    """Smooths a raw atlas loaded by `_load_raw_atlas` and writes `out`."""
    log(f"=== FFT Gaussian Smoothing (sigma = {sigma} bp) ===")
    if not chr_data:
        raise SystemExit(f"No allowed chromosomes found "
                         f"(include_sex_chromosomes={include_sex})")
    names = list(chr_data)
    log(f"Chromosome filter: include_sex_chromosomes={include_sex} "
        f"({len(names)} chroms retained)")
    n = sum(len(v["START"]) for v in chr_data.values())
    log(f"Loaded {n} CpG sites.")
    kernel = gaussian_kernel(sigma)
    span = {c: _span(chr_data[c]) for c in names}
    work = sorted(names, key=lambda c: -span[c])
    # The independent FFTs of one chromosome run in 3 threads (same numbers as
    # 1): for every chromosome when they all fit in the thread budget, else
    # for the two longest, whose FFTs set the wall time. The 2 extra threads of
    # each of those two come out of the pool, so the busy threads never exceed
    # `threads`.
    fit = 3 * len(names) <= threads
    big = set(work[:2]) if threads >= 8 and not fit else set()

    def inner(c):
        return 3 if fit or c in big else 1

    nproc = max(1, min(threads - 2 * len(big), len(names)))
    log(f"Processing {len(names)} chromosomes using {nproc} cores...")
    gate = None
    inner_room = inner
    if mem_gb is not None:
        # arrays, results, parent base
        held = (
            sum(a.nbytes for v in chr_data.values() for a in v.values()) / 1e9 +
            16 * n / 1e9 + runtime.SMOOTH_PARENT_BASE_GB)
        room = mem_gb - held
        need = max(runtime.smooth_gb(span[c]) for c in names)
        if need > room:
            raise RuntimeError(
                f"smoothing needs about {need:.1f} GB for the largest "
                f"chromosome plus {held:.1f} GB "
                f"of loaded data, more than the {mem_gb:.1f} GB available")
        gate = MemGate(room)

        def inner_room(c):
            # 3 FFT threads need about 1.8 times the memory of 1 (same
            # numbers); fall back to 1 when they do not fit
            if inner(c) > 1 and runtime.smooth_gb(span[c], 3) <= room:
                return 3
            return 1

    order = {c: i for i, c in enumerate(work)}

    def task(c, ex):
        return _smooth_one(c, chr_data[c], kernel, ex, inner_room(c), gate,
                           order[c], fft_child)

    _write_chrom_ordered(out, names, work, task, nproc,
                         chr_data if index else None)
    log(f"Written {n} sites to {out} ({time.time() - t0:.0f} s)")


def smooth_atlas(raw_atlas: str,
                 sigma: int,
                 out: str,
                 include_sex: bool = True,
                 threads: int = 1,
                 log: Callable[[str], None] = _log,
                 index: bool = True,
                 mem_gb: float | None = None,
                 fft_child: bool = True) -> None:
    """Smooths a raw atlas into the reference atlas (v2.2 SMOOTH_ATLAS).

    The output is BGZF (threaded compression) and, with `index`, gets its tabix
    index `out`.tbi in the same step. A raw atlas the fast reader cannot handle
    (not BGZF, odd tokens, DMRDETECTOR_NO_JIT=1) goes through the pandas path
    (`_smooth_atlas_pool`). With `fft_child` the FFT of every chromosome runs in
    a short-lived process, which frees the scipy.fft plan cache of its length
    when it exits (the plans of all chromosomes together are about 8 bytes per
    padded base, 23 GB for a human genome); `fft_child=False` runs them in this
    process, with the same numbers.

    Args:
        raw_atlas: Raw atlas (7 columns).
        sigma: Gaussian sigma in bp.
        out: Output atlas (9 columns, bgzip).
        include_sex: Keep chrX and chrY.
        threads: Worker threads.
        log: Receives progress messages.
        index: Write the tabix index.
        mem_gb: Memory budget in GB; chromosomes are smoothed together only
            while
            their summed estimates (`runtime.smooth_gb`, on top of the loaded
            arrays) fit in it. The numbers do not depend on it.
        fft_child: Run each FFT in a child process.
    """
    t0 = time.time()
    chr_data = None
    if not bgzf.NO_JIT:
        try:
            chr_data = _load_raw_atlas(raw_atlas, include_sex, threads)
        except bgzf.ParseError as e:
            log(f"fast reader not used ({e}); reading with pandas")
    if chr_data is not None:
        _smooth_loaded(chr_data, sigma, out, include_sex, threads, log, t0,
                       mem_gb, fft_child, index)
        return
    _smooth_atlas_pool(raw_atlas, sigma, out, include_sex, threads, log, t0,
                       mem_gb, fft_child)
    if index:
        bed.tabix_bed(out)


def _write_chrom_ordered(out, names, work, fn, nproc, chr_data=None):
    """Submits fn(chrom, pool) in `work` order and writes the chunks in order.

    `work` is largest first; fn returns the chunk futures of compress_rows
    results. Chunks are written in `names` order. With `chr_data` the tabix
    index `out`.tbi is built from the virtual offsets the writer reports for
    every row (no second read).
    """
    from concurrent.futures import ThreadPoolExecutor
    contigs = []
    with bgzf.BgzfWriter(out) as w:
        with ThreadPoolExecutor(nproc) as ex:
            futs = {c: ex.submit(fn, c, ex) for c in work}
            for c in names:
                voff = []
                for f in futs[c].result():
                    voff.append(w.add(*f.result()))
                if chr_data is not None:
                    contigs.append((chr_data[c]["START"], chr_data[c]["END"],
                                    np.concatenate(voff)))
                futs[c] = None
        final = w.end
    if chr_data is not None and not bgzf.write_tbi(out, names, contigs, final):
        bed.tabix_bed(out)


def _smooth_atlas_pool(raw_atlas,
                       sigma,
                       out,
                       include_sex,
                       threads,
                       log,
                       t0,
                       mem_gb=None,
                       fft_child=True):
    """First v3 release: pandas read, forked workers, serial BGZF write.

    With `fft_child` the FFT of every chromosome runs in a short-lived process
    of its own (as in `_smooth_one`), so the scipy.fft plan cache of one length
    is freed when that process exits; a ProcessPoolExecutor cannot replace its
    forked workers after each task.
    """
    log(f"=== FFT Gaussian Smoothing (sigma = {sigma} bp) ===")
    df = pd.read_csv(raw_atlas,
                     sep="\t",
                     header=None,
                     names=bed.RAW_ATLAS_COLS[:7],
                     usecols=range(7),
                     dtype={
                         "CHROM": str,
                         "START": np.int64,
                         "END": np.int64
                     },
                     float_precision="round_trip")
    df = prepare_atlas(df, include_sex)
    chroms = bed.c_order(df["CHROM"])
    log(f"Chromosome filter: include_sex_chromosomes={include_sex} "
        f"({len(chroms)} chroms retained)")
    log(f"Loaded {len(df)} CpG sites.")
    codes, names = bed.chrom_codes(df["CHROM"].to_numpy(dtype=object))
    bounds = np.searchsorted(codes, np.arange(len(names) + 1))
    chr_data = {}
    for i, c in enumerate(names):
        sl = slice(bounds[i], bounds[i + 1])
        chr_data[c] = {
            k: df[k].to_numpy()[sl]
            for k in ("START", "END", "N", "COV", "PHI", "RAW_METH")
        }
        chr_data[c]["COV"] = chr_data[c]["COV"].astype(np.float64)
    del df
    _G.update({"chr": chr_data, "kernel": gaussian_kernel(sigma)})
    # largest chromosomes first; at most `threads` run at once
    work = sorted(names, key=lambda c: -_span(chr_data[c]))
    nproc = max(1, min(threads, len(names)))
    # the largest estimates of the chromosomes smoothed together must fit
    if mem_gb is not None:
        est = [runtime.smooth_gb(_span(chr_data[c])) for c in work]
        room = mem_gb - sum(v["START"].nbytes * 6
                            for v in chr_data.values()) / 1e9
        if est[0] > room:
            raise RuntimeError(
                f"smoothing needs about {est[0]:.1f} GB for the largest "
                "chromosome, more than the "
                f"{room:.1f} GB available")
        while nproc > 1 and sum(est[:nproc]) > room:
            nproc -= 1
    log(f"Processing {len(names)} chromosomes using {nproc} cores...")
    with runtime.fork_pool(nproc, "smoothing") as ex:
        futs = {c: ex.submit(_smooth_task_pool, c, fft_child) for c in work}
        # a result is dropped once it is written
        bed.write_bgzip_text(out, (futs.pop(c).result() for c in names))
    n = sum(len(v["START"]) for v in chr_data.values())
    log(f"Written {n} sites to {out} ({time.time() - t0:.0f} s)")


def _main_build(args) -> int:
    """Runs `atlas build` for the parsed arguments; returns 0."""
    out = args.out or "raw_atlas.bed.gz"
    build_atlas(args.beds,
                out,
                args.min_samples_per_site,
                args.shrinkage_weight,
                args.threads,
                args.tmpdir,
                window=args.window,
                mem_gb=args.mem_gb)
    return 0


def _main_smooth(args) -> int:
    """Runs `atlas smooth` for the parsed arguments; returns 0."""
    out = args.out or (f"ref_atlas_{args.sigma}.bed.gz"
                       if args.scope == "autosomal" else
                       f"ref_atlas_{args.sigma}_{args.scope}.bed.gz")
    smooth_atlas(args.raw_atlas,
                 args.sigma,
                 out,
                 args.include_sex_chromosomes,
                 args.threads,
                 index=args.index)
    return 0


def main(args: argparse.Namespace) -> int:
    """Runs the `atlas build` or `atlas smooth` subcommand."""
    return args.subfunc(args)


def add_parser(subparsers: argparse._SubParsersAction) -> None:
    """Adds the `atlas` command with `build` and `smooth` to `subparsers`."""
    p = subparsers.add_parser("atlas", help="build the raw atlas and smooth it")
    p.set_defaults(func=main)
    sub = p.add_subparsers(dest="atlas_cmd", required=True)

    b = sub.add_parser(
        "build",
        help="aggregate CTRL BEDs into the raw atlas (empirical-Bayes PHI)")
    b.add_argument(
        "beds",
        nargs="+",
        help="per-sample 5-column .bed.gz of one scope (.tbi used if present)")
    b.add_argument("--out",
                   help="default raw_atlas.bed.gz (7 columns, headerless)")
    b.add_argument("--min-samples-per-site",
                   type=int,
                   default=settings.ATLAS_MIN_SAMPLES_PER_SITE,
                   help="atlas_min_samples_per_site")
    b.add_argument("--shrinkage-weight",
                   type=float,
                   default=settings.ATLAS_SHRINKAGE_WEIGHT,
                   help="atlas_shrinkage_weight")
    b.add_argument("--threads", type=int, default=1)
    b.add_argument("--tmpdir",
                   help="scratch for indexes of inputs lacking a .tbi")
    b.add_argument(
        "--mem-gb",
        type=float,
        help="memory budget (GB): sets the window size and the number of "
        "windows "
        "aggregated at once; the output does not depend on it (default: no "
        "limit)")
    b.add_argument(
        "--window",
        type=int,
        help="window size in bp instead of the one derived from --mem-gb")
    b.set_defaults(subfunc=_main_build)

    s = sub.add_parser(
        "smooth",
        help=
        "Gaussian-smooth the raw atlas (FFT), write ref_atlas_<sigma>.bed.gz")
    s.add_argument("--raw-atlas", required=True)
    s.add_argument(
        "--smooth-sigma",
        "--sigma",
        dest="sigma",
        type=int,
        default=settings.SMOOTH_SIGMA,
        help="Gaussian sigma in bp (v2.2 smooth_sigma_narrow; the wide sigma "
        "is "
        "not built in v3)")
    s.add_argument("--scope",
                   default="autosomal",
                   choices=["autosomal", "male_X", "female_X", "Y"])
    s.add_argument("--out", help="default ref_atlas_<sigma>[_<scope>].bed.gz")
    s.add_argument("--include-sex-chromosomes",
                   action=argparse.BooleanOptionalAction,
                   default=True)
    s.add_argument("--threads", type=int, default=1)
    s.add_argument("--index",
                   action=argparse.BooleanOptionalAction,
                   default=True,
                   help="write the tabix index <out>.tbi")
    s.set_defaults(subfunc=_main_smooth)
