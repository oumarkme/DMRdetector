"""Atlas files: table readers and the fast panel reader of CpG island and ICR.

The panel reader tokenises only the columns the panels need (CHROM, START, END,
N, PHI, REF_METH), filters rows while it parses and inflates BGZF blocks in a
thread pool. Input it cannot read raises `Unsupported`; callers use pandas.
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
import mmap
import os
import zlib

import numpy as np
import pandas as pd

from dmrdetector.fileio import bed
from dmrdetector.fileio import bgzf
from dmrdetector.fileio import numparse
from dmrdetector.util.jit import jit


ATLAS_COLS = [
    "CHROM", "START", "END", "N", "COV", "PHI", "RAW_METH", "REF_METH",
    "REF_VAR"
]


def read_atlas_pandas(path: str) -> pd.DataFrame:
    """Reads a headerless atlas BED (gzip or plain) with pandas.

    It has 9 columns, or the first ones of a raw atlas. Decimals are parsed
    with `float_precision="round_trip"`: the default parser of pandas can lose
    the last digits, and the byte-level reader of the fast path is correctly
    rounded, so only the exact parse makes both routes give the same bits.
    """
    with bed.open_text(path) as fh:
        ncols = len(fh.readline().rstrip("\n").split("\t"))
    return pd.read_csv(path,
                       sep="\t",
                       header=None,
                       names=ATLAS_COLS[:ncols],
                       dtype={"CHROM": str},
                       na_values=["NA"],
                       float_precision="round_trip")


class LazyAtlas:
    """Path of an atlas that `detect` reads with the fast reader.

    It behaves like the DataFrame (attributes, items, len) for callers that need
    the table, which is loaded on first use.
    """

    def __init__(self, path: str) -> None:
        """Stores the path; nothing is read yet."""
        self.path = os.fspath(path)
        self._df = None

    def frame(self) -> pd.DataFrame:
        """Returns the atlas DataFrame, reading it on first use."""
        if self._df is None:
            self._df = read_atlas_pandas(self.path)
        return self._df

    def __getattr__(self, name):
        """Passes attribute access to the DataFrame."""
        if name in ("path", "_df"):
            raise AttributeError(name)
        return getattr(self.frame(), name)

    def __getitem__(self, key):
        """Passes item access to the DataFrame."""
        return self.frame()[key]

    def __len__(self):
        """Returns the number of atlas rows."""
        return len(self.frame())


def read_atlas(path: str) -> LazyAtlas:
    """Returns a LazyAtlas for `path`."""
    return LazyAtlas(path)


def unread_path(atlas: object) -> str | None:
    """Returns the file of an atlas that is not loaded, else None.

    That is a path, or a LazyAtlas whose table has not been read; a DataFrame
    gives None.
    """
    if isinstance(atlas, LazyAtlas) and atlas._df is None:
        return atlas.path
    return os.fspath(atlas) if isinstance(atlas, (str, os.PathLike)) else None


CPG, ICR = 0, 1  # kinds
_NAN = np.nan


class Unsupported(Exception):
    """The file is not in a form the fast reader handles."""


def tune_malloc() -> object | None:
    """Raises the glibc mmap and trim thresholds to keep buffers in the heap.

    The decompression buffers are several MB. With the default mmap threshold
    every buffer is an mmap/munmap pair, and page faults contend on the process
    mmap lock (12 threads were slower than 4). Best effort; other C libraries
    are left alone.

    Returns:
        The ctypes handle of libc (the caller may reset the thresholds), or None
        when they could not be set.
    """
    try:
        import ctypes
        libc = ctypes.CDLL("libc.so.6")
        libc.mallopt(-3, 1 << 30)  # M_MMAP_THRESHOLD
        libc.mallopt(-1, 1 << 30)  # M_TRIM_THRESHOLD
    except (OSError, AttributeError):
        return None
    return libc


@jit
def _numeric_token(buf, i, j):
    """Returns True when buf[i:j] is a number Python float() reads.

    That is [+-] digits [. digits] [e [+-] digits], or . digits ...
    """
    p = i
    if p < j and (buf[p] == 45 or buf[p] == 43):
        p += 1
    nd = 0
    while p < j and buf[p] >= 48 and buf[p] <= 57:
        p += 1
        nd += 1
    if p < j and buf[p] == 46:
        p += 1
        while p < j and buf[p] >= 48 and buf[p] <= 57:
            p += 1
            nd += 1
    if nd == 0:
        return False
    if p < j and (buf[p] == 101 or buf[p] == 69):
        p += 1
        if p < j and (buf[p] == 45 or buf[p] == 43):
            p += 1
        ne = 0
        while p < j and buf[p] >= 48 and buf[p] <= 57:
            p += 1
            ne += 1
        if ne == 0:
            return False
    return p == j


@jit
def _num(buf, i, j):
    """Returns (value, status) of buf[i:j].

    Status 0 is ok (NA tokens give NaN), 1 means Python float() must read it, 2
    is an error.
    """
    if i >= j:
        return _NAN, 0
    st, v = numparse.parse_decimal(buf, i, j)
    if st == 0:
        return v, 0
    L = j - i
    # NA, NaN or nan
    if ((L == 2 and buf[i] == 78 and buf[i + 1] == 65) or
        (L == 3 and buf[i + 1] == 97 and
         ((buf[i] == 78 and buf[i + 2] == 78) or
          (buf[i] == 110 and buf[i + 2] == 110)))):
        return _NAN, 0
    if _numeric_token(buf, i, j):
        return _NAN, 1
    return _NAN, 2


@jit
def _int(buf, i, j):
    st, v = numparse.parse_int(buf, i, j)
    return v, st == 0


@jit
def _chrom_code(buf, i, j, kind):
    """Returns the chromosome code of buf[i:j], or -1.

    CPG: chr1..chr22 are 0..21, chrX 22, chrY 23 (exact names). ICR: an optional
    case-insensitive `chr` prefix, then 1..22 without a leading zero, gives
    0..21 (v2.2 and pandas: normalize_primary_chrom, then isin chr1..chr22).
    """
    p = i
    if kind == 0:
        if j - i < 4 or buf[i] != 99 or buf[i + 1] != 104 or buf[i + 2] != 114:
            return -1
        p = i + 3
    else:
        if j - i >= 3 and (buf[i] | 32) == 99 and (buf[i + 1] | 32) == 104 and (
                buf[i + 2] | 32) == 114:
            p = i + 3
    L = j - p
    if L == 1:
        c = np.int64(buf[p])
        if c >= 49 and c <= 57:
            return c - 49
        if kind == 0:
            if c == 88:
                return 22
            if c == 89:
                return 23
        return -1
    if L == 2:
        a = np.int64(buf[p])
        b = np.int64(buf[p + 1])
        if a >= 49 and a <= 57 and b >= 48 and b <= 57:
            v = (a - 48) * 10 + (b - 48)
            if v <= 22:
                return v - 1
    return -1


@jit
def parse_lines(buf, kind, max_code, min_n, strict, m_lo, m_hi, m_hi_strict, cA,
                cPM, cOff, o_code, o_start, o_end, o_phi, o_meth, slow, counts):
    """Tokenises the complete lines of `buf` and keeps the wanted rows.

    A row is kept with a wanted chromosome (code <= max_code) and N > min_n
    (strict) or N >= min_n, and REF_METH in [m_lo, m_hi) (m_hi_strict) or [m_lo,
    m_hi] (NaN is dropped; a number that needs Python float() is kept and
    decided by the caller). START is parsed for kept rows only. For ICR the
    START/END must also closed-overlap a candidate interval (per chromosome
    code: sorted starts cA, prefix maximum of ends cPM, offsets cOff).

    `counts` receives lines, chromosome ok, kept after N, errors, slow numbers.

    Returns:
        The number of rows written.
    """
    n = buf.size
    pos = 0
    k = 0
    fs = np.empty(8, np.int64)
    fe = np.empty(8, np.int64)
    ncap = slow.shape[0]
    while pos < n:
        counts[0] += 1
        nf = 0
        p = pos
        while nf < 8:  # fields end at tab (9) or newline (10)
            fs[nf] = p
            while p < n and (buf[p] > 10 or buf[p] < 9):
                p += 1
            fe[nf] = p
            nf += 1
            if nf < 8:
                if p >= n or buf[p] == 10:
                    break
                p += 1
        q = p
        while q < n and buf[q] != 10:
            q += 1
        if nf < 8:
            counts[3] += 1
            pos = q + 1
            continue
        code = _chrom_code(buf, fs[0], fe[0], kind)
        if code < 0 or code > max_code:
            pos = q + 1
            continue
        counts[1] += 1
        nv, sn = _num(buf, fs[3], fe[3])
        if sn != 0:
            counts[3] += 1
            pos = q + 1
            continue
        if strict:
            keep = nv > min_n
        else:
            keep = nv >= min_n
        if not keep:
            pos = q + 1
            continue
        counts[2] += 1
        mv, sm = _num(buf, fs[7], fe[7])
        if sm == 2:
            counts[3] += 1
            pos = q + 1
            continue
        if sm == 0:
            if m_hi_strict:
                keep = mv >= m_lo and mv < m_hi
            else:
                keep = mv >= m_lo and mv <= m_hi
            if not keep:
                pos = q + 1
                continue
        s, ok = _int(buf, fs[1], fe[1])
        if not ok:
            counts[3] += 1
            pos = q + 1
            continue
        e = np.int64(0)
        if kind == 1:
            e, ok2 = _int(buf, fs[2], fe[2])
            if not ok2:
                counts[3] += 1
                pos = q + 1
                continue
            a0 = cOff[code]
            a1 = cOff[code + 1]
            # number of candidate starts <= e (binary search in cA[a0:a1])
            lo = a0
            hi = a1
            while lo < hi:
                mid = (lo + hi) >> 1
                if cA[mid] <= e:
                    lo = mid + 1
                else:
                    hi = mid
            if lo == a0 or cPM[lo - 1] < s:
                pos = q + 1
                continue
        if sm == 1:
            if counts[4] < ncap:
                slow[counts[4], 0] = k
                slow[counts[4], 1] = 1
                slow[counts[4], 2] = fs[7]
                slow[counts[4], 3] = fe[7]
            counts[4] += 1
            mv = _NAN
        pv = _NAN
        if kind == 0:
            pv, sp = _num(buf, fs[5], fe[5])
            if sp == 2:
                counts[3] += 1
                pos = q + 1
                continue
            if sp == 1:
                if counts[4] < ncap:
                    slow[counts[4], 0] = k
                    slow[counts[4], 1] = 0
                    slow[counts[4], 2] = fs[5]
                    slow[counts[4], 3] = fe[5]
                counts[4] += 1
                pv = _NAN
        o_code[k] = code
        o_start[k] = s
        o_end[k] = e
        o_phi[k] = pv
        o_meth[k] = mv
        k += 1
        pos = q + 1
    return k


def _parse_bytes(data: bytes, kind, max_code, min_n, strict, mfilt, win):
    """Parses the complete lines of `data`.

    Returns:
        (code, start, end, phi, meth, counts), or None without a newline.
    """
    nl = data.count(b"\n")
    if nl == 0:
        return None
    buf = np.frombuffer(data, np.uint8)
    cap = nl
    o_code = np.empty(cap, np.int8)
    o_start = np.empty(cap, np.int64)
    o_end = np.empty(cap, np.int64)
    o_phi = np.empty(cap, np.float64)
    o_meth = np.empty(cap, np.float64)
    scap = 64
    while True:
        slow = np.zeros((scap, 4), np.int64)
        counts = np.zeros(5, np.int64)
        k = parse_lines(buf, kind, max_code, float(min_n), strict, mfilt[0],
                        mfilt[1], mfilt[2], win[0], win[1], win[2], o_code,
                        o_start, o_end, o_phi, o_meth, slow, counts)
        if counts[4] <= scap:
            break
        scap = int(counts[4])
    drop = None
    for r, col, a, b in slow[:counts[4]].tolist():
        v = float(data[a:b])
        (o_meth if col == 1 else o_phi)[r] = v
        # decided here: the kernel kept the row without testing REF_METH
        if col == 1:
            ok = v >= mfilt[0] and (v < mfilt[1] if mfilt[2] else v <= mfilt[1])
            if not ok:
                if drop is None:
                    drop = np.ones(k, np.bool_)
                drop[r] = False
    if drop is not None:
        return (o_code[:k][drop], o_start[:k][drop], o_end[:k][drop],
                o_phi[:k][drop], o_meth[:k][drop], counts)
    if k < cap:
        return (o_code[:k].copy(), o_start[:k].copy(), o_end[:k].copy(),
                o_phi[:k].copy(), o_meth[:k].copy(), counts)
    return o_code, o_start, o_end, o_phi, o_meth, counts


def read_columns(path: str,
                 kind: int,
                 *,
                 max_code: int,
                 min_n: float,
                 strict: bool,
                 meth_lo: float = -np.inf,
                 meth_hi: float = np.inf,
                 meth_hi_strict: bool = False,
                 windows: tuple | None = None,
                 threads: int | None = None,
                 batch_bytes: int = 4 << 20) -> tuple:
    """Reads the parsed, filtered rows of the atlas at `path`.

    Args:
        path: Plain, gzip or BGZF atlas.
        kind: CPG or ICR.
        max_code: Largest chromosome code kept.
        min_n: N filter.
        strict: N > min_n instead of N >= min_n.
        meth_lo: Smallest REF_METH kept.
        meth_hi: Largest REF_METH kept.
        meth_hi_strict: REF_METH < meth_hi instead of <=.
        windows: (cA, cPM, cOff) from `candidate_windows` (ICR).
        threads: Reader threads.
        batch_bytes: Bytes per parse task.

    Returns:
        (code int8, start, end, phi, meth, counts) in file order; counts is
        [lines, chromosome ok, kept after N, errors, slow numbers].

    Raises:
        Unsupported: The file is empty or has a broken block chain.
    """
    if windows is None:
        z = np.zeros(1, np.int64)
        windows = (z, z, np.zeros(25, np.int64))
    threads = threads or bed.default_threads(8)
    import multiprocessing as mp
    import threading
    # fork only from a single-threaded, non-daemonic process (a daemonic worker
    # may not have children; forking next to other threads can inherit held
    # locks); otherwise use threads
    use_threads = (os.environ.get("DMRDETECTOR_EXECUTOR") == "thread" or
                   mp.current_process().daemon or threading.active_count() > 1)
    if use_threads:
        tune_malloc()
    args = (kind, max_code, min_n, strict, (float(meth_lo), float(meth_hi),
                                            bool(meth_hi_strict)), windows)
    with open(path, "rb") as fh:
        if fh.read(2) == b"\x1f\x8b":
            gz = True
        else:
            gz = False
        fh.seek(0)
        size = os.fstat(fh.fileno()).st_size
        if size == 0:
            raise Unsupported("empty file")
        mm = mmap.mmap(fh.fileno(), 0, access=mmap.ACCESS_READ)
        try:
            std = gz and os.pread(fh.fileno(), 16, 0) == bgzf.MAGIC
            blocks = None
            if gz and not std:
                t = bgzf.scan_blocks(np.frombuffer(mm, np.uint8))
                if t is not None:
                    blocks = list(
                        zip(t[3].tolist(), (t[0] + t[1] - 8 - t[3]).tolist(),
                            t[2].tolist()))
            if std or blocks is not None:
                # standard header: workers find the blocks of their own byte
                # range
                if std:
                    groups = [(a, min(a + _CHUNK, size))
                              for a in range(0, size, _CHUNK)]
                    work = _bgzf_chunk_work
                else:
                    groups, cur, acc = [], [], 0
                    for off, ln, isz in blocks:
                        cur.append((off, ln, isz))
                        acc += 65536
                        if acc >= batch_bytes:
                            groups.append(cur)
                            cur, acc = [], 0
                    if cur:
                        groups.append(cur)
                    work = _bgzf_work

                global _MM, _ARGS, _FD, _SIZE
                _MM, _ARGS, _FD, _SIZE = mm, args, fh.fileno(), size
                if use_threads:
                    with ThreadPoolExecutor(threads) as ex:
                        results = list(ex.map(work, groups))
                else:  # fork: the mmap and the compiled kernel are inherited
                    from dmrdetector import runtime  # RunError, pool helper
                    with runtime.fork_pool(threads, "atlas reader") as ex:
                        results = list(ex.map(work, groups))
            else:

                def chunks():
                    if gz:
                        d = zlib.decompressobj(31)
                        pos = 0
                        while pos < size:
                            b = mm[pos:pos + (8 << 20)]
                            pos += len(b)
                            while b:
                                yield d.decompress(b)
                                if d.eof:
                                    b = d.unused_data
                                    d = zlib.decompressobj(31)
                                else:
                                    b = b""
                    else:
                        pos = 0
                        while pos < size:
                            yield mm[pos:pos + batch_bytes]
                            pos += batch_bytes

                with ThreadPoolExecutor(threads) as ex:
                    futs, pend = [], b""
                    for ch in chunks():
                        pend += ch
                        if len(pend) >= batch_bytes:
                            futs.append(ex.submit(_split_parse, pend, args))
                            pend = b""
                    if pend:
                        futs.append(ex.submit(_split_parse, pend, args))
                    results = [f.result() for f in futs]
        finally:
            mm.close()
    if std:
        pos = 0  # every byte of the file must belong to exactly one block
        for r, first, last_end in results:
            if first is None:
                continue
            if first != pos:
                raise Unsupported("BGZF block chain broken")
            pos = last_end
        if pos != size:
            raise Unsupported("BGZF block chain broken")
        results = [r for r, _f, _l in results]
    return _stitch(results, args)


_MM = None
_FD = -1
_SIZE = 0
_ARGS = None


def _bgzf_work(g):
    a, b = g[0][0], g[-1][0] + g[-1][1]
    # one read per batch (page faults of a mmap'ed NAS file were slow)
    raw = os.pread(_FD, b - a, a)
    return _split_parse(bgzf.inflate_group(raw, a, g), _ARGS)


_CHUNK = 3 << 19  # compressed bytes per worker task (about 64 blocks)


def _bgzf_chunk_work(ab):
    a, b = ab
    raw = os.pread(_FD, min(b + 65536, _SIZE) - a, a)
    g = bgzf.blocks_starting_in(raw, a, b, _SIZE)
    if not g:
        return ("nonl", b""), None, None
    first = g[0][0] - 18
    last_end = g[-1][0] + g[-1][1] + 8
    return _split_parse(bgzf.inflate_group(raw, a, g), _ARGS), first, last_end


def _split_parse(data: bytes, args):
    """Cuts `data` into a head line, whole lines and a tail.

    The head continues the previous batch; the whole lines are parsed.
    """
    first = data.find(b"\n")
    if first < 0:
        return ("nonl", data)
    last = data.rfind(b"\n")
    body = _parse_bytes(data[first + 1:last +
                             1], *args) if last > first else None
    return ("ok", data[:first + 1], body, data[last + 1:])


def _stitch(results, args):
    pieces = []
    carry = b""
    for r in results:
        if r[0] == "nonl":
            carry += r[1]
            continue
        _, head, body, tail = r
        pieces.append(_parse_bytes(carry + head, *args))
        pieces.append(body)
        carry = tail
    if carry:
        pieces.append(_parse_bytes(carry + b"\n", *args))
    pieces = [p for p in pieces if p is not None]
    counts = np.zeros(5, np.int64)
    for p in pieces:
        counts += p[5]
    if counts[3]:
        raise Unsupported("unparsable lines")
    if not pieces:
        e = np.empty(0)
        return (np.empty(0, np.int8), np.empty(0, np.int64),
                np.empty(0, np.int64), e, e, counts)
    cols = [np.concatenate([p[i] for p in pieces]) for i in range(5)]
    return (*cols, counts)


def candidate_windows(
        codes: np.ndarray, a: np.ndarray,
        b: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Returns the candidate windows per chromosome code (0..23).

    These are the sorted candidate starts, the prefix maximum of the ends and
    the offsets. A site [s, e] can overlap a candidate (closed) only if some A_i
    <= e and B_i >= s.
    """
    order = np.lexsort((a, codes))
    codes, a, b = codes[order], a[order], b[order]
    off = np.searchsorted(codes, np.arange(25), side="left").astype(np.int64)
    pm = np.empty_like(b)
    for c in range(24):
        lo, hi = off[c], off[c + 1]
        if hi > lo:
            pm[lo:hi] = np.maximum.accumulate(b[lo:hi])
    if a.size == 0:
        z = np.zeros(1, np.int64)
        return z, z, off
    return a.astype(np.int64), pm.astype(np.int64), off
