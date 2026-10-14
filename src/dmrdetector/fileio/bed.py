"""Text helpers: chromosome names, R fwrite number text, gzip, tabix.

Number text follows R `data.table::fwrite` (15 significant digits, fixed or
scientific by width). `fmt_doubles` rounds exactly; fwrite can differ in the
15th digit for about 3 % of random ratios (relative size 1e-15).
"""

from __future__ import annotations

from collections.abc import Callable
import gzip
import os
import re
import sys
from typing import Iterable, Sequence, TextIO

import numpy as np
import pandas as pd


PRIMARY = {str(i) for i in range(1, 23)} | {"X", "Y"}
AUTOSOMES = [f"chr{i}" for i in range(1, 23)]
RAW_ATLAS_COLS = ["CHROM", "START", "END", "N", "COV", "PHI", "RAW_METH"]


def c_order(chroms: Iterable[str]) -> list[str]:
    """Returns the unique names in C-locale (byte) order.

    This is the order of data.table setkey() and `sort -k1,1` under LC_ALL=C:
    chr1, chr10, ..., chr19, chr2, ..., chrX.
    """
    return sorted(set(chroms))


def normalize_primary_chrom(
        x: Iterable[str] | pd.Series | np.ndarray) -> pd.Series:
    """Maps chr1-22, X and Y (any case, with or without chr) to `chrN`.

    Other names are returned unchanged, as the v2.2 function does.

    Returns:
        A Series with the index of `x` (a RangeIndex for lists and arrays).
    """
    s = pd.Series(x, dtype=object).astype(str)
    base = s.str.replace(r"^[cC][hH][rR]", "", regex=True).str.upper()
    primary = base.isin([str(i) for i in range(1, 23)] + ["X", "Y"])
    out = s.copy()
    out[primary] = "chr" + base[primary]
    return out


def chrom_codes(chrom: np.ndarray) -> tuple[np.ndarray, list[str]]:
    """Returns integer codes in C-locale order, plus the ordered names."""
    names = c_order(chrom)
    cat = pd.Categorical(chrom, categories=names, ordered=True)
    return cat.codes.astype(np.int32), names


def _fwrite_one(v: float) -> str:
    """Returns one double as fwrite writes it.

    That is 15 significant digits, the shorter of fixed and scientific, fixed on
    ties.
    """
    if v != v:
        return "NaN"
    if v in (float("inf"), float("-inf")):
        return "Inf" if v > 0 else "-Inf"
    # digits and decimal exponent from scientific form
    m = re.match(r"^(-?)(\d)(?:\.(\d+))?e([+-]\d+)$", "%.14e" % v)
    sign, d0 = m.group(1), m.group(2)
    rest, ex = m.group(3) or "", int(m.group(4))
    digits = (d0 + rest).rstrip("0") or "0"
    if digits == "0":
        return "0"
    nd = len(digits)
    # scientific string
    sci = (sign + digits[0] + ("." + digits[1:] if nd > 1 else "") + "e" +
           ("-" if ex < 0 else "+") + ("%02d" % abs(ex)))
    # fixed string
    if ex >= 0:
        if nd <= ex + 1:
            fixed = sign + digits + "0" * (ex + 1 - nd)
        else:
            fixed = sign + digits[:ex + 1] + "." + digits[ex + 1:]
    else:
        fixed = sign + "0." + "0" * (-ex - 1) + digits
    return fixed if len(fixed) <= len(sci) else sci


def fmt_doubles(x: np.ndarray) -> list[str]:
    """Formats a float array like fwrite.

    Fast path through %.15g; only values whose fixed/scientific choice can
    differ from %g (a single significant digit with exponent -4 or 5..14, or
    magnitude >= 1e15) take the slow path.
    """
    x = np.asarray(x, dtype=np.float64)
    out = ["%.15g" % v for v in x.tolist()]
    ax = np.abs(x)
    cand = np.flatnonzero(((ax >= 1e-4) & (ax < 1e-3)) | ((ax >= 1e5) &
                                                          (ax < 1e15)) |
                          (ax >= 1e15) | ~np.isfinite(x))
    for i in cand:
        out[i] = _fwrite_one(float(x[i]))
    return out


def fmt_int(x: np.ndarray) -> list[str]:
    """Formats an integer array as text."""
    return np.asarray(x).astype(np.int64).astype(str).tolist()


def awk_number(x: np.ndarray) -> list[str]:
    """Formats numbers as awk `print x` does.

    Integral values are written as integers, others with OFMT %.6g (v2.2 prints
    coverage and N_mod this way after `$10+0`).
    """
    x = np.asarray(x, dtype=np.float64)
    integral = x == np.floor(x)
    out = np.where(integral, x.astype(np.int64).astype(str), "")
    nz = np.flatnonzero(~integral)
    if nz.size:
        out = out.astype(object)
        for i in nz:
            out[i] = "%.6g" % x[i]
    return out.tolist()


def join_columns(cols: Sequence[Sequence[str]]) -> str:
    """Joins parallel string columns into tab-separated lines."""
    if not len(cols[0]):
        return ""
    return "\n".join("\t".join(r) for r in zip(*cols)) + "\n"


def fmt_double(x: float) -> str:
    """Formats a double like data.table::fwrite; NaN becomes an empty cell."""
    x = float(x)
    return "" if x != x else _fwrite_one(x)


def fwrite_table(
    df: pd.DataFrame,
    path: str,
    int_cols: Iterable[str] = (),
    str_cols: Iterable[str] = (),
    float_cols: Iterable[str] = ()
) -> None:
    """Writes `df` tab-separated and unquoted, like fwrite.

    The file is gzip when `path` ends in .gz. Columns in `int_cols` and
    `float_cols` are written as integers and as fwrite doubles; the others as
    text.
    """
    cols = list(df.columns)
    parts = []
    for c in cols:
        s = df[c]
        if c in float_cols:
            parts.append([fmt_double(v) for v in s.to_numpy(dtype=float)])
        elif c in int_cols:
            parts.append([("" if pd.isna(v) else str(int(v))) for v in s])
        else:
            parts.append([("" if
                           (v is None or
                            (isinstance(v, float) and v != v)) else str(v))
                          for v in s])
    lines = ["\t".join(cols)] + ["\t".join(r) for r in zip(*parts)
                                ] if parts else ["\t".join(cols)]
    data = ("\n".join(lines) + "\n").encode()
    if path.endswith(".gz"):
        with gzip.open(path, "wb", compresslevel=6) as fh:
            fh.write(data)
    else:
        with open(path, "wb") as fh:
            fh.write(data)


def write_atlas_text(chrom: np.ndarray, start: np.ndarray, end: np.ndarray,
                     n: np.ndarray, cov: np.ndarray, *float_cols:
                     np.ndarray) -> str:
    """Returns an atlas block as text: integers, then fwrite doubles."""
    cols = [
        chrom.tolist() if isinstance(chrom, np.ndarray) else list(chrom),
        fmt_int(start),
        fmt_int(end),
        fmt_int(n),
        fmt_int(cov)
    ]
    cols += [fmt_doubles(c) for c in float_cols]
    return join_columns(cols)


def is_gzip(path: str) -> bool:
    """Returns True when the file starts with the gzip magic bytes."""
    with open(path, "rb") as fh:
        return fh.read(2) == b"\x1f\x8b"


def open_text(path: str) -> TextIO:
    """Opens a plain or gzip text file, as `zcat-or-cat` does in v2.2."""
    return gzip.open(path, "rt") if is_gzip(path) else open(path, "rt")


def write_bgzip_text(path: str, text_chunks: Iterable[str]) -> None:
    """Writes text chunks to a BGZF file with pysam."""
    import pysam
    w = pysam.BGZFile(path, "wb")
    try:
        for t in text_chunks:
            if t:
                w.write(t.encode())
    finally:
        w.close()


def tabix_bed(path: str) -> None:
    """Creates `path`.tbi (preset bed) for a bgzip file, overwriting it."""
    import pysam
    pysam.tabix_index(path, preset="bed", force=True, keep_original=True)


def write_bed_bgzip_tabix(path: str, text: str) -> None:
    """Writes sorted BED text as bgzip plus tabix index (`tabix -p bed -f`).

    An empty text gives a valid empty bgzip file and no index.
    """
    write_bgzip_text(path, [text])
    if text:
        tabix_bed(path)
    elif os.path.exists(path + ".tbi"):
        os.remove(path + ".tbi")


def ensure_tabix(path: str, scratch: str, check_age: bool = False) -> str:
    """Returns a path with a usable .tbi.

    That is `path` itself if indexed, else a symlink that is indexed in a folder
    of `scratch` named by a hash of the absolute path, so that files with one
    name in different folders do not share a link or an index (inputs stay
    untouched). With `check_age` an index older than its data counts as
    missing.
    """
    import hashlib

    def fresh(p):
        return os.path.exists(p + ".tbi") and (not check_age or
                                               os.path.getmtime(p + ".tbi")
                                               >= os.path.getmtime(p))

    if fresh(path):
        return path
    sub = os.path.join(
        scratch,
        hashlib.sha1(os.path.abspath(path).encode()).hexdigest()[:16])
    os.makedirs(sub, exist_ok=True)
    link = os.path.join(sub, os.path.basename(path))
    if not os.path.lexists(link):
        os.symlink(os.path.abspath(path), link)
    if not fresh(link):
        tabix_bed(link)
    return link


def default_threads(default: int = 1) -> int:
    """Returns the threads for the bulk readers and writers.

    DMRDETECTOR_THREADS if set, else `default`. `build` and `run` set the
    variable in every worker to that task's thread share (runtime.pool_map), so
    a run never exceeds --threads.
    """
    v = os.environ.get("DMRDETECTOR_THREADS")
    return max(1, int(v)) if v else max(1, default)


def map_text(path: str,
             fn: Callable[[np.ndarray, bool], object],
             threads: int | None = None,
             piece: int = 4 << 20) -> list:
    """Applies fn(uint8 array, first) to the text of a plain, gzip or BGZF file.

    fn sees whole lines only: the text is cut at newlines and the lines that
    straddle two pieces are passed to fn on their own, so fn must work on any
    set of lines. It is called from several threads and must release the GIL to
    scale. `first` is True for the piece at the start of the file.

    Returns:
        The results of fn in file order.
    """
    import mmap
    threads = threads or default_threads()
    if os.path.getsize(path) == 0:
        return []
    with open(path, "rb") as fh:
        mm = mmap.mmap(fh.fileno(), 0, access=mmap.ACCESS_READ)
        try:
            blocks = None
            if mm[:2] == b"\x1f\x8b":
                from dmrdetector.fileio import bgzf
                blocks = bgzf.scan_blocks(np.frombuffer(mm, np.uint8))
            if blocks is not None:
                return _map_bgzf(mm, blocks, fn, threads, piece)
            return _map_stream(mm, fn, threads, piece)
        finally:
            try:
                mm.close()
            # a numpy view of the file is still alive (exception in flight)
            except BufferError:
                pass


def _call(fn, b: bytes, first: bool = False):
    """Calls fn on the bytes `b`, or returns None for no bytes."""
    return fn(np.frombuffer(b, dtype=np.uint8), first) if b else None


def _map_bgzf(mm, blocks, fn, threads, piece):
    from concurrent.futures import ThreadPoolExecutor

    from dmrdetector.fileio import bgzf
    coff, csize, usize, doff = blocks
    triples = list(
        zip(doff.tolist(), (coff + csize - 8 - doff).tolist(), usize.tolist()))
    groups, cur, size = [], [], 0
    for b in triples:
        cur.append(b)
        size += b[2]
        if size >= piece:
            groups.append(cur)
            cur, size = [], 0
    if cur:
        groups.append(cur)

    def work(i):
        data = bgzf.inflate_group(mm, 0, groups[i])
        if i == 0:
            head, rest = b"", data
        else:
            k = data.find(b"\n")
            head, rest = (data, b"") if k < 0 else (data[:k + 1], data[k + 1:])
        k = rest.rfind(b"\n")
        return head, _call(fn, rest[:k + 1], i == 0), rest[k + 1:]

    with ThreadPoolExecutor(threads) as ex:
        parts = list(ex.map(work, range(len(groups))))
    out, carry = [], b""
    for i, (head, res, tail) in enumerate(parts):
        if i > 0:
            carry += head
            # no newline in this piece: keep carrying
            if not carry.endswith(b"\n"):
                continue
            out.append(_call(fn, carry))
            carry = b""
        out.append(res)
        carry += tail
    if carry:
        out.append(_call(fn, carry))
    return [r for r in out if r is not None]


def _map_stream(mm, fn, threads, piece):
    from concurrent.futures import ThreadPoolExecutor
    import zlib
    gz = mm[:2] == b"\x1f\x8b"

    def pieces():
        d = zlib.decompressobj(31) if gz else None
        pos, carry, n = 0, b"", len(mm)
        while pos < n:
            raw = mm[pos:pos + (1 << 20 if gz else piece)]
            pos += len(raw)
            if gz:
                data = b""
                while raw:
                    data += d.decompress(raw)
                    if d.eof:
                        raw = d.unused_data
                        d = zlib.decompressobj(31)
                    else:
                        raw = b""
            else:
                data = raw
            carry += data
            if len(carry) >= piece:
                k = carry.rfind(b"\n")
                if k >= 0:
                    yield carry[:k + 1]
                    carry = carry[k + 1:]
        if carry:
            yield carry

    out = []
    with ThreadPoolExecutor(threads) as ex:
        futs = []
        for b in pieces():
            futs.append(ex.submit(_call, fn, b, not out and not futs))
            if len(futs) > 2 * threads:
                out.append(futs.pop(0).result())
        out += [f.result() for f in futs]
    return [r for r in out if r is not None]


def overlap_counts(start: np.ndarray, end: np.ndarray, b_start: np.ndarray,
                   b_end: np.ndarray) -> np.ndarray:
    """Counts the B intervals that overlap each A interval.

    Half-open BED semantics (a.start < b.end and b.start < a.end), as bedtools
    intersect. B intervals may overlap each other.
    """
    cnt = np.zeros(start.shape[0], dtype=np.int64)
    for bs, be in zip(np.asarray(b_start).tolist(), np.asarray(b_end).tolist()):
        cnt += (start < be) & (bs < end) & (end > start)
    return cnt


def sorted_point_overlap(start: np.ndarray, end: np.ndarray,
                         p0: np.ndarray) -> np.ndarray:
    """Returns True where [start, end) holds a point of `p0` (sorted)."""
    if p0.size == 0 or start.size == 0:
        return np.zeros(start.shape[0], dtype=bool)
    # exists p0 with start-1 < p0 < end  i.e.  start <= p0 <= end-1
    lo = np.searchsorted(p0, start, side="left")
    hi = np.searchsorted(p0, end, side="left")
    return (hi > lo) & (end > start)


def bed_sort_order(chrom: np.ndarray, start: np.ndarray) -> np.ndarray:
    """Row order of `sort -k1,1 -k2,2n` (C locale; stable among ties)."""
    codes, _ = chrom_codes(chrom)
    return np.lexsort((start, codes))


def warn(msg: str) -> None:
    """Prints a message to stderr."""
    print(msg, file=sys.stderr)
