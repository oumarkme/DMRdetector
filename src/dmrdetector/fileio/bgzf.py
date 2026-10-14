"""BGZF blocks and tabix indexes: block scan, inflate, compress, write, index.

Inflate uses libdeflate through ctypes (zlib when it is missing or
DMRDETECTOR_NO_LIBDEFLATE=1). Writers deflate in a thread pool and report the
virtual offset of every row, which `write_tbi` turns into a tabix index without
a second pass. Design notes: doc/ALGORITHM.md, section Fast readers.
"""

from __future__ import annotations

from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor
import ctypes
import gzip
import os
import struct
import sys
import threading
import zlib

import numpy as np

from dmrdetector.util.jit import jit


BLOCK = 0xFF00  # uncompressed bytes per BGZF block, as written by htslib
EOF_BLOCK = bytes.fromhex(
    "1f8b08040000000000ff0600424302001b0003000000000000000000")
# header of an htslib/bgzip block up to BSIZE
MAGIC = bytes.fromhex("1f8b08040000000000ff060042430200")
_HDR = MAGIC[:16]
NO_JIT = os.environ.get("DMRDETECTOR_NO_JIT") == "1"
MAX_TBI_POS = 1 << 29  # largest END of a tabix (not csi) index


class ParseError(Exception):
    """The fast reader cannot handle this file; use the pandas path."""


_LD = None
_TLS = threading.local()


def libdeflate() -> ctypes.CDLL | bool:
    """Returns the ctypes handle of libdeflate, or False.

    libdeflate ships with htslib and is several times faster than zlib.
    DMRDETECTOR_NO_LIBDEFLATE=1 forces zlib.
    """
    global _LD
    if _LD is None:
        _LD = False
        if os.environ.get("DMRDETECTOR_NO_LIBDEFLATE") != "1":
            import ctypes.util
            for c in (os.path.join(sys.prefix, "lib", "libdeflate.so.0"),
                      ctypes.util.find_library("deflate")):
                if not c:
                    continue
                try:
                    lib = ctypes.CDLL(c)
                    vp, sz = ctypes.c_void_p, ctypes.c_size_t
                    lib.libdeflate_alloc_decompressor.restype = vp
                    lib.libdeflate_free_decompressor.argtypes = [vp]
                    lib.libdeflate_deflate_decompress.argtypes = [
                        vp, vp, sz, vp, sz, vp
                    ]
                    lib.libdeflate_deflate_decompress.restype = ctypes.c_int
                    lib.libdeflate_alloc_compressor.restype = vp
                    lib.libdeflate_alloc_compressor.argtypes = [ctypes.c_int]
                    lib.libdeflate_free_compressor.argtypes = [vp]
                    lib.libdeflate_deflate_compress.argtypes = [
                        vp, vp, sz, vp, sz
                    ]
                    lib.libdeflate_deflate_compress.restype = sz
                    _LD = lib
                    break
                except (OSError, AttributeError):
                    continue
    return _LD


# level of the BGZF files the stages write: preprocessed beds (read a few times)
# take the fast level, atlases (kept as the reference) level 6. libdeflate level
# 1 makes 35 MB per 150 MB of text in 0.8 s, level 6 31 MB in 1.9 s.
DEFAULT_LEVEL = int(
    os.environ.get("DMRDETECTOR_BGZF_LEVEL", "1" if libdeflate() else "4"))
ARCHIVE_LEVEL = int(os.environ.get("DMRDETECTOR_BGZF_LEVEL", "6"))


@jit
def _scan(raw, coff, csize, usize, doff):
    """Walk the block chain of the file bytes `raw` from offset 0.

    Fills the offset and size of every block, the uncompressed size from its
    trailer and the offset of its deflate data; returns the number of blocks, or
    -1 when the chain is not a BGZF file.
    """
    n = 0
    p = 0
    size = raw.shape[0]
    while p < size:
        if p + 18 > size:
            return -1
        if raw[p] != 31 or raw[p + 1] != 139 or raw[p + 2] != 8 or (raw[p + 3] &
                                                                    4) == 0:
            return -1
        xlen = np.int64(raw[p + 10]) | (np.int64(raw[p + 11]) << 8)
        q = p + 12
        end = q + xlen
        bsize = -1
        while q + 4 <= end:
            slen = np.int64(raw[q + 2]) | (np.int64(raw[q + 3]) << 8)
            if raw[q] == 66 and raw[q + 1] == 67 and slen == 2:
                bsize = (np.int64(raw[q + 4]) | (np.int64(raw[q + 5]) << 8)) + 1
            q += 4 + slen
        if bsize < 26 + xlen - 6 or p + bsize > size:
            return -1
        coff[n] = p
        csize[n] = bsize
        doff[n] = end
        t = p + bsize - 4
        usize[n] = (np.int64(raw[t]) | (np.int64(raw[t + 1]) << 8) |
                    (np.int64(raw[t + 2]) << 16) | (np.int64(raw[t + 3]) << 24))
        n += 1
        p += bsize
    return n


def scan_blocks(
    raw: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray] | None:
    """Scans the blocks of a BGZF file held in `raw`.

    Args:
        raw: uint8 array (bytes through np.frombuffer, or an mmap).

    Returns:
        (coff, csize, usize, doff) of every block, or None when it is not BGZF.
        The deflate data of block k is raw[doff[k]:coff[k] + csize[k] - 8].
    """
    if raw.size < 28 or raw[0] != 31 or raw[1] != 139:
        return None
    cap = raw.size // 28 + 1
    coff = np.empty(cap, np.int64)
    csize = np.empty(cap, np.int64)
    usize = np.empty(cap, np.int64)
    doff = np.empty(cap, np.int64)
    n = _scan(raw, coff, csize, usize, doff)
    if n < 0:
        return None
    return coff[:n].copy(), csize[:n].copy(), usize[:n].copy(), doff[:n].copy()


def blocks_starting_in(raw: bytes, a: int, b: int,
                       size: int) -> list[tuple[int, int, int]]:
    """Finds the blocks of a standard-header BGZF file that start in [a, b).

    A block start is found by the header pattern (see MAGIC) and accepted when
    the next block, or the end of the file, follows where BSIZE says. Workers
    can split a file by byte ranges without a serial scan; the caller checks
    that the ranges tile the file.

    Args:
        raw: The file bytes from offset `a`, long enough to hold the last block.
        a: First byte of the range.
        b: End of the range.
        size: File size.

    Returns:
        [(offset of deflate data, its length, uncompressed size)].
    """
    out = []
    s = 0 if a == 0 else raw.find(MAGIC)
    while 0 <= s < b - a:
        if len(raw) < s + 18:
            break
        end = s + (raw[s + 16] | (raw[s + 17] << 8)) + 1
        nxt_ok = a + end == size or raw[end:end + 16] == MAGIC
        if end - s >= 26 and nxt_ok:
            if end > len(raw):
                break
            out.append((a + s + 18, end - 8 - s - 18,
                        int.from_bytes(raw[end - 4:end], "little")))
            s = end
        else:
            s = raw.find(MAGIC, s + 1)
    return out


def inflate_blocks(src: np.ndarray, off: np.ndarray, ln: np.ndarray,
                   isz: np.ndarray, dst: np.ndarray,
                   dst_off: np.ndarray) -> None:
    """Inflates raw-deflate payloads.

    Block k is src[off[k]:off[k] + ln[k]] and fills dst[dst_off[k]:+isz[k]],
    with libdeflate through ctypes, else zlib.

    Raises:
        ValueError: A block is corrupt or its size differs from the trailer.
    """
    off, ln, isz, dst_off = (
        np.asarray(x).tolist() for x in (off, ln, isz, dst_off))
    lib = libdeflate()
    if lib:
        sp, dp = src.ctypes.data, dst.ctypes.data
        d = lib.libdeflate_alloc_decompressor()
        try:
            for o, n, m, p in zip(off, ln, isz, dst_off):
                if m and lib.libdeflate_deflate_decompress(
                        d, sp + o, n, dp + p, m, None) != 0:
                    raise ValueError("bad deflate block")
        finally:
            lib.libdeflate_free_decompressor(d)
        return
    mv = memoryview(src)
    for o, n, m, p in zip(off, ln, isz, dst_off):
        if m:
            try:
                data = zlib.decompress(mv[o:o + n], -15)
            except zlib.error as e:
                raise ValueError(f"bad deflate block ({e})") from None
            if len(data) != m:
                raise ValueError("bad deflate block")
            dst[p:p + m] = np.frombuffer(data, np.uint8)


def inflate_group(raw: bytes | np.ndarray, base: int,
                  group: list[tuple[int, int, int]]) -> bytearray:
    """Returns the concatenated payloads of a group of BGZF blocks.

    Args:
        raw: File bytes from offset `base` (bytes, mmap or uint8 array).
        base: File offset of raw[0].
        group: [(offset of deflate data, length, uncompressed size)] in file
            offsets.
    """
    isz = [g[2] for g in group]
    out = bytearray(sum(isz))
    if out:
        pos = np.concatenate(([0], np.cumsum(isz)[:-1]))
        inflate_blocks(np.frombuffer(raw, np.uint8),
                       [g[0] - base for g in group], [g[1] for g in group], isz,
                       np.frombuffer(out, np.uint8), pos)
    return out


def read_blocks(
    path: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray] | None:
    """Returns (raw bytes as uint8 array, coff, csize, usize, doff), or None.

    None means the file is not BGZF.
    """
    raw = np.fromfile(path, dtype=np.uint8)
    t = scan_blocks(raw)
    return None if t is None else (raw, *t)


def count_newlines(path: str, threads: int = 1) -> int | None:
    """Counts the newline bytes of a BGZF file.

    The file is inflated in `threads` threads and nothing is kept. Returns None
    when the file is not BGZF.
    """
    rb = read_blocks(path)
    if rb is None:
        return None
    raw, coff, csize, usize, doff = rb
    triples = list(
        zip(doff.tolist(), (coff + csize - 8 - doff).tolist(), usize.tolist()))
    groups = [triples[i:i + 64] for i in range(0, len(triples), 64)]

    def count(g):
        """Returns the newlines of one group of blocks (GIL released)."""
        d = np.frombuffer(inflate_group(raw, 0, g), np.uint8)
        return count_nl(d, 0, d.size)

    with ThreadPoolExecutor(max(1, threads)) as ex:
        return int(sum(ex.map(count, groups)))


def inflate(path: str, threads: int = 1) -> np.ndarray | None:
    """Returns the decompressed text of a BGZF file as a uint8 array.

    A trailing newline is guaranteed. Returns None when the file is not BGZF.
    """
    rb = read_blocks(path)
    if rb is None:
        return None
    raw, coff, csize, usize, doff = rb
    uoff = np.concatenate(([0], np.cumsum(usize)))[:-1]
    total = int(usize.sum())
    out = np.empty(total + 1, np.uint8)
    nb = len(coff)
    dlen = coff + csize - 8 - doff
    t = max(1, min(threads, nb // 4 or 1))
    cut = np.linspace(0, nb, t + 1).astype(int)

    def run(i):
        """Inflates the blocks of part `i` into `out`."""
        s = slice(cut[i], cut[i + 1])
        inflate_blocks(raw, doff[s], dlen[s], usize[s], out, uoff[s])

    try:
        if t == 1:
            run(0)
        else:
            with ThreadPoolExecutor(t) as ex:
                for f in [ex.submit(run, i) for i in range(t)]:
                    f.result()
    except ValueError as e:
        raise ParseError(f"corrupt BGZF block in file ({e})") from None
    if total == 0:
        return out[:0]
    if out[total - 1] != 10:
        out[total] = 10
        return out[:total + 1]
    return out[:total]


def whole_blocks_len(raw: bytes) -> int:
    """Returns the length of the longest prefix of whole blocks in `raw`.

    `raw` starts at a block boundary and holds htslib-style BGZF blocks.

    Raises:
        ParseError: The first block does not have the standard header.
    """
    n = len(raw)
    p = 0
    while p + 18 <= n:
        if raw[p:p + 16] != _HDR:
            raise ParseError("BGZF block with a non-standard header")
        b = (raw[p + 16] | (raw[p + 17] << 8)) + 1
        if p + b > n:
            break
        p += b
    return p


def inflate_chunk(raw: bytes) -> np.ndarray:
    """Returns the text of the whole BGZF blocks in `raw` as a uint8 array.

    One thread, nothing appended.

    Raises:
        ParseError: The blocks are not a valid chain.
    """
    arr = np.frombuffer(raw, np.uint8)
    t = scan_blocks(arr) if arr.size >= 28 else None
    if t is None:
        raise ParseError("not a BGZF block chain")
    coff, csize, usize, doff = t
    uoff = np.concatenate(([0], np.cumsum(usize)))[:-1]
    out = np.empty(int(usize.sum()), np.uint8)
    try:
        inflate_blocks(arr, doff, coff + csize - 8 - doff, usize, out, uoff)
    except ValueError as e:
        raise ParseError(f"corrupt BGZF block in file ({e})") from None
    return out


@jit
def count_nl(buf, lo, hi):
    """Counts the newline bytes of buf[lo:hi]."""
    c = 0
    for i in range(lo, hi):
        if buf[i] == 10:
            c += 1
    return c


class _Compressors(dict):
    """Per-thread libdeflate compressors by level, freed with the thread."""

    lib = None

    def __del__(self):
        try:
            for c in self.values():
                self.lib.libdeflate_free_compressor(c)
        except Exception:  # interpreter shutdown
            pass


def _tls_compressor(lib, level):
    d = getattr(_TLS, "c", None)
    if d is None:
        d = _TLS.c = _Compressors()
        d.lib = lib
    if level not in d:
        d[level] = lib.libdeflate_alloc_compressor(level)
    return d[level]


def deflate_raw(addr: int, n: int, level: int) -> bytes:
    """Returns the raw deflate of n bytes at address addr.

    libdeflate is used when available, else zlib.
    """
    lib = libdeflate()
    if lib:
        buf = getattr(_TLS, "out", None)
        if buf is None:
            buf = _TLS.out = ctypes.create_string_buffer(70000)
        k = lib.libdeflate_deflate_compress(_tls_compressor(lib, level), addr,
                                            n, ctypes.addressof(buf), 70000)
        if k:
            return buf.raw[:k]
    c = zlib.compressobj(min(level, 9), zlib.DEFLATED, -15)
    return c.compress(ctypes.string_at(addr, n)) + c.flush()


def _block(addr: int, n: int, level: int, crc: int) -> bytes:
    data = deflate_raw(addr, n, level)
    return b"".join(
        (_HDR, struct.pack("<H",
                           len(data) + 25), data, struct.pack("<II", crc, n)))


@jit
def _cuts(rowoff, nrows, limit):
    """Returns greedy block cuts at row boundaries.

    Block k holds rows cut[k]..cut[k+1] with at most `limit` bytes.
    """
    cut = np.empty(nrows + 2, dtype=np.int64)
    cut[0] = 0
    nc = 1
    r0 = 0
    r = 0
    while r < nrows:
        if rowoff[r + 1] - rowoff[r0] > limit:
            cut[nc] = r
            nc += 1
            r0 = r
        r += 1
    cut[nc] = nrows
    nc += 1
    return cut[:nc]


def compress_bytes(data: bytes | memoryview,
                   level: int | None = None) -> list[bytes]:
    """Returns the BGZF blocks of `data` (no EOF block).

    The blocks are cut every BLOCK uncompressed bytes.
    """
    level = DEFAULT_LEVEL if level is None else level
    arr = np.frombuffer(data, np.uint8)
    base, mv = arr.ctypes.data, memoryview(arr)
    return [
        _block(base + a,
               min(arr.size, a + BLOCK) - a, level, zlib.crc32(mv[a:a + BLOCK]))
        for a in range(0, arr.size, BLOCK)
    ]


def compress_rows(
        buf: np.ndarray,
        rowoff: np.ndarray,
        nrows: int,
        level: int | None = None) -> tuple[list[bytes], np.ndarray, np.ndarray]:
    """Compresses a text segment of `nrows` lines into line-aligned blocks.

    Args:
        buf: The text as uint8 array.
        rowoff: rowoff[k] is the start of line k, rowoff[nrows] the end.
        nrows: Number of lines.
        level: Compression level (default: DEFAULT_LEVEL).

    Returns:
        (blocks, rows per block, in-block byte offset of every row).
    """
    level = DEFAULT_LEVEL if level is None else level
    if nrows == 0:
        return [], np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64)
    cut = _cuts(rowoff, nrows, BLOCK)
    blocks = []
    mv = memoryview(buf)
    base = buf.ctypes.data
    for k in range(len(cut) - 1):
        a, b = int(rowoff[cut[k]]), int(rowoff[cut[k + 1]])
        blocks.append(_block(base + a, b - a, level, zlib.crc32(mv[a:b])))
    brows = np.diff(cut)
    inoff = rowoff[:nrows] - np.repeat(rowoff[cut[:-1]], brows)
    return blocks, brows, inoff.astype(np.int64)


class BgzfWriter:
    """A BGZF file written in order.

    `add` takes the output of `compress_rows` and returns the virtual offset
    (block offset << 16 | offset in the block) of every row start. `close`
    writes the 28-byte EOF block.
    """

    def __init__(self, path: str) -> None:
        """Opens `path` for writing."""
        self.fh = open(path, "wb")
        self.pos = 0

    def write_blocks(self, blocks: list[bytes]) -> None:
        """Appends compressed blocks."""
        self.fh.writelines(blocks)
        self.pos += sum(map(len, blocks))

    def add(self, blocks: list[bytes], brows: np.ndarray,
            inoff: np.ndarray) -> np.ndarray:
        """Appends compressed blocks; returns the virtual row offsets."""
        if not blocks:
            return np.zeros(0, dtype=np.int64)
        sizes = np.fromiter(map(len, blocks), dtype=np.int64, count=len(blocks))
        coff = self.pos + np.concatenate(([0], np.cumsum(sizes)[:-1]))
        self.fh.writelines(blocks)
        self.pos += int(sizes.sum())
        return (np.repeat(coff, brows) << 16) | inoff

    @property
    def end(self) -> int:
        """Virtual offset of the end of the data written so far."""
        return self.pos << 16

    def close(self) -> None:
        """Writes the EOF block and closes the file."""
        if not self.fh.closed:
            self.fh.write(EOF_BLOCK)
            self.fh.close()

    def __enter__(self):
        """Returns the writer."""
        return self

    def __exit__(self, *exc):
        """Closes the writer."""
        self.close()


def write_segments(path: str,
                   jobs: list,
                   make_seg: Callable[[object], tuple],
                   threads: int,
                   level: int | None = None,
                   order: list[int] | None = None) -> tuple[np.ndarray, int]:
    """Writes the text of make_seg(job) for every job, in job order.

    Deflation runs in `threads` threads.

    Args:
        path: Output file.
        jobs: Job descriptors.
        make_seg: Takes a job, returns (buf, rowoff, nrows).
        threads: Deflate threads.
        level: Compression level.
        order: Indices of the jobs in the order they are started, when it
            differs
            from the output order (a job may wait for a result of earlier ones).

    Returns:
        (virtual offset of every row start in job order, virtual offset of the
            end
        of the data).
    """
    voffs = []

    def work(job):
        """Returns the compressed blocks of one job."""
        buf, rowoff, nrows = make_seg(job)
        return compress_rows(buf, rowoff, nrows, level)

    with BgzfWriter(path) as w:
        if jobs:
            with ThreadPoolExecutor(max(1, threads)) as ex:
                if order is None:
                    results = ex.map(work, jobs)
                else:
                    futs = {i: ex.submit(work, jobs[i]) for i in order}
                    results = (futs.pop(i).result() for i in range(len(jobs)))
                for blocks, brows, inoff in results:
                    voffs.append(w.add(blocks, brows, inoff))
        final = w.end
    return (np.concatenate(voffs)
            if voffs else np.zeros(0, dtype=np.int64)), final


@jit
def _index_contig(beg, end, vstart, vlast, run_bin, run_u, run_v, lin):
    """One pass over the records of a contig (sorted by START).

    run_*: one entry per run of consecutive records in the same UCSC bin (min
    shift 14, 5 levels): the bin, the virtual offset where the run starts and
    where it ends (the start of the next record; `vlast` for the last run).
    lin[w] = start offset of the first record that touches the 16 kb window w
    (lin is -1 where none does). Returns the number of runs.
    """
    nr = 0
    cur = -1
    for i in range(beg.shape[0]):
        b = beg[i]
        e = end[i] - 1
        w0 = b >> 14
        w1 = e >> 14
        if w0 == w1:
            bn = 4681 + w0
        elif (b >> 17) == (e >> 17):
            bn = 585 + (b >> 17)
        elif (b >> 20) == (e >> 20):
            bn = 73 + (b >> 20)
        elif (b >> 23) == (e >> 23):
            bn = 9 + (b >> 23)
        elif (b >> 26) == (e >> 26):
            bn = 1 + (b >> 26)
        else:
            bn = 0
        if bn != cur:
            if nr > 0:
                run_v[nr - 1] = vstart[i]
            run_bin[nr] = bn
            run_u[nr] = vstart[i]
            nr += 1
            cur = bn
        for w in range(w0, max(w0, w1) + 1):
            if lin[w] < 0:
                lin[w] = vstart[i]
    if nr > 0:
        run_v[nr - 1] = vlast
    return nr


def _compress_binning(chunks: dict) -> None:
    """Applies htslib compress_binning to {bin: [[start, end], ...]}.

    The dict is the bins of one contig. Bins are visited from the finest level
    up; a bin whose chunks cover less than 65,536 bytes of compressed data
    (block address of the last end minus block address of the first start) moves
    into its parent bin if that bin exists. Then the chunks of each bin that
    share a BGZF block are joined.
    """
    for lvl in range(5, 0, -1):
        lo, hi = ((1 << (3 * lvl)) - 1) // 7, ((1 << (3 * lvl + 3)) - 1) // 7
        for b in sorted(k for k in chunks if lo <= k < hi):
            lst = chunks[b]
            if lvl < 5 and len(lst) > 1:
                lst.sort()
            if (lst[-1][1] >> 16) - (lst[0][0] >> 16) < 0x10000 and (
                (b - 1) >> 3) in chunks:
                chunks[(b - 1) >> 3].extend(lst)
                del chunks[b]
    for lst in chunks.values():
        m = 0
        for k in range(1, len(lst)):
            if lst[m][1] >> 16 >= lst[k][0] >> 16:
                lst[m][1] = max(lst[m][1], lst[k][1])
            else:
                m += 1
                lst[m] = lst[k]
        del lst[m + 1:]


def tabix_index_bytes(names: list[str], contigs: list, vfinal: int) -> bytes:
    """Bytes of the (uncompressed) tabix index of a sorted BED-like file.

    names: contig names in file order; contigs[i] = (START, END, virtual offset
    of the start of each record) arrays of the records of contig i (the offsets
    are what `BgzfWriter.add` returns; an empty or missing entry for a contig
    without records); vfinal: virtual offset of the end of the data. Follows
    htslib hts_idx_push: one chunk per run of records in the same bin, bins that
    span less than 64 KiB of compressed data move into an existing parent bin,
    chunks of a bin that share a BGZF block are merged, the linear index holds
    the start offset of the first record touching each 16 kb window, and the
    pseudo-bin 37450 carries the offset range and the record count of the
    contig. Records may span several windows.
    """
    out = [
        struct.pack("<4s8i", b"TBI\x01", len(names), 0x10000, 1, 2, 3, 35, 0,
                    sum(len(x.encode()) + 1 for x in names)),
        b"".join(x.encode() + b"\0" for x in names)
    ]
    # the end of a contig's last record is the start of the first record of the
    # next contig that has one
    nxt = [vfinal] * (len(names) + 1)
    for t in range(len(names) - 1, -1, -1):
        nxt[t] = int(contigs[t][2][0]
                    ) if t < len(contigs) and contigs[t] is not None and len(
                        contigs[t][2]) else nxt[t + 1]
    off_beg = 0
    for t in range(len(names)):
        c = contigs[t] if t < len(contigs) else None
        chunks = {}
        if c is not None and len(c[2]):
            beg, end, vstart = (np.ascontiguousarray(x, np.int64) for x in c)
            m = len(vstart)
            vlast = nxt[t + 1]
            nw = (max(int(end.max()) - 1, int(beg.max())) >> 14) + 1
            lin = np.full(nw, -1, np.int64)
            rb, ru, rv = (np.empty(m, np.int64) for _ in range(3))
            nr = _index_contig(beg, end, vstart, vlast, rb, ru, rv, lin)
            for bn, u, v in zip(rb[:nr].tolist(), ru[:nr].tolist(),
                                rv[:nr].tolist()):
                chunks.setdefault(bn, []).append([u, v])
            _compress_binning(chunks)
            chunks[37450] = [[off_beg, vlast], [m, 0]]
            off_beg = vlast
            have = lin >= 0
            # windows without a record take the offset of the next window that
            # has one (htslib)
            if not have.all():
                nx = np.where(have, np.arange(nw), nw)
                nx = np.minimum.accumulate(nx[::-1])[::-1]
                ok = nx < nw
                lin[~have & ok] = lin[nx[~have & ok]]
            ioff = lin.astype(np.uint64)
        else:
            ioff = np.zeros(0, np.uint64)
        out.append(struct.pack("<i", len(chunks)))
        for bid in sorted(chunks):
            out.append(struct.pack("<Ii", bid, len(chunks[bid])))
            out.append(
                np.asarray(chunks[bid], np.uint64).astype("<u8").tobytes())
        out.append(struct.pack("<i", len(ioff)))
        out.append(ioff.astype("<u8").tobytes())
    out.append(struct.pack("<Q", 0))  # n_no_coor, as hts_idx_save writes it
    return b"".join(out)


def write_tbi(path: str,
              names: list[str],
              contigs: list,
              vfinal: int,
              level: int | None = None) -> bool:
    """Writes path.tbi for the sorted, bgzip-compressed BED `path`.

    See tabix_index_bytes for the arguments.

    Returns:
        False, writing nothing, when there is no record or an END is beyond what
        a tabix (not csi) index holds; True otherwise.
    """
    ends = [int(np.max(c[1])) for c in contigs if c is not None and len(c[1])]
    if not ends or max(ends) > MAX_TBI_POS:
        return False
    with open(path + ".tbi", "wb") as fh:
        fh.writelines(
            compress_bytes(tabix_index_bytes(names, contigs, vfinal), level))
        fh.write(EOF_BLOCK)
    return True


def _read_tbi(path: str):
    """Returns (inflated index bytes, contig names, offset of the first ref)."""
    with gzip.open(path + ".tbi", "rb") as fh:
        buf = fh.read()
    magic, n_ref, _fmt, _cs, _cb, _ce, _meta, _skip, l_nm = struct.unpack_from(
        "<4s8i", buf, 0)
    if magic != b"TBI\x01":
        raise ValueError(f"not a tabix index: {path}.tbi")
    names = buf[36:36 + l_nm].split(b"\0")[:n_ref]
    return buf, names, 36 + l_nm


def tbi_extents(path: str) -> dict[str, int]:
    """Returns the upper bound of the covered coordinate per contig.

    It is the number of 16 kb linear-index windows of `path`.tbi times 16384.
    """
    buf, names, off = _read_tbi(path)
    n_ref = len(names)
    out = {}
    for i in range(n_ref):
        (n_bin,) = struct.unpack_from("<i", buf, off)
        off += 4
        for _ in range(n_bin):
            _bin, n_chunk = struct.unpack_from("<Ii", buf, off)
            off += 8 + 16 * n_chunk
        (n_intv,) = struct.unpack_from("<i", buf, off)
        off += 4 + 8 * n_intv
        out[names[i].decode()] = n_intv * 16384
    return out


@jit
def _scan_bins(buf, off, n_bin):
    """Walks the `n_bin` bins of one tabix reference at byte `off` of `buf`.

    Args:
        buf: The inflated index as uint8 array.
        off: Byte offset of the first bin.
        n_bin: Number of bins.

    Returns:
        (offset after the bins, smallest chunk begin, largest chunk end) over
            all
        bins except the pseudo-bin 37450.
    """
    first = np.int64(1) << np.int64(62)
    last = np.int64(0)
    for _ in range(n_bin):
        b = np.int64(buf[off]) | (np.int64(buf[off + 1]) << 8) | (
            np.int64(buf[off + 2]) << 16) | (np.int64(buf[off + 3]) << 24)
        nc = np.int64(buf[off + 4]) | (np.int64(buf[off + 5]) << 8) | (
            np.int64(buf[off + 6]) << 16) | (np.int64(buf[off + 7]) << 24)
        off += 8
        for _c in range(nc):
            if b != 37450:
                beg = np.int64(0)
                end = np.int64(0)
                for k in range(8):
                    beg |= np.int64(buf[off + k]) << np.int64(8 * k)
                    end |= np.int64(buf[off + 8 + k]) << np.int64(8 * k)
                if beg < first:
                    first = beg
                if end > last:
                    last = end
            off += 16
    return off, first, last


def tbi_index(path: str) -> dict:
    """Read `path`.tbi: {contig: (lin_lo, lin_hi, v_first, v_last)}.

    `lin_lo[w]` is a lower bound and `lin_hi[w]` an upper bound of the virtual
    offset of the first record that overlaps the 16 kb window w (the linear
    index with its empty windows filled from the left and from the right;
    `lin_hi` stays 2^64 - 1 after the last record); `v_first` and `v_last` are
    the first and the end virtual offset of all records of the contig (from the
    bin chunks, pseudo-bin excluded). Contigs without records are left out.
    """
    buf, names, off = _read_tbi(path)
    n_ref = len(names)
    arr = np.frombuffer(buf, np.uint8)
    out = {}
    big = np.iinfo(np.uint64).max
    for i in range(n_ref):
        (n_bin,) = struct.unpack_from("<i", buf, off)
        off, first, last = _scan_bins(arr, off + 4, n_bin)
        (n_intv,) = struct.unpack_from("<i", buf, off)
        off += 4
        lin = np.frombuffer(buf, "<u8", n_intv, off)
        off += 8 * n_intv
        if last == 0:
            continue
        lo = np.maximum.accumulate(lin)
        hi = np.minimum.accumulate(np.where(lin == 0, big, lin)[::-1])[::-1]
        out[names[i].decode()] = (lo, hi, int(first), int(last))
    return out
