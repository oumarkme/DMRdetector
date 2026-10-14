"""Methylation BED input: header detection, byte-level parse, write, split.

The fast path takes ordinary files only (tab or space separated, digit-only
START, END, coverage and N_mod); for anything else `read_filtered` returns None
and the caller runs the pandas code. Output rows follow `sort -k1,1 -k2,2n`.
"""

from __future__ import annotations

import re
import threading

import numpy as np
import pandas as pd

from dmrdetector.fileio import bed
from dmrdetector.fileio import bgzf
from dmrdetector.fileio import numparse
from dmrdetector.util.jit import jit


_TAB, _NL, _SP = 9, 10, 32
# ASCII codes of "#" and of "chrom" (the start of a header line)
_HASH = 35
_C, _H, _R, _O, _M = 99, 104, 114, 111, 109
SEG_ROWS = 200_000
CHUNK_NAMES = 65536

_NEED_COLS = {"5col": 5, "modkit": 12, "pbcpg": 7}


def scan_header(path: str, fmt: str = "auto") -> tuple[str, int, int]:
    """Detects the format from the first data line, as v2.2 does.

    Lines that start with `#` or `chrom`, and blank lines, are header lines.
    Fields split on whitespace runs (awk NF): 18 or more columns is modkit, 7 to
    9 is pbcpg, else 5col.

    Args:
        path: Methylation BED, plain or gzip.
        fmt: "auto" or a format that overrides the detected one.

    Returns:
        (format, number of columns of the first data line, number of leading
        header lines). The column count is 0 for a file without data lines.

    Raises:
        ValueError: `fmt` is unknown, or the first data line has too few
            columns for it.
    """
    nhead, ncols, detected = 0, 0, "5col"
    with bed.open_text(path) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("chrom"):
                nhead += 1
                continue
            ncols = len(line.split())
            if ncols == 0:
                nhead += 1
                continue
            if ncols >= 18:
                detected = "modkit"
            elif 7 <= ncols <= 9:
                detected = "pbcpg"
            break
    if fmt == "auto":
        return detected, ncols, nhead
    if fmt not in _NEED_COLS:
        raise ValueError(
            f"unknown input format {fmt!r} (auto, 5col, modkit, pbcpg)")
    if ncols and ncols < _NEED_COLS[fmt]:
        raise ValueError(
            f"{path}: input_format {fmt} needs at least "
            f"{_NEED_COLS[fmt]} columns, the first data line has {ncols}")
    return fmt, ncols, nhead


def read_methylation_bed(path: str,
                         fmt: str = "auto") -> tuple[pd.DataFrame, str]:
    """Reads a methylation BED in any of the three formats with pandas.

    Leading header lines (`#...`, `chrom...`) are dropped; columns split on
    whitespace like awk.

    Args:
        path: Plain or gzip file.
        fmt: "auto" or a forced format (see scan_header).

    Returns:
        (frame with CHROM (str), START, END (int64), COVERAGE, N_MOD (float64),
        format name).
    """
    fmt, ncols, nhead = scan_header(path, fmt)
    empty = pd.DataFrame({
        "CHROM": pd.Series(dtype=object),
        "START": pd.Series(dtype=np.int64),
        "END": pd.Series(dtype=np.int64),
        "COVERAGE": pd.Series(dtype=np.float64),
        "N_MOD": pd.Series(dtype=np.float64)
    })
    if ncols == 0:
        return empty, fmt
    cov_i, mod_i = {"modkit": (9, 11), "pbcpg": (5, 6), "5col": (3, 4)}[fmt]
    df = pd.read_csv(path,
                     sep=r"\s+",
                     header=None,
                     usecols=[0, 1, 2, cov_i, mod_i],
                     skiprows=nhead,
                     comment="#",
                     dtype={
                         0: str,
                         1: np.int64,
                         2: np.int64,
                         cov_i: np.float64,
                         mod_i: np.float64
                     },
                     float_precision="round_trip",
                     compression="gzip" if bed.is_gzip(path) else None)
    df = df[[0, 1, 2, cov_i, mod_i]]
    df.columns = ["CHROM", "START", "END", "COVERAGE", "N_MOD"]
    df["CHROM"] = df["CHROM"].astype(object)
    return df.reset_index(drop=True), fmt


@jit
def _parse(buf, first, cov_i, mod_i, ncols, strict, name_off, name_len):
    """Parses whole lines of `buf`.

    `first` skips leading header lines (`#...`, `chrom...`, blank). `strict`
    accepts tab separators only, exactly one tab between fields and canonical
    integers (no leading zeros), as the filter step writes them.

    Returns:
        (flag, n, code, start, end, cov, mod, nnames, nskip). A nonzero flag
            means
        the bytes need the generic code.
    """
    n = buf.shape[0]
    cap = 1
    for i in range(n):
        if buf[i] == _NL:
            cap += 1
    code = np.empty(cap, dtype=np.int32)
    start = np.empty(cap, dtype=np.int64)
    end = np.empty(cap, dtype=np.int64)
    cov = np.empty(cap, dtype=np.float64)
    mod = np.empty(cap, dtype=np.float64)
    nrow = 0
    nnames = 0
    nskip = 0
    last = -1
    in_head = first and not strict
    tokS = np.empty(5, dtype=np.int64)
    tokE = np.empty(5, dtype=np.int64)
    v = np.empty(4, dtype=np.int64)
    p = 0
    while p < n:
        e = p
        while e < n and buf[e] != _NL:
            e += 1
        if in_head:
            # header / blank lines of the file start
            q = p
            while q < e and (buf[q] == _SP or buf[q] == _TAB):
                q += 1
            skip = False
            if q == e:
                skip = True
            elif buf[p] == _HASH:
                skip = True
            elif e - p >= 5 and buf[p] == _C and buf[p + 1] == _H and buf[
                    p + 2] == _R and buf[p + 3] == _O and buf[p + 4] == _M:
                skip = True
            if skip:
                nskip += 1
                p = e + 1
                continue
            in_head = False
        if e == p or buf[p] == _SP or buf[p] == _TAB:
            return (1, nrow, code, start, end, cov, mod, nnames, nskip)
        f = 0
        q = p
        while q < e:
            fs = q
            while q < e:
                c = buf[q]
                if c == _SP or c == _TAB:
                    break
                if c < 33 or c > 126 or c == _HASH or c == 34:
                    return (1, nrow, code, start, end, cov, mod, nnames, nskip)
                q += 1
            if f == 0:
                tokS[0] = fs
                tokE[0] = q
            elif f == 1:
                tokS[1] = fs
                tokE[1] = q
            elif f == 2:
                tokS[2] = fs
                tokE[2] = q
            if f == cov_i:
                tokS[3] = fs
                tokE[3] = q
            if f == mod_i:
                tokS[4] = fs
                tokE[4] = q
            f += 1
            if q < e:
                # separator run
                if strict:
                    if buf[q] != _TAB:
                        return (1, nrow, code, start, end, cov, mod, nnames,
                                nskip)
                    q += 1
                    if q >= e or buf[q] == _TAB or buf[q] == _SP:
                        return (1, nrow, code, start, end, cov, mod, nnames,
                                nskip)
                else:
                    while q < e and (buf[q] == _SP or buf[q] == _TAB):
                        q += 1
                    if q >= e:  # trailing whitespace
                        return (1, nrow, code, start, end, cov, mod, nnames,
                                nskip)
        if f != ncols:
            return (1, nrow, code, start, end, cov, mod, nnames, nskip)
        # numbers
        for t in range(1, 5):
            s = tokS[t]
            x = tokE[t]
            ln = x - s
            maxd = 18 if t < 3 else 15
            if ln < 1 or ln > maxd:
                return (1, nrow, code, start, end, cov, mod, nnames, nskip)
            if (strict and ln > 1 and buf[s] == 48) or buf[s] == 45:
                return (1, nrow, code, start, end, cov, mod, nnames, nskip)
            # digits only: a sign or any other byte is refused
            st, acc = numparse.parse_int(buf, s, x)
            if st != 0:
                return (1, nrow, code, start, end, cov, mod, nnames, nskip)
            v[t - 1] = acc
        start[nrow] = v[0]
        end[nrow] = v[1]
        cov[nrow] = float(v[2])
        mod[nrow] = float(v[3])
        # contig id (chunk-local)
        cs = tokS[0]
        ce = tokE[0]
        cid = -1
        if last >= 0 and name_len[last] == ce - cs:
            same = True
            for x in range(ce - cs):
                if buf[name_off[last] + x] != buf[cs + x]:
                    same = False
                    break
            if same:
                cid = last
        if cid < 0:
            for m in range(nnames):
                if name_len[m] == ce - cs:
                    same = True
                    for x in range(ce - cs):
                        if buf[name_off[m] + x] != buf[cs + x]:
                            same = False
                            break
                    if same:
                        cid = m
                        break
        if cid < 0:
            if nnames >= name_off.shape[0]:
                return (1, nrow, code, start, end, cov, mod, nnames, nskip)
            name_off[nnames] = cs
            name_len[nnames] = ce - cs
            cid = nnames
            nnames += 1
        last = cid
        code[nrow] = cid
        nrow += 1
        p = e + 1
    return 0, nrow, code, start, end, cov, mod, nnames, nskip


@jit
def _drop_mask(start, end, code, hid, het_off, het_pos, ndrop_out):
    """Flags the rows that hold a het-SNP position p (start <= p < end).

    This is sorted_point_overlap of the pandas code: the rows where
    searchsorted(p, end) > searchsorted(p, start) and end > start.
    """
    n = start.shape[0]
    drop = np.zeros(n, dtype=np.bool_)
    cnt = 0
    for i in range(n):
        h = hid[code[i]]
        if h < 0:
            continue
        a = het_off[h]
        b = het_off[h + 1]
        if a == b:
            continue
        s = start[i]
        t = end[i]
        if t <= s or het_pos[b - 1] < s or het_pos[a] >= t:
            continue
        lo = a
        hi = b
        while lo < hi:  # first p >= s
            mid = (lo + hi) >> 1
            if het_pos[mid] < s:
                lo = mid + 1
            else:
                hi = mid
        if lo < b and het_pos[lo] < t:
            drop[i] = True
            cnt += 1
    ndrop_out[0] = cnt
    return drop


@jit
def _blocks(gid, start):
    """Contig blocks of rows.

    Returns (ok, nblocks, block_gid, block_a, block_b). ok is False when a
    contig appears in two blocks or starts are not strictly increasing inside a
    block (the generic sort then applies).
    """
    n = gid.shape[0]
    bg = np.empty(n + 1, dtype=np.int32)
    ba = np.empty(n + 1, dtype=np.int64)
    bb = np.empty(n + 1, dtype=np.int64)
    nb = 0
    if n == 0:
        return True, 0, bg, ba, bb
    ba[0] = 0
    bg[0] = gid[0]
    nb = 1
    for i in range(1, n):
        if gid[i] != gid[i - 1]:
            bb[nb - 1] = i
            bg[nb] = gid[i]
            ba[nb] = i
            nb += 1
            if nb > 100000:
                return False, nb, bg, ba, bb
        elif start[i] <= start[i - 1]:
            return False, nb, bg, ba, bb
    bb[nb - 1] = n
    return True, nb, bg, ba, bb


@jit
def _format(gid, start, end, cov, nmod, a, b, nm_buf, nm_off, nm_len, out,
            rowoff):
    """Writes rows a..b as `CHROM START END COV NMOD` lines.

    Integral values use the awk number format. rowoff[k] is the byte offset of
    row a+k and rowoff[b-a] the total.

    Returns:
        The total bytes, or -1 when a value is not a non-negative integer below
        1e15 (the caller formats that segment with the generic code).
    """
    p = 0
    for r in range(a, b):
        x = cov[r]
        y = nmod[r]
        if not (x >= 0.0 and x < 1e15 and y >= 0.0 and y < 1e15 and
                x == np.floor(x) and y == np.floor(y)):
            return -1
        rowoff[r - a] = p
        g = gid[r]
        o = nm_off[g]
        for i in range(nm_len[g]):
            out[p + i] = nm_buf[o + i]
        p += nm_len[g]
        out[p] = _TAB
        p = numparse.put_int(out, p + 1, start[r])
        out[p] = _TAB
        p = numparse.put_int(out, p + 1, end[r])
        out[p] = _TAB
        p = numparse.put_int(out, p + 1, np.int64(x))
        out[p] = _TAB
        p = numparse.put_int(out, p + 1, np.int64(y))
        out[p] = _NL
        p += 1
    rowoff[b - a] = p
    return p


class Table:
    """Rows of a methylation BED as parallel arrays.

    `blocks` is [(gid, a, b)] in output order (C-locale contig names).
    """

    def __init__(self,
                 names: list[str],
                 gid: np.ndarray,
                 start: np.ndarray,
                 end: np.ndarray,
                 cov: np.ndarray,
                 nmod: np.ndarray,
                 n_in: int = 0,
                 n_drop: int = 0) -> None:
        """Stores the arrays; `names[gid]` is the contig name of a row."""
        self.names, self.gid = names, gid
        self.start, self.end, self.cov, self.nmod = start, end, cov, nmod
        self.n_in, self.n_drop = n_in, n_drop
        self.blocks = None

    def __len__(self):
        """Returns the number of rows."""
        return len(self.gid)

    def order_blocks(self) -> bool:
        """Sets `blocks`; returns False when the generic sort is needed.

        The rows of every contig must be contiguous and strictly increasing in
        START.
        """
        ok, nb, bg, ba, bb = _blocks(self.gid, self.start)
        if not ok:
            return False
        gids = bg[:nb].tolist()
        if len(set(gids)) != len(gids):
            return False
        self.blocks = sorted(zip(gids, ba[:nb].tolist(), bb[:nb].tolist()),
                             key=lambda t: self.names[t[0]])
        return True

    def add_name(self, name: str) -> int:
        """Returns the gid of a contig name, adding it when new."""
        if name in self.names:
            return self.names.index(name)
        self.names = list(self.names) + [name]
        return len(self.names) - 1


def read_filtered(path: str,
                  fmt: str,
                  ncols: int,
                  min_cov: float,
                  het: dict,
                  threads: int,
                  strict: bool = False,
                  cols: tuple | None = None) -> Table | None:
    """Parses a methylation BED; drops het-SNP rows and low coverage.

    Args:
        path: Methylation BED, plain or gzip.
        fmt: "5col", "modkit" or "pbcpg".
        ncols: Number of columns of the data lines.
        min_cov: Keep rows with coverage above this.
        het: {contig: sorted het-SNP positions}.
        threads: Reader threads.
        strict: Accept the output form of the filter step only.
        cols: (coverage column, N_mod column), 0-based; default from `fmt`.

    Returns:
        A Table (`n_in` counts every parsed row, `n_drop` the rows on a het SNP
        whatever their coverage), or None when the file needs the generic code.
    """
    cov_i, mod_i = cols or {
        "modkit": (9, 11),
        "pbcpg": (5, 6),
        "5col": (3, 4)
    }[fmt]
    hnames = list(het)
    hetid = {n: i for i, n in enumerate(hnames)}
    het_off = np.zeros(len(hnames) + 1, dtype=np.int64)
    for i, n in enumerate(hnames):
        het_off[i + 1] = het_off[i] + het[n].size
    het_pos = np.concatenate([het[n] for n in hnames]).astype(
        np.int64) if hnames else np.zeros(0, dtype=np.int64)
    ids: dict[str, int] = {}
    lock = threading.Lock()
    min_cov = float(min_cov)

    def task(arr, first):
        off = np.zeros(CHUNK_NAMES, dtype=np.int64)
        ln = np.zeros(CHUNK_NAMES, dtype=np.int64)
        flag, n, code, start, end, cov, mod, nn, _ = _parse(
            arr, first, cov_i, mod_i, ncols, strict, off, ln)
        if flag:
            return ()  # anything the byte reader does not accept
        names = [bytes(arr[off[i]:off[i] + ln[i]]).decode() for i in range(nn)]
        with lock:
            gl = np.array([ids.setdefault(nm, len(ids)) for nm in names],
                          dtype=np.int32)
        hid = np.array([hetid.get(nm, -1) for nm in names], dtype=np.int32)
        code, start, end, cov, mod = code[:n], start[:n], end[:
                                                              n], cov[:
                                                                      n], mod[:
                                                                              n]
        nd = np.zeros(1, dtype=np.int64)
        drop = _drop_mask(start, end, code, hid, het_off, het_pos,
                          nd) if het_pos.size else np.zeros(n, dtype=bool)
        keep = ~drop & (cov > min_cov)
        return n, int(nd[0]), gl[code[keep]] if nn else code[keep], start[
            keep], end[keep], cov[keep], mod[keep]

    res = bed.map_text(path, task, threads)
    if any(len(r) == 0 for r in res):
        return None
    names = [None] * len(ids)
    for k, v in ids.items():
        names[v] = k
    if not res:
        z = np.zeros(0)
        return Table(names, z.astype(np.int32), z.astype(np.int64),
                     z.astype(np.int64), z, z)
    cat = [np.concatenate([r[i] for r in res]) for i in range(2, 7)]
    n_in, n_drop = sum(r[0] for r in res), sum(r[1] for r in res)
    del res
    return Table(names, cat[0].astype(np.int32, copy=False), cat[1], cat[2],
                 cat[3], cat[4], n_in, n_drop)


def _generic_text(tab_arrays, names, a, b):
    """Returns the text of rows a..b by the generic awk-number formatter."""
    gid, start, end, cov, nmod = tab_arrays
    cols = [[names[g] for g in gid[a:b].tolist()],
            start[a:b].astype(str).tolist(), end[a:b].astype(str).tolist(),
            bed.awk_number(cov[a:b]),
            bed.awk_number(nmod[a:b])]
    raw = bed.join_columns(cols).encode()
    arr = np.frombuffer(raw, dtype=np.uint8)
    nl = np.flatnonzero(arr == 10)
    rowoff = np.concatenate(([0], nl + 1)).astype(np.int64)
    return arr, rowoff, b - a


def _make_seg_fn(names, nm):
    nm_buf, nm_off, nm_len = nm
    maxname = int(nm_len.max()) if len(nm_len) else 0

    def make(job):
        src, a, b = job
        n = b - a
        gid, start, end, cov, nmod = src
        out = np.empty(n * (maxname + 76) + 8, dtype=np.uint8)
        rowoff = np.empty(n + 1, dtype=np.int64)
        r = _format(gid, start, end, cov, nmod, a, b, nm_buf, nm_off, nm_len,
                    out, rowoff)
        if r < 0:
            return _generic_text(src, names, a, b)
        return out[:r], rowoff, n

    return make


def write_blocks(path: str,
                 names: list,
                 blocks: list,
                 threads: int,
                 level: int,
                 nm: tuple | None = None) -> int:
    """Writes blocks as bgzip text plus a tabix index.

    Args:
        path: Output file.
        names: Contig names by gid.
        blocks: [(gid, src, a, b)], one contig each, in output order; src is the
            (gid, start, end, cov, nmod) arrays.
        threads: Writer threads.
        level: Compression level.
        nm: Name table (default from `names`).

    Returns:
        The number of rows.
    """
    nm = nm or _names_table(names)
    jobs, spans = [], []
    for g, src, a, b in blocks:
        for s in range(a, b, SEG_ROWS):
            jobs.append((src, s, min(b, s + SEG_ROWS)))
        spans.append((g, src, a, b))
    nrows = sum(b - a for _, _, a, b in spans)
    voff, final = bgzf.write_segments(path, jobs, _make_seg_fn(names, nm),
                                      threads, level)
    if nrows == 0:
        return 0
    contigs, o = [], 0
    for _, src, a, b in spans:
        contigs.append((src[1][a:b], src[2][a:b], voff[o:o + b - a]))
        o += b - a
    if not bgzf.write_tbi(path, [names[g] for g, _, _, _ in spans], contigs,
                          final, level):
        bed.tabix_bed(path)  # an END beyond the 2^29 limit of a tabix index
    return nrows


def _names_table(names):
    enc = [n.encode() for n in names]
    ln = np.array([len(e) for e in enc], dtype=np.int64)
    off = np.concatenate(([0], np.cumsum(ln)[:-1])).astype(
        np.int64) if len(enc) else np.zeros(0, dtype=np.int64)
    return np.frombuffer(b"".join(enc), dtype=np.uint8), off, ln


def write_table(path: str,
                tab: Table,
                threads: int,
                level: int = bgzf.DEFAULT_LEVEL) -> int:
    """Writes a Table as bgzip plus tabix BED in `sort -k1,1 -k2,2n` order.

    Returns:
        The number of rows.
    """
    src = (tab.gid, tab.start, tab.end, tab.cov, tab.nmod)
    return write_blocks(path, tab.names,
                        [(g, src, a, b) for g, a, b in tab.blocks], threads,
                        level)


def split_table(tab: Table,
                par: dict,
                out_prefix: str,
                auto_re: re.Pattern,
                x_re: re.Pattern,
                y_re: re.Pattern,
                par2_shift: int,
                threads: int,
                level: int = bgzf.DEFAULT_LEVEL) -> dict | None:
    """Writes {prefix}_autosomal, _chrX and _chrY .bed.gz from a Table.

    The selection, PAR folding and sort rules are those of the generic split
    (steps.preprocess._split_generic).

    Returns:
        {part: rows}, or None when the generic code is needed.
    """
    names = list(tab.names)
    src0 = (tab.gid, tab.start, tab.end, tab.cov, tab.nmod)
    core, xblocks, yblocks = [], [], []
    for g, a, b in tab.blocks:
        nm = names[g]
        if auto_re.match(nm):
            core.append((g, src0, a, b))
        elif x_re.match(nm):
            xblocks.append((g, a, b))
        elif y_re.match(nm):
            yblocks.append((g, a, b))

    def split_par(blks):
        """Returns (non-PAR rows, PAR rows repeated by multiplicity)."""
        npar_sel, par_rows = [], []
        for g, a, b in blks:
            s, e = tab.start[a:b], tab.end[a:b]
            nm = names[g]
            cnt = bed.overlap_counts(s, e, *par[nm]) if nm in par else np.zeros(
                b - a, dtype=np.int64)
            npar_sel.append(a + np.flatnonzero(cnt == 0))
            par_rows.append(np.repeat(np.arange(a, b), cnt))
        z = np.zeros(0, dtype=np.int64)
        return (np.concatenate(npar_sel) if npar_sel else
                z), (np.concatenate(par_rows) if par_rows else z)

    xnp, xpar = split_par(xblocks)
    ynp, ypar = split_par(yblocks)

    # PAR block: chrX PAR rows then chrY PAR rows (PAR2 shifted), sorted by
    # START, ties by the generic rule
    parblock = None
    if xpar.size or ypar.size:
        rows = np.concatenate((xpar, ypar))
        st = tab.start[rows].copy()
        en = tab.end[rows].copy()
        p2 = np.concatenate((np.zeros(xpar.size, dtype=bool), tab.start[ypar]
                             >= 10_000_000))
        st[p2] += par2_shift
        en[p2] += par2_shift
        part = pd.DataFrame({
            "CHROM": "chrX",
            "START": st,
            "END": en,
            "COV": bed.awk_number(tab.cov[rows]),
            "NMOD": bed.awk_number(tab.nmod[rows]),
            "ROW": np.arange(rows.size)
        })
        part = part.iloc[bed.bed_sort_order(
            part["CHROM"].to_numpy(dtype=object), part["START"].to_numpy())]
        part = tie_sort(part)
        order = part["ROW"].to_numpy().astype(np.int64)
        gx = tab.add_name("chrX")
        names = list(tab.names)
        parblock = (gx, (np.full(rows.size, gx,
                                 dtype=np.int32)[order], st[order], en[order],
                         tab.cov[rows][order], tab.nmod[rows][order]))

    def sel_src(sel):
        return (tab.gid[sel], tab.start[sel], tab.end[sel], tab.cov[sel],
                tab.nmod[sel])

    nm = _names_table(names)
    # autosomal: core contigs and the chrX PAR block in contig-name order
    auto = list(core)
    if parblock is not None:
        auto.append((parblock[0], parblock[1], 0, len(parblock[1][0])))
    auto.sort(key=lambda t: names[t[0]])
    counts = {}
    counts["autosomal"] = write_blocks(f"{out_prefix}_autosomal.bed.gz", names,
                                       auto, threads, level, nm)
    for part_name, sel in (("chrX", xnp), ("chrY", ynp)):
        blocks = []
        if sel.size:
            s = sel_src(sel)
            gs = s[0]
            bnd = np.flatnonzero(
                np.concatenate(([True], gs[1:] != gs[:-1], [True])))
            for i in range(len(bnd) - 1):
                blocks.append(
                    (int(gs[bnd[i]]), s, int(bnd[i]), int(bnd[i + 1])))
            blocks.sort(key=lambda t: names[t[0]])
        counts[part_name] = write_blocks(f"{out_prefix}_{part_name}.bed.gz",
                                         names, blocks, threads, level, nm)
    return counts


def tie_sort(part: pd.DataFrame) -> pd.DataFrame:
    """Orders tied (CHROM, START) rows by the rest of the line text.

    This is the last resort of `sort`.
    """
    part = part.reset_index(drop=True)
    key = part["CHROM"] + "\t" + part["START"].astype(str)
    dup = key.duplicated(keep=False).to_numpy()
    if not dup.any():
        return part
    pos = np.flatnonzero(dup)
    sub = part.iloc[pos]
    gid = (key.iloc[pos] != key.iloc[pos].shift()).cumsum().to_numpy()
    tup = list(
        zip(gid.tolist(), sub["END"].astype(str), sub["COV"], sub["NMOD"]))
    idx = sorted(range(len(tup)), key=tup.__getitem__)
    out = part.copy()
    out.iloc[pos] = sub.iloc[idx].to_numpy()
    return out
