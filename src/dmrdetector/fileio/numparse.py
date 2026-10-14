"""Number parsing and formatting kernels (numba): decimals, R fwrite text.

Parsing is correctly rounded: up to 15 digits is one IEEE operation, up to 19 a
double-double product, and a tie goes to Python `float()`. A token a kernel
cannot read sets a flag and the caller falls back to pandas. Formatting writes
15 significant digits as R `fwrite` does. Notes: doc/ALGORITHM.md, Fast readers.
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction

import numpy as np

from dmrdetector.fileio import bgzf
from dmrdetector.util.jit import jit


_P10_MIN, _P10_MAX = -345, 308
# 1e0 .. 1e22 are exact doubles
P10_EXACT = np.array([float(10**k) for k in range(23)])


def _dd_table():
    hi = np.empty(_P10_MAX - _P10_MIN + 1)
    lo = np.empty_like(hi)
    for i, k in enumerate(range(_P10_MIN, _P10_MAX + 1)):
        f = Fraction(10)**k
        h = float(f)
        hi[i] = h
        lo[i] = float(f - Fraction(h)) if h != 0.0 else 0.0
    return hi, lo


P10_HI, P10_LO = _dd_table()


def _map(fn, ranges, threads):
    if threads <= 1 or len(ranges) == 1:
        return [fn(*r) for r in ranges]
    with ThreadPoolExecutor(threads) as ex:
        futs = [ex.submit(fn, *r) for r in ranges]
        return [f.result() for f in futs]


@jit
def _parse_float(buf, p, q, p10x, p10hi, p10lo):
    """Returns (status, value) of a correctly rounded decimal.

    Status 0 is ok, 1 means the token is not handled.
    """
    neg = False
    if p < q and buf[p] == 45:
        neg = True
        p += 1
    elif p < q and buf[p] == 43:
        p += 1
    m = np.uint64(0)
    nd = 0
    fd = 0
    anyd = False
    while p < q and buf[p] >= 48 and buf[p] <= 57:
        d = np.uint64(buf[p] - 48)
        anyd = True
        if m != 0 or d != 0:
            if nd >= 19:
                return 1, 0.0
            m = m * np.uint64(10) + d
            nd += 1
        p += 1
    if p < q and buf[p] == 46:
        p += 1
        while p < q and buf[p] >= 48 and buf[p] <= 57:
            d = np.uint64(buf[p] - 48)
            anyd = True
            if m != 0 or d != 0:
                if nd >= 19:
                    return 1, 0.0
                m = m * np.uint64(10) + d
                nd += 1
            fd += 1
            p += 1
    if not anyd:
        return 1, 0.0
    e10 = 0
    if p < q and (buf[p] == 101 or buf[p] == 69):
        p += 1
        eneg = False
        if p < q and buf[p] == 45:
            eneg = True
            p += 1
        elif p < q and buf[p] == 43:
            p += 1
        if p >= q:
            return 1, 0.0
        while p < q and buf[p] >= 48 and buf[p] <= 57:
            e10 = e10 * 10 + np.int64(buf[p] - 48)
            if e10 > 100000:
                return 1, 0.0
            p += 1
        if eneg:
            e10 = -e10
    if p != q:
        return 1, 0.0
    if m == 0:
        return 0, (-0.0 if neg else 0.0)
    ex = e10 - fd
    val = 0.0
    if m < np.uint64(9007199254740992) and ex >= -22 and ex <= 22:
        v = np.float64(np.int64(m))
        if ex >= 0:
            val = v * p10x[ex]
        else:
            val = v / p10x[-ex]
    else:
        if m >= np.uint64(4000000000000000000) or ex < -290 or ex > 280:
            return 1, 0.0
        mi = np.int64(m)
        mh = np.float64(mi)
        ml = np.float64(mi - np.int64(mh))
        th = p10hi[ex + 345]
        tl = p10lo[ex + 345]
        # Dekker product mh * th
        c = 134217729.0 * mh
        ah = c - (c - mh)
        al = mh - ah
        c = 134217729.0 * th
        bh = c - (c - th)
        bl = th - bh
        pr = mh * th
        er = ((ah * bh - pr) + ah * bl + al * bh) + al * bl
        er += mh * tl + ml * th
        val = pr + er
        # the product is exact to about 2^-100 relative; when pr + er lies
        # within 1e-9 of half a unit of the last place from a tie the rounding
        # is not certain: Python float() decides
        lo = er - (val - pr)
        if lo >= 0.0:
            half = (np.nextafter(val, np.inf) - val) * 0.5
        else:
            half = (val - np.nextafter(val, 0.0)) * 0.5
        if abs(abs(lo) - half) <= half * 1e-9:
            return 1, 0.0
    return 0, (-val if neg else val)


@jit
def parse_decimal(buf, p, q):
    """Returns (status, value) of buf[p:q] with the shared 10^k tables."""
    return _parse_float(buf, p, q, P10_EXACT, P10_HI, P10_LO)


@jit
def parse_int(buf, p, q):
    """Returns (status, value) of the digits buf[p:q]; status 0 is ok."""
    neg = False
    if p < q and buf[p] == 45:
        neg = True
        p += 1
    if p >= q or q - p > 18:
        return 1, 0
    v = 0
    while p < q:
        c = buf[p]
        if c < 48 or c > 57:
            return 1, 0
        v = v * 10 + np.int64(c - 48)
        p += 1
    return 0, (-v if neg else v)


@jit
def _parse_stream(buf, lo, hi, row0, cols, kinds, slots, outi, outf, bad,
                  change, ls, p10x, p10hi, p10lo):
    """Parses the lines of buf[lo:hi] (lo at a line start, hi after a newline).

    Line starts go to `ls` and the columns `cols` to outi / outf as in
    _parse_fields. change[row] is 1 when the first field differs from the
    previous line's (not set for the first row of the segment).
    """
    nc = cols.shape[0]
    p = lo
    r = row0
    prev_s = -1
    prev_l = 0
    while p < hi:
        ls[r] = p
        flag = 0
        if buf[p] == 10 or buf[p] == 35:
            flag = 2
            q = p
        else:
            f = 0
            j = 0
            tok = p
            q = p
            c0l = 0
            while True:
                q = tok
                while buf[q] != 9 and buf[q] != 10:
                    q += 1
                if f == 0:
                    c0l = q - tok
                if f == cols[j]:
                    if kinds[j] == 1:
                        st, v = parse_int(buf, tok, q)
                        if st == 0:
                            outi[slots[j], r] = v
                        else:
                            flag = 1
                    else:
                        st, v = _parse_float(buf, tok, q, p10x, p10hi, p10lo)
                        if st == 0:
                            outf[slots[j], r] = v
                        else:
                            flag = 1
                    j += 1
                    if j == nc:
                        break
                if buf[q] == 10:
                    flag = 2
                    break
                f += 1
                tok = q + 1
            if f == 0 and nc > 0 and cols[0] != 0:
                pass
            # chromosome run change (first field)
            if prev_s >= 0:
                if c0l != prev_l:
                    change[r] = 1
                else:
                    for k in range(c0l):
                        if buf[p + k] != buf[prev_s + k]:
                            change[r] = 1
                            break
            if f >= 0:
                prev_s = p
                prev_l = c0l
        while buf[q] != 10:
            q += 1
        bad[r] = flag
        r += 1
        p = q + 1
    return r


def _first_token_len(buf, p):
    q = p
    while buf[q] != 9 and buf[q] != 10:
        q += 1
    return q - p


def read_table(path: str,
               spec: list[tuple[int, str]],
               threads: int = 1) -> tuple | None:
    """Inflates and parses a BGZF text table in two passes.

    Args:
        path: BGZF file.
        spec: [(column, "i" or "f")]; the first column (chromosome) is not
            needed.
        threads: Threads for inflating and parsing.

    Returns:
        (ints, floats, (run_starts, run_names)), empty arrays for an empty file;
        None when the file is not BGZF.

    Raises:
        bgzf.ParseError: A line cannot be handled.
    """
    if bgzf.NO_JIT:
        return None
    buf = bgzf.inflate(path, threads)
    if buf is None:
        return None
    return parse_buffer(buf, spec, threads)


def parse_buffer(buf: np.ndarray,
                 spec: list[tuple[int, str]],
                 threads: int = 1) -> tuple:
    """Parses an inflated table (arguments and result as in read_table)."""
    spec = sorted(spec)
    cols = np.array([c for c, _ in spec], np.int64)
    kinds = np.array([1 if k == "i" else 2 for _, k in spec], np.int64)
    slots = np.zeros(len(spec), np.int64)
    ni = nf = 0
    for j, k in enumerate(kinds):
        slots[j] = ni if k == 1 else nf
        ni += k == 1
        nf += k != 1
    if buf.size == 0:
        return ([np.zeros(0, np.int64) for _ in range(ni)],
                [np.zeros(0) for _ in range(nf)], (np.zeros(0, np.int64), []))
    nseg = threads * 2 if buf.size > 1 << 22 else 1
    cut = [0]
    for i in range(1, nseg):
        c = max(cut[-1], buf.size * i // nseg)
        while c < buf.size and buf[c - 1] != 10:
            c += 1
        cut.append(c)
    cut.append(buf.size)
    segs = [(cut[i], cut[i + 1]) for i in range(nseg) if cut[i + 1] > cut[i]]
    counts = _map(lambda a, b: bgzf.count_nl(buf, a, b), segs, threads)
    n = int(sum(counts))
    row0 = np.concatenate(([0], np.cumsum(counts)))
    outi = np.zeros((max(ni, 1), n), np.int64)
    outf = np.zeros((max(nf, 1), n), np.float64)
    bad = np.zeros(n, np.uint8)
    change = np.zeros(n, np.uint8)
    ls = np.empty(n, np.int64)
    args = (cols, kinds, slots, outi, outf, bad, change, ls, P10_EXACT, P10_HI,
            P10_LO)
    _map(
        lambda i: _parse_stream(buf, segs[i][0], segs[i][1], int(row0[i]), *args
                               ), [(i,) for i in range(len(segs))], threads)
    change[0] = 1
    for i in range(1, len(segs)):
        r = int(row0[i])
        if r < n:
            la, lb = _first_token_len(buf,
                                      ls[r - 1]), _first_token_len(buf, ls[r])
            change[r] = 1 if la != lb or not np.array_equal(
                buf[ls[r - 1]:ls[r - 1] + la], buf[ls[r]:ls[r] + lb]) else 0
    if (bad == 2).any():
        raise bgzf.ParseError(
            "comment, empty line or line with too few columns")
    for r in np.flatnonzero(bad == 1):
        a = int(ls[r])
        e = a
        while buf[e] != 10:
            e += 1
        toks = bytes(buf[a:e]).split(b"\t")
        for j, (c, k) in enumerate(spec):
            try:
                s = toks[c].decode()
                if k == "i":
                    if not s.lstrip("-").isdigit():
                        raise ValueError(s)
                    outi[slots[j], r] = int(s)
                else:
                    if s.strip() != s or s == "" or s.lower().lstrip("+-") in (
                            "nan", "inf", "infinity"):
                        raise ValueError(s)
                    outf[slots[j], r] = float(s)
            except ValueError as ex:
                raise bgzf.ParseError(f"cannot parse {ex}") from None
    starts = np.flatnonzero(change)
    names = [
        bytes(buf[ls[r]:ls[r] + _first_token_len(buf, ls[r])]).decode()
        for r in starts
    ]
    ints = [outi[slots[j]] for j, k in enumerate(kinds) if k == 1]
    floats = [outf[slots[j]] for j, k in enumerate(kinds) if k == 2]
    return ints, floats, (starts, names)


@jit
def _digits15(x, p10hi, p10lo, digs):
    """Correctly rounded 15 significant digits of finite x != 0 (positive).

    Fills digs[0:15], returns the decimal exponent of the first digit, or -9999
    when the rounding is too close to call (caller uses Python).
    """
    # exponent estimate
    if x < 1e-280 or x > 1e280:
        return -9999
    e = int(np.floor(np.log10(x)))
    for _ in range(3):
        k = 14 - e
        if k < -300 or k > 300:
            return -9999
        th = p10hi[k + 345]
        tl = p10lo[k + 345]
        c = 134217729.0 * x
        ah = c - (c - x)
        al = x - ah
        c = 134217729.0 * th
        bh = c - (c - th)
        bl = th - bh
        pr = x * th
        er = ((ah * bh - pr) + ah * bl + al * bh) + al * bl
        er += x * tl
        hi = pr + er
        lo = er - (hi - pr)
        if hi >= 1e15:
            e += 1
            continue
        if hi < 1e14:
            e -= 1
            continue
        fl = np.floor(hi)
        frac = (hi - fl) + lo
        # fl is exact (hi < 2^53 * ... hi < 1e15 < 2^50 so hi - fl exact)
        if frac < 0.0:
            fl -= 1.0
            frac += 1.0
        elif frac >= 1.0:
            fl += 1.0
            frac -= 1.0
        if abs(frac - 0.5) < 1e-6:
            return -9999
        d = np.int64(fl)
        if frac > 0.5:
            d += 1
        if d >= 1000000000000000:
            if d == 1000000000000000 and frac > 0.5:
                # rounded up to 1e15: 15 digits 100000000000000 with exponent +
                # 1
                d = 100000000000000
                e += 1
            else:
                e += 1
                continue
        for i in range(14, -1, -1):
            digs[i] = d % 10
            d //= 10
        return e
    return -9999


@jit
def put_int(out, o, v):
    """Writes the integer `v` as text at out[o:]; returns the new offset."""
    if v == 0:
        out[o] = 48
        return o + 1
    if v < 0:
        out[o] = 45
        o += 1
        v = -v
    n = 0
    t = v
    while t > 0:
        n += 1
        t //= 10
    for i in range(n - 1, -1, -1):
        out[o + i] = 48 + v % 10
        v //= 10
    return o + n


@jit
def _put_double(out, o, x, p10hi, p10lo, digs):
    """Writes x as fwrite does; returns the new offset, or -1 for Python."""
    if x != x or x == np.inf or x == -np.inf:
        return -1
    if x == 0.0:
        if np.signbit(x):
            out[o] = 45
            o += 1
        out[o] = 48
        return o + 1
    neg = x < 0
    ax = -x if neg else x
    ex = _digits15(ax, p10hi, p10lo, digs)
    if ex == -9999:
        return -1
    nd = 15
    while nd > 1 and digs[nd - 1] == 0:
        nd -= 1
    sci_len = (1 if nd == 1 else nd +
               1) + 2 + (3 if (ex >= 100 or ex <= -100) else 2)
    if ex >= 0:
        fix_len = (ex + 1) if nd <= ex + 1 else nd + 1
    else:
        fix_len = 2 - ex - 1 + nd
    if neg:
        out[o] = 45
        o += 1
    # fmt_double: values in these ranges take the fwrite rule (shorter text),
    #   all others the printf %g rule
    cand = (ax >= 1e-4 and ax < 1e-3) or (ax >= 1e5)
    if cand:
        use_fixed = fix_len <= sci_len
    else:
        use_fixed = ex >= -4 and ex < 15
    if use_fixed:
        if ex >= 0:
            if nd <= ex + 1:
                for i in range(nd):
                    out[o] = 48 + digs[i]
                    o += 1
                for i in range(ex + 1 - nd):
                    out[o] = 48
                    o += 1
            else:
                for i in range(ex + 1):
                    out[o] = 48 + digs[i]
                    o += 1
                out[o] = 46
                o += 1
                for i in range(ex + 1, nd):
                    out[o] = 48 + digs[i]
                    o += 1
        else:
            out[o] = 48
            out[o + 1] = 46
            o += 2
            for i in range(-ex - 1):
                out[o] = 48
                o += 1
            for i in range(nd):
                out[o] = 48 + digs[i]
                o += 1
    else:
        out[o] = 48 + digs[0]
        o += 1
        if nd > 1:
            out[o] = 46
            o += 1
            for i in range(1, nd):
                out[o] = 48 + digs[i]
                o += 1
        out[o] = 101
        o += 1
        if ex < 0:
            out[o] = 45
            ae = -ex
        else:
            out[o] = 43
            ae = ex
        o += 1
        if ae < 10:
            out[o] = 48
            o += 1
        o = put_int(out, o, ae)
    return o


@jit
def _format_rows(chrom_bytes, chrom_len, ints, floats, out, rowend, p10hi,
                 p10lo, badrow):
    """Writes rows `chrom<TAB>ints<TAB>floats` as fwrite text.

    ints is (ni, n) int64, written as integers; floats is (nf, n).

    Returns:
        The bytes used, or -1 - r when row r needs Python (badrow[r] is set; the
        caller formats that row).
    """
    n = ints.shape[1]
    digs = np.zeros(16, np.int64)
    o = 0
    for r in range(n):
        start = o
        for k in range(chrom_len):
            out[o] = chrom_bytes[k]
            o += 1
        for c in range(ints.shape[0]):
            out[o] = 9
            o += 1
            o = put_int(out, o, ints[c, r])
        for c in range(floats.shape[0]):
            out[o] = 9
            o += 1
            o2 = _put_double(out, o, floats[c, r], p10hi, p10lo, digs)
            if o2 < 0:
                badrow[r] = 1
                o = start
                break
            o = o2
        if badrow[r] == 0:
            out[o] = 10
            o += 1
        rowend[r] = o
    return o


def format_rows(chrom: str, ints: np.ndarray, floats: np.ndarray) -> bytes:
    """Returns the text of rows `chrom<TAB>ints...<TAB>floats...` as bytes.

    Args:
        chrom: Contig name of all rows.
        ints: (ni, n) int64 columns, written as integers.
        floats: (nf, n) float64 columns, written as fwrite writes doubles.
    """
    from dmrdetector.fileio import bed
    n = ints.shape[1]
    if n == 0:
        return b""
    cb = np.frombuffer(chrom.encode(), np.uint8)
    out = np.empty(
        n * (len(cb) + 1 + 12 * ints.shape[0] + 26 * floats.shape[0] + 2),
        np.uint8)
    rowend = np.zeros(n, np.int64)
    badrow = np.zeros(n, np.uint8)
    used = _format_rows(cb, len(cb), np.ascontiguousarray(ints),
                        np.ascontiguousarray(floats), out, rowend, P10_HI,
                        P10_LO, badrow)
    if not badrow.any():
        return out[:used].tobytes()
    # rare rows (exact ties, zero-adjacent, non-finite): format with the Python
    # formatter, keep order
    rs = np.concatenate(([0], rowend[:-1]))
    pieces = []
    prev = 0
    for r in np.flatnonzero(badrow):
        pieces.append(out[rs[prev]:rs[r]].tobytes() if r > prev else b"")
        cols = [chrom] + [str(int(v)) for v in ints[:, r]] + [
            bed.fmt_doubles(floats[c, r:r + 1])[0]
            for c in range(floats.shape[0])
        ]
        pieces.append(("\t".join(cols) + "\n").encode())
        prev = r + 1
    pieces.append(out[rs[prev]:used].tobytes() if prev < n else b"")
    return b"".join(pieces)
