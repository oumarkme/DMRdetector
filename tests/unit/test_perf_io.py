"""Tests of the fast paths against the plain Python, pandas and tabix paths.

The fast paths are the numba readers and formatter, in-memory BGZF reads and
threaded writers. Compiled and interpreted kernels must agree on corner cases.
"""
import gzip
import os
from pathlib import Path
import time

import numpy as np
import pandas as pd
import pysam
import pytest

from dmrdetector.fileio import bed as dio
from dmrdetector.fileio import bgzf as bz
from dmrdetector.fileio import numparse
from dmrdetector.steps import atlas  # noqa: E402
from dmrdetector.steps import call_dmr  # noqa: E402
from dmrdetector.steps import target  # noqa: E402
from dmrdetector.util import rstats

from helpers import parse_tbi
from helpers import write_bgzf
from helpers import write_call_dmr_inputs


pytestmark = pytest.mark.skipif(
    bz.NO_JIT, reason="fast paths are off with DMRDETECTOR_NO_JIT=1")


def _parse_one(tokens, kind="f", use_py=False):
    """Parses column 1 of one-column lines with the stream kernel.

    This is the kernel `parse_buffer` uses.

    Returns:
        (values, flags).
    """
    text = ("\n".join("x\t" + t for t in tokens) + "\n").encode()
    buf = np.frombuffer(text, np.uint8).copy()
    n = len(tokens)
    outi = np.zeros((1, n), np.int64)
    outf = np.zeros((1, n), np.float64)
    bad = np.zeros(n, np.uint8)
    change = np.zeros(n, np.uint8)
    ls = np.empty(n, np.int64)
    cols = np.array([1], np.int64)
    kinds = np.array([1 if kind == "i" else 2], np.int64)
    slots = np.zeros(1, np.int64)
    fn = numparse._parse_stream.py_func if use_py else numparse._parse_stream
    assert fn(buf, 0, buf.size, 0, cols, kinds, slots, outi, outf, bad, change,
              ls, numparse.P10_EXACT, numparse.P10_HI, numparse.P10_LO) == n
    return (outi[0] if kind == "i" else outf[0]), bad


TOKENS = [
    "0", "0.5", "1", "-0.25", "1e-5", "1E3", "00012.5", "0.000352488493104247",
    "1.10346917362394e-05", "123456789012345678", "0.1234567890123456789",
    "9007199254740993", "4.94065645841247e-324", "1.7976931348623157e308",
    "1e400", ".5", "5.", "-0", "+3", "1e", "abc", "", "1.2.3", "NA", "nan",
    "1234567890123456789012", "9.999999999999999e22", "1.5e-30",
    "2.2250738585072014e-308"
]


def test_parse_float_equals_python_and_interpreted():
    rng = np.random.default_rng(0)
    toks = TOKENS + [
        "%.15g" % v
        for v in rng.random(3000) * 10.0**rng.integers(-40, 40, 3000)
    ]
    toks += ["%.17g" % v for v in rng.random(2000)
            ] + ["%.12f" % v for v in rng.random(1000)]
    got, bad = _parse_one(toks)
    for t, v, b in zip(toks, got, bad):
        if b:
            continue  # handed to the Python fallback by parse_buffer
        assert v == float(t), t
    assert bad[TOKENS.index("abc")] == 1 and bad[TOKENS.index("")] == 1 and bad[
        TOKENS.index("1e")] == 1
    assert not bad[:len(TOKENS)][[
        TOKENS.index(t)
        for t in ("0", "0.5", "1e-5", "00012.5", ".5", "5.", "-0")
    ]].any()
    small = toks[:400]
    a, ba = _parse_one(small)
    p, bp = _parse_one(small, use_py=True)
    assert (a.view(np.int64) == p.view(np.int64)).all() and (ba == bp).all()


def test_parse_int_and_flags():
    got, bad = _parse_one(
        ["0", "17", "-3", "12.5", "x", "99999999999999999999"], "i")
    assert list(got[:3]) == [0, 17, -3] and list(bad) == [0, 0, 0, 1, 1, 1]
    got, bad = _parse_one(["7"], "i", use_py=True)
    assert got[0] == 7 and bad[0] == 0


def test_format_equals_fmt_double_and_interpreted():
    rng = np.random.default_rng(1)
    x = np.concatenate([
        rng.random(20000),
        rng.random(20000) * 10.0**rng.integers(-30, 30, 20000),
        np.round(rng.random(5000), 15),
        rng.integers(0, 2000, 5000) / rng.integers(1, 2000, 5000),
        np.array([
            0.0, -0.0, 1.0, 0.5, 100000.0, 123456.0, 1e-4, 1e-5, 1.5e-4, 1e15,
            1e14, 99999.99999999999, 0.1, 1 / 3, 2.0**-22,
            2.384185791015625e-07, 5e-324, 1.7976931348623157e308,
            999999999999999.5, 0.99999999999999989, 0.0000999999999999999,
            np.inf, -np.inf, np.nan, 123456789012345678.0
        ])
    ])
    x = np.concatenate([x, -x[:500]])
    txt = numparse.format_rows("chr1", np.zeros((1, x.size), np.int64),
                               x[None, :]).decode().split("\n")[:-1]
    assert [t.split("\t")[2] for t in txt
           ] == [("" if False else s) for s in dio.fmt_doubles(x)]
    # interpreted kernel on a small sample (slow) gives the same bytes as the
    # compiled one
    xs = x[:300]
    ints = np.zeros((1, xs.size), np.int64)
    out1 = np.empty(xs.size * 60, np.uint8)
    out2 = np.empty(xs.size * 60, np.uint8)
    cb = np.frombuffer(b"chr1", np.uint8)
    args = (cb, 4, ints, xs[None, :].copy())
    r1, r2 = np.zeros(xs.size, np.int64), np.zeros(xs.size, np.int64)
    b1, b2 = np.zeros(xs.size, np.uint8), np.zeros(xs.size, np.uint8)
    u1 = numparse._format_rows(*args, out1, r1, numparse.P10_HI,
                               numparse.P10_LO, b1)
    u2 = numparse._format_rows.py_func(*args, out2, r2, numparse.P10_HI,
                                       numparse.P10_LO, b2)
    assert u1 == u2 and (out1[:u1] == out2[:u2]).all() and (b1 == b2).all()


def test_bgzf_roundtrip_and_pysam(tmp_path):
    rng = np.random.default_rng(2)
    text = "".join(
        f"chr1\t{i}\t{i + 1}\t{rng.integers(5, 50)}\t{rng.integers(0, 5)}\n"
        for i in range(40000)).encode()
    p = str(tmp_path / "x.bed.gz")
    with bz.BgzfWriter(p) as w:
        for i in range(0, len(text), 100000):
            w.write_blocks(bz.compress_bytes(text[i:i + 100000]))
    assert gzip.open(p).read() == text
    buf = bz.inflate(p, 3)
    assert buf.tobytes() == text
    dio.tabix_bed(p)
    with pysam.TabixFile(p) as tbx:
        assert len(list(tbx.fetch("chr1", 100, 200))) == 100
    ints, floats, (rs, names) = numparse.read_table(p, [(1, "i"), (3, "f")], 3)
    assert ints[0][5] == 5 and floats[0].size == 40000 and names == ["chr1"]


def test_read_table_corner_cases(tmp_path):

    def w(name, data):
        p = str(tmp_path / name)
        write_bgzf(p, data)
        return p

    ints, floats, runs = numparse.read_table(w("a.gz", b""), [(1, "i")], 1)
    assert ints[0].size == 0 and runs[1] == []
    ints, floats, runs = numparse.read_table(
        w("b.gz", b"chr1\t5\t6\nchr2\t7\t8"), [(1, "i")], 1)  # no final newline
    assert list(ints[0]) == [5, 7] and runs[1] == ["chr1", "chr2"]
    with pytest.raises(bz.ParseError):
        numparse.read_table(w("c.gz", b"#h\nchr1\t5\t6\n"), [(1, "i")], 1)
    with pytest.raises(bz.ParseError):
        numparse.read_table(w("d.gz", b"chr1\t5\t6\nchr1\tNA\t6\n"), [(1, "i")],
                            1)
    with pytest.raises(bz.ParseError):
        numparse.read_table(w("e.gz", b"chr1\t5\n"), [(2, "i")], 1)
    plain = tmp_path / "f.bed"
    plain.write_bytes(b"chr1\t5\t6\n")
    assert numparse.read_table(str(plain), [(1, "i")], 1) is None  # not BGZF
    gz = tmp_path / "g.gz"
    with gzip.open(gz, "wb") as fh:
        fh.write(b"chr1\t5\t6\n")
    assert numparse.read_table(str(gz), [(1, "i")],
                               1) is None  # plain gzip, not BGZF
    ints, _, _ = numparse.read_table(w("h.gz", b"chr1\t1e3\t6\n"), [(2, "i")],
                                     1)  # fallback parse of an int token
    assert ints[0][0] == 6


# ---- atlas ----------------------------------------------------------------


def _bed5(path, rows):
    rows = sorted(rows, key=lambda r: (r[0], r[1]))
    dio.write_bed_bgzip_tabix(
        str(path), "".join("\t".join(map(str, r)) + "\n" for r in rows))
    return str(path)


def _same_content(a, b):
    return gzip.open(a).read() == gzip.open(b).read()


@pytest.fixture(scope="module")
def beds(tmp_path_factory):
    d = tmp_path_factory.mktemp("beds")
    rng = np.random.default_rng(5)
    out = []
    for k in range(5):
        rows = []
        for chrom in ("chr1", "chr10", "chr2"):
            pos = np.unique(rng.integers(100, 300000, 4000))
            for p in pos:
                if rng.random() < 0.85:
                    c = int(rng.integers(6, 60))
                    rows.append((chrom, int(p), int(p) + 1, c,
                                 int(rng.integers(0, c + 1))))
        # duplicate start with a different END, and a duplicated row
        rows += [("chr1", 5000, 5001, 20, 10), ("chr1", 5000, 5002, 20, 5),
                 ("chr1", 5000, 5001, 30, 3)]
        out.append(_bed5(d / f"S{k}_autosomal.bed.gz", rows))
    return out


def test_build_fast_equals_legacy(beds, tmp_path):
    a, b = str(tmp_path / "a.gz"), str(tmp_path / "b.gz")
    ta = atlas.build_atlas(beds,
                           a,
                           3,
                           0.1,
                           threads=3,
                           window=50_000,
                           log=lambda m: None)
    import time as _t
    tb = atlas._build_atlas_tabix(sorted(beds), b, 3, 0.1, 2,
                                  str(tmp_path / "s"), 50_000, lambda m: None,
                                  _t.time())
    assert ta == tb and _same_content(a, b)
    # thread count does not change the output
    c = str(tmp_path / "c.gz")
    atlas.build_atlas(beds,
                      c,
                      3,
                      0.1,
                      threads=1,
                      window=10_000_000,
                      log=lambda m: None)
    assert _same_content(a, c)


def test_smooth_fast_equals_legacy_and_index(beds, tmp_path):
    raw = str(tmp_path / "raw.gz")
    atlas.build_atlas(beds, raw, 3, 0.1, threads=2, log=lambda m: None)
    a, b = str(tmp_path / "a.gz"), str(tmp_path / "b.gz")
    atlas.smooth_atlas(raw, 40, a, True, threads=3, log=lambda m: None)
    atlas._smooth_atlas_pool(raw, 40, b, True, 2, lambda m: None, time.time())
    assert _same_content(a, b)
    assert os.path.exists(a + ".tbi")
    # the index written in the same step equals the one htslib builds from the
    # finished file (same content)
    import shutil
    shutil.copy(a, str(tmp_path / "c.gz"))
    pysam.tabix_index(str(tmp_path / "c.gz"),
                      preset="bed",
                      force=True,
                      keep_original=True)
    assert parse_tbi(a + ".tbi") == parse_tbi(str(tmp_path / "c.gz") + ".tbi")
    with pysam.TabixFile(a) as tbx:
        assert sorted(tbx.contigs) == ["chr1", "chr10", "chr2"]
        assert len(list(tbx.fetch("chr2", 100000, 200000))) > 0
    atlas.smooth_atlas(raw,
                       40,
                       str(tmp_path / "n.gz"),
                       True,
                       threads=1,
                       log=lambda m: None,
                       index=False)
    assert not os.path.exists(str(tmp_path / "n.gz") + ".tbi")
    # shuffled / plain-gzip raw atlas takes the slow path and gives the same
    # text
    txt = gzip.open(raw).read().decode().splitlines()
    rng = np.random.default_rng(0)
    rng.shuffle(txt)
    sh = str(tmp_path / "shuf.gz")
    with gzip.open(sh, "wt") as fh:
        fh.write("\n".join(txt) + "\n")
    s2 = str(tmp_path / "s2.gz")
    atlas.smooth_atlas(sh, 40, s2, True, threads=2, log=lambda m: None)
    assert _same_content(a, s2)


def test_smooth_chromosome_matches_scipy_reference():
    from scipy import fft as sfft
    rng = np.random.default_rng(3)
    start = np.unique(rng.integers(1000, 90000, 1500))
    end = start + 1
    cov = rng.integers(5, 60, start.size).astype(float)
    raw = rng.random(start.size)
    phi = rng.random(start.size) * 0.2
    k = atlas.gaussian_kernel(60)
    ref, var = atlas.smooth_chromosome(start, end, cov, raw, phi, k)
    span = int(end.max() - start.min() + 1)
    idx = start - start.min()
    n = rstats.nextn(span + k.size - 1)
    kf = sfft.rfft(k, n=n)

    def cmul(a, b):
        # every product and sum rounded on its own, as the program does; the
        # `*` of numpy fuses them on a CPU with AVX2 and FMA3
        out = np.empty(a.shape, np.complex128)
        out.real = a.real * b.real - a.imag * b.imag
        out.imag = a.real * b.imag + a.imag * b.real
        return out

    def conv(vals):
        v = np.zeros(span)
        v[idx] = vals
        return sfft.irfft(cmul(sfft.rfft(v, n=n), kf),
                          n=n)[k.size // 2:k.size // 2 + span][idx]

    sm, sc = conv(raw * cov), conv(cov)
    exp = np.zeros(idx.size)
    ok = sc > 1e-9
    exp[ok] = sm[ok] / sc[ok]
    assert (ref == exp).all() and (var == exp *
                                   (1 - exp) * phi).all()  # bitwise
    r3, v3 = atlas.smooth_chromosome(start, end, cov, raw, phi, k, inner=3)
    assert (r3 == ref).all() and (v3 == var).all()


# ---- call-dmr ---------------------------------------------------------------


def _raw_text(df, path):
    """Returns the text that `call_dmr.write_raw` makes of a raw DMR table."""
    call_dmr.write_raw(df, str(path))
    return Path(path).read_text()


# The pandas reader (the route of NO_JIT, of unusual files and of a one-thread
# read of an indexed file) parses decimals with float_precision="round_trip",
# as correctly rounded as the byte-level reader. The two tables are therefore
# equal in every column, bit for bit and as text.


def _assert_same_dmr_table(fast, slow, tmp_path, tag):
    assert _raw_text(fast, tmp_path / f"{tag}_f.tsv") == _raw_text(
        slow, tmp_path / f"{tag}_s.tsv")
    pd.testing.assert_frame_equal(fast, slow, check_exact=True)


def test_pandas_round_trip_parse_equals_python_float_and_numparse():
    """Checks the parser that the pandas fallbacks must use.

    Tokens of 15 to 17 significant digits, like those of the atlas and of the
    raw DMR tables, are read by pandas with float_precision="round_trip", by
    Python float() and by the byte-level kernel. The default pandas parser can
    lose the last digits of such a token (it reads 0.000535864718615665 as
    0.0005358647186156), which is why the fallbacks do not use it.
    """
    from io import StringIO
    rng = np.random.default_rng(11)
    toks = []
    for _ in range(30000):
        nd = int(rng.integers(15, 18))
        d = "".join(map(str, rng.integers(0, 10, nd)))
        toks.append(f"0.{d}" if rng.random() < 0.5 else
                    f"{d[0]}.{d[1:]}e{int(rng.integers(-30, 3))}")
    toks += ["0.000535864718615665", "0.000111476322925955"]
    want = np.array([float(t) for t in toks])
    got = pd.read_csv(StringIO("\n".join(toks)),
                      header=None,
                      float_precision="round_trip")[0].to_numpy()
    assert (got == want).all()
    vals, flags = _parse_one(toks)
    ok = flags == 0  # a rare token is handed over to float()
    assert ok.mean() > 0.99 and (vals[ok] == want[ok]).all()


def test_call_dmr_fast_equals_pandas(tmp_path, monkeypatch):
    d = str(tmp_path / "syn")
    write_call_dmr_inputs(d)
    for scope in ("autosomal", "male_X"):
        args = (f"{d}/{scope}.sample.bed.gz", f"{d}/{scope}.atlas.bed.gz",
                f"{d}/{scope}.cpg.bed",
                f"{d}/icr.bed" if scope == "autosomal" else None)
        fast = call_dmr.call_dmr(*args,
                                 scope=scope,
                                 log=lambda m: None,
                                 threads=3)
        with monkeypatch.context() as mp:
            mp.setattr(call_dmr, "_read_fast", lambda *a, **k: None)
            slow = call_dmr.call_dmr(*args, scope=scope, log=lambda m: None)
        _assert_same_dmr_table(fast, slow, tmp_path, scope)
        assert len(fast) > 0
        # atlas copy without .tbi: same result (the legacy path would read it
        # with pandas)
        import shutil
        noidx = str(tmp_path / f"{scope}.noidx.bed.gz")
        shutil.copy(args[1], noidx)
        again = call_dmr.call_dmr(args[0],
                                  noidx,
                                  *args[2:],
                                  scope=scope,
                                  log=lambda m: None,
                                  threads=2)
        pd.testing.assert_frame_equal(fast, again)
        dedupe = call_dmr.call_dmr(*args,
                                   scope=scope,
                                   log=lambda m: None,
                                   tabix_duplicates=False)
        with monkeypatch.context() as mp:
            mp.setattr(call_dmr, "_read_fast", lambda *a, **k: None)
            dedupe2 = call_dmr.call_dmr(*args,
                                        scope=scope,
                                        log=lambda m: None,
                                        tabix_duplicates=False)
        _assert_same_dmr_table(dedupe, dedupe2, tmp_path, f"{scope}_dedupe")


# ---- target
# -------------------------------------------------------------------


def _tgt_inputs(tmp_path):
    from test_target import _build
    a, s, tb, pos = _build(tmp_path, sample_shift=-0.45)
    return a, s, tb, pos


def test_target_unindexed_inputs_and_het_variants(tmp_path):
    a, s, tb, pos = _tgt_inputs(tmp_path)
    hap = [("chr1", int(p), int(p) + 1, 5, 1) for p in pos]
    h1 = _bed5(tmp_path / "S1.hap1.bed.gz", hap)
    vcf = tmp_path / "s.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="g">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        f"chr1\t{int(pos[3]) + 1}\t.\tA\tG\t50\tPASS\t.\tGT\t0|1\n"
        f"chr1\t{int(pos[9]) + 1}\t.\tA\tG\t50\tPASS\t.\tGT\t1|1\n"
        f"chr1\t{int(pos[12])}\t.\tC\tT\t50\tPASS\t.\tGT\t0/1\n"
        f"chr1\t{int(pos[20]) + 1}\t.\tAC\tA\t50\tPASS\t.\tGT\t0/1\n")
    gz = str(tmp_path / "s.vcf.gz")
    pysam.tabix_compress(str(vcf), gz, force=True)
    pysam.tabix_index(gz, preset="vcf", force=True)
    from dmrdetector.fileio.vcf import het_snp_sites

    def smp(**k):
        return [{
            "sample_id": "S1",
            "type": "CASE",
            "merged": s,
            "hap1": h1,
            **k
        }]

    base = target.target_report(smp(vcf=str(vcf)), a, tb,
                                log=lambda m: None)  # full scan (no index)
    r_idx = target.target_report(smp(vcf=gz), a, tb,
                                 log=lambda m: None)  # region fetch (indexed)
    r_het = target.target_report(smp(het=het_snp_sites(str(vcf))),
                                 a,
                                 tb,
                                 log=lambda m: None)  # precomputed sites
    pd.testing.assert_frame_equal(base, r_idx)
    pd.testing.assert_frame_equal(base, r_het)
    assert base[(base.REGION_NAME == "R1") &
                (base.HAPLOTYPE == "merged")].iloc[0].NUM_CPG == 39
    # sample bed and atlas without index: whole-file readers
    for f in (a, s, h1):
        os.remove(f + ".tbi")
    noidx = target.target_report(smp(vcf=gz), a, tb, log=lambda m: None)
    pd.testing.assert_frame_equal(base, noidx)


def test_merge_order_equals_sort_grouping():
    rng = np.random.default_rng(9)
    for trial in range(30):
        nf = int(rng.integers(1, 6))
        ss, ee, fo = [], [], [0]
        for f in range(nf):
            n = int(rng.integers(0, 60))
            s = np.sort(rng.integers(0, 40, n))
            e = s + rng.choice([1, 1, 1, 2, 3], n)
            ss.append(s)
            ee.append(e)
            fo.append(fo[-1] + n)
        s, e = np.concatenate(ss), np.concatenate(ee)
        if s.size == 0:
            continue
        cov = rng.integers(5, 50, s.size).astype(float)
        nm = np.floor(cov * rng.random(s.size))
        ref = atlas.aggregate_sites(s, e, cov, nm)
        got = atlas.aggregate_sites(s, e, cov, nm, fo)
        for x, y in zip(ref, got):
            assert (x == y).all()
        order = np.empty(s.size, np.int64)
        gb = np.empty(s.size + 1, np.int64)
        order2 = order.copy()
        gb2 = gb.copy()
        g1 = atlas._merge_order(s, e, np.array(fo, np.int64), order, gb)
        g2 = atlas._merge_order.py_func(s, e, np.array(fo, np.int64), order2,
                                        gb2)
        assert g1 == g2 and (order == order2).all() and (gb[:g1 + 1]
                                                         == gb2[:g2 + 1]).all()


@pytest.fixture(scope="module")
def wide_beds(tmp_path_factory):
    """Returns 5 BEDs whose sites lie in several 10 Mbp cells.

    The cells are the grid of the output blocks. The sites are on both sides of
    a cell boundary, with equal START and different END, plus one record much
    wider than a tabix window.
    """
    d = tmp_path_factory.mktemp("wide")
    rng = np.random.default_rng(11)
    out = []
    for k in range(5):
        rows = []
        for chrom in ("chr1", "chr10", "chr2"):
            for lo, hi, n in ((100, 300_000, 3000), (9_960_000, 10_040_000,
                                                     1500), (29_990_000,
                                                             30_010_000, 600)):
                for p in np.unique(rng.integers(lo, hi, n)):
                    if rng.random() < 0.85:
                        c = int(rng.integers(6, 60))
                        rows.append((chrom, int(p), int(p) + 1, c,
                                     int(rng.integers(0, c + 1))))
        rows += [("chr1", 5000, 5001, 20, 10), ("chr1", 5000, 5002, 20, 5),
                 ("chr1", 10_000_000, 10_000_001, 25, 5),
                 ("chr1", 9_999_999, 10_200_000, 30, 3)
                ]  # a 200 kb record starting before the cell boundary
        out.append(_bed5(d / f"S{k}_autosomal.bed.gz", rows))
    return out


def _md5(path):
    import hashlib
    return hashlib.md5(open(path, "rb").read()).hexdigest()


def test_raw_atlas_bytes_do_not_depend_on_window_or_budget(wide_beds, tmp_path):
    """Checks the raw atlas bytes against window size, budget and threads.

    Three window sizes (one that does not divide the 10 Mbp grid), three memory
    budgets and 1 or 3 threads give a byte-identical compressed file, and its
    text equals the tabix path.
    """
    ref = str(tmp_path / "ref.gz")
    t_ref = atlas.build_atlas(wide_beds,
                              ref,
                              3,
                              0.1,
                              threads=1,
                              window=10_000_000,
                              log=lambda m: None)
    want = _md5(ref)
    for i, kw in enumerate(
        (dict(window=50_000, threads=3), dict(window=7_000_000, threads=2),
         dict(window=156_250, threads=1), dict(mem_gb=0.7, threads=3),
         dict(mem_gb=2, threads=2), dict(mem_gb=500, threads=4))):
        o = str(tmp_path / f"o{i}.gz")
        assert atlas.build_atlas(wide_beds, o, 3, 0.1, log=lambda m: None,
                                 **kw) == t_ref
        assert _md5(o) == want, kw
    leg = str(tmp_path / "leg.gz")
    import time as _t
    atlas._build_atlas_tabix(sorted(wide_beds), leg, 3, 0.1, 2,
                             str(tmp_path / "s"), 10_000_000, lambda m: None,
                             _t.time())
    assert _same_content(ref, leg)
    cells = {(ln.split(b"\t")[0], int(ln.split(b"\t")[1]) // 10_000_000)
             for ln in gzip.open(ref).read().splitlines()}
    assert len(cells) >= 6  # the output spans several cells of the grid


def test_indexed_window_reader_matches_tabix(tmp_path):
    """Checks _IndexedBed.rows against the rows with START in [w0, w1).

    The file has many BGZF blocks, wide records and a second contig; the windows
    start and end inside, between and beyond the data.
    """
    rng = np.random.default_rng(3)
    rows = []
    for chrom, n in (("chr1", 150_000), ("chr2", 20_000)):
        pos = np.unique(rng.integers(1, 20_000_000, n))
        for p in pos:
            rows.append(
                (chrom, int(p),
                 int(p) + (1 if rng.random() < 0.999 else 40_000), 7, 3))
    rows.sort(key=lambda r: (r[0], r[1]))
    path = _bed5(tmp_path / "w.bed.gz", rows)
    ib = atlas._IndexedBed(path)
    for chrom, w0, w1 in (("chr1", 0, 10_000_000), ("chr1", 10_000_000,
                                                    20_000_000),
                          ("chr1", 123_456, 123_457), ("chr1", 7_000_000,
                                                       7_050_000),
                          ("chr1", 19_990_000, 30_000_000), ("chr2", 0, 78_125),
                          ("chr2", 5_000_000, 5_078_125),
                          ("chr1", 25_000_000, 26_000_000), ("chrX", 0, 10)):
        got = ib.rows(chrom, w0, w1)
        want = [r for r in rows if r[0] == chrom and w0 <= r[1] < w1]
        if not want:
            assert got is None, (chrom, w0, w1)
            continue
        assert got[0].tolist() == [r[1] for r in want] and got[1].tolist() == [
            r[2] for r in want
        ], (chrom, w0, w1)


def test_raw_atlas_plan_follows_the_budget():
    from dmrdetector import runtime
    size, ext = 2.51e9, runtime.SCOPE_EXTENT["autosomal"]  # 12 samples
    prev = None
    for mem in (16, 24, 48, 96, 288):
        w, k = runtime.plan_raw_atlas(size, ext, 3, 24, mem)
        floor, one = runtime.raw_atlas_floor_gb(size, ext,
                                                3), runtime.raw_atlas_window_gb(
                                                    size, ext, w)
        assert floor + k * one <= mem
        assert runtime.RAW_WINDOW_MIN <= w <= runtime.RAW_WINDOW_MAX
        assert 1 <= k <= 24
        assert prev is None or (
            w, k
        ) >= prev  # more memory never gives a smaller window or fewer workers
        prev = (w, k)
    assert runtime.plan_raw_atlas(size, ext, 3, 24,
                                  288) == (runtime.RAW_WINDOW_MAX, 24)
    for n in (12, 40, 130):  # the least memory does not grow with the cohort
        assert runtime.raw_atlas_min_gb(size * n / 12, ext, 3) < 6.5
    assert runtime.plan_raw_atlas(size, ext, 3, 4,
                                  runtime.raw_atlas_min_gb(size, ext, 3) -
                                  0.1) is None


def test_raw_atlas_stops_when_not_even_one_small_window_fits(
        wide_beds, tmp_path):
    with pytest.raises(
            SystemExit,
            match=
            r"raw atlas of 5 samples: needs at least .* GB .* more than the "
            r"0\.1 GB available"):
        atlas.build_atlas(wide_beds,
                          str(tmp_path / "o.gz"),
                          3,
                          0.1,
                          threads=2,
                          mem_gb=0.1,
                          log=lambda m: None)


def test_raw_atlas_rebuilds_an_index_older_than_its_data(wide_beds, tmp_path):
    """Checks that an index older than its data is rebuilt in scratch.

    The windows are read through the tabix index, so a stale index is not
    trusted: the input and its .tbi stay untouched and the output does not
    change.
    """
    import shutil
    d = tmp_path / "in"
    d.mkdir()
    files = []
    for f in wide_beds:
        g = str(d / os.path.basename(f))
        shutil.copy(f, g)
        shutil.copy(f + ".tbi", g + ".tbi")
        files.append(g)
    ref = str(tmp_path / "ref.gz")
    atlas.build_atlas(files, ref, 3, 0.1, threads=2, log=lambda m: None)
    # a stale index: the data are newer than the index and the index is wrong
    # (taken from another file)
    shutil.copy(files[1] + ".tbi", files[0] + ".tbi")
    os.utime(files[0] + ".tbi", (1, 1))
    before = open(files[0] + ".tbi", "rb").read()
    out = str(tmp_path / "out.gz")
    atlas.build_atlas(files,
                      out,
                      3,
                      0.1,
                      threads=2,
                      scratch=str(tmp_path / "scr"),
                      log=lambda m: None)
    assert _md5(out) == _md5(ref) and open(files[0] + ".tbi",
                                           "rb").read() == before
    # the new index is in a folder of scratch named by a hash of the path
    assert len(
        list((tmp_path / "scr").glob("*/" + os.path.basename(files[0]) +
                                     ".tbi"))) == 1


def test_raw_atlas_removes_the_scratch_folder_it_made(wide_beds, tmp_path,
                                                      monkeypatch):
    """Checks that a call without `scratch` leaves no folder in the temp dir."""
    import tempfile
    tmpd = tmp_path / "tmp"
    tmpd.mkdir()
    monkeypatch.setattr(tempfile, "tempdir", str(tmpd))
    atlas.build_atlas(list(wide_beds),
                      str(tmp_path / "out.gz"),
                      3,
                      0.1,
                      threads=2,
                      log=lambda m: None)
    assert list(tmpd.iterdir()) == []
    mine = tmp_path / "mine"
    mine.mkdir()
    atlas.build_atlas(list(wide_beds),
                      str(tmp_path / "out2.gz"),
                      3,
                      0.1,
                      threads=2,
                      scratch=str(mine),
                      log=lambda m: None)
    assert mine.is_dir()  # a folder given by the caller is theirs to remove
