import gzip
import os
import struct
import subprocess
import sys
import tempfile
import zlib

import numpy as np
import oracle_cpg_islands_v3_orig as oracle_cpg
import oracle_icr_v3_orig as oracle_icr
import pandas as pd
import pytest

from dmrdetector import settings
from dmrdetector.fileio import atlas_io
from dmrdetector.fileio import bed
from dmrdetector.steps import cpg_islands
from dmrdetector.steps import icr


def _atlas(rows):
    return pd.DataFrame(rows,
                        columns=[
                            "CHROM", "START", "END", "N", "COV", "PHI",
                            "RAW_METH", "REF_METH", "REF_VAR"
                        ])


def test_unread_path_names_the_file_only_while_it_is_not_loaded(tmp_path):
    path = tmp_path / "a.bed"
    path.write_text("chr1\t1\t2\t8\t10\t0.1\t0.2\t0.2\t0.0\n")
    lazy = atlas_io.read_atlas(str(path))
    assert atlas_io.unread_path(str(path)) == str(path)
    assert atlas_io.unread_path(path) == str(path)
    assert atlas_io.unread_path(lazy) == str(path)
    lazy.frame()
    assert atlas_io.unread_path(lazy) is None
    assert atlas_io.unread_path(_atlas([])) is None


def test_cpg_island_with_shore_extension():
    # 20 CORE sites 10 bp apart (190 bp, density ok), shore 100 bp right, far
    # site not used
    core = [("chr1", 1000 + 10 * i, 1001 + 10 * i, 8, 10, 0.01, 0.02, 0.02, 0.0)
            for i in range(20)]
    shore = [("chr1", 1290, 1291, 8, 10, 0.01, 0.2, 0.2, 0.0)
            ]  # DIST 100 to last core (1190)
    far = [("chr1", 5000, 5001, 8, 10, 0.01, 0.9, 0.9, 0.0)]
    res = cpg_islands.detect(_atlas(core + shore + far), log=lambda m: None)
    assert res.values.tolist() == [["chr1", 1000, 1190, 1000, 1290, 191, 20]]


def test_cpg_n_filter_is_strict():
    rows = [("chr1", 1000 + 10 * i, 1001 + 10 * i, 5, 10, 0.01, 0.02, 0.02, 0.0)
            for i in range(20)]
    assert cpg_islands.detect(_atlas(rows), log=lambda m: None).empty


def test_female_x_hemi():
    rows = [("chrX", 1000 + 10 * i, 1001 + 10 * i, 8, 10, 0.5, 0.5, 0.5, 0.0)
            for i in range(20)]
    res = cpg_islands.detect(_atlas(rows), "female_X", log=lambda m: None)
    assert res.values.tolist() == [["chrX", 1000, 1190, 1000, 1190, 191, 20]]


def test_icr_overlapping_candidates_count_twice_like_v1():
    cand = pd.DataFrame({
        "CHROM": ["chr1", "chr1"],
        "START": [100, 150],
        "END": [300, 400],
        "GENE": ["A", "B"]
    })
    rows = [("chr1", 160 + 10 * i, 161 + 10 * i, 6, 10, 0.01, 0.5, 0.5, 0.0)
            for i in range(3)]
    res = icr.detect(_atlas(rows), cand, log=lambda m: None)
    assert dict(zip(res["GENE"], res["NUM_CPG"])) == {"A": 6, "B": 6}
    res1 = icr.detect(_atlas(rows),
                      cand,
                      v1_duplicates=False,
                      log=lambda m: None)
    assert res1.empty  # 3 sites < min_sites


def _quiet(msg):
    return None


CH = [
    "chr1", "chr2", "chr10", "chr22", "chrX", "chrY", "chrM", "chrUn_gl000220"
]
GAPS = [1, 2, 7, 10, 149, 150, 151, 152, 400]
METH = [
    0.0, 0.05, 0.0999999999, 0.1, 0.2999999, 0.3, 0.30000001, 0.4, 0.5, 0.6,
    0.6000001, 0.9, np.nan
]
PHI = [0.0, 0.04999999, 0.05, 0.2, np.nan]
NS = [3, 5, 6, 12]


def bgzf_bytes(data: bytes, block: int = 65280) -> bytes:
    """Returns `data` as a BGZF container (standard header, one BC field).

    `block` is the number of uncompressed bytes per block.
    """
    out = []
    for i in range(0, len(data), block):
        chunk = data[i:i + block]
        co = zlib.compressobj(6, zlib.DEFLATED, -15)
        z = co.compress(chunk) + co.flush()
        bsize = len(z) + 25
        out.append(
            b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00" +
            struct.pack("<H", bsize) + z +
            struct.pack("<II", zlib.crc32(chunk), len(chunk)))
    out.append(
        bytes.fromhex(
            "1f8b08040000000000ff0600424302001b0003000000000000000000"))
    return b"".join(out)


def fmt(x, rng):
    if isinstance(x, float) and np.isnan(x):
        return "NA"
    k = rng.integers(0, 4)
    if k == 0:
        return repr(float(x))
    if k == 1:
        return f"{x:.15g}"
    if k == 2:
        return f"{x:.8e}"
    return f"{x:.17g}"


def make_atlas(rng, n_rows=400, chroms=CH):
    rows = []
    for ch in chroms:
        pos = int(rng.integers(1, 5000))
        for _ in range(int(rng.integers(0, max(1, n_rows // len(chroms) * 2)))):
            pos += int(rng.choice(
                GAPS)) if rng.random() < 0.97 else 0  # some duplicate STARTs
            end = pos + 1 if rng.random() < 0.9 else pos + int(
                rng.integers(0, 4))
            rows.append([
                ch, pos, end,
                int(rng.choice(NS)),
                int(rng.integers(5, 300)),
                float(rng.choice(PHI)),
                float(rng.random()),
                float(rng.choice(METH)),
                float(rng.random() * 0.01)
            ])
    order = rng.permutation(len(rows)) if rng.random() < 0.3 else np.arange(
        len(rows))
    rows = [rows[i] for i in order]
    lines = []
    for r in rows:
        lines.append("\t".join([
            r[0],
            str(r[1]),
            str(r[2]),
            str(r[3]),
            str(r[4]),
            fmt(r[5], rng),
            fmt(r[6], rng),
            fmt(r[7], rng),
            fmt(r[8], rng)
        ]))
    return "\n".join(lines) + ("\n" if lines else "")


def write_variants(d, text, block=300):
    """Writes plain, gzip and bgzf copies of one atlas text.

    The bgzf copy has small blocks, so lines straddle blocks.
    """
    raw = text.encode()
    paths = {
        "plain": os.path.join(d, "a.bed"),
        "gz": os.path.join(d, "a.bed.gz"),
        "bgzf": os.path.join(d, "b.bed.gz")
    }
    open(paths["plain"], "wb").write(raw)
    open(paths["gz"], "wb").write(gzip.compress(raw))
    open(paths["bgzf"], "wb").write(bgzf_bytes(raw, block))
    return paths


def _same(a, b):
    assert list(a.columns) == list(b.columns)
    assert len(a) == len(b)
    for c in a.columns:
        assert a[c].tolist() == b[c].tolist(), c


def teardown_module(module=None):
    atlas_io._CHUNK = 3 << 19


def test_cpg_matches_oracle_fuzz():
    rng = np.random.default_rng(1)
    with tempfile.TemporaryDirectory() as d:
        for seed in range(60):
            text = make_atlas(rng, n_rows=int(rng.choice([0, 20, 400, 1500])))
            atlas_io._CHUNK = int(rng.choice([150, 700, 3 << 19
                                             ]))  # byte range per reader task
            paths = write_variants(d,
                                   text,
                                   block=int(rng.choice([97, 300, 5000])))
            ref_atlas = oracle_cpg.read_atlas(paths["plain"]) if text else None
            for scope in ("autosomal", "female_X", "male_X", "Y"):
                for include_sex in (True, False):
                    kw = dict(include_sex=include_sex, log=_quiet)
                    if seed % 3 == 0:
                        kw.update(max_gap=int(rng.choice([149, 150, 151])),
                                  min_length=int(rng.choice([100, 150, 151])))
                    if text:
                        want = oracle_cpg.detect(ref_atlas, scope, **kw)
                    else:
                        continue
                    _same(cpg_islands.detect(ref_atlas, scope, **kw),
                          want)  # DataFrame input
                    for k in ("plain", "gz", "bgzf"):
                        _same(cpg_islands.detect(paths[k], scope, **kw),
                              want)  # fast reader
                        _same(
                            cpg_islands.detect(atlas_io.read_atlas(paths[k]),
                                               scope, **kw), want)  # LazyAtlas


def test_cpg_corner_gaps_and_edges():

    def isl(gap):
        # two dense blocks (20 CpG, 10 bp apart, 190 bp long) separated by `gap`
        # bp
        starts = [1000 + 10 * i for i in range(20)
                 ] + [1190 + gap + 10 * i for i in range(20)]
        rows = [("chr1", s, s + 1, 8, 10, 0.01, 0.02, 0.02, 0.0) for s in starts
               ]
        got = cpg_islands.detect(_atlas(rows), log=_quiet)
        _same(got, oracle_cpg.detect(_atlas(rows), log=_quiet))
        return got.values.tolist()

    # gap exactly 150 joins, 151 breaks
    assert len(isl(150)) == 1 and len(isl(151)) == 2 and len(isl(149)) == 1
    # shore at exactly 150 from the island end is included, at 151 not; window
    # edge s0 - 150
    core = [("chr1", 1000 + 10 * i, 1001 + 10 * i, 8, 10, 0.01, 0.02, 0.02, 0.0)
            for i in range(20)]
    for off, inside in ((150, True), (151, False)):
        sh = [("chr1", 1190 + off, 1191 + off, 8, 10, 0.01, 0.2, 0.2, 0.0),
              ("chr1", 1000 - off, 1001 - off, 8, 10, 0.01, 0.2, 0.2, 0.0)]
        r = cpg_islands.detect(_atlas(core + sh), log=_quiet).values.tolist()[0]
        assert (r[3], r[4]) == ((1000 - off, 1190 + off) if inside else
                                (1000, 1190))
        _same(cpg_islands.detect(_atlas(core + sh), log=_quiet),
              oracle_cpg.detect(_atlas(core + sh), log=_quiet))
    # LENGTH 149 / 150 / 151 (10 CpG are enough for the density rule)
    for last, keep in ((148, False), (149, True), (150, True)):
        rows = [("chr1", s, s + 1, 8, 10, 0.01, 0.02, 0.02, 0.0)
                for s in (0, 20, 40, 60, 80, 100, 120, 140, last, last)]
        rows = sorted(set(rows))
        got = cpg_islands.detect(_atlas(rows), log=_quiet)
        _same(got, oracle_cpg.detect(_atlas(rows), log=_quiet))
        assert (len(got) == 1) == keep
    # exact thresholds are not CORE (strict <), N == 5 dropped, NaN never CORE
    for meth, phi in ((0.1, 0.01), (0.01, 0.05), (np.nan, 0.01), (0.01,
                                                                  np.nan)):
        rows = [("chr1", 1000 + 10 * i, 1001 + 10 * i, 8, 10, phi, 0.02, meth,
                 0.0) for i in range(20)]
        assert cpg_islands.detect(_atlas(rows), log=_quiet).empty
    # empty results keep the header-only form
    e = cpg_islands.detect(_atlas([("chr1", 5, 6, 8, 10, 0.5, 0.9, 0.9, 0.0)]),
                           log=_quiet)
    assert e.empty and list(e.columns) == cpg_islands.OUT_COLS
    assert e.to_csv(sep="\t",
                    index=False) == oracle_cpg.detect(
                        _atlas([("chr1", 5, 6, 8, 10, 0.5, 0.9, 0.9, 0.0)]),
                        log=_quiet).to_csv(sep="\t", index=False)


def test_kernel_jit_equals_python():
    rng = np.random.default_rng(5)
    kern = cpg_islands._chrom_kernel
    py = getattr(kern, "py_func", kern)
    for _ in range(80):
        n = int(rng.integers(0, 300))
        start = np.cumsum(rng.choice(GAPS + [0], n)).astype(np.int64) + 1000
        meth = rng.choice(METH, n)
        phi = rng.choice(PHI, n)
        for female in (False, True):
            a = kern(start, meth, phi, female, 0.1, 0.05, 0.3, 0.4, 0.6, 150,
                     150)
            b = py(start, meth, phi, female, 0.1, 0.05, 0.3, 0.4, 0.6, 150, 150)
            assert a[0] == b[0] and a[6] == b[6] and a[7] == b[7]
            for x, y in zip(a[1:6], b[1:6]):
                assert x[:a[0]].tolist() == y[:b[0]].tolist()


_NO_JIT_SCRIPT = r"""
import sys, hashlib, numpy as np
sys.path.insert(0, sys.argv[1])
from dmrdetector.fileio import atlas_io
from dmrdetector.steps import cpg_islands
txt = open(sys.argv[2], "rb").read()
h = hashlib.md5()
for kind, mh, strict in ((0, 0.3, True), (1, 0.6, False)):
    z = np.zeros(1, np.int64)
    win = (np.array([0, 40], np.int64), np.array([5000, 60000], np.int64),
           np.arange(25).clip(0, 1) * 2)
    n = txt.count(b"\n")
    o = [np.empty(n, np.int8), np.empty(n, np.int64), np.empty(n, np.int64),
         np.empty(n), np.empty(n)]
    cnt = np.zeros(5, np.int64)
    slow = np.zeros((10000, 4), np.int64)
    k = atlas_io.parse_lines(
        np.frombuffer(txt, np.uint8), kind, 23 if kind == 0 else 21, 5.0,
        strict, -np.inf if kind == 0 else 0.4, mh, kind == 0, win[0], win[1],
        win[2], *o, slow, cnt)
    for x in o:
        h.update(np.nan_to_num(x[:k], nan=-7.5).tobytes())
    h.update(cnt.tobytes())
    h.update(slow.tobytes())
rng = np.random.default_rng(3)
n = 500
start = np.cumsum(rng.choice([1, 10, 149, 150, 151, 400], n)).astype(np.int64)
meth = rng.choice([0.0, 0.05, 0.1, 0.3, 0.5, np.nan], n)
phi = rng.choice([0.0, 0.05, 0.01, np.nan], n)
for female in (False, True):
    r = cpg_islands._chrom_kernel(start, meth, phi, female, 0.1, 0.05, 0.3,
                                  0.4, 0.6, 150, 150)
    h.update(np.array([r[0], r[6], r[7]]).tobytes())
    for x in r[1:6]:
        h.update(x[:r[0]].tobytes())
print(h.hexdigest())
"""


def test_interpreted_equals_jitted():
    """Checks that DMRDETECTOR_NO_JIT=1 gives the same results as the jit.

    Every kernel then runs as plain Python; the check runs in a separate process
    because the variable is read at import.
    """
    rng = np.random.default_rng(9)
    src = os.path.join(
        os.path.dirname(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "src")
    with tempfile.TemporaryDirectory() as d:
        txt = os.path.join(d, "t.bed")
        open(txt,
             "w").write(make_atlas(rng, 300, chroms=["chr1", "chr2", "chrX"]))
        script = os.path.join(d, "s.py")
        open(script, "w").write(_NO_JIT_SCRIPT)
        res = []
        for nojit in ("0", "1"):
            env = dict(os.environ, DMRDETECTOR_NO_JIT=nojit)
            p = subprocess.run([sys.executable, script, src, txt],
                               env=env,
                               capture_output=True,
                               text=True)
            assert p.returncode == 0, p.stderr
            res.append(p.stdout.strip())
        assert res[0] == res[1], res


def test_reader_numbers_and_stitching_match_pandas():
    rng = np.random.default_rng(11)
    vals = [
        0.0, 1.0, 0.1, 0.05552484840187, 7.73208042839979e-05, 1e-30, 1.5e-310,
        123456789012345678.0, 0.30000000000000004, 0.1 + 0.2, 5e-324,
        1.7976931348623157e308, 0.99999999999999989, 2.2250738585072014e-308,
        1e22, 1e23, 12345.678901234567
    ]
    vals += [
        float(x) for x in rng.random(300) * 10.0**rng.integers(-12, 3, 300)
    ]
    lines = []
    for i, v in enumerate(vals):
        lines.append("\t".join([
            "chr1",
            str(100 + 10 * i),
            str(101 + 10 * i), "9", "20",
            fmt(v, rng), "0.5",
            fmt(v, rng), "0.0"
        ]))
    text = "\n".join(lines) + "\n"
    want = pd.read_csv(__import__("io").StringIO(text),
                       sep="\t",
                       header=None,
                       names=oracle_cpg.ATLAS_COLS)
    # exact values: Python float() is correctly rounded; the pandas C parser is
    # up to 1 ulp off for 15 to 17 digits
    tok = [ln.split("\t") for ln in lines]
    want["PHI"] = [float(t[5]) for t in tok]
    want["REF_METH"] = [float(t[7]) for t in tok]
    with tempfile.TemporaryDirectory() as d:
        for block, batch, chunk in ((50, 1, 120), (50, 1 << 30, 400),
                                    (1000, 1, 3000), (65280, 4 << 20, 3 << 19)):
            atlas_io._CHUNK = chunk
            p = os.path.join(d, "x.gz")
            open(p, "wb").write(bgzf_bytes(text.encode(), block))
            code, start, end, phi, meth, cnt = atlas_io.read_columns(
                p,
                atlas_io.CPG,
                max_code=23,
                min_n=0,
                strict=True,
                threads=3,
                batch_bytes=batch)
            assert cnt[3] == 0 and start.tolist() == want["START"].tolist()
            assert phi.view(
                np.int64).tolist() == want["PHI"].to_numpy(float).view(
                    np.int64).tolist()
            assert meth.view(
                np.int64).tolist() == want["REF_METH"].to_numpy(float).view(
                    np.int64).tolist()
        # not understood -> Unsupported (callers fall back to pandas)
        for bad in ("chr1\t1\t2\t9\n", "chr1\t1\t2\t9\t20\t0.1\t0.2\t0.3x\t0\n",
                    "chr1\t1\t2\t9\t20\tinf\t0.2\t0.3\t0\n",
                    "chr1\t1\t2\t9\t20\t0.1\t0.2\t0.3\r\n"):
            p = os.path.join(d, "bad.bed")
            open(p, "w").write(bad)
            with pytest.raises(atlas_io.Unsupported):
                atlas_io.read_columns(p,
                                      atlas_io.CPG,
                                      max_code=23,
                                      min_n=0,
                                      strict=True)


def make_cands(rng):
    chroms = ["chr1", "chr2", "chr10", "chrX", "1", "Chr2"]
    rows = []
    for g in range(int(rng.integers(1, 7))):
        c = rng.choice(chroms)
        a = int(rng.integers(100, 3000))
        rows.append((c, a, a + int(rng.integers(20, 1500)),
                     f"G{int(rng.integers(0, 4))}"))
    return pd.DataFrame(rows,
                        columns=["CHROM", "START", "END", "GENE"
                                ]).assign(CHROM=lambda x: IC_norm(x["CHROM"]))


def IC_norm(s):
    return oracle_icr.normalize_primary_chrom(s)


def test_icr_matches_oracle_fuzz_all_routes():
    rng = np.random.default_rng(21)
    with tempfile.TemporaryDirectory() as d:
        for seed in range(80):
            text = make_atlas(rng,
                              n_rows=int(rng.choice([30, 400, 1500])),
                              chroms=["chr1", "chr2", "chr10", "chrX", "chr22"])
            if not text:
                continue
            cand = make_cands(rng)
            cand = cand[cand["CHROM"].isin(
                bed.AUTOSOMES)].reset_index(drop=True)
            atlas_io._CHUNK = int(rng.choice([150, 700, 3 << 19]))
            paths = write_variants(d, text, block=int(rng.choice([97, 5000])))
            for v1d in (True, False):
                kw = dict(v1_duplicates=v1d,
                          min_sites=int(rng.choice([1, 3, 5])),
                          log=_quiet)
                if seed % 4 == 0:
                    kw.update(min_samples=6,
                              merge_distance=int(rng.choice([149, 150, 151])))
                want = oracle_icr.detect(oracle_cpg.read_atlas(paths["plain"]),
                                         cand, **kw)
                _same(
                    icr.detect(oracle_cpg.read_atlas(paths["plain"]), cand,
                               **kw), want)
                for k in ("plain", "gz", "bgzf"):
                    _same(icr.detect(paths[k], cand, index="scan", **kw), want)
                    _same(icr.detect(atlas_io.read_atlas(paths[k]), cand, **kw),
                          want)


def test_icr_tabix_route_matches_oracle():
    pysam = pytest.importorskip("pysam")
    rng = np.random.default_rng(31)
    with tempfile.TemporaryDirectory() as d:
        for seed in range(25):
            text = make_atlas(rng,
                              n_rows=800,
                              chroms=["chr1", "chr2", "chr10", "chr22"])
            # tabix needs sorted input
            rows = sorted((ln.split("\t") for ln in text.splitlines()),
                          key=lambda r: (r[0], int(r[1])))
            text = "\n".join("\t".join(r) for r in rows) + "\n"
            cand = make_cands(rng)
            cand = cand[cand["CHROM"].isin(
                bed.AUTOSOMES)].reset_index(drop=True)
            if cand.empty:
                continue
            plain = os.path.join(d, "t.bed")
            open(plain, "w").write(text)
            gz = plain + ".gz"
            pysam.tabix_compress(plain, gz, force=True)
            want = oracle_icr.detect(oracle_cpg.read_atlas(plain),
                                     cand,
                                     log=_quiet)
            _same(icr.detect(gz, cand, index="build", log=_quiet),
                  want)  # index built in scratch
            assert not os.path.exists(
                gz + ".tbi")  # the input directory is untouched
            pysam.tabix_index(gz, preset="bed", force=True, keep_original=True)
            _same(icr.detect(gz, cand, index="auto", log=_quiet),
                  want)  # existing .tbi
            _same(icr.detect(gz, cand, index="scan", log=_quiet), want)
            os.remove(gz + ".tbi")


def test_icr_n_ge_and_window_edges():
    cand = pd.DataFrame({
        "CHROM": ["chr1"],
        "START": [1000],
        "END": [1100],
        "GENE": ["A"]
    })

    def run(rows, **kw):
        return icr.detect(_atlas(rows), cand, min_sites=1, log=_quiet, **kw)

    # N == 5 kept (>=), N == 4 dropped; REF_METH exactly 0.4 and 0.6 kept
    keep = [("chr1", 1000 + i, 1001 + i, 5, 10, 0.01, 0.5, m, 0.0)
            for i, m in enumerate((0.4, 0.6))]
    assert run(keep).values.tolist()[0][5] == 2
    assert run([("chr1", 1010, 1011, 4, 10, 0.01, 0.5, 0.5, 0.0)]).empty
    # site ending exactly at START (END == A) is not half-open overlapping: no
    # hit; START == B is a closed hit only
    assert run([("chr1", 999, 1000, 6, 10, 0.01, 0.5, 0.5, 0.0)]).empty
    r = run([("chr1", 1100, 1101, 6, 10, 0.01, 0.5, 0.5, 0.0)])
    _same(
        r,
        oracle_icr.detect(_atlas([("chr1", 1100, 1101, 6, 10, 0.01, 0.5, 0.5,
                                   0.0)]),
                          cand,
                          min_sites=1,
                          log=_quiet))


def test_cli_defaults_come_from_defaults_module():
    import argparse
    root = argparse.ArgumentParser()
    sub = root.add_subparsers()
    cpg_islands.add_parser(sub)
    icr.add_parser(sub)
    c = root.parse_args(["cpg-islands", "--ref-atlas", "x", "--out", "y"])
    assert (c.include_sex_chromosomes, c.min_samples_per_site,
            c.core_max_methylation, c.core_max_dispersion,
            c.shore_max_methylation, c.sex_f_ref_meth_min, c.sex_f_ref_meth_max,
            c.max_gap_distance, c.min_island_length) == (True, 5, 0.1, 0.05,
                                                         0.3, 0.40, 0.60, 150,
                                                         150)
    assert (c.min_samples_per_site, c.core_max_methylation,
            c.max_gap_distance) == (settings.CPG_MIN_SAMPLES_PER_SITE,
                                    settings.CPG_CORE_MAX_METHYLATION,
                                    settings.CPG_MAX_GAP_DISTANCE)
    i = root.parse_args(["icr", "--ref-atlas", "x", "--icr-gene-list", "g"])
    assert (i.icr_ref_meth_min, i.icr_ref_meth_max, i.icr_min_sites,
            i.dmr_merge_distance, i.min_samples_per_site, i.no_v1_duplicates,
            i.index) == (0.40, 0.60, 5, 150, 5, False, "auto")
    assert root.parse_args([
        "icr", "--ref-atlas", "x", "--icr-gene-list", "g", "--no-v1-duplicates"
    ]).no_v1_duplicates
    assert settings.ICR_V1_DUPLICATES is True
