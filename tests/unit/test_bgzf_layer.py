"""Tests of the BGZF layer: scan, inflate, writer, tabix index, number parsing.

Each part is checked against an independent reference (gzip, pysam/htslib,
Python float). Run it with the jit, with DMRDETECTOR_NO_JIT=1 (plain Python
kernels) and with DMRDETECTOR_NO_LIBDEFLATE=1 (zlib).
"""
import gzip
import os
import struct
import zlib

import numpy as np
import pysam
import pytest

from dmrdetector.fileio import bgzf as bz  # noqa: E402
from dmrdetector.fileio import numparse  # noqa: E402

from helpers import parse_tbi


def _blocks_with_extra_field(data: bytes, step: int = 5000) -> bytes:
    """Returns BGZF with a subfield before BC in the block headers (XLEN 12).

    It is valid but not the htslib header, so the byte-range split cannot be
    used and the chain walk must read XLEN.
    """
    out = b""
    for i in range(0, len(data), step):
        chunk = data[i:i + step]
        c = zlib.compressobj(6, zlib.DEFLATED, -15)
        body = c.compress(chunk) + c.flush()
        extra = b"XX\x02\x00ab" + b"BC\x02\x00" + struct.pack(
            "<H",
            len(body) + 12 + 12 + 8 - 1)
        header = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff"
        trailer = struct.pack("<II", zlib.crc32(chunk), len(chunk))
        out += (header + struct.pack("<H", len(extra)) + extra + body + trailer)
    return out + bz.EOF_BLOCK


def test_scan_and_inflate_standard_and_extra_field(tmp_path):
    rng = np.random.default_rng(0)
    data = rng.integers(0, 256, 300_000,
                        dtype=np.uint8).tobytes() + b"abc\n" * 50_000
    std = b"".join(bz.compress_bytes(data, 3)) + bz.EOF_BLOCK
    for name, raw in (("std", std), ("extra", _blocks_with_extra_field(data))):
        p = tmp_path / f"{name}.gz"
        p.write_bytes(raw)
        assert gzip.open(p).read() == data
        arr = np.frombuffer(raw, np.uint8)
        coff, csize, usize, doff = bz.scan_blocks(arr)
        assert int(
            usize.sum()) == len(data) and int(coff[-1] + csize[-1]) == len(raw)
        assert int(doff[0]) == (18 if name == "std" else 24)
        triples = list(
            zip(doff.tolist(), (coff + csize - 8 - doff).tolist(),
                usize.tolist()))
        assert bytes(bz.inflate_group(raw, 0, triples)) == data
        buf = bz.inflate(str(p), 3)
        assert buf.tobytes() == data + (b"" if data.endswith(b"\n") else b"\n")
        assert bz.count_newlines(str(p), 3) == data.count(b"\n")
    # byte-range split of a standard file: every block is found exactly once
    ranges = [(a, min(a + 4001, len(std))) for a in range(0, len(std), 4001)]
    got = []
    for a, b in ranges:
        got += bz.blocks_starting_in(std[a:min(b + 65536, len(std))], a, b,
                                     len(std))
    assert sum(g[2] for g in got) == len(data) and len(got) == len(
        bz.scan_blocks(np.frombuffer(std, np.uint8))[2])
    assert bz.scan_blocks(np.frombuffer(gzip.compress(data),
                                        np.uint8)) is None  # plain gzip
    assert bz.scan_blocks(np.frombuffer(b"chr1\t1\t2\n" * 10,
                                        np.uint8)) is None  # text


def test_corrupt_block_is_an_error(tmp_path):
    raw = bytearray(b"".join(bz.compress_bytes(b"chr1\t5\t6\n" * 20000)) +
                    bz.EOF_BLOCK)
    raw[100] ^= 0xFF
    p = tmp_path / "bad.gz"
    p.write_bytes(bytes(raw))
    with pytest.raises(bz.ParseError):
        bz.inflate(str(p), 2)


@pytest.mark.parametrize("level", [1, 4, 6])
def test_writer_offsets_and_order(tmp_path, level):
    rng = np.random.default_rng(1)
    lines = [
        f"chr1\t{i}\t{i + 1}\t{'x' * int(rng.integers(0, 90))}\n".encode()
        for i in range(30_000)
    ]
    jobs = [(a, min(a + 4000, len(lines))) for a in range(0, len(lines), 4000)]

    def make(job):
        txt = b"".join(lines[job[0]:job[1]])
        rowoff = np.concatenate(
            ([0], np.cumsum([len(x) for x in lines[job[0]:job[1]]
                            ]))).astype(np.int64)
        return np.frombuffer(txt, np.uint8), rowoff, job[1] - job[0]

    outs = []
    for kw in ({}, {"order": list(range(len(jobs)))[::-1]}):
        p = str(tmp_path / f"w{len(outs)}.gz")
        voff, final = bz.write_segments(p, jobs, make, 4, level, **kw)
        outs.append(open(p, "rb").read())
        assert gzip.open(p).read() == b"".join(lines) and outs[-1].endswith(
            bz.EOF_BLOCK)
        assert final == (len(outs[-1]) - len(bz.EOF_BLOCK)) << 16
    assert outs[0] == outs[1]  # start order does not change the file


def _write_bed(path, rows, level=1):
    """Writes a BED of `rows` [(chrom, start, end)], sorted.

    Returns:
        (names, contigs, final) as the writer reports them.
    """
    names = list(dict.fromkeys(r[0] for r in rows))
    lines = [f"{c}\t{s}\t{e}\t1\t2\n".encode() for c, s, e in rows]
    jobs = [(a, min(a + 3000, len(lines))) for a in range(0, len(lines), 3000)]

    def make(job):
        seg = lines[job[0]:job[1]]
        rowoff = np.concatenate(
            ([0], np.cumsum([len(x) for x in seg]))).astype(np.int64)
        return np.frombuffer(b"".join(seg), np.uint8), rowoff, len(seg)

    voff, final = bz.write_segments(path, jobs, make, 3, level)
    contigs, o = [], 0
    for c in names:
        sel = [r for r in rows if r[0] == c]
        contigs.append(
            (np.array([r[1] for r in sel],
                      np.int64), np.array([r[2] for r in sel],
                                          np.int64), voff[o:o + len(sel)]))
        o += len(sel)
    return names, contigs, final


def _rows(rng, n, wide):
    out = []
    for chrom in ("chr1", "chr10", "chr2", "chrX"):
        k = int(rng.integers(0, n))
        starts = np.sort(rng.integers(0, 3_000_000, k))
        for s in starts.tolist():
            w = 1
            if wide and rng.random() < 0.02:
                w = int(rng.choice([2, 150, 16384, 16385, 40_000, 300_000]))
            out.append((chrom, s, s + w))
    return out


@pytest.mark.parametrize("wide", [False, True])
def test_tbi_equals_htslib(tmp_path, wide):
    """Checks that the index from the writer's offsets equals htslib's.

    The records are 1 bp (every BED of the pipeline) or span several 16 kb
    windows.
    """
    rng = np.random.default_rng(3 if wide else 4)
    bad = []
    for it in range(6):
        rows = _rows(rng, 25_000 if it % 2 == 0 else 300, wide)
        if not rows:
            continue
        p = str(tmp_path / f"t{it}.bed.gz")
        names, contigs, final = _write_bed(p, rows)
        assert bz.write_tbi(p, names, contigs, final)
        ref = str(tmp_path / f"r{it}.bed.gz")
        open(ref, "wb").write(open(p, "rb").read())
        pysam.tabix_index(ref, preset="bed", force=True, keep_original=True)
        if parse_tbi(p + ".tbi") != parse_tbi(ref + ".tbi"):
            bad.append(it)
        with pysam.TabixFile(p) as tbx:  # and the index works
            for c in names:
                assert sum(1 for _ in tbx.fetch(c)) == sum(
                    1 for r in rows if r[0] == c)
    assert not bad, f"tbi differs from htslib for iterations {bad}"


def test_tbi_refuses_end_beyond_tabix_limit(tmp_path):
    p = str(tmp_path / "big.bed.gz")
    names, contigs, final = _write_bed(p, [("chr1", 5, 6),
                                           ("chr1", 600_000_000, 600_000_001)])
    assert bz.write_tbi(p, names, contigs, final) is False
    assert not os.path.exists(p + ".tbi")


def test_parse_decimal_never_off_by_an_ulp():
    rng = np.random.default_rng(5)
    toks = []
    for _ in range(20000):  # 16 to 19 digits: the double-double path
        nd = int(rng.integers(16, 20))
        d = "".join(map(str, rng.integers(0, 10, nd)))
        toks.append(f"0.{d}" if rng.random() <
                    0.5 else f"{d[:3]}.{d[3:]}e{int(rng.integers(-30, 30))}")
    # tokens next to a tie: the exact midpoint of two doubles, cut to 19 digits,
    # and its neighbors
    for v in rng.random(3000) * 10.0**rng.integers(-5, 5, 3000):
        lo, hi = float(v), float(np.nextafter(v, np.inf))
        from fractions import Fraction
        mid = (Fraction(lo) + Fraction(hi)) / 2
        s = f"{float(mid):.18e}"
        toks += [s]
    toks += [
        "9007199254740993", "9007199254740995", "18014398509481985",
        "1.00000000000000011102230246251565404"
    ]
    text = ("\n".join(toks) + "\n").encode()
    buf = np.frombuffer(text, np.uint8)
    pos = np.concatenate(([0], np.flatnonzero(buf == 10) + 1))
    n_slow = n_short = 0
    for t, a, b in zip(toks, pos[:-1], pos[1:]):
        st, v = numparse.parse_decimal(buf, int(a), int(b) - 1)
        short = sum(c.isdigit() for c in t.split("e")[0].lstrip("0.")) <= 18
        n_short += short
        if st:
            n_slow += short
        else:
            assert v == float(t), t
    # the hand-over to float() is rare below 19 digits
    assert n_slow < n_short // 100


def test_one_fwrite_formatter():
    """Checks that the scalar and the array fwrite formatters follow one rule.

    bed.fmt_double (scalar, annotate and target tables) and bed.fmt_doubles
    (arrays, atlas text) differ only for NaN (an empty cell in the scalar
    formatter) and negative zero (scalar "0", array "-0"), as before.
    """
    from dmrdetector.fileio import bed as dio
    rng = np.random.default_rng(7)
    x = np.concatenate([
        rng.random(50_000),
        rng.random(50_000) * 10.0**rng.integers(-30, 30, 50_000),
        np.round(rng.random(5000), 3),
        np.array([
            0.0, 1e-4, 1e-5, 1e5, 123456.0, 1e15, 1e-300, np.inf, -np.inf,
            5e-324
        ])
    ])
    x = np.concatenate([x, -x[:1000]])
    assert [dio.fmt_double(v) for v in x] == dio.fmt_doubles(x)
    assert dio.fmt_double(float("nan")) == ""
