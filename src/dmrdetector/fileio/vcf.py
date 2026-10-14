"""Heterozygous-SNP scan of a VCF over raw bytes (numba).

A record counts when some ALT is a SNP by the htslib typing of `alt_is_snp` and
some sample has two or more distinct called alleles in GT. The kernel reads
ordinary VCF lines only; on anything unusual it sets a flag and the caller uses
the pysam loop for the whole file, so the result never depends on how the kernel
treats odd input.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from dmrdetector.util.jit import jit


if TYPE_CHECKING:
    import pysam

_TAB, _NL, _COMMA, _COLON, _SLASH, _PIPE, _DOT = 9, 10, 44, 58, 47, 124, 46


@jit
def _up(c):
    return c - 32 if 97 <= c <= 122 else c


@jit
def _snp(buf, rb, re_, ab, ae):
    """Returns alt_is_snp of REF = buf[rb:re_] and ALT = buf[ab:ae]."""
    lr = re_ - rb
    la = ae - ab
    c0 = buf[ab]
    if la == 1 and c0 == 42:  # '*'
        return False
    if lr == 1 and la == 1:
        return c0 != 46 and c0 != 88 and buf[rb] != c0
    if c0 == 60 or c0 == 91 or c0 == 93:  # '<' '[' ']'
        return False
    i = 0
    while i < lr and i < la and _up(buf[rb + i]) == _up(buf[ab + i]):
        i += 1
    if i == lr or i == la:
        return False
    r = lr - 1
    a = la - 1
    ca = buf[ab + a]
    if ca == 91 or ca == 93:
        return False
    if _up(buf[rb + r]) == _up(ca) and r != i and a != i:
        while r > i and a > i and _up(buf[rb + r]) == _up(buf[ab + a]):
            r -= 1
            a -= 1
    return r == i and a == i


@jit
def _het_gt(buf, s, e, k, nalt):
    """Returns the genotype of the sample field buf[s:e], subfield index k.

    Returns:
        1 for heterozygous, 0 for not, -1 for an anomaly.
    """
    # locate subfield k
    p = s
    j = 0
    while j < k:
        while p < e and buf[p] != _COLON:
            p += 1
        if p >= e:
            return 0  # trailing subfield absent: GT missing
        p += 1
        j += 1
    first = -1
    het = 0
    nall = 0
    while p < e and buf[p] != _COLON:
        c = buf[p]
        if c == _DOT:
            p += 1
            nall += 1
        elif 48 <= c <= 57:
            v = 0
            while p < e and 48 <= buf[p] <= 57:
                v = v * 10 + (int(buf[p]) - 48)
                if v > 1000000:
                    return -1
                p += 1
            if v > nalt:
                return -1
            nall += 1
            if first < 0:
                first = v
            elif v != first:
                het = 1
        else:
            return -1
        if p < e and buf[p] != _COLON:
            if buf[p] != _SLASH and buf[p] != _PIPE:
                return -1
            p += 1
            if p >= e or buf[p] == _COLON:
                return -1
    if nall == 0:
        return -1
    return het


@jit
def _scan(buf, name_off, name_len):
    """Scan whole lines.

    Returns (flag, nhit, hit_cid, hit_pos, nnames); name_off/name_len receive
    the contig names (chunk-local ids) as byte ranges of buf.
    """
    n = buf.shape[0]
    cap = 0
    for i in range(n):
        if buf[i] == _NL:
            cap += 1
    cap += 1
    hit_cid = np.empty(cap, dtype=np.int32)
    hit_pos = np.empty(cap, dtype=np.int64)
    nhit = 0
    nnames = 0
    last_cid = -1
    p = 0
    while p < n:
        if buf[p] == _NL:  # blank line
            return 1, nhit, hit_cid, hit_pos, nnames
        if buf[p] == 35:  # '#'
            while p < n and buf[p] != _NL:
                p += 1
            p += 1
            continue
        # field boundaries of this line
        e = p
        while e < n and buf[e] != _NL:
            e += 1
        # fields 0..8: start and end offsets
        fs = np.empty(9, dtype=np.int64)
        fe = np.empty(9, dtype=np.int64)
        nf = 0
        fs[0] = p
        q = p
        while q <= e and nf < 9:
            if q == e or buf[q] == _TAB:
                fe[nf] = q
                nf += 1
                if nf < 9:
                    fs[nf] = q + 1
                if q == e:
                    break
            q += 1
        if nf < 8:
            return 1, nhit, hit_cid, hit_pos, nnames
        # POS
        pos = 0
        if fe[1] <= fs[1]:
            return 1, nhit, hit_cid, hit_pos, nnames
        for x in range(fs[1], fe[1]):
            c = buf[x]
            if c < 48 or c > 57:
                return 1, nhit, hit_cid, hit_pos, nnames
            pos = pos * 10 + (int(c) - 48)
        rb = fs[3]
        re_ = fe[3]
        ab = fs[4]
        ae = fe[4]
        if re_ <= rb or ae <= ab:
            return 1, nhit, hit_cid, hit_pos, nnames
        for x in range(rb, re_):
            if buf[x] >= 128 or buf[x] == 13:
                return 1, nhit, hit_cid, hit_pos, nnames
        for x in range(ab, ae):
            if buf[x] >= 128 or buf[x] == 13:
                return 1, nhit, hit_cid, hit_pos, nnames
        if buf[e - 1] == 13:
            return 1, nhit, hit_cid, hit_pos, nnames
        if ae - ab == 1 and buf[ab] == _DOT:  # no ALT
            p = e + 1
            continue
        # any ALT a SNP? (comma-separated)
        snp = False
        nalt = 1
        s0 = ab
        for x in range(ab, ae + 1):
            if x == ae or buf[x] == _COMMA:
                if x == s0:
                    return 1, nhit, hit_cid, hit_pos, nnames
                if not snp and _snp(buf, rb, re_, s0, x):
                    snp = True
                s0 = x + 1
                if x < ae:
                    nalt += 1
        if not snp:
            p = e + 1
            continue
        # FORMAT: index of GT
        if nf < 9 or fe[8] <= fs[8]:
            p = e + 1
            continue
        # find the GT key
        k = -1
        j = 0
        ks = fs[8]
        t = fs[8]
        while t <= fe[8]:
            if t == fe[8] or buf[t] == _COLON:
                if t - ks == 2 and buf[ks] == 71 and buf[ks + 1] == 84:
                    k = j
                    break
                j += 1
                ks = t + 1
            t += 1
        if k < 0:
            p = e + 1
            continue
        # htslib does not expose a GT that is not the first FORMAT key: leave it
        # to the pysam loop
        if k > 0:
            return 1, nhit, hit_cid, hit_pos, nnames
        # samples
        het = 0
        s = fe[8] + 1
        while s <= e and het == 0:
            f = s
            while f < e and buf[f] != _TAB:
                f += 1
            r = _het_gt(buf, s, f, k, nalt)
            if r < 0:
                return 1, nhit, hit_cid, hit_pos, nnames
            het = r
            s = f + 1
        if het:
            # contig id (chunk-local)
            cs = fs[0]
            ce = fe[0]
            cid = -1
            if last_cid >= 0 and name_len[last_cid] == ce - cs:
                same = True
                for x in range(ce - cs):
                    if buf[name_off[last_cid] + x] != buf[cs + x]:
                        same = False
                        break
                if same:
                    cid = last_cid
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
                    return 1, nhit, hit_cid, hit_pos, nnames
                name_off[nnames] = cs
                name_len[nnames] = ce - cs
                cid = nnames
                nnames += 1
            last_cid = cid
            hit_cid[nhit] = cid
            hit_pos[nhit] = pos - 1
            nhit += 1
        p = e + 1
    return 0, nhit, hit_cid, hit_pos, nnames


def _chunk(buf: np.ndarray, first: bool = False):
    off = np.zeros(65536, dtype=np.int64)
    ln = np.zeros(65536, dtype=np.int64)
    flag, nhit, cid, pos, nn = _scan(buf, off, ln)
    names = [bytes(buf[off[i]:off[i] + ln[i]]).decode() for i in range(nn)]
    return flag, names, cid[:nhit].copy(), pos[:nhit].copy()


def scan_het_sites(vcf: str,
                   threads: int | None = None) -> dict[str, np.ndarray] | None:
    """Returns {contig: sorted unique 0-based het-SNP positions}.

    Returns None when the file needs the pysam loop.
    """
    from dmrdetector.fileio import bed
    res = bed.map_text(vcf, _chunk, threads)
    sites: dict[str, list[np.ndarray]] = {}
    for flag, names, cid, pos in res:
        if flag:
            return None
        for i, nm in enumerate(names):
            sites.setdefault(nm, []).append(pos[cid == i])
    return {
        c: np.unique(np.concatenate(v))
        for c, v in sites.items()
        if any(x.size for x in v)
    }


def sample_count(vcf: str) -> int:
    """Returns the number of sample columns of the #CHROM line (0 if none).

    Args:
        vcf: VCF, plain or gzip.
    """
    import gzip
    with open(vcf, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    with (gzip.open(vcf, "rt") if gz else open(vcf)) as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                return max(0, len(line.rstrip("\r\n").split("\t")) - 9)
            if not line.startswith("#"):
                break
    return 0


def alt_is_snp(ref: str, alt: str) -> bool:
    """Returns True when htslib types REF>ALT as a SNP.

    This is VCF_SNP in vcf.c bcf_set_variant_type, the test behind `bcftools
    view -v snps`. Besides a single-base REF and ALT it includes equal-length
    alleles that differ at one base once the common prefix and suffix are
    trimmed (REF CGTGTGT, ALT TGTGTGT).
    """
    if alt == "*":
        return False
    if len(ref) == 1 and len(alt) == 1:
        return alt != "." and alt != "X" and ref != alt
    if alt[0] == "<" or alt[0] in "[]":
        return False
    i = 0
    while i < len(ref) and i < len(alt) and ref[i].upper() == alt[i].upper():
        i += 1
    if i == len(ref) or i == len(alt):  # pure insertion or deletion
        return False
    re_, ae = len(ref) - 1, len(alt) - 1
    if alt[ae] in "[]":
        return False
    if ref[re_].upper() == alt[ae].upper() and re_ != i and ae != i:
        while re_ > i and ae > i and ref[re_].upper() == alt[ae].upper():
            re_ -= 1
            ae -= 1
    return re_ == i and ae == i


def het_position(rec: "pysam.VariantRecord") -> int | None:
    """0-based position of a pysam record that is a heterozygous SNP, else None.

    A record counts when some ALT is a SNP by `alt_is_snp` and some sample has
    two or more distinct called alleles in GT (phasing is not needed).
    """
    ref, alts = rec.ref, rec.alts
    if alts is None or ref is None:
        return None
    if not any(alt_is_snp(ref, alt) for alt in alts):
        return None
    for call in rec.samples.values():
        gt = call.get("GT")
        if gt is not None and len({a for a in gt if a is not None}) > 1:
            return rec.pos - 1
    return None


def _het_snp_sites_pysam(vcf: str) -> dict[str, np.ndarray]:
    """Reference implementation of het_snp_sites (pysam record loop)."""
    import pysam
    sites: dict[str, list[int]] = {}
    with pysam.VariantFile(vcf) as vf:
        for rec in vf:
            pos = het_position(rec)
            if pos is not None:
                sites.setdefault(rec.chrom, []).append(pos)
    return {
        c: np.unique(np.asarray(p, dtype=np.int64)) for c, p in sites.items()
    }


def het_snp_sites(vcf: str,
                  threads: int | None = None) -> dict[str, np.ndarray]:
    """Sorted 0-based positions of heterozygous SNPs per contig (exact name).

    The byte-level scan runs first; a file with anything unusual goes through
    the pysam loop, the reference rule.
    """
    sites = scan_het_sites(vcf, threads)
    return sites if sites is not None else _het_snp_sites_pysam(vcf)
