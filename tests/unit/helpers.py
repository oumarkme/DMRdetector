"""Synthetic inputs and small readers shared by the unit tests."""

from __future__ import annotations

import gzip
import os
import struct

import numpy as np
import pandas as pd
import pysam

from dmrdetector.fileio import bgzf


CONTIGS = {  # name: (n_sites, [(lo, hi), ...])
    "chr21": (3000, [(5_010_000, 5_400_000)]),
    "chr22": (2000, [(16_050_000, 16_300_000)]),
    "chrX": (1500, [(9_000, 2_800_000), (60_000_000, 60_300_000),
                    (155_690_000, 156_040_000)]),
    "chrY": (600, [(9_000, 2_800_000), (20_000_000, 20_100_000),
                   (56_880_000, 57_215_000)]),
}


def write_bgzf(path, data):
    """Writes `data` (bytes) as one BGZF file."""
    with bgzf.BgzfWriter(path) as writer:
        writer.write_blocks(bgzf.compress_bytes(data))


def read_bed5(path):
    """Reads a 5-column BED (CHROM START END COVERAGE N_METH), plain or gz."""
    return pd.read_csv(path,
                       sep="\t",
                       header=None,
                       usecols=[0, 1, 2, 3, 4],
                       names=["CHROM", "START", "END", "COVERAGE", "N_METH"],
                       dtype={0: str},
                       comment="#")


def parse_tbi(path):
    """Parses a tabix index.

    Returns:
        (header fields, names, [(bins, intervals)] per reference, tail bytes).
    """
    with gzip.open(path, "rb") as fh:
        data = fh.read()
    head = struct.unpack_from("<4s8i", data, 0)
    n_ref, name_len = head[1], head[8]
    off = 36
    names = data[off:off + name_len].split(b"\0")[:n_ref]
    off += name_len
    refs = []
    for _ in range(n_ref):
        (n_bin,) = struct.unpack_from("<i", data, off)
        off += 4
        bins = {}
        for _ in range(n_bin):
            bin_id, n_chunk = struct.unpack_from("<Ii", data, off)
            off += 8
            bins[bin_id] = [
                struct.unpack_from("<QQ", data, off + 16 * k)
                for k in range(n_chunk)
            ]
            off += 16 * n_chunk
        (n_intv,) = struct.unpack_from("<i", data, off)
        off += 4
        intervals = list(struct.unpack_from(f"<{n_intv}Q", data, off))
        off += 8 * n_intv
        refs.append((bins, intervals))
    return head, names, refs, data[off:]


def _sites(rng, name):
    n, regions = CONTIGS[name]
    per = n // len(regions)
    return np.unique(
        np.concatenate([rng.integers(lo, hi, per) for lo, hi in regions]))


def _write_lines(path, lines):
    with gzip.open(path, "wt") as fh:
        fh.write("\n".join(lines) + "\n")


def generate(outdir, n=6, seed=1):
    """Writes n synthetic samples S01.. and returns their ids.

    Per sample: S0k.{5col,modkit,pbcpg}.bed.gz (the same CpG table in three
    formats) and S0k.vcf (het SNPs on CpG positions, hom SNPs, indels, a contig
    name that is not in the BED). Coordinates are 0-based BED.
    """
    os.makedirs(outdir, exist_ok=True)
    rng = np.random.default_rng(seed)
    base = {c: _sites(rng, c) for c in CONTIGS}
    ids = []
    for k in range(1, n + 1):
        rows = []  # chrom, start, end, cov, nmod
        for c, pos in base.items():
            p = pos[rng.random(pos.size) < 0.9]
            cov = np.clip(
                rng.poisson(20, p.size) + rng.integers(0, 4, p.size), 1, None)
            pm = np.where((p // 5000) % 3 == 0, 0.05, 0.85)
            nm = rng.binomial(cov, pm)
            width = np.where(rng.random(p.size) < 0.05, 2, 1)
            for a, w, cv, m in zip(p, width, cov, nm):
                rows.append((c, int(a), int(a + w), int(cv), int(m)))
            # a duplicated start, as modkit writes with two mod codes
            for i in rng.choice(p.size, 3, replace=False):
                rows.append(
                    (c, int(p[i]), int(p[i]) + 1, int(cov[i]) + 2, int(nm[i])))
        rows.sort(key=lambda r: (r[0], r[1]))
        sid = f"S{k:02d}"
        _write_lines(os.path.join(outdir, f"{sid}.5col.bed.gz"),
                     ["\t".join(map(str, r)) for r in rows])
        modkit = []
        for c, a, b, cv, m in rows:
            head = [
                c,
                str(a),
                str(b), "m",
                str(min(cv, 1000)), "+",
                str(a),
                str(b), "255,0,0"
            ]
            tail = [
                str(cv), f"{100 * m / cv:.2f}",
                str(m),
                str(cv - m), "0", "0", "0", "0", "0"
            ]
            modkit.append("\t".join(head + tail))
        _write_lines(os.path.join(outdir, f"{sid}.modkit.bed.gz"), modkit)
        pbcpg = [
            "#chrom\tbegin\tend\tmod_score\ttype\tcov\test_mod_count\t"
            "est_unmod_count\tdiscretized_mod_score"
        ]
        for c, a, b, cv, m in rows:
            count = f"{m:.2f}" if m % 5 == 0 else str(m)
            pbcpg.append("\t".join([
                c,
                str(a),
                str(b), f"{100 * m / cv:.1f}", "Total",
                str(cv), count,
                str(cv - m), "0"
            ]))
        _write_lines(os.path.join(outdir, f"{sid}.pbcpg.bed.gz"), pbcpg)
        _write_vcf(os.path.join(outdir, f"{sid}.vcf"), sid, rows, rng)
        ids.append(sid)
    return ids


def _write_vcf(path, sid, rows, rng):
    vcf = [
        "##fileformat=VCFv4.2", "##contig=<ID=chr21>", "##contig=<ID=chrX>",
        "##contig=<ID=chrY>", "##contig=<ID=21>",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sid
    ]
    pick = {c: [r for r in rows if r[0] == c] for c in CONTIGS}
    recs = []
    for c in ("chr21", "chrX", "chrY"):
        rr = pick[c]
        for j, i in enumerate(rng.choice(len(rr), 40, replace=False)):
            r = rr[i]
            pos1 = r[1] + 1 if j % 3 else r[2]  # first or last base of a CpG
            gt = ["0/1", "1/1", "0|1", "1/2", "./1", "0/0"][j % 6]
            ref, alt = ("C", "T") if j % 5 else ("CG", "C")  # a few indels
            if j % 11 == 0:
                alt = "T,G"
            if gt == "1/2" and "," not in alt:
                alt = "T,G"
            recs.append(
                (c, pos1, ref, alt, gt, "q10" if j % 4 == 0 else "PASS"))
    r = pick["chr21"][5]
    recs.append(("21", r[1] + 1, "C", "T", "0/1", "PASS"))  # not in the BED
    recs.sort(key=lambda x: (x[0], x[1]))
    for c, p1, ref, alt, gt, flt in recs:
        vcf.append(f"{c}\t{p1}\t.\t{ref}\t{alt}\t50\t{flt}\t.\tGT\t{gt}")
    with open(path, "w") as fh:
        fh.write("\n".join(vcf) + "\n")


def write_call_dmr_inputs(out, gtf=None, heavy=False):
    """Writes sample and atlas BEDs, CpG panels and an ICR panel into `out`.

    Per scope (autosomal, female_X, male_X): <scope>.sample.bed.gz,
    <scope>.atlas.bed.gz (both with .tbi) and <scope>.cpg.bed; plus icr.bed and
    icr_genes.bed. Blocks of shifted methylation sit at random gene positions of
    `gtf` when given, so annotation finds overlaps; several panel regions
    overlap or repeat. `heavy` raises the overlap.
    """
    os.makedirs(out, exist_ok=True)
    rng = np.random.default_rng(7)
    genes = None
    if gtf:
        from dmrdetector.steps.annotate import load_genes
        genes = load_genes(gtf)
        genes = genes[genes["seq"].isin(["1", "2", "10", "X"]) &
                      (genes["end"] - genes["start"] > 3000)]

    def block_positions(chrom_ens, n_blocks):
        if genes is None:
            return [(int(rng.integers(1_000_000, 50_000_000)), 6000)
                    for _ in range(n_blocks)]
        g = genes[genes["seq"] == chrom_ens].sample(n_blocks,
                                                    random_state=int(
                                                        rng.integers(1e6)))
        return [(int(s) + int(rng.integers(0, 500)),
                 int(rng.integers(1500, 4000))) for s in g["start"]]

    def make_scope(scope, chroms, shift_scale):
        regions, atlas_rows, sample_rows = [], [], []
        for chrom, ens in chroms:
            for bs, bl in block_positions(ens, 40):
                pos = np.unique(bs +
                                np.cumsum(rng.integers(8, 90, size=bl // 40)))
                ref = np.clip(rng.beta(2, 2, size=pos.size), 0.02, 0.98)
                # three regimes: strong shift, partial shift, none
                kind = rng.integers(0, 3)
                shift = (shift_scale if rng.random() < 0.5 else
                         -shift_scale) if kind < 2 else 0.0
                cov = rng.integers(8, 45, size=pos.size)
                scale = 1.0 if kind == 0 else rng.random(pos.size)
                frac = np.clip(
                    ref + shift * scale + rng.normal(0, 0.05, pos.size), 0, 1)
                nm = np.round(frac * cov).astype(int)
                for p, r, c, n in zip(pos, ref, cov, nm):
                    atlas_rows.append(
                        (chrom, int(p), int(p) + 1, 5, 120, 0.05, r, r, 0.01))
                    sample_rows.append(
                        (chrom, int(p), int(p) + 1, int(c), int(n)))
                lo, hi = int(pos[0]), int(pos[-1]) + 2
                regions.append((chrom, lo + 5, hi - 5, lo, hi))
                if rng.random() < (0.3 if heavy else 0.04):  # overlapping
                    mid = (lo + hi) // 2
                    regions.append((chrom, mid, hi + 400, mid - 20, hi + 420))
                if rng.random() < (0.15 if heavy else 0.04):  # duplicate
                    regions.append((chrom, lo + 5, hi - 5, lo, hi))
            # background sites outside any panel region
            bg = np.unique(
                rng.integers(100_000_000,
                             120_000_000,
                             size=3000 if heavy else 20000))
            for p in bg:
                r = float(np.clip(rng.beta(2, 2), 0.02, 0.98))
                atlas_rows.append(
                    (chrom, int(p), int(p) + 1, 5, 120, 0.05, r, r, 0.01))
                sample_rows.append((chrom, int(p), int(p) + 1, 20,
                                    int(round(20 * rng.random()))))
        for name, rows in (("sample", sample_rows), ("atlas", atlas_rows)):
            df = pd.DataFrame(rows).drop_duplicates([0, 1]).sort_values(
                [0, 1], kind="stable")
            # drop some sites: the inner join is not trivial
            if name == "sample":
                df = df.sample(frac=0.97,
                               random_state=1).sort_values([0, 1],
                                                           kind="stable")
            path = os.path.join(out, f"{scope}.{name}.bed")
            df.to_csv(path, sep="\t", header=False, index=False)
            pysam.tabix_compress(path, path + ".gz", force=True)
            pysam.tabix_index(path + ".gz", preset="bed", force=True)
            os.remove(path)
        cpg = pd.DataFrame(
            regions,
            columns=["CHROM", "START", "END", "START_EXTEND", "END_EXTEND"])
        cpg["LENGTH"] = cpg["END"] - cpg["START"]
        cpg["NUM_CPG"] = 10
        cpg = cpg.sort_values(["CHROM", "START"], kind="stable")
        cpg.to_csv(os.path.join(out, f"{scope}.cpg.bed"), sep="\t", index=False)

    make_scope("autosomal", [("chr1", "1"), ("chr2", "2"), ("chr10", "10")],
               0.55)
    make_scope("female_X", [("chrX", "X")], 0.35)
    make_scope("male_X", [("chrX", "X")], 0.55)

    # ICR panel: shrunken CpG panel regions, overlapping and repeated
    cpg = pd.read_csv(os.path.join(out, "autosomal.cpg.bed"), sep="\t")
    sub = cpg.sample(frac=0.5, random_state=3)
    icr = pd.DataFrame({
        "CHROM": sub["CHROM"],
        "START": sub["START"] - 200,
        "END": sub["END"] + 50
    })
    icr["GENE"] = [f"G{i % 17}" for i in range(len(icr))]
    icr["ICR_ID"] = [f"icr{i}" for i in range(len(icr))]
    icr["NUM_CPG"] = 6
    icr = pd.concat(
        [icr, icr.iloc[:3], icr.iloc[:2].assign(GENE="GALT", ICR_ID="dup")])
    icr = icr.sort_values(["CHROM", "START"], kind="stable")
    icr.to_csv(os.path.join(out, "icr.bed"), sep="\t", index=False)
    genes_bed = pd.DataFrame({
        "CHROM": icr["CHROM"],
        "START": icr["START"] - 1000,
        "END": icr["END"] + rng.integers(0, 5000, len(icr)),
        "GENE": icr["GENE"]
    })
    genes_bed.drop_duplicates().to_csv(os.path.join(out, "icr_genes.bed"),
                                       sep="\t",
                                       index=False)
