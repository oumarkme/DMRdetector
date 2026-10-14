"""Gene annotation of raw DMR tables (`annotate`) and `get-genes`.

v2.2 read the R package EnsDb.Hsapiens.v86; this version reads an Ensembl GRCh38
GTF (`dmrdetector get-genes` downloads and verifies one). The gene of a DMR is
the best by biotype rank, then largest overlap, then longest gene, then (start,
gene_id) order as in EnsDb; ICR-panel DMRs take the curated panel gene symbol.
"""

from __future__ import annotations

import argparse
import gzip
import os
import re
import shutil
import ssl
import subprocess
import sys
import urllib.request

import numpy as np
import pandas as pd

from dmrdetector.fileio import bed
from dmrdetector.util import rstats


# v2.2 biotype priority (lower = preferred); unlisted biotypes rank 99
_IG_TR_GENES = ("IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene", "TR_C_gene",
                "TR_D_gene", "TR_J_gene", "TR_V_gene")
_PSEUDOGENES = ("pseudogene", "processed_pseudogene", "unprocessed_pseudogene",
                "transcribed_processed_pseudogene",
                "transcribed_unprocessed_pseudogene", "IG_pseudogene",
                "TR_pseudogene")
BIOTYPE_RANK = {
    "protein_coding": 1,
    **dict.fromkeys(_IG_TR_GENES, 2),
    "lncRNA": 3,
    "lincRNA": 3,
    "antisense": 3,
    "miRNA": 4,
    "snoRNA": 5,
    "snRNA": 5,
    "rRNA": 6,
    "Mt_tRNA": 7,
    "Mt_rRNA": 7,
    "misc_RNA": 8,
    "processed_transcript": 9,
    **dict.fromkeys(_PSEUDOGENES, 10),
}
UNRANKED = 99

OUTPUT_COLS = [
    "CHROM", "START", "END", "GENE", "PANEL", "NUM_CPG", "NUM_READS",
    "MEAN_SAMPLE", "MEAN_REF", "MEAN_DIFF", "P_VALUE", "P_FDR", "DIRECTION"
]
OUTPUT_NAMES = [
    "CHROM", "START", "END", "GENE", "PANEL", "NUM_CPG", "NUM_READS",
    "SAMPLE_METH (%)", "REF_METH (%)", "DIFF (%)", "P_VALUE", "P_FDR",
    "DIRECTION"
]

_ATTR = {
    k: re.compile(rf'{k} "([^"]*)"')
    for k in ("gene_id", "gene_name", "gene_biotype")
}


def ensembl_style(name: str) -> str:
    """Returns an Ensembl contig name for a UCSC name.

    chr1-22, X and Y drop `chr` and chrM becomes MT; other names are unchanged
    (v2.2 leaves e.g. chrMT unmatched, so it gets no gene).
    """
    if name[:3] == "chr":
        n = name[3:]
        if n in bed.PRIMARY:
            return n
        if n == "M":
            return "MT"
    return name


def load_genes(gtf: str) -> pd.DataFrame:
    """Returns the gene table of a GTF in EnsDb order.

    Args:
        gtf: Ensembl GRCh38 GTF, plain or .gz.

    Returns:
        One row per `gene` feature (1-based closed coordinates) with symbol and
        biotype rank.
    """
    op = gzip.open if gtf.endswith(".gz") else open
    rows = []
    with op(gtf, "rt") as fh:
        for line in fh:
            if line[0] == "#" or "\tgene\t" not in line:
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            a = f[8]
            gid = _ATTR["gene_id"].search(a)
            gname = _ATTR["gene_name"].search(a)
            gbio = _ATTR["gene_biotype"].search(a)
            gid = gid.group(1) if gid else ""
            gname = gname.group(1) if gname else ""
            rows.append((ensembl_style(f[0]), int(f[3]), int(f[4]), gid, gname,
                         gbio.group(1) if gbio else ""))
    g = pd.DataFrame(
        rows,
        columns=["seq", "start", "end", "gene_id", "gene_name", "biotype"])
    # EnsDb genes() order is (seqname, start, gene_id); v2.2 breaks remaining
    # ties by this hit order
    g = g.sort_values(["seq", "start", "gene_id"],
                      kind="stable").reset_index(drop=True)
    g["symbol"] = np.where(g["gene_name"] != "", g["gene_name"], g["gene_id"])
    g["rank"] = g["biotype"].map(BIOTYPE_RANK).fillna(UNRANKED).astype(int)
    return g


def assign_genes(dmrs: pd.DataFrame, genes: pd.DataFrame) -> np.ndarray:
    """Returns the best gene symbol per DMR (None when no gene overlaps).

    Coordinates are closed, as in GRanges.
    """
    out = np.full(len(dmrs), None, dtype=object)
    qseq = np.array([ensembl_style(c) for c in dmrs["CHROM"]], dtype=object)
    qs = dmrs["START"].to_numpy(np.int64)
    qe = dmrs["END"].to_numpy(np.int64)
    bychr = {
        s: (np.flatnonzero((genes["seq"] == s).to_numpy()),)
        for s in genes["seq"].unique()
    }
    gs, ge = genes["start"].to_numpy(np.int64), genes["end"].to_numpy(np.int64)
    gr, sym = genes["rank"].to_numpy(), genes["symbol"].to_numpy(dtype=object)
    for i in range(len(dmrs)):
        if qseq[i] not in bychr:
            continue
        idx = bychr[qseq[i]][0]
        hit = idx[(gs[idx] <= qe[i]) & (ge[idx] >= qs[i])]
        if hit.size == 0:
            continue
        ov = np.minimum(qe[i], ge[hit]) - np.maximum(qs[i], gs[hit]) + 1
        ln = ge[hit] - gs[hit] + 1
        # subject order breaks remaining ties
        best = np.lexsort((hit, -ln, -ov, gr[hit]))[0]
        out[i] = sym[hit[best]]
    return out


def override_icr_genes(dmrs: pd.DataFrame, gene_col: np.ndarray, icr_panel: str,
                       icr_gene_list: str | None) -> None:
    """Sets the panel gene symbol of the ICR rows of `dmrs` in `gene_col`.

    The symbol is that of the largest-overlap ICR panel entry; ties go to the
    smallest gene-list span, then the name.
    """
    icr_rows = np.flatnonzero((dmrs["PANEL"] == "ICR").to_numpy())
    if icr_rows.size == 0 or not icr_panel or not os.path.exists(icr_panel):
        return
    pan = pd.read_csv(icr_panel, sep="\t", dtype={0: str})
    if len(pan) == 0 or "GENE" not in pan.columns:
        return
    pan = pd.DataFrame({
        "CHROM": bed.normalize_primary_chrom(pan["CHROM"]).to_numpy(),
        "START": pan["START"].astype(np.int64),
        "END": pan["END"].astype(np.int64),
        "PANEL_GENE": pan["GENE"].astype(str)
    })
    span = None
    if icr_gene_list and os.path.exists(icr_gene_list):
        gl = pd.read_csv(icr_gene_list, sep="\t", dtype={0: str}).iloc[:, :4]
        gl.columns = ["CHROM", "START", "END", "GENE"]
        gl["w"] = gl["END"].astype(np.int64) - gl["START"].astype(np.int64)
        span = gl.groupby("GENE")["w"].min().to_dict()
    big = np.iinfo(np.int32).max
    for r in icr_rows:
        c, s, e = dmrs["CHROM"].iat[r], int(dmrs["START"].iat[r]), int(
            dmrs["END"].iat[r])
        p = pan[(pan["CHROM"] == c) & (pan["START"] <= e) & (pan["END"] >= s)]
        if len(p) == 0:
            continue
        ovw = np.minimum(e, p["END"].to_numpy()) - np.maximum(
            s, p["START"].to_numpy()) + 1
        names = p["PANEL_GENE"].tolist()
        ms = [(span.get(n, big) if span is not None else 0) for n in names]
        best = min(range(len(names)),
                   key=lambda k: (-int(ovw[k]), ms[k], names[k]))
        gene_col[r] = names[best]


def annotate(raw_files: list[str],
             gtf: str,
             icr_panel: str | None = None,
             icr_gene_list: str | None = None,
             *,
             genes: pd.DataFrame | None = None) -> pd.DataFrame:
    """Returns the annotated DMR table with the v2.2 output columns and names.

    Values are still unrounded decimals, as in the v2.2 output.

    Args:
        raw_files: Raw DMR tables (all scopes of one sample).
        gtf: Ensembl GTF; not opened when `genes` is given.
        icr_panel: ICR panel BED for the ICR gene override.
        icr_gene_list: Imprinting gene BED for the tie-break by span.
        genes: Gene table from `load_genes`; pass it when annotating many
            samples so the GTF is read once.
    """
    frames = [
        pd.read_csv(f,
                    sep="\t",
                    dtype={"CHROM": str},
                    float_precision="round_trip") for f in raw_files
    ]
    dmrs = pd.concat(frames, ignore_index=True)
    if len(dmrs) == 0:
        return pd.DataFrame({c: pd.Series(dtype=object) for c in OUTPUT_NAMES})
    dmrs["CHROM"] = bed.normalize_primary_chrom(dmrs["CHROM"]).to_numpy()
    gene = assign_genes(dmrs, genes if genes is not None else load_genes(gtf))
    override_icr_genes(dmrs, gene, icr_panel, icr_gene_list)
    dmrs["GENE"] = gene
    dmrs = dmrs[OUTPUT_COLS].copy()
    dmrs["NUM_READS"] = rstats.r_round(dmrs["NUM_READS"].to_numpy(float))
    for c in ("MEAN_SAMPLE", "MEAN_REF", "MEAN_DIFF"):
        dmrs[c] = rstats.r_round(dmrs[c].to_numpy(float) * 100, 1)
    dmrs.columns = OUTPUT_NAMES
    return dmrs.sort_values(["CHROM", "START"],
                            kind="stable").reset_index(drop=True)


def write_annotated(df: pd.DataFrame, path: str) -> None:
    """Writes an annotated table as tab-separated text, as fwrite does."""
    bed.fwrite_table(df,
                     path,
                     int_cols=("START", "END", "NUM_CPG"),
                     float_cols=("NUM_READS", "SAMPLE_METH (%)", "REF_METH (%)",
                                 "DIFF (%)", "P_VALUE", "P_FDR"))


ENSEMBL_FTP = "https://ftp.ensembl.org/pub"


def _bsd_sum_py(path: str) -> tuple[int, int]:
    """BSD `sum` in Python (exact, slow: about a minute for 141 MB).

    The rotate-right and the add modulo 2**16 do not combine into a linear map,
    so a vectorized or prefix-scan form is not exact.
    """
    s = n = 0
    with open(path, "rb") as fh:
        while chunk := fh.read(1 << 20):
            for b in chunk:
                s = (((s >> 1) + ((s & 1) << 15)) + b) & 0xFFFF
            n += len(chunk)
    return s, -(-n // 1024)


def bsd_sum(path: str) -> tuple[int, int]:
    """Returns (checksum, 1 KiB block count) as `sum` prints for CHECKSUMS.

    It uses `sum -r` (coreutils, a second for a GTF) when available, else the
    Python loop.
    """
    exe = shutil.which("sum")
    if exe:
        try:
            r = subprocess.run([exe, "-r", path],
                               capture_output=True,
                               text=True,
                               timeout=600)
            f = r.stdout.split()
            if r.returncode == 0 and len(f) >= 2:
                return int(f[0]), int(f[1])
        except (OSError, ValueError, subprocess.SubprocessError):
            pass
    return _bsd_sum_py(path)


def _ssl_fallback() -> None:
    """Point OpenSSL at a system CA bundle when the Python default has none."""
    cafile = ssl.get_default_verify_paths().cafile
    if os.environ.get("SSL_CERT_FILE") or (cafile and os.path.exists(cafile)):
        return
    for f in ("/etc/ssl/certs/ca-certificates.crt",
              "/etc/pki/tls/certs/ca-bundle.crt", "/etc/ssl/cert.pem"):
        if os.path.exists(f):
            os.environ["SSL_CERT_FILE"] = f
            return


def _get(url: str) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": "dmrdetector"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.read().decode("utf-8", "replace")


def _download(url: str, dest: str) -> None:
    req = urllib.request.Request(url, headers={"User-Agent": "dmrdetector"})
    with urllib.request.urlopen(req, timeout=120) as r, open(dest, "wb") as fh:
        while chunk := r.read(1 << 20):
            fh.write(chunk)


def parse_checksums(text: str, name: str) -> tuple[int, int] | None:
    """Returns (checksum, blocks) of `name` in an Ensembl CHECKSUMS text.

    Each line is `<sum> <blocks> <file>`; the result is None when `name` is not
    listed.
    """
    expect = None
    for line in text.splitlines():
        f = line.split()
        if len(f) == 3 and f[2] == name:
            expect = (int(f[0]), int(f[1]))
    return expect


def latest_release() -> int:
    """Returns the newest Ensembl release whose human GTF is in CHECKSUMS.

    It tries the newest three releases, since a new release folder can be
    half-populated.

    Raises:
        SystemExit: The server is unreachable or no release qualifies.
    """
    _ssl_fallback()
    try:
        nums = sorted(
            {
                int(n) for n in re.findall(r'href="release-(\d+)/"',
                                           _get(f"{ENSEMBL_FTP}/"))
            },
            reverse=True)
    except OSError as e:
        raise SystemExit(f"cannot reach {ENSEMBL_FTP} ({e}); retry later or "
                         "give a release number") from None
    for n in nums[:3]:
        try:
            ck = _get(f"{ENSEMBL_FTP}/release-{n}/gtf/homo_sapiens/CHECKSUMS")
        except OSError:
            continue
        if parse_checksums(ck, f"Homo_sapiens.GRCh38.{n}.gtf.gz"):
            return n
    raise SystemExit(
        "could not resolve the latest Ensembl release; give a release number "
        "with --release N")


def get_genes(release: int | str, outdir: str) -> str:
    """Downloads an Ensembl GRCh38 GTF and verifies it against CHECKSUMS.

    Args:
        release: Release number or `latest`.
        outdir: Folder for the GTF.

    Returns:
        The path of the verified GTF.

    Raises:
        SystemExit: The release is not listed or the checksum differs.
    """
    _ssl_fallback()
    release = latest_release() if str(release) == "latest" else int(release)
    base = f"{ENSEMBL_FTP}/release-{release}/gtf/homo_sapiens"
    name = f"Homo_sapiens.GRCh38.{release}.gtf.gz"
    os.makedirs(outdir, exist_ok=True)
    dest = os.path.join(outdir, name)
    ck_path = os.path.join(outdir, f"CHECKSUMS.release-{release}")
    _download(f"{base}/CHECKSUMS", ck_path)
    with open(ck_path) as fh:
        expect = parse_checksums(fh.read(), name)
    if expect is None:
        raise SystemExit(f"{name} not listed in {base}/CHECKSUMS")
    if not (os.path.exists(dest) and bsd_sum(dest) == expect):
        part = dest + ".part"
        _download(f"{base}/{name}", part)
        got = bsd_sum(part)
        if got != expect:
            os.remove(part)
            raise SystemExit(
                f"checksum mismatch for {name}: expected {expect}, got {got}")
        os.replace(part, dest)
    print(
        f"Ensembl release {release}: {dest} (checksum {expect[0]} blocks "
        f"{expect[1]} verified)",
        file=sys.stderr)
    return dest


def _release(v: str) -> str:
    if v != "latest" and not re.fullmatch(r"\d+", v):
        raise argparse.ArgumentTypeError("must be a release number or 'latest'")
    return v


def add_parser(sub: argparse._SubParsersAction) -> None:
    """Adds the `annotate` and `get-genes` commands to `sub`."""
    p = sub.add_parser("annotate",
                       help="ANNOTATE_DMRS: add gene symbols to raw DMR tables",
                       description="ANNOTATE_DMRS: gene annotation of raw DMR "
                       "tables, plus `get-genes`.")
    p.add_argument(
        "--raw",
        required=True,
        nargs="+",
        help=
        "dmr_raw_<id>_<type>_<scope>.tsv.gz files of one sample (all scopes)")
    p.add_argument("--gtf",
                   required=True,
                   help="Ensembl GRCh38 GTF (see get-genes)")
    p.add_argument("--icr-panel",
                   help="ICR panel BED (GENE column) for the ICR gene override")
    p.add_argument(
        "--icr-gene-list",
        help="imprinting gene BED (CHROM START END GENE), tie-break by span")
    p.add_argument("--out", required=True, help="output dmr_<id>_<type>.tsv")
    p.set_defaults(func=main)

    g = sub.add_parser(
        "get-genes",
        help="download and verify an Ensembl GRCh38 GTF",
        description="Download an Ensembl GRCh38 GTF (release N or `latest`) "
        "and verify it against the release CHECKSUMS.")
    g.add_argument(
        "--release",
        required=True,
        type=_release,
        help=
        "Ensembl release number (EnsDb.Hsapiens.v86 = 86) or `latest` (newest "
        "release with the GTF "
        "in its CHECKSUMS)")
    g.add_argument("--outdir",
                   required=True,
                   help="directory for the GTF (no default)")
    g.set_defaults(func=main_get_genes)


def main(args: argparse.Namespace) -> int:
    """Runs `annotate` on the raw tables given with --raw."""
    files = sorted(args.raw)  # v2.2 list.files() order
    df = annotate(files, args.gtf, args.icr_panel, args.icr_gene_list)
    write_annotated(df, args.out)
    print(f"Wrote {len(df)} annotated DMRs to {args.out}", file=sys.stderr)
    return 0


def main_get_genes(args: argparse.Namespace) -> int:
    # stdout: path of the verified GTF
    """Runs `get-genes` and prints the path of the verified GTF."""
    try:
        print(get_genes(args.release, args.outdir))
    except SystemExit as e:
        print(f"ERROR: get-genes: {e.code}", file=sys.stderr)
        return 1
    except OSError as e:
        print(f"ERROR: get-genes: cannot download from {ENSEMBL_FTP} ({e}); "
              "check the release number and the network (compute nodes often "
              "have no internet access: run it on a login node, or download "
              "the GTF elsewhere and pass it with --gtf)",
              file=sys.stderr)
        return 1
    return 0
