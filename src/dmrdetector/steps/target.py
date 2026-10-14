"""Target-region report: methylation of disease-critical regions.

For each sample and autosomal target region the CpGs of the sample are joined to
the autosomal reference atlas and summarized for the merged track and, with hap1
and hap2 BEDs, for each haplotype. DIFF, P_FDR and DIRECTION are defined in
doc/target_region_methods.md and doc/ALGORITHM.md.
"""

from __future__ import annotations

import argparse
from collections.abc import Callable
from collections.abc import Iterable
import csv
import hashlib
import io
import os
import re
import tempfile

import numpy as np
import pandas as pd

from dmrdetector import data
from dmrdetector import runtime
from dmrdetector import settings
from dmrdetector.fileio import bed
from dmrdetector.fileio import methbed
from dmrdetector.fileio import vcf as vcf_io
from dmrdetector.steps import call_dmr
from dmrdetector.steps import preprocess
from dmrdetector.util import rstats


def _no_jit() -> bool:
    return os.environ.get("DMRDETECTOR_NO_JIT") == "1"


DEFAULT_ICR_GENES = data.asset("icr_gene_list.bed")
TRACKS = ("merged", "hap1", "hap2")
OUT_COLS = [
    "SAMPLE_ID", "SAMPLE_TYPE", "CHROM", "START", "END", "REGION_NAME",
    "GENE_SYMBOL", "DISEASE_NAME", "HAPLOTYPE", "NUM_CPG", "SAMPLE_METH (%)",
    "REF_METH (%)", "REF_SD (%)", "COVERAGE", "DIFF", "P_FDR", "DIRECTION"
]
_FLOAT_OUT = ("SAMPLE_METH (%)", "REF_METH (%)", "REF_SD (%)", "COVERAGE",
              "DIFF", "P_FDR")


def read_target_bed(path: str,
                    log: Callable[[str], None] = print) -> pd.DataFrame:
    """Reads the target BED.

    Header and comment lines are ignored (a row needs an integer START).
    Non-autosomal regions are dropped with a message.

    Returns:
        CHROM, START, END, REGION_NAME, GENE_SYMBOL, DISEASE_NAME and REGION_ID.
    """
    rows = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\r\n").split("\t")
            if len(f) < 6 and line.strip() and re.fullmatch(
                    r"\d+", f[1] if len(f) > 1 else ""):
                raise ValueError(
                    "target BED needs 6 columns: CHROM START END REGION_NAME "
                    "GENE_SYMBOL DISEASE_NAME")
            if len(f) < 6 or f[0].startswith("#") or not re.fullmatch(
                    r"\d+", f[1]):
                continue
            rows.append(f[:6])
    tgt = pd.DataFrame(rows,
                       columns=[
                           "CHROM", "START", "END", "REGION_NAME",
                           "GENE_SYMBOL", "DISEASE_NAME"
                       ])
    tgt["CHROM"] = bed.normalize_primary_chrom(tgt["CHROM"]).to_numpy()
    tgt[["START", "END"]] = tgt[["START", "END"]].astype(np.int64)
    skip = ~tgt["CHROM"].isin(bed.AUTOSOMES)
    if skip.any():
        log("Skipping non-autosomal target regions (autosomal atlas only): " +
            ", ".join(tgt.loc[skip, "REGION_NAME"]))
    tgt = tgt[~skip].reset_index(drop=True)
    tgt["REGION_ID"] = np.arange(len(tgt))
    return tgt


def _read_intervals(path: str | None) -> pd.DataFrame:
    """Returns the first three columns of a BED-like file as CHROM, START, END.

    Header and comment lines are skipped.
    """
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return pd.DataFrame({"CHROM": [], "START": [], "END": []})
    df = pd.read_csv(path,
                     sep="\t",
                     header=None,
                     usecols=[0, 1, 2],
                     dtype=str,
                     comment="#")
    df.columns = ["CHROM", "START", "END"]
    df = df[df["START"].str.fullmatch(r"\d+", na=False)].copy()
    df["CHROM"] = bed.normalize_primary_chrom(df["CHROM"]).to_numpy()
    df[["START", "END"]] = df[["START", "END"]].astype(np.int64)
    return df


def icr_overlap(tgt: pd.DataFrame, *panels: str | None) -> np.ndarray:
    """Returns True where a target region overlaps (half-open) a panel interval.

    The intervals come from the given BED files.
    """
    hit = np.zeros(len(tgt), dtype=bool)
    for p in panels:
        iv = _read_intervals(p)
        for chrom, grp in iv.groupby("CHROM"):
            m = (tgt["CHROM"] == chrom).to_numpy()
            if not m.any():
                continue
            s, e = tgt["START"].to_numpy()[m], tgt["END"].to_numpy()[m]
            gs, ge = grp["START"].to_numpy(), grp["END"].to_numpy()
            hit[m] |= ((s[:, None] < ge[None, :]) &
                       (gs[None, :] < e[:, None])).any(axis=1)
    return hit


def _cov_cols(nf: int, fmt: str = "auto") -> tuple[int, int]:
    """Returns the 0-based (coverage, modified reads) columns of a layout.

    With `fmt` auto the layout is given by its column count (the v2 cov_cols
    rule); a forced format ("5col", "modkit" or "pbcpg") fixes the columns.
    """
    if fmt != "auto":
        return {"modkit": (9, 11), "pbcpg": (5, 6), "5col": (3, 4)}[fmt]
    if nf >= 18:
        return 9, 11
    if 7 <= nf <= 9:
        return 5, 6
    return 3, 4


def _regions_from_arrays(names, bounds, start, end, regions):
    """Finds the rows overlapping each region [START, END), as tabix fetch.

    The regions are taken one after another in file order. `names` and `bounds`
    are the chromosome runs of the arrays (contig name as in the file).

    Returns:
        [(REGION_ID, name, row index)].
    """
    by_norm = dict(zip(bed.normalize_primary_chrom(names), names))
    run_of = {n: i for i, n in enumerate(names)}
    maxlen = {}
    out = []
    for rid, chrom, s, e in zip(regions["REGION_ID"], regions["CHROM"],
                                regions["START"], regions["END"]):
        name = by_norm.get(chrom)
        if name is None:
            continue
        i = run_of[name]
        a, b = int(bounds[i]), int(bounds[i + 1])
        st, en = start[a:b], end[a:b]
        if i not in maxlen:
            maxlen[i] = int((en - st).max()) if b > a else 0
        lo = a + int(
            np.searchsorted(st, int(s) - max(maxlen[i], 0), side="left"))
        hi = a + int(np.searchsorted(st, int(e), side="left"))
        if hi > lo:
            ix = np.flatnonzero(end[lo:hi] > int(s)) + lo
            if ix.size:
                out.append((int(rid), name, ix))
    return out


def _load_whole(path: str, cols, threads: int):
    """Reads a BGZF bed in memory.

    `cols` is the ("s", "i" or "f") spec.

    Returns:
        ({column: array}, names, bounds), or None.
    """
    spec = [(c, k) for c, k, _ in cols]
    # numba is imported only when a file is read here
    from dmrdetector.fileio import bgzf
    from dmrdetector.fileio import numparse
    try:
        buf = None if bgzf.NO_JIT else bgzf.inflate(path, threads)
        if buf is None:
            return None
        r = numparse.parse_buffer(buf[_header_len(buf):], spec, threads)
    except bgzf.ParseError:
        return None
    ints, floats, (rs, names) = r
    it, ft = iter(ints), iter(floats)
    arrs = {
        nm: (next(it) if k == "i" else next(ft)) for _, k, nm in sorted(cols)
    }
    n = next(iter(arrs.values())).size
    if len(set(names)) != len(names):
        return None
    bounds = np.append(rs, n)
    for i in range(len(names)):
        a, b = int(bounds[i]), int(bounds[i + 1])
        if b - a > 1 and (np.diff(arrs["START"][a:b]) < 0).any():
            return None
    return arrs, names, bounds


def _is_header(line: bytes | str) -> bool:
    """Tells a header line: `#...`, `chrom...` or blank (as scan_header)."""
    if isinstance(line, bytes):
        return line.startswith((b"#", b"chrom")) or not line.strip()
    return line.startswith(("#", "chrom")) or not line.strip()


def _header_len(buf: np.ndarray) -> int:
    """Returns the byte length of the leading header lines of a table."""
    pos, n = 0, buf.size
    while pos < n:
        end = pos
        while end < n and buf[end] != 10:
            end += 1
        if not _is_header(bytes(buf[pos:end])):
            break
        pos = end + 1
    return min(pos, n)


def _head(path: str) -> tuple[int, str | None]:
    """Returns the number of leading header lines and the first data line."""
    nhead = 0
    with bed.open_text(path) as fh:
        for ln in fh:
            if _is_header(ln):
                nhead += 1
                continue
            return nhead, ln.rstrip("\n")
    return nhead, None


def _first_line_fields(path: str) -> int | None:
    """Returns the columns (split on whitespace) of the first data line.

    None means no data line, or a data line that is not tab-separated.
    """
    line = _head(path)[1]
    if line is None or "\t".join(line.split()) != line:
        return None
    return len(line.split())


def _is_bgzf(path: str) -> bool:
    """Returns True when the file starts with a BGZF block header."""
    with open(path, "rb") as fh:
        h = fh.read(18)
    return h[:4] == b"\x1f\x8b\x08\x04" and h[12:14] == b"BC"


def _tab_chunks(path: str, nhead: int, rows: int = 100_000):
    """Yields the data lines of a table as tab-separated text chunks."""
    buf = []
    with bed.open_text(path) as fh:
        for i, ln in enumerate(fh):
            if i < nhead:
                continue
            f = ln.split()
            if f:
                buf.append("\t".join(f))
            if len(buf) >= rows:
                yield "\n".join(buf) + "\n"
                buf = []
    if buf:
        yield "\n".join(buf) + "\n"


def _index_is_fresh(path: str) -> bool:
    """Returns True when `path`.tbi exists and is not older than the data."""
    return (os.path.exists(path + ".tbi") and
            os.path.getmtime(path + ".tbi") >= os.path.getmtime(path))


def _indexable(path: str, scratch: str) -> str | None:
    """Returns a path that tabix can index and read, or None without data.

    An indexed file is used as it is, unless its index is older than the data:
    that index is not trusted and a new one is made in `scratch`. An unindexed
    BGZF file of tab-separated data lines is indexed through a symlink. Any
    other file (header lines, spaces between columns as in modkit, gzip that is
    not BGZF, plain text) is rewritten as a BGZF copy without the header lines
    and with tabs. Each input gets a folder of its own in `scratch`, so files
    with one name in different folders do not share an index.
    """
    if _index_is_fresh(path):
        return path
    sub = os.path.join(
        scratch,
        hashlib.sha1(os.path.abspath(path).encode()).hexdigest()[:16])
    nhead, first = _head(path)
    if first is None:
        return None
    try:
        if nhead == 0 and _is_bgzf(path) and _first_line_fields(path):
            return bed.ensure_tabix(path, scratch, check_age=True)
        os.makedirs(sub, exist_ok=True)
        clean = os.path.join(sub, "track.bed.gz")
        bed.write_bgzip_text(clean, _tab_chunks(path, nhead))
        bed.tabix_bed(clean)
    except OSError as e:
        raise runtime.RunError(
            f"{path}: cannot index the track ({e}); it must be sorted by "
            "contig and START, with each contig in one block (for example "
            "`sort -k1,1 -k2,2n`)") from None
    return clean


def _check_indexed_hg38(tbx, path: str) -> None:
    """Stops when an indexed track has a record beyond the hg38 contig length.

    It reads the records that end after the length plus the margin of
    `preprocess.check_hg38`, with one seek per primary contig, so the test is
    the one made on the whole file by `preprocess`.
    """
    beyond = {}
    for name in tbx.contigs:
        m = preprocess.PRIMARY_CONTIG.match(name)
        length = preprocess.HG38_LENGTHS.get("chr" + m.group(2)) if m else None
        if length is None:
            continue
        ends = [
            int(ln.split()[2])
            for ln in tbx.fetch(name, length + preprocess.END_MARGIN)
            if not ln.startswith(("#", "chrom"))
        ]
        if ends:
            beyond[name] = max(ends)
    preprocess.check_hg38(beyond, path)


def _fetch_regions(
        path: str,
        regions: pd.DataFrame,
        scratch: str,
        check_build: bool = False) -> list[tuple[int, str, list[str]]]:
    """Reads the target regions with tabix.

    The contig may be named with or without `chr` in the file. With
    `check_build` the file is first checked for records beyond the hg38 contig
    lengths (`_check_indexed_hg38`).

    Returns:
        [(REGION_ID, contig name in the file, lines)].
    """
    import pysam
    path = _indexable(path, scratch)
    if path is None:
        return []
    out = []
    with pysam.TabixFile(path) as tbx:
        if check_build:
            _check_indexed_hg38(tbx, path)
        names = list(tbx.contigs)
        by_norm = dict(zip(bed.normalize_primary_chrom(names), names))
        for rid, chrom, s, e in zip(regions["REGION_ID"], regions["CHROM"],
                                    regions["START"], regions["END"]):
            name = by_norm.get(chrom)
            if name is None:
                continue
            lines = [
                ln for ln in tbx.fetch(name, int(s), int(e))
                if not ln.startswith(("#", "chrom"))
            ]
            if lines and "\t".join(lines[0].split()) != lines[0]:
                lines = ["\t".join(ln.split()) for ln in lines]
            if lines:
                out.append((int(rid), name, lines))
    return out


_TRACK_COLS = [
    "REGION_ID", "RAW_CHROM", "CHROM", "START", "END", "COVERAGE", "N_METH"
]


def _empty_track() -> pd.DataFrame:
    """Returns a table of the track columns without rows."""
    return pd.DataFrame({
        c: pd.Series(dtype=object if c in ("RAW_CHROM", "CHROM") else float)
        for c in _TRACK_COLS
    })


def _read_track_fast(path: str,
                     tgt: pd.DataFrame,
                     threads: int,
                     input_format: str = "auto"):
    """Returns the CpGs of a track in the target regions, or None.

    The track is an unindexed BGZF bed, inflated and parsed in memory; the rows
    and their order equal the tabix region reads. None means the fast reader
    does not apply.

    Raises:
        runtime.RunError: A contig ends beyond its hg38 length.
    """
    if _no_jit() or bed.is_gzip(path) is False:
        return None
    nf = _first_line_fields(path)
    if nf is None:
        return None
    ci, mi = _cov_cols(nf, input_format)
    got = _load_whole(path, [(1, "i", "START"), (2, "i", "END"),
                             (ci, "f", "COVERAGE"), (mi, "f", "N_METH")],
                      threads)
    if got is None:
        return None
    arrs, names, bounds = got
    preprocess.check_hg38(
        {
            names[i]: int(arrs["END"][bounds[i]:bounds[i + 1]].max())
            for i in range(len(names))
            if bounds[i + 1] > bounds[i]
        }, path)
    hits = _regions_from_arrays(names, bounds, arrs["START"], arrs["END"], tgt)
    if not hits:
        return _empty_track()
    ix = np.concatenate([h[2] for h in hits])
    rid = np.concatenate([np.full(h[2].size, h[0]) for h in hits])
    raw = np.concatenate([np.full(h[2].size, h[1], dtype=object) for h in hits])
    df = pd.DataFrame({
        "REGION_ID": rid,
        "RAW_CHROM": raw,
        "CHROM": bed.normalize_primary_chrom(raw).to_numpy(),
        "START": arrs["START"][ix],
        "END": arrs["END"][ix],
        "COVERAGE": arrs["COVERAGE"][ix],
        "N_METH": arrs["N_METH"][ix]
    })
    return df[_TRACK_COLS]


def _check_track_format(path: str, input_format: str) -> str:
    """Detects the format of a track and stops on a file `preprocess` rejects.

    Haplotype tracks and tracks of the stand-alone `target` command get the
    checks that `run` makes on the merged track: a modkit file with h or a rows,
    or with one row per strand, is refused.

    Returns:
        The format: 5col, modkit or pbcpg.

    Raises:
        runtime.RunError: The forced format does not fit the file, or the file
            would count a CpG twice.
    """
    try:
        fmt = methbed.scan_header(path, input_format)[0]
    except ValueError as e:
        raise runtime.RunError(f"{path}: cannot read as input_format "
                               f"{input_format} ({e}); check input_format or "
                               "use auto") from None
    if fmt == "modkit":
        preprocess.check_modkit(path)
    return fmt


def _read_track_sites(path: str,
                      tgt: pd.DataFrame,
                      scratch: str,
                      threads: int = 4,
                      input_format: str = "auto") -> pd.DataFrame:
    """Returns the CpGs of one track inside the target regions.

    Columns: REGION_ID, RAW_CHROM, CHROM, START, END, COVERAGE, N_METH. The
    track is checked as `preprocess` checks a sample (`_check_track_format`,
    `preprocess.check_hg38`).

    Raises:
        runtime.RunError: The format or the build check fails.
    """
    _check_track_format(path, input_format)
    if not _index_is_fresh(path):  # no index, or one that is not trusted
        fast = _read_track_fast(path, tgt, threads, input_format)
        if fast is not None:
            return fast
    parts = _fetch_regions(path, tgt, scratch, check_build=True)
    if not parts:
        return _empty_track()
    ci, mi = _cov_cols(len(parts[0][2][0].split("\t")), input_format)
    frames = []
    for rid, name, lines in parts:
        f = pd.read_csv(io.StringIO("\n".join(lines)),
                        sep="\t",
                        header=None,
                        usecols=[1, 2, ci, mi],
                        dtype=np.float64,
                        float_precision="round_trip")
        f.columns = ["START", "END", "COVERAGE", "N_METH"]
        f[["START", "END"]] = f[["START", "END"]].astype(np.int64)
        f.insert(0, "REGION_ID", rid)
        f.insert(1, "RAW_CHROM", name)
        frames.append(f)
    df = pd.concat(frames, ignore_index=True)
    df["CHROM"] = bed.normalize_primary_chrom(df["RAW_CHROM"]).to_numpy()
    return df[_TRACK_COLS]


class _HetSites:
    """Heterozygous SNP positions of a VCF by exact contig name.

    Positions are read only where needed. With an indexed VCF (.tbi or .csi)
    only the windows of the target regions (+-2 bp) are fetched; the positions
    equal those of the whole-VCF scan `vcf.het_snp_sites` for every CpG inside a
    target region. A VCF without an index is scanned once. Sites computed
    already can be given as `sites`.
    """

    def __init__(self, vcf=None, tgt=None, sites=None):
        self.vcf, self.tgt, self.sites = vcf, tgt, sites
        self.vf = None
        self.cache = {}

    def get(self, contig):
        if self.sites is not None:
            return self.sites.get(contig)
        if self.vcf is None:
            return None
        if self.vf is None:
            self.vf = False
            if os.path.exists(self.vcf + ".tbi") or os.path.exists(self.vcf +
                                                                   ".csi"):
                try:
                    import pysam
                    self.vf = pysam.VariantFile(self.vcf)
                except (OSError, ValueError):
                    self.vf = False
            if self.vf is False:
                self.sites = vcf_io.het_snp_sites(self.vcf)
                return self.sites.get(contig)
        return self._fetch(contig)

    def _fetch(self, contig):
        if contig in self.cache:
            return self.cache[contig]
        vf = self.vf
        norm = bed.normalize_primary_chrom([contig])[0]
        m = (self.tgt["CHROM"] == norm).to_numpy()
        pos = []
        for s, e in sorted(
                zip(self.tgt["START"].to_numpy()[m].tolist(),
                    self.tgt["END"].to_numpy()[m].tolist())):
            try:
                recs = list(vf.fetch(contig, max(0, s - 2), e + 2))
            except (ValueError, OSError):
                # contig not in the index: no records, as in the full scan
                self.cache[contig] = None
                return None
            for rec in recs:
                het = vcf_io.het_position(rec)
                if het is not None:
                    pos.append(het)
        out = np.unique(np.asarray(pos, dtype=np.int64)) if pos else None
        self.cache[contig] = out
        return out

    def close(self):
        if self.vf:
            self.vf.close()


def _drop_het(df: pd.DataFrame, het) -> pd.DataFrame:
    """Removes the CpGs whose [START, END) contains a heterozygous SNP.

    This follows v2 and preprocess; contigs compare as exact strings. `het` is a
    dict contig -> sorted positions, or a _HetSites.
    """
    if het is None or df.empty:
        return df
    if isinstance(het, dict):
        if not het:
            return df
        het = _HetSites(sites=het)
    drop = np.zeros(len(df), dtype=bool)
    raw = df["RAW_CHROM"].to_numpy()
    start, end = df["START"].to_numpy(), df["END"].to_numpy()
    for c in pd.unique(raw):
        p0 = het.get(c)
        if p0 is not None and len(p0):
            m = raw == c
            drop[m] = bed.sorted_point_overlap(start[m], end[m], p0)
    return df[~drop]


def _summarise(sites: pd.DataFrame, atlas: pd.DataFrame, min_cov: float, het,
               tgt: pd.DataFrame, want_diff: bool):
    """Summarizes one track per region.

    Returns:
        (DataFrame by REGION_ID, {REGION_ID: per-CpG DIFF array}).
    """
    s = _drop_het(sites[sites["COVERAGE"] > min_cov], het)
    s = s.assign(SAMPLE_METH=s["N_METH"] / s["COVERAGE"])
    # like-for-like CpG set
    ov = s.merge(atlas, on=["CHROM", "START"], how="inner")
    ov["DIFF"] = ov["SAMPLE_METH"] - ov["REF_METH"]
    g = ov.groupby("REGION_ID", sort=True)
    res = pd.DataFrame({
        "NUM_CPG": g.size(),
        "SAMPLE_METH": g["SAMPLE_METH"].mean(),
        "REF_METH": g["REF_METH"].mean(),
        "REF_SD": np.sqrt(g["REF_VAR"].mean()),
        "COVERAGE": g["COVERAGE"].mean(),
        "DIFF": g["DIFF"].mean()
    })
    res = res.reindex(tgt["REGION_ID"])
    res["NUM_CPG"] = res["NUM_CPG"].fillna(0).astype(np.int64)
    diffs = {}
    if want_diff:
        diffs = {int(r): grp["DIFF"].to_numpy() for r, grp in g}
    return res, diffs


def _atlas_over_regions(atlas: str, tgt: pd.DataFrame, scratch: str,
                        threads: int) -> pd.DataFrame:
    """Returns REF_METH (column 8) and REF_VAR (column 9) of the atlas CpGs.

    The CpGs are those in the target regions; columns are CHROM, START, REF_METH
    and REF_VAR.
    """
    if not _index_is_fresh(atlas) and not _no_jit() and bed.is_gzip(atlas):
        got = _load_whole(atlas, [(1, "i", "START"), (2, "i", "END"),
                                  (7, "f", "REF_METH"), (8, "f", "REF_VAR")],
                          threads)
        if got is not None:
            arrs, names, bounds = got
            hits = _regions_from_arrays(names, bounds, arrs["START"],
                                        arrs["END"], tgt)
            frames = []
            for _rid, name, ix in hits:
                f = pd.DataFrame({
                    "START": arrs["START"][ix],
                    "REF_METH": arrs["REF_METH"][ix],
                    "REF_VAR": arrs["REF_VAR"][ix]
                })
                f["CHROM"] = bed.normalize_primary_chrom([name] *
                                                         len(f)).to_numpy()
                frames.append(f)
            return _unique_atlas(frames)
    frames = []
    for _rid, name, lines in _fetch_regions(atlas, tgt, scratch):
        f = pd.read_csv(io.StringIO("\n".join(lines)),
                        sep="\t",
                        header=None,
                        usecols=[1, 7, 8],
                        dtype=np.float64,
                        float_precision="round_trip")
        f.columns = ["START", "REF_METH", "REF_VAR"]
        f["START"] = f["START"].astype(np.int64)
        f["CHROM"] = bed.normalize_primary_chrom([name] * len(f)).to_numpy()
        frames.append(f)
    return _unique_atlas(frames)


def _unique_atlas(frames):
    """Joins the per-region atlas frames and drops repeated CpGs.

    With no frame the result is empty but typed like the joined one, so that it
    can still be merged with the CpGs of a sample.
    """
    return (pd.concat(frames, ignore_index=True).drop_duplicates(
        ["CHROM", "START"]) if frames else pd.DataFrame({
            "CHROM": pd.Series(dtype=object),
            "START": pd.Series(dtype=np.int64),
            "REF_METH": pd.Series(dtype=np.float64),
            "REF_VAR": pd.Series(dtype=np.float64)
        }))


def target_report(
        samples: Iterable[dict],
        atlas: str,
        target_bed: str,
        out: str | None = None,
        *,
        icr_region: str | None = None,
        icr_gene_list: str | None = None,
        haplotypes: bool = True,
        min_coverage: float = settings.TARGET_MIN_COVERAGE,
        hap_min_coverage: float = settings.TARGET_HAP_MIN_COVERAGE,
        fdr_threshold: float = call_dmr.DEFAULTS["fdr_threshold"],
        cpg_diff_threshold: float = call_dmr.DEFAULTS["cpg_diff_threshold"],
        icr_diff_threshold: float = call_dmr.DEFAULTS["icr_diff_threshold"],
        input_format: str = settings.INPUT_FORMAT,
        scratch: str | None = None,
        log: Callable[[str], None] = print,
        threads: int = 4) -> pd.DataFrame:
    """Builds the target-region table for all samples.

    `out` (target_region.tsv) is written when given. Tracks and atlas with a
    .tbi are read by region; BGZF files without one are read in memory, other
    files get an index in `scratch`. An indexed VCF is read in the target
    regions only; an unindexed one is scanned once.

    Args:
        samples: Iterable of dicts with sample_id, type ("CASE" or "CTRL"),
            merged (BED.gz; 5-column, modkit or pb-CpG-tools), optional hap1 and
            hap2 (BED.gz), vcf (for het-SNP removal; required whenever hap1 or
            hap2 are given, optional for a preprocessed merged BED) and het
            (heterozygous SNP sites already computed, the dict of
            `vcf.het_snp_sites(vcf)`; replaces the VCF read).
        atlas: Autosomal narrow-bandwidth atlas BED.gz (9 columns);
            tabix-indexed
            (else indexed in scratch).
        target_bed: Target BED (see read_target_bed).
        out: Output table, or None.
        icr_region: ICR panel; regions overlapping it use the ICR threshold.
        icr_gene_list: Imprinting gene BED, used the same way (default:
            bundled).
        haplotypes: Report the hap1 and hap2 tracks.
        min_coverage: Merged track: keep a CpG with coverage above this.
        hap_min_coverage: Haplotype tracks: the same.
        fdr_threshold: A region is called when P_FDR is below this.
        cpg_diff_threshold: Minimum |DIFF| of a non-ICR region.
        icr_diff_threshold: Minimum |DIFF| of an ICR region.
        input_format: Format of the tracks: auto (detect), 5col, modkit or
            pbcpg (the `input_format` setting).
        scratch: Folder for indexes (default: a temporary folder).
        log: Receives progress messages.
        threads: Threads for reading.

    Returns:
        The table with SAMPLE_ID, SAMPLE_TYPE, the BED columns, HAPLOTYPE,
            NUM_CPG,
        SAMPLE_METH (%), REF_METH (%), REF_SD (%), COVERAGE, DIFF, P_FDR and
        DIRECTION.
    """
    owns = scratch is None
    scratch_ctx = tempfile.TemporaryDirectory(
        prefix="dmrdetector_target_") if owns else None
    scratch = scratch_ctx.name if owns else scratch
    try:
        tgt = read_target_bed(target_bed, log)
        if tgt.empty:
            raise ValueError("no autosomal target regions left")
        if icr_region is None and icr_gene_list is None:
            icr_gene_list = DEFAULT_ICR_GENES
        is_icr = icr_overlap(tgt, icr_region, icr_gene_list)
        log("ICR thresholds used for: " +
            (", ".join(tgt.loc[is_icr, "REGION_NAME"]) or "none"))
        thr = np.where(is_icr, icr_diff_threshold, cpg_diff_threshold)

        atl = _atlas_over_regions(atlas, tgt, scratch, threads)

        samples = sorted(
            samples,
            key=lambda d:
            (0 if str(d["type"]).upper() == "CASE" else 1, d["sample_id"]))
        tables = []
        for smp in samples:
            het = None
            if smp.get("het") is not None:
                het = smp["het"]
            elif smp.get("vcf"):
                het = _HetSites(smp["vcf"], tgt)
            tracks = [("merged", smp["merged"])]
            if haplotypes:
                tracks += [(t, smp.get(t)) for t in ("hap1", "hap2")]
            rows, pvals = [], np.full(len(tgt), np.nan)
            for label, path in tracks:
                if not path:
                    continue
                cov = min_coverage if label == "merged" else hap_min_coverage
                res, diffs = _summarise(
                    _read_track_sites(path, tgt, scratch, threads,
                                      input_format), atl, cov, het, tgt,
                    label == "merged")
                if label == "merged":
                    pvals = np.array([
                        rstats.t_test_one_sample_p(diffs[r])
                        if r in diffs else np.nan for r in tgt["REGION_ID"]
                    ])
                out_t = tgt.drop(columns="REGION_ID").copy()
                out_t["HAPLOTYPE"] = label
                out_t["NUM_CPG"] = res["NUM_CPG"].to_numpy()
                out_t["SAMPLE_METH (%)"] = rstats.r_round(
                    res["SAMPLE_METH"].to_numpy() * 100, 1)
                out_t["REF_METH (%)"] = rstats.r_round(
                    res["REF_METH"].to_numpy() * 100, 1)
                out_t["REF_SD (%)"] = rstats.r_round(
                    res["REF_SD"].to_numpy() * 100, 1)
                out_t["COVERAGE"] = rstats.r_round(res["COVERAGE"].to_numpy(),
                                                   1)
                out_t["DIFF"] = res["DIFF"].to_numpy(
                ) if label == "merged" else np.nan
                out_t["P_FDR"] = np.nan
                out_t["DIRECTION"] = None
                if label == "merged":
                    evaluable = ~np.isnan(pvals)
                    padj = np.full(len(tgt), np.nan)
                    padj[evaluable] = rstats.p_adjust_bh(pvals[evaluable])
                    d = out_t["DIFF"].to_numpy()
                    dirn = np.where(
                        padj < fdr_threshold,
                        np.where(d > thr, "hyper",
                                 np.where(d < -thr, "hypo", "normal")),
                        "normal")
                    out_t["P_FDR"] = padj
                    out_t["DIRECTION"] = np.where(evaluable, dirn, None)
                out_t.insert(0, "SAMPLE_TYPE", str(smp["type"]).upper())
                out_t.insert(0, "SAMPLE_ID", smp["sample_id"])
                rows.append(out_t)
            if not rows:
                continue
            df = pd.concat(rows, ignore_index=True)
            if not (df.loc[df["HAPLOTYPE"] == "merged", "NUM_CPG"] > 0).any():
                log(f"WARNING: no target region of {smp['sample_id']} has a "
                    "CpG in the sample that is also in the reference atlas "
                    "(contig names or coverage); the report lists every "
                    "region with NUM_CPG 0 and no P_FDR or DIRECTION")
            df["_h"] = df["HAPLOTYPE"].map({t: i for i, t in enumerate(TRACKS)})
            df = df.sort_values(["CHROM", "START", "END", "_h"],
                                kind="stable").drop(columns="_h")
            tables.append(df)
            if isinstance(het, _HetSites):
                het.close()
            log(f"Wrote {len(df)} target rows for {smp['sample_id']}")
        res_all = pd.concat(
            tables, ignore_index=True)[OUT_COLS] if tables else pd.DataFrame(
                columns=OUT_COLS)
    finally:
        if scratch_ctx is not None:
            scratch_ctx.cleanup()
    if out:
        write_target_table(res_all, out)
    return res_all


def write_target_table(df: pd.DataFrame, path: str) -> None:
    """Writes the table tab-separated and unquoted, with NA as `NA`.

    This is v2 fwrite(na = "NA"); the file is gzip when `path` ends in .gz.
    """
    cols = []
    for c in OUT_COLS:
        s = df[c]
        if c in _FLOAT_OUT:
            cols.append(
                [bed.fmt_double(v) or "NA" for v in s.to_numpy(dtype=float)])
        elif c in ("START", "END", "NUM_CPG"):
            cols.append([str(int(v)) for v in s])
        else:
            cols.append([
                "NA" if (v is None or
                         (isinstance(v, float) and v != v)) else str(v)
                for v in s
            ])
    text = "\n".join(["\t".join(OUT_COLS)] + ["\t".join(r)
                                              for r in zip(*cols)]) + "\n"
    if path.endswith(".gz"):
        import gzip
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        with open(path, "w") as fh:
            fh.write(text)


def read_samplesheet(path: str, haplotypes: bool = True) -> list[dict]:
    """Reads a CSV of sample_id,type,meth_bed[,vcf,meth_bed_hap1,meth_bed_hap2].

    A sex column is ignored. The haplotype BEDs come from the two columns only;
    file names are not interpreted.
    """
    out = []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            row = {k.strip(): (v or "").strip() for k, v in row.items() if k}
            d = {
                "sample_id": row["sample_id"],
                "type": row["type"].upper(),
                "merged": row["meth_bed"],
                "vcf": row.get("vcf") or None
            }
            if haplotypes:
                for tag, col in (("hap1", "meth_bed_hap1"), ("hap2",
                                                             "meth_bed_hap2")):
                    p = row.get(col)
                    if p:
                        d[tag] = p
            out.append(d)
    return out


def main(args: argparse.Namespace) -> int:
    """Runs the `target` command."""
    samples = read_samplesheet(args.samplesheet, not args.no_haplotypes)
    log = bed.warn
    try:
        target_report(samples,
                      args.atlas,
                      args.target_bed,
                      args.out,
                      icr_region=args.icr_region,
                      icr_gene_list=args.icr_gene_list,
                      haplotypes=not args.no_haplotypes,
                      min_coverage=args.target_min_coverage,
                      hap_min_coverage=args.target_hap_min_coverage,
                      fdr_threshold=args.fdr_threshold,
                      cpg_diff_threshold=args.cpg_diff_threshold,
                      icr_diff_threshold=args.icr_diff_threshold,
                      input_format=args.input_format,
                      log=log,
                      threads=args.threads)
    except runtime.RunError as e:
        log(f"ERROR: {e}")
        return 1
    return 0


def add_parser(sub: argparse._SubParsersAction) -> None:
    """Adds the `target` command to `sub`."""
    p = sub.add_parser("target",
                       help="target-region report (target_region.tsv)",
                       description="Target-region report: per-sample "
                       "methylation "
                       "of disease-critical regions (v2 module, ported).")
    p.set_defaults(func=main)
    p.add_argument(
        "--samplesheet",
        required=True,
        help="CSV: sample_id,type,meth_bed[,vcf,meth_bed_hap1,meth_bed_hap2]")
    p.add_argument("--atlas",
                   required=True,
                   help="autosomal narrow-bandwidth atlas BED.gz")
    p.add_argument(
        "--target-bed",
        required=True,
        help=
        "CHROM START END REGION_NAME GENE_SYMBOL DISEASE_NAME (header allowed)")
    p.add_argument("--out", default="target_region.tsv")
    p.add_argument(
        "--icr-region",
        help=
        "icr_region.bed in use (regions overlapping it use the ICR thresholds)")
    p.add_argument(
        "--icr-gene-list",
        help="icr_gene_list.bed (imprinted genes; same use); default: bundled "
        "assets/icr_gene_list.bed "
        "when neither ICR file is given")
    p.add_argument("--no-haplotypes",
                   action="store_true",
                   help="merged track only")
    p.add_argument("--threads",
                   type=int,
                   default=4,
                   help="threads for BEDs without .tbi (read in memory)")
    p.add_argument("--target-min-coverage",
                   type=float,
                   default=settings.TARGET_MIN_COVERAGE,
                   help="merged track: coverage must exceed this")
    p.add_argument("--target-hap-min-coverage",
                   type=float,
                   default=settings.TARGET_HAP_MIN_COVERAGE,
                   help="hap1/hap2: coverage must exceed this")
    p.add_argument("--input-format",
                   choices=settings.INPUT_FORMATS,
                   default=settings.INPUT_FORMAT,
                   help="format of the methylation BEDs [%(default)s]")
    p.add_argument("--cpg-diff-threshold",
                   type=float,
                   default=call_dmr.DEFAULTS["cpg_diff_threshold"])
    p.add_argument("--icr-diff-threshold",
                   type=float,
                   default=call_dmr.DEFAULTS["icr_diff_threshold"])
    p.add_argument("--fdr-threshold",
                   "--dmr-fdr-threshold",
                   dest="fdr_threshold",
                   type=float,
                   default=call_dmr.DEFAULTS["fdr_threshold"])
