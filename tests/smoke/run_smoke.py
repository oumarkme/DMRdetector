#!/usr/bin/env python3
"""End-to-end smoke test of `dmrdetector build` and `dmrdetector run`.

It uses synthetic data written in pure Python. Scenarios: build (sex-stratified,
autosomal only, too few controls of one sex, `type` column, failure and
resume, --keep-tmp), one `run` per case (output folder, temporary folder, log,
no --sex, v1-layout reference, target region with and without haplotype BEDs),
the three BED formats, settings file, --set and flag precedence, bad inputs.

  PYTHONPATH=src python tests/smoke/run_smoke.py [--workdir DIR] [--threads N]
  [--keep]

Needs numpy, pandas and pysam. Writes only under --workdir (default: a new
directory in $TMPDIR, removed at the end unless --keep).
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import pysam


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DATA = Path(os.environ.get("DMRDETECTOR_TEST_DATA",
                           ROOT / "tests" / "data"))  # test data folder
GTF_FULL = Path(
    os.environ.get("DMRDETECTOR_TEST_GTF",
                   DATA / "Homo_sapiens.GRCh38.86.gtf.gz"))

CTRL = [(f"CTRL{i:02d}", "M" if i <= 4 else "F") for i in range(1, 9)
       ]  # 4 M + 4 F
CASE = [("CASE01", "F"), ("CASE02", "M")]
ISLANDS22 = [16_050_000, 20_000_000, 30_000_000, 40_000_000, 45_000_000]
ICR_ISLAND = 20_309_238  # inside DGCR6L of the bundled icr_gene_list.bed
ISLAND_LEN = 1500
STEP = 5


# relaxed thresholds for the tiny data (algorithms unchanged, as in the v1 smoke
# test)
def _set(**kv) -> list[str]:
    return [x for k, v in kv.items() for x in ("--set", f"{k}={v}")]


BUILD_OPTS = _set(cpg_min_samples_per_site=2,
                  cpg_min_island_length=100,
                  icr_min_sites=5)
RUN_OPTS = _set(dmr_cpg_min_sites=10,
                dmr_cpg_diff_threshold=0.3,
                icr_min_sites=5,
                icr_diff_threshold=0.25,
                dmr_fdr_threshold=0.05)

results: list[tuple[str, bool, str]] = []


def check(name: str, ok: bool, detail: str = "") -> bool:
    results.append((name, bool(ok), detail))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" +
          (f"  ({detail})" if detail and not ok else ""),
          flush=True)
    return bool(ok)


# ---------------------------------------------------------------------------
# synthetic data
# ---------------------------------------------------------------------------


def site_list(rng) -> list[tuple[str, int]]:
    sites = []
    for s in ISLANDS22 + [ICR_ISLAND]:
        sites += [("chr22", p) for p in range(s, s + ISLAND_LEN, STEP)]
    sites += [("chr22", int(p)) for p in rng.choice(
        np.arange(1_000_000, 49_000_000), 3000, replace=False)]
    sites += [
        ("chrX", p) for p in range(30_000_000, 30_000_000 + ISLAND_LEN, STEP)
    ]
    sites += [("chrX", int(p)) for p in rng.choice(
        np.arange(5_000_000, 150_000_000), 400, replace=False)]
    sites += [
        ("chrY", p) for p in range(20_000_000, 20_000_000 + ISLAND_LEN, STEP)
    ]
    sites += [("chrY", int(p)) for p in rng.choice(
        np.arange(5_000_000, 55_000_000), 100, replace=False)]
    return sorted(set(sites))


def p_meth(chrom: str, pos: int, kind: str, sex: str) -> float:
    if chrom == "chr22":
        if any(s <= pos < s + ISLAND_LEN for s in ISLANDS22):
            return 0.02 if kind == "CTRL" else 0.8
        if ICR_ISLAND <= pos < ICR_ISLAND + ISLAND_LEN:
            return 0.5 if kind == "CTRL" else 0.95
    elif chrom == "chrX" and 30_000_000 <= pos < 30_000_000 + ISLAND_LEN:
        if kind == "CTRL":
            return 0.5 if sex == "F" else 0.02
        return 0.95 if sex == "F" else 0.8
    elif chrom == "chrY" and 20_000_000 <= pos < 20_000_000 + ISLAND_LEN:
        return 0.02 if kind == "CTRL" else 0.8
    return 0.9


def to_modkit(c, s, e, n,
              k) -> str:  # 18 columns: Nvalid_cov = col 10, Nmod = col 12
    return "\t".join(
        map(str, [
            c, s, e, "m", n, "+", s, e, "255,0,0", n,
            round(100 * k / n, 2), k, n - k, 0, 0, 0, 0, 0
        ])) + "\n"


def to_pbcpg(
        c, s, e, n,
        k) -> str:  # 9 columns with header: cov = col 6, est_mod_count = col 7
    return "\t".join(
        map(str, [
            c, s, e,
            round(100 * k / n, 2), "Total", n, k, n - k,
            round(k / n, 3)
        ])) + "\n"


FORMATS = {
    "5col": (lambda c, s, e, n, k: f"{c}\t{s}\t{e}\t{n}\t{k}\n", ""),
    "modkit": (to_modkit, ""),
    "pbcpg": (
        to_pbcpg,
        "chrom\tbegin\tend\tmod_score\ttype\tcov\test_mod_count\test_unmod_count\tdiscretized_mod_score\n"
    ),
}


def write_data(outdir: Path) -> None:
    rng = np.random.default_rng(20261006)
    sites = site_list(rng)
    het = [s for s in sites if s[0] == "chr22" and s[1] % 7 == 0][:20]
    for fmt, (fn, header) in FORMATS.items():
        d = outdir / "data" / fmt
        d.mkdir(parents=True, exist_ok=True)
        rng2 = np.random.default_rng(7)  # same counts in every format
        for sid, kind, sex in [(s, "CTRL", x) for s, x in CTRL
                              ] + [(s, "CASE", x) for s, x in CASE]:
            cov = np.maximum(12, rng2.poisson(25, len(sites)))
            plain = d / f"{sid}.bed"
            with open(plain, "w") as fh:
                fh.write(header)
                for (c, p), n in zip(sites, cov):
                    fh.write(
                        fn(c, p, p + 1, int(n),
                           int(rng2.binomial(int(n), p_meth(c, p, kind, sex)))))
            pysam.tabix_compress(str(plain),
                                 str(d / f"{sid}.bed.gz"),
                                 force=True)  # bgzip, as real inputs
            plain.unlink()
            if fmt == "5col":
                with open(d / f"{sid}.vcf", "w") as fh:
                    fh.write("##fileformat=VCFv4.2\n##contig=<ID=chr22>\n")
                    fh.write(
                        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                    )
                    fh.write(
                        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sid}\n"
                    )
                    for c, p in het:
                        fh.write(
                            f"{c}\t{p + 1}\t.\tC\tT\t50\tPASS\t.\tGT\t0/1\n")


def sheet(path: Path,
          data: Path,
          fmt="5col",
          subset=None,
          extra_col=False) -> Path:
    """Writes the control samplesheet of `build` and returns its path.

    Columns: sample_id,sex,meth_bed,vcf. With `extra_col` the v2.2 layout
    sample_id,type,sex,meth_bed,vcf: the controls are CTRL rows, followed by
    the two CASE rows of the data set (`subset` does not filter these).
    """
    cols = ["sample_id"
           ] + (["type"] if extra_col else []) + ["sex", "meth_bed", "vcf"]
    rows = [",".join(cols)]
    for s, x in CTRL:
        if subset is None or s in subset:
            cells = [s] + (["CTRL"] if extra_col else []) + [
                x, f"{data}/{fmt}/{s}.bed.gz", f"{data}/5col/{s}.vcf"
            ]
            rows.append(",".join(cells))
    if extra_col:
        for s, x in CASE:
            rows.append(",".join([
                s, "CASE", x, f"{data}/{fmt}/{s}.bed.gz", f"{data}/5col/{s}.vcf"
            ]))
    path.write_text("\n".join(rows) + "\n")
    return path


def sample_args(data: Path,
                sid: str,
                sex: str | None = None,
                fmt="5col",
                hap=False) -> list[str]:
    """Per-sample options of `run`."""
    a = [
        "--meth_bed_combined", f"{data}/{fmt}/{sid}.bed.gz", "--vcf",
        f"{data}/5col/{sid}.vcf", "--sample-id", sid
    ]
    if sex:
        a += ["--sex", sex]
    if hap:
        a += [
            "--meth_bed_hap1", f"{data}/{fmt}/{sid}.bed.gz", "--meth_bed_hap2",
            f"{data}/{fmt}/{sid}.bed.gz"
        ]
    return a


def small_gtf(dest: Path) -> Path:
    """Writes the test GTF restricted to chr22, X and Y and returns its path.

    The synthetic data lie there, so the file is fast to load.
    """
    with gzip.open(GTF_FULL, "rt") as fi, gzip.open(dest, "wt") as fo:
        for line in fi:
            if line[0] == "#" or line.split("\t", 1)[0] in ("22", "X", "Y"):
                fo.write(line)
    return dest


# ---------------------------------------------------------------------------
# running
# ---------------------------------------------------------------------------

# empty folder for commands without their own cwd: ./settings.toml of the
# repository root must not be read
NEUTRAL_CWD = None


def cli(args: list[str], cwd=None, env=None) -> subprocess.CompletedProcess:
    cwd = cwd or NEUTRAL_CWD
    e = dict(os.environ,
             PYTHONPATH=os.environ.get("DMRDETECTOR_SRC", str(ROOT / "src")),
             **(env or {}))
    return subprocess.run([sys.executable, "-m", "dmrdetector.cli", *args],
                          capture_output=True,
                          text=True,
                          env=e,
                          cwd=cwd)


def snapshot(d: Path) -> dict:
    return {
        str(p.relative_to(d)): p.stat().st_mtime_ns
        for p in d.rglob("*")
        if p.is_file()
    }


def dmr_rows(path: Path) -> list[list[str]]:
    with open(path) as fh:
        return [ln.rstrip("\n").split("\t") for ln in fh][1:]


def gz_bytes(path: Path) -> bytes:
    with gzip.open(path, "rb") as fh:
        return fh.read()


def atlas_index_problems(d: Path) -> list[str]:
    """Returns the atlases of a reference directory that lack a usable `.tbi`.

    Every smoothed atlas must have its `.tbi`, and tabix must return the rows of
    every contig.
    """
    from collections import Counter
    atlases = sorted(d.glob("ref_atlas_*.bed.gz"))
    bad = [] if atlases else ["no atlas"]
    for a in atlases:
        if not Path(str(a) + ".tbi").is_file():
            bad.append(f"{a.name}: no .tbi")
            continue
        want = Counter(
            ln.split(b"\t", 1)[0].decode() for ln in gz_bytes(a).splitlines())
        with pysam.TabixFile(str(a)) as tbx:
            got = {c: sum(1 for _ in tbx.fetch(c)) for c in tbx.contigs}
        if got != dict(want):
            bad.append(f"{a.name}: tabix rows {got} != {dict(want)}")
    return bad


def same_files(a: Path, b: Path, names) -> bool:
    return all((a / n).is_file() and (b / n).is_file() and
               (a / n).read_bytes() == (b / n).read_bytes() for n in names)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--workdir")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--keep", action="store_true")
    a = ap.parse_args()
    if not GTF_FULL.is_file():
        print(
            f"Ensembl GTF not found: {GTF_FULL}\nDownload one with "
            "`dmrdetector get-genes --release 86 --outdir DIR` and set "
            "DMRDETECTOR_TEST_GTF to the file.",
            file=sys.stderr)
        return 2
    base = Path(a.workdir) if a.workdir else Path(
        tempfile.mkdtemp(prefix="dmrdetector_smoke_"))
    base.mkdir(parents=True, exist_ok=True)
    global NEUTRAL_CWD
    NEUTRAL_CWD = base / "neutral"
    NEUTRAL_CWD.mkdir(exist_ok=True)
    try:
        return scenarios(base, a.threads)
    finally:
        if not a.keep:
            shutil.rmtree(base, ignore_errors=True)
        else:
            print(f"kept {base}")


def listing(d: Path) -> list[str]:
    return sorted(p.name for p in d.iterdir()) if d.is_dir() else []


def scenarios(base: Path, threads: int) -> int:
    write_data(base)
    data = base / "data"
    gtf = small_gtf(base / "genes.gtf.gz")
    full = sheet(base / "full.csv", data)
    res = ["--threads", str(threads)]
    scratch = base / "tmp"  # --tmp-dir of every run: must be empty afterwards
    scratch.mkdir()

    def build_cmd(sh, out, *x):
        return cli([
            "build", "--input",
            str(sh), "--outdir",
            str(out), *res, *BUILD_OPTS, *x
        ])

    def run_cmd(sid,
                sex,
                out,
                *x,
                ref=None,
                fmt="5col",
                hap=False,
                opts=RUN_OPTS,
                g=True,
                tmp=True):
        a = [
            "run", *sample_args(data, sid, sex, fmt, hap), "--outdir",
            str(out), *res, *opts
        ]
        if ref is not None:
            a += ["--reference", str(ref)]
        if g:
            a += ["--gtf", str(gtf)]
        if tmp:
            a += ["--tmp-dir", str(scratch)]
        return cli([*a, *x])

    # ---- 1. build: sex-stratified reference, controls only ----
    ref = base / "ref_sex"
    r = build_cmd(full, ref)
    check("build (sex-stratified): exit 0", r.returncode == 0, r.stderr[-1500:])
    atl = [
        f"ref_atlas_150{s}.bed.gz" for s in ("", "_male_X", "_female_X", "_Y")
    ]
    want = sorted(
        atl + [f"{a}.tbi" for a in atl] +
        [f"cpg_signal{s}.bed" for s in ("", "_male_X", "_female_X", "_Y")] +
        ["icr_region.bed", "reference.json", "build.log"])
    check(
        "build: the reference directory holds exactly the reference files and "
        "build.log (no work folder, no run_info.json)",
        listing(ref) == want, str(listing(ref)))
    probs = atlas_index_problems(ref)
    check("build: every smoothed atlas has its .tbi and tabix returns all rows",
          not probs and len(atl) == 4, "; ".join(probs))
    if (ref / "reference.json").is_file():
        rj = json.loads((ref / "reference.json").read_text())
        import hashlib
        ok = (rj["samplesheet_md5"] == hashlib.md5(
            full.read_bytes()).hexdigest() and len(rj["samples"]) == 8 and
              set(rj["scopes"]) == {"autosomal", "male_X", "female_X", "Y"} and
              rj["scopes"]["autosomal"]["n_sites"] > 1000 and
              rj["parameters"]["smooth_sigma"] == 150 and
              rj["scopes"]["male_X"]["n_ctrl"] == 4 and
              rj["samplesheet"] == full.name and str(base) not in str(rj) and
              rj["icr_gene_list_md5"] == hashlib.md5(
                  (Path(os.environ.get("DMRDETECTOR_SRC", ROOT / "src")) /
                   "dmrdetector/data/assets/icr_gene_list.bed").read_bytes()
              ).hexdigest())
        check(
            "build: reference.json (md5, samples, scopes, atlas stats, "
            "parameters, md5 of the ICR gene list; the samplesheet by file "
            "name, no absolute path)", ok, str(rj)[:300])
    log = (ref / "build.log").read_text() if (ref /
                                              "build.log").is_file() else ""
    check(
        "build.log: version, command line, parameters with source, step "
        "timings, status", "DMRdetector 3." in log and "command:" in log and
        "cpg_min_island_length = 100  # flag" in log and
        "smooth_sigma = 150  # default" in log and "[steps]" in log and
        "status:   completed" in log and str(base) not in log, log[:400])
    vl = next((x for x in log.splitlines() if x.startswith("versions:")), "")
    check(
        "build.log: versions of python, numpy, scipy, pandas, pysam and numba",
        all(f" {n} " in vl + " "
            for n in ("python", "numpy", "scipy", "pandas", "pysam", "numba"))
        and "not installed" not in vl, vl)
    r = build_cmd(sheet(base / "legacy.csv", data, extra_col=True),
                  base / "ref_legacy")
    n_max, n_ids = 0, []
    if (base / "ref_legacy/ref_atlas_150.bed.gz").is_file():
        with gzip.open(base / "ref_legacy/ref_atlas_150.bed.gz", "rt") as fh:
            n_max = max(int(x.split("\t")[3]) for x in fh)
        n_ids = [
            s["sample_id"]
            for s in json.loads((base / "ref_legacy/reference.json").read_text(
            ))["samples"]
        ]
    check(
        "build: a `type` column keeps the CTRL rows (8 CTRL and 2 CASE rows "
        "give an atlas with N = 8), and a NOTICE gives the number of rows "
        "dropped", r.returncode == 0 and n_max == len(CTRL) and
        n_ids == [s for s, _ in CTRL] and
        "2 of 10 rows are dropped" in r.stderr, f"N max {n_max}; " +
        r.stderr[-400:])
    only_case = base / "only_case.csv"
    only_case.write_text("sample_id,type,sex,meth_bed,vcf\n" + "".join(
        f"{s},CASE,{x},{data}/5col/{s}.bed.gz,{data}/5col/{s}.vcf\n"
        for s, x in CASE))
    r = build_cmd(only_case, base / "ref_only_case")
    check(
        "build: a `type` column without a CTRL row stops with an error",
        r.returncode == 1 and "no row with type CTRL" in r.stderr and
        not (base / "ref_only_case/ref_atlas_150.bed.gz").exists(),
        r.stderr[-300:])

    typo = base / "typo.csv"
    typo.write_text(
        (base / "legacy.csv").read_text().replace(",CTRL,", ",contol,", 1))
    r = build_cmd(typo, base / "ref_typo")
    check(
        "build: a `type` value that is neither a control nor CASE stops with an "
        "error that names the value, as in v2.2", r.returncode == 1 and
        "neither a control nor CASE" in r.stderr and
        "'contol' (row 2)" in r.stderr and
        not (base / "ref_typo/ref_atlas_150.bed.gz").exists(), r.stderr[-300:])

    # ---- 2. build: --keep-tmp, resume after a failure ----
    wk = base / "ref_keep_work"
    r = build_cmd(full, base / "ref_keep", "--tmp-dir", str(wk), "--keep-tmp")
    kept = sorted(wk.glob("dmrdetector-build-*"))
    check(
        "build --keep-tmp: work files kept in "
        "<tmp-dir>/dmrdetector-build-<hash> and the path printed",
        r.returncode == 0 and len(kept) == 1 and
        (kept[0] / "raw_atlas").is_dir() and
        (kept[0] / ".dmrdetector-build").is_file() and str(kept[0]) in r.stderr,
        r.stderr[-400:])
    shutil.rmtree(wk, ignore_errors=True)
    junk = base / "junk.bed.gz"
    junk.write_text("not a bed file\n")
    badsheet = base / "bad_ctrl.csv"
    badsheet.write_text("\n".join(
        ln.replace(f"{data}/5col/CTRL03.bed.gz", str(junk))
        for ln in full.read_text().splitlines()) + "\n")
    outf = base / "ref_fail"
    r = build_cmd(badsheet, outf)
    check(
        "build with an unreadable control: exit 1, work folder kept, "
        "build.log says FAILED", r.returncode == 1 and
        (outf / "work/stamps").is_dir() and
        "status:   FAILED" in (outf / "build.log").read_text() and
        not (outf / "ref_atlas_150.bed.gz").exists(), r.stderr[-500:])
    r = build_cmd(full, outf)
    lg = (outf / "build.log").read_text() if (outf /
                                              "build.log").is_file() else ""
    check(
        "build rerun resumes: finished preprocessing is skipped, work folder "
        "removed after success", r.returncode == 0 and
        "preprocess_CTRL01: skipped" in lg and
        "preprocess_CTRL03: done" in lg and not (outf / "work").exists() and
        (outf / "ref_atlas_150.bed.gz").is_file(), lg[-600:] + r.stderr[-300:])
    check(
        "build: the rerun gives the same reference as the first build",
        gz_bytes(outf / "ref_atlas_150.bed.gz") == gz_bytes(
            ref / "ref_atlas_150.bed.gz") and
        (outf / "cpg_signal.bed").read_bytes()
        == (ref / "cpg_signal.bed").read_bytes())

    # ---- 3. run, one sample per invocation ----
    o2 = base / "run_F"
    r = run_cmd("CASE01", "F", o2, ref=ref)
    check("run (CASE01, F): exit 0", r.returncode == 0, r.stderr[-1500:])
    check("run: the output folder holds exactly <id>.dmr.tsv and <id>.log",
          listing(o2) == ["CASE01.dmr.tsv", "CASE01.log"], str(listing(o2)))
    check("run: the temporary folder is removed (--tmp-dir is empty)",
          listing(scratch) == [], str(listing(scratch)))
    o3 = base / "run_M"
    r = run_cmd("CASE02", "M", o3, ref=ref)
    chroms = {
        s: ({x[0] for x in dmr_rows(o / f"{s}.dmr.tsv")} if
            (o / f"{s}.dmr.tsv").is_file() else set()
           ) for s, o in (("CASE01", o2), ("CASE02", o3))
    }
    check("run: CASE01 (F) chr22 + chrX, no chrY",
          {"chr22", "chrX"} <= chroms["CASE01"] and
          "chrY" not in chroms["CASE01"], str(chroms))
    check("run: CASE02 (M) chr22 + chrX + chrY", {"chr22", "chrX", "chrY"}
          <= chroms["CASE02"], str(chroms))
    if (o2 / "CASE01.dmr.tsv").is_file():
        rows = dmr_rows(o2 / "CASE01.dmr.tsv")
        check("run: ICR panel DMR found", any(x[4] == "ICR" for x in rows),
              str(rows[:2]))
    lg = (o2 / "CASE01.log").read_text() if (o2 /
                                             "CASE01.log").is_file() else ""
    check(
        "log: version, command line, inputs, parameters with source, step "
        "timings, status", "DMRdetector 3." in lg and "command:" in lg and
        "sample_id = CASE01" in lg and
        "dmr_fdr_threshold = 0.05  # flag" in lg and
        "smooth_sigma = 150  # default" in lg and
        "call_dmr_CASE01_autosomal: done" in lg and
        "status:   completed" in lg and str(base) not in lg and
        f"threads: {threads}, memory budget " in lg, lg[:600])
    r2 = run_cmd("CASE01", "F", base / "run_F2", ref=ref)
    check(
        "run twice: identical DMR table", r2.returncode == 0 and
        same_files(o2, base / "run_F2", ["CASE01.dmr.tsv"]))

    # ---- 4. --keep-tmp, no --sex, sample ID ----
    o4 = base / "run_keep"
    kt = base / "keep_tmp"
    r_kt = run_cmd("CASE01",
                   "F",
                   o4,
                   "--keep-tmp",
                   "--tmp-dir",
                   str(kt),
                   ref=ref,
                   tmp=False)
    kept = [p for p in kt.glob("dmrdetector-*")] if kt.is_dir() else []
    check(
        "run --keep-tmp: the temporary folder is kept in --tmp-dir and its "
        "path is printed", r_kt.returncode == 0 and len(kept) == 1 and
        (kept[0] / "preprocessed").is_dir() and str(kept[0]) in r_kt.stderr and
        listing(o4) == ["CASE01.dmr.tsv", "CASE01.log"], r_kt.stderr[-400:])
    o5 = base / "run_nosex"
    r = run_cmd("CASE02", None, o5, ref=ref)
    ch = {x[0] for x in dmr_rows(o5 / "CASE02.dmr.tsv")
         } if (o5 / "CASE02.dmr.tsv").is_file() else set()
    check(
        "run without --sex: one-line notice, autosomal scope only (chr22 "
        "DMRs, no chrX/chrY)", r.returncode == 0 and
        r.stderr.count("--sex not given") == 1 and "chr22" in ch and
        not (ch & {"chrX", "chrY"}), r.stderr[-400:] + str(ch))
    o6 = base / "run_defid"
    a = [
        "run", "--meth_bed_combined", f"{data}/5col/CASE01.bed.gz", "--vcf",
        f"{data}/5col/CASE01.vcf", "--sex", "F", "--reference",
        str(ref), "--gtf",
        str(gtf), "--outdir",
        str(o6), *res, *RUN_OPTS
    ]
    r = cli(a)
    check("run: the default sample ID is the BED file name up to its first '.'",
          r.returncode == 0 and listing(o6) == ["CASE01.dmr.tsv", "CASE01.log"],
          r.stderr[-300:])
    r = cli([
        "run", "--meth-bed-combined", f"{data}/5col/CASE01.bed.gz", "--vcf",
        f"{data}/5col/CASE01.vcf", "--sex", "F", "--sample-id", "X.1",
        "--reference",
        str(ref), "--gtf",
        str(gtf), "--outdir",
        str(base / "run_hy"), *res, *RUN_OPTS
    ])
    check(
        "run: hyphen spelling --meth-bed-combined and --sample-id",
        r.returncode == 0 and
        listing(base / "run_hy") == ["X.1.dmr.tsv", "X.1.log"], r.stderr[-300:])
    wd = base / "cwd"
    wd.mkdir()
    r = cli([
        "run", *sample_args(data, "CASE01", "F"), "--reference",
        str(ref), "--gtf",
        str(gtf), *res, *RUN_OPTS
    ],
            cwd=wd)
    check(
        "run: the default --outdir is ./dmrdetector", r.returncode == 0 and
        listing(wd / "dmrdetector") == ["CASE01.dmr.tsv", "CASE01.log"],
        r.stderr[-300:])

    # ---- 5. reference as v1 layout; autosomal-only reference ----
    v1 = base / "ref_v1layout"
    for sub in ("atlas", "cpg_islands", "icr_regions"):
        (v1 / sub).mkdir(parents=True, exist_ok=True)
    for sfx in ("", "_male_X", "_female_X", "_Y"):
        shutil.copy(ref / f"ref_atlas_150{sfx}.bed.gz",
                    v1 / "atlas" / f"ref_atlas_150{sfx}.tsv.gz")
        shutil.copy(ref / f"cpg_signal{sfx}.bed",
                    v1 / "cpg_islands" / f"cpg_signal{sfx}.bed")
    shutil.copy(ref / "icr_region.bed", v1 / "icr_regions" / "icr_region.bed")
    o7 = base / "run_v1ref"
    r = run_cmd("CASE01", "F", o7, ref=v1)
    check(
        "run vs v1-layout reference (atlas/, cpg_islands/, icr_regions/; "
        ".tsv.gz atlases): same DMR table", r.returncode == 0 and
        same_files(o2, o7, ["CASE01.dmr.tsv"]), r.stderr[-800:])
    refa = base / "ref_auto"
    r = build_cmd(full, refa, "--set", "include_sex_chromosomes=false")
    check(
        "build with include_sex_chromosomes=false: autosomal files only",
        r.returncode == 0 and (refa / "ref_atlas_150.bed.gz").is_file() and
        not list(refa.glob("*male_X*")) and not list(refa.glob("*_Y.*")),
        r.stderr[-800:])
    check("build (autosomal only): the atlas has its .tbi",
          not atlas_index_problems(refa))
    r = run_cmd("CASE01",
                "F",
                base / "run_auto",
                "--set",
                "include_sex_chromosomes=false",
                ref=refa)
    ch = {x[0] for x in dmr_rows(base / "run_auto/CASE01.dmr.tsv")
         } if r.returncode == 0 else set()
    check(
        "run on it: chr22 DMRs only, no warning", r.returncode == 0 and
        "chr22" in ch and not (ch & {"chrX", "chrY"}) and
        "Falling back" not in r.stderr,
        str(ch) + r.stderr[-400:])
    r = run_cmd("CASE01", "F", base / "run_fallback", ref=refa)
    check(
        "run with --sex on an autosomal-only reference: warning and autosomal "
        "scope only", r.returncode == 0 and
        "Falling back to autosomal-only" in r.stderr, r.stderr[-400:])
    few = base / "few.csv"
    few.write_text("\n".join(
        ln for ln in full.read_text().splitlines()
        if not ln.startswith(("CTRL02,", "CTRL03,", "CTRL04,"))) + "\n")
    r = build_cmd(few, base / "ref_few")
    check(
        "build: one male control -> male_X and Y scopes dropped with warnings",
        r.returncode == 0 and "scope male_X" in r.stderr and
        not (base / "ref_few/ref_atlas_150_male_X.bed.gz").exists() and
        (base / "ref_few/ref_atlas_150_female_X.bed.gz").is_file(),
        r.stderr[-600:])

    # ---- 6. settings file, --set and flags ----
    cfgf = base / "my.toml"
    cfgf.write_text(f'[dmr]\ndmr_fdr_threshold = 0.02\ndmr_cpg_min_sites = '
                    f'10\n[reference]\nreference = "{ref}"\ngtf = "{gtf}"\n'
                    f'[resources]\nthreads = 2\n[paths]\noutdir = '
                    f'"{base / "run_cfg"}"\ntmp_dir = "{scratch}"\n')
    r = cli([
        "run", *sample_args(data, "CASE01", "F"), "--config",
        str(cfgf), "--set", "dmr_cpg_diff_threshold=0.3", "--set",
        "dmr_fdr_threshold=0.03"
    ],
            cwd=wd)
    lg = (base / "run_cfg/CASE01.log").read_text() if (
        base / "run_cfg/CASE01.log").is_file() else ""
    check(
        "config: reference, gtf, outdir, tmp_dir and threads come from the "
        "settings file; --set > file > default in the log",
        r.returncode == 0 and "dmr_fdr_threshold = 0.03  # flag" in lg and
        "dmr_cpg_min_sites = 10  # file" in lg and
        "threads = 2  # file" in lg and
        "dmr_merge_distance = 150  # default" in lg and
        f"reference = {ref.name}  # file" in lg and str(base) not in lg,
        r.stderr[-400:] + lg[:900])
    r = cli([
        "run", *sample_args(data, "CASE01", "F"), "--reference",
        str(ref), "--gtf",
        str(gtf), "--outdir",
        str(base / "run_t1"), "--threads", "1", *RUN_OPTS
    ])
    lg = (base / "run_t1/CASE01.log").read_text() if (
        base / "run_t1/CASE01.log").is_file() else ""
    check(
        "threads: --threads 1 is raised to 2 with a warning; threads and "
        "memory budget in the log", r.returncode == 0 and
        "--threads 1 raised to 2" in r.stderr and
        "threads: 2, memory budget " in lg and
        "raised to 2" in lg.split("[warnings]")[1].split("[messages]")[0],
        r.stderr[-400:])
    (wd / "settings.toml").write_text("[dmr]\ndmr_fdr_threshold = 0.04\n")
    r = cli(["config", "show", "run"], cwd=wd)
    check(
        "config show: ./settings.toml is picked up; value, source and meaning "
        "per key",
        "dmr_fdr_threshold = 0.04  # file: report DMRs with FDR" in r.stdout and
        "[dmr]" in r.stdout and "reference = (unset)  # default" in r.stdout and
        "atlas_min_samples_per_site" not in r.stdout and
        "deprecated" not in r.stderr, r.stdout[:300] + r.stderr[-300:])
    (wd / "dmrdetector.toml"
    ).write_text("[dmr]\ndmr_fdr_threshold = 0.07\ndmr_merge_distance = 150\n")
    r = cli(["config", "show", "run"], cwd=wd)
    check(
        "config show: with both files, settings.toml wins and the warning "
        "says so", "dmr_fdr_threshold = 0.04  # file" in r.stdout and
        "dmrdetector.toml is deprecated; rename it to settings.toml" in r.stderr
        and "settings.toml is used" in r.stderr,
        r.stdout[:200] + r.stderr[-300:])
    (wd / "settings.toml").unlink()
    r = cli(["config", "show", "run"], cwd=wd)
    check(
        "config show: ./dmrdetector.toml is still read, with a deprecation "
        "warning on stderr", "dmr_fdr_threshold = 0.07  # file" in r.stdout and
        "dmrdetector.toml is deprecated; rename it to settings.toml"
        in r.stderr, r.stdout[:200] + r.stderr[-300:])
    r = cli([
        "run", *sample_args(data, "CASE01", "F"), "--reference",
        str(ref), "--gtf",
        str(gtf), "--outdir",
        str(base / "run_legacy"), "--threads", "2", *RUN_OPTS
    ],
            cwd=wd)
    lg = (base / "run_legacy/CASE01.log").read_text() if (
        base / "run_legacy/CASE01.log").is_file() else ""
    check(
        "run: the deprecation warning is on stderr and in the [warnings] "
        "section of the log", r.returncode == 0 and
        "dmrdetector.toml is deprecated; rename it to settings.toml" in r.stderr
        and "dmrdetector.toml is deprecated"
        in lg.split("[warnings]")[1].split("[messages]")[0] and
        "dmr_merge_distance = 150  # file" in lg, r.stderr[-400:])
    (wd / "dmrdetector.toml").unlink()
    r = cli(["config", "show", "build", "--set", "dmr_fdr_threshold=0.1"],
            cwd=wd)
    check("config show build: --set of a run-only key is rejected",
          r.returncode == 1 and "applies to run only" in r.stderr,
          r.stderr[-300:])
    r = cli(["config", "init", "--out", str(base / "new.toml")])
    tmpl = ROOT / "settings.toml"
    check(
        "config init writes settings.toml, the same text as the template in "
        "the repository root", r.returncode == 0 and
        (base / "new.toml").read_text() == tmpl.read_text(), r.stderr[-300:])
    for name, text, msg in [
        ("unknown key", "[dmr]\ndmr_fdr_treshold = 0.1\n",
         "unknown key 'dmr_fdr_treshold' in [dmr] (did you mean "
         "'dmr_fdr_threshold'?)"),
        ("bad type", "[dmr]\ndmr_fdr_threshold = \"high\"\n", "must be float"),
        ("out of range", "[dmr]\ndmr_fdr_threshold = 3\n",
         "must be between 0 and 1"),
        ("wrong section", "[atlas]\ndmr_fdr_threshold = 0.1\n",
         "belongs to [dmr]"),
        ("unknown section", "[dmrs]\nx = 1\n", "unknown section [dmrs]"),
        ("bad format name", "[input]\ninput_format = \"bam\"\n",
         "must be one of auto")
    ]:
        f = base / f"bad_{name.replace(' ', '_')}.toml"
        f.write_text(text)
        r = run_cmd("CASE01", "F", base / "x", "--config", str(f), ref=ref)
        check(
            f"config error ({name}): exit 1 and message", r.returncode == 1 and
            msg in r.stderr and "Traceback" not in r.stderr, r.stderr[-300:])
    for name, item, msg in [
        ("unknown key", "dmr_fdr_treshold=0.1",
         "did you mean 'dmr_fdr_threshold'"),
        ("wrong command", "atlas_shrinkage_weight=0.2",
         "applies to build only"),
        ("bad value", "dmr_fdr_threshold=high", "must be float"),
        ("range", "dmr_cpg_diff_threshold=1.5", "must be between 0 and 1"),
        ("no value", "dmr_fdr_threshold", "expected KEY=VALUE")
    ]:
        r = run_cmd("CASE01", "F", base / "x", "--set", item, ref=ref, opts=[])
        check(
            f"--set error ({name}): exit 1 and message", r.returncode == 1 and
            msg in r.stderr and "Traceback" not in r.stderr, r.stderr[-300:])
    for flag_ in (["--dmr-fdr-threshold",
                   "0.05"], ["--input",
                             "x.csv"], ["--include-ctrl"], ["--ref-atlas", "a"],
                  ["--cpg-panel", "a"], ["--icr-panel", "a"], ["--force"],
                  ["--input-format",
                   "5col"], ["--no-icr-v1-duplicates"], ["--cpus", "4"]):
        r = run_cmd("CASE01", "F", base / "x", *flag_, ref=ref)
        check(f"removed flag {flag_[0]}: rejected", r.returncode == 2 and
              "unrecognized arguments" in r.stderr, r.stderr[-200:])
    r = cli([
        "build", "--input",
        str(full), "--outdir",
        str(base / "x"), "--dry-run"
    ])
    check("removed flag build --dry-run: rejected", r.returncode == 2 and
          "unrecognized arguments" in r.stderr, r.stderr[-200:])
    r = cli(["--help"])
    check(
        "--help lists build, run, config, get-genes and not the stage commands",
        all(w in r.stdout for w in ("build", "run", "config", "get-genes")) and
        not any(w in r.stdout
                for w in ("preprocess", "cpg-islands", "call-dmr")) and
        "SUPPRESS" not in r.stdout, r.stdout[:600])
    r = cli(["call-dmr", "--help"])
    check("stage commands stay available (call-dmr --help)",
          r.returncode == 0 and "--atlas" in r.stdout, r.stderr[-200:])

    # ---- 7. input formats ----
    for fmt in ("modkit", "pbcpg"):
        rf = base / f"ref_{fmt}"
        r = build_cmd(sheet(base / f"full_{fmt}.csv", data, fmt), rf)
        same = ((rf / "ref_atlas_150.bed.gz").is_file() and
                gz_bytes(rf / "ref_atlas_150.bed.gz") == gz_bytes(
                    ref / "ref_atlas_150.bed.gz") and
                gz_bytes(rf / "ref_atlas_150_Y.bed.gz") == gz_bytes(
                    ref / "ref_atlas_150_Y.bed.gz") and
                (rf / "cpg_signal.bed").read_bytes()
                == (ref / "cpg_signal.bed").read_bytes())
        check(
            f"build from {fmt} files (auto-detected): atlases and panels "
            "identical to 5-column", r.returncode == 0 and same,
            r.stderr[-800:])
        o8 = base / f"run_{fmt}"
        r = run_cmd("CASE01", "F", o8, ref=ref, fmt=fmt)
        check(
            f"run on a {fmt} file (auto-detected): DMR table identical to "
            "5-column", r.returncode == 0 and
            same_files(o2, o8, ["CASE01.dmr.tsv"]), r.stderr[-800:])
        o8b = base / f"run_{fmt}_forced"
        r = run_cmd("CASE01",
                    "F",
                    o8b,
                    "--set",
                    f"input_format={fmt}",
                    ref=ref,
                    fmt=fmt)
        check(f"run --set input_format={fmt}: same DMR table",
              r.returncode == 0 and same_files(o8, o8b, ["CASE01.dmr.tsv"]),
              r.stderr[-800:])
    r = run_cmd("CASE01",
                "F",
                base / "run_wrongfmt",
                "--set",
                "input_format=pbcpg",
                ref=ref,
                fmt="modkit")
    check(
        "run with input_format=pbcpg on a modkit file: clear error, no "
        "output, temporary folder removed", r.returncode == 1 and
        "cannot read as input_format pbcpg" in r.stderr and
        "Traceback" not in r.stderr and
        listing(base / "run_wrongfmt") == ["CASE01.log"] and
        listing(scratch) == [],
        r.stderr[-400:] + str(listing(base / "run_wrongfmt")))
    check(
        "run failure: the log says FAILED and names the error",
        "status:   FAILED" in (base / "run_wrongfmt/CASE01.log").read_text() and
        "cannot read as input_format"
        in (base / "run_wrongfmt/CASE01.log").read_text())

    # ---- 8. target region ----
    tb = base / "target.bed"
    tb.write_text(
        f"chr22\t{ICR_ISLAND}\t{ICR_ISLAND + ISLAND_LEN}\tT1\tDGCR6L\tTest "
        "syndrome (hypo)\n")
    o9 = base / "run_target"
    r = run_cmd("CASE01",
                "F",
                o9,
                "--target-region",
                str(tb),
                ref=ref,
                hap=True)
    tr = o9 / "CASE01.target_region.tsv"
    txt = tr.read_text() if tr.is_file() else ""
    check(
        "run --target-region FILE: <id>.target_region.tsv with merged, hap1 "
        "and hap2 rows; the folder holds exactly 3 files", r.returncode == 0 and
        listing(o9)
        == ["CASE01.dmr.tsv", "CASE01.log", "CASE01.target_region.tsv"] and
        all(t in txt for t in ("\tmerged\t", "\thap1\t", "\thap2\t")) and
        listing(scratch) == [], r.stderr[-1200:] + txt[:300])
    check("run --target-region: the DMR table equals the run without it",
          same_files(o2, o9, ["CASE01.dmr.tsv"]))
    o9b = base / "run_target_nohap"
    r = run_cmd("CASE01", "F", o9b, "--target-region", str(tb), ref=ref)
    txt = (o9b / "CASE01.target_region.tsv").read_text() if (
        o9b / "CASE01.target_region.tsv").is_file() else ""
    check("run --target-region without hap BEDs: merged track only",
          r.returncode == 0 and "\tmerged\t" in txt and "\thap1\t" not in txt,
          r.stderr[-400:])
    for fmt in ("modkit", "pbcpg"):
        of = base / f"run_target_{fmt}"
        r = run_cmd("CASE01",
                    "F",
                    of,
                    "--target-region",
                    str(tb),
                    ref=ref,
                    fmt=fmt,
                    hap=True)
        check(
            f"run --target-region on {fmt} files (combined, hap1, hap2): "
            "target_region.tsv and DMR table identical to 5-column",
            r.returncode == 0 and listing(scratch) == [] and
            same_files(o9, of, ["CASE01.target_region.tsv", "CASE01.dmr.tsv"]),
            r.stderr[-800:])
    om = base / "run_target_mixed"
    r = run_cmd("CASE01",
                "F",
                om,
                "--target-region",
                str(tb),
                "--meth_bed_hap1",
                f"{data}/pbcpg/CASE01.bed.gz",
                "--meth_bed_hap2",
                f"{data}/modkit/CASE01.bed.gz",
                ref=ref)
    check(
        "run --target-region with a 5-column combined BED and pb-CpG-tools "
        "and modkit haplotype BEDs: target_region.tsv identical to 5-column",
        r.returncode == 0 and listing(scratch) == [] and
        same_files(o9, om, ["CASE01.target_region.tsv"]), r.stderr[-800:])
    r = run_cmd("CASE01",
                "F",
                base / "x",
                "--target-region",
                str(base / "nope.bed"),
                ref=ref)
    check("run --target-region with a missing file: clear error",
          r.returncode == 1 and "BED not found" in r.stderr, r.stderr[-300:])
    r = run_cmd("CASE01",
                "F",
                base / "run_tflag",
                "--target-region",
                "--reference",
                str(ref),
                ref=None)
    check(
        "run --target-region without a path: rejected (the flag takes a file)",
        r.returncode == 2, r.stderr[-300:])

    # ---- 9. bad input ----
    bad = base / "bad"
    bad.mkdir()
    d5 = data / "5col"
    good = f"CTRL01,M,{d5}/CTRL01.bed.gz,{d5}/CTRL01.vcf"
    hdr = "sample_id,sex,meth_bed,vcf\n"
    cases = {
        "empty": ("", "is empty"),
        "missing column": ("sample_id,sex,meth_bed\nCTRL01,M,x\n",
                           "missing required column(s): vcf"),
        "empty sample_id": (hdr + f",M,{d5}/CTRL01.bed.gz,{d5}/CTRL01.vcf\n",
                            "row 2: sample_id is empty"),
        "duplicate id":
            (hdr + f"{good}\n{good}\n", "duplicate sample_id 'CTRL01'"),
        "bad sex": (hdr + f"CTRL01,X,{d5}/CTRL01.bed.gz,{d5}/CTRL01.vcf\n",
                    "sex must be M or F, got 'X'"),
        "missing file": (hdr + f"CTRL01,M,{d5}/nope.bed.gz,{d5}/CTRL01.vcf\n",
                         "meth_bed file not found"),
    }
    for name, (text, msg) in cases.items():
        f = bad / (name.replace(" ", "_") + ".csv")
        f.write_text(text)
        r = cli(["build", "--input", str(f), "--outdir", str(bad / "out")])
        check(
            f"bad samplesheet ({name}) in build: exit 1 and message",
            r.returncode == 1 and msg in r.stderr and
            "Traceback" not in r.stderr, r.stderr[-300:])
    r = run_cmd("CASE01", "F", bad / "o0", ref=ref, g=False)
    check(
        "run without --gtf: clear error before any work, nothing created",
        r.returncode == 1 and "--gtf FILE" in r.stderr and
        "dmrdetector get-genes --release N --outdir DIR" in r.stderr and
        not (bad / "o0").exists(), r.stderr[-300:])
    r = run_cmd("CASE01",
                "F",
                bad / "o0",
                "--gtf",
                str(base / "none.gtf"),
                ref=ref,
                g=False)
    check("run with a missing GTF file: clear error", r.returncode == 1 and
          "gtf file not found" in r.stderr, r.stderr[-300:])
    r = cli([
        "run", *sample_args(data, "CASE01", "F"), "--gtf",
        str(gtf), "--outdir",
        str(bad / "o1"), *res
    ],
            cwd=bad)
    check(
        "run without a reference: clear error", r.returncode == 1 and
        "no reference" in r.stderr and not (bad / "o1").exists(),
        r.stderr[-300:])
    r = run_cmd("CASE01", "F", bad / "o3", ref=base / "nodir")
    check("run with a missing reference directory: clear error",
          r.returncode == 1 and "is not a directory" in r.stderr,
          r.stderr[-300:])
    (bad / "emptyref").mkdir()
    r = run_cmd("CASE01", "F", bad / "o4", ref=bad / "emptyref")
    check("run with an empty reference directory: clear error",
          r.returncode == 1 and "no autosomal reference atlas" in r.stderr,
          r.stderr[-300:])
    r = cli([
        "run", "--vcf",
        str(d5 / "CASE01.vcf"), "--reference",
        str(ref), "--gtf",
        str(gtf)
    ])
    check(
        "run without --meth_bed_combined: clear error", r.returncode == 1 and
        "--meth_bed_combined FILE is required" in r.stderr, r.stderr[-300:])
    r = run_cmd("CASE01", "F", bad / "o5", "--sample-id", "a/b", ref=ref)
    r2 = run_cmd("CASE01", "Q", bad / "o6", ref=ref)
    check(
        "run with a bad --sample-id or --sex: clear error",
        r.returncode == 1 and "may contain only" in r.stderr and
        r2.returncode == 1 and "--sex must be M or F" in r2.stderr,
        r.stderr[-200:] + r2.stderr[-200:])
    r = cli([
        "build", "--input",
        str(full), "--outdir",
        str(bad / "o7"), "--set", "include_sex_chromosomes=false", "--threads",
        "0"
    ])
    check("build with --threads 0: clear error", r.returncode == 1 and
          "at least 2" in r.stderr, r.stderr[-300:])

    hg19 = bad / "hg19.bed.gz"
    with gzip.open(hg19, "wt") as fh:
        fh.write("".join(
            f"chr22\t{p}\t{p + 1}\t20\t10\n" for p in (20_000_000, 30_000_000,
                                                       60_000_000)))
    r = run_cmd("CASE01",
                "F",
                bad / "o8",
                "--meth_bed_combined",
                str(hg19),
                ref=ref)
    check(
        "run with hg19-sized coordinates: ERROR naming hg38, no traceback",
        r.returncode == 1 and "ERROR:" in r.stderr and
        "hg38/GRCh38" in r.stderr and "Traceback" not in r.stderr,
        r.stderr[-300:])
    mk = bad / "modkit_h.bed.gz"

    def mrow(p, code):
        return "\t".join([
            "chr22",
            str(p),
            str(p + 1), code, "20", ".",
            str(p),
            str(p + 1), "255,0,0", "20", "50.00", "10", "10", "0", "0", "0",
            "0", "0"
        ]) + "\n"

    with gzip.open(mk, "wt") as fh:
        fh.write("".join(
            mrow(p, c)
            for p in range(20_000_000, 20_000_300, 6)
            for c in ("m", "h")))
    r = run_cmd("CASE01",
                "F",
                bad / "o9",
                "--meth_bed_combined",
                str(mk),
                ref=ref)
    check(
        "run with 5hmC rows in a modkit file: ERROR with --combine-strands, "
        "no traceback", r.returncode == 1 and "ERROR:" in r.stderr and
        "--combine-strands" in r.stderr and "Traceback" not in r.stderr,
        r.stderr[-300:])

    n_fail = sum(not ok for _, ok, _ in results)
    print(f"\n{len(results) - n_fail}/{len(results)} checks passed")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
