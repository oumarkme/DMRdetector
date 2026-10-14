"""`dmrdetector run`: DMRs of one sample against a reference from `build`.

Steps: preprocess and split, call DMRs in each scope the reference has, annotate
genes, and with --target-region the target-region report. The output folder gets
<sample_id>.dmr.tsv, <sample_id>.log and the optional report; the temporary
folder is removed at the end unless --keep-tmp (doc/output.md).
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from dataclasses import field
import functools
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
import time

from dmrdetector import __version__
from dmrdetector import data
from dmrdetector import runtime
from dmrdetector import settings


REF_SCOPE_ARGS = {
    "autosomal": "",
    "male_X": "_male_X",
    "female_X": "_female_X",
    "Y": "_Y"
}

HELP = (
    "Detect DMRs of one sample against a reference made by `dmrdetector "
    "build`. Writes <sample_id>.dmr.tsv, "
    "<sample_id>.log and, with --target-region, <sample_id>.target_region.tsv "
    "into --outdir.")
EPILOG = ("Thresholds and the other options are set in the settings file "
          "(`dmrdetector config init`) or for one run "
          "with --set; `dmrdetector config show run` lists them. `dmrdetector "
          "get-genes` downloads a GTF.")


@dataclass
class Reference:
    """Files of a reference: atlases and panels by scope, ICR panel, notes."""

    atlas: dict = field(default_factory=dict)  # scope -> path
    panel: dict = field(default_factory=dict)  # scope -> path
    icr_panel: str | None = None
    sex_branches: bool = False
    info: dict = field(default_factory=dict)  # reference.json (if any)
    warnings: list = field(default_factory=list)


def _find(dirs, names) -> str | None:
    """Returns the first existing dirs[i]/names[j] (directories first)."""
    for d in dirs:
        for n in names:
            p = os.path.join(d, n)
            if os.path.isfile(p):
                return p
    return None


def resolve_reference(refdir: str, sigma: int, include_sex: bool) -> Reference:
    """Finds the atlases and panels of a reference directory.

    Args:
        refdir: Directory written by `build`, or a v2.2 results directory
            (atlas/, cpg_islands/, icr_regions/).
        sigma: Smoothing sigma of the atlas names; reference.json wins.
        include_sex: Whether the sex-specific scopes are wanted.

    Returns:
        The Reference. The sex scopes are used only when all three atlases
        and all three CpG panels exist.

    Raises:
        runtime.RunError: No autosomal atlas or CpG panel, or a listed file is
            missing.
    """
    ref = Reference()
    if not os.path.isdir(refdir):
        raise runtime.RunError(f"--reference {refdir} is not a directory"
                               f"{runtime.container_hint()}")
    rj = os.path.join(refdir, "reference.json")
    if os.path.isfile(rj):
        with open(rj) as fh:
            ref.info = json.load(fh)
        sigma = ref.info.get("smooth_sigma", sigma)
    adirs = [refdir, os.path.join(refdir, "atlas")]
    cdirs = [refdir, os.path.join(refdir, "cpg_islands")]
    idirs = [refdir, os.path.join(refdir, "icr_regions")]
    missing = []
    for sc, suf in REF_SCOPE_ARGS.items():
        a = _find(adirs, [
            f"ref_atlas_{sigma}{suf}.bed.gz", f"ref_atlas_{sigma}{suf}.tsv.gz"
        ])
        c = _find(cdirs, [f"cpg_signal{suf}.bed"])
        if a:
            ref.atlas[sc] = a
        if c:
            ref.panel[sc] = c
        if sc != "autosomal" and not (a and c):
            missing.append(f"{sc} ({'atlas' if not a else 'CpG panel'})")
    if "autosomal" not in ref.atlas:
        raise runtime.RunError(
            f"no autosomal reference atlas: looked for ref_atlas_{sigma}.bed.gz"
            f" or .tsv.gz in {', '.join(adirs)}; pass --reference DIR (made by "
            "`dmrdetector build`)")
    if "autosomal" not in ref.panel:
        raise runtime.RunError("no autosomal CpG island panel "
                               f"(cpg_signal.bed) in {', '.join(cdirs)}")
    ref.icr_panel = _find(idirs, ["icr_region.bed"])
    if not ref.icr_panel:
        ref.warnings.append(
            "no ICR panel (icr_region.bed): ICR regions are not tested")
    full = not missing
    ref.sex_branches = include_sex and full
    if include_sex and not full:
        ref.warnings.append(f"the reference lacks: {', '.join(missing)}. "
                            "Falling back to autosomal-only DMR calling")
    if not ref.sex_branches:
        for sc in ("male_X", "female_X", "Y"):
            ref.atlas.pop(sc, None)
            ref.panel.pop(sc, None)
    for sc in ref.atlas:
        for f in (ref.atlas[sc], ref.panel.get(sc)):
            if f and not os.path.isfile(f):
                raise runtime.RunError(f"reference file not found: {f}")
    return ref


def expected_scopes(s: runtime.Sample, sex_branches: bool) -> list[str]:
    """Returns the scopes a sample is tested in."""
    if not sex_branches:
        return ["autosomal"]
    return [
        sc for sc in runtime.SCOPES if runtime.SCOPE_SEX[sc] in (None, s.sex)
    ]


def default_sample_id(meth_bed: str) -> str:
    """Returns the file name of the methylation BED up to its first '.'."""
    return os.path.basename(meth_bed).split(".", 1)[0]


NO_MATCH_CAUSES = (
    "likely causes: contig names other than chr1-chr22 and chrX/chrY (for "
    "example RefSeq NC_000001.11 instead of chr1), coordinates that are not "
    "hg38/GRCh38, or no CpG left after the coverage filter (dmr_min_coverage)")


def check_matched(sample_id: str, scope: str, stats: dict) -> str | None:
    """Judges the CpGs of a sample that match the reference atlas in a scope.

    Args:
        sample_id: Sample ID, for the message.
        scope: Scope name.
        stats: sample_sites, atlas_sites and merged, as `call_dmr` counts them
            (all inside the panel regions).

    Returns:
        None if some CpG matched, else the message with the counts.
    """
    if stats["merged"] > 0:
        return None
    return (f"{sample_id} {scope}: 0 sample CpGs match the reference atlas "
            f"(sample CpGs in the panel regions: {stats['sample_sites']}; "
            f"atlas sites in the panel regions: {stats['atlas_sites']})")


def _t_call_dmr(t: dict):
    """Pool task: DMR calling of one sample and scope.

    No match with the atlas is an error in the autosomal scope and a warning
    in a sex scope (a sex scope can be empty for a valid sample).
    """
    from dmrdetector.steps import call_dmr
    out = t["out"]
    stats: dict = {}

    def work():
        with runtime.temp_output(out) as tmp:
            df = call_dmr.call_dmr(
                t["bed"],
                t["atlas"],
                t["panel"],
                t["icr_panel"],
                scope=t["scope"],
                include_sex=t["include_sex"],
                tabix_duplicates=t["tabix_duplicates"],
                threads=t["threads"],
                log=lambda m: runtime.log(f"  {t['id']} {t['scope']}: {m}"),
                stats=stats,
                **t["thr"])
            call_dmr.write_raw(df, tmp)
            runtime.commit(tmp, out)

    key = dict(step="call_dmr", scope=t["scope"])
    done = runtime.guarded(t["stamps"], f"call_dmr_{t['id']}_{t['scope']}", key,
                           [out], work)
    msg = check_matched(t["id"], t["scope"], stats) if stats else None
    if msg and t["scope"] == "autosomal":
        raise runtime.RunError(f"{msg}; {NO_MATCH_CAUSES}")
    if msg:
        runtime.log(f"WARNING: {msg}; no DMR can be called in this scope")
    return [done]


# gene table loaded once in the parent before the fork; workers inherit it (not
# pickled per task)
_GENES = None


def _t_annotate(t: dict):
    """Pool task: gene annotation of one sample."""
    from dmrdetector.steps import annotate
    out = t["out"]

    def work():
        with runtime.temp_output(out) as tmp:
            df = annotate.annotate(sorted(t["raw"]),
                                   t["gtf"],
                                   t["icr_panel"],
                                   t["icr_gene_list"],
                                   genes=_GENES)
            annotate.write_annotated(df, tmp)
            os.replace(tmp, out)
        runtime.log(f"  {t['id']}: {len(df)} annotated DMRs")

    return [
        runtime.guarded(t["stamps"], f"annotate_{t['id']}",
                        dict(step="annotate"), [out], work)
    ]


def _publish(src: Path, dest: Path) -> None:
    """Copies a finished file to `dest` through a temporary name."""
    part = dest.with_name(f".{dest.name}.part")
    try:
        shutil.copyfile(src, part)
        os.replace(part, dest)
    finally:
        try:
            part.unlink()
        except OSError:
            pass


def _inputs(args) -> runtime.Sample:
    """Returns the validated per-sample inputs of the command line."""
    for flag, v in (("--meth_bed_combined", args.meth_bed_combined),
                    ("--vcf", args.vcf)):
        if not v:
            raise runtime.RunError(
                f"{flag} FILE is required (the merged-haplotype methylation "
                "BED and the VCF of the sample)")
    runtime.check_files([("--meth_bed_combined", args.meth_bed_combined),
                         ("--vcf", args.vcf),
                         ("--meth_bed_hap1", args.meth_bed_hap1),
                         ("--meth_bed_hap2", args.meth_bed_hap2)])
    sid = args.sample_id or default_sample_id(args.meth_bed_combined)
    if not args.sample_id and not sid:
        raise runtime.RunError(
            "cannot derive a sample ID from the file name; pass --sample-id "
            "NAME")
    runtime.check_sample_id(sid, "" if args.sample_id else "default ")
    sex = (args.sex or "").upper()
    if args.sex and sex not in ("M", "F"):
        raise runtime.RunError(f"--sex must be M or F, got '{args.sex}'")

    def ab(v):
        return os.path.abspath(v) if v else None

    return runtime.Sample(sid, sex, ab(args.meth_bed_combined), ab(args.vcf),
                          ab(args.meth_bed_hap1), ab(args.meth_bed_hap2))


def run(args: argparse.Namespace) -> int:
    """Detects the DMRs of one sample.

    Returns:
        0 on success. Errors propagate; the log is written in every case and
        the temporary folder is removed unless --keep-tmp.
    """
    global _GENES
    t_start = time.time()
    cfg, src, cfg_path = settings.effective(
        args, "run", warn=lambda m: runtime.log(f"WARNING: {m}"))
    s = _inputs(args)
    # a bad value stops here, before any folder is made
    runtime.resolve_threads(cfg["threads"])
    if not cfg["gtf"]:
        raise runtime.RunError(
            "no gene annotation: pass --gtf FILE (an Ensembl GRCh38 GTF, "
            "plain or .gz), or set `gtf` in the "
            "settings file; download one with `dmrdetector get-genes "
            "--release N --outdir DIR`")
    if not os.path.isfile(cfg["gtf"]):
        raise runtime.RunError(
            f"gtf file not found: {cfg['gtf']}{runtime.container_hint()}")
    if not cfg["reference"]:
        raise runtime.RunError(
            "no reference: pass --reference DIR (written by `dmrdetector "
            "build`), or set `reference` in the settings file")
    notices = []
    include_sex = cfg["include_sex_chromosomes"]
    if include_sex and not s.sex:
        notices.append(
            "NOTICE: --sex not given: only the autosomal scope (chr1-22 and "
            "the PARs) is analyzed; "
            "the sex-specific scopes (chrX, chrY) are skipped")
    ref = resolve_reference(cfg["reference"], cfg["smooth_sigma"],
                            include_sex and bool(s.sex))
    par_bed = cfg["par_bed"] or data.asset("par_hg38.bed")
    gene_list = cfg["icr_gene_list"] or data.asset("icr_gene_list.bed")
    runtime.check_files([("par_bed", par_bed if ref.sex_branches else None),
                         ("icr_gene_list", gene_list)])
    ref.warnings += gene_list_warnings(ref.info, gene_list)
    gtf = os.path.abspath(cfg["gtf"])
    target_bed = None
    if cfg["target_region"]:
        target_bed = (data.asset("target_region.hg38.bed")
                      if cfg["target_region"] == settings.BUILTIN_TARGET else
                      os.path.abspath(cfg["target_region"]))
        if not os.path.isfile(target_bed):
            raise runtime.RunError(
                f"--target-region BED not found: {cfg['target_region']}"
                f"{runtime.container_hint()}")
    elif s.meth_bed_hap1 or s.meth_bed_hap2:
        notices.append(
            "NOTICE: --meth_bed_hap1/--meth_bed_hap2 are used only with "
            "--target-region; ignored")

    outdir = Path(cfg["outdir"]).resolve()
    try:
        outdir.mkdir(parents=True, exist_ok=True)
        if cfg["tmp_dir"]:
            os.makedirs(cfg["tmp_dir"], exist_ok=True)
        workdir = Path(
            tempfile.mkdtemp(prefix="dmrdetector-", dir=cfg["tmp_dir"] or
                             None)).resolve()
    except OSError as e:
        raise runtime.RunError(
            f"cannot create the output or temporary folder: {e}") from None
    ctx = runtime.make_ctx("run", cfg, src, cfg_path, outdir, workdir)
    error = None
    try:
        runtime.prepare_dirs(ctx, "preprocessed", "split", "dmr_raw",
                             "annotated")
        runtime.log(f"DMRdetector {__version__} run: sample {s.sample_id}, "
                    f"scopes={','.join(ref.atlas)}, "
                    f"{runtime.budget_text(ctx)}")
        runtime.log_settings(ctx)
        for w in notices:
            runtime.log(w)
        for w in ref.warnings:
            runtime.log(f"WARNING: {w}")
        runtime.remove_stale_tmp(outdir)
        for w in replaced_outputs(outdir, s.sample_id, bool(target_bed)):
            runtime.log(f"WARNING: {w}")
        _steps(ctx, cfg, s, ref, par_bed, gene_list, gtf)
        # the DMR table is final here: a failure of the target report must not
        # discard it
        _publish(workdir / "annotated" / f"{s.sample_id}.dmr.tsv",
                 outdir / f"{s.sample_id}.dmr.tsv")
        if target_bed:
            _target_report(ctx, cfg, s, ref, gene_list, target_bed)
            _publish(workdir / "target_region.tsv",
                     outdir / f"{s.sample_id}.target_region.tsv")
        runtime.log(
            f"finished in {time.time() - t_start:.0f} s; outputs in {outdir}")
        return 0
    except BaseException as e:
        error = e
        raise
    finally:
        _GENES = None
        inputs = dict(sample_id=s.sample_id,
                      sex=s.sex or None,
                      meth_bed_combined=s.meth_bed,
                      vcf=s.vcf,
                      meth_bed_hap1=s.meth_bed_hap1,
                      meth_bed_hap2=s.meth_bed_hap2,
                      reference=cfg["reference"],
                      scopes=",".join(ref.atlas),
                      gtf=gtf,
                      target_region=target_bed,
                      outdir=outdir,
                      temporary_folder=workdir)
        try:
            runtime.write_log(outdir / f"{s.sample_id}.log", ctx, "run",
                              t_start, inputs, error)
        except OSError as e:
            print(f"WARNING: could not write the log: {e}", file=sys.stderr)
        runtime.restore_tmp()
        runtime.log.set_path(None)
        if cfg["keep_tmp"]:
            print(f"temporary files kept in {workdir}", file=sys.stderr)
        else:
            # the mkdtemp folder of this run, never --tmp-dir itself
            shutil.rmtree(workdir, ignore_errors=True)


def gene_list_warnings(info: dict, gene_list: str) -> list[str]:
    """Returns a warning when the ICR gene list is not the one of the build.

    `build` records the md5 of its gene list in reference.json
    (`icr_gene_list_md5`). The list gives the GENE labels of the ICR DMRs (the
    ICR panel itself is read from the reference). A reference without the key
    (a v2.2 layout, a download, a build before 3.0.8) is not checked.

    Args:
        info: Parsed reference.json (empty when there is none).
        gene_list: Path of the gene list this run uses.
    """
    recorded = info.get("icr_gene_list_md5")
    if not recorded or runtime.md5_file(gene_list) == recorded:
        return []
    return ["the ICR gene list differs from the one used to build the "
            "reference (GENE labels may differ)"]


def replaced_outputs(outdir: Path, sample_id: str,
                     with_target: bool) -> list[str]:
    """Returns warnings about files of `outdir` that this run replaces or keeps.

    Nothing is deleted: a target report of an earlier run stays when this run
    makes none.
    """
    out = []
    for suffix, what in ((".dmr.tsv", "DMR table"), (".log", "log")):
        if (outdir / f"{sample_id}{suffix}").exists():
            out.append(f"{sample_id}{suffix} (the {what} of an earlier run) "
                       "will be replaced")
    tr = outdir / f"{sample_id}.target_region.tsv"
    if tr.exists():
        out.append(f"{tr.name} (the target-region report of an earlier run) " +
                   ("will be replaced" if with_target else
                    "is left in place; this run writes no target-region "
                    "report, so the file does not belong to it"))
    return out


def _steps(ctx, cfg, s, ref, par_bed, gene_list, gtf) -> None:
    """Runs the steps of one sample: preprocess, call DMRs, annotate."""
    from dmrdetector.steps import preprocess
    global _GENES
    st, out = ctx.stamps, ctx.workdir
    include_sex = cfg["include_sex_chromosomes"]
    preprocess.preprocess_samples(ctx, [s], cfg["dmr_min_coverage"],
                                  ref.sex_branches, par_bed)

    dtasks, raws = [], []
    for sc in expected_scopes(s, ref.sex_branches):
        female_x = sc == "female_X"
        # the v2.2 name is dmr_raw_<id>_<type>_<scope>.tsv.gz, and the sorted
        # names fix the order in which the annotation reads the tables
        raw = str(ctx.workdir / "dmr_raw" /
                  f"dmr_raw_{s.sample_id}_case_{sc}.tsv.gz")
        raws.append(raw)
        dtasks.append(
            dict(id=s.sample_id,
                 scope=sc,
                 bed=preprocess.scope_bed(ctx, s.sample_id, sc,
                                          ref.sex_branches),
                 atlas=ref.atlas[sc],
                 panel=ref.panel[sc],
                 icr_panel=ref.icr_panel,
                 out=raw,
                 include_sex=include_sex,
                 tabix_duplicates=not cfg["dmr_dedupe_sites"],
                 stamps=st,
                 thr=dict(cpg_diff_threshold=None
                          if female_x else cfg["dmr_cpg_diff_threshold"],
                          cpg_min_sites=None
                          if female_x else cfg["dmr_cpg_min_sites"],
                          sex_f_diff_threshold=cfg["sex_f_diff_threshold"],
                          sex_f_min_sites=cfg["sex_f_min_sites"],
                          icr_diff_threshold=cfg["icr_diff_threshold"],
                          icr_min_sites=cfg["icr_min_sites"],
                          merge_distance=cfg["dmr_merge_distance"],
                          fdr_threshold=cfg["dmr_fdr_threshold"])))
    runtime.log(f"call DMRs: {len(dtasks)} scope tasks")
    t0 = time.time()
    runtime.pool_map(
        ctx,
        _t_call_dmr,
        dtasks,
        runtime.TASK_GB["call_dmr"],
        "call_dmr",
        # the autosomal tables hold 99% of the sites
        [6.0 if t["scope"] == "autosomal" else 1.0 for t in dtasks])
    runtime.log(f"DMR calling finished ({time.time() - t0:.0f} s)")

    from dmrdetector.steps import annotate
    runtime.check_fits("gene table", runtime.TASK_GB["gene_table"], ctx.mem_gb,
                       ctx.threads)
    runtime.log(f"annotate: loading genes from {gtf}")
    t0 = time.time()
    _GENES = annotate.load_genes(gtf)
    ctx.timings["load_genes"] = {
        "status": "loaded",
        "seconds": round(time.time() - t0, 2)
    }
    atasks = [
        dict(id=s.sample_id,
             raw=raws,
             gtf=gtf,
             icr_panel=ref.icr_panel,
             icr_gene_list=gene_list,
             out=str(out / "annotated" / f"{s.sample_id}.dmr.tsv"),
             stamps=st)
    ]
    runtime.pool_map(ctx, _t_annotate, atasks, runtime.TASK_GB["annotate"],
                     "annotate")
    _GENES = None


def _target_report(ctx, cfg, s, ref, gene_list, target_bed) -> None:
    """Runs the target-region report of one sample.

    A failure is logged as a warning that says the DMR table is published, and
    raised again, so that the run ends with status 1.
    """
    from dmrdetector.steps import target
    st, out = ctx.stamps, ctx.workdir
    tout = str(out / "target_region.tsv")
    rec = dict(sample_id=s.sample_id, type="CASE", merged=s.meth_bed, vcf=s.vcf)
    if cfg["target_haplotypes"]:
        for tag, p in (("hap1", s.meth_bed_hap1), ("hap2", s.meth_bed_hap2)):
            if p:
                rec[tag] = p

    def work():
        tdir = ctx.scratch / "target"
        tdir.mkdir(exist_ok=True)
        with runtime.temp_output(tout) as tmp:
            target.target_report(
                [rec],
                ref.atlas["autosomal"],
                target_bed,
                tmp,
                icr_region=ref.icr_panel,
                icr_gene_list=gene_list,
                haplotypes=cfg["target_haplotypes"],
                min_coverage=cfg["target_min_coverage"],
                hap_min_coverage=cfg["target_hap_min_coverage"],
                fdr_threshold=cfg["dmr_fdr_threshold"],
                cpg_diff_threshold=cfg["dmr_cpg_diff_threshold"],
                icr_diff_threshold=cfg["icr_diff_threshold"],
                input_format=cfg["input_format"],
                scratch=str(tdir),
                log=runtime.log,
                threads=ctx.threads)
            os.replace(tmp, tout)

    try:
        runtime.check_fits("target region report",
                           runtime.TASK_GB["target_report"], ctx.mem_gb,
                           ctx.threads)
        runtime.log("target region report")
        runtime.record(
            ctx,
            runtime.safe("target region", runtime.guarded, st, "target_region",
                         {"step": "target"}, [tout], work))
    except runtime.RunError as e:
        ctx.status_note = ("the DMR table was published; the target-region "
                           "report failed")
        old = ctx.outdir / f"{s.sample_id}.target_region.tsv"
        first = str(e).strip().splitlines()[0]
        runtime.log(
            f"WARNING: target-region report failed: {first}. "
            f"{s.sample_id}.dmr.tsv was published and is complete; the run is "
            "incomplete" + (f", and the {old.name} in the folder is from an "
                            "earlier run" if old.exists() else ""))
        raise


def add_parser(
        subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """Adds the `run` command to `subparsers` and returns its parser."""
    p = subparsers.add_parser(
        "run",
        help="detect DMRs of one sample against a reference (made by `build`)",
        description=HELP,
        epilog=EPILOG,
        formatter_class=settings.HelpFormatter,
        allow_abbrev=False,
        add_help=False,
        usage="dmrdetector run --meth_bed_combined FILE --vcf FILE [options]")
    g = p.add_argument_group("Inputs (one sample per run)")
    g.add_argument("--meth_bed_combined",
                   metavar="FILE",
                   help="merged-haplotype methylation BED (required)")
    g.add_argument(
        "--vcf",
        metavar="FILE",
        help="VCF for the het-SNP filter; need not be phased (required)")
    g.add_argument("--sex",
                   metavar="M|F",
                   help="sample sex; without it, autosomal scope only")
    g.add_argument("--sample-id",
                   metavar="NAME",
                   help="output name [default: BED name up to first '.']")
    g.add_argument("--meth_bed_hap1",
                   metavar="FILE",
                   help="haplotype 1 BED (target-region report)")
    g.add_argument("--meth_bed_hap2",
                   metavar="FILE",
                   help="haplotype 2 BED (target-region report)")
    # same options with hyphens, hidden from --help
    for k in ("combined", "hap1", "hap2"):
        p.add_argument(f"--meth-bed-{k}",
                       dest=f"meth_bed_{k}",
                       help=argparse.SUPPRESS)
    g = p.add_argument_group("Reference and annotation")
    settings.add_flag(g, "reference")
    settings.add_flag(g, "gtf")
    settings.add_flag(g, "target_region")
    g = p.add_argument_group("Output")
    settings.add_flag(g, "outdir")
    settings.add_flag(g, "tmp_dir")
    settings.add_flag(g, "keep_tmp")
    g = p.add_argument_group("Resources")
    settings.add_flag(g, "threads")
    g = p.add_argument_group("Other")
    settings.add_common(g)
    g.add_argument("-h",
                   "--help",
                   action="help",
                   help="show this help and exit")
    p.set_defaults(func=functools.partial(runtime.run_command, run))
    return p
