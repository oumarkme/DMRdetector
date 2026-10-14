"""`dmrdetector build`: a reference (atlas and panels) from control samples.

Per sample: preprocess and split by chromosome component. Per scope: raw atlas,
smoothed atlas, CpG island panel; the autosomal scope also gets the ICR panel.
Files and options are described in doc/output.md and the README.
"""

from __future__ import annotations

import argparse
from datetime import datetime
import functools
import os
from pathlib import Path
import sys
import time

from dmrdetector import __version__
from dmrdetector import data
from dmrdetector import runtime
from dmrdetector import settings


HELP = ("Build a reference (atlas, CpG island and ICR panels) from control "
        "samples. The samplesheet lists the controls: "
        "sample_id,sex,meth_bed,vcf (sex M or F). Work files go to --tmp-dir "
        "(default <outdir>/work); build makes the folder itself and refuses an "
        "existing one it did not make. They are kept after a failure so that "
        "a rerun resumes, and removed after success.")
EPILOG = ("Thresholds and the other options are set in the settings file "
          "(`dmrdetector config init`) or per run with --set.")

# below this the scopes run one after another, each with the whole budget
CONCURRENT_SCOPES_MIN_THREADS = 16


def plan_scopes(cfg: dict,
                ctrl: list[runtime.Sample]) -> tuple[dict, list[str]]:
    """Plans the scopes of the reference.

    Args:
        cfg: Effective options.
        ctrl: Control samples.

    Returns:
        ({scope: [sample_id]}, warnings). A sex scope needs
        sex_atlas_min_controls samples of that sex.
    """
    warns: list[str] = []
    scopes = {}
    for sc in (runtime.SCOPES if cfg["include_sex_chromosomes"] else
               ("autosomal",)):
        sx = runtime.SCOPE_SEX[sc]
        ids = [s.sample_id for s in ctrl if sx is None or s.sex == sx]
        if sc != "autosomal" and len(ids) < cfg["sex_atlas_min_controls"]:
            warns.append(
                f"scope {sc}: {len(ids)} CTRL samples (< "
                "sex_atlas_min_controls "
                f"{cfg['sex_atlas_min_controls']}); no reference for this scope"
            )
            continue
        scopes[sc] = ids
    return scopes, warns


def scope_plan(threads: int,
               scopes: list[str]) -> tuple[list[float], int | None]:
    """Returns (thread weights, max concurrent scopes).

    The autosomal build and smoothing (about 850 CPU s against about 45 for a
    sex scope) are the critical path. With 16 or more threads the scopes run
    side by side: the sex scopes get one thread each (their single-chromosome
    FFT gains nothing from 2) and the autosomal scope the rest of the budget.
    With fewer threads a sex scope running beside it would take threads the
    autosomal scope needs for most of its run, so the scopes run in turn,
    autosomal first, each with all threads.
    """
    if threads >= CONCURRENT_SCOPES_MIN_THREADS:
        return [100.0 if sc == "autosomal" else 1.0 for sc in scopes], None
    return [1.0] * len(scopes), 1


def scope_need_gb(
        beds: list[str],
        scope: str,
        min_samples: int = settings.ATLAS_MIN_SAMPLES_PER_SITE) -> float:
    """Returns the least memory (GB) with which one scope worker can run.

    The largest of: the raw atlas on the smallest window
    (`runtime.raw_atlas_min_gb`), the smoothing of the longest chromosome with
    the loaded raw atlas arrays, and a panel. The arrays follow the number of
    sites (at most one row per `min_samples` input rows and one per
    RAW_BP_PER_SITE_MIN bp), not the number of samples.
    """
    size = sum(os.path.getsize(b) for b in beds)
    extent = runtime.SCOPE_EXTENT[scope]
    sites = min(size / runtime.RAW_ROW_FILE_BYTES / max(1, min_samples),
                extent / runtime.RAW_BP_PER_SITE_MIN)
    held_gb = runtime.SMOOTH_HELD_BYTES_PER_SITE * sites / 1e9
    return max(
        runtime.raw_atlas_min_gb(size, extent, min_samples),
        runtime.smooth_gb(runtime.SCOPE_LARGEST_SPAN[scope]) + held_gb +
        runtime.SMOOTH_PARENT_BASE_GB, runtime.TASK_GB["panel"])


def scope_memory(needs: list[float], workers: int,
                 mem_gb: float) -> list[float]:
    """Returns the memory share of every scope worker.

    Each gets its estimate; the one with the largest estimate (the autosomal
    scope, whose chromosomes are smoothed side by side) also gets whatever the
    other `workers - 1` scopes that can run beside it leave of the budget.
    """
    out = list(needs)
    top = max(range(len(needs)), key=needs.__getitem__)
    others = sorted((n for i, n in enumerate(needs) if i != top),
                    reverse=True)[:max(0, workers - 1)]
    out[top] = max(needs[top], mem_gb - sum(others))
    return out


def _panel_steps(t: dict) -> list:
    """Runs the CpG island panel of one scope, or the ICR panel."""
    from dmrdetector.fileio import atlas_io
    out = t["out"]
    if t["kind"] == "cpg":
        from dmrdetector.steps import cpg_islands

        def work():
            with runtime.temp_output(out) as tmp:
                res = cpg_islands.detect(
                    atlas_io.read_atlas(t["atlas"]),
                    t["scope"],
                    include_sex=t["include_sex"],
                    min_samples=t["min_samples"],
                    core_max_meth=t["core_max_meth"],
                    core_max_phi=t["core_max_phi"],
                    shore_max_meth=t["shore_max_meth"],
                    female_min_meth=t["female_min"],
                    female_max_meth=t["female_max"],
                    max_gap=t["max_gap"],
                    min_length=t["min_length"],
                    log=lambda m: runtime.log(f"  cpg {t['scope']}: {m.strip()}"
                                             ))
                res.to_csv(tmp, sep="\t", index=False)
                runtime.commit(tmp, out)

        name = f"cpg_islands_{t['scope']}"
    else:
        from dmrdetector.steps import icr

        def work():
            with runtime.temp_output(out) as tmp:
                cand = icr.read_candidates(t["gene_list"])
                res = icr.detect(
                    atlas_io.read_atlas(t["atlas"]),
                    cand,
                    ref_lo=t["ref_lo"],
                    ref_hi=t["ref_hi"],
                    min_sites=t["min_sites"],
                    merge_distance=t["merge_distance"],
                    min_samples=t["min_samples"],
                    v1_duplicates=t["v1_duplicates"],
                    log=lambda m: runtime.log(f"  icr: {m.strip()}"))
                res.to_csv(tmp, sep="\t", index=False)
                runtime.commit(tmp, out)

        name = "icr_regions"
    key = dict(t["params"], atlas=runtime.stat_of(t["atlas"]))
    if t["kind"] == "icr":
        key["gene_list"] = runtime.stat_of(t["gene_list"])
    return [runtime.guarded(t["stamps"], name, key, [out], work)]


def _t_scope(t: dict) -> list:
    """Runs everything one scope needs after preprocessing, in one worker.

    Raw atlas, smoothed atlas with its tabix index, CpG island panel and, for
    the autosomal scope, the ICR panel. Scopes are independent, so the small sex
    scopes run next to the autosomal FFTs (bound by chr1).
    """
    from dmrdetector.steps import atlas
    sc, st, threads = t["scope"], t["stamps"], t["threads"]
    raw, final, include_sex = t["raw"], t["final"], t["include_sex"]
    steps = []

    def build():
        """Builds the raw atlas of the scope into `raw`."""
        with runtime.temp_output(raw) as tmp:
            atlas.build_atlas(t["beds"],
                              tmp,
                              t["min_samples_atlas"],
                              t["shrinkage"],
                              threads=threads,
                              scratch=t["scratch"],
                              log=lambda m: runtime.log(f"  atlas {sc}: {m}"),
                              mem_gb=t.get("mem_gb"))
            os.replace(tmp, raw)

    key = dict(step="build_atlas",
               scope=sc,
               min_samples=t["min_samples_atlas"],
               w=t["shrinkage"],
               beds=[runtime.stat_of(b) for b in t["beds"]])
    runtime.log(f"atlas {sc}: {len(t['beds'])} CTRL samples, {threads} threads")
    steps.append(runtime.guarded(st, f"build_atlas_{sc}", key, [raw], build))

    def smooth():
        """Smooths the raw atlas into the reference atlas `final`."""
        with runtime.temp_output(final) as tmp:
            atlas.smooth_atlas(raw,
                               t["sigma"],
                               tmp,
                               include_sex,
                               threads,
                               log=lambda m: runtime.log(f"  smooth {sc}: {m}"),
                               mem_gb=t.get("mem_gb"))
            # the atlas and its .tbi move together; the index is never rebuilt
            runtime.commit(tmp, final)

    key = dict(step="smooth",
               scope=sc,
               sigma=t["sigma"],
               sex=include_sex,
               raw=runtime.stat_of(raw))
    steps.append(
        runtime.guarded(st, f"smooth_atlas_{sc}", key, [final, final + ".tbi"],
                        smooth))
    for pt in t["panels"]:
        steps += _panel_steps(dict(pt, atlas=final, stamps=st))
    return steps


def run(args: argparse.Namespace) -> int:
    """Builds the reference of `args.input` into `args.outdir`.

    Returns:
        0 on success. Errors propagate; the log is written in every case.
    """
    cfg, src, cfg_path = settings.effective(
        args, "build", warn=lambda m: runtime.log(f"WARNING: {m}"))
    ctrl = runtime.parse_samplesheet(args.input)
    include_sex = cfg["include_sex_chromosomes"]
    if include_sex:
        runtime.require_sex(ctrl, "include_sex_chromosomes is true")
    par_bed = cfg["par_bed"] or data.asset("par_hg38.bed")
    gene_list = cfg["icr_gene_list"] or data.asset("icr_gene_list.bed")
    runtime.check_files([("par_bed", par_bed if include_sex else None),
                         ("icr_gene_list", gene_list)])
    scopes, warns = plan_scopes(cfg, ctrl)
    outdir = Path(args.outdir).resolve()
    workdir = runtime.build_workdir(cfg["tmp_dir"], outdir)
    ctx = runtime.make_ctx("build", cfg, src, cfg_path, outdir, workdir)
    outdir.mkdir(parents=True, exist_ok=True)
    # a folder that this build did not make is refused before it is written to
    claim = runtime.claim_workdir(workdir, outdir)
    overwrite, error = False, None
    try:
        overwrite = _holds_reference(outdir) and not claim.resumed
        runtime.prepare_dirs(ctx, "preprocessed", "split", "raw_atlas")
        runtime.remove_stale_tmp(outdir)
    except BaseException:
        runtime.restore_tmp()
        runtime.log.set_path(None)
        runtime.release_workdir(claim, False, True)
        raise
    t_start = time.time()
    runtime.log(f"DMRdetector {__version__} build: {len(ctrl)} control "
                f"samples, scopes={','.join(scopes)}, "
                f"{runtime.budget_text(ctx)}")
    runtime.log_settings(ctx)
    for w in warns:
        runtime.log(f"WARNING: {w}")
    if overwrite:
        runtime.log(f"WARNING: {outdir} already holds a reference "
                    "(reference.json or ref_atlas files); its files are "
                    "overwritten")
    try:
        _steps(ctx, args, cfg, ctrl, scopes, par_bed, gene_list, cfg_path)
        runtime.log(
            f"finished in {time.time() - t_start:.0f} s; reference in {outdir}")
        return 0
    except BaseException as e:
        error = e
        raise
    finally:
        runtime.write_log(outdir / "build.log",
                          ctx,
                          "build",
                          t_start,
                          dict(input=os.path.abspath(args.input),
                               outdir=outdir,
                               work_dir=workdir,
                               control_samples=len(ctrl),
                               scopes=",".join(scopes)),
                          error,
                          files=[
                              f for s in ctrl
                              for f in (s.meth_bed, s.vcf, s.meth_bed_hap1,
                                        s.meth_bed_hap2)
                          ])
        runtime.restore_tmp()
        runtime.log.set_path(None)
        # the folder is removed only if its marker names this outdir
        removed = runtime.release_workdir(claim, error is None, cfg["keep_tmp"])
        if not removed:
            print(f"work files kept in {workdir}", file=sys.stderr)


def _reference_parameters(cfg: dict) -> dict:
    """Returns the options recorded in reference.json.

    With log_input_paths = name (the default) a path option is recorded by its
    last component, as in the logs; with full, as given.
    """
    by_name = cfg.get("log_input_paths", "name") == "name"
    out = {}
    for k, v in cfg.items():
        if k in ("threads", "tmp_dir", "keep_tmp", "log_input_paths"):
            continue
        if by_name and v and k in runtime.PATH_OPTIONS:
            v = os.path.basename(str(v).rstrip(os.sep))
        out[k] = v
    return out


def _holds_reference(outdir: Path) -> bool:
    """Returns True when `outdir` already holds a reference of an earlier build.

    That is a reference.json or an atlas file of any smoothing sigma.
    """
    return (outdir / "reference.json").is_file() or any(
        outdir.glob("ref_atlas_*.bed.gz"))


def _steps(ctx, args, cfg, ctrl, scopes, par_bed, gene_list, cfg_path) -> None:
    from dmrdetector.steps import preprocess
    out, st = ctx.outdir, ctx.stamps
    include_sex = cfg["include_sex_chromosomes"]
    min_cov = cfg["atlas_min_read_coverage"]
    preprocess.preprocess_samples(ctx, ctrl, min_cov, include_sex, par_bed)

    sigma = cfg["smooth_sigma"]
    cpg_params = dict(min_samples=cfg["cpg_min_samples_per_site"],
                      core_max_meth=cfg["cpg_core_max_methylation"],
                      core_max_phi=cfg["cpg_core_max_dispersion"],
                      shore_max_meth=cfg["cpg_shore_max_methylation"],
                      female_min=cfg["sex_f_ref_meth_min"],
                      female_max=cfg["sex_f_ref_meth_max"],
                      max_gap=cfg["cpg_max_gap_distance"],
                      min_length=cfg["cpg_min_island_length"])
    icr_params = dict(ref_lo=cfg["icr_ref_meth_min"],
                      ref_hi=cfg["icr_ref_meth_max"],
                      min_sites=cfg["icr_min_sites"],
                      merge_distance=cfg["dmr_merge_distance"],
                      min_samples=cfg["cpg_min_samples_per_site"],
                      v1_duplicates=cfg["icr_v1_duplicates"])
    atlas: dict[str, str] = {}
    tasks = []
    for sc, ids in scopes.items():
        atlas[sc] = str(out / f"ref_atlas_{sigma}{runtime.suffix(sc)}.bed.gz")
        scratch = ctx.scratch / f"atlas{runtime.suffix(sc)}"
        scratch.mkdir(exist_ok=True)
        panels = [
            dict(kind="cpg",
                 scope=sc,
                 out=str(out / f"cpg_signal{runtime.suffix(sc)}.bed"),
                 include_sex=include_sex,
                 params=dict(cpg_params, scope=sc, include_sex=include_sex),
                 **cpg_params)
        ]
        if sc == "autosomal":
            panels.append(
                dict(kind="icr",
                     out=str(out / "icr_region.bed"),
                     gene_list=gene_list,
                     params=dict(icr_params, step="icr"),
                     **icr_params))
        tasks.append(
            dict(scope=sc,
                 beds=[
                     preprocess.scope_bed(ctx, i, sc, include_sex) for i in ids
                 ],
                 raw=str(ctx.workdir / "raw_atlas" /
                         f"raw_atlas{runtime.suffix(sc)}.bed.gz"),
                 final=atlas[sc],
                 scratch=str(scratch),
                 sigma=sigma,
                 include_sex=include_sex,
                 stamps=st,
                 min_samples_atlas=cfg["atlas_min_samples_per_site"],
                 shrinkage=cfg["atlas_shrinkage_weight"],
                 panels=panels))
    needs = [
        scope_need_gb(t["beds"], t["scope"], cfg["atlas_min_samples_per_site"])
        for t in tasks
    ]
    weights, max_workers = scope_plan(ctx.threads, [t["scope"] for t in tasks])
    workers, _ = runtime.thread_plan(ctx, needs, weights, max_workers, "atlas")
    tasks = [
        dict(t, mem_gb=m)
        for t, m in zip(tasks, scope_memory(needs, workers, ctx.mem_gb))
    ]
    runtime.log(f"atlas, smoothing and panels: {len(tasks)} scopes " +
                ("in parallel" if max_workers is None else "one after another"))
    runtime.pool_map(ctx, _t_scope, tasks, needs, "atlas", weights, max_workers)

    # reference.json
    ref_scopes = {}
    for sc, ids in scopes.items():
        cpg = out / f"cpg_signal{runtime.suffix(sc)}.bed"
        ref_scopes[sc] = dict(atlas=os.path.basename(atlas[sc]),
                              cpg_panel=cpg.name,
                              n_ctrl=len(ids),
                              n_sites=runtime.count_lines(
                                  atlas[sc], ctx.threads),
                              n_islands=max(0,
                                            runtime.count_lines(str(cpg)) - 1))
    ref = dict(version=__version__,
               created=datetime.now().isoformat(timespec="seconds"),
               smooth_sigma=sigma,
               include_sex_chromosomes=include_sex,
               scopes=ref_scopes,
               icr_panel="icr_region.bed",
               n_icr_regions=max(
                   0,
                   runtime.count_lines(str(out / "icr_region.bed")) - 1),
               icr_gene_list_md5=runtime.md5_file(gene_list),
               samplesheet=os.path.basename(args.input),
               samplesheet_md5=runtime.md5_file(args.input),
               samples=[dict(sample_id=s.sample_id, sex=s.sex) for s in ctrl],
               input_format=cfg["input_format"],
               parameters=_reference_parameters(cfg),
               config_file=os.path.basename(cfg_path) if cfg_path else None)
    runtime.write_json(out / "reference.json", ref)
    runtime.log(f"reference written to {out}: " + ", ".join(
        f"{k} {v['n_sites']} sites" for k, v in ref_scopes.items()))


def add_parser(
        subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """Adds the `build` command to `subparsers` and returns its parser."""
    p = subparsers.add_parser(
        "build",
        help="build a reference (atlas and panels) from control samples",
        description=HELP,
        epilog=EPILOG,
        formatter_class=settings.HelpFormatter,
        allow_abbrev=False,
        add_help=False,
        usage="dmrdetector build --input FILE --outdir DIR [options]")
    g = p.add_argument_group("Inputs")
    g.add_argument("--input",
                   required=True,
                   metavar="FILE",
                   help="control samplesheet: sample_id,sex,meth_bed,vcf")
    g = p.add_argument_group("Output")
    g.add_argument("--outdir",
                   required=True,
                   metavar="DIR",
                   help="reference directory to write")
    settings.add_flag(g, "tmp_dir",
                      "work files folder [default: <outdir>/work]")
    settings.add_flag(g, "keep_tmp", "keep the work files after success")
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
