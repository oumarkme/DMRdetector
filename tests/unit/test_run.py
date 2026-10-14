"""Unit tests of the support code of `build` and `run`.

They cover the settings file, --set, dedicated flags, samplesheet, resources,
scopes, reference resolution, and the output and temporary folders of `run` (the
steps are stubbed). tests/smoke/run_smoke.py covers the end to end behavior.
"""

import fcntl
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

import pytest

from dmrdetector import build
from dmrdetector import cli
from dmrdetector import data
from dmrdetector import run
from dmrdetector import runtime
from dmrdetector import settings


ROOT = Path(__file__).resolve().parents[2]
DEPRECATED = "dmrdetector.toml is deprecated; rename it to settings.toml"
# v2.2 parameter name (nextflow.config) -> option key
V1_PARAMS = {
    o.key: o.key
    for o in settings.OPTIONS
    if o.key not in ("reference", "gtf", "target_region", "outdir", "keep_tmp",
                     "icr_gene_list", "par_bed", "tmp_dir", "input_format",
                     "smooth_sigma", "target_haplotypes", "target_min_coverage",
                     "target_hap_min_coverage", "threads", "icr_v1_duplicates",
                     "dmr_dedupe_sites", "log_input_paths")
}
V1_PARAMS["smooth_sigma_narrow"] = "smooth_sigma"


def _nextflow_params() -> dict:
    """Returns the numeric and boolean `params` of the v1 nextflow.config.

    The values are frozen in a fixture (v1 tag v2.2-manuscript).
    """
    return json.loads(
        (Path(__file__).parent / "fixtures" / "v1_params.json").read_text())


def test_defaults_match_v1_nextflow_config():
    nf = _nextflow_params()
    assert len(nf) > 20
    for v1_name, key in V1_PARAMS.items():
        if v1_name in nf:
            assert settings.BY_KEY[key].default == nf[v1_name], v1_name
    # every numeric v1 parameter is represented (smooth_sigma_wide: the wide
    # atlas is not built in v3)
    assert set(nf) - set(V1_PARAMS) <= {"smooth_sigma_wide", "help"}


def _write(tmp_path, text, name="c.toml"):
    p = tmp_path / name
    p.write_text(text)
    return str(p)


@pytest.mark.parametrize("text,msg", [
    ("[dmr]\ndmr_fdr_treshold = 0.1\n", "did you mean 'dmr_fdr_threshold'"),
    ("[dmr]\ndmr_fdr_threshold = 'x'\n", "must be float"),
    ("[dmr]\ndmr_fdr_threshold = 2\n", "must be between 0 and 1"),
    ("[dmr]\ndmr_cpg_min_sites = 1.5\n", "must be int"),
    ("[sex]\ninclude_sex_chromosomes = 1\n", "must be bool"),
    ("[atlas]\ndmr_fdr_threshold = 0.1\n", "belongs to [dmr]"),
    ("[nope]\nx = 1\n", "unknown section"),
    ("[input]\ninput_format = 'bam'\n", "must be one of"),
    ("dmr = 3\n", "must be a [section] table"),
    ("[dmr\n", "not valid TOML"),
])
def test_config_errors(tmp_path, text, msg):
    with pytest.raises(settings.ConfigError, match=re.escape(msg)):
        settings.load_file(_write(tmp_path, text))


def test_config_int_accepted_for_float(tmp_path):
    assert settings.load_file(_write(tmp_path,
                                     "[dmr]\ndmr_fdr_threshold = 1\n")) == {
                                         "dmr_fdr_threshold": 1.0
                                     }


def _args(command, argv):
    """Returns the parsed arguments of `dmrdetector <command> ...`.

    The required options are filled in.
    """
    req = {"build": ["--input", "s.csv", "--outdir", "o"], "run": []}[command]
    return cli.build_parser().parse_args([command, *req, *argv])


def test_precedence_flag_file_default(tmp_path):
    f = _write(
        tmp_path, "[dmr]\ndmr_fdr_threshold = 0.02\ndmr_cpg_min_sites = "
        "12\n[resources]\nthreads = 6\n")
    vals, src, path = settings.effective(
        _args(
            "run",
            ["--config", f, "--set", "dmr_fdr_threshold=0.5", "--threads", "8"
            ]), "run")
    assert (vals["dmr_fdr_threshold"], src["dmr_fdr_threshold"]) == (0.5,
                                                                     "flag")
    assert (vals["threads"], src["threads"]) == (8, "flag"
                                                )  # dedicated flag > file
    assert (vals["dmr_cpg_min_sites"], src["dmr_cpg_min_sites"]) == (12, "file")
    assert (vals["dmr_merge_distance"],
            src["dmr_merge_distance"]) == (150, "default")
    assert path == f
    assert "atlas_min_samples_per_site" not in vals  # build-only option


def test_file_sets_the_per_installation_defaults(tmp_path):
    f = _write(
        tmp_path, '[reference]\nreference = "ref"\ngtf = '
        '"g.gtf.gz"\n[target]\ntarget_region = "t.bed"\n'
        '[paths]\noutdir = "o"\ntmp_dir = "/fast"\nicr_gene_list = '
        '"l.bed"\npar_bed = "p.bed"\n'
        '[resources]\nthreads = 6\nkeep_tmp = true\n[input]\ninput_format = '
        '"modkit"\n')
    vals, src, _ = settings.effective(_args("run", ["--config", f]), "run")
    assert {
        k: vals[k] for k in ("reference", "gtf", "target_region", "outdir",
                            "tmp_dir", "icr_gene_list", "par_bed", "threads",
                            "keep_tmp", "input_format")
    } == dict(reference="ref",
              gtf="g.gtf.gz",
              target_region="t.bed",
              outdir="o",
              tmp_dir="/fast",
              icr_gene_list="l.bed",
              par_bed="p.bed",
              threads=6,
              keep_tmp=True,
              input_format="modkit")
    assert all(src[k] == "file" for k in ("reference", "gtf", "target_region",
                                          "threads", "keep_tmp"))
    vals, src, _ = settings.effective(
        _args("run", [
            "--config", f, "--reference", "other", "--gtf", "x.gtf",
            "--target-region", "u.bed"
        ]), "run")
    assert (vals["reference"], vals["gtf"], vals["target_region"],
            src["gtf"]) == ("other", "x.gtf", "u.bed", "flag")
    vals, src, _ = settings.effective(
        _args("run",
              ["--config", _write(tmp_path, "", "empty.toml")]), "run")
    assert vals["outdir"] == "dmrdetector" and vals[
        "reference"] is None and vals["keep_tmp"] is False


def test_keys_of_the_other_command_in_the_file_are_ignored(tmp_path):
    f = _write(
        tmp_path,
        "[atlas]\natlas_shrinkage_weight = 0.2\n[dmr]\ndmr_fdr_threshold = "
        "0.05\n")
    vals, _, _ = settings.effective(_args("build", ["--config", f]), "build")
    assert vals[
        "atlas_shrinkage_weight"] == 0.2 and "dmr_fdr_threshold" not in vals
    vals, _, _ = settings.effective(_args("run", ["--config", f]), "run")
    assert vals[
        "dmr_fdr_threshold"] == 0.05 and "atlas_shrinkage_weight" not in vals


def test_default_file_in_cwd(tmp_path):
    (tmp_path / "settings.toml").write_text("[atlas]\nsmooth_sigma = 200\n")
    warnings = []
    vals, src, path = settings.effective(_args("build", []),
                                         "build",
                                         cwd=str(tmp_path),
                                         warn=warnings.append)
    assert vals["smooth_sigma"] == 200 and src[
        "smooth_sigma"] == "file" and path.endswith(
            "settings.toml") and not warnings
    vals, src, path = settings.effective(_args("build", []),
                                         "build",
                                         cwd=str(tmp_path / "nowhere"),
                                         warn=warnings.append)
    assert vals["smooth_sigma"] == 150 and path is None and not warnings


def test_legacy_file_name_is_read_with_a_warning(tmp_path, capsys):
    (tmp_path / "dmrdetector.toml").write_text("[atlas]\nsmooth_sigma = 200\n")
    warnings = []
    vals, src, path = settings.effective(_args("build", []),
                                         "build",
                                         cwd=str(tmp_path),
                                         warn=warnings.append)
    assert vals["smooth_sigma"] == 200 and src[
        "smooth_sigma"] == "file" and path.endswith("dmrdetector.toml")
    assert len(warnings) == 1
    assert DEPRECATED in warnings[0]
    settings.effective(_args("build", []), "build",
                       cwd=str(tmp_path))  # default: stderr
    assert "WARNING: " + DEPRECATED in capsys.readouterr().err


def test_settings_toml_wins_over_the_legacy_file(tmp_path):
    (tmp_path / "dmrdetector.toml").write_text("[atlas]\nsmooth_sigma = 200\n")
    (tmp_path / "settings.toml").write_text("[atlas]\nsmooth_sigma = 300\n")
    warnings = []
    vals, _, path = settings.effective(_args("build", []),
                                       "build",
                                       cwd=str(tmp_path),
                                       warn=warnings.append)
    assert vals["smooth_sigma"] == 300 and path.endswith("settings.toml")
    assert len(warnings) == 1
    assert DEPRECATED in warnings[0]
    assert "settings.toml is used" in warnings[0]
    assert "dmrdetector.toml is ignored" in warnings[0]


def test_explicit_config_does_not_look_at_the_working_directory(tmp_path):
    (tmp_path / "dmrdetector.toml").write_text("[atlas]\nsmooth_sigma = 200\n")
    f = _write(tmp_path, "[atlas]\nsmooth_sigma = 250\n", "mine.toml")
    warnings = []
    vals, _, path = settings.effective(_args("build", ["--config", f]),
                                       "build",
                                       cwd=str(tmp_path),
                                       warn=warnings.append)
    assert vals["smooth_sigma"] == 250 and path == f and not warnings


def test_root_settings_toml_is_the_generated_text():
    """Checks that the root settings.toml is the generated text.

    The template is never edited by hand. Regenerate it with `python -m
    dmrdetector config init --out settings.toml --force` (repository root,
    PYTHONPATH=src).
    """
    assert (ROOT / "settings.toml").read_text() == settings.default_toml(
    ), "stale: regenerate settings.toml with `config init`"


def test_settings_toml_lists_every_default_and_comments_the_rest(tmp_path):
    """Checks the layout and the content of settings.toml.

    Every default is a live line with one comment line above it; keys without a
    default are commented out. The shipped copy loads, has the types of the
    defaults and the v1 values, and changes nothing when read.
    """
    import tomllib
    text = settings.default_toml()
    parsed = tomllib.loads(text)
    live = {k: v for section in parsed.values() for k, v in section.items()}
    with_default = {
        o.key: o.default for o in settings.OPTIONS if o.default is not None
    }
    assert live == with_default and list(parsed) == settings.SECTIONS
    for o in settings.OPTIONS:
        if o.default is None:
            assert f"# {o.key} = " in text and o.key not in live
        else:
            settings._check_value(o.key, o, o.default)
            lines = text.splitlines()
            idx = next(
                i for i, ln in enumerate(lines) if ln.startswith(f"{o.key} = "))
            assert lines[idx - 1].startswith("# ")
    assert 'flag > --set' in text
    assert 'Delete a line' in text
    assert 'settings.toml' in text
    shipped = tomllib.loads((ROOT / "settings.toml").read_text())
    live = {k: v for section in shipped.values() for k, v in section.items()}
    assert all(
        type(settings.BY_KEY[k].default) is type(v) for k, v in live.items())
    nf = _nextflow_params()
    assert all(live[key] == nf[v1] for v1, key in V1_PARAMS.items() if v1 in nf)
    for command in ("build", "run"):
        vals, _, path = settings.effective(
            _args(command, ["--config", str(ROOT / "settings.toml")]), command)
        base, _, _ = settings.effective(_args(command, []),
                                        command,
                                        cwd=str(tmp_path))
        assert vals == base and path == str(ROOT / "settings.toml")


@pytest.mark.parametrize(
    "item,key,value",
    [
        ("dmr_fdr_threshold=0.05", "dmr_fdr_threshold", 0.05),
        ("dmr_fdr_threshold=1e-3", "dmr_fdr_threshold", 0.001),
        ("dmr_cpg_min_sites=12", "dmr_cpg_min_sites", 12),
        ("dmr_cpg_diff_threshold=1", "dmr_cpg_diff_threshold",
         1.0),  # an int is accepted for a float
        ("dmr_dedupe_sites=true", "dmr_dedupe_sites", True),
        ("include_sex_chromosomes=false", "include_sex_chromosomes", False),
        ("input_format=modkit", "input_format",
         "modkit"),  # bare word for a text option
        ('input_format="pbcpg"', "input_format", "pbcpg"),
        ("gtf=/some/dir/genes.gtf.gz", "gtf", "/some/dir/genes.gtf.gz"),
        ("dmr-fdr-threshold=0.2", "dmr_fdr_threshold", 0.2),  # hyphens accepted
        (" dmr_fdr_threshold = 0.2 ", "dmr_fdr_threshold", 0.2),
    ])
def test_set_parsing(item, key, value):
    got = settings.parse_set([item], "run")
    assert got == {key: value} and type(got[key]) is type(value)


@pytest.mark.parametrize(
    "item,command,msg",
    [
        ("dmr_fdr_treshold=0.1", "run", "did you mean 'dmr_fdr_threshold'"),
        ("nonsense=1", "run", "unknown option"),
        ("atlas_shrinkage_weight=0.2", "run", "applies to build only"),
        ("dmr_fdr_threshold=0.1", "build", "applies to run only"),
        ("dmr_fdr_threshold=high", "run", "must be float"),
        ("dmr_cpg_min_sites=1.5", "run", "must be int"),
        ("dmr_dedupe_sites=maybe", "run", "must be bool"),
        ("dmr_fdr_threshold=1.5", "run", "must be between 0 and 1"),
        ("dmr_cpg_min_sites=-1", "run", "must be between 0 and inf"),
        ("input_format=bam", "run", "must be one of"),
        ("dmr_fdr_threshold", "run", "expected KEY=VALUE"),
        ("=0.1", "run", "expected KEY=VALUE"),
        ("sex=M", "run", "unknown option"),  # per-sample inputs are flags only
    ])
def test_set_errors(item, command, msg):
    with pytest.raises(settings.ConfigError, match=re.escape(msg)):
        settings.parse_set([item], command)


def test_set_twice_and_with_the_dedicated_flag_is_an_error():
    with pytest.raises(settings.ConfigError, match="given twice"):
        settings.parse_set(["dmr_fdr_threshold=0.1", "dmr_fdr_threshold=0.2"],
                           "run")
    with pytest.raises(settings.ConfigError, match="both given"):
        settings.effective(
            _args("run", ["--threads", "4", "--set", "threads=6"]), "run")


def test_set_reaches_every_option_of_the_command():
    for command in ("build", "run"):
        for o in settings.OPTIONS:
            if command in o.commands and o.default is not None:
                value = o.default
                if isinstance(value, bool):
                    value = "true" if value else "false"
                text = f"{o.key}={value}"
                assert settings.parse_set([text], command)[o.key] == o.default


def test_removed_flags_are_rejected(capsys):
    p = cli.build_parser()
    for argv in (
        ["run", "--dmr-fdr-threshold", "0.1"],
        ["run", "--input", "s.csv"],
        ["run", "--include-ctrl"],
        ["run", "--ref-atlas", "a"],
        ["run", "--cpg-panel", "a"],
        ["run", "--icr-panel", "a"],
        ["run", "--force"],
        ["run", "--dry-run"],
        ["run", "--input-format", "5col"],
        ["run", "--no-icr-v1-duplicates"],
        ["run", "--smooth-sigma", "150"],
        ["run", "--include-sex-chromosomes"],
        ["run", "--dmr-dedupe-sites"],
        ["run", "--ref"],
        ["run", "--target"],  # no abbreviations
        ["build", "--input", "s", "--outdir", "o", "--force"],
        ["build", "--input", "s", "--outdir", "o", "--dry-run"],
        ["build", "--input", "s", "--outdir", "o", "--icr-min-sites", "5"],
        [
            "build", "--input", "s", "--outdir", "o",
            "--no-include-sex-chromosomes"
        ]):
        with pytest.raises(SystemExit) as e:
            p.parse_args(argv)
        assert e.value.code == 2, argv
    capsys.readouterr()


def test_both_spellings_of_the_per_sample_flags():
    a = _args("run", [
        "--meth_bed_combined", "a", "--meth-bed-hap1", "b", "--meth_bed_hap2",
        "c", "--vcf", "v", "--sex", "F"
    ])
    assert (a.meth_bed_combined, a.meth_bed_hap1, a.meth_bed_hap2, a.vcf,
            a.sex) == ("a", "b", "c", "v", "F")
    assert _args("run", ["--meth-bed-combined", "x"]).meth_bed_combined == "x"
    assert _args("run", ["--target-region", "t.bed"]).target_region == "t.bed"
    with pytest.raises(SystemExit):
        _args("run", ["--target-region"])


def _help(capsys, *argv):
    with pytest.raises(SystemExit):
        cli.build_parser().parse_args([*argv, "--help"])
    return capsys.readouterr().out


def test_stage_commands_are_marked_as_advanced(capsys):
    assert "Advanced command for scripting" in _help(capsys, "call-dmr")


def test_run_and_build_help_are_grouped_and_short(capsys):
    for cmd, groups in (("run", ("Inputs", "Reference and annotation", "Output",
                                 "Resources", "Other")),
                        ("build", ("Inputs", "Output", "Resources", "Other"))):
        out = _help(capsys, cmd)
        for g in groups:
            assert re.search(rf"^{g}", out, re.M), (cmd, g)
        flags = [ln for ln in out.splitlines() if re.match(r"\s+(--|-h)", ln)]
        assert len(flags) <= 25, (cmd, len(flags))
        assert max(len(ln) for ln in out.splitlines()
                  ) <= 80, cmd  # fits an 80-column terminal, one line per flag
        assert '--dmr-fdr-threshold' not in out
        assert '--smooth-sigma' not in out
        assert '--set KEY=VALUE' in out


def test_config_show_lists_value_source_and_meaning(tmp_path, capsys):
    f = _write(tmp_path, "[dmr]\ndmr_fdr_threshold = 0.02\n")
    assert cli.main([
        "config", "show", "run", "--config", f, "--set", "dmr_cpg_min_sites=7"
    ]) == 0
    out = capsys.readouterr().out
    assert '[dmr]' in out
    assert 'dmr_fdr_threshold = 0.02  # file: report DMRs with FDR' in out
    assert 'dmr_cpg_min_sites = 7  # flag:' in out
    assert 'dmr_merge_distance = 150  # default:' in out
    assert 'reference = (unset)  # default:' in out
    assert 'atlas_min_samples_per_site' not in out
    assert cli.main(
        ["config", "show", "build", "--set", "dmr_fdr_threshold=0.1"]) == 1


def _sheet(tmp_path, rows, header="sample_id,sex,meth_bed,vcf"):
    for n in ("a.bed.gz", "a.vcf"):
        (tmp_path / n).write_text("x")
    p = tmp_path / "s.csv"
    p.write_text(header + "\n" + "\n".join(
        r.replace("BED", str(tmp_path /
                             "a.bed.gz")).replace("VCF", str(tmp_path /
                                                             "a.vcf"))
        for r in rows) + "\n")
    return str(p)


def test_samplesheet_ok(tmp_path):
    s = runtime.parse_samplesheet(
        _sheet(tmp_path, ["A,m,BED,VCF", "B,F,BED,VCF"]))
    assert [(x.sample_id, x.sex) for x in s] == [("A", "M"), ("B", "F")]


def test_samplesheet_type_column_keeps_the_control_rows(tmp_path):
    # v2.2: only CTRL rows (any case; 'control' too) build the atlas
    s = runtime.parse_samplesheet(
        _sheet(tmp_path, [
            "A,CTRL,M,BED,VCF", "B,case,F,BED,VCF", "C,Control,F,BED,VCF",
            "D,ctrl,M,BED,VCF", "E,CaSe,M,BED,VCF"
        ],
               header="sample_id,type,sex,meth_bed,vcf"))
    assert [x.sample_id for x in s] == ["A", "C", "D"]


def test_samplesheet_unknown_type_values_stop_the_build_as_in_v2_2(tmp_path):
    header = "sample_id,type,sex,meth_bed,vcf"
    # an empty type, a typo and a repeated typo; the fix is named, and each
    # distinct value is listed once with its first row
    path = _sheet(tmp_path, [
        "A,CTRL,M,BED,VCF", "B,,F,BED,VCF", "C,contol,F,BED,VCF",
        "D,contol,M,BED,VCF", "E,CASE,M,BED,VCF"
    ],
                  header=header)
    with pytest.raises(runtime.RunError) as e:
        runtime.parse_samplesheet(path)
    msg = str(e.value)
    assert "2 value(s) that are neither a control nor CASE" in msg
    assert "empty (row 3)" in msg and "'contol' (row 4)" in msg
    assert msg.count("contol") == 1
    assert "Write CTRL (or control) for a control and CASE for a row" in msg
    # only the first five distinct values are shown
    rows = [f"S{i},T{i},M,BED,VCF" for i in range(8)]
    with pytest.raises(runtime.RunError, match=r"8 value\(s\).* and 3 more"):
        runtime.parse_samplesheet(_sheet(tmp_path, rows, header=header))
    # the types are checked before the files of any row
    path = _sheet(tmp_path, ["A,CTRL,M,/nope.bed.gz,VCF", "B,x,M,BED,VCF"],
                  header=header)
    with pytest.raises(runtime.RunError, match="neither a control nor CASE"):
        runtime.parse_samplesheet(path)


def test_samplesheet_optional_sex_column(tmp_path):
    s = runtime.parse_samplesheet(
        _sheet(tmp_path, ["A,BED,VCF"], header="sample_id,meth_bed,vcf"))
    assert s[0].sex == ""


@pytest.mark.parametrize("rows,header,msg", [
    ([], "sample_id,sex,meth_bed", "missing required column(s): vcf"),
    (["A,M,BED,VCF", "A,M,BED,VCF"], None, "row 3: duplicate sample_id 'A'"),
    ([",M,BED,VCF"], None, "row 2: sample_id is empty"),
    (["A,X,BED,VCF"], None, "sex must be M or F, got 'X'"),
    (["A,M,/nope.bed.gz,VCF"], None, "meth_bed file not found"),
    (["A,M,BED,"], None, "vcf is empty"),
    (["A/b,M,BED,VCF"], None, "may contain only"),
    ([], None, "has a header but no samples"),
])
def test_samplesheet_errors(tmp_path, rows, header, msg):
    path = _sheet(tmp_path, rows, header) if header else _sheet(tmp_path, rows)
    with pytest.raises(runtime.RunError, match=re.escape(msg)):
        runtime.parse_samplesheet(path)


def test_samplesheet_empty_and_missing(tmp_path):
    p = tmp_path / "e.csv"
    p.write_text("")
    with pytest.raises(runtime.RunError, match="is empty"):
        runtime.parse_samplesheet(str(p))
    with pytest.raises(runtime.RunError, match="not found"):
        runtime.parse_samplesheet(str(tmp_path / "none.csv"))


def test_require_sex():
    s = [runtime.Sample("A", "", "b", "v")]
    with pytest.raises(runtime.RunError, match="sex .* is required"):
        runtime.require_sex(s, "why")


def _res(**kw):
    cfg = {o.key: o.default for o in settings.OPTIONS}
    cfg.update(kw)
    return cfg


def test_threads_default_is_half_of_the_available_threads_at_least_two():
    no_cap = math.inf  # the memory of the host does not matter here
    assert runtime.resolve_threads(None, 96, no_cap) == (48, 384.0, None)
    assert runtime.resolve_threads(None, 7, no_cap)[:2] == (3, 24.0)
    assert runtime.resolve_threads(None, 2, no_cap)[:2] == (2, 16.0)
    assert runtime.resolve_threads(None, 1, no_cap)[:2] == (2, 16.0)


def test_threads_memory_is_eight_gb_per_thread_and_one_is_raised_to_two():
    no_cap = math.inf
    assert runtime.resolve_threads(30, 96, no_cap) == (30, 240.0, None)
    n, mem, warn = runtime.resolve_threads(1, 96, no_cap)
    assert (n, mem) == (2, 16.0) and warn and "raised to 2" in warn
    with pytest.raises(runtime.RunError, match="at least 2"):
        runtime.resolve_threads(0, 96, no_cap)
    with pytest.raises(runtime.RunError):
        runtime.resolve_threads(-3, 96, no_cap)


def test_available_threads_follows_affinity_and_cgroup_quota(
        tmp_path, monkeypatch):
    monkeypatch.setattr(runtime.os, "sched_getaffinity",
                        lambda _: set(range(40)))
    monkeypatch.setattr(runtime, "cgroup_cpu_limit", lambda *a, **k: None)
    assert runtime.available_threads() == 40
    monkeypatch.setattr(runtime, "cgroup_cpu_limit", lambda *a, **k: 12.5)
    assert runtime.available_threads() == 12
    monkeypatch.setattr(runtime, "cgroup_cpu_limit", lambda *a, **k: 100.0)
    assert runtime.available_threads() == 40


def test_cgroup_quota_v2_v1_and_none(tmp_path):
    v2 = tmp_path / "v2"
    v2.mkdir()
    (v2 / "cpu.max").write_text("max 100000\n")
    assert runtime.cgroup_cpu_limit(str(v2), str(tmp_path / "none")) is None
    (v2 / "cpu.max").write_text("250000 100000\n")
    assert runtime.cgroup_cpu_limit(str(v2), str(tmp_path / "none")) == 2.5
    v1 = tmp_path / "v1" / "cpu"
    v1.mkdir(parents=True)
    (v1 / "cpu.cfs_quota_us").write_text("-1\n")
    (v1 / "cpu.cfs_period_us").write_text("100000\n")
    assert runtime.cgroup_cpu_limit(str(tmp_path / "v1"),
                                    str(tmp_path / "none")) is None
    (v1 / "cpu.cfs_quota_us").write_text("800000\n")
    assert runtime.cgroup_cpu_limit(str(tmp_path / "v1"),
                                    str(tmp_path / "none")) == 8.0


def test_n_workers_threads_and_memory_limits():
    assert runtime.n_workers(12, 100, [4] * 100) == 12  # threads limit
    assert runtime.n_workers(100, 16, [4] * 100) == 4  # memory limit: 4 x 4 GB
    assert runtime.n_workers(12, 1000, [4] * 3) == 3  # tasks
    assert runtime.n_workers(12, 16, [20, 1, 1, 1]) == 1
    assert runtime.n_workers(12, 30, [20, 5, 5, 5, 5]) == 3  # 20 + 5 + 5 = 30
    assert runtime.n_workers(12, 30, [10, 10, 10, 10], max_workers=2) == 2


def test_a_task_larger_than_the_budget_is_an_error():
    ctx = runtime.Ctx("run", {}, Path("o"), Path("o/work"), 2, 16.0, False)
    with pytest.raises(
            runtime.RunError,
            match=r"estimated peak memory 20.0 GB exceeds the budget of 16 GB "
            r".*--threads 3"):
        runtime.thread_plan(ctx, [20.0, 1.0], [1.0, 1.0], label="atlas")


def test_no_duplicate_resource_options_remain():
    assert [o.key for o in settings.OPTIONS if o.section == "resources"
           ] == ["threads", "keep_tmp"]
    for gone in ("cpus", "max_cpus", "max_memory", "gb_per_task",
                 "core_fraction", "memory_fraction"):
        assert gone not in settings.BY_KEY
    for command in ("build", "run"):
        flags = {
            a.option_strings[0] for a in cli.build_parser()._subparsers.
            _group_actions[0].choices[command]._actions if a.option_strings
        }
        assert "--threads" in flags and not flags & {
            "--cpus", "--max-cpus", "--max-memory"
        }


def _ctrl(n_m, n_f):
    return [runtime.Sample(f"M{i}", "M", "b", "v") for i in range(n_m)
           ] + [runtime.Sample(f"F{i}", "F", "b", "v") for i in range(n_f)]


def test_build_scopes_min_controls():
    cfg = _res()
    sc, w = build.plan_scopes(cfg, _ctrl(4, 4))
    assert list(sc) == ["autosomal", "male_X", "female_X", "Y"] and not w
    sc, w = build.plan_scopes(cfg, _ctrl(2, 5))
    assert list(sc) == ["autosomal", "female_X"
                       ] and len(w) == 2 and sc["autosomal"].__len__() == 7
    sc, _ = build.plan_scopes(_res(include_sex_chromosomes=False), _ctrl(4, 4))
    assert list(sc) == ["autosomal"]


def test_expected_scopes():
    f, m = runtime.Sample("a", "F", "b",
                          "v"), runtime.Sample("b", "M", "b", "v")
    assert run.expected_scopes(f, True) == ["autosomal", "female_X"]
    assert run.expected_scopes(m, True) == ["autosomal", "male_X", "Y"]
    assert run.expected_scopes(m, False) == ["autosomal"]


def test_default_sample_id_is_the_file_name_up_to_the_first_dot():
    assert run.default_sample_id(
        "/x/y/SAMPLE01-run7.combined.bed.gz") == "SAMPLE01-run7"
    assert run.default_sample_id("S1.bed") == "S1" and run.default_sample_id(
        "S1") == "S1"


def _touch(*paths):
    for p in paths:
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text("x")


def _flat(d, sex=True, sigma=150):
    names = ["ref_atlas_%d.bed.gz" % sigma, "cpg_signal.bed", "icr_region.bed"]
    if sex:
        for s in ("male_X", "female_X", "Y"):
            names += [f"ref_atlas_{sigma}_{s}.bed.gz", f"cpg_signal_{s}.bed"]
    _touch(*[d / n for n in names])


def test_reference_flat_dir(tmp_path):
    _flat(tmp_path)
    r = run.resolve_reference(str(tmp_path), 150, True)
    assert r.sex_branches and set(r.atlas) == set(
        runtime.SCOPES) and r.icr_panel and not r.warnings


def test_reference_v1_layout_tsv(tmp_path):
    for s in ("", "_male_X", "_female_X", "_Y"):
        _touch(tmp_path / "atlas" / f"ref_atlas_150{s}.tsv.gz",
               tmp_path / "cpg_islands" / f"cpg_signal{s}.bed")
    _touch(tmp_path / "icr_regions" / "icr_region.bed")
    r = run.resolve_reference(str(tmp_path), 150, True)
    assert r.sex_branches and r.atlas["Y"].endswith(
        "ref_atlas_150_Y.tsv.gz") and r.icr_panel.endswith("icr_region.bed")


def test_reference_prefers_bed_gz_over_tsv_gz(tmp_path):
    _touch(tmp_path / "atlas" / "ref_atlas_150.bed.gz",
           tmp_path / "atlas" / "ref_atlas_150.tsv.gz",
           tmp_path / "cpg_signal.bed")
    assert run.resolve_reference(str(tmp_path), 150,
                                 False).atlas["autosomal"].endswith(".bed.gz")


def test_reference_fallback_when_sex_set_incomplete(tmp_path):
    _flat(tmp_path, sex=False)
    _touch(tmp_path / "ref_atlas_150_Y.bed.gz")
    r = run.resolve_reference(str(tmp_path), 150, True)
    assert not r.sex_branches and list(r.atlas) == ["autosomal"] and any(
        "Falling back" in w for w in r.warnings)
    assert not run.resolve_reference(str(tmp_path), 150, False).warnings


def test_reference_without_icr_panel_warns(tmp_path):
    _touch(tmp_path / "ref_atlas_150.bed.gz", tmp_path / "cpg_signal.bed")
    r = run.resolve_reference(str(tmp_path), 150, False)
    assert r.icr_panel is None and any("ICR" in w for w in r.warnings)


def test_reference_errors(tmp_path):
    with pytest.raises(runtime.RunError, match="not a directory"):
        run.resolve_reference(str(tmp_path / "x"), 150, True)
    with pytest.raises(runtime.RunError, match="no autosomal reference atlas"):
        run.resolve_reference(str(tmp_path), 150, True)
    _touch(tmp_path / "ref_atlas_150.bed.gz")
    with pytest.raises(runtime.RunError, match="CpG island panel"):
        run.resolve_reference(str(tmp_path), 150, True)


def test_reference_json_sigma(tmp_path):
    _flat(tmp_path, sigma=300)
    (tmp_path / "reference.json").write_text('{"smooth_sigma": 300}')
    assert run.resolve_reference(
        str(tmp_path), 150,
        True).atlas["autosomal"].endswith("ref_atlas_300.bed.gz")


def test_guarded_skip_rerun(tmp_path):
    out = tmp_path / "o.txt"
    calls = []

    def fn():
        calls.append(1)
        out.write_text("x")

    st = str(tmp_path)
    assert runtime.guarded(st, "s", {"a": 1}, [str(out)], fn)[1] == "done"
    assert runtime.guarded(st, "s", {"a": 1}, [str(out)], fn)[1] == "skipped"
    assert runtime.guarded(st, "s", {"a": 2}, [str(out)],
                           fn)[1] == "done"  # parameter changed
    out.unlink()
    assert runtime.guarded(st, "s", {"a": 2}, [str(out)],
                           fn)[1] == "done"  # output missing
    assert len(calls) == 3


def test_guarded_failure_leaves_no_stamp(tmp_path):

    def boom():
        raise ValueError("x")

    with pytest.raises(ValueError):
        runtime.guarded(str(tmp_path), "s", {}, [str(tmp_path / "o")], boom)
    assert not (tmp_path / "s.json").exists()


def test_commit_moves_index_first(tmp_path):
    final = str(tmp_path / "f.bed.gz")
    tmp = runtime.tmp_name(final)
    Path(tmp).write_text("d")
    Path(tmp + ".tbi").write_text("i")
    runtime.commit(tmp, final)
    assert Path(final).read_text() == "d" and Path(
        final + ".tbi").read_text() == "i" and not Path(tmp).exists()


def test_temp_output_removes_what_was_not_moved(tmp_path):
    final = str(tmp_path / "f.bed.gz")
    extra = str(tmp_path / "x.bed.gz")
    with runtime.temp_output(final, extra) as tmp:
        Path(tmp).write_text("d")
        runtime.commit(tmp, final)
    assert Path(final).read_text() == "d" and not Path(tmp).exists()
    with pytest.raises(ValueError):
        with runtime.temp_output(final, extra) as tmp:
            for f in (tmp, tmp + ".tbi", extra, extra + ".tbi"):
                Path(f).write_text("partial")
            raise ValueError("boom")
    assert Path(final).read_text() == "d"  # the earlier output stays
    assert not any(p.name != "f.bed.gz" for p in tmp_path.iterdir())


def test_write_atomic_failure_leaves_no_temporary_file_and_the_old_file(
        tmp_path):
    final = tmp_path / "state.json"
    runtime.write_atomic(final, "old")
    assert final.read_text() == "old"
    with pytest.raises(UnicodeEncodeError):
        runtime.write_atomic(final, "\ud800")  # a lone surrogate fails to write
    assert final.read_text() == "old"
    assert [p.name for p in tmp_path.iterdir()] == ["state.json"]
    runtime.write_atomic(final, "new")
    assert final.read_text() == "new"
    assert [p.name for p in tmp_path.iterdir()] == ["state.json"]


def test_forced_format_checks_columns(tmp_path):
    from dmrdetector.fileio import methbed
    f = tmp_path / "a.bed"
    f.write_text("chr1\t1\t2\t10\t5\n")
    with pytest.raises(ValueError, match="needs at least 12"):
        methbed.read_methylation_bed(str(f), "modkit")
    df, fmt = methbed.read_methylation_bed(str(f), "5col")
    assert fmt == "5col" and len(df) == 1


class Env:
    """Inputs on disk and a stub for run._steps.

    The stub records what it was given and writes the step outputs.
    """

    def __init__(self, tmp_path, monkeypatch):
        self.tmp = tmp_path
        self.cwd = tmp_path / "cwd"
        self.cwd.mkdir()
        monkeypatch.chdir(self.cwd)
        self.ref = tmp_path / "ref"
        _flat(self.ref)
        for n in ("S1.combined.bed.gz", "S1.vcf", "g.gtf.gz", "t.bed",
                  "h1.bed.gz", "h2.bed.gz"):
            (tmp_path / n).write_text("x")
        self.seen, self.fail, self.target_fail = {}, None, None
        from dmrdetector.steps import target
        monkeypatch.setattr(run, "_steps", self._steps)
        monkeypatch.setattr(target, "target_report", self._target_report)

    def _steps(self, ctx, cfg, s, ref, par_bed, gene_list, gtf):
        self.seen = dict(workdir=ctx.workdir,
                         sample=s,
                         sex_branches=ref.sex_branches,
                         scopes=list(ref.atlas),
                         gtf=gtf,
                         target_bed=None,
                         cfg=cfg,
                         threads=ctx.threads)
        (ctx.workdir / "annotated" /
         f"{s.sample_id}.dmr.tsv").write_text("CHROM\tSTART\n")
        (ctx.workdir / "scratch" / "intermediate.txt").write_text("x")
        runtime.log("stub step done")
        ctx.timings["stub"] = {"status": "done", "seconds": 0.1}
        if self.fail:
            raise self.fail

    def _target_report(self, samples, atlas, target_bed, out, **kwargs):
        self.seen["target_bed"] = target_bed
        self.seen["target_kwargs"] = kwargs
        if self.target_fail:
            raise self.target_fail
        Path(out).write_text("REGION_NAME\n")

    def args(self, *extra, sex="F", gtf=True, reference=True, sample_id=None):
        a = [
            "run", "--meth_bed_combined",
            str(self.tmp / "S1.combined.bed.gz"), "--vcf",
            str(self.tmp / "S1.vcf"), "--threads", "2"
        ]
        if sex:
            a += ["--sex", sex]
        if gtf:
            a += ["--gtf", str(self.tmp / "g.gtf.gz")]
        if reference:
            a += ["--reference", str(self.ref)]
        if sample_id:
            a += ["--sample-id", sample_id]
        return [*a, *extra]

    def go(self, *extra, **kw):
        return cli.main(self.args(*extra, **kw))


@pytest.fixture
def env(tmp_path, monkeypatch):
    return Env(tmp_path, monkeypatch)


def _names(d):
    return sorted(p.name for p in Path(d).iterdir())


def test_run_temporary_folder_is_made_in_tmp_dir_and_removed(env, tmp_path):
    tmpd = tmp_path / "scratch"
    assert env.go("--tmp-dir", str(tmpd), "--outdir", str(tmp_path / "o")) == 0
    assert env.seen["workdir"].parent == tmpd.resolve(
    ) and env.seen["workdir"].name.startswith("dmrdetector-")
    assert not env.seen["workdir"].exists() and _names(tmpd) == [
    ]  # the folder of the run is gone, --tmp-dir stays


def test_run_never_removes_what_the_user_put_in_tmp_dir(env, tmp_path):
    tmpd = tmp_path / "scratch"
    tmpd.mkdir()
    (tmpd / "mine.txt").write_text("keep")
    (tmpd / "dmrdetector-mine").mkdir()
    env.go("--tmp-dir", str(tmpd), "--outdir", str(tmp_path / "o"))
    assert _names(tmpd) == ["dmrdetector-mine", "mine.txt"]


def test_run_temporary_folder_removed_after_failure_and_log_says_failed(
        env, tmp_path, capsys):
    out = tmp_path / "o"
    env.fail = runtime.RunError("step exploded")
    assert env.go("--outdir", str(out), "--tmp-dir",
                  str(tmp_path / "scratch")) == 1
    assert "ERROR: step exploded" in capsys.readouterr().err
    assert not env.seen["workdir"].exists()
    assert _names(out) == ["S1.log"]  # no partial DMR table
    log = (out / "S1.log").read_text()
    assert 'status:   FAILED' in log
    assert 'step exploded' in log
    assert 'stub step done' in log
    env.fail = RuntimeError("unexpected")
    with pytest.raises(
            RuntimeError
    ):  # not a user error: the traceback is kept, the cleanup still happens
        env.go("--outdir", str(out), "--tmp-dir", str(tmp_path / "scratch"))
    assert not env.seen["workdir"].exists() and _names(
        tmp_path / "scratch") == []


def test_run_restores_the_temporary_directory_of_the_process(env, tmp_path):
    before = (os.environ.get("TMPDIR"), tempfile.tempdir, tempfile.gettempdir())
    env.go("--outdir", str(tmp_path / "o"))
    assert (os.environ.get("TMPDIR"), tempfile.tempdir
           ) == before[:2] and tempfile.gettempdir() == before[2]


def test_run_without_reference_or_inputs(env, tmp_path, capsys):
    assert env.go("--outdir", str(tmp_path / "o"), reference=False) == 1
    assert "no reference: pass --reference DIR" in capsys.readouterr().err
    assert cli.main([
        "run", "--vcf",
        str(tmp_path / "S1.vcf"), "--gtf",
        str(tmp_path / "g.gtf.gz")
    ]) == 1
    assert "--meth_bed_combined FILE is required" in capsys.readouterr().err
    assert cli.main(
        ["run", "--meth_bed_combined",
         str(tmp_path / "S1.combined.bed.gz")]) == 1
    assert "--vcf FILE is required" in capsys.readouterr().err
    assert env.go(
        sex="Q") == 1 and "--sex must be M or F" in capsys.readouterr().err
    assert not (tmp_path / "o").exists() and not (env.cwd /
                                                  "dmrdetector").exists()


def test_run_sample_id_default_and_rules(env, tmp_path, capsys):
    env.go("--outdir", str(tmp_path / "o1"))
    assert env.seen["sample"].sample_id == "S1"
    env.go("--outdir", str(tmp_path / "o2"), sample_id="Case_1.a-b")
    assert _names(tmp_path / "o2") == ["Case_1.a-b.dmr.tsv", "Case_1.a-b.log"]
    for bad in ("a/b", "a b", "ä"):
        assert env.go("--outdir", str(tmp_path / "o3"), sample_id=bad) == 1
        assert "may contain only letters, digits" in capsys.readouterr().err
    (tmp_path / "A B.bed.gz").write_text("x")
    assert cli.main([
        "run", "--meth_bed_combined",
        str(tmp_path / "A B.bed.gz"), "--vcf",
        str(tmp_path / "S1.vcf"), "--gtf",
        str(tmp_path / "g.gtf.gz"), "--reference",
        str(env.ref)
    ]) == 1
    assert "default sample_id 'A B' may contain only" in capsys.readouterr().err


def test_run_without_sex_gives_one_notice_and_the_autosomal_scope(
        env, tmp_path, capsys):
    out = tmp_path / "o"
    assert env.go("--outdir", str(out), sex=None) == 0
    err = capsys.readouterr().err
    assert err.count("--sex not given") == 1 and err.count("NOTICE") == 1
    assert env.seen["sex_branches"] is False and env.seen["scopes"] == [
        "autosomal"
    ]
    assert "NOTICE: --sex not given" in (
        out /
        "S1.log").read_text().split("[warnings]")[1].split("[messages]")[0]
    assert env.go("--outdir", str(out), sex="M") == 0
    assert "--sex not given" not in capsys.readouterr().err
    assert env.seen["sex_branches"] is True and env.seen["scopes"] == list(
        runtime.SCOPES)
    assert env.go("--outdir",
                  str(out),
                  "--set",
                  "include_sex_chromosomes=false",
                  sex=None) == 0
    assert "--sex not given" not in capsys.readouterr(
    ).err  # nothing to announce when chrX and chrY are off


def test_run_target_region_is_a_path_and_haplotype_files_come_from_the_options(
        env, tmp_path, capsys):
    t = str(tmp_path / "t.bed")
    assert env.go("--outdir", str(tmp_path / "o"),
                  "--target-region", t, "--meth_bed_hap1",
                  str(tmp_path / "h1.bed.gz"), "--meth-bed-hap2",
                  str(tmp_path / "h2.bed.gz")) == 0
    s = env.seen["sample"]
    assert env.seen["target_bed"] == t and s.meth_bed_hap1 == str(
        tmp_path / "h1.bed.gz") and s.meth_bed_hap2 == str(
            tmp_path / "h2.bed.gz")
    assert env.go("--outdir", str(tmp_path / "o"), "--target-region",
                  str(tmp_path / "missing.bed")) == 1
    assert "--target-region BED not found" in capsys.readouterr().err
    assert env.go("--outdir", str(tmp_path / "o"), "--meth_bed_hap1",
                  str(tmp_path / "none.bed.gz")) == 1
    assert "--meth_bed_hap1 file not found" in capsys.readouterr().err
    env.seen = {}
    assert env.go("--outdir", str(tmp_path / "o2"), "--meth_bed_hap1",
                  str(tmp_path / "h1.bed.gz")) == 0
    assert "used only with --target-region" in capsys.readouterr(
    ).err and env.seen["target_bed"] is None
    # the haplotype files are not found by file-name convention
    (tmp_path / "S1.hap1.bed.gz").write_text("x")
    env.go("--outdir", str(tmp_path / "o3"), "--target-region", t)
    assert env.seen["sample"].meth_bed_hap1 is None


def test_run_reports_option_errors_before_creating_anything(
        env, tmp_path, capsys):
    assert env.go("--outdir", str(tmp_path / "o"), "--set",
                  "dmr_fdr_threshold=7") == 1
    assert "must be between 0 and 1" in capsys.readouterr().err and not (
        tmp_path / "o").exists()
    assert env.go("--outdir", str(tmp_path / "o"), "--config",
                  str(tmp_path / "none.toml")) == 1
    assert "config file not found" in capsys.readouterr().err


def test_build_removes_only_its_own_work_folder(tmp_path, monkeypatch):
    """Checks that `build` removes its own folder after success only.

    The folder name below --tmp-dir follows the outdir; --tmp-dir itself stays,
    and a failure keeps the folder for the rerun.
    """
    sheet = _sheet(tmp_path, ["A,M,BED,VCF", "B,F,BED,VCF", "C,F,BED,VCF"])
    tmpd, out = tmp_path / "scratch", tmp_path / "ref"
    tmpd.mkdir()
    (tmpd / "mine.txt").write_text("keep")
    calls = []

    def fake_steps(ctx, *a, **k):
        calls.append(ctx.workdir)
        (ctx.workdir / "stamps" / "x.json").write_text("{}")
        if len(calls) == 1:
            raise runtime.RunError("boom")

    monkeypatch.setattr(build, "_steps", fake_steps)
    args = [
        "build", "--input", sheet, "--outdir",
        str(out), "--tmp-dir",
        str(tmpd), "--threads", "2", "--set", "include_sex_chromosomes=false"
    ]
    assert cli.main(args) == 1
    assert calls[0] == runtime.build_workdir(str(tmpd), out.resolve()) and (
        calls[0] / "stamps" / "x.json").is_file()  # kept: a rerun resumes
    assert "status:   FAILED" in (out / "build.log").read_text()
    assert cli.main(args) == 0
    assert _names(tmpd) == [
        "mine.txt"
    ] and "status:   completed" in (out / "build.log").read_text()
    assert cli.main([*args, "--keep-tmp"
                    ]) == 0 and _names(tmpd) == [calls[0].name, "mine.txt"]
    out2 = tmp_path / "ref2"
    assert cli.main([
        "build", "--input", sheet, "--outdir",
        str(out2), "--threads", "2", "--set", "include_sex_chromosomes=false"
    ]) == 0
    assert not (out2 / "work").exists(
    )  # default work folder <outdir>/work is removed after success


def test_package_versions_line_names_python_and_the_numerical_packages():
    line = runtime.package_versions()
    for name in ("python", "numpy", "scipy", "pandas", "pysam", "numba"):
        assert re.search(rf"(^|, ){name} \d", line), line


def test_run_target_region_builtin_selects_the_packaged_bed(
        env, tmp_path, capsys):
    from dmrdetector import data
    packaged = data.asset("target_region.hg38.bed")
    assert os.path.isfile(packaged)
    assert env.go("--outdir", str(tmp_path / "o"), "--target-region",
                  "builtin") == 0
    assert env.seen["target_bed"] == packaged
    assert (tmp_path / "o" / "S1.target_region.tsv").is_file()
    # the settings file and --set give the same
    f = _write(tmp_path, '[target]\ntarget_region = "builtin"\n')
    env.seen = {}
    assert env.go("--outdir", str(tmp_path / "o2"), "--config", f) == 0
    assert env.seen["target_bed"] == packaged
    env.seen = {}
    assert env.go("--outdir", str(tmp_path / "o3"), "--set",
                  "target_region=builtin") == 0
    assert env.seen["target_bed"] == packaged
    # a path keeps working, relative to the working directory
    (tmp_path / "mine.bed").write_text("x")
    env.seen = {}
    assert env.go("--outdir", str(tmp_path / "o4"), "--target-region",
                  str(tmp_path / "mine.bed")) == 0
    assert env.seen["target_bed"] == str(tmp_path / "mine.bed")


def test_builtin_is_documented_in_help_and_settings(capsys):
    out = _help(capsys, "run")
    assert "builtin" in out
    assert "builtin" in settings.default_toml()


def test_set_suggests_the_nearest_key_for_a_short_name():
    with pytest.raises(settings.ConfigError,
                       match=re.escape("did you mean 'dmr_fdr_threshold'")):
        settings.parse_set(["dmr_fdr=0.05"], "run")
    with pytest.raises(settings.ConfigError) as e:
        settings.parse_set(["zzzzzz=1"], "run")
    assert "did you mean" not in str(e.value)


def test_wrong_command_lists_only_the_visible_commands(capsys):
    with pytest.raises(SystemExit):
        cli.build_parser().parse_args(["zzz"])
    err = capsys.readouterr().err
    assert "invalid choice: 'zzz'" in err
    for hidden in ("preprocess", "atlas", "cpg-islands", "call-dmr", "target"):
        assert f"'{hidden}'" not in err
    for shown in ("build", "run", "config", "get-genes"):
        assert f"'{shown}'" in err
    # a hidden command is still accepted
    a = cli.build_parser().parse_args([
        "target", "--samplesheet", "s.csv", "--atlas", "a", "--target-bed", "t"
    ])
    assert a.command == "target"
    assert "(v2.2)" not in _help(capsys, "atlas") + _help(capsys, "preprocess")


def test_samplesheet_type_column_gives_a_notice_with_the_count(
        tmp_path, capsys):
    header = "sample_id,type,sex,meth_bed,vcf"
    runtime.parse_samplesheet(
        _sheet(tmp_path, ["A,CTRL,M,BED,VCF", "B,CASE,F,BED,VCF"],
               header=header))
    err = capsys.readouterr().err
    assert "NOTICE" in err and "1 of 2 rows are dropped" in err
    assert "1 kept as controls" in err
    # nothing to drop, or no type column: silent
    runtime.parse_samplesheet(
        _sheet(tmp_path, ["A,CTRL,M,BED,VCF"], header=header))
    runtime.parse_samplesheet(_sheet(tmp_path, ["A,M,BED,VCF"]))
    assert "'type' column" not in capsys.readouterr().err


def test_samplesheet_type_column_without_a_control_is_an_error(tmp_path):
    with pytest.raises(runtime.RunError, match="no row with type CTRL"):
        runtime.parse_samplesheet(
            _sheet(tmp_path, ["A,CASE,M,BED,VCF"],
                   header="sample_id,type,sex,meth_bed,vcf"))


def test_samplesheet_dropped_rows_are_not_validated(tmp_path):
    # a case row may name files that do not exist and repeat an id; a control
    # row is still checked
    path = _sheet(tmp_path, [
        "A,CTRL,M,BED,VCF", "A,CASE,F,/nope.bed.gz,/nope.vcf",
        "B,CTRL,X,BED,VCF"
    ],
                  header="sample_id,type,sex,meth_bed,vcf")
    with pytest.raises(runtime.RunError, match="row 4 .*sex must be M or F"):
        runtime.parse_samplesheet(path)
    s = runtime.parse_samplesheet(
        _sheet(tmp_path, ["A,CTRL,M,BED,VCF", "A,CASE,F,/nope.bed.gz,VCF"],
               header="sample_id,type,sex,meth_bed,vcf"))
    assert [x.sample_id for x in s] == ["A"]


def test_samplesheet_relative_path_hint(tmp_path, monkeypatch):
    sheet_dir = tmp_path / "sheets"
    sheet_dir.mkdir()
    (sheet_dir / "a.bed.gz").write_text("x")
    (sheet_dir / "a.vcf").write_text("x")
    sheet = sheet_dir / "s.csv"
    sheet.write_text("sample_id,sex,meth_bed,vcf\nA,M,a.bed.gz,a.vcf\n")
    monkeypatch.chdir(tmp_path)
    with pytest.raises(runtime.RunError) as e:
        runtime.parse_samplesheet(str(sheet))
    assert "file not found: a.bed.gz" in str(e.value)
    assert "counts from the working directory" in str(e.value)
    # no hint when the folders agree
    monkeypatch.chdir(sheet_dir)
    assert runtime.parse_samplesheet("s.csv")[0].sample_id == "A"
    (sheet_dir / "a.vcf").unlink()
    with pytest.raises(runtime.RunError) as e:
        runtime.parse_samplesheet("s.csv")
    assert "working directory" not in str(e.value)


def test_input_that_is_not_a_regular_file_is_named(tmp_path):
    with pytest.raises(runtime.RunError, match="not a regular file"):
        runtime.parse_samplesheet(os.devnull)
    with pytest.raises(runtime.RunError, match="not a regular file"):
        runtime.check_files([("--vcf", os.devnull)])


def test_run_logs_the_settings_file(env, tmp_path, capsys, monkeypatch):
    monkeypatch.chdir(tmp_path)
    assert env.go("--outdir", str(tmp_path / "o")) == 0
    err = capsys.readouterr().err
    assert "settings: built-in defaults" in err and "NOTICE" not in err
    log = (tmp_path / "o" / "S1.log").read_text()
    assert "settings file: (none)" in log
    f = _write(tmp_path, "[dmr]\ndmr_fdr_threshold = 0.02\n", "settings.toml")
    assert env.go("--outdir", str(tmp_path / "o2")) == 0
    err = capsys.readouterr().err
    assert f"settings: {f}" in err
    assert "differs from the packaged defaults in 1 setting(s): " \
        "dmr_fdr_threshold" in err
    # a file that holds the defaults only gives no NOTICE
    Path(f).write_text(settings.default_toml())
    assert env.go("--outdir", str(tmp_path / "o3")) == 0
    err = capsys.readouterr().err
    assert f"settings: {f}" in err and "differs from the packaged" not in err


def test_build_warns_when_it_overwrites_a_reference_but_not_on_resume(
        tmp_path, monkeypatch, capsys):
    sheet = _sheet(tmp_path, ["A,M,BED,VCF", "B,F,BED,VCF", "C,F,BED,VCF"])
    out = tmp_path / "ref"
    calls = []

    def fake_steps(ctx, *a, **k):
        calls.append(1)
        (ctx.workdir / "stamps" / "x.json").write_text("{}")
        if len(calls) == 2:
            raise runtime.RunError("boom")

    monkeypatch.setattr(build, "_steps", fake_steps)
    args = [
        "build", "--input", sheet, "--outdir",
        str(out), "--threads", "2", "--set", "include_sex_chromosomes=false"
    ]
    assert cli.main(args) == 0
    assert "already holds a reference" not in capsys.readouterr().err
    (out / "reference.json").write_text("{}")  # what a finished build leaves
    assert cli.main(args) == 1  # overwrite, then the failure keeps stamps
    assert "already holds a reference" in capsys.readouterr().err
    assert cli.main(args) == 0  # resume of the interrupted build
    assert "already holds a reference" not in capsys.readouterr().err


def test_empty_methylation_bed_warns(tmp_path, capsys):
    from dmrdetector.steps import preprocess
    bed = tmp_path / "empty.bed"
    bed.write_text("")
    vcf = tmp_path / "s.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\t"
                   "FILTER\tINFO\tFORMAT\tS\n")
    (tmp_path / "split").mkdir()
    (tmp_path / "stamps").mkdir()
    t = dict(id="E",
             meth_bed=str(bed),
             vcf=str(vcf),
             min_cov=5,
             fmt="auto",
             out=str(tmp_path / "E.bed.gz"),
             stamps=str(tmp_path / "stamps"),
             split=False,
             split_dir=str(tmp_path / "split"),
             par_bed=None,
             threads=1)
    preprocess._t_preprocess(t)
    err = capsys.readouterr().err
    assert "WARNING: E: no CpG records in" in err


# ---- work folder of `build` (marker, lock, shared --tmp-dir) -----------------


def _build_cli(sheet, out, *extra):
    return [
        "build", "--input", sheet, "--outdir",
        str(out), "--threads", "2", "--set", "include_sex_chromosomes=false",
        *extra
    ]


def _fake_steps(monkeypatch, each=None):
    seen = []

    def steps(ctx, *a, **k):
        seen.append(ctx.workdir)
        (ctx.workdir / "stamps" / "x.json").write_text("{}")
        if each:
            each(ctx, len(seen))

    monkeypatch.setattr(build, "_steps", steps)
    return seen


def test_build_never_touches_a_work_folder_it_did_not_make(
        tmp_path, monkeypatch, capsys):
    sheet = _sheet(tmp_path, ["A,M,BED,VCF", "B,F,BED,VCF", "C,F,BED,VCF"])
    out = tmp_path / "ref"
    work = out / "work"
    (work / "mine").mkdir(parents=True)
    (work / "mine" / "data.txt").write_text("precious")
    (work / "notes.txt").write_text("keep")
    seen = _fake_steps(monkeypatch)
    assert cli.main(_build_cli(sheet, out)) == 1
    err = capsys.readouterr().err
    assert err.startswith("ERROR: work folder ") and str(work) in err
    assert "--tmp-dir" in err and "not made by dmrdetector build" in err
    assert seen == [] and _names(work) == ["mine", "notes.txt"]
    assert (work / "mine" / "data.txt").read_text() == "precious"
    assert not (out / "build.log").exists()


def test_build_adopts_an_empty_work_folder_and_removes_it_after_success(
        tmp_path, monkeypatch):
    sheet = _sheet(tmp_path, ["A,M,BED,VCF", "B,F,BED,VCF", "C,F,BED,VCF"])
    out = tmp_path / "ref"
    (out / "work").mkdir(parents=True)
    _fake_steps(monkeypatch)
    assert cli.main(_build_cli(sheet, out)) == 0
    assert not (out / "work").exists()


def test_build_resumes_its_interrupted_work_folder(tmp_path, monkeypatch):
    sheet = _sheet(tmp_path, ["A,M,BED,VCF", "B,F,BED,VCF", "C,F,BED,VCF"])
    out = tmp_path / "ref"
    found = []

    def each(ctx, n):
        found.append(n)
        if n == 1:
            raise runtime.RunError("boom")

    seen = _fake_steps(monkeypatch, each)
    # the stamp of the first attempt exists when the second starts
    first = (out / "work" / "stamps" / "x.json")
    assert cli.main(_build_cli(sheet, out)) == 1 and first.is_file()
    marker = json.loads((out / "work" / runtime.WORK_MARKER).read_text())
    assert marker["outdir"] == str(
        out.resolve()) and marker["status"] == "running"
    assert marker["version"] and marker["pid"] == os.getpid()
    assert cli.main(_build_cli(sheet, out)) == 0 and len(seen) == 2
    assert not (out / "work").exists()


def test_work_folder_of_another_outdir_is_refused(tmp_path):
    wd = tmp_path / "w"
    claim = runtime.claim_workdir(wd, tmp_path / "ref1")
    runtime.release_workdir(claim, False, True)
    with pytest.raises(runtime.RunError, match="belongs to a build of .*ref1"):
        runtime.claim_workdir(wd, tmp_path / "ref2")


def test_a_second_build_cannot_lock_the_same_work_folder(tmp_path):
    wd, out = tmp_path / "w", tmp_path / "ref"
    claim = runtime.claim_workdir(wd, out)
    assert not claim.resumed
    with pytest.raises(runtime.RunError, match="another dmrdetector build"):
        runtime.claim_workdir(wd, out)
    assert (wd / runtime.WORK_MARKER).is_file()  # the refused build left it
    runtime.release_workdir(claim, False, True)
    again = runtime.claim_workdir(wd, out)  # the lock is gone: a resume
    assert again.resumed
    runtime.release_workdir(again, True, False)
    assert not wd.exists()


def test_a_lock_error_of_the_file_system_is_a_warning(tmp_path, monkeypatch,
                                                      capsys):

    def no_lock(fd, op):
        raise OSError(37, "No locks available")

    monkeypatch.setattr(fcntl, "flock", no_lock)
    claim = runtime.claim_workdir(tmp_path / "w", tmp_path / "ref")
    assert "WARNING: cannot lock the work folder" in capsys.readouterr().err
    runtime.release_workdir(claim, True, False)


def test_the_work_folder_is_removed_only_with_its_own_marker(tmp_path):
    wd, out = tmp_path / "w", tmp_path / "ref"
    claim = runtime.claim_workdir(wd, out)
    (wd / runtime.WORK_MARKER).write_text(json.dumps({"outdir": "/elsewhere"}))
    assert runtime.release_workdir(claim, True, False) is False and wd.is_dir()
    claim = runtime.WorkClaim(wd, out, False)
    (wd / runtime.WORK_MARKER).unlink()
    assert runtime.release_workdir(claim, True, False) is False and wd.is_dir()


def test_builds_sharing_a_tmp_dir_keep_their_own_folders(
        tmp_path, monkeypatch, capsys):
    sheet = _sheet(tmp_path, ["A,M,BED,VCF", "B,F,BED,VCF", "C,F,BED,VCF"])
    tmpd = tmp_path / "scratch"
    out1, out2 = tmp_path / "ref1", tmp_path / "ref2"
    inner = []

    def each(ctx, n):
        if n != 1:
            return
        # while build 1 runs: a build of another outdir, then one of the same
        inner.append(cli.main(_build_cli(sheet, out2, "--tmp-dir", str(tmpd))))
        inner.append(cli.main(_build_cli(sheet, out1, "--tmp-dir", str(tmpd))))
        assert (ctx.workdir / "stamps" / "x.json").is_file()
        assert (ctx.workdir / runtime.WORK_MARKER).is_file()

    seen = _fake_steps(monkeypatch, each)
    assert cli.main(_build_cli(sheet, out1, "--tmp-dir", str(tmpd))) == 0
    assert inner == [0, 1]  # the same outdir is locked by build 1
    assert "another dmrdetector build" in capsys.readouterr().err
    assert seen[0] != seen[1] and seen[0].parent == seen[1].parent
    assert seen[0] == runtime.build_workdir(str(tmpd), out1.resolve())
    assert _names(tmpd) == []  # each build removed its own folder only


def test_the_overwrite_warning_sees_atlas_files_without_reference_json(
        tmp_path):
    assert not build._holds_reference(tmp_path)
    (tmp_path / "other.txt").write_text("x")
    assert not build._holds_reference(tmp_path)
    (tmp_path / "ref_atlas_150.bed.gz").write_text("x")
    assert build._holds_reference(tmp_path)


def test_stale_temporary_files_of_dead_processes_are_removed(tmp_path):
    child = subprocess.run(
        [sys.executable, "-c", "import os; print(os.getpid())"],
        capture_output=True,
        text=True)
    dead = int(child.stdout)  # finished: no process has this id
    old = 1.0e9
    now = old + 7200

    def make(name, mtime=old, folder=tmp_path):
        p = folder / name
        p.write_text("x")
        os.utime(p, (mtime, mtime))

    make(f".tmp{dead}.ref_atlas_150.bed.gz")  # dead and old: removed
    make(f".tmp{dead}.recent.tsv", mtime=now - 60)  # dead but recent: kept
    make(f".tmp{os.getpid()}.alive.tsv")  # a live process: kept
    make("results.tsv")  # not a temporary name: kept
    make(f".tmp{dead}")  # no file name part: kept
    (tmp_path / "sub").mkdir()
    make(f".tmp{dead}.nested.tsv", folder=tmp_path / "sub")  # below: kept
    removed = runtime.remove_stale_tmp(tmp_path, now=now)
    assert removed == [f".tmp{dead}.ref_atlas_150.bed.gz"]
    assert _names(tmp_path) == sorted([
        f".tmp{dead}.recent.tsv", f".tmp{os.getpid()}.alive.tsv", "results.tsv",
        f".tmp{dead}", "sub"
    ])
    assert (tmp_path / "sub" / f".tmp{dead}.nested.tsv").exists()


def test_a_stamp_of_another_version_is_not_resumed(tmp_path, monkeypatch):
    stamps = tmp_path / "stamps"
    stamps.mkdir()
    out = tmp_path / "out.txt"
    calls = []

    def step():
        calls.append(1)
        out.write_text("x")

    def go():
        return runtime.guarded(str(stamps), "s", {"a": 1}, [str(out)], step)[1]

    assert go() == "done" and go() == "skipped"
    monkeypatch.setattr(runtime, "__version__", "0.0.0-other")
    assert go() == "done" and go() == "skipped" and len(calls) == 2
    monkeypatch.setattr(runtime, "STAMP_FORMAT", runtime.STAMP_FORMAT + 1)
    assert go() == "done"


# ---- threads, memory, cgroups ------------------------------------------------


def test_threads_above_the_available_threads_are_lowered_with_a_warning():
    n, mem, warn = runtime.resolve_threads(30, available=8, mem_cap_gb=math.inf)
    assert (n, mem) == (8, 64.0) and "lowered to 8" in warn
    assert runtime.resolve_threads(8, 8, math.inf) == (8, 64.0, None)
    n, mem, warn = runtime.resolve_threads(5, available=1, mem_cap_gb=math.inf)
    assert (n, mem) == (2, 16.0) and "lowered to 2" in warn  # the minimum


def test_memory_budget_is_capped_by_the_memory_the_process_may_use():
    n, mem, warn = runtime.resolve_threads(10, 10, 32.0)
    assert (n, mem) == (10, 32.0)
    assert "lowered from 80 GB to 32 GB" in warn
    n, mem, warn = runtime.resolve_threads(None, 10, 100.0)
    assert (n, mem, warn) == (5, 40.0, None)  # 8 GB x 5 is below the cap


def test_memory_cap_is_the_smaller_of_meminfo_and_cgroup(tmp_path, monkeypatch):
    info = tmp_path / "meminfo"
    info.write_text("MemTotal:  200000000 kB\nMemFree: 1 kB\n"
                    "MemAvailable:   50000000 kB\n")
    assert runtime.meminfo_available_gb(str(info)) == pytest.approx(51.2)
    assert runtime.meminfo_available_gb(str(tmp_path / "none")) is None
    monkeypatch.setattr(runtime, "meminfo_available_gb", lambda *a: 100.0)
    monkeypatch.setattr(runtime, "cgroup_memory_limit_gb", lambda *a: 40.0)
    assert runtime.memory_cap_gb() == 40.0
    monkeypatch.setattr(runtime, "cgroup_memory_limit_gb", lambda *a: None)
    assert runtime.memory_cap_gb() == 100.0
    monkeypatch.setattr(runtime, "meminfo_available_gb", lambda *a: None)
    assert runtime.memory_cap_gb() is None


def test_cgroup_memory_limit_v2_v1_and_unlimited(tmp_path):
    none = str(tmp_path / "none")
    v2 = tmp_path / "v2"
    v2.mkdir()
    (v2 / "memory.max").write_text("max\n")
    assert runtime.cgroup_memory_limit_gb(str(v2), none) is None
    (v2 / "memory.max").write_text("32000000000\n")
    assert runtime.cgroup_memory_limit_gb(str(v2), none) == 32.0
    v1 = tmp_path / "v1" / "memory"
    v1.mkdir(parents=True)
    (v1 / "memory.limit_in_bytes").write_text("9223372036854771712\n")
    assert runtime.cgroup_memory_limit_gb(str(tmp_path / "v1"), none) is None
    (v1 / "memory.limit_in_bytes").write_text("16000000000\n")
    assert runtime.cgroup_memory_limit_gb(str(tmp_path / "v1"), none) == 16.0


def test_cgroup_limits_of_the_parents_apply(tmp_path):
    root = tmp_path / "cg"
    (root / "a" / "b").mkdir(parents=True)
    me = tmp_path / "self"
    me.write_text("0::/a/b\n")
    (root / "cpu.max").write_text("max 100000\n")
    (root / "a" / "cpu.max").write_text("200000 100000\n")
    (root / "a" / "b" / "cpu.max").write_text("max 100000\n")
    (root / "a" / "memory.max").write_text("64000000000\n")
    assert runtime.cgroup_cpu_limit(str(root), str(me)) == 2.0
    assert runtime.cgroup_memory_limit_gb(str(root), str(me)) == 64.0
    (root / "a" / "b" / "cpu.max").write_text("100000 100000\n")  # leaf wins
    assert runtime.cgroup_cpu_limit(str(root), str(me)) == 1.0
    # an unrelated line of the file must not hide the group of the controller
    v1 = tmp_path / "v1"
    (v1 / "cpu,cpuacct" / "g1" / "inner").mkdir(parents=True)
    (v1 / "cpu,cpuacct" / "g1" / "cpu.cfs_quota_us").write_text("300000\n")
    (v1 / "cpu,cpuacct" / "g1" / "cpu.cfs_period_us").write_text("100000\n")
    (v1 / "cpu,cpuacct" / "g1" / "inner" /
     "cpu.cfs_quota_us").write_text("-1\n")
    (v1 / "cpu,cpuacct" / "g1" / "inner" /
     "cpu.cfs_period_us").write_text("100000\n")
    me1 = tmp_path / "self1"
    me1.write_text("12:pids:/user.slice\n4:cpu,cpuacct:/g1/inner\n")
    assert runtime.cgroup_cpu_limit(str(v1), str(me1)) == 3.0


def test_a_budget_below_eight_gb_per_thread_is_reported_as_the_machine():
    with pytest.raises(runtime.RunError) as e:
        runtime.check_fits("x", 20.0, 16.0, threads=4)
    assert "has 16 GB" not in str(e.value) and "available" in str(e.value)
    assert "--threads" not in str(e.value)
    with pytest.raises(runtime.RunError, match="use --threads 3 or more"):
        runtime.check_fits("x", 20.0, 16.0, threads=2)


def _die(task):
    os._exit(9)


def test_a_killed_worker_gives_a_clear_error(tmp_path):
    (tmp_path / "scratch").mkdir()
    ctx = runtime.Ctx("run", {}, tmp_path, tmp_path, 4, 100.0)
    with pytest.raises(runtime.RunError) as e:
        runtime.pool_map(ctx, _die, [{}, {}], 1.0, "kill")
    assert "a worker process was killed, probably out of memory" in str(e.value)
    assert "use fewer --threads" in str(e.value)


# ---- target-region failure, replaced outputs, paths in the log ---------------


def test_a_failed_target_report_keeps_the_dmr_table(env, tmp_path, capsys):
    out = tmp_path / "o"
    out.mkdir()
    (out / "S1.target_region.tsv").write_text("older run\n")
    env.target_fail = ValueError("no autosomal target regions left")
    assert env.go("--outdir", str(out), "--target-region",
                  str(tmp_path / "t.bed")) == 1
    err = capsys.readouterr().err
    assert ("WARNING: target-region report failed: target region: "
            "ValueError") in err
    assert ("S1.dmr.tsv was published and is complete; the run is "
            "incomplete") in err
    assert "S1.target_region.tsv in the folder is from an earlier run" in err
    assert "ERROR: target region" in err
    assert (out / "S1.dmr.tsv").read_text() == "CHROM\tSTART\n"
    assert (out / "S1.target_region.tsv").read_text() == "older run\n"
    log = (out / "S1.log").read_text()
    assert "status:   FAILED (the DMR table was published; the target-region " \
        "report failed)" in log
    assert not env.seen["workdir"].exists()


def test_replaced_outputs_are_announced_and_never_deleted(
        env, tmp_path, capsys):
    out = tmp_path / "o"
    assert env.go("--outdir", str(out), "--target-region",
                  str(tmp_path / "t.bed")) == 0
    assert "will be replaced" not in capsys.readouterr().err
    assert env.go("--outdir", str(out), "--target-region",
                  str(tmp_path / "t.bed")) == 0
    err = capsys.readouterr().err
    for name in ("S1.dmr.tsv", "S1.log", "S1.target_region.tsv"):
        assert f"WARNING: {name} (" in err
    # a run without a target report leaves the old one and says so
    assert env.go("--outdir", str(out)) == 0
    err = capsys.readouterr().err
    assert ("S1.target_region.tsv (the target-region report of an earlier "
            "run) is left in place") in err
    assert (out / "S1.target_region.tsv").is_file()


def test_the_log_names_files_not_folders_unless_full_paths_are_asked(
        env, tmp_path):
    out = tmp_path / "o"
    assert env.go("--outdir", str(out)) == 0
    log = (out / "S1.log").read_text()
    assert str(tmp_path) not in log
    assert "meth_bed_combined = S1.combined.bed.gz (1 bytes)" in log
    assert "vcf = S1.vcf (1 bytes)" in log
    assert "outdir = o" in log and "log_input_paths = 'name'" not in log
    assert "log_input_paths = name  # default" in log
    assert "--reference ref" in log  # the command line is shortened too
    out2 = tmp_path / "o2"
    assert env.go("--outdir", str(out2), "--set", "log_input_paths=full") == 0
    full = (out2 / "S1.log").read_text()
    assert f"meth_bed_combined = {tmp_path / 'S1.combined.bed.gz'}" in full
    assert f"vcf = {tmp_path / 'S1.vcf'}" in full


def test_log_input_paths_accepts_name_and_full_only(tmp_path):
    f = _write(tmp_path, '[paths]\nlog_input_paths = "relative"\n')
    with pytest.raises(settings.ConfigError, match="one of name, full"):
        settings.load_file(f)
    assert settings.BY_KEY["log_input_paths"].default == "name"


def test_build_log_and_reference_json_show_file_names(tmp_path, monkeypatch):
    sheet = _sheet(tmp_path, ["A,M,BED,VCF", "B,F,BED,VCF", "C,F,BED,VCF"])
    out = tmp_path / "ref"
    _fake_steps(monkeypatch)
    assert cli.main(_build_cli(sheet, out)) == 0
    log = (out / "build.log").read_text()
    assert str(tmp_path) not in log and "input = s.csv (" in log


def test_path_redaction_replaces_whole_paths_only():
    names = {"/data/run1": "run1", "/data": "data"}
    text = ("in /data/run1/x.bed.gz, /data/other and /mnt/data/run1/y, "
            "/data2/z, /data/run1.bak")
    got = runtime.redact_paths(text, names)
    # /data/run1.bak is a longer name than /data/run1, but it lies in /data
    assert got == ("in run1/x.bed.gz, data/other and /mnt/data/run1/y, "
                   "/data2/z, data/run1.bak")


# ---- the ICR gene list is bound to the reference ---------------------------

GENE_LIST_WARNING = ("WARNING: the ICR gene list differs from the one used to "
                     "build the reference (GENE labels may differ)")


def _write_reference_json(ref, **kv):
    (Path(ref) / "reference.json").write_text(json.dumps(kv))


def test_gene_list_warning_only_when_the_reference_records_another_md5(
        tmp_path):
    lst = tmp_path / "genes.bed"
    lst.write_text("chr1\t1\t2\tG1\n")
    md5 = runtime.md5_file(str(lst))
    assert run.gene_list_warnings({"icr_gene_list_md5": md5}, str(lst)) == []
    assert run.gene_list_warnings({}, str(lst)) == []  # v2.2 layout, Zenodo
    lst.write_text("chr1\t1\t2\tG2\n")
    got = run.gene_list_warnings({"icr_gene_list_md5": md5}, str(lst))
    assert [f"WARNING: {w}" for w in got] == [GENE_LIST_WARNING]


def test_run_warns_when_its_gene_list_is_not_the_one_of_the_build(
        env, tmp_path, capsys):
    bundled = runtime.md5_file(str(data.asset("icr_gene_list.bed")))
    _write_reference_json(env.ref, icr_gene_list_md5=bundled)
    assert env.go("--outdir", str(tmp_path / "o1")) == 0
    assert GENE_LIST_WARNING not in capsys.readouterr().err
    other = tmp_path / "other.bed"
    other.write_text("chr1\t1\t2\tG1\n")
    assert env.go("--outdir", str(tmp_path / "o2"), "--set",
                  f"icr_gene_list={other}") == 0
    err = capsys.readouterr().err
    assert err.count(GENE_LIST_WARNING) == 1
    assert GENE_LIST_WARNING in (tmp_path / "o2" /
                                 "S1.log").read_text().split("[messages]")[0]
    # a reference without the key stays silent with any list
    _write_reference_json(env.ref, smooth_sigma=150)
    assert env.go("--outdir", str(tmp_path / "o3"), "--set",
                  f"icr_gene_list={other}") == 0
    assert "gene list differs" not in capsys.readouterr().err
