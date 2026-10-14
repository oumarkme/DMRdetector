"""Options of `build` and `run`: defaults, settings file, `--set`, flags.

OPTIONS is the single source of every default (the values of v2.2
`nextflow.config`). `effective` merges flag or --set, settings file and default;
`default_toml` writes settings.toml; `config show` lists the result.
"""

from __future__ import annotations

import argparse
from collections.abc import Callable
from dataclasses import dataclass
import difflib
import os
import sys
import tomllib


BUILD = ("build",)
RUN = ("run",)
BOTH = ("build", "run")

INPUT_FORMATS = ("auto", "5col", "modkit", "pbcpg")
LOG_INPUT_PATH_CHOICES = ("name", "full")


@dataclass(frozen=True)
class Opt:
    """One option: its place in settings.toml, default and limits."""

    section: str
    key: str  # v2.2 parameter name where one exists
    default: object  # None = unset (commented out in the TOML file)
    type: type  # bool, int, float or str
    # one-line meaning (`config show`, the commented TOML file, --help of a
    # dedicated flag)
    help: str
    commands: tuple = BOTH
    # True: the option also has a dedicated command-line flag (the others:
    # settings file or --set)
    flag: bool = False
    lo: float | None = None  # smallest accepted value (numbers)
    hi: float | None = None  # largest accepted value (numbers)


FRAC = dict(lo=0.0, hi=1.0)
COUNT = dict(lo=0)

# yapf: disable
OPTIONS: list[Opt] = [
    # [input]
    Opt("input", "input_format", "auto", str,
        "methylation BED format: auto (detect), 5col, modkit (bedMethyl, 18 "
        "columns) or pbcpg (pb-CpG-tools, 9 columns)"),
    # [reference]
    Opt("reference", "reference", None, str,
        "reference directory written by `build` (v2.2 layout accepted)",
        RUN, flag=True),
    Opt("reference", "gtf", None, str,
        "Ensembl GRCh38 GTF (plain or .gz); `dmrdetector get-genes` "
        "downloads one",
        RUN, flag=True),
    # [atlas]
    Opt("atlas", "atlas_min_samples_per_site", 3, int,
        "keep a CpG only when at least this many CTRL samples cover it",
        BUILD, lo=1),
    Opt("atlas", "atlas_min_read_coverage", 5, int,
        "keep a CpG of a CTRL sample only when its read coverage is above "
        "this",
        BUILD, **COUNT),
    Opt("atlas", "atlas_shrinkage_weight", 0.1, float,
        "weight of the global trend in the shrunken per-site dispersion (0 "
        "to 1)",
        BUILD, **FRAC),
    Opt("atlas", "smooth_sigma", 150, int,
        "Gaussian smoothing sigma of the reference atlas (bp); also names "
        "ref_atlas_<sigma>.bed.gz",
        lo=1),
    # [cpg_islands]
    Opt("cpg_islands", "cpg_min_samples_per_site", 5, int,
        "number of CTRL samples per atlas site: CpG island panels use sites "
        "with more than this many (ICR panel: at least this many); not the "
        "CpG count of icr_min_sites",
        BUILD, **COUNT),
    Opt("cpg_islands", "cpg_core_max_methylation", 0.1, float,
        "island core: reference methylation below this",
        BUILD, **FRAC),
    Opt("cpg_islands", "cpg_core_max_dispersion", 0.05, float,
        "island core: reference dispersion (PHI) below this",
        BUILD, **FRAC),
    Opt("cpg_islands", "cpg_shore_max_methylation", 0.3, float,
        "island shore: reference methylation below this",
        BUILD, **FRAC),
    Opt("cpg_islands", "cpg_max_gap_distance", 150, int,
        "largest gap (bp) between CpGs inside one island",
        BUILD, lo=1),
    Opt("cpg_islands", "cpg_min_island_length", 150, int,
        "shortest island (bp)",
        BUILD, lo=1),
    # [icr]
    Opt("icr", "icr_ref_meth_min", 0.4, float,
        "ICR: reference methylation at least this (hemi-methylated)",
        BUILD, **FRAC),
    Opt("icr", "icr_ref_meth_max", 0.6, float,
        "ICR: reference methylation at most this",
        BUILD, **FRAC),
    Opt("icr", "icr_min_sites", 5, int,
        "CpGs per ICR: the ICR panel needs at least this many, an ICR DMR "
        "more than this many (not a number of controls; see "
        "cpg_min_samples_per_site)",
        **COUNT),
    Opt("icr", "icr_diff_threshold", 0.25, float,
        "ICR DMR calls: the absolute methylation difference to the "
        "reference must exceed this",
        RUN, **FRAC),
    Opt("icr", "icr_v1_duplicates", True, bool,
        "count a CpG once per overlapping candidate gene as v2.2 did "
        "(false: once)",
        BUILD),
    # [dmr]
    Opt("dmr", "dmr_min_coverage", 10, int,
        "keep a CpG of a tested sample only when its read coverage is above "
        "this",
        RUN, **COUNT),
    Opt("dmr", "dmr_cpg_diff_threshold", 0.4, float,
        "CpG island DMR calls: the absolute methylation difference must "
        "exceed this",
        RUN, **FRAC),
    Opt("dmr", "dmr_merge_distance", 150, int,
        "CpGs at most this far apart (bp) form one DMR (also groups ICR "
        "CpGs)",
        **COUNT),
    Opt("dmr", "dmr_fdr_threshold", 0.01, float,
        "report DMRs with FDR (Benjamini-Hochberg) below this",
        RUN, **FRAC),
    Opt("dmr", "dmr_cpg_min_sites", 30, int,
        "CpG island DMR calls: a DMR needs more than this many CpGs",
        RUN, **COUNT),
    Opt("dmr", "dmr_dedupe_sites", False, bool,
        "read each CpG once when panel regions overlap (false reproduces "
        "v2.2, which repeats CpGs and inflates NUM_CPG)",
        RUN),
    # [sex]
    Opt("sex", "include_sex_chromosomes", True, bool,
        "analyze chrX and chrY in sex-specific scopes (male_X, female_X, "
        "Y); false restricts to chr1-22"),
    Opt("sex", "sex_atlas_min_controls", 3, int,
        "minimum CTRL samples of a sex to build a sex-specific scope",
        BUILD, lo=1),
    Opt("sex", "sex_f_ref_meth_min", 0.4, float,
        "female_X islands: reference methylation at least this "
        "(hemi-methylated)",
        BUILD, **FRAC),
    Opt("sex", "sex_f_ref_meth_max", 0.6, float,
        "female_X islands: reference methylation at most this",
        BUILD, **FRAC),
    Opt("sex", "sex_f_diff_threshold", 0.25, float,
        "female_X DMR calls: the absolute methylation difference must "
        "exceed this",
        RUN, **FRAC),
    Opt("sex", "sex_f_min_sites", 5, int,
        "female_X DMR calls: a DMR needs more than this many CpGs",
        RUN, **COUNT),
    # [target]
    Opt("target", "target_region", None, str,
        "BED of target regions, or \"builtin\" for the 21 regions packaged "
        "with the program; `run` then also writes "
        "<sample_id>.target_region.tsv",
        RUN, flag=True),
    Opt("target", "target_haplotypes", True, bool,
        "include hap1/hap2 tracks when --meth_bed_hap1/--meth_bed_hap2 are "
        "given",
        RUN),
    Opt("target", "target_min_coverage", 8, int,
        "merged track: keep a CpG only when coverage is above this",
        RUN, **COUNT),
    Opt("target", "target_hap_min_coverage", 4, int,
        "hap1/hap2 tracks: keep a CpG only when coverage is above this",
        RUN, **COUNT),
    # [paths]
    Opt("paths", "outdir", "dmrdetector", str,
        "folder of the output files of `run` (`build` takes --outdir on the "
        "command line)",
        RUN, flag=True),
    Opt("paths", "tmp_dir", None, str,
        "folder for the temporary files (fast local disk; default: the "
        "system temporary directory, for `build` <outdir>/work)",
        flag=True),
    Opt("paths", "log_input_paths", "name", str,
        "paths in <sample_id>.log and build.log: name (the last component, "
        "and the size of an input file; safe to share) or full (absolute "
        "paths)"),
    Opt("paths", "icr_gene_list", None, str,
        "imprinting gene BED (default: bundled assets/icr_gene_list.bed); "
        "`build` stores its md5 in reference.json and `run` warns when its "
        "list differs"),
    Opt("paths", "par_bed", None, str,
        "hg38 PAR coordinates BED (default: bundled assets/par_hg38.bed)"),
    # [resources]
    Opt("resources", "threads", None, int,
        "maximum threads the pipeline may use at any moment, over all "
        "worker processes and their threads; memory is budgeted as 8 GB per "
        "thread. Minimum 2 (1 is raised to 2). Default: 50% of the "
        "available threads, at least 2",
        flag=True),
    Opt("resources", "keep_tmp", False, bool,
        "keep the temporary folder and print its path (`build`: keep the "
        "work folder after success)",
        flag=True),
]
# yapf: enable

SECTIONS = list(dict.fromkeys(o.section for o in OPTIONS))
BY_KEY = {o.key: o for o in OPTIONS}
assert len(BY_KEY) == len(OPTIONS), "duplicate option key"

# module constants: ATLAS_MIN_SAMPLES_PER_SITE, SMOOTH_SIGMA, DMR_FDR_THRESHOLD,
# ...
for _o in OPTIONS:
    globals()[_o.key.upper()] = _o.default

DEFAULT_FILE = "settings.toml"
# value of `target_region` that selects the BED packaged with the program
BUILTIN_TARGET = "builtin"
# read with a warning in 3.0.2, to be removed in the next release
LEGACY_FILE = "dmrdetector.toml"


class ConfigError(Exception):
    """Invalid settings file or value; shown as `ERROR: ...` by the commands."""


def _check_value(where: str, o: Opt, v):
    """Checks type, choice and range of one value; returns it (float cast).

    Args:
        where: Source of the value, for the error message.
        o: The option.
        v: The value.

    Raises:
        ConfigError: The value is invalid.
    """
    if o.type is bool:
        ok = isinstance(v, bool)
    elif o.type is int:
        ok = isinstance(v, int) and not isinstance(v, bool)
    elif o.type is float:
        ok = isinstance(v, (int, float)) and not isinstance(v, bool)
        v = float(v) if ok else v
    else:
        ok = isinstance(v, str)
    if not ok:
        raise ConfigError(
            f"{where} must be {o.type.__name__}, got {v!r} ({type(v).__name__})"
            + (" (use true or false)" if o.type is bool else ""))
    check_choice(o, v, where)
    check_range(o, v, where)
    return v


def check_range(o: Opt, v: float, where: str = "") -> None:
    """Raises ConfigError when a number is outside the limits of `o`."""
    if o.type in (int, float) and ((o.lo is not None and v < o.lo) or
                                   (o.hi is not None and v > o.hi)):
        lo = "-inf" if o.lo is None else f"{o.lo:g}"
        hi = "inf" if o.hi is None else f"{o.hi:g}"
        raise ConfigError(
            f"{where or o.key} must be between {lo} and {hi}, got {v!r}")


def check_choice(o: Opt, v: str, where: str = "") -> None:
    """Raises ConfigError when a text option has a value outside its choices.

    The options with choices are `input_format` and `log_input_paths`.
    """
    choices = {
        "input_format": INPUT_FORMATS,
        "log_input_paths": LOG_INPUT_PATH_CHOICES
    }.get(o.key)
    if choices is not None and v not in choices:
        raise ConfigError(f"{where or o.key} must be one of "
                          f"{', '.join(choices)}, got {v!r}")


def load_file(path: str) -> dict:
    """Reads and validates a settings file.

    Returns:
        {key: value} for the keys the file sets.

    Raises:
        ConfigError: The file is missing, not valid TOML, or has an unknown
            section or key, or a wrong type or range.
    """
    try:
        with open(path, "rb") as fh:
            raw = tomllib.load(fh)
    except FileNotFoundError:
        raise ConfigError(f"config file not found: {path}") from None
    except tomllib.TOMLDecodeError as e:
        raise ConfigError(f"{path}: not valid TOML: {e}") from None
    out = {}
    for section, table in raw.items():
        if section not in SECTIONS:
            hint = difflib.get_close_matches(section, SECTIONS, 1)
            raise ConfigError(f"{path}: unknown section [{section}]" +
                              (f" (did you mean [{hint[0]}]?)" if hint else
                               f" (sections: {', '.join(SECTIONS)})"))
        if not isinstance(table, dict):
            raise ConfigError(f"{path}: '{section}' must be a [section] table")
        in_section = {o.key: o for o in OPTIONS if o.section == section}
        for key, v in table.items():
            if key not in in_section:
                if key in BY_KEY:
                    raise ConfigError(
                        f"{path}: key '{key}' belongs to "
                        f"[{BY_KEY[key].section}], not [{section}]")
                hint = difflib.get_close_matches(key, in_section, 1)
                raise ConfigError(
                    f"{path}: unknown key '{key}' in [{section}]" +
                    (f" (did you mean '{hint[0]}'?)" if hint else ""))
            out[key] = _check_value(f"{path}: [{section}] {key}",
                                    in_section[key], v)
    return out


METAVAR = {
    "reference": "DIR",
    "gtf": "FILE",
    "target_region": "FILE.bed",
    "outdir": "DIR",
    "tmp_dir": "DIR",
    "threads": "N"
}
FLAG_HELP = {  # one line each in --help; the full meaning is in `config show`
    "reference": "reference directory from `build` (v2.2 layout ok)",
    "gtf": "Ensembl GRCh38 GTF, plain or .gz",
    "target_region": "target BED, or builtin (packaged 21 regions): also "
                     "write <id>.target_region.tsv",
    "outdir": "output folder [default: dmrdetector]",
    "tmp_dir": "temporary files folder [default: system temp]",
    "threads": "max threads (8 GB each, min 2) [default: half]",
    "keep_tmp": "keep the temporary folder, print its path",
}


class HelpFormatter(argparse.HelpFormatter):
    """Wide help, so that each flag fits on one line."""

    def __init__(self, prog):
        """Sets a wide help position so that each flag fits on one line."""
        super().__init__(prog, max_help_position=30, width=80)


def flag(key: str) -> str:
    """Returns the dedicated flag of an option: tmp_dir -> --tmp-dir."""
    return "--" + key.replace("_", "-")


def add_flag(group: argparse._ArgumentGroup,
             key: str,
             help: str | None = None) -> None:
    """Adds the dedicated flag of option `key` to an argparse group.

    Defaults are None so that the settings file and the built-in default can be
    told apart.
    """
    o = BY_KEY[key]
    help_ = (help or FLAG_HELP.get(key, o.help)).replace("%", "%%")
    if o.type is bool:
        group.add_argument(flag(key),
                           dest=key,
                           action="store_true",
                           default=None,
                           help=help_)
    else:
        group.add_argument(flag(key),
                           dest=key,
                           type=o.type,
                           default=None,
                           metavar=METAVAR.get(key, o.type.__name__.upper()),
                           help=help_)


def add_common(group: argparse._ArgumentGroup) -> None:
    """Adds `--config` and `--set` (`build`, `run`, `config show`)."""
    group.add_argument("--config",
                       metavar="FILE",
                       help=f"settings file [default: ./{DEFAULT_FILE}]")
    group.add_argument("--set",
                       metavar="KEY=VALUE",
                       action="append",
                       default=[],
                       help="one-off change, e.g. dmr_fdr_threshold=0.05")


def _toml_literal(text: str):
    """Returns the TOML value of `text`, or `text` when it is not TOML."""
    try:
        return tomllib.loads(f"v = {text}")["v"]
    except tomllib.TOMLDecodeError:
        return text


def parse_set(items: list[str], command: str) -> dict:
    """Parses the `--set KEY=VALUE` items.

    The value is a TOML literal (0.05, true, "text"; a bare word is a string for
    a text option), checked for key, type and range.

    Returns:
        {key: value}.

    Raises:
        ConfigError: Bad item, unknown key or key of the other command, a key
            given twice, or an invalid value.
    """
    out: dict = {}
    for item in items or []:
        key, eq, text = item.partition("=")
        key, text = key.strip().replace("-", "_"), text.strip()
        if not eq or not key:
            raise ConfigError(f"--set {item!r}: expected KEY=VALUE, "
                              "e.g. --set dmr_fdr_threshold=0.05")
        valid = [o.key for o in OPTIONS if command in o.commands]
        if key not in valid:
            if key in BY_KEY:
                raise ConfigError(f"--set {key}: this option applies to "
                                  f"{' and '.join(BY_KEY[key].commands)} only, "
                                  f"not to `{command}`")
            # low cutoff: dmr_fdr must find dmr_fdr_threshold (ratio 0.58)
            hint = difflib.get_close_matches(key, valid, 1, cutoff=0.5)
            raise ConfigError(f"--set {key}: unknown option" +
                              (f" (did you mean '{hint[0]}'?)" if hint else
                               "; `dmrdetector config show` lists the keys"))
        o = BY_KEY[key]
        if key in out:
            raise ConfigError(f"--set {key} is given twice")
        v = _toml_literal(text)
        if o.type is str and not isinstance(v, str):
            v = text
        out[key] = _check_value(f"--set {key}:", o, v)
    return out


def _warn_stderr(msg: str) -> None:
    print(f"WARNING: {msg}", file=sys.stderr)


def default_file(cwd: str | None = None,
                 warn: Callable[[str], None] = _warn_stderr) -> str | None:
    """Returns the settings file of a folder, or None.

    ./settings.toml, else the deprecated ./dmrdetector.toml (with a warning).
    """
    cwd = cwd or os.getcwd()
    new, old = os.path.join(cwd, DEFAULT_FILE), os.path.join(cwd, LEGACY_FILE)
    if os.path.isfile(new):
        if os.path.isfile(old):
            warn(f"{LEGACY_FILE} is deprecated; rename it to {DEFAULT_FILE} "
                 f"or delete it. Both exist in {cwd}: {DEFAULT_FILE} is used "
                 f"and {LEGACY_FILE} is ignored")
        return new
    if os.path.isfile(old):
        warn(f"{LEGACY_FILE} is deprecated; rename it to {DEFAULT_FILE}. "
             "It is read in this release only")
        return old
    return None


def effective(
    args: argparse.Namespace,
    command: str,
    cwd: str | None = None,
    warn: Callable[[str],
                   None] = _warn_stderr) -> tuple[dict, dict, str | None]:
    """Merges flag or --set, settings file and default for one command.

    Args:
        args: Parsed arguments of the command.
        command: "build" or "run".
        cwd: Folder searched for the settings file (default: the current one).
        warn: Receives the deprecation warning (default: stderr); `build` and
            `run` pass a function that also writes it to the log.

    Returns:
        (values, sources, settings file). A source is "flag" (dedicated flag
        or --set), "file" or "default".

    Raises:
        ConfigError: A flag and --set give the same key, or a value is invalid.
    """
    path = getattr(args, "config", None)
    if path is None:
        path = default_file(cwd, warn)
    file_vals = load_file(path) if path else {}
    set_vals = parse_set(getattr(args, "set", None), command)
    vals, src = {}, {}
    for o in OPTIONS:
        if command not in o.commands:
            continue
        cli = getattr(args, o.key, None) if o.flag else None
        if cli is not None and o.key in set_vals:
            raise ConfigError(
                f"{flag(o.key)} and --set {o.key} both given; use one of them")
        if cli is None and o.key in set_vals:
            cli = set_vals[o.key]
        if cli is not None:
            v, s = cli, "flag"
            check_choice(o, v, flag(o.key))
            check_range(o, v, flag(o.key))
        elif o.key in file_vals:
            v, s = file_vals[o.key], "file"
        else:
            v, s = o.default, "default"
        vals[o.key], src[o.key] = v, s
    return vals, src, path


def _toml_value(v) -> str:
    """Returns the TOML text of a default."""
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, str):
        return '"' + v.replace("\\", "\\\\").replace('"', '\\"') + '"'
    return repr(v)


def default_toml() -> str:
    """Returns the text of settings.toml.

    Every option has its default and a one-line comment. This is the only
    generator of the file (`dmrdetector config init`); a unit test checks the
    copy in the repository root against it.
    """
    lines = [
        "# DMRdetector settings: every option of `dmrdetector build` and "
        "`dmrdetector run` with its default value.",
        f"# The file is read from ./{DEFAULT_FILE} in the working directory, "
        "or from the file given by --config FILE.",
        "# Precedence: flag > --set KEY=VALUE > this file > built-in default "
        "(a flag and --set for the same key is an error).",
        "# Delete a line to restore the default of its key. Keys without a "
        "default are commented out; remove the # to set them.",
        "# `dmrdetector config show run` lists the values in force. Create "
        "this file with `dmrdetector config init`.",
        "# Key names are the parameter names of the earlier R/Nextflow version "
        "(v2.2) where one exists.",
    ]
    for section in SECTIONS:
        lines += ["", f"[{section}]"]
        for o in OPTIONS:
            if o.section != section:
                continue
            used = "build and run" if len(o.commands) == 2 else o.commands[0]
            lines.append(f"# {o.help} (used by: {used})")
            if o.default is None:
                blank = '""' if o.type is str else "0"
                lines.append(f"# {o.key} = {blank}")
            else:
                lines.append(f"{o.key} = {_toml_value(o.default)}")
    return "\n".join(lines) + "\n"


def main_init(args: argparse.Namespace) -> int:
    """Runs `config init`: writes settings.toml (or stdout with --out -)."""
    text = default_toml()
    if args.out == "-":
        sys.stdout.write(text)
        return 0
    if os.path.exists(args.out) and not args.force:
        print(
            f"ERROR: {args.out} exists; use --force to overwrite or --out "
            "to choose another file",
            file=sys.stderr)
        return 1
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    from dmrdetector import runtime  # runtime imports this module
    runtime.write_atomic(args.out, text)
    print(f"wrote {args.out}", file=sys.stderr)
    return 0


def main_show(args: argparse.Namespace) -> int:
    """Runs `config show`: every option with its value, source and meaning."""
    try:
        vals, src, path = effective(args, args.for_command)
    except ConfigError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    print(f"# `{args.for_command}`: value, source (default, file or flag) "
          "and meaning of every option")
    print(f"# settings file: {path or '(none; built-in defaults)'}")
    for section in SECTIONS:
        keys = [o for o in OPTIONS if o.section == section and o.key in vals]
        if not keys:
            continue
        print(f"\n[{section}]")
        for o in keys:
            v = "(unset)" if vals[o.key] is None else repr(vals[o.key])
            print(f"{o.key} = {v}  # {src[o.key]}: {o.help}")
    return 0


def add_parser(
        subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """Adds the `config` command to `subparsers` and returns its parser."""
    p = subparsers.add_parser(
        "config",
        help="write or inspect the settings file (settings.toml)",
        description=
        "Settings file (settings.toml), dedicated flags and `--set KEY=VALUE` "
        "of `build` and `run`.")
    sub = p.add_subparsers(dest="config_cmd", required=True)
    i = sub.add_parser(
        "init", help="write settings.toml with every option and its default")
    i.add_argument("--out",
                   default=DEFAULT_FILE,
                   help="output file, or - for stdout [default: ./%s]" %
                   DEFAULT_FILE)
    i.add_argument("--force",
                   action="store_true",
                   help="overwrite an existing file")
    i.set_defaults(func=main_init)
    s = sub.add_parser(
        "show",
        help="print every option of `build` or `run` with its value, source "
        "and "
        "meaning")
    s.add_argument("for_command",
                   nargs="?",
                   default="run",
                   choices=["build", "run"],
                   metavar="{build,run}",
                   help="command to show [default: run]")
    add_common(s)
    s.set_defaults(func=main_show)
    return p
