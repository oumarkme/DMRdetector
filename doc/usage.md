# Usage reference

Details that the README leaves out: building the image, the demo data, mounting data into the container, cluster jobs, input checks, the settings file, building a reference, memory and run time, and the error messages with their fixes. Every `dmrdetector` command in this page is checked against the real flags by `tests/unit/test_docs_commands.py`.

## Building the image

Supported platform: Linux x86_64 and WSL2. The image is built for x86_64 only. macOS (Intel or Apple Silicon), ARM Linux and native Windows are not supported. On Apple Silicon, Docker Desktop can run an x86_64 image under emulation; we did not test this, and it is not suitable for real data. On WSL2 we checked only that the image starts (`dmrdetector --version`, Apptainer from conda); the other commands were run on Linux.

The build needs Apptainer 1.x or Docker and a network connection. About 1 GB of free space in `/tmp` should be enough (an estimate; the unpacked image is 464 MB). Set `APPTAINER_TMPDIR` if `/tmp` is small. On one server (Apptainer 1.4.5) the build took 3 min 39 s, including the compilation of the numerical kernels, and gave a 168 MB file. It runs the package smoke test at the end, which loads the compiled kernels from the image and checks that none is compiled again. Without the cache, the smoke test prints `absent` and the build still succeeds: a machine with less memory than `build` budgets (about 16 GB) cannot compile the kernels, and the image then compiles on first use. A release build fails in that case when it is given `--build-arg REQUIRE_CACHE=1` (Apptainer) or `--build-arg DMRDETECTOR_REQUIRE_CACHE=1` (Docker); the commands are in the README.

- `--fakeroot` refused (a FATAL message about fakeroot or `/etc/subuid`): try the same command without `--fakeroot`. It worked on our server (4 min 9 s), which allows unprivileged user namespaces. Otherwise ask your administrator for a subuid range, or build on any machine where you have root or fakeroot and copy `dmrdetector.sif` to the cluster. The `.sif` is one file.
- Docker: `docker build -f env/Dockerfile -t dmrdetector:3.2.0 .`, then `docker run --rm dmrdetector:3.2.0 --version`. We could not test this recipe because no Docker daemon was available to us. It repeats the steps of the Apptainer recipe, and a unit test checks that both recipes agree. Neither recipe sets `NUMBA_CACHE_DIR`; the program then uses a private folder per user in the temporary directory and copies the compiled kernels of the image into it (see "Numba cache" below). Report problems in an issue.
- Singularity: the `singularity` command of our Apptainer 1.4.5 installation is an alias of `apptainer`. We did not test SingularityCE itself, which reads the same `.def` file in principle.
- Same versions on every machine: the base image is pinned by its sha256 digest (Python 3.12.15, Debian trixie), and every Python package by exact version and wheel sha256 in `env/requirements.lock` (licenses and audit in [dependencies.md](dependencies.md)). The log of every `build` and `run` records the versions.

## Demo data

The Quick start of the README runs on the synthetic data set that `tests/smoke/generate_data.py` writes (8 controls and 2 cases on chr22, chrX and chrY). The thresholds in the demo `settings.toml` are relaxed because the data are tiny. Do not analyze real data from inside `demo/`, because that file is read automatically. With a pip installation, run `generate_data.py` with `venv/bin/python` and call `venv/bin/dmrdetector`; outside the image the numerical kernels compile on first use. In a fresh container with 8 threads, `build` took 62 s (median of 3) and `run` 4 s (3.0.8: 80 s and 12 s in a fresh cache folder, see "Numba cache").

`run` writes `dmrdetector/CASE01.dmr.tsv`, `CASE01.target_region.tsv` and `CASE01.log`. The DMR table has seven rows; the README shows two and leaves out `P_VALUE`. Two warnings are expected:

- The female demo samples (four controls and CASE01) have no heterozygous SNP on chrX in their synthetic VCF, so `build` and `run` print `WARNING: ... sex F, but the VCF has 0 heterozygous SNPs on chrX outside the PARs`. Real whole-genome VCFs of females have more than 40,000 ([sex_chromosomes.md](sex_chromosomes.md)).
- `--target-region builtin` selects the 21 regions packaged with the program. The demo reads lie on chr22 and none of the 21 regions does, so `run` prints `no target region of CASE01 has a CpG in the sample that is also in the reference atlas`, and the report lists every region with `NUM_CPG` 0 and `NA` in `P_FDR` and `DIRECTION`. Use `--target-region target_region.bed` (the chr22 region written to `demo/`) to see filled rows.

## Running the image

**Apptainer** mounts your home directory, `/tmp` and the current directory. Other folders are not visible until you add them with `--bind`. A file that exists on the host but is not mounted gives `file not found` (the message then ends with a mounting hint).

```bash
dmrdetector() { apptainer exec --bind /data,/refs /path/to/dmrdetector.sif dmrdetector "$@"; }
```

We checked these cases on one server (Apptainer 1.4.5). Paths with spaces work when the `--bind` argument is quoted. Relative paths work and count from the current directory, also in the samplesheet. A read-only data folder works (`--bind /data:/data:ro`), because the pipeline writes only to `--outdir` and `--tmp-dir`. A current directory that is not writable gives `cannot create the output or temporary folder` until you pass `--outdir` and `--tmp-dir` on a writable disk.

**Docker.** The image entry point is `dmrdetector`, so the arguments after the image name are the command:

```bash
docker run --rm -u "$(id -u):$(id -g)" -v /data:/data -w /data dmrdetector:3.2.0 run --help
```

This recipe has not been built by us; see the README.

**Numba cache.** The numerical kernels are compiled code that numba keeps in `NUMBA_CACHE_DIR`. The image carries the kernels, compiled at build time for a generic x86-64 CPU (`NUMBA_CPU_NAME=generic` in both recipes) and stored in `/opt/dmrdetector/numba-cache`. Numba cannot use a read-only cache folder, so on its first call `dmrdetector` copies the kernels (about 3 MB, files with mode 0600) into the cache folder in use, once per folder and image, and compiles only what the copy lacks. The folder in use is `NUMBA_CACHE_DIR` when you set it. Otherwise it is `dmrdetector-numba-<uid>-<id>` in the temporary directory, created with mode 0700 for the current user; `<id>` is a hash of the file `SEED` that the image build writes next to the kernels (it names the version, the CPU name and a hash of the kernel files), so two builds of the image never share a folder. Numba unpickles the files of its cache, so whoever can write there can run code as you. `dmrdetector` therefore uses the folder only when it belongs to you, is not a link and has no group or world write bit, and when the temporary directory is not writable by others without the sticky bit (like `/tmp`). Otherwise it makes a new private folder with mode 0700, copies the kernels into it, and removes it when the process ends; the per-user folder is never removed. If the copy fails, numba compiles as in 3.0.8. `DMRDETECTOR_NUMBA_SEED` names another folder to copy from; set it empty to skip the copy. On the demo data in a fresh cache folder, the first `run` took 3.9 s instead of 11.6 s with 3.0.8 and `build` took 62 s instead of 80 s. Apptainer shares the host `/tmp`, so the cache survives between runs until that folder is cleaned. If the node `/tmp` is wiped, or you want the cache elsewhere, set a folder you own:

```bash
mkdir -p "$HOME/.cache/dmrdetector"
export APPTAINERENV_NUMBA_CACHE_DIR="$HOME/.cache/dmrdetector"
```

In Docker, add `-e NUMBA_CACHE_DIR=/cache -v $HOME/.cache/dmrdetector:/cache`; without it the cache lives in the container and is lost with `--rm`. Because the kernels are compiled for a generic CPU, a cache folder shared by nodes with different CPU models is valid on all of them. The kernels give the same output as kernels compiled for the host CPU (see [validation.md](validation.md)).

## Cluster jobs (SLURM example)

`--threads` caps the threads of the whole pipeline; the default is half of the threads the job may use (CPU affinity is respected, and at least 2). On a shared server always set it, because the default counts the threads you may use, not the threads that are idle. Memory is budgeted as 8 GB per thread, which is a scheduling rule: the measured peaks are lower (see below). Put the temporary files on node-local scratch with `--tmp-dir` (by default they go to `TMPDIR` or `/tmp`) and bind that folder.

```bash
#!/bin/bash
#SBATCH --job-name=dmrdetector
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
set -euo pipefail

SIF=/path/to/dmrdetector.sif
SCRATCH="${TMPDIR:-/tmp}/dmr_${SLURM_JOB_ID:-test}"
mkdir -p "$SCRATCH" "$HOME/.cache/dmrdetector"
export APPTAINERENV_NUMBA_CACHE_DIR="$HOME/.cache/dmrdetector"

apptainer exec --bind /data,/refs,"$SCRATCH" "$SIF" dmrdetector run \
  --meth_bed_combined /data/P001.combined.bed.gz --vcf /data/P001.vcf --sex F \
  --reference /refs/reference --gtf /refs/Homo_sapiens.GRCh38.86.gtf.gz \
  --threads "${SLURM_CPUS_PER_TASK:-8}" --tmp-dir "$SCRATCH" --outdir /data/results
rm -rf "$SCRATCH"
```

We ran the commands of this script on a workstation with the variables unset; we did not run it under SLURM. `build` of a large atlas needs more memory (see Memory); give it `--cpus-per-task` and `--mem` of 8 GB per thread. Some sites forbid `apptainer build` on compute nodes and allow only `apptainer exec`; build the image once on a login or build host.

## Input checks

- **Genome build.** The run stops when a record ends beyond the length of its hg38 contig (chr1 to chr22, chrX, chrY; contigs named `chr1` or `1`). The check sees only the CpGs in the file. It catches an hg19 file that has CpGs where a contig is longer in hg19 than in hg38, and it cannot prove that a file is hg38.
- **modkit.** Run `modkit pileup --cpg --combine-strands` against an hg38 reference, so that each CpG is one row. The reader takes read depth (column 10) and modified reads (column 12) of every row and ignores the modification code (column 4) and the strand (column 6). A file with 5hmC (`h`) or 6mA (`a`) rows, or with one row per strand, would count CpGs twice; the run stops when such rows appear among the first 20,000 data lines. Remove other codes with `awk '$4 == "m"'`.
- **pb-CpG-tools.** Pass the combined BED (`<prefix>.combined.bed.gz`, 9 columns with a header line) to `--meth_bed_combined`. The haplotype BEDs of the same run go to `--meth_bed_hap1` and `--meth_bed_hap2`.
- **5-column BED.** `CHROM START END coverage modified_reads`, no header.
- **VCF.** A record counts when its genotype has two or more different called alleles (`0/1`, `0|1`, `1/2`). The genotype need not be phased. Indels, multiallelic records without a SNP allele, and homozygous or missing genotypes are ignored.
- **Sample ID.** Letters, digits, `.`, `_` and `-`. The default is the BED file name up to its first `.`; set `--sample-id` when that is not a valid ID.
- **Contig names of the VCF.** The program matches the VCF and the BED by exact contig name, as v2.2 did. When the VCF has heterozygous SNPs and none of its contig names is a contig of the BED (for example `1` against `chr1`), the program warns and removes no CpG. The result for that input does not change.
- **VCF with several samples.** Every sample column is read, and a CpG is removed when any sample is heterozygous there. The program warns when the VCF has more than one sample column.
- **Haplotype tracks.** The modkit and hg38 checks above also apply to `--meth_bed_hap1`, `--meth_bed_hap2` and to the tracks of the stand-alone `target` command, and `input_format` (`--input-format` for `target`) sets their columns. An unindexed track must be sorted by contig and START, with each contig in one block (`sort -k1,1 -k2,2n`); otherwise the report stops with `cannot index the track`. A `.tbi` that is older than its track is not used; the program makes a new index in the temporary folder.
- **Nothing to compare.** `run` and `build` warn when no CpG is left after filtering, when no CpG lies on chr1 to chr22 (for example RefSeq contig names such as `NC_000001.11`), and when a sex scope of the sample is empty. `run` stops with `0 sample CpGs match the reference atlas` when no CpG of the sample matches the autosomal atlas; the message gives the counts and the likely causes. In a sex scope it only warns.
- **Sex.** `run` and `build` warn when a sample labeled F has as many chrY CpGs as a male; see [sex_chromosomes.md](sex_chromosomes.md).

## Work folder and temporary files

`build` makes its work folder itself: `<outdir>/work`, or with `--tmp-dir` a folder `dmrdetector-build-<12 hex digits>` below it, where the digits are a hash of the resolved `--outdir`, so builds of different references can share one `--tmp-dir`. The folder holds a marker file `.dmrdetector-build` (version, outdir, process ID, time) and a lock file. An existing folder is used only when it holds the marker of a build of the same `--outdir`, or is empty. Otherwise `build` stops with an error that names the folder, and it never writes to or deletes that folder. A second `build` that finds the folder locked also stops. An interrupted `build` of the same `--outdir` resumes; a resume under another program version reruns every step, because each step stamp includes the version. After a success the folder is removed, only when its marker names this `--outdir`; `--keep-tmp` keeps it. A file system without `flock` gives a warning instead of a lock.

At the start, `build` and `run` remove the files `.tmp<pid>.<name>` that a killed process left in the output folder, only when no process with that ID exists on this host and the file is at least one hour old. They look nowhere else and do not enter subfolders. The `.<name>.part` copies of `run` are not matched. `run` warns when it replaces `<id>.dmr.tsv`, `<id>.log` or `<id>.target_region.tsv`, and when it leaves a `<id>.target_region.tsv` of an earlier run because it makes no target-region report; it deletes none of them.

## Settings file

`dmrdetector config init` writes `./settings.toml` with every option, its default value and a one-line comment. Edit the values you need and delete the lines you do not; a deleted line restores the default. `build` and `run` read `./settings.toml` when it exists, or the file given by `--config`. Precedence: flag, then `--set KEY=VALUE`, then file, then built-in default (a flag and `--set` for the same key is an error). Keys without a default (`reference`, `gtf`, `target_region`, `tmp_dir`, `icr_gene_list`, `par_bed`, `threads`) are commented out. `dmrdetector config show run` lists every key with its value, source and meaning. Unknown keys, wrong types and out-of-range values stop the run with a message that suggests the nearest valid key. The first lines of the log of `build` and `run` give the file in use (`settings: ` and its path, or `settings: built-in defaults`); a file that changes a default also gets a NOTICE that lists the changed keys. Set `target_region = "builtin"` (or pass `--target-region builtin`) to use the 21 regions packaged with the program; inside the image the packaged file is at a path that exists only there, so it cannot be passed from outside. A file named `builtin` in the working directory must be written as `./builtin`.

A file named `dmrdetector.toml` (the name in 3.0.0 and 3.0.1) is still read in this release, with the warning `dmrdetector.toml is deprecated; rename it to settings.toml`; when both files exist, `settings.toml` is used.

The imprinting gene list (`icr_gene_list`, default: the bundled `assets/icr_gene_list.bed`) is used by both commands. `build` takes the ICR panel from it and writes its md5 to `reference.json` as `icr_gene_list_md5`. `run` takes the `GENE` label of an ICR DMR from the list, so it prints `WARNING: the ICR gene list differs from the one used to build the reference (GENE labels may differ)` when the md5 of its list is another one. Set the same `icr_gene_list` for `build` and `run`. A reference without `icr_gene_list_md5` (a download, a v2.2 layout, a build before 3.0.8) is not checked.

```toml
[reference]
reference = "/refs/reference"
gtf       = "/refs/Homo_sapiens.GRCh38.86.gtf.gz"

[target]
target_region = "/refs/target_region.hg38.bed"

[resources]
threads = 36
```

```bash
dmrdetector run --meth_bed_combined P001.bed.gz --vcf P001.vcf --sex F --set dmr_fdr_threshold=0.05
```

## Memory and run time

The budget is 8 GB per thread; the minimum is 2 threads (16 GB), and `--threads 1` is raised to 2 with a warning. A `--threads` value above the available threads is lowered with a warning. The available threads are the CPU affinity, lowered to the cgroup quota of the process and of its parent groups. The budget is lowered to the smaller of MemAvailable in `/proc/meminfo` and the cgroup memory limit (`memory.max` of cgroup v2, `memory.limit_in_bytes` of v1). When the budget comes from the machine and a task does not fit, the message says so instead of asking for more threads. A worker that the kernel kills for lack of memory gives `a worker process was killed, probably out of memory; use fewer --threads`. A task starts only when the estimated peak memory of the tasks running together fits in the budget. The estimates are in `TASK_GB` of `src/dmrdetector/runtime.py` (preprocessing 2.5 GB, call-dmr 5 GB per sample and scope, gene table 2 GB, annotation 1 GB, target report 4 GB). They come from single-task peak measurements on one server.

Measured on a shared server (96 cores, 2 TB RAM) with one real ONT sample (pb-CpG-tools BED, 29,528,733 CpG rows) and a reference of four scopes. Memory is the peak of the summed resident set size of all processes, sampled every 0.5 s; shared pages are counted once per process, so it is an upper bound.

| `run` | Wall time | Peak memory |
|---|---|---|
| `--threads 2`, cold numba cache | 79 s | 2.4 GB |
| `--threads 8`, warm numba cache | 44 s | 6.6 GB |

`build` with 12 controls (four scopes and panels) took 143 to 146 s at 36 threads and peaked at 115 to 120 GB; at `--threads 2` it peaked at 13 GB (not timed). These `build` figures come from an earlier measurement whose sampling method is not recorded. Preprocessing of one sample with 4 threads took 15 s (memory not measured). `--version`, `--help` and `config` start in about 1 s. A large `build` can stop early with `raw atlas of N samples: needs at least X GB ... raise --threads`.

## Error messages

All user errors print one line that starts with `ERROR:` and exit with status 1. Command-line syntax errors come from the argument parser and exit with status 2.

| Message | Cause and fix |
|---|---|
| `--meth_bed_combined FILE is required` | The sample files are missing. Pass `--meth_bed_combined` and `--vcf`. |
| `--meth_bed_combined file not found: PATH` (also `--vcf`, `gtf`, `--target-region BED`, `samplesheet`) | The path is wrong, or the folder is not mounted in the container: add `--bind DIR` (Apptainer) or `-v DIR:DIR` (Docker). A relative path in a samplesheet counts from the working directory, not from the folder of the CSV; the message says so when the two folders differ. A device such as `/dev/null` gives `not a regular file`. |
| `--sex must be M or F, got 'X'` | Use `M` or `F` (lowercase is accepted). |
| `threads must be at least 2, got 0` | Use `--threads 2` or more, or delete the option for the default. |
| `no reference: pass --reference DIR` | Give the folder written by `build`, or the files of the Zenodo atlas. |
| `--reference DIR is not a directory` | The path is a file or does not exist. |
| `no autosomal reference atlas: looked for ref_atlas_150.bed.gz ...` | The folder is not a reference, or `smooth_sigma` differs from the one used by `build`. |
| `no gene annotation: pass --gtf FILE` | Pass an Ensembl GRCh38 GTF, or download one with `get-genes`. |
| `FILE: chrN has a record ending at ..., but chrN is ... bp in hg38/GRCh38` | The file is not hg38 (for example hg19). Align to GRCh38 and call the methylation again. |
| `FILE: this modkit bedMethyl file has rows for more than one modification code` | Keep only the 5mC rows (`awk '$4 == "m"'`) and use `--combine-strands`. |
| `preprocess: ValueError: FILE: input_format 5col needs at least 5 columns` | The BED has too few columns, or the format is mis-detected. Check the file, or set `input_format` in `settings.toml`. |
| `preprocess: an input file is truncated or not valid gzip` | An input file is damaged. Test each file with `gzip -t FILE` and copy it again. |
| `preprocess: ValueError: invalid file ... is it VCF/BCF format?` | The `--vcf` file is not a VCF. |
| `samplesheet row N (sample_id=ID): sex must be M or F` | Fix the `sex` column of `controls.csv`. |
| `sex (M or F) is required for the sex-specific scopes ...` | Add a `sex` column, or set `include_sex_chromosomes = false`. |
| `samplesheet row N: duplicate sample_id` | Every control needs its own ID. |
| `samplesheet FILE: the 'type' column has N value(s) that are neither a control nor CASE` | Every value of `type` must be `CTRL` or `control` (a control) or `CASE` (a row that `build` leaves out), in any case; an empty value is an error too, as in v2.2. The message lists the values and their first row. |
| `samplesheet FILE has a 'type' column but no row with type CTRL (or control)` | `build` reads only control rows of a sheet with a `type` column. Fix the `type` values, or delete the column when every row is a control. |
| `atlas: no CpG site passed the atlas filters` | Fewer than `atlas_min_samples_per_site` controls cover any site, or the files are not hg38 CpG tables. Add controls or check the inputs. |
| `--set KEY: unknown option (did you mean ...)` | A mistyped key. `dmrdetector config show run` lists the keys. |
| `--set KEY: must be float, got ...` / `must be between 0 and 1` | The value has the wrong type or range. |
| `estimated peak memory ... exceeds the budget` | Raise `--threads` (8 GB per thread). |
| `estimated peak memory ... exceeds the N GB that this machine (or its cgroup) has available` | The machine, not `--threads`, limits the budget. Free memory or use a larger machine. |
| `STEP: a worker process was killed, probably out of memory; use fewer --threads` | The kernel killed a worker of the step named first (`preprocess`, `atlas`, `smoothing`, `atlas reader`, ...). Lower `--threads`; the budget then admits fewer tasks at once. Since 3.0.10 every worker pool reports this; before, a killed worker of the atlas or smoothing steps could leave the program waiting for ever. |
| `work folder DIR exists and was not made by dmrdetector build` | `DIR` has files and no `.dmrdetector-build` marker. `build` leaves it alone: move or remove it, or choose another folder with `--tmp-dir`. |
| `another dmrdetector build is using the work folder DIR` | A second `build` with the same `--outdir` and work folder is running. Wait for it, or use another `--tmp-dir`. |
| `ID autosomal: 0 sample CpGs match the reference atlas` | No CpG of the sample is in the atlas. Check the contig names (`chr1` or `1`, not RefSeq names), that the sample is hg38, and `dmr_min_coverage`. |
| `FILE: cannot index the track` | A target-region track is not sorted by contig and START. Sort it with `sort -k1,1 -k2,2n`. |
| `target-region report failed` (WARNING, then `ERROR`, exit 1) | The DMR table `ID.dmr.tsv` was written and is complete; the target-region report was not. The log says `status: FAILED (the DMR table was published; ...)`. Fix the track and run again. |
| `cannot create the output or temporary folder` | The folder is not writable. Pass `--outdir` and `--tmp-dir` on a writable disk. |
| `get-genes: cannot download from ...` | No network, or a wrong release. Compute nodes often have no internet access: run it on a login node, or download the GTF elsewhere. |

Warnings go to the screen and to the log (`[warnings]` section). The common ones: `no target region of ID has a CpG in the sample that is also in the reference atlas` (the target-region report is still written, with `NUM_CPG` 0 in every row; check that the BED and the sample use the same contigs), `no CpG records in FILE` (the methylation BED is empty), `the 'type' column keeps the CTRL rows only; N of M rows are dropped` (`build`, a NOTICE: the samplesheet has a `type` column and N rows are not `CTRL` or `control`), `the ICR gene list differs from the one used to build the reference (GENE labels may differ)` (`run`: the list in use has another md5 than the list recorded in `reference.json`), `N of the kept methylation records have END - START != 1` (a CpG record is not one base long; see [ALGORITHM.md](ALGORITHM.md)), `DIR already holds a reference; its files are overwritten` (`build`), `scope male_X: 2 CTRL samples (< sex_atlas_min_controls 3); no reference for this scope` (`build` skips that scope), `the reference lacks: male_X (atlas), Y (atlas). Falling back to autosomal-only DMR calling` (`run` uses the autosomes only), and `--threads 1 raised to 2`, `--threads N lowered to M`, `memory budget lowered from X GB to Y GB`, `no CpG is left after filtering`, `no CpG on chr1-chr22 or the PARs`, `no CpG after filtering in the SCOPE scope`, `the VCF has heterozygous SNPs on contigs ..., but none of these names is a contig of the methylation BED`, `the VCF has N sample columns`, `the VCF has no sample column (sites-only?)` and `the VCF has no heterozygous SNP call` (the heterozygous-SNP filter removes no CpG), `sex F, but N CpGs remain on chrY`, `sex M, but the VCF has N heterozygous SNPs on chrX outside the PARs` and `sex F, but the VCF has N heterozygous SNPs on chrX outside the PARs` (the label disagrees with the chrX SNP count; see [sex_chromosomes.md](sex_chromosomes.md)), and `FILE (the DMR table of an earlier run) will be replaced`.

## Building a reference

Build the reference from controls of the same technology, base caller and tissue as the patients: `dmrdetector build --input controls.csv --outdir reference`. `build` overwrites an existing reference in `--outdir`. The columns of `controls.csv` are `sample_id,sex,meth_bed,vcf` with `sex` M or F; the handling of a `type` column is in the error table below. Relative paths count from the working directory, not from the folder of the CSV.

A sex-specific scope for chrX and chrY needs at least 3 controls of that sex (`sex_atlas_min_controls`); a scope with fewer is skipped with a warning, and `run` then analyzes the autosomes only. Three control counts are easy to confuse. The atlas keeps a site covered by at least 3 controls (`atlas_min_samples_per_site`). `cpg_min_samples_per_site` (5) is also a number of controls: CpG island panels use sites with more than 5 controls (so at least 6 per scope) and the ICR panel uses sites with at least 5. `icr_min_sites` (5) is a number of CpGs, the minimum per ICR region, and is unrelated to the number of controls. With fewer controls the panels are empty or `build` stops with `no CpG site passed the atlas filters`. We did not test how many controls give stable calls; the manuscript atlas has 69.

## Reference of the Zenodo record

The Zenodo record 10.5281/zenodo.18978867 holds the atlases of the manuscript as separate files (`ref_atlas_150*.bed.gz`, `ref_atlas_500*.bed.gz`, `cpg_signal*.bed`, `icr_region.bed`), with no `.tbi` index and no `reference.json`. The record has no description. A folder with `ref_atlas_150.bed.gz`, `ref_atlas_150_{male_X,female_X,Y}.bed.gz`, `cpg_signal*.bed` and `icr_region.bed` works as `--reference`; we tested this layout by removing the indexes and `reference.json` of a reference built with this version. The `ref_atlas_500*` files are for visualization and are not read by default. The files hold 25,719 CpG islands and 445 ICRs. `ref_atlas_150.bed.gz` has 28,624,561 sites on chr1 to chr22 (plus 59,926 pseudoautosomal sites) and a maximum `N` of 69, in line with the 69 controls of the manuscript.

## FAQ

- What coverage is needed? A CpG of the patient counts when its coverage is above 10 (`dmr_min_coverage`), and a control CpG when it is above 5 (`atlas_min_read_coverage`). Both are settings. We did not test how DMR calls change with sequencing depth.
- Must the data be phased? No. Phased hap1 and hap2 BEDs are used only for the optional target-region report.
- Other species or genome builds? Not supported. The panels, the contig lengths and the target regions are hg38.
- Why is the sample sex needed? chrX and chrY are tested in sex-specific scopes ([sex_chromosomes.md](sex_chromosomes.md)). Without `--sex`, only the autosomes are analyzed.
