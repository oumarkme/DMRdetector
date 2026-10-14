#!/bin/bash
#
# Compiles the numerical kernels of DMRdetector into the cache folder $1 by
# running the demo of the README (generated data, three BED formats, both
# sexes, target-region report with haplotype tracks) and one large file per
# format, so that the image can carry the compiled kernels (see
# src/dmrdetector/util/jit.py).
#
# usage: warm_cache.sh CACHE_DIR SMOKE_DIR
#   CACHE_DIR  folder to fill (NUMBA_CPU_NAME of the caller decides the CPU
#              the kernels are compiled for)
#   SMOKE_DIR  folder with generate_data.py and run_smoke.py (tests/smoke)
#              and warm_big.py (env)
#
# The script exits with 0 and writes no SEED file when a step fails (for
# example on a build machine with less memory than `build` budgets), so that
# an image without the cache still builds; it then works, but compiles on
# first use. With DMRDETECTOR_REQUIRE_CACHE=1 (release builds) it exits with 1
# instead, which stops the build.
#
# The SEED file names the program version, the CPU name and a hash of the
# cache files; the program names its per-user cache folder from a hash of it
# (src/dmrdetector/util/jit.py), so two builds never share a folder.

readonly CACHE=$1
readonly SMOKE=$2
readonly PY=/opt/venv/bin/python
readonly CLI=/opt/venv/bin/dmrdetector

WORK=$(mktemp -d) || exit 1
trap 'rm -rf "${WORK}"' EXIT
mkdir -p "${CACHE}"
export NUMBA_CACHE_DIR="${CACHE}" DMRDETECTOR_NUMBA_SEED= TMPDIR="${WORK}"

# Leaves the cache folder empty and exits with 0, or with 1 when the cache is
# required.
fail() {
  rm -rf "${CACHE:?}"/*
  if [[ "${DMRDETECTOR_REQUIRE_CACHE:-0}" == "1" ]]; then
    echo "ERROR: kernel cache not built (DMRDETECTOR_REQUIRE_CACHE=1): $*" >&2
    exit 1
  fi
  echo "WARNING: kernel cache not built: $*" >&2
  exit 0
}

"${PY}" "${SMOKE}/generate_data.py" "${WORK}/demo" || fail "demo data"
cd "${WORK}/demo" || fail "demo folder"
"${CLI}" build --input controls.csv --outdir reference --threads 2 \
  || fail "build"
for fmt in 5col modkit pbcpg; do
  for case_sex in CASE01:F CASE02:M; do
    id=${case_sex%%:*}
    sex=${case_sex##*:}
    bed="data/${fmt}/${id}.bed.gz"
    "${CLI}" run --meth_bed_combined "${bed}" --vcf "data/5col/${id}.vcf" \
      --sex "${sex}" --sample-id "${id}" --reference reference \
      --gtf genes.gtf --target-region target_region.bed \
      --meth_bed_hap1 "${bed}" --meth_bed_hap2 "${bed}" \
      --outdir "out_${fmt}_${id}" --threads 2 || fail "run ${fmt} ${id}"
  done
done
# Files of several 4 MiB pieces: the kernels also see read-only buffers (see
# warm_big.py).
"${PY}" "${SMOKE}/warm_big.py" "${SMOKE}" big || fail "large files"
for fmt in 5col modkit pbcpg; do
  bed="big/BIG.${fmt}.bed.gz"
  "${CLI}" run --meth_bed_combined "${bed}" --vcf data/5col/CASE01.vcf \
    --sex F --sample-id BIG --reference reference --gtf genes.gtf \
    --target-region target_region.bed \
    --meth_bed_hap1 "${bed}" --meth_bed_hap2 "${bed}" \
    --outdir "out_big_${fmt}" --threads 2 || fail "run large ${fmt}"
done
version=$("${CLI}" --version) || fail "version"
[[ -n "$(find "${CACHE}" -name '*.nbi' -print -quit)" ]] \
  || fail "no kernel was compiled"
# hash of the names and contents of the cache files (names sorted in the C
# locale), so that the SEED of two builds differs
manifest=$(set -o pipefail; cd "${CACHE}" \
  && find . -type f \( -name '*.nbi' -o -name '*.nbc' \) -print0 \
  | LC_ALL=C sort -z | xargs -0 sha256sum | sha256sum | cut -d' ' -f1) \
  || fail "manifest"
printf '%s numba cpu %s\nmanifest sha256 %s\n' "${version}" \
  "${NUMBA_CPU_NAME:-host}" "${manifest}" >"${CACHE}/SEED"
chmod -R a+rX "${CACHE}"
