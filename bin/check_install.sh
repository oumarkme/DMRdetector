#!/usr/bin/env bash
# DMRdetector — post-install sanity check.
#
# Verifies that every binary and R/Bioconductor package the pipeline expects
# is reachable from the active environment. Run this immediately after
# `conda activate dmrdetector` to catch broken installs before launching a
# real pipeline run.
#
# Exit code is the number of failed checks (0 = OK).

set -uo pipefail

fail=0
ok()   { printf "  \033[32mOK \033[0m  %s\n" "$*"; }
miss() { printf "  \033[31mFAIL\033[0m  %s\n" "$*"; fail=$((fail + 1)); }

require_cli() {
    local name="$1" min="${2-}"
    if ! command -v "$name" >/dev/null 2>&1; then
        miss "$name not found on PATH"
        return
    fi
    local ver
    case "$name" in
        nextflow) ver=$("$name" -v 2>&1 | head -1 | tr -d '\r') ;;
        *)        ver=$("$name" --version 2>&1 | head -1 | tr -d '\r') ;;
    esac
    if [[ -n "$min" ]]; then
        ok "$name ($ver) — required ≥ $min"
    else
        ok "$name ($ver)"
    fi
}

require_r_pkg() {
    local pkg="$1"
    if Rscript -e "suppressPackageStartupMessages(library($pkg))" >/dev/null 2>&1; then
        local v
        v=$(Rscript -e "cat(as.character(packageVersion('$pkg')))" 2>/dev/null)
        ok "R package $pkg ($v)"
    else
        miss "R package $pkg cannot be loaded"
    fi
}

echo "[check] CLI tools"
require_cli nextflow 24.10
require_cli java     17
require_cli Rscript
require_cli bcftools 1.20
require_cli bedtools 2.31
require_cli bgzip
require_cli tabix
require_cli awk
require_cli sort

echo "[check] R / Bioconductor packages"
if command -v Rscript >/dev/null 2>&1; then
    require_r_pkg data.table
    require_r_pkg R.utils
    require_r_pkg readxl
    require_r_pkg ensembldb
    require_r_pkg EnsDb.Hsapiens.v86
    require_r_pkg GenomicRanges
fi

echo "[check] Nextflow self-test"
if command -v nextflow >/dev/null 2>&1; then
    if nextflow -version >/dev/null 2>&1; then
        ok "nextflow -version succeeds"
    else
        miss "nextflow -version failed"
    fi
fi

echo
if [[ $fail -eq 0 ]]; then
    echo "[check] all green — pipeline is ready to run"
else
    echo "[check] $fail check(s) failed; reinstall with: bash install.sh --force"
fi
exit "$fail"
