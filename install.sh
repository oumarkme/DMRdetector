#!/usr/bin/env bash
# DMRdetector installer — creates the `dmrdetector` conda environment.
#
# Usage:
#   bash install.sh              # create the env (skip if it exists)
#   bash install.sh --update     # update an existing env
#   bash install.sh --force      # remove and recreate
#
# Auto-detects mamba, micromamba, or conda. After install, prints next steps.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_FILE="$SCRIPT_DIR/environment.yml"
ENV_NAME=dmrdetector

mode="create"
case "${1:-}" in
    --update)  mode="update" ;;
    --force)   mode="force"  ;;
    --help|-h)
        sed -n '2,11p' "$0" | sed 's/^# *//'
        exit 0
        ;;
    "") ;;
    *)  echo "Unknown argument: $1 (use --help)" >&2; exit 2 ;;
esac

# Pick the fastest available solver.
if   command -v mamba       >/dev/null 2>&1; then solver=mamba
elif command -v micromamba  >/dev/null 2>&1; then solver=micromamba
elif command -v conda       >/dev/null 2>&1; then solver=conda
else
    echo "ERROR: need conda, mamba, or micromamba on PATH." >&2
    echo "Install miniforge: https://github.com/conda-forge/miniforge#install" >&2
    exit 1
fi
echo "[install] solver: $solver"

env_exists() {
    "$solver" env list 2>/dev/null | awk '{print $1}' | grep -qx "$ENV_NAME"
}

case "$mode" in
    force)
        if env_exists; then
            echo "[install] removing existing $ENV_NAME"
            "$solver" env remove -y -n "$ENV_NAME"
        fi
        echo "[install] creating $ENV_NAME from $ENV_FILE"
        "$solver" env create -y -f "$ENV_FILE"
        ;;
    update)
        if ! env_exists; then
            echo "ERROR: --update requested but env '$ENV_NAME' not found." >&2
            exit 1
        fi
        echo "[install] updating $ENV_NAME from $ENV_FILE"
        "$solver" env update -f "$ENV_FILE" --prune
        ;;
    create)
        if env_exists; then
            echo "[install] env '$ENV_NAME' already exists — skipping create"
            echo "[install] (use --update to refresh, --force to recreate)"
        else
            echo "[install] creating $ENV_NAME from $ENV_FILE (5–15 min)"
            "$solver" env create -y -f "$ENV_FILE"
        fi
        ;;
esac

echo
echo "[install] done. Next steps:"
echo "  conda activate $ENV_NAME"
echo "  bash $SCRIPT_DIR/bin/check_install.sh"
echo "  nextflow run $SCRIPT_DIR/main.nf --help"
