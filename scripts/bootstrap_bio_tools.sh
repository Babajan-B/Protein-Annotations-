#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONDA_BIN="${CONDA_BIN:-/Users/jaan/miniconda3/bin/conda}"
ENV_PREFIX="${ENV_PREFIX:-$ROOT_DIR/envs/bio-tools}"
CONDA_PKGS_DIRS="${CONDA_PKGS_DIRS:-$ROOT_DIR/.conda/pkgs}"
CONDA_ENVS_DIRS="${CONDA_ENVS_DIRS:-$ROOT_DIR/.conda/envs}"

mkdir -p "$CONDA_PKGS_DIRS" "$CONDA_ENVS_DIRS" "$ROOT_DIR/resources/tool-downloads" "$ROOT_DIR/resources/databases"

if [[ ! -x "$CONDA_BIN" ]]; then
  echo "conda not found at $CONDA_BIN" >&2
  exit 1
fi

echo "Creating/updating conda env at $ENV_PREFIX"
CONDA_PKGS_DIRS="$CONDA_PKGS_DIRS" CONDA_ENVS_DIRS="$CONDA_ENVS_DIRS" \
  "$CONDA_BIN" create -y -p "$ENV_PREFIX" -c conda-forge -c bioconda \
  python=3.11 hmmer hhsuite psipred blast clustalo mafft cmake numba llvmlite

echo "Installing EVcouplings into the same env when possible"
if ! PATH="$ENV_PREFIX/bin:$PATH" "$ENV_PREFIX/bin/pip" install --no-build-isolation evcouplings; then
  echo "EVcouplings could not be installed automatically in this environment." >&2
  echo "The rest of the toolchain is still usable; EVcouplings remains optional." >&2
fi

cat <<EOF

Openly installable tools have been prepared in:
  $ENV_PREFIX

Next:
1. Review resources/manual_downloads.md
2. Download the academic/manual tools into resources/tool-downloads
3. Copy .env.predictors.example to .env.predictors and set the *_COMMAND_TEMPLATE values
EOF
