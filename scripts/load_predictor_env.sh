#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_FILE="${1:-$ROOT_DIR/.env.predictors}"

if [[ ! -f "$ENV_FILE" ]]; then
  echo "Predictor env file not found: $ENV_FILE" >&2
  exit 1
fi

set -a
source "$ENV_FILE"
set +a

echo "Loaded predictor environment from $ENV_FILE"
