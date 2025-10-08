#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="trying"
THREADS="${THREADS:-1}"
TARGETS="${*:-all}"

# Initialize micromamba shell hook if needed
if ! command -v micromamba >/dev/null 2>&1; then
  echo "micromamba not found in PATH" >&2
  exit 1
fi

# Activate env (works in non-interactive script)
eval "$(micromamba shell hook -s bash)"
micromamba activate "${ENV_NAME}" || { echo "Failed to activate env ${ENV_NAME}"; exit 1; }

echo "[run] Using env ${ENV_NAME}" 
which python || true
snakemake --version || { echo "snakemake not installed in env ${ENV_NAME}"; exit 2; }

# Dry run first for feedback
echo "[run] Dry run: snakemake -n -j${THREADS} ${TARGETS}" >&2
snakemake -n -j"${THREADS}" ${TARGETS}

# Real run
echo "[run] Executing: snakemake -j${THREADS} ${TARGETS}" >&2
snakemake -j"${THREADS}" ${TARGETS}
