#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="goodpan"
ENV_PATH="/home/alex/micromamba/envs/goodpan"
THREADS="${THREADS:-24}"
TARGETS="${*:-all}"  # Default to 'all' rule, not --rerun-incomplete

# Simply add the environment bin directory to PATH
export PATH="${ENV_PATH}/bin:$PATH"

# Verify the environment is accessible
if [ ! -f "${ENV_PATH}/bin/python" ]; then
    echo "Error: goodpan environment not found at ${ENV_PATH}" >&2
    exit 1
fi

echo "[run] Using env ${ENV_NAME}" 
echo "[run] Python: $(which python)"
echo "[run] Snakemake: $(snakemake --version)"
echo "[run] ISS: $(which iss)"

# Real run with rerun-incomplete and increased latency wait
# The build_all_indexes rule will be built automatically as a dependency
echo "[run] Executing: snakemake -j${THREADS} --latency-wait 60 --rerun-incomplete ${TARGETS}" >&2
snakemake -j"${THREADS}" --latency-wait 60 --rerun-incomplete ${TARGETS}
