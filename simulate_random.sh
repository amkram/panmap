#!/bin/bash

# Script to simulate <n> reads from a random node in a panman file
# with mutations applied (5 SNPs, 1 insertion, 2 deletions by default)
# Usage: ./simulate_random.sh <panman_file> <n_reads> [seed] [threads]
#
# Arguments:
#   panman_file: Path to the .panman file
#   n_reads: Number of reads to simulate
#   seed: (optional) Random seed for reproducibility
#   threads: (optional) Number of threads for simulation (default: 4)

set -e

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <panman_file> <n_reads> [seed] [threads]"
    echo "  panman_file: Path to the .panman file"
    echo "  n_reads: Number of reads to simulate"
    echo "  seed: (optional) Random seed for reproducibility"
    echo "  threads: (optional) Number of threads for simulation (default: 4)"
    exit 1
fi

PANMAN_FILE="$(realpath "$1")"
N_READS="$2"
SEED="${3:-$RANDOM}"
CPUS="${4:-4}"

# Check if panman file exists
if [ ! -f "$PANMAN_FILE" ]; then
    echo "Error: Panman file not found: $PANMAN_FILE"
    exit 1
fi

# Get script directory and binary paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ORIGINAL_DIR="$(pwd)"
PANMAP_BIN="$SCRIPT_DIR/build/bin/panmap"
SIMULATE_BIN="$SCRIPT_DIR/build/bin/simulate"

# Check if binaries exist
if [ ! -f "$PANMAP_BIN" ]; then
    echo "Error: panmap binary not found: $PANMAP_BIN"
    echo "Please build the project first (./make.sh or cmake + make)"
    exit 1
fi

if [ ! -f "$SIMULATE_BIN" ]; then
    echo "Error: simulate binary not found: $SIMULATE_BIN"
    echo "Please build the project first (./make.sh or cmake + make)"
    exit 1
fi

rm -f ./panmap_*.map
rm -f ./simulated_*

# Create temp directory for output
TEMP_DIR=$(mktemp -d -t simulate_random_XXXXXX)
trap "rm -rf $TEMP_DIR" EXIT

# Determine output file prefix early
PANMAN_BASENAME=$(basename "$PANMAN_FILE" .panman)
OUTPUT_PREFIX="simulated_${PANMAN_BASENAME}_${N_READS}reads_seed${SEED}"

echo "[simulate_random] =========================================="
echo "[simulate_random] Panman file: $PANMAN_FILE"
echo "[simulate_random] Number of reads: $N_READS"
echo "[simulate_random] Seed: $SEED"
echo "[simulate_random] Threads: $CPUS"
echo "[simulate_random] Mutations: 5 SNPs, 1 insertion, 2 deletions"
echo "[simulate_random] =========================================="
echo "[simulate_random] Output files will be created in current directory:"
echo "  $ORIGINAL_DIR/${OUTPUT_PREFIX}_R1.fastq"
echo "  $ORIGINAL_DIR/${OUTPUT_PREFIX}_R2.fastq"
echo "  $ORIGINAL_DIR/${OUTPUT_PREFIX}_node.fasta"
echo "[simulate_random] =========================================="
echo ""

# Step 1: Dump a random node sequence from the panman file
echo "[simulate_random] Step 1: Dumping random node from panman..."

# Clean up any old random node files to avoid confusion
PANMAN_DIR=$(dirname "$PANMAN_FILE")
rm -f "$PANMAN_DIR"/*.random.*.fa 2>/dev/null

"$PANMAP_BIN" "$PANMAN_FILE" --dump-random-node --seed "$SEED" -t 8 2>&1 | tee "$TEMP_DIR/dump_log.txt"

# panmap writes the random node file next to the .panman file, so look there
# The filename pattern is: {panman_file}.random.{node_name}.fa
RANDOM_NODE_FA=$(ls "$PANMAN_DIR"/"$(basename "$PANMAN_FILE")".random.*.fa 2>/dev/null | head -n1)

if [ -z "$RANDOM_NODE_FA" ]; then
    echo "Error: Failed to generate random node fasta"
    echo "Checked directory: $PANMAN_DIR"
    echo "Pattern: $(basename "$PANMAN_FILE").random.*.fa"
    ls -la "$PANMAN_DIR"/*.random.*.fa 2>/dev/null || echo "No random files found"
    exit 1
fi

# Copy to temp directory for simulate to use
cp "$RANDOM_NODE_FA" "$TEMP_DIR/"
RANDOM_NODE_FA="$TEMP_DIR/$(basename "$RANDOM_NODE_FA")"

echo "[simulate_random] Random node fasta: $RANDOM_NODE_FA"

# Extract node name from the fasta header
NODE_NAME=$(head -n1 "$RANDOM_NODE_FA" | sed 's/^>//')
echo "[simulate_random] Selected node: $NODE_NAME"

# Step 2: Simulate reads from the random node
echo "[simulate_random] Step 2: Simulating $N_READS reads with $CPUS threads..."
echo "[simulate_random] Mutations: 5 SNPs, 1 insertion, 2 deletions"

# Default parameters (matching Snakefile defaults)
MODEL="NovaSeq"

# Run simulate with mutations: 5 SNPs, 1 insertion, 2 deletions (3 indels total)
cd "$TEMP_DIR"
"$SIMULATE_BIN" --panmat "$PANMAN_FILE" --ref "$NODE_NAME" --out_dir "$TEMP_DIR" \
    --n_reads "$N_READS" --rep 1 --model "$MODEL" --cpus "$CPUS" --seed "$SEED" \
    --mutnum 5 1 2 --indel_len 1 9

# Find generated read files
R1_FASTQ=$(find "$TEMP_DIR" -name "*_R1.fastq" 2>/dev/null | head -n1)
R2_FASTQ=$(find "$TEMP_DIR" -name "*_R2.fastq" 2>/dev/null | head -n1)

if [ -z "$R1_FASTQ" ] || [ -z "$R2_FASTQ" ]; then
    echo "Error: Failed to generate read files"
    exit 1
fi

# Copy output files to original working directory
cp "$R1_FASTQ" "$ORIGINAL_DIR/${OUTPUT_PREFIX}_R1.fastq"
cp "$R2_FASTQ" "$ORIGINAL_DIR/${OUTPUT_PREFIX}_R2.fastq"
cp "$RANDOM_NODE_FA" "$ORIGINAL_DIR/${OUTPUT_PREFIX}_node.fasta"

# Get absolute paths of created files
FINAL_R1="$ORIGINAL_DIR/${OUTPUT_PREFIX}_R1.fastq"
FINAL_R2="$ORIGINAL_DIR/${OUTPUT_PREFIX}_R2.fastq"
FINAL_NODE="$ORIGINAL_DIR/${OUTPUT_PREFIX}_node.fasta"

echo ""
echo "[simulate_random] =========================================="
echo "[simulate_random] SUCCESS! Files created:"
echo "  $FINAL_R1"
echo "  $FINAL_R2"
echo "  $FINAL_NODE"
echo "[simulate_random] =========================================="
echo "[simulate_random] Node: $NODE_NAME"
echo "[simulate_random] Seed: $SEED (use this to reproduce)"
echo ""
echo "[simulate_random] Running placement..."
echo "[simulate_random] Command: build/bin/panmap $PANMAN_FILE ${OUTPUT_PREFIX}_R1.fastq ${OUTPUT_PREFIX}_R2.fastq --stop align -t 1 -k 27 -s 8 -l 1 -f"
echo "[simulate_random] =========================================="
echo ""

# Run placement
cd "$ORIGINAL_DIR"
"$PANMAP_BIN" "$PANMAN_FILE" "${OUTPUT_PREFIX}_R1.fastq" "${OUTPUT_PREFIX}_R2.fastq" --stop genotype -t 1 -k 27 -s 8 -l 1 -f

echo ""
echo "[simulate_random] =========================================="
echo "[simulate_random] Placement completed!"
echo "[simulate_random] Expected node: $NODE_NAME"
echo ""
echo "[simulate_random] Checking placement results..."

# Check if the placement file exists
PLACEMENT_FILE="$PANMAN_FILE.placement.tsv"
if [ -f "$PLACEMENT_FILE" ]; then
    echo "[simulate_random] Placement results from $PLACEMENT_FILE:"
    echo "[simulate_random] ----------------------------------------"
    cat "$PLACEMENT_FILE" | column -t
    echo "[simulate_random] ----------------------------------------"
    
    # Extract the best node(s) for each metric (may be comma-separated)
    # Use -F'\t' to properly handle tab-separated fields with empty values
    echo ""
    echo "[simulate_random] Summary:"
    RAW_NODES=$(grep "^raw" "$PLACEMENT_FILE" | awk -F'\t' '{print $4}')
    JACCARD_NODES=$(grep "^jaccard" "$PLACEMENT_FILE" | awk -F'\t' '{print $4}')
    COSINE_NODES=$(grep "^cosine" "$PLACEMENT_FILE" | awk -F'\t' '{print $4}')
    WEIGHTED_NODES=$(grep "^weighted_jaccard" "$PLACEMENT_FILE" | awk -F'\t' '{print $4}')
    
    # Helper function to check if expected node is in comma-separated list
    check_node_in_list() {
        local expected="$1"
        local node_list="$2"
        # Convert comma-separated list to space-separated and check each
        echo "$node_list" | tr ',' '\n' | grep -Fxq "$expected"
        return $?
    }
    
    if check_node_in_list "$NODE_NAME" "$RAW_NODES"; then
        echo "[simulate_random] ✓ Raw matches correctly placed at $NODE_NAME"
        [ "$RAW_NODES" != "$NODE_NAME" ] && echo "    (tied with: $RAW_NODES)"
    else
        echo "[simulate_random] ✗ Raw matches placed at $RAW_NODES (expected: $NODE_NAME)"
    fi
    
    if check_node_in_list "$NODE_NAME" "$JACCARD_NODES"; then
        echo "[simulate_random] ✓ Jaccard correctly placed at $NODE_NAME"
        [ "$JACCARD_NODES" != "$NODE_NAME" ] && echo "    (tied with: $JACCARD_NODES)"
    else
        echo "[simulate_random] ✗ Jaccard placed at $JACCARD_NODES (expected: $NODE_NAME)"
    fi
    
    if check_node_in_list "$NODE_NAME" "$COSINE_NODES"; then
        echo "[simulate_random] ✓ Cosine correctly placed at $NODE_NAME"
        [ "$COSINE_NODES" != "$NODE_NAME" ] && echo "    (tied with: $COSINE_NODES)"
    else
        echo "[simulate_random] ✗ Cosine placed at $COSINE_NODES (expected: $NODE_NAME)"
    fi
    
    if check_node_in_list "$NODE_NAME" "$WEIGHTED_NODES"; then
        echo "[simulate_random] ✓ Weighted Jaccard correctly placed at $NODE_NAME"
        [ "$WEIGHTED_NODES" != "$NODE_NAME" ] && echo "    (tied with: $WEIGHTED_NODES)"
    else
        echo "[simulate_random] ✗ Weighted Jaccard placed at $WEIGHTED_NODES (expected: $NODE_NAME)"
    fi
else
    echo "[simulate_random] Warning: Placement file not found: $PLACEMENT_FILE"
fi

echo "[simulate_random] =========================================="
