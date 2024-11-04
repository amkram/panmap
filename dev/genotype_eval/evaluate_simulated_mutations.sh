#! /bin/bash

# Command line options
TREE_PATH=$1
NUM_MUTATIONS=$2
COVERAGES=$3
NUM_REPLICATES=$4
OUTPUT_PREFIX=$5

# Parse coverage string and convert to array of coverages and numbers of reads
IFS=',' read -r -a COVERAGE_ARRAY <<< "$COVERAGES"
READS_PER_COVERAGE=200
NUM_READS_ARRAY=()

for COVERAGE in "${COVERAGE_ARRAY[@]}"; do
    NUM_READS=$((COVERAGE * READS_PER_COVERAGE))
    NUM_READS_ARRAY+=($NUM_READS)
done

# Create output directory
OUTPUT_DIR=$(dirname "$OUTPUT_PREFIX")
OUTPUT_PREFIX=$(basename "$OUTPUT_PREFIX")
mkdir -p "$OUTPUT_DIR"

# Set up for running simulated data
PANMAT_PATH="../panmans/sars_20000.panman"
REF_STRAIN="USA/CA-CDPH-500076119/2022|OP769690.1|2022-08-27"

for i in "${!NUM_READS_ARRAY[@]}"; do
  NUM_READS=${NUM_READS_ARRAY[$i]}
  COVERAGE=${COVERAGE_ARRAY[$i]}
  simulate \
    --panmat "$PANMAT_PATH" \
    --ref "$REF_STRAIN" \
    --out_dir "$OUTPUT_DIR" \
    --prefix "$OUTPUT_PREFIX"_"$COVERAGE"x \
    --mutnum "$NUM_MUTATIONS" 0 0 \
    --mut_spec "$PANMAT_PATH.mm" \
    --rep "$NUM_REPLICATES" \
    --n_reads "$NUM_READS" \
    --cpus 4
done
