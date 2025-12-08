#!/bin/bash
set -e

# Define paths
PANMAP_BIN="./build/bin/panmap"
TEST_DATA_DIR="src/test/data"
PANMAN_FILE="${TEST_DATA_DIR}/test.pmat"
READS_FILE="${TEST_DATA_DIR}/test.fastq"
OUTPUT_DIR="verify_out"
PREFIX="${OUTPUT_DIR}/test"

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Run panmap with minimap2 alignment
echo "Running panmap with minimap2..."
${PANMAP_BIN} --aligner minimap2 --prefix ${PREFIX} ${PANMAN_FILE} ${READS_FILE}

# Check outputs
echo "Checking outputs..."

if [ -f "${PREFIX}.sam" ]; then
    echo "SAM file created: ${PREFIX}.sam"
    head -n 5 "${PREFIX}.sam"
else
    echo "Error: SAM file not created"
    exit 1
fi

if [ -f "${PREFIX}.bam" ]; then
    echo "BAM file created: ${PREFIX}.bam"
else
    echo "Error: BAM file not created"
    exit 1
fi

if [ -f "${PREFIX}.vcf" ]; then
    echo "VCF file created: ${PREFIX}.vcf"
    head -n 5 "${PREFIX}.vcf"
else
    echo "Error: VCF file not created"
    # Genotyping might be skipped if no variants found or configured to skip?
    # But default is to run genotyping.
    # If placement fails, alignment won't run.
    # If alignment fails (no mapped reads), genotyping might fail or produce empty VCF.
    exit 1
fi

echo "Verification successful!"
