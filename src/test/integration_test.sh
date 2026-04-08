#!/bin/bash
# Integration tests for panmap
# Runs end-to-end pipeline and verifies outputs
set -euo pipefail

PANMAP="${1:?Usage: $0 <panmap-binary>}"
TESTDATA="$(dirname "$0")/data"
TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

PASS=0
FAIL=0

check() {
    local name="$1"; shift
    if "$@" >/dev/null 2>&1; then
        echo "  PASS  $name"
        PASS=$((PASS + 1))
    else
        echo "  FAIL  $name"
        FAIL=$((FAIL + 1))
    fi
}

echo "=== panmap integration tests ==="
echo "Binary: $PANMAP"
echo "Data:   $TESTDATA"
echo ""

# Test 1: Index build
echo "[1] Index build"
$PANMAP "$TESTDATA/rsv_4K.panman" --stop index -o "$TMPDIR/idx_test" -t 2 >/dev/null 2>&1
check "index file created" test -f "$TMPDIR/idx_test.idx"
check "index file non-empty" test -s "$TMPDIR/idx_test.idx"

# Test 2: Placement (known leaf node)
echo "[2] Placement - leaf node MZ515733.1"
$PANMAP "$TESTDATA/rsv_4K.panman" "$TESTDATA/MZ515733.1.fa" \
    --stop place -o "$TMPDIR/place_leaf" -t 2 -i "$TMPDIR/idx_test.idx" >/dev/null 2>&1
check "placement file created" test -f "$TMPDIR/place_leaf.placement.tsv"
check "placed to MZ515733.1" grep -q "MZ515733.1" "$TMPDIR/place_leaf.placement.tsv"
SCORE=$(awk -F'\t' '/^log_raw/ {print $2}' "$TMPDIR/place_leaf.placement.tsv")
check "log_raw score > 50" awk "BEGIN {exit ($SCORE > 50) ? 0 : 1}"

# Test 3: Full pipeline (placement + alignment + genotyping)
echo "[3] Full pipeline - leaf node"
$PANMAP "$TESTDATA/rsv_4K.panman" "$TESTDATA/MZ515733.1.fa" \
    --stop genotype -o "$TMPDIR/full_leaf" -t 2 -i "$TMPDIR/idx_test.idx" >/dev/null 2>&1
check "placement.tsv exists" test -f "$TMPDIR/full_leaf.placement.tsv"
check "bam exists" test -f "$TMPDIR/full_leaf.bam"
check "vcf exists" test -f "$TMPDIR/full_leaf.vcf"
check "ref.fa exists" test -f "$TMPDIR/full_leaf.ref.fa"
# Leaf node input = exact match, so 0 variant calls expected
NVARS=$(grep -cv '^#' "$TMPDIR/full_leaf.vcf" || true)
check "0 variants for exact match" test "$NVARS" -eq 0

# Test 4: Placement with internal node sequence
echo "[4] Placement - internal node"
# Pick a random node fasta from test data
NODE_FA=$(ls "$TESTDATA"/rsv_4K.panman.random.node_*.fa | head -1)
$PANMAP "$TESTDATA/rsv_4K.panman" "$NODE_FA" \
    --stop place -o "$TMPDIR/place_node" -t 2 -i "$TMPDIR/idx_test.idx" >/dev/null 2>&1
check "placement file created" test -f "$TMPDIR/place_node.placement.tsv"
NODE_SCORE=$(awk -F'\t' '/^log_raw/ {print $2}' "$TMPDIR/place_node.placement.tsv")
check "log_raw score > 0" awk "BEGIN {exit ($NODE_SCORE > 0) ? 0 : 1}"

# Test 5: Full pipeline with internal node (should produce variants)
echo "[5] Full pipeline - internal node"
$PANMAP "$TESTDATA/rsv_4K.panman" "$NODE_FA" \
    --stop genotype -o "$TMPDIR/full_node" -t 2 -i "$TMPDIR/idx_test.idx" >/dev/null 2>&1
check "vcf exists" test -f "$TMPDIR/full_node.vcf"
check "bam exists" test -f "$TMPDIR/full_node.bam"

# Test 6: Placement with FASTQ input (quality scores present)
echo "[6] Placement - fastq input"
$PANMAP "$TESTDATA/rsv_4K.panman" "$TESTDATA/MZ515733.1.fastq" \
    --stop place -o "$TMPDIR/place_fq" -t 2 -i "$TMPDIR/idx_test.idx" >/dev/null 2>&1
check "placement file created" test -f "$TMPDIR/place_fq.placement.tsv"
check "placed to MZ515733.1" grep -q "MZ515733.1" "$TMPDIR/place_fq.placement.tsv"
FQ_SCORE=$(awk -F'\t' '/^log_raw/ {print $2}' "$TMPDIR/place_fq.placement.tsv")
check "log_raw score > 50" awk "BEGIN {exit ($FQ_SCORE > 50) ? 0 : 1}"

# Test 7: Full pipeline with FASTQ (placement + alignment + genotyping)
echo "[7] Full pipeline - fastq"
$PANMAP "$TESTDATA/rsv_4K.panman" "$TESTDATA/MZ515733.1.fastq" \
    --stop genotype -o "$TMPDIR/full_fq" -t 2 -i "$TMPDIR/idx_test.idx" >/dev/null 2>&1
check "vcf exists" test -f "$TMPDIR/full_fq.vcf"
check "bam exists" test -f "$TMPDIR/full_fq.bam"

# Test 8: Version flag
echo "[8] CLI flags"
check "--version works" $PANMAP --version
check "--help works" $PANMAP --help

# Summary
echo ""
echo "=== Results: $PASS passed, $FAIL failed ==="
exit $FAIL
