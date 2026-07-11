#!/bin/bash
# End-to-end tests for panmap: drive the real binary through CLI pipelines and verify
# outputs. Run via: ctest -R e2e   (or directly: run_e2e.sh <panmap-binary>)
# No `set -e`: a crashing step must not abort the suite; we run every test, print a
# full summary, and surface the failing step's output (which otherwise goes to a log).
set -uo pipefail

PANMAP="${1:?Usage: $0 <panmap-binary>}"
TESTDATA="$(cd "$(dirname "$0")/.." && pwd)/data"
TMPDIR=$(mktemp -d)
cp "$TESTDATA/rsv_4K.panman" "$TMPDIR/"   # writable: index is built next to the panman
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

# Run a panmap pipeline step. On failure (including a crash) print its output, record
# a failure, and let the suite continue.
drive() {
    local out="$TMPDIR/_step.log"
    if "$@" >"$out" 2>&1; then return 0; fi
    local rc=$?
    echo "  !! step failed (exit $rc): $*"
    tail -20 "$out" | sed 's/^/     | /'
    FAIL=$((FAIL + 1))
    return 0
}

echo "=== panmap e2e tests ==="
echo "Binary: $PANMAP"
echo "Data:   $TESTDATA"
echo ""

echo "[1] Index build"
drive $PANMAP "$TMPDIR/rsv_4K.panman" --stop index -t 2
check "index file created" test -f "$TMPDIR/rsv_4K.panman.idx"
check "index file non-empty" test -s "$TMPDIR/rsv_4K.panman.idx"

echo "[2] Placement - leaf node MZ515733.1"
drive $PANMAP "$TMPDIR/rsv_4K.panman" "$TESTDATA/MZ515733.1.fa" \
    --stop place -o "$TMPDIR/place_leaf" -t 2
check "placement file created" test -f "$TMPDIR/place_leaf.placement.tsv"
check "placed to MZ515733.1" grep -q "MZ515733.1" "$TMPDIR/place_leaf.placement.tsv"
SCORE=$(awk -F'\t' '/^log_raw/ {print $2}' "$TMPDIR/place_leaf.placement.tsv")
check "log_raw score > 50" awk "BEGIN {exit ($SCORE > 50) ? 0 : 1}"

echo "[3] Full pipeline - leaf node"
drive $PANMAP "$TMPDIR/rsv_4K.panman" "$TESTDATA/MZ515733.1.fa" \
    --stop genotype -o "$TMPDIR/full_leaf" -t 2
check "placement.tsv exists" test -f "$TMPDIR/full_leaf.placement.tsv"
check "bam exists" test -f "$TMPDIR/full_leaf.bam"
check "vcf exists" test -f "$TMPDIR/full_leaf.vcf"
check "ref.fa exists" test -f "$TMPDIR/full_leaf.ref.fa"
# Leaf node input = exact match, so 0 variant calls expected
NVARS=$(grep -cv '^#' "$TMPDIR/full_leaf.vcf" || true)
check "0 variants for exact match" test "$NVARS" -eq 0

echo "[4] Placement - internal node"
NODE_FA=$(ls "$TESTDATA"/rsv_4K.panman.random.node_*.fa | head -1)
NODE_ID=$(basename "$NODE_FA" | sed -E 's/.*\.random\.(node_[0-9]+)\.fa/\1/')
drive $PANMAP "$TMPDIR/rsv_4K.panman" "$NODE_FA" \
    --stop place -o "$TMPDIR/place_node" -t 2
check "placement file created" test -f "$TMPDIR/place_node.placement.tsv"
# Feeding a node's own sequence should place it back on that node.
check "placed to $NODE_ID" grep -q "$NODE_ID" "$TMPDIR/place_node.placement.tsv"
NODE_SCORE=$(awk -F'\t' '/^log_raw/ {print $2}' "$TMPDIR/place_node.placement.tsv")
check "log_raw score > 50" awk "BEGIN {exit ($NODE_SCORE > 50) ? 0 : 1}"

# --force-leaf sends the alignment to the closest leaf, which differs from the internal
# node, so variants are expected. (Placement now allows internal nodes past --stop place,
# so without --force-leaf this input would place back on the node itself and call 0
# variants; the flag is required to exercise the genome-vs-leaf path.) The input is a
# single genome sequence (depth 1), so disable the consensus depth/QUAL gate.
echo "[5] Full pipeline - internal node (force-leaf)"
drive $PANMAP "$TMPDIR/rsv_4K.panman" "$NODE_FA" \
    --stop genotype --force-leaf --min-depth 1 --min-qual 0 -o "$TMPDIR/full_node" -t 2
check "vcf exists" test -f "$TMPDIR/full_node.vcf"
check "bam exists" test -f "$TMPDIR/full_node.bam"
NODE_NVARS=$(grep -cv '^#' "$TMPDIR/full_node.vcf" || true)
check "internal node produces variants" test "$NODE_NVARS" -gt 0

echo "[6] Placement - fastq input"
drive $PANMAP "$TMPDIR/rsv_4K.panman" "$TESTDATA/MZ515733.1.fastq" \
    --stop place -o "$TMPDIR/place_fq" -t 2
check "placement file created" test -f "$TMPDIR/place_fq.placement.tsv"
check "placed to MZ515733.1" grep -q "MZ515733.1" "$TMPDIR/place_fq.placement.tsv"
FQ_SCORE=$(awk -F'\t' '/^log_raw/ {print $2}' "$TMPDIR/place_fq.placement.tsv")
check "log_raw score > 50" awk "BEGIN {exit ($FQ_SCORE > 50) ? 0 : 1}"

echo "[7] Full pipeline - fastq"
drive $PANMAP "$TMPDIR/rsv_4K.panman" "$TESTDATA/MZ515733.1.fastq" \
    --stop genotype -o "$TMPDIR/full_fq" -t 2
check "vcf exists" test -f "$TMPDIR/full_fq.vcf"
check "bam exists" test -f "$TMPDIR/full_fq.bam"

# --meta auto-builds the .midx index and exercises the EM path.
echo "[8] Metagenomic abundance"
drive $PANMAP "$TMPDIR/rsv_4K.panman" --meta --stop index -t 2
check "mgsr index created" test -s "$TMPDIR/rsv_4K.panman.midx"
drive $PANMAP "$TMPDIR/rsv_4K.panman" "$TESTDATA/MZ515733.1.fastq" \
    --meta -o "$TMPDIR/meta" -t 2
ABUND="$TMPDIR/meta.mgsr.abundance.out"
check "abundance output created" test -s "$ABUND"
# All reads are from MZ515733.1, so it must be the dominant haplotype.
check "MZ515733.1 is dominant" bash -c "sort -t\$'\t' -k2 -g -r '$ABUND' | head -1 | grep -q 'MZ515733.1'"
# Proportions must sum to ~1.0.
SUM=$(awk -F'\t' '{s+=$2} END {print s}' "$ABUND")
check "proportions sum to ~1.0" awk "BEGIN {exit ($SUM > 0.99 && $SUM < 1.01) ? 0 : 1}"

echo "[9] Consensus stage - self-match"
drive $PANMAP "$TMPDIR/rsv_4K.panman" "$TESTDATA/MZ515733.1.fa" \
    --stop consensus -o "$TMPDIR/cons" -t 2
check "consensus.fa created" test -f "$TMPDIR/cons.consensus.fa"
check "ref.fa created" test -f "$TMPDIR/cons.ref.fa"
# Self-match => 0 variants => consensus must equal the placement reference.
check "consensus equals reference" python3 -c '
import sys
def rd(p): return "".join(l.strip() for l in open(p) if not l.startswith(">")).upper()
sys.exit(0 if rd(sys.argv[1]) == rd(sys.argv[2]) else 1)' "$TMPDIR/cons.consensus.fa" "$TMPDIR/cons.ref.fa"

echo "[10] Variant calling correctness - known SNPs"
python3 - "$TESTDATA/MZ515733.1.fa" "$TMPDIR/snp_reads.fastq" "$TMPDIR/truth.txt" <<'PY'
import sys
inp, reads_out, truth_out = sys.argv[1], sys.argv[2], sys.argv[3]
seq = [l.strip() for l in open(inp) if not l.startswith(">")]
s = list("".join(seq).upper())
truth = []
for p in (4000, 7000, 10000):           # 0-based, mid-genome
    ref = s[p]; alt = next(b for b in "ACGT" if b != ref); s[p] = alt
    truth.append((p + 1, ref, alt))      # VCF POS is 1-based
g = "".join(s)
L, step = 150, 5                          # ~30x tiling so each SNP has read support
with open(reads_out, "w") as o:
    n = 0
    for i in range(0, len(g) - L, step):
        o.write(f"@r{n}\n{g[i:i+L]}\n+\n{'I'*L}\n"); n += 1
with open(truth_out, "w") as t:
    for p, r, a in truth: t.write(f"{p}\t{r}\t{a}\n")
PY
drive $PANMAP "$TMPDIR/rsv_4K.panman" "$TMPDIR/snp_reads.fastq" \
    --stop genotype -o "$TMPDIR/snp" -t 2
check "vcf created" test -f "$TMPDIR/snp.vcf"
while IFS=$'\t' read -r POS REF ALT; do
    check "SNP ${REF}${POS}${ALT} called with correct POS/REF/ALT" \
        awk -v p="$POS" -v r="$REF" -v a="$ALT" -F'\t' \
        '!/^#/ && $2==p && $4==r && $5==a {ok=1} END{exit ok?0:1}' "$TMPDIR/snp.vcf"
done < "$TMPDIR/truth.txt"

# The two files are independent halves of the tiling, not true mate pairs, so this
# exercises the two-file input path rather than insert-size / paired-end handling.
echo "[11] Two FASTQ inputs (R1 + R2)"
python3 - "$TESTDATA/MZ515733.1.fa" "$TMPDIR/pe_R1.fastq" "$TMPDIR/pe_R2.fastq" <<'PY'
import sys
inp, r1p, r2p = sys.argv[1], sys.argv[2], sys.argv[3]
g = "".join(l.strip() for l in open(inp) if not l.startswith(">")).upper()
L, step = 150, 15
reads = [g[i:i+L] for i in range(0, len(g) - L, step)]
half = len(reads) // 2
with open(r1p, "w") as o:
    for n, r in enumerate(reads[:half]): o.write(f"@p{n}/1\n{r}\n+\n{'I'*L}\n")
with open(r2p, "w") as o:
    for n, r in enumerate(reads[half:2*half]): o.write(f"@p{n}/2\n{r}\n+\n{'I'*L}\n")
PY
drive $PANMAP "$TMPDIR/rsv_4K.panman" "$TMPDIR/pe_R1.fastq" "$TMPDIR/pe_R2.fastq" \
    --stop genotype -o "$TMPDIR/pe" -t 2
check "placement.tsv created" test -f "$TMPDIR/pe.placement.tsv"
check "placed to MZ515733.1" grep -q "MZ515733.1" "$TMPDIR/pe.placement.tsv"
check "bam created" test -f "$TMPDIR/pe.bam"
check "vcf created" test -f "$TMPDIR/pe.vcf"

echo "[12] Metagenomic mixture - 70/30 recovery"
python3 - "$TESTDATA/MZ515733.1.fa" "$TESTDATA/rsv_4K.panman.random.node_1330.fa" "$TMPDIR/mix.fastq" <<'PY'
import sys
def rd(p): return "".join(l.strip() for l in open(p) if not l.startswith(">")).upper()
a, b = rd(sys.argv[1]), rd(sys.argv[2])
out = open(sys.argv[3], "w")
def emit(g, n, pre):
    L = 150; step = max(1, (len(g) - L) // n); c = i = 0
    while c < n and i + L <= len(g):
        out.write(f"@{pre}{c}\n{g[i:i+L]}\n+\n{'I'*L}\n"); c += 1; i += step
emit(a, 700, "A"); emit(b, 300, "B"); out.close()
PY
drive $PANMAP "$TMPDIR/rsv_4K.panman" "$TMPDIR/mix.fastq" \
    --meta -o "$TMPDIR/mix" -t 2
MIX="$TMPDIR/mix.mgsr.abundance.out"
check "mixture abundance created" test -s "$MIX"
check "exactly 2 haplotypes" test "$(grep -c . "$MIX")" -eq 2
check "MZ515733.1 ~0.70 (majority)" \
    awk -F'\t' '$1=="MZ515733.1" && $2>0.55 && $2<0.82 {ok=1} END{exit ok?0:1}' "$MIX"
check "node_1330 ~0.30 (minority)" \
    awk -F'\t' '$1=="node_1330" && $2>0.18 && $2<0.45 {ok=1} END{exit ok?0:1}' "$MIX"
MSUM=$(awk -F'\t' '{s+=$2} END{print s}' "$MIX")
check "proportions sum to ~1.0" awk "BEGIN {exit ($MSUM > 0.99 && $MSUM < 1.01) ? 0 : 1}"

echo "[13] CLI flags"
check "--version works" $PANMAP --version
check "--help works" $PANMAP --help

echo ""
echo "=== Results: $PASS passed, $FAIL failed ==="
exit $FAIL
