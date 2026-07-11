#!/usr/bin/env bash
# Run the three README demos and diff against examples/expected/.
# Usage: examples/check_examples.sh [path/to/panmap]   (default: ./build/bin/panmap, else panmap on PATH)
# Exit 0 iff all match. Normalizes fields that vary between runs (see below).

set -uo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

# ---- locate panmap -----------------------------------------------------------
PANMAP="${1:-}"
if [[ -z "$PANMAP" ]]; then
  if [[ -x ./build/bin/panmap ]]; then PANMAP=./build/bin/panmap
  elif command -v panmap >/dev/null 2>&1; then PANMAP=panmap
  else echo "ERROR: panmap not found (build it or pass its path as \$1)" >&2; exit 2
  fi
fi
echo "Using panmap: $PANMAP"
"$PANMAP" --version >/dev/null 2>&1 || { echo "ERROR: '$PANMAP' is not runnable" >&2; exit 2; }

data_dir=examples/data
exp_dir=examples/expected
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

pass=0 fail=0
ok()   { echo "  PASS  $1"; pass=$((pass+1)); }
bad()  { echo "  FAIL  $1"; fail=$((fail+1)); }

# VCF: compare variant identity (CHROM/POS/REF/ALT + GT). bcftools QUAL/INFO stats aren't
# reproducible across htslib versions / CPU arch; consensus.fa pins that they're applied.
norm_vcf() { awk -F'\t' '!/^#/ { split($10, g, ":"); print $1"\t"$2"\t"$4"\t"$5"\t"g[1] }' "$1"; }

# Resolve the index column to read names (FASTQ write order is thread-dependent),
# emit sorted readname/node/taxon triples.
resolve_assignments() {
  local fq=$1 out=$2
  awk -F'\t' -v OFS='\t' '
    FNR==NR { if (FNR%4==1) { nm=$0; sub(/^@/,"",nm); names[c++]=nm } next }
    { node=$1; taxon=$2; n=split($4,idx,","); for (i=1;i<=n;i++) print names[idx[i]], node, taxon }
  ' "$fq" "$out" | sort
}

# ---- Demo 1: single-sample pipeline -----------------------------------------
echo
echo "[1/3] single-sample pipeline"
o="$work/isolate"
if "$PANMAP" "$data_dir/sars_20000_twilight_dipper.panman" \
     "$data_dir/isolate_R1.fastq.gz" "$data_dir/isolate_R2.fastq.gz" \
     -o "$o" >"$work/d1.log" 2>&1; then
  e="$exp_dir/single_sample"
  diff -q "$e/isolate.placement.tsv" "$o.placement.tsv" >/dev/null && ok "placement.tsv" || bad "placement.tsv"
  diff -q "$e/isolate.ref.fa"        "$o.ref.fa"        >/dev/null && ok "ref.fa"        || bad "ref.fa"
  if diff -q "$e/isolate.consensus.fa" "$o.consensus.fa" >/dev/null; then ok "consensus.fa"; else
    bad "consensus.fa"; echo "    consensus diff (< expected / > got):"; diff "$e/isolate.consensus.fa" "$o.consensus.fa" | head -12; fi
  if diff <(norm_vcf "$e/isolate.vcf") <(norm_vcf "$o.vcf") >/dev/null; then ok "vcf (variant records)"; else
    bad "vcf (variant records)"; echo "    vcf variant diff (< expected / > got):"; diff <(norm_vcf "$e/isolate.vcf") <(norm_vcf "$o.vcf") | sed 's/^/      /'; fi
else
  bad "single-sample run failed (see $work/d1.log)"; cat "$work/d1.log"
fi

# ---- Demo 2: metagenomic abundance ------------------------------------------
echo
echo "[2/3] metagenomic abundance (--meta)"
o="$work/example"
if "$PANMAP" "$data_dir/sars_20000_twilight_dipper.panman" $data_dir/sars20000_5hap_*.fastq.gz \
     --meta --threads 4 --em-delta-threshold 0.00001 --output "$o" >"$work/d2.log" 2>&1; then
  diff <(sort "$exp_dir/meta_abundance/example.mgsr.abundance.out") \
       <(sort "$o.mgsr.abundance.out") >/dev/null \
    && ok "abundance.out (haplotype fractions)" || bad "abundance.out (haplotype fractions)"
else
  bad "metagenomic run failed (see $work/d2.log)"; cat "$work/d2.log"
fi

# ---- Demo 3: filter and assign ----------------------------------------------
echo
echo "[3/3] filter and assign (--filter-and-assign)"
[[ -f "$data_dir/v_mtdna.panman" ]] || cp data/vertebrate_mtdna/v_mtdna.panman "$data_dir/"
o="$work/subsampled"
if "$PANMAP" "$data_dir/v_mtdna.panman" "$data_dir/subsampled.fastq.gz" --meta --filter-and-assign \
     -k 15 -s 8 -l 1 --discard 0.6 --dust 5 --taxonomic-metadata "$data_dir/v_mtdna.meta.tsv" \
     -t 4 --breadth-ratio --output "$o" >"$work/d3.log" 2>&1; then
  e="$exp_dir/filter_assign"
  # Assigned reads as a set (order is thread-dependent).
  diff <(paste - - - - < "$e/subsampled.mgsr.assignedReads.fastq" | sort) \
       <(paste - - - - < "$o.mgsr.assignedReads.fastq"           | sort) >/dev/null \
    && ok "assignedReads.fastq (read set)" || bad "assignedReads.fastq (read set)"
  # Read -> node assignment (indices resolved to names).
  diff <(resolve_assignments "$e/subsampled.mgsr.assignedReads.fastq" "$e/subsampled.mgsr.assignedReads.out") \
       <(resolve_assignments "$o.mgsr.assignedReads.fastq"           "$o.mgsr.assignedReads.out") >/dev/null \
    && ok "assignedReads.out (read->node)" || bad "assignedReads.out (read->node)"
  # Read -> LCA node assignment.
  diff <(resolve_assignments "$e/subsampled.mgsr.assignedReads.fastq" "$e/subsampled.mgsr.assignedReadsLCANode.out") \
       <(resolve_assignments "$o.mgsr.assignedReads.fastq"           "$o.mgsr.assignedReadsLCANode.out") >/dev/null \
    && ok "assignedReadsLCANode.out (read->LCA)" || bad "assignedReadsLCANode.out (read->LCA)"
else
  bad "filter-and-assign run failed (see $work/d3.log)"; cat "$work/d3.log"
fi

# ---- summary -----------------------------------------------------------------
echo
echo "-------------------------------------------"
echo "  $pass passed, $fail failed"
echo "-------------------------------------------"
[[ $fail -eq 0 ]] && echo "All README demos reproduced the expected output." || echo "Some checks failed (see logs in the temp dir printed above)."
exit $(( fail > 0 ? 1 : 0 ))
