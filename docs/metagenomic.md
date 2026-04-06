# Metagenomic Mode

Metagenomic mode (`--meta`) scores reads from a mixture sample against every node in the PanMAN to estimate haplotype abundance or assign reads directly to nodes.

---

## Abundance estimation

Estimate lineage abundances from a mixed sample (e.g., SARS-CoV-2 in wastewater).

!!! note "Prerequisites"
    This walkthrough requires minimap2, samtools, ivar, and pangolin in addition to panmap.

### 1. Preprocess reads

```bash
mkdir example_run && cd example_run

# Fetch reads
prefetch SRR19707934 && fasterq-dump SRR19707934

# Align to reference
minimap2 -a --sam-hit-only --MD -2 -x sr \
  ../examples/data/NC_045512v2.fa \
  SRR19707934_1.fastq SRR19707934_2.fastq \
  | samtools sort -@ 8 -o SRR19707934.sorted.bam

# Trim primers
ivar trim -e \
  -b ../examples/data/SNAPaddtlCov.bed \
  -p SRR19707934.trimmed.bam \
  -i SRR19707934.sorted.bam \
  -q 1 -m 80 -x 3
samtools sort -@ 8 -o SRR19707934.trimmed.sorted.bam SRR19707934.trimmed.bam

# Prepare amplicon stacks and trimmed reads
python3 ../scripts/trim_and_stack_amplicon.py stack \
  SRR19707934.sorted.bam ../examples/data/SNAPaddtlCov.bed \
  -o SRR19707934.amplicon_stacks.tsv
samtools fastq --no-sc SRR19707934.trimmed.sorted.bam \
  > SRR19707934.trimmed.fastq
```

### 2. Run panmap

```bash
# Build index
panmap ../examples/data/sars_20000_twilight_dipper.panman \
  --index-mgsr sars_20000_twilight_dipper.idx

# Estimate abundances
panmap ../examples/data/sars_20000_twilight_dipper.panman \
  SRR19707934.trimmed.fastq \
  --meta \
  --index sars_20000_twilight_dipper.idx \
  --amplicon-depth SRR19707934.amplicon_stacks.tsv \
  --mask-reads-relative-frequency 0.01 \
  em-delta-threshold 0.00001 \
  --output SRR19707934 -t 8
```

Output: `SRR19707934.mgsr.abundance.out` -- estimated node abundances.

### 3. Identify lineages

```bash
panmap ../examples/data/sars_20000_twilight_dipper.panman \
  --dump-sequences "$(cut -f1 SRR19707934.mgsr.abundance.out \
    | tr ',' '\n' | grep '\S' | tr '\n' ' ' | sed 's/ $//g')" \
  --output SRR19707934

pangolin SRR19707934.dump-sequences.fa \
  --outfile SRR19707934.pangolin.csv

python3 ../scripts/assign_lineage.py \
  SRR19707934.mgsr.abundance.out \
  SRR19707934.pangolin.csv \
  > SRR19707934.mgsr.lineage.abundance.out
```

---

## Filter and assign reads

Assign reads directly to pangenome nodes. Useful for ancient or environmental DNA (aeDNA) where the goal is species identification rather than abundance estimation.

### 1. Build index

For aeDNA reads, use `-k 15 -s 8 -l 1`:

```bash
mkdir example_run && cd example_run

panmap ../examples/data/v_mtdna.panman \
  --index-mgsr v_mtdna.idx \
  -k 15 -s 8 -l 1
```

### 2. Run filter and assign

```bash
panmap ../examples/data/v_mtdna.panman \
  ../examples/data/subsampled.fastq.gz \
  --meta \
  -i v_mtdna.idx \
  --filter-and-assign \
  --discard 0.6 --dust 5 \
  --taxonomic-metadata ../examples/data/v_mtdna.meta.tsv \
  -t 4 --breadth-ratio \
  --output subsampled
```

### Output files

| File | Contents |
|------|----------|
| `.mgsr.assignedReads.fastq` | Assigned reads |
| `.mgsr.assignedReads.out` | Read count and indices per node |
| `.mgsr.assignedReadsLCANode.out` | Read count and indices for the LCA of each read's assigned nodes |

---

## Metagenomic options reference

### Indexing

| Option | Description | Default |
|--------|-------------|---------|
| `--index-mgsr <file>` | Build/rebuild MGSR index at this path | -- |
| `--index-full` | Build full index (default is lite) | off |
| `--index-packed` | Build packed Cap'n Proto message | off |
| `--read-packed` | Read packed Cap'n Proto message | off |
| `--no-progress` | Disable progress bars | off |

### EM algorithm

| Option | Description | Default |
|--------|-------------|---------|
| `--top-oc` | Top N nodes by overlap coefficient to send to EM | `1000` |
| `--mask-reads` | Mask reads with k-min-mer total occurrence <= threshold | `0` |
| `--mask-seeds` | Mask k-min-mer seeds with total occurrence <= threshold | `0` |
| `--amplicon-depth <file>` | Amplicon depth TSV for frequency-based masking | -- |
| `--mask-reads-relative-frequency` | Mask reads with k-min-mer relative frequency < threshold * amplicon_depth | `0.0` |
| `--mask-seeds-relative-frequency` | Mask seeds with relative frequency < threshold * amplicon_depth | `0.0` |
| `--em-convergence-threshold` | Converge when likelihood difference < threshold | `0.00001` |
| `--em-delta-threshold` | Converge when max proportion change < threshold | `0.0` |
| `--em-maximum-rounds` | Maximum EM rounds | `5` |
| `--em-maximum-iterations` | Maximum EM iterations per round | `1000` |
| `--em-leaves-only` | Only run EM on leaf (sample) nodes | off |

### Filter and assign

| Option | Description | Default |
|--------|-------------|---------|
| `--filter-and-assign` | Assign reads to nodes without EM | off |
| `--dust` | Discard reads with DUST score > threshold | `100.0` (no filtering) |
| `--discard` | Discard reads with parsimony score < threshold * total seeds | `0.0` (no discard) |
| `--mask-read-ends` | Mask N bases from read ends (for aeDNA damage) | `0` |
| `--taxonomic-metadata <file>` | TSV with taxonomic metadata per node | -- |
| `--maximum-families` | Discard reads spanning > N taxonomic families | `1` |
| `--ambiguous-score-threshold-ratio` | Discard reads scoring outside max families by ratio | `0.0` |
| `--ambiguous-score-threshold` | Discard reads scoring outside max families by absolute value | `0` |
| `--breadth-ratio` | Calculate observed/expected breadth ratio | off |
| `--pseudochain` | Use pseudo-chains for read scoring | off |
| `--batch-files-path <file>` | TSV file containing batch file paths | -- |
| `--batch-size` | Batch size for filtering and assigning | `1000000` |
