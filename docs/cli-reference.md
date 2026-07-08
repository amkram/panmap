# CLI Reference

```
panmap <panman> [reads1.fq] [reads2.fq] [options]
```

Use `--help` for common options or `--help-all` for the full list. This page
mirrors `panmap --help-all` for version 0.1.2.

## General

| Option | Description | Default |
|--------|-------------|---------|
| `-h, --help` | Show common options | -- |
| `--help-all` | Show all options | -- |
| `-V, --version` | Show version | -- |
| `-o, --output` | Output file prefix | derived from reads filename |
| `-t, --threads` | Number of threads | `1` |
| `--stop` | Stop after: `index`, `place`, `align`, `genotype`, `consensus` | `consensus` |
| `--meta` | Enable metagenomic mode | off |
| `-a, --aligner` | `minimap2` or `bwa` | `minimap2` |
| `-v, --verbose` | Verbose output | off |
| `-q, --quiet` | Errors only | off |
| `--no-color` | Disable colored output | off |

## Index

| Option | Description | Default |
|--------|-------------|---------|
| `-i, --index` | Load a pre-built index from this path | auto-built at `<panman>.idx` (`.midx` under `--meta`) |
| `--index-out` | Write the built index to this path | next to the panman |
| `-f, --reindex` | Force rebuild index | off |
| `--index-packed` | Build packed Cap'n Proto message | off |
| `--read-packed` | Read packed Cap'n Proto message | off |
| `--zstd-level` | ZSTD compression level for index (1-22) | `7` |

## Seed parameters

| Option | Description | Default |
|--------|-------------|---------|
| `-k, --kmer` | Syncmer k-mer length | `19` |
| `-s, --syncmer` | Syncmer s parameter | `8` |
| `-l, --lmer` | Syncmers per seed | `3` |
| `--offset` | Syncmer offset | `0` |
| `--open-syncmer` | Use open syncmers | off |
| `--hpc` | Homopolymer-compressed seeds | off |
| `--flank-mask` | Mask bp at ends of genomes | `250` |
| `--seed-mask-fraction` | Mask top seed fraction | `0` |
| `--min-seed-quality` | Minimum Phred score for seed region | `0` |
| `--min-read-support` | Minimum reads for a seed (2 = filter singletons) | `1` |
| `--extent-guard` | Guard seed deletions at genome extent boundaries | off |

## Read processing

| Option | Description | Default |
|--------|-------------|---------|
| `--trim-start` | Trim bases from read start | `0` |
| `--trim-end` | Trim bases from read end | `0` |
| `--dedup` | Deduplicate reads | off |

## Placement

| Option | Description | Default |
|--------|-------------|---------|
| `--force-leaf` | Restrict placement to leaf nodes (auto-enabled unless `--stop place`) | off |
| `--refine` | Alignment-based refinement of top candidates | off |
| `--refine-top-pct` | Top fraction of nodes to refine | `0.01` |
| `--refine-max-top-n` | Max nodes to align against | `150` |
| `--refine-neighbor-radius` | Expand to neighbors within N branches | `2` |
| `--refine-max-neighbor-n` | Max additional nodes from neighbor expansion | `150` |

## Genotyping

| Option | Description | Default |
|--------|-------------|---------|
| `--impute` | Impute N's from parent (skip `_->N` mutations) | off |
| `--no-mutation-spectrum` | Disable mutation-spectrum filtering in VCF | off |
| `--baq` | Enable Base Alignment Quality in mpileup | off |

## Metagenomic: EM algorithm

| Option | Description | Default |
|--------|-------------|---------|
| `--top-oc` | Top N nodes by overlap coefficient to send to EM | `1000` |
| `--mask-reads` | Mask reads whose k-min-mers have total occurrence ≤ threshold | `0` |
| `--mask-seeds` | Mask k-min-mer seeds with total occurrence ≤ threshold | `0` |
| `--amplicon-depth` | Amplicon-depth TSV for frequency-based masking | -- |
| `--mask-reads-relative-frequency` | Mask reads with relative frequency < threshold × amplicon depth | `0` |
| `--mask-seeds-relative-frequency` | Mask seeds with relative frequency < threshold × amplicon depth | `0` |
| `--em-convergence-threshold` | Converge when the likelihood difference < threshold | `1e-05` |
| `--em-delta-threshold` | Converge when the max proportion change < threshold | `0` |
| `--em-maximum-rounds` | Maximum EM rounds | `5` |
| `--em-maximum-iterations` | Maximum EM iterations per round | `1000` |
| `--em-leaves-only` | Only run EM on leaf (sample) nodes | off |
| `--no-progress` | Disable progress bars | off |

## Metagenomic: filter and assign

| Option | Description | Default |
|--------|-------------|---------|
| `--filter-and-assign` | Assign reads to nodes without running EM | off |
| `--dust` | Discard reads with a PRINSEQ DUST score > threshold | `100` (no filtering) |
| `--discard` | Discard reads with max parsimony score < threshold × total seeds | `0` (no discard) |
| `--mask-read-ends` | Mask N bases from both read ends (for aeDNA damage) | `0` |
| `--taxonomic-metadata` | TSV with taxonomic metadata per node | -- |
| `--taxonomic-rank` | Taxonomic rank (column in the metadata TSV) to filter/assign on | `Family` |
| `--maximum-taxon-number` | Discard reads spanning more than N distinct taxa at that rank | `1` |
| `--ambiguous-score-threshold-ratio` | Discard reads scoring outside the max-scoring taxa by this ratio | `0` |
| `--ambiguous-score-threshold` | Discard reads scoring outside the max-scoring taxa by this absolute value | `0` |
| `--breadth-ratio` | Compute observed/expected breadth ratio | off |
| `--pseudochain` | Use pseudo-chains for read scoring | off |
| `--batch-size` | Batch size for filtering and assigning | `1000000` |

## Batch mode

| Option | Description | Default |
|--------|-------------|---------|
| `--batch` | File listing samples, one per line: `reads1 [reads2] [output_prefix]` (works in normal and `--meta` modes) | -- |

## Developer

Low-level flags for debugging and reproducibility.

| Option | Description | Default |
|--------|-------------|---------|
| `--seed` | Random seed | `42` |
| `--random-seed` | Legacy string seed for the dump utilities (hashed) | -- |
| `--reference-node` | Skip placement and use this node as the reference | -- |
| `--dump-sequence` | Write the FASTA for a single node and exit | -- |
| `--dump-sequences` | Write FASTAs for a list of node IDs to `<output>.dump-sequences.fa` | -- |
| `--simulate-snps` | SNPs to inject per `--dump-sequences` node (positional) | -- |
| `--dump-random-nodeIDs` | Write N random leaf node IDs to `<output>.randomNodeIDs.txt` | `0` |
| `--dump-all-scores` | Write all node placement scores to a TSV | -- |
| `--write-meta-read-scores-filtered` | Write filtered per-read meta scores to a TSV | off |
| `--write-meta-read-scores-unfiltered` | Write unfiltered per-read meta scores to a TSV | off |
| `--write-ocranks` | Write overlap-coefficient ranks to a TSV | off |
