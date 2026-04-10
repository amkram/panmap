# CLI Reference

```
panmap <panman> [reads1.fq] [reads2.fq] [options]
```

Use `--help` for common options or `--help-all` for the full list.

## General

| Option | Description | Default |
|--------|-------------|---------|
| `-h, --help` | Show common options | -- |
| `--help-all` | Show all options | -- |
| `-V, --version` | Show version | -- |
| `-o, --output` | Output file prefix | derived from reads filename |
| `-t, --threads` | Number of threads | `1` |
| `-a, --aligner` | `minimap2` or `bwa` | `minimap2` |
| `--stop` | Stop after: `index`, `place`, `align`, `genotype`, `consensus` | `consensus` |
| `--meta` | Enable metagenomic mode | off |
| `-v, --verbose` | Verbose output | off |
| `-q, --quiet` | Errors only | off |
| `--no-color` | Disable colored output | off |

## Index

| Option | Description | Default |
|--------|-------------|---------|
| `-i, --index` | Path to pre-built index file | derived from output/panman |
| `-f, --reindex` | Force rebuild index | off |
| `--index-mgsr` | Build/rebuild MGSR index (metagenomic) | -- |
| `--index-full` | Build full index instead of lite | off |
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
| `--force-leaf` | Restrict placement to leaf nodes (auto-enabled with `--stop genotype`) | off |
| `--refine` | Alignment-based refinement of top candidates | off |
| `--refine-top-pct` | Top fraction of nodes to refine | `0.01` |
| `--refine-max-top-n` | Max nodes to align against | `150` |
| `--refine-neighbor-radius` | Expand to neighbors within N branches | `2` |
| `--refine-max-neighbor-n` | Max additional nodes from neighbor expansion | `150` |

## Genotyping

| Option | Description | Default |
|--------|-------------|---------|
| `--impute` | Impute N's from parent (skip _->N mutations) | off |
| `--no-mutation-spectrum` | Disable mutation spectrum filtering in VCF | off |
| `--baq` | Enable Base Alignment Quality in mpileup | off |

## Metagenomic: EM algorithm

| Option | Description | Default |
|--------|-------------|---------|
| `--top-oc` | Top N nodes by overlap coefficient to send to EM | `1000` |
| `--mask-reads` | Mask reads with k-min-mer total occurrence <= threshold | `0` |
| `--mask-seeds` | Mask k-min-mer seeds with total occurrence <= threshold | `0` |
| `--amplicon-depth` | Amplicon depth TSV for frequency-based masking | -- |
| `--mask-reads-relative-frequency` | Mask reads with relative frequency < threshold * amplicon_depth | `0.0` |
| `--mask-seeds-relative-frequency` | Mask seeds with relative frequency < threshold * amplicon_depth | `0.0` |
| `--em-convergence-threshold` | Converge when likelihood difference < threshold | `0.00001` |
| `--em-delta-threshold` | Converge when max proportion change < threshold | `0.0` |
| `--em-maximum-rounds` | Maximum EM rounds | `5` |
| `--em-maximum-iterations` | Maximum EM iterations per round | `1000` |
| `--em-leaves-only` | Only run EM on leaf (sample) nodes | off |
| `--no-progress` | Disable progress bars | off |

## Metagenomic: filter and assign

| Option | Description | Default |
|--------|-------------|---------|
| `--filter-and-assign` | Assign reads to nodes without EM | off |
| `--dust` | Discard reads with DUST score > threshold | `100.0` (no filtering) |
| `--discard` | Discard reads with parsimony score < threshold * total seeds | `0.0` (no discard) |
| `--mask-read-ends` | Mask N bases from read ends (for aeDNA damage) | `0` |
| `--taxonomic-metadata` | TSV with taxonomic metadata per node | -- |
| `--maximum-families` | Discard reads spanning > N taxonomic families | `1` |
| `--ambiguous-score-threshold-ratio` | Discard reads scoring outside max families by ratio | `0.0` |
| `--ambiguous-score-threshold` | Discard reads scoring outside max families by absolute value | `0` |
| `--breadth-ratio` | Calculate observed/expected breadth ratio | off |
| `--pseudochain` | Use pseudo-chains for read scoring | off |
| `--batch-files-path` | TSV file containing batch file paths | -- |
| `--batch-size` | Batch size for filtering and assigning | `1000000` |

## Batch mode

| Option | Description | Default |
|--------|-------------|---------|
| `--batch` | File listing samples, one per line: `reads1 [reads2] [output_prefix]` | -- |

## Utility

| Option | Description |
|--------|-------------|
| `--dump-sequence <node>` | Extract FASTA for a single node |
| `--dump-sequences <nodes...>` | Extract FASTA for multiple node IDs |
