# Single-Sample Mode

The default mode. Places reads from one sample onto the pangenome tree, aligns to the best-matching reference, and genotypes variants.

## Basic usage

```bash
# Full pipeline (default: through consensus)
panmap ref.panman reads_R1.fq reads_R2.fq -t 8 -o sample

# Stop early
panmap ref.panman reads_R1.fq reads_R2.fq --stop genotype -t 8 -o sample
```

## Pipeline stages

| Stage | `--stop` value | Output | Description |
|-------|---------------|--------|-------------|
| Index | `index` | `<panman>.idx` | Builds seed index from the PanMAN (next to the PanMAN) |
| Place | `place` | `<prefix>.placement.tsv` | Places reads onto the pangenome tree |
| Align | `align` | `<prefix>.bam` | Aligns reads to closest reference |
| Genotype | `genotype` | `<prefix>.vcf` | Calls variants from alignments |
| Consensus | `consensus` (default) | `<prefix>.consensus.fa` | Generates consensus FASTA from variants |

!!! note
    `--force-leaf` is enabled automatically whenever the pipeline runs past
    placement (the default run, `--stop align`, `--stop genotype`, etc.). Only
    `--stop place` allows reads to be placed on internal nodes.

## Example: SARS-CoV-2 genotyping from Illumina reads

This example places paired-end reads onto a SARS-CoV-2 pangenome and calls variants.

### 1. Get a PanMAN

Download or build a PanMAN for your organism. For SARS-CoV-2, a pre-built PanMAN with 20,000 samples is included:

```bash
ls examples/data/sars_20000_twilight_dipper.panman
```

### 2. Run the full pipeline

```bash
panmap examples/data/sars_20000_twilight_dipper.panman \
  reads_R1.fq.gz reads_R2.fq.gz \
  --stop genotype -t 8 -o my_sample
```

### 3. Output files

| File | Contents |
|------|----------|
| `my_sample.placement.tsv` | Placement on the pangenome tree |
| `my_sample.bam` | Reads aligned to the closest reference |
| `my_sample.vcf` | Called variants |
| `my_sample.consensus.fa` | Consensus sequence |

The seed index is built once next to the PanMAN (`sars_20000_twilight_dipper.panman.idx`) and reused automatically on later runs.

### 4. Reuse the index

Use `--reindex` to force a rebuild, or `--index <path>` to load a pre-built index from another location:

```bash
panmap examples/data/sars_20000_twilight_dipper.panman \
  sample2_R1.fq.gz sample2_R2.fq.gz \
  --stop consensus -t 8 -o sample2
```

## Common options

| Option | Description | Default |
|--------|-------------|---------|
| `-o, --output` | Output file prefix | derived from reads filename |
| `-t, --threads` | Number of threads | `1` |
| `-a, --aligner` | `minimap2` or `bwa` | `minimap2` |
| `--stop` | Stop after: `index`, `place`, `align`, `genotype`, `consensus` | `consensus` |
| `--force-leaf` | Restrict placement to leaf nodes | off (auto-enabled unless `--stop place`) |
| `--refine` | Alignment-based refinement of top candidates | off |

## Choosing an aligner

- **minimap2** (default): fast, good for most cases
- **bwa**: higher sensitivity, may be preferable for short reads

```bash
panmap ref.panman reads.fq -a bwa -o sample
```

## Refining placement

`--refine` uses alignment scores to re-rank placement candidates. Slower but can improve accuracy when seed-based placement is ambiguous.

```bash
panmap ref.panman reads.fq --refine -o sample
```

Refinement parameters:

| Option | Description | Default |
|--------|-------------|---------|
| `--refine-top-pct` | Top fraction of nodes to refine | `0.01` (1%) |
| `--refine-max-top-n` | Max nodes to align against | `150` |
| `--refine-neighbor-radius` | Expand to neighbors within N branches | `2` |
| `--refine-max-neighbor-n` | Max additional nodes from neighbor expansion | `150` |

## Advanced options

| Option | Description | Default |
|--------|-------------|---------|
| `-i, --index` | Path to pre-built index file | derived from output/panman |
| `-f, --reindex` | Force rebuild index | off |
| `--dedup` | Deduplicate reads | off |
| `--impute` | Impute N's from parent (skip _->N mutations) | off |
| `--no-mutation-spectrum` | Disable mutation spectrum filtering in VCF genotyping | off |
| `--baq` | Enable Base Alignment Quality in mpileup | off |
| `--batch` | Batch file listing samples (one per line: `reads1 [reads2] [output_prefix]`) | -- |

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
| `--trim-start` | Trim bases from read start | `0` |
| `--trim-end` | Trim bases from read end | `0` |
| `--extent-guard` | Guard seed deletions at genome extent boundaries | off |
