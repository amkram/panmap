# Single-Sample Mode

The default mode. Places reads from one sample onto the pangenome tree, aligns to the best-matching reference, and genotypes variants.

## Basic usage

```bash
# Run full pipeline (default -- through consensus)
panmap ref.panman reads_R1.fq reads_R2.fq -t 8 -o sample

# Stop at an earlier stage
panmap ref.panman reads_R1.fq reads_R2.fq --stop genotype -t 8 -o sample
```

## Pipeline stages

By default, panmap runs through **consensus**. Use `--stop` to stop at an earlier stage.

| Stage | `--stop` value | Output | Description |
|-------|---------------|--------|-------------|
| Index | `index` | `<prefix>.idx` | Builds seed index from the PanMAN |
| Place | `place` | `<prefix>.placement.tsv` | Places reads onto the pangenome tree |
| Align | `align` | `<prefix>.bam` | Aligns reads to closest reference |
| Genotype | `genotype` | `<prefix>.vcf` | Calls variants from alignments |
| Consensus | `consensus` (default) | `<prefix>.consensus.fa` | Generates consensus FASTA from variants |

!!! note
    When `--stop genotype` is used, `--force-leaf` is enabled automatically.

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
| `my_sample.idx` | Seed index (reusable for future runs with `--index`) |
| `my_sample.placement.tsv` | Placement on the pangenome tree |
| `my_sample.bam` | Reads aligned to the closest reference |
| `my_sample.vcf` | Called variants |
| `my_sample.consensus.fa` | Consensus sequence |

### 4. Reuse the index

Once built, the index can be reused across samples:

```bash
panmap examples/data/sars_20000_twilight_dipper.panman \
  sample2_R1.fq.gz sample2_R2.fq.gz \
  --index my_sample.idx --stop consensus -t 8 -o sample2
```

## Common options

| Option | Description | Default |
|--------|-------------|---------|
| `-o, --output` | Output file prefix | derived from reads filename |
| `-t, --threads` | Number of threads | `1` |
| `-a, --aligner` | `minimap2` or `bwa` | `minimap2` |
| `--stop` | Stop after: `index`, `place`, `align`, `genotype`, `consensus` | `consensus` |
| `--force-leaf` | Restrict placement to leaf nodes | off (auto-enabled with `--stop genotype`) |
| `--refine` | Alignment-based refinement of top candidates | off |

## Choosing an aligner

- **minimap2** (default) -- fast, good for most use cases
- **bwa** -- may be preferable for short reads or when higher sensitivity is needed

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
