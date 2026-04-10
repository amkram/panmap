# panmap

Pangenome-based sequence placement, alignment, and genotyping.

[Documentation](https://amkram.github.io/panmap/) | [Preprint](https://www.biorxiv.org/content/10.64898/2026.03.29.711974v1)

## Install

```bash
conda install -c bioconda panmap
```

Or with Docker:

```bash
docker pull alanalohaucsc/panmap:latest
```

## Quick start

```bash
# Place and genotype paired-end reads
panmap ref.panman reads_R1.fq reads_R2.fq --stop genotype -t 8 -o sample

# Metagenomic abundance estimation
panmap ref.panman reads.fq --meta --index ref.idx -t 8 -o sample
```

## Pipeline

```
index  -->  place  -->  align  -->  genotype
 .idx    .placement.tsv  .bam       .vcf
```

By default, panmap stops after placement. Use `--stop` to run further stages.

## Modes

- **Single-sample** (default): Place reads, align to closest reference, call variants (BAM + VCF)
- **Metagenomic** (`--meta`): Estimate haplotype abundance or assign reads to pangenome nodes

## Links

- [Full documentation](https://amkram.github.io/panmap/)
- [Installation options](https://amkram.github.io/panmap/installation/)
- [CLI reference](https://amkram.github.io/panmap/cli-reference/)
- [PanMAN format](https://github.com/TurakhiaLab/panman)
