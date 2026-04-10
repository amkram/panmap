# panmap

**Pangenome-based sequence placement, alignment, and genotyping.**

panmap takes sequencing reads and a pangenome in [PanMAN](https://github.com/TurakhiaLab/panman) format, places the reads onto the pangenome tree, aligns them to the closest reference, and calls variants.

## At a glance

```bash
# Install
conda install -c bioconda panmap

# Place and genotype paired-end reads
panmap ref.panman reads_R1.fq reads_R2.fq --stop genotype -t 8 -o sample

# Metagenomic abundance estimation
panmap ref.panman reads.fq --meta --index ref.idx -t 8 -o sample
```

## Modes

**Single-sample** (default)
:   Places reads from a single sample, aligns to the best-matching reference, and genotypes variants. Outputs BAM and VCF files.

**Metagenomic** (`--meta`)
:   Scores reads against every node in the PanMAN to estimate haplotype abundance or assign reads directly to nodes.

## Documentation

- [Installation](installation.md) -- Bioconda, Docker, building from source
- [Quick Start](quickstart.md) -- Pipeline overview and basic examples
- [Single-Sample Mode](single-sample.md) -- Genotyping walkthrough
- [Metagenomic Mode](metagenomic.md) -- Wastewater and aeDNA workflows
- [CLI Reference](cli-reference.md) -- All options and flags
