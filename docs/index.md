# panmap

**Pangenome-based sequence placement, alignment, and genotyping.**

panmap takes sequencing reads and a pangenome in [PanMAN](https://github.com/TurakhiaLab/panman) format, places the reads onto the pangenome tree, aligns them to the closest reference, and calls variants.

## Modes

**Single-sample** (default)
:   Places reads from a single sample, aligns to the best-matching reference, and genotypes variants. Outputs BAM and VCF files.

**Metagenomic** (`--meta`)
:   Scores reads against every node in the PanMAN to estimate haplotype abundance or assign reads directly to nodes.

## At a glance

```bash
# Install
conda install -c bioconda panmap

# Place reads (stops after placement by default)
panmap ref.panman reads_R1.fq reads_R2.fq -t 8 -o sample

# Run full pipeline through genotyping
panmap ref.panman reads_R1.fq reads_R2.fq --stop genotype -t 8 -o sample
```

## Documentation

- [Installation](installation.md) -- Docker and building from source
- [Quick Start](quickstart.md) -- First analysis in minutes
- [Single-Sample Mode](single-sample.md) -- Default pipeline walkthrough
- [Metagenomic Mode](metagenomic.md) -- Abundance estimation and read assignment
- [CLI Reference](cli-reference.md) -- All options and flags
