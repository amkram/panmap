# panmap

**Pangenome-based sequence placement, alignment, and genotyping.**

[![Bioconda](https://img.shields.io/conda/vn/bioconda/panmap?label=bioconda)](https://anaconda.org/bioconda/panmap)
[![biocontainer](https://img.shields.io/badge/biocontainer-quay.io-blue)](https://quay.io/repository/biocontainers/panmap?tab=tags)
[![License: MIT](https://img.shields.io/github/license/amkram/panmap)](https://github.com/amkram/panmap/blob/main/LICENSE)
[![bioRxiv](https://img.shields.io/badge/bioRxiv-10.64898%2F2026.03.29.711974-b31b1b)](https://www.biorxiv.org/content/10.64898/2026.03.29.711974v1)

panmap takes sequencing reads and a pangenome in [PanMAN](https://github.com/TurakhiaLab/panman)
format, places each sample onto the pangenome, and then either aligns and
genotypes it against the closest reference (single-sample mode) or estimates
haplotype abundances and assigns reads to pangenome nodes (metagenomic mode).

## At a glance

```bash
# Install
conda install -c conda-forge -c bioconda panmap

# Place, genotype, and build a consensus (single-sample, default)
panmap ref.panman reads_R1.fq reads_R2.fq -t 8 -o sample

# Estimate haplotype abundances from a mixture (metagenomic)
panmap ref.panman reads.fq --meta -t 8 -o sample
```

## Modes

**Single-sample** (default)
:   Places reads from a single sample, aligns to the best-matching reference, and
    genotypes variants. Produces a placement, BAM, VCF, and consensus FASTA.

**Metagenomic** (`--meta`)
:   Scores reads against every node in the PanMAN to estimate haplotype abundance,
    or assigns reads directly to nodes for taxonomic identification.

## Documentation

<div class="grid cards" markdown>

- :material-download: **[Installation](installation.md)** — Bioconda, Docker, building from source
- :material-rocket-launch: **[Quick Start](quickstart.md)** — pipeline overview and basic examples
- :material-dna: **[Single-Sample Mode](single-sample.md)** — genotyping walkthrough
- :material-flask: **[Metagenomic Mode](metagenomic.md)** — wastewater and aeDNA workflows
- :material-file-document: **[Outputs](outputs.md)** — every file panmap writes
- :material-console: **[CLI Reference](cli-reference.md)** — all options and flags

</div>

## Citing panmap

If you use panmap in your work, please cite the preprint:

> panmap: pangenome-based sequence placement, alignment, and genotyping.
> *bioRxiv* [10.64898/2026.03.29.711974](https://www.biorxiv.org/content/10.64898/2026.03.29.711974v1).
