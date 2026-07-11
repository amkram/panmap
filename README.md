# Panmap

[![CI](https://github.com/amkram/panmap/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/amkram/panmap/actions/workflows/ci.yml)
[![bioRxiv](https://img.shields.io/badge/bioRxiv-10.64898%2F2026.03.29.711974-b31b1b)](https://www.biorxiv.org/content/10.64898/2026.03.29.711974v1)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/panmap?label=bioconda)](https://anaconda.org/bioconda/panmap)
[![Conda downloads](https://img.shields.io/conda/dn/bioconda/panmap?label=conda%20downloads)](https://anaconda.org/bioconda/panmap)
[![biocontainer](https://img.shields.io/badge/biocontainer-quay.io-blue)](https://quay.io/repository/biocontainers/panmap?tab=tags)
[![License: MIT](https://img.shields.io/github/license/amkram/panmap)](LICENSE)

Pangenome-based sequence placement, alignment, and genotyping.

[Documentation](https://amkram.github.io/panmap/) · [Preprint](https://www.biorxiv.org/content/10.64898/2026.03.29.711974v1)

Given reads and a [PanMAN](https://github.com/TurakhiaLab/panman) pangenome, panmap places each sample on the pangenome, then either aligns and genotypes it against the closest reference (single-sample mode) or estimates haplotype abundances and assigns reads to pangenome nodes (`--meta`).

## Install

With conda (~2 min):

```bash
conda install -c conda-forge -c bioconda panmap
panmap -h
```

With Docker (pulled from registry):

```bash
docker pull quay.io/biocontainers/panmap:0.1.2--0
docker run --rm quay.io/biocontainers/panmap:0.1.2--0 panmap -h
```

Or build locally with Docker (~6 min):

```bash
docker build -t panmap .
docker run --rm panmap panmap -h
```

Or build from source with conda:

```bash
conda env create -f environment.yml && conda activate panmap
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j
build/bin/panmap -h
```

## System requirements

Any modern Linux/macOS machine with an x86-64 CPU or Apple Silicon and ≥8 GB RAM.

Developed and benchmarked on a Linux server (dual Intel Xeon Gold 6338, 128 threads, Ubuntu 22.04 LTS), and also tested on macOS 14.3 (MacBook Air, Apple M3, 8 cores, 16 GB RAM, arm64).


## Single-sample mode (default)

**Usage**:

```bash
panmap [options] <pangenome.panman> <reads1.fastq> [reads2.fastq]
```

**Example**:

```bash
panmap examples/data/panmans/sars_20000_twilight_dipper.panman \
       examples/data/reads/isolate_R1.fastq.gz examples/data/reads/isolate_R2.fastq.gz
```
This places, genotypes, and assembles a small set of SARS-CoV-2 reads against a 20,000 genome PanMAN (~0.6s).

Outputs (prefix defaults to the reads filename, `isolate_R1`):

| File | Contents |
|---|---|
| `.placement.tsv` | best-match node per placement metric |
| `.ref.fa` | closest reference genome |
| `.bam`, `.bam.bai` | reads aligned to the closest reference |
| `.vcf` | variants vs the closest reference |
| `.consensus.fa` | sample consensus |


By default, panmap runs the full pipeline (`index → place → align → genotype → consensus`). Use `--stop <stage>` to run fewer stages (single-sample only; ignored with `--meta`).


## Metagenomic mode (`--meta`)

**Estimate haplotype abundances**:

```bash
panmap --meta [options] <pangenome.panman> <reads...>
```

**Example** (~2 min):

```bash
mkdir -p examples/output
panmap --meta -t 4 --em-delta-threshold 0.00001 -o examples/output/example \
       examples/data/panmans/sars_20000_twilight_dipper.panman \
       examples/data/reads/sars20000_5hap_*.fastq.gz
```

Writes `examples/output/example.mgsr.abundance.out` (haplotype abundance per sample). These are simulated shotgun-sequencing reads of SARS-CoV-2 mixtures. For wastewater samples see the [full documentation](https://amkram.github.io/panmap/) or [examples/wastewater](examples/wastewater).

**Filter and assign reads**:

```bash
panmap --meta --filter-and-assign [options] <pangenome.panman> <reads...>
```

**Example** (~2.5 min):

```bash
mkdir -p examples/output
panmap examples/data/panmans/v_mtdna.panman examples/data/reads/subsampled.fastq.gz \
       --meta --filter-and-assign -k 15 -s 8 -l 1 --discard 0.6 --dust 5 \
       --breadth-ratio -t 4 --taxonomic-metadata examples/data/metadata/v_mtdna.meta.tsv -o examples/output/subsampled
```

Writes (prefix `examples/output/subsampled`):

- `.mgsr.assignedReads.fastq`: the reads that were assigned
- `.mgsr.assignedReads.out`: per node, number of reads assigned and their indices into the fastq
- `.mgsr.assignedReadsLCANode.out`: per LCA node, number of reads assigned and their indices; a read's LCA node is the LCA of all nodes it was assigned to
- `.mgsr.breadths.out`: observed/expected breadth ratio per node (from `--breadth-ratio`)

## Verify

Reference outputs for the demos are in [examples/expected/](examples/expected/).

## Links

[Full documentation](https://amkram.github.io/panmap/) · [Installation options](https://amkram.github.io/panmap/installation/) · [CLI reference](https://amkram.github.io/panmap/cli-reference/) · [PanMAN format](https://github.com/TurakhiaLab/panman)
