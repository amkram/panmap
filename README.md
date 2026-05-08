# panmap

Pangenome-based sequence placement, alignment, and genotyping.

[Documentation](https://amkram.github.io/panmap/) | [Preprint](https://www.biorxiv.org/content/10.64898/2026.03.29.711974v1)

## System requirements

Panmap was developed and benchmarked on a Linux compute server (dual Intel Xeon Gold 6338, 128 threads, Ubuntu 22.04 LTS)
and also tested on macOS 14.3 (MacBook Air, Apple M3, 8 cores, 16 GB RAM, arm64).

For typical use, Panmap runs on any Linux machine with a modern x86-64 CPU and at least 8 GB of RAM. macOS
(Apple Silicon) is also supported through conda and Docker. No non-standard hardware is required.

## Installation

Panmap can be installed with conda, which takes about 2 minutes

```bash
conda install -c conda-forge -c bioconda panmap
panmap -h
```


Or built with Docker, which takes about 6 minutes:

```bash
docker build -t panmap .
docker run --rm panmap panmap -h
```

Or pulled with Docker:

```bash
docker pull quay.io/biocontainers/panmap:0.1.1--0
docker run --rm panmap panmap -h
```

## Quick start (demo)

### Single-sample Pipeline

```
index  -->  place  -->  align  -->  genotype  -->  consensus
 .idx    .placement.tsv  .bam       .vcf        .consensus.fa
```

By default, panmap runs the full pipeline. Use `--stop <stage>` to run fewer stages.

**Example:**

```bash
panmap examples/data/sars_20000_twilight_dipper.panman examples/data/isolate_R1.fastq.gz examples/data/isolate_R2.fastq.gz
```

This places, genotypes, and assembles a small set of SARS-CoV-2 reads against a 20,000 genome PanMAN (should complete in ~0.8s).

**Outputs**:
- Placement:  `isolate_R1.placement.tsv`
- Closest reference:  `isolate_R1.ref.fa`
- Read alignments to closest reference:  `isolate_R1.bam` + `isolate_R1.bam.bai`
- Variants from closest reference:   `isolate_R1.vcf`
- Consensus:  `isolate_R1.consensus.fa`

### For metagenomic samples

**Esimating haplotype abundance**

This demo run should take about 2 minutes to complete.

```bash
data_dir=examples/data; output_dir=examples/output; mkdir -p $output_dir

# Build an index for metagenomics mode
panmap $data_dir/sars_20000_twilight_dipper.panman --index-mgsr $output_dir/sars_20000_twilight_dipper.idx 

# Run panmap with --meta option
panmap  $data_dir/sars_20000_twilight_dipper.panman $data_dir/*.fastq.gz --meta --index $output_dir/sars_20000_twilight_dipper.idx --threads 4 --em-delta-threshold 0.00001 --output $output_dir/example
```

This outputs a `.mgsr.abundance.out` file containing the haplotype abundance for each sample.

Reads used above were simulated shotgun-sequencing reads of SARS-CoV-2 mixtures. For wastewater samples, refer to the
full documentation linked below or README in examples/wastewater for more details.

**Filter and assign reads**

This demo run should take about 2.5 minutes to complete.

```bash
data_dir=examples/data; output_dir=examples/output; mkdir -p $output_dir
cp data/vertebrate_mtdna/v_mtdna.panman $data_dir/

# Build an index for metagenomics mode with ancient dna parameters
panmap $data_dir/v_mtdna.panman --index-mgsr $data_dir/v_mtdna.idx -k 15 -s 8 -l 1

# Run Panmap with --filter-and-assign option
panmap $data_dir/v_mtdna.panman $data_dir/subsampled.fastq.gz --meta -i $data_dir/v_mtdna.idx --filter-and-assign \
  --discard 0.6 --dust 5 --taxonomic-metadata $data_dir/v_mtdna.meta.tsv -t 4 --breadth-ratio --output $output_dir/subsampled
```

This outputs 3 files:

- `.mgsr.assignedReads.fastq` file containing the reads that were assigned

- `.mgsr.assignedReads.out` file containing the number of reads assigned to each node and the indices of the reads assigned, with respect to the the .mgsr.assignedReads.fastq file

- `.mgsr.assignedReadsLCANode.out` file containg the number of reads assigned to the LCA node and the indices of the reads assigned. As reads may be assigned to multiple nodes, the LCA node of a read if the LCA of all the nodes it was assigned to

## Modes

- **Single-sample** (default): Place reads, align to closest reference, call variants (BAM + VCF)
- **Metagenomic** (`--meta`): Estimate haplotype abundance or assign reads to pangenome nodes

## Links

- [Full documentation](https://amkram.github.io/panmap/)
- [Installation options](https://amkram.github.io/panmap/installation/)
- [CLI reference](https://amkram.github.io/panmap/cli-reference/)
- [PanMAN format](https://github.com/TurakhiaLab/panman)
