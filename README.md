# panmap

Pangenome-based sequence placement, alignment, and genotyping.

[Documentation](https://amkram.github.io/panmap/) | [Preprint](https://www.biorxiv.org/content/10.64898/2026.03.29.711974v1)

## Install

```bash
conda install -c conda-forge -c bioconda panmap
panmap -h
```

Or build with Docker:

```bash
docker build -t panmap .
docker run --rm panmap panmap -h
```

Or pull with Docker:

```bash
docker pull quay.io/biocontainers/panmap:0.1.1--0
docker run --rm panmap panmap -h
```

## Quick start

### For isolated samples 

```bash
# Place and genotype paired-end reads
panmap ref.panman reads_R1.fq reads_R2.fq -t 8 -o sample
```

**Pipeline**

```
index  -->  place  -->  align  -->  genotype  -->  consensus
 .idx    .placement.tsv  .bam       .vcf        .consensus.fa
```

By default, panmap stops after placement. Use `--stop` to run further stages.

### For metagenomic samples

**Esimating haplotype abundance**

```bash
data_dir=examples/data
# Build an index for metagenomics mode
panmap $data_dir/sars_20000_twilight_dipper.panman --index-mgsr $data_dir/sars_20000_twilight_dipper.idx 

# Run panmap with --meta option
panmap  $data_dir/sars_20000_twilight_dipper.panman $data_dir/*.fastq.gz --meta --index $data_dir/sars_20000_twilight_dipper.idx --threads 8 --em-delta-threshold 0.00001
```

Reads used above were simulated shotgun-sequencing reads of SARS-CoV-2 mixtures. For wastewater samples, refer to the
full documentation linked below or README in examples/wastewater for more details.

**Filter and assign reads**

```bash
data_dir=examples/data
cp data/vertebrate_mtdna/v_mtdna.panman $data_dir/

# Build an index for metagenomics mode with ancient dna parameters
panmap $data_dir/v_mtdna.panman --index-mgsr $data_dir/v_mtdna.idx -k 15 -s 8 -l 1

# Run Panmap with --filter-and-assign option
panmap $data_dir/v_mtdna.panman $data_dir/subsampled.fastq.gz --meta -i $data_dir/v_mtdna.idx --filter-and-assign \
  --discard 0.6 --dust 5 --taxonomic-metadata $data_dir/v_mtdna.meta.tsv -t 4 --breadth-ratio --output $data_dir/subsampled

```

## Modes

- **Single-sample** (default): Place reads, align to closest reference, call variants (BAM + VCF)
- **Metagenomic** (`--meta`): Estimate haplotype abundance or assign reads to pangenome nodes

## Links

- [Full documentation](https://amkram.github.io/panmap/)
- [Installation options](https://amkram.github.io/panmap/installation/)
- [CLI reference](https://amkram.github.io/panmap/cli-reference/)
- [PanMAN format](https://github.com/TurakhiaLab/panman)
