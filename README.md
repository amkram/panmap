# panmap

Pangenome-based sequence placement, alignment, and genotyping.

Given a pangenome (in [PanMAN](https://github.com/TurakhiaLab/panman) format) and sequencing reads, panmap places reads onto the pangenome tree, aligns them to the closest reference, and calls variants.

### Modes

- **Single-sample** (default): Places reads from a single sample, aligns to the best-matching reference, and genotypes variants (BAM + VCF output).
- **Metagenomic** (`--meta`): Scores reads from a mixture sample against every node in the PanMAN, and uses the scoring informtion to estimate haplotype abundance or directly assign reads to nodes.



### Installation with Docker (recommended)

```bash
docker build -t panmap .
```

See below for building from source.

## Usage

```
panmap <panman> [reads1.fq] [reads2.fq] [options]
```

### Pipeline stages

By default, panmap runs through genotyping. Use `--stop` to control how far the pipeline runs:

| Stage      | Output                  |
|------------|-------------------------|
| `index`    | `.idx` (seed index)     |
| `place`    | `.placement.tsv`        |
| `align`    | `.bam`                  |
| `genotype` | `.vcf`                  |

### Key options

```
-o, --output <prefix>     Output file prefix
-t, --threads <N>         Number of threads (default: 1)
-a, --aligner <str>       minimap2 (default) or bwa
--stop <stage>            Stop after: index, place, align, genotype
--meta                    Metagenomic mode
-k, --kmer <19>           Syncmer k
-s, --syncmer <8>         Syncmer s
--refine                  Alignment-based refinement of top candidates
--force-leaf              Restrict placement to leaf nodes
-v, --verbose             Verbose output
-q, --quiet               Errors only
```

Run `panmap --help` for the full option list.

### Example usage

**Place and genotype paired-end reads:**
```bash
panmap ref.panman reads_R1.fq reads_R2.fq --stop genotype -t 8 -o sample
```

**Metagenomic mode, estimating SARS-CoV-2 lineage abundances:**

For this, you need to have installed minimap2, samtools, ivar, and pangolin.

First step is to fetch some reads and preprocess them:

```bash
mkdir example_run && cd example_run
prefetch SRR19707934 && fasterq-dump SRR19707934
minimap2 -a --sam-hit-only --MD -2 -x sr ../examples/data/NC_045512v2.fa SRR19707934_1.fastq SRR19707934_2.fastq > SRR19707934.sam
samtools sort -o SRR19707934.sorted.bam SRR19707934.sam -@ 8
ivar trim -e -b ../examples/data/SNAPaddtlCov.bed -p SRR19707934.trimmed.bam -i SRR19707934.sorted.bam -q 1 -m 80 -x 3
samtools sort -o SRR19707934.trimmed.sorted.bam SRR19707934.trimmed.bam -@ 8
python3 ../scripts/trim_and_stack_amplicon.py stack SRR19707934.sorted.bam ../examples/data/SNAPaddtlCov.bed -o SRR19707934.amplicon_stacks.tsv
samtools fastq --no-sc SRR19707934.trimmed.sorted.bam > SRR19707934.trimmed.fastq
```

Then build an index and run panmap:

```bash
../build/bin/panmap ../examples/data/sars_20000_twilight_dipper.panman --index-mgsr sars_20000_twilight_dipper.idx 
../build/bin/panmap  ../examples/data/sars_20000_twilight_dipper.panman  SRR19707934.trimmed.fastq --meta --index sars_20000_twilight_dipper.idx --amplicon-depth SRR19707934.amplicon_stacks.tsv --mask-reads-relative-frequency 0.01 em-delta-threshold  0.00001 --output SRR19707934 --threads 8
```

Panmap outputs an abundance file `SRR19707934.mgsr.abundance.out` containing the estimated abundances of nodes.

To get the lineage proportions, we can retrieve the sequences from the PanMAN and use pangolin to identify their lineages:

```bash
../build/bin/panmap ../examples/data/sars_20000_twilight_dipper.panman --dump-sequences "$(cut -f1 SRR19707934.mgsr.abundance.out | tr ',' '\n' | grep '\S' | tr '\n' ' ' | sed 's/ $//g')" --output SRR19707934
pangolin pangolin SRR19707934.dump-sequences.fa --outfile SRR19707934.pangolin.csv
python3 ../scripts/assign_lineage.py SRR19707934.mgsr.abundance.out SRR19707934.pangolin.csv > SRR19707934.mgsr.lineage.abundance.out
```

**Metagenomic mode, filter and assign reads:**

We first build an index for the vertebrate mitochondrial PanMAN. We recommend using the `-k 15 -s 8 -l 1` seed parameters for aeDNA reads.

```bash
mkdir example_run && cd example_run
../build/bin/panmap ../examples/data/v_mtdna.panman --index-mgsr v_mtdna.idx -k 15 -s 8 -l 1
```

Then we run panmap with the `--filter-and-assign` option to filter and assign reads to the PanMAN.

```bash
../build/bin/panmap ../examples/data/v_mtdna.panman ../examples/data/subsampled.fastq.gz --meta -i v_mtdna.idx --filter-and-assign --discard 0.6 --dust 5 --taxonomic-metadata ../examples/data/v_mtdna.meta.tsv -t 4 --breadth-ratio --output subsampled
```

This outputs 3 files:

`.mgsr.assignedReads.fastq` file containing the reads that were assigned

`.mgsr.assignedReads.out` file containing the number of reads assigned to each node and the indices of the reads assigned, with respect to the the `.mgsr.assignedReads.fastq` file

`.mgsr.assignedReadsLCANode.out` file containg the number of reads assigned to the LCA node and the indices of the reads assigned. *As reads may be assigned to multiple nodes, the LCA node of a read if the LCA of all the nodes it was assigned to.*

### Advanced

#### Building from source:
Dependencies:

- `CMake >= 3.12`
- `C++17 compiler`
- `Protobuf (protobuf-compiler, libprotobuf-dev)`
- `Boost (program_options, iostreams, filesystem, system, date_time)`
- `zlib`
- `Cap'n Proto (capnproto, libcapnp-dev)`
- `Eigen3 (libeigen3-dev)`
- `htslib 1.20`
