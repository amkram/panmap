# panmap

Pangenome-based sequence placement, alignment, and genotyping.

Given a pangenome (in [PanMAN](https://github.com/TurakhiaLab/panman) format) and sequencing reads, panmap places reads onto the pangenome tree, aligns them to the closest reference, and calls variants.

### Modes

- **Single-sample** (default): Places reads from a single sample, aligns to the best-matching reference, and genotypes variants (BAM + VCF output).
- **Metagenomic** (`--meta`): Assigns reads from a metagenome to multiple references with EM-based abundance estimation.



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

### Examples

Place and genotype paired-end reads:
```bash
panmap ref.panman reads_R1.fq reads_R2.fq --stop genotype -t 8 -o sample
```

Metagenomic assignment:
```bash
panmap ref.panman sample.fq --meta --top-oc 1000 -t 24
```

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
