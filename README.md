# panmap

Pangenome-based sequence placement, alignment, and genotyping.

Given a pangenome (in [PanMAN](https://github.com/TurakhiaLab/panman) format) and sequencing reads, panmap places reads onto the pangenome tree, aligns them to the closest reference, and calls variants.

### Modes

- **Single-sample** (default): Places reads from a single sample, aligns to the best-matching reference, and genotypes variants (BAM + VCF output).
- **Metagenomic** (`--meta`): Scores reads against every node in the PanMAN to estimate haplotype abundance or assign reads to nodes.

### Install

**From source:**
```bash
git clone https://github.com/amkram/panmap.git
cd panmap && mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

See the [installation docs](https://amkram.github.io/panmap/installation/) for dependencies and platform-specific notes.

**Docker:**
```bash
docker build -t panmap .
docker run --rm panmap --help
```

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

### Examples

**Place and genotype paired-end reads:**
```bash
panmap sars.panman reads_R1.fq reads_R2.fq -t 8 -o sample
# Output: sample.placement.tsv, sample.bam, sample.vcf
```

**Index only (reuse for multiple samples):**
```bash
panmap sars.panman --stop index -o sars
panmap sars.panman reads.fq -i sars.idx -o sample1
panmap sars.panman reads2.fq -i sars.idx -o sample2
```

**Metagenomic mode:**
```bash
panmap sars.panman --index-mgsr sars_mgsr.idx
panmap sars.panman reads_R1.fq reads_R2.fq \
  --meta -i sars_mgsr.idx -t 8
```

See the [documentation](https://amkram.github.io/panmap/) for detailed usage guides.

### Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for build instructions, testing, and code style.

### Citation

If you use panmap in your research, please cite:

> Kramer, A. et al. (2024). Pangenome-based sequence analysis with panmap. *bioRxiv*. [doi:10.1101/2024.XX.XX](https://doi.org/10.1101/2024.XX.XX)
