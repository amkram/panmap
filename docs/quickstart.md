# Quick Start

## Synopsis

```
panmap <panman> [reads1.fq] [reads2.fq] [options]
```

## Pipeline

panmap runs five stages in sequence. By default it runs through **consensus**. Use `--stop` to stop earlier.

```
index  -->  place  -->  align  -->  genotype  -->  consensus
 .idx    .placement.tsv  .bam       .vcf        .consensus.fa
```

## Single-sample genotyping

Place reads onto the pangenome, align to the closest reference, call variants, and generate a consensus:

```bash
panmap ref.panman reads_R1.fq reads_R2.fq -t 8 -o sample
```

This produces `sample.bam`, `sample.vcf`, and `sample.consensus.fa`.

## Metagenomic abundance estimation

Estimate which lineages are present in a mixed sample:

```bash
panmap ref.panman reads.fq --meta -t 8 -o sample
```

Output: `sample.mgsr.abundance.out`

## Partial pipelines

```bash
# Place reads (default)
panmap ref.panman reads.fq -o sample

# Place and align, skip genotyping
panmap ref.panman reads.fq --stop align -o sample
```

## Managing indexes

panmap builds the index automatically on first run and reuses it after. The path is
derived from the panman (not `-o`): a placement index at `<panman>.idx` and an MGSR
index at `<panman>.midx` under `--meta`.

```bash
# Build the index only, next to the panman (no reads needed)
panmap ref.panman --stop index          # -> ref.panman.idx   (placement)
panmap ref.panman --meta --stop index   # -> ref.panman.midx  (metagenomic)

# Build to a custom path with --index-out (output), then load it with --index (input)
panmap ref.panman --stop index --index-out /data/ref.idx
panmap ref.panman reads.fq --index /data/ref.idx -o sample

# Force a rebuild
panmap ref.panman reads.fq --reindex -o sample
```

## Next steps

- [Single-Sample Mode](single-sample.md) -- full walkthrough with examples
- [Metagenomic Mode](metagenomic.md) -- wastewater and aeDNA workflows
- [CLI Reference](cli-reference.md) -- all options
