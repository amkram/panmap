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
# The metagenomic index (.midx) is built automatically on first run
panmap ref.panman reads.fq --meta -t 8 -o sample
```

Output: `sample.mgsr.abundance.out`

## Partial pipelines

```bash
# Place reads (default; the index is built automatically if absent)
panmap ref.panman reads.fq -o sample

# Place and align, skip genotyping
panmap ref.panman reads.fq --stop align -o sample
```

## Managing indexes

panmap builds the index automatically on first run: a placement index (`.idx`) in
the default pipeline and an MGSR index (`.midx`) under `--meta`. To build one
explicitly or reuse a prebuilt index:

```bash
# Build the index only (no reads)
panmap ref.panman --stop index -o ref          # -> ref.idx  (placement)
panmap ref.panman --meta --stop index -o ref   # -> ref.midx (metagenomic)

# Reuse a prebuilt index with --index
panmap ref.panman reads.fq --index ref.idx -o sample
panmap ref.panman reads.fq --meta --index ref.midx -o sample

# Force a rebuild
panmap ref.panman reads.fq --reindex -o sample
```

## Next steps

- [Single-Sample Mode](single-sample.md) -- full walkthrough with examples
- [Metagenomic Mode](metagenomic.md) -- wastewater and aeDNA workflows
- [CLI Reference](cli-reference.md) -- all options
