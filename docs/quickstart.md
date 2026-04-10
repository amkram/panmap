# Quick Start

## Synopsis

```
panmap <panman> [reads1.fq] [reads2.fq] [options]
```

## Pipeline

panmap runs five stages in sequence. By default it stops after **placement**. Use `--stop` to run further.

```
index  -->  place  -->  align  -->  genotype  -->  consensus
 .idx    .placement.tsv  .bam       .vcf        .consensus.fa
```

## Single-sample genotyping

Place reads onto the pangenome, align to the closest reference, call variants, and generate a consensus:

```bash
panmap ref.panman reads_R1.fq reads_R2.fq \
  --stop consensus -t 8 -o sample
```

This produces `sample.bam`, `sample.vcf`, and `sample.consensus.fa`.

## Metagenomic abundance estimation

Estimate which lineages are present in a mixed sample:

```bash
# Build metagenomic index (once per pangenome)
panmap ref.panman --index-mgsr ref.idx

# Estimate abundances
panmap ref.panman reads.fq \
  --meta --index ref.idx -t 8 -o sample
```

Output: `sample.mgsr.abundance.out`

## Partial pipelines

```bash
# Build index only
panmap ref.panman --stop index -o ref

# Place reads (default)
panmap ref.panman reads.fq -o sample

# Place and align, skip genotyping
panmap ref.panman reads.fq --stop align -o sample
```

## Next steps

- [Single-Sample Mode](single-sample.md) -- full walkthrough with examples
- [Metagenomic Mode](metagenomic.md) -- wastewater and aeDNA workflows
- [CLI Reference](cli-reference.md) -- all options
