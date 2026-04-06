# Quick Start

## Synopsis

```
panmap <panman> [reads1.fq] [reads2.fq] [options]
```

## Pipeline

panmap runs four stages in sequence. By default, it stops after **placement**. Use `--stop` to run further.

```
index  ──>  place  ──>  align  ──>  genotype
 .idx    .placement.tsv  .bam       .vcf
```

## Example: paired-end genotyping

```bash
panmap ref.panman reads_R1.fq reads_R2.fq --stop genotype -t 8 -o sample
```

This runs the full pipeline and produces:

| File | Contents |
|------|----------|
| `sample.idx` | Seed index (reusable) |
| `sample.placement.tsv` | Tree placement |
| `sample.bam` | Aligned reads |
| `sample.vcf` | Called variants |

## Running partial pipelines

```bash
# Build index only
panmap ref.panman --stop index -o ref

# Place reads (default behavior)
panmap ref.panman reads.fq -o sample

# Place and align, but skip genotyping
panmap ref.panman reads.fq --stop align -o sample
```

## Next steps

- [Single-Sample Mode](single-sample.md) -- pipeline details and options
- [Metagenomic Mode](metagenomic.md) -- abundance estimation and read assignment
- [CLI Reference](cli-reference.md) -- full option list
