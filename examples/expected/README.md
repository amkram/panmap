# Expected outputs for the README demos

Reference results for the three quick-start demos. Reproduce and diff a fresh run:

```bash
examples/check_examples.sh
```

| Directory | Demo | Files |
|-----------|------|-------|
| `single_sample/` | default pipeline | `isolate.placement.tsv`, `isolate.ref.fa`, `isolate.vcf`, `isolate.consensus.fa` |
| `meta_abundance/` | `--meta` abundance | `example.mgsr.abundance.out` |
| `filter_assign/` | `--filter-and-assign` | `subsampled.mgsr.assignedReads.fastq`, `.assignedReads.out`, `.assignedReadsLCANode.out` |

BAM and index files are not stored.
