# Expected outputs for the README demos

Reference results for the three quick-start demos. Reproduce and diff a fresh run:

```bash
examples/check_examples.sh
```

```
examples/expected/
├── single_sample/                              default pipeline: place → align → genotype → consensus
│   ├── isolate.placement.tsv                   best-match node per placement metric
│   ├── isolate.ref.fa                          placement node genome (alignment reference)
│   ├── isolate.vcf                             variants called vs the reference
│   └── isolate.consensus.fa                    sample consensus (reference + variants)
├── meta_abundance/                             --meta: haplotype abundance
│   └── example.mgsr.abundance.out              node → estimated fraction
└── filter_assign/                              --filter-and-assign: filter reads, assign to nodes
    ├── subsampled.mgsr.assignedReads.fastq         reads passing the filter
    ├── subsampled.mgsr.assignedReads.out           node → read count + read indices
    └── subsampled.mgsr.assignedReadsLCANode.out    LCA node → read count + read indices
```

BAM and index files are not stored.
