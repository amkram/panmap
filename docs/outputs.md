# Outputs

Every file panmap writes is prefixed with the output prefix (`-o`, or the reads
filename by default). Indexes are the exception: they are written next to the
PanMAN so they can be shared across samples.

## Single-sample mode

Running the full pipeline (`panmap ref.panman reads.fq -o sample`) produces:

| File | Contents |
|------|----------|
| `sample.placement.tsv` | Best-matching node for each placement metric |
| `sample.ref.fa` | Closest reference genome (with `.ref.fa.fai` index) |
| `sample.bam` | Reads aligned to the closest reference (with `.bam.bai` index) |
| `sample.vcf` | Variants called against the closest reference |
| `sample.consensus.fa` | Sample consensus sequence |

Stopping early with `--stop <stage>` produces only the files up to that stage.

### `.placement.tsv`

One row per placement metric, giving the metric's best-scoring node:

```
metric                score       nodes
log_raw               74.238025   England/CAMC-141DBEA/2021|OU190639.1|2021-03-19
log_cosine            0.787143    England/CAMC-141DBEA/2021|OU190639.1|2021-03-19
containment           0.196456    England/CAMC-141DBEA/2021|OU190639.1|2021-03-19
weighted_containment  1.065240    England/CAMC-141DBEA/2021|OU190639.1|2021-03-19
log_containment       0.454292    England/CAMC-141DBEA/2021|OU190639.1|2021-03-19
```

The metrics score how well the sample's seeds match each node's seed set:

| Metric | Meaning |
|--------|---------|
| `log_raw` | Log of the raw number of shared seeds |
| `log_cosine` | Log cosine similarity of the seed-count vectors |
| `containment` | Fraction of the sample's seeds contained in the node |
| `weighted_containment` | Containment weighted by seed multiplicity |
| `log_containment` | Log-scaled containment (the default placement score) |

The reference for alignment and genotyping is the node chosen by
`log_containment`.

### `.vcf`, `.bam`, `.consensus.fa`

Standard formats: the BAM holds reads aligned to `sample.ref.fa`, the VCF holds
variants relative to that reference, and the consensus FASTA applies those
variants to the reference. All are viewable in IGV, `samtools`, `bcftools`, etc.

## Metagenomic mode

### Abundance (`--meta`)

| File | Contents |
|------|----------|
| `sample.mgsr.abundance.out` | Estimated proportion per haplotype |

Tab-separated `node_id` and proportion, sorted by proportion (proportions sum to
~1.0):

```
USA/MA-CDCBI-.../2021|OK176722.1|2021-09-07     0.50052
USA/SEARCH-6918/2021|PP046785.1|2021-01-21      0.20476
Germany/IMS-10036-.../2021|OU071802.1|2021-02-03 0.15086
```

### Filter and assign (`--meta --filter-and-assign`)

| File | Contents |
|------|----------|
| `sample.mgsr.assignedReads.fastq` | The reads that were assigned |
| `sample.mgsr.assignedReads.out` | Per node: number of reads assigned and their indices into the FASTQ |
| `sample.mgsr.assignedReadsLCANode.out` | Per LCA node: assigned-read counts and indices (a read's LCA node is the LCA of all nodes it was assigned to) |

## Indexes

Indexes are written next to the PanMAN, not under the output prefix, and are
reused automatically on later runs.

| File | Built by | Purpose |
|------|----------|---------|
| `<panman>.idx` | default / `--stop index` | Placement index (single-sample) |
| `<panman>.midx` | `--meta` | Metagenomic (MGSR) index |

Use `--index-out <path>` to write to a custom location and `--index <path>` to
load one. `--reindex` forces a rebuild. See the
[Quick Start](quickstart.md#managing-indexes) for details.
