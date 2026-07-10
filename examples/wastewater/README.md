# Estimating SARS-CoV-2 lineage abundances from wastewater samples

Steps used to process the Point Loma, San Diego wastewater samples and produce the manuscript results.

Requires minimap2, samtools, ivar, and pangolin.

Fetch reads and preprocess them:

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

Run panmap in metagenomic mode. It derives and caches the metagenomic index (`sars_20000_twilight_dipper.panman.midx`) next to the PanMAN on first run, so there is no separate build step:

```bash
../build/bin/panmap ../examples/data/sars_20000_twilight_dipper.panman SRR19707934.trimmed.fastq --meta --amplicon-depth SRR19707934.amplicon_stacks.tsv --mask-reads-relative-frequency 0.01 --em-delta-threshold 0.00001 --output SRR19707934 --threads 8
```

Panmap writes per-node abundance estimates to `SRR19707934.mgsr.abundance.out`.

For lineage proportions, retrieve the sequences from the PanMAN and identify their lineages with pangolin:

```bash
../build/bin/panmap ../examples/data/sars_20000_twilight_dipper.panman --dump-sequences "$(cut -f1 SRR19707934.mgsr.abundance.out | tr ',' '\n' | grep '\S' | tr '\n' ' ' | sed 's/ $//g')" --output SRR19707934
pangolin SRR19707934.dump-sequences.fa --outfile SRR19707934.pangolin.csv
python3 ../scripts/assign_lineage.py SRR19707934.mgsr.abundance.out SRR19707934.pangolin.csv > SRR19707934.mgsr.lineage.abundance.out
```