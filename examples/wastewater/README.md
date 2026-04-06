# Estimating SARS-CoV-2 lineage abundances from wastewater samples using Panmap

Procedures below were used to process wastewater samples collected from Point Loma, San Diego and generate the results
as shown in the manuscript


For this, you need to have installed minimap2, samtools, ivar, and pangolin.

First step is to fetch some reads and preprocess them:

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

Then build an index and run panmap:

```bash
../build/bin/panmap ../examples/data/sars_20000_twilight_dipper.panman --index-mgsr sars_20000_twilight_dipper.idx 
../build/bin/panmap  ../examples/data/sars_20000_twilight_dipper.panman  SRR19707934.trimmed.fastq --meta --index sars_20000_twilight_dipper.idx --amplicon-depth SRR19707934.amplicon_stacks.tsv --mask-reads-relative-frequency 0.01 em-delta-threshold  0.00001 --output SRR19707934 --threads 8
```

Panmap outputs an abundance file `SRR19707934.mgsr.abundance.out` containing the estimated abundances of nodes.

To get the lineage proportions, we can retrieve the sequences from the PanMAN and use pangolin to identify their lineages:

```bash
../build/bin/panmap ../examples/data/sars_20000_twilight_dipper.panman --dump-sequences "$(cut -f1 SRR19707934.mgsr.abundance.out | tr ',' '\n' | grep '\S' | tr '\n' ' ' | sed 's/ $//g')" --output SRR19707934
pangolin pangolin SRR19707934.dump-sequences.fa --outfile SRR19707934.pangolin.csv
python3 ../scripts/assign_lineage.py SRR19707934.mgsr.abundance.out SRR19707934.pangolin.csv > SRR19707934.mgsr.lineage.abundance.out
```