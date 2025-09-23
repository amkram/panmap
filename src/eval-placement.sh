#!/bin/bash

prefix=/private/groups/corbettlab/alex/panmap-testing-docker/LATEST_RESULTS/simulated_data/sim_generaltest

rm rsv_batch.tsv sars_batch.tsv tb_batch.tsv
awk -F'\t' '{print $1,$2,$3,$8}' ../../eval_data/rsv_pre.tsv | while read id num_mutations num_reads true_placement; do
    echo -e "rsv_4000_${id}_${num_mutations}_muts_${num_reads}_reads\t${prefix}_rsv_4000_${id}_R1.fastq\t${prefix}_rsv_4000_${id}_R2.fastq\t${true_placement}" >> rsv_batch.tsv
done

awk -F'\t' '{print $1,$2,$3,$7}' ../../eval_data/sars_pre.tsv | while read id num_mutations num_reads true_placement; do
    echo -e "sars_20000_${id}_${num_mutations}_muts_${num_reads}_reads\t${prefix}_sars_20000_${id}_R1.fastq\t${prefix}_sars_20000_${id}_R2.fastq\t${true_placement}" >> sars_batch.tsv
done

awk -F'\t' '{print $1,$2,$3,$7}' ../../eval_data/tb_pre.tsv | while read id num_mutations num_reads true_placement; do
    if [ "$true_placement" != "" ]; then
        echo -e "tb_400_${id}_${num_mutations}_muts_${num_reads}_reads\t${prefix}_tb_400_${id}_R1.fastq\t${prefix}_tb_400_${id}_R2.fastq\t${true_placement}" >> tb_batch.tsv
    fi
done

set -x

# while read -r line; do
#     sample_id=$(echo "$line" | cut -f 1)
#     r1=$(echo "$line" | cut -f 2)
#     r2=$(echo "$line" | cut -f 3)
#     true_node=$(echo "$line" | cut -f 4)

#     ../build/bin/panmap ../data/tb_400.panman $r1 $r2 -p "tb_placement_${sample_id}"
# done < tb_batch.tsv

# ./create-references.sh

# ./process-placements.sh

# align each placement node fasta with the matching true node fasta
# with minimap2 then call variants (no minimum size)
# collapsing indels into one variant
# report the number of variants, the number of aligned bases, and the proportion (variants / aligned bases)

set +x
