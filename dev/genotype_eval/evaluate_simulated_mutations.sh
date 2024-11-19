#! /bin/bash

# Command line options
TREE_PATH=$1
COVERAGES=$2
NUM_REPLICATES=$3
OUTPUT_PREFIX=$4

# Parse coverage string and convert to array of coverages and numbers of reads
IFS=',' read -r -a COVERAGE_ARRAY <<< "$COVERAGES"
READS_PER_COVERAGE=200
NUM_READS_ARRAY=()

for COVERAGE in "${COVERAGE_ARRAY[@]}"; do
    NUM_READS=$((COVERAGE * READS_PER_COVERAGE))
    NUM_READS_ARRAY+=($NUM_READS)
done

# Create output directory
OUTPUT_DIR=$(dirname "$OUTPUT_PREFIX")
OUTPUT_PREFIX=$(basename "$OUTPUT_PREFIX")
mkdir -p \
  "$OUTPUT_DIR" \
  "$OUTPUT_DIR/${OUTPUT_PREFIX}_bam" \
  "$OUTPUT_DIR/${OUTPUT_PREFIX}_bcf_vcf" \
  "$OUTPUT_DIR/${OUTPUT_PREFIX}_panmap_output" \
  "$OUTPUT_DIR/${OUTPUT_PREFIX}_roc_bcf_vcf" \
  "$OUTPUT_DIR/${OUTPUT_PREFIX}_roc_panmap"

# Set up for running simulated data
PANMAT_PATH="../panmans/sars_20000.panman"
REF_STRAIN="USA/CA-CDPH-500076119/2022|OP769690.1|2022-08-27"
REF_STRAIN_MODIFIED=$(echo "$REF_STRAIN" | tr '/' '_')
REF_STRAIN_MODIFIED=$(echo "$REF_STRAIN_MODIFIED" | tr ' ' '_')


simulate \
  --panmat "$PANMAT_PATH" \
  --ref "$REF_STRAIN" \
  --out_dir "$OUTPUT_DIR" \
  --prefix "$OUTPUT_PREFIX" \
  --mut_spec "$PANMAT_PATH.mm" \
  --mut_spec_type snp \
  --rep "$NUM_REPLICATES" \
  --model NovaSeq \
  --no-reads

for ((i=0; i<NUM_REPLICATES; i++)); do
  sed "s/ref/${REF_STRAIN_MODIFIED}/g" "${OUTPUT_DIR}/${OUTPUT_PREFIX}_vcfTrue/${REF_STRAIN_MODIFIED}.var.${i}.vcf" > "${OUTPUT_DIR}/${OUTPUT_PREFIX}_vcfTrue/${REF_STRAIN_MODIFIED}.var.${i}.modified.vcf"
  mv "${OUTPUT_DIR}/${OUTPUT_PREFIX}_vcfTrue/${REF_STRAIN_MODIFIED}.var.${i}.modified.vcf" "${OUTPUT_DIR}/${OUTPUT_PREFIX}_vcfTrue/${REF_STRAIN_MODIFIED}.var.${i}.vcf"
  for j in "${!NUM_READS_ARRAY[@]}"; do
    NUM_READS=${NUM_READS_ARRAY[$j]}
    COVERAGE=${COVERAGE_ARRAY[$j]}
    READS_PREFIX="$OUTPUT_DIR/${OUTPUT_PREFIX}_reads/${REF_STRAIN_MODIFIED}.reads.${i}.${COVERAGE}x"
    echo "Simulating with $NUM_READS reads for replicate $i at "$COVERAGE"x coverage..."
    
    iss generate \
      --model NovaSeq \
      --genomes "$OUTPUT_DIR/${OUTPUT_PREFIX}_varFasta/${REF_STRAIN_MODIFIED}.var.${i}.fa" \
      --output "$READS_PREFIX" \
      --n_reads "$NUM_READS" \
      --cpus 4 \
      --quiet

    READ_R1="${READS_PREFIX}_R1.fastq"
    READ_R2="${READS_PREFIX}_R2.fastq"

    echo "Running bwa aln and sampe to generate SAM files..."
    # Index the reference fasta file using BWA
    bwa index "$OUTPUT_DIR/${OUTPUT_PREFIX}_refFasta/${REF_STRAIN_MODIFIED}.fa" > /dev/null 2>&1

    # Align reads using BWA
    bwa aln -t 4 "$OUTPUT_DIR/${OUTPUT_PREFIX}_refFasta/${REF_STRAIN_MODIFIED}.fa" "$READ_R1" > "${READS_PREFIX}_R1.sai" 2> /dev/null
    bwa aln -t 4 "$OUTPUT_DIR/${OUTPUT_PREFIX}_refFasta/${REF_STRAIN_MODIFIED}.fa" "$READ_R2" > "${READS_PREFIX}_R2.sai" 2> /dev/null

    # Generate SAM file using BWA sampe
    bwa sampe "$OUTPUT_DIR/${OUTPUT_PREFIX}_refFasta/${REF_STRAIN_MODIFIED}.fa" "${READS_PREFIX}_R1.sai" "${READS_PREFIX}_R2.sai" "$READ_R1" "$READ_R2" 2> /dev/null | samtools sort -o "${OUTPUT_DIR}/${OUTPUT_PREFIX}_bam/${REF_STRAIN_MODIFIED}.${i}.${COVERAGE}x.sorted.sam" > /dev/null 2>&1

    echo "Calling variants with bcftools to generate VCF files..."
    
    # Convert SAM to BAM
    samtools view -Sb "${OUTPUT_DIR}/${OUTPUT_PREFIX}_bam/${REF_STRAIN_MODIFIED}.${i}.${COVERAGE}x.sorted.sam" > "${OUTPUT_DIR}/${OUTPUT_PREFIX}_bam/${REF_STRAIN_MODIFIED}.${i}.${COVERAGE}x.sorted.bam" 2> /dev/null
    
    # Sort BAM file
    samtools sort "${OUTPUT_DIR}/${OUTPUT_PREFIX}_bam/${REF_STRAIN_MODIFIED}.${i}.${COVERAGE}x.sorted.bam" -o "${OUTPUT_DIR}/${OUTPUT_PREFIX}_bam/${REF_STRAIN_MODIFIED}.${i}.${COVERAGE}x.sorted.bam" 2> /dev/null
    
    # Index the sorted BAM file
    samtools index "${OUTPUT_DIR}/${OUTPUT_PREFIX}_bam/${REF_STRAIN_MODIFIED}.${i}.${COVERAGE}x.sorted.bam" 2> /dev/null
    
    # Call variants using bcftools
    bcftools mpileup -Ou -f "${OUTPUT_DIR}/${OUTPUT_PREFIX}_refFasta/${REF_STRAIN_MODIFIED}.fa" "${OUTPUT_DIR}/${OUTPUT_PREFIX}_bam/${REF_STRAIN_MODIFIED}.${i}.${COVERAGE}x.sorted.bam" 2> /dev/null > "${OUTPUT_DIR}/${OUTPUT_PREFIX}_bcf_vcf/${REF_STRAIN_MODIFIED}.variants.${i}.${COVERAGE}x.bcf"
    
    bcftools call --ploidy 1 -V indels -cv -Ov -o "${OUTPUT_DIR}/${OUTPUT_PREFIX}_bcf_vcf/${REF_STRAIN_MODIFIED}.variants.${i}.${COVERAGE}x.varonly.vcf" "${OUTPUT_DIR}/${OUTPUT_PREFIX}_bcf_vcf/${REF_STRAIN_MODIFIED}.variants.${i}.${COVERAGE}x.bcf" > /dev/null 2>&1

    bcftools call --ploidy 1 -V indels -cA -Ov -o "${OUTPUT_DIR}/${OUTPUT_PREFIX}_bcf_vcf/${REF_STRAIN_MODIFIED}.variants.${i}.${COVERAGE}x.vcf" "${OUTPUT_DIR}/${OUTPUT_PREFIX}_bcf_vcf/${REF_STRAIN_MODIFIED}.variants.${i}.${COVERAGE}x.bcf" > /dev/null 2>&1

    echo "Applying mutation spectrum to bcf vcf..."

    echo "Calling variants with panmap mutation spectrum..."
    python3 apply_mut_spectrum.py \
      "${OUTPUT_DIR}/${OUTPUT_PREFIX}_bcf_vcf/${REF_STRAIN_MODIFIED}.variants.${i}.${COVERAGE}x.vcf" \
      "${PANMAT_PATH}.mm" \
      "${OUTPUT_DIR}/${OUTPUT_PREFIX}_bcf_vcf/${REF_STRAIN_MODIFIED}.variants.${i}.${COVERAGE}x.mutmat_s1.vcf" \
      1 \
      no


    echo "Creating roc files for bcf vcf..."
    python3 filter_genotype_var.py \
      "${OUTPUT_DIR}/${OUTPUT_PREFIX}_bcf_vcf/${REF_STRAIN_MODIFIED}.variants.${i}.${COVERAGE}x.varonly.vcf" \
      "${OUTPUT_DIR}/${OUTPUT_PREFIX}_bcf_vcf/${REF_STRAIN_MODIFIED}.variants.${i}.${COVERAGE}x.varonly.filtered.vcf"
    vcfroc \
      -t "${OUTPUT_DIR}/${OUTPUT_PREFIX}_vcfTrue/${REF_STRAIN_MODIFIED}.var.${i}.vcf" \
      -r "${OUTPUT_DIR}/${OUTPUT_PREFIX}_refFasta/${REF_STRAIN_MODIFIED}.fa" \
      "${OUTPUT_DIR}/${OUTPUT_PREFIX}_bcf_vcf/${REF_STRAIN_MODIFIED}.variants.${i}.${COVERAGE}x.varonly.filtered.vcf" \
      > "${OUTPUT_DIR}/${OUTPUT_PREFIX}_roc_bcf_vcf/${REF_STRAIN_MODIFIED}.${i}.${COVERAGE}x.varonly.roc"

    echo "Creating roc files for bcf vcf with mutmat applied..."
    python3 filter_genotype_var.py \
      "${OUTPUT_DIR}/${OUTPUT_PREFIX}_bcf_vcf/${REF_STRAIN_MODIFIED}.variants.${i}.${COVERAGE}x.mutmat_s1.vcf" \
      "${OUTPUT_DIR}/${OUTPUT_PREFIX}_bcf_vcf/${REF_STRAIN_MODIFIED}.variants.${i}.${COVERAGE}x.mutmat_s1.filtered.vcf"
    vcfroc \
      -t "${OUTPUT_DIR}/${OUTPUT_PREFIX}_vcfTrue/${REF_STRAIN_MODIFIED}.var.${i}.vcf" \
      -r "${OUTPUT_DIR}/${OUTPUT_PREFIX}_refFasta/${REF_STRAIN_MODIFIED}.fa" \
      "${OUTPUT_DIR}/${OUTPUT_PREFIX}_bcf_vcf/${REF_STRAIN_MODIFIED}.variants.${i}.${COVERAGE}x.mutmat_s1.filtered.vcf" \
      > "${OUTPUT_DIR}/${OUTPUT_PREFIX}_roc_bcf_vcf/${REF_STRAIN_MODIFIED}.${i}.${COVERAGE}x.mutmat_s1.roc"

  done
done

SCALEING_FACTORS=(16.5)
for scaling_factor in "${SCALEING_FACTORS[@]}"; do  
  for ((i=0; i<NUM_REPLICATES; i++)); do
    for j in "${!NUM_READS_ARRAY[@]}"; do
      NUM_READS=${NUM_READS_ARRAY[$j]}
      COVERAGE=${COVERAGE_ARRAY[$j]}
      READS_PREFIX="$OUTPUT_DIR/${OUTPUT_PREFIX}_reads/${REF_STRAIN_MODIFIED}.reads.${i}.${COVERAGE}x"
      READ_R1="${READS_PREFIX}_R1.fastq"
      READ_R2="${READS_PREFIX}_R2.fastq"

      echo "Applying scaled mutation spectrum scaled by $scaling_factor for replicate $i at ${COVERAGE}x coverage..."
      
      echo "Calling variants with panmap mutation spectrum..."
      python3 apply_mut_spectrum.py \
        "${OUTPUT_DIR}/${OUTPUT_PREFIX}_bcf_vcf/${REF_STRAIN_MODIFIED}.variants.${i}.${COVERAGE}x.vcf" \
        "${PANMAT_PATH}.mm" \
        "${OUTPUT_DIR}/${OUTPUT_PREFIX}_bcf_vcf/${REF_STRAIN_MODIFIED}.variants.${i}.${COVERAGE}x.mutmat_s${scaling_factor}.vcf" \
        "$scaling_factor" \
        no

      echo "Creating roc files for bcf vcf with mutmat applied..."
      python3 filter_genotype_var.py \
        "${OUTPUT_DIR}/${OUTPUT_PREFIX}_bcf_vcf/${REF_STRAIN_MODIFIED}.variants.${i}.${COVERAGE}x.mutmat_s${scaling_factor}.vcf" \
        "${OUTPUT_DIR}/${OUTPUT_PREFIX}_bcf_vcf/${REF_STRAIN_MODIFIED}.variants.${i}.${COVERAGE}x.mutmat_s${scaling_factor}.filtered.vcf"
      vcfroc \
        -t "${OUTPUT_DIR}/${OUTPUT_PREFIX}_vcfTrue/${REF_STRAIN_MODIFIED}.var.${i}.vcf" \
        -r "${OUTPUT_DIR}/${OUTPUT_PREFIX}_refFasta/${REF_STRAIN_MODIFIED}.fa" \
        "${OUTPUT_DIR}/${OUTPUT_PREFIX}_bcf_vcf/${REF_STRAIN_MODIFIED}.variants.${i}.${COVERAGE}x.mutmat_s${scaling_factor}.filtered.vcf" \
        > "${OUTPUT_DIR}/${OUTPUT_PREFIX}_roc_bcf_vcf/${REF_STRAIN_MODIFIED}.${i}.${COVERAGE}x.mutmat_s${scaling_factor}.roc"
    done
  done
done

