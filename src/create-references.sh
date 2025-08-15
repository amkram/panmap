#!/bin/bash

# Script to create reference FASTA files for the alignment process
# This extracts reference sequences from panmap files

# Paths to binaries and files
panmap_binary="../build/bin/panmap"
tb_panman="../data/tb_400.panman"

tb_mask_bed="R00000039_repregions.bed"

tb_ref="GCF_000195955.2_ASM19595v2_genomic.fna"
# Define paths for reference FASTA and GFF used by liftoff
ref_fa="${tb_ref}"
ref_gff="gff_files/$(basename "${tb_mask_bed}" .bed).gff" # Use the GFF file created from BED conversion

sars_ref="GCA_009858895.3_ASM985889v3_genomic.fna"

# Create output directories
mkdir -p reference annotations intermediate_files gff_files

# Convert BED to GFF using bed2gtf
echo "Converting BED to GFF format..."
set -x
# First convert Windows line endings (CRLF) to Unix line endings (LF)
tr -d '\r' <"${tb_mask_bed}" >"${tb_mask_bed}.unix"
# Simple BED to GFF conversion for repetitive regions
# Format: convert BED (chr start end) to GFF (chr source feature start end score strand phase attributes)
awk -v OFS="\t" '{print $1, "R00000039_repregions", "masked_region", $2+1, $3, ".", ".", ".", "ID=masked_"NR";"}' "${tb_mask_bed}.unix" >"gff_files/$(basename "${tb_mask_bed}" .bed).gff"
rm "${tb_mask_bed}.unix"

echo "Converted ${tb_mask_bed} to GFF format"

set +x

# Function to calculate the percentage of bases masked
calculate_masked_percentage() {
	local gff_file="$1"
	local fasta_file="$2"

	# Sum the lengths of all masked_region features in the GFF file
	echo "Calculating masked regions for ${gff_file}..."
	local masked_sum=$(awk -F'\t' '{if ($3 == "masked_region") sum += $5 - $4 + 1} END {print sum}' "${gff_file}")

	# Get the reference sequence length
	echo "Calculating reference length for ${fasta_file}..."
	local ref_length=$(grep -v "^>" "${fasta_file}" | tr -d '\n' | wc -c)

	# Calculate percentage
	local percentage=$(awk "BEGIN {printf \"%.2f\", (${masked_sum} / ${ref_length}) * 100}")

	echo "Results for $(basename "${gff_file}"):"
	echo "  Total masked bases: ${masked_sum}"
	echo "  Reference length: ${ref_length}"
	echo "  Percentage masked: ${percentage}%"
	echo ""
}

echo "masked_region" >features.txt

while read -r line; do
	# Extract node ID from each line of the batch file
	id=$(echo "${line}" | cut -f 4)
	echo "Processing node: ${id}"

	# Step 1: Dump the node's FASTA sequence using panmap
	echo "Extracting FASTA for node ${id}..."
	${panmap_binary} "$tb_panman" --dump-sequence "${id}"
	node_fa="${tb_panman}.${id}.fa"

	if [[ -f "${node_fa}" ]]; then
		echo "Successfully extracted FASTA for ${id}"

		# Create symlink with cleaner name
		ln -sf "../${node_fa}" "reference/${id}.fa"

		# Step 2: Use liftoff to transfer annotations from reference to this node
		echo "Running liftoff for ${id}..."
		liftoff "${node_fa}" "${ref_fa}" -g "${ref_gff}" -o "annotations/${id}.gff" \
			-f features.txt -dir "intermediate_files/${id}" -p 4

		# Calculate masking percentage
		if [[ -f "annotations/${id}.gff" ]]; then
			calculate_masked_percentage "annotations/${id}.gff" "${node_fa}"
		fi

		echo "Completed processing for ${id}"
	else
		echo "Failed to extract sequence for ${id}"
	fi
done <tb_batch.tsv

# Also calculate masking for the main mask GFF if it exists
if [[ -f "gff_files/tb_mask.gff" ]] && [[ -f "${ref_fa}" ]]; then
	calculate_masked_percentage "gff_files/tb_mask.gff" "${ref_fa}"
fi
