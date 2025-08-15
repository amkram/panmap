#!/bin/bash

# Main script to run the entire workflow for processing placements and creating chain files

set -e # Exit on error

# Make scripts executable
chmod +x create_references.sh
chmod +x process_placements.sh

echo "============================================="
echo "STEP 1: Creating reference FASTA files"
echo "============================================="
./create_references.sh

echo "============================================="
echo "STEP 2: Processing batch files and creating chain files"
echo "============================================="
./process_placements.sh

echo "============================================="
echo "WORKFLOW COMPLETED"
echo "============================================="
echo "Chain files have been created in the chain_files directory."
echo "Check for any warnings or errors in the output above."
