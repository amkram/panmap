
#!/bin/bash

# Function to simulate reads for a single random node
simulate_single() {
    local num_reads=$1
    local iteration=$2
    
    echo "=== Simulation $iteration: Generating $num_reads reads ===" >&2
    
    # Get random node and generate FASTA
    build/bin/panmap rsv_4K.panman --dump-random-node > tmpOut_${iteration}.txt
    
    random_node_id=$(grep -oP 'Random node \K\S+' tmpOut_${iteration}.txt)
    fasta_filename=$(grep -oP 'sequence written to \K\S+' tmpOut_${iteration}.txt)
    
    echo "  Node: $random_node_id" >&2
    echo "  FASTA: $fasta_filename" >&2
    
    # Generate reads using Docker
    docker run -v $(pwd):/data  \
        quay.io/biocontainers/insilicoseq:2.0.0--pyh7cba7a3_0  \
        iss generate -n $num_reads --genomes /data/$fasta_filename --model novaseq --output /data/"sim_${iteration}_${random_node_id}_reads" >&2
    
    # Store the paths for batch file
    r1="sim_${iteration}_${random_node_id}_reads_R1.fastq"
    r2="sim_${iteration}_${random_node_id}_reads_R2.fastq"
    
    # Output ONLY the batch entry to stdout
    echo "sim_${iteration}_${random_node_id},$(pwd)/${r1},$(pwd)/${r2}"
    
    # Clean up temp file
    rm -f tmpOut_${iteration}.txt
}

# Main wrapper function
run_batch_simulation() {
    local num_reads=${1:-6}  # Default to 6 reads if not specified
    local num_iterations=${2:-10}  # Default to 10 iterations if not specified
    
    echo "Starting batch simulation: $num_iterations iterations with $num_reads reads each"
    
    
    # Create batch file
    batch_file="batch_placement.txt"
    > $batch_file  # Clear/create the file
    
    # Run simulations and collect batch entries
    for i in $(seq 1 $num_iterations); do
        batch_entry=$(simulate_single $num_reads $i)
        echo "$batch_entry" >> $batch_file
    done
    
    echo "=== Batch file created: $batch_file ==="
    echo "Contents:"
    cat $batch_file
    
    echo "=== Starting batch placement ==="
    # Run batch placement
    build/bin/panmap rsv_4K.panman --batch $batch_file
    
    echo "=== Batch placement completed ==="
    echo "Results should be in files matching: batch_results_sim_*_placement"
}

# If script is run directly, execute the wrapper
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    # Check if arguments provided
    if [[ $# -eq 0 ]]; then
        echo "Usage: $0 [num_reads] [num_iterations]"
        echo "  num_reads: Number of reads to simulate per iteration (default: 6)"
        echo "  num_iterations: Number of simulation iterations (default: 10)"
        echo ""
        echo "Running with defaults..."
        run_batch_simulation
    else
        run_batch_simulation "$@"
    fi
fi



