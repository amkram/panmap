NUM_READS=${1:-20}

cp ../manuscript/data/ref-panmans-new-format/rsv_4K.panman .

# build/bin/panmap -f -k 32 -s 8 rsv_4K.panman

build/bin/panmap rsv_4K.panman --dump-random-node > tmpOut.txt

random_node_id=$(grep -oP 'Random node \K\S+' tmpOut.txt)
fasta_filename=$(grep -oP 'sequence written to \K\S+' tmpOut.txt)

echo "Parsed node: $random_node_id"
echo "FASTA filename: $fasta_filename"

docker run --rm \
  -v "$(pwd):/data" \
  --entrypoint iss \
  quay.io/biocontainers/insilicoseq:2.0.1--pyh7cba7a3_0 \
  generate -n "${NUM_READS:?NUM_READS not set}" \
    --genomes "/data/${fasta_filename:?fasta_filename not set}" \
    --mode perfect \
    --output "/data/reads/${random_node_id:?random_node_id not set}_reads"

r1="reads/${random_node_id}_reads_R1.fastq"
r2="reads/${random_node_id}_reads_R2.fastq"

echo "Simulated $NUM_READS read-pairs for $random_node_id"
echo "R1: $r1"
echo "R2: $r2"



