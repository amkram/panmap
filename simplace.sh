
cp ../manuscript/data/ref-panmans-new-format/rsv_4K.panman .

# build/bin/panmap -f -k 32 -s 8 rsv_4K.panman

build/bin/panmap rsv_4K.panman --dump-random-node > tmpOut.txt

random_node_id=$(grep -oP 'Random node \K\S+' tmpOut.txt)
fasta_filename=$(grep -oP 'sequence written to \K\S+' tmpOut.txt)

echo "Parsed node: $random_node_id"
echo "FASTA filename: $fasta_filename"


docker run -v $(pwd):/data  \
    quay.io/biocontainers/insilicoseq:2.0.0--pyh7cba7a3_0  \
    iss generate -n 6 --genomes /data/$fasta_filename --model novaseq --output /data/"${random_node_id}_reads"

r1="${random_node_id}_R1.fastq"
r2="${random_node_id}_R2.fastq"

echo "Simulated 3 read-pairs for $random_node_id"
echo "R1: $r1"
echo "R2: $r2"



