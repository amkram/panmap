import sys
from collections import defaultdict

abundance_file = sys.argv[1]
lineage_report = sys.argv[2]

node_to_lineage = {}
with open(lineage_report, 'r') as f:
  f.readline()
  for line in f:
    node, lineage = line.strip().split(',')[0:2]
    node_to_lineage[node.strip()] = lineage.strip()


lineage_abundance = defaultdict(float)
with open(abundance_file, 'r') as f:
  for line in f:
    fields = line.strip().split('\t')
    nodes = fields[0].split(',')
    num_nodes = len(nodes)
    abundance = float(fields[1])
    for node in nodes:
      if node not in node_to_lineage:
        print(f"Node {node} not found in lineage report")
        exit(1)
      lineage_abundance[node_to_lineage[node]] += abundance / num_nodes
  
sorted_lineage_abundance = sorted(lineage_abundance.items(), key=lambda x: x[1], reverse=True)
for lineage, abundance in sorted_lineage_abundance:
  print(f"{lineage}\t{abundance}")