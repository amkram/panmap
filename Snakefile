"""
Panmap benchmarking workflow - simplified version
"""
import hashlib, pathlib, os, re, subprocess, shutil
from pathlib import Path

OUTPUT_DIR = "workflow_output"

# ============================================================================
# CONFIGURATION
# ============================================================================
PANGENOMES = {
    'rsv_4K': {'path': 'panmans/rsv_4K.panman', 'genome_size': 15000},
    'tb_400': {'path': 'panmans/tb_400.panman', 'genome_size': 4400000}
}

CONFIG = {
    'pangenomes': ['rsv_4K'],
    'k_values': range(11, 60, 8),
    's_values': [4, 8],
    'l_values': [1, 3],
    'include_internal': [True],
    'coverage_levels': [.1, 0.5, 1, 10, 100, 1000],
    'mutation_rates': [0.0001, 0.0005, 0.001],
    'replicates': 100,
    'model': 'NovaSeq',
    'read_length': 150
}

# ============================================================================
# EXPERIMENT GENERATION
# ============================================================================
def generate_experiments():
    experiments = []
    for i, (pan, k, s, l, inc_int, mut_rate) in enumerate(
        (pan, k, s, l, inc, mut) 
        for pan in CONFIG['pangenomes']
        for k in CONFIG['k_values']
        for s in CONFIG['s_values']
        for l in CONFIG['l_values']
        for inc in CONFIG['include_internal']
        for mut in CONFIG['mutation_rates']
    ):
        genome_size = PANGENOMES[pan]['genome_size']
        # index_tag excludes mutation_rate since index doesn't depend on it
        index_tag = f"k{k}_s{s}_l{l}_{'int' if inc_int else 'noint'}"
        experiments.append({
            'id': f'exp{i}',
            'pangenome_name': pan,
            'panman_path': PANGENOMES[pan]['path'],
            'genome_size': genome_size,
            'pan_stem': pan,
            'k': k, 's': s, 'l': l,
            'include_internal': inc_int,
            'mutation_rate': mut_rate,
            'model': CONFIG['model'],
            'coverage_levels': CONFIG['coverage_levels'],
            'num_reads_values': [int((genome_size * cov) / CONFIG['read_length']) for cov in CONFIG['coverage_levels']],
            'replicates': CONFIG['replicates'],
            'tag': f"k{k}_s{s}_l{l}_{'int' if inc_int else 'noint'}" + (f"_mut{mut_rate}" if mut_rate > 0 else ""),
            'index_tag': index_tag  # for shared index lookup
        })
    return experiments

EXPERIMENTS = generate_experiments()
EXP_BY_ID = {e['id']: e for e in EXPERIMENTS}

def exp_root(eid, pan_stem, tag):
    return f"{OUTPUT_DIR}/experiments/{eid}/{pan_stem}/{tag}"

# Generate unique index combinations (pangenome + k/s/l, independent of mutation_rate)
INDEX_COMBOS = list({(e['pan_stem'], e['k'], e['s'], e['l'], e['include_internal'], e['index_tag'], e['panman_path']) 
                     for e in EXPERIMENTS})

def index_root(pan_stem, index_tag):
    return f"{OUTPUT_DIR}/indexes/{pan_stem}/{index_tag}"

# Generate all read combinations
READS = [(e['id'], e['pan_stem'], e['tag'], cov, n, rep) 
         for e in EXPERIMENTS 
         for cov, n in zip(e['coverage_levels'], e['num_reads_values'])
         for rep in range(e['replicates'])]

print(f"[config] {len(EXPERIMENTS)} experiments, {len(INDEX_COMBOS)} unique indexes, {len(READS)} total read sets")

# ============================================================================
# FILE PATHS
# ============================================================================
def paths(template):
    return [template.format(eid=eid, pan=pan, tag=tag, cov=cov, n=n, rep=rep) 
            for (eid, pan, tag, cov, n, rep) in READS]

EXP_ROOT = f"{OUTPUT_DIR}/experiments/{{eid}}/{{pan}}/{{tag}}"
GENOME_FASTA = paths(f"{EXP_ROOT}/genomes/reads/cov{{cov}}_{{n}}_rep{{rep}}/random_node.fasta")
MUT_FASTA = paths(f"{EXP_ROOT}/mutgenomes/reads/cov{{cov}}_{{n}}_rep{{rep}}/mutated.fasta")
FASTQ_R1 = paths(f"{EXP_ROOT}/reads/cov{{cov}}_{{n}}_rep{{rep}}_R1.fastq")
FASTQ_R2 = paths(f"{EXP_ROOT}/reads/cov{{cov}}_{{n}}_rep{{rep}}_R2.fastq")
META_FILES = paths(f"{EXP_ROOT}/results/reads/cov{{cov}}_{{n}}_rep{{rep}}.txt")
PLACEMENTS = paths(f"{EXP_ROOT}/placements/reads/cov{{cov}}_{{n}}_rep{{rep}}/placements.tsv")
DETAILED = paths(f"{EXP_ROOT}/placements/reads/cov{{cov}}_{{n}}_rep{{rep}}/detailed.tsv")
ACCURACY = paths(f"{EXP_ROOT}/alignments/reads/cov{{cov}}_{{n}}_rep{{rep}}/accuracy.tsv")
# Shared indexes - one per (pangenome, k, s, l) combination
INDEX_FILES = [f"{index_root(pan, idx_tag)}/index.pmi" for (pan, k, s, l, inc, idx_tag, path) in INDEX_COMBOS]
MM_FILES = [f"{index_root(pan, idx_tag)}/{pan}.panman.mm" for (pan, k, s, l, inc, idx_tag, path) in INDEX_COMBOS]
INDEX_TIME_LOGS = [f"{index_root(pan, idx_tag)}/index_time.log" for (pan, k, s, l, inc, idx_tag, path) in INDEX_COMBOS]
PLACEMENT_TIME_LOGS = paths(f"{EXP_ROOT}/placements/reads/cov{{cov}}_{{n}}_rep{{rep}}/time.log")

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================
def sanitize_filename(name):
    return re.sub(r'[\/\\:\*\?"<>|]', '_', name)

def parse_time_log(log_file):
    """Parse /usr/bin/time -v output"""
    stats = {}
    try:
        with open(log_file) as f:
            for line in f:
                if 'Elapsed (wall clock) time' in line:
                    m = re.search(r'(\d+):(\d+):(\d+\.\d+)', line) or re.search(r'(\d+):(\d+\.\d+)', line)
                    if m:
                        g = m.groups()
                        stats['wall_time_sec'] = sum(float(x) * (60 ** (len(g)-1-i)) for i, x in enumerate(g))
                elif 'Maximum resident set size' in line:
                    m = re.search(r':\s*(\d+)', line)
                    if m: stats['max_rss_mb'] = int(m.group(1)) / 1024.0
                elif 'User time (seconds)' in line:
                    m = re.search(r':\s*([\d.]+)', line)
                    if m: stats['user_time_sec'] = float(m.group(1))
                elif 'System time (seconds)' in line:
                    m = re.search(r':\s*([\d.]+)', line)
                    if m: stats['system_time_sec'] = float(m.group(1))
    except FileNotFoundError:
        pass
    return stats

def count_variants_with_exclusion(cs_field, target_start, target_end, exclude_bp):
    """Count SNPs/indels from CS string, excluding bp from ends"""
    snps = indels = 0
    if target_end - target_start < 2 * exclude_bp:
        return snps, indels
    exclude_start, exclude_end = target_start + exclude_bp, target_end - exclude_bp
    pos, i = target_start, 0
    while i < len(cs_field):
        c = cs_field[i]
        if c == ':':
            i += 1
            num = ""
            while i < len(cs_field) and cs_field[i].isdigit():
                num += cs_field[i]; i += 1
            if num: pos += int(num)
        elif c == '*':
            if exclude_start <= pos < exclude_end: snps += 1
            pos += 1; i += 3 if i + 2 < len(cs_field) else len(cs_field)
        elif c == '+':
            if exclude_start <= pos < exclude_end: indels += 1
            i += 1
            while i < len(cs_field) and cs_field[i] in 'acgtACGT': i += 1
        elif c == '-':
            if exclude_start <= pos < exclude_end: indels += 1
            i += 1; del_len = 0
            while i < len(cs_field) and cs_field[i] in 'acgtACGT': del_len += 1; i += 1
            pos += del_len
        else: i += 1
    return snps, indels

# ============================================================================
# RULES
# ============================================================================
rule all:
    input:
        f"{OUTPUT_DIR}/reports/experiments_summary.tsv",
        f"{OUTPUT_DIR}/reports/placements_summary.tsv",
        f"{OUTPUT_DIR}/reports/placements_by_metric.tsv",
        f"{OUTPUT_DIR}/reports/alignment_accuracy_summary.tsv",
        f"{OUTPUT_DIR}/plots/placement_accuracy_plots.done",
        # Per-pangenome accuracy plots
        expand(f"{OUTPUT_DIR}/plots/{{pan}}/k_analysis/accuracy_by_k_l1.png", pan=PANGENOMES.keys()),
        expand(f"{OUTPUT_DIR}/plots/{{pan}}/k_analysis/accuracy_by_k_l3.png", pan=PANGENOMES.keys()),
        f"{OUTPUT_DIR}/reports/index_performance_summary.tsv",
        f"{OUTPUT_DIR}/reports/placement_performance_summary.tsv",
        # Comprehensive results tables
        f"{OUTPUT_DIR}/reports/index_stats.tsv",
        f"{OUTPUT_DIR}/reports/placement_stats.tsv",
        f"{OUTPUT_DIR}/reports/combined_results.tsv",
        # Per-pangenome performance plots
        expand(f"{OUTPUT_DIR}/plots/{{pan}}/performance/index_time_by_k.png", pan=PANGENOMES.keys()),
        expand(f"{OUTPUT_DIR}/plots/{{pan}}/performance/placement_time_by_k.png", pan=PANGENOMES.keys()),
        *GENOME_FASTA, *MUT_FASTA, *INDEX_FILES, *MM_FILES, *FASTQ_R1, *FASTQ_R2,
        *META_FILES, *PLACEMENTS, *DETAILED, *ACCURACY

rule experiments_summary:
    output: tsv=f"{OUTPUT_DIR}/reports/experiments_summary.tsv"
    run:
        Path(output.tsv).parent.mkdir(parents=True, exist_ok=True)
        with open(output.tsv, 'w') as f:
            f.write('\t'.join(['id','tag','pan_stem','panman_path','include_internal','k','s','l','mutation_rate','genome_size','num_reads_values','replicates','model'])+'\n')
            for e in EXPERIMENTS:
                f.write('\t'.join(map(str, [e['id'], e['tag'], e['pan_stem'], e['panman_path'], str(e['include_internal']).lower(), e['k'], e['s'], e['l'], e['mutation_rate'], e['genome_size'], ';'.join(map(str,e['num_reads_values'])), e['replicates'], e['model']]))+'\n')

rule dump_random_node:
    input:
        panmap_bin="build/bin/panmap",
        panman=lambda wc: EXP_BY_ID[wc.eid]['panman_path']
    output: f"{OUTPUT_DIR}/experiments/{{eid}}/{{pan_stem}}/{{tag}}/genomes/reads/cov{{cov}}_{{n}}_rep{{rep}}/random_node.fasta"
    params:
        seed=lambda wc: int(hashlib.md5(f"{wc.eid}_random_node_cov{wc.cov}_{wc.n}_rep{wc.rep}".encode()).hexdigest()[:8], 16) % (2**31 - 1)
    shell:
        r'''
        set -e; mkdir -p $(dirname {output})
        PANMAP=$(realpath {input.panmap_bin}); PANMAN=$(realpath {input.panman}); OUT=$(realpath -m {output})
        TMPDIR=$(mktemp -d); trap "rm -rf $TMPDIR" EXIT
        cp "$PANMAN" "$TMPDIR/temp.panman"; cd "$TMPDIR"
        "$PANMAP" temp.panman --dump-random-node --seed {params.seed} | tail -1
        mv temp.panman.random.*.fa "$OUT" 2>/dev/null || {{ echo "ERROR: No file generated"; exit 1; }}
        '''

rule index_experiment:
    input:
        bin="build/bin/panmap",
        pan=lambda wc: PANGENOMES[wc.pan_stem]['path']
    output:
        index=f"{OUTPUT_DIR}/indexes/{{pan_stem}}/{{index_tag}}/index.pmi",
        mm=f"{OUTPUT_DIR}/indexes/{{pan_stem}}/{{index_tag}}/{{pan_stem}}.panman.mm",
        time_log=f"{OUTPUT_DIR}/indexes/{{pan_stem}}/{{index_tag}}/index_time.log"
    params:
        pan=lambda wc: PANGENOMES[wc.pan_stem]['path'],
        k=lambda wc: int(re.search(r'k(\d+)', wc.index_tag).group(1)),
        s=lambda wc: int(re.search(r's(\d+)', wc.index_tag).group(1)),
        l=lambda wc: int(re.search(r'l(\d+)', wc.index_tag).group(1))
    shell:
        r'''
        set -e; mkdir -p $(dirname {output.index})
        [[ -f {output.index} && -s {output.index} ]] && {{ [[ ! -f {output.mm} ]] && cp {params.pan}.mm {output.mm}; exit 0; }}
        TMP=$(mktemp -d); trap "rm -rf $TMP" EXIT
        /usr/bin/time -v {input.bin} {params.pan} -k {params.k} -s {params.s} -l {params.l} -i "$TMP/idx.pmi" --stop index -f 2> {output.time_log}
        mv "$TMP/idx.pmi" {output.index}; cp {params.pan}.mm {output.mm}
        '''

rule simulate_reads:
    threads: 4
    input:
        fasta=f"{OUTPUT_DIR}/experiments/{{eid}}/{{pan_stem}}/{{tag}}/genomes/reads/cov{{cov}}_{{n}}_rep{{rep}}/random_node.fasta",
        panman=lambda wc: EXP_BY_ID[wc.eid]['panman_path'],
        simulate_bin="build/bin/simulate",
        mm=lambda wc: f"{index_root(wc.pan_stem, EXP_BY_ID[wc.eid]['index_tag'])}/{wc.pan_stem}.panman.mm"
    output:
        r1=f"{OUTPUT_DIR}/experiments/{{eid}}/{{pan_stem}}/{{tag}}/reads/cov{{cov}}_{{n}}_rep{{rep}}_R1.fastq",
        r2=f"{OUTPUT_DIR}/experiments/{{eid}}/{{pan_stem}}/{{tag}}/reads/cov{{cov}}_{{n}}_rep{{rep}}_R2.fastq",
        meta=f"{OUTPUT_DIR}/experiments/{{eid}}/{{pan_stem}}/{{tag}}/results/reads/cov{{cov}}_{{n}}_rep{{rep}}.txt",
        mut_fa=f"{OUTPUT_DIR}/experiments/{{eid}}/{{pan_stem}}/{{tag}}/mutgenomes/reads/cov{{cov}}_{{n}}_rep{{rep}}/mutated.fasta"
    params:
        out_dir=lambda wc: f"{OUTPUT_DIR}/experiments/{wc.eid}/{wc.pan_stem}/{wc.tag}/mutgenomes/reads/cov{wc.cov}_{wc.n}_rep{wc.rep}",
        seed=lambda wc: int(hashlib.md5(f"{wc.eid}_{wc.cov}_{wc.n}_{wc.rep}".encode()).hexdigest()[:8], 16) % (2**31 - 1),
        model=lambda wc: EXP_BY_ID[wc.eid]['model'],
        mut_rate=lambda wc: EXP_BY_ID[wc.eid].get('mutation_rate', 0.0)
    shell:
        r'''
        set -e; mkdir -p {params.out_dir} $(dirname {output.r1}) $(dirname {output.meta})
        node=$(head -1 {input.fasta} | sed 's/^>//')
        # Calculate actual genome length from FASTA (excluding header lines)
        genome_size=$(grep -v '^>' {input.fasta} | tr -d '\n' | wc -c)
        # Calculate mutations based on actual genome size and mutation rate
        snps=$(python3 -c "print(int({params.mut_rate} * $genome_size))")
        indels=$(python3 -c "print(int({params.mut_rate} * $genome_size / 4))")
        ins=$(((indels + 1) / 2)); del=$((indels / 2))
        mm_flag=""; [ -f {input.mm} ] && mm_flag="--mut_spec {input.mm} --mut_spec_type snp"
        {input.simulate_bin} --panmat {input.panman} --ref "$node" --out_dir {params.out_dir} --n_reads {wildcards.n} --rep 1 --model {params.model} --cpus 4 --seed {params.seed} --mutnum $snps $ins $del $mm_flag 2>&1 | tee {params.out_dir}/simulate.log
        find {params.out_dir} -name "*_R1.fastq" -exec cp {{}} {output.r1} \;
        find {params.out_dir} -name "*_R2.fastq" -exec cp {{}} {output.r2} \;
        find {params.out_dir} -name "*.fa" ! -name "*.fai" -exec cp {{}} {output.mut_fa} \;
        true_node=$(grep -oP '(curNodeID|Applied.*to node)\K[^\s]+$' {params.out_dir}/simulate.log | tail -1 || echo "$node")
        echo "id={wildcards.eid}	coverage={wildcards.cov}	reads={wildcards.n}	{wildcards.tag}	mutation_rate={params.mut_rate}	genome_size=$genome_size	requested_snps=$snps	requested_indels=$indels	true_node=$true_node" > {output.meta}
        '''

rule place_reads:
    input:
        panmap_bin="build/bin/panmap",
        panman=lambda wc: EXP_BY_ID[wc.eid]['panman_path'],
        index=lambda wc: f"{index_root(wc.pan_stem, EXP_BY_ID[wc.eid]['index_tag'])}/index.pmi",
        r1=f"{OUTPUT_DIR}/experiments/{{eid}}/{{pan_stem}}/{{tag}}/reads/cov{{cov}}_{{n}}_rep{{rep}}_R1.fastq",
        r2=f"{OUTPUT_DIR}/experiments/{{eid}}/{{pan_stem}}/{{tag}}/reads/cov{{cov}}_{{n}}_rep{{rep}}_R2.fastq",
        meta=f"{OUTPUT_DIR}/experiments/{{eid}}/{{pan_stem}}/{{tag}}/results/reads/cov{{cov}}_{{n}}_rep{{rep}}.txt"
    output:
        placement=f"{OUTPUT_DIR}/experiments/{{eid}}/{{pan_stem}}/{{tag}}/placements/reads/cov{{cov}}_{{n}}_rep{{rep}}/placements.tsv",
        detailed=f"{OUTPUT_DIR}/experiments/{{eid}}/{{pan_stem}}/{{tag}}/placements/reads/cov{{cov}}_{{n}}_rep{{rep}}/detailed.tsv",
        time_log=f"{OUTPUT_DIR}/experiments/{{eid}}/{{pan_stem}}/{{tag}}/placements/reads/cov{{cov}}_{{n}}_rep{{rep}}/time.log"
    params:
        outprefix=lambda wc: f"{OUTPUT_DIR}/experiments/{wc.eid}/{wc.pan_stem}/{wc.tag}/placements/reads/cov{wc.cov}_{wc.n}_rep{wc.rep}/result",
        k=lambda wc: EXP_BY_ID[wc.eid]['k'],
        s=lambda wc: EXP_BY_ID[wc.eid]['s'],
        l=lambda wc: EXP_BY_ID[wc.eid]['l']
    shell:
        r'''
        set -e; mkdir -p $(dirname {output.placement})
        true_node=$(grep -oP 'true_node=\K[^\t]+' {input.meta} || echo "unknown")
        /usr/bin/time -v {input.panmap_bin} {input.panman} {input.r1} {input.r2} --output {params.outprefix} -t 1 --stop place -i {input.index} 2> {output.time_log}
        mv {params.outprefix}.placement.tsv {output.placement}
        echo -e "experiment\treads\treplicate\ttrue_node\tmetric\tscore\tbest_node\tnode_id\ttimestamp\tk\ts\tl" > {output.detailed}
        tail -n +2 {output.placement} | while IFS=$'\t' read metric score hits nodes; do
            echo -e "{wildcards.eid}\t{wildcards.n}\t{wildcards.rep}\t$true_node\t$metric\t$score\t$nodes\t$nodes\t$(date '+%Y-%m-%d %H:%M:%S')\t{params.k}\t{params.s}\t{params.l}"
        done >> {output.detailed}
        '''

rule align_placement_accuracy:
    input:
        panmap_bin="build/bin/panmap",
        panman=lambda wc: EXP_BY_ID[wc.eid]['panman_path'],
        placement=f"{OUTPUT_DIR}/experiments/{{eid}}/{{pan_stem}}/{{tag}}/placements/reads/cov{{cov}}_{{n}}_rep{{rep}}/placements.tsv",
        meta=f"{OUTPUT_DIR}/experiments/{{eid}}/{{pan_stem}}/{{tag}}/results/reads/cov{{cov}}_{{n}}_rep{{rep}}.txt"
    output:
        accuracy_0bp=f"{OUTPUT_DIR}/experiments/{{eid}}/{{pan_stem}}/{{tag}}/alignments/reads/cov{{cov}}_{{n}}_rep{{rep}}/accuracy_0bp.tsv",
        accuracy_50bp=f"{OUTPUT_DIR}/experiments/{{eid}}/{{pan_stem}}/{{tag}}/alignments/reads/cov{{cov}}_{{n}}_rep{{rep}}/accuracy_50bp.tsv",
        accuracy_150bp=f"{OUTPUT_DIR}/experiments/{{eid}}/{{pan_stem}}/{{tag}}/alignments/reads/cov{{cov}}_{{n}}_rep{{rep}}/accuracy_150bp.tsv",
        accuracy_500bp=f"{OUTPUT_DIR}/experiments/{{eid}}/{{pan_stem}}/{{tag}}/alignments/reads/cov{{cov}}_{{n}}_rep{{rep}}/accuracy_500bp.tsv",
        alignment_results=f"{OUTPUT_DIR}/experiments/{{eid}}/{{pan_stem}}/{{tag}}/alignments/reads/cov{{cov}}_{{n}}_rep{{rep}}/accuracy.tsv",
        log_dir=directory(f"{OUTPUT_DIR}/experiments/{{eid}}/{{pan_stem}}/{{tag}}/alignments/reads/cov{{cov}}_{{n}}_rep{{rep}}/logs")
    run:
        log_dir = Path(output.log_dir)
        log_dir.mkdir(parents=True, exist_ok=True)
        
        # Read metadata
        true_node = "unknown"
        try:
            with open(input.meta) as f:
                fields = dict(p.split('=',1) for p in f.read().strip().split('\t') if '=' in p)
                true_node = fields.get("true_node", "unknown")
        except: pass
        
        # Read placements
        placement_nodes = {}
        try:
            with open(input.placement) as f:
                for line in f.readlines()[1:]:
                    parts = line.strip().split('\t')
                    if len(parts) >= 4 and parts[3]:
                        placement_nodes[parts[0]] = parts[3]
        except: pass
        
        exclusion_configs = [(0, output.accuracy_0bp), (50, output.accuracy_50bp), (150, output.accuracy_150bp), (500, output.accuracy_500bp)]
        header = "metric\ttrue_node\tplacement_node\tgenome_length\tsnps\tindels\ttotal_variants\talignment_length\tidentity\texcluded_bp\n"
        for _, f in exclusion_configs:
            with open(f, 'w') as out: out.write(header)
        
        if true_node == "unknown":
            for _, f in exclusion_configs:
                with open(f, 'a') as out:
                    for m in ["raw","jaccard","cosine","weighted_jaccard"]:
                        out.write(f"{m}\tunknown\t\t0\t0\t0\t0\t0\t0.0\t0\n")
            return
        
        for metric in ["raw", "jaccard", "cosine", "weighted_jaccard"]:
            placement_node = placement_nodes.get(metric, "")
            genome_length = snps = indels = 0
            identity = 1.0 if placement_node == true_node else 0.0
            
            if not placement_node or placement_node == true_node:
                # Get genome length from true node
                try:
                    true_fa = str(log_dir / f"{sanitize_filename(true_node)}.fa")
                    subprocess.run([input.panmap_bin, input.panman, "--dump-sequence", true_node, "-o", true_fa], capture_output=True, timeout=60)
                    with open(true_fa) as f:
                        genome_length = sum(len(l.strip()) for l in f if not l.startswith('>'))
                except: pass
                for excl, f in exclusion_configs:
                    with open(f, 'a') as out:
                        out.write(f"{metric}\t{true_node}\t{placement_node}\t{genome_length}\t0\t0\t0\t{genome_length}\t{identity:.6f}\t{excl}\n")
            else:
                # Align sequences
                try:
                    true_fa = str(log_dir / f"{sanitize_filename(true_node)}.fa")
                    place_fa = str(log_dir / f"{sanitize_filename(placement_node)}.fa")
                    subprocess.run([input.panmap_bin, input.panman, "--dump-sequence", true_node, "-o", true_fa], capture_output=True, timeout=60)
                    subprocess.run([input.panmap_bin, input.panman, "--dump-sequence", placement_node, "-o", place_fa], capture_output=True, timeout=60)
                    with open(true_fa) as f:
                        genome_length = sum(len(l.strip()) for l in f if not l.startswith('>'))
                    
                    paf = log_dir / f"{metric}.paf"
                    with open(paf, 'w') as pf:
                        subprocess.run(["minimap2", "-cx", "asm20", "--cs", true_fa, place_fa], stdout=pf, timeout=60)
                    
                    results_by_excl = {}
                    aln_len = 0
                    with open(paf) as f:
                        for line in f:
                            parts = line.strip().split('\t')
                            if len(parts) >= 12:
                                aln_len, matches = int(parts[10]), int(parts[9])
                                cs = next((p[5:] for p in parts[12:] if p.startswith('cs:Z:')), None)
                                if cs:
                                    for excl, _ in exclusion_configs:
                                        results_by_excl[excl] = count_variants_with_exclusion(cs, int(parts[7]), int(parts[8]), excl)
                                    identity = matches / aln_len if aln_len else 0
                                break
                    
                    for excl, f in exclusion_configs:
                        s, i = results_by_excl.get(excl, (0, 0))
                        with open(f, 'a') as out:
                            out.write(f"{metric}\t{true_node}\t{placement_node}\t{genome_length}\t{s}\t{i}\t{s+i}\t{aln_len}\t{identity:.6f}\t{excl}\n")
                except Exception as e:
                    print(f"[align] Error for {metric}: {e}")
                    for excl, f in exclusion_configs:
                        with open(f, 'a') as out:
                            out.write(f"{metric}\t{true_node}\t{placement_node}\t0\t0\t0\t0\t0\t0.0\t{excl}\n")
        
        # Create combined accuracy.tsv from individual files
        with open(output.alignment_results, 'w') as out:
            out.write(header)
            for excl, f in exclusion_configs:
                with open(f) as inp:
                    for line in inp.readlines()[1:]:  # Skip header
                        out.write(line)

# ============================================================================
# SUMMARY RULES
# ============================================================================
rule placements_summary:
    input: placements=PLACEMENTS, metas=META_FILES
    output: tsv=f"{OUTPUT_DIR}/reports/placements_summary.tsv"
    run:
        Path(output.tsv).parent.mkdir(parents=True, exist_ok=True)
        meta_map = {}
        for m in input.metas:
            try:
                with open(m) as f:
                    meta_map[m] = dict(p.split('=',1) for p in f.read().strip().split('\t') if '=' in p)
            except: meta_map[m] = {}
        
        with open(output.tsv, 'w') as out:
            out.write("id\tpan_stem\ttag\treads\treplicate\ttrue_node\ttop_score\tmutation_rate\tgenome_size\n")
            for (eid, pan, tag, cov, n, rep) in READS:
                m = f"{exp_root(eid, pan, tag)}/results/reads/cov{cov}_{n}_rep{rep}.txt"
                meta = meta_map.get(m, {})
                out.write(f"{eid}\t{pan}\t{tag}\t{n}\t{rep}\t{meta.get('true_node','')}\t0\t{meta.get('mutation_rate','0')}\t{meta.get('genome_size','0')}\n")

rule placements_by_metric:
    input: placements=PLACEMENTS, metas=META_FILES
    output: tsv=f"{OUTPUT_DIR}/reports/placements_by_metric.tsv"
    run:
        Path(output.tsv).parent.mkdir(parents=True, exist_ok=True)
        meta_map = {m: dict(p.split('=',1) for p in open(m).read().strip().split('\t') if '=' in p) if os.path.exists(m) else {} for m in input.metas}
        with open(output.tsv, 'w') as out:
            out.write("id\tpan_stem\ttag\treads\treplicate\ttrue_node\tmetric\tscore\tnodes\n")
            for (eid, pan, tag, cov, n, rep) in READS:
                p = f"{exp_root(eid, pan, tag)}/placements/reads/cov{cov}_{n}_rep{rep}/placements.tsv"
                m = f"{exp_root(eid, pan, tag)}/results/reads/cov{cov}_{n}_rep{rep}.txt"
                true_node = meta_map.get(m, {}).get('true_node', '')
                try:
                    with open(p) as f:
                        for line in f.readlines()[1:]:
                            parts = line.strip().split('\t')
                            if len(parts) >= 4:
                                out.write(f"{eid}\t{pan}\t{tag}\t{n}\t{rep}\t{true_node}\t{parts[0]}\t{parts[1]}\t{parts[3]}\n")
                except: pass

rule alignment_accuracy_summary:
    input: metas=META_FILES
    output: tsv=f"{OUTPUT_DIR}/reports/alignment_accuracy_summary.tsv"
    run:
        Path(output.tsv).parent.mkdir(parents=True, exist_ok=True)
        with open(output.tsv, 'w') as out:
            out.write("experiment_id\tpangenome\ttag\tcoverage\treads\treplicate\tmetric\ttrue_node\tplacement_node\tgenome_length\tsnps\tindels\ttotal_variants\talignment_length\tidentity\tmutation_rate\texcluded_bp\n")
            for (eid, pan, tag, cov, n, rep) in READS:
                meta = {}
                try:
                    with open(f"{exp_root(eid, pan, tag)}/results/reads/cov{cov}_{n}_rep{rep}.txt") as f:
                        meta = dict(p.split('=',1) for p in f.read().strip().split('\t') if '=' in p)
                except: pass
                mut_rate = meta.get("mutation_rate", 0)
                for excl in [0, 50, 150, 500]:
                    acc_file = f"{exp_root(eid, pan, tag)}/alignments/reads/cov{cov}_{n}_rep{rep}/accuracy_{excl}bp.tsv"
                    try:
                        with open(acc_file) as f:
                            for line in f.readlines()[1:]:
                                parts = line.strip().split('\t')
                                if len(parts) >= 10:
                                    out.write(f"{eid}\t{pan}\t{tag}\t{cov}\t{n}\t{rep}\t{parts[0]}\t{parts[1]}\t{parts[2]}\t{parts[3]}\t{parts[4]}\t{parts[5]}\t{parts[6]}\t{parts[7]}\t{parts[8]}\t{mut_rate}\t{parts[9]}\n")
                    except: pass

# ============================================================================
# PERFORMANCE SUMMARIES
# ============================================================================
rule index_performance_summary:
    input: logs=INDEX_TIME_LOGS
    output:
        summary=f"{OUTPUT_DIR}/reports/index_performance_summary.tsv",
        per_pan=expand(f"{OUTPUT_DIR}/reports/{{pan}}/index_performance_summary.tsv", pan=PANGENOMES.keys())
    run:
        Path(output.summary).parent.mkdir(parents=True, exist_ok=True)
        with open(output.summary, 'w') as out:
            out.write('experiment_id\tpangenome\ttag\tk\ts\tl\twall_time_sec\tuser_time_sec\tsystem_time_sec\tmax_rss_mb\n')
            for e in EXPERIMENTS:
                stats = parse_time_log(f"{exp_root(e['id'], e['pan_stem'], e['tag'])}/indexes/index_time.log")
                if 'wall_time_sec' in stats:
                    out.write(f"{e['id']}\t{e['pan_stem']}\t{e['tag']}\t{e['k']}\t{e['s']}\t{e['l']}\t{stats.get('wall_time_sec',0)}\t{stats.get('user_time_sec',0)}\t{stats.get('system_time_sec',0)}\t{stats.get('max_rss_mb',0)}\n")
        for pan in PANGENOMES:
            Path(f"{OUTPUT_DIR}/reports/{pan}").mkdir(parents=True, exist_ok=True)
            with open(f"{OUTPUT_DIR}/reports/{pan}/index_performance_summary.tsv", 'w') as out:
                out.write('experiment_id\tpangenome\ttag\tk\ts\tl\twall_time_sec\tmax_rss_mb\n')
                for e in EXPERIMENTS:
                    if e['pan_stem'] == pan:
                        stats = parse_time_log(f"{exp_root(e['id'], e['pan_stem'], e['tag'])}/indexes/index_time.log")
                        if 'wall_time_sec' in stats:
                            out.write(f"{e['id']}\t{pan}\t{e['tag']}\t{e['k']}\t{e['s']}\t{e['l']}\t{stats.get('wall_time_sec',0)}\t{stats.get('max_rss_mb',0)}\n")

rule placement_performance_summary:
    input: logs=PLACEMENT_TIME_LOGS
    output:
        summary=f"{OUTPUT_DIR}/reports/placement_performance_summary.tsv",
        per_pan=expand(f"{OUTPUT_DIR}/reports/{{pan}}/placement_performance_summary.tsv", pan=PANGENOMES.keys())
    run:
        Path(output.summary).parent.mkdir(parents=True, exist_ok=True)
        with open(output.summary, 'w') as out:
            out.write('experiment_id\tpangenome\ttag\tcoverage\treads\treplicate\tk\ts\tl\twall_time_sec\tmax_rss_mb\n')
            for (eid, pan, tag, cov, n, rep) in READS:
                e = EXP_BY_ID[eid]
                stats = parse_time_log(f"{exp_root(eid, pan, tag)}/placements/reads/cov{cov}_{n}_rep{rep}/time.log")
                if 'wall_time_sec' in stats:
                    out.write(f"{eid}\t{pan}\t{tag}\t{cov}\t{n}\t{rep}\t{e['k']}\t{e['s']}\t{e['l']}\t{stats.get('wall_time_sec',0)}\t{stats.get('max_rss_mb',0)}\n")
        for pan_name in PANGENOMES:
            Path(f"{OUTPUT_DIR}/reports/{pan_name}").mkdir(parents=True, exist_ok=True)
            with open(f"{OUTPUT_DIR}/reports/{pan_name}/placement_performance_summary.tsv", 'w') as out:
                out.write('experiment_id\tpangenome\ttag\tcoverage\treads\treplicate\tk\ts\tl\twall_time_sec\tmax_rss_mb\n')
                for (eid, pan, tag, cov, n, rep) in READS:
                    if pan == pan_name:
                        e = EXP_BY_ID[eid]
                        stats = parse_time_log(f"{exp_root(eid, pan, tag)}/placements/reads/cov{cov}_{n}_rep{rep}/time.log")
                        if 'wall_time_sec' in stats:
                            out.write(f"{eid}\t{pan}\t{tag}\t{cov}\t{n}\t{rep}\t{e['k']}\t{e['s']}\t{e['l']}\t{stats.get('wall_time_sec',0)}\t{stats.get('max_rss_mb',0)}\n")

# ============================================================================
# PLOTTING
# ============================================================================
rule plot_performance:
    input:
        index_perf=f"{OUTPUT_DIR}/reports/index_performance_summary.tsv",
        placement_perf=f"{OUTPUT_DIR}/reports/placement_performance_summary.tsv"
    output:
        # Per-pangenome performance plots
        expand(f"{OUTPUT_DIR}/plots/{{pan}}/performance/index_time_by_k.png", pan=PANGENOMES.keys()),
        expand(f"{OUTPUT_DIR}/plots/{{pan}}/performance/index_memory_by_k.png", pan=PANGENOMES.keys()),
        expand(f"{OUTPUT_DIR}/plots/{{pan}}/performance/placement_time_by_k.png", pan=PANGENOMES.keys()),
        expand(f"{OUTPUT_DIR}/plots/{{pan}}/performance/placement_memory_by_k.png", pan=PANGENOMES.keys()),
        expand(f"{OUTPUT_DIR}/plots/{{pan}}/performance/index_time_by_k_l1.png", pan=PANGENOMES.keys()),
        expand(f"{OUTPUT_DIR}/plots/{{pan}}/performance/index_time_by_k_l3.png", pan=PANGENOMES.keys()),
        expand(f"{OUTPUT_DIR}/plots/{{pan}}/performance/placement_time_by_k_l1.png", pan=PANGENOMES.keys()),
        expand(f"{OUTPUT_DIR}/plots/{{pan}}/performance/placement_time_by_k_l3.png", pan=PANGENOMES.keys())
    run:
        import pandas as pd
        import matplotlib.pyplot as plt
        import seaborn as sns
        
        index_df = pd.read_csv(input.index_perf, sep='\t')
        placement_df = pd.read_csv(input.placement_perf, sep='\t')
        
        def boxplot(data, x, y, title, ylabel, outfile, hue=None):
            plt.figure(figsize=(10, 6))
            if len(data) == 0:
                plt.text(0.5, 0.5, 'No data', ha='center', va='center', transform=plt.gca().transAxes)
            else:
                sns.boxplot(data=data, x=x, y=y, hue=hue, palette='Set2')
            plt.xlabel('k-mer size (k)'); plt.ylabel(ylabel); plt.title(title)
            plt.tight_layout(); plt.savefig(outfile, dpi=300); plt.close()
        
        # Generate per-pangenome plots
        for pan in PANGENOMES:
            out_dir = Path(f"{OUTPUT_DIR}/plots/{pan}/performance")
            out_dir.mkdir(parents=True, exist_ok=True)
            
            pan_index = index_df[index_df['pangenome']==pan] if len(index_df) else index_df
            pan_placement = placement_df[placement_df['pangenome']==pan] if len(placement_df) else placement_df
            
            boxplot(pan_index, 'k', 'wall_time_sec', f'Index Build Time - {pan}', 'Wall time (s)', out_dir / 'index_time_by_k.png', 'l')
            boxplot(pan_index, 'k', 'max_rss_mb', f'Index Memory - {pan}', 'Peak memory (MB)', out_dir / 'index_memory_by_k.png', 'l')
            boxplot(pan_index[pan_index['l']==1] if len(pan_index) else pan_index, 'k', 'wall_time_sec', f'Index Time (l=1) - {pan}', 'Wall time (s)', out_dir / 'index_time_by_k_l1.png')
            boxplot(pan_index[pan_index['l']==3] if len(pan_index) else pan_index, 'k', 'wall_time_sec', f'Index Time (l=3) - {pan}', 'Wall time (s)', out_dir / 'index_time_by_k_l3.png')
            boxplot(pan_placement, 'k', 'wall_time_sec', f'Placement Time - {pan}', 'Wall time (s)', out_dir / 'placement_time_by_k.png', 'l')
            boxplot(pan_placement, 'k', 'max_rss_mb', f'Placement Memory - {pan}', 'Peak memory (MB)', out_dir / 'placement_memory_by_k.png', 'l')
            boxplot(pan_placement[pan_placement['l']==1] if len(pan_placement) else pan_placement, 'k', 'wall_time_sec', f'Placement Time (l=1) - {pan}', 'Wall time (s)', out_dir / 'placement_time_by_k_l1.png')
            boxplot(pan_placement[pan_placement['l']==3] if len(pan_placement) else pan_placement, 'k', 'wall_time_sec', f'Placement Time (l=3) - {pan}', 'Wall time (s)', out_dir / 'placement_time_by_k_l3.png')

rule plot_placement_accuracy:
    input: accuracy_summary=f"{OUTPUT_DIR}/reports/alignment_accuracy_summary.tsv"
    output: done=touch(f"{OUTPUT_DIR}/plots/placement_accuracy_plots.done")
    run:
        import pandas as pd
        import matplotlib.pyplot as plt
        import numpy as np
        
        if not Path(input.accuracy_summary).exists():
            return
        
        df = pd.read_csv(input.accuracy_summary, sep='\t')
        if len(df) == 0:
            return
        
        def extract_ksl(tag):
            k = int(re.search(r'k(\d+)', tag).group(1)) if re.search(r'k(\d+)', tag) else 0
            s = int(re.search(r's(\d+)', tag).group(1)) if re.search(r's(\d+)', tag) else 0
            l = int(re.search(r'l(\d+)', tag).group(1)) if re.search(r'l(\d+)', tag) else 1
            return k, s, l
        
        df[['k', 's', 'l']] = df['tag'].apply(lambda x: pd.Series(extract_ksl(x)))
        df['distance'] = df['total_variants']
        df['k_s_l'] = df.apply(lambda r: f"k{r['k']}_s{r['s']}_l{r['l']}", axis=1)
        df['category'] = df['reads'].astype(str) + ' reads, μ=' + df['mutation_rate'].astype(str)
        
        # Generate plots per pangenome
        for pan in df['pangenome'].unique():
            pan_df = df[df['pangenome'] == pan]
            plots_dir = Path(f"{OUTPUT_DIR}/plots/{pan}/accuracy")
            plots_dir.mkdir(parents=True, exist_ok=True)
            
            # Save per-pangenome summary
            pan_df.to_csv(plots_dir / "accuracy_data.tsv", sep='\t', index=False)
            
            # Create plots per k_s_l and metric for this pangenome
            for ksl in pan_df['k_s_l'].unique():
                ksl_data = pan_df[pan_df['k_s_l'] == ksl]
                for metric in ksl_data['metric'].unique():
                    metric_data = ksl_data[ksl_data['metric'] == metric]
                    if len(metric_data) == 0: continue
                    
                    categories = sorted(metric_data['category'].unique())
                    fig, ax = plt.subplots(figsize=(12, 8))
                    
                    distance_bins = list(range(6)) + ['>5']
                    plot_data = []
                    for cat in categories:
                        cat_data = metric_data[metric_data['category'] == cat]
                        total = len(cat_data)
                        if total == 0: continue
                        counts = [len(cat_data[cat_data['distance']==d])/total for d in range(6)]
                        counts.append(len(cat_data[cat_data['distance']>5])/total)
                        plot_data.append(counts)
                    
                    x = np.arange(len(distance_bins))
                    width = 0.8 / max(len(categories), 1)
                    colors = plt.cm.viridis(np.linspace(1, 0, len(categories)))
                    
                    for i, (cat, data) in enumerate(zip(categories, plot_data)):
                        ax.bar(x + (i - len(categories)/2 + 0.5) * width, data, width, label=cat, color=colors[i], alpha=0.8)
                    
                    ax.set_xlabel('Distance from expected (SNPs + Indels)')
                    ax.set_ylabel('Proportion')
                    ax.set_title(f'{pan}: Placement accuracy, {ksl}, {metric}')
                    ax.set_xticks(x)
                    ax.set_xticklabels([str(d) for d in distance_bins])
                    ax.set_ylim(0, 1)
                    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
                    ax.grid(alpha=0.3)
                    
                    plt.tight_layout()
                    plt.savefig(plots_dir / f"placement_accuracy_{ksl}_{metric}.png", dpi=300, bbox_inches='tight')
                    plt.close()

rule plot_accuracy_by_k:
    input: summary=f"{OUTPUT_DIR}/reports/alignment_accuracy_summary.tsv"
    output:
        expand(f"{OUTPUT_DIR}/plots/{{pan}}/k_analysis/accuracy_by_k_l1.png", pan=PANGENOMES.keys()),
        expand(f"{OUTPUT_DIR}/plots/{{pan}}/k_analysis/accuracy_by_k_l1.pdf", pan=PANGENOMES.keys()),
        expand(f"{OUTPUT_DIR}/plots/{{pan}}/k_analysis/accuracy_by_k_l3.png", pan=PANGENOMES.keys()),
        expand(f"{OUTPUT_DIR}/plots/{{pan}}/k_analysis/accuracy_by_k_l3.pdf", pan=PANGENOMES.keys())
    run:
        import pandas as pd
        import matplotlib.pyplot as plt
        import seaborn as sns
        
        df = pd.read_csv(input.summary, sep='\t')
        
        if 'total_variants' not in df.columns or len(df) == 0:
            for pan in PANGENOMES:
                out_dir = Path(f"{OUTPUT_DIR}/plots/{pan}/k_analysis")
                out_dir.mkdir(parents=True, exist_ok=True)
                for f in ['accuracy_by_k_l1.png', 'accuracy_by_k_l1.pdf', 'accuracy_by_k_l3.png', 'accuracy_by_k_l3.pdf']:
                    Path(out_dir / f).touch()
            return
        
        df['correct'] = (df['total_variants'] == 0).astype(int)
        df['k'] = df['tag'].apply(lambda t: int(re.search(r'k(\d+)', t).group(1)) if re.search(r'k(\d+)', t) else 0)
        df['l'] = df['tag'].apply(lambda t: int(re.search(r'l(\d+)', t).group(1)) if re.search(r'l(\d+)', t) else 1)
        
        agg = df.groupby(['pangenome','k','l','metric','reads'], as_index=False).agg({'correct':'mean'})
        agg['accuracy'] = agg['correct'] * 100
        
        def plot_l(pan_data, pan, l_val, png, pdf):
            sub = pan_data[pan_data['l'] == l_val]
            if sub.empty:
                Path(png).touch(); Path(pdf).touch()
                return
            plt.figure(figsize=(10, 6))
            sns.set_style('whitegrid')
            for cat in sorted(sub.apply(lambda r: f"{r['metric']}:{r['reads']}", axis=1).unique()):
                m, r = cat.split(':')
                d = sub[(sub['metric']==m) & (sub['reads']==int(r))].sort_values('k')
                if not d.empty:
                    plt.plot(d['k'], d['accuracy'], marker='o', label=f"{m} ({r} reads)")
            plt.xlabel('k'); plt.ylabel('Accuracy (%)'); plt.title(f'{pan}: Accuracy vs k (l={l_val})')
            plt.xticks(sorted(sub['k'].unique())); plt.ylim(0, 100)
            plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
            plt.tight_layout(); plt.savefig(png, dpi=300); plt.savefig(pdf); plt.close()
        
        # Generate per-pangenome plots
        for pan in PANGENOMES:
            out_dir = Path(f"{OUTPUT_DIR}/plots/{pan}/k_analysis")
            out_dir.mkdir(parents=True, exist_ok=True)
            pan_agg = agg[agg['pangenome'] == pan]
            plot_l(pan_agg, pan, 1, out_dir / 'accuracy_by_k_l1.png', out_dir / 'accuracy_by_k_l1.pdf')
            plot_l(pan_agg, pan, 3, out_dir / 'accuracy_by_k_l3.png', out_dir / 'accuracy_by_k_l3.pdf')

# ============================================================================
# COMPREHENSIVE STATS TABLES
# ============================================================================

def get_file_size_mb(path):
    """Get file size in MB, returns 0 if file doesn't exist."""
    try:
        return os.path.getsize(path) / (1024 * 1024)
    except:
        return 0

def get_gzip_size_mb(path):
    """Get gzip-compressed size in MB without modifying the file."""
    import subprocess
    try:
        result = subprocess.run(
            ['gzip', '-c', path],
            capture_output=True,
            timeout=300
        )
        return len(result.stdout) / (1024 * 1024)
    except:
        return 0

rule index_stats:
    """Comprehensive index statistics including file sizes."""
    input: 
        indexes=INDEX_FILES,
        time_logs=INDEX_TIME_LOGS
    output:
        tsv=f"{OUTPUT_DIR}/reports/index_stats.tsv"
    run:
        Path(output.tsv).parent.mkdir(parents=True, exist_ok=True)
        with open(output.tsv, 'w') as out:
            out.write('experiment_id\tpangenome\ttag\tk\ts\tl\tindex_size_mb\tindex_size_gzip_mb\twall_time_sec\tuser_time_sec\tsystem_time_sec\tmax_rss_mb\tpct_cpu\n')
            for e in EXPERIMENTS:
                exp_dir = exp_root(e['id'], e['pan_stem'], e['tag'])
                index_path = f"{exp_dir}/indexes/index.pmi"
                time_log = f"{exp_dir}/indexes/index_time.log"
                
                size_mb = get_file_size_mb(index_path)
                size_gzip_mb = get_gzip_size_mb(index_path)
                stats = parse_time_log(time_log)
                
                out.write(f"{e['id']}\t{e['pan_stem']}\t{e['tag']}\t{e['k']}\t{e['s']}\t{e['l']}\t")
                out.write(f"{size_mb:.2f}\t{size_gzip_mb:.2f}\t{stats.get('wall_time_sec',0)}\t{stats.get('user_time_sec',0)}\t")
                out.write(f"{stats.get('system_time_sec',0)}\t{stats.get('max_rss_mb',0)}\t{stats.get('pct_cpu','0')}\n")

rule placement_stats:
    """Comprehensive placement statistics with accuracy."""
    input:
        placements=PLACEMENTS,
        time_logs=PLACEMENT_TIME_LOGS,
        metas=META_FILES,
        accuracy=ACCURACY
    output:
        tsv=f"{OUTPUT_DIR}/reports/placement_stats.tsv"
    run:
        Path(output.tsv).parent.mkdir(parents=True, exist_ok=True)
        
        # Build meta lookup
        meta_map = {}
        for m in input.metas:
            if os.path.exists(m):
                try:
                    meta_map[m] = dict(p.split('=',1) for p in open(m).read().strip().split('\t') if '=' in p)
                except:
                    meta_map[m] = {}
        
        # Build accuracy lookup (from accuracy files)
        acc_map = {}
        for acc_file in input.accuracy:
            if os.path.exists(acc_file):
                try:
                    with open(acc_file) as f:
                        for line in f.readlines()[1:]:
                            parts = line.strip().split('\t')
                            if len(parts) >= 7:
                                # key = (file, metric)
                                key = (acc_file, parts[0])
                                acc_map[key] = {
                                    'snps': parts[4] if len(parts) > 4 else '0',
                                    'indels': parts[5] if len(parts) > 5 else '0',
                                    'total_variants': parts[6] if len(parts) > 6 else '0'
                                }
                except:
                    pass
        
        with open(output.tsv, 'w') as out:
            out.write('experiment_id\tpangenome\ttag\tcoverage\treads\treplicate\tk\ts\tl\t')
            out.write('wall_time_sec\tmax_rss_mb\ttrue_node\tmutation_rate\t')
            out.write('jaccard_distance\tweighted_distance\tcosine_distance\n')
            
            for (eid, pan, tag, cov, n, rep) in READS:
                e = EXP_BY_ID[eid]
                exp_dir = exp_root(eid, pan, tag)
                time_log = f"{exp_dir}/placements/reads/cov{cov}_{n}_rep{rep}/time.log"
                meta_file = f"{exp_dir}/results/reads/cov{cov}_{n}_rep{rep}.txt"
                acc_file = f"{exp_dir}/alignments/reads/cov{cov}_{n}_rep{rep}/accuracy_0bp.tsv"
                
                stats = parse_time_log(time_log)
                meta = meta_map.get(meta_file, {})
                
                # Get distances for each metric
                jacc_dist = acc_map.get((acc_file, 'jaccard'), {}).get('total_variants', '')
                weight_dist = acc_map.get((acc_file, 'weighted_jaccard'), {}).get('total_variants', '')
                cos_dist = acc_map.get((acc_file, 'cosine'), {}).get('total_variants', '')
                
                out.write(f"{eid}\t{pan}\t{tag}\t{cov}\t{n}\t{rep}\t{e['k']}\t{e['s']}\t{e['l']}\t")
                out.write(f"{stats.get('wall_time_sec',0)}\t{stats.get('max_rss_mb',0)}\t")
                out.write(f"{meta.get('true_node','')}\t{meta.get('mutation_rate','')}\t")
                out.write(f"{jacc_dist}\t{weight_dist}\t{cos_dist}\n")

rule combined_results:
    """Combined table with all index and placement data for analysis."""
    input:
        index_stats=f"{OUTPUT_DIR}/reports/index_stats.tsv",
        placement_stats=f"{OUTPUT_DIR}/reports/placement_stats.tsv"
    output:
        tsv=f"{OUTPUT_DIR}/reports/combined_results.tsv"
    run:
        import pandas as pd
        Path(output.tsv).parent.mkdir(parents=True, exist_ok=True)
        
        # Read index stats
        if os.path.exists(input.index_stats):
            idx_df = pd.read_csv(input.index_stats, sep='\t')
            idx_df = idx_df.rename(columns={
                'wall_time_sec': 'index_wall_time_sec',
                'user_time_sec': 'index_user_time_sec',
                'system_time_sec': 'index_system_time_sec',
                'max_rss_mb': 'index_max_rss_mb',
                'pct_cpu': 'index_pct_cpu',
                'index_size_gzip_mb': 'index_size_gzip_mb'
            })
        else:
            idx_df = pd.DataFrame()
        
        # Read placement stats
        if os.path.exists(input.placement_stats):
            place_df = pd.read_csv(input.placement_stats, sep='\t')
            place_df = place_df.rename(columns={
                'wall_time_sec': 'placement_wall_time_sec',
                'max_rss_mb': 'placement_max_rss_mb'
            })
        else:
            place_df = pd.DataFrame()
        
        if idx_df.empty and place_df.empty:
            pd.DataFrame().to_csv(output.tsv, sep='\t', index=False)
            return
        
        if idx_df.empty:
            place_df.to_csv(output.tsv, sep='\t', index=False)
            return
            
        if place_df.empty:
            idx_df.to_csv(output.tsv, sep='\t', index=False)
            return
        
        # Merge on experiment_id, pangenome, tag
        merged = pd.merge(
            place_df,
            idx_df[['experiment_id', 'pangenome', 'tag', 'index_size_mb', 'index_size_gzip_mb',
                    'index_wall_time_sec', 'index_user_time_sec', 
                    'index_system_time_sec', 'index_max_rss_mb', 'index_pct_cpu']],
            on=['experiment_id', 'pangenome', 'tag'],
            how='left'
        )
        
        # Reorder columns for readability
        col_order = [
            'experiment_id', 'pangenome', 'tag', 'k', 's', 'l',
            'coverage', 'reads', 'replicate', 'mutation_rate', 'true_node',
            'index_size_mb', 'index_size_gzip_mb', 'index_wall_time_sec', 'index_max_rss_mb',
            'placement_wall_time_sec', 'placement_max_rss_mb',
            'jaccard_distance', 'weighted_distance', 'cosine_distance',
            'index_user_time_sec', 'index_system_time_sec', 'index_pct_cpu'
        ]
        # Only include columns that exist
        col_order = [c for c in col_order if c in merged.columns]
        merged = merged[col_order]
        
        merged.to_csv(output.tsv, sep='\t', index=False)

# ============================================================================
# UTILITY RULES
# ============================================================================
rule clean:
    run:
        import shutil, glob
        for d in [Path(OUTPUT_DIR), Path('.snakemake')]:
            if d.exists(): shutil.rmtree(d, ignore_errors=True)
        for e in EXPERIMENTS:
            pan = Path(e['panman_path']).resolve()
            for pat in [f"{pan}.placement.tsv", f"{pan}.random.*.fa", f"{pan}.pmi", str(pan.parent / 'read_seeds_debug.log')]:
                for f in glob.glob(pat):
                    try: os.remove(f)
                    except: pass
