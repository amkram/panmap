OUTPUT_DIR = "workflow_output"

# Available pangenome datasets with genome size estimates
PANGENOMES = {
    'rsv_4K': {
        'path': 'panmans/rsv_4K.panman',
        'genome_size': 15000  # RSV genome ~15kb
    },
    'sars_20K': {
        'path': 'panmans/sars_20000.panman', 
        'genome_size': 30000  # SARS-CoV-2 genome ~30kb
    },
    'tb_400': {
        'path': 'panmans/tb_400.panman',
        'genome_size': 4000000  # Mycobacterium tuberculosis genome ~4Mb
    }
}

# Experiment parameters - modify these to customize experiments
CONFIG = {
    'pangenomes': ['rsv_4K'],   # Which datasets to test - multiple pangenomes
    'k_values': [18,31],               # k-mer sizes - test different sensitivity levels (reduced from 4 to 2)
    's_values': [8],                    # Minimizer spacing - test density vs speed tradeoff (reduced to 1)
    'include_internal': [False],            # Include internal nodes? [False, True]
    'coverage_levels': [10,100],        # Coverage levels (X coverage) (removed 1x)
    'mutation_rates': [0.0001], # Mutation rates per base pair (reduced to 1)
    'replicates': 10,                        # Replicates per condition (reduced from 100 to 10)
    'model': 'NovaSeq',                     # Sequencing model
    'read_length': 150                      # Average read length for coverage calculation
}

def generate_experiments():
    """Generate all parameter combinations"""
    experiments = []
    exp_id = 0
    
    for pangenome in CONFIG['pangenomes']:
        for k in CONFIG['k_values']:
            for s in CONFIG['s_values']:
                for include_internal in CONFIG['include_internal']:
                    for mutation_rate in CONFIG['mutation_rates']:
                        tag = f"k{k}_s{s}_{'int' if include_internal else 'noint'}"
                        if mutation_rate > 0:
                            tag += f"_mut{mutation_rate}"
                        
                        # Calculate read counts for each coverage level
                        genome_size = PANGENOMES[pangenome]['genome_size']
                        read_length = CONFIG['read_length']
                        read_counts = []
                        for coverage in CONFIG['coverage_levels']:
                            num_reads = int((genome_size * coverage) / read_length)
                            read_counts.append(num_reads)
                        
                        experiments.append({
                            'id': f'exp{exp_id}',
                            'pangenome_name': pangenome,
                            'panman_path': PANGENOMES[pangenome]['path'],
                            'genome_size': genome_size,
                            'pan_stem': pangenome,
                            'k': k, 's': s,
                            'include_internal': include_internal,
                            'mutation_rate': mutation_rate,
                            'model': CONFIG['model'],
                            'coverage_levels': CONFIG['coverage_levels'],
                            'num_reads_values': read_counts,
                            'replicates': CONFIG['replicates'],
                            'tag': tag
                        })
                        exp_id += 1
    
    return experiments

EXPERIMENTS = generate_experiments()
EXP_BY_ID = {exp['id']: exp for exp in EXPERIMENTS}

def _exp_root(eid, pan_stem, tag):
    return f"{OUTPUT_DIR}/experiments/{eid}/{pan_stem}/{tag}"

def _tag(exp):
    return exp['tag']

# Print experiment summary
print(f"[config] Generated {len(EXPERIMENTS)} experiments")
print(f"[config] Coverage levels: {CONFIG['coverage_levels']} with read length {CONFIG['read_length']}")
print(f"[config] Read counts per coverage:")
genome_size = list(PANGENOMES.values())[0]['genome_size']  # Use first genome size as reference
for cov in CONFIG['coverage_levels']:
    reads = int((genome_size * cov) / CONFIG['read_length'])
    print(f"  {cov}x coverage: {reads} reads")
print(f"[config] Experiments:")
for exp in EXPERIMENTS[:3]:  # Show first 3
    print(f"  {exp['id']}: {exp['pangenome_name']} - {exp['tag']} (mutation_rate: {exp['mutation_rate']})")
if len(EXPERIMENTS) > 3:
    print(f"  ... and {len(EXPERIMENTS)-3} more")

import pathlib, os, re

# =============================================================================
# GENERATED FILE PATHS
# =============================================================================

READS = [(e['id'], e['pan_stem'], e['tag'], cov, n, rep) 
         for e in EXPERIMENTS 
         for cov, n in zip(e['coverage_levels'], e['num_reads_values'])
         for rep in range(e['replicates'])]

FASTQS = {
    'r1': [f"{_exp_root(eid, pan_stem, tag)}/reads/cov{cov}_{n}_rep{rep}_R1.fastq" for (eid, pan_stem, tag, cov, n, rep) in READS],
    'r2': [f"{_exp_root(eid, pan_stem, tag)}/reads/cov{cov}_{n}_rep{rep}_R2.fastq" for (eid, pan_stem, tag, cov, n, rep) in READS]
}
RESULT_FILES = [f"{_exp_root(eid, pan_stem, tag)}/results/reads/cov{cov}_{n}_rep{rep}.txt" for (eid, pan_stem, tag, cov, n, rep) in READS]
INDEX_FILES = [f"{_exp_root(e['id'], e['pangenome_name'], e['tag'])}/indexes/index.pmi" for e in EXPERIMENTS]
MM_FILES = [f"{_exp_root(e['id'], e['pangenome_name'], e['tag'])}/indexes/{e['pangenome_name']}.panman.mm" for e in EXPERIMENTS]
GENOME_FASTA = [f"{_exp_root(eid, pan_stem, tag)}/genomes/reads/cov{cov}_{n}_rep{rep}/random_node.fasta" for (eid, pan_stem, tag, cov, n, rep) in READS]
MUT_FASTA = [f"{_exp_root(eid, pan_stem, tag)}/mutgenomes/reads/cov{cov}_{n}_rep{rep}/mutated.fasta" for (eid, pan_stem, tag, cov, n, rep) in READS]
PLACEMENTS = [f"{_exp_root(eid, pan_stem, tag)}/placements/reads/cov{cov}_{n}_rep{rep}/placements.tsv" for (eid, pan_stem, tag, cov, n, rep) in READS]
DETAILED_PLACEMENTS = [f"{_exp_root(eid, pan_stem, tag)}/placements/reads/cov{cov}_{n}_rep{rep}/detailed.tsv" for (eid, pan_stem, tag, cov, n, rep) in READS]
PLACEMENT_LOGS = [f"{_exp_root(eid, pan_stem, tag)}/placements/reads/cov{cov}_{n}_rep{rep}/panmap.log" for (eid, pan_stem, tag, cov, n, rep) in READS]
ALIGNMENT_ACCURACY = [f"{_exp_root(eid, pan_stem, tag)}/alignments/reads/cov{cov}_{n}_rep{rep}/accuracy.tsv" for (eid, pan_stem, tag, cov, n, rep) in READS]
AGGREGATION_MARKERS = [f"{OUTPUT_DIR}/results/aggregated/{eid}/{pan_stem}/{tag}/cov{cov}_{n}_rep{rep}.done" for (eid, pan_stem, tag, cov, n, rep) in READS]

rule all:
    input:
        f"{OUTPUT_DIR}/reports/experiments_summary.tsv",
        f"{OUTPUT_DIR}/reports/placements_summary.tsv",
        f"{OUTPUT_DIR}/reports/placements_by_metric.tsv",
        f"{OUTPUT_DIR}/reports/alignment_accuracy_summary.tsv",
        f"{OUTPUT_DIR}/plots/placement_accuracy_plots.done",
        *GENOME_FASTA,
        *MUT_FASTA,
        *INDEX_FILES,
        *MM_FILES,
        *FASTQS['r1'],
        *FASTQS['r2'],
        *RESULT_FILES,
        *PLACEMENTS,
        *DETAILED_PLACEMENTS,
        *ALIGNMENT_ACCURACY,
        *AGGREGATION_MARKERS

rule experiments_summary:
    output:
        tsv=f"{OUTPUT_DIR}/reports/experiments_summary.tsv"
    run:
        # Ensure output directory's reports path exists
        pathlib.Path(output.tsv).parent.mkdir(parents=True, exist_ok=True)
        header = ['id','tag','pan_stem','panman_path','include_internal','k','s','mutation_rate','genome_size','num_reads_values','replicates','model']
        with open(output.tsv, 'w') as tf:
            tf.write('\t'.join(header)+'\n')
            for e in EXPERIMENTS:
                tf.write('\t'.join([
                    e['id'], e['tag'], e['pan_stem'], str(e['panman_path']), str(e['include_internal']).lower(), str(e['k']), str(e['s']),
                    str(e['mutation_rate']), str(e['genome_size']), ';'.join(map(str,e['num_reads_values'])), str(e['replicates']), e['model']
                ])+'\n')

rule dump_random_node:
    input:
        panmap_bin="build/bin/panmap",
        panman=lambda wc: EXP_BY_ID[wc.eid]['panman_path']
    output:
        f"{OUTPUT_DIR}/experiments/{{eid}}/{{pan_stem}}/{{tag}}/genomes/reads/cov{{cov}}_{{n}}_rep{{rep}}/random_node.fasta"
    run:
        ex = EXP_BY_ID[wildcards.eid]
        print(f"[dump_random_node] eid={wildcards.eid} pan_stem={wildcards.pan_stem} cov={wildcards.cov} n={wildcards.n} rep={wildcards.rep}")
        outp = pathlib.Path(output[0]).resolve()  # Use absolute path
        print(f"[dump_random_node] output path: {outp}")
        outp.parent.mkdir(parents=True, exist_ok=True)
        
        print(f"[dump_random_node] running panmap --dump-random-node for {ex['pan_stem']} replicate {wildcards.rep}")
        import glob, shutil, os, tempfile, time, hashlib
        
        # Create deterministic seed based on experiment ID AND replicate to ensure different random node per replicate
        seed_string = f"{ex['id']}_random_node_cov{wildcards.cov}_{wildcards.n}_rep{wildcards.rep}"
        seed_hash = hashlib.md5(seed_string.encode()).hexdigest()[:8]
        seed_int = int(seed_hash, 16) % (2**31 - 1)
        print(f"[dump_random_node] using deterministic seed: {seed_int} for experiment {ex['id']} replicate {wildcards.rep}")
        
        # Use a temporary working directory to avoid race conditions
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = pathlib.Path(tmpdir)
            panmap_bin_abs = str(pathlib.Path(input.panmap_bin).resolve())
            panman_abs = str(pathlib.Path(input.panman).resolve())
            
            # Copy panman file to temp directory to avoid conflicts
            temp_panman = tmpdir_path / f"{ex['pan_stem']}_{wildcards.eid}.panman"
            shutil.copy2(panman_abs, temp_panman)
            
            # Run panmap in temp directory with deterministic seed
            original_cwd = os.getcwd()
            try:
                os.chdir(tmpdir)
                result = shell(f"{panmap_bin_abs} {temp_panman} --dump-random-node --seed {seed_int}", check=False)
                
                # Find generated random node fasta in temp directory
                candidates = list(tmpdir_path.glob("*.random.*.fa"))
                
                if candidates:
                    # Sort by modification time and take the newest
                    newest = max(candidates, key=lambda x: x.stat().st_mtime)
                    print(f"[dump_random_node] copying {newest} -> {outp}")
                    print(f"[dump_random_node] source file size: {newest.stat().st_size} bytes")
                    
                    # Ensure output directory exists
                    outp.parent.mkdir(parents=True, exist_ok=True)
                    
                    # Remove existing target if present
                    if outp.exists() or outp.is_symlink():
                        try:
                            outp.unlink()
                        except Exception:
                            pass
                    
                    # Copy to final location and verify
                    shutil.copyfile(newest, outp)
                    
                    # Force filesystem sync
                    import os
                    try:
                        if hasattr(os, 'fsync'):
                            with open(outp, 'r+b') as f:
                                os.fsync(f.fileno())
                    except Exception:
                        pass
                    
                    # Verify the file was copied successfully
                    if outp.exists() and outp.stat().st_size > 0:
                        print(f"[dump_random_node] successfully copied to {outp}, size: {outp.stat().st_size} bytes")
                    else:
                        print(f"[dump_random_node] copy failed - output file doesn't exist or is empty")
                        raise RuntimeError(f"Failed to copy random node sequence for {ex['pan_stem']}")
                        
                else:
                    print(f"[dump_random_node] panmap failed (exit code {result}) or no candidates found for {ex['pan_stem']}")
                    print(f"[dump_random_node] candidates found: {candidates}")
                    raise RuntimeError(f"Failed to generate random node sequence for {ex['pan_stem']}")
                    
            finally:
                os.chdir(original_cwd)

rule index_experiment:
        input:
                bin="build/bin/panmap",
                pan=lambda wc: EXP_BY_ID[wc.eid]['panman_path']
        output:
                index=f"{OUTPUT_DIR}/experiments/{{eid}}/{{pan_stem}}/{{tag}}/indexes/index.pmi",
                mm=f"{OUTPUT_DIR}/experiments/{{eid}}/{{pan_stem}}/{{tag}}/indexes/{{pan_stem}}.panman.mm"
        params:
                pan=lambda wc: EXP_BY_ID[wc.eid]['panman_path'],
                k=lambda wc: EXP_BY_ID[wc.eid]['k'],
                s=lambda wc: EXP_BY_ID[wc.eid]['s']
        shell:
                r'''
                echo "[index_experiment] eid={wildcards.eid} pan_stem={wildcards.pan_stem} tag={wildcards.tag}" >&2
                mkdir -p $(dirname {output.index})
                {input.bin} -f -k {params.k} -s {params.s} {params.pan} --index-output {output.index}
                cp {params.pan}.mm {output.mm}
                echo "[index_experiment] Index built successfully: {output.index}" >&2
                '''

rule simulate_reads:
    threads: 4
    input:
        fasta=f"{OUTPUT_DIR}/experiments/{{eid}}/{{pan_stem}}/{{tag}}/genomes/reads/cov{{cov}}_{{n}}_rep{{rep}}/random_node.fasta",
        panman=lambda wc: EXP_BY_ID[wc.eid]['panman_path'],
        simulate_bin="build/bin/simulate",
        mm=f"{OUTPUT_DIR}/experiments/{{eid}}/{{pan_stem}}/{{tag}}/indexes/{{pan_stem}}.panman.mm"
    output:
        r1=f"{OUTPUT_DIR}/experiments/{{eid}}/{{pan_stem}}/{{tag}}/reads/cov{{cov}}_{{n}}_rep{{rep}}_R1.fastq",
        r2=f"{OUTPUT_DIR}/experiments/{{eid}}/{{pan_stem}}/{{tag}}/reads/cov{{cov}}_{{n}}_rep{{rep}}_R2.fastq",
        meta=f"{OUTPUT_DIR}/experiments/{{eid}}/{{pan_stem}}/{{tag}}/results/reads/cov{{cov}}_{{n}}_rep{{rep}}.txt",
        mut_fa=f"{OUTPUT_DIR}/experiments/{{eid}}/{{pan_stem}}/{{tag}}/mutgenomes/reads/cov{{cov}}_{{n}}_rep{{rep}}/mutated.fasta"
    params:
        model='NovaSeq',
        cpus=1,
        mm_type='both'  # Use both SNP and indel matrices
    run:
        ex = EXP_BY_ID[wildcards.eid]
        print(f"[simulate_reads] eid={wildcards.eid} pan_stem={wildcards.pan_stem} tag={wildcards.tag} cov={wildcards.cov} n={wildcards.n} rep={wildcards.rep}")
        out_dir = pathlib.Path(output.mut_fa).parent
        out_dir.mkdir(parents=True, exist_ok=True)
        reads_dir = pathlib.Path(output.r1).parent
        reads_dir.mkdir(parents=True, exist_ok=True)
        
        if pathlib.Path(input.simulate_bin).exists():
            print(f"[simulate_reads] running simulate for {ex['pangenome_name']} cov={wildcards.cov} n={wildcards.n} rep={wildcards.rep}")
            include_internal_flag = "--include_internal" if ex['include_internal'] else ""
            
            # Create deterministic seed based on experiment parameters 
            import hashlib
            seed_string = f"{ex['id']}_{wildcards.cov}_{wildcards.n}_{wildcards.rep}"
            seed_hash = hashlib.md5(seed_string.encode()).hexdigest()[:8]
            seed_int = int(seed_hash, 16) % (2**31 - 1)
            
            # Calculate mutations based on mutation rate and genome size
            mutation_rate = ex.get('mutation_rate', 0.0)
            genome_size = ex.get('genome_size', 15000)  # Default to RSV size
            
            # Calculate total number of mutations across the genome
            total_snps = int(mutation_rate * genome_size) if mutation_rate > 0 else 0
            total_indels = int(mutation_rate * genome_size / 4) if mutation_rate > 0 else 0  # Indels = mutation_rate * genome_size / 4
            
            print(f"[simulate_reads] Mutation rate: {mutation_rate}, Genome size: {genome_size}")
            print(f"[simulate_reads] Generating {total_snps} SNPs and {total_indels} indels")
            
            # Use mutation matrix if available
            mm_flag = f"--in_mm {input.mm}" if pathlib.Path(input.mm).exists() else ""
            
            # Capture simulate output to parse mutation counts
            import subprocess
            simulate_cmd = (
                f"{input.simulate_bin} --in_panman {input.panman} --in_fasta {input.fasta} --out_dir {out_dir} "
                f"--n_reads {wildcards.n} --n_replicates 1 --model {ex['model']} --cpus 4 "
                f"--seed {seed_int} --mutnum {total_snps} {total_indels} 0 --mm_type {params.mm_type} "
                f"{mm_flag} {include_internal_flag}"
            )
            
            print(f"[simulate_reads] Running: {simulate_cmd}")
            result = subprocess.run(simulate_cmd, shell=True, capture_output=True, text=True)
            simulate_output = result.stdout + result.stderr
            print(f"[simulate_reads] Simulate output:\n{simulate_output}")
            
            # Parse actual mutation counts from simulate output
            applied_snps = applied_insertions = applied_deletions = 0
            insertion_lengths = deletion_lengths = ""
            selected_node = "unknown"  # prefer original (unsanitized) if available
            selected_node_sanitized = "unknown"
            
            for line in simulate_output.split('\n'):
                if "Applied" in line and "SNPs" in line:
                    # Parse lines like: "Applied 15 SNPs to node KJ643480.1"
                    # or "Applied 2 SNPs, 1 insertions (lengths: 3), 1 deletions (lengths: 2) to node KJ643480.1"
                    import re
                    
                    # Extract SNP count
                    snp_match = re.search(r'(\d+) SNPs', line)
                    if snp_match:
                        applied_snps = int(snp_match.group(1))
                    
                    # Extract insertion count and lengths
                    ins_match = re.search(r'(\d+) insertions', line)
                    if ins_match:
                        applied_insertions = int(ins_match.group(1))
                        # Extract insertion lengths if present
                        ins_len_match = re.search(r'insertions \(lengths: ([0-9,]+)\)', line)
                        if ins_len_match:
                            insertion_lengths = ins_len_match.group(1)
                    
                    # Extract deletion count and lengths
                    del_match = re.search(r'(\d+) deletions', line)
                    if del_match:
                        applied_deletions = int(del_match.group(1))
                        # Extract deletion lengths if present
                        del_len_match = re.search(r'deletions \(lengths: ([0-9,]+)\)', line)
                        if del_len_match:
                            deletion_lengths = del_len_match.group(1)
                    
                    # Extract node name
                    if "to node" in line:
                        selected_node = line.split("to node")[-1].strip()
                elif line.strip().startswith("curNodeID:"):
                    # Prefer the original (unsanitized) node identifier when available
                    selected_node = line.split(":", 1)[-1].strip()
                elif "curNodeID (from FASTA):" in line:
                    selected_node = line.split(":", 1)[-1].strip()
            
            print(f"[simulate_reads] Parsed mutations: SNPs={applied_snps}, Insertions={applied_insertions} (lengths: {insertion_lengths}), Deletions={applied_deletions} (lengths: {deletion_lengths}), Node={selected_node}")
            
            import glob, shutil
            produced_fa = glob.glob(str(out_dir / '*varFasta*/*.fa')) or glob.glob(str(out_dir / '*.fa'))
            shutil.copyfile(produced_fa[0], output.mut_fa) if produced_fa else pathlib.Path(output.mut_fa).write_text(f">mut_{ex['pangenome_name']}\nN\n")
            produced_r1 = glob.glob(str(out_dir / '*reads*/*_R1.fastq'))
            produced_r2 = glob.glob(str(out_dir / '*reads*/*_R2.fastq'))
            shutil.copyfile(produced_r1[0], output.r1) if produced_r1 else pathlib.Path(output.r1).write_text("@r/1\nA\n+\nI\n")
            shutil.copyfile(produced_r2[0], output.r2) if produced_r2 else pathlib.Path(output.r2).write_text("@r/2\nT\n+\nI\n")
            
            # Extract true_node from the generated FASTQ read header
            true_node = selected_node  # default to original from simulate output
            try:
                with open(output.r1, 'r') as f:
                    first_line = f.readline().strip()
                    if first_line.startswith('@'):
                        # Extract node ID from FASTQ header (format: @NODE_ID_... or @NODE_ID.var_...)
                        header_node = first_line[1:]  # Remove @
                        if '.var_' in header_node:
                            selected_node_sanitized = header_node.split('.var_')[0]
                        elif '_' in header_node:
                            # ISS format: @NODE_ID_read_number/pair -> extract NODE_ID
                            selected_node_sanitized = '_'.join(header_node.split('_')[:-2]) if header_node.count('_') >= 2 else header_node.split('_')[0]
                        else:
                            selected_node_sanitized = header_node.split('/')[0]  # Remove /1 or /2 if present
                        
                        # Keep true_node as the original if we found it earlier; otherwise fall back to sanitized
                        if true_node == "unknown" or not true_node:
                            true_node = selected_node_sanitized
            except Exception:
                pass
        else:
            print(f"[simulate_reads] simulate_bin missing, writing synthetic outputs for {ex['pangenome_name']} n={wildcards.n}")
            pathlib.Path(output.mut_fa).write_text(f">mut_{ex['pangenome_name']}\n{'N'*100}\n")
            pathlib.Path(output.r1).write_text(f"@r/1\n{'A'*50}\n+\n{'I'*50}\n")
            pathlib.Path(output.r2).write_text(f"@r/2\n{'T'*50}\n+\n{'I'*50}\n")
            true_node = "synthetic"
            applied_snps = applied_insertions = applied_deletions = 0
            insertion_lengths = deletion_lengths = ""
        
        # Write comprehensive metadata including mutation counts
        pathlib.Path(output.meta).parent.mkdir(parents=True, exist_ok=True)
        with open(output.meta,'w') as meta_file:
            meta_file.write(f"id={ex['id']}\t")
            meta_file.write(f"coverage={wildcards.cov}\t")
            meta_file.write(f"reads={wildcards.n}\t")
            meta_file.write(f"{wildcards.tag}\t")
            meta_file.write(f"mutation_rate={ex['mutation_rate']}\t")
            meta_file.write(f"genome_size={ex['genome_size']}\t")
            meta_file.write(f"requested_snps={total_snps}\t")
            meta_file.write(f"requested_indels={total_indels}\t")
            meta_file.write(f"applied_snps={applied_snps}\t")
            meta_file.write(f"applied_insertions={applied_insertions}\t")
            meta_file.write(f"applied_deletions={applied_deletions}\t")
            meta_file.write(f"insertion_lengths={insertion_lengths}\t")
            meta_file.write(f"deletion_lengths={deletion_lengths}\t")
            # Record both forms for downstream consumers
            meta_file.write(f"true_node={true_node}\t")
            meta_file.write(f"true_node_sanitized={selected_node_sanitized}\n")

rule place_reads:
    input:
        panmap_bin="build/bin/panmap",
        panman=lambda wc: EXP_BY_ID[wc.eid]['panman_path'],
        index=f"{OUTPUT_DIR}/experiments/{{eid}}/{{pan_stem}}/{{tag}}/indexes/index.pmi",
        r1=f"{OUTPUT_DIR}/experiments/{{eid}}/{{pan_stem}}/{{tag}}/reads/cov{{cov}}_{{n}}_rep{{rep}}_R1.fastq",
        r2=f"{OUTPUT_DIR}/experiments/{{eid}}/{{pan_stem}}/{{tag}}/reads/cov{{cov}}_{{n}}_rep{{rep}}_R2.fastq"
    output:
        placement=f"{OUTPUT_DIR}/experiments/{{eid}}/{{pan_stem}}/{{tag}}/placements/reads/cov{{cov}}_{{n}}_rep{{rep}}/placements.tsv",
        detailed=f"{OUTPUT_DIR}/experiments/{{eid}}/{{pan_stem}}/{{tag}}/placements/reads/cov{{cov}}_{{n}}_rep{{rep}}/detailed.tsv",
        raw_log=f"{OUTPUT_DIR}/experiments/{{eid}}/{{pan_stem}}/{{tag}}/placements/reads/cov{{cov}}_{{n}}_rep{{rep}}/panmap.log"
    params:
        prefix=lambda wc: f"{OUTPUT_DIR}/experiments/{wc.eid}/{wc.pan_stem}/{wc.tag}/placements/reads/cov{wc.cov}_{wc.n}_rep{wc.rep}/result"
    run:
        ex = EXP_BY_ID[wildcards.eid]
        print(f"[place_reads] eid={wildcards.eid} pan_stem={wildcards.pan_stem} tag={wildcards.tag} n={wildcards.n} rep={wildcards.rep}")
        out_dir = pathlib.Path(output.placement).parent
        out_dir.mkdir(parents=True, exist_ok=True)
        if not pathlib.Path(input.panmap_bin).exists():
            print(f"[place_reads] panmap_bin missing, writing placeholder for {ex['pan_stem']} n={wildcards.n}")
            pathlib.Path(output.placement).write_text(f"PLACEMENTS placeholder\nindex={input.index}\nr1={input.r1}\nr2={input.r2}\n")
        else:
            print(f"[place_reads] running panmap placement for {ex['pan_stem']} n={wildcards.n}")
            import subprocess, re
            # Capture panmap output directly to parse placement results
            try:
                result = subprocess.run([
                    input.panmap_bin, input.panman, 
                    "-i", input.index,
                    "-k", str(ex['k']), "-s", str(ex['s']),
                    input.r1, input.r2
                ], capture_output=True, text=True, timeout=300)
                log_output = result.stderr + result.stdout
                
                # Save raw panmap log for debugging
                with open(output.raw_log, 'w') as f:
                    f.write(log_output)
                
                # Parse placement results from log output - enhanced to capture more details
                raw_hits = jaccard_sim = cosine_sim = weighted_score = 0
                raw_node = jaccard_node = cosine_node = weighted_node = ""
                all_candidates = []
                placement_details = {}
                
                for line in log_output.split('\n'):
                    line = line.strip()
                    if 'Raw Hits:' in line:
                        match = re.search(r'Raw Hits:\s*(\d+)\s*\(best:\s*([^)]+)\)', line)
                        if match:
                            raw_hits = int(match.group(1))
                            raw_node = match.group(2)
                            placement_details["raw"] = {"score": raw_hits, "best_node": raw_node}
                    elif 'Jaccard Similarity:' in line:
                        match = re.search(r'Jaccard Similarity:\s*([\d.]+)\s*\(best:\s*([^)]+)\)', line)
                        if match:
                            jaccard_sim = float(match.group(1))
                            jaccard_node = match.group(2)
                            placement_details["jaccard"] = {"score": jaccard_sim, "best_node": jaccard_node}
                    elif 'Cosine Similarity:' in line:
                        match = re.search(r'Cosine Similarity:\s*([\d.]+)\s*\(best:\s*([^)]+)\)', line)
                        if match:
                            cosine_sim = float(match.group(1))
                            cosine_node = match.group(2)
                            placement_details["cosine"] = {"score": cosine_sim, "best_node": cosine_node}
                    elif 'Weighted Score:' in line:
                        match = re.search(r'Weighted Score:\s*([\d.]+)\s*\(best:\s*([^)]+)\)', line)
                        if match:
                            weighted_score = float(match.group(1))
                            weighted_node = match.group(2)
                            placement_details["weighted"] = {"score": weighted_score, "best_node": weighted_node}
                    # Look for potential candidate lines or additional placement details
                    elif re.search(r'node_\d+|[A-Z]{2}\d+\.\d+', line) and ('score' in line.lower() or 'hit' in line.lower()):
                        all_candidates.append(line)
                
                # Write structured TSV with parsed results
                with open(output.placement, 'w') as f:
                    f.write("metric\tscore\thits\tnodes\n")
                    f.write(f"raw\t{raw_hits}\t{raw_hits}\t{raw_node}\n")
                    f.write(f"jaccard\t{jaccard_sim}\t\t{jaccard_node}\n")
                    f.write(f"cosine\t{cosine_sim}\t\t{cosine_node}\n")
                    f.write(f"weighted\t{weighted_score}\t\t{weighted_node}\n")
                
                # Write detailed placement information
                with open(output.detailed, 'w') as f:
                    f.write("experiment\treads\treplicate\ttrue_node\tmetric\tscore\tbest_node\ttimestamp\n")
                    import time
                    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
                    # Get true_node from metadata
                    true_node = "unknown"
                    try:
                        meta_file = f"{OUTPUT_DIR}/experiments/{wildcards.eid}/{wildcards.pan_stem}/{wildcards.tag}/results/reads/cov{wildcards.cov}_{wildcards.n}_rep{wildcards.rep}.txt"
                        if pathlib.Path(meta_file).exists():
                            with open(meta_file) as mf:
                                meta_line = mf.read().strip()
                                fields = dict(part.split('=',1) for part in meta_line.split('\t') if '=' in part)
                                true_node = fields.get("true_node", "unknown")
                    except Exception:
                        pass
                    
                    for metric, details in placement_details.items():
                        score = details["score"]
                        best_node = details["best_node"]
                        f.write(f"{wildcards.eid}\t{wildcards.n}\t{wildcards.rep}\t{true_node}\t{metric}\t{score}\t{best_node}\t{timestamp}\n")
                
                print(f"[place_reads] parsed: raw={raw_hits}, jaccard={jaccard_sim}, cosine={cosine_sim}, weighted={weighted_score}")
                
            except Exception as e:
                print(f"[place_reads] panmap failed or parsing error: {e}")
                # Write empty results on failure
                with open(output.placement, 'w') as f:
                    f.write("metric\tscore\thits\tnodes\n")
                    f.write("raw\t0\t0\t\n")
                    f.write("jaccard\t0\t\t\n")
                    f.write("cosine\t0\t\t\n")
                    f.write("weighted\t0\t\t\n")
                with open(output.detailed, 'w') as f:
                    f.write("experiment\treads\treplicate\ttrue_node\tmetric\tscore\tbest_node\ttimestamp\n")
                with open(output.raw_log, 'w') as f:
                    f.write(f"ERROR: {e}\n")

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
        alignment_results=f"{OUTPUT_DIR}/experiments/{{eid}}/{{pan_stem}}/{{tag}}/alignments/reads/cov{{cov}}_{{n}}_rep{{rep}}/accuracy.tsv",  # Keep original for compatibility
        log_dir=directory(f"{OUTPUT_DIR}/experiments/{{eid}}/{{pan_stem}}/{{tag}}/alignments/reads/cov{{cov}}_{{n}}_rep{{rep}}/logs")
    run:
        import subprocess, pathlib, tempfile, os, re

        out_dir = pathlib.Path(output.alignment_results).parent
        out_dir.mkdir(parents=True, exist_ok=True)

        log_dir = pathlib.Path(output.log_dir)
        log_dir.mkdir(parents=True, exist_ok=True)

        # Read metadata to get true_node
        true_node = "unknown"
        try:
            with open(input.meta) as mf:
                meta_line = mf.read().strip()
                fields = dict(part.split('=',1) for part in meta_line.split('\t') if '=' in part)
                true_node = fields.get("true_node", "unknown")
        except Exception:
            pass

        if true_node == "unknown":
            # Write empty results for all versions
            for output_file in [output.accuracy_0bp, output.accuracy_50bp, output.accuracy_150bp, output.accuracy_500bp, output.alignment_results]:
                with open(output_file, 'w') as f:
                    f.write("metric\ttrue_node\tplacement_node\tgenome_length\tsnps\tindels\ttotal_variants\talignment_length\tidentity\texcluded_bp\n")
            return

        # Read placement results for each metric
        placement_nodes = {}
        try:
            with open(input.placement) as pf:
                lines = pf.readlines()
                for line in lines[1:]:  # Skip header
                    parts = line.strip().split('\t')
                    if len(parts) >= 4:
                        metric = parts[0]
                        node = parts[3]
                        if node and node != "":
                            placement_nodes[metric] = node
        except Exception:
            pass
        
        # Helper to mirror C++ sanitizeFilename for output files
        def sanitize_filename(name: str) -> str:
            # Replace characters: / \ : * ? " < > |
            return re.sub(r'[\/\\:\*\?"<>|]', '_', name)

        # Function to process variants with different exclusion levels
        def count_variants_with_exclusion(cs_field, target_start, target_end, exclude_bp):
            """Count variants excluding specified bp from both ends"""
            snps = indels = 0
            
            # Calculate exclusion boundaries
            exclude_start = target_start + exclude_bp
            exclude_end = target_end - exclude_bp
            
            # Skip if alignment is too short for filtering
            if target_end - target_start < 2 * exclude_bp:
                return snps, indels
            
            position = target_start  # Current position in target sequence
            
            # Parse CS string element by element
            i = 0
            while i < len(cs_field):
                if cs_field[i] == ':':
                    # Match segment (:12)
                    i += 1
                    num_str = ""
                    while i < len(cs_field) and cs_field[i].isdigit():
                        num_str += cs_field[i]
                        i += 1
                    if num_str:
                        position += int(num_str)
                elif cs_field[i] == '*':
                    # SNP (*ag)
                    if exclude_start <= position < exclude_end:
                        snps += 1
                    position += 1
                    i += 3 if i + 2 < len(cs_field) else len(cs_field)
                elif cs_field[i] == '+':
                    # Insertion (+acgt)
                    if exclude_start <= position < exclude_end:
                        indels += 1
                    i += 1
                    while i < len(cs_field) and cs_field[i] in 'acgtACGT':
                        i += 1
                elif cs_field[i] == '-':
                    # Deletion (-acgt)
                    if exclude_start <= position < exclude_end:
                        indels += 1
                    i += 1
                    del_len = 0
                    while i < len(cs_field) and cs_field[i] in 'acgtACGT':
                        del_len += 1
                        i += 1
                    position += del_len
                else:
                    i += 1
            
            return snps, indels

        # Define exclusion levels and corresponding output files
        exclusion_configs = [
            (0, output.accuracy_0bp),
            (50, output.accuracy_50bp), 
            (150, output.accuracy_150bp),
            (500, output.accuracy_500bp)
        ]

        # Initialize all output files with headers
        for exclude_bp, output_file in exclusion_configs:
            with open(output_file, 'w') as f:
                f.write("metric\ttrue_node\tplacement_node\tgenome_length\tsnps\tindels\ttotal_variants\talignment_length\tidentity\texcluded_bp\n")
        
        
        for metric in ["raw", "jaccard", "cosine", "weighted"]:
            placement_node = placement_nodes.get(metric, "")
            
            if not placement_node or placement_node == true_node:
                # Perfect match or no placement - get genome length from true node
                genome_length = 0
                try:
                    # Extract true node sequence (panmap writes to file: {panman}.{sanitized(node)}.fa)
                    true_fa_path = f"{input.panman}.{sanitize_filename(true_node)}.fa"
                    # Remove existing file if present
                    if os.path.exists(true_fa_path):
                        os.unlink(true_fa_path)
                    # Run extraction
                    result = subprocess.run([
                        input.panmap_bin, input.panman,
                        "--dump-sequence", true_node
                    ], capture_output=True, text=True, timeout=60)
                    if result.returncode != 0:
                        raise Exception(f"dump-sequence failed for {true_node}: {result.stderr}")
                    # Verify file exists and has content
                    if not os.path.exists(true_fa_path) or os.path.getsize(true_fa_path) == 0:
                        raise Exception(f"FASTA not created: {true_fa_path}")
                    # Get genome length (count only sequence lines)
                    with open(true_fa_path) as tf:
                        for line in tf:
                            if not line.startswith('>'):
                                genome_length += len(line.strip())
                except Exception:
                    pass
                
                snps = indels = total_variants = 0
                alignment_length = genome_length
                identity = 1.0 if placement_node == true_node else 0.0
                
                # Write to all output files with different exclusion levels
                for exclude_bp, output_file in exclusion_configs:
                    with open(output_file, 'a') as f:
                        f.write(f"{metric}\t{true_node}\t{placement_node}\t{genome_length}\t{snps}\t{indels}\t{total_variants}\t{alignment_length}\t{identity:.6f}\t{exclude_bp}\n")
                
                # Write to original file for compatibility
                with open(output.alignment_results, 'a') as f:
                    f.write(f"{metric}\t{true_node}\t{placement_node}\t{genome_length}\t{snps}\t{indels}\t{total_variants}\t{alignment_length}\t{identity:.6f}\n")
            
            else:
                # Different nodes - need alignment
                print(f"[align_placement_accuracy] {metric}: aligning {true_node} vs {placement_node}")
                try:
                    # Extract both node sequences - panmap creates files with sanitized names
                    true_fa_path = f"{input.panman}.{sanitize_filename(true_node)}.fa"
                    place_fa_path = f"{input.panman}.{sanitize_filename(placement_node)}.fa"

                    # Remove existing files if they exist
                    if os.path.exists(true_fa_path):
                        os.unlink(true_fa_path)
                    if os.path.exists(place_fa_path):
                        os.unlink(place_fa_path)

                    # Extract true node
                    result1 = subprocess.run([
                        input.panmap_bin, input.panman,
                        "--dump-sequence", true_node
                    ], capture_output=True, text=True, timeout=60)
                    if result1.returncode != 0:
                        print(f"[align_placement_accuracy] Failed to extract {true_node}: {result1.stderr}")
                        raise Exception(f"Failed to extract {true_node}")

                    # Extract placement node
                    result2 = subprocess.run([
                        input.panmap_bin, input.panman,
                        "--dump-sequence", placement_node
                    ], capture_output=True, text=True, timeout=60)
                    if result2.returncode != 0:
                        print(f"[align_placement_accuracy] Failed to extract {placement_node}: {result2.stderr}")
                        raise Exception(f"Failed to extract {placement_node}")

                    if not (os.path.exists(true_fa_path) and os.path.exists(place_fa_path)):
                        raise Exception("FASTA files not created")

                    genome_length = 0
                    with open(true_fa_path) as tf:
                        for line in tf:
                            if not line.startswith('>'):
                                genome_length += len(line.strip())

                    paf_path = log_dir / f"{metric}_{true_node}_vs_{placement_node}.paf"
                    minimap_log_path = log_dir / f"{metric}_{true_node}_vs_{placement_node}_minimap.log"
                    with open(paf_path, 'w') as paf_file, open(minimap_log_path, 'w') as log_file:
                        result3 = subprocess.run([
                            "minimap2", "-cx", "asm20", "--cs",
                            true_fa_path, place_fa_path
                        ], stdout=paf_file, stderr=log_file, timeout=60, text=True)
                    if result3.returncode != 0:
                        with open(minimap_log_path) as log_f:
                            stderr_content = log_f.read()
                        raise Exception(f"minimap2 failed: {stderr_content}")

                    # Parse PAF
                    snps = indels = 0
                    alignment_length = 0
                    identity = 0.0
                    results_by_exclusion = {}
                    with open(paf_path) as paf:
                        for line in paf:
                            if line.startswith('#'):
                                continue
                            parts = line.strip().split('\t')
                            if len(parts) < 12:
                                continue
                            alignment_length = int(parts[10])
                            matches = int(parts[9])
                            cs_field = None
                            for field in parts[12:]:
                                if field.startswith('cs:Z:'):
                                    cs_field = field[5:]
                                    break
                            if cs_field:
                                target_start = int(parts[7])
                                target_end = int(parts[8])
                                for exclude_bp, _out in exclusion_configs:
                                    snps_excl, indels_excl = count_variants_with_exclusion(cs_field, target_start, target_end, exclude_bp)
                                    results_by_exclusion[exclude_bp] = (snps_excl, indels_excl)
                                if alignment_length > 0:
                                    identity = matches / alignment_length
                            break  # first alignment only

                    # Write per-exclusion outputs
                    for exclude_bp, out_file in exclusion_configs:
                        s_excl, i_excl = results_by_exclusion.get(exclude_bp, (0,0))
                        total_excl = s_excl + i_excl
                        with open(out_file, 'a') as f_out:
                            f_out.write(f"{metric}\t{true_node}\t{placement_node}\t{genome_length}\t{s_excl}\t{i_excl}\t{total_excl}\t{alignment_length}\t{identity:.6f}\t{exclude_bp}\n")

                    # Append legacy summary
                    total_variants = snps + indels
                    with open(output.alignment_results, 'a') as legacy:
                        legacy.write(f"{metric}\t{true_node}\t{placement_node}\t{genome_length}\t{snps}\t{indels}\t{total_variants}\t{alignment_length}\t{identity:.6f}\n")

                    # Save copies of fasta (debug)
                    import shutil
                    saved_true_fa = log_dir / f"{metric}_{true_node}.fa"
                    saved_place_fa = log_dir / f"{metric}_{placement_node}.fa"
                    shutil.copy2(true_fa_path, saved_true_fa)
                    shutil.copy2(place_fa_path, saved_place_fa)
                    os.unlink(true_fa_path)
                    os.unlink(place_fa_path)
                except Exception as e:
                    print(f"[align_placement_accuracy] Error for {metric}: {e}")
                    for exclude_bp, out_file in exclusion_configs:
                        with open(out_file, 'a') as f_out:
                            f_out.write(f"{metric}\t{true_node}\t{placement_node}\t0\t0\t0\t0\t0\t0.0\t{exclude_bp}\n")
                    with open(output.alignment_results, 'a') as legacy:
                        legacy.write(f"{metric}\t{true_node}\t{placement_node}\t0\t0\t0\t0\t0\t0.0\n")

rule aggregate_results:
    input:
        detailed=f"{OUTPUT_DIR}/experiments/{{eid}}/{{pan_stem}}/{{tag}}/placements/reads/cov{{cov}}_{{n}}_rep{{rep}}/detailed.tsv"
    output:
        touch(f"{OUTPUT_DIR}/results/aggregated/{{eid}}/{{pan_stem}}/{{tag}}/cov{{cov}}_{{n}}_rep{{rep}}.done")
    run:
        import time
        # Aggregate results incrementally
        results_dir = pathlib.Path(f"{OUTPUT_DIR}/results")
        results_dir.mkdir(parents=True, exist_ok=True)
        
        # Append to master results file
        master_file = results_dir / "all_placements.tsv"
        header_written = master_file.exists()
        
        with open(master_file, 'a') as master:
            if not header_written:
                master.write("experiment\treads\treplicate\ttrue_node\tmetric\tscore\tbest_node\ttimestamp\n")
            
            # Copy detailed results to master file
            try:
                with open(input.detailed, 'r') as detailed:
                    lines = detailed.readlines()
                    if len(lines) > 1:  # Skip header
                        for line in lines[1:]:
                            master.write(line)
                print(f"[aggregate_results] Added results for {wildcards.eid} cov={wildcards.cov} n={wildcards.n} rep={wildcards.rep}")
            except Exception as e:
                print(f"[aggregate_results] Error aggregating {input.detailed}: {e}")
        
        # Also create per-experiment summary
        exp_summary_file = results_dir / f"{wildcards.eid}_summary.tsv"
        exp_header_written = exp_summary_file.exists()
        
        with open(exp_summary_file, 'a') as exp_file:
            if not exp_header_written:
                exp_file.write("reads\treplicate\ttrue_node\tmetric\tscore\tbest_node\ttimestamp\n")
            
            try:
                with open(input.detailed, 'r') as detailed:
                    lines = detailed.readlines()
                    if len(lines) > 1:
                        for line in lines[1:]:
                            # Remove experiment column since it's in filename
                            parts = line.strip().split('\t')
                            if len(parts) >= 7:
                                exp_file.write('\t'.join(parts[1:]) + '\n')
            except Exception as e:
                print(f"[aggregate_results] Error writing experiment summary: {e}")

rule placements_summary:
    input:
        placements=PLACEMENTS,
        metas=RESULT_FILES
    output:
        tsv=f"{OUTPUT_DIR}/reports/placements_summary.tsv"
    run:
        import re
        pathlib.Path(output.tsv).parent.mkdir(parents=True, exist_ok=True)
        # Map meta path to metadata fields including mutation counts
        meta_map = {}
        for m in input.metas:
            try:
                with open(m) as fh:
                    line = fh.read().strip()
                # expecting key=val pairs separated by tabs
                fields = dict(part.split('=',1) for part in line.split('\t') if '=' in part)
                meta_map[m] = fields
            except Exception:
                meta_map[m] = {}
        
        header = ["id","pan_stem","tag","reads","replicate","true_node","top_node","top_score","rows",
                 "mutation_rate","genome_size","requested_snps","requested_indels",
                 "applied_snps","applied_insertions","applied_deletions","insertion_lengths","deletion_lengths"]
        with open(output.tsv, 'w') as tf:
            tf.write('\t'.join(header)+"\n")
            for (eid, pan_stem, tag, cov, n, rep) in READS:
                p = f"{OUTPUT_DIR}/experiments/{eid}/{pan_stem}/{tag}/placements/reads/cov{cov}_{n}_rep{rep}/placements.tsv"
                m = f"{OUTPUT_DIR}/experiments/{eid}/{pan_stem}/{tag}/results/reads/cov{cov}_{n}_rep{rep}.txt"
                genome_fa = f"{OUTPUT_DIR}/experiments/{eid}/{pan_stem}/{tag}/genomes/reads/cov{cov}_{n}_rep{rep}/random_node.fasta"
                
                # Determine true_node with fallback to random_node.fasta header
                meta_data = meta_map.get(m, {})
                true_node = meta_data.get('true_node', '')
                if not true_node or true_node in ("random_node", "unknown"):
                    try:
                        with open(genome_fa, 'r') as gf:
                            for ln in gf:
                                if ln.startswith('>'):
                                    hdr = ln[1:].strip()
                                    true_node = hdr.split()[0]
                                    break
                    except Exception:
                        true_node = "random_node"
                
                # Extract mutation information from metadata
                mutation_rate = meta_data.get('mutation_rate', '0')
                genome_size = meta_data.get('genome_size', '0')
                requested_snps = meta_data.get('requested_snps', '0')
                requested_indels = meta_data.get('requested_indels', '0')
                applied_snps = meta_data.get('applied_snps', '0')
                applied_insertions = meta_data.get('applied_insertions', '0')
                applied_deletions = meta_data.get('applied_deletions', '0')
                insertion_lengths = meta_data.get('insertion_lengths', '')
                deletion_lengths = meta_data.get('deletion_lengths', '')
                
                top_node = "n/a"
                top_score = "0"
                rows = n  # Use the read count as rows since we know it
                
                try:
                    with open(p, 'r') as pf:
                        # Skip header line
                        next(pf)
                        weighted_score = "0"
                        for line in pf:
                            parts = line.strip().split('\t')
                            if len(parts) >= 2:
                                metric = parts[0]
                                score = parts[1]
                                if metric == "weighted" and score != "0":
                                    weighted_score = score
                                elif metric == "cosine" and score != "0" and weighted_score == "0":
                                    weighted_score = score
                                elif metric == "jaccard" and score != "0" and weighted_score == "0":
                                    weighted_score = score
                        top_score = weighted_score
                except FileNotFoundError:
                    pass
                
                tf.write('\t'.join(map(str,[eid, pan_stem, tag, n, rep, true_node, top_node, top_score, rows,
                                          mutation_rate, genome_size, requested_snps, requested_indels,
                                          applied_snps, applied_insertions, applied_deletions, 
                                          insertion_lengths, deletion_lengths]))+"\n")

rule placements_by_metric:
    input:
        placements=PLACEMENTS,
        metas=RESULT_FILES
    output:
        tsv=f"{OUTPUT_DIR}/reports/placements_by_metric.tsv"
    run:
        import re
        pathlib.Path(output.tsv).parent.mkdir(parents=True, exist_ok=True)
        # Map meta path to true_node
        true_map = {}
        for m in input.metas:
            try:
                with open(m) as fh:
                    line = fh.read().strip()
                fields = dict(part.split('=',1) for part in line.split('\t') if '=' in part)
                true_map[m] = fields.get('true_node') or ''
            except Exception:
                true_map[m] = ''
        header = ["id","pan_stem","tag","reads","replicate","true_node","metric","score","nodes","hits","rows"]
        with open(output.tsv, 'w') as tf:
            tf.write('\t'.join(header)+"\n")
            for (eid, pan_stem, tag, cov, n, rep) in READS:
                p = f"{OUTPUT_DIR}/experiments/{eid}/{pan_stem}/{tag}/placements/reads/cov{cov}_{n}_rep{rep}/placements.tsv"
                m = f"{OUTPUT_DIR}/experiments/{eid}/{pan_stem}/{tag}/results/reads/cov{cov}_{n}_rep{rep}.txt"
                genome_fa = f"{OUTPUT_DIR}/experiments/{eid}/{pan_stem}/{tag}/genomes/reads/cov{cov}_{n}_rep{rep}/random_node.fasta"
                
                # Determine true_node with fallback to random_node.fasta header
                true_node = true_map.get(m, '')
                if not true_node or true_node in ("random_node", "unknown"):
                    try:
                        with open(genome_fa, 'r') as gf:
                            for ln in gf:
                                if ln.startswith('>'):
                                    hdr = ln[1:].strip()
                                    true_node = hdr.split()[0]
                                    break
                    except Exception:
                        true_node = "random_node"
                
                rows = 0
                raw_hits = None
                cosine_score = None
                jaccard_score = None
                raw_nodes = []
                cosine_nodes = []
                jaccard_nodes = []
                
                try:
                    with open(p, 'r') as pf:
                        # Skip header line
                        next(pf)
                        for line in pf:
                            parts = line.strip().split('\t')
                            if len(parts) >= 4:
                                metric = parts[0]
                                score = parts[1] if parts[1] else "0"
                                hits = parts[2] if parts[2] else ""
                                nodes = parts[3] if parts[3] else ""
                                
                                tf.write('\t'.join(map(str,[eid, pan_stem, tag, n, rep, true_node, metric, score, nodes, hits, n]))+"\n")
                except FileNotFoundError:
                    # Write empty rows for missing files
                    for metric in ["raw", "jaccard", "cosine", "weighted"]:
                        tf.write('\t'.join(map(str,[eid, pan_stem, tag, n, rep, true_node, metric, "0", "", "", n]))+"\n")

rule alignment_accuracy_summary:
    # Relaxed inputs: we intentionally do NOT require every accuracy.tsv file as an input
    # so that we can summarize partially generated accuracy variant files without
    # forcing creation of placeholder legacy files for missing replicates.
    # The rule will iterate over READS and gracefully handle missing variant files.
    input:
        metas=RESULT_FILES
    output:
        tsv=f"{OUTPUT_DIR}/reports/alignment_accuracy_summary.tsv"
    run:
        import pathlib
        pathlib.Path(output.tsv).parent.mkdir(parents=True, exist_ok=True)
        
        # Create comprehensive header with experiment info
        header = [
            "experiment_id", "pangenome", "tag", "coverage", "reads", "replicate",
            "metric", "true_node", "placement_node", "genome_length", 
            "snps", "indels", "total_variants", "alignment_length", "identity",
            "mutation_rate", "applied_snps", "applied_insertions", "applied_deletions",
            "excluded_bp"
        ]
        
        with open(output.tsv, 'w') as tf:
            tf.write('\t'.join(header) + "\n")
            
            for (eid, pan_stem, tag, cov, n, rep) in READS:
                # Get experiment info
                ex = EXP_BY_ID[eid]
                
                # Get metadata for mutation info
                meta_file = f"{OUTPUT_DIR}/experiments/{eid}/{pan_stem}/{tag}/results/reads/cov{cov}_{n}_rep{rep}.txt"
                mutation_rate = applied_snps = applied_insertions = applied_deletions = 0
                
                try:
                    with open(meta_file) as mf:
                        meta_line = mf.read().strip()
                        fields = dict(part.split('=',1) for part in meta_line.split('\t') if '=' in part)
                        mutation_rate = float(fields.get("mutation_rate", 0))
                        applied_snps = int(fields.get("applied_snps", 0))
                        applied_insertions = int(fields.get("applied_insertions", 0))
                        applied_deletions = int(fields.get("applied_deletions", 0))
                except Exception:
                    pass
                
                # Read all exclusion versions (0,50,150,500) plus legacy accuracy.tsv if present
                exclusion_levels = [0,50,150,500]
                base_path = f"{OUTPUT_DIR}/experiments/{eid}/{pan_stem}/{tag}/alignments/reads/cov{cov}_{n}_rep{rep}"
                for excl in exclusion_levels:
                    variant_file = f"{base_path}/accuracy_{excl}bp.tsv"
                    try:
                        with open(variant_file) as af:
                            lines = af.readlines()
                            for line in lines[1:]:  # skip header
                                parts = line.strip().split('\t')
                                # Expected columns: metric true_node placement_node genome_length snps indels total_variants alignment_length identity excluded_bp
                                if len(parts) >= 10:
                                    metric, true_node, placement_node, genome_length, snps, indels, total_variants, alignment_length, identity, excluded_bp = parts[:10]
                                    row = [
                                        eid, pan_stem, tag, cov, n, rep,
                                        metric, true_node, placement_node, genome_length,
                                        snps, indels, total_variants, alignment_length, identity,
                                        mutation_rate, applied_snps, applied_insertions, applied_deletions,
                                        excluded_bp
                                    ]
                                    tf.write('\t'.join(map(str,row)) + "\n")
                    except FileNotFoundError:
                        # Write placeholder rows if variant file missing
                        for metric in ["raw","jaccard","cosine","weighted"]:
                            row = [
                                eid, pan_stem, tag, cov, n, rep,
                                metric, "unknown", "", 0,
                                0, 0, 0, 0, 0.0,
                                mutation_rate, applied_snps, applied_insertions, applied_deletions,
                                excl
                            ]
                            tf.write('\t'.join(map(str,row)) + "\n")



rule compare_variants:
    input:
        panmap_bin="build/bin/panmap",
        panman=lambda wc: EXP_BY_ID[wc.eid]['panman_path'],
        true_node_fasta=f"{OUTPUT_DIR}/experiments/{{eid}}/{{pan_stem}}/{{tag}}/mutgenomes/reads/cov{{cov}}_{{n}}_rep{{rep}}/mutated.fasta",
        placement_results=f"{OUTPUT_DIR}/experiments/{{eid}}/{{pan_stem}}/{{tag}}/placements/reads/cov{{cov}}_{{n}}_rep{{rep}}/placements.tsv",
        meta=f"{OUTPUT_DIR}/experiments/{{eid}}/{{pan_stem}}/{{tag}}/results/reads/cov{{cov}}_{{n}}_rep{{rep}}.txt"
    output:
        variant_analysis=f"{OUTPUT_DIR}/experiments/{{eid}}/{{pan_stem}}/{{tag}}/variants/reads/cov{{cov}}_{{n}}_rep{{rep}}/variants.tsv"
    run:
        import subprocess
        ex = EXP_BY_ID[wildcards.eid]
        print(f"[compare_variants] eid={wildcards.eid} pan_stem={wildcards.pan_stem} tag={wildcards.tag} cov={wildcards.cov} n={wildcards.n} rep={wildcards.rep}")
        
        out_dir = pathlib.Path(output.variant_analysis).parent
        out_dir.mkdir(parents=True, exist_ok=True)
        
        # Get the true node from metadata
        true_node = "unknown"
        try:
            with open(input.meta) as mf:
                meta_line = mf.read().strip()
                fields = dict(part.split('=',1) for part in meta_line.split('\t') if '=' in part)
                true_node = fields.get("true_node", "unknown")
        except Exception:
            pass
        
        # Get the best predicted node from placement results
        predicted_node = "unknown"
        best_score = 0
        try:
            with open(input.placement_results) as pf:
                for line in pf:
                    if line.startswith('metric'):
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) >= 4 and parts[0] == 'weighted':
                        score = float(parts[1])
                        if score > best_score:
                            best_score = score
                            predicted_node = parts[3] if len(parts) > 3 else "unknown"
        except Exception as e:
            print(f"[compare_variants] Error reading placement results: {e}")
        
        # Create predicted node FASTA by extracting from pangenome
        predicted_fasta = out_dir / "predicted_node.fasta"
        if predicted_node != "unknown" and pathlib.Path(input.panmap_bin).exists():
            try:
                # Use panmap to dump the specific node sequence
                result = subprocess.run([
                    input.panmap_bin, input.panman, 
                    "--dump-node", predicted_node
                ], capture_output=True, text=True, cwd=out_dir)
                
                # Find the generated FASTA file
                node_files = list(out_dir.glob(f"{predicted_node}*.fa")) + list(out_dir.glob(f"{predicted_node}*.fasta"))
                if node_files:
                    shutil.move(str(node_files[0]), str(predicted_fasta))
                else:
                    # Fallback: create empty file
                    predicted_fasta.write_text(f">{predicted_node}\\nN\\n")
            except Exception as e:
                print(f"[compare_variants] Error extracting predicted node: {e}")
                predicted_fasta.write_text(f">{predicted_node}\\nN\\n")
        else:
            predicted_fasta.write_text(f">{predicted_node}\\nN\\n")
        
        # Align true and predicted sequences with minimap2
        paf_file = out_dir / "alignment.paf"
        try:
            subprocess.run([
                "minimap2", "-c", "-x", "asm5",
                str(predicted_fasta), str(input.true_node_fasta)
            ], stdout=open(paf_file, 'w'), stderr=subprocess.DEVNULL, check=True)
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            print(f"[compare_variants] minimap2 failed: {e}")
            # Create empty PAF file
            paf_file.write_text("")
        
        # Parse PAF and count variants
        total_snps = 0
        total_indels = 0
        true_node_length = 0
        predicted_node_length = 0
        alignment_length = 0
        
        # Get sequence lengths from FASTA files
        try:
            with open(input.true_node_fasta) as f:
                seq = ""
                for line in f:
                    if not line.startswith('>'):
                        seq += line.strip()
                true_node_length = len(seq)
        except Exception:
            pass
            
        try:
            with open(predicted_fasta) as f:
                seq = ""
                for line in f:
                    if not line.startswith('>'):
                        seq += line.strip()
                predicted_node_length = len(seq)
        except Exception:
            pass
        
        # Parse PAF for variants (simplified - counts mismatches and gaps)
        try:
            with open(paf_file) as paf:
                for line in paf:
                    fields = line.strip().split('\\t')
                    if len(fields) >= 12:
                        # PAF format: query_length, target_start, target_end, matches, alignment_length
                        matches = int(fields[9])
                        alignment_length = int(fields[10])
                        if alignment_length > 0:
                            mismatches = alignment_length - matches
                            total_snps += mismatches
                            
                        # Look for CIGAR string in optional fields for indels
                        for field in fields[12:]:
                            if field.startswith('cg:Z:'):
                                cigar = field[5:]
                                # Simple CIGAR parsing for indels
                                import re
                                insertions = sum(int(x) for x in re.findall(r'(\\d+)I', cigar))
                                deletions = sum(int(x) for x in re.findall(r'(\\d+)D', cigar))
                                total_indels += insertions + deletions
                                break
        except Exception as e:
            print(f"[compare_variants] Error parsing PAF: {e}")
        
        # Write results
        with open(output.variant_analysis, 'w') as f:
            f.write("experiment\\ttrue_node\\tpredicted_node\\ttrue_length\\tpredicted_length\\talignment_length\\tsnps\\tindels\\ttotal_variants\\tcoverage\\treads\\treplicate\\n")
            f.write(f"{wildcards.eid}\\t{true_node}\\t{predicted_node}\\t{true_node_length}\\t{predicted_node_length}\\t{alignment_length}\\t{total_snps}\\t{total_indels}\\t{total_snps + total_indels}\\t{wildcards.cov}\\t{wildcards.n}\\t{wildcards.rep}\\n")

rule plot_placement_accuracy:
    input:
        accuracy_files=ALIGNMENT_ACCURACY
    output:
        done=touch(f"{OUTPUT_DIR}/plots/placement_accuracy_plots.done"),
        plots_dir=directory(f"{OUTPUT_DIR}/plots")
    run:
        import pandas as pd
        import numpy as np
        from pathlib import Path
        import warnings
        warnings.filterwarnings('ignore')
        
        # Create plots directory
        plots_dir = Path(output.plots_dir)
        plots_dir.mkdir(parents=True, exist_ok=True)
        
        # Filter for files that actually exist
        existing_files = [f for f in input.accuracy_files if Path(f).exists()]
        
        print(f"[plot_placement_accuracy] Input specified {len(input.accuracy_files)} accuracy files")
        print(f"[plot_placement_accuracy] Found {len(existing_files)} existing accuracy files")
        
        if len(existing_files) == 0:
            print(f"[plot_placement_accuracy] No accuracy files found, creating empty plots directory")
            # Create a summary file indicating no data
            summary_file = plots_dir / "summary.txt"
            with open(summary_file, 'w') as f:
                f.write("No accuracy data available for plotting.\n")
            return
        
        # Load and combine all available accuracy data
        all_data = []
        for accuracy_file in existing_files:
            try:
                print(f"[plot_placement_accuracy] Loading {accuracy_file}")
                df = pd.read_csv(accuracy_file, sep='\t')
                
                # Extract experiment info from file path
                # Format: workflow_output/experiments/{eid}/{pan_stem}/{tag}/alignments/reads/cov{cov}_{n}_rep{rep}/accuracy.tsv
                path_parts = Path(accuracy_file).parts
                if len(path_parts) >= 8:
                    eid = path_parts[2]
                    pan_stem = path_parts[3] 
                    tag = path_parts[4]
                    reads_info = path_parts[7]  # cov{cov}_{n}_rep{rep}
                    
                    # Parse cov{cov}_{n}_rep{rep}
                    import re
                    match = re.match(r'cov([\d.]+)_(\d+)_rep(\d+)', reads_info)
                    if match:
                        coverage = float(match.group(1))
                        reads = int(match.group(2))
                        replicate = int(match.group(3))
                        
                        # Add metadata to dataframe
                        df['experiment'] = eid
                        df['pan_stem'] = pan_stem
                        df['tag'] = tag
                        df['coverage'] = coverage
                        df['reads'] = reads
                        df['replicate'] = replicate
                        
                        # Extract k, s, and mutation rate from tag
                        # Tag format: k{k}_s{s}_noint_mut{rate}
                        k_match = re.search(r'k(\d+)', tag)
                        s_match = re.search(r's(\d+)', tag)
                        mut_match = re.search(r'mut([\d.]+)', tag)
                        
                        df['k'] = int(k_match.group(1)) if k_match else 0
                        df['s'] = int(s_match.group(1)) if s_match else 0
                        df['k_s_params'] = f"k{df['k'].iloc[0]}_s{df['s'].iloc[0]}"
                        
                        if mut_match:
                            df['mutation_rate'] = float(mut_match.group(1))
                        else:
                            df['mutation_rate'] = 0.0
                        
                        # Try to read detailed mutation info from metadata file
                        # Format coverage to match the original file naming convention
                        coverage_str = str(int(coverage)) if coverage == int(coverage) else str(coverage)
                        meta_file = f"workflow_output/experiments/{eid}/{pan_stem}/{tag}/results/reads/cov{coverage_str}_{reads}_rep{replicate}.txt"
                        mutation_breakdown = "unknown"
                        calculated_mutation_rate = 0.0
                        calculated_genome_size = 0
                        try:
                            if pathlib.Path(meta_file).exists():
                                with open(meta_file, 'r') as mf:
                                    meta_line = mf.read().strip()
                                    # Split on both tabs and spaces, then filter for parts with '='
                                    all_parts = re.split(r'[\t\s]+', meta_line)
                                    fields = dict(part.split('=',1) for part in all_parts if '=' in part)
                                    
                                    # Extract mutation rate and genome size for consistent calculation
                                    meta_mutation_rate = float(fields.get('mutation_rate', 0.0))
                                    meta_genome_size = int(fields.get('genome_size', 0))
                                    
                                    if meta_mutation_rate > 0 and meta_genome_size > 0:
                                        # Use the same calculation as in simulation
                                        calculated_snps = int(meta_mutation_rate * meta_genome_size)
                                        calculated_indels = int(meta_mutation_rate * meta_genome_size / 4)
                                        calculated_mutation_rate = meta_mutation_rate
                                        calculated_genome_size = meta_genome_size
                                        
                                        # Create breakdown based on calculated values for consistency
                                        parts = []
                                        if calculated_snps > 0:
                                            parts.append(f"{calculated_snps}snp")
                                        if calculated_indels > 0:
                                            parts.append(f"{calculated_indels}indel")
                                        
                                        if parts:
                                            mutation_breakdown = '+'.join(parts)
                                        else:
                                            mutation_breakdown = "0mut"
                                    else:
                                        # Fallback to applied values if no mutation rate/genome size
                                        applied_snps = int(fields.get('applied_snps', 0))
                                        applied_insertions = int(fields.get('applied_insertions', 0))
                                        applied_deletions = int(fields.get('applied_deletions', 0))
                                        
                                        parts = []
                                        if applied_snps > 0:
                                            parts.append(f"{applied_snps}snp")
                                        if applied_insertions > 0:
                                            parts.append(f"{applied_insertions}ins")
                                        if applied_deletions > 0:
                                            parts.append(f"{applied_deletions}del")
                                        
                                        if parts:
                                            mutation_breakdown = '+'.join(parts)
                                        else:
                                            mutation_breakdown = "0mut"
                                        
                        except Exception as e:
                            print(f"[plot_placement_accuracy] Could not parse metadata {meta_file}: {e}")
                            print(f"[plot_placement_accuracy] Meta line: {meta_line if 'meta_line' in locals() else 'NO_LINE'}")
                            mutation_breakdown = "unknown"
                        
                        df['mutation_breakdown'] = mutation_breakdown
                        df['calculated_mutation_rate'] = calculated_mutation_rate
                        df['calculated_genome_size'] = calculated_genome_size
                        
                        all_data.append(df)
                        
            except Exception as e:
                print(f"[plot_placement_accuracy] Error loading {accuracy_file}: {e}")
                continue
        
        if len(all_data) == 0:
            print(f"[plot_placement_accuracy] No valid data loaded")
            summary_file = plots_dir / "summary.txt"
            with open(summary_file, 'w') as f:
                f.write("No valid accuracy data could be loaded for plotting.\n")
            return
        
        # Combine all data
        df = pd.concat(all_data, ignore_index=True)
        print(f"[plot_placement_accuracy] Combined {len(df)} records from {len(all_data)} files")
        print(f"[plot_placement_accuracy] Columns: {list(df.columns)}")
        
        # Calculate distance from expected genome (total variants)
        if 'total_variants' in df.columns:
            df['distance'] = df['total_variants']
        else:
            df['distance'] = df.get('snps', 0) + df.get('indels', 0)
        
        # Count replicates for each combination to add to legend
        # Group by reads and mutation rate to count replicates
        if 'calculated_mutation_rate' in df.columns and 'reads' in df.columns and 'replicate' in df.columns:
            replicate_counts = df.groupby(['reads', 'calculated_mutation_rate'])['replicate'].nunique().reset_index()
            replicate_counts.columns = ['reads', 'calculated_mutation_rate', 'n_replicates']
            df = df.merge(replicate_counts, on=['reads', 'calculated_mutation_rate'], how='left')
        elif 'mutation_breakdown' in df.columns and 'reads' in df.columns and 'replicate' in df.columns:
            replicate_counts = df.groupby(['reads', 'mutation_breakdown'])['replicate'].nunique().reset_index()
            replicate_counts.columns = ['reads', 'mutation_breakdown', 'n_replicates']
            df = df.merge(replicate_counts, on=['reads', 'mutation_breakdown'], how='left')
        else:
            # Default to counting unique replicates per reads
            if 'replicate' in df.columns and 'reads' in df.columns:
                replicate_counts = df.groupby(['reads'])['replicate'].nunique().reset_index()
                replicate_counts.columns = ['reads', 'n_replicates']
                df = df.merge(replicate_counts, on=['reads'], how='left')
            else:
                df['n_replicates'] = 1
        
        # Create read count + mutation rate categories for legend with replicate count
        if 'calculated_mutation_rate' in df.columns and 'reads' in df.columns:
            # Format mutation rate for legend (e.g., "1e-4" or "0.0001")
            df['mutation_rate_str'] = df['calculated_mutation_rate'].apply(
                lambda x: f"{x:.0e}" if x > 0 and x < 0.001 else f"{x:.4f}" if x > 0 else "0"
            )
            df['category'] = df['reads'].astype(str) + ' reads, μ=' + df['mutation_rate_str'].astype(str) + ' (n=' + df['n_replicates'].astype(str) + ')'
            # Add separate columns for proper sorting
            df['reads_sort'] = df['reads']
            df['mutations_sort'] = df['calculated_mutation_rate'] * 10000  # Use mutation rate for sorting
        elif 'mutation_breakdown' in df.columns and 'reads' in df.columns:
            # Fallback to mutation breakdown if mutation rate not available
            df['category'] = df['reads'].astype(str) + ' reads, ' + df['mutation_breakdown'].astype(str) + ' (n=' + df['n_replicates'].astype(str) + ')'
            # Add separate columns for proper sorting
            df['reads_sort'] = df['reads']
            df['mutations_sort'] = df.get('mutation_rate', 0) * 10000  # Use mutation rate for sorting
        else:
            df['category'] = df.get('reads', 'unknown').astype(str) + ' reads (n=' + df['n_replicates'].astype(str) + ')'
            df['reads_sort'] = df.get('reads', 0)
            df['mutations_sort'] = 0
        
        # Create a detailed summary TSV file instead of plots for now
        summary_file = plots_dir / "placement_accuracy_summary.tsv"
        
        print(f"[plot_placement_accuracy] Creating summary table")
        
        # Calculate summary statistics
        summary_data = []
        metrics = df['metric'].unique() if 'metric' in df.columns else ['all']
        
        for metric in metrics:
            if metric != 'all':
                metric_data = df[df['metric'] == metric].copy()
            else:
                metric_data = df.copy()
            
            if len(metric_data) == 0:
                continue
            
            # Sort categories by reads first, then by mutations
            if 'reads_sort' in metric_data.columns and 'mutations_sort' in metric_data.columns:
                category_order = metric_data[['category', 'reads_sort', 'mutations_sort']].drop_duplicates()
                category_order = category_order.sort_values(['reads_sort', 'mutations_sort'])
                categories = category_order['category'].tolist()
            else:
                categories = sorted(metric_data['category'].unique())
            
            for category in categories:
                cat_data = metric_data[metric_data['category'] == category]
                total_samples = len(cat_data)
                
                if total_samples == 0:
                    continue
                
                # Calculate distance distribution
                for distance in range(6):  # 0-5+ variants
                    if distance < 5:
                        count = len(cat_data[cat_data['distance'] == distance])
                        distance_label = str(distance)
                    else:
                        count = len(cat_data[cat_data['distance'] >= 5])
                        distance_label = "5+"
                    
                    proportion = count / total_samples if total_samples > 0 else 0
                    
                    summary_data.append({
                        'metric': metric,
                        'category': category,
                        'distance': distance_label,
                        'count': count,
                        'total_samples': total_samples,
                        'proportion': proportion
                    })
        
        # Save summary table
        summary_df = pd.DataFrame(summary_data)
        summary_df.to_csv(summary_file, sep='\t', index=False)
        print(f"[plot_placement_accuracy] Saved summary table: {summary_file}")
        
        # Also save the raw combined data for future plotting
        raw_data_file = plots_dir / "combined_accuracy_data.tsv"
        df.to_csv(raw_data_file, sep='\t', index=False)
        print(f"[plot_placement_accuracy] Saved raw data: {raw_data_file}")
        
        # Create a simple text report
        report_file = plots_dir / "placement_accuracy_report.txt"
        with open(report_file, 'w') as f:
            f.write("Placement Accuracy Analysis Report\n")
            f.write("==================================\n\n")
            f.write(f"Total accuracy files processed: {len(all_data)}\n")
            f.write(f"Total records: {len(df)}\n")
            f.write(f"Metrics analyzed: {', '.join(metrics)}\n")
            f.write(f"Categories found: {len(df['category'].unique())}\n\n")
            
            f.write("Summary by metric:\n")
            for metric in metrics:
                if metric != 'all':
                    metric_data = df[df['metric'] == metric]
                else:
                    metric_data = df
                    
                perfect_matches = len(metric_data[metric_data['distance'] == 0])
                total = len(metric_data)
                accuracy = perfect_matches / total * 100 if total > 0 else 0
                
                f.write(f"  {metric}: {perfect_matches}/{total} perfect matches ({accuracy:.1f}%)\n")
        
        print(f"[plot_placement_accuracy] Created report: {report_file}")
        
        # Create visual plots with matplotlib - separate by (k,s) parameters
        try:
            import matplotlib.pyplot as plt
            import seaborn as sns
            import numpy as np
            
            plt.style.use('default')
            sns.set_palette("viridis")
            
            # Calculate distance from expected genome (total variants)
            if 'total_variants' in df.columns:
                df['distance'] = df['total_variants']
            else:
                df['distance'] = df.get('snps', 0) + df.get('indels', 0)
            
            # Get unique (k,s) parameter combinations
            k_s_combinations = df['k_s_params'].unique() if 'k_s_params' in df.columns else ['all']
            
            for k_s_params in k_s_combinations:
                print(f"[plot_placement_accuracy] Creating plots for {k_s_params}")
                
                # Filter data for this (k,s) combination
                if k_s_params != 'all':
                    k_s_data = df[df['k_s_params'] == k_s_params].copy()
                else:
                    k_s_data = df.copy()
                
                if len(k_s_data) == 0:
                    print(f"[plot_placement_accuracy] No data for {k_s_params}")
                    continue
                
                # Create plots for each metric within this (k,s) combination
                metrics = k_s_data['metric'].unique() if 'metric' in k_s_data.columns else ['combined']
                
                for metric in metrics:
                    print(f"[plot_placement_accuracy] Creating plot for {k_s_params}, metric: {metric}")
                    
                    # Filter data for this metric
                    if metric != 'combined':
                        metric_data = k_s_data[k_s_data['metric'] == metric].copy()
                    else:
                        metric_data = k_s_data.copy()
                    
                    if len(metric_data) == 0:
                        print(f"[plot_placement_accuracy] No data for {k_s_params}, metric {metric}")
                        continue
                    
                    # Create the plot
                    fig, ax = plt.subplots(figsize=(12, 8))
                    
                    # Calculate proportions for each category and distance
                    # Sort categories by reads first, then by mutations
                    if 'reads_sort' in metric_data.columns and 'mutations_sort' in metric_data.columns:
                        category_order = metric_data[['category', 'reads_sort', 'mutations_sort']].drop_duplicates()
                        category_order = category_order.sort_values(['reads_sort', 'mutations_sort'])
                        categories = category_order['category'].tolist()
                    else:
                        categories = sorted(metric_data['category'].unique())
                        
                    max_distance = min(int(metric_data['distance'].max()), 10) if len(metric_data) > 0 else 5
                    
                    # Create distance bins (0, 1, 2, 3, 4, >=5)
                    distance_bins = list(range(min(max_distance + 1, 6))) + ['>5']
                    
                    # Prepare data for plotting
                    plot_data = []
                    for cat in categories:
                        cat_data = metric_data[metric_data['category'] == cat]
                        total_samples = len(cat_data)
                        
                        if total_samples == 0:
                            continue
                        
                        # Count samples in each distance bin
                        distance_counts = []
                        for d in range(min(max_distance + 1, 6)):
                            count = len(cat_data[cat_data['distance'] == d])
                            distance_counts.append(count / total_samples)
                        
                        # Count samples with distance > 5
                        count_high = len(cat_data[cat_data['distance'] > 5])
                        distance_counts.append(count_high / total_samples)
                        
                        plot_data.append(distance_counts)
                    
                    # Create the bar plot
                    x = np.arange(len(distance_bins))
                    width = 0.8 / len(categories) if len(categories) > 0 else 0.8
                    colors = plt.cm.viridis(np.linspace(1, 0, len(categories)))  # Inverted gradient
                    
                    for i, (cat, data) in enumerate(zip(categories, plot_data)):
                        offset = (i - len(categories)/2 + 0.5) * width
                        bars = ax.bar(x + offset, data, width, label=cat, alpha=0.8, color=colors[i])
                    
                    # Customize the plot
                    ax.set_xlabel('Dist. from expected genome (# SNPs + Indels)', fontsize=12)
                    ax.set_ylabel('Proportion of samples', fontsize=12)
                    ax.set_title(f'Placement accuracy, {k_s_params}, {metric} metric', fontsize=14, fontweight='bold')
                    ax.set_xticks(x)
                    ax.set_xticklabels([str(d) for d in distance_bins])
                    ax.set_ylim(0, 1.0)
                    
                    # Add reference line at y=1.0
                    ax.axhline(y=1.0, color='black', linestyle='--', alpha=0.3)
                    
                    # Add legend
                    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
                    
                    # Add grid
                    ax.grid(True, alpha=0.3)
                    
                    # Tight layout
                    plt.tight_layout()
                    
                    # Save the plot with (k,s) parameter in filename
                    plot_file = plots_dir / f"placement_accuracy_{k_s_params}_{metric}.png"
                    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
                    print(f"[plot_placement_accuracy] Saved plot: {plot_file}")
                    
                    # Also save as PDF
                    plot_file_pdf = plots_dir / f"placement_accuracy_{k_s_params}_{metric}.pdf"
                    plt.savefig(plot_file_pdf, bbox_inches='tight')
                    
                    plt.close()
            
            print(f"[plot_placement_accuracy] Created visual plots for {len(k_s_combinations)} (k,s) combinations")
            
            # Create comprehensive PDF with all plots in a grid
            print(f"[plot_placement_accuracy] Creating comprehensive PDF with all plots...")
            
            metrics = ['raw', 'jaccard', 'cosine', 'weighted']
            n_metrics = len(metrics)
            n_k_s_combinations = len(k_s_combinations)
            
            # Calculate grid dimensions (rows=k_s combinations, cols=metrics)
            n_rows = n_k_s_combinations
            n_cols = n_metrics
            
            # Create large figure for comprehensive view
            fig_width = n_cols * 5  # 5 inches per plot
            fig_height = n_rows * 4  # 4 inches per plot
            fig, axes = plt.subplots(n_rows, n_cols, figsize=(fig_width, fig_height))
            
            # Handle case where we have only one row or column
            if n_rows == 1 and n_cols == 1:
                axes = [[axes]]
            elif n_rows == 1:
                axes = [axes]
            elif n_cols == 1:
                axes = [[ax] for ax in axes]
            
            # Plot each combination
            for row_idx, k_s_params in enumerate(k_s_combinations):
                k_s_data = df[df['k_s_params'] == k_s_params]
                
                for col_idx, metric in enumerate(metrics):
                    ax = axes[row_idx][col_idx]
                    
                    if metric != 'combined':
                        metric_data = k_s_data[k_s_data['metric'] == metric].copy()
                    else:
                        metric_data = k_s_data.copy()
                    
                    if len(metric_data) == 0:
                        ax.text(0.5, 0.5, f'No data\n{k_s_params}\n{metric}', 
                               ha='center', va='center', transform=ax.transAxes)
                        ax.set_xticks([])
                        ax.set_yticks([])
                        continue
                    
                    # Sort categories by reads first, then by mutations
                    if 'reads_sort' in metric_data.columns and 'mutations_sort' in metric_data.columns:
                        category_order = metric_data[['category', 'reads_sort', 'mutations_sort']].drop_duplicates()
                        category_order = category_order.sort_values(['reads_sort', 'mutations_sort'])
                        categories = category_order['category'].tolist()
                    else:
                        categories = sorted(metric_data['category'].unique())
                        
                    max_distance = min(int(metric_data['distance'].max()), 10) if len(metric_data) > 0 else 5
                    distance_bins = list(range(min(max_distance + 1, 6))) + ['>5']
                    
                    # Prepare data for plotting
                    plot_data = []
                    for cat in categories:
                        cat_data = metric_data[metric_data['category'] == cat]
                        total_samples = len(cat_data)
                        
                        if total_samples == 0:
                            continue
                        
                        # Count samples in each distance bin
                        distance_counts = []
                        for d in range(min(max_distance + 1, 6)):
                            count = len(cat_data[cat_data['distance'] == d])
                            distance_counts.append(count / total_samples)
                        
                        # Count samples with distance > 5
                        count_high = len(cat_data[cat_data['distance'] > 5])
                        distance_counts.append(count_high / total_samples)
                        
                        plot_data.append(distance_counts)
                    
                    if len(plot_data) > 0:
                        # Create the bar plot
                        x = np.arange(len(distance_bins))
                        width = 0.8 / len(categories) if len(categories) > 0 else 0.8
                        colors = plt.cm.viridis(np.linspace(1, 0, len(categories)))  # Inverted gradient
                        
                        for i, (cat, data) in enumerate(zip(categories, plot_data)):
                            offset = (i - len(categories)/2 + 0.5) * width
                            bars = ax.bar(x + offset, data, width, alpha=0.8, color=colors[i])
                        
                        # Customize the subplot
                        ax.set_xticks(x)
                        ax.set_xticklabels([str(d) for d in distance_bins], fontsize=8)
                        ax.set_ylim(0, 1.0)
                        ax.grid(True, alpha=0.3)
                        
                        # Only add labels on edge plots to save space
                        if row_idx == n_rows - 1:  # Bottom row
                            ax.set_xlabel('Dist. from expected genome (# SNPs + Indels)', fontsize=10)
                        if col_idx == 0:  # Left column
                            ax.set_ylabel('Proportion of samples', fontsize=10)
                        
                        # Add title for each subplot
                        ax.set_title(f'{k_s_params}, {metric}', fontsize=11, fontweight='bold')
                        
                        # Add reference line at y=1.0
                        ax.axhline(y=1.0, color='black', linestyle='--', alpha=0.3)
                    else:
                        ax.text(0.5, 0.5, f'No data\n{k_s_params}\n{metric}', 
                               ha='center', va='center', transform=ax.transAxes)
                        ax.set_xticks([])
                        ax.set_yticks([])
            
            # Add overall title
            fig.suptitle('Placement Accuracy Analysis by (k,s) Parameters and Metrics', 
                        fontsize=16, fontweight='bold', y=0.98)
            
            # Create a legend for the entire figure (use data from first non-empty plot)
            legend_categories = None
            legend_colors = None
            for k_s_params in k_s_combinations:
                k_s_data = df[df['k_s_params'] == k_s_params]
                if len(k_s_data) > 0:
                    metric_data = k_s_data[k_s_data['metric'] == 'raw']  # Use raw metric for legend
                    if len(metric_data) > 0 and 'category' in metric_data.columns:
                        if 'reads_sort' in metric_data.columns and 'mutations_sort' in metric_data.columns:
                            category_order = metric_data[['category', 'reads_sort', 'mutations_sort']].drop_duplicates()
                            category_order = category_order.sort_values(['reads_sort', 'mutations_sort'])
                            legend_categories = category_order['category'].tolist()
                        else:
                            legend_categories = sorted(metric_data['category'].unique())
                        legend_colors = plt.cm.viridis(np.linspace(1, 0, len(legend_categories)))
                        break
            
            if legend_categories and legend_colors is not None:
                # Create legend handles
                legend_handles = [plt.Rectangle((0,0),1,1, color=color, alpha=0.8) 
                                for color in legend_colors]
                fig.legend(legend_handles, legend_categories, 
                          loc='center right', bbox_to_anchor=(0.98, 0.5), 
                          title='Read count & mutation rate (μ)', fontsize=10)
            
            # Adjust layout to make room for legend
            plt.tight_layout()
            plt.subplots_adjust(right=0.85, top=0.94)
            
            # Save comprehensive PDF
            comprehensive_pdf = plots_dir / "placement_accuracy_comprehensive.pdf"
            plt.savefig(comprehensive_pdf, bbox_inches='tight', dpi=300)
            print(f"[plot_placement_accuracy] Saved comprehensive PDF: {comprehensive_pdf}")
            
            # Also save as PNG for quick viewing
            comprehensive_png = plots_dir / "placement_accuracy_comprehensive.png"
            plt.savefig(comprehensive_png, bbox_inches='tight', dpi=300)
            print(f"[plot_placement_accuracy] Saved comprehensive PNG: {comprehensive_png}")
            
            plt.close()
            
        except ImportError:
            print(f"[plot_placement_accuracy] Note: Install matplotlib/seaborn for visual plots")
        
        print(f"[plot_placement_accuracy] Completed analysis of {len(df)} records")


rule clean:
    run:
        import shutil, pathlib, os, glob
        outdir = pathlib.Path(OUTPUT_DIR)
        print(f"[clean] Removing {outdir} if it exists...")
        if outdir.exists():
            shutil.rmtree(outdir, ignore_errors=True)
        # Remove Snakemake working directory
        sm = pathlib.Path('.snakemake')
        if sm.exists():
            print(f"[clean] Removing {sm}...")
            shutil.rmtree(sm, ignore_errors=True)
        # Remove per-panman byproducts in their directories
        for e in EXPERIMENTS:
            pan = pathlib.Path(e['panman_path']).resolve()
            pan_dir = pan.parent
            patterns = [
                str(pan_dir / '*.placement.tsv'),
                f"{str(pan)}.placement.tsv",
                str(pan_dir / '*.random.*.fa'),
                f"{str(pan)}.random.*.fa",
                f"{str(pan)}.pmi",
                str(pan_dir / 'read_seeds_debug.log'),
            ]
            for pat in patterns:
                for f in glob.glob(pat):
                    try:
                        print(f"[clean] Removing {f}")
                        os.remove(f)
                    except FileNotFoundError:
                        pass
