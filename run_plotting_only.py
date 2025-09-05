#!/usr/bin/env python3

import pandas as pd
import numpy as np
from pathlib import Path
import warnings
import re
warnings.filterwarnings('ignore')

def main():
    # Create plots directory
    plots_dir = Path("workflow_output/plots")
    plots_dir.mkdir(parents=True, exist_ok=True)
    
    # Find all existing accuracy files (legacy + exclusion variants)
    accuracy_files = list(Path("workflow_output").glob("**/accuracy.tsv"))
    # New exclusion-aware files
    for pattern in ["**/accuracy_0bp.tsv","**/accuracy_50bp.tsv","**/accuracy_150bp.tsv","**/accuracy_500bp.tsv"]:
        accuracy_files.extend(Path("workflow_output").glob(pattern))
    existing_files = [str(f) for f in accuracy_files if f.exists()]
    
    print(f"Found {len(existing_files)} existing accuracy files")
    
    if len(existing_files) == 0:
        print("No accuracy files found")
        return
    
    # Load and combine all available accuracy data
    all_data = []
    for accuracy_file in existing_files:
        try:
            print(f"Loading {accuracy_file}")
            df = pd.read_csv(accuracy_file, sep='\t')
            # Add excluded_bp column if present/absent
            if 'excluded_bp' not in df.columns:
                # Infer from filename if possible
                if accuracy_file.endswith('accuracy_0bp.tsv'):
                    df['excluded_bp'] = 0
                elif accuracy_file.endswith('accuracy_50bp.tsv'):
                    df['excluded_bp'] = 50
                elif accuracy_file.endswith('accuracy_150bp.tsv'):
                    df['excluded_bp'] = 150
                elif accuracy_file.endswith('accuracy_500bp.tsv'):
                    df['excluded_bp'] = 500
                else:
                    df['excluded_bp'] = -1  # legacy/unknown
            
            # Extract experiment info from file path
            path_parts = Path(accuracy_file).parts
            if len(path_parts) >= 8:
                eid = path_parts[2]
                pan_stem = path_parts[3] 
                tag = path_parts[4]
                reads_info = path_parts[7]  # cov{cov}_{n}_rep{rep}
                
                # Parse cov{cov}_{n}_rep{rep}
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
                    coverage_str = str(int(coverage)) if coverage == int(coverage) else str(coverage)
                    meta_file = f"workflow_output/experiments/{eid}/{pan_stem}/{tag}/results/reads/cov{coverage_str}_{reads}_rep{replicate}.txt"
                    mutation_breakdown = "unknown"
                    calculated_mutation_rate = 0.0
                    calculated_genome_size = 0
                    try:
                        if Path(meta_file).exists():
                            with open(meta_file, 'r') as mf:
                                meta_line = mf.read().strip()
                                all_parts = re.split(r'[\t\s]+', meta_line)
                                fields = dict(part.split('=',1) for part in all_parts if '=' in part)
                                
                                meta_mutation_rate = float(fields.get('mutation_rate', 0.0))
                                meta_genome_size = int(fields.get('genome_size', 0))
                                
                                if meta_mutation_rate > 0 and meta_genome_size > 0:
                                    calculated_snps = int(meta_mutation_rate * meta_genome_size)
                                    calculated_indels = int(meta_mutation_rate * meta_genome_size / 4)
                                    calculated_mutation_rate = meta_mutation_rate
                                    calculated_genome_size = meta_genome_size
                                    
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
                        print(f"Could not parse metadata {meta_file}: {e}")
                        mutation_breakdown = "unknown"
                    
                    df['mutation_breakdown'] = mutation_breakdown
                    df['calculated_mutation_rate'] = calculated_mutation_rate
                    df['calculated_genome_size'] = calculated_genome_size
                    
                    all_data.append(df)
                    
        except Exception as e:
            print(f"Error loading {accuracy_file}: {e}")
            continue
    
    if len(all_data) == 0:
        print("No valid data loaded")
        return
    
    # Combine all data
    df = pd.concat(all_data, ignore_index=True)
    print(f"Combined {len(df)} records from {len(all_data)} files")
    print(f"Columns: {list(df.columns)}")
    
    # Calculate distance from expected genome
    if 'total_variants' in df.columns:
        df['distance'] = df['total_variants']
    else:
        df['distance'] = df.get('snps', 0) + df.get('indels', 0)
    
    # Count replicates for each combination to add to legend
    if 'calculated_mutation_rate' in df.columns and 'reads' in df.columns and 'replicate' in df.columns:
        replicate_counts = df.groupby(['reads', 'calculated_mutation_rate'])['replicate'].nunique().reset_index()
        replicate_counts.columns = ['reads', 'calculated_mutation_rate', 'n_replicates']
        df = df.merge(replicate_counts, on=['reads', 'calculated_mutation_rate'], how='left')
    elif 'mutation_breakdown' in df.columns and 'reads' in df.columns and 'replicate' in df.columns:
        replicate_counts = df.groupby(['reads', 'mutation_breakdown'])['replicate'].nunique().reset_index()
        replicate_counts.columns = ['reads', 'mutation_breakdown', 'n_replicates']
        df = df.merge(replicate_counts, on=['reads', 'mutation_breakdown'], how='left')
    else:
        if 'replicate' in df.columns and 'reads' in df.columns:
            replicate_counts = df.groupby(['reads'])['replicate'].nunique().reset_index()
            replicate_counts.columns = ['reads', 'n_replicates']
            df = df.merge(replicate_counts, on=['reads'], how='left')
        else:
            df['n_replicates'] = 1
    
    # Create read count + mutation rate categories for legend with replicate count
    if 'calculated_mutation_rate' in df.columns and 'reads' in df.columns:
        df['mutation_rate_str'] = df['calculated_mutation_rate'].apply(
            lambda x: f"{x:.0e}" if x > 0 and x < 0.001 else f"{x:.4f}" if x > 0 else "0"
        )
        df['category'] = df['reads'].astype(str) + ' reads, μ=' + df['mutation_rate_str'].astype(str) + ' (n=' + df['n_replicates'].astype(str) + ')'
        df['reads_sort'] = df['reads']
        df['mutations_sort'] = df['calculated_mutation_rate'] * 10000
    elif 'mutation_breakdown' in df.columns and 'reads' in df.columns:
        df['category'] = df['reads'].astype(str) + ' reads, ' + df['mutation_breakdown'].astype(str) + ' (n=' + df['n_replicates'].astype(str) + ')'
        df['reads_sort'] = df['reads']
        df['mutations_sort'] = df.get('mutation_rate', 0) * 10000
    else:
        df['category'] = df.get('reads', 'unknown').astype(str) + ' reads (n=' + df['n_replicates'].astype(str) + ')'
        df['reads_sort'] = df.get('reads', 0)
        df['mutations_sort'] = 0
    
    # Save raw combined data
    raw_data_file = plots_dir / "combined_accuracy_data.tsv"
    df.to_csv(raw_data_file, sep='\t', index=False)
    print(f"Saved raw data: {raw_data_file}")
    
    # Create simple text report
    report_file = plots_dir / "placement_accuracy_report.txt"
    with open(report_file, 'w') as f:
        f.write("Placement Accuracy Analysis Report\n")
        f.write("==================================\n\n")
        f.write(f"Total accuracy files processed: {len(all_data)}\n")
        f.write(f"Total records: {len(df)}\n")
        f.write(f"Categories found: {len(df['category'].unique())}\n\n")
        
        metrics = df['metric'].unique() if 'metric' in df.columns else ['all']
        f.write(f"Metrics analyzed: {', '.join(metrics)}\n\n")
        
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
    
    print(f"Created report: {report_file}")
    
    # Create visual plots with matplotlib
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
        
        plt.style.use('default')
        sns.set_palette("viridis")
        
        # Get unique (k,s) parameter combinations
        k_s_combinations = df['k_s_params'].unique() if 'k_s_params' in df.columns else ['all']
        
        for k_s_params in k_s_combinations:
            print(f"Creating plots for {k_s_params}")
            
            # Filter data for this (k,s) combination
            if k_s_params != 'all':
                k_s_data = df[df['k_s_params'] == k_s_params].copy()
            else:
                k_s_data = df.copy()
            
            if len(k_s_data) == 0:
                continue
            
            # Create plots for each metric
            metrics = k_s_data['metric'].unique() if 'metric' in k_s_data.columns else ['combined']
            
            for metric in metrics:
                print(f"Creating plot for {k_s_params}, metric: {metric}")
                
                # Filter data for this metric
                if metric != 'combined':
                    metric_data = k_s_data[k_s_data['metric'] == metric].copy()
                else:
                    metric_data = k_s_data.copy()
                
                if len(metric_data) == 0:
                    continue
                
                # Create the plot
                fig, ax = plt.subplots(figsize=(12, 8))
                
                # Sort categories by reads first, then by mutations
                if 'reads_sort' in metric_data.columns and 'mutations_sort' in metric_data.columns:
                    category_order = metric_data[['category', 'reads_sort', 'mutations_sort']].drop_duplicates()
                    category_order = category_order.sort_values(['reads_sort', 'mutations_sort'])
                    categories = category_order['category'].tolist()
                else:
                    categories = sorted(metric_data['category'].unique())
                    
                max_distance = min(int(metric_data['distance'].max()), 10) if len(metric_data) > 0 else 5
                
                # Create distance bins (0, 1, 2, 3, 4, >=5)
                distance_bins = list(range(min(max_distance + 1, 6)))
                if max_distance >= 5:
                    distance_bins.append('>5')
                
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
                    
                    # Count samples with distance > 5 if applicable
                    if max_distance >= 5:
                        count_high = len(cat_data[cat_data['distance'] > 5])
                        distance_counts.append(count_high / total_samples)
                    
                    plot_data.append(distance_counts)
                
                if len(plot_data) == 0:
                    continue
                
                # Create the bar plot
                x = np.arange(len(distance_bins))
                width = 0.8 / len(categories) if len(categories) > 0 else 0.8
                colors = plt.cm.viridis(np.linspace(1, 0, len(categories)))
                
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
                
                # Save the plot
                plot_file = plots_dir / f"placement_accuracy_{k_s_params}_{metric}.png"
                plt.savefig(plot_file, dpi=300, bbox_inches='tight')
                print(f"Saved plot: {plot_file}")
                
                # Also save as PDF
                plot_file_pdf = plots_dir / f"placement_accuracy_{k_s_params}_{metric}.pdf"
                plt.savefig(plot_file_pdf, bbox_inches='tight')
                
                plt.close()
        
        print(f"Created visual plots for {len(k_s_combinations)} (k,s) combinations")
        
    except ImportError:
        print("Note: Install matplotlib/seaborn for visual plots")
    
    print(f"Completed analysis of {len(df)} records")
    
    # Create touch file to indicate completion
    touch_file = plots_dir / "placement_accuracy_plots.done"
    touch_file.touch()
    print(f"Created completion marker: {touch_file}")

if __name__ == "__main__":
    main()
