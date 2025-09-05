#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import re
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

def parse_paf_file(paf_file):
    """Parse a PAF file and extract alignment differences by position"""
    mutations_by_position = []
    
    try:
        with open(paf_file, 'r') as f:
            for line in f:
                if line.strip():
                    fields = line.strip().split('\t')
                    if len(fields) >= 12:
                        # PAF format: query_name, query_len, query_start, query_end, strand, 
                        # target_name, target_len, target_start, target_end, matches, alignment_len, mapq
                        query_name = fields[0]
                        target_name = fields[5]
                        target_start = int(fields[7])
                        target_end = int(fields[8])
                        
                        # Define exclusion zones (first and last 100bp of alignment)
                        exclusion_size = 100
                        alignment_start_exclude = target_start + exclusion_size
                        alignment_end_exclude = target_end - exclusion_size
                        
                        # Skip alignments that are too short (less than 200bp total)
                        if target_end - target_start < 2 * exclusion_size:
                            continue
                        
                        # Look for cs string (coordinate string) in optional fields
                        cs_string = None
                        cigar_string = None
                        for field in fields[12:]:
                            if field.startswith('cs:Z:'):
                                cs_string = field[5:]  # Remove 'cs:Z:' prefix
                            elif field.startswith('cg:Z:'):
                                cigar_string = field[5:]  # Remove 'cg:Z:' prefix
                        
                        if cs_string:
                            # Parse cs string to extract mutation positions
                            mutations = parse_cs_string(cs_string, target_start, alignment_start_exclude, alignment_end_exclude)
                            mutations_by_position.extend(mutations)
                        elif cigar_string:
                            # Fallback to parsing CIGAR if no cs string
                            mutations = parse_cigar_string(cigar_string, target_start, alignment_start_exclude, alignment_end_exclude)
                            mutations_by_position.extend(mutations)
    
    except Exception as e:
        print(f"Error parsing {paf_file}: {e}")
    
    return mutations_by_position

def parse_cs_string(cs_string, start_pos, exclude_start, exclude_end):
    """Parse minimap2 cs string to extract mutation positions, excluding first and last 100bp"""
    mutations = []
    current_pos = start_pos
    
    # Split cs string into operations
    # cs format: :6*ag:3+gat:10-cc:5
    # : = match, * = substitution, + = insertion, - = deletion
    
    i = 0
    while i < len(cs_string):
        if cs_string[i] == ':':
            # Match - advance position by the number
            i += 1
            num_str = ''
            while i < len(cs_string) and cs_string[i].isdigit():
                num_str += cs_string[i]
                i += 1
            if num_str:
                current_pos += int(num_str)
        
        elif cs_string[i] == '*':
            # Substitution - record mutation at current position if within allowed region
            if exclude_start <= current_pos < exclude_end:
                mutations.append({
                    'type': 'substitution',
                    'position': current_pos,
                    'reference': cs_string[i+1] if i+1 < len(cs_string) else '',
                    'alternate': cs_string[i+2] if i+2 < len(cs_string) else ''
                })
            current_pos += 1
            i += 3  # Skip *XY
        
        elif cs_string[i] == '+':
            # Insertion - record at current position if within allowed region
            i += 1
            insert_seq = ''
            while i < len(cs_string) and cs_string[i] not in ':*+-':
                insert_seq += cs_string[i]
                i += 1
            if exclude_start <= current_pos < exclude_end:
                mutations.append({
                    'type': 'insertion',
                    'position': current_pos,
                    'sequence': insert_seq
                })
        
        elif cs_string[i] == '-':
            # Deletion - record at current position and advance if within allowed region
            i += 1
            delete_seq = ''
            while i < len(cs_string) and cs_string[i] not in ':*+-':
                delete_seq += cs_string[i]
                i += 1
            if exclude_start <= current_pos < exclude_end:
                mutations.append({
                    'type': 'deletion',
                    'position': current_pos,
                    'sequence': delete_seq
                })
            current_pos += len(delete_seq)
        
        else:
            i += 1
    
    return mutations

def parse_cigar_string(cigar_string, start_pos, exclude_start, exclude_end):
    """Parse CIGAR string to extract approximate mutation positions, excluding first and last 100bp"""
    mutations = []
    current_pos = start_pos
    
    # Parse CIGAR: 8518M17I6736M
    cigar_ops = re.findall(r'(\d+)([MIDNSHPX=])', cigar_string)
    
    for length, op in cigar_ops:
        length = int(length)
        
        if op == 'M':  # Match/mismatch
            current_pos += length
        elif op == 'I':  # Insertion
            # Record insertion(s) at current position if within allowed region
            if exclude_start <= current_pos < exclude_end:
                for _ in range(length):
                    mutations.append({
                        'type': 'insertion',
                        'position': current_pos
                    })
        elif op == 'D':  # Deletion
            # Record deletion(s) and advance position if within allowed region
            for i in range(length):
                pos = current_pos + i
                if exclude_start <= pos < exclude_end:
                    mutations.append({
                        'type': 'deletion',
                        'position': pos
                    })
            current_pos += length
        elif op in 'SH':  # Soft/hard clipping
            pass  # Don't advance position for clipping
        elif op == 'N':  # Skipped region
            current_pos += length
    
    return mutations

def create_genomic_mutation_plot():
    """Create plot showing distribution of alignment differences across genome"""
    
    output_dir = Path("workflow_output")
    plots_dir = output_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    
    # Find all PAF files
    paf_files = list(output_dir.glob("**/alignments/**/logs/*.paf"))
    
    print(f"Found {len(paf_files)} PAF files to analyze")
    
    if len(paf_files) == 0:
        print("No PAF files found!")
        return
    
    # Parse all PAF files and collect mutations by position
    all_mutations = []
    processed_files = 0
    
    for paf_file in paf_files[:100]:  # Limit to first 100 files for faster processing
        print(f"Processing {paf_file.name}...")
        mutations = parse_paf_file(paf_file)
        all_mutations.extend(mutations)
        processed_files += 1
        
        if processed_files % 10 == 0:
            print(f"Processed {processed_files} files, found {len(all_mutations)} mutations so far")
    
    print(f"Total mutations found: {len(all_mutations)}")
    
    if len(all_mutations) == 0:
        print("No mutations found in PAF files!")
        return
    
    # Convert to DataFrame
    mutations_df = pd.DataFrame(all_mutations)
    
    # Determine genome size (approximate from max position)
    max_position = mutations_df['position'].max()
    genome_size = max_position + 1000  # Add some buffer
    
    print(f"Estimated genome size: {genome_size} bp")
    
    # Create position bins (0-25, 26-50, etc.)
    bin_size = 25
    num_bins = int(np.ceil(genome_size / bin_size))
    bins = np.arange(0, num_bins * bin_size + 1, bin_size)
    
    # Assign each mutation to a bin
    mutations_df['bin'] = pd.cut(mutations_df['position'], bins=bins, 
                                labels=[f"{i*bin_size}-{(i+1)*bin_size-1}" for i in range(num_bins)],
                                include_lowest=True)
    
    # Count mutations per bin
    mutations_per_bin = mutations_df.groupby('bin').size().reset_index()
    mutations_per_bin.columns = ['bin', 'mutation_count']
    
    # Create bin centers for plotting
    bin_centers = []
    bin_labels = []
    for i in range(len(mutations_per_bin)):
        bin_label = mutations_per_bin.iloc[i]['bin']
        start_pos = i * bin_size + bin_size/2  # Center of bin
        bin_centers.append(start_pos)
        bin_labels.append(bin_label)
    
    mutations_per_bin['bin_center'] = bin_centers
    
    # Create the plot
    fig, ax = plt.subplots(figsize=(16, 8))
    
    # Plot points with lines connecting them
    ax.plot(mutations_per_bin['bin_center'], mutations_per_bin['mutation_count'], 
            marker='o', linewidth=2, markersize=6, alpha=0.8, color='steelblue')
    
    # Fill area under the curve
    ax.fill_between(mutations_per_bin['bin_center'], mutations_per_bin['mutation_count'], 
                   alpha=0.3, color='steelblue')
    
    # Customize the plot
    ax.set_xlabel('Genomic Position (bp)', fontsize=12)
    ax.set_ylabel('Number of Alignment Differences', fontsize=12)
    ax.set_title('Distribution of Alignment Differences Across Genome\n(25bp bins, excluding first/last 100bp of alignments)', 
                fontsize=14, fontweight='bold')
    
    # Add grid
    ax.grid(True, alpha=0.3)
    
    # Format x-axis to show every nth bin label
    n_labels = min(20, len(bin_centers))  # Show at most 20 labels
    step = max(1, len(bin_centers) // n_labels)
    ax.set_xticks(bin_centers[::step])
    ax.set_xticklabels([f"{int(x)}" for x in bin_centers[::step]], rotation=45)
    
    plt.tight_layout()
    
    # Save the plot
    plot_file = plots_dir / "genomic_mutation_distribution.pdf"
    plt.savefig(plot_file, bbox_inches='tight', dpi=300)
    print(f"Saved plot: {plot_file}")
    
    plot_file_png = plots_dir / "genomic_mutation_distribution.png"
    plt.savefig(plot_file_png, bbox_inches='tight', dpi=300)
    print(f"Saved plot: {plot_file_png}")
    
    plt.close()
    
    # Also create a plot by mutation type
    if 'type' in mutations_df.columns:
        fig, ax = plt.subplots(figsize=(16, 8))
        
        # Count mutations by type and bin
        type_counts = mutations_df.groupby(['bin', 'type']).size().reset_index()
        type_counts.columns = ['bin', 'type', 'count']
        
        # Pivot to get mutation types as columns
        type_pivot = type_counts.pivot(index='bin', columns='type', values='count').fillna(0)
        
        # Add bin centers
        type_pivot = type_pivot.reset_index()
        type_pivot['bin_center'] = [i * bin_size + bin_size/2 for i in range(len(type_pivot))]
        
        # Plot each mutation type
        mutation_types = ['substitution', 'insertion', 'deletion']
        colors = ['red', 'green', 'blue']
        
        for mut_type, color in zip(mutation_types, colors):
            if mut_type in type_pivot.columns:
                ax.plot(type_pivot['bin_center'], type_pivot[mut_type], 
                       marker='o', linewidth=2, markersize=4, alpha=0.8, 
                       color=color, label=mut_type.capitalize())
        
        # Customize the plot
        ax.set_xlabel('Genomic Position (bp)', fontsize=12)
        ax.set_ylabel('Number of Mutations', fontsize=12)
        ax.set_title('Distribution of Mutation Types Across Genome\n(25bp bins, excluding first/last 100bp of alignments)', 
                    fontsize=14, fontweight='bold')
        
        # Add legend and grid
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Format x-axis
        ax.set_xticks(type_pivot['bin_center'][::step])
        ax.set_xticklabels([f"{int(x)}" for x in type_pivot['bin_center'][::step]], rotation=45)
        
        plt.tight_layout()
        
        # Save the mutation type plot
        type_plot_file = plots_dir / "genomic_mutation_types_distribution.pdf"
        plt.savefig(type_plot_file, bbox_inches='tight', dpi=300)
        print(f"Saved mutation types plot: {type_plot_file}")
        
        type_plot_file_png = plots_dir / "genomic_mutation_types_distribution.png"
        plt.savefig(type_plot_file_png, bbox_inches='tight', dpi=300)
        print(f"Saved mutation types plot: {type_plot_file_png}")
        
        plt.close()
    
    # Save summary data
    summary_file = plots_dir / "genomic_mutation_summary.tsv"
    mutations_per_bin.to_csv(summary_file, sep='\t', index=False)
    print(f"Saved summary data: {summary_file}")
    
    # Print summary statistics
    print(f"\nSummary Statistics:")
    print(f"Total mutations analyzed: {len(all_mutations)}")
    print(f"Mutations per 25bp bin - Mean: {mutations_per_bin['mutation_count'].mean():.2f}")
    print(f"Mutations per 25bp bin - Std: {mutations_per_bin['mutation_count'].std():.2f}")
    print(f"Mutations per 25bp bin - Max: {mutations_per_bin['mutation_count'].max()}")
    print(f"Number of genomic bins: {len(mutations_per_bin)}")
    
    if 'type' in mutations_df.columns:
        print(f"\nMutation type distribution:")
        type_counts = mutations_df['type'].value_counts()
        for mut_type, count in type_counts.items():
            print(f"  {mut_type}: {count} ({count/len(mutations_df)*100:.1f}%)")

if __name__ == "__main__":
    create_genomic_mutation_plot()
