#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

def create_comprehensive_grid_plot():
    """Create the comprehensive grid plot with all k,s combinations and metrics"""
    
    output_dir = Path("workflow_output")
    plots_dir = output_dir / "plots"
    
    # Load the combined data that was already created
    combined_data_file = plots_dir / "combined_accuracy_data.tsv"
    
    if not combined_data_file.exists():
        print("No combined data file found. Run the main plotting script first.")
        return
    
    print("Loading combined accuracy data...")
    df = pd.read_csv(combined_data_file, sep='\t')
    
    print(f"Loaded {len(df)} records")
    
    # Calculate distance from expected genome (total variants)
    if 'total_variants' in df.columns:
        df['distance'] = df['total_variants']
    else:
        df['distance'] = df.get('snps', 0) + df.get('indels', 0)
    
    # Get unique (k,s) parameter combinations and metrics
    k_s_combinations = sorted(df['k_s_params'].unique()) if 'k_s_params' in df.columns else ['all']
    metrics = ['raw', 'jaccard', 'cosine', 'weighted']
    exclusion_levels = sorted(df['excluded_bp'].unique()) if 'excluded_bp' in df.columns else [0]
    
    print(f"Creating comprehensive grid for {len(k_s_combinations)} k,s combinations, {len(metrics)} metrics, {len(exclusion_levels)} exclusion levels")
    
    # PERFORMANCE EVALUATION: Find best performing metric/k-s combinations
    print("\n" + "="*60)
    print("PERFORMANCE EVALUATION")
    print("="*60)
    
    # Calculate average distance for each metric/k-s combination
    performance_data = []
    
    for k_s_params in k_s_combinations:
        k_s_data = df[df['k_s_params'] == k_s_params]
        for excl in exclusion_levels:
            excl_data = k_s_data[k_s_data.get('excluded_bp',0) == excl]
            for metric in metrics:
                metric_data = excl_data[excl_data['metric'] == metric]
                if len(metric_data) > 0:
                    avg_distance = metric_data['distance'].mean()
                    median_distance = metric_data['distance'].median()
                    std_distance = metric_data['distance'].std()
                    perfect_matches = len(metric_data[metric_data['distance'] == 0])
                    total_samples = len(metric_data)
                    perfect_match_rate = perfect_matches / total_samples * 100
                    performance_data.append({
                        'k_s_params': k_s_params,
                        'metric': metric,
                        'excluded_bp': excl,
                        'avg_distance': avg_distance,
                        'median_distance': median_distance,
                        'std_distance': std_distance,
                        'perfect_matches': perfect_matches,
                        'total_samples': total_samples,
                        'perfect_match_rate': perfect_match_rate
                    })
    
    # Create performance dataframe and sort by average distance
    perf_df = pd.DataFrame(performance_data)
    if 'excluded_bp' in perf_df.columns:
        print("\nTOP PERFORMERS BY EXCLUSION LEVEL:")
        for excl in exclusion_levels:
            sub = perf_df[perf_df['excluded_bp']==excl]
            if len(sub) == 0:
                continue
            best = sub.sort_values(['avg_distance','median_distance']).head(3)
            print(f"  Excluded {excl}bp:")
            for _, row in best.iterrows():
                print(f"    {row['k_s_params']} {row['metric']} avg_dist={row['avg_distance']:.3f} perfect%={row['perfect_match_rate']:.1f}")
    perf_df_sorted = perf_df.sort_values('avg_distance')
    
    print(f"\nBEST PERFORMING COMBINATIONS (by lowest average distance):")
    print("-" * 80)
    print(f"{'Rank':<4} {'K-S Params':<12} {'Metric':<10} {'Avg Dist':<9} {'Med Dist':<9} {'Perfect %':<9} {'Samples':<8}")
    print("-" * 80)
    
    for idx, row in perf_df_sorted.head(10).iterrows():
        print(f"{idx+1:<4} {row['k_s_params']:<12} {row['metric']:<10} "
              f"{row['avg_distance']:<9.3f} {row['median_distance']:<9.1f} "
              f"{row['perfect_match_rate']:<9.1f} {row['total_samples']:<8.0f}")
    
    print("\nWORST PERFORMING COMBINATIONS (by highest average distance):")
    print("-" * 80)
    print(f"{'Rank':<4} {'K-S Params':<12} {'Metric':<10} {'Avg Dist':<9} {'Med Dist':<9} {'Perfect %':<9} {'Samples':<8}")
    print("-" * 80)
    
    for idx, row in perf_df_sorted.tail(10).iterrows():
        rank = len(perf_df_sorted) - len(perf_df_sorted.tail(10)) + (idx - perf_df_sorted.tail(10).index[0]) + 1
        print(f"{rank:<4} {row['k_s_params']:<12} {row['metric']:<10} "
              f"{row['avg_distance']:<9.3f} {row['median_distance']:<9.1f} "
              f"{row['perfect_match_rate']:<9.1f} {row['total_samples']:<8.0f}")
    
    # Summary by metric
    print(f"\nPERFORMANCE SUMMARY BY METRIC:")
    print("-" * 50)
    metric_summary = perf_df.groupby('metric').agg({
        'avg_distance': ['mean', 'std', 'min', 'max'],
        'perfect_match_rate': ['mean', 'std', 'min', 'max']
    }).round(3)
    
    for metric in metrics:
        metric_data = perf_df[perf_df['metric'] == metric]
        if len(metric_data) > 0:
            print(f"\n{metric.upper()} metric:")
            print(f"  Average distance: {metric_data['avg_distance'].mean():.3f} ± {metric_data['avg_distance'].std():.3f}")
            print(f"  Best performance: {metric_data['avg_distance'].min():.3f} (k_s: {metric_data.loc[metric_data['avg_distance'].idxmin(), 'k_s_params']})")
            print(f"  Perfect match rate: {metric_data['perfect_match_rate'].mean():.1f}% ± {metric_data['perfect_match_rate'].std():.1f}%")
    
    # Summary by k-s parameters
    print(f"\nPERFORMANCE SUMMARY BY K-S PARAMETERS:")
    print("-" * 50)
    
    for k_s in k_s_combinations:
        k_s_data = perf_df[perf_df['k_s_params'] == k_s]
        if len(k_s_data) > 0:
            best_metric = k_s_data.loc[k_s_data['avg_distance'].idxmin(), 'metric']
            best_distance = k_s_data['avg_distance'].min()
            print(f"\n{k_s}:")
            print(f"  Best metric: {best_metric} (avg distance: {best_distance:.3f})")
            print(f"  Average across metrics: {k_s_data['avg_distance'].mean():.3f} ± {k_s_data['avg_distance'].std():.3f}")
            print(f"  Perfect match rate: {k_s_data['perfect_match_rate'].mean():.1f}% ± {k_s_data['perfect_match_rate'].std():.1f}%")
    
    # Save detailed performance report
    performance_file = plots_dir / "performance_evaluation.tsv"
    perf_df_sorted.to_csv(performance_file, sep='\t', index=False)
    print(f"\nDetailed performance data saved to: {performance_file}")
    
    print("="*60)
    print("")
    
    n_k_s_combinations = len(k_s_combinations)
    n_metrics = len(metrics)
    
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
    
    print("Creating individual subplots...")
    
    # Plot each combination
    for row_idx, k_s_params in enumerate(k_s_combinations):
        print(f"  Processing {k_s_params} (row {row_idx + 1}/{n_rows})")
        k_s_data = df[df['k_s_params'] == k_s_params]
        
        for col_idx, metric in enumerate(metrics):
            ax = axes[row_idx][col_idx]
            
            metric_data = k_s_data[k_s_data['metric'] == metric].copy()
            
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
    
    print("Adding title and legend...")
    
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
    
    print("Saving comprehensive plots...")
    
    # Save comprehensive PDF
    comprehensive_pdf = plots_dir / "placement_accuracy_comprehensive.pdf"
    plt.savefig(comprehensive_pdf, bbox_inches='tight', dpi=300)
    print(f"Saved comprehensive PDF: {comprehensive_pdf}")
    
    # Also save as PNG for quick viewing
    comprehensive_png = plots_dir / "placement_accuracy_comprehensive.png"
    plt.savefig(comprehensive_png, bbox_inches='tight', dpi=300)
    print(f"Saved comprehensive PNG: {comprehensive_png}")
    
    plt.close()
    
    # CREATE PERFORMANCE VISUALIZATION
    print("Creating performance visualization...")
    
    # Create performance heatmap
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    # Pivot data for heatmap - Average Distance
    perf_pivot_avg = perf_df.pivot(index='k_s_params', columns='metric', values='avg_distance')
    
    # Plot average distance heatmap
    sns.heatmap(perf_pivot_avg, annot=True, fmt='.3f', cmap='RdYlGn_r', 
                ax=ax1, cbar_kws={'label': 'Average Distance'})
    ax1.set_title('Average Distance from Target Genome\n(Lower is Better)', fontweight='bold')
    ax1.set_xlabel('Metric')
    ax1.set_ylabel('K-S Parameters')
    
    # Pivot data for heatmap - Perfect Match Rate
    perf_pivot_perfect = perf_df.pivot(index='k_s_params', columns='metric', values='perfect_match_rate')
    
    # Plot perfect match rate heatmap
    sns.heatmap(perf_pivot_perfect, annot=True, fmt='.1f', cmap='RdYlGn', 
                ax=ax2, cbar_kws={'label': 'Perfect Match Rate (%)'})
    ax2.set_title('Perfect Match Rate\n(Higher is Better)', fontweight='bold')
    ax2.set_xlabel('Metric')
    ax2.set_ylabel('K-S Parameters')
    
    plt.tight_layout()
    
    # Save performance visualization
    performance_viz_pdf = plots_dir / "performance_evaluation.pdf"
    plt.savefig(performance_viz_pdf, bbox_inches='tight', dpi=300)
    print(f"Saved performance visualization PDF: {performance_viz_pdf}")
    
    performance_viz_png = plots_dir / "performance_evaluation.png"
    plt.savefig(performance_viz_png, bbox_inches='tight', dpi=300)
    print(f"Saved performance visualization PNG: {performance_viz_png}")
    
    plt.close()
    
    # Create bar plot showing top 10 best performers
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # Get top 10 best performers
    top_10 = perf_df_sorted.head(10).copy()
    top_10['label'] = top_10['k_s_params'] + '_' + top_10['metric']
    
    # Create bar plot
    bars = ax.bar(range(len(top_10)), top_10['avg_distance'], 
                  color=plt.cm.viridis(np.linspace(0, 1, len(top_10))))
    
    # Customize plot
    ax.set_xlabel('K-S Parameter & Metric Combination')
    ax.set_ylabel('Average Distance from Target')
    ax.set_title('Top 10 Best Performing Combinations\n(by Average Distance from Target Genome)', 
                fontweight='bold', fontsize=14)
    ax.set_xticks(range(len(top_10)))
    ax.set_xticklabels(top_10['label'], rotation=45, ha='right')
    ax.grid(True, alpha=0.3)
    
    # Add value labels on bars
    for i, (bar, value) in enumerate(zip(bars, top_10['avg_distance'])):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                f'{value:.3f}', ha='center', va='bottom', fontsize=10)
    
    plt.tight_layout()
    
    # Save top performers plot
    top_performers_pdf = plots_dir / "top_performers.pdf"
    plt.savefig(top_performers_pdf, bbox_inches='tight', dpi=300)
    print(f"Saved top performers PDF: {top_performers_pdf}")
    
    top_performers_png = plots_dir / "top_performers.png"
    plt.savefig(top_performers_png, bbox_inches='tight', dpi=300)
    print(f"Saved top performers PNG: {top_performers_png}")
    
    plt.close()
    
    # CREATE PERFORMANCE LINE PLOT WITH METRIC × CATEGORY COMBINATIONS
    print("\nCreating performance line plot with metric × category combinations...")
    
    # Calculate average distance for each metric × category × k,s combination
    line_plot_data = []
    all_categories = sorted(df['category'].unique())
    
    # Exclude the 0mut category with only n=2 samples
    categories = [cat for cat in all_categories if '0mut' not in cat]
    print(f"Excluded categories with 0mut: {[cat for cat in all_categories if '0mut' in cat]}")
    print(f"Using {len(categories)} categories for line plot")
    
    metrics = ['raw', 'jaccard', 'cosine', 'weighted']
    
    # Create expanded categories that include metric information
    for category in categories:
        for metric in metrics:
            # Filter data for this category and metric
            cat_metric_data = df[(df['category'] == category) & (df['metric'] == metric)]
            
            if len(cat_metric_data) > 0:
                for k_s_params in k_s_combinations:
                    k_s_cat_metric_data = cat_metric_data[cat_metric_data['k_s_params'] == k_s_params]
                    if len(k_s_cat_metric_data) > 0:
                        avg_distance = k_s_cat_metric_data['distance'].mean()
                        
                        # Create combined category label with metric
                        combined_category = f"{category} ({metric})"
                        
                        line_plot_data.append({
                            'category': category,
                            'metric': metric,
                            'combined_category': combined_category,
                            'k_s_params': k_s_params,
                            'avg_distance': avg_distance,
                            'n_samples': len(k_s_cat_metric_data)
                        })
    
    line_df = pd.DataFrame(line_plot_data)
    
    if len(line_df) > 0:
        # Create the performance line plot
        fig, ax = plt.subplots(figsize=(16, 10))
        
        # Sort k,s combinations for proper ordering on x-axis
        # Extract k and s values for sorting
        k_s_sort_data = []
        for k_s in k_s_combinations:
            if k_s != 'all':
                # Extract k and s from k{k}_s{s} format
                import re
                k_match = re.search(r'k(\d+)', k_s)
                s_match = re.search(r's(\d+)', k_s)
                if k_match and s_match:
                    k_val = int(k_match.group(1))
                    s_val = int(s_match.group(1))
                    k_s_sort_data.append((k_s, k_val, s_val))
        
        # Sort by k first, then by s
        k_s_sort_data.sort(key=lambda x: (x[1], x[2]))
        k_s_ordered = [x[0] for x in k_s_sort_data]
        
        # Create x-axis positions
        x_positions = range(len(k_s_ordered))
        x_labels = [f"k{x[1]},s{x[2]}" for x in k_s_sort_data]
        
        # Get unique combined categories and sort them
        combined_categories = sorted(line_df['combined_category'].unique())
        
        # Create color palette - use different base colors for each original category
        # and different line styles for each metric
        base_categories = sorted(df['category'].unique())
        metric_styles = {'raw': '-', 'jaccard': '--', 'cosine': '-.', 'weighted': ':'}
        metric_markers = {'raw': 'o', 'jaccard': 's', 'cosine': '^', 'weighted': 'D'}
        
        # Generate colors for base categories
        base_colors = plt.cm.tab10(np.linspace(0, 1, len(base_categories)))
        
        # Plot lines for each metric × category combination
        for combined_cat in combined_categories:
            # Extract original category and metric
            for base_cat_idx, base_cat in enumerate(base_categories):
                for metric in metrics:
                    expected_combined = f"{base_cat} ({metric})"
                    if combined_cat == expected_combined:
                        # Get data for this combination
                        cat_metric_data = line_df[line_df['combined_category'] == combined_cat]
                        
                        # Get y values for this combination in the correct k,s order
                        y_values = []
                        for k_s in k_s_ordered:
                            k_s_data = cat_metric_data[cat_metric_data['k_s_params'] == k_s]
                            if len(k_s_data) > 0:
                                y_values.append(k_s_data['avg_distance'].iloc[0])
                            else:
                                y_values.append(np.nan)  # Missing data point
                        
                        # Plot the line, handling missing data
                        mask = ~np.isnan(y_values)
                        if np.any(mask):
                            x_plot = np.array(x_positions)[mask]
                            y_plot = np.array(y_values)[mask]
                            
                            ax.plot(x_plot, y_plot, 
                                   linestyle=metric_styles[metric],
                                   marker=metric_markers[metric], 
                                   linewidth=2, markersize=4, 
                                   label=combined_cat, 
                                   color=base_colors[base_cat_idx], 
                                   alpha=0.8)
        
        # Customize the plot
        ax.set_xlabel('(k,s) Parameter Combinations', fontsize=12)
        ax.set_ylabel('Average Distance from Target Genome', fontsize=12)
        ax.set_title('Performance by (k,s) Parameters, Read/Mutation Categories, and Metrics', 
                    fontsize=14, fontweight='bold')
        
        # Set x-axis
        ax.set_xticks(x_positions)
        ax.set_xticklabels(x_labels, rotation=45, ha='right')
        
        # Add grid
        ax.grid(True, alpha=0.3)
        
        # Add legend with smaller font to accommodate more entries
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8, ncol=1)
        
        # Add horizontal line at y=0 for reference (perfect placement)
        ax.axhline(y=0, color='red', linestyle='-', alpha=0.5, linewidth=1)
        
        # Add text box explaining line styles
        metric_legend_text = "Line styles: raw (—), jaccard (- -), cosine (-·), weighted (···)"
        ax.text(0.02, 0.98, metric_legend_text, transform=ax.transAxes, 
               fontsize=10, verticalalignment='top', 
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        
        plt.tight_layout()
        
        # Save performance line plot
        performance_line_pdf = plots_dir / "performance_by_k_s_parameters_with_metrics.pdf"
        plt.savefig(performance_line_pdf, bbox_inches='tight', dpi=300)
        print(f"Saved performance line plot PDF: {performance_line_pdf}")
        
        performance_line_png = plots_dir / "performance_by_k_s_parameters_with_metrics.png"
        plt.savefig(performance_line_png, bbox_inches='tight', dpi=300)
        print(f"Saved performance line plot PNG: {performance_line_png}")
        
        plt.close()
        
        # Print summary of best performing combinations for each category × metric
        print("\nBest performing (k,s) for each category × metric combination:")
        print("-" * 70)
        for combined_cat in sorted(combined_categories):
            cat_metric_data = line_df[line_df['combined_category'] == combined_cat]
            if len(cat_metric_data) > 0:
                best_row = cat_metric_data.loc[cat_metric_data['avg_distance'].idxmin()]
                print(f"{combined_cat}: {best_row['k_s_params']} (avg distance: {best_row['avg_distance']:.3f})")
    
    print("\nPerformance evaluation completed!")
    print("Comprehensive grid plot completed!")

if __name__ == "__main__":
    create_comprehensive_grid_plot()
