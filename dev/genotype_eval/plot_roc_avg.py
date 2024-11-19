import os
import numpy as np
import scipy.interpolate
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
import json
from collections import defaultdict
import sys

def plotPRC(ax, prc_data_wp, prc_data_wo, coverage):
    aucs_wp = []
    aucs_wo = []

    # Define common recall points for shading (used only for variance visualization)
    recall_points = np.linspace(0, 1, num=100)

    def compute_mean_std(prc_data):
        interpolated_precisions = []
        for recall, precision in prc_data:
            recall = np.array(recall)
            precision = np.array(precision)
            if len(recall) == 0 or len(precision) == 0:
                interpolated_precisions.append(np.full_like(recall_points, np.nan))
            else:
                # Sort recall and precision in ascending order of recall
                sorted_indices = np.argsort(recall)
                recall_sorted = recall[sorted_indices]
                precision_sorted = precision[sorted_indices]
                # Interpolate precision at common recall points
                precision_interp = np.interp(recall_points, recall_sorted, precision_sorted, left=np.nan, right=np.nan)
                interpolated_precisions.append(precision_interp)
        interpolated_precisions = np.array(interpolated_precisions)
        mean_precision = np.nanmean(interpolated_precisions, axis=0)
        std_precision = np.nanstd(interpolated_precisions, axis=0)
        return mean_precision, std_precision

    def calculate_auc(prc_data):
        aucs = []
        for recall, precision in prc_data:
            recall = np.array(recall)
            precision = np.array(precision)
            if len(recall) == 0 or len(precision) == 0:
                aucs.append(0)
            else:
                sorted_indices = np.argsort(recall)
                recall_sorted = recall[sorted_indices]
                precision_sorted = precision[sorted_indices]
                auc = np.trapz(precision_sorted, recall_sorted)
                aucs.append(auc)
        return aucs

    # Compute mean and std precision for 'with prior' and 'without prior'
    mean_precision_wp, std_precision_wp = compute_mean_std(prc_data_wp)
    mean_precision_wo, std_precision_wo = compute_mean_std(prc_data_wo)

    # Calculate AUC for each replicate
    aucs_wp = calculate_auc(prc_data_wp)
    aucs_wo = calculate_auc(prc_data_wo)

    # Plot mean precision-recall curve with shaded regions
    ax.plot(recall_points, mean_precision_wp, color='blue', label='Tree-specific prior')
    ax.fill_between(recall_points, mean_precision_wp - std_precision_wp, mean_precision_wp + std_precision_wp, color='blue', alpha=0.3)
    ax.plot(recall_points, mean_precision_wo, color='orange', label='BCFtools flat prior')
    ax.fill_between(recall_points, mean_precision_wo - std_precision_wo, mean_precision_wo + std_precision_wo, color='orange', alpha=0.3)

    # Compute mean and std of AUCs
    mean_auc_wp = np.mean(aucs_wp)
    std_auc_wp = np.std(aucs_wp)
    mean_auc_wo = np.mean(aucs_wo)
    std_auc_wo = np.std(aucs_wo)

    # Display the AUCs on the plot
    textstr = '\n'.join((
        r'Tree-specific prior: Mean AUC = %.3f, STD = %.3f' % (mean_auc_wp, std_auc_wp),
        r'BCFtools flat prior: Mean AUC = %.3f, STD = %.3f' % (mean_auc_wo, std_auc_wo)
    ))

    ax.text(0.95, 0.05, textstr, fontsize=8,
            verticalalignment='bottom', horizontalalignment='right',
            transform=ax.transAxes, bbox=dict(facecolor='white', alpha=0.3))
    ax.set_xlim(-0.1, 1.1)
    ax.set_ylim(-0.1, 1.1)
    ax.set_title(f'{coverage}X Coverage')
    ax.set_xlabel('Recall')
    ax.set_ylabel('Precision')
    ax.grid(True)

    # Return mean and std of AUCs
    return mean_auc_wp, std_auc_wp, mean_auc_wo, std_auc_wo

def get_stats_from_roc(ref_name, prefix, cov, rep, scaling_factor, prior):
    # dir_path = f'{prefix}_roc_bcf_vcf'
    if prior:
        dir_path = f'{prefix}_roc_panmap'
    else:
        dir_path = f'{prefix}_roc_bcf_vcf'
    TPRs = []
    PRCs = []
    
    file_path = ''
    if prior:
        # file_path = f'{dir_path}/{ref_name}.{rep}.{cov}x.mutmat_s{scaling_factor}.roc'
        file_path = f'{dir_path}/{ref_name}.{rep}.{cov}x.roc'
    else:
        file_path = f'{dir_path}/{ref_name}.{rep}.{cov}x.varonly.roc'
    fh = open(file_path, 'r')
    for line in fh.readlines()[1:]:
        fields = list(map(float, line.strip().split()))
        TP = (fields[1] - fields[2]) + (fields[4] - fields[5]) + (fields[7] - fields[8])
        FP = fields[2] + fields[5] + fields[8]
        FN = fields[3] + fields[6] + fields[9]
        assert(TP + FN > 0)
        TPR = round(TP / (TP + FN), 3)
        PRC = round(TP / (TP + FP), 3) if (TP + FP > 0) else 1
        TPRs.append(TPR)
        PRCs.append(PRC)
    fh.close()
    return TPRs, PRCs

############################
## Precision Recall Curve ##
############################
ref_name = sys.argv[1]
prefix = sys.argv[2]
num_reps = int(sys.argv[3])
scaling_factor = float(sys.argv[4])
# Ensure scaling_factor is an integer if it is a whole number
if scaling_factor.is_integer():
    scaling_factor = int(scaling_factor)

specific_rep = int(sys.argv[5]) if len(sys.argv) > 5 else None

# Assuming there are 10 coverage levels (adjust nrows and ncols as needed)
fig, axs = plt.subplots(nrows=2, ncols=5, figsize=(20, 8)) # Adjust the layout as necessary
axs = axs.flatten() # Flatten the array of axes for easy indexing
coverages = [i+1 for i in range(10)]
mean_aucs_msp = []
std_aucs_msp = []
mean_aucs_bcf = []
std_aucs_bcf = []
reps = [specific_rep] if specific_rep is not None else range(num_reps)
for i, cov in enumerate(coverages):
    prc_data_msp = []
    prc_data_bcf = []
    for rep in reps:
        TPRs_msp, PRCs_msp = get_stats_from_roc(ref_name, prefix, cov, rep, scaling_factor, prior=True)
        TPRs_bcf, PRCs_bcf = get_stats_from_roc(ref_name, prefix, cov, rep, scaling_factor, prior=False)
        prc_data_msp.append((TPRs_msp, PRCs_msp))
        prc_data_bcf.append((TPRs_bcf, PRCs_bcf))

    mean_auc_msp, std_auc_msp, mean_auc_bcf, std_auc_bcf = plotPRC(axs[i], prc_data_msp, prc_data_bcf, cov)
    mean_aucs_msp.append(mean_auc_msp)
    std_aucs_msp.append(std_auc_msp)
    mean_aucs_bcf.append(mean_auc_bcf)
    std_aucs_bcf.append(std_auc_bcf)

# Creating legend labels
line1, = axs[0].plot([], [], color='blue', label='Tree-specific prior')
line2, = axs[0].plot([], [], color='orange', label='BCFtools flat prior')

# Adding a figure-level legend
fig.legend(handles=[line1, line2], loc='lower center', ncol=2)

# Adjust layout to make room for the figure-level legend
plt.tight_layout(rect=[0, 0.05, 1, 0.95]) 

fig.suptitle('Precision-Recall Curves for Different Read Coverages, SNP (MiSeq, bwa)', fontsize=16)
plt.savefig(f'{prefix}_precision_recall_curve_s{scaling_factor}.png')

coverages = np.array(coverages)  # Convert coverages to numpy array if it isn't already

# Create the plot for Mean AUC comparison
plt.figure(figsize=(10, 6))
plt.errorbar(coverages, mean_aucs_msp, yerr=std_aucs_msp, fmt='o-', color='blue', label='MSP', capsize=5)
plt.errorbar(coverages, mean_aucs_bcf, yerr=std_aucs_bcf, fmt='s-', color='orange', label='BCF', capsize=5)

# Add labels, legend, and grid
plt.xlabel('Coverage')
plt.xticks(coverages, [f'{i}X' for i in coverages])
plt.ylabel('Mean AUC')
plt.ylim(0, 1)
plt.title('Comparison of Mean AUCs Between MSP and BCF Across Coverages')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.7)

# Display the plot
plt.tight_layout()
plt.savefig(f'{prefix}_mean_auc_comparison_s{scaling_factor}.png')