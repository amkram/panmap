from collections import defaultdict
import matplotlib.patches as mplpatches
import matplotlib.pyplot as plt
import numpy as np
import json
import os


##########
## Data ##
##########
'''
groups:
    haplotype N:
        prop 0 (others):
            5x: 0.05, 0.02
            10x: 0 
        prop 0.5: 
            5x: 0.51, 0.45, 0.52...
            10x: ...
        prop 0.2: ...
        prop 0.1: ...
        ..
'''

'''
files:
    1 haplotype:
        replicate 0:
            true abundance: 
            estimated: 
                5x cov
                10x cov
        replicate 1:
        ...
    3 haplotypes:
        ...
'''
files = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
for file in os.listdir('sim_results_pure_reads/'):
    fields = file.split('_')
    n_haplotype = int(fields[0][1:])
    replicate = int(fields[1][1:])
    if file.endswith('abundance.txt'):
        files[n_haplotype][replicate]['true'] = [file]
    else:
        files[n_haplotype][replicate]['estimates'].append(file)

groups = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
for n_haplotype in (1, 3, 5, 10):
    for replicate in range(10):
        true_abundance = {}
        with open('sim_results_pure_reads/' + files[n_haplotype][replicate]['true'][0]) as fh:
            for line in fh.readlines():
                fields = line.strip().split()
                true_abundance[fields[0]] = float(fields[1])
        for out_abundance_file in files[n_haplotype][replicate]['estimates']:
            coverage = (int(out_abundance_file.split('_')[2]) / 2) / 200
            estimate_abundance = {}
            with open('sim_results_pure_reads/' + out_abundance_file) as fh:
                for line in fh.readlines():
                    if line.startswith('@'): continue
                    fields = line.strip().split()
                    estimate_abundance[fields[0]] = float(fields[1])
            
            seen = {}
            for hap in true_abundance:
                seen[hap] = True
                if hap in estimate_abundance:
                    groups[n_haplotype][true_abundance[hap]][coverage].append(estimate_abundance[hap])
                else:
                    undetected = True
                    for haplotypes_str in estimate_abundance:
                        if hap in haplotypes_str:
                            groups[n_haplotype][true_abundance[hap]][coverage].append(estimate_abundance[haplotypes_str])
                            undetected = False
                            seen[haplotypes_str] = True
                            break
                    if undetected:
                        if n_haplotype == 10 and coverage == 200 and true_abundance[hap] == 0.15: print(replicate)
                        groups[n_haplotype][true_abundance[hap]][coverage].append(0)
            
            others = 0
            for hap in estimate_abundance:
                if hap in seen: continue
                others += estimate_abundance[hap]
            groups[n_haplotype][0][coverage].append(others)
            
quartiles_data = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
for n_haplotype in groups:
    for true_abun in groups[n_haplotype]:
        for cov in groups[n_haplotype][true_abun]:
            quartiles_data[n_haplotype][true_abun][cov] = list(np.quantile(groups[n_haplotype][true_abun][cov], [0, 0.25, 0.5, 0.75, 1.0]))

# print(json.dumps(quartiles_data, indent=2))

###############
## Functions ##
###############
def make_box(panel, xcoor, width, quartiles, facecolor):
    rectangle = mplpatches.Rectangle(
        (xcoor, quartiles[1]), width, quartiles[3] - quartiles[1],
        facecolor=facecolor,
        edgecolor='black',
        linewidth=0.3,
        alpha=0.5
    )

    # box
    panel.add_patch(rectangle)

    # middle line
    panel.plot(
        [xcoor, xcoor + width],
        [quartiles[2], quartiles[2]],
        color='black',
        linewidth=0.3
    )

    # whiskers
    panel.plot(
        [xcoor + width / 2, xcoor + width / 2],
        [quartiles[3], quartiles[4]],
        color='black',
        linewidth=0.3
    )
    panel.plot(
        [xcoor + width / 2, xcoor + width / 2],
        [quartiles[0], quartiles[1]],
        color='black',
        linewidth=0.3
    )


#################
## figure size ##
#################

# Define style
plt.style.use('BME163')

figure_size = (6.5, 5)
plt.figure(figsize=figure_size)

##################
## Legend panel ##
##################
colors = ['red', 'blue', 'purple', 'yellow', 'orange']
legend_panel_size = (1.5, 0.25)
legend_panel_coor = (6.5 / 2 - 1.5 / 2, 0.2)
legend_panel = plt.axes([
   legend_panel_coor[0] / figure_size[0],
   legend_panel_coor[1] / figure_size[1],
   legend_panel_size[0] / figure_size[0],
   legend_panel_size[1] / figure_size[1]],
   frameon=False)
legend_panel.tick_params(
    left=False, labelleft=False)
legend_panel.set_xlim(-0.2, 1.0)
legend_panel.set_ylim(-0.2, 1.0)
legend_gap = 0.05
legend_box_width = (1 - legend_gap * 4) / 5
legend_panel.set_xticks([(legend_box_width + legend_gap) * i + legend_box_width / 2 for i in range(5)])
legend_panel.set_xticklabels(['5x', '10x', '50x', '100x', '200x'], fontsize=6)
legend_panel.tick_params(axis='x', which='both', length=0)
legend_panel.text(
    -0.1, 0.5, 'Average depth', fontsize=6,
    horizontalalignment='right',
    verticalalignment='center')
legned_quartiles = [0, 0.25 , 0.5, 0.75, 1.0]
for i, cov in enumerate((5.0, 10.0, 50.0, 100.0, 200.0)):
    cur_legend_box_xcoor = i * (legend_box_width + legend_gap)
    make_box(legend_panel, cur_legend_box_xcoor, legend_box_width, legned_quartiles, colors[i])


#######################
## 1 haplotype panel ##
#######################

hap_1_x_size = 1.2
hap_1_panel_size = (hap_1_x_size, 1.5)
hap_1_panel_coor = (0.5, 3)
hap_1_panel = plt.axes([
    hap_1_panel_coor[0] / figure_size[0],
    hap_1_panel_coor[1] / figure_size[1],
    hap_1_panel_size[0] / figure_size[0],
    hap_1_panel_size[1] / figure_size[1]])
hap_1_panel.set_xlim(0, 4)
hap_1_panel.set_ylim(-0.1, 1.1)
hap_1_panel.set_xticks([1, 3])
hap_1_panel.set_xticklabels(['haplotype', 'others'])
hap_1_panel.set_yticks([0.2 * i for i in range(6)])
hap_1_panel.set_ylabel('Abundance')
hap_1_panel.set_title(
    '1 haplotype (1.0)',
    loc='left', fontsize=9)
plt.grid(color = 'gray', linewidth = 0.2, alpha=0.3)


#######################
## 3 haplotype panel ##
#######################

hap_3_panel_size = (hap_1_x_size * 2, 1.5)
hap_3_panel_coor = (0.5, 0.8)
hap_3_panel = plt.axes([
    hap_3_panel_coor[0] / figure_size[0],
    hap_3_panel_coor[1] / figure_size[1],
    hap_3_panel_size[0] / figure_size[0],
    hap_3_panel_size[1] / figure_size[1]])
hap_3_panel.set_xlim(0, 8)
hap_3_panel.set_ylim(-0.1, 0.65)
hap_3_panel.set_xticks([1, 3, 5, 7])
hap_3_panel.set_xticklabels(['', 'haplotypes', '', 'others'])
hap_3_panel.set_yticks([0.1 * i for i in range(7)])
hap_3_panel.set_ylabel('Abundance')
hap_3_panel.set_title(
    '3 haplotypes (0.5, 0.3, 0.2)',
    loc='left', fontsize=9)
plt.grid(color = 'gray', linewidth = 0.2, alpha=0.3)

#######################
## 5 haplotype panel ##
#######################

hap_5_panel_size = (hap_1_x_size * 3, 1.5)
hap_5_panel_coor = (2.5, 3)
hap_5_panel = plt.axes([
    hap_5_panel_coor[0] / figure_size[0],
    hap_5_panel_coor[1] / figure_size[1],
    hap_5_panel_size[0] / figure_size[0],
    hap_5_panel_size[1] / figure_size[1]])
hap_5_panel.set_xlim(0, 12)
hap_5_panel.set_ylim(-0.1, 0.65)
hap_5_panel.set_xticks([1, 3, 5, 7, 9, 11])
hap_5_panel.set_xticklabels(['', '', 'haplotypes', '', '', 'others'])
hap_5_panel.set_yticks([0.1 * i for i in range(7)])
hap_5_panel.set_ylabel('Abundance')
hap_5_panel.set_title(
    '5 haplotypes (0.5, 0.2, 0.15, 0.1, 0.05)',
    loc='left', fontsize=9)
plt.grid(color = 'gray', linewidth = 0.2, alpha=0.3)

########################
## 10 haplotype panel ##
########################

hap_10_panel_size = (hap_1_x_size * 2, 1.5)
hap_10_panel_coor = (3.7, 0.8)
hap_10_panel = plt.axes([
    hap_10_panel_coor[0] / figure_size[0],
    hap_10_panel_coor[1] / figure_size[1],
    hap_10_panel_size[0] / figure_size[0],
    hap_10_panel_size[1] / figure_size[1]])
hap_10_panel.set_xlim(0, 8)
hap_10_panel.set_ylim(-0.1, 0.65)
hap_10_panel.set_xticks([1, 3, 5, 7])
hap_10_panel.set_xticklabels(['', 'haplotypes', '', 'others'])
hap_10_panel.set_yticks([0.1 * i for i in range(7)])
hap_10_panel.set_ylabel('Abundance')
hap_10_panel.set_title(
    '10 haplotypes (4 x 0.15, 2 x 0.1, 4 x 0.05)',
    loc='left', fontsize=9)
plt.grid(color = 'gray', linewidth = 0.2, alpha=0.3)

##############
## Whiskeys ##
##############
group_range = 1.8
box_width = group_range / 5
offset = 0.1
true_abundances = {
    1: [1.0, 0],
    3: [0.5, 0.3, 0.2, 0],
    5: [0.5, 0.2, 0.15, 0.1, 0.05, 0],
    10: [0.15, 0.1, 0.05, 0]
}
num_haplotypes = [1, 3, 5, 10]
panels = [hap_1_panel, hap_3_panel, hap_5_panel, hap_10_panel]
for num_haplotype, panel in zip(num_haplotypes, panels):
    for i, abun in enumerate(true_abundances[num_haplotype]):
        panel.plot(
            [offset + 2 * i, offset + 2 * i + group_range],
            [abun, abun],
            color='red',
            linestyle='--',
            linewidth=0.5
        )
        for j, cov in enumerate((5.0, 10.0, 50.0, 100.0, 200.0)):
            xcoor = offset + i * 2 + j * box_width
            points    = groups[num_haplotype][abun][cov]
            quartiles = quartiles_data[num_haplotype][abun][cov]
            panel.plot(
                [xcoor + box_width / 2 for i in range(len(points))],
                points,
                marker='o',
                markersize=1,
                color='blue',
                alpha=0.7,
                markeredgewidth=0,
                linewidth=0)
            make_box(panel, xcoor, box_width, quartiles, colors[j])

plt.savefig('sim_results_pure_reads.png', dpi=1200)