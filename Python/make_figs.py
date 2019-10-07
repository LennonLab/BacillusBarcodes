from __future__ import division
import os
import pandas as pd
import numpy as np
import  matplotlib.pyplot as plt
from matplotlib.lines import Line2D

path = os.path.expanduser("~/GitHub/BacillusBarcodes")

bc1 = 'TCAGACGTGTGA'
bc2 = 'CAGCCTAGTACG'

bc1_freq_0 = {'GSF2346-F1-A':0.502632411,'GSF2346-F1-B':0.515677676,
            'GSF2346-F1-C':0.59174401, 'GSF2346-F2-A':0.415243398,
            'GSF2346-F2-B':0.43089695, 'GSF2346-F2-C':0.486277276,
            'GSF2346-F3-A':0.425032152, 'GSF2346-F3-B':0.459793851,
            'GSF2346-F3-C':0.49423082}




# assume two fold dilution
# log_2 (100) ~= 6.64
time = 6.64
# assume f_0 = 0.5
# fitness = e ^ ((ln(f_t) - ln(f_0))/ 6.64 ) - 1

# do this for each barcode, calculate mean, get relative fitness for strain of interest


#F1 – BC1-WT & BC2-R1
#F2 – BC2-WT & BC1-R1
#F3 – BC1-WT & BC2-WT

fit_dic = {}
with open(path + '/data/barcode_counts.txt', 'r') as fh:
    for line in fh:
        line_split = line.strip().split('\t')
        if line_split[0] == 'fastq_file':
            continue
        name = line_split[0].split('_')[0][:-1]
        bc1 = int(line_split[1])
        bc2 = int(line_split[2])

        f1 = bc1 / (bc1+bc2)
        f2 = 1 - f1

        print(name, f1)

        f1_0 = bc1_freq_0[name]
        f2_0 = 1- f1_0

        x1 = np.exp( (np.log(f1) - np.log(f1_0)) / time ) - 1
        x2 = np.exp( (np.log(f2) - np.log(f2_0)) / time ) - 1
        mean_x = (x1*f1) + (x2*f2)

        #if 'F2' in name:
        #    fitness = x1
        #else:
        #    fitness = x2

        name_bc1 = name + '-bc1'
        name_bc2 = name + '-bc2'

        if name_bc1 not in fit_dic:
            fit_dic[name_bc1] = [x1]
        else:
            fit_dic[name_bc1].append(x1)

        if name_bc2 not in fit_dic:
            fit_dic[name_bc2] = [x2]
        else:
            fit_dic[name_bc2].append(x2)


mean_fit_dict = {}
colow_dict = []
for key, value in fit_dic.items():
    key_split = key.split('-')
    new_key = key_split[1] + '-' + key_split[3]
    if new_key not in mean_fit_dict:
        mean_fit_dict[new_key] = [np.mean(value)]
    else:
        mean_fit_dict[new_key].append(np.mean(value))

fig = plt.figure()
plt.axhline(y=0, color = 'dimgrey', lw = 2, ls = '--', zorder=1)
for key, value in mean_fit_dict.items():

    i = int(key.split('-')[0][1])
    if 'bc1' in key:
        c = 'b'
    else:
        c = 'r'
    x = np.random.normal(i, 0.04, size=len(value))
    plt.scatter(x, value, color = c, alpha=0.8, zorder=2)


x1 = [1,1.5,2,2.5,3]
squad = ['BC1-WT vs. BC2-R1','','BC2-WT vs BC1-R1','', 'BC1-WT vs BC2-WT']

plt.xticks(x1, squad)#, rotation=45)
plt.ylabel('Fitness', fontsize=18)


legend_elements = [Line2D([0], [0], marker='o', color='b', label='Barcode 1',
                        markersize=10, linestyle="None"),
                   Line2D([0], [0], marker='o', color='r', label='Barcode 2',
                        markersize=10, linestyle="None")]

plt.legend(handles=legend_elements, loc='upper right')


fig_name = path + '/figs/fitness.png'
fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
