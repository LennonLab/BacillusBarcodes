from __future__ import division
import os
import pandas as pd
import numpy as np
import  matplotlib.pyplot as plt

path = os.path.expanduser("~/GitHub/BacillusBarcodes")

bc1 = 'TCAGACGTGTGA'
bc2 = 'CAGCCTAGTACG'

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

        x1 = np.exp( (np.log(f1) - np.log(0.5)) / time ) - 1
        x2 = np.exp( (np.log(f2) - np.log(0.5)) / time ) - 1

        print(f1)

        mean_x = (x1*f1) + (x2*f2)

        if 'F2' in name:
            fitness = x1
        else:
            fitness = x2


        if name not in fit_dic:
            fit_dic[name] = [fitness]
        else:
            fit_dic[name].append(fitness)

mean_fit_dict = {}

for key, value in fit_dic.items():
    fit_mean = np.mean(value)
    flask = key.split('-')[1]
    print(flask)
    if flask not in mean_fit_dict:
        mean_fit_dict[flask] = [fit_mean]
    else:
        mean_fit_dict[flask].append(fit_mean)

print(mean_fit_dict)
