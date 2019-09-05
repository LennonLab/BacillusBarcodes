from __future__ import division
import sys
import os
import gzip
from itertools import groupby

bc1 = 'TCAGACGTGTGA'
bc2 = 'CAGCCTAGTACG'

directory = '/N/dc2/projects/muri2/Task2/BacillusBarcodes/data/fastq'

quality_dict = {'!':0, '"':1, '#':2, '$':3, '%':4, '&':5, "'":6, '(':7, ')':8,
                '*':9, '#':10, ',':11, '-':12, '.':13, '/':14, '0':15, '1':16,
                '2':17, '3':18, '4':19, '5':20, '6':21, '7':22, '8':23, '9':24,
                ':':25, ';':26, '<':27, '=':28, '>':29, '?':30, '@':31, 'A':32,
                'B':33, 'C':34, 'D':35, 'E':36, 'F':37, 'G':38, 'H':39, 'I':40,
                'J':41, 'K':42}



def process(lines=None):
    ks = ['name', 'sequence', 'optional', 'quality']
    return {k: v for k, v in zip(ks, lines)}




def run_parse_fastq():
    df_out = open('/N/dc2/projects/muri2/Task2/BacillusBarcodes/data/barcode_counts.txt', 'w')
    df_out.write('\t'.join(['fastq_file', 'bc_1', 'bc_2', 'min_lib_bc', 'mean_lib_bc', 'max_lib_bc']) + '\n')

    for filename in os.listdir(directory):
        if filename.endswith("_001.fastq.gz"):
            file_no_ext = filename.split('.')[0]
            #sys.stderr.write("Sample %s analyzed" % (str(file_no_ext)))
            print(os.path.join(directory, filename), file_no_ext)
            library_bcs = []
            n=4
            n_bc1 = 0
            n_bc2 = 0
            n_no_match = 0
            with gzip.open(os.path.join(directory, filename), 'r') as fh:
                lines = []
                for line in fh:
                    lines.append(line.rstrip())
                    if len(lines) == n:
                        record = process(lines)
                        strain_bc = record['sequence'][34:34+12]
                        strain_bc = strain_bc.decode("utf-8")
                        strain_bc_quality = record['quality'][34:34+12]
                        mean_strain_bc_quality = sum([quality_dict[i] for i in list(strain_bc_quality.decode("utf-8"))]) / len(strain_bc_quality)
                        if mean_strain_bc_quality > 30:
                            if strain_bc == bc1:
                                n_bc1 += 1
                            elif strain_bc == bc2:
                                n_bc2 += 1
                            else:
                                n_no_match += 1
                        library_bc = record['sequence'][:8]
                        library_bcs.append(library_bc.decode("utf-8"))

                        lines = []


            counts = [(i, len(list(c))) for i,c in groupby(library_bcs)]
            just_counts = [i[1] for i in counts]

            df_out.write('\t'.join([file_no_ext, str(n_bc1), str(n_bc2), str(min(just_counts)), str(sum(just_counts)/len(just_counts)), str(max(just_counts))]) + '\n')




    df_out.close()


run_parse_fastq()
