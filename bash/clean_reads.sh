#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=10gb,walltime=4:00:00
#PBS -M wrshoema@iu.edu
#PBS -m abe
#PBS -j oe


module load cutadapt


for file in "/N/dc2/projects/muri2/Task2/BacillusBarcodes/data/fastq/GSF2346-"*"_R1_001.fastq.gz";
do
  file_clean="${file/.fastq.gz/_clean.fastq.gz}"
  cutadapt -q 30 --minimum-length 150 -o $file_clean $file
done
