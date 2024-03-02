#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=3
#SBATCH --time=48:00:00

# mamba activate snps 
#submit as: for RUN in $(awk '{print $1}' ~/merondun/cuculus_host/gene_hunt/Coverage_SamplesN2.txt | sed '1d'); do sbatch -J ~/merondun/cuculus_host/gene_hunt/2.Calculate_Coverage_for_log2FM.sh ${RUN}; done 
RUN=$1

mkdir work run 
bamdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/bams/Illumina_Alignments_Merged
#calculate coverage over 500bp windows
mosdepth --threads 3 --no-per-base work/${RUN}_500 --mapq 20 --by 500 --use-median ${bamdir}/${RUN}.bam

zcat work/${RUN}_500*gz | awk -v s="${RUN}" '{OFS="\t"}{print $0, s}' > raw/${RUN}_500bpMQ20.cov