#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00

# mamba activate snps
#submit as: for RUN in $(cat CNV_Samples.list); do sbatch -J COV_${RUN} ~/merondun/cuculus_host/gene_hunt/2.Calculate_Coverage_for_log2FM.sh ${RUN}; done
RUN=$1
CHR=$2

mkdir work raw
bamdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/bams/Illumina_Alignments_Merged

#calculate coverage over 500BP  windows
mosdepth --threads 1 --chrom ${CHR} --no-per-base work/${RUN}_${CHR}_500 --mapq 20 --by 500 --use-median ${bamdir}/${RUN}.bam

zcat work/${RUN}_${CHR}_500*gz | awk -v s="${RUN}" '{OFS="\t"}{print $0, s}' > raw/${RUN}_${CHR}_500MQ20.cov