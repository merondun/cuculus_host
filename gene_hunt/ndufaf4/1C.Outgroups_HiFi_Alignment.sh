#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00

# Align pacbio HiFi data for outgroups 

RUN=$1

echo "Aligning sample: ${RUN}"

SCRATCH=/tmp/$SLURM_JOB_ID
genome=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa
outdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/cuculiformes_NDUFAF4_blasting/Alignment_Based_Approach/Whole_Genome_Alignment/10gb/bams
rawdata=/dss/dsslegfs01/pr53da/pr53da-dss-0021/rawdata/External_Data_SRA_GenomeArk/SRA_Cuculiformes_Outgroups
subdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/cuculiformes_NDUFAF4_blasting/Alignment_Based_Approach/subset_10gb

# Subset 10gb, and then align, discarding unaligned reads to save space (-F 4)
bbduk.sh in=${rawdata}/${RUN}_1.fastq.gz out=${subdir}/${RUN}.10gb.fastq.gz maxbasesout=10000000000
minimap2 -ax map-pb -t 8 ${genome} ${subdir}/${RUN}.10gb.fastq.gz > ${SCRATCH}/${RUN}.sam
samtools sort ${SCRATCH}/${RUN}.sam | samtools view -F 4 -b > ${outdir}/${RUN}.bam
samtools index -b ${outdir}/${RUN}.bam;  