#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00

# Align cuckoo illumina data 

RUN=$1

echo "Aligning sample: ${RUN}"

SCRATCH=/tmp/$SLURM_JOB_ID
genome=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa
qcdata=/dss/dsslegfs01/pr53da/pr53da-dss-0021/shared_resources/trimmed_fastq
outdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/cuculiformes_NDUFAF4_blasting/Alignment_Based_Approach/Whole_Genome_Alignment/bams
subdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/cuculiformes_NDUFAF4_blasting/Alignment_Based_Approach/subset_5gb

# Subset 5gb, and then align, discarding unaligned reads to save space (-F 4)
#bbduk.sh in=${qcdata}/${RUN}.trim.fastq.gz out=${subdir}/${RUN}.5gb.fastq.gz maxbasesout=5000000000
bwa mem -M -p -t 10 ${genome} ${subdir}/${RUN}.5gb.fastq.gz | samtools sort -@10 -o ${SCRATCH}/${RUN}.bam -;
samtools view -F 4 -b ${SCRATCH}/${RUN}.bam > ${outdir}/${RUN}.bam
samtools index -b ${outdir}/${RUN}.bam;       