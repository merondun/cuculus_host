#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=12:00:00

SCRATCH=tmp/$SLURM_JOB_ID

#Bam cleaning script. Merge by sample, Mark duplicates, remove reads overhanging scaffolds, and summarize

RUN=$1

RGdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/trimmed_fastq/RGs
bamdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/bams/2021_04
genome=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa

mkdir tmp
mkdir $SCRATCH

#Make final bam and stats folder
mkdir ${bamdir}/merged
mkdir ${bamdir}/merged/stats
mkdir ${bamdir}/merged/stats/alignments
mkdir ${bamdir}/merged/stats/coverage

#merge all the libraries together
samtools merge - ${bamdir}/${RUN}*.bam | samtools sort -@3 -o $SCRATCH/${RUN}.bam
samtools index $SCRATCH/${RUN}.bam

sambamba markdup --tmpdir $SCRATCH/ $SCRATCH/${RUN}.bam $SCRATCH/${RUN}.dup.bam
samtools index $SCRATCH/${RUN}.dup.bam

#Remove reads over scaffold end
gatk CleanSam -I $SCRATCH/${RUN}.dup.bam -O ${bamdir}/merged/${RUN}.bam -R ${genome} --CREATE_INDEX true

Summary Stats
gatk CollectAlignmentSummaryMetrics -R ${genome} -I ${bamdir}/merged/${RUN}.bam -O ${bamdir}/merged/stats/alignments/${RUN}.alignments.txt -LEVEL SAMPLE -LEVEL READ_GROUP

#calculate coverage
mosdepth --threads 3 --use-median --by 100000 --fast-mode --no-per-base ${bamdir}/merged/stats/coverage/${RUN} ${bamdir}/merged/${RUN}.bam
