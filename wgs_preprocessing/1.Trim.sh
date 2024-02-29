#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=12:00:00

#Submit at library level 
RUN=$1
wd=/dss/dsslegfs01/pr53da/pr53da-dss-0021/rawdata/WGS/ILLUMINA/
outdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/trimmed_fastq

RUN=$1

file1=$(ls ${wd}/${RUN}__R1__*)
file2=$(ls ${wd}/${RUN}__R2__*)

SCRATCH=/tmp/$SLURM_JOB_ID

#adapter trim
bbduk.sh t=6 -Xmx24g in=${file1} in2=${file2} out=${outdir}/${RUN}.trim.fastq.gz minlen=25 qtrim=rl trimq=2 ktrim=r k=23 mink=11 ref=~/modules/bbmap/adapters.fa hdist=1 tpe tbo