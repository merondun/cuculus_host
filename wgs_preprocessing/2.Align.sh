#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=12:00:00

#example submit: for sample in $(cat Libraries.list); do sbatch -J ${sample} Align.sh ${sample}; done Plot Alignment:
SCRATCH=tmp/$SLURM_JOB_ID

wd=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/trimmed_fastq
qcdata=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/trimmed_fastq
outdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/bams/2021_04
genome=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa

RUN=$1

bwa index ${genome}
#Mapping script, will map reads with the prefix identified in the list ${lists}/${RUN}

        ID="$(cat ${wd}/RGs/${RUN}.ID | tail -1)";
        PU="$(cat ${wd}/RGs/${RUN}.PU | tail -1)";
        SM="$(cat ${wd}/RGs/${RUN}.SM | tail -1)";
        LB="$(cat ${wd}/RGs/${RUN}.LB | tail -1)";

        #library mapping
        bwa mem -M -p -t 10 -R "@RG\tID:${ID}\tSM:${SM}\tPL:ILLUMINA\tPU:${PU}\tLB:${LB}" ${genome} ${qcdata}/${RUN}.trim.fastq.gz | samtools sort -@10 -o ${SCRATCH}/${RUN}.bam -;
        samtools view -b -f 1 -F 524 ${SCRATCH}/${RUN}.bam > ${outdir}/${RUN}.bam
        samtools index -b ${outdir}/${RUN}.bam;