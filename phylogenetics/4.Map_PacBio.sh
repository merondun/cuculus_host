#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=20
#SBATCH --time=220:00:00

# mamba activate pacbio
# for i in $(egrep -v 'Neomo|Taper' RUNS.list); do sbatch -J PACBIO_${i} ~/merondun/cuculus_host/phylogenetics/4.Map_PacBio.sh ${i}; done
RUN=$1 

genomedir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/
bamdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/bams/PacBio_Alignments
covdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/bams/coverage

#index
# minimap2 -t 20 -d $genomedir/GCA_017976375.1_bCucCan1.pri_genomic.CHR.mmi $genomedir/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa

#map and index 
ID=$(grep ${RUN} README.txt | awk '{print $1}')
minimap2 -t 20 -ax map-pb $genomedir/GCA_017976375.1_bCucCan1.pri_genomic.CHR.mmi ${RUN}_1.fastq.gz | samtools view -@20 -bS - | samtools sort -@ 20 -o $bamdir/${ID}.bam -
samtools index -@20 $bamdir/${ID}.bam

#calculate coverage over 5KB windows
mosdepth --threads 20 --no-per-base $covdir/${ID}_5KB --mapq 30 --by 5000 --use-median $bamdir/${ID}.bam
zcat $covdir/${ID}_5KB*gz | awk -v s="${ID}" '{OFS="\t"}{print $0, s}' > $covdir/${ID}_5KB-MQ30.cov