#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --time=48:00:00

bam=$1

#this was first assessed from external study, this SNP file (of unpruned SNPs) was further filtered to leave approx 50K SNPs 
#bcftools +prune -m 0.1 -w 30kb --nsites-per-win 1  autos_canorus_LD.vcf.gz -Ob | bcftools view -G -Oz -o ../with_optatus/somalier/Sites.vcf.gz
genome=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa

somalier extract -d extracted/ --sites Sites.vcf.gz -f ${genome} ${bam}

#afterwards 
somalier relate --min-depth 5 extracted/*.somalier