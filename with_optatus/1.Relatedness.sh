#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00

GROUP=$1

for CHR in $(cat somalier/MacroChrs.list); do 
bcftools view --threads 10 --force-samples --samples-file ${GROUP}.list ../../merged/${CHR}.vcf.gz -Ob | \
        bcftools view --types snps --threads 10 --min-af 0.3 --max-af 0.7 --min-alleles 2 --max-alleles 2 -i 'F_MISSING<0.1' -Ob | \
        bcftools +prune -m 0.1 -w 10kb --nsites-per-win 1 -o related/${CHR}_${GROUP}_MAF50-LDr1.vcf.gz
bcftools index --threads 10 related/${CHR}_${GROUP}_MAF50-LDr1.vcf.gz
bcftools index -n --threads 10 related/${CHR}_${GROUP}_MAF50-LDr1.vcf.gz

done 

bcftools concat related/*${GROUP}*MAF50-LDr1*gz -Oz -o related/${GROUP}_MAF50-LDr1_Macros.vcf.gz
bcftools index related/${GROUP}_MAF50-LDr1_Macros.vcf.gz

genome=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa

cd somalier
mkdir ${GROUP}
somalier extract -d ${GROUP} --sites ../related/${GROUP}_MAF50-LDr1_Macros.vcf.gz -f ${genome} ../related/${GROUP}_MAF50-LDr1_Macros.vcf.gz
somalier relate --min-depth 5 ${GROUP}/*.somalier -o ${GROUP}