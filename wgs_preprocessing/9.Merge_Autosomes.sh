#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=20
#SBATCH --time=200:00:00

mkdir merged_unrelfull
#merge autosomes
bcftools concat --file-list Autosomes.list --threads 20 -Oz -o merged_unrelfull/autos.vcf.gz
bcftools index --threads 20 merged_unrelfull/autos.vcf.gz

#LD prune
~/modules/plink2 --threads 20 --vcf merged_unrelfull/autos.vcf.gz --allow-extra-chr --set-missing-var-ids @:# \
        --rm-dup --indep-pairwise 50 5 0.1 --maf 0.05 --hwe 1e-10 --max-alleles 2 --min-alleles 2 --out merged_unrelfull/autos_canorus_LD
        
#extract, also a vcf and run PCA 
~/modules/plink2 --threads 20 --vcf merged_unrelfull/autos.vcf.gz --allow-extra-chr --set-missing-var-ids @:# \
        --extract merged_unrelfull/autos_canorus_LD.prune.in \
        --make-bed --recode vcf bgz --pca --out merged_unrelfull/autos_canorus_LD
bcftools index --threads 20 merged_unrelfull/autos_canorus_LD.vcf.gz
sed -i 's/chr_//g' merged_unrelfull/autos_canorus_LD.bim