#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=5
#SBATCH --time=200:00:00

# mamba activate snps
# for CHR in $(cat Chromosomes.list); do sbatch -J FST_${CHR} ~/merondun/cuculus_host/gene_hunt/1A.Subset_Groups_AllinOne_FST.sh ${CHR}; done
CHR=$1

#genotypes BELOW this will be set to missing
MINDP=3

if [[ $CHR = 'chr_W' || $CHR = 'chr_MT' ]]; then
        PLOIDY=1

        echo "FILTERING AND MERGING VARIANT SITES FOR ${CHR}, PLOIDY: ${PLOIDY}"
        bcftools view --threads 5 --samples-file AllSamples.list --force-samples -Ou ../../../merged/snps_only/${CHR}.SNPS.vcf.gz | \
                bcftools view --types snps --min-alleles 2 --max-alleles 2 --threads 5 | \
                #set genotypes below MINDP to missing
                bcftools +setGT -- -t q -i "FMT/DP < ${MINDP}" -n "./." | \
                #set het genotypes to missing based on binomial test
                bcftools +setGT -Ou -- -t "b:AD<1e-5" -n "./." | \
                #set weakly het genotypes to major allele
                bcftools +setGT -Ou -- --target-gt q --new-gt M -i 'GT=="het"' | \
                #set to haploid, can skip this for most purposes
                bcftools +fixploidy -Ou - -- -f ${PLOIDY} | \
                #update AC fields
                bcftools +fill-tags -Ou -- -t AC,AN | \
                bcftools view --min-alleles 2 --max-alleles 2 --types snps --min-ac 1 -i 'MQ > 40 & F_MISSING < 0.1' -Oz -o vcfs/${CHR}.SNP.DP3.vcf.gz
        bcftools index --threads 5 vcfs/${CHR}.SNP.DP3.vcf.gz

else
        PLOIDY=2

        echo "FILTERING AND MERGING VARIANT SITES FOR ${CHR}, PLOIDY: ${PLOIDY}"
        bcftools view --threads 5 --samples-file AllSamples.list -Ou ../../../merged/snps_only/${CHR}.SNPS.vcf.gz | \
                bcftools view --types snps --min-alleles 2 --max-alleles 2 --threads 5 | \
                #set genotypes below MINDP to missing
                bcftools +setGT -- -t q -i "FMT/DP < ${MINDP}" -n "./." | \
                #update AC fields
                bcftools +fill-tags -Ou -- -t AC,AN | \
                bcftools view --min-alleles 2 --max-alleles 2 --types snps --min-ac 1 -i 'MQ > 40 & F_MISSING < 0.1' -Oz -o vcfs/${CHR}.SNP.DP3.vcf.gz
        bcftools index --threads 5 vcfs/${CHR}.SNP.DP3.vcf.gz

fi

