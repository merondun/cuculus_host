#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=5
#SBATCH --time=48:00:00

CHR=$1
GENOME=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa
mask=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/Coverage_Masks/GCA_017976375.1_bCucCan1.pri_genomic.CHR.N75-DoubleCoverage.mask.bed

echo "WORKING ON ${CHR}"

# Use bcftools to perform sample-wise filtering on the vcf file.
bcftools view --types snps --force-samples --samples-file Full_Samples_n302.list -Oz ../merged/${CHR}.vcf.gz -o snp_full/${CHR}.vcf.gz
bcftools index --threads 5  snp_full/${CHR}.vcf.gz

# Calculate some statistics for the next round of filtering
SAMPN=$(bcftools query -l snp_full/${CHR}.vcf.gz | wc -l)
AVGDP=$(bcftools query -f '%DP\n' snp_full/${CHR}.vcf.gz | datamash median 1 | datamash round 1) # Compute the median coverage depth across all sites in the vcf file
DPHI=$(($AVGDP*2)) # Calculate thresholds for filtering based on depth of coverage

# Apply filters
bcftools view --types snps --threads 5 -e "QUAL < 20 || INFO/DP > $DPHI || INFO/DP < $SAMPN || MQ < 30 || RPBZ < -3 || RPBZ > 3" -Oz -o snp_full/${CHR}.IF.vcf.gz snp_full/${CHR}.vcf.gz
bcftools index --threads 5 snp_full/${CHR}.IF.vcf.gz

# For diploid chromosomes, convert low-coverage sites to missing
MINDP=3
bcftools +setGT -Oz -o snp_full/${CHR}.IF-GF.vcf.gz snp_full/${CHR}.IF.vcf.gz -- -t q -i "FMT/DP < ${MINDP}" -n "./."
bcftools index --threads 5 snp_full/${CHR}.IF-GF.vcf.gz

#for W and MT, set allele imbalance violations (heterozygosity) to missing for FEMALE, no W for MALE, diploid Z for MALE
if [[ $CHR == 'chr_MT' ]]
then
        echo "HAPLOID, SETTING HETEROZYGOUS SITES TO MISSING FOR BOTH SEXES"
        #All samples
        bcftools view --samples-file Full_Samples_n302.list -Ou snp_full/${CHR}.IF-GF.vcf.gz | \
                bcftools view --min-alleles 2 --max-alleles 2 --min-ac 1 --threads 5 -i 'F_MISSING<0.1' | \
                bedtools subtract -header -a - -b ${mask} | \
                bcftools +setGT -Ou -- -t "b:AD<1e-5" -n "./." | \
            bcftools +setGT -Oz -o snp_full/${CHR}.IF-GF-MM1.vcf.gz -- --target-gt q --new-gt M -i 'GT=="het"'
        bcftools index --threads 5 snp_full/${CHR}.IF-GF-MM1.vcf.gz

        echo "HAPLOID, SETTING HETEROZYGOUS SITES TO MISSING FOR BOTH SEXES"
        #Females only
        bcftools view --samples-file FFull_Samples_n302.list -Ou snp_full/${CHR}.IF-GF.vcf.gz | \
                bcftools view --min-alleles 2 --max-alleles 2 --min-ac 1 --threads 5 -i 'F_MISSING<0.1' | \
                bedtools subtract -header -a - -b ${mask} | \
                bcftools +setGT -Ou -- -t "b:AD<1e-5" -n "./." | \
            bcftools +setGT -Oz -o snp_full/${CHR}.F_IF-GF-MM1.vcf.gz -- --target-gt q --new-gt M -i 'GT=="het"'
        bcftools index --threads 5 snp_full/${CHR}.F_IF-GF-MM1.vcf.gz

elif [[ $CHR == 'chr_W' ]]
then
        #ONLY FEMALES
        echo "HAPLOID, SETTING HETEROZYGOUS SITES TO MISSING AND ONLY GRAB FEMALES "
        bcftools view --samples-file FFull_Samples_n302.list -Ou snp_full/${CHR}.IF-GF.vcf.gz | \
                bcftools view --min-alleles 2 --max-alleles 2 --min-ac 1 --threads 5 -i 'F_MISSING<0.1' | \
                bedtools subtract -header -a - -b ${mask} | \
                bcftools +setGT -Ou -- -t "b:AD<1e-5" -n "./." | \
            bcftools +setGT -Oz -o snp_full/${CHR}.IF-GF-MM1.vcf.gz -- --target-gt q --new-gt M -i 'GT=="het"'
        bcftools index --threads 5 snp_full/${CHR}.IF-GF-MM1.vcf.gz

elif [[ $CHR == 'chr_Z' ]]
then
        echo "DIPLOID AUTOSOME"
        #all together
        bcftools view --samples-file MFull_Samples_n302.list -Ou snp_full/${CHR}.IF-GF.vcf.gz | \
                bedtools subtract -header -a - -b ${mask} | \
                bcftools view -Oz -o snp_full/${CHR}.M_IF-GF-MM1.vcf.gz --min-alleles 2 --max-alleles 2 --min-ac 1 --threads 5 -i 'F_MISSING<0.1'
        bcftools index --threads 5 snp_full/${CHR}.M_IF-GF-MM1.vcf.gz

        #Females second
        echo "HAPLOID, SETTING HETEROZYGOUS SITES TO MISSING AND ONLY GRAB FEMALES "
        bcftools view --samples-file FFull_Samples_n302.list -Ou snp_full/${CHR}.IF-GF.vcf.gz | \
        bcftools view --min-alleles 2 --max-alleles 2 --min-ac 1 --threads 5 -i 'F_MISSING<0.1' | \
                bedtools subtract -header -a - -b ${mask} | \
                bcftools +setGT -Ou -- -t "b:AD<1e-5" -n "./." | \
                bcftools +setGT -Oz -o snp_full/${CHR}.F_IF-GF-MM1.vcf.gz -- --target-gt q --new-gt M -i 'GT=="het"'
        bcftools index --threads 5 snp_full/${CHR}.F_IF-GF-MM1.vcf.gz

        echo "DIPLOID AUTOSOME"
        #all together
        bcftools view --samples-file Full_Samples_n302.list -Ou snp_full/${CHR}.IF-GF.vcf.gz | \
                bedtools subtract -header -a - -b ${mask} | \
                bcftools view -Oz -o snp_full/${CHR}.IF-GF-MM1.vcf.gz --min-alleles 2 --max-alleles 2 --min-ac 1 --threads 5 -i 'F_MISSING<0.1'
        bcftools index --threads 5 snp_full/${CHR}.IF-GF-MM1.vcf.gz

else
        echo "DIPLOID AUTOSOME"
        #all together
        bcftools view --samples-file Full_Samples_n302.list -Ou snp_full/${CHR}.IF-GF.vcf.gz | \
                bcftools view -Oz -o snp_full/${CHR}.IF-GF-MM1.vcf.gz --min-alleles 2 --max-alleles 2 --min-ac 1 --threads 5 -i 'F_MISSING<0.1'
        bcftools index --threads 5 snp_full/${CHR}.IF-GF-MM1.vcf.gz
fi
