#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00

#to submit: for i in $(cat ~/merondun/cuculus_host/Chromosomes.list); do sbatch -J FILT_${i} ~/merondun/cuculus_host/wgs_preprocessing/8.Unrelated_SNP_Filter.sh ${i}; done 
CHR=$1
GENOME=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa
mask=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/Coverage_Masks/GCA_017976375.1_bCucCan1.pri_genomic.CHR.N75-DoubleCoverage.mask.bed

mkdir snp_unrelfull 

echo "WORKING ON ${CHR}"

# Use bcftools to perform sample-wise filtering on the vcf file.
bcftools view --types snps --force-samples --samples-file Unrelated_2023OCT27.list -Oz ../merged/${CHR}.vcf.gz -o snp_unrelfull/${CHR}.vcf.gz
bcftools index --threads 10  snp_unrelfull/${CHR}.vcf.gz

# Calculate some statistics for the next round of filtering
SAMPN=$(bcftools query -l snp_unrelfull/${CHR}.vcf.gz | wc -l)
AVGDP=$(bcftools query -f '%DP\n' snp_unrelfull/${CHR}.vcf.gz | datamash median 1 | datamash round 1) # Compute the median coverage depth across all sites in the vcf file
DPHI=$(($AVGDP*2)) # Calculate thresholds for filtering based on depth of coverage

# Apply filters
bcftools view --types snps --threads 10 -e "QUAL < 20 || INFO/DP > $DPHI || INFO/DP < $SAMPN || MQ < 30 || RPBZ < -3 || RPBZ > 3" -Oz -o snp_unrelfull/${CHR}.IF.vcf.gz snp_unrelfull/${CHR}.vcf.gz
bcftools index --threads 10 snp_unrelfull/${CHR}.IF.vcf.gz

#Minimum 2X
MINDP=3

bcftools +setGT -Oz -o snp_unrelfull/${CHR}.IF-GF.vcf.gz snp_unrelfull/${CHR}.IF.vcf.gz -- -t q -i "FMT/DP < ${MINDP}" -n "./."
bcftools index --threads 10 snp_unrelfull/${CHR}.IF-GF.vcf.gz

#for W and MT, set allele imbalance violations (heterozygosity) to missing for FEMALE, no W for MALE, diploid Z for MALE
if [[ $CHR == 'chr_MT' || $CHR == 'chr_W' ]]
then
        echo "HAPLOID, SETTING HETEROZYGOUS SITES TO MISSING FOR BOTH SEXES"
        #All samples
        bcftools view --samples-file Unrelated_2023OCT27.list -Ob snp_unrelfull/${CHR}.IF-GF.vcf.gz | \
                bcftools view --types snps --min-alleles 2 --max-alleles 2 --threads 10 | \
                #remove SNPs in bad coverage regions 
                bedtools subtract -header -a - -b ${mask} | \
                #set het genotypes to missing based on binomial test 
                bcftools +setGT -Ob -- -t "b:AD<1e-5" -n "./." | \
                #set weakly het genotypes to major allele 
                bcftools +setGT -Ob -- --target-gt q --new-gt M -i 'GT=="het"' | \
                #set to haploid 
                bcftools +fixploidy -Ob - -- -f 1 | \
                #update AC fields 
                bcftools +fill-tags -Ob -- -t AC,AN | \
                bcftools view --min-alleles 2 --max-alleles 2 --types snps -i 'F_MISSING<0.1' --min-ac 1 --max-af 0.9999 -Oz -o snp_unrelfull/${CHR}.IF-GF-MM1.vcf.gz
        bcftools index --threads 10 snp_unrelfull/${CHR}.IF-GF-MM1.vcf.gz

elif [[ $CHR == 'chr_Z' ]]
then
        echo "DIPLOID AUTOSOME"
        #all together
        bcftools view --samples-file Unrelated_2023OCT27.list -Ov snp_unrelfull/${CHR}.IF-GF.vcf.gz | \
                bedtools subtract -header -a - -b ${mask} | \
                bcftools +fill-tags -Ob -- -t AC,AN | \
                bcftools view -Oz -o snp_unrelfull/${CHR}.IF-GF-MM1.vcf.gz --types snps --min-alleles 2 --max-alleles 2 --min-ac 1 --max-af 0.9999 --threads 10 -i 'F_MISSING<0.1'
        bcftools index --threads 10 snp_unrelfull/${CHR}.IF-GF-MM1.vcf.gz

        #Females second, apply same filtering as chrW/MT
        grep '_F$' Unrelated_2023OCT27.list > Females.tmp 
        echo "HAPLOID, SETTING HETEROZYGOUS SITES TO MISSING AND ONLY GRAB FEMALES "
        bcftools view --samples-file Females.tmp -Ob snp_unrelfull/${CHR}.IF-GF.vcf.gz | \
                bcftools view --types snps --min-alleles 2 --max-alleles 2 --threads 10 | \
                #remove SNPs in bad coverage regions 
                bedtools subtract -header -a - -b ${mask} | \
                #set het genotypes to missing based on binomial test 
                bcftools +setGT -Ob -- -t "b:AD<1e-5" -n "./." | \
                #set weakly het genotypes to major allele 
                bcftools +setGT -Ob -- --target-gt q --new-gt M -i 'GT=="het"' | \
                #set to haploid 
                bcftools +fixploidy -Ob - -- -f 1 | \
                #update AC fields 
                bcftools +fill-tags -Ob -- -t AC,AN | \
                bcftools view --min-alleles 2 --max-alleles 2 --types snps -i 'F_MISSING<0.1' --min-ac 1 --max-af 0.9999 -Oz -o snp_unrelfull/${CHR}.F_IF-GF-MM1.vcf.gz
        bcftools index --threads 10 snp_unrelfull/${CHR}.F_IF-GF-MM1.vcf.gz

        #Males last, diploid 
        echo "DIPLOID AUTOSOME, MALES ONLY "
        grep '_M$' Unrelated_2023OCT27.list > Males.tmp 
        bcftools view --samples-file Males.tmp -Ov snp_unrelfull/${CHR}.IF-GF.vcf.gz | \
                bedtools subtract -header -a - -b ${mask} | \
                bcftools +fill-tags -Ob -- -t AC,AN | \
                bcftools view -Oz -o snp_unrelfull/${CHR}.M_IF-GF-MM1.vcf.gz --types snps --min-alleles 2 --max-alleles 2 --min-ac 1 --max-af 0.9999 --threads 10 -i 'F_MISSING<0.1'
        bcftools index --threads 10 snp_unrelfull/${CHR}.M_IF-GF-MM1.vcf.gz

        rm Males.tmp Females.tmp

else
        echo "DIPLOID AUTOSOME"
        #all together
        bcftools view --samples-file Unrelated_2023OCT27.list -Ob snp_unrelfull/${CHR}.IF-GF.vcf.gz | \
                bcftools +fill-tags -Ob -- -t AC,AN | \
                bcftools view -Oz -o snp_unrelfull/${CHR}.IF-GF-MM1.vcf.gz --types snps --min-alleles 2 --max-alleles 2 --min-ac 1 --max-af 0.9999 --threads 10 -i 'F_MISSING<0.1'
        bcftools index --threads 10 snp_unrelfull/${CHR}.IF-GF-MM1.vcf.gz
fi
