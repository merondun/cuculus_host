#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=5
#SBATCH --time=48:00:00

#mamba activate snps
#submit with: for i in $(cat Chromosomes.list); do for j in $(cat Groups.list); do sbatch -J FILT_${i}_${j} ~/merondun/cuculus_host/with_optatus/3.Filter_SNP_Unrelated.sh ${i} ${j}; done ; done

CHR=$1
GENOME=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa
mask=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/Coverage_Masks/GCA_017976375.1_bCucCan1.pri_genomic.CHR.N75-DoubleCoverage.mask.bed

mkdir both_related

echo "WORKING ON ${CHR}"

# Use bcftools to perform sample-wise filtering on the vcf file.
bcftools view --types snps --force-samples --samples-file Female_N5_EggClade_E1-E2-E4-E5-E6-E11.txt -Oz ../../merged/${CHR}.vcf.gz -o both_related/${CHR}_BothN5.vcf.gz
bcftools index --threads 5  both_related/${CHR}_BothN5.vcf.gz

# Calculate some statistics for the next round of filtering
SAMPN=$(bcftools query -l both_related/${CHR}_BothN5.vcf.gz | wc -l)
AVGDP=$(bcftools query -f '%DP\n' both_related/${CHR}_BothN5.vcf.gz | datamash median 1 | datamash round 1) # Compute the median coverage depth across all sites in the vcf file
DPHI=$(($AVGDP*2)) # Calculate thresholds for filtering based on depth of coverage

# Apply filters
bcftools view --types snps --threads 5 -e "QUAL < 20 || INFO/DP > $DPHI || INFO/DP < $SAMPN || MQ < 30 || RPBZ < -3 || RPBZ > 3" -Oz -o both_related/${CHR}_BothN5.IF.vcf.gz both_related/${CHR}_BothN5.vcf.gz
bcftools index --threads 5 both_related/${CHR}_BothN5.IF.vcf.gz

# For diploid chromosomes, convert low-coverage sites to missing
MINDP=5
bcftools +setGT -Oz -o both_related/${CHR}_BothN5.IF-GF.vcf.gz both_related/${CHR}_BothN5.IF.vcf.gz -- -t q -i "FMT/DP < ${MINDP}" -n "./."
bcftools index --threads 5 both_related/${CHR}_BothN5.IF-GF.vcf.gz

#Filter, with MAF 
bcftools view -Ou both_related/${CHR}_BothN5.IF-GF.vcf.gz | \
        bcftools view --threads 5 --min-alleles 2 --max-alleles 2 --types snps \
            -i 'F_MISSING<0.1'  --min-af 0.05 --max-af 0.95 -Oz -o both_related/${CHR}_BothN5.IF-GF-MM1-MAF05.vcf.gz
bcftools index --threads 5 both_related/${CHR}_BothN5.IF-GF-MM1-MAF05.vcf.gz
