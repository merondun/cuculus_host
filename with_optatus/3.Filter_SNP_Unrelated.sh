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
GROUP=$2
GENOME=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa
mask=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/Coverage_Masks/GCA_017976375.1_bCucCan1.pri_genomic.CHR.N75-DoubleCoverage.mask.bed

mkdir unrelated

echo "WORKING ON ${CHR}"

# Use bcftools to perform sample-wise filtering on the vcf file.
bcftools view --types snps --force-samples --samples-file Unrel_${GROUP}.list -Oz ../../merged/${CHR}.vcf.gz -o unrelated/${CHR}_${GROUP}.vcf.gz
bcftools index --threads 5  unrelated/${CHR}_${GROUP}.vcf.gz

# Calculate some statistics for the next round of filtering
SAMPN=$(bcftools query -l unrelated/${CHR}_${GROUP}.vcf.gz | wc -l)
AVGDP=$(bcftools query -f '%DP\n' unrelated/${CHR}_${GROUP}.vcf.gz | datamash median 1 | datamash round 1) # Compute the median coverage depth across all sites in the vcf file
DPHI=$(($AVGDP*2)) # Calculate thresholds for filtering based on depth of coverage

# Apply filters
bcftools view --types snps --threads 5 -e "QUAL < 20 || INFO/DP > $DPHI || INFO/DP < $SAMPN || MQ < 30 || RPBZ < -3 || RPBZ > 3" -Oz -o unrelated/${CHR}_${GROUP}.IF.vcf.gz unrelated/${CHR}_${GROUP}.vcf.gz
bcftools index --threads 5 unrelated/${CHR}_${GROUP}.IF.vcf.gz

# For diploid chromosomes, convert low-coverage sites to missing
MINDP=5
bcftools +setGT -Oz -o unrelated/${CHR}_${GROUP}.IF-GF.vcf.gz unrelated/${CHR}_${GROUP}.IF.vcf.gz -- -t q -i "FMT/DP < ${MINDP}" -n "./."
bcftools index --threads 5 unrelated/${CHR}_${GROUP}.IF-GF.vcf.gz

#for W and MT, set allele imbalance violations (heterozygosity) to missing for FEMALE, no W for MALE, diploid Z for MALE
if [[ $CHR == 'chr_MT' || $CHR == 'chr_W' ]]
then
        echo "HAPLOID, SETTING HETEROZYGOUS SITES TO MISSING FOR BOTH SEXES"
        #All samples
        bcftools view --force-samples --samples-file Unrel_${GROUP}.list -Ou unrelated/${CHR}_${GROUP}.IF-GF.vcf.gz | \
                bcftools view --min-alleles 2 --max-alleles 2 --min-ac 1 --threads 5 -i 'F_MISSING<0.1' | \
                bedtools subtract -header -a - -b ${mask} | \
                bcftools +setGT -Ou -- -t "b:AD<1e-5" -n "./." | \
                bcftools +setGT  -- --target-gt q --new-gt M -i 'GT=="het"' | \
                bcftools +fixploidy -Ou -- -f 1 | \
                bcftools view --types snps -i 'F_MISSING<0.1' --threads 5 --min-ac 1 --max-af 0.9999 -Oz -o unrelated/${CHR}_${GROUP}.IF-GF-MM1.vcf.gz
        bcftools index --threads 5 unrelated/${CHR}_${GROUP}.IF-GF-MM1.vcf.gz

        #Separate invariant sites too
        echo "FILTERING AND MERGING INVARIANT SITES FOR ${CHR}"
        bcftools view --threads 5 -Ou --force-samples --samples-file Unrel_${GROUP}.list ../../merged/${CHR}.vcf.gz | \
                bcftools view --max-ac 0 -i 'F_MISSING<0.1 & MQ > 30' --threads 5 -Ou | \
                bcftools +fixploidy -Oz -o unrelated/${CHR}_${GROUP}.1N.vcf.gz - -- -f 1
        bcftools index --threads 5 unrelated/${CHR}_${GROUP}.1N.vcf.gz

        #re-merge the invariant and filtered SNPs
        bcftools concat --threads 5 -Ob unrelated/${CHR}_${GROUP}.1N.vcf.gz unrelated/${CHR}_${GROUP}.IF-GF-MM1.vcf.gz | \
                bcftools sort -Oz -o unrelated/${CHR}_${GROUP}.AllSites.vcf.gz
        bcftools index --threads 5 unrelated/${CHR}_${GROUP}.AllSites.vcf.gz

        #create simon's divergence input file from the filtered vcf and the raw vcf
        echo "CREATE SIMON MARTIN DIVERGENCE OUTPUT "
        python ~/modules/genomics_general/VCF_processing/parseVCF.py --ploidyMismatchToMissing --ploidy 1 --skipIndels -i unrelated/${CHR}_${GROUP}.AllSites.vcf.gz | \
                bgzip > unrelated/${CHR}_${GROUP}.AllSites.geno.gz

elif [[ $CHR == 'chr_Z' ]]
then
        echo "DIPLOID AUTOSOME"
        #males first 
        grep '_M$' Unrel_${GROUP}.list > ${GROUP}_M.tmp
        bcftools view --force-samples --samples-file ${GROUP}_M.tmp -Ov unrelated/${CHR}_${GROUP}.IF-GF.vcf.gz | \
                bedtools subtract -header -a - -b ${mask} | \
                bcftools view -Oz -o unrelated/${CHR}_${GROUP}.M_IF-GF-MM1.vcf.gz --min-alleles 2 --max-alleles 2 --min-ac 1 --threads 5 -i 'F_MISSING<0.1'
        bcftools index --threads 5 unrelated/${CHR}_${GROUP}.M_IF-GF-MM1.vcf.gz

        #Females second
        grep '_F$' Unrel_${GROUP}.list > ${GROUP}_F.tmp
        echo "HAPLOID, SETTING HETEROZYGOUS SITES TO MISSING AND ONLY GRAB FEMALES "
        bcftools view --force-samples --samples-file ${GROUP}_F.tmp -Ou unrelated/${CHR}_${GROUP}.IF-GF.vcf.gz | \
        bcftools view --min-alleles 2 --max-alleles 2 --min-ac 1 --threads 5 -i 'F_MISSING<0.1' | \
                bedtools subtract -header -a - -b ${mask} | \
                bcftools +setGT -Ou -- -t "b:AD<1e-5" -n "./." | \
                bcftools +setGT  -- --target-gt q --new-gt M -i 'GT=="het"' | \
                bcftools +fixploidy -Ou -- -f 1 | \
                bcftools view --types snps -i 'F_MISSING<0.1' --threads 5 --min-ac 1 --max-af 0.9999 -Oz -o unrelated/${CHR}_${GROUP}.IF-GF-MM1.vcf.gz
        bcftools index --threads 5 unrelated/${CHR}_${GROUP}.F_IF-GF-MM1.vcf.gz

        echo "DIPLOID AUTOSOME"
        #all together
        bcftools view --force-samples --samples-file Unrel_${GROUP}.list -Ov unrelated/${CHR}_${GROUP}.IF-GF.vcf.gz | \
                bedtools subtract -header -a - -b ${mask} | \
                bcftools view -Oz -o unrelated/${CHR}_${GROUP}.IF-GF-MM1.vcf.gz --min-alleles 2 --max-alleles 2 --min-ac 1 --threads 5 -i 'F_MISSING<0.1'
        bcftools index --threads 5 unrelated/${CHR}_${GROUP}.IF-GF-MM1.vcf.gz

        #Separate invariant sites too
        echo "FILTERING AND MERGING INVARIANT SITES FOR ${CHR}"
        bcftools view --threads 5 -Ou --force-samples --samples-file ${GROUP}_F.tmp ../../merged/${CHR}.vcf.gz | \
                bcftools view --max-ac 0 -i 'F_MISSING<0.1 & MQ > 30' --threads 5 -Ou | \
                bcftools +fixploidy -Oz -o unrelated/${CHR}_${GROUP}.F_1N.vcf.gz - -- -f 1
        bcftools index --threads 5 unrelated/${CHR}_${GROUP}.F_1N.vcf.gz

        #re-merge the invariant and filtered SNPs
        bcftools concat --threads 5 -Ob unrelated/${CHR}_${GROUP}.F_1N.vcf.gz unrelated/${CHR}_${GROUP}.F_IF-GF-MM1.vcf.gz | \
                bcftools sort -Oz -o unrelated/${CHR}_${GROUP}.F_AllSites.vcf.gz
        bcftools index --threads 5 unrelated/${CHR}_${GROUP}.F_AllSites.vcf.gz

        #create simon's divergence input file from the filtered vcf and the raw vcf
        echo "CREATE SIMON MARTIN DIVERGENCE OUTPUT "
        python ~/modules/genomics_general/VCF_processing/parseVCF.py --ploidyMismatchToMissing --ploidy 1 --skipIndels -i unrelated/${CHR}_${GROUP}.F_AllSites.vcf.gz | \
                bgzip > unrelated/${CHR}_${GROUP}.F_AllSites.geno.gz

        rm ${GROUP}_M.tmp ${GROUP}_F.tmp

else
        echo "DIPLOID AUTOSOME"
        #all together
        bcftools view --force-samples --samples-file Unrel_${GROUP}.list -Ou unrelated/${CHR}_${GROUP}.IF-GF.vcf.gz | \
                bcftools view -Oz -o unrelated/${CHR}_${GROUP}.IF-GF-MM1.vcf.gz --min-alleles 2 --max-alleles 2 --min-ac 1 --threads 5 -i 'F_MISSING<0.1'
        bcftools index --threads 5 unrelated/${CHR}_${GROUP}.IF-GF-MM1.vcf.gz

        #Separate invariant sites too
        echo "FILTERING AND MERGING INVARIANT SITES FOR ${CHR}"
        bcftools view --threads 5 -Ou --force-samples --samples-file Unrel_${GROUP}.list ../../merged/${CHR}.vcf.gz | \
                bcftools view --max-ac 0 -i 'F_MISSING<0.1 & MQ > 30' --threads 5 -Oz -o unrelated/${CHR}_${GROUP}.1N.vcf.gz
        bcftools index --threads 5 unrelated/${CHR}_${GROUP}.1N.vcf.gz

        #re-merge the invariant and filtered SNPs
        bcftools concat --threads 5 -Ob unrelated/${CHR}_${GROUP}.1N.vcf.gz unrelated/${CHR}_${GROUP}.IF-GF-MM1.vcf.gz | \
                bcftools sort -Oz -o unrelated/${CHR}_${GROUP}.AllSites.vcf.gz
        bcftools index --threads 5 unrelated/${CHR}_${GROUP}.AllSites.vcf.gz

        #create simon's divergence input file from the filtered vcf and the raw vcf
        echo "CREATE SIMON MARTIN DIVERGENCE OUTPUT "
        python ~/modules/genomics_general/VCF_processing/parseVCF.py --ploidyMismatchToMissing --ploidy 2 --skipIndels -i unrelated/${CHR}_${GROUP}.AllSites.vcf.gz | \
                bgzip > unrelated/${CHR}_${GROUP}.AllSites.geno.gz

fi

