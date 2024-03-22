#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=10
#SBATCH --time=200:00:00

# mamba activate snps
# for CHR in $(cat Chromosomes.list); do sbatch -J FST_${CHR} ~/merondun/cuculus_host/gene_hunt/1A.Subset_Groups_AllinOne_FST.sh ${CHR}; done 
CHR=$1

#genotypes BELOW this will be set to missing 
MINDP=3

if [[ $CHR = 'chr_W' || $CHR = 'chr_Z' || $CHR = 'chr_MT' ]]; then      
        PLOIDY=1
else
        PLOIDY=2
fi 

mkdir vcfs divergence

echo "FILTERING AND MERGING VARIANT SITES FOR ${CHR}, PLOIDY: ${PLOIDY}"
bcftools view --threads 10 --samples-file ManyHost.list -Ou ../../merged/${CHR}.vcf.gz | \
        bcftools view --types snps --min-alleles 2 --max-alleles 2 --threads 10 | \
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
        bcftools view --min-alleles 2 --max-alleles 2 --types snps --min-ac 1 -Oz -o vcfs/${CHR}.SNP.DP3.vcf.gz 
bcftools index --threads 10 vcfs/${CHR}.SNP.DP3.vcf.gz 

#also filter on DP 
bcftools query -f '%CHROM\t%POS\t%INFO/DP\n' vcfs/${CHR}.SNP.DP3.vcf.gz > ${CHR}_dp_stats.txt
mean=$(cat ${CHR}_dp_stats.txt | datamash mean 3 | xargs printf "%.0f\n")
sd=$(cat ${CHR}_dp_stats.txt | datamash sstdev 3 | xargs printf "%.0f\n")
low=$(echo "$mean - 2*$sd" | bc)
high=$(echo "$mean + 2*$sd" | bc)

rm ${CHR}_dp_stats.txt

#filter, INCLUDE singletons
bcftools view vcfs/${CHR}.SNP.DP3.vcf.gz --min-alleles 2 --max-alleles 2 --min-ac 1 --max-af 0.999 --types snps -i "MQ>40 & INFO/DP > ${low} & INFO/DP < ${high} & F_MISSING < 0.1" -Oz -o vcfs/${CHR}.SNP.DP3-AC1-MQ40.vcf.gz
bcftools index --threads 10 vcfs/${CHR}.SNP.DP3-AC1-MQ40.vcf.gz

#Separate invariant sites too
echo "FILTERING AND MERGING INVARIANT SITES FOR ${CHR}"
bcftools view --threads 10 -Ou --samples-file ManyHost.list ../../merged/${CHR}.vcf.gz | \
        bcftools view --max-ac 0 -i 'F_MISSING < 0.1 & MQ > 40' --threads 10 -Ou | \
        bcftools +fixploidy -Oz -o vcfs/${CHR}.1N.vcf.gz - -- -f ${PLOIDY}
bcftools index --threads 10 vcfs/${CHR}.1N.vcf.gz

#re-merge the invariant and filtered SNPs
bcftools concat --threads 10 -Ob vcfs/${CHR}.1N.vcf.gz vcfs/${CHR}.SNP.DP3-AC1-MQ40.vcf.gz | \
        bcftools sort -Oz -o vcfs/${CHR}.AllSites.vcf.gz
bcftools index --threads 10 vcfs/${CHR}.AllSites.vcf.gz
tabix vcfs/${CHR}.AllSites.vcf.gz

#pixy
TMPDIR=/dss/dsslegfs01/pr53da/pr53da-dss-0021/tmp
pixy --stats fst --vcf vcfs/${CHR}.AllSites.vcf.gz --populations ManyHost.pop --window_size 20000 --n_cores 10 --output_folder divergence --output_prefix ${CHR}.20KB
pixy --stats fst --vcf vcfs/${CHR}.AllSites.vcf.gz --populations ManyHost.pop --window_size 1 --n_cores 10 --output_folder divergence --output_prefix ${CHR}.BP