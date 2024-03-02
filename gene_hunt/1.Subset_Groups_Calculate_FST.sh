#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=5
#SBATCH --time=200:00:00

# mamba activate snps
# for CHR in $(cat Chromosomes.list); do sbatch -J FST_${CHR} ~/merondun/cuculus_host/gene_hunt/1.Subset_Groups_Calculate_FST.sh ${CHR}; done 
CHR=$1
#window size for FST calculations, and maximum number of missing sites 
winsize=500
missize=200
#genotypes BELOW this will be set to missing 
MINDP=3

if [[ $CHR = 'chr_W' || $CHR = 'chr_Z' || $CHR = 'chr_MT' ]]; then      
        PLOIDY=1
else
        PLOIDY=2
fi 

mkdir vcfs fst_bp fst_bp/work fst_bp/out divergence divergence/work divergence/out dnds

echo "FILTERING AND MERGING VARIANT SITES FOR ${CHR}, PLOIDY: ${PLOIDY}"
bcftools view --threads 5 --samples-file GeneHunt_Samples.list -Ou ../../merged/${CHR}.vcf.gz | \
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
        bcftools view --min-alleles 2 --max-alleles 2 --types snps --min-ac 1 -Oz -o vcfs/${CHR}.SNP.DP3.vcf.gz 
bcftools index --threads 5 vcfs/${CHR}.SNP.DP3.vcf.gz 

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
bcftools view --threads 5 -Ou --samples-file GeneHunt_Samples.list ../../merged/${CHR}.vcf.gz | \
        bcftools view --max-ac 0 -i 'F_MISSING < 0.1 & MQ > 40' --threads 5 -Ou | \
        bcftools +fixploidy -Oz -o vcfs/${CHR}.1N.vcf.gz - -- -f ${PLOIDY}
bcftools index --threads 5 vcfs/${CHR}.1N.vcf.gz

#re-merge the invariant and filtered SNPs
bcftools concat --threads 5 -Ob vcfs/${CHR}.1N.vcf.gz vcfs/${CHR}.SNP.DP3-AC1-MQ40.vcf.gz | \
        bcftools sort -Oz -o vcfs/${CHR}.AllSites.vcf.gz
bcftools index --threads 5 vcfs/${CHR}.AllSites.vcf.gz

#create simon's divergence input file from the filtered vcf and the raw vcf
echo "CREATE SIMON MARTIN DIVERGENCE OUTPUT "
python ~/modules/genomics_general/VCF_processing/parseVCF.py --ploidyMismatchToMissing --ploidy ${PLOIDY} --skipIndels -i vcfs/${CHR}.AllSites.vcf.gz | \
        bgzip > vcfs/${CHR}.AllSites.geno.gz

#calculate fst within BOTH blues vs other reversion
for GROUP in $(cat GROUPS.list); do 

        p1=$(echo ${GROUP} | sed 's/__.*//g')
        p2=$(echo ${GROUP} | sed 's/.*__//g')

        echo "CALCULATING FST/DXY/Tajima/PI FOR ${GROUP}"
        popgenWindows.py -w $winsize -m $missize -g vcfs/${CHR}.AllSites.geno.gz \
                -o divergence/work/${CHR}_${GROUP}.csv.gz --analysis popFreq popDist popPairDist \
                -f phased --windType coordinate --ploidy ${PLOIDY} -T 5 -p ${p1} -p ${p2} --popsFile populations/${GROUP}.pop
        zcat divergence/work/${CHR}_${GROUP}.csv.gz | tr ',' '\t' > divergence/out/${CHR}_${GROUP}.txt

        #and also for bp-FST using vcftools 
        if [[ $CHR = 'chr_W' || $CHR = 'chr_Z' || $CHR = 'chr_MT' ]]; then      
                #calculate fst HAPLOID
                ~/modules/vcftools/bin/vcftools --haploid --gzvcf vcfs/${CHR}.SNP.DP3-AC1-MQ40.vcf.gz --weir-fst-pop populations/${p1}.list --weir-fst-pop populations/${p2}.list --out fst_bp/work/${CHR}_${GROUP}
        else
                #calculate fst DIPLOID
                ~/modules/vcftools/bin/vcftools --gzvcf vcfs/${CHR}.SNP.DP3-AC1-MQ40.vcf.gz --weir-fst-pop populations/${p1}.list --weir-fst-pop populations/${p2}.list --out fst_bp/work/${CHR}_${GROUP}
        fi 

        #reformat it nicely
        awk -v h=${p1} -v g=${p2} '{OFS="\t"}{print $1, $2, $3, h, g}' fst_bp/work/${CHR}_${GROUP}.weir.fst | grep -v 'nan' > fst_bp/out/${CHR}_${GROUP}.fst

done 

#and finally, calculate dNdS on the vcf
Rscript ~/merondun/cuculus_host/gene_hunt/3.Calculate_dNdS.R ${CHR}