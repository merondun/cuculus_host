#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=2
#SBATCH --time=48:00:00

# mamba activate snps
# for i in $(cat PairwiseComparisons.list); do sbatch -J FST_${i} ~/merondun/cuculus_host/correlations_geography_genetics/2.Pairwise_Distance_FST.sh ${i}; done 
GROUP=$1

for CHR in $(cat Chromosomes.list); do

echo "WORKING ON CHR: ${GROUP} and ${CHR}"

p1=$(echo ${GROUP} | sed 's/__.*//g')
p2=$(echo ${GROUP} | sed 's/.*__//g')

mkdir work out

if [[ $CHR == 'chr_Z' ]]
then

    calculate fst
    ~/modules/vcftools/bin/vcftools --gzvcf ../../snp_unrelfull/${CHR}.IF-GF-MM1.vcf.gz --fst-window-size 100000 --fst-window-step 100000 --weir-fst-pop populations/${p1}.pop --weir-fst-pop populations/${p2}.pop --out work/${CHR}_${GROUP}
    awk -v g=${GROUP} '{OFS="\t"}{print $1, $2, $3,$4,$5, g}' work/${CHR}_${GROUP}.windowed.weir.fst > out/${CHR}_${GROUP}.fst

    #calculate also for only males 
    ~/modules/vcftools/bin/vcftools --gzvcf ../../snp_unrelfull/${CHR}.M_IF-GF-MM1.vcf.gz --fst-window-size 100000 --fst-window-step 100000 --weir-fst-pop populations/${p1}.pop --weir-fst-pop populations/${p2}.pop --out work/${CHR}_${GROUP}.M
    awk -v g=${GROUP} '{OFS="\t"}{print $1, $2, $3,$4,$5, g"__M"}' work/${CHR}_${GROUP}.M.windowed.weir.fst > out/${CHR}_${GROUP}.M.fst

else

    calculate fst
    ~/modules/vcftools/bin/vcftools --gzvcf ../../snp_unrelfull/${CHR}.IF-GF-MM1.vcf.gz --fst-window-size 100000 --fst-window-step 100000 --weir-fst-pop populations/${p1}.pop --weir-fst-pop populations/${p2}.pop --out work/${CHR}_${GROUP}
    awk -v g=${GROUP} '{OFS="\t"}{print $1, $2, $3,$4,$5, g}' work/${CHR}_${GROUP}.windowed.weir.fst > out/${CHR}_${GROUP}.fst

fi

done