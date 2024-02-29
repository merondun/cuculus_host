#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=3
#SBATCH --time=80:00:00

GROUP=$1

P1=$(echo ${GROUP} | sed 's/__.*//g')
P2=$(echo ${GROUP} | sed 's/.*__//g')
popdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/with_optatus/gwas/pops
vcfdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/with_optatus/both_related

for CHR in $(cat Chromosomes.list); do 

mkdir fst fst/out fst/work

echo "WORKING ON CHR: ${CHR}"

#calculate fst 
~/modules/vcftools/bin/vcftools --gzvcf $vcfdir/${CHR}_BothN5.IF-GF-MM1-MAF05.vcf.gz --weir-fst-pop $popdir/${P1}.list --weir-fst-pop $popdir/${P1}.list --out fst/work/${GROUP}_${CHR}
awk -v p1=${P1} -v p2=${P2} '{OFS="\t"}{print $1, $2, $3, p1, p2}' fst/work/${GROUP}_${CHR}.weir.fst | grep -v 'nan' > fst/out/${GROUP}_${CHR}.fst

#done 