#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=3
#SBATCH --time=48:00:00

pair=$1
p1=$(echo $pair | sed 's/@.*//g')
p2=$(echo $pair | sed 's/.*@//g')

mkdir fst_bp fst_bp/work fst_bp/out

for CHR in $(cat Chromosomes.list); do

if [[ $CHR = 'chr_W' || $CHR = 'chr_Z' || $CHR = 'chr_MT' ]]; then
        vcftools --gzvcf vcfs/${CHR}.SNP.DP3-AC1-MQ40.vcf.gz --out fst_bp/work/${CHR}_${p1}__${p2} --weir-fst-pop pops/${p1}.list --weir-fst-pop pops/${p2}.list --fst-window-size 1 --haploid

else
        vcftools --gzvcf vcfs/${CHR}.SNP.DP3-AC1-MQ40.vcf.gz --out fst_bp/work/${CHR}_${p1}__${p2} --weir-fst-pop pops/${p1}.list --weir-fst-pop pops/${p2}.list --fst-window-size 1

fi

awk -v p1=${p1} -v p2=${p2} '{OFS="\t"}{print p1, p2, $1, $2, $5}' fst_bp/work/${CHR}_${p1}__${p2}.windowed.weir.fst | sed '1d' > fst_bp/out/${CHR}_${p1}__${p2}.fst.txt

done
