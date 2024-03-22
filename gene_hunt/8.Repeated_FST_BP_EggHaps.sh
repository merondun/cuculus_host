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

mkdir mirror_bp mirror_bp/work mirror_bp/out

for CHR in $(cat Chromosomes.list); do

vcftools --gzvcf vcfs/${CHR}.vcf.gz --out mirror_bp/work/${CHR}_${p1}__${p2} --weir-fst-pop pops/${p1}.list --weir-fst-pop pops/${p2}.list --fst-window-size 1
awk -v p1=${p1} -v p2=${p2} '{OFS="\t"}{print $1, $2, $2, p1, p2, $5}' mirror_bp/work/${CHR}_${p1}__${p2}.windowed.weir.fst | sed '1d' > mirror_bp/out/${CHR}_${p1}__${p2}.fst.txt

done
