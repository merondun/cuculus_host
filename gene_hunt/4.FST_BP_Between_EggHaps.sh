#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --mem=5000mb
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00

#mamba activate snps
# for CAT in $(cat Pairwise_Contrasts.list); do sbatch -J FST_${CAT} FST.sh ${CAT}; done

CAT=$1

p1=$(echo ${CAT} | sed 's/@.*//g')
p2=$(echo ${CAT} | sed 's/.*@//g')

mkdir fst fst/work fst/out

for CHR in $(cat Chromosomes.list); do

if [[ $CHR = 'chr_W' || $CHR = 'chr_MT' || $CHR = 'chr_Z' ]]; then

        echo "WORKING ON FEMALE HAPLOID"
        vcftools --haploid --gzvcf female_vcfs/${CHR}.FEM.DP3.vcf.gz --out fst/work/${CHR}_${p1}_${p2} --weir-fst-pop pops/${p1}.list --weir-fst-pop pops/${p2}.list --fst-window-size 1 --max-missing 0.1

else

        echo "WORKING ON FEMALE DIPLOID"
        vcftools --gzvcf female_vcfs/${CHR}.FEM.DP3.vcf.gz --out fst/work/${CHR}_${p1}_${p2} --weir-fst-pop pops/${p1}.list --weir-fst-pop pops/${p2}.list --fst-window-size 1 --max-missing 0.1

fi

    awk -v p1=${p1} -v p2=${p2} '{OFS="\t"}{print $1, $2, $2, p1, p2, $5}' fst/work/${CHR}_${p1}_${p2}.windowed.weir.fst | \
            sed '1d' | bedtools intersect -a - -b Gene_Lookup.bed -wao  | \
            bedtools intersect -a - -b female_vcfs/Annotated_Variants_2024APR2.bed -wao | \
        awk '{OFS="\t"}{print $1, $2, $4, $5,$6,$11,$16,$17}' | \
        sed 's/ID=//g' | sed 's/gene-//g' > fst/out/${CHR}_${p1}_${p2}.fst.txt
awk '$5 == 1' fst/out/${CHR}_${p1}_${p2}.fst.txt > fst/out/${CHR}_${p1}_${p2}.fst1.txt

done

