#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --time=48:00:00

RUN=$1

popdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/manyhost_hunt/males/randomizing_fst
vcfdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/manyhost_hunt/males/female_vcfs
outdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/manyhost_hunt/males/randomizing_fst

popfile=$popdir/${RUN}.popfile
p1="E10"
p2="E6"

mkdir fst fst/work fst/out

# Separate the randomized pops 
grep -w ${p1} ${popfile} > fst/work/${RUN}_${p1}.list
grep -w ${p2} ${popfile} > fst/work/${RUN}_${p2}.list


for CHR in $(cat Chromosomes.list); do

if [[ $CHR = 'chr_W' || $CHR = 'chr_MT' || $CHR = 'chr_Z' ]]; then

        echo "WORKING ON FEMALE HAPLOID"
        vcftools --haploid --gzvcf ${vcfdir}/${CHR}.FEM.DP3.vcf.gz --out fst/work/${CHR}_${RUN} \
                --weir-fst-pop fst/work/${RUN}_${p1}.list --weir-fst-pop fst/work/${RUN}_${p2}.list --fst-window-size 1 --max-missing 0.1

else

        echo "WORKING ON FEMALE DIPLOID"
        vcftools --gzvcf ${vcfdir}/${CHR}.FEM.DP3.vcf.gz --out fst/work/${CHR}_${RUN} \
                --weir-fst-pop fst/work/${RUN}_${p1}.list --weir-fst-pop fst/work/${RUN}_${p2}.list --fst-window-size 1 --max-missing 0.1

fi

        # Format the output and merge it with the gene and dN / dS data 
        awk -v run=${RUN} '{OFS="\t"}{print $1, $2, $2, run, $5}' fst/work/${CHR}_${RUN}.windowed.weir.fst | \
                sed '1d' | bedtools intersect -a - -b Gene_Lookup.bed -wao  | \
                bedtools intersect -a - -b ${vcfdir}/Annotated_Variants_2024APR2.bed -wao | \
                awk '{OFS="\t"}{print $1, $2, $4, $5,$10,$15,$16}' | \
                sed 's/ID=//g' | sed 's/gene-//g' > fst/out/${CHR}_${RUN}.fst.txt

        # Subset fixed snps 
        awk '$4 == 1' fst/out/${CHR}_${RUN}.fst.txt > fst/out/${CHR}_${RUN}.fst1.txt

done

