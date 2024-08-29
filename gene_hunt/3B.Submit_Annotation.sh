#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_highmem
#SBATCH --cpus-per-task=10
#SBATCH --time=8:00:00

CHR=$1

#Rscript DNDS.R ${CHR}

sed '1d' raw_dnds/Annotated_Variants_${CHR}__2024APR1.txt | awk '{OFS="\t"}{print $1, $2, $2, $5, $6}' | sed 's/gene-//g' | sed 's/ID=//g' > raw_dnds/${CHR}.bed

