#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_highmem
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00

# mamba activate R
# for CHR in $(cat Chromosomes.list); do sbatch -J dnds_${CHR} ~/merondun/cuculus_host/gene_hunt/3B.Submit_dNdS.sh ${CHR}; done 
CHR=$1

Rscript ~/merondun/cuculus_host/gene_hunt/3.Calculate_dNdS.R ${CHR}