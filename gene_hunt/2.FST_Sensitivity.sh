#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00

# mamba activate snps
# for CHR in $(cat Chromosomes.list); do sbatch -J FST_${CHR} ~/merondun/cuculus_host/gene_hunt/2.FST_Sensitivity.sh ${CHR}; done 
CHR=$1

#genotypes BELOW this will be set to missing 
for WIN in 1000 5000 25000 50000 250000; do

#pixy
TMPDIR=/dss/dsslegfs01/pr53da/pr53da-dss-0021/tmp
pixy --stats fst --vcf vcfs/${CHR}.AllSites.vcf.gz --populations ManyHost.pop --window_size ${WIN} --n_cores 10 --output_folder divergence --output_prefix ${CHR}_${WIN}

done 
