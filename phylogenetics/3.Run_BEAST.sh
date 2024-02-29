#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=20
#SBATCH --time=220:00:00

#default base env, or snps
#sbatch -J chr_W_SNAPP ~/merondun/cuculus_host/phylogenetics/3.Run_BEAST.sh chr_W
RUN=$1
~/modules/beast/bin/beast -threads 20 -overwrite -beagle_SSE -seed 777 -java ${RUN}.xml