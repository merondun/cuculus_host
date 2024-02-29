#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_highmem
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00

# mamba activate R
# sbatch ~/merondun/cuculus_host/correlations_geography_genetics/3.Pairwise_Distance_AMOVA.sh
Rscript ~/merondun/cuculus_host/correlations_geography_genetics/4.AMOVA.R