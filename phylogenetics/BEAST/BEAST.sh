#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=10
#SBATCH --time=2-00:00:00

# mamba activate beast

RUN=$1
beast -threads 10 -overwrite -beagle_SSE -seed 777 -java ${RUN}.xml
