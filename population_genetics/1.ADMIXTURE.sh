#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=10
#SBATCH --time=200:00:00

#for i in {1..10}; do sbatch -J BAD_BOY_SERGIO_${i} ~/merondun/cuculus_host/population_genetics/2.ADMIXTURE.sh ${i}; done
K=$1

cd admixture

#Run Admixture
admixture -j7 --cv=5 ../merged_unrelfull/autos_canorus_LD.bed ${K} > autos_canorus_LD.log${K}.out

#Evaluate Runs
~/modules/evalAdmix/evalAdmix -plink ../merged_unrelfull/autos_canorus_LD -fname autos_canorus_LD.${K}.P -qname autos_canorus_LD.${K}.Q -P 10 -o eval_${K}