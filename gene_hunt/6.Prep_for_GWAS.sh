#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=5
#SBATCH --time=200:00:00

CHR=$1

#Do first
plink --threads 10 --vcf vcfs/${CHR}.SNP.DP3.vcf.gz --allow-extra-chr --const-fid --set-missing-var-ids @:# --make-bed --out beds/${CHR} --chr-set 39
awk '{print $2, $2, $3, $4, $5, "1"}' beds/${CHR}.fam > beds/${CHR}.tmp
mv beds/${CHR}.tmp beds/${CHR}.fam

