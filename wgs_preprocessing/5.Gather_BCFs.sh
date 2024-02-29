#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=5
#SBATCH --time=48:00:00

genome=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa

mkdir ../filtered
mkdir ../merged
RUN=$1
CHR=$(echo ${RUN} | sed 's/__.*//g')
COMPARTMENT=$(echo ${RUN} | sed 's/.*__//g')

echo "${CHR} in compartment ${COMPARTMENT}"
mask=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/${CHR}.consmask
bcftools concat --threads 5 -a -r ${CHR} -Ou ../raw/*${COMPARTMENT}.bcf | \
        bcftools sort -Oz -o ../merged/${CHR}.vcf.gz -
bcftools index --threads 5 ../merged/${CHR}.vcf.gz