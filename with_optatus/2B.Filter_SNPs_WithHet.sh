#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00

CHR=$1
GENOME=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa

mkdir related

echo "WORKING ON ${CHR}"

if [[ $CHR == 'chr_W' || $CHR == 'chr_Z' ]]; then
    MINDP=5
else
    MINDP=10
fi

# More filtereing, make tree 
bcftools view --threads 10 -e 'F_MISSING > 0.1' --genotype ^het related/${CHR}.IF-GF.vcf.gz -Oz -o related/${CHR}.IF-GF-MM1-NoHet.vcf.gz
bcftools index --threads 10 related/${CHR}.IF-GF-MM1-NoHet.vcf.gz
bcftools index -n related/${CHR}.IF-GF-MM1-NoHet.vcf.gz

python ~/modules/vcf2phylip.py -i related/${CHR}.IF-GF-MM1-NoHet.vcf.gz -f --output-folder related
iqtree --redo -keep-ident -T 10 -s related/${CHR}.IF-GF-MM1-NoHet.min4.phy --seqtype DNA -m MFP+ASC -alrt 1000 -B 1000
iqtree --redo -keep-ident -T 10 -s related/${CHR}.IF-GF-MM1-NoHet.min4.phy.varsites.phy --seqtype DNA -m MFP+ASC -alrt 1000 -B 1000