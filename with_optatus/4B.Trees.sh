#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00

CHR=$1

python ~/modules/vcf2phylip.py -i both_related/${CHR}_BothN5.IF-GF-MM1-MAF05.vcf.gz -f --output-folder both_related
iqtree --redo -keep-ident -T 10 -s both_related/${CHR}_BothN5.IF-GF-MM1-MAF05.min4.phy --seqtype DNA -m MFP+ASC -alrt 1000 -B 1000
iqtree --redo -keep-ident -T 10 -s both_related/${CHR}_BothN5.IF-GF-MM1-MAF05.min4.phy.varsites.phy --seqtype DNA -m MFP+ASC -alrt 1000 -B 1000