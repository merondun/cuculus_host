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

# Use bcftools to perform sample-wise filtering on the vcf file.
bcftools view --threads 10 --types snps --force-samples --samples-file AllSamplesHost.list -Oz ../../merged/${CHR}.vcf.gz -o related/${CHR}.vcf.gz
bcftools index --threads 10  related/${CHR}.vcf.gz

# Calculate some statistics for the next round of filtering
SAMPN=$(bcftools query -l related/${CHR}.vcf.gz | wc -l)
AVGDP=$(bcftools query -f '%DP\n' related/${CHR}.vcf.gz | datamash median 1 | datamash round 1) # Compute the median coverage depth across all sites in the vcf file
DPHI=$(($AVGDP*2)) # Calculate thresholds for filtering based on depth of coverage

# Apply filters
bcftools view --types snps --threads 10 -e "QUAL < 20 || INFO/DP > $DPHI || INFO/DP < $SAMPN || MQ < 30 || RPBZ < -3 || RPBZ > 3" -Oz -o related/${CHR}.IF.vcf.gz related/${CHR}.vcf.gz
bcftools index --threads 10 related/${CHR}.IF.vcf.gz

# For diploid chromosomes, convert low-coverage sites to missing
bcftools +setGT -Oz -o related/${CHR}.IF-GF.vcf.gz related/${CHR}.IF.vcf.gz -- -t q -i "FMT/DP < ${MINDP}" -n "./."
bcftools index --threads 10 related/${CHR}.IF-GF.vcf.gz

# More filtereing, make tree 
bcftools view --threads 10 -e 'F_MISSING > 0.1' --genotype ^het related/${CHR}.IF-GF.vcf.gz -Oz -o related/${CHR}.IF-GF-MM1-NoHet.vcf.gz
bcftools index --threads 10 related/${CHR}.IF-GF-MM1-NoHet.vcf.gz
bcftools index -n related/${CHR}.IF-GF-MM1-NoHet.vcf.gz

python ~/modules/vcf2phylip.py -i related/${CHR}.IF-GF-MM1-NoHet.vcf.gz -f --output-folder related
iqtree --redo -keep-ident -T 10 -s related/${CHR}.IF-GF-MM1-NoHet.min4.phy --seqtype DNA -m MFP+ASC -alrt 1000 -B 1000
iqtree --redo -keep-ident -T 10 -s related/${CHR}.IF-GF-MM1-NoHet.min4.phy.varsites.phy --seqtype DNA -m MFP+ASC -alrt 1000 -B 1000