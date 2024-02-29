#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=5
#SBATCH --time=200:00:00

mask=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/Coverage_Masks/GCA_017976375.1_bCucCan1.pri_genomic.CHR.N75-DoubleCoverage.mask.bed
CHR=$1

if [[ $CHR = 'chr_W' || $CHR = 'chr_Z' || $CHR = 'chr_MT' ]]; then

        #Separate invariant sites too
        echo "FILTERING AND MERGING INVARIANT SITES FOR ${CHR}"
        bcftools view --threads 5 -Ou --samples-file AllSamples.list ../../merged/${CHR}.vcf.gz | \
                bcftools view --max-ac 0 -i 'F_MISSING < 0.5 & MQ > 20' --threads 5 -Ou | \
                bcftools +fixploidy -Oz -o vcfs_raw/${CHR}.1N.vcf.gz - -- -f 1
        bcftools index --threads 5 vcfs_raw/${CHR}.1N.vcf.gz

        #re-merge the invariant and filtered SNPs
        bcftools concat --threads 5 -Ob vcfs_raw/${CHR}.1N.vcf.gz vcfs_raw/${CHR}.vcf.gz | \
                bcftools sort -Oz -o vcfs_raw/${CHR}.AllSites.vcf.gz
        bcftools index --threads 5 vcfs_raw/${CHR}.AllSites.vcf.gz

        #create simon's divergence input file from the filtered vcf and the raw vcf
        echo "CREATE SIMON MARTIN DIVERGENCE OUTPUT "
        python ~/modules/genomics_general/VCF_processing/parseVCF.py --ploidyMismatchToMissing --ploidy 1 --skipIndels -i vcfs_raw/${CHR}.AllSites.vcf.gz | \
                bgzip > vcfs_raw/${CHR}.AllSites.geno.gz

else
        #Separate invariant sites too
        echo "FILTERING AND MERGING INVARIANT SITES FOR ${CHR}"
        bcftools view --threads 5 -Ou --samples-file AllSamples.list ../../merged/${CHR}.vcf.gz | \
                bcftools view --max-ac 0 -i 'F_MISSING < 0.5 & MQ > 20' --threads 5 -Oz -o vcfs_raw/${CHR}.1N.vcf.gz 
        bcftools index --threads 5 vcfs_raw/${CHR}.1N.vcf.gz

        #re-merge the invariant and filtered SNPs
        bcftools concat --threads 5 -Ob vcfs_raw/${CHR}.1N.vcf.gz vcfs_raw/${CHR}.vcf.gz | \
                bcftools sort -Oz -o vcfs_raw/${CHR}.AllSites.vcf.gz
        bcftools index --threads 5 vcfs_raw/${CHR}.AllSites.vcf.gz

        #create simon's divergence input file from the filtered vcf and the raw vcf
        echo "CREATE SIMON MARTIN DIVERGENCE OUTPUT "
        python ~/modules/genomics_general/VCF_processing/parseVCF.py --ploidyMismatchToMissing --ploidy 2 --skipIndels -i vcfs_raw/${CHR}.AllSites.vcf.gz | \
                bgzip > vcfs_raw/${CHR}.AllSites.geno.gz

fi