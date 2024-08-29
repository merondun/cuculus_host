#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --time=24:00:00

RUN=$1

# Define paths to your files
genome=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa
regions=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/cuculiformes_NDUFAF4_blasting/Alignment_Based_Approach/WholeGene_Regions.bed

mkdir ${RUN}
cd ${RUN}

# Generate consensus sequence in VCF format
bcftools mpileup -f $genome -Ou -a AD -R $regions ../${RUN}.bam | \
		bcftools call -m -Ov | \
			bcftools norm -f $genome -Oz -o ${RUN}.vcf.gz
bcftools index ${RUN}.vcf.gz

# If coverage is below 1x, or MQ < 30 - exclude! 
MINDP=1
bcftools view -e "DP < ${MINDP} || MQ < 30 || F_MISSING > 0.1" -Oz -o ${RUN}.Filtered.vcf.gz ${RUN}.vcf.gz
bcftools index ${RUN}.Filtered.vcf.gz

# Create FASTA file with '-' for regions with no coverage
bcftools consensus -f $genome -o ${RUN}.fa -H 1 --absent - ${RUN}.Filtered.vcf.gz

# Extract region
bedtools getfasta -fi ${RUN}.fa -bed ${regions} -fo ${RUN}_regions.fa -nameOnly
awk -v prefix="${RUN}" '/^>/ {split($0,a," "); sub(">","",a[1]); print ">" prefix "_" a[1]} !/^>/ {print}' ${RUN}_regions.fa > ../gene_alignment/${RUN}.fa
