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
fastadir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/cuculiformes_NDUFAF4_blasting/Alignment_Based_Approach/Whole_Genome_Alignment/10gb/fastas
bamdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/cuculiformes_NDUFAF4_blasting/Alignment_Based_Approach/Whole_Genome_Alignment/10gb/bams

cd $fastadir
mkdir results scratch scratch/${RUN}
cd scratch/${RUN}

# Generate consensus sequence in VCF format
bcftools mpileup -f $genome -Ou -a AD -R $regions ${bamdir}/${RUN}.bam | \
	        bcftools call -m -Ov | \
		        bcftools norm -f $genome -Oz -o ${RUN}.vcf.gz
bcftools index ${RUN}.vcf.gz

# If coverage is below 1x, or MQ < 50 - exclude!
MINDP=1
bcftools view -V indels -e "DP < ${MINDP} || MQ < 50 || F_MISSING > 0.1" -Oz -o ${RUN}.Filtered.vcf.gz ${RUN}.vcf.gz
bcftools index ${RUN}.Filtered.vcf.gz

# Create FASTA file for chrW with '-' for regions with no coverage
samtools faidx $genome chr_W:21149910-21177468 | \
	    bcftools consensus ${RUN}.Filtered.vcf.gz --absent - | \
	        awk -v run=${RUN}_W '/^>/{print ">" run; next} {print}' > ${fastadir}/results/${RUN}.W.fa

# and for chr3
samtools faidx $genome chr_3:25511697-25515677 | \
	    bcftools consensus ${RUN}.Filtered.vcf.gz --absent - | \
	        awk -v run=${RUN}_3 '/^>/{print ">" run; next} {print}' > ${fastadir}/results/${RUN}.3.fa
