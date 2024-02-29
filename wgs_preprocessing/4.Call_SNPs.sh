### FIRST! Split the genome into chunks, I do 100 for 100 parallel jobs 

genome=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa

#create dict
gatk CreateSequenceDictionary -R $genome
#split genome by Ns 
gatk ScatterIntervalsByNs -R $genome -O scatter.interval_list -OT ACGT

#split intervals
gatk SplitIntervals -R $genome -L scatter.interval_list --scatter-count 100 -O interval_files/

#Now, call SNPs on each chunk: 


#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=5
#SBATCH --time=240:00:00

vcfdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/2022_11/raw_vcf
genome=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa
ints=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/2022_11/snpcalling/intervals

CHR=$1
GROUP=$2

mkdir $vcfdir
if [[ $CHR == 0102 ]]
then
	maxdp=20000
else
	maxdp=100
fi
echo "WORKING ON INTERVAL ${CHR} WITH MAX DEPTH ${maxdp}"
#bcftools 1.16
bcftools mpileup --max-depth ${maxdp} -C 50 -threads 5 -f ${genome} -b ${GROUP}.bam -R $ints/$CHR.bed -a "AD,ADF,ADR,DP,SP" | \
	bcftools call --ploidy 2 --threads 5 -m -Oz -o ${vcfdir}/${CHR}_bcft.vcf.gz