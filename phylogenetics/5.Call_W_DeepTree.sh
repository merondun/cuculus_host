#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=20
#SBATCH --time=24:00:00

#mamba activate pacbio
#to submit all regions: for i in $(cat REGIONS.list); do sbatch -J PACBIO_${i} ~/merondun/cuculus_host/phylogenetics/5.Call_W_DeepTree.sh ${i} ; done
mkdir vcfs

CHR=$1 

#mask with male-biased coverage 
mask=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/Coverage_Masks/GCA_017976375.1_bCucCan1.pri_genomic.CHR.N75-DoubleCoverage.mask.bed
genome=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa

if [[ $CHR == 'MT' ]]
then
	maxdp=2000000
    ploidy=1
elif [[ $CHR == 'W' ]]
then 
	maxdp=150
    ploidy=1
else
	maxdp=150
    ploidy=2
fi

echo "WORKING ON INTERVAL ${CHR} WITH MAX DEPTH ${maxdp} and PLOIDY: ${ploidy}"
#bcftools 1.16
bcftools mpileup -Ob --min-MQ 30 --skip-indels --max-depth ${maxdp} -C 50 --threads 20 -f ${genome} --bam-list deep_cuculiformes.bams -R region_${CHR}.bed -a "AD,DP" | \
	bcftools call -v --ploidy ${ploidy} --threads 20 -m -Oz -o vcfs/chr_${CHR}_bcft.vcf.gz

#minimum coverage, LESS than this set to missing 
MINDP=3
echo "FILTERING AND MERGING VARIANT SITES FOR ${CHR}"
bcftools view --threads 20 -Ou vcfs/chr_${CHR}_bcft.vcf.gz | \
        bcftools view --types snps --min-alleles 2 --max-alleles 2 --threads 10 | \
        #remove SNPs in bad coverage regions 
        bedtools subtract -header -a - -b ${mask} | \
        #set genotypes below MINDP to missing 
        bcftools +setGT -- -t q -i "FMT/DP < ${MINDP}" -n "./." | \
        #update AC fields 
        bcftools +fill-tags -Ou -- -t AC,AN | \
        bcftools view --min-alleles 2 --max-alleles 2 --types snps --min-ac 1 -Oz -o vcfs/${CHR}.SNP.DP3.vcf.gz 
bcftools index --threads 10 vcfs/${CHR}.SNP.DP3.vcf.gz 

#also filter on DP 
bcftools query -f '%CHROM\t%POS\t%INFO/DP\t%INFO/MQ\n' vcfs/${CHR}.SNP.DP3.vcf.gz > ${CHR}_dp_stats.txt
mean=$(cat ${CHR}_dp_stats.txt | datamash mean 3 | xargs printf "%.0f\n")
sd=$(cat ${CHR}_dp_stats.txt | datamash sstdev 3 | xargs printf "%.0f\n")
low=$(echo "$mean - 2*$sd" | bc)
high=$(echo "$mean + 2*$sd" | bc)

rm ${CHR}_dp_stats.txt

#filter, include singletons 
bcftools view vcfs/${CHR}.SNP.DP3.vcf.gz --min-alleles 2 --max-alleles 2 --min-ac 1 --max-af 0.999 --types snps -i "MQ>40 & INFO/DP > ${low} & INFO/DP < ${high} & F_MISSING < 0.1" -Oz -o vcfs/${CHR}.SNP.DP3-AC1-MQ40.vcf.gz
bcftools index --threads 10 vcfs/${CHR}.SNP.DP3-AC1-MQ40.vcf.gz 

#create tree  
python ~/modules/vcf2phylip.py -i vcfs/${CHR}.SNP.DP3-AC1-MQ40.vcf.gz -f --output-folder ml_trees
iqtree --redo -keep-ident -T 20 -s ml_trees/${CHR}.SNP.DP3-AC1-MQ40.min4.phy --seqtype DNA -m MFP+ASC -alrt 1000 -B 1000

#filter with NO SINGLETONS 
bcftools view vcfs/${CHR}.SNP.DP3.vcf.gz --min-alleles 2 --max-alleles 2 --min-ac 2 --max-af 0.999 --types snps -i "MQ>40 & INFO/DP > ${low} & INFO/DP < ${high} & F_MISSING < 0.1" -Oz -o vcfs/${CHR}.SNP.DP3-AC2-MQ40.vcf.gz
bcftools index --threads 10 vcfs/${CHR}.SNP.DP3-AC2-MQ40.vcf.gz 

#create tree  
python ~/modules/vcf2phylip.py -i vcfs/${CHR}.SNP.DP3-AC2-MQ40.vcf.gz -f --output-folder ml_trees
iqtree --redo -keep-ident -T 20 -s ml_trees/${CHR}.SNP.DP3-AC2-MQ40.min4.phy --seqtype DNA -m MFP+ASC -alrt 1000 -B 1000