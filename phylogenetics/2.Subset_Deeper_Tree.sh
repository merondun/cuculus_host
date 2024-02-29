#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=24:00:00

#mamba activate snps
#to submit all iterations: for j in $(cat SUBSETS.list); do sbatch -J TREE_${j} ~/merondun/cuculus_host/phylogenetics/2.Subset_Deeper_Tree.sh chr_MT ${j} ; done 

#mask with male-biased coverage 
mask=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/Coverage_Masks/GCA_017976375.1_bCucCan1.pri_genomic.CHR.N75-DoubleCoverage.mask.bed

CHR=$1
SAMP=$2

#minimum coverage, LESS than this set to missing 
MINDP=3
echo "FILTERING AND MERGING VARIANT SITES FOR ${CHR}"
bcftools view --threads 10 --samples-file ${CHR}_${SAMP}.list -Ou ../../../../merged/${CHR}.vcf.gz | \
        bcftools view --types snps --min-alleles 2 --max-alleles 2 --threads 10 | \
        #remove SNPs in bad coverage regions 
        bedtools subtract -header -a - -b ${mask} | \
        #set genotypes below MINDP to missing 
        bcftools +setGT -- -t q -i "FMT/DP < ${MINDP}" -n "./." | \
        #set het genotypes to missing based on binomial test 
        bcftools +setGT -Ou -- -t "b:AD<1e-5" -n "./." | \
        #set weakly het genotypes to major allele 
        bcftools +setGT -Ou -- --target-gt q --new-gt M -i 'GT=="het"' | \
        #set to haploid, can skip this for most purposes 
        #bcftools +fixploidy -Ou - -- -f 1 | \
        #update AC fields 
        bcftools +fill-tags -Ou -- -t AC,AN | \
        bcftools view --min-alleles 2 --max-alleles 2 --types snps --min-ac 1 -Oz -o vcfs/${CHR}_${SAMP}.SNP.DP3.vcf.gz 
bcftools index --threads 10 vcfs/${CHR}_${SAMP}.SNP.DP3.vcf.gz 

#also filter on DP 
bcftools query -f '%CHROM\t%POS\t%INFO/DP\t%INFO/MQ\n' vcfs/${CHR}_${SAMP}.SNP.DP3.vcf.gz > ${CHR}_${SAMP}_dp_stats.txt
mean=$(cat ${CHR}_${SAMP}_dp_stats.txt | datamash mean 3 | xargs printf "%.0f\n")
sd=$(cat ${CHR}_${SAMP}_dp_stats.txt | datamash sstdev 3 | xargs printf "%.0f\n")
low=$(echo "$mean - 2*$sd" | bc)
high=$(echo "$mean + 2*$sd" | bc)

rm ${CHR}_${SAMP}_dp_stats.txt

#filter 
bcftools view vcfs/${CHR}_${SAMP}.SNP.DP3.vcf.gz --min-alleles 2 --max-alleles 2 --min-ac 1 --max-af 0.999 --types snps -i "MQ>40 & INFO/DP > ${low} & INFO/DP < ${high} & F_MISSING < 0.1" -Oz -o vcfs/${CHR}_${SAMP}.SNP.DP3-AC1-MQ40.vcf.gz
bcftools index --threads 10 vcfs/${CHR}_${SAMP}.SNP.DP3-AC1-MQ40.vcf.gz 

#create tree  
python ~/modules/vcf2phylip.py -i vcfs/${CHR}_${SAMP}.SNP.DP3-AC1-MQ40.vcf.gz -f --output-folder ml_trees
iqtree --redo -keep-ident -T 20 -s ml_trees/${CHR}_${SAMP}.SNP.DP3-AC1-MQ40.min4.phy --seqtype DNA -m MFP+ASC -alrt 1000 -B 1000

#filter with NO SINGLETONS 
bcftools view vcfs/${CHR}_${SAMP}.SNP.DP3.vcf.gz --min-alleles 2 --max-alleles 2 --min-ac 2 --max-af 0.999 --types snps -i "MQ>40 & INFO/DP > ${low} & INFO/DP < ${high} & F_MISSING < 0.1" -Oz -o vcfs/${CHR}_${SAMP}.SNP.DP3-AC2-MQ40.vcf.gz
bcftools index --threads 10 vcfs/${CHR}_${SAMP}.SNP.DP3-AC2-MQ40.vcf.gz 

#create tree  
python ~/modules/vcf2phylip.py -i vcfs/${CHR}_${SAMP}.SNP.DP3-AC2-MQ40.vcf.gz -f --output-folder ml_trees
iqtree --redo -keep-ident -T 20 -s ml_trees/${CHR}_${SAMP}.SNP.DP3-AC2-MQ40.min4.phy --seqtype DNA -m MFP+ASC -alrt 1000 -B 1000