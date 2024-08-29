#  Cuculus Evolutionary Genetics

Here we examine population genetic variation across Eurasian *Cuculus*. Herein lies our primary target, the diversification of host-specificity and egg morphology.

```bash
title: "Cuculus Evolutionary Genetics"
author: "Justin Merondun"
date: "`r Sys.Date()`"
output:
  html_document:
    theme: journal
    toc: true
    toc_float:
      collapsed: true
```



# Processing

## Trimming

Business as usual; have 658 individual SRRs to trim and align. BBDUK + BWA. 

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=12:00:00
RUN=$1
wd=/dss/dsslegfs01/pr53da/pr53da-dss-0021/rawdata/WGS/ILLUMINA/
outdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/trimmed_fastq

RUN=$1

file1=$(ls ${wd}/${RUN}__R1__*)
file2=$(ls ${wd}/${RUN}__R2__*)

SCRATCH=/tmp/$SLURM_JOB_ID

#adapter trim
bbduk.sh t=6 -Xmx24g in=${file1} in2=${file2} out=${outdir}/${RUN}.trim.fastq.gz minlen=25 qtrim=rl trimq=2 ktrim=r k=23 mink=11 ref=~/modules/bbmap/adapters.fa hdist=1 tpe tbo
```

and submit:

```bash
for i in $(cat Libraries.list); do sbatch -J $i Trim.sh $i; done 
```

## Identify Read Group

This script will find read groups from the fastq files, provided you have normal fastq, and your files are named in our convention (SM__SRR.something.fastq.gz).

```bash
#!/bin/bash

#In a directory with the cleaned .fq.gz files, this will create ID (for us just SRR ID), SM (most important, for VCF files), LB (in case same library across multiple lanes), and PU (run data) fields. It will also count the number of SRRs per SM for merging later. 

mkdir RGs

for i in $( ls *.trim.fastq.gz | sed 's/\..*//g' ); do
        zcat $i.trim.fastq.gz | head -n 1 > ./RGs/$i.txt
        awk 'BEGIN {FS=":"}{print FILENAME, $3, $4, $NF}' ./RGs/$i.txt > ./RGs/$i.tmp
        ex -sc '%s/.txt//g' -c x ./RGs/$i.tmp
        ex -sc '%s:./RGs/::g' -c x ./RGs/$i.tmp

        #okay, now subet the parts we want
        sed 's/__.*//g' ./RGs/$i.tmp > ./RGs/$i.SM
        sed 's/.*__//g' ./RGs/$i.tmp | cut -d ' ' -f 1 > ./RGs/$i.ID
        awk '{print $4}' ./RGs/$i.tmp > ./RGs/$i.tmp2
        paste ./RGs/$i.SM ./RGs/$i.tmp2 | sed 's/\t/./g' > ./RGs/$i.LB

        #run details for PU
        zcat $i.trim.fastq.gz |head -n 1 |cut -f 3-5 -d ":" > ./RGs/$i.PU

        #Determine number of runs per sample, and output to .numlibs file
        ID="$(cat ./RGs/${i}.ID | tail -1)";
        SM="$(cat ./RGs/${i}.SM | tail -1)";
        echo ${ID} >> ./RGs/${SM}.libco
        cat ./RGs/${SM}.libco | sort -u > ./RGs/${SM}.numlibs

done

rm ./RGs/*tmp*
rm ./RGs/*txt*
```



## Alignment

Still on the 658 individual SRR accessions. Submit positional library.

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=12:00:00

SCRATCH=tmp/$SLURM_JOB_ID

wd=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/trimmed_fastq
qcdata=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/trimmed_fastq
outdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/bams/2021_04
genome=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa

RUN=$1

bwa index ${genome}
#Mapping script, will map reads with the prefix identified in the list ${lists}/${RUN}

        ID="$(cat ${wd}/RGs/${RUN}.ID | tail -1)";
        PU="$(cat ${wd}/RGs/${RUN}.PU | tail -1)";
        SM="$(cat ${wd}/RGs/${RUN}.SM | tail -1)";
        LB="$(cat ${wd}/RGs/${RUN}.LB | tail -1)";

        #library mapping
        bwa mem -M -p -t 10 -R "@RG\tID:${ID}\tSM:${SM}\tPL:ILLUMINA\tPU:${PU}\tLB:${LB}" ${genome} ${qcdata}/${RUN}.trim.fastq.gz | samtools sort -@10 -o ${SCRATCH}/${RUN}.bam -;
        samtools view -b -f 1 -F 524 ${SCRATCH}/${RUN}.bam > ${outdir}/${RUN}.bam
        samtools index -b ${outdir}/${RUN}.bam;
```

And submit by sample:

```bash
for sample in $(cat Libraries.list); do sbatch -J ${sample} Align.sh ${sample}; done Plot Alignment:
```


## Merge by Sample

This script will merge, mark duplicates, remove reads overhanging scaffold ends, and summarize (coverage + alignments). Submit positional SAMPLE (n=390). 

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=12:00:00

SCRATCH=tmp/$SLURM_JOB_ID

#Bam cleaning script. Merge by sample, Mark duplicates, remove reads overhanging scaffolds, and summarize

RUN=$1

RGdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/trimmed_fastq/RGs
bamdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/bams/2021_04
genome=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa

mkdir tmp
mkdir $SCRATCH

#Make final bam and stats folder
mkdir ${bamdir}/merged
mkdir ${bamdir}/merged/stats
mkdir ${bamdir}/merged/stats/alignments
mkdir ${bamdir}/merged/stats/coverage

#merge all the libraries together
samtools merge - ${bamdir}/${RUN}*.bam | samtools sort -@3 -o $SCRATCH/${RUN}.bam
samtools index $SCRATCH/${RUN}.bam

sambamba markdup --tmpdir $SCRATCH/ $SCRATCH/${RUN}.bam $SCRATCH/${RUN}.dup.bam
samtools index $SCRATCH/${RUN}.dup.bam

#Remove reads over scaffold end
gatk CleanSam -I $SCRATCH/${RUN}.dup.bam -O ${bamdir}/merged/${RUN}.bam -R ${genome} --CREATE_INDEX true

Summary Stats
gatk CollectAlignmentSummaryMetrics -R ${genome} -I ${bamdir}/merged/${RUN}.bam -O ${bamdir}/merged/stats/alignments/${RUN}.alignments.txt -LEVEL SAMPLE -LEVEL READ_GROUP

#calculate coverage
mosdepth --threads 3 --use-median --by 100000 --fast-mode --no-per-base ${bamdir}/merged/stats/coverage/${RUN} ${bamdir}/merged/${RUN}.bam

```

**Sample Stats**

**Alignments**

Check out alignment stats within /stats/alignments/

```bash
mkdir output
grep 'TOTAL' 001_CB_ORW_CHN_F.alignments.txt > output/Alignments.txt

for i in $(ls *.alignments.txt | sed 's/\..*//g'); do 
awk 'NR==13' ${i}.alignments.txt >> output/Alignments.txt
done
```

**Coverage**

Within /coverage/

```bash
#take the unzipped bed.gz file from mosdepth, and change it into the below
for i in $(ls *.regions.bed | rev | cut -c13- | rev); do awk -v s=${i} 'BEGIN{print s}; {print $4}' ${i}.regions.bed > ${i}.cov; done

#create header file with chr / start / end, since all files on the same genome, they all have the same length and this is consistent across samples 
awk '{OFS="\t"}BEGIN{print "CHR","START","END"};{print $1, $2, $3}' 053_CC_RED_FIN_F.regions.bed | sed 's/\ /\t/g' > 0001Header.cov

#combine them all together..
paste *.cov > Master.coverage
```

Also print the GW mean of median:

```bash
for i in $(ls *.regions.bed | rev | cut -c13- | rev); do awk -v s=${i} '{ total += $4 } END { print s, total/NR }' ${i}.regions.bed >> gw_median.coverage ; done
```

And by chromosome:

```bash
for i in $(ls *.summary.txt | sed 's/\..*//g'); do
grep -i -w 'CM030676.1\|CM030677.1\|total\|CM030679.1' ${i}.mosdepth.summary.txt | awk -v s=${i} 'BEGIN{OFS="\t"}{print s, $1, $4}' >> summarize/Compartment.coverage.txt ; done
```

# Variant Calling

Split Genome

```bash
genome=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa

#create dict
gatk CreateSequenceDictionary -R $genome
#split genome by Ns 
gatk ScatterIntervalsByNs -R $genome -O scatter.interval_list -OT ACGT

#split intervals
gatk SplitIntervals -R $genome -L scatter.interval_list --scatter-count 100 -O interval_files/
```

post-hoc, split only scaffold regions into 50 for fast processing:

```bash
genome=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.SCAFFOLDS.fa

#create dict
gatk CreateSequenceDictionary -R $genome
#split genome by Ns 
gatk ScatterIntervalsByNs -R $genome -O scatter.interval_list -OT ACGT

#split intervals
gatk SplitIntervals -R $genome -L scatter.interval_list --scatter-count 50 -O .
```



## SNP calling

```bash
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
bcftools mpileup --max-depth ${maxdp} -C 50 -threads 5 -f ${genome} -b ${GROUP}.bam -R $ints/$CHR.bed -a "AD,DP,GQ" | \
	bcftools call --ploidy 2 --threads 5 -m -Oz -o ${vcfdir}/${CHR}_bcft.vcf.gz
```

Scaffolds:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=5
#SBATCH --time=240:00:00

genome=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa

REGION=$1
maxdp=150

#bcftools 1.16
echo "WORKING ON INTERVAL ${ints} WITH MAX DEPTH ${maxdp} AND SAMPLES FILE: All.bams"
bcftools mpileup -Ou -I --max-depth ${maxdp} -C 50 --threads 5 -f ${genome} -b All.bams -R ${REGION}.bed -a "AD,ADF,ADR,DP,SP" | \
                bcftools call -a GQ --ploidy 2 --threads 5 -m -Ob -o ${REGION}.bcf
bcftools index ${REGION}.bcf
```

merge scaffolds: 

```bash
bcftools concat --threads 10 -a -Ou *.bcf |  bcftools sort -Oz -o ../../merged/scaffolds.vcf.gz -
bcftools index --threads 10 ../../merged/scaffolds.vcf.gz
```

Merge by chromosome

```bash
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
```

Also grab all 'scaffold' chromosomes and merge:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_highmem
#SBATCH --cpus-per-task=10
#SBATCH --time=80:00:00

bcftools concat --threads 10 -a -r scaffolds.list -Ou ../raw/*.bcf |  bcftools sort -Oz -o ../merged/scaffolds.vcf.gz -
bcftools index --threads 10 ../merged/scaffolds.vcf.gz
```

## SNP Filtering; N = 299

First on the full dataset.

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=5
#SBATCH --time=48:00:00

CHR=$1
GENOME=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa
mask=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/Coverage_Masks/GCA_017976375.1_bCucCan1.pri_genomic.CHR.N75-DoubleCoverage.mask.bed

echo "WORKING ON ${CHR}"

# Use bcftools to perform sample-wise filtering on the vcf file.
bcftools view --types snps --force-samples --samples-file Full_Samples_n302.list -Oz ../merged/${CHR}.vcf.gz -o snp_full/${CHR}.vcf.gz
bcftools index --threads 5  snp_full/${CHR}.vcf.gz

# Calculate some statistics for the next round of filtering
SAMPN=$(bcftools query -l snp_full/${CHR}.vcf.gz | wc -l)
AVGDP=$(bcftools query -f '%DP\n' snp_full/${CHR}.vcf.gz | datamash median 1 | datamash round 1) # Compute the median coverage depth across all sites in the vcf file
DPHI=$(($AVGDP*2)) # Calculate thresholds for filtering based on depth of coverage

# Apply filters
bcftools view --types snps --threads 5 -e "QUAL < 20 || INFO/DP > $DPHI || INFO/DP < $SAMPN || MQ < 30 || RPBZ < -3 || RPBZ > 3" -Oz -o snp_full/${CHR}.IF.vcf.gz snp_full/${CHR}.vcf.gz
bcftools index --threads 5 snp_full/${CHR}.IF.vcf.gz

# For diploid chromosomes, convert low-coverage sites to missing
MINDP=3
bcftools +setGT -Oz -o snp_full/${CHR}.IF-GF.vcf.gz snp_full/${CHR}.IF.vcf.gz -- -t q -i "FMT/DP < ${MINDP}" -n "./."
bcftools index --threads 5 snp_full/${CHR}.IF-GF.vcf.gz

#for W and MT, set allele imbalance violations (heterozygosity) to missing for FEMALE, no W for MALE, diploid Z for MALE
if [[ $CHR == 'chr_MT' ]]
then
        echo "HAPLOID, SETTING HETEROZYGOUS SITES TO MISSING FOR BOTH SEXES"
        #All samples
        bcftools view --samples-file Full_Samples_n302.list -Ou snp_full/${CHR}.IF-GF.vcf.gz | \
                bcftools view --min-alleles 2 --max-alleles 2 --min-ac 1 --threads 5 -i 'F_MISSING<0.1' | \
                bedtools subtract -header -a - -b ${mask} | \
                bcftools +setGT -Ou -- -t "b:AD<1e-5" -n "./." | \
            bcftools +setGT -Oz -o snp_full/${CHR}.IF-GF-MM1.vcf.gz -- --target-gt q --new-gt M -i 'GT=="het"'
        bcftools index --threads 5 snp_full/${CHR}.IF-GF-MM1.vcf.gz

        echo "HAPLOID, SETTING HETEROZYGOUS SITES TO MISSING FOR BOTH SEXES"
        #Females only
        bcftools view --samples-file FFull_Samples_n302.list -Ou snp_full/${CHR}.IF-GF.vcf.gz | \
                bcftools view --min-alleles 2 --max-alleles 2 --min-ac 1 --threads 5 -i 'F_MISSING<0.1' | \
                bedtools subtract -header -a - -b ${mask} | \
                bcftools +setGT -Ou -- -t "b:AD<1e-5" -n "./." | \
            bcftools +setGT -Oz -o snp_full/${CHR}.F_IF-GF-MM1.vcf.gz -- --target-gt q --new-gt M -i 'GT=="het"'
        bcftools index --threads 5 snp_full/${CHR}.F_IF-GF-MM1.vcf.gz

elif [[ $CHR == 'chr_W' ]]
then
        #ONLY FEMALES
        echo "HAPLOID, SETTING HETEROZYGOUS SITES TO MISSING AND ONLY GRAB FEMALES "
        bcftools view --samples-file FFull_Samples_n302.list -Ou snp_full/${CHR}.IF-GF.vcf.gz | \
                bcftools view --min-alleles 2 --max-alleles 2 --min-ac 1 --threads 5 -i 'F_MISSING<0.1' | \
                bedtools subtract -header -a - -b ${mask} | \
                bcftools +setGT -Ou -- -t "b:AD<1e-5" -n "./." | \
            bcftools +setGT -Oz -o snp_full/${CHR}.IF-GF-MM1.vcf.gz -- --target-gt q --new-gt M -i 'GT=="het"'
        bcftools index --threads 5 snp_full/${CHR}.IF-GF-MM1.vcf.gz

elif [[ $CHR == 'chr_Z' ]]
then
        echo "DIPLOID AUTOSOME"
        #all together
        bcftools view --samples-file MFull_Samples_n302.list -Ou snp_full/${CHR}.IF-GF.vcf.gz | \
                bedtools subtract -header -a - -b ${mask} | \
                bcftools view -Oz -o snp_full/${CHR}.M_IF-GF-MM1.vcf.gz --min-alleles 2 --max-alleles 2 --min-ac 1 --threads 5 -i 'F_MISSING<0.1'
        bcftools index --threads 5 snp_full/${CHR}.M_IF-GF-MM1.vcf.gz

        #Females second
        echo "HAPLOID, SETTING HETEROZYGOUS SITES TO MISSING AND ONLY GRAB FEMALES "
        bcftools view --samples-file FFull_Samples_n302.list -Ou snp_full/${CHR}.IF-GF.vcf.gz | \
        bcftools view --min-alleles 2 --max-alleles 2 --min-ac 1 --threads 5 -i 'F_MISSING<0.1' | \
                bedtools subtract -header -a - -b ${mask} | \
                bcftools +setGT -Ou -- -t "b:AD<1e-5" -n "./." | \
                bcftools +setGT -Oz -o snp_full/${CHR}.F_IF-GF-MM1.vcf.gz -- --target-gt q --new-gt M -i 'GT=="het"'
        bcftools index --threads 5 snp_full/${CHR}.F_IF-GF-MM1.vcf.gz

        echo "DIPLOID AUTOSOME"
        #all together
        bcftools view --samples-file Full_Samples_n302.list -Ou snp_full/${CHR}.IF-GF.vcf.gz | \
                bedtools subtract -header -a - -b ${mask} | \
                bcftools view -Oz -o snp_full/${CHR}.IF-GF-MM1.vcf.gz --min-alleles 2 --max-alleles 2 --min-ac 1 --threads 5 -i 'F_MISSING<0.1'
        bcftools index --threads 5 snp_full/${CHR}.IF-GF-MM1.vcf.gz

else
        echo "DIPLOID AUTOSOME"
        #all together
        bcftools view --samples-file Full_Samples_n302.list -Ou snp_full/${CHR}.IF-GF.vcf.gz | \
                bcftools view -Oz -o snp_full/${CHR}.IF-GF-MM1.vcf.gz --min-alleles 2 --max-alleles 2 --min-ac 1 --threads 5 -i 'F_MISSING<0.1'
        bcftools index --threads 5 snp_full/${CHR}.IF-GF-MM1.vcf.gz
fi

```

## Remove Relatives

Calculate relatedness on the full set 

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=5
#SBATCH --time=48:00:00

#merge autosomes
bcftools concat --file-list Autosomes.list --threads 10 -Oz -o merged_full/autos.vcf.gz
bcftools index --threads 10 merged_full/autos.vcf.gz

#LD prune
~/modules/plink2 --threads 20 --vcf merged_full/autos.vcf.gz --keep Canorus_Samples_n300.list --allow-extra-chr --set-missing-var-ids @:# \
        --rm-dup --indep-pairwise 50 5 0.1 --maf 0.05 --hwe 1e-10 --max-alleles 2 --min-alleles 2 --out merged_full/autos_canorus_LD
#extract, also a vcf
~/modules/plink2 --threads 20 --vcf merged_full/autos.vcf.gz --keep Canorus_Samples_n300.list --allow-extra-chr --set-missing-var-ids @:# \
        --extract merged_full/autos_canorus_LD.prune.in \
        --make-bed --recode vcf bgz --pca --out merged_full/autos_canorus_LD
bcftools index --threads 5 merged_full/autos_canorus_LD.vcf.gz
sed -i 's/chr_//g' merged_full/autos_canorus_LD.bim

#run admixture 
for K in {2..10}; do
echo "Runnign admixture at K: ${K}"
admixture -j7 --cv=5 ../merged_full/autos_canorus_LD.bed ${K} > autos_canorus_LD.log${K}.out
done

#calculated relatedness 
zcat merged_full/autos_canorus_LD.vcf.gz | sed 's/VCFv4.3/VCFv4.2/g' | vcftools --maf 0.05 --vcf - --relatedness2 --out merged_full/autos_canorus_LD

#calculate missingness
zcat merged_full/autos_canorus_LD.vcf.gz | sed 's/VCFv4.3/VCFv4.2/g' | vcftools --maf 0.05 --vcf - --missing-indv --out merged_full/autos_canorus_LD

#also calculate IBS0
plink --vcf merged_full/autos_canorus_LD.vcf.gz --genome full --const-fid --allow-extra-chr
awk '{print $2, $4, $15}' plink.genome > merged_full/autos_canorus_LD.ibs
```

Calculate the number of SNPs between individuals to design maternal / paternal relatives. This will output double the amount of SNPs because it considers both genotypes if it is encoded as diploid.

```python
import allel
import itertools
import numpy as np
import gzip

# Function to read a VCF file and return a GenotypeArray
def read_vcf(filename):
    callset = allel.read_vcf(filename)
    genotypes = allel.GenotypeArray(callset['calldata/GT'])
    return genotypes, callset['samples']

# Function to compute pairwise comparisons
def pairwise_comparisons(genotypes, samples):
    results = []
    for id1, id2 in itertools.combinations(range(len(samples)), 2):
        g1 = genotypes[:, id1]
        g2 = genotypes[:, id2]

        # Filter out missing data
        non_missing = ~(g1.is_missing() | g2.is_missing())
        g1 = g1.compress(non_missing, axis=0)
        g2 = g2.compress(non_missing, axis=0)

        # Count the number of sites
        number_sites = non_missing.sum()

        # Count the number of SNPs
        number_snps = np.count_nonzero(g1 != g2)

        results.append((samples[id1], samples[id2], number_sites, number_snps))
    return results

# Main script
if __name__ == '__main__':
    # Path to the VCF file
    vcf_path = 'chr_MT.IF-GF-MM1-AC2.vcf.gz'

    # Read the VCF file
    genotypes, samples = read_vcf(vcf_path)

    # Perform pairwise comparisons
    comparisons = pairwise_comparisons(genotypes, samples)

    # Output the results
    with open('pairwise_comparisons.txt', 'w') as f:
        for comparison in comparisons:
            f.write(f'{comparison[0]} {comparison[1]} {comparison[2]} {comparison[3]}\n')

    print('Pairwise comparisons have been saved to pairwise_comparisons.txt')

```

Identify relatives: 

```bash
#### Determine relatedness with PHI statistic 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/merged_full')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(viridis)
library(adegenet)
library(vcfR)
library(ggpubr)
library(scales)

# Analyze LD-pruned relatedness from vcftools --relatedness2
rel <- read.table('autos_canorus_LD.relatedness2',header=T) %>% select(1,2,7) %>% as_tibble
names(rel) <- c('ID_A','ID_B','PHI')
rel = rel %>% mutate(PHI = pmax(0, pmin(0.5, PHI)))
rel = rel %>% filter(ID_A != ID_B)

# Compare with plink IBS0
pk = read.table('autos_canorus_LD.ibs',sep=' ',header=TRUE); names(pk) = c('ID_A','ID_B','IBS0')
relpk = left_join(pk,rel)

# Values from here for designating relationships https://www.kingrelatedness.com/manual.shtml
fam = relpk %>% mutate(Relationship = ifelse(PHI > 0.354, 'First Degree',
                                             ifelse(PHI > 0.177, 'First Degree',
                                                    ifelse(PHI > 0.0884,'Second Degree',
                                                           ifelse(PHI > 0.0442, 'Third Degree',
                                                                  ifelse(PHI <= 0.0442,'Unrelated','Unassigned'))))))
fam %>% filter(PHI > 0.354) #usually > 0.354 is MZ twin / duplicate, but since our highest value is 0.393 and most around 0.37, seems more likely they are just first degree
phi_ibs = fam %>% ggplot(aes(x=PHI,y=IBS0,col=Relationship))+
  geom_point()+
  scale_color_viridis(discrete=TRUE)+
  theme_bw()
pdf('../figures/Relatedness__PhivIBS_2023OCT27.pdf',height=5,width=6)
phi_ibs
dev.off()

# Add mtDNA differences, calculated externally which identifies pairwise SNPs between samples for mtDNA
mt = read.table('chrMT_Pairwise_Differences.txt')
names(mt) = c('ID_A','ID_B','Sites','SNPs')
# The script counted diploid genotypes as 2 SNPs, so divide by 2
mt <- mt %>% mutate(SNPs = SNPs / 2)
fam2 = left_join(fam,mt) %>% drop_na(SNPs)

# Remove redundant comparisons, e.g. ID_A vs ID_B or pairwise redundancies
famrm = fam2 %>% 
  select(-IBS0) %>% 
  rowwise() %>% 
  mutate(pair = sort(c(ID_A,ID_B)) %>% paste(collapse = ",")) %>%
  group_by(pair,Sites,SNPs) %>%
  distinct(pair, .keep_all = T) %>% 
  separate(pair,into=c('ID_A','ID_B'),remove=F,sep=',') %>% ungroup() %>% select(-pair)

# Remove unrelated individuals ... first merge with metadata
md = read_tsv('~/merondun/cuculus_host/Metadata_Host.txt')
md = md %>% select(c(ID,HostParentShort,Habitat = Habitat_Optatus,Egg,KDist,Sampling_Year,Hap,Sex,Age))
mda = md

# Add an '_A' and '_B' to each metadata field 
names(mda) = paste0(names(mda),'_A')
mdb = md
names(mdb) = paste0(names(mdb),'_B')
fam3 = left_join(famrm,mda) %>% 
  left_join(.,mdb)

fam3 <- fam3 %>% filter(!grepl('148_',ID_A) & !grepl('148_',ID_B))

# Inspect, look to see which haps don't align, etc 
fam3 %>% filter(Relationship == 'First Degree' & Hap_A != Hap_B) %>% data.frame
fam3 %>% filter(Relationship == 'First Degree' & KDist_A != KDist_B) %>% data.frame
fam3 %>% filter(Relationship == 'First Degree' & Egg_A != Egg_B) %>% data.frame

# Only retain those with known egg morph, which are NESTLINGS, and which are caught within 2 years ( no parent / sibling)
famd = fam3 %>% drop_na(Egg_A) %>% drop_na(Egg_B)
famd = famd %>% filter(KDist_A == KDist_B) # Only compare within the same distance clade
#famd = famd %>% filter(Age_A == 'Young' & Age_B == 'Young') # Only compare within YOUNG NESTLINGS!
famd %>% count(Age_A,Age_B)
#famd = famd %>% filter(abs(Sampling_Year_A - Sampling_Year_B) <= 2) # Only compare within the same sampling year (e.g. SIBLINGS)
famd %>% count(Age_A,Age_B)

# Within first degree relatives, what's the distribution of the # of mtDNA snps? 
famd %>% filter(Relationship == 'First Degree') %>% 
  ggplot(aes(x=SNPs))+geom_histogram(show.legend = F)+theme_bw()+
  scale_x_continuous(breaks=pretty_breaks(n=14))

# If mtDNA haplotypes are identical, assign as maternal. 
famd = famd %>% mutate(Line = ifelse(SNPs <= 0,'Maternal','Paternal'))
famd %>% filter(Relationship == 'First Degree') %>%  count(Line)

#simply assign unrelated as paternal, move pie chart in final plot 
famd = famd %>% mutate(Line = ifelse(Relationship == 'Unrelated','Paternal',Line))

nrow(famd)
#[1] 2806 or [1] 2948 with no age and year restrictions
#how many unique individuals?
length(unique(c(famd$ID_A,famd$ID_B)))

#function to get matched data
get_matched_data <- function(df) {
  df %>% 
    mutate(
      Host = ifelse(HostParentShort_A == HostParentShort_B, 'Matched', 'Unmatched'),
      Habitat = ifelse(Habitat_A == Habitat_B, 'Matched', 'Unmatched'),
      Year = ifelse(abs(Sampling_Year_A - Sampling_Year_B) <= 2, 'Matched', 'Unmatched'),
      Egg = ifelse(Egg_A == Egg_B, 'Matched', 'Unmatched'),
      Distance = ifelse(KDist_A == KDist_B, 'Matched', 'Unmatched'),
      Haplogroup = ifelse(Hap_A == Hap_B, 'Matched', 'Unmatched')
    ) %>% 
    gather(key = "Variable", value = "Matched", Host, Habitat, Year, Egg, Distance,Haplogroup) %>% 
    group_by(Relationship,Line,Variable) %>% 
    count(Matched)
}

#Divide data into frames and store in a list
data_frames <- list(
  fem = famd %>% filter(Sex_A == 'F' & Sex_B == 'F'),
  mal = famd %>% filter(Sex_A == 'M' & Sex_B == 'M'),
  fm = famd %>% filter(Sex_A != Sex_B)
)

#function to each data frame and store results in a list
results_list <- lapply(names(data_frames), function(frame_name) {
  result_df <- get_matched_data(data_frames[[frame_name]])
  result_df$frame <- frame_name
  result_df
})

#combine
mm <- bind_rows(results_list)
mm$Relationship <- factor(mm$Relationship,levels=c('First Degree','Second Degree','Third Degree','Unrelated'))
mm <- mm %>%
  mutate(frame = case_when(
    frame == "fem" ~ "Female-Female",
    frame == "mal" ~ "Male-Male",
    frame == "fm" ~ "Intersexual",
    TRUE ~ frame)) %>% dplyr::rename(Count = n)

#proportions for pie charts
piece = mm %>% select(-frame) %>% group_by(Relationship,Line,Variable,Matched) %>% na.omit %>% summarize(Count = sum(Count)) %>% ungroup %>% group_by(Relationship,Line,Variable) %>% 
  mutate(Total = sum(Count),
         Proportion = Count/Total,
         Percent = paste0(round(Count/Total,3)*100,'% (',Total,')'))
piece = piece %>% filter(!grepl('Year|Distance|Hap',Variable))

#plot
relpie = piece %>% 
  ggplot(aes(x="",y=Proportion,fill=Matched))+
  geom_bar(stat='identity')+
  coord_polar("y", start=0)+xlab('')+ylab('')+
  facet_grid(Relationship~Variable+Line)+
  scale_fill_manual(values=c('forestgreen','grey95','black'))+
  theme_bw(base_size=7)+labs(fill='Phenotype Matching')+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
#add labels showing the % matched, and then the total sample size 
labs = piece %>% 
  filter(Matched == 'Matched')
piesR = relpie + 
  geom_label(data = labs, 
             aes(y=Inf,x=-Inf,label=Percent),fill='white',
             size=2,vjust=.3,hjust=.5,alpha=0.8)
piesR

pdf('../figures/20240717-Relatives_Phenotype_Matching-0Mismatch.pdf',height=4.5,width=6)
piesR
dev.off()

# 3RD DEGREE paternal relatives that don't match eggs: which egg types? 
# Egg_A Egg_B     n
# <chr> <chr> <int>
#   1 E1    E1       12
#   2 E10   E10       1
#   3 E4    E4        1
#   4 E6    E10      10
#   5 E6    E4        6
#   6 E6    E6       48

#also add missingness, we will remove individual with more missing data if both are females or both are males 
miss <- read.table('autos_canorus_LD.imiss',head=T) %>% select(c(1,5))
names(miss) = c('ID_A','Missing_A')
miss$Missing_A = as.numeric(miss$Missing_A)
miss = miss %>% na.omit
fam4 = fam3 %>% left_join(.,miss) %>% 
  merge(.,miss %>% dplyr::rename(ID_B=ID_A,Missing_B=Missing_A)) 

#Remove unrelateds
reli = fam4 %>% filter(Relationship != 'Unrelated')
relatives = fam4

#Unique samples to test 
samples <- unique(c(reli$ID_A,reli$ID_B))
rm <- NULL 
choice <- reli
#Re-run this multiple times, until it stops removing any individuals. 
for (run in seq(1,5,1)) { 
  
  cat('Running ',run,'\n')
  #grab one sample at a time 
  for (samp in samples) {
    #grab all the records with this sample 
    sf.a <- choice[grepl(samp,choice$ID_A),] 
    sf.b <- choice[grepl(samp,choice$ID_B),] 
    sf <- rbind(sf.a,sf.b) %>% arrange(desc(PHI))
    if(nrow(sf) < 1) next
    #if one if male remove, otherwise, if sampled clade is higher, remove, otherwise if missing data is higher, remove, otherwise random (sample by ID# higher#)
    sf1 <- sf %>% mutate(Remove = ifelse(Sex_A == 'F' & Sex_B == 'M', ID_B,
                                         ifelse(Sex_A == 'M' & Sex_B == 'F', ID_A,
                                                ifelse(Missing_A > Missing_B, ID_A,
                                                       ifelse(Missing_A < Missing_B, ID_B,
                                                              'Problem')))))
    sf1 %>% select(contains(c('ID')))
    #take the ID from the sample to remove
    bad <- head(na.omit(sf1$Remove),1) 
    #unless there are no samples to remove.. 
    if(length(bad) < 1) next
    #remove that bad sample from the pool
    cat('Removing bad sample: ',bad,'\n')
    choice <- choice[!grepl(bad,choice$ID_A),]
    choice <- choice[!grepl(bad,choice$ID_B),]
    #and then restart the whole process 
    rm <- rbind(rm,bad)
    rm(bad)
  }
}

#these are the bad individuals we will remove
rms <- as.data.frame(rm)
names(rms) <- 'Remove'
row.names(rms) <- NULL

#remove those bad IDs from the full dataset, make sure to filter form both ID_A and ID_B columns 
keep <- relatives %>% filter(!ID_A %in% rms$Remove) %>% filter(!ID_B %in% rms$Remove)
#see how many individuals were in the full dataset (sanity check), and then how many are retained after filtering 
length(unique(c(relatives$ID_A,relatives$ID_B)))
length(unique(c(keep$ID_A,keep$ID_B)))
length(unique(c(relatives$ID_A,relatives$ID_B))) - length(unique(rms$Remove))

write.table(unique(c(keep$ID_A,keep$ID_B)),'Unrelated_2023OCT27.list',quote=F,row.names=F,sep='/t',col.names=F)
samps = read_tsv('Unrelated_2023OCT27.list',col_names = F)

#calculate pairwise geographic distance
md = read_tsv('../Cuckoo_Full_Metadata_2023OCT3.txt')
ll = md %>% select(ID,Latitude,Longitude)
ll1 = geosphere::distm(ll %>% select(Longitude,Latitude)) %>% as.data.frame
names(ll1) = ll$ID
ll1$ID = ll$ID
ll2 = ll1 %>% pivot_longer(!(ID),names_to='ID_B',values_to = 'GDistance') %>% dplyr::rename(ID_A=ID)
ll2 = ll2 %>% mutate(GDistance = GDistance/1000)

#Add genetic distance
dp = left_join(fam4,ll2)
```

## SNP Filtering Unrelated: N = 202

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=5
#SBATCH --time=48:00:00

CHR=$1
GENOME=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa
mask=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/Coverage_Masks/GCA_017976375.1_bCucCan1.pri_genomic.CHR.N75-DoubleCoverage.mask.bed

echo "WORKING ON ${CHR}"

# Use bcftools to perform sample-wise filtering on the vcf file.
bcftools view --types snps --force-samples --samples-file Unrelated_2023OCT27.list -Oz ../merged/${CHR}.vcf.gz -o snp_unrelfull/${CHR}.vcf.gz
bcftools index --threads 5  snp_unrelfull/${CHR}.vcf.gz

# Calculate some statistics for the next round of filtering
SAMPN=$(bcftools query -l snp_unrelfull/${CHR}.vcf.gz | wc -l)
AVGDP=$(bcftools query -f '%DP\n' snp_unrelfull/${CHR}.vcf.gz | datamash median 1 | datamash round 1) # Compute the median coverage depth across all sites in the vcf file
DPHI=$(($AVGDP*2)) # Calculate thresholds for filtering based on depth of coverage

# Apply filters
bcftools view --types snps --threads 5 -e "QUAL < 20 || INFO/DP > $DPHI || INFO/DP < $SAMPN || MQ < 30 || RPBZ < -3 || RPBZ > 3" -Oz -o snp_unrelfull/${CHR}.IF.vcf.gz snp_unrelfull/${CHR}.vcf.gz
bcftools index --threads 5 snp_unrelfull/${CHR}.IF.vcf.gz

# For diploid chromosomes, convert low-coverage sites to missing
MINDP=3
bcftools +setGT -Oz -o snp_unrelfull/${CHR}.IF-GF.vcf.gz snp_unrelfull/${CHR}.IF.vcf.gz -- -t q -i "FMT/DP < ${MINDP}" -n "./."
bcftools index --threads 5 snp_unrelfull/${CHR}.IF-GF.vcf.gz

#for W and MT, set allele imbalance violations (heterozygosity) to missing for FEMALE, no W for MALE, diploid Z for MALE
if [[ $CHR == 'chr_MT' ]]
then
        echo "HAPLOID, SETTING HETEROZYGOUS SITES TO MISSING FOR BOTH SEXES"
        #All samples
        bcftools view --samples-file Unrelated_2023OCT27.list -Ou snp_unrelfull/${CHR}.IF-GF.vcf.gz | \
                bcftools view --min-alleles 2 --max-alleles 2 --min-ac 1 --threads 5 -i 'F_MISSING<0.1' | \
                bedtools subtract -header -a - -b ${mask} | \
                bcftools +setGT -Ou -- -t "b:AD<1e-5" -n "./." | \
            bcftools +setGT -Oz -o snp_unrelfull/${CHR}.IF-GF-MM1.vcf.gz -- --target-gt q --new-gt M -i 'GT=="het"'
        bcftools index --threads 5 snp_unrelfull/${CHR}.IF-GF-MM1.vcf.gz

        echo "HAPLOID, SETTING HETEROZYGOUS SITES TO MISSING FOR BOTH SEXES"
        #Females only
        bcftools view --samples-file FUnrelated_2023OCT27.list -Ou snp_unrelfull/${CHR}.IF-GF.vcf.gz | \
                bcftools view --min-alleles 2 --max-alleles 2 --min-ac 1 --threads 5 -i 'F_MISSING<0.1' | \
                bedtools subtract -header -a - -b ${mask} | \
                bcftools +setGT -Ou -- -t "b:AD<1e-5" -n "./." | \
            bcftools +setGT -Oz -o snp_unrelfull/${CHR}.F_IF-GF-MM1.vcf.gz -- --target-gt q --new-gt M -i 'GT=="het"'
        bcftools index --threads 5 snp_unrelfull/${CHR}.F_IF-GF-MM1.vcf.gz

elif [[ $CHR == 'chr_W' ]]
then
        #ONLY FEMALES
        echo "HAPLOID, SETTING HETEROZYGOUS SITES TO MISSING AND ONLY GRAB FEMALES "
        bcftools view --samples-file FUnrelated_2023OCT27.list -Ou snp_unrelfull/${CHR}.IF-GF.vcf.gz | \
                bcftools view --min-alleles 2 --max-alleles 2 --min-ac 1 --threads 5 -i 'F_MISSING<0.1' | \
                bedtools subtract -header -a - -b ${mask} | \
                bcftools +setGT -Ou -- -t "b:AD<1e-5" -n "./." | \
            bcftools +setGT -Oz -o snp_unrelfull/${CHR}.IF-GF-MM1.vcf.gz -- --target-gt q --new-gt M -i 'GT=="het"'
        bcftools index --threads 5 snp_unrelfull/${CHR}.IF-GF-MM1.vcf.gz

elif [[ $CHR == 'chr_Z' ]]
then
        echo "DIPLOID AUTOSOME"
        #all together
        bcftools view --samples-file MUnrelated_2023OCT27.list -Ov snp_unrelfull/${CHR}.IF-GF.vcf.gz | \
                bedtools subtract -header -a - -b ${mask} | \
                bcftools view -Oz -o snp_unrelfull/${CHR}.M_IF-GF-MM1.vcf.gz --min-alleles 2 --max-alleles 2 --min-ac 1 --threads 5 -i 'F_MISSING<0.1'
        bcftools index --threads 5 snp_unrelfull/${CHR}.M_IF-GF-MM1.vcf.gz

        #Females second
        echo "HAPLOID, SETTING HETEROZYGOUS SITES TO MISSING AND ONLY GRAB FEMALES "
        bcftools view --samples-file FUnrelated_2023OCT27.list -Ou snp_unrelfull/${CHR}.IF-GF.vcf.gz | \
        bcftools view --min-alleles 2 --max-alleles 2 --min-ac 1 --threads 5 -i 'F_MISSING<0.1' | \
                bedtools subtract -header -a - -b ${mask} | \
                bcftools +setGT -Ou -- -t "b:AD<1e-5" -n "./." | \
                bcftools +setGT -Oz -o snp_unrelfull/${CHR}.F_IF-GF-MM1.vcf.gz -- --target-gt q --new-gt M -i 'GT=="het"'
        bcftools index --threads 5 snp_unrelfull/${CHR}.F_IF-GF-MM1.vcf.gz

        echo "DIPLOID AUTOSOME"
        #all together
        bcftools view --samples-file Unrelated_2023OCT27.list -Ov snp_unrelfull/${CHR}.IF-GF.vcf.gz | \
                bedtools subtract -header -a - -b ${mask} | \
                bcftools view -Oz -o snp_unrelfull/${CHR}.IF-GF-MM1.vcf.gz --min-alleles 2 --max-alleles 2 --min-ac 1 --threads 5 -i 'F_MISSING<0.1'
        bcftools index --threads 5 snp_unrelfull/${CHR}.IF-GF-MM1.vcf.gz

else
        echo "DIPLOID AUTOSOME"
        #all together
        bcftools view --samples-file Unrelated_2023OCT27.list -Ou snp_unrelfull/${CHR}.IF-GF.vcf.gz | \
                bcftools view -Oz -o snp_unrelfull/${CHR}.IF-GF-MM1.vcf.gz --min-alleles 2 --max-alleles 2 --min-ac 1 --threads 5 -i 'F_MISSING<0.1'
        bcftools index --threads 5 snp_unrelfull/${CHR}.IF-GF-MM1.vcf.gz
fi

```

## Merge Autosomes & LD Pruning

Merge entire dataset: 

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=10
#SBATCH --time=200:00:00

#merge autosomes
bcftools concat --file-list Autosomes.list --threads 10 -Oz -o merged_unrelfull/autos.vcf.gz
bcftools index --threads 10 merged_unrelfull/autos.vcf.gz

#LD prune
~/modules/plink2 --threads 10 --vcf merged_unrelfull/autos.vcf.gz --keep CanorusUnrelated_2023OCT27.list --allow-extra-chr --set-missing-var-ids @:# \
        --rm-dup --indep-pairwise 50 5 0.1 --maf 0.05 --hwe 1e-10 --max-alleles 2 --min-alleles 2 --out merged_unrelfull/autos_canorus_LD
        
#extract, also a vcf and run PCA 
~/modules/plink2 --threads 10 --vcf merged_unrelfull/autos.vcf.gz --keep CanorusUnrelated_2023OCT27.list --allow-extra-chr --set-missing-var-ids @:# \
        --extract merged_unrelfull/autos_canorus_LD.prune.in \
        --make-bed --recode vcf bgz --pca --out merged_unrelfull/autos_canorus_LD
bcftools index --threads 10 merged_unrelfull/autos_canorus_LD.vcf.gz
sed -i 's/chr_//g' merged_unrelfull/autos_canorus_LD.bim
```

 

# Analyses: Matrilineal Phylogenetics

## Output Maternal Fastas 

Subset samples, create tree:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=20
#SBATCH --time=24:00:00

#mamba activate snps
#to submit all iterations: for i in $(cat CHRS.list); do for j in $(cat SUBSETS.list); do sbatch -J TREE_${i}_${j} ~/merondun/cuculus_host/phylogenetics/1.Subset_Samples_Filter.sh ${i} ${j} ; done ; done
mkdir vcfs ml_trees

#mask with male-biased coverage 
mask=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/Coverage_Masks/GCA_017976375.1_bCucCan1.pri_genomic.CHR.N75-DoubleCoverage.mask.bed

CHR=$1
SAMP=$2

#minimum coverage, LESS than this set to missing 
MINDP=3
echo "FILTERING AND MERGING VARIANT SITES FOR ${CHR}"
bcftools view --threads 10 --samples-file ${CHR}_${SAMP}.list -Ou ../../../merged/${CHR}.vcf.gz | \
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

#filter, include singletons 
bcftools view vcfs/${CHR}_${SAMP}.SNP.DP3.vcf.gz --min-alleles 2 --max-alleles 2 --min-ac 1 --max-af 0.999 --types snps -i "MQ>40 & INFO/DP > ${low} & INFO/DP < ${high} & F_MISSING < 0.1" -Oz -o vcfs/${CHR}_${SAMP}.SNP.DP3-AC1-MQ40.vcf.gz
bcftools index --threads 10 vcfs/${CHR}_${SAMP}.SNP.DP3-AC1-MQ40.vcf.gz 

#create tree  
python ~/modules/vcf2phylip.py -i vcfs/${CHR}_${SAMP}.SNP.DP3-AC1-MQ40.vcf.gz -f --output-folder ml_trees
iqtree --redo -keep-ident -T 20 -s ml_trees/${CHR}_${SAMP}.SNP.DP3-AC1-MQ40.min4.phy --seqtype DNA -m "MFP+ASC" -alrt 1000 -B 1000
iqtree --redo -keep-ident -T 20 -s ml_trees/${CHR}_${SAMP}.SNP.DP3-AC1-MQ40.min4.phy.varsites.phy --seqtype DNA -m "MFP+ASC" -alrt 1000 -B 1000

#filter with NO SINGLETONS 
bcftools view vcfs/${CHR}_${SAMP}.SNP.DP3.vcf.gz --min-alleles 2 --max-alleles 2 --min-ac 2 --max-af 0.999 --types snps -i "MQ>40 & INFO/DP > ${low} & INFO/DP < ${high} & F_MISSING < 0.1" -Oz -o vcfs/${CHR}_${SAMP}.SNP.DP3-AC2-MQ40.vcf.gz
bcftools index --threads 10 vcfs/${CHR}_${SAMP}.SNP.DP3-AC2-MQ40.vcf.gz 

#create tree NO SINGLETONS   
python ~/modules/vcf2phylip.py -i vcfs/${CHR}_${SAMP}.SNP.DP3-AC2-MQ40.vcf.gz -f --output-folder ml_trees
iqtree --redo -keep-ident -T 20 -s ml_trees/${CHR}_${SAMP}.SNP.DP3-AC2-MQ40.min4.phy --seqtype DNA -m "MFP+ASC" -alrt 1000 -B 1000
iqtree --redo -keep-ident -T 20 -s ml_trees/${CHR}_${SAMP}.SNP.DP3-AC2-MQ40.min4.phy.varsites.phy --seqtype DNA -m "MFP+ASC" -alrt 1000 -B 1000
```

## Tree Sample Sensitivity

Plot many trees:

```R
#### Plot many W and MT trees with different sample subsets 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host')
.libPaths('~/mambaforge/envs/R/lib/R/library')
library(ggtree)
library(ape)
library(tidyverse)
library(treeio)
library(viridis)
library(ggpubr)
library(RColorBrewer)
library(phytools)
library(ggnewscale)

md = read_tsv('~/merondun/cuculus_host/Metadata_Host.txt')

###### CIRCULAR W TREE  #######
#plot tree 
trees = list.files('beast_dating/variant_only/ml_trees/',pattern='.*AC2.*contree',full.names = TRUE)
counter = 0 
for (tree in trees) { counter = counter + 1;
iqtree = read.iqtree(tree)
iqtr = midpoint.root(as.phylo(iqtree))
lab = gsub('-MQ40.*','',gsub('.*chr_','',tree))
cat ('Making tree: ',lab,'\n')

#plot with outgroups 
circ = ggtree(iqtr, layout = "circular",branch.length='none') %<+% md  +
  geom_nodepoint(mapping=aes(subset=(as.numeric(label) >= 95)),col='black',fill='grey90',pch=23,size=0.75,show.legend=F)+
  geom_tippoint(aes(fill = Hap,shape=Species),size=2)+
  scale_fill_manual('W Haplogroup',values=md$HapCol,breaks=md$Hap)+
  scale_shape_manual(values=c(21,4,8))+
  guides(fill=guide_legend(nrow=4,override.aes=list(shape=21)),
         shape=guide_legend(nrow=5))+
  ggtitle(lab)+
  theme(legend.text = element_text(size = 8),legend.title = element_text(size = 10),
        legend.key.size = unit(0.2, "cm"),    legend.position = 'top')
circ

apetree = ggtree(iqtr, layout = "ape") %<+% md  +
  geom_nodepoint(mapping=aes(subset=(as.numeric(label) >= 95)),col='black',fill='grey90',pch=23,size=0.75,show.legend=F)+
  geom_tippoint(aes(fill = Hap,shape=Species),size=1.25)+
  scale_fill_manual('W Haplogroup',values=md$HapCol,breaks=md$Hap)+
  scale_shape_manual(values=c(21,4,8))+
  guides(fill=guide_legend(nrow=4,override.aes=list(shape=21)),
         shape=guide_legend(nrow=5))+
  theme(legend.text = element_text(size = 8),legend.title = element_text(size = 10),
        legend.key.size = unit(0.2, "cm"),    legend.position = 'top')
apetree

ca = ggarrange(circ,apetree,nrow=2,heights=c(0.6,0.4),legend = 'none') 
assign(paste0('p',counter),ca)

}

pdf('figures/ML_Trees-VariableSampleSubsets-All-UnrelMNLR-BEAST_2024MAR14.pdf',height=14,width=7)
ggarrange(p1,p5,p2,p6,p4,p8,common.legend = TRUE,nrow=3,ncol=2)
dev.off()


```



## Assign Haplogroups

Using fasta files from VCF, clsuter using tree:

```R
#### Plot Trees, Assign Haplogroups 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host')
.libPaths('~/mambaforge/envs/R/lib/R/library')
library(ggtree)
library(ape)
library(tidyverse)
library(treeio)
library(viridis)
library(adegenet)
library(ggpubr)
library(vcfR)
library(adegenet)
library(RColorBrewer)
library(phytools)
library(ggnewscale)

md = read_tsv('Cuckoo_Full_Metadata_2023OCT3.txt')

###### W #######
## Plot Tree 
iqtree = read.iqtree('beast_dating/variant_only/ml_trees/chr_W_All.SNP.DP3-AC2-MQ40.min4.phy.contree')
iqtr = midpoint.root(as.phylo(iqtree))

#plot cladogram, first find nodes to collapse 
insp = ggtree(iqtr, layout = "dendrogram") %<+% md  + 
  geom_nodelab(aes(label=node),geom = 'label',size=3)+  geom_nodepoint(mapping=aes(subset=(as.numeric(label) >= 95)),col='black',fill='grey90',pch=23,size=1.5,show.legend=F)+
  geom_tippoint(aes(fill = KDist,shape=Plumage),size=5)+
  scale_shape_manual(values=c(24,21,25))+
  #scale_fill_manual('W Haplogroup',values=md$Wcol,breaks=md$W)+
  scale_fill_manual('Distance Group',values=md$KDCol,breaks=md$KDist)+
  geom_tiplab()+
  guides(fill=guide_legend(nrow=2,override.aes=list(shape=21)))
insp
pdf('figures/W_TREE-Inspection-AC2_2024FEB27.pdf',height=25,width=55)
insp
dev.off()

#find the nodes in which all descedents will be assigned a haplogroup 
labs = data.frame(Nodes = c(219,201,165,308,290,277,231))
ds = NULL; counter = 0
for (node in labs$Nodes){ counter = counter + 1 
d = data.frame(ID = get_taxa_name(insp, node = node), W = paste0('W',counter))
ds = rbind(ds,d)
}

#re-inspect 
inspW = ggtree(iqtr, layout = "ape") %<+% ds2  + 
  geom_tippoint(aes(col = W),size=2,pch=16)+
  scale_color_manual('W Haplogroup',values=ds2$Wcol,breaks=ds2$W)+
  geom_treescale(offset = 0.01)
inspW
insp$data %>% filter(is.na(W) & isTip == TRUE)

#add colors
cols = ds %>% select(W) %>% unique() %>% mutate(Wcol = brewer.pal(7,'Paired'))
ds2 = left_join(ds,cols)
write_tsv(ds2,file='figures/W_Haplogroups_VCF_2023DEC19.txt')
```

## Plot Collapsed Tree

Collapse chrW tree on supported nodes for visualization (fig 1)

```R
#### Collapse chrW Tree 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host')
.libPaths('~/mambaforge/envs/R/lib/R/library')
library(ggtree)
library(ape)
library(tidyverse)
library(treeio)
library(viridis)
library(adegenet)
library(ggpubr)
library(vcfR)
library(adegenet)
library(RColorBrewer)
library(phytools)
library(ggnewscale)

md = read_tsv('~/merondun/cuculus_host/Metadata_Host.txt')

##### COLLAPSE W TREE #####
#midpoint root 
iqtree = read.iqtree('beast_dating/variant_only/ml_trees/chr_W_Host.SNP.DP3-AC2-MQ40.min4.phy.contree')
iqtr = midpoint.root(as.phylo(iqtree))
p = ggtree(iqtr, layout = "dendrogram") %>% ggtree::rotate(65) %>% ggtree::rotate(74) %>% ggtree::rotate(80) %<+% md 
pdf('figures/W_TREE-COLLAPSE_Inspection_2024FEB27.pdf',height=75,width=55)
p + geom_tippoint(aes(fill = Hap),pch=21,size=4)+scale_fill_manual(values=md$HapCol,breaks=md$Hap) + 
  geom_nodelab(aes(label=node),geom = 'label',size=2) 
dev.off()
hapcols = md %>% select(Hap,HapCol) %>% unique %>% arrange(Hap) %>% filter(!grepl('CM|CP',Hap)) %>% na.omit
#collapse
p = ggtree(iqtr, layout = "dendrogram") %>% ggtree::rotate(65) %>% ggtree::rotate(74) %>% ggtree::rotate(80) %<+% md + 
  geom_tippoint(aes(fill=Hap),size=3,pch=21) +
  scale_fill_manual(values=md$HapCol,breaks=md$Hap) + 
  geom_nodepoint(mapping=aes(subset=(as.numeric(label) >= 95)),col='black',fill='grey90',pch=23,size=1.5,show.legend=F)
p2 = p %>% collapse(node=75) + geom_point2(aes(subset=(node==75)), shape=21, size=3, fill=hapcols[1,2]); p2
p2 = collapse(p2,node=72) + geom_point2(aes(subset=(node==72)), shape=21, size=3, fill=hapcols[[2,2]]); p2
p2 = collapse(p2,node=108) + geom_point2(aes(subset=(node==108)), shape=21, size=3, fill=hapcols[[3,2]]); p2
#p2 = collapse(p2,node=) + geom_point2(aes(subset=(node==292)), shape=21, size=3, fill=hapcols[[4,2]]);p2
p2 = collapse(p2,node=105) + geom_point2(aes(subset=(node==105)), shape=21, size=3, fill=hapcols[[5,2]]);p2
p2 = collapse(p2,node=103) + geom_point2(aes(subset=(node==103)), shape=21, size=3, fill=hapcols[[6,2]]);p2
p2 = collapse(p2,node=83) + geom_point2(aes(subset=(node==83)), shape=21, size=3, fill=hapcols[[7,2]]);p2
p2 = p2 + geom_treescale(x=0.05,offset = 0.01) 

pdf('~/merondun/cuculus_host/phylogenetics/ML_Tree_chr_W-Collapse_Haplogroups_2024FEB27.pdf',height=1,width=6)
p2
dev.off()

gp = ggtree(iqtr, layout = "dendrogram") %>% ggtree::rotate(65) %>% ggtree::rotate(74) %>% ggtree::rotate(80) %<+% md+ 
  #geom_nodelab(aes(label=node),geom = 'label',size=3)+
  geom_nodepoint(mapping=aes(subset=(as.numeric(label) >= 95)),col='black',fill='grey90',pch=23,size=1.5,show.legend=F) +
  geom_tippoint(aes(fill = Hap,shape=SpeciesShort),size=2.5,alpha=0.95)+
  scale_fill_manual(values=md$HapCol,breaks=md$Hap)+
  scale_shape_manual(values=c(21,4,8))+
  geom_treescale(x=0.05,offset = 0.01)+
  theme(legend.position='none')
gp

pdf('~/merondun/cuculus_host/phylogenetics/W_TREE-Ape_AllSamples_2024FEB27.pdf',height=1.75,width=5.5)
gp
dev.off()

```

## BEAST

Following this [tutorial](https://beast2-dev.github.io/beast-docs/beast2/DivergenceDating/DivergenceDatingTutorial.html).

From the host-type fastas:

```bash
seqret -sequence ml_trees/chr_MT_HostOG.SNP.DP3-AC2-MQ40.min4.fasta -outseq nexus/chr_MT_HostOG.SNP.DP3-AC2-MQ40.nex -osformat nexus
seqret -sequence ml_trees/chr_W_HostOG.SNP.DP3-AC2-MQ40.min4.fasta -outseq nexus/chr_W_HostOG.SNP.DP3-AC2-MQ40.nex -osformat nexus
```

For BEAUTI:

```bash
#for beauti
4 gamma categories, HKY, empirical sub rate 
clock rate strict at 5.05E-9
coalescent exponential, priors with lognormal CM / CP ancestor at M = 1.7  S=0.2
100M chains, log every 1k 
```

for SNAPP from this [tutorial](https://github.com/mmatschiner/tutorials/blob/master/divergence_time_estimation_with_snp_data/README.md) 

```bash
#chr_W, run 2 iterations of this with a different SNP sample
bcftools view ../vcfs/chr_W_HostOG.SNP.DP3-AC2-MQ40.vcf.gz -Ov -o chr_W_HostOG.SNP.DP3-AC2-MQ40.vcf
ruby snapp_prep.rb -o chr_W -v chr_W_HostOG.SNP.DP3-AC2-MQ40.vcf -t chr_W_samples.txt -c constraints.txt -m 1000 -l 100000
mv snapp.xml chr_W.xml

#also re-run with 5K SNPs 
ruby snapp_prep.rb --xml chr_W-BC.xml -o chr_W-BC -v chr_W_HostOG.SNP.DP3-AC2-MQ40.vcf -t chr_W_samples.txt -c constraints_beast.txt -m 1000 -l 500000

#and second replicate
ruby snapp_prep.rb --xml chr_W-BC2.xml -o chr_W-BC2 -v chr_W_HostOG.SNP.DP3-AC2-MQ40.vcf -t chr_W_samples.txt -c constraints_beast.txt -m 1000 -l 500000

#and with double SNPs 
ruby snapp_prep.rb --xml chr_W-2KBC2.xml -o chr_W-2KBC2 -v chr_W_HostOG.SNP.DP3-AC2-MQ40.vcf -t chr_W_samples.txt -c constraints_beast.txt -m 2000 -l 500000
```

For the SNAPP constraints file I will specify the crown according to the whole tree, with a lenient prior corresponding to a 95% CI between 3.06 - 8.16 MYA

```bash
#Constraints give this CI: 3.06 - 8.16 MYA
meanlog <- log(5)  # This is an approximation; see the explanation above
sdlog <- 0.25

# Calculate the 95% CI on the log scale
ci_lower <- qlnorm(0.025, meanlog=meanlog, sdlog=sdlog)
ci_upper <- qlnorm(0.975, meanlog=meanlog, sdlog=sdlog)

ci_lower
# [1] 3.06316
ci_upper
# [1] 8.161508

#file
cat constraints.txt
lognormal(0,5,0.25)  crown CP,CM,W1,W2,W3,W4,W5,W6,W7
```

Run BEAST & Annotate:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=20
#SBATCH --time=200:00:00

RUN=$1
beast -threads 20 -overwrite -java ${RUN}.xml
```

For adding invariant sites: https://groups.google.com/g/beast-users/c/QfBHMOqImFE

Count AGCT:

```bash
grep -v "^>" chr_W.fa | fold -w1 | sort | uniq -c

chr_W:
T	G	C	A
6108268	4562746	4578162	6088177

chr_MT:
T	G	C	A
4749	2471	5775	6703

e.g. 
   <data id='chr_W_HostOG' spec='FilteredAlignment' filter='-' data='@chr_W_HostOGOriginal' constantSiteWeights='6088177 4578162 4562746 6108268'/>
```

Info from the BEAST google groups help forum: 

```bash
You can use a FilteredAlignment to insert constant sites and set the constantSiteWeights attribute. Say, your original alignment is called xyz, so the XML produced by BEAUti contains something like

    <data id="xyz" name="alignment">

It is easiest to rename this to say xyzOriginal,

    <data id="xyzOriginal" name="alignment">

then add another data element, just after the closing </data> element of the alignment would be a good spot, that say

   <data id='xyz' spec='FilteredAlignment' filter='-' data='@xyzOriginal' constantSiteWeights='100 200 300 400'/>

Note id='xyz' and data='@xyzOriginal' should match what you have in the XML.

The constant weights at the end add weights for DNA in order A,C,G,T, so it adds 100 constant sites with all As, 200 with all Cs etc.

In the output to screen, it should report statistics of the xyzOriginal as something like:

6 taxa
768 sites
69 patterns

followed by statistics of the filtered alignment

Filter -
6 taxa
768 sites + 1000 constant sites
69 patterns

where the total number of constant sites added are reported as well.
```

## Plot: SNAPP

```bash
#### Plot SNAPP annotated trees 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/beast_dating/variant_only/snapp/1Mchains/')
.libPaths('~/mambaforge/envs/R/lib/R/library')
library(ggtree)
library(phytools)
library(ape)
library(treeio)
library(viridis)
library(ggpubr)
library(RColorBrewer)
library(tidyverse)

#Constraints give this CI: 3.06 - 8.16 MYA
meanlog <- log(5)  # This is an approximation; see the explanation above
sdlog <- 0.25

# Calculate the 95% CI on the log scale
ci_lower <- qlnorm(0.025, meanlog=meanlog, sdlog=sdlog)
ci_upper <- qlnorm(0.975, meanlog=meanlog, sdlog=sdlog)

ci_lower
ci_upper

#Read in metadata
md = read_tsv('~/merondun/cuculus_host/Metadata_Host.txt')

files = list.files('.',paste0('.*ann'))
file = 'chr_W.trees.ann'
counter = 0
for (file in files){
  counter = counter +  1 
  iqtree = read.beast(file)
  if (counter == 1) { 
    gg = ggtree(iqtree,layout='dendrogram') %<+% md 
    } 
  else {
    #Rotate the nodes for consistency 
    gg = ggtree(iqtree,layout='dendrogram') %<+% md %>% ggtree::rotate(10) %>% ggtree::rotate(15) %>% ggtree::rotate(13)
  }
  
  #for rotating clades
  #ggtree(iqtree,layout='dendrogram') + geom_nodelab(aes(label=node),size=20)
      
  #add label for 95% CIs
  lab = gsub('.trees.*','',file)
  heights = gg$data$height_0.95_HPD
  df = as.data.frame(do.call(rbind, heights)) #convert the list to a data frame
  df$node_value = 1:nrow(df) # Add node values as a new column
  colnames(df) = c("value1", "value2", "node")
  df = df[, c("node", "value1", "value2")]
  df = df %>% mutate(lab = paste0(round(value1,1),' - ',round(value2,1))) %>% select(!c(value1,value2))
  
  leg = md %>% select(Hap,HapCol,SpeciesShort) %>% unique %>% drop_na(Hap) %>% mutate(shape = c(rep(21,6),4,21,8))
  gg$data = left_join(gg$data,df)
  ggp = gg  +
    geom_range(range='height_0.95_HPD', color='red', alpha=.6, size=1) +
    geom_tippoint(aes(fill = label,shape=label),size=1.5)+
    geom_nodelab(aes(label=lab),size=1.5,vjust=1) +
    ggtitle(lab)+
    #geom_tiplab(size=2)+
    #geom_nodelab(aes(x=branch, label=round(posterior, 2)), vjust=-.5, size=3) +
    scale_color_continuous(low="darkgreen", high="red") +
    scale_fill_manual(values=leg$HapCol,breaks=leg$Hap)+
    scale_shape_manual(values=leg$shape,breaks=leg$Hap)+
    theme(legend.position=c(.1, .8))+
    geom_treescale(x = 1)+
    guides(fill=guide_legend(override.aes=list(shape=21)))+
    theme(legend.position='right')
  ggp
  assign(paste0('p',counter),ggp)
} 

ggarrange(p1,p2,common.legend = TRUE)

pdf('../../../../figures/SNAPP_Divergence_Dating_All_2024MAR16.pdf',height=3.5,width=7)
ggarrange(p1,p2,common.legend = TRUE)
dev.off()

```

## Plot: BEAST

```R
#### Plot BEAST annotated trees 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/beast_dating/variant_only/nexus/')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(ggtree)
library(phytools)
library(ape)
library(treeio)
library(viridis)
library(ggpubr)
library(RColorBrewer)
library(tidyverse)

#Read in metadata
md = read_tsv('~/merondun/cuculus_host/Metadata_Host.txt')

files = list.files('.',paste0('.*ann'))

counter = 0
for (file in files){
  counter = counter +  1 
  iqtree = read.beast(file) 
  gg = ggtree(iqtree,layout='rectangular') %<+% md
  
  #add label for 95% CIs
  lab = gsub('.trees.*','',file)
  heights = gg$data$height_0.95_HPD
  df = as.data.frame(do.call(rbind, heights)) #convert the list to a data frame
  df$node_value = 1:nrow(df) # Add node values as a new column
  colnames(df) = c("value1", "value2", "node")
  df = df[, c("node", "value1", "value2")]
  df = df %>% mutate(lab = paste0(round(value1,1),' - ',round(value2,1))) %>% select(!c(value1,value2))
  
  leg = md %>% select(Hap,HapCol,SpeciesShort) %>% unique %>% drop_na(Hap)
  gg$data = left_join(gg$data,df)
  ggp = gg  +
    geom_range(range='height_0.95_HPD', color='red', alpha=.6, size=1) +
    geom_tippoint(aes(fill = Hap,shape=SpeciesShort),size=1.5)+
    geom_nodelab(aes(label=lab),size=1.5,vjust=1) +
    ggtitle(lab)+
    #geom_tiplab(size=2)+
    #geom_nodelab(aes(x=branch, label=round(posterior, 2)), vjust=-.5, size=3) +
    scale_color_continuous(low="darkgreen", high="red") +
    scale_fill_manual(values=leg$HapCol,breaks=leg$Hap)+
    scale_shape_manual(values=c(21,4,8))+
    theme(legend.position=c(.1, .8))+
    geom_treescale(x = 5)+
    guides(fill=guide_legend(override.aes=list(shape=21)))+
    theme(legend.position='right')
  ggp
  assign(paste0('p',counter),ggp)
} 

ggarrange(p1,p2,p3,p4,common.legend = TRUE)

pdf('../../../figures/BEAST_Divergence_Dating_All_2024MAR14.pdf',height=9,width=7)
ggarrange(p1,p2,p3,p4,common.legend = TRUE)
dev.off()

```







To get the sequence for analysis, take the individuals from above subset across the blue egg groups, then:

* Grab W chromosome genes between 10-75KB
* Find open reading frames in those genes shared among all samples
* Extract those ORFS, recombine them into a single fasta for codon positions in beauti
* Randomly subset 1 individual from each population as a representative, create 4 iterations to assess sensitivity 

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=24:00:00

#mamba activate beast
#Run this in an empty directory. Specify the paths to the standard files and it will output a file for each gene (indicated by name in field 4 of the bed annotation), with the sample names from the vcf as the sequence header. This is primarily useful for haploidized VCFs e.g. chrW or mtDNA

mkdir vcfs work nex
#subset vcf
mask=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/Coverage_Masks/GCA_017976375.1_bCucCan1.pri_genomic.CHR.N75-DoubleCoverage.mask.bed

MINDP=3
echo "FILTERING AND MERGING VARIANT SITES FOR chr_W"
bcftools view --threads 10 --samples-file BEAST_AllSamples.list -Ou ../../../merged/chr_W.vcf.gz | \
        bcftools view --types snps --min-alleles 2 --max-alleles 2 --threads 10 | \ #grab snps
        bedtools subtract -header -a - -b ${mask} | \ #remove male coverage sites 
        bcftools +setGT -- -t q -i "FMT/DP < ${MINDP}" -n "./." | \ #set sites below 2x to missing
        bcftools +setGT -Ou -- -t "b:AD<1e-5" -n "./." | \ #set heterozygous genotypes to missing
        bcftools +setGT -Ou -- --target-gt q --new-gt M -i 'GT=="het"' | \ #set weakly het to major allele 
        bcftools +fixploidy -Ou - -- -f 1 | \ #set to haplid
        bcftools +fill-tags -Ou -- -t AC,AN | \ #update AC fields 
        bcftools view --min-alleles 2 --max-alleles 2 --types snps --min-ac 1 -Oz -o vcfs/chr_W.SNP.DP5.vcf.gz 
bcftools index --threads 10 vcfs/chr_W.SNP.DP5.vcf.gz 

#also filter on DP 
bcftools query -f '%CHROM\t%POS\t%INFO/DP\t%INFO/MQ\n' vcfs/chr_W.SNP.DP5.vcf.gz > dp_stats.txt

#IN R:
library(tidyverse); library(meRo)
d = read_tsv('dp_stats.txt',col_names=F)
#DP
ds = d %>% sum_stats(X3)
   mean    sd    se median   iqr conf_low conf_high
1 1414.  354. 0.706   1465   300    1413.     1416.
ds$mean-ds$sd*2
[1] 706.2658
ds$mean+ds$sd*2
[1] 2122.033

#filter 
bcftools view vcfs/chr_W.SNP.DP5.vcf.gz --min-alleles 2 --max-alleles 2 --min-ac 2 --max-af 0.999 --types snps -i 'MQ>40 & INFO/DP > 706 & INFO/DP < 2122 & F_MISSING < 0.1' -Oz -o vcfs/chr_W.SNP.DP5-AC2-MQ40.vcf.gz

#with no missing data 
bcftools view vcfs/chr_W.SNP.vcf.gz -i 'F_MISSING=0' | bcftools view --min-alleles 2 --max-alleles 2 --min-ac 2 --max-af 0.999 --types snps -i 'MQ>40 & INFO/DP > 706 & INFO/DP < 2122' -Oz -o vcfs/chr_W.SNP.DP5-AC2-MQ40-MM0.vcf.gz
```

Full chrW:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00

#mamba activate beast
#Run this in an empty directory. Specify the paths to the standard files and it will output a file for each gene (indicated by name in field 4 of the bed annotation), with the sample names from the vcf as the sequence header. This is primarily useful for haploidized VCFs e.g. chrW or mtDNA

#ref fasta, gene in bed format, sample list, and vcf file
fasta=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/beast_dating/allsamples/chr_W.fa
annotation=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/beast_dating/allsamples/Gene_Coordinates.bed
samples=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/beast_dating/allsamples/BEAST_AllSamples.list
vcf=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/beast_dating/allsamples/vcfs/chr_W.AC2.AllSites.vcf.gz

#loop through samples
sample=$1
echo ${sample}
bcftools consensus --haplotype A ${vcf} --sample ${sample} -f ${fasta} | sed "1 s/^>.*/>${sample}/" > fastas/chrW_${sample}.AC2.fa
```

Extract: `for i in $(cat BEAST_AllSamples.list); do sbatch -J EXTRACT_${i} 2.B_ExtractSamples_Singletons.sh ${i}; done` 

Merge:

```bash
cat chrW*.AC2.fa > chr_W.AC2.fa
seqkit subseq -r 1000000:1250000 chr_W.fa > chr_W_1-125Mb.fa
seqkit subseq -r 2000000:2100000 chr_W.fa > chr_W_2-21Mb.fa
seqkit subseq -r 2000000:2250000 chr_W.fa > chr_W_2-225Mb.fa
seqkit subseq -r 4000000:4250000 chr_W.fa > chr_W_4-425Mb.fa
seqkit subseq -r 5000000:5250000 chr_W.fa > chr_W_5-525Mb.fa
seqkit subseq -r 10000000:10250000 chr_W.fa > chr_W_10-10.25Mb.fa
seqkit subseq -r 14000000:14250000 chr_W.fa > chr_W_14-14.25Mb.fa

for i in $(ls *chr_W_*fa); do 
seqret -sequence ${i} -outseq ../nex/${i}.nex -osformat nexus
done 
```

For BEAUTI:

```bash
#for beauti
4 gamma categories, HKY, empirical sub rate 
clock rate strict at 5.05E-9
coalescent exponential, priors with lognormal CM / CP ancestor at M = 1.7  S=0.2
100M chains, log every 1k 
```





Plot:

```R
#### Plot BEAST annotated trees 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/beast_dating/allsamples/no_singletons')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(ggtree)
library(phytools)
library(ape)
library(treeio)
library(viridis)
library(ggpubr)
library(RColorBrewer)
library(tidyverse)

#Read in metadata
md = read_tsv('../../../Cuckoo_Full_Metadata_2023OCT3.txt')

files = list.files('.',paste0('.*ann'))
files = list.files('.',paste0('.*EXP.trees.ann'))
files = 'chr_W_5-525.trees.ann'
counter = 0
for (file in files){
  counter = counter +  1 
  iqtree = read.beast(file) 
  gg = ggtree(iqtree,layout='dendrogram') %<+% md
  
  #add label for 95% CIs
  lab = gsub('.trees.*','',file)
  heights = gg$data$height_0.95_HPD
  df = as.data.frame(do.call(rbind, heights)) #convert the list to a data frame
  df$node_value = 1:nrow(df) # Add node values as a new column
  colnames(df) = c("value1", "value2", "node")
  df = df[, c("node", "value1", "value2")]
  df = df %>% mutate(lab = paste0(round(value1,1),' - ',round(value2,1))) %>% select(!c(value1,value2))
  
  leg = md %>% select(Hap,HapCol,SpeciesShort) %>% unique
  gg$data = left_join(gg$data,df)
  ggp = gg  +
    geom_range(range='height_0.95_HPD', color='red', alpha=.6, size=1) +
    geom_tippoint(aes(fill = Hap,shape=SpeciesShort),size=1.5)+
    geom_nodelab(aes(label=lab),size=1) +
    ggtitle(lab)+
    #geom_tiplab(size=2)+
    #geom_nodelab(aes(x=branch, label=round(posterior, 2)), vjust=-.5, size=3) +
    scale_color_continuous(low="darkgreen", high="red") +
    scale_fill_manual(values=md$HapCol,breaks=md$Hap)+
    scale_shape_manual(values=c(21,4,8,2,2,2))+
    theme(legend.position=c(.1, .8))+
    geom_treescale(x = 10)+
    guides(fill=guide_legend(override.aes=list(shape=21)))+
    theme(legend.position='right')
  ggp
  assign(paste0('p',counter),ggp)
} 
p1
p2
pdf('../../../figures/BEAST_Divergence_Dating_All_2024FEB22.pdf',height=11,width=9)
ggarrange(p1,p2,p3,p4,p5,p6,common.legend = TRUE,nrow=3,ncol=2)
dev.off()

pdf('../../../figures/BEAST_Divergence_Dating_Sensitivity_2024FEB18.pdf',height=10,width=14)
ggarrange(p1,p2,p3,p4,common.legend = TRUE)
dev.off()

pdf('../../../figures/BEAST_Divergence_Dating_HighestESS_2024FEB18.pdf',height=3,width=8)
p1
dev.off()


#read in ML tree
ml = read.iqtree('chrW_AllORFS.fa.contree')
gg = ggtree(ml,layout='dendrogram') %<+% md
gg  +
  geom_tippoint(aes(fill = Hap,shape=SpeciesShort),size=3)+
  #geom_nodelab(aes(x=branch, label=round(posterior, 2)), vjust=-.5, size=3) +
  scale_fill_manual(values=md$HapCol,breaks=md$Hap)+
  scale_shape_manual(values=c(21,4,8))+
  guides(fill=guide_legend(override.aes=list(shape=21)))
```



Subset ORFs corresponding to those genes, then subset n=1 individual from each population, repeating this 4 times: 

```R
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=24:00:00

#mamba activate beast
#Run this in an empty directory. Specify the paths to the standard files and it will output a file for each gene (indicated by name in field 4 of the bed annotation), with the sample names from the vcf as the sequence header. This is primarily useful for haploidized VCFs e.g. chrW or mtDNA

#ref fasta, gene in bed format, sample list, and vcf file
fasta=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/beast_dating/n7_haplogroups/chr_W.fa
annotation=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/beast_dating/n7_haplogroups/Gene_Coordinates.bed
samples=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/beast_dating/n7_haplogroups/BEAST_AllSamples.list
vcf=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/beast_dating/n7_haplogroups/vcfs/chr_W.AllSites.vcf.gz

#loop through genes
for gene in $(cat Retained_Genes_Between10k-75k_BSP.txt); do

        echo ${gene}
        coords=$(grep -w ${gene} ${annotation} | awk '{print $1":"$2"-"$3}')

        #loop through samples
        for sample in $(cat ${samples}); do
        echo ${sample}
        samtools faidx ${fasta} ${coords} | bcftools consensus --haplotype A ${vcf} --sample ${sample} | sed "1 s/^>.*/>${sample}/" >> fastas/${gene}.fa
#
        done
done

#merge them into a megafile
seqkit concat fastas/*fa > Retained_Genes_Between10k-75k_BSP.fa

#get orfs from the megafile starting with ATG and stopping at a stop codon
getorf -methionine -minsize 300 -sequence Retained_Genes_Between10k-75k_BSP.fa -outseq ORFS.fa -table 0 -find 1

#grab the unique IDs, grab the one shared by 25 samples, convert it to stranded bed format  #switch 25 to the number of samples!
grep '>' ORFS.fa | sed 's/.*\[/[/g' | \
        sed 's/\].*/]/g' | sort | uniq -c | \
        awk '$1 == 25' | sed 's/\[/\t/g' | \
        sed 's/\]//g' | sed 's/ - /\t/g' | \
        awk '{OFS="\t"}{print "W", $2, $3, $4}' | \
        awk '{OFS="\t"}{if ($2 > $3) print $1, $3, $2,$3":"$2,$2-$3,"\t-"; else print $1, $2, $3,$2":"$3,$3-$2,"\t+"}' | \
        bedtools sort -i - > ORFS.bed

#since we have overlapping ORFS, just grab a random ORF when they overlap, we will analyze that
bedtools merge -i ORFS.bed -c 4 -o distinct | awk '{print $4}' | sed 's/,.*//g'  > ORFS_USE.txt

mkdir ORFS
#loop through the ORFs, if it is negative sense take the reverse complement. Otherwise just add the normal reading frame
for REGION in $(cat ORFS_USE.txt); do
        strand=$(grep -w ${REGION} ORFS.bed | awk '{print $6}')
        echo "WORKING ON REGION: ${REGION} WHICH IS STRAND: ${strand}"

        if [[ ${strand} = '-' ]]
                then
                seqkit subseq --seq-type dna --region ${REGION} Retained_Genes_Between10k-75k_BSP.fa | \
                        seqkit seq --reverse --complement --seq-type dna > ORFS/${REGION}.fa
                else
                seqkit subseq --seq-type dna --region ${REGION} Retained_Genes_Between10k-75k_BSP.fa > ORFS/${REGION}.fa
        fi

done

seqkit concat ORFS/*fa > chrW_AllORFS.fa

#And then just subset those individuals, rename them as 'W7','W1','CM',etc
for IT in {1..4}; do
seqtk subseq chrW_AllORFS.fa BEAST_Samples_IT${IT}.list > work/IT${IT}.fa
awk '{print $1, $2}' BEAST_Samples_IT${IT}.pop | sed -f <(awk '{print "s/>" $1 "/>" $2 "/"}' -) work/IT${IT}.fa > work/IT${IT}_nm.fa
seqret -sequence work/IT${IT}_nm.fa -outseq nex/Full_${IT}.nex -osformat nexus
done

```

For beauti:

```bash
#for beauti
codon divided, linked tree and clock
4 gamma categories, HKY, empirical sub rate 
clock rate strict at 5.05E-9 
calibrated yule, priors with lognormal CM / CP ancestor at M = 1.7  S=0.2
birthrate Y with gamma distribution alpha = 0.001 and beta 1000
100M chains, log every 1k 
```

After creating the xml in beauti, run:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=10
#SBATCH --time=200:00:00

RUN=$1

beast -threads 10 -overwrite -java ${RUN}.xml
```

Run treeannotator: mean heights, 10% burn-in, save as .boot somewhere. 

```R
#### Plot BEAST annotated trees 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/beast_dating/n7_haplogroups/mu_5.05e-09/')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(ggtree)
library(phytools)
library(ape)
library(treeio)
library(viridis)
library(ggpubr)
library(RColorBrewer)
library(tidyverse)

#Read in metadata
md = read_tsv('../../../Cuckoo_Full_Metadata_2023OCT3.txt')

files = list.files('.',paste0('.*ann'))
counter = 0
for (file in files){
  counter = counter +  1 
  iqtree = read.beast(file) 
  gg = ggtree(iqtree) %<+% md
  
  #add label for 95% CIs
  heights = gg$data$height_0.95_HPD
  df = as.data.frame(do.call(rbind, heights)) #convert the list to a data frame
  df$node_value = 1:nrow(df) # Add node values as a new column
  colnames(df) = c("value1", "value2", "node")
  df = df[, c("node", "value1", "value2")]
  df = df %>% mutate(lab = paste0(round(value1,2),' - ',round(value2,2))) %>% select(!c(value1,value2))
  
  leg = md %>% select(Hap,HapCol,SpeciesShort) %>% unique
  gg$data = gg$data %>% mutate(Species = ifelse(label == 'CM','CM',ifelse(label == 'CP','CP',
                                                                          ifelse(grepl('^W',label),'CC',NA))),
                               Hap = ifelse(label == 'CM','CM',ifelse(label == 'CP','CP',
                                                                      ifelse(grepl('^W',label),label,NA))))
  gg$data = left_join(gg$data,df)
  ggp = gg  +
    geom_range(range='height_0.95_HPD', color='red', alpha=.6, size=3) +
    geom_tippoint(aes(fill = Hap,shape=Species),size=3)+
    geom_nodelab(aes(label=lab), vjust=-1, size=1.5) +
    #geom_nodelab(aes(x=branch, label=round(posterior, 2)), vjust=-.5, size=3) +
    scale_color_continuous(low="darkgreen", high="red") +
    scale_fill_manual(values=md$HapCol,breaks=md$Hap)+
    scale_shape_manual(values=c(21,4,8))+
    theme(legend.position=c(.1, .8))+
    geom_treescale(x = 0.15)+
    guides(fill=guide_legend(override.aes=list(shape=21)))
  ggp
  assign(paste0('p',counter),ggp)
} 

pdf('../../../figures/BEAST_Divergence_Dating_2024JAN16.pdf',height=6,width=9)
ggarrange(p1,p2,p3,p4,common.legend = TRUE)
dev.off()
```

# Analyses: Population Genetic Differentiation

## Assign Geographic Clades 

Use k-means clustering to assign samples into discrete geographic 'populations' for analysis. 

```R
#### Determine relatedness with PHI statistic 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/relatedness')
.libPaths('~/mambaforge/envs/r/lib/R/library')
#Igraph approach
library(tidyverse)
library(RColorBrewer)
library(geosphere)
library(igraph)
library(spThin)
library(sf)
library(ggspatial)
library(factoextra)
library(ggpubr)

# Read metadata and filter for necessary columns
md = read_tsv('../Cuckoo_Full_Metadata_2023OCT3.txt') 
all = md %>% filter(Retained_Full == 1 & SpeciesShort == 'CC')

ds = all %>% mutate(INDEX = row_number())
thinned <- thin(ds,lat.col = 'Latitude',long.col = 'Longitude',spec.col = 'SpeciesShort',
                thin.par=25,reps=1,out.dir='.',out.base=sp,locs.thinned.list.return=T,write.files = F)
takeindex <- thinned %>% data.frame() %>% row.names() %>% as.vector()
thindat <- ds[grepl(paste0('^',takeindex,'$',collapse='|'),ds$INDEX),]


##### K-Means Clustering #####
sites <- st_as_sf(thindat, coords = c("Longitude", "Latitude"), 
                  crs = 4326, agr = "constant")

dist = st_distance(sites$geometry,  by_element = F)
d1 <- as.data.frame(dist)
d2 <- as.data.frame(lapply(d1, function(y) gsub(" [m]", "", y)))
names(d2) <- thindat$ID
row.names(d2) <- thindat$ID
d3 <- as.dist(as.matrix(d2))
cs1 <- fviz_nbclust(as.matrix(d3),kmeans,method = 'gap_stat',k.max = 25)
opt = cs1$data %>% ggplot(aes(group='Hi',x=clusters,y=gap,ymin=ymin,ymax=ymax))+
  geom_errorbar(col='cadetblue3')+ylab('Gap Statistic')+xlab('Number of Clusters')+
  geom_vline(xintercept=14,lty=2)+
  geom_line(col='cadetblue3')+
  geom_point(pch=16,size=3,col='cadetblue3')+
  theme_bw()
opt
pdf('../figures/OptimalKDistance_Clusters_2023OCT31.pdf',height=5,width=6)
opt
dev.off()

##### Clade Designations #####
#create new clade designation based on species and geographic distances 
world <- map_data("world")
set.seed(9999)

final = kmeans(d2, 14, nstart = 25)
d4 = thindat %>% mutate(Cluster = final$cluster)
dsite = st_as_sf(d4, coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant") 
ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),col='grey50',fill='grey95') +
  geom_sf(data = dsite, 
          aes(fill=as.factor(Cluster),shape=as.factor(Cluster)),
          size=5,show.legend = T) +
  scale_shape_manual(values=rep(c(21,24,25,22),8))+
  scale_fill_manual(values=c(brewer.pal(12,'Paired'),rev(brewer.pal(12,'Paired'))))+
  xlab('Longitude')+ylab('Latitude')+
  coord_sf(xlim = c(min(d4$Longitude)-5, max(d4$Longitude)+5), 
           ylim = c(min(d4$Latitude)-5, max(d4$Latitude)+5), expand = FALSE)+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"),
        legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))+
  annotation_scale(line_width=0.5)

#add names / colors
avgs = d4 %>% group_by(Cluster) %>% summarize(long = mean(Longitude)) %>% arrange(long) %>% 
  mutate(KDist = paste0('D',Cluster), KDCol = rep_len(c(viridis(12, option = 'turbo'),rev(viridis(4, option = 'turbo'))), 14), KDShape = rep_len(c(21, 24, 25, 22), 14)) %>% select(-long)

d5 = left_join(d4,avgs)
dsites = st_as_sf(d5, coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant") 
ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),col='grey50',fill='grey95') +
  geom_sf(data = dsites, 
          aes(fill=KDist,shape=KDist),
          size=5,show.legend = T) +
  xlab('Longitude')+ylab('Latitude')+
  scale_fill_manual(values=d5$KDCol,breaks=d5$KDist)+
  scale_shape_manual(values=d5$KDShape,breaks=d5$KDist)+
  coord_sf(xlim = c(min(d4$Longitude)-5, max(d4$Longitude)+5), 
           ylim = c(min(d4$Latitude)-5, max(d4$Latitude)+5), expand = FALSE)+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"),
        legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))+
  annotation_scale(line_width=0.5)

write.table(d5,'../figures/KDistance_Thinned_2023OCT31.txt',quote=F,sep='\t',row.names=F)

##### Re-assign thinned samples #####
#now, since we thinned the dataset, we need to go back to the full dataset and reassign clades based on the closest point
Bfull <- d5 %>% select(Longitude,Latitude,KDist,KDCol,KDShape,ID) %>% data.frame() 
Bfull <- Bfull %>% mutate(INDEX_B = row_number())
B <- Bfull %>% select(Longitude,Latitude)
ad <- all
Afull <- ad %>% mutate(INDEX_A = row_number())
A <- Afull %>% select(Longitude,Latitude)
for(i in 1:nrow(A)){
  #calucate distance against all of B
  distances<-geosphere::distGeo(A[i,], B)/1000
  #rank the calculated distances
  ranking<-rank(distances, ties.method = "first")
  
  #find the 3 shortest and store the indexes of B back in A
  A$shortest[i]<-which(ranking ==1) #Same as which.min()
  A$shorter[i]<-which(ranking==2)
  A$short[i]<-which(ranking ==3)
  
  #store the distances back in A
  A$shortestD[i]<-distances[A$shortest[i]] #Same as min()
  A$shorterD[i]<-distances[A$shorter[i]]
  A$shortD[i]<-distances[A$short[i]]
}
A <- A %>% mutate(INDEX_A = row_number())
A <- A %>% select(INDEX_A,shortest,Longitude,Latitude)
Aall <- merge(Afull,A,by=Reduce(intersect, list(names(Afull),names(A))))
Aall <- Aall %>% rename(INDEX_B = shortest)
Btake <- Bfull %>% select(INDEX_B,KDist,KDCol,KDShape)
as <- merge(Aall,Btake,by='INDEX_B')

#plot to confirm
as1 = as %>% mutate(LatJit = jitter(Latitude,amount = 1),
                    LonJit = jitter(Longitude,amount = 1)) 
ksites = st_as_sf(as1, coords = c("LonJit", "LatJit"), crs = 4326, agr = "constant") 
allsamp = ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),col='grey50',fill='grey95') +
  geom_sf(data = ksites, 
          aes(fill=KDist,shape=KDist),
          size=5,show.legend = T) +
  xlab('Longitude')+ylab('Latitude')+
  scale_fill_manual(values=d5$KDCol,breaks=d5$KDist)+
  scale_shape_manual(values=d5$KDShape,breaks=d5$KDist)+
  coord_sf(xlim = c(min(d4$Longitude)-5, max(d4$Longitude)+5), 
           ylim = c(min(d4$Latitude)-5, max(d4$Latitude)+5), expand = FALSE)+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"),
        legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))+
  annotation_scale(line_width=0.5)
allsamp

pdf('../figures/GeographicClusters_2023OCT31.pdf',height=8,width=12)
allsamp
dev.off()

write.table(as,'../figures/KDistance_All_2023OCT31.txt',quote=F,sep='\t',row.names=F)

```

## PCA

Plot PCA from plink on LD-pruned data:

```bash
setwd('~/merondun/cuculus_host/population_genetics/')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(viridis)
library(ggpubr)
library(vegan)
library(ggord)
library(recipes)
library(ggtext)
library(sf)
library(maps)
library(ggspatial)

md = read_tsv('~/merondun/cuculus_host/Metadata_Host.txt')
autos_canorus_LD-PCA.eigenvalues.txt
#read plink
prefix = 'autos_canorus_LD'
vec = read.table('autos_canorus_LD-PCA.eigenvectors.txt',header=F)
names(vec) = c('ID',paste0('PC',seq(1,10,1)))
vec = vec %>% mutate(ID = gsub('_F_.*','_F',ID),ID = gsub('_M_.*','_M',ID))
val = read.table('autos_canorus_LD-PCA.eigenvalues.txt',header=F)
val = val %>% mutate(VE = V1/sum(V1))

#merge with metadata, add a $distance vector based on lat/long PC1
vc = left_join(vec,md)
summary(prcomp(vc[, c("Longitude", "Latitude")], scale. = TRUE)) #66% of variance 
vc = vc %>% 
  mutate(Longitude = as.numeric(Longitude), Latitude = as.numeric(Latitude)) %>%
  { cbind(., pca = prcomp(.[, c("Longitude", "Latitude")], scale. = TRUE)$x[,1]) } %>%
  rename(Distance = pca) 

#superimpose
world <- map_data("world")

# Convert your data to a simple feature object
sites <- st_as_sf(vc, coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant")

#add same custom grey color scale 
kcols = data.frame(Kcols = gray.colors(5,start=0.05,end=0.95)[c(1,4,3,5,2)], AncestryK5 = paste0('K',seq(1,5,1)))

# spatil coordinates
geo_coords <- as.matrix(vc[, c("Longitude", "Latitude")])

# PCA coordinates
vc$PC1_flip = vc$PC1*-1
pca_coords <- as.matrix(vc[, c("PC1", "PC2")])

# Procrustes transformation
proc_transform <- vegan::procrustes(geo_coords, pca_coords,  symmetric = FALSE, scale = TRUE, reflection = 'best')

# Procrustes similarity statistic
procrustes_statistic <- proc_transform$ss / sum(geo_coords^2)

# Extract rotation angle in degrees
rotation_matrix <- proc_transform$rotation
rotation_angle <- atan2(rotation_matrix[2, 1], rotation_matrix[1, 1]) * (180 / pi)

# Permutation test
set.seed(123) # For reproducibility
perm_test <- protest(geo_coords, pca_coords, permutations = 100000)
cat("Procrustes similarity statistic:", procrustes_statistic, "\n")
cat("Rotation angle (degrees counterclockwise):", rotation_angle, "\n")
cat("Permutation test p-value:", perm_test$signif, "\n")

#Procrustes similarity statistic: 0.1368981 
#Rotation angle (degrees counterclockwise): -6.212215 
#Permutation test p-value: 9.9999e-06 

#extract frame 
new_pc = data.frame(nLong = proc_transform$X[,1], nLat = proc_transform$X[,2], nPC1 = proc_transform$Yrot[,1], nPC2 = proc_transform$Yrot[,2], ID = vc$ID, AncestryK5 = vc$AncestryK5, PC1 = vc$PC1, PC2 = vc$PC2, Lat = vc$Latitude, Long = vc$Longitude)

rp = new_pc %>% ggplot(aes(x=PC1,y=-PC2,fill=AncestryK5))+geom_point(pch=21)+scale_fill_manual(values=kcols$Kcols,breaks=kcols$AncestryK5)+
  xlab(paste0('PC1: ',round(val[1,2],3)*100,'%'))+  ylab(paste0('PC2: ',round(val[2,2],3)*100,'%'))+
  theme_bw(base_size=8)+theme(legend.position='top')
rp
new_pc %>% ggplot(aes(x=nPC1,y=nPC2,fill=AncestryK5))+geom_point(pch=21)+scale_fill_manual(values=kcols$Kcols,breaks=kcols$AncestryK5)+theme_bw(base_size=8)+ggtitle('Raw PCA')

#plot 
new_pc %>% ggplot()+
  geom_segment(aes(x=nLong,xend=nPC1,y=nLat,yend=nPC2),alpha=0.1)+
  geom_point(aes(x=nPC1,y=nPC2,fill=AncestryK5),pch=21)+
  scale_fill_manual(values=kcols$Kcols,breaks=kcols$AncestryK5)+
  theme_classic()

# Convert to simple feature object
library(scales)
no_rotation = vc %>% mutate(PC1_scaled = rescale(PC1, to = c(min(Longitude), max(Longitude))),
                            PC2_scaled = rescale(PC2*-1, to = c(min(Latitude), max(Latitude))))
sitesG <- st_as_sf(no_rotation, coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant")
sitesP <- st_as_sf(no_rotation, coords = c("PC1_scaled", "PC2_scaled"), crs = 4326, agr = "constant")

# Get world map data
world <- map_data("world")

# Plotting
pp = no_rotation %>% 
  ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), color = 'grey90', fill = 'white') +
  geom_segment(aes(x=PC1_scaled,xend=Longitude,y=PC2_scaled,yend=Latitude),alpha=0.1,col='darkred')+
  geom_point(aes(x=PC1_scaled,y=PC2_scaled,fill = AncestryK5), size = 1, pch = 21, show.legend = TRUE, color = 'grey20', alpha = 0.9) +
  geom_point(aes(x=Longitude,y=Latitude,fill = AncestryK5), size = 0.5, pch = 21, show.legend = TRUE, color = 'grey20', alpha = 0.2) +
  #geom_sf(data = sitesG, aes(fill = AncestryK5), size = 2, pch = 21, show.legend = TRUE, color = 'grey20', alpha = 0.2) +
  #geom_sf(data = sitesP, aes(fill = AncestryK5), size = 2, pch = 21, show.legend = TRUE, color = 'grey20', alpha = 0.9) +
  #geom_point(data = vc, aes(x = PC1_rot, y = PC2_rot), color = 'red', size = 2) +  # Add rotated PCA points
  scale_fill_manual(values = kcols$Kcols, breaks = kcols$AncestryK5) +  # Use unique Haplotypes for colors
  coord_cartesian(xlim = c(min(no_rotation$Longitude)-5, max(no_rotation$Longitude)+5),
                  ylim = c(min(no_rotation$Latitude)-5, max(no_rotation$Latitude)+5)) +
  theme_classic(base_size=8) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.background = element_rect(fill = "aliceblue")) +
  theme(legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.key.size = unit(0.1, 'cm'),legend.position='none') +
  annotate('text',x=-Inf,y=Inf,label=paste0('Procrustes SS: ',signif(procrustes_statistic,3),' angle = ',signif(rotation_angle,3), ' p = ',signif(perm_test$signif,3)),
           vjust=2,hjust=-0.1,size=2)+
  xlab(paste0('PC1 (',signif(val[1,2],3)*100,'%) / Longitude'))+
  ylab(paste0('PC2 (',signif(val[2,2],3)*100,'%) / Latitude'))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))

pdf('PCA_Superimposed-SpatialDistribution_2024FEB29.pdf',height=3,width=6)
ggarrange(rp,pp,widths=c(0.4,0.6))
dev.off()

#also slightly jitter points by 1degree for visualization 
# Convert your data to a simple feature object
sitesjit <- st_as_sf(vc %>% mutate(Longitude = jitter(Longitude,amount=1),
                                   Latitude = jitter(Latitude,amount=1)), 
                     coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant")

#add same custom grey color scale 
kcols = data.frame(Kcols = gray.colors(5,start=0.05,end=0.95)[c(1,4,3,5,2)], AncestryK5 = paste0('K',seq(1,5,1)))

# Plotting the map with PCA points and lines
kmap = ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),col='grey90',fill='white') + #add basemap
  geom_sf(data = sitesjit , 
          aes(fill=AncestryK5),
          size=2,pch=21,show.legend = T,col='grey20',alpha=0.9) +  #add the geometry for the points, I use pch 21-25 so the points have outlines
  scale_fill_manual(values=kcols$Kcols,breaks=kcols$AncestryK5)+ #custom fill encoded from metadata
  xlab('')+ylab('')+
  coord_sf(xlim = c(min(vc$Longitude)-5, max(vc$Longitude)+5), 
           ylim = c(min(vc$Latitude)-5, max(vc$Latitude)+5), expand = FALSE)+ #expand the map boundaries so that the labels show up, probably an easier way but I just add a gray rectangle in inkscape afterwards in between the labels and the map 
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"))+ #make the ocean blue
  theme(legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))+ 
  annotation_scale(line_width=0.5,location='bl')+  #add scale bar (requires ggspatial)
  guides(fill=guide_legend(override.aes=list(shape=21))) #ensure the legend has a fill-able shape
kmap
hapmap = ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),col='grey90',fill='white') + #add basemap
  geom_sf(data = sitesjit , 
          aes(fill=Hap),
          size=2,pch=21,show.legend = T,col='grey20',alpha=0.9) +  #add the geometry for the points, I use pch 21-25 so the points have outlines
  scale_fill_manual(values=md$HapCol,breaks=md$Hap)+ #custom fill encoded from metadata
  xlab('')+ylab('')+
  coord_sf(xlim = c(min(vc$Longitude)-5, max(vc$Longitude)+5), 
           ylim = c(min(vc$Latitude)-5, max(vc$Latitude)+5), expand = FALSE)+ #expand the map boundaries so that the labels show up, probably an easier way but I just add a gray rectangle in inkscape afterwards in between the labels and the map 
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"))+ #make the ocean blue
  theme(legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))+ 
  annotation_scale(line_width=0.5,location='bl')+  #add scale bar (requires ggspatial)
  guides(fill=guide_legend(override.aes=list(shape=21))) #ensure the legend has a fill-able shape
hapmap

pdf('SpatialDistribution_Haplogroups-AncestryK5_2024FEB29.pdf',height=7,width=7)
ggarrange(kmap,hapmap,nrow=2)
dev.off()


```



## ADMIXTURE

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=10
#SBATCH --time=200:00:00

K=$1

cd admixture

#Run Admixture
admixture -j7 --cv=5 ../merged_unrelfull/autos_canorus_LD.bed ${K} > autos_canorus_LD.log${K}.out

#Evaluate Runs
~/modules/evalAdmix/evalAdmix -plink ../merged_unrelfull/autos_canorus_LD -fname autos_canorus_LD.${K}.P -qname autos_canorus_LD.${K}.Q -P 10 -o eval_${K}
```

```bash
for i in {2..10}; do sbatch -J BAD_BOY_SERGIO_${i} Admixture_Eval.sh ${i}; done
```

CV Error: combine first `grep "CV" autos_canorus_LD*out | sed 's/.*log//g' | sed 's/.out:CV//g' > ~/merondun/cuculus_host/population_genetics/autos_canorus_LD-ADMIXTURE.CVs.txt` 

```R
setwd('~/merondun/cuculus_host/population_genetics/')
library(tidyverse)
library(ggplot2)
library(scales)
library(viridis)

cv <- read.table('autos_canorus_LD-ADMIXTURE.CVs.txt',header=FALSE)
names(cv) <- c('K','d1','d2','Error')
cvs = 
  cv %>% ggplot(aes(x=K,y=Error))+
  geom_line(show.legend = F,col='black')+
  geom_point(show.legend = F,col='black',size=2)+
  scale_color_manual(values=viridis(3))+
  scale_x_continuous(breaks=function(x) pretty(x,n=10))+
  ylab('C-V Error')+
  theme_classic()
cvs
pdf('~/merondun/cuculus_host/population_genetics/ADMIXTURE_CV-Error_2024FEB28.pdf',height=2,width=3)
cvs
dev.off()
```

Plot:

```R
### Plot ADMIXTURE
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/admixture')
.libPaths('~/mambaforge/envs/R/lib/R/library')
library(tidyverse)
library(viridis)
library(ggpubr)
library(RColorBrewer)
library(meRo) #devtools::install_github('merondun/meRo')

prefix = 'autos_canorus_LD' #1546013 SNPs 
qdir = '.' #directory with Q files

#read in admixture q's 
admix = melt_admixture(prefix = prefix, qdir = qdir)

#read in metadata
md = read_tsv('../Cuckoo_Full_Metadata_2023OCT3.txt')
admixmd = left_join(admix,md)

#visualize corelation of residuals between individuals 0 = good fit, if same population, positively correlated, different histories = negative corrleated
source('~/modules/evalAdmix/visFuns.R')

#this is only if you want to visualize correlation structure of a single run, recommend only for final one
Kval=5
pop <- admixmd %>% select(ID,KDist) %>% unique %>% select(KDist) # N length character vector with each individual population assignment
q <- as.matrix(read.table(paste0(prefix,'.',Kval,'.Q'))) # admixture porpotions q is optional for visualization but if used for ordering plot might look better
r <- as.matrix(read.table(paste0("eval_",Kval)))
plotAdmix(q=q, pop=pop$KDist)
plotCorRes(cor_mat = r, pop = pop$KDist, title = "Admixture evaluation as correlation of residuals", max_z=0.25, min_z=-0.25)

#loop through all admixture runs and extract the average correlation values, we want to MINIMIZE this! (closest to 0)
evaldat = NULL; for (Kval in seq(2,10,1)){
  r <- as.matrix(read.table(paste0("eval_",Kval)))
  mean_value <- mean(r,na.rm=TRUE)
  median_value <- median(r,na.rm=TRUE)
  sd_value <- sd(r,na.rm=TRUE)
  iqr_value <- IQR(r,na.rm=TRUE)
  valdat = data.frame(K = Kval,mean = mean_value,median=median_value,sd=sd_value,iqr=iqr_value)
  evaldat = rbind(valdat,evaldat)
}

#plot, for main figure show the n=3 lowest median
targs = evaldat %>% slice_min(median,n=3)
ep = evaldat %>% 
  ggplot(aes(x=K,y=median,ymin=median-iqr,ymax=median+iqr))+
  geom_rect(data=targs,aes(xmin=K-0.25,xmax=K+0.25,ymin=-Inf,ymax=Inf),fill='darkseagreen3')+
  geom_text(aes(y = 0.014, label = format(signif(median, 2), scientific = TRUE)),size=2) +  ylim(c(-0.015,0.015))+
  geom_hline(yintercept=0,lty=2)+
  geom_point()+ylab('Median +/- IQR Correlation of Residuals') +
  geom_errorbar()+
  theme_bw() + 
  scale_x_continuous(breaks = seq(min(evaldat$K), max(evaldat$K), by = 1)) +
  coord_flip()
ep
pdf('~/merondun/cuculus_host/population_genetics//ADMIXTURE_evalAdmix_MinimizeCorrelation_2024FEB28.pdf',height=5,width=4)
ep
dev.off()

#add a distance vector based on lat/long PC1
admixmd = admixmd %>% 
  mutate(Longitude = as.numeric(Longitude), Latitude = as.numeric(Latitude)) %>%
  { cbind(., pca = prcomp(.[, c("Longitude", "Latitude")], scale. = TRUE)$x[,1]) } %>%
  rename(Distance = pca) 

#if you want to add CV error directly on the label
cv = read.table('~/merondun/cuculus_host/population_genetics/autos_canorus_LD-ADMIXTURE.CVs.txt',header=FALSE)
names(cv) = c('Specified_K','d1','d2','Error')
cv = cv %>% mutate(label = paste0('K',Specified_K,' (',round(Error,2),')')) %>% select(!c(d1,d2,Error))

admixmd = admixmd %>% mutate(ID = fct_reorder(ID,desc(Distance))) 
admixmd = admixmd %>% mutate(Klab = paste0('K',Specified_K))
ad = left_join(admixmd,cv)
#save input
write.table(ad,file='~/merondun/cuculus_host/population_genetics/autos_canorus_LD-ADMIXTURE.Qmatrix.txt',quote=F,sep='\t',row.names=F)
ad = read_tsv('~/merondun/cuculus_host/population_genetics/autos_canorus_LD-ADMIXTURE.Qmatrix.txt')
k_order = ad %>% select(Klab,Specified_K) %>% arrange(Specified_K) %>% unique
ad$Hap = factor(ad$Hap,levels=c('W1','W2','W3','W4','W5','W6','W7'))
ad$Klab = factor(ad$Klab,levels=k_order$Klab)
ad = ad %>% group_by(Hap) %>% mutate(ID = reorder(factor(ID),Distance))
k_order = ad %>% ungroup %>% select(label,Specified_K) %>% arrange(Specified_K) %>% unique
ad$label = factor(ad$label,levels=k_order$label)
ad$KDist = factor(ad$KDist,levels=c('D6','D9','D12','D4','D2','D5','D1','D8','D13','D11','D14','D7','D10','D3'))

#add custom grey color scale 
kcols = data.frame(Kcols = gray.colors(5,start=0.05,end=0.95)[c(1,4,3,5,2)], Kcluster = paste0('K',seq(1,5,1)))

#Now plot
adplot = ad %>% 
  #filter(grepl('^2$|^3|^5', Specified_K)) %>%  # Specify the levels you want 
  #filter(Specified_K < 10) %>%  # Specify the levels you want 
  ggplot(aes(x = ID, y = Q, fill = factor(K))) +
  geom_col(color = NA, width=1) +
  scale_fill_manual(values = viridis(10, option = 'turbo')) +
  #scale_fill_manual(values = kcols$Kcols, breaks = kcols$Kcluster)+
  facet_grid(label ~ KDist, switch = "y", scales = "free", space = "free") +  # Switch facet labels to the left
  #facet_grid(label ~ Hap, switch = "y", scales = "free", space = "free") +  # Switch facet labels to the left
  theme_minimal() + labs(x = "", y = "Autosomal Ancestry Coefficient (Q)") +
  scale_y_continuous(expand = c(0, 0), n.breaks = 3, position = "right") +  # Move y-axis to the right
  scale_x_discrete(expand = expand_scale(add = 1)) +
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 7),
    panel.grid = element_blank(),
    strip.text.y.left = element_text(size = 8, angle = 0),
    strip.placement = "outside",  # Ensure facet labels are outside the plot area
    legend.position = 'bottom'
  ) 

adplot

pdf('~/merondun/cuculus_host/population_genetics/ADMIXTURE_K2K3K5-HAPS_2024FEB28.pdf',height=2.5,width=7.25)
adplot
dev.off()

pdf('~/merondun/cuculus_host/population_genetics/ADMIXTURE_K2-K10-HAPS_2024FEB28.pdf',height=8,width=7.25)
pdf('~/merondun/cuculus_host/population_genetics/ADMIXTURE_K2-K10-DIST_2024FEB28.pdf',height=8,width=7.25)
adplot
dev.off()

#assign individuals a K, and then plot PCA to see how those Ks fit 
#read PCA data
vec = read.table('~/merondun/cuculus_host/population_genetics/autos_canorus_LD-PCA.eigenvectors.txt',header=TRUE,comment.char = '')
vec = vec %>% dplyr::rename(ID = X.IID)
val = read.table('~/merondun/cuculus_host/population_genetics/autos_canorus_LD-PCA.eigenvalues.txt',header=FALSE,comment.char = '')
val = val %>% mutate(VE = paste0(round(V1/sum(V1),2)*100,'%'))

for (k in c(2,3,5)){ 
  #ks = ad %>% filter(Specified_K == k) %>% group_by(ID) %>% 
  #  slice_max(Q) %>% select(ID,K) %>% dplyr::rename(KCluster=K)
  
  ks = ad %>% 
    drop_na(HostParentShort) %>% 
    filter(Specified_K == k) %>%
    group_by(ID) %>%
    arrange(desc(Q)) %>%
    mutate(Rank = row_number()) %>%
    summarise(
      KCluster = case_when(
        max(Q) > 0.55 ~ first(K),
        TRUE ~ paste(K[1], K[2], sep = "_")
      )
    )
  cat('Number of samples intermediate for K',k,': ',ks %>% filter(grepl('_',KCluster)) %>% nrow,'\n')
  
  #merge with metadata, add a $distance vector based on lat/long PC1
  vc = left_join(ks,vec)
  p1 = vc %>% ggplot(aes(x=PC1,y=PC2,fill=KCluster))+
    geom_point(pch=21)+xlab(paste0('PC1 VE: ',val[1,2]))+ylab(paste0('PC2 VE: ',val[2,2]))+
    scale_fill_viridis(discrete=TRUE,option='turbo')+
    theme_bw()
  p2 = vc %>% ggplot(aes(x=PC3,y=PC4,fill=KCluster))+
    geom_point(pch=21)+xlab(paste0('PC3 VE: ',val[3,2]))+ylab(paste0('PC4 VE: ',val[4,2]))+
    scale_fill_viridis(discrete=TRUE,option='turbo')+
    theme_bw()
  ps = ggarrange(p1,p2,common.legend = TRUE,nrow=1,ncol=2)
  assign(paste0('k',k),ps)
}
pdf('~/merondun/cuculus_host/population_genetics/PCA_AssignIndividualsCluster_2024FEB28.pdf',height=7,width=6)
ggarrange(k2,k3,k5,nrow=3,ncol=1)
dev.off()

k5

#simply assign based on K5 
ks = ad %>% filter(Specified_K == 5) %>% group_by(ID) %>% slice_max(Q) %>% dplyr::rename(KCluster=K) %>% select(ID,Hap,KDist,Gens,Habitat,KCluster) %>% drop_na(Gens)
write.table(ks,file='~/merondun/cuculus_host/population_genetics/output_Admixture_AssignedAncestry_K5_Gens_2024FEB28.txt',quote=F,sep='\t',row.names=F)

#assign all samples 
ks = ad %>% filter(Specified_K == 5) %>% group_by(ID) %>% slice_max(Q) %>% dplyr::rename(KCluster=K) %>% select(ID,Hap,KDist,Gens,Habitat,KCluster)
write.table(ks,file='~/merondun/cuculus_host/population_genetics/output_Admixture_AssignedAncestry_K5_AllSamples_2024FEB28.txt',quote=F,sep='\t',row.names=F)
```

## Plot Tesselation + mtDNA Pies

```R
### Plot ADMIXTURE across landscape (tesselation)
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/admixture')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(terra)
library(tidyverse)
library(viridis)
library(stringr)
library(RColorBrewer)
library(meRo) #devtools::install_github('merondun/meRo')
library(LEA) #install with bioconductor, you don't actually need this if you impot your own q-matrix
library(tess3r) #install with github devtools
library(rworldmap) #for ggplot mapping 
library(sf) #for spatial plotting of distance vectors to confirm
library(ggspatial) #to add scale bars onto maps
library(ggpubr)

md = read_tsv('~/merondun/cuculus_host/Metadata_Host.txt')
admixmd = read_tsv('~/merondun/cuculus_host/population_genetics/autos_canorus_LD-ADMIXTURE.Qmatrix.txt')

### Tesselation
for (k in c(2,3,5)){
  cat('Tesselating K : ',k,'\n')
  #load in Q matrix 
  show_q = read.table(paste0('~/merondun/cuculus_host/population_genetics/admixture_q_files/autos_canorus_LD.',k,'.Q')) #read in the specific file 
  show_q_mat = as.matrix(show_q) #convert it to a matrix
  class(show_q_mat) = c('tess3Q','matrix','array') #make sure tess3r thinks that it's actually a tess object
  #grab coordinates 
  coords = admixmd %>% select(ID,Longitude,Latitude) %>% unique %>% select(Longitude,Latitude) #convert lat and long 
  coords_mat = as.matrix(coords) #convert coordinates to matrix
  #get map 
  map.polygon <- getMap(resolution = "low")
  pl = ggtess3Q(show_q_mat, coords_mat, map.polygon = map.polygon,col.palette = brewer.pal(k,'RdYlBu'))
  
  #plot 
  kp = pl +
    geom_path(data = map.polygon, aes(x = long, y = lat, group = group)) +
    xlim(min(coords$Longitude)-5,max(coords$Longitude)+5) +
    ylim(min(coords$Latitude)-5,max(coords$Latitude)+5) +
    coord_equal() +
    geom_point(data = coords, aes(x = Longitude, y = Latitude), size = 2,pch=21,fill='white') +
    xlab("Longitude") + ylab("Latitude") + ggtitle(paste0('K = ',k))+
    theme_classic()
  assign(paste0('k',k),kp)
}

pdf('~/merondun/cuculus_host/population_genetics/TESSELATION_K2-K5_2024FEB28.pdf',height=10,width=6)
png('~/merondun/cuculus_host/population_genetics/TESSELATION_K2-K5_2024FEB28.png',units='in',res=600,height=10,width=6)
ggarrange(k2,k3,k5,nrow=3)
dev.off()

#create  color/shape groups...
md = md %>% 
  filter(Analysis_PopulationGenetics == 1) %>% 
  mutate(HostParentShort = ifelse(Species_Latin == 'C. poliocephalus','Outgroup',HostParentShort)) %>% 
  arrange(Hap,HapCol)

#count the proportion of each haplogroup within each distance group 
mdc = md %>% count(KDist,Hap) %>% ungroup %>% group_by(KDist) %>% 
  mutate(Total = sum(n),
         Proportion = n/Total,
         Percent = paste0(round(n/Total,3)*100,'% (',n,')'))
#join that with lat/long 
mdcc = left_join(mdc,md %>% group_by(KDist) %>% summarize(Latitude=mean(Latitude),Longitude=mean(Longitude))) %>% 
  left_join(. , md %>% select(Hap,HapCol) %>% unique) %>% 
  arrange(KDist,Hap) %>% filter(KDist != 'CP')

#we ALSO need to create a scaling factor, based on how many haps are in each $Distance region 
scal = mdcc %>%
  select(KDist, Total) %>%
  unique() %>%
  ungroup() %>%
  mutate(Min = min(Total),
         Max = max(Total),
         Scaling_factor = ((Total - Min) / (Max - Min) * 10) + 2)
#make sure the scaling factor is linear
scal %>% ggplot(aes(x=Total,y=Scaling_factor))+geom_point()+theme_bw()
#add that scaling factor back to the haplogroups 
mdcc = left_join(mdcc,scal %>% select(KDist,Scaling_factor))

#add custom grey color scale 
kcols = data.frame(Kcols = gray.colors(5,start=0.05,end=0.95)[c(1,4,3,5,2)], Kcluster = paste0('K',seq(1,5,1)))

##### Plot pies across the world 
plot_pie <- function(data) {
  ggplot(data, aes(x = "", y = n, fill = Hap,label=paste0(KDist,': ',Total))) +
    geom_bar(col='black',lwd=0.5,width = 1, stat = "identity") +
    #geom_label(aes(x=Inf,y=Inf),fill='white',vjust=1.5,size=2)+
    #geom_text(size=2,vjust=-1)+
    coord_polar("y") +
    scale_fill_manual(values=md$HapCol,breaks=md$Hap)+
    theme_void() +
    theme(legend.position = "none")
}

#set up map and convert df to coordinate frame
world = map_data("world")
sites = st_as_sf(mdcc, coords = c("Longitude", "Latitude"), 
                 crs = 4326, agr = "constant") 

# Main map plot
#display.brewer.all(5,colorblindFriendly = TRUE)
p = 
  #ggplot()+ #for showing labels only 
  #geom_polygon(data = world, aes(x = long, y = lat, group = group),col='grey90',fill='white') +
  ggtess3Q(show_q_mat, coords_mat, map.polygon = map.polygon,col.palette = kcols$Kcols) + 
  geom_path(data = map.polygon, aes(x = long, y = lat, group = group)) +
  geom_sf(data = sites, aes(geometry = geometry), size = 0.1, alpha = 0.1) +
  coord_sf() +
  xlab('')+ylab('')+
  coord_sf(xlim = c(min(mdcc$Longitude)-5, max(mdcc$Longitude)+5), 
           ylim = c(min(mdcc$Latitude)-5, max(mdcc$Latitude)+5), expand = FALSE)+
  theme_classic(base_size = 8)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"))+
  theme(legend.position = 'top',legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))+
  annotation_scale(line_width=0.5,location='br')
#p

# Add pies
for (i in unique(mdcc$KDist)) {
  subset_data = mdcc %>% filter(KDist == i)
  lon = unique(subset_data$Longitude)
  lat = unique(subset_data$Latitude)
  scale_factor = unique(subset_data$Scaling_factor)
  cat('Scaling factor is : ',scale_factor,'\n')
  pie = plot_pie(subset_data)
  p <- p + annotation_custom(ggplotGrob(pie), 
                             xmin = lon - scale_factor, 
                             xmax = lon + scale_factor, 
                             ymin = lat - scale_factor, 
                             ymax = lat + scale_factor)
}

#p

#pdf('~/merondun/cuculus_host/population_genetics/TESSELATION_mtDNAHapPies_2024FEB28.pdf',height=6,width=9)
png('~/merondun/cuculus_host/population_genetics/TESSELATION_mtDNAHapPies_2024FEB28.png',units='in',res=600,height=3,width=5)
p
dev.off()

#re-run with labels to add the KDistance labels with sample size 
pdf('../figures/ADMIXTURE_TESSELATION_mtDNAHaps_2024FEB19_labs.pdf',height=3,width=5)
p
dev.off()

```

## Distance Correlations 

```R
#### Calculate pairwise geographic distance between distance groups 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/fst/2024feb_correlations/')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(viridis)
library(geosphere)
library(sf)
library(ggspatial)
library(RColorBrewer)

#Read in metadata
md = read_tsv('../../Cuckoo_Full_Metadata_2023OCT3.txt')

#groups for fst: distance ~ mtDNA
analyze_these = md %>% 
  filter(Analysis_PopulationGenetics == 1 & SpeciesShort == 'CC' & Sex == 'F') %>% select(KDist) %>%
  group_by(KDist) %>% mutate(DistanceHaps = n()) %>% unique %>% filter(DistanceHaps >= 3) %>% pull(KDist)

#calculate geographic distance between the groups
dists = md %>% filter(Analysis_PopulationGenetics == 1 & SpeciesShort == 'CC' & KDist %in% analyze_these) %>% group_by(KDist) %>% summarize(meanLat = mean(Latitude),meanLong = mean(Longitude))

#plot 
world = map_data("world")
sites = st_as_sf(dists, coords = c("meanLong", "meanLat"), 
                 crs = 4326, agr = "constant") 
fst_compars = ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),col='grey90',fill='white') +
  geom_sf(data = sites, 
          aes(fill=KDist,shape=KDist),
          size=4,alpha=0.9,show.legend = T,stroke=0.5) +
  scale_shape_manual(values=md$KDShape,breaks=md$KDist)+
  scale_fill_manual(values=md$KDCol,breaks=md$KDist)+
  xlab('')+ylab('')+
  coord_sf(xlim = c(min(dists$meanLong)-5, max(dists$meanLong)+5), 
           ylim = c(min(dists$meanLat)-5, max(dists$meanLat)+5), expand = FALSE)+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"))+
  theme(legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))+
  annotation_scale(line_width=0.5)

pdf('~/merondun/cuculus_host/correlations_geography_genetics/Geographic_Clusters_FST_Comparisons_2024FEB28.pdf',height=6,width=9)
fst_compars
dev.off()

# initialize an empty data frame to store pairwise distances
pairwise_distances <- data.frame()
sub_data = dists

# calculate pairwise distances for the current 'W' KDist
for(i in 1:(nrow(sub_data) - 1)) {
  for(j in (i + 1):nrow(sub_data)) {
    
    point1 <- c(sub_data$meanLong[i], sub_data$meanLat[i])
    point2 <- c(sub_data$meanLong[j], sub_data$meanLat[j])
    
    distance_km <- distHaversine(point1, point2) / 1000  # convert to km
    
    # append the result to the pairwise_distances data frame
    pairwise_distances <- rbind(pairwise_distances, 
                                data.frame(P1 = sub_data$KDist[i], 
                                           P2 = sub_data$KDist[j],
                                           Distance_km = distance_km))
  }
}

pairwise_distances = pairwise_distances %>% mutate(Group = paste0(P1,'__',P2))

# show the calculated pairwise distances
write_tsv(pairwise_distances,file='~/merondun/cuculus_host/correlations_geography_genetics/Pairwise_GeographicDistance_Km_2024FEB28.txt')

#write out populations
for (pop in unique(analyze_these)) {
  su = md %>% filter(Retained_Full_Unrelated & KDist == pop)
  write.table(su$ID,file=paste0('populations/',pop,'.pop'),quote=F,sep='\t',row.names=F,col.names=F)
}

```

Comparisons:

```bash
awk '{print $4}' ~/merondun/cuculus_host/correlations_geography_genetics/Pairwise_GeographicDistance_Km_2024FEB28.txt | sed '1d' > PairwiseComparisons.list
cat PairwiseComparisons.list 
D1__D10
D1__D12
D1__D13
D1__D14
D1__D3
D1__D4
D1__D5
D1__D7
D1__D9
D10__D12
```

Calculate FST:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=2
#SBATCH --time=48:00:00

# mamba activate snps
# for i in $(cat PairwiseComparisons.list); do sbatch -J FST_${i} ~/merondun/cuculus_host/correlations_geography_genetics/2.Pairwise_Distance_FST.sh ${i}; done 
GROUP=$1

for CHR in $(cat Chromosomes.list); do

echo "WORKING ON CHR: ${GROUP} and ${CHR}"

p1=$(echo ${GROUP} | sed 's/__.*//g')
p2=$(echo ${GROUP} | sed 's/.*__//g')

mkdir work out

if [[ $CHR == 'chr_Z' ]]
then

    calculate fst
    ~/modules/vcftools/bin/vcftools --gzvcf ../../snp_unrelfull/${CHR}.IF-GF-MM1.vcf.gz --fst-window-size 100000 --fst-window-step 100000 --weir-fst-pop populations/${p1}.pop --weir-fst-pop populations/${p2}.pop --out work/${CHR}_${GROUP}
    awk -v g=${GROUP} '{OFS="\t"}{print $1, $2, $3,$4,$5, g}' work/${CHR}_${GROUP}.windowed.weir.fst > out/${CHR}_${GROUP}.fst

    #calculate also for only males 
    ~/modules/vcftools/bin/vcftools --gzvcf ../../snp_unrelfull/${CHR}.M_IF-GF-MM1.vcf.gz --fst-window-size 100000 --fst-window-step 100000 --weir-fst-pop populations/${p1}.pop --weir-fst-pop populations/${p2}.pop --out work/${CHR}_${GROUP}.M
    awk -v g=${GROUP} '{OFS="\t"}{print $1, $2, $3,$4,$5, g"__M"}' work/${CHR}_${GROUP}.M.windowed.weir.fst > out/${CHR}_${GROUP}.M.fst

else

    calculate fst
    ~/modules/vcftools/bin/vcftools --gzvcf ../../snp_unrelfull/${CHR}.IF-GF-MM1.vcf.gz --fst-window-size 100000 --fst-window-step 100000 --weir-fst-pop populations/${p1}.pop --weir-fst-pop populations/${p2}.pop --out work/${CHR}_${GROUP}
    awk -v g=${GROUP} '{OFS="\t"}{print $1, $2, $3,$4,$5, g}' work/${CHR}_${GROUP}.windowed.weir.fst > out/${CHR}_${GROUP}.fst

fi

done
```

For mtDNA, calculate PHI using an AMOVA, subset pops:

```bash
for i in $(cat Populations.list); do samtools faidx chr_MT_All.SNP.DP3-AC2-MQ40.min4.fasta $(cat populations/${i}.pop) > fastas/${i}.fa; done
```

And repeat only for W:

```bash
for i in $(cat Populations.list); do samtools faidx chr_W_All.SNP.DP3-AC2-MQ40.min4.fasta $(grep '_F$' populations/${i}.pop) > fastas/${i}.W.fa; done
```

And then in R to calculate phi for mtDNA/W:

```R
library(tidyverse)
library(ape)
library(pegas)
library(adegenet)
library(mmod)
library(poppr)

comps = read_tsv('PairwiseComparisons.list',col_names = F)

#save results here 
fdM = NULL

#for chr_MT
for (comp in comps$X1) {
  cat('Working on comparison: ',comp,'\n')
  
  #p1 and p2
  p1 = gsub('__.*','',comp)
  p2 = gsub('.*__','',comp)
  
  seqs1 = read.dna(paste0('fastas/',p1,'.fa'), format = "fasta")
  seqs2 = read.dna(paste0('fastas/',p2,'.fa'), format = "fasta")
  
  combined_seqs = rbind(seqs1,seqs2)
  pops = c(rep(p1,nrow(seqs1)),rep(p2,nrow(seqs2)))
  
  #import to genind, then convert to genclone 
  binary_matrix <- as.matrix(combined_seqs)
  genind_obj <- DNAbin2genind(combined_seqs, pop = pops)
  
  gc = as.genclone(genind_obj)

  #add strata (populations)
  stratadf = data.frame(ind = indNames(genind_obj), stratum = pops)
  strata(gc) = stratadf
  res = poppr.amova(gc,hier=~stratum,clonecorrect = FALSE)
  resd = data.frame(Group = paste0(p1,'__',p2),PHI = res$statphi[[1]])
  
  fdM = rbind(fdM,resd)
}

###for W
fd = NULL
for (comp in comps$X1) {
  cat('Working on comparison: ',comp,'\n')
  
  p1 = gsub('__.*','',comp)
  p2 = gsub('.*__','',comp)
  
  seqs1 = read.dna(paste0('fastas/',p1,'.W.fa'), format = "fasta")
  seqs2 = read.dna(paste0('fastas/',p2,'.W.fa'), format = "fasta")
  
  combined_seqs = rbind(seqs1,seqs2)
  pops = c(rep(p1,nrow(seqs1)),rep(p2,nrow(seqs2)))
  
  binary_matrix <- as.matrix(combined_seqs)
  genind_obj <- DNAbin2genind(combined_seqs, pop = pops)

  gc = as.genclone(genind_obj)
  stratadf = data.frame(ind = indNames(genind_obj), stratum = pops)
  strata(gc) = stratadf
  res = poppr.amova(gc,hier=~stratum,clonecorrect = FALSE)
  resd = data.frame(Group = paste0(p1,'__',p2),PHI = res$statphi[[1]])
  
  fd = rbind(fd,resd)
}

dats = cbind(fd %>% dplyr::rename(PHI_W = PHI),fdM %>% 
        select(-Group)) 
#dats %>% ggplot(aes(x=PHI_W,y=PHI))+
  # geom_smooth(method='lm')+
  # geom_point()+
  # theme_bw()

write_tsv(dats,file='PHI_chr_MT-W_2024FEB28.txt')
```

Plot final correlations:

```R
#### Plot Distance ~ FST / PHIst
setwd('~/merondun/cuculus_host/correlations_geography_genetics/')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(viridis)
library(geosphere)
library(broom)
library(ggpubr)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(meRo)
library(RColorBrewer)

#Read in metadata
md = read_tsv('~/merondun/cuculus_host/Metadata_Host.txt')

d = read_tsv('Pairwise_GeographicDistance_Km_2024FEB28.txt')
f = read_tsv('Pairwise_FST_DistanceGroups_Autosomes_2024FEB28.txt')
p = read_tsv('Pairwise_PHIst_DistanceGroups_W-MT_2024FEB28.txt')

#correlation between W/MT
p %>% ggplot(aes(x=PHI_W,y=PHI))+geom_point()+theme_bw()
cor.test(p$PHI_W,p$PHI)
# 
# Pearson's product-moment correlation
# 
# data:  p$PHI_W and p$PHI
# t = 8.8912, df = 43, p-value = 2.695e-11
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.6693965 0.8884474
# sample estimates:
#       cor 
# 0.8047956 

#prep and merge frames 
p = p %>% dplyr::rename(chr_W = PHI_W, chr_MT = PHI) %>% pivot_longer(!Group) %>% 
  dplyr::rename(FST = value, chr = name) %>% mutate(FST = pmax(0, pmin(1, FST)))
names(f) = c('chr','start','end','snps','FST','Group')
z_comparison = f %>% filter(chr == 'chr_Z') %>% mutate(Subset = ifelse(grepl('__M$',Group),'Males_Only','All_Samples'),Group = gsub('__M$','',Group)) %>% 
  select(-snps) %>% pivot_wider(names_from = Subset,values_from=FST) %>% na.omit 
z_comparison %>% ggplot(aes(x=All_Samples,y=Males_Only))+geom_point()+theme_bw()
cor.test(z_comparison$All_Samples,z_comparison$Males_Only)
# Pearson's product-moment correlation
# 
# data:  z_comparison$All_Samples and z_comparison$Males_Only
# t = 356.72, df = 34922, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.8835307 0.8880480
# sample estimates:
#       cor 
# 0.8858103 

zc = f %>% filter(grepl('__M$',Group)) #only grab males for Z 
amc = f %>% filter(chr != 'chr_Z') #exclude the male only comparison for autosomes 
f2 = rbind(zc,amc) %>% mutate(FST = pmax(0, pmin(1, FST))) %>% select(!c(start,end,snps)) %>% mutate(Group = gsub('__M','',Group))
pf = rbind(p,f2)

#merge with geographic distance 
df = left_join(pf,d) 

#assign AvZ 
df = df %>% mutate(AvZ = ifelse(chr == 'chr_Z','Z',ifelse(chr == 'chr_MT','mtDNA',ifelse(chr == 'chr_W','W','Autosome'))))

dfs = df %>% group_by(AvZ,Group,P1,P2,Distance_km) %>% 
  sum_stats(FST)
dfs %>% ggplot(aes(x=log(Distance_km),y=log(mean),col=AvZ))+
  geom_point()+
  geom_smooth()+
  theme_bw()

# Spearman's rho and cor test within each AvZ level and gather coefficients
cors = dfs %>% group_by(AvZ) %>%
  summarize(
    rho = cor(mean, Distance_km,method='spearman'),
    p_value = cor.test(mean, Distance_km,method='spearman')$p.value) %>% 
  mutate(padj = p.adjust(p_value,method='bonferroni'),
         signif = case_when(
           padj < 0.05 ~ "*",
           TRUE ~ " (n.s.)"),
         label = paste0('rs: ',round(rho, 2), signif))

cols = brewer.pal(4,'Dark2')
dfs$AvZ = factor(dfs$AvZ,levels=c('Autosome','Z','W','mtDNA'))

#lmm, account for P1 and P2 
model_summaries <- dfs %>%
  group_by(AvZ) %>%
  do({
    model <- lmer(log(pmax(mean, 0.005)) ~ log(pmax(Distance_km, 0.005)) + (1|P1) + (1|P2), data = .)
    summary_df <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE)
  }) %>%
  dplyr::bind_rows()
#output fixed effect of distance, bonferroni correction 
model_summaries %>% filter(grepl('Distance',term)) %>% select(-effect,-term) %>% ungroup %>% 
  mutate(padj = p.adjust(p.value,method='bonferroni'),
         signif = ifelse(padj < 0.05,'*','n.s.'))
# A tibble: 4 × 10
# AvZ      estimate std.error statistic    df  p.value conf.low conf.high     padj signif
# <fct>       <dbl>     <dbl>     <dbl> <dbl>    <dbl>    <dbl>     <dbl>    <dbl> <chr> 
#   1 Autosome    1.22     0.0614     19.9   30.3 6.01e-19   1.10       1.35  2.41e-18 *     
#   2 Z           1.09     0.0747     14.6   34.1 3.07e-16   0.941      1.24  1.23e-15 *     
#   3 W           0.614    0.331       1.86  41.4 7.05e- 2  -0.0538     1.28  2.82e- 1 n.s.  
# 4 mtDNA       0.499    0.183       2.73  42.3 9.30e- 3   0.130      0.868 3.72e- 2 *   
  

# Create the plot
pp1 = ggplot(dfs, aes(x = log(pmax(Distance_km,0.005)), y = log(pmax(mean,0.005)), color = AvZ,shape=AvZ)) +
  geom_point(size=1) + #0.25 for main plot 
  geom_smooth(method = "lm", se = TRUE,lwd=0.75) +  #0.5 for main plot 
  # geom_text(data = cors[1,], aes(x = Inf, y = -Inf, label = label), vjust = -6, hjust = 1.73, col = cols[1]) +
  # geom_text(data = cors[2,], aes(x = Inf, y = -Inf, label = label), vjust = -4.5, hjust = 1.73, col = cols[2]) +
  # geom_text(data = cors[3,], aes(x = Inf, y = -Inf, label = label), vjust = -3, hjust = 1.1, col = cols[3]) +
  # geom_text(data = cors[4,], aes(x = Inf, y = -Inf, label = label), vjust = -1.5, hjust = 2.02, col = cols[4]) +
  labs(title = "Relationship between Mean and Distance_km",
       x = "log(Geographic Distance)",
       y = "log(FST | ΦST)",
       color = "AvZ") +
  scale_color_brewer(palette='Dark2')+
  facet_wrap(.~AvZ,scales='free',nrow=4,ncol=1)+ #for dual plot sensitivity
  scale_shape_manual(values=c(15,16,17,3))+
  theme_bw(base_size=6) + theme(legend.position = 'none')
pp1

pdf('~/merondun/cuculus_host/correlations_geography_genetics/Correlations_logFST-logDIST_2024FEB29.pdf',height=1.5,width=1.5)
pp1
dev.off()

# Spearman's rho and cor test within each AvZ level and gather coefficients
dfs = dfs %>% mutate(onephi = mean / (1 - mean))
cors2 = dfs %>% group_by(AvZ) %>%
  summarize(
    rho = cor(onephi, Distance_km,method='spearman'),
    p_value = cor.test(onephi, Distance_km,method='spearman')$p.value) %>% 
  mutate(padj = p.adjust(p_value,method='bonferroni'),
         signif = case_when(
           padj < 0.05 ~ "*",
           TRUE ~ " (n.s.)"),
         label = paste0('rs: ',round(rho, 2), signif))

pp2 = ggplot(dfs, aes(x = Distance_km, y = pmin(onephi,1), color = AvZ,shape=AvZ)) +
  geom_point(size=1) +
  geom_smooth(method = "lm", se = TRUE,lwd=0.75) +
  # geom_text(data = cors[1,], aes(x = Inf, y = -Inf, label = label), vjust = -6, hjust = 1.73, col = cols[1]) + 
  # geom_text(data = cors[2,], aes(x = Inf, y = -Inf, label = label), vjust = -4.5, hjust = 1.73, col = cols[2]) + 
  # geom_text(data = cors[3,], aes(x = Inf, y = -Inf, label = label), vjust = -3, hjust = 1.1, col = cols[3]) + 
  # geom_text(data = cors[4,], aes(x = Inf, y = -Inf, label = label), vjust = -1.5, hjust = 2.02, col = cols[4]) + 
  labs(title = "Relationship between FST | ΦST and Distance",
       x = "Geographic Distance (km)",
       y = "FST | ΦST / 1 - FST | ΦST)",
       color = "AvZ") +
  scale_color_manual('Compartment',values=brewer.pal(4,'Dark2'))+
  scale_shape_manual('Compartment',values=c(16,17,15,18))+
  theme_bw(base_size=6)+
  facet_wrap(.~AvZ,scales='free',nrow=4,ncol=1)+
  theme(legend.position='top')
pp2
pdf('~/merondun/cuculus_host/correlations_geography_genetics/Correlations_1-logFST-logDIST_2024FEB29.pdf',height=1.5,width=1.5)
pp
dev.off()

write_tsv(cors,file='~/merondun/cuculus_host/correlations_geography_genetics/Correlations_logFST-logDIST_cortest-results.txt')

pdf('~/merondun/cuculus_host/correlations_geography_genetics/Correlations-BothMethods_2024FEB29.pdf',height=6,width=4)
ggarrange(pp1,pp2,common.legend = TRUE,nrow=1,ncol=2)
dev.off()

```

# Analyses: Host & Habitat Associations 

## Extract Land Class 

First, extract land class from the coordinates:

```R
#Extract land class 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/spatial')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(terra)
library(sf)
library(tidyverse)
library(viridis)

file_list = list.files(path = ".", pattern = "\\.hdf$", full.names = TRUE)

# Initialize an empty list to store SpatRasters
rasters = list()

for (file in file_list) {
  raster_layer = rast(file, lyrs = "LC_Type1")
  
  # Append to the list
  rasters[[length(rasters) + 1]] <- raster_layer
}

# Merge all rasters into one
merged_raster = do.call(merge, rasters)
#plot(merged_raster)

writeRaster(merged_raster, 'LandClass_2010_MODIS.tif',overwrite=TRUE)
merged_raster = rast('LandClass_2010_MODIS.tif')

#Extract points
md = read_tsv('~/merondun/cuculus_host/Metadata_Host.txt')
md = md %>% drop_na(Egg)

# Convert to an sf object
points_sf <- st_as_sf(md, coords = c("Longitude", "Latitude"), crs = 4326)
points_sf_transformed <- st_transform(points_sf, crs = crs(merged_raster))

# Extract the land class values at the given coordinates
land_class_values <- terra::extract(merged_raster, points_sf_transformed)
lc_legend = read_tsv('LC.legend',col_names = F)
names(lc_legend) = c('LC_Type1','Habitat')
lc = left_join(land_class_values,lc_legend)
lcp = lc %>% group_by(Habitat) %>% summarize(total = n()) %>% 
  ggplot(aes(x=Habitat,fill=Habitat,y=total,label=total))+
  geom_bar(stat='identity')+geom_text(vjust=-1)+theme_bw()
pdf('../figures/LandClass_ValuesCounts_SINUproj_2024MAR05.pdf',height=4,width=6)
lcp
dev.off()
lc2 = cbind(md %>% select(-Habitat),lc %>% select(-ID))

#save it 
write.table(lc2,file='../admixture/Admixture_KClusters_K5_LandClass_Input_2024MAR05.txt',quote=F,sep='\t',row.names=F)

#plot to confirm
# Calculate the extent of the points
xmin <- min(st_coordinates(points_sf_transformed)[, "X"]) - 20
xmax <- max(st_coordinates(points_sf_transformed)[, "X"]) + 20
ymin <- min(st_coordinates(points_sf_transformed)[, "Y"]) - 20
ymax <- max(st_coordinates(points_sf_transformed)[, "Y"]) + 20

# Create a SpatExtent object
new_extent <- ext(xmin, xmax, ymin, ymax)

# Clip the raster to this new extent
clipped_raster <- crop(merged_raster, new_extent)

#save the points
st_write(points_sf_transformed, "Cuckoo_Locations_sinu.shp",delete_layer = TRUE)

# Plot to verify
pdf('../figures/LandClass_Values_SINUproj_2024MAR05.pdf',height=5,width=9)
plot(clipped_raster)
plot(points_sf_transformed, add = TRUE, col = "red", pch = 20)
dev.off()

```



## dbRDA

```R
#### dbRDA for e.g. egg type ~ maternal haplogroup associations 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/mantel')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(viridis)
library(geosphere)
library(sf)
library(ggspatial)
library(RColorBrewer)
library(cluster)
library(vegan)
library(ape)
library(ggpubr)
library(data.table)
library(ecodist)

#calculate geographic distance matrix
mdf = read_tsv('~/merondun/cuculus_host/Metadata_Host.txt') 

#only grab samples with known egg (hash out filter for related individuals)
md_egg = mdf %>%
  filter(Analysis_PopulationGenetics == 1) %>%
  drop_na(Egg) %>%
  select(ID = ID, Host = HostParentShort, Egg, Habitat, Haplogroup = Hap, Ancestry = AncestryK5, Geography = KDist, Latitude, Longitude, CountryFull, HapCol)

##### Initialize, summary stats on raw data #####
#ensure they are all factor variables
md_egg = md_egg %>% mutate_at(c('Host','Egg','Habitat','Haplogroup','Ancestry','Geography'),as.factor)

#only retain Host where we have at least 2 cuckoos 
md_egg = md_egg %>% group_by(Host) %>% mutate(TotalHost = n()) %>% ungroup %>% group_by(Habitat) %>% mutate(TotalHabitat = n()) %>% ungroup %>% group_by(Egg) %>% mutate(TotalEgg = n())  %>% ungroup
minobs=2
md_egg = md_egg %>% filter(TotalHost >= minobs & TotalHabitat >= minobs & TotalEgg >= minobs) 

# If you want to exclude the blue clades W1, W2, W3! 
#md_egg <-  md_egg %>% filter(!grepl('W1|W2|W3',Haplogroup))
md <- md_egg

kept_order_samples = md$ID

#Calculate pairwise geographic distance 
md_sf = st_as_sf(md, coords = c("Longitude", "Latitude"))
st_crs(md_sf) = 4326 # Set the projection as ESPG 4326 (long_lat)
#Apply st_distance 
geo_m = as.data.frame(st_distance(md_sf))
geo_df = geo_m %>% mutate(across(everything(), ~ as.numeric(gsub(" \\[m\\]$", "", .))/1000))
names(geo_df) = md$ID
rownames(geo_df) = md$ID
geo_mat = as.matrix(geo_df)

#visualize the geographic distance matrix
geo_mds = as.data.frame(cmdscale(geo_mat,2))
geo_mds$ID = rownames(geo_mds)
geo_mds_input = left_join(geo_mds,md)
geo_p = ggplot(geo_mds_input,aes(x=V1,y=V2,color=CountryFull)) + 
  geom_point() + theme_bw(base_size=6) + xlab('MDS1')+ylab('MDS2')+
  scale_color_viridis(discrete=TRUE,option='turbo')+
  theme(legend.position = "top", # Moves legend to top
        legend.text = element_text(size = 6), # Adjusts text size
        legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))
geo_p

#load in autosomal distance, calculated: ~/modules/VCF2Dis-1.50/bin/VCF2Dis -InPut autos_canorus_LD.vcf.gz -OutPut autos_canorus_LD.pdist 
auto = read.table('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/merged_unrelfull/autos_canorus_LD.pdist',header=F)
#add proper names, because VCF2DIS truncates them
auto = auto %>% mutate(IDNum = gsub('_.*','',V1))
auto_id = left_join(auto,mdf %>% select(ID) %>% mutate(IDNum = gsub('_.*','',ID)))  
ord = auto_id$ID
autos = auto_id %>% select(!c(IDNum,ID,V1))
names(autos) = ord
rownames(autos) = ord

#extract columns in same order as geographic distance
common_names = intersect(rownames(autos), rownames(geo_mat))
auto_aligned = autos[common_names, common_names] #subset 'auto' matrix to match the 'geo_mat' matrix rows and columns
data.frame(mat1 = rownames(auto_aligned), mat2 = rownames(geo_mat)) %>% filter(mat1 != mat2) #sanity check for alignment 
auto_mat = as.matrix(auto_aligned)

#visualize autosomal distance matrix in terms of ancestry K = 5
kcols = data.frame(Kcols = gray.colors(5,start=0.05,end=0.95)[c(1,4,3,5,2)], Kcluster = paste0('K',seq(1,5,1)))
auto_mds = as.data.frame(cmdscale(auto_mat,2))
auto_mds$ID = rownames(auto_mds)
auto_mds_input = left_join(auto_mds,md)
auto_p = ggplot(auto_mds_input,aes(x=V1,y=V2,fill=Ancestry,shape=Ancestry)) + 
  geom_point(size=1.5) + theme_bw(base_size=6) + xlab('MDS1')+ylab('MDS2')+
  scale_fill_manual('K',values=kcols$Kcols,breaks=kcols$Kcluster)+
  scale_shape_manual(values=c(21,22,23,24,25))+
  guides(fill=guide_legend(override.aes=list(shape=21)))+
  theme(legend.position = "top", # Moves legend to top
        legend.text = element_text(size = 6), # Adjusts text size
        legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))
auto_p

#calculate mtDNA distance
seq = read.dna('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/beast_dating/variant_only/ml_trees/chr_MT_All.SNP.DP3-AC1-MQ40.min4.fasta',format='fasta')
dna_dist = as.data.frame(as.matrix(dist.dna(seq,model='JC69')))

#extract columns in same order as geographic distance
common_names2 = intersect(rownames(dna_dist), rownames(geo_mat))
dna_aligned = dna_dist[common_names2, common_names2] #subset 'auto' matrix to match the 'geo_mat' matrix rows and columns
data.frame(mat1 = rownames(dna_aligned), mat2 = rownames(geo_mat)) %>% filter(mat1 != mat2) #sanity check for alignment 
dna_mat = as.matrix(dna_aligned)

#visualize the mtDNA distance matrix in terms of haplogroup 
dna_mds = as.data.frame(cmdscale(dna_mat,2))
dna_mds$ID = rownames(dna_mds)
dna_mds_input = left_join(dna_mds,md)
dna_p = ggplot(dna_mds_input,aes(x=V1,y=V2,color=Haplogroup)) + 
  geom_point() + theme_bw(base_size=6) + xlab('MDS1')+ylab('MDS2')+
  scale_color_manual(values=md$HapCol,breaks=md$Haplogroup)+
  theme(legend.position = "top", # Moves legend to top
        legend.text = element_text(size = 6), # Adjusts text size
        legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))
dna_p

pdf('../figures/20240828_dbRDA_DistanceMatrix_Inputs_n80.pdf',height=5,width=2.5)
ggarrange(geo_p,auto_p,dna_p,nrow=3)
dev.off()

geoscat = as.data.frame(geo_mat) %>% mutate(ID = rownames(.)) %>% pivot_longer(!ID,names_to = 'ID2',values_to = 'EDistance') %>% mutate(Explanatory = 'Biogeographical')
autoscat = as.data.frame(auto_mat) %>% mutate(ID = rownames(.)) %>% pivot_longer(!ID,names_to = 'ID2',values_to = 'EDistance') %>% mutate(Explanatory = 'Autosomal')
mtscat = as.data.frame(dna_mat) %>% mutate(ID = rownames(.)) %>% pivot_longer(!ID,names_to = 'ID2',values_to = 'EDistance') %>% mutate(Explanatory = 'Matrilineal')

##### dbRDA: reverse #####
covars = md %>% select(ID,Egg,Host,Habitat) %>% mutate_all(as.factor)
inputs = c('dna_mat','auto_mat','geo_mat')
dbr = list()
for (inp in inputs) {
  
  #constrained ordination with step selection 
  null_formula_str = as.formula(paste(inp, "~ 1"))
  m1f = as.formula(paste(inp, "~ Egg"))
  m2f = as.formula(paste(inp, "~ Host"))
  m3f = as.formula(paste(inp, "~ Habitat"))
  mff= as.formula(paste(inp, "~ Egg + Host + Habitat"))
  null = dbrda(null_formula_str, covars, dist="gow",scaling=TRUE)  # Model with intercept only
  m1 = dbrda(m1f, md, dist="gow",scaling=TRUE) # Model with all explanatory variables
  m2 = dbrda(m2f, md, dist="gow",scaling=TRUE) # Model with all explanatory variables
  m3 = dbrda(m3f, md, dist="gow",scaling=TRUE) # Model with all explanatory variables
  mf = dbrda(mff, md, dist="gow",scaling=TRUE) # Model with all explanatory variables
  
  ## With scope present
  step = ordistep(null, scope = formula(mff), perm.max = 200)
  stepres = as.data.frame(step$anova)
  stepres = stepres %>% mutate(Response = gsub('\\+ ','',rownames(.))) %>% select(Response,stepF = F)
  
  #by terms
  r1 = as.data.frame(anova(m1, by="terms", permu=10000)) # test for sign. environ. variables
  r2 = as.data.frame(anova(m2, by="terms", permu=10000)) # test for sign. environ. variables
  r3 = as.data.frame(anova(m3, by="terms", permu=10000)) # test for sign. environ. variables
  
  # adjusted R^2
  a1 = round(RsquareAdj(m1)$adj.r.squared,3)
  p1 = round(anova(m1)[1,4],3) # overall test of the significant of the analysis
  a2 = round(RsquareAdj(m2)$adj.r.squared,3)
  p2 = round(anova(m2)[1,4],3) # overall test of the significant of the analysis
  a3 = round(RsquareAdj(m3)$adj.r.squared,3)
  p3 = round(anova(m3)[1,4],3) # overall test of the significant of the analysis
  
  #save results
  lab = ifelse(inp == 'dna_mat','Haplogroup',ifelse(inp == 'auto_mat','Ancestry','Geography'))
  
  dbrda_results = rbind(r1,r2,r3) %>% drop_na(F) %>% mutate(adjR2 = c(a1,a2,a3), p = c(p1,p2,p3)) %>% dplyr::rename(anova_p = 'Pr(>F)') %>% 
    mutate(Response = rownames(.),Test = lab) %>% left_join(.,stepres)
  
  dbr[[lab]] = dbrda_results
}

dbrf = rbindlist(dbr) %>% as_tibble

db_save =  dbrf %>% select(-Df,SumOfSqs) %>% arrange(Response) %>% mutate(padj = p.adjust(p,method='bonferroni'))
dbrda_p = db_save %>% 
  ggplot(aes(y=Response,x=adjR2,fill=Test,label=signif(F,3)))+
  geom_text(hjust=-0.5,position=position_dodge(width=0.9),size=0.5)+
  geom_bar(col='black',stat='identity',position=position_dodge(width=0.9),lwd=0.25)+
  scale_fill_manual(values=c(brewer.pal(3,'Set2')))+
  xlab('Adjusted R2')+ylab('')+
  xlim(c(0,1.05))+
  theme_bw(base_size=6)+
  theme(legend.position='top')

db_save$Response <- factor(db_save$Response,levels=c('Host','Habitat','Egg'))
dbrda_pnt = db_save %>% 
  ggplot(aes(y=Response,x=adjR2,fill=Test))+
  geom_bar(stat='identity',position=position_dodge(width=0.9))+
  scale_fill_manual(values=c(brewer.pal(3,'Set2')))+
  xlab('Adjusted R2')+ylab('')+
  xlim(c(0,1.05))+
  theme_bw(base_size=6)+
  theme(legend.position='top')
dbrda_pnt

pdf('../figures/20240828_dbRDA.pdf',height=2,width=1.5)
dbrda_pnt
dev.off()

write.table(db_save,'../figures/20240828_dbRDA_Results.txt',quote=F,sep='\t',row.names=F)


```

## MNLR

Create a MNLR with egg or host or habitat as the response variable and Geography + Autosomal K + Haplogroups as the predictors. In short:

* Only retain response variables where there are at least n=2 observations
* Downsample all response classes so that all classes have n=2 observations
* Fit 7 multinomial logistic regression models, each with n=100 bootstraps using all combinations of predictors
* Extract AUC, and use the model to predict response variable on the full dataset again (too small for unseen data prediction)
* Repeat the above procedure 100 times so that different downsampled observations are included 
* Determine which classes are predicted correctly (% correct) from the confusion matrix on real / predicted responses across bootstraps

```R
# Egg, Host, Habitat associations, MNLR, Model Selection
#### Find associations between MT/W haplotypes and features 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host')
.libPaths('~/mambaforge/envs/rfs/lib/R/library')
library(tidyverse)
library(viridis)
library(RColorBrewer)
library(vegan)
library(FactoMineR)
library(factoextra)
library(ggrepel)
library(ggpubr)
library(ggdist)
library(nnet)
library(purrr)
library(broom)
library(meRo) #devtools::install_github('merondun/meRo',force=TRUE,upgrade='never')
library(caret)

set.seed(123)
md = read_tsv('~/merondun/cuculus_host/Metadata_Host.txt')

#only grab samples with known egg (hash out filter for related individuals)
md_egg = md %>%
  filter(Analysis_PopulationGenetics == 1) %>%
  drop_na(Egg) %>%
  select(ID = ID, Host = HostParentShort, Egg, Habitat, Haplogroup = Hap, Ancestry = AncestryK5, Geography = KDist)

##### Initialize, summary stats on raw data #####
#ensure they are all factor variables
md_egg = md_egg %>% mutate_at(c('Host','Egg','Habitat','Haplogroup','Ancestry','Geography'),as.factor)

#only retain Host where we have at least 2 cuckoos 
md_egg = md_egg %>% group_by(Host) %>% mutate(TotalHost = n()) %>% ungroup %>% group_by(Habitat) %>% mutate(TotalHabitat = n()) %>% ungroup %>% group_by(Egg) %>% mutate(TotalEgg = n())  %>% ungroup
minobs=2
md_egg = md_egg %>% filter(TotalHost >= minobs & TotalHabitat >= minobs & TotalEgg >= minobs) 
md_egg %>% filter(!ID %in% md$ID)
#write.table(md$ID,file='randomforest/Samples_Retained_Unrelated_MNLR_2024MAR14.txt',quote=F,sep='\t',row.names=F,col.names=F)

# If you want to exclude the blue clades W1, W2, W3! 
md_egg <-  md_egg %>% filter(!grepl('W1|W2|W3',Haplogroup))

#### Count proportions first, count proportions for host and habitat and egg  
hp = md_egg %>% group_by(Host) %>% mutate(Total = n()) %>% 
  pivot_longer(c(Ancestry,Geography,Haplogroup)) %>% ungroup %>% 
  group_by(Host,name,value) %>% 
  summarize(Parts = n(),
            Proportion = Parts / Total,
            Percent = paste0(round( Parts / Total,3)*100,'% (',Total,')')) %>% unique
tp = md_egg %>% group_by(Habitat) %>% mutate(Total = n()) %>% 
  pivot_longer(c(Ancestry,Geography,Haplogroup)) %>% ungroup %>% 
  group_by(Habitat,name,value) %>% 
  summarize(Parts = n(),
            Proportion = Parts / Total,
            Percent = paste0(round( Parts / Total,3)*100,'% (',Total,')')) %>% unique
ep = md_egg %>% group_by(Egg) %>% mutate(Total = n()) %>% 
  pivot_longer(c(Ancestry,Geography,Haplogroup)) %>% ungroup %>% 
  group_by(Egg,name,value) %>% 
  summarize(Parts = n(),
            Proportion = Parts / Total,
            Percent = paste0(round( Parts / Total,3)*100,'% (',Total,')')) %>% unique
egglev = ep %>% ungroup %>% select(Egg) %>% unique %>% mutate(ord = as.numeric(gsub('E','',Egg))) %>% arrange(ord)
ep$Egg = factor(ep$Egg,levels=egglev$Egg)

# Bind them together 
ap = rbind(hp %>% ungroup %>% mutate(Response = Host, variable = 'Host') %>% select(-Host), 
           tp %>% ungroup %>% mutate(Response = Habitat, variable = 'Habitat') %>% select(-Habitat),
           ep %>% ungroup %>% mutate(Response = Egg, variable = 'Egg') %>% select(-Egg))

# Just for ordering the covariates nicely 
ord = ap %>% ungroup %>% select(name,value) %>% unique %>% mutate(ord = as.numeric(str_sub(value,2))) %>% group_by(name) %>% arrange(ord)
ap$value = factor(ap$value,levels=ord$value)

# Plot proportions
app = ap %>% 
  ggplot(aes(y=Response,x=value,fill=Proportion,label=Parts))+
  facet_grid(variable~name,scales='free',space='free')+
  scale_fill_gradient('Proportion Observed',low='yellow',high='red')+
  geom_tile()+
  geom_text(size=1)+
  ylab('')+xlab('')+
  theme_bw(base_size=6)+
  theme(legend.position = 'top')
app

#pdf('figures/Proportions_HostHabitat_2024MAR14.pdf',height=6.5,width=6)
app
#dev.off()

##### Model Selection #####
#assess covariate importance with model selection, using MNLR
vars = c('Habitat','Host','Egg')
vars = 'Egg'

# With bootstrapping
ctrl = trainControl(method = "boot",   # Use bootstrapping
                    number = 100,      # Replicates 
                    summaryFunction = multiClassSummary,  #caret now attempts to stratify levels
                    classProbs = TRUE, #For classification models, to save class probabilities
                    savePredictions = TRUE) #save predictions for each fold

# Change punctuation e.g. 'A. pal' to A_pal'
md_cv = md_egg %>% mutate(Host = gsub('\\. ','_',Host))

# For determining which predictors improve model:
adat = NULL; conf_dat = NULL; class_dat = NULL; new_preds = NULL; counter = 0
for (minnum in c(2)) {  # Loop through 2, 3, or 4 minimum observations per category
  
  # Filter at the very base level, ensuring that across egg / host / habitat we have the same individuals with representative sampling
  md_subbed = md_cv %>% filter(TotalHost >= minnum & TotalHabitat >= minnum & TotalEgg >= minnum) %>%
    mutate_at(c('Host','Egg','Habitat','Haplogroup','Ancestry','Geography'),as.factor)
  
  for (rep in seq(1,10,1)){  # Create 10 replicate models
    for (var in vars) { counter = counter + 1;
    
    # Ensure that we have adequate levels, only a sanity since it is already filtered
    retained = md_subbed %>% ungroup %>% count(!!sym(var)) %>% filter(n >= minnum) %>% pull(var)
    length(retained)
    # Number of samples to subsample
    subsamp = md_subbed %>% ungroup %>% count(!!sym(var)) %>% filter(n >= minnum) %>% summarize(min = min(n)) %>% pull(min)
    
    cat('Downsampling to n = ',subsamp,', requiring min = ',minnum,' for variable: ',var,', ','replicate: ',rep,'\n')
    # Subsampling
    mdi = md_subbed %>%
      filter(!!sym(var) %in% retained) %>%
      group_by(!!sym(var)) %>%
      sample_n(min(n(), subsamp),replace = TRUE) %>%
      ungroup()
    mdi %>% count(!!sym(var))
    
    # Ensure the factor levels are dropped which aren't in the dataset after filtering
    mdi = droplevels(mdi)
    
    # First MNLR on combinations
    formula_1 = as.formula(paste(var, "~ Haplogroup + Ancestry + Geography"))
    m1 = train(formula_1, data = mdi, method = "multinom", trControl = ctrl, metric = "AUC", trace = FALSE)
    
    formula_2 = as.formula(paste(var, "~ Haplogroup + Geography"))
    m2 = train(formula_2, data = mdi, method = "multinom", trControl = ctrl, metric = "AUC", trace = FALSE)
    
    formula_3 = as.formula(paste(var, "~ Haplogroup "))
    m3 = train(formula_3, data = mdi, method = "multinom", trControl = ctrl, metric = "AUC", trace = FALSE)
    
    formula_4 = as.formula(paste(var, "~ Haplogroup + Ancestry"))
    m4 = train(formula_4, data = mdi, method = "multinom", trControl = ctrl, metric = "AUC", trace = FALSE)
    
    formula_5 = as.formula(paste(var, "~ Ancestry"))
    m5 = train(formula_5, data = mdi, method = "multinom", trControl = ctrl, metric = "AUC", trace = FALSE)
    
    formula_6 = as.formula(paste(var, "~ Ancestry + Geography"))
    m6 = train(formula_6, data = mdi, method = "multinom", trControl = ctrl, metric = "AUC", trace = FALSE)
    
    formula_7 = as.formula(paste(var, "~ Geography"))
    m7 = train(formula_7, data = mdi, method = "multinom", trControl = ctrl, metric = "AUC", trace = FALSE)
    
    models = c('m1','m2','m3','m4','m5','m6','m7')
    
    # Extract model fit 
    for (model in models) {
      # Output model fit from confusion matrix
      mo = get(model)
      
      # Get AIC
      final_model = mo$finalModel;
      AIC = AIC(final_model)
      
      # Save the model results
      dat = data.frame(Model = model, Iteration = rep, Variable = var, Subsampled = subsamp, MinObs = minnum,
                       decay = mo$results$decay, logLoss = mo$results$logLoss, logLossSD = mo$results$logLossSD, AIC = AIC, AUC = mo$results$AUC, Accuracy = mo$results$Accuracy, AccuracySD = mo$results$AccuracySD, Kappa = mo$results$Kappa)
      dat_best = dat %>% slice_min(logLoss)
      adat = rbind(adat,dat_best)
      
      # Also save training confusion matrix
      pred = mo$pred
      pred_best = pred %>% filter(decay == dat_best$decay)
      
      # Predict against real data
      predicted_classes = predict(mo, newdata = mdi, type = "raw")
      new = mdi %>% select(reference = !!sym(var)) %>% mutate(predicted = predicted_classes)
      
      conf_new = confusionMatrix(new$predicted, new$reference)
      conf_real = as.data.frame(conf_new$table) %>% mutate(
        Prediction = gsub('_','\\. ',Prediction),
        Reference = gsub('_','\\. ',Reference)) %>%
        mutate(Model = model, Iteration = counter, Variable = var, Subsample = subsamp, MinObs = minnum, AUC=dat_best$AUC,logloss=dat_best$logLoss,Accuracy = dat_best$Accuracy,AccuracySD=dat_best$AccuracySD)
      new_preds = rbind(new_preds,conf_real)
      rm(conf_real,dat,dat_best)
      
    } # Exit model loop
    } # Exit response variable loop
  } # Exit iteration loop
} # Exit minimum samples loop

# Write the no-blue data
write.table(adat,'/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/randomforest/Model_Selection_NoW1W2W3_Boot-2Obs_2024AUG06.txt',quote=F,sep='\t',row.names=F)
write.table(new_preds,'/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/randomforest/ConfusionMatrix_NoW1W2W3_Boot-2Obs_2024AUG06.txt',quote=F,sep='\t',row.names=F)

write.table(adat,'/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/randomforest/Model_Selection_Boot-2Obs_2024AUG06.txt',quote=F,sep='\t',row.names=F)
write.table(new_preds,'/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/randomforest/ConfusionMatrix_Boot-2Obs_2024AUG06.txt',quote=F,sep='\t',row.names=F)

# Or start here and read in saved data 
adat = read_tsv('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/randomforest/Model_Selection_Boot-2Obs_2024AUG06.txt')
conf_mats = read_tsv('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/randomforest/ConfusionMatrix_Boot-2Obs_2024AUG06.txt')

# Plot results 
leg = data.frame(Model = c('m1','m2','m3','m4','m5','m6','m7'),
                 Label = c('W+K+B','W+B','W','W+K','K','K+B','B'))
model_dat = adat %>% left_join(.,leg)

# Plot AUC for each model across all the iterations
model_dat$Label = factor(model_dat$Label, levels=c('W+K+B','W+B','W+K','K+B','B','W','K'))
app = model_dat %>% ggplot(aes(y=Label,x=AUC,fill=Variable))+
  geom_boxplot()+
  scale_fill_manual(values=brewer.pal(3,'RdBu'))+
  xlab('AUC')+ylab('Model Covariates')+
  theme_bw(base_size=6)+theme(legend.position = 'top',strip.text.y.right = element_text(angle = 0))
app 

pdf('figures/2024AUG06_MNLR_ModelSelection_AUC.pdf',height=4,width=6)
app
dev.off()

cols <- brewer.pal(3,'Set2')[c(1,2,3)]
model_dat$Label <- factor(model_dat$Label,levels=c('K','B','W'))
model_dat$Variable <- factor(model_dat$Variable,levels=c('Host','Habitat','Egg'))
auc_plot_input <- model_dat %>%
  filter(MinObs == 2) %>% 
  filter(Label == 'W' | Label == 'K' | Label == 'B') %>% 
  group_by(Label,Variable) %>% 
  sum_stats(AUC)

auc_plot <- auc_plot_input %>% 
  ggplot(aes(y=Variable,fill=Label,x=mean,xmin=conf_low,xmax=conf_high))+
  geom_bar(stat='identity',position=position_dodge(width=0.9))+
  scale_fill_manual(values=cols)+
  xlab('')+ylab('AUC (Mean & 95% CI)')+
  geom_errorbar(width=0.25,position=position_dodge(width=0.9))+
  theme_bw(base_size=6)+
  coord_cartesian(xlim=c(0.7,1))+
  theme(legend.position='top')
auc_plot

pdf('figures/2024AUG06_MNLR_ModelSelection_AUC.pdf',height=2,width=1.5)
auc_plot
dev.off()

write.table(auc_plot_input,file='figures/20240806_AUC_Results_Boot-2Obs.txt',quote=F,sep='\t',row.names=F)

#order full, single plot, make sure the 3 variables are in order 
egglevs = conf_mats %>% filter(Variable == 'Egg') %>% select(Prediction,Variable) %>% mutate(ord = as.numeric(gsub('E','',Prediction)),Prediction) %>% unique %>% arrange(Variable,ord) %>% select(-ord)
hostlevs = conf_mats %>% filter(Variable == 'Host') %>% select(Prediction,Variable) %>% unique %>% arrange(Variable,Prediction) 
hablevs = conf_mats %>% filter(Variable == 'Habitat') %>% select(Prediction,Variable) %>% unique %>% arrange(Variable,Prediction) 
all_levs = rbind(hostlevs,egglevs,hablevs)
conf_mats = conf_mats %>% mutate(Prediction = factor(Prediction,levels=all_levs$Prediction),
                                 Reference = factor(Reference,levels=all_levs$Prediction))

### Plot how the addition of haplogroup improves predictions show K+B (m6) vs W+K+B (m1)
# ancestry + geogrpahy
repredictions_ancgeo = conf_mats %>% 
  filter(Model == 'm6') %>%  #main figure, plot model 2 (W+B) 
  group_by(Variable,MinObs,Iteration,Reference) %>% 
  mutate(TotalIteration = sum(Freq), #within each of the 100 replicates, calculate the frequency of correct classifications
         Proportion = Freq / TotalIteration) %>% ungroup %>% 
  group_by(Prediction,Reference,Variable) %>% 
  sum_stats(Proportion)

#plot 
ancgeo_plot = repredictions_ancgeo %>% 
  ggplot(aes(x=Prediction,y=Reference,fill=mean,col=sd))+
  geom_tile(col=NA)+
  #geom_point(size=0.2)+
  facet_wrap(Variable~.,scales='free')+
  scale_color_continuous(low='white',high='black')+
  scale_fill_continuous(low='white',high='darkblue')+  theme_bw(base_size=6)+
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1),legend.position='top')
ancgeo_plot

# ancestry + geogrpahy + hap 
repredictions_ancgeohap = conf_mats %>% 
  filter(Model == 'm1') %>%  #main figure, plot model 2 (W+B) 
  group_by(Variable,MinObs,Iteration,Reference) %>% 
  mutate(TotalIteration = sum(Freq), #within each of the 100 replicates, calculate the frequency of correct classifications
         Proportion = Freq / TotalIteration) %>% ungroup %>% 
  group_by(Prediction,Reference,Variable) %>% 
  sum_stats(Proportion)

#plot 
ancgeohap_plot = repredictions_ancgeohap %>% 
  ggplot(aes(x=Prediction,y=Reference,fill=mean,col=sd))+
  geom_tile(col=NA)+
  #geom_point(size=0.2)+
  scale_color_continuous(low='white',high='black')+
  facet_wrap(Variable~.,scales='free')+
  scale_fill_continuous(low='white',high='darkblue')+  theme_bw(base_size=6)+
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1),legend.position='top')
ancgeohap_plot

pdf('figures/2024AUG08_MNLR_ConfusionMatrix-Repredictions_M1vsM6.pdf',height=4.5,width=5)
ggarrange(ancgeo_plot,ancgeohap_plot,nrow = 2,common.legend = TRUE)
dev.off()

```

## Discordance Analysis: Egg

Determine the number of shifts for each egg type, binarizing each egg and comparing mtDNA and autosomal trees:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=5
#SBATCH --time=48:00:00

auto=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/merged/snps_only/chr_1.SNPS.vcf.gz

mtdna=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/merged/snps_only/chr_MT.SNPS.vcf.gz

#Subset VCFS
bcftools view --samples-file Samples.list -Ou ${auto} | \
                        bcftools view --types snps --min-alleles 2 --max-alleles 2 --min-af 0.05 --max-af 0.95 -Ou | \
                                        bcftools +prune -m 0.2 --window 5kb -Oz -o autos_LD_n89.vcf.gz
bcftools index autos_LD_n89.vcf.gz

# and mtDNA
bcftools view --samples-file Samples.list -Ou ${mtdna} | \
                        bcftools view --types snps --min-alleles 2 --max-alleles 2 --min-af 0.05 --max-af 0.95 -Oz -o chr_MT_n89.vcf.gz
bcftools index chr_MT_n89.vcf.gz

#create tree, autosomes
iqtree --redo -keep-ident -T 5 -s ml_trees/autos_LD_n89.min4.phy --seqtype DNA -m "GTR+ASC" -B 1000
iqtree --redo -keep-ident -T 5 -s ml_trees/autos_LD_n89.min4.phy.varsites.phy --seqtype DNA -m "GTR+ASC" -B 1000

#create tree, mtdna
iqtree --redo -keep-ident -T 5 -s ml_trees/chr_MT_n89.min4.phy --seqtype DNA -m "GTR+ASC" -B 1000
iqtree --redo -keep-ident -T 5 -s ml_trees/chr_MT_n89.min4.phy.varsites.phy --seqtype DNA -m "GTR+ASC" -B 1000
```

And then compare egg shifts across both mtDNA and autosomal tree:

```R
 

# Read in trees and metadata 
md <- read_tsv('~/merondun/cuculus_host/Metadata_Host.txt') %>% filter(Analysis_PopulationGenetics == 1) %>% drop_na(Egg)
md$dummy <- 1 # Assign fake dummy var so that all connections are plotted 

### mtDNA tree
m = read.iqtree('ml_trees/chr_MT_n89.min4.phy.varsites.phy.contree')
m1 <- root(as.phylo(m),outgroup = '387_CP_MBW_RUS_F')
mt_tree <- drop.tip(m1,c('387_CP_MBW_RUS_F','386_CP_MBW_RUS_M'))

### AUTOSOME tree 
a = read.iqtree('ml_trees/autos_LD_n89.min4.phy.varsites.phy.contree')
a1 <- root(as.phylo(a),outgroup = '387_CP_MBW_RUS_F')
a_tree <- drop.tip(a1,c('387_CP_MBW_RUS_F','386_CP_MBW_RUS_M'))

# Store results
results <- list()
for (egg in unique(md$Egg)) { 
  for (tree in c('mt_tree','a_tree')) { 
    
    cat('Working on egg type: ',egg,' for tree: ',tree,'\n')
    
    # Change egg to binary trait, only target egg is 1 all else is 0 
    md_mod <- md %>% mutate(Egg = ifelse(Egg == egg,egg,'E0'))
    
    # Generate base tree 
    targ_tree <- get(tree)
    ggt <- ggtree(targ_tree, layout = "circular",branch.length='none') %<+% md_mod
    
    #grab only egg
    phenos <- as.data.frame(ggt$data %>% filter(isTip == TRUE))
    egg_mat <- as.matrix(phenos %>% select(Egg))
    phenotypes <- setNames(as.vector(egg_mat[, 1]), targ_tree$tip.label)
    
    #inspect tree
    ggtree(targ_tree, layout = "rectangular",branch.length='none') %<+% md_mod + 
      geom_tippoint(aes(fill=Egg),pch=21,size=2)+
      scale_fill_brewer(palette='Set2')
    
    #Plot probabilities 
    t2 <- multi2di(targ_tree)
    
    # Quick ape method 
    fitER <- ape::ace(phenotypes,t2,model="ER",type="discrete")
    
    # Determine ancestral state likelihood number of shifts using MCMC
    sim <- make.simmap(t2, phenotypes, model="ER", nsim=100,Q='mcmc')
    
    # Output stats on the number of shifts across the 100 replicates 
    vec <- data.frame(shifts = countSimmap(sim)$Tr[,1])
    shifts <- vec %>% mutate(Egg = egg, Tree = tree)
    
    results[[paste0(egg,'_',tree)]] <- shifts
    
  }
}

full_results <- rbindlist(results)

write.table(full_results,file='Results_Binary_EggComparison_2024JULY16.txt',quote=F,sep='\t',row.names=F)

# Read in 
full_results <- read_tsv('Results_Binary_EggComparison_2024JULY16.txt')

# Add egg colors
md = read_tsv('~/merondun/cuculus_host/Metadata_Host.txt')
egglevs = md %>% filter(Analysis_Mantel == 1) %>% select(Egg) %>% unique %>% mutate(ord = as.numeric(gsub('E','',Egg))) %>% arrange(ord) %>% mutate(col = viridis(11, option='turbo'))
full_results$Egg = factor(full_results$Egg,levels=egglevs$Egg)

# Each egg / tree has 100 bootstrapps, so add an identifier for each
fr <- full_results %>% 
  group_by(Egg, Tree) %>% 
  mutate(rep = row_number()) 

# And for each bootstrap compare mtDNA vs autosomal shifts 
mt_vs_aut <- fr %>% 
  pivot_wider(names_from = Tree,values_from = shifts) %>% 
  mutate(shifts = mt_tree - a_tree)

# Calculate summary stats incl. 95% CIs
confs <- mt_vs_aut %>% group_by(Egg) %>% sum_stats(shifts)

plot_mtaut <- mt_vs_aut %>% ungroup %>% 
  ggplot(aes(x=Egg,y=shifts,fill=Egg))+
  geom_hline(yintercept=0,lty=2)+
  geom_hline(yintercept=c(-1,1),lty=3,col='darkgray')+
  geom_errorbar(data=confs,width=0.2,col='black',aes(x=Egg,ymin=conf_low,ymax=conf_high),inherit.aes=FALSE,position=position_nudge(x=-0.15))+
  ggdist::stat_halfeye(adjust = 1,width = .3,.width = 0,justification = -.2, point_colour = NA,alpha = 0.95,normalize='groups')+
  #geom_boxplot(width = .15,outlier.shape = NA, alpha = 0.3,position=position_nudge(x=0.15)) +
  #gghalves::geom_half_point(aes(col=Egg),side='l',range_scale = .4,alpha = .3)+
  scale_fill_manual(values=egglevs$col,breaks=egglevs$Egg)+
  scale_color_manual(values=egglevs$col,breaks=egglevs$Egg)+
  theme_bw()+
  ylab('')+xlab('')+
  coord_flip()
plot_mtaut

png('../../figures/20240717_mtDNA-vs-Autosomal-shifts-MCMC.png',units='in',res=300,height=4,width=4)
pdf('../../figures/20240717_mtDNA-vs-Autosomal-shifts-MCMC.pdf',height=4,width=4)
plot_mtaut
dev.off()

#### Plot tree discordance ####
# This will also estimate transitions from each egg type to other eggs,  but doesn't seem as reliable as binary above approach

full_tree_results <- list()
full_egg_results <- list()

for (tree in c('mt_tree','a_tree')) {
  
  cat('Running full ML reconstruction on tree: ',tree,'\n')
  
  # Generate base tree 
  targ_tree <- get(tree)
  ggt <- ggtree(targ_tree, layout = "circular",branch.length='none') %<+% md 
  
  #grab only egg
  phenos <- as.data.frame(ggt$data %>% filter(isTip == TRUE))
  egg_mat <- as.matrix(phenos %>% select(Egg))
  phenotypes <- setNames(as.vector(egg_mat[, 1]), targ_tree$tip.label)
  
  #Plot probabilities 
  t2 <- multi2di(targ_tree)
  
  # Quick ape method 
  fitER <- ape::ace(phenotypes,t2,model="ER",type="discrete")
  
  # Determine ancestral state likelihood number of shifts using MCMC
  simfull <- make.simmap(t2, phenotypes, model="ER", nsim=100,Q='mcmc')
  
  # Output stats on the number of shifts across the 100 replicates 
  vec <- data.frame(shifts = countSimmap(simfull)$Tr[,1])
  shifts <- vec %>% mutate(Tree = tree, MLloglik = fitER$loglik, MLest = fitER$rates, MLse = fitER$se) 
  
  # Extract number straight from the MCMC approach
  mat <- countSimmap(simfull)$Tr
  col_names <- colnames(mat)
  
  # Calculate average transitions for each Egg
  calculate_statistics <- function(mat, prefix) {
    # Identify columns that contain the specific egg prefix
    relevant_cols <- grep(paste0("^", prefix, ",E|E,", prefix, "$"), col_names, value = TRUE)
    
    # Sum the values in these columns
    total_transitions <- rowSums(mat[, relevant_cols])
    
    # Calculate statistics
    stats <- list(
      Mean = mean(total_transitions),
      Median = median(total_transitions),
      SD = sd(total_transitions)
    )
    
    return(stats)
  }
  
  # Apply to each E group
  E_groups <- paste("E", 1:11, sep="")  # Adjust based on your specific E groups
  averages <- sapply(E_groups, function(E) calculate_statistics(mat, E))
  
  egg_shifts <- as.data.frame(t(averages)) %>% mutate(Egg = rownames(.), Tree = tree)
  rownames(egg_shifts) <- NULL
  
  full_tree_results[[tree]] <- shifts
  full_egg_results[[tree]] <- egg_shifts
  
  # Extract nodes and the proportions for pies
  nodes <- data.frame(
    node=1:t2$Nnode+Ntip(t2),
    fitER$lik.anc)
  
  ## cols parameter indicate which columns store stats
  pies <- nodepie(nodes, cols=2:12,outline.color='black',outline.size = 0.1)
  pies <- lapply(pies, function(g) g+scale_fill_manual(values = egglevs$col,breaks=egglevs$Egg))
  
  t3 <- full_join(t2, data.frame(label = names(phenotypes), stat = phenotypes ), by = 'label')
  tp <- ggtree(t3,layout='rectangular',branch.length = 'none') %<+% md
  tp$data$dummy <- 1
  tp_final <- tp + geom_inset(pies, width = .09, height = .09)+ 
    geom_tippoint(aes(fill=Egg),pch=21,size=1)+
    scale_fill_manual(values=egglevs$col,breaks=egglevs$Egg) 
  assign(paste0(tree,'_pies'),tp_final)
  
}

# Grab the model results
full_search_res <- rbindlist(full_tree_results)
full_search_res %>% group_by(Tree,MLloglik,MLest,MLse) %>% sum_stats(shifts)

# Also bind the egg data 
full_search_eggs <- rbindlist(full_egg_results)
full_search_eggs
full_search_eggs$Mean <- unlist(full_search_eggs$Mean)
full_search_eggs$Median <- unlist(full_search_eggs$Median)
full_search_eggs$SD <- unlist(full_search_eggs$SD)

write.table(full_search_res,file='Full_Search_Results_20240806.txt',quote=F,sep='\t',row.names=F)
write.table(full_search_eggs,file='Full_Search_EggResults_20240806.txt',quote=F,sep='\t',row.names=F)

full_search_res <- read_tsv('Full_Search_Results_20240806.txt')

# Plot all connections 
discord <- simple.tanglegram(tree1=mt_tree_pies,tree2=a_tree_pies,column=dummy,value=1,t2_pad=2,l_color='black',tiplab=F)

pdf('../../figures/20240717_TreeCompare-mtDNA-Auto-NodePiesML.pdf',height=5,width=7)
discord
dev.off()

# Swap 
discord2 <- simple.tanglegram(tree1=a_tree_pies,tree2=mt_tree_pies,column=dummy,value=1,t2_pad=2,l_color='black',tiplab=F)

pdf('../../figures/20240717_TreeCompare-Auto-mtDNA-NodePiesML.pdf',height=5,width=7)
discord2
dev.off()

# And as before, add the shift data 
frf <- full_search_res %>% select(Tree,shifts) %>% 
  group_by(Tree) %>% 
  mutate(rep = row_number()) 

# And for each bootstrap compare mtDNA vs autosomal shifts 
mt_vs_aut_full <- frf %>%
  pivot_wider(names_from = Tree,values_from = shifts) %>% 
  mutate(shifts = mt_tree - a_tree)

# Calculate summary stats incl. 95% CIs
confs_full <- mt_vs_aut_full %>% sum_stats(shifts)

## For binding with the overall reconstruction! ~~ load in data above FIRST ~~
mt_vs_autboth <- rbind(mt_vs_aut,mt_vs_aut_full %>% mutate(Egg = 'Multiclass'))
confs_all <- rbind(confs,confs_full %>% mutate(Egg = 'Multiclass'))
eggboth <- rbind(egglevs,data.frame(Egg='Multiclass',ord=12,col='grey80')) %>% arrange(desc(ord))
eggboth$Egg <- factor(eggboth$Egg,levels=eggboth$Egg)
mt_vs_autboth$Egg <- factor(mt_vs_autboth$Egg,levels=eggboth$Egg)
confs_all$Egg <- factor(confs_all$Egg,levels=eggboth$Egg)

# Plot
plot_mtaut <- mt_vs_autboth %>% ungroup %>% 
  ggplot(aes(x=Egg,y=shifts,fill=Egg))+
  geom_hline(yintercept=0,lty=2)+
  geom_hline(yintercept=c(-1,1),lty=3,col='darkgray')+
  geom_errorbar(data=confs_all,width=0.2,col='black',aes(x=Egg,ymin=conf_low,ymax=conf_high),inherit.aes=FALSE,position=position_nudge(x=-0.15))+
  geom_boxplot(width = .15,outlier.shape = NA, alpha = 0.9,lwd=0.25,position=position_nudge(x=0.15)) +
  #ggdist::stat_halfeye(width = .3,.width = 0,justification = -.2, point_colour = NA,alpha = 0.95,normalize='groups')+
  scale_fill_manual(values=eggboth$col,breaks=eggboth$Egg)+
  scale_color_manual(values=eggboth$col,breaks=eggboth$Egg)+
  theme_bw()+
  ylab('')+xlab('')+
  coord_flip(ylim=c(-15,15))
plot_mtaut

png('../../figures/20240807_mtDNA-vs-Autosomal-shifts-MCMC-Multiclass.png',units='in',res=300,height=4,width=4)
pdf('../../figures/20240807_mtDNA-vs-Autosomal-shifts-MCMC-Multiclass-Boxes.pdf',height=4.5,width=4)
plot_mtaut
dev.off()


# Plot full search egg shifts
eggdf <- as.data.frame(full_search_eggs)
mt_vs_aut <- eggdf %>% 
  select(Mean,Tree,Egg) %>% 
  pivot_wider(names_from = Tree,values_from = Mean) %>% 
  mutate(mt_tree = map_dbl(mt_tree, unlist),
         a_tree = map_dbl(a_tree, unlist),
         shifts = mt_tree - a_tree)
mt_vs_aut$Egg = factor(mt_vs_aut$Egg,levels=egglevs$Egg)


plot_mtaut_alt <- mt_vs_aut %>% ungroup %>% 
  ggplot(aes(x=Egg,y=shifts,fill=Egg))+
  geom_hline(yintercept=0,lty=2)+
  geom_hline(yintercept=c(-1,1),lty=3,col='darkgray')+
  geom_point(pch=21,size=3)+
  scale_fill_manual(values=egglevs$col,breaks=egglevs$Egg)+
  scale_color_manual(values=egglevs$col,breaks=egglevs$Egg)+
  theme_bw()+
  ylab('')+xlab('')+
  coord_flip()
plot_mtaut_alt

#inspect tree
for (egg in unique(egglevs$Egg)) { 
  
  cat('Making individual plot for egg: ',egg,'\n')
  md_inspect <- md %>% mutate(Egg = ifelse(Egg == egg,'1','0'))
  t1 <- ggtree(mt_tree, layout = "rectangular",branch.length='none') %<+% md_inspect + 
    geom_tippoint(aes(fill=Egg),pch=21,size=1)+
    scale_fill_manual(values=rainbow(3))
  t2 <- ggtree(a_tree, layout = "rectangular",branch.length='none') %<+% md_inspect + 
    geom_tippoint(aes(fill=Egg),pch=21,size=1)+
    scale_fill_manual(values=rainbow(3))
  ts <- simple.tanglegram(tree1=t1,tree2=t2,column=Egg,value=1,t2_pad=2,l_color='black',tiplab=F) + ggtitle(egg)
  assign(paste0('p',egg),ts)
  
}

full_plots <- ggarrange(pE1,pE2,pE3,pE4,pE5,pE6,pE7,pE8,pE9,pE10,pE11,common.legend = TRUE,nrow=4,ncol=3)

pdf('../../figures/20240722_mtDNA-vs-Auto_Individual_Eggs.pdf',height=7,width=5)
full_plots
dev.off()
```



# Analyses: Egg Hunt

## GWAS

### Sample Selection

```R
#### Identify samples to use for GWAS, unrelated individuals - but include males 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/manyhost_hunt/all_samples')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(viridis)
library(ggpubr)
library(sf)
library(ggspatial)
md = read_tsv('~/merondun/cuculus_host/Metadata_Host.txt')
eggcols = md %>% ungroup %>% select(Egg) %>% unique %>% mutate(ord = as.numeric(gsub('E','',Egg))) %>% arrange(ord) %>% 
  drop_na(Egg) %>% mutate(col = viridis(12,option='turbo')[1:11]) 
eggcols = rbind(eggcols,data.frame(Egg='E6W3',ord=6,col='white'))

##### Identify Samples   #####
md = read_tsv('~/merondun/cuculus_host/Metadata_Host.txt')
md = md %>% drop_na(Egg) %>%
  mutate(Group = ifelse(Egg == 'E6' & Hap == 'W3','E6W3',Egg))

#so while we will analyze for those groups: we will still include other samples as the background
bg = md %>% drop_na(Egg) %>% filter(Analysis_PopulationGenetics == 1)  %>% select(ID)
write.table(bg$ID,file='AllSamples.list',quote=F,sep='\t',row.names=F,col.names=F)
write.table(mds %>% count(Group) %>% pull(Group),file='GROUPS.list',quote=F,sep='\t',row.names=F,col.names=F)

#save for all v 1 
allv1 = md %>% select(ID,Group)
write.table(allv1,file='AllSamples.pop',quote=F,sep='\t',row.names=F,col.names=F)

#jitter points up to 1 lat/long for viewing
md = md %>% mutate(LatJit = jitter(Latitude,amount =2),
                   LonJit = jitter(Longitude,amount=2))
sites = st_as_sf(md, coords = c("LonJit", "LatJit"), 
                 crs = 4326, agr = "constant")

#set up map and convert df to coordinate frame
world = map_data("world")

gwas_sample_plot = ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), col='grey90', fill='white') +
  geom_sf(data = sites, 
          aes(fill=as.factor(Group)),
          size=3,show.legend = T,pch=21) +
  xlab('')+ylab('')+
  coord_sf(xlim = c(min(md$Longitude)-5, max(md$Longitude)+5), 
           ylim = c(min(md$Latitude)-5, max(md$Latitude)+5), expand = FALSE)+
  scale_fill_manual(values=eggcols$col,breaks=eggcols$Egg)+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"))+
  theme(legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))+
  annotation_scale(line_width=0.5)+
  guides(fill=guide_legend(override.aes=list(shape=21)))
gwas_sample_plot

pdf('../../figures/ManyHost_SpatialDistribution-GWAS_2024APR2.pdf',height=4,width=7)
gwas_sample_plot
dev.off()
```

![image-20240403160833844](C:\Users\herit\AppData\Roaming\Typora\typora-user-images\image-20240403160833844.png)

### Subset VCF

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=5
#SBATCH --time=200:00:00

# mamba activate snps
# for CHR in $(cat Chromosomes.list); do sbatch -J FST_${CHR} ~/merondun/cuculus_host/gene_hunt/1A.Subset_Groups_AllinOne_FST.sh ${CHR}; done
CHR=$1

#genotypes BELOW this will be set to missing
MINDP=3

if [[ $CHR = 'chr_W' || $CHR = 'chr_MT' ]]; then
        PLOIDY=1

        echo "FILTERING AND MERGING VARIANT SITES FOR ${CHR}, PLOIDY: ${PLOIDY}"
        bcftools view --threads 5 --samples-file AllSamples.list --force-samples -Ou ../../../merged/snps_only/${CHR}.SNPS.vcf.gz | \
                bcftools view --types snps --min-alleles 2 --max-alleles 2 --threads 5 | \
                #set genotypes below MINDP to missing
                bcftools +setGT -- -t q -i "FMT/DP < ${MINDP}" -n "./." | \
                #set het genotypes to missing based on binomial test
                bcftools +setGT -Ou -- -t "b:AD<1e-5" -n "./." | \
                #set weakly het genotypes to major allele
                bcftools +setGT -Ou -- --target-gt q --new-gt M -i 'GT=="het"' | \
                #set to haploid, can skip this for most purposes
                bcftools +fixploidy -Ou - -- -f ${PLOIDY} | \
                #update AC fields
                bcftools +fill-tags -Ou -- -t AC,AN | \
                bcftools view --min-alleles 2 --max-alleles 2 --types snps --min-ac 1 -i 'MQ > 40 & F_MISSING < 0.1' -Oz -o vcfs/${CHR}.SNP.DP3.vcf.gz
        bcftools index --threads 5 vcfs/${CHR}.SNP.DP3.vcf.gz

else
        PLOIDY=2

        echo "FILTERING AND MERGING VARIANT SITES FOR ${CHR}, PLOIDY: ${PLOIDY}"
        bcftools view --threads 5 --samples-file AllSamples.list -Ou ../../../merged/snps_only/${CHR}.SNPS.vcf.gz | \
                bcftools view --types snps --min-alleles 2 --max-alleles 2 --threads 5 | \
                #set genotypes below MINDP to missing
                bcftools +setGT -- -t q -i "FMT/DP < ${MINDP}" -n "./." | \
                #update AC fields
                bcftools +fill-tags -Ou -- -t AC,AN | \
                bcftools view --min-alleles 2 --max-alleles 2 --types snps --min-ac 1 -i 'MQ > 40 & F_MISSING < 0.1' -Oz -o vcfs/${CHR}.SNP.DP3.vcf.gz
        bcftools index --threads 5 vcfs/${CHR}.SNP.DP3.vcf.gz

fi

```

Annotate variants:

```R
#!/usr/bin/env rscript
args <- commandArgs(trailingOnly = TRUE)

library(tidyverse)
library(data.table)
library(VariantAnnotation)
library(GenomicFeatures)

chr = args[1]
txdb = makeTxDbFromGFF("/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR-chr_MT.gff")
fasta_seq <- readDNAStringSet("/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa")

#load in vcf
vcf = readVcf(paste0("vcfs/",chr,'.SNP.DP3.vcf.gz'), genome = "cuckoo")

#annotate
if (chr == 'chr_MT') {
  code = getGeneticCode("SGC1")
  coding_predictions = predictCoding(query = vcf, subject = txdb, seqSource = fasta_seq, genetic.code= code) #for chr_MT
  dat = as_tibble(coding_predictions) %>% dplyr::select(seqnames,start,REF,varAllele,CONSEQUENCE,GENEID)
} else {
  code = getGeneticCode("SGC0")
  coding_predictions = predictCoding(query = vcf, subject = txdb, seqSource = fasta_seq, genetic.code= code,) #for other chrs
  dat = as_tibble(coding_predictions) %>% dplyr::select(seqnames,start,REF,varAllele,CONSEQUENCE,GENEID)

}

write.table(dat %>% unique %>% arrange(seqnames,start),file=paste0('raw_dnds/Annotated_Variants_',chr,'__2024APR1.txt'),quote=F,sep='\t',row.names=F)

```

and submit:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_highmem
#SBATCH --cpus-per-task=10
#SBATCH --time=8:00:00

CHR=$1

#Rscript DNDS.R ${CHR}

sed '1d' raw_dnds/Annotated_Variants_${CHR}__2024APR1.txt | awk '{OFS="\t"}{print $1, $2, $2, $5, $6}' | sed 's/gene-//g' | sed 's/ID=//g' > raw_dnds/${CHR}.bed

```



### GEMMA

Prep beds:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=5
#SBATCH --time=200:00:00

CHR=$1

#Do first
plink --threads 10 --vcf vcfs/${CHR}.SNP.DP3.vcf.gz --allow-extra-chr --const-fid --set-missing-var-ids @:# --make-bed --out beds/${CHR} --chr-set 39
awk '{print $2, $2, $3, $4, $5, "1"}' beds/${CHR}.fam > beds/${CHR}.tmp
mv beds/${CHR}.tmp beds/${CHR}.fam

```

From the same autosomal p-distance matrix used for mantel tests, subset the individuals in the same order as the bed file to use as the covariance matrix:

```bash
#### Plot the p-distance matrix on the subset GWAS individuals 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/manyhost_hunt/males')
.libPaths('~/mambaforge/envs/R/lib/R/library')
library(tidyverse)
library(viridis)
library(ggpubr)
library(meRo)
library(RColorBrewer)
library(ggpubr)
library(scales)

md = read_tsv('~/merondun/cuculus_host/Metadata_Host.txt')
eggcols = md %>% ungroup %>% select(Egg) %>% unique %>% mutate(ord = as.numeric(gsub('E','',Egg))) %>% arrange(ord) %>% 
  drop_na(Egg) %>% mutate(col = viridis(11,option='turbo'))
eggcols = rbind(eggcols,data.frame(Egg='E6W3',ord=12,col='white'))

#load in autosomal distance, calculated: ~/modules/VCF2Dis-1.50/bin/VCF2Dis -InPut autos_canorus_LD.vcf.gz -OutPut autos_canorus_LD.pdist 
auto = read.table('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/merged_unrelfull/autos_canorus_LD.pdist',header=F)
#add proper names, because VCF2DIS truncates them
auto = auto %>% mutate(IDNum = gsub('_.*','',V1))
auto_id = left_join(auto,md %>% select(ID) %>% mutate(IDNum = gsub('_.*','',ID)))  
ord = auto_id$ID
autos = auto_id %>% select(!c(IDNum,ID,V1))
names(autos) = ord
rownames(autos) = ord

#the order we want, from the bed file 
ids = read.table('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/manyhost_hunt/males/beds/chr_1.fam')

#extract columns in same order as geographic distance
common_names = intersect(rownames(autos),ids$V1)
auto_aligned = autos[common_names, common_names] #subset 'auto' matrix to match the 'geo_mat' matrix rows and columns
auto_mat = as.matrix(auto_aligned)

#visualize autosomal distance matrix in terms of ancestry K = 5
kcols = data.frame(Kcols = gray.colors(5,start=0.05,end=0.95)[c(1,4,3,5,2)], AncestryK5 = paste0('K',seq(1,5,1)))
auto_mds = as.data.frame(cmdscale(auto_mat,2))
auto_mds$ID = rownames(auto_mds)
auto_mds_input = left_join(auto_mds,md)
auto_p = ggplot(auto_mds_input,aes(x=-V1,y=V2,fill=AncestryK5,shape=AncestryK5)) + 
  geom_point(size=1.5) + theme_bw() + xlab('MDS1')+ylab('MDS2')+ggtitle('Autosomal Distance')+
  scale_fill_manual(name = "Ancestry K5", values = setNames(kcols$Kcols, kcols$Kcluster)) +
  scale_shape_manual(name = "Ancestry K5", values = setNames(c(21, 22, 23, 24, 25), kcols$Kcluster)) +
  theme(legend.position = "top", # Moves legend to top
        legend.text = element_text(size = 6), # Adjusts text size
        legend.title = element_text(size = 6))
auto_p

pdf('../../figures/GWAS_Covariance_Matrix_2024APR3.pdf',height=3,width=3)
auto_p
dev.off()

write.table(auto_mat,file='beds/Covariance_Matrix.txt',quote=F,sep=' ',row.names=F,col.names=F)

##### Females
#the order we want, from the bed file 
ids = read.table('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/manyhost_hunt/males/beds/chr_W.fam')

#extract columns in same order as geographic distance
common_names = intersect(rownames(autos),ids$V1)
auto_aligned = autos[common_names, common_names] #subset 'auto' matrix to match the 'geo_mat' matrix rows and columns
auto_mat = as.matrix(auto_aligned)

#visualize autosomal distance matrix in terms of ancestry K = 5
kcols = data.frame(Kcols = gray.colors(5,start=0.05,end=0.95)[c(1,4,3,5,2)], Kcluster = paste0('K',seq(1,5,1)))
auto_mds = as.data.frame(cmdscale(auto_mat,2))
auto_mds$ID = rownames(auto_mds)
auto_mds_input = left_join(auto_mds,md)
auto_p = ggplot(auto_mds_input,aes(x=-V1,y=V2,fill=AncestryK5,shape=AncestryK5)) + 
  geom_point(size=1.5) + theme_bw() + xlab('MDS1')+ylab('MDS2')+ggtitle('Autosomal Distance')+
  scale_fill_manual('K',values=kcols$Kcols,breaks=kcols$Kcluster)+
  scale_shape_manual(values=c(21,22,23,24,25))+
  guides(fill=guide_legend(override.aes=list(shape=21)))+
  theme(legend.position = "top", # Moves legend to top
        legend.text = element_text(size = 6), # Adjusts text size
        legend.title = element_text(size = 6))
auto_p

pdf('../../figures/GWAS_Covariance_Matrix-Females_2024APR3.pdf',height=3,width=3)
auto_p
dev.off()

write.table(auto_mat,file='beds/Covariance_Matrix_N60.txt',quote=F,sep=' ',row.names=F,col.names=F)
```

![image-20240403111915618](C:\Users\herit\AppData\Roaming\Typora\typora-user-images\image-20240403111915618.png)

And then submit GWAS, using an all v one strategy:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=5
#SBATCH --time=200:00:00

TARGET=$1

mkdir gemma gemma/${TARGET}
cd gemma/${TARGET}

for CHR in $(cat ../../Chromosomes.list); do

if [[ $CHR = 'chr_W' ]]; then
        #grab phenotypes, females only
        awk -v target="$TARGET" '{if ($2 == target) print 1; else print 0}' ../../Females.pop > ${CHR}.phenotypes

        #run gemma
        ../../gemma-0.98.6-pre1 -bfile ../../beds/${CHR} -k ../../beds/Covariance_Matrix_N60.txt -p ${CHR}.phenotypes -lm 1 -o ${CHR} -miss 0.1 -maf 0.01

else

        #grab phenotypes
        awk -v target="$TARGET" '{if ($2 == target) print 1; else print 0}' ../../AllSamples.pop > ${CHR}.phenotypes

        #run gemma
        ../../gemma-0.98.6-pre1 -bfile ../../beds/${CHR} -k ../../beds/Covariance_Matrix.txt -p ${CHR}.phenotypes -lm 1 -o ${CHR} -miss 0.1 -maf 0.01

fi

#save output
awk -v t=${TARGET} '{OFS="\t"}{print $1, $3, $3, $11, t}' output/${CHR}.assoc.txt | \
        sed '1d' | \
        bedtools intersect -a - -b ../../raw_dnds/${CHR}.bed -wao | \
        awk '{OFS="\t"}{print $1, $2, $4, $5, $9, $10}' | \
        sed 's/gene-//g' | sed 's/ID=//g' > ../../gwas/${TARGET}_${CHR}.gwas
done
```

### Plot

```R
#### Identify candidates from GWAS results
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/manyhost_hunt/males')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(viridis)
library(ggpubr)
library(meRo)
library(RColorBrewer)
library(zoo)
library(ggpubr)
library(scales)
library(data.table)

md = read_tsv('~/merondun/cuculus_host/Metadata_Host.txt')
eggcols = md %>% ungroup %>% select(Egg) %>% unique %>% mutate(ord = as.numeric(gsub('E','',Egg))) %>% arrange(ord) %>% 
  drop_na(Egg) %>% mutate(col = viridis(12,option='turbo')[1:11])
eggcols = rbind(eggcols,data.frame(Egg='E6W3',ord=5,col='white')) %>% arrange(ord)

pops = read_tsv('AllSamples.pop',col_names=F)
names(pops) = c('ID','EggType')
left_join(pops,md) %>% count(Egg,Hap)

##### Load and prep GWAS SNP data   ######
tidy_bp = fread('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/manyhost_hunt/males/GWAS_Results_2024APR2-MultiGene.txt.gz')
names(tidy_bp) = c('chr','start','end','p','Egg','effect','gene','intersect_gene')
setDT(tidy_bp)
genes = unique(tidy_bp, by = c("chr", "start", "p", "Egg", "intersect_gene")) # For sup fig all v one GWAS 
#genes = unique(tidy_bp, by = c("chr", "start", "p", "Egg", "effect", "gene"))

#strip chr_, and assign AvZ
genes[, chr := gsub('chr_', '', chr)][
  , AvZ := fifelse(chr == 'W', 'W',
                   fifelse(chr == 'Z', 'Z',
                           fifelse(chr == 'MT', 'MT', 'Autosome')))]

##### Proportions of Significant SNPs #####
td = genes %>% group_by(Egg) %>% mutate(padj = p.adjust(p,method='bonferroni')) %>% ungroup

candidates = td %>% filter(padj < 0.05)

#total windows, count the proportion of sites which are FST = 1.0 in ALL contrasts for a certain mutation level
chr_windows = td %>% group_by(Egg,AvZ) %>% count() %>% ungroup %>% group_by(Egg,AvZ) %>% slice_max(n) %>% select(Egg,AvZ,NumberWindows = n) %>% distinct
prop = candidates  %>%
  count(AvZ,Egg) %>% unique %>%  #count number of SNPs by phylo 
  left_join(chr_windows,.) %>% ungroup %>% 
  replace_na(list(n=0)) %>% 
  ungroup %>%
  mutate(Proportion = n/NumberWindows)

#plot, divide into top and bottom plots for visualization 
prop$Egg = factor(prop$Egg,levels=eggcols$Egg)
prop$AvZ = factor(prop$AvZ,levels=c('Autosome','Z','W','MT'))
props = prop %>% ggplot(aes(y=Egg,x=Proportion,fill=Egg,label=n))+
  geom_bar(stat='identity',position=position_dodge(width=0.9),col='black')+
  #geom_text(position=position_dodge(width=0.9),hjust=-.25,size=1)+
  scale_fill_manual(values=eggcols$col,breaks=eggcols$Egg)+
  scale_x_continuous(n.breaks = 4)+
  facet_grid(.~AvZ,scales='free')+
  theme_bw(base_size=6)+xlab('Proportion of FST = 1.0 SNPs Present in At Least 1 Comparison')+ylab('')+
  theme(legend.position = "top", # Moves legend to bottom
        legend.text = element_text(size = 6), # Adjusts text size
        legend.title = element_text(size = 6)) + # Adjusts title size
  guides(fill = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 1)))
props

pdf('../../figures/Proportions_GWAS_2024APR3.pdf',height=3,width=6)
props
dev.off()

##### Candidate Genes, here use the nonsynonymnous SNPs after bonferonni correction #####
#filter for bonferonni corrected significant sites 
cands = nonsyn[p <= 0.05/21210834]

#replace LOC IDs with known homologs
locs = read_tsv('../LOC_Lookup.txt')
cands = left_join(cands,locs %>% dplyr::rename(gene=Label)) %>% 
  mutate(gene = ifelse(is.na(ID),gene,paste0(ID,'-L'))) %>% select(-ID) %>% unique

#strip chr_, and assign AvZ
cands[, chr := gsub('chr_', '', chr)][
  , AvZ := fifelse(chr == 'W', 'W',
                   fifelse(chr == 'Z', 'Z',
                           fifelse(chr == 'MT', 'MT', 'Autosome')))]

#out of those, which are overrepresented genes 
length(unique(cands$Egg))
nonsyn_cands = cands %>% filter(gene != '.' & effect != '.' & effect != 'synonymous') 
targets = nonsyn_cands %>% 
  filter(gene != '.' & effect != '.' & effect != 'synonymous') %>% 
  #first count the number of mutations for each gene, for each egg
  group_by(AvZ,Egg,gene) %>% 
  summarize(snps = n()) %>% ungroup %>% 
  #now, find the genes which have multiple eggs with dN 
  group_by(AvZ,gene) %>% 
  summarize(Group = paste0(Egg,collapse='__'), total_snps = sum(snps),Groups = n()) %>% ungroup %>% 
  filter(Groups >= 4) 
targets

#save those genes 
write.table(targets$gene,file='../../figures/Candidate_Genes_GWAS-2Contrast_2024APR2.txt',quote=F,sep='\t',row.names=F,col.names=F)

#plot those 
gene_candidates = nonsyn_cands %>% filter(gene %in% targets$gene) %>% 
  #make all genes uppercase, and replace some of the mtDNA genes 
  mutate(gene = toupper(gsub('ID=','',gene)),
         gene = gsub('COX1','CO1',gene),
         gene = gsub('NAD1','ND1',gene),
         gene = gsub('NAD2','ND2',gene),
         gene = gsub('NAD3','ND3',gene),
         gene = gsub('NAD4','ND4',gene),
         gene = gsub('NAD5','ND5',gene),
         gene = gsub('NAD6','ND6',gene)) %>% 
  #arrange genes according to the highest number of eggs present 
  select(chr,Egg,gene,AvZ) %>% distinct %>% 
  group_by(gene) %>% 
  mutate(total = n()) %>% ungroup %>% 
  arrange(desc(total)) %>% 
  mutate(gene = factor(gene, levels = unique(gene))) %>%  
  ggplot(aes(y=gene,fill=Egg))+
  geom_bar(position='stack',col='black')+
  facet_grid(AvZ~.,scales='free',space='free')+
  scale_x_continuous(breaks= pretty_breaks(n=3))+
  scale_fill_manual(values=eggcols$col,breaks=eggcols$Egg)+
  theme_bw(base_size=6)+
  theme(legend.position='top')
gene_candidates

pdf('../../figures/Fixed_dN-SNPS_GWAS_2024APR2.pdf',height=3.5,width=2)
gene_candidates
dev.off()

#what are those SNPs?
snp_sites = nonsyn_cands %>% filter(gene %in% targets$gene) %>% arrange(chr,start,gene)
write.table(snp_sites,file='../../figures/GWAS_SNP-Sites_2024APR4.txt',quote=F,sep='\t',row.names=F)

##### Plot GWAS #####
library(karyoploteR)
#prep karyoplot
genome = read.table('/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa.bed',header=FALSE) %>% filter(str_detect(V1,'scaff',negate=T)) %>% arrange(desc(V2))
names(genome) = c('chr','start','end'); genome$chr = gsub('chr_','',genome$chr); 
genome = genome %>% mutate(chr = factor(chr, levels = c(as.character(1:39), "Z", "W", "MT"))) %>%  arrange(chr)
G = makeGRangesFromDataFrame(genome,seqnames.field = 'chr',start.field = 'start',end.field = 'end')

#add colors 
chrs = genome %>% select(chr) %>% mutate(Color = if_else(row_number() %% 2 == 0, 'grey60', 'black'))
candsgr = makeGRangesFromDataFrame(cands,keep.extra.columns = TRUE,seqnames.field = 'chr',start.field = 'start',end.field = 'end')

# General karyoplot parameters
pdf('../../figures/KARYOPLOT_GWAS_2024APR3.pdf',height=5,width=7.5)
pp = getDefaultPlotParams(plot.type=4)
pp$leftmargin = 0.1
kp = plotKaryotype(plot.type=4, genome = G,
                   labels.plotter = NULL,
                   plot.params = pp)
kpAddChromosomeNames(kp, yoffset = -5,cex=0.7)
kpAddBaseNumbers(kp,tick.dist = 50000000,minor.ticks = FALSE,cex=0.4)

#Loop through each track (all species together, CC, and CO), and plot them on their own layer 
counter = 0; tracks = 12
for (egg in eggcols$Egg) { 
  col = eggcols %>% filter(Egg == egg) %>% pull(col)
  counter = counter + 1; at = autotrack(current.track = counter, total.tracks = tracks);at$r1 <- at$r1-0.02
  subgr = candsgr[mcols(candsgr)$Egg == egg]
  kpPlotDensity(kp,data=subgr,r0=at$r0,r1=at$r1,col=col)
  kpAxis(kp,r0=at$r0,r1=at$r1,cex=0.3,numticks = 2)
  kpAddLabels(kp,cex=0.25,labels = egg,r0=at$r0+.01, r1=at$r1,col="black",srt=0,label.margin = 0.02)
  
}

dev.off()
```

## Contrast FST

Hierarchical contrasts:

```R
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(sf)
library(ggspatial)
md = read_tsv('~/merondun/cuculus_host/Metadata_Host.txt')
md = md %>% drop_na(Egg) %>% filter(Sex == 'F' & Analysis_PopulationGenetics == 1) %>% mutate(Group = paste0(Egg,'_',Hap))
keep = md %>% count(Group) %>% filter(n >= 4)
mds = md %>% filter(Group %in% keep$Group)
mds %>% count(Group)

#save metadata, and a file .list for each pop listing the individuals 
for (group in unique(mds$Group)){
  d = mds %>% filter(Group == group)
  write.table(d$ID,file=paste0('manyhost_hunt/males/pops/',group,'.list'),quote=F,sep='\t',row.names=F,col.names=F)
}
groups = mds %>% select(Group) %>% unique %>% pull(Group)
pairwise_combinations <- combn(groups, 2)

# Convert the combinations into a dataframe
pairwise_combinations_df <- data.frame(
  p1 = pairwise_combinations[1,],
  p2 = pairwise_combinations[2,]
)

#and add a phylo category
phylo = pairwise_combinations_df %>% mutate(Group = paste0(p1,'__',p2)) %>% 
  separate(p1,into=c('E1','H1'),remove=F) %>% 
  separate(p2,into=c('E2','H2'),remove=F) %>% 
  mutate(
    Depth = ifelse(grepl('W1|W2|W3',p1) & !grepl('W1|W2|W3',p2),'Ancient',
                   ifelse(grepl('W1|W2|W3',p2) & !grepl('W1|W2|W3',p1),'Ancient','Contemporary')),
    Comparison = ifelse(grepl('W1|W2',Group) & (grepl('E1_',p1) | grepl('E1_',p2)) & Depth == 'Ancient','Blue',
                        ifelse(Depth == 'Ancient','Reversion',
                               ifelse(Depth == 'Contemporary' & E1 == E2,'Control',
                                      ifelse((grepl('W1|W2|W3',p1) & grepl('W1|W2|W3',p2)) & Depth == 'Contemporary','Reversion','Diversification')))),
    Phylo = paste0(Depth,'_',Comparison)
  ) %>% arrange(Phylo) %>% 
  filter(Phylo != 'Ancient_Reversion')
phylo 

write.table(phylo,'manyhost_hunt/males/Contrast_Metadata.txt',quote=F,sep='\t',row.names=F)

dfs = phylo %>% mutate(Group = paste0(p1,'@',p2))
write.table(dfs$Group,file='manyhost_hunt/males/Pairwise_Contrasts.list',quote=F,sep='\t',row.names=F,col.names=F)

#jitter points up to 1 lat/long for viewing
mds = mds %>% mutate(LatJit = jitter(Latitude,amount =2),
                     LonJit = jitter(Longitude,amount=2))
sites = st_as_sf(mds, coords = c("LonJit", "LatJit"), 
                 crs = 4326, agr = "constant")

#set up map and convert df to coordinate frame
world = map_data("world")
cols = brewer.pal(4,'Paired')

imm_spat = ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), col='grey90', fill='white') +
  geom_sf(data = sites, 
          aes(fill=as.factor(Group)),
          size=3,show.legend = T,pch=21) +
  xlab('')+ylab('')+
  coord_sf(xlim = c(min(md$Longitude)-5, max(md$Longitude)+5), 
           ylim = c(min(md$Latitude)-5, max(md$Latitude)+5), expand = FALSE)+
  scale_fill_viridis(discrete=TRUE,option='turbo')+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"))+
  theme(legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))+
  annotation_scale(line_width=0.5)
imm_spat

pdf('figures/Contrasts_Haplotype-Females_SpatialDistribution_2024APR4.pdf',height=3,width=7)
imm_spat
dev.off()
```

![image-20240404092834008](C:\Users\herit\AppData\Roaming\Typora\typora-user-images\image-20240404092834008.png)

### Estimate FST

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --mem=5000mb
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00

#mamba activate snps
# for CAT in $(cat Pairwise_Contrasts.list); do sbatch -J FST_${CAT} FST.sh ${CAT}; done

CAT=$1

p1=$(echo ${CAT} | sed 's/@.*//g')
p2=$(echo ${CAT} | sed 's/.*@//g')

mkdir fst fst/work fst/out

for CHR in $(cat Chromosomes.list); do

if [[ $CHR = 'chr_W' || $CHR = 'chr_MT' || $CHR = 'chr_Z' ]]; then

        echo "WORKING ON FEMALE HAPLOID"
        vcftools --haploid --gzvcf female_vcfs/${CHR}.FEM.DP3.vcf.gz --out fst/work/${CHR}_${p1}_${p2} --weir-fst-pop pops/${p1}.list --weir-fst-pop pops/${p2}.list --fst-window-size 1 --max-missing 0.1

else

        echo "WORKING ON FEMALE DIPLOID"
        vcftools --gzvcf female_vcfs/${CHR}.FEM.DP3.vcf.gz --out fst/work/${CHR}_${p1}_${p2} --weir-fst-pop pops/${p1}.list --weir-fst-pop pops/${p2}.list --fst-window-size 1 --max-missing 0.1

fi

    awk -v p1=${p1} -v p2=${p2} '{OFS="\t"}{print $1, $2, $2, p1, p2, $5}' fst/work/${CHR}_${p1}_${p2}.windowed.weir.fst | \
            sed '1d' | bedtools intersect -a - -b Gene_Lookup.bed -wao  | \
            bedtools intersect -a - -b female_vcfs/Annotated_Variants_2024APR2.bed -wao | \
        awk '{OFS="\t"}{print $1, $2, $4, $5,$6,$11,$16,$17}' | \
        sed 's/ID=//g' | sed 's/gene-//g' > fst/out/${CHR}_${p1}_${p2}.fst.txt
awk '$5 == 1' fst/out/${CHR}_${p1}_${p2}.fst.txt > fst/out/${CHR}_${p1}_${p2}.fst1.txt

done

```

Identify candidates

```R
#### Identify candidates from female-only FST scan 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/manyhost_hunt/males')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(viridis)
library(ggpubr)
library(meRo)
library(RColorBrewer)
library(zoo)
library(ggpubr)
library(scales)
library(data.table)

md = read_tsv('~/merondun/cuculus_host/Metadata_Host.txt')
eggcols = md %>% ungroup %>% select(Egg) %>% unique %>% mutate(ord = as.numeric(gsub('E','',Egg))) %>% arrange(ord) %>% 
  drop_na(Egg) %>% mutate(col = viridis(12,option='turbo')[1:11])
eggcols = rbind(eggcols,data.frame(Egg='E6W3',ord=5,col='white')) %>% arrange(ord)

##### Load and prep FST BP data  ######
tidy_bp = fread('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/manyhost_hunt/males/FST_Contrasts_2024APR4.txt.gz')
names(tidy_bp) = c('chr','start','p1','p2','FST','intgene','effect','gene')

#Assign AvZ
tidy_bp = tidy_bp %>% 
  mutate(
    chr = gsub('chr_','',chr),
    AvZ = case_when(
      chr == 'W' ~ 'W',
      chr == 'Z' ~ 'Z',
      chr == 'MT' ~ 'MT',
      TRUE ~ 'Autosome'),
    #add floor and ceiling to FST between 0 - 1
    FST = pmax(0, pmin(1, FST)),
    site = paste0(chr,'_',start)
  ) 

phylo = read_tsv('Contrast_Metadata.txt') %>% select(p1,p2,Group,Phylo)
tidy_dat = left_join(tidy_bp,phylo) %>% mutate(gene = toupper(gene),
                                               intgene = toupper(intgene)) 

#grab only the fixed sites!
tidy_fixed = tidy_dat %>% filter(FST == 1) 
write.table(tidy_fixed,file='FST1_TIDY_2024APR4.txt',quote=F,sep='\t',row.names=F)
tidy_fixed = read_tsv('FST1_TIDY_2024APR4.txt')

fwrite(tidy_dat,'FST_TIDY_2024APR4.txt.gz',scipen = getOption('scipen', 999L),sep='\t')
tidy_dat = fread('FST_TIDY_2024APR4.txt.gz')

###### Identify Candidate SNPS ###### 
#these comparisons do not differ by egg, so remove those sites
tidy_control = tidy_fixed %>% filter(Phylo == 'Contemporary_Control')

#remove those
tidy_fixed_targ = tidy_fixed %>%
  filter(Phylo != 'Contemporary_Control') %>% #13 comparisons left
  filter(!site %in% tidy_control$site)

#count the number of contrasts, for now, since we may have duplicate sites (overlapping CDS) assigned, ensure that we only have single
tidy_prop_input = tidy_fixed_targ %>% select(chr,start,AvZ,p1,p2,Group,Phylo,FST) %>% distinct
tidy_prop_input %>% select(Group,Phylo) %>% distinct %>% count(Phylo)
tidy_prop_input %>% count(Phylo,Group,chr,start) %>% filter(n > 1) #this should be zero

#requiring at least one contrast...
tidy_prop_input_filtered = tidy_prop_input %>% group_by(Phylo,AvZ,chr,start) %>% count %>%
  filter(
    ifelse(Phylo == 'Ancient_Blue', n >= 3, 
           ifelse(Phylo == 'Contemporary_Diversification', n >= 1,
                  ifelse(Phylo == 'Contemporary_Reversion', n >= 1, n >= 999))))

#requiring all contrasts for the same SNP (used in primary figure!) 
tidy_prop_input_filtered = tidy_prop_input %>% group_by(Phylo,AvZ,chr,start) %>% count %>%
  filter(
    ifelse(Phylo == 'Ancient_Blue', n >= 6, 
           ifelse(Phylo == 'Contemporary_Diversification', n >= 2,
                  ifelse(Phylo == 'Contemporary_Reversion', n >= 2, n >= 999))))

##### Proportions of FST=1 SNPs #####
#total windows, count the proportion of sites which are FST = 1.0 in ALL contrasts for a certain mutation level
chr_windows = tidy_dat %>% select(chr,start,AvZ,Group,Phylo) %>% distinct %>% count(AvZ,Group,Phylo) %>% ungroup %>% group_by(Phylo,AvZ) %>% slice_max(n) %>% select(Phylo,AvZ,NumberWindows = n) %>% unique %>% filter(!grepl('Control',Phylo))
tidy_fixed_prop = tidy_prop_input_filtered  %>% ungroup %>% 
  count(AvZ,Phylo) %>% distinct %>%  #count number of SNPs by phylo 
  left_join(chr_windows,.) %>% ungroup %>% 
  replace_na(list(n=0)) %>% 
  ungroup %>%
  mutate(Proportion = n/NumberWindows)
tidy_fixed_prop$AvZ = factor(tidy_fixed_prop$AvZ,levels=c('Autosome','Z','W','MT'))
tidy_fixed_prop$Phylo = factor(tidy_fixed_prop$Phylo,levels=c('Ancient_Blue','Contemporary_Reversion','Contemporary_Diversification'))

#check colorblind
# library(colorblindcheck)
# pal <- c('cadetblue3', 'plum', 'indianred1')
# pal_hex <- col2rgb(pal) / 255
# pal_hex <- apply(pal_hex, 2, function(x) rgb(x[1], x[2], x[3]))
# show_col(pal_hex)
# palette_check(pal_hex,plot=TRUE)

#plot, divide into top and bottom plots for visualization 
props = tidy_fixed_prop %>% ggplot(aes(y=AvZ,x=Proportion,fill=Phylo,label=n))+
  geom_bar(stat='identity',position=position_dodge(width=0.9))+
  geom_text(position=position_dodge(width=0.9),hjust=-.25,size=3)+
  scale_fill_manual(values=rev(brewer.pal(3,'Set2')))+
  #coord_cartesian(xlim=c(0,1.0))+
  scale_x_continuous(n.breaks = 4)+
  #facet_grid(AvZ~.,scales='free',nrow=4,ncol=1)+
  theme_bw(base_size=6)+xlab('Proportion of FST = 1.0 SNPs Present in At Least 1 Comparison')+ylab('')+
  theme(legend.position = "top", # Moves legend to bottom
        legend.text = element_text(size = 6), # Adjusts text size
        legend.title = element_text(size = 6)) + # Adjusts title size
  guides(fill = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 1)))
props

pdf('../../figures/20240718_FST-BP_Females_Proportions-Half.pdf',height=2,width=4)
props
dev.off()

#are there regions of MT/W which contain more differentiated windows?
dens = tidy_prop_input_filtered %>% filter((chr == 'W' | chr == 'MT')) %>%
  ggplot(aes(x=start,fill=Phylo))+
  scale_fill_manual(values=c('cadetblue3','goldenrod1','indianred1'))+
  geom_density(alpha=0.8)+facet_wrap(.~chr,scales='free',nrow=2,ncol=1)+theme_bw(base_size=6) +theme(legend.position='top')
dens

pdf('../../figures/FST-BP_Females_Density_2024APR4.pdf',height=4,width=2.5)
dens
dev.off()

##### dN and dS SNPs #####

#these are our strongest candidate SNPs  
candidates = left_join(tidy_prop_input_filtered %>% select(-n),tidy_fixed_targ)

#replace LOC IDs with known homologs
locs = read_tsv('../LOC_Lookup.txt')

#based only on overlapping genes 
cands = left_join(candidates,locs %>% dplyr::rename(intgene=Label)) %>% 
  mutate(gene = ifelse(is.na(ID),intgene,paste0(ID,'-L'))) %>% select(-ID) %>% unique %>% 
  mutate(gene = gsub('NAD','ND',gene), gene = gsub('_A','',gene)) 
cand_genes = cands %>% filter(gene != '.') %>% ungroup %>% count(gene,Phylo,AvZ,chr) 
all_targets = cand_genes %>% ungroup %>% group_by(AvZ,gene) %>% 
  summarize(Group = paste0(Phylo,collapse='__'), total_snps = sum(n),Groups = n()) %>% ungroup 
targets = all_targets %>% filter(Groups >= 1)

#based on nonsynonymous/nonsense (shown in primary figure!) 
cands = left_join(candidates,locs %>% dplyr::rename(gene=Label)) %>% 
  mutate(gene = ifelse(is.na(ID),gene,paste0(ID,'-L'))) %>% distinct %>% 
  mutate(gene = gsub('NAD','ND',gene), gene = gsub('_A','',gene)) 
cand_genes = cands %>% filter(effect != '.' & effect != 'synonymous') %>% ungroup %>% count(gene,Phylo,AvZ,chr) 
cands %>% filter(effect == 'nonsense') 

#plot how many dN SNPs are in each gene, across each phylo category 
cand_genes$Phylo = factor(cand_genes$Phylo,levels=c('Ancient_Blue', 'Contemporary_Reversion','Contemporary_Diversification'))

#Vertical 
p = cand_genes %>% 
  filter(gene %in% targets$gene) %>%   #hash out for dN 
  select(Phylo,gene,AvZ) %>% distinct %>% 
  group_by(gene) %>% 
  mutate(total = n()) %>% ungroup %>% 
  arrange(desc(total)) %>% 
  mutate(gene = factor(gene, levels = unique(gene))) %>%  
  ggplot(aes(y=gene,fill=Phylo))+
  geom_bar(position='stack')+
  facet_grid(AvZ~.,scales='free',space='free')+
  scale_x_continuous(breaks= pretty_breaks(n=3))+
  #scale_shape_manual(values=rev(c(16,15,18,17)))+
  scale_fill_manual(values=rev(brewer.pal(3,'Set2')))+
  theme_bw(base_size=6)+
  theme(legend.position='top')
p

# Horizontal  
cand_genes$AvZ <- factor(cand_genes$AvZ,levels=c('MT','W','Z','Autosome'))
p = cand_genes %>% 
  filter(gene %in% targets$gene) %>%   #hash out for dN 
  select(Phylo,gene,AvZ) %>% distinct %>% 
  group_by(gene) %>% 
  mutate(total = n()) %>% ungroup %>% 
  arrange(desc(total)) %>% 
  mutate(gene = factor(gene, levels = unique(gene))) %>%  
  ggplot(aes(x=gene,fill=Phylo))+
  geom_bar(position='stack')+
  facet_grid(.~AvZ,scales='free',space='free')+
  scale_y_continuous(breaks= pretty_breaks(n=3))+
  scale_fill_manual(values=rev(brewer.pal(3,'Set2')))+
  theme_bw(base_size=6)+xlab('')+ylab('')+
  theme(legend.position='top',
        axis.text.x=element_text(angle=45,hjust=1,vjust=1))
p

pdf('../../figures/20240722_FST-BP_Females_dN_Mutations_Genes_Horizontal.pdf',height=2.5,width=6)
p
dev.off()

#save those genes
genes = targets %>% mutate(gene = gsub('-L','',gene))
write.table(genes$gene,file='../../figures/FST-BP_Females_CandidatesGenesOverlap_2024APR4.txt',quote=F,sep='\t',row.names=F,col.names=F)

#save nonsynonymous genes
dn_genes = cand_genes %>% mutate(gene = gsub('-L','',gene))
write.table(dn_genes$gene,file='../../figures/FST-BP_Females_CandidatesGenesdN_2024APR4.txt',quote=F,sep='\t',row.names=F,col.names=F)

#grab the actual SNPs from the main positions 
inspection = cands %>% filter(effect == 'nonsynonymous' & gene %in% cand_genes$gene) %>% ungroup %>% select(chr,start,gene,Phylo) %>% unique %>% group_by(chr,gene,Phylo) %>% 
  #slice_max(start) %>% 
  data.frame
inspect_bed = inspection %>% filter(grepl('ND2|ND5|ND4|TLDC2|NDUFAF4',gene)) %>% mutate(end = start,chr=paste0('chr_',chr)) %>% select(chr,start,end,Phylo,gene)
write.table(inspect_bed,file='Fixed_SNP_Candidates_PositionsForGenotypePlot_2024JULY18.txt',quote=F,sep='\t',row.names=F,col.names=F)

#in bash:
for chr in $(awk '{print $1}' Fixed_SNP_Candidates_PositionsForGenotypePlot_2024JULY18.txt | sort | uniq); do
echo "WORKING ON CHR: ${chr}"
bedtools intersect -header -a female_vcfs/${chr}.FEM.DP3.vcf.gz -b Fixed_SNP_Candidates_PositionsForGenotypePlot_2024JULY18.txt | bcftools view -Oz -o plot_genos/${chr}.vcf.gz
done
bcftools concat -Oz -o plot_genos/Targets.vcf.gz plot_genos/chr_*.vcf.gz

##### Plot Genotypes #####
library(vcfR)
library(tidyverse)
library(data.table)
library(RColorBrewer)
#plot genotypes 
vcfr = read.vcfR('plot_genos/Targets.vcf.gz')
gt_data <- vcfR2tidy(vcfr, format_fields = 'GT')
gt = gt_data$gt
#remember vcfs are 1-based
names(gt) = c('Key','site','ID','Genotype','Allele')

#add metadata since I want to group by morph
md = read_tsv('~/merondun/cuculus_host/Metadata_Host.txt')
md = md %>% drop_na(Egg)
gtm = left_join(md %>% dplyr::select(ID,Hap,Egg),gt) #ensure that the 'ID' from the vcf matches the 'ID' from your metadata

#now join with the annotated snps
agt = gtm %>% drop_na(site) %>% replace_na(list(Genotype = 'Missing'))
agt %>% count(Key,site)
aga = left_join(agt,inspection %>% dplyr::rename(site=start)) %>% drop_na(gene) %>% mutate(gene = toupper(gsub('ID=','',gene)))
pg = aga %>% arrange(chr,site)
sites = pg %>% select(chr,site) %>% unique
pg$site = factor(pg$site,levels=sites$site)

samples = read_tsv('Female_Hunt.list',col_names=F)
names(samples) = 'ID'

pg = pg %>% filter(ID %in% samples$ID)
pg$Phylo = factor(pg$Phylo,levels=c('Ancient_Blue','Contemporary_Reversion','Contemporary_Diversification'))
pg$gene = factor(pg$gene,levels=c('NDUFAF4-L','ND2','ND4','ND5','TLDC2'))
pg$Egg = factor(pg$Egg,levels=eggcols$Egg)
gtp = pg %>% 
  group_by(chr,gene,site,Genotype,Allele,ID,Hap,Egg) %>% 
  mutate(Genotype = ifelse(Genotype == '0','0/0',
                           ifelse(Genotype == '1','1/1',Genotype))) %>% 
  #summarize(Group = paste0(Phylo,collapse='__')) %>% 
  ggplot(aes(y = ID, x = as.factor(site), fill = Genotype)) +
  geom_tile() +
  ylab("Individuals") +
  facet_grid(Hap+Egg~gene+chr+Phylo,space='free',scales='free')+
  scale_fill_manual(values=c('grey70','grey40','grey10','white'))+
  #scale_fill_manual(values=c('grey70','grey20','red'))+
  theme_minimal(base_size=6) + ylab('')+xlab('') +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.x = element_text(size=4,angle = 90, vjust = 0.5, hjust = 1),
        strip.text.y.right = element_text(angle = 0),legend.position='top',
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
gtp

pdf('../../figures/Genotypes-Targets_2024APR4.pdf',height=3.5,width=1.75)
gtp
dev.off()
```

## Blast NDUFAF4 other Cuculiformes

### Blast-Based

Only one other Cuculiformes species with a W assembly: [Phaenicophaeus curvirostris](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_032191515.2/). 

Download the genome, prep for blast:

```bash
# Make blast db
makeblastdb -in GCA_032191515.2_BPBGC_Pcur_1.0_genomic.fna -parse_seqids -dbtype nucl

#gene to blast
ndufaf4=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/cuculiformes_NDUFAF4_blasting/NDUFAF4.fa

#also blast against cuckoo 
cuckoo_blastdb=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/blast_db/GCA_017976375.1_bCucCan1.pri_genomic.CHR

#first blast against  cuckoo
blastn -query $ndufaf4 -db $cuckoo_blastdb -evalue 1e-6 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand' -num_threads 5
```

As expected, with common cuckoo it blasts to chr3 and chrW:

| qseqid                         | sseqid | pident | length | mismatch | gapopen | qstart | qend | sstart   | send     | evalue   | bitscore | sstrand |
| ------------------------------ | ------ | ------ | ------ | -------- | ------- | ------ | ---- | -------- | -------- | -------- | -------- | ------- |
| NC_071403.1:c25515677-25511697 | chr_3  | 100    | 3981   | 0        | 0       | 1      | 3981 | 25515677 | 25511697 | 0        | 7352     | minus   |
| NC_071403.1:c25515677-25511697 | chr_3  | 85.891 | 645    | 59       | 6       | 3032   | 3675 | 25516398 | 25515785 | 0        | 658      | minus   |
| NC_071403.1:c25515677-25511697 | chr_W  | 95.745 | 141    | 4        | 2       | 3285   | 3424 | 21177317 | 21177456 | 3.25E-56 | 226      | plus    |

And against the chestnut-breasted malkoha:

```bash
phaen_blastdb=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/cuculiformes_NDUFAF4_blasting/Phaenicophaeus_curvirostris/GCA_032191515.2_BPBGC_Pcur_1.0_genomic.fna

blastn -query $ndufaf4 -db $phaen_blastdb -evalue 1e-6 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand' -num_threads 5
```

| qseqid                         | sseqid           | pident | length | mismatch | gapopen | qstart | qend | sstart  | send    | evalue   | bitscore | sstrand |
| ------------------------------ | ---------------- | ------ | ------ | -------- | ------- | ------ | ---- | ------- | ------- | -------- | -------- | ------- |
| NC_071403.1:c25515677-25511697 | gb\|CM063363.1\| | 87.021 | 2612   | 244      | 40      | 1396   | 3980 | 5656863 | 5659406 | 0        | 2857     | plus    |
| NC_071403.1:c25515677-25511697 | gb\|CM063363.1\| | 87.5   | 336    | 33       | 5       | 3032   | 3365 | 5653399 | 5653727 | ######## | 379      | plus    |
| NC_071403.1:c25515677-25511697 | gb\|CM063363.1\| | 84.737 | 190    | 22       | 7       | 1      | 186  | 5654507 | 5654693 | 2.21E-43 | 183      | plus    |
| NC_071403.1:c25515677-25511697 | gb\|CM063363.1\| | 83.871 | 124    | 20       | 0       | 3558   | 3681 | 5654302 | 5654425 | 6.32E-24 | 119      | plus    |

Interesting, how about if we go deeper into the phylogeny:

| Species                                                      | Common               | Accession       | Link                                                         |
| ------------------------------------------------------------ | -------------------- | --------------- | ------------------------------------------------------------ |
| [Clamator_glandarius](https://www.ncbi.nlm.nih.gov/datasets/taxonomy/78203) | Great_Spotted_Cuckoo | GCA_033459285.1 | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/033/459/285/GCA_033459285.1_BPBGC_Cgla_1.0/GCA_033459285.1_BPBGC_Cgla_1.0_genomic.fna.gz |
| [Ceuthmochares_aereus](https://www.ncbi.nlm.nih.gov/datasets/taxonomy/1961834) | Blue_Malkoha         | GCA_013398935.1 | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/013/398/935/GCA_013398935.1_ASM1339893v1/GCA_013398935.1_ASM1339893v1_genomic.fna.gz |
| Tapera_naevia                                                | Striped_Cuckoo       | GCA_033556795.1 | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/033/556/795/GCA_033556795.1_BPBGC_Tnae_1.0/GCA_033556795.1_BPBGC_Tnae_1.0_genomic.fna.gz |
| [Coccyzus_lansbergi](https://ebird.org/species/gyccuc)       | Grey_capped_Cuckoo   | GCA_033558105.1 | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/033/558/105/GCA_033558105.1_BPBGC_Clan_1.0/GCA_033558105.1_BPBGC_Clan_1.0_genomic.fna.gz |
| [Geococcyx_californianus](https://www.ncbi.nlm.nih.gov/datasets/taxonomy/8947) | Greater_Roadrunner   | GCA_013389885.1 | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/013/389/885/GCA_013389885.1_ASM1338988v1/GCA_013389885.1_ASM1338988v1_genomic.fna.gz |
| [Piaya_cayana](https://www.ncbi.nlm.nih.gov/datasets/taxonomy/33601) | Squirrel_Cuckoo      | GCA_013389865.1 | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/013/389/865/GCA_013389865.1_ASM1338986v1/GCA_013389865.1_ASM1338986v1_genomic.fna.gz |
| [Dromococcyx_pavoninus](https://www.ncbi.nlm.nih.gov/datasets/taxonomy/3054325) | Pavonine_Cuckoo      | GCA_033558185.1 | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/033/558/185/GCA_033558185.1_BPBGC_Dpav_1.0/GCA_033558185.1_BPBGC_Dpav_1.0_genomic.fna.gz |
| Cuculus_canorus                                              | Common_Cuckoo        | GCA_017976375.1 | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/017/976/375/GCF_017976375.1_bCucCan1.pri/GCF_017976375.1_bCucCan1.pri_genomic.fna.gz |

Loop it:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --time=48:00:00

#gene to blast
ndufaf4=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/cuculiformes_NDUFAF4_blasting/NDUFAF4.fa
results_dir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/cuculiformes_NDUFAF4_blasting/results
base_dir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/cuculiformes_NDUFAF4_blasting/
RUN=$1

#make directory
mkdir ${RUN}

#grab the genome link for wget
file_path=$(grep ${RUN} All_Species.txt | awk '{print $4}')

#go into directory and download file
cd ${RUN}
wget ${file_path}
gunzip *gz

#make blastdb
makeblastdb -in *fna -parse_seqids -dbtype nucl

#blast
blastn -query $ndufaf4 -db *fna -evalue 1e-6 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand' -num_threads 2 > $results_dir/${RUN}.blast.txt

#grab those sequences as well to make a tree
awk -v r=${RUN} '{
    OFS="\t";
    start=$9;
    end=$10;
    if (start > end) {
        print $2, end, start, r"_"NR, "0", "-";
    } else {
        print $2, start, end, r"_"NR, "0", "+";
    }
}' $results_dir/${RUN}.blast.txt | sed 's/gb|//g' | sed 's/|//g' | sed 's/ref//g' > ${RUN}.bed

# Merge any hits found within 3KB of one another
bedtools sort -i ${RUN}.bed | bedtools merge -i - -d 3000 -c 4,5,6 -o first,distinct,distinct > ${RUN}_merged.bed

bedtools getfasta -fi *fna -bed ${RUN}_merged.bed -fo $results_dir/${RUN}.fa -nameOnly

```

How many sequences per species:

```bash
seqkit stats *fa
processed files:  10 / 10 [======================================] ETA: 0s. done
file                        format  type  num_seqs  sum_len  min_len  avg_len  max_len
Clamator_glandarius.fa      FASTA   DNA          1    5,313    5,313    5,313    5,313
Cuculus_canorus.fa          FASTA   DNA          2    4,840      139    2,420    4,701
Dromococcyx_pavoninus.fa    FASTA   DNA          1    3,650    3,650    3,650    3,650
Geococcyx_californianus.fa  FASTA   DNA          1    3,760    3,760    3,760    3,760
Piaya_cayana.fa             FASTA   DNA          1    5,032    5,032    5,032    5,032
Tapera_naevia.fa            FASTA   DNA          1    2,528    2,528    2,528    2,528
```

### Fastq Alignment-Based

From the merged blast hits on the cuckoo genome, extract the sequence corresponding to the 2 NDUFAF4 hits for inspection:

```bash
cat Cuculus_canorus_merged.bed
NC_071403.1     25511697        25516398        Cuculus_canorus_1       0       -
NC_071440.1     21177317        21177456        Cuculus_canorus_3       0       +
```

Align all reads from each subsampled 5gb library to the whole cuckoo genome. This will take 4 runs, corresponding to illumina / HiFi for cuckoos vs others. 

Cuckoo HiFi data:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00

# Cuckoo HiFi genome individual

RUN=cuckoo_HiFi

echo "Aligning sample: ${RUN}"

SCRATCH=/tmp/$SLURM_JOB_ID
genome=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa
outdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/cuculiformes_NDUFAF4_blasting/Alignment_Based_Approach/Whole_Genome_Alignment/bams
cuckoo_hifi=/dss/dsslegfs01/pr53da/pr53da-dss-0021/rawdata/External_Data_SRA_GenomeArk/PacBio_GenomeArk/Cuculus_canorus/bCucCan1_pacbio.fastq.gz
subdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/cuculiformes_NDUFAF4_blasting/Alignment_Based_Approach/subset_5gb

# Subset 5gb, and then align, discarding unaligned reads to save space (-F 4)
bbduk.sh in=${cuckoo_hifi} out=${subdir}/${RUN}.5gb.fastq.gz maxbasesout=5000000000
minimap2 -ax map-pb -t 8 ${genome} ${subdir}/${RUN}.5gb.fastq.gz > ${SCRATCH}/${RUN}.sam
samtools sort ${SCRATCH}/${RUN}.sam | samtools view -F 4 -b > ${outdir}/${RUN}.bam
samtools index -b ${outdir}/${RUN}.bam;
```

Cuckoo illumina data:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00

# Align cuckoo illumina data

RUN=$1

echo "Aligning sample: ${RUN}"

SCRATCH=/tmp/$SLURM_JOB_ID
genome=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa
qcdata=/dss/dsslegfs01/pr53da/pr53da-dss-0021/shared_resources/trimmed_fastq
outdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/cuculiformes_NDUFAF4_blasting/Alignment_Based_Approach/Whole_Genome_Alignment/bams
subdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/cuculiformes_NDUFAF4_blasting/Alignment_Based_Approach/subset_5gb

# Subset 5gb, and then align, discarding unaligned reads to save space (-F 4)
bbduk.sh in=${qcdata}/${RUN}.trim.fastq.gz out=${subdir}/${RUN}.5gb.fastq.gz maxbasesout=5000000000
bwa mem -M -p -t 10 ${genome} ${subdir}/${RUN}.5gb.fastq.gz | samtools sort -@10 -o ${SCRATCH}/${RUN}.bam -;
samtools view -F 4 -b ${SCRATCH}/${RUN}.bam > ${outdir}/${RUN}.bam
samtools index -b ${outdir}/${RUN}.bam
```

Outgroup HiFi data:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00

# Align pacbio HiFi data for outgroups

RUN=$1

echo "Aligning sample: ${RUN}"

SCRATCH=/tmp/$SLURM_JOB_ID
genome=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa
outdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/cuculiformes_NDUFAF4_blasting/Alignment_Based_Approach/Whole_Genome_Alignment/bams
rawdata=/dss/dsslegfs01/pr53da/pr53da-dss-0021/rawdata/External_Data_SRA_GenomeArk/SRA_Cuculiformes_Outgroups
subdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/cuculiformes_NDUFAF4_blasting/Alignment_Based_Approach/subset_5gb

# Subset 5gb, and then align, discarding unaligned reads to save space (-F 4)
bbduk.sh in=${rawdata}/${RUN}_1.fastq.gz out=${subdir}/${RUN}.5gb.fastq.gz maxbasesout=5000000000
minimap2 -ax map-pb -t 8 ${genome} ${subdir}/${RUN}.5gb.fastq.gz > ${SCRATCH}/${RUN}.sam
samtools sort ${SCRATCH}/${RUN}.sam | samtools view -F 4 -b > ${outdir}/${RUN}.bam
samtools index -b ${outdir}/${RUN}.bam
```

And outgroup Illumina:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00

# Align other species, illumina data

RUN=$1

echo "Aligning sample: ${RUN}"

SCRATCH=/tmp/$SLURM_JOB_ID
genome=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa
rawdata=/dss/dsslegfs01/pr53da/pr53da-dss-0021/rawdata/External_Data_SRA_GenomeArk/SRA_Cuculiformes_Outgroups
outdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/cuculiformes_NDUFAF4_blasting/Alignment_Based_Approach/Whole_Genome_Alignment/bams
subdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/cuculiformes_NDUFAF4_blasting/Alignment_Based_Approach/subset_5gb

# Subset 5gb, and then align, discarding unaligned reads to save space (-F 4)
bbduk.sh in=${rawdata}/${RUN}_1.fastq.gz out=${subdir}/${RUN}.5gb.fastq.gz maxbasesout=5000000000
bwa mem -M -t 10 ${genome} ${subdir}/${RUN}.5gb.fastq.gz | samtools sort -@10 -o ${SCRATCH}/${RUN}.bam -;
samtools view -F 4 -b ${SCRATCH}/${RUN}.bam > ${outdir}/${RUN}.bam
samtools index -b ${outdir}/${RUN}.bam;
```

Calculate coverage:

```bash
for RUN in $(ls *bam | sed 's/.bam//g'); do 

cds=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/cuculiformes_NDUFAF4_blasting/Alignment_Based_Approach/CDS_Regions_N3_BothParalogs.bed

# Calculate coverage MQ60, --by 50 reports coverage in 50bp regions 
mosdepth --threads 3 --mapq 60 --by ${cds} --fast-mode --no-per-base ../coverage/${RUN}_cds ${RUN}.bam

done 

# After, in /coverage/, merge into single file:
for i in $(ls *_cds.regions.bed | sed 's/_cds.regions.bed//g'); do awk -v i=${i} '{OFS="\t"}{print $1, $2, $3, $4, i}' ${i}_cds.regions.bed > ${i}.cds.cov ; done
cat *.cds.cov > NDUFAF4_MQ60_Coverage_2024JULY25_CDS.txt
```

### Extract Consensus

```bash
awk '{OFS="\t"}{print $1, $4, $5, $1"_"$4"_"$5, ".", $7}' NDUFAF4_Cuckoo_BothParalogs.gff > WholeGene_Regions.bed
```

Extract the fasta for that region from all samples:

```bash
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
```

Ensuring that we only trim and create a tree for the species which have coverage for chrW / chr3: 

```bash
cat *fa > All_Except_064.fa

# Align and keep the coordinates of the chr3 autosome 
mafft --thread 10 --auto --addfragments All_Except_064.fa --keeplength --reorder 064_CC_GRW_BGR_M__SRR11531726.fasta > NDUFAF4_chr3anchor.fa

# Align freely
mafft --thread 10 --auto All.fa > NDUFAF4_freealign.fa

# With python, split the MSA
from Bio import SeqIO

# Path to your MSA file in FASTA format
msa_file = 'NDUFAF4_chr3anchor.fa'

# Read the MSA file and create individual FASTA files
with open(msa_file, 'r') as msa:
    for record in SeqIO.parse(msa, 'fasta'):
        filename = f"{record.id}.fa"
        with open(filename, 'w') as output_file:
            SeqIO.write(record, output_file, 'fasta')

# And then identify gaps 
seqkit stats -a *

# Identify the samples with < 50% gaps, and create a tree:
for i in $(cat Keep.list); do cat ${i} >> NDUFAF4_Gapless50.fa; done

# Tree
iqtree --redo -keep-ident -T 20 -s NDUFAF4_Gapless50.fa --seqtype DNA -m "GTR" -B 1000
```

Plot Coverage and Tree:

```R
#### Plot NDUFAF4 Coverage 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/cuculiformes_NDUFAF4_blasting/Alignment_Based_Approach/Whole_Genome_Alignment/coverage/')
.libPaths('~/mambaforge/envs/R/lib/R/library')
library(tidyverse)
library(viridis)

# Load in coverage data 
c <- read_tsv('NDUFAF4_MQ60_Coverage_2024JULY25_CDS_WGA.txt')
c$Species <- factor(c$Species,levels=c('Cuculus_canorus','Cuculus_micropterus','Cuculus_poliocephalus','Clamator_glandarius','Coccyzus_lansbergi','Piaya_cayana','Dromococcyx_pavoninus','Tapera_naevia','Geococcyx_californianus'))
c <- c %>% arrange(Species,Coverage)
c$ID <- factor(c$ID,levels=unique(c$ID))


c %>% 
  group_by(Species,Sex,chr,Object,start,end) %>% 
  summarize(Coverage = mean(Coverage)) %>% 
  filter(Species == 'Cuculus_canorus' & Sex == 'Female') %>% 
  ggplot(aes(y=Object,x=Coverage,col=chr))+
  geom_point()+
  theme_bw()

cp <- c %>% 
  filter(Object == 'Gene') %>% 
  group_by(Species,Sex,chr) %>% 
  summarize(Coverage = mean(Coverage)) %>% 
  ggplot(aes(y=Species,x=Coverage,fill=chr,shape=chr))+
  geom_point(size=2,position=position_jitter(height=0.15))+
  facet_grid(.~Sex,scales='free')+
  scale_shape_manual(values=c(21,24))+
  theme_bw(base_size=6)+
  theme(legend.position='top')
cp

pdf('../../../../figures/20240807_NDUFAF4_Coverage.pdf',height=2.75,width=2.25)
cp
dev.off()

# All objects
allcds <- c %>% 
  group_by(Species,Object,Sex,chr) %>% 
  summarize(Coverage = mean(Coverage)) %>% 
  ggplot(aes(y=Species,x=Coverage,fill=chr))+
  geom_point(size=3,pch=21,position=position_jitter(height=0.15))+
  facet_grid(Object~Sex,scales='free')+
  theme_bw()

png('NDUFAF4_Coverage-IncludingCDS.png',units='in',res=300,height=3.5,width=7)
cp
dev.off()

# Plot Tree afterwards
library(ggtree)
library(treeio)

# Free alignment
iqtree = read.iqtree('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/cuculiformes_NDUFAF4_blasting/Alignment_Based_Approach/Whole_Genome_Alignment/bams/gene_alignment/free_align/NDUFAF4_Gapless50.fa.contree')

iqtr = root(as.phylo(iqtree),'Geococcyx_californianus_SRR9994302_chr_3_25511697_25515677')

md <- read_tsv('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/cuculiformes_NDUFAF4_blasting/Alignment_Based_Approach/Whole_Genome_Alignment/bams/Metadata.txt')
g <- ggtree(iqtr,layout='rectangular')
#g <- ggtree(iqtr,layout='ape')
g$data <- g$data %>% mutate(Chromosome = ifelse(grepl('chr_W',label),'chrW','chr3'),
                            label = gsub('_chr_W.*','',label),
                            label = gsub('_chr_3.*','',label))
g$data <- left_join(g$data,md %>% dplyr::rename(label = ID))

#plot with outgroups 
gtree <- g +
  geom_tippoint(aes(col=Chromosome,shape=Species),size=1,stroke=1)+
  geom_nodelab(aes(label=node),geom = 'text',size=1.5)+ 
  scale_shape_manual(values=c(1,2,3,4,5,6,7,8,9))+
  xlim(c(0,1))+
  #geom_tiplab(aes(label=Species),size=2,offset = 0.05)+
  geom_tiplab(size=2,offset = 0.05)+
  theme(legend.position='none')
gtree

pdf('../../../../figures/20240806_NDUFAF4-ParalogTree.pdf',height=2.75,width=2.5)
gtree
dev.off()

```





# Sensitivities

## PAR ID

```R
##### CNV #####
#what's happening with log2FM in the PAR area? 
cnv = read_tsv('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/manyhost_hunt/coverage/CNV_log2_EggSex-ZW_500BP_2024MAR19.txt')
cnv = cnv %>% mutate(chr = gsub('chr_','',chr),
                     #chr = ifelse(chr == 'Z' & start > 77e6,'W',chr),
                     #start = ifelse(chr == 'W' & start > 77e6,start-77e6+22312972,start)
                     ) %>% 
  drop_na(log2)

#assign phylo to each egg and sex comparison 
egg_cnv = cnv %>% filter(comparison == 'Egg') %>% mutate(Group = paste0(p1,'__',p2)) %>% drop_na(log2)
egg_cnv = egg_cnv %>% mutate(Depth = ifelse(grepl('W1|W2|W3',p1) & !grepl('W1|W2|W3',p2),'Ancient',
                                            ifelse(grepl('W1|W2|W3',p2) & !grepl('W1|W2|W3',p1),'Ancient','Contemporary')),
                             Comparison = ifelse(grepl('W1|W2',Group) & Depth == 'Ancient','Blue',
                                                 ifelse(Depth == 'Ancient','Reversion',
                                                        ifelse(Group == 'W2__E1__W1__E1' | Group == 'W7__E6__W5__E6','Control',
                                                               ifelse((grepl('W1|W2|W3',p1) & grepl('W1|W2|W3',p2)) & Depth == 'Contemporary','Reversion','Diversification')))),
                             Phylo = paste0(Depth,'_',Comparison))
sex_cnv = cnv %>% filter(comparison != 'Egg') %>% mutate(Group = comparison) %>% drop_na(log2) %>% 
  mutate(Phylo = Group)

#region elevated in blue:
targ_start = 7778296-5000
targ_end = 7783688+5000

#region between E10 
targ_start = 7739546-5000
targ_end = 7739546+5000

#region FECH
targ_start = 77714931-5e4
targ_end = 77729710+5e4
PAR = sex_cnv %>% filter(chr == 'Z') %>% group_by(chr,start) %>% sum_stats(log2) 
parID = PAR %>% filter(conf_low >= -0.5 & conf_high <= 0.5) %>% ggplot(aes(x=start))+geom_histogram()+theme_bw()
PAR_regions = PAR %>% mutate(PAR = ifelse(conf_low >= -0.5 & conf_high <= 0.5 & start >= 77294000,'PAR','ZW'))

PAR_boundary = PAR_regions %>%
  group_by(chr) %>%
  mutate(consecutive = with(rle(PAR), rep(seq_along(lengths), lengths))) %>%
  group_by(chr, consecutive) %>%
  filter(PAR == "PAR" & n() >= 5) %>%
  summarise(start_range = min(start), end_range = max(start)) 

PAR_boundary %>% ungroup %>%  summarize(chr=unique(chr),start=min(start_range),end=max(end_range))
# chr      start      end
# <chr>    <dbl>    <dbl>
#   1 Z     77659000 78141000

parplot = PAR_regions %>% filter(start >= 77e6) %>% 
  ggplot(aes(x=start,y=mean,col=PAR))+
  geom_point()+
  geom_hline(yintercept=0,lty=2)+
  geom_vline(xintercept=c(min(PAR_boundary$start_range),max(PAR_boundary$end_range)),lty=2)+
  theme_bw()

pdf('../figures/PAR_Plot_2024MAR22.pdf',height=3.5,width=7)
parplot
dev.off()
```

![image-20240322112728027](C:\Users\herit\AppData\Roaming\Typora\typora-user-images\image-20240322112728027.png)

