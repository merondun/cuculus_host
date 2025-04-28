#  Cuculus Evolutionary Genetics

Here we examine population genetic variation across Eurasian *Cuculus*. Herein lies our primary target, the diversification of host-specific egg morphology.

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

Metadata summaries:

```bash
#### metadata summaries 
setwd('~/EvoBioWolf/CUCKOO_gentes/')
.libPaths('~/mambaforge/envs/R/lib/R/library')
library(tidyverse)

md = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt')

# Species count
md %>% count(SpeciesShort)

# 90% of sampling of canorus is within...
md %>% 
  filter(SpeciesShort == 'CC') %>% 
  pull(Sampling_Year) %>% 
  quantile(probs = c(0.05, 0.95), na.rm = TRUE)

# Eggs (also unrelated)
md %>% drop_na(Egg) %>% count(SpeciesShort)
md %>% filter(Analysis_PopulationGenetics == 1) %>% drop_na(Egg) %>% count(SpeciesShort)

# Tissue
md %>% filter(grepl('CO|CC',SpeciesShort)) %>% mutate(Source = ifelse(grepl('Tissue|Muscle|Heart|Embryo',Source),'Muscle',
                                                                      ifelse(grepl('Feather',Source),'Feather','Blood'))) %>% count(Source)

# Reads
summary(md$Raw_Bases_Gb)
sum(md$Raw_Bases_Gb)
sum(md$Raw_Reads_M)
summary(md$Mapped_Reads_M)
summary(md$MeanCoverage)

# Gene Hunt
md %>% filter(Analysis_GeneHunt == 1) %>% count(SpeciesShort)

```



# Processing

## Trimming

Business as usual; have hundreds of individual SRRs to trim and align. BBDUK + BWA. 

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

Still on the individual SRR accessions. Submit positional library.

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

Split Genome into chunks:

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

Calls SNPs with BCFTools

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

## SNP Filtering

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

## Relative Analyses

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

Strict mapping and missingness, since we want to be incredibly accurate:

```bash
bcftools view chr_MT_All.SNP.DP3.vcf.gz --samples-file Samples.list | bcftools view --min-ac 1 -e 'F_MISSING > 0.05 || MQ < 40' -Oz -o chr_MT.DP3.MM05.MQ40.vcf.gz
```

Count snps: `python 02_Count_mtDNA_Mutations.py --vcf chr_MT.DP3.MM05.MQ40.vcf.gz --out chrMT_Pairwise_Differences.txt` 

```python
import allel
import itertools
import numpy as np
import gzip
import argparse

# Read VCF > genotypeArray
def read_vcf(filename):
    callset = allel.read_vcf(filename)
    genotypes = allel.GenotypeArray(callset['calldata/GT'])
    return genotypes, callset['samples']

# Compute pairwise comparisons
def pairwise_comparisons(genotypes, samples):
    results = []
    ploidy = genotypes.ploidy  #ploidy 

    for id1, id2 in itertools.combinations(range(len(samples)), 2):
        g1 = genotypes[:, id1]
        g2 = genotypes[:, id2]

        # Filter out missing data
        non_missing = ~(g1.is_missing() | g2.is_missing())
        g1 = g1.compress(non_missing, axis=0)
        g2 = g2.compress(non_missing, axis=0)

        # Count the number of valid sites
        number_sites = non_missing.sum()

        #Count the number of SNPs correctly accounting for ploidy
        number_snps = np.count_nonzero(g1.to_n_alt() != g2.to_n_alt())

        results.append((samples[id1], samples[id2], number_sites, number_snps))
    return results

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compute pairwise SNP differences from a VCF file.')
    parser.add_argument('--vcf', required=True, help='Input VCF file (gzipped)')
    parser.add_argument('--out', required=True, help='Output file for pairwise comparisons')
    args = parser.parse_args()

    # Read VCF
    genotypes, samples = read_vcf(args.vcf)

    # Perform pairwise comparisons
    comparisons = pairwise_comparisons(genotypes, samples)

    # Output results
    with open(args.out, 'w') as f:
        for comparison in comparisons:
            f.write(f'{comparison[0]} {comparison[1]} {comparison[2]} {comparison[3]}\n')

    print(f'Pairwise comparisons saved to {args.out}')

```

### Plot Shared Relative Phenotypes

Identify relatives, shared phenotypes among them:

```bash
#### Determine relatedness  
setwd('~/EvoBioWolf/CUCKOO_gentes/relatedness/')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(viridis)
library(scales)

# Analyze LD-pruned relatedness from vcftools --relatedness2, KING coefficient 
rel <- read.table('autos_canorus_LD.relatedness2',header=T) %>% select(1,2,7) %>% as_tibble
names(rel) <- c('ID_A','ID_B','PHI')
rel = rel %>% mutate(PHI = pmax(0, pmin(0.5, PHI)))
rel = rel %>% filter(ID_A != ID_B)

# Values from here for designating relationships https://www.kingrelatedness.com/manual.shtml
fam = rel %>% mutate(Relationship = ifelse(PHI > 0.354, 'First Degree',
                                           ifelse(PHI > 0.177, 'First Degree',
                                                  ifelse(PHI > 0.0884,'Second Degree',
                                                         ifelse(PHI > 0.0442, 'Third Degree',
                                                                ifelse(PHI <= 0.0442,'Unrelated','Unassigned'))))))

fam %>% filter(PHI > 0.354) #usually > 0.354 is MZ twin / duplicate, but since our highest value is 0.393 and most around 0.37, seems more likely they are just first degree

# Add mtDNA differences, calculated externally which identifies pairwise SNPs between samples for mtDNA
mt = read.table('chrMT_Pairwise_Differences.txt')

names(mt) = c('ID_A','ID_B','Sites','SNPs')
# The script counted diploid genotypes as 2 SNPs, so divide by 2
#mt <- mt %>% mutate(SNPs = SNPs / 2)
fam2 = left_join(fam,mt) %>% drop_na(SNPs)

# Remove redundant comparisons, e.g. ID_A vs ID_B or pairwise redundancies
famrm = fam2 %>% 
  rowwise() %>% 
  mutate(pair = sort(c(ID_A,ID_B)) %>% paste(collapse = ",")) %>%
  group_by(pair,Sites,SNPs) %>%
  distinct(pair, .keep_all = T) %>% 
  separate(pair,into=c('ID_A','ID_B'),remove=F,sep=',') %>% ungroup() %>% select(-pair)

# Merge with metadata
md = read_tsv('Full_Metadata.txt')v
md <- read_tsv('../Metadata_Host.txt')
md = md %>% select(c(ID,Egg,GeographicGroup,Sampling_Year,Haplogroup,Sex,Age))
mda = md

# Add an '_A' and '_B' to each metadata field 
names(mda) = paste0(names(mda),'_A')
mdb = md
names(mdb) = paste0(names(mdb),'_B')
fam3 = left_join(famrm,mda) %>% 
  left_join(.,mdb)

# Inspect, look to see which Haplogroups don't align, etc 
fam3 %>% filter(Relationship == 'First Degree' & Haplogroup_A != Haplogroup_B) %>% data.frame
fam3 %>% filter(Relationship == 'First Degree' & GeographicGroup_A != GeographicGroup_B) %>% data.frame
fam3 %>% filter(Relationship == 'First Degree' & Egg_A != Egg_B) %>% data.frame

# Only retain those with known egg morph, which are NESTLINGS, and which are caught within 2 years ( no parent / sibling)
famd = fam3 %>% drop_na(Egg_A) %>% drop_na(Egg_B)
famd = famd %>% filter(GeographicGroup_A == GeographicGroup_B) # Only compare within the same distance clade
famd = famd %>% filter(Age_A == 'Young' & Age_B == 'Young') # Only compare within YOUNG NESTLINGS!
famd %>% count(Age_A,Age_B)
famd = famd %>% filter(abs(Sampling_Year_A - Sampling_Year_B) <= 2) # Only compare within the same sampling year (e.g. SIBLINGS)
famd %>% count(Age_A,Age_B) #n = 2806

# Within first degree relatives, what's the distribution of the # of mtDNA snps? 
mtdna_mutations <- famd %>% filter(Relationship == 'First Degree') %>% 
  ggplot(aes(x=SNPs))+geom_histogram(show.legend = F)+theme_bw()+
  scale_x_continuous(breaks=pretty_breaks(n=14))

png('20250207_mtDNAMutations.png',units='in',res=300,height=3,width=5)
mtdna_mutations
dev.off()

# Unhash this for making the plots with difference mtDNA difference thresholds
#for (snp_count in seq(0,3,1)){

# snp_count = 2 shows a similar number of unrelated paternal / maternal, and first degree more maternal (n=81 compared to 34 paternal), which we expect for cuckoos 
snp_count <- 0

# If mtDNA Haplogrouplotypes are identical, assign as maternal. 
famd = famd %>% mutate(Line = ifelse(SNPs <= snp_count,'Maternal','Paternal'))
famd %>% filter(Relationship == 'First Degree') %>%  count(Line)
# Line         n
#   1 Maternal    56
# 2 Paternal    59

famd %>% filter(Relationship != 'Unrelated') %>% count(Line)
# Line         n
#   1 Maternal   112
# 2 Paternal   163

#simply assign unrelated as paternal, move pie chart in final plot 
famd = famd %>% mutate(Line = ifelse(Relationship == 'Unrelated','Paternal',Line))

# and remove first degree relatives, OR NOT! 
#famd = famd %>% filter(Relationship != 'First Degree')

famd %>% filter(Relationship != 'Unrelated') %>% count(Relationship,Line)
famd %>% filter(Relationship != 'Unrelated') %>% nrow
# Relationship  Line         n
#   1 First Degree  Maternal    56
# 2 First Degree  Paternal    59
# 3 Second Degree Maternal    29
# 4 Second Degree Paternal    38
# 5 Third Degree  Maternal    27
# 6 Third Degree  Paternal    66

nrow(famd)
#[1] 2806 with first degree, or [1] 2948 with no age and year restrictions
#how many unique individuals?
length(unique(c(famd$ID_A,famd$ID_B)))

#function to get matched data
get_matched_data <- function(df) {
  df %>% 
    mutate(
      Year = ifelse(abs(Sampling_Year_A - Sampling_Year_B) <= 2, 'Matched', 'Unmatched'),
      Egg = ifelse(Egg_A == Egg_B, 'Matched', 'Unmatched'),
      Distance = ifelse(GeographicGroup_A == GeographicGroup_B, 'Matched', 'Unmatched'),
      Haplogrouplogroup = ifelse(Haplogroup_A == Haplogroup_B, 'Matched', 'Unmatched')
    ) %>% 
    gather(key = "Variable", value = "Matched", Year, Egg, Distance,Haplogrouplogroup) %>% 
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
piece = piece %>% filter(!grepl('Year|Distance|Haplogroup',Variable))

#plot
relpie = piece %>% 
  ggplot(aes(x="",y=Proportion,fill=Matched))+
  geom_bar(stat='identity')+
  coord_polar("y", start=0)+xlab('')+ylab('')+
  facet_grid(Variable+Line~Relationship)+
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

pdf(paste0('20250207-Relatives_Phenotype_Matching-',snp_count,'Mismatch-FIRSTDEGREE.pdf'),height=4.5,width=6)
print(piesR)
dev.off()

# png(paste0('20250207-Relatives_Phenotype_Matching-',snp_count,'Mismatch-FIRSTDEGREE.png'),units='in',res=300,height=4.5,width=6)
# print(piesR)
# dev.off()

#} # Unhash this for doing the loop based on mtDNA differences 

famd %>% filter(Relationship != 'Unrelated' & Line == 'Paternal') %>%
  count(Relationship,Egg_A,Egg_B)


famd %>% filter(Relationship != 'Unrelated' & Line == 'Maternal') %>% 
  count(Relationship,Egg_A,Egg_B)


famd %>% count(Relationship,Egg_A,Egg_B,Line)

## Aggregate
df_agg <- famd %>% 
  count(Relationship,Egg_A,Egg_B,Line) %>% 
  mutate(
    Egg = pmin(Egg_A, Egg_B),  # ensures consistent grouping
    Match = ifelse(Egg_A == Egg_B, n, 0),
    Unmatch = ifelse(Egg_A != Egg_B, n, 0)
  ) %>%
  group_by(Relationship, Egg, Line) %>%
  summarise(
    Match = sum(Match),
    Unmatch = sum(Unmatch),
    .groups = "drop"
  )

# Supplementary table 
write.table(df_agg,file='20250207_EggMismatch_Counts.txt',quote=F,sep='\t',row.names=F)

# Plot matched 
counts_mismatch <- df_agg %>% 
  filter(Relationship != 'Unrelated') %>% 
  pivot_longer(c(Match,Unmatch)) %>% 
  ggplot(aes(y=Egg,x=value,fill=name))+
  xlab('Count')+
  scale_fill_manual('Matched',values=c('forestgreen','grey95'))+
  geom_bar(col='black',stat='identity',position=position_dodge(width=0.5),lwd=0.25)+
  facet_grid(Relationship~Line,scales='free',space='free_y')+
  theme_bw(base_size=8)+
  scale_x_continuous(breaks = pretty_breaks(n = 3))
counts_mismatch

pdf('~/symlinks/host/figures/20250404-Relatives_EggMismatch.pdf',height=4.5,width=3)
counts_mismatch
dev.off()

# Perform Fisher's test for each Relationship category
results <- df_agg %>%
  filter(Relationship != 'Unrelated') %>% 
  group_by(Relationship) %>%
  summarise(
    p_value = list(
      fisher.test(matrix(c(sum(Match[Line == "Maternal"]),
                           sum(Unmatch[Line == "Maternal"]),
                           sum(Match[Line == "Paternal"]),
                           sum(Unmatch[Line == "Paternal"])),
                         nrow = 2))$p.value
    )
  ) %>%
  unnest(p_value) %>% 
  mutate(padj = p.adjust(p_value,method='bonferroni'))
options(scipen = 9999)
print(results)
# Relationship  p_value   padj
#   1 First Degree  1       1     
# 2 Second Degree 0.502   1     
# 3 Third Degree  0.00460 0.0138


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
```



## SNP Filtering Unrelated Samples

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

## Estimate Maternal Phylogenies

Subset samples. In general we have 2 sample sets: 

* one with females (n=87) and outgroups (n=6) together for BEAST (n=93) .
* one with all samples for mtDNA (n=136) and females for W (n=87). 

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --mem=30000mb
#SBATCH --cpus-per-task=10
#SBATCH --time=24:00:00

WD=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/phylogenetics/202503_fulln136
raw_vcfs=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/merged/snps_only
#mask with male-biased coverage
mask=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/Coverage_Masks/GCA_017976375.1_bCucCan1.pri_genomic.CHR.N75-DoubleCoverage.mask.bed

cd ${WD}

#mamba activate snps
# sbatch 1.Subset_Samples_Filter.sh
mkdir -p vcfs ml_trees

for CHR in chr_MT chr_W; do

#minimum coverage, LESS than this set to missing
MINDP=3
echo "FILTERING AND MERGING VARIANT SITES FOR ${CHR}"
bcftools view --force-samples --threads 10 --samples-file Samples.list -Ou ${raw_vcfs}/${CHR}.SNPS.vcf.gz | \
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
        bcftools view --min-alleles 2 --max-alleles 2 --types snps --min-ac 1 -Oz -o vcfs/${CHR}.SNP.DP3.vcf.gz
bcftools index --threads 10 vcfs/${CHR}.SNP.DP3.vcf.gz

#filter, include singletons
bcftools view vcfs/${CHR}.SNP.DP3.vcf.gz --min-alleles 2 --max-alleles 2 --min-ac 1 --max-af 0.999 --types snps -i "MQ>40 & F_MISSING < 0.1" -Oz -o vcfs/${CHR}.SNP.DP3-AC1-MQ40-MM1.vcf.gz
bcftools index --threads 10 vcfs/${CHR}.SNP.DP3-AC1-MQ40-MM1.vcf.gz

#filter, exclude singletons, note that this must be --MIN-AC 3, BECAUSE ITS PSEUDO-HAPLOID
bcftools view vcfs/${CHR}.SNP.DP3.vcf.gz --min-alleles 2 --max-alleles 2 --min-ac 3 --max-af 0.999 --types snps -i "MQ>40 & F_MISSING < 0.1" -Oz -o vcfs/${CHR}.SNP.DP3-AC2-MQ40-MM1.vcf.gz
bcftools index --threads 10 vcfs/${CHR}.SNP.DP3-AC2-MQ40-MM1.vcf.gz

#create tree
python ~/modules/vcf2phylip.py -i vcfs/${CHR}.SNP.DP3-AC1-MQ40-MM1.vcf.gz -f --output-folder ml_trees
iqtree --redo -keep-ident -T 10 -s ml_trees/${CHR}.SNP.DP3-AC1-MQ40-MM1.min4.phy --seqtype DNA -m "MFP+ASC" -alrt 1000 -B 1000
iqtree --redo -keep-ident -T 10 -s ml_trees/${CHR}.SNP.DP3-AC1-MQ40-MM1.min4.phy.varsites.phy --seqtype DNA -m "MFP+ASC" -alrt 1000 -B 1000

#create tree NO SINGLETONS
python ~/modules/vcf2phylip.py -i vcfs/${CHR}.SNP.DP3-AC2-MQ40-MM1.vcf.gz -f --output-folder ml_trees
iqtree --redo -keep-ident -T 10 -s ml_trees/${CHR}.SNP.DP3-AC2-MQ40-MM1.min4.phy --seqtype DNA -m "MFP+ASC" -alrt 1000 -B 1000
iqtree --redo -keep-ident -T 10 -s ml_trees/${CHR}.SNP.DP3-AC2-MQ40-MM1.min4.phy.varsites.phy --seqtype DNA -m "MFP+ASC" -alrt 1000 -B 1000

done

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

## Plot Phylogenies

### Plot W and mtDNA Trees

Trees were estimated here: 

```bash
/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/phylogenetics/202503_fulln136

#full susbset with n=252
/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/phylogenetics/202503_fulln252
```

Plot many trees:

```R
#### Plot many W and MT trees with different sample subsets 
setwd('~/EvoBioWolf/CUCKOO_gentes/phylogenetics/trees/')
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

md = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt')

###### CIRCULAR W TREE  #######
#plot tree 
trees = list.files('.',pattern='.*AC1.*contree',full.names = TRUE)
counter = 0 
for (tree in trees) { counter = counter + 1;
iqtree = read.iqtree(tree)
iqtr = midpoint.root(as.phylo(iqtree))
lab = gsub('.SNP.*','',gsub('.*chr_','',tree))
cat ('Making tree: ',lab,'\n')

# Haplogroups 
circ = ggtree(iqtr, layout = "circular",branch.length='none') %<+% md  +
  geom_nodepoint(mapping=aes(subset=(as.numeric(label) >= 95)),col='black',fill='grey90',pch=23,size=0.75,show.legend=F)+
  geom_tippoint(aes(fill = Haplogroup,shape=Species),size=1)+
  scale_fill_manual('Haplogroup',values=md$HaplogroupColor,breaks=md$Haplogroup)+
  scale_shape_manual(values=md$Shape,breaks=md$Species)+
  guides(fill=guide_legend(nrow=4,override.aes=list(shape=21)),
         shape=guide_legend(nrow=5))+
  ggtitle(lab)+
  theme(legend.text = element_text(size = 8),legend.title = element_text(size = 10),
        legend.key.size = unit(0.2, "cm"),    legend.position = 'top')
circ

apetree = ggtree(iqtr, layout = "ape") %<+% md  +
  geom_nodepoint(mapping=aes(subset=(as.numeric(label) >= 95)),col='black',fill='grey90',pch=23,size=0.75,show.legend=F)+
  geom_tippoint(aes(fill = Haplogroup,shape=Species),size=1)+
  scale_fill_manual('Haplogroup',values=md$HaplogroupColor,breaks=md$Haplogroup)+
  scale_shape_manual(values=md$Shape,breaks=md$Species)+
  guides(fill=guide_legend(nrow=4,override.aes=list(shape=21)),
         shape=guide_legend(nrow=5))+
  theme(legend.text = element_text(size = 8),legend.title = element_text(size = 10),
        legend.key.size = unit(0.2, "cm"),    legend.position = 'top')
apetree


# Egg
circ_egg = ggtree(iqtr, layout = "circular",branch.length='none') %<+% md  +
  geom_nodepoint(mapping=aes(subset=(as.numeric(label) >= 95)),col='black',fill='grey90',pch=23,size=0.75,show.legend=F)+
  geom_tippoint(aes(fill = Egg,shape=Egg),size=1)+
  scale_fill_manual('Egg',values=md$EggCol,breaks=md$Egg)+
  scale_shape_manual(values=md$EggShape,breaks=md$Egg)+
  ggtitle(lab)+
  theme(legend.text = element_text(size = 8),legend.title = element_text(size = 10),
        legend.key.size = unit(0.2, "cm"),    legend.position = 'top')
circ_egg

apetree_egg = ggtree(iqtr, layout = "ape") %<+% md  +
  geom_nodepoint(mapping=aes(subset=(as.numeric(label) >= 95)),col='black',fill='grey90',pch=23,size=0.75,show.legend=F)+
  geom_tippoint(aes(fill = Egg,shape=Egg),size=1)+
  scale_fill_manual('Egg',values=md$EggCol,breaks=md$Egg)+
  scale_shape_manual(values=md$EggShape,breaks=md$Egg)+
  guides(fill=guide_legend(nrow=4,override.aes=list(shape=21)),
         shape=guide_legend(nrow=5))+
  theme(legend.text = element_text(size = 8),legend.title = element_text(size = 10),
        legend.key.size = unit(0.2, "cm"),    legend.position = 'top')
apetree_egg
ca = ggarrange(circ,apetree,nrow=2,heights=c(0.6,0.4),legend = 'none') 
ca_egg = ggarrange(circ_egg,apetree_egg,nrow=2,heights=c(0.6,0.4),legend = 'none') 
assign(paste0('p',counter),ca)
assign(paste0('e',counter),ca_egg)

}

ggarrange(e1,e2,p1,p2)
ggsave('~/symlinks/host/figures/20250330_MT-W-Trees-All.pdf',
       ggarrange(e1,e2,p1,p2),
       height=9,width=6,dpi=600)

```



### Plot Collapsed Tree

Collapse chrW tree on supported nodes for visualization (fig 1)

```R
#### Collapse chrW Tree 
setwd('~/EvoBioWolf/CUCKOO_gentes/phylogenetics')
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

md = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt')

##### COLLAPSE W TREE #####
w = read.iqtree('trees/chr_W.SNP.DP3-AC1-MQ40-MM1.min4.phy.contree')
w <- root(as.phylo(w),outgroup = '387_CP_MBW_RUS_F')
outs <- w$tip.label[grepl('_CM_|_CP_',w$tip.label)]
w_tree <- drop.tip(w,outs)

p = ggtree(w_tree, layout = "dendrogram") %>% ggtree::rotate(166) %>% ggtree::rotate(168) %>% ggtree::rotate(123) %<+% md 
inspect_tree <- p + geom_tippoint(aes(fill = Haplogroup),pch=21,size=4)+scale_fill_manual(values=md$HaplogroupColor,breaks=md$Haplogroup) + 
  geom_nodelab(aes(label=node),geom = 'label',size=2) 
ggsave('~/symlinks/host/figures/20250319_W_TREE-COLLAPSE_Inspection.pdf',
       inspect_tree,height=15,width=15,dpi=300,limitsize = FALSE)

hapcols = md %>% select(Haplogroup,HaplogroupColor) %>% unique %>% arrange(Haplogroup) %>% filter(!grepl('CM|CP',Haplogroup)) %>% na.omit

#collapse
p = ggtree(w_tree, layout = "dendrogram") %>% 
  ggtree::rotate(166) %>% ggtree::rotate(168) %>% ggtree::rotate(123) %<+% md + 
  geom_tippoint(aes(fill=Haplogroup,shape=SpeciesShort),size=3) +
  scale_fill_manual(values=md$HaplogroupColor,breaks=md$Haplogroup) + 
  scale_shape_manual(values=md$Shape,breaks=md$SpeciesShort)+
  geom_nodepoint(mapping=aes(subset=(as.numeric(label) >= 95)),col='black',fill='grey90',pch=23,size=1.5,show.legend=F)
#MCC1
p2 = p %>% collapse(node=117) + geom_point2(aes(subset=(node==117)), shape=21, size=3, fill=hapcols$HaplogroupColor[hapcols$Haplogroup == 'MCC1']); p2
#MCC2
p2 = collapse(p2,node=91) + geom_point2(aes(subset=(node==91)), shape=21, size=3, fill=hapcols$HaplogroupColor[hapcols$Haplogroup == 'MCC2']); p2
#MCC3
p2 = collapse(p2,node=99) + geom_point2(aes(subset=(node==99)), shape=21, size=3, fill=hapcols$HaplogroupColor[hapcols$Haplogroup == 'MCC3']); p2
#BIN MCC4 MCO1 and MCO2 (rufous)
p2 = collapse(p2,node=167) + geom_point2(aes(subset=(node==167)), shape=21, size=3, fill=hapcols$HaplogroupColor[hapcols$Haplogroup == 'MCC4']); p2
#MCC5
p2 = collapse(p2,node=162) + geom_point2(aes(subset=(node==162)), shape=21, size=3, fill=hapcols$HaplogroupColor[hapcols$Haplogroup == 'MCC5']); p2
#MCC6
p2 = collapse(p2,node=160) + geom_point2(aes(subset=(node==160)), shape=21, size=3, fill=hapcols$HaplogroupColor[hapcols$Haplogroup == 'MCC6']); p2
#MCC7
p2 = collapse(p2,node=125) + geom_point2(aes(subset=(node==125)), shape=21, size=3, fill=hapcols$HaplogroupColor[hapcols$Haplogroup == 'MCC7']); p2
#MCO3
p2 = collapse(p2,node=143) + geom_point2(aes(subset=(node==143)), shape=24, size=3, fill=hapcols$HaplogroupColor[hapcols$Haplogroup == 'MCO3']); p2
#MCO4
p2 = collapse(p2,node=145) + geom_point2(aes(subset=(node==145)), shape=24, size=3, fill=hapcols$HaplogroupColor[hapcols$Haplogroup == 'MCO4']); p2
p2 = p2 + geom_treescale(x=0.05,offset = 0.01) 

ggsave('~/symlinks/host/figures/20250319_W_TREE-COLLAPSED.pdf',
       p2,height=1,width=6,dpi=300)

```

### Plot Dendrogram with Phenotype Bars

```bash
#### Plot Trees, Assign Haplogroups 
setwd('~/EvoBioWolf/CUCKOO_gentes/phylogenetics/')
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
library(rcartocolor)

md = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt') 

#only grab samples with known egg (hash out filter for related individuals)
md_egg = md %>%
  filter(Analysis_FullGensAssociations == 1) %>%
  select(ID = ID, Egg, Haplogroup, Geography = GeographicGroup,Ancestry = AncestryA5)
md_egg %>% count(Egg)

md_egg = md_egg %>% group_by(Egg) %>% mutate(TotalEgg = n()) %>% ungroup
minobs=2
md_egg = md_egg %>% filter(TotalEgg >= minobs)
md_egg %>% count(Egg)

## Plot Trees
## mtDNA
mt = read.iqtree('trees/chr_MT.SNP.DP3-AC1-MQ40-MM1.min4.phy.contree')
mt <- midpoint.root(as.phylo(mt))

# voila 
mtp = ggtree(mt, layout = "dendrogram") %<+% md  + 
  #geom_nodelab(aes(label=node),geom = 'label',size=3)+  
  geom_nodepoint(mapping=aes(subset=(as.numeric(label) >= 95)),col='black',fill='grey90',pch=23,size=1.5,show.legend=F)+
  geom_tippoint(aes(fill = Haplogroup,shape=SpeciesShort),size=2)+
  scale_shape_manual(values=c(21,24))+
  scale_fill_manual(values=md$HaplogroupColor,breaks=md$Haplogroup)+
  #scale_fill_manual('Egg',values=md$HostCol,breaks=md$Egg)+
  #geom_tiplab()+
  guides(fill=guide_legend(nrow=2,override.aes=list(shape=21)))+
  theme(legend.position='top')
mtp

## W chromosome
w = read.iqtree('trees/chr_W.SNP.DP3-AC1-MQ40-MM1.min4.phy.contree')
w <- midpoint.root(as.phylo(w))

#plot cladogram, first find nodes to collapse 
wp = ggtree(w, layout = "dendrogram") %<+% md  + 
  geom_nodepoint(mapping=aes(subset=(as.numeric(label) >= 95)),col='black',fill='grey90',pch=23,size=1.5,show.legend=F)+
  geom_tippoint(aes(fill = Haplogroup,shape=SpeciesShort),size=2)+
  scale_shape_manual(values=c(21,24))+
  scale_fill_manual(values=md$HaplogroupColor,breaks=md$Haplogroup)+
  guides(fill=guide_legend(nrow=2,override.aes=list(shape=21)))+
  theme(legend.position='top')
wp

ggsave('~/symlinks/host/figures/20250330_N136-All-MT-W.pdf',
       ggarrange(mtp,wp,nrow=2,common.legend = TRUE),width=6,height=5,dpi=300)


#### W tree with egg: ancestry+geography+haplogroup
md <- md %>% arrange(Haplogroup,HaplogroupColor)
p = ggtree(w, layout = 'dendrogram') %<+% md + 
  geom_nodepoint(mapping=aes(subset=(as.numeric(label) >= 95)),col='black',fill='grey90',pch=23,size=0.75,show.legend=F)+
  geom_tippoint(aes(fill=Haplogroup,shape=SpeciesShort),size=1.25)+
  geom_tiplab(align=TRUE,size=1)+
  scale_shape_manual(values=c(21,24))+
  scale_fill_manual(values=md$HaplogroupColor,breaks=md$Haplogroup)+
  geom_treescale()+
  guides(fill=guide_legend(override.aes=list(shape=21)))
p

# Alt, with egg shape
p = ggtree(w, layout = 'dendrogram',branch.length = 'none') %<+% md + 
  geom_nodepoint(mapping=aes(subset=(as.numeric(label) >= 95)),col='black',fill='grey90',pch=23,size=0.75,show.legend=F)+
  geom_tippoint(aes(fill=Egg,shape=Egg),size=1.25)+
  geom_tiplab(align=TRUE,size=1)+
  scale_shape_manual(values=md$EggShape,breaks=md$Egg)+
  scale_fill_manual(values=md$EggCol,breaks=md$Egg)+
  geom_treescale()
p

phenos = as.data.frame(p$data %>% filter(isTip == TRUE))
rownames(phenos) = phenos$label

#second response variable: egg
p1 = p + new_scale_fill()
egg = phenos %>% select(E = Egg)
md <- md %>% arrange(EggOrder)
p1 = gheatmap(p1, egg,offset=.015,width=0.05) +
  scale_fill_manual('E',values=md$EggCol,breaks=md$Egg,na.value='white')
p1

#third variable: ancestry
p2 = p1 + new_scale_fill()
anc = phenos %>% select(N = AncestryA5)
md <- md %>% arrange(AncestryA5)
p2 = gheatmap(p2, anc,offset=.03,width=0.05) +
  scale_fill_manual('N',values=md$AncestryColor,breaks=md$AncestryA5)
p2

# fourth: geography
p3 = p2 + new_scale_fill()
md <- md %>% arrange(GeographicGroup)
geo = phenos %>% select(G = GeographicGroup)
p3 = gheatmap(p3, geo,offset=.045,width=0.05) +
  scale_fill_manual('G',values=md$GeoColor,breaks=md$GeographicGroup)
p3

# # Save first as haplogroup
# ggsave('~/symlinks/host/figures/20250330_Dendrogram-W-Eggs-Phenotypes-HAPLOGROUP.pdf',
#        p3,
#        height=5,width=8,dpi=300)

# And then save with egg 
ggsave('~/symlinks/host/figures/20250330_Dendrogram-W-Eggs-Phenotypes-EGG.pdf',
       p3,
       height=5,width=8,dpi=300)

```



## TWISST

To go deeper than C. poliocephalus, download Clamator glandarius female data (HiFi) from the interwebs. Align it to the genome, and call SNPs (only on the W for brevity).

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=20
#SBATCH --time=48:00:00

WD=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/phylogenetics/202503_twisst/deeper
READS=/dss/dsslegfs01/pr53da/pr53da-dss-0021/rawdata/External_Data_SRA_GenomeArk/SRA_Cuculiformes_Outgroups/Clamator_glandarius_SRR26807982_1.fastq.gz
GENOME=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa

mkdir -p coverage

# Align reads
pbmm2 align \
    --preset HIFI \
    --sort \
    --num-threads 20 \
    --rg "@RG\tID:1\tSM:999_CG_UNK_UNK_F\tPL:HiFi\tLB:999_CG_UNK_UNK_F" \
    ${GENOME} ${READS} Clamator_glandarius.bam
samtools index Clamator_glandarius.bam

# call SNPs on chrW
bcftools mpileup --threads 20 -Ou -r chr_W -f ${GENOME} Clamator_glandarius.bam | \
        bcftools call --threads 20 -mv -Oz -o Clamator_glandarius.chrW.vcf.gz
bcftools index Clamator_glandarius.chrW.vcf.gz

# coverage
mosdepth --threads 20 --mapq 30 --by 1000000 --fast-mode --no-per-base coverage/Clamator_glandarius Clamator_glandarius.bam
```

Prepare and run:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=10
#SBATCH --time=200:00:00

#mamba activate twisst

WD=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/phylogenetics/202503_twisst
#VCF=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/phylogenetics/202503_snapp/vcfs/chr_W.SNP.DP3-AC1-MQ40-MM1.vcf.gz
VCF=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/phylogenetics/202503_twisst/deeper/vcfs/chr_W.SNP.DP3-AC2-MQ40-MM1.vcf.gz

cd ${WD}

mkdir -p input output

# Force haploid
bcftools view ${VCF} | \
  bcftools +fixploidy -Ou - -- -f 1 | \
  bcftools view -Oz -o input/chr_W.SNP.DP3-AC1-MQ40-MM1-HAPLOID.vcf.gz

# Create geno file
python ~/modules/genomics_general/VCF_processing/parseVCF.py --ploidyMismatchToMissing --ploidy 1 --skipIndels -i input/chr_W.SNP.DP3-AC1-MQ40-MM1-HAPLOID.vcf.gz | \
        bgzip > input/chr_W.SNP.DP3-AC1-MQ40-MM1-HAPLOID.geno.gz

#only retain SNPs with MAC >= 2
python ~/modules/genomics_general/filterGenotypes.py --threads 10 -i input/chr_W.SNP.DP3-AC1-MQ40-MM1-HAPLOID.geno.gz -o input/chr_W.SNP.DP3-AC2-MQ40-MM1-HAPLOID.geno.gz --minVarCount 2 --minAlleles 2 --maxAlleles 2

for win in 50 100 500 1000; do
  miss=$(echo ${win} | awk -v w=${win} '{print w/5}')
  missind=$(echo ${win} | awk -v w=${win} '{print w/10}')

  #calculate trees in windows
  python ~/modules/genomics_general/phylo/phyml_sliding_windows.py -T 10 \
          -g input/chr_W.SNP.DP3-AC2-MQ40-MM1-HAPLOID.geno.gz --prefix output/chr_W.SNP.DP3-AC2-MQ40-MM1-HAPLOID.phyml_bionj.Sw${win}-m${miss}-ind${missind} \
         --windSize ${win} --minSites ${miss} --minPerInd ${missind} --windType sites --model GTR --optimise n

  #run twisst
  ~/modules/twisst/twisst.py -t output/chr_W.SNP.DP3-AC2-MQ40-MM1-HAPLOID.phyml_bionj.Sw${win}-m${miss}-ind${missind}.trees.gz \
    -w output/chr_W.SNP.DP3-AC2-MQ40-MM1-HAPLOID.weights.Sw${win}-m${miss}-ind${missind}.csv.gz --outgroup CG \
    --outputTopos output/chr_W.SNP.DP3-AC2-MQ40-MM1-HAPLOID.topologies.Sw${win}-m${miss}-ind${missind}.trees -g CG -g CP -g CM -g BLUE -g RUFOUS -g OTHER --method complete --groupsFile Females_N94-Simple.pop

done

```

### Plot

```R
#### Plot twisst
setwd('~/EvoBioWolf/CUCKOO_gentes/phylogenetics/twisst/')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(viridis)
library(RColorBrewer)
library(ggpubr)
library(meRo)
library(ggdist)
source('~/modules/twisst/plot_twisst.R')

# Import data 
weights_files = list.files('.',pattern='.*weights.*.csv.gz')
window_data_files = list.files('.',pattern='.*data.tsv')

windat <- list()
for (size in seq(1,4,1)) {
# Read in the data
win <- gsub('-.*','',gsub('.*weights.Sw','',weights_files[size]))
twisst_data <- import.twisst(weights_files=weights_files[size],window_data_files=window_data_files[size])
data <- twisst_data$weights_raw$chr_W
data$win <- as.numeric(win)
data <- data %>% pivot_longer(!win)
windat[[size]] <- data
}
full_dat <- rbindlist(windat) %>% as_tibble

topord <- full_dat %>% select(name) %>% distinct %>% mutate(ord=as.numeric(gsub('topo','',name))) %>% arrange(ord)
full_dat$name <- factor(full_dat$name,levels=topord$name)

# Summarize topos
topo_sums <- full_dat %>% 
  group_by(name,win) %>% 
  sum_stats(value)
top_topos <-  topo_sums %>% ungroup  %>%  filter(win == 500) %>% slice_max(mean,n=10) 
topo_sums$name <- factor(topo_sums$name,levels=top_topos$name)
extract_topos <- top_topos %>% head(n=10) %>% pull(name)
plot <- topo_sums %>% 
  filter(name %in% extract_topos) %>% 
  ggplot(aes(x=name,y=mean,ymin=conf_low,ymax=conf_high,col=name))+
  geom_point()+
  facet_grid(win~.,scales='free')+
  scale_color_manual(values=topo_cols)+
  geom_errorbar()+
  theme_bw(base_size=6)+
  coord_cartesian(ylim=c(0,1))
plot

ggsave('~/symlinks/host/figures/20250325_TWISST_chrW_Sensitivity.pdf',
       plot, height = 5,width=7,dpi=300)

##### Plot topologies using Ape #####
pdf('~/symlinks/host/figures/20250320_Twisst_Plotted_Topologies.pdf',height=1.25,width=15)

par(mfrow = c(1,length(twisst_data$topos)), mar = c(1,1,2,1), xpd=NA)
for (n in 1:length(twisst_data$topos)){
  plot.phylo(twisst_data$topos[[n]], type = "clad", edge.color=topo_cols[n], edge.width=3, label.offset=.2, cex = 1)
  mtext(side=3,text=paste0("topo",n))
}

dev.off()

topo_ids <- as.integer(gsub('topo','',extract_topos))
pdf('~/symlinks/host/figures/20250325_Twisst_Top_Topologies.pdf',height=1.25,width=15)

par(mfrow = c(1,length(topo_ids)), mar = c(1,1,2,1), xpd=NA)
for (i in seq_along(topo_ids)) {
  t <- topo_ids[i]
  plot.phylo(twisst_data$topos[[t]], type = "clad", edge.color=topo_cols[i], edge.width=3, label.offset=.2, cex = 1)
  mtext(side=3, text=paste0("topo", i))
}

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
# with CG: 3.5 0.25 gives 95% 20.3 - 54.1
100M chains, log every 1k 
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
T       G       C       A
6108268 4562746 4578162 6088177

chr_MT:
T       G       C       A
4749    2471    5775    6703

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

Treeannotator:

```bash
for i in $(ls *trees); do 
/dss/dsshome1/lxc07/di39dux/modules/beast/bin/treeannotator -b 10 -height mean -file ${i} ${i}.ann
done 
```



### Plot

```bash
#### Plot BEAST annotated trees 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/phylogenetics/202503_beast/inputs')
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
md = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt') %>% select(ID,Haplogroup,HaplogroupColor,SpeciesShort,Shape) 
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
  
  leg = md %>% select(Haplogroup,HaplogroupColor,SpeciesShort,Shape) %>% unique %>% drop_na(HaplogroupColor)
  gg$data = left_join(gg$data,df)
  ggp = gg  +
    geom_range(range='height_0.95_HPD', color='red', alpha=.6, size=1) +
    geom_tippoint(aes(fill = Haplogroup,shape=SpeciesShort),size=1.5)+
    geom_nodelab(aes(label=lab),size=1.5,vjust=1) +
    ggtitle(lab)+
    #geom_tiplab(size=2)+
    #geom_nodelab(aes(x=branch, label=round(posterior, 2)), vjust=-.5, size=3) +
    scale_color_continuous(low="darkgreen", high="red") +
    scale_fill_manual(values=leg$HaplogroupColor,breaks=leg$Haplogroup)+
    scale_shape_manual(values=leg$Shape,breaks=leg$SpeciesShort)+
    theme(legend.position=c(.1, .8))+
    geom_treescale(x = 5)+
    guides(fill=guide_legend(override.aes=list(shape=21)))+
    theme(legend.position='right')
  ggp
  assign(paste0('p',counter),ggp)
} 

ggarrange(p1,p2,p3,p4,common.legend = TRUE)

pdf('~/symlinks/host/figures/20250424_BEAST_Divergence_Dating_All-Dual.pdf',height=9,width=7)
ggarrange(p1,p2,p3,p4,common.legend = TRUE)
dev.off()

```



# Analyses: Population Genetic Differentiation

ADMIXTURE, PCA, Distance correlations: performed separately on C. canorus and C. optatus.

## Assign Geographic Clades 

Use k-means clustering to assign samples into discrete geographic 'populations' for analysis. 

```R
#### Assign geographic distance groups with k-means
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

### C. canorus

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



### C. optatus

```bash
# Plot PCA for C. optatus 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/manyhost_hunt/optatus')
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

md = read_tsv('metadata_optatus.txt')
#read plink
prefix = 'autos_optatus_LD'
vec = read.table('autosomal_files/autos_optatus_N50_LD.eigenvec',header=F)
names(vec) = c('ID',paste0('PC',seq(1,10,1)))
val = read.table('autosomal_files/autos_optatus_N50_LD.eigenval',header=F)
val = val %>% mutate(VE = V1/sum(V1))

#merge with metadata, add a $distance vector based on lat/long PC1
vc = left_join(vec,md)
summary(prcomp(vc[, c("Longitude", "Latitude")], scale. = TRUE)) #87% of variance 
vc = vc %>% 
  mutate(Longitude = as.numeric(Longitude), Latitude = as.numeric(Latitude)) %>%
  { cbind(., pca = prcomp(.[, c("Longitude", "Latitude")], scale. = TRUE)$x[,1]) } %>%
  rename(Distance = pca) 

#superimpose
world <- map_data("world")

# Convert your data to a simple feature object
sites <- st_as_sf(vc, coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant")

#add same custom grey color scale 
kcols = data.frame(Kcols = gray.colors(5,start=0.05,end=0.95)[c(1,4,3,5,2)], AncestryA5 = paste0('K',seq(1,5,1)))

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
new_pc = data.frame(nLong = proc_transform$X[,1], nLat = proc_transform$X[,2], nPC1 = proc_transform$Yrot[,1], nPC2 = proc_transform$Yrot[,2], ID = vc$ID, AncestryA5 = vc$AncestryA5, PC1 = vc$PC1, PC2 = vc$PC2, Lat = vc$Latitude, Long = vc$Longitude)

rp = new_pc %>% ggplot(aes(x=PC1,y=-PC2,fill=AncestryA5))+geom_point(pch=21)+scale_fill_manual(values=kcols$Kcols,breaks=kcols$AncestryA5)+
  xlab(paste0('PC1: ',round(val[1,2],3)*100,'%'))+  ylab(paste0('PC2: ',round(val[2,2],3)*100,'%'))+
  theme_bw(base_size=8)+theme(legend.position='top')
rp
new_pc %>% ggplot(aes(x=nPC1,y=nPC2,fill=AncestryA5))+geom_point(pch=21)+scale_fill_manual(values=kcols$Kcols,breaks=kcols$AncestryA5)+theme_bw(base_size=8)+ggtitle('Raw PCA')

#plot 
new_pc %>% ggplot()+
  geom_segment(aes(x=nLong,xend=nPC1,y=nLat,yend=nPC2),alpha=0.1)+
  geom_point(aes(x=nPC1,y=nPC2,fill=AncestryA5),pch=21)+
  scale_fill_manual(values=kcols$Kcols,breaks=kcols$AncestryA5)+
  theme_classic()

# Convert to simple feature object
library(scales)
no_rotation = vc %>% mutate(PC1_scaled = rescale(PC1*-1, to = c(min(Longitude), max(Longitude))),
                            PC2_scaled = rescale(PC2, to = c(min(Latitude), max(Latitude))))
sitesG <- st_as_sf(no_rotation, coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant")
sitesP <- st_as_sf(no_rotation, coords = c("PC1_scaled", "PC2_scaled"), crs = 4326, agr = "constant")

# Get world map data
world <- map_data("world")

# Plotting
pp = no_rotation %>% 
  ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), color = 'grey90', fill = 'white') +
  geom_segment(aes(x=PC1_scaled,xend=Longitude,y=PC2_scaled,yend=Latitude),alpha=0.1,col='darkred')+
  geom_point(aes(x=PC1_scaled,y=PC2_scaled,fill = AncestryA5), size = 1, pch = 21, show.legend = TRUE, color = 'grey20', alpha = 0.9) +
  geom_point(aes(x=Longitude,y=Latitude,fill = AncestryA5), size = 0.5, pch = 21, show.legend = TRUE, color = 'grey20', alpha = 0.2) +
  #geom_sf(data = sitesG, aes(fill = AncestryA5), size = 2, pch = 21, show.legend = TRUE, color = 'grey20', alpha = 0.2) +
  #geom_sf(data = sitesP, aes(fill = AncestryA5), size = 2, pch = 21, show.legend = TRUE, color = 'grey20', alpha = 0.9) +
  #geom_point(data = vc, aes(x = PC1_rot, y = PC2_rot), color = 'red', size = 2) +  # Add rotated PCA points
  scale_fill_manual(values = kcols$Kcols, breaks = kcols$AncestryA5) +  # Use unique Haplotypes for colors
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
pp

pdf('../../figures/20250127_PCA_Superimposed-SpatialDistribution_Coptatus.pdf',height=3,width=6)
ggarrange(rp,pp,widths=c(0.4,0.6))
dev.off()

#also slightly jitter points by 1degree for visualization 
# Convert your data to a simple feature object
sitesjit <- st_as_sf(vc %>% mutate(Longitude = jitter(Longitude,amount=1),
                                   Latitude = jitter(Latitude,amount=1)), 
                     coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant")

#add same custom grey color scale 
kcols = data.frame(Kcols = gray.colors(5,start=0.05,end=0.95)[c(1,4,3,5,2)], AncestryA5 = paste0('K',seq(1,5,1)))

# Plotting the map with PCA points and lines
kmap = ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),col='grey90',fill='white') + #add basemap
  geom_sf(data = sitesjit , 
          aes(fill=AncestryA5),
          size=2,pch=21,show.legend = T,col='grey20',alpha=0.9) +  #add the geometry for the points, I use pch 21-25 so the points have outlines
  scale_fill_manual(values=kcols$Kcols,breaks=kcols$AncestryA5)+ #custom fill encoded from metadata
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
          aes(fill=Haplogroup),
          size=2,pch=21,show.legend = T,col='grey20',alpha=0.9) +  #add the geometry for the points, I use pch 21-25 so the points have outlines
  scale_fill_manual(values=md$HaplogroupColor,breaks=md$Haplogroup)+ #custom fill encoded from metadata
  xlab('')+ylab('')+
  coord_sf(xlim = c(min(vc$Longitude)-5, max(vc$Longitude)+5), 
           ylim = c(min(vc$Latitude)-5, max(vc$Latitude)+5), expand = FALSE)+ #expand the map boundaries so that the labels show up, probably an easier way but I just add a gray rectangle in inkscape afterwards in between the labels and the map 
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"))+ #make the ocean blue
  theme(legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))+ 
  annotation_scale(line_width=0.5,location='bl')+  #add scale bar (requires ggspatial)
  guides(fill=guide_legend(override.aes=list(shape=21))) #ensure the legend has a fill-able shape
hapmap

pdf('../../figures/20250127_SpatialDistribution_Haplogroups-AncestryA4.pdf',height=7,width=7)
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

Also run for C. optatus:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=20
#SBATCH --time=48:00:00

# mamba activate snps

wd=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/manyhost_hunt/optatus
cd $wd

mkdir -p autosomal_files 

#merge autosomes
bcftools concat --file-list Autosomes.list --threads 20 -Oz -o autosomal_files/autos.vcf.gz
bcftools index --threads 20 autosomal_files/autos.vcf.gz

#LD prune
~/modules/plink2 --threads 20 --vcf autosomal_files/autos.vcf.gz --allow-extra-chr --set-missing-var-ids @:# \
        --rm-dup --indep-pairwise 50 5 0.1 --maf 0.05 --hwe 1e-10 --max-alleles 2 --min-alleles 2 --out autosomal_files/autos_optatus_LD
        
#extract, also a vcf and run PCA 
~/modules/plink2 --threads 20 --vcf autosomal_files/autos.vcf.gz --allow-extra-chr --set-missing-var-ids @:# \
        --extract autosomal_files/autos_optatus_LD.prune.in \
        --make-bed --recode vcf bgz --pca --out autosomal_files/autos_optatus_LD
bcftools index --threads 20 autosomal_files/autos_optatus_LD.vcf.gz
sed -i 's/chr_//g' autosomal_files/autos_optatus_LD.bim
```

### Plot, Fig 1B

Plot:

```R
### Plot ADMIXTURE across landscape (tesselation)
setwd('~/EvoBioWolf/CUCKOO_gentes/population_genetics/')
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
library(ggnewscale)

# Hash out which species to run, entire script will run afterwards 
admix_run = 'autos_optatus_N50_LD'
sp='CO'
admix_run = 'autos_canorus_LD'
sp='CC'

qdir = 'admixture_q_files' #directory with Q files

admix = melt_admixture(prefix = admix_run, qdir = qdir)

#read in metadata
md = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt')  
admixmd = left_join(admix,md) 

# Reorder individuals baseed on longitude
admixmd = admixmd %>% mutate(ID = fct_reorder(ID,Longitude))
if (sp == "CO") {
  kclust <- 'ACO'
  admixmd <- admixmd %>% mutate(K = gsub('K', kclust, K))
  kcols <- brewer.pal(7,'Purples')[c(1,3,5,6,7)]
  spshape=24
  filt="Cuculus optatus"
} else {
  kclust <- 'ACC'
  admixmd <- admixmd %>% mutate(K = gsub('K', kclust, K))
  kcols <- brewer.pal(6,'Greys')[c(2,3,4,5,6)]
  spshape=21
  filt="Cuculus canorus"
}

# Plot CV error 
cv = read.table(paste0(admix_run,'-ADMIXTURE.CVs.txt'),header=FALSE)
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
ggsave(paste0('~/symlinks/host/figures/20250318_ADMIXTURE-CV-Error_',sp,'.pdf'),cvs,height=2,width=3,dpi=600)

#if you want to add CV error directly on the label
names(cv) = c('Specified_K','d1','d2','Error')
cv = cv %>% mutate(label = paste0('K',Specified_K,' (',round(Error,2),')')) %>% select(!c(d1,d2,Error)) %>% arrange(Specified_K)
cv

# K = 2-5, by haplogroup  
adplot =
  admixmd %>% filter(Specified_K == 5 | Specified_K == 2) %>%  #specify the levels you want 
  mutate(Specified_K = paste0('K',Specified_K)) %>% 
  ggplot(aes(x = factor(ID), y = Q, fill = factor(K))) +
  geom_col(size = 0.1) +
  facet_grid(Specified_K~Haplogroup, scales = "free", space = "free") +
  theme_minimal(base_size=6) + labs(x = "",y = "") +
  scale_y_continuous(expand = c(0, 0),n.breaks = 3) +
  scale_fill_manual(values=kcols)+
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    #axis.text.x = element_blank(),
    axis.text.x=element_text(angle=90,size=5),
    axis.text.y = element_text(size=3),
    panel.grid = element_blank(),
    legend.position = 'bottom',
    plot.title = element_text(size=6)
  )
adplot

ggsave(paste0('~/symlinks/host/figures/20250330_Admixture-Fig1-K2-5_',sp,'.pdf'),adplot,height=2,width=6,dpi=600)

# K = 2-10, by geographic group 
geo_ord <- md %>% select(GeographicGroup) %>% distinct %>% mutate(geoord = as.numeric(gsub('GCC|GCO','',GeographicGroup))) %>% 
  arrange(geoord) %>% data.frame
admixmd <- admixmd %>% mutate(Specified_K = paste0('K',Specified_K))
k_ord <- admixmd %>%  select(Specified_K,K) %>% distinct %>% mutate(sp_k = as.numeric(gsub('K','',Specified_K))) %>% 
  arrange(sp_k) %>% data.frame
admixmd$GeographicGroup <- factor(admixmd$GeographicGroup,levels=geo_ord$GeographicGroup)
admixmd$Specified_K <- factor(admixmd$Specified_K,levels=unique(k_ord$Specified_K))
admixmd$K <- factor(admixmd$K,levels=unique(k_ord$K))
adplot =
  admixmd %>%  #specify the levels you want 
  ggplot(aes(x = factor(ID), y = Q, fill = factor(K))) +
  geom_col(size = 0.1) +
  facet_grid(Specified_K~GeographicGroup, scales = "free", space = "free") +
  theme_minimal(base_size=6) + labs(x = "",y = "") +
  scale_y_continuous(expand = c(0, 0),n.breaks = 3) +
  #scale_fill_manual(values=viridis(10,option='turbo'))+
  scale_fill_manual(values=brewer.pal(10,'Spectral'))+
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    #axis.text.x = element_blank(),
    axis.text.x=element_text(angle=90,size=5),
    axis.text.y = element_text(size=3),
    panel.grid = element_blank(),
    legend.position = 'bottom',
    plot.title = element_text(size=6)
  )
adplot

ggsave(paste0('~/symlinks/host/figures/20250318_Admixture-GeographicGroupK10_',sp,'.pdf'),adplot,height=8,width=6,dpi=600)


# K = 2-10, by geographic group 
egg_ord <- md %>% select(Egg,EggOrder) %>% distinct %>% arrange(EggOrder)
admixmd$Egg <- factor(admixmd$Egg,levels=egg_ord$Egg)
adplot =
  admixmd %>%  
  drop_na(Egg) %>% #specify the levels you want 
  ggplot(aes(x = factor(ID), y = Q, fill = factor(K))) +
  geom_col(size = 0.1) +
  facet_grid(Specified_K~Egg, scales = "free", space = "free") +
  theme_minimal(base_size=6) + labs(x = "",y = "") +
  scale_y_continuous(expand = c(0, 0),n.breaks = 3) +
  #scale_fill_manual(values=viridis(10,option='turbo'))+
  scale_fill_manual(values=brewer.pal(10,'Spectral'))+
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    #axis.text.x = element_blank(),
    axis.text.x=element_text(angle=90,size=5),
    axis.text.y = element_text(size=3),
    panel.grid = element_blank(),
    legend.position = 'bottom',
    plot.title = element_text(size=6)
  )
adplot

ggsave(paste0('~/symlinks/host/figures/20250318_Admixture-EggK10_',sp,'.pdf'),adplot,height=8,width=4,dpi=600)
                     
# Just K = 5, sorted by ancestry in each geo group
adplot =
  admixmd %>% filter(Specified_K == 5 | Specified_K == 2) %>%  #specify the levels you want 
  mutate(Specified_K = paste0('K',Specified_K)) %>% 
  ggplot(aes(x = factor(ID), y = Q, fill = factor(K))) +
  geom_col(size = 0.1) +
  facet_grid(AncestryA5+GeographicGroup~Specified_K, scales = "free", space = "free") +
  theme_minimal(base_size=6) + labs(x = "",y = "") +
  scale_y_continuous(expand = c(0, 0),n.breaks = 3) +
  scale_fill_manual(values=kcols)+
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    #axis.text.x = element_blank(),
    axis.text.x=element_text(angle=90,size=5),
    axis.text.y = element_text(size=3),
    panel.grid = element_blank(),
    legend.position = 'bottom',
    plot.title = element_text(size=6)
  )+
  coord_flip()
adplot

ggsave(paste0('~/symlinks/host/figures/20250408_Admixture-AncGeoK2-5_',sp,'.pdf'),adplot,height=12,width=6,dpi=600)
```

## Plot Tesselation + mtDNA Pies

```R
### Plot ADMIXTURE across landscape (tesselation)
setwd('~/EvoBioWolf/CUCKOO_gentes/population_genetics/')
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
library(ggnewscale)

# Hash out which species to run, entire script will run afterwards 
admix_run = 'autos_optatus_N50_LD'
sp='CO'
admix_run = 'autos_canorus_LD'
sp='CC'

qdir = 'admixture_q_files' #directory with Q files

admix = melt_admixture(prefix = admix_run, qdir = qdir)

#read in metadata
md = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt') %>% 
  select(ID,Haplogroup,Longitude,Latitude,Egg,GeographicGroup,HaplogroupColor) 
admixmd = left_join(admix,md) 

# Reorder individuals baseed on longitude
admixmd = admixmd %>% mutate(ID = fct_reorder(ID,Longitude))
if (sp == "CO") {
  kclust <- 'ACO'
  admixmd <- admixmd %>% mutate(K = gsub('K', kclust, K))
  kcols <- brewer.pal(7,'Purples')[c(1,3,5,6,7)]
  spshape=24
  filt="Cuculus optatus"
} else {
  kclust <- 'ACC'
  admixmd <- admixmd %>% mutate(K = gsub('K', kclust, K))
  kcols <- brewer.pal(6,'Greys')[c(2,3,4,5,6)]
  spshape=21
  filt="Cuculus canorus"
}


# We need to manually re-specify our K orders to reflect geography
admixmd <- admixmd %>% 
  mutate(K = case_when(
    # for optatus, swap 
    Specified_K == 5 & K == "ACO3" ~ "ACO5",
    Specified_K == 5 & K == "ACO5" ~ "ACO3",
    # for canorus, swap 
    Specified_K == 5 & K == "ACC5" ~ "ACC2",
    Specified_K == 5 & K == "ACC4" ~ "ACC5",
    Specified_K == 5 & K == "ACC2" ~ "ACC4",
    TRUE ~ K
  ))


#loop through all admixture runs and extract the average correlation values from evalAdmix, we want to MINIMIZE this! (closest to 0)
evaldat = NULL; for (Kval in seq(2,10,1)){
  r <- as.matrix(read.table(paste0("evalAdmix/",sp,"_eval_",Kval)))
  mean_value <- mean(r,na.rm=TRUE)
  median_value <- median(r,na.rm=TRUE)
  sd_value <- sd(r,na.rm=TRUE)
  iqr_value <- IQR(r,na.rm=TRUE)
  valdat = data.frame(K = Kval,mean = mean_value,median=median_value,sd=sd_value,iqr=iqr_value)
  evaldat = rbind(valdat,evaldat)
}

#plot, for main figure show the n=3 lowest median
targs = evaldat %>% slice_min(abs(median),n=3)
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

ggsave(paste0('~/symlinks/host/figures/20250318_',sp,'_evalAdmix.pdf'),ep,height=5,width=6,dpi=600)

# K = 5 
adplot =
  admixmd %>% filter(Specified_K == 5 ) %>%  #specify the levels you want 
  mutate(Specified_K = paste0('K',Specified_K)) %>% 
  ggplot(aes(x = factor(ID), y = Q, fill = factor(K))) +
  geom_col(size = 0.1) +
  facet_grid(Specified_K~., scales = "free", space = "free") +
  theme_minimal(base_size=6) + labs(x = "",y = "") +
  scale_y_continuous(expand = c(0, 0),n.breaks = 3) +
  scale_fill_manual(values=kcols)+
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    #axis.text.x = element_blank(),
    axis.text.x=element_text(angle=90,size=5),
    axis.text.y = element_text(size=3),
    panel.grid = element_blank(),
    legend.position = 'bottom',
    plot.title = element_text(size=6)
  )
adplot

ggsave(paste0('~/symlinks/host/figures/20250318_Admixture_',sp,'.pdf'),adplot,height=2.5,width=6,dpi=600)

# Assign K1 - K5 
ks = admixmd %>% filter(Specified_K == 5) %>% group_by(ID) %>% slice_max(Q) %>% dplyr::rename(KCluster=K) %>% select(ID,KCluster)

### Tesselation
# Import birdlife shapefiles, breeding == 2, presence ==1 means extant
bg <-  st_read('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/spatial/birdlife_breeding_distributions/SppDataRequest.shp')
filtered_data <- bg[bg$PRESENCE == 1 & bg$SEASONAL == 2 & bg$SCI_NAME == filt, ]

#which K value to plot 
show_k = 5

# Extract the q values and the long/lat 
tesselation_qs <- admixmd %>% filter(Specified_K == show_k) %>% select(ID,K,Q,Longitude,Latitude,Haplogroup) %>% 
  pivot_wider(names_from = K,values_from = Q) %>% select(ID, Longitude, Latitude, Haplogroup, paste0(kclust,seq(1,show_k,1)))
cluster_centroids <- left_join(tesselation_qs,ks) %>% group_by(KCluster) %>% 
  summarize(lat = mean(Latitude),lon=mean(Longitude))

#plot using ggplot
map.polygon <- getMap(resolution = "high")

# Jitter the lat/long and add haplogroup colors 
sitesp = st_as_sf(tesselation_qs %>% mutate(loj = jitter(Longitude,amount=1),laj = jitter(Latitude,amount=1)),remove = F, coords = c("loj", "laj"), crs = 4326, agr = "constant") 
max_col <- 5+show_k-1
pl = ggtess3Q(tesselation_qs[5:max_col], as.matrix(tesselation_qs[2:3]), map.polygon = map.polygon,col.palette = kcols)

# Plot the K1-5, look at color matching, etc 
check_cols = pl +
  geom_path(data = map.polygon, aes(x = long, y = lat, group = group),col='white',lwd=0.2) +
  coord_sf(xlim = c(min(admixmd$Longitude)-5, max(admixmd$Longitude)+5), 
           ylim = c(min(admixmd$Latitude)-5, max(admixmd$Latitude)+5), expand = FALSE)+    
  new_scale_fill()+
  new_scale_color()+
  geom_label(data = cluster_centroids, aes(x = lon, y = lat,label=KCluster,fill=KCluster))+
  scale_fill_manual(values=kcols)+
  xlab("Longitude") + ylab("Latitude") + 
  theme_classic(base_size=8)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"))+
  theme(legend.position = 'top',legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))
check_cols 

# Plot haplogroups 
k5p1 = pl +
  geom_path(data = map.polygon, aes(x = long, y = lat, group = group),col='white',lwd=0.2) +
  coord_sf(xlim = c(min(admixmd$Longitude)-5, max(admixmd$Longitude)+5), 
           ylim = c(min(admixmd$Latitude)-5, max(admixmd$Latitude)+5), expand = FALSE)+    
  new_scale_fill()+
  new_scale_color()+
  geom_point(data = sitesp, aes(x = loj, y = laj,fill=Haplogroup),shape=spshape, size = 2,col='black') +
  scale_fill_manual(values=md$HaplogroupColor,breaks=md$Haplogroup)+
  xlab("Longitude") + ylab("Latitude") + 
  theme_classic(base_size=8)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"))+
  theme(legend.position = 'top',legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))
k5p1

ggsave(paste0('~/symlinks/host/figures/20250318_Tesselation_',sp,'.pdf'),ggarrange(check_cols,k5p1,nrow=2),height=7,width=7,dpi=300)

##### plot all K... this will ANCESTRY Q WITHIN PIES!  #####
maximum_k = 5

for (kval in seq(2,maximum_k,1)) {
  
  cat('Working on K = ',kval,'\n')
  tesselation_qs <- admixmd %>% filter(Specified_K == kval) %>% select(ID,K,Q,Longitude,Latitude,Haplogroup) %>% 
    pivot_wider(names_from = K,values_from = Q) %>% select(ID, Longitude, Latitude, Haplogroup, paste0(kclust,seq(1,kval,1)))
  max_col <- 5 + kval - 1

  # First, grab the individuals and calculate the mean Q values within each cluster. Cluster will be geographic reigon, also calculate mean lat/long for plotting
  kept = admixmd %>% select(ID,Specified_K,K,Q,Latitude,Longitude,Group = GeographicGroup)
  group_summaries = kept %>% 
    filter(Specified_K == kval) %>%
    group_by(Group,K) %>%  #within each group and K, average lat/long/q and count number of individuals 
    summarize(Lat = mean(Latitude),
              Long = mean(Longitude),
              Q = mean(Q),
              N = n_distinct(ID)) %>% 
    ungroup() %>% #
    #Calculate scaling factors for the pies based on num samples
    mutate(MinN = min(N),
           MaxN = max(N)) %>%
    group_by(Group) %>%
    mutate(Scaling_factor = ((N - MinN) / (MaxN - MinN) * 10) + 2) %>%
    select(-MinN, -MaxN) 
  
  ##### Plot pies across the world 
  plot_pie <- function(data) {
    ggplot(data, aes(x = "", y = Q, fill = K,)) +
      geom_bar(col='white',lwd=0.5,width = 1, stat = "identity") +
      coord_polar("y") +
      scale_fill_manual(values=kcols)+
      theme_void() +
      theme(legend.position = "none")
  }
  
  #set up map and make a sf object from the summaries 
  sites = st_as_sf(group_summaries, coords = c("Long", "Lat"), crs = 4326, agr = "constant") 
  
  # Main map plot
  p = 
    ggtess3Q(tesselation_qs[5:max_col], as.matrix(tesselation_qs[2:3]), map.polygon = filtered_data,col.palette = kcols) + 
    geom_path(data = map.polygon, aes(x = long, y = lat, group = group)) +
    geom_sf(data = sites, aes(geometry = geometry), size = 0.1, alpha = 0.1, pch=26) +
    xlab('')+ylab('')+
    coord_sf(xlim = c(min(admixmd$Longitude)-5, max(admixmd$Longitude)+5), 
             ylim = c(min(admixmd$Latitude)-5, max(admixmd$Latitude)+5), expand = FALSE)+    
    theme_classic(base_size = 8)+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"))+
    theme(legend.position = 'top',legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))
  p
  
  # Add pies
  for (i in unique(group_summaries$Group)) {
    subset_data = group_summaries %>% filter(Group == i)
    lon = unique(subset_data$Long)
    lat = unique(subset_data$Lat)
    scale_factor = unique(subset_data$Scaling_factor)
    cat('Scaling factor is : ',scale_factor,' for group : ',i,'\n')
    pie = plot_pie(subset_data)
    p <- p + annotation_custom(ggplotGrob(pie), 
                               xmin = lon - scale_factor, 
                               xmax = lon + scale_factor, 
                               ymin = lat - scale_factor, 
                               ymax = lat + scale_factor)
  }
  p
  assign(paste0('t',kval),p)
  
}

t5
png(paste0('~/symlinks/host/figures/20250318_Tesselation_',sp,'_AllK.png'),res=600,units='in',height=7,width=10)
print(ggarrange(t2,t3,t4,t5,ncol=2,nrow=2))
dev.off()

##### plot all K... this will show HAPLOGROUPS within PIES!  #####
show_hap_k = 5 
tesselation_qs <- admixmd %>% filter(Specified_K == show_hap_k) %>% select(ID,K,Q,Longitude,Latitude,Haplogroup,GeographicGroup) %>% 
  pivot_wider(names_from = K,values_from = Q) %>% select(ID, Longitude, Latitude, Haplogroup, GeographicGroup, paste0(kclust,seq(1,show_hap_k,1)))
max_col <- 6 + show_hap_k - 1

# First, grab the individuals and calculate the mean Q values within each cluster. Cluster will be geographic reigon, also calculate mean lat/long for plotting
#count the proportion of each haplogroup within each distance group 
haps = tesselation_qs %>% count(GeographicGroup,Haplogroup) %>% ungroup %>% group_by(GeographicGroup) %>% 
  mutate(Total = sum(n),
         Proportion = n/Total,
         Percent = paste0(round(n/Total,3)*100,'% (',n,')'))
coords <- tesselation_qs %>% group_by(GeographicGroup) %>% 
  summarize(Latitude=mean(Latitude),
            Longitude=mean(Longitude))
group_summaries <- left_join(haps,coords) %>% left_join(., md %>% select(Haplogroup,HaplogroupColor) %>% unique) %>% 
  arrange(GeographicGroup,Haplogroup)

#we ALSO need to create a scaling factor, based on how many haps are in each $Distance region 
scal = group_summaries %>%
  select(GeographicGroup, Total) %>%
  unique() %>%
  ungroup() %>%
  mutate(Min = min(Total),
         Max = max(Total),
         Scaling_factor = ((Total - Min) / (Max - Min) * 10) + 2)
#make sure the scaling factor is linear
scal %>% ggplot(aes(x=Total,y=Scaling_factor))+geom_point()+theme_bw()
#add that scaling factor back to the haplogroups 
group_summaries = left_join(group_summaries,scal %>% select(GeographicGroup,Scaling_factor))

##### Plot pies across the world 
#set up map and make a sf object from the summaries 
sites = st_as_sf(group_summaries, coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant") 

# Main map plot
p = 
  ggtess3Q(tesselation_qs[6:max_col], as.matrix(tesselation_qs[2:3]), map.polygon = map.polygon,col.palette = kcols) + 
  geom_path(data = map.polygon, aes(x = long, y = lat, group = group)) +
  geom_sf(data = sites, aes(geometry = geometry), size = 0.1, alpha = 0.1, pch=26) +
  xlab('')+ylab('')+
  coord_sf(xlim = c(min(admixmd$Longitude)-5, max(admixmd$Longitude)+5), 
           ylim = c(min(admixmd$Latitude)-5, max(admixmd$Latitude)+5), expand = FALSE)+    
  theme_classic(base_size = 8)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"))+
  theme(legend.position = 'top',legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))
p

# Pie function 
plot_pie <- function(data) {
  ggplot(data, aes(x = "", y = n, fill = Haplogroup,label=paste0(GeographicGroup,': ',Total))) +
    geom_bar(col='black',lwd=0.5,width = 1, stat = "identity") +
    coord_polar("y") +
    scale_fill_manual(values=md$HaplogroupColor,breaks=md$Haplogroup)+
    theme_void() +
    theme(legend.position = "none")
}

# Add pies
for (i in unique(group_summaries$GeographicGroup)) {
  subset_data = group_summaries %>% filter(GeographicGroup == i)
  lon = unique(subset_data$Longitude)
  lat = unique(subset_data$Latitude)
  scale_factor = unique(subset_data$Scaling_factor)
  cat('Scaling factor is : ',scale_factor,' for group : ',i,'\n')
  pie = plot_pie(subset_data)
  p <- p + annotation_custom(ggplotGrob(pie), 
                             xmin = lon - scale_factor, 
                             xmax = lon + scale_factor, 
                             ymin = lat - scale_factor, 
                             ymax = lat + scale_factor)
  }
p

png(paste0('~/symlinks/host/figures/20250318_Tesselation-mtDNA-Pies_',sp,'_K5.png'),height=4,width=7,units='in',res=300)
p
dev.off()

final_cols <- data.frame(KCluster=paste0(kclust,seq(1,5,1)),Kcols=kcols)
final_cols %>% ggplot(aes(x = KCluster, y = 1, fill = Kcols)) +geom_tile(color = "black") + scale_fill_identity() + theme(legend.position = "none")

#simply assign based on K5 
write.table(left_join(ks,final_cols),file=paste0('~/symlinks/host/figures/20250318_',sp,'_assigned_k5.txt'),quote=F,sep='\t',row.names=F)

```

## Plot Spatial: Geography, Eggs, Ancestry

```R
#### Assign geographic distance groups with k-means
setwd('~/EvoBioWolf/CUCKOO_gentes/correlations_geography_genetics')
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
md = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt') %>% drop_na(GeographicGroup) %>% arrange(GeoOrder)
world <- map_data("world")

# Plot Geographic Groups 
as1 = md %>% mutate(LatJit = jitter(Latitude,amount = 1),
                    LonJit = jitter(Longitude,amount = 1)) 
ksites = st_as_sf(as1, coords = c("LonJit", "LatJit"), crs = 4326, agr = "constant") 
allsamp = ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),col='grey50',fill='grey95') +
  geom_sf(data = ksites, 
          aes(fill=GeographicGroup,shape=GeographicGroup),
          size=3,show.legend = T) +
  xlab('Longitude')+ylab('Latitude')+
  scale_fill_manual(values=md$GeoColor,breaks=md$GeographicGroup)+
  scale_shape_manual(values=md$GeoShape,breaks=md$GeographicGroup)+
  coord_sf(xlim = c(min(as1$Longitude)-5, max(as1$Longitude)+5), 
           ylim = c(min(as1$Latitude)-5, max(as1$Latitude)+5), expand = FALSE)+
  theme_classic()+
  facet_grid(SpeciesShort~.)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"),
        legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))+
  annotation_scale(line_width=0.5)+
  theme(legend.position='top')
allsamp

ggsave('~/symlinks/host/figures/20250319_GeographicGroups.pdf',allsamp,dpi=300,height=8,width=9)

# Plot haplogroups 
as1 = md %>% mutate(LatJit = jitter(Latitude,amount = 2),
                    LonJit = jitter(Longitude,amount = 2)) 

hap_colors <- md %>% drop_na(HaplogroupColor) %>% distinct(Haplogroup, HaplogroupColor) %>% deframe()
hap_shapes <- data.frame(names = names(hap_colors),shape = rep_len(c(21,22,23,24),11)) %>% deframe()

# re-plot
ksites = st_as_sf(as1, coords = c("LonJit", "LatJit"), crs = 4326)

plothaps = ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), col = 'grey50', fill = 'grey95') +
  geom_sf(data = ksites, aes(fill = Haplogroup, shape = Haplogroup), size = 2, show.legend = TRUE) +
  scale_fill_manual(values = hap_colors) +
  scale_shape_manual(values = hap_shapes) +
  coord_sf(xlim = c(min(as1$Longitude) - 5, max(as1$Longitude) + 5), 
           ylim = c(min(as1$Latitude) - 5, max(as1$Latitude) + 5), expand = FALSE) +
  facet_grid(SpeciesShort ~ .) +
  labs(x = 'Longitude', y = 'Latitude') +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.background = element_rect(fill = "aliceblue"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.key.size = unit(0.1, 'cm'),
        legend.position = 'top') +
  annotation_scale(line_width = 0.5)
plothaps
ggsave('~/symlinks/host/figures/20250408_HaplogroupLocations.pdf',plothaps,dpi=300,height=8,width=9)



# Plot Egg
md = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt') %>% drop_na(Egg) %>% arrange(EggOrder)
egg = md %>% mutate(LatJit = jitter(Latitude,amount = 2),
                    LonJit = jitter(Longitude,amount = 2)) 
eggsites = st_as_sf(egg, coords = c("LonJit", "LatJit"), crs = 4326, agr = "constant") 
allsamp = ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),col='grey50',fill='grey95') +
  geom_sf(data = eggsites, 
          aes(fill=Egg,shape=Egg),
          size=3,show.legend = T) +
  xlab('Longitude')+ylab('Latitude')+
  scale_fill_manual(values=md$EggCol,breaks=md$Egg)+
  scale_shape_manual(values=md$EggShape,breaks=md$Egg)+
  coord_sf(xlim = c(min(egg$Longitude)-5, max(egg$Longitude)+5), 
           ylim = c(min(egg$Latitude)-5, max(egg$Latitude)+5), expand = FALSE)+
  theme_classic()+
  facet_grid(SpeciesShort~.)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"),
        legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))+
  annotation_scale(line_width=0.5)+
  theme(legend.position='top')
allsamp

ggsave('~/symlinks/host/figures/20250319_EggLocations.pdf',allsamp,dpi=300,height=8,width=9)

# Plot Ancestry
md = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt') %>% drop_na(AncestryA5) %>% arrange(AncestryA5)
anc = md %>% mutate(LatJit = jitter(Latitude,amount = 2),
                    LonJit = jitter(Longitude,amount = 2)) 
ancsites = st_as_sf(anc, coords = c("LonJit", "LatJit"), crs = 4326, agr = "constant") 
plot_anc = ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),col='grey50',fill='grey95') +
  geom_sf(data = ancsites, 
          aes(fill=AncestryA5,shape=SpeciesShort),
          size=3,show.legend = T) +
  xlab('Longitude')+ylab('Latitude')+
  scale_fill_manual(values=md$AncestryColor,breaks=md$AncestryA5)+
  scale_shape_manual(values=md$Shape,breaks=md$SpeciesShort)+
  coord_sf(xlim = c(min(anc$Longitude)-5, max(anc$Longitude)+5), 
           ylim = c(min(anc$Latitude)-5, max(anc$Latitude)+5), expand = FALSE)+
  theme_classic()+
  facet_grid(SpeciesShort~.)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"),
        legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))+
  annotation_scale(line_width=0.5)+
  theme(legend.position='top')+
  guides(fill=guide_legend(override.aes=list(shape=21)))
plot_anc

ggsave('~/symlinks/host/figures/20250319_AncestryLocations.pdf',plot_anc,dpi=300,height=8,width=9)

```



## Distance Correlations 

First, identify the samples and groups for comparison (C. canorus first):

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

### Estimate FST

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

### Plot Results

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



### Repeat C. optatus

```R
#### Calculate pairwise geographic distance between distance groups 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/fst/2025apr_correlations/')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(viridis)
library(geosphere)
library(sf)
library(ggspatial)
library(RColorBrewer)

#Read in metadata
md = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt')

#groups for fst: distance ~ mtDNA
analyze_these = md %>% 
  filter(Analysis_PopulationGenetics == 1 & SpeciesShort == 'CO' & Sex == 'F') %>% select(GeographicGroup) %>%
  group_by(GeographicGroup) %>% mutate(DistanceHaps = n()) %>% unique %>% filter(DistanceHaps >= 3) %>% pull(GeographicGroup)

#calculate geographic distance between the groups
dists = md %>% filter(Analysis_PopulationGenetics == 1 & SpeciesShort == 'CO' & GeographicGroup %in% analyze_these) %>% group_by(GeographicGroup) %>% summarize(meanLat = mean(Latitude),meanLong = mean(Longitude))

#plot 
world = map_data("world")
sites = st_as_sf(dists, coords = c("meanLong", "meanLat"), 
                 crs = 4326, agr = "constant") 
fst_compars = ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),col='grey90',fill='white') +
  geom_sf(data = sites, 
          aes(fill=GeographicGroup,shape=GeographicGroup),
          size=4,alpha=0.9,show.legend = T,stroke=0.5) +
  scale_shape_manual(values=md$GeoShape,breaks=md$GeographicGroup)+
  scale_fill_manual(values=md$GeoColor,breaks=md$GeographicGroup)+
  xlab('')+ylab('')+
  coord_sf(xlim = c(min(dists$meanLong)-5, max(dists$meanLong)+5), 
           ylim = c(min(dists$meanLat)-5, max(dists$meanLat)+5), expand = FALSE)+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"))+
  theme(legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))+
  annotation_scale(line_width=0.5)

pdf('~/symlinks/host/figures/20250405_Geographic_Clusters_FST_Comparisons.pdf',height=6,width=9)
fst_compars
dev.off()

# initialize an empty data frame to store pairwise distances
pairwise_distances <- data.frame()
sub_data = dists

# calculate pairwise distances for the current 'W' GeographicGroup
for(i in 1:(nrow(sub_data) - 1)) {
  for(j in (i + 1):nrow(sub_data)) {
    
    point1 <- c(sub_data$meanLong[i], sub_data$meanLat[i])
    point2 <- c(sub_data$meanLong[j], sub_data$meanLat[j])
    
    distance_km <- distHaversine(point1, point2) / 1000  # convert to km
    
    # append the result to the pairwise_distances data frame
    pairwise_distances <- rbind(pairwise_distances, 
                                data.frame(P1 = sub_data$GeographicGroup[i], 
                                           P2 = sub_data$GeographicGroup[j],
                                           Distance_km = distance_km))
  }
}

pairwise_distances = pairwise_distances %>% mutate(Group = paste0(P1,'__',P2))

# show the calculated pairwise distances
write_tsv(pairwise_distances,file='~/EvoBioWolf/CUCKOO_gentes/correlations_geography_genetics/Pairwise_GeographicDistance_Km_CO_2025APR05.txt')

#write out populations
for (pop in unique(analyze_these)) {
  su = md %>% filter(Analysis_PopulationGenetics ==1 & GeographicGroup == pop)
  write.table(su$ID,file=paste0('populations/',pop,'.pop'),quote=F,sep='\t',row.names=F,col.names=F)
}

```

As above:

Comparisons:

```bash
awk '{print $4}' ~/EvoBioWolf/CUCKOO_gentes/correlations_geography_genetics/Pairwise_GeographicDistance_Km_CO_2025APR05.txt | sed '1d' > PairwiseComparisons.list
 cat PairwiseComparisons.list
GCO2__GCO3
GCO2__GCO5
GCO2__GCO6
GCO3__GCO5
GCO3__GCO6
GCO5__GCO6
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
# cat PairwiseComparisons.list | xargs -I {} sbatch -J FST_{} 2B.Pairwise_Distance_FST_CO.sh {} 
GROUP=$1

VCFS=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/manyhost_hunt/optatus/vcfs

for CHR in $(cat Chromosomes.list); do

echo "WORKING ON CHR: ${GROUP} and ${CHR}"

p1=$(echo ${GROUP} | sed 's/__.*//g')
p2=$(echo ${GROUP} | sed 's/.*__//g')

mkdir -p work out

# Calculate FST
~/modules/vcftools/bin/vcftools --gzvcf ${VCFS}/${CHR}.SNP.DP3.vcf.gz --max-missing 0.1 --fst-window-size 100000 --fst-window-step 100000 --weir-fst-pop populations/${p1}.pop --weir-fst-pop populations/${p2}.pop --out work/${CHR}_${GROUP}
awk -v g=${GROUP} '{OFS="\t"}{print $1, $2, $3,$4,$5, g}' work/${CHR}_${GROUP}.windowed.weir.fst > out/${CHR}_${GROUP}.fst

done
```

For mtDNA, calculate PHI using an AMOVA, subset pops:

```bash
for i in $(cat Populations.list); do samtools faidx chr_MT.SNP.DP3-AC1-MQ40-MM1.min4.fasta $(cat populations/${i}.pop) > fastas/${i}.fa; done
```

And repeat only for W:

```bash
for i in $(cat Populations.list); do samtools faidx chr_W.SNP.DP3-AC1-MQ40-MM1.min4.fasta $(grep '_F$' populations/${i}.pop) > fastas/${i}.W.fa; done
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

write_tsv(dats,file='PHI_chr_MT-W_2025APR05.txt')
```

Plot final correlations:

```R
#### Plot Distance ~ FST / PHIst for C. optatus
setwd('~/EvoBioWolf/CUCKOO_gentes/correlations_geography_genetics/')
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
md = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt')

d = read_tsv('Pairwise_GeographicDistance_Km_CO_2025APR05.txt')
f = read_tsv('Pairwise_FST_DistanceGroups_Autosomes_CO_2025APR05.txt')
p = read_tsv('Pairwise_PHIst_DistanceGroups_W-MT_CO_2025APR05.txt')

#correlation between W/MT
p %>% ggplot(aes(x=PHI_W,y=PHI))+geom_point()+theme_bw()
cor.test(p$PHI_W,p$PHI)

# data:  p$PHI_W and p$PHI
# t = 7.2041, df = 4, p-value = 0.001968
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.6971846 0.9961463
# sample estimates:
#   cor 
# 0.9635574 

#prep and merge frames 
p = p %>% dplyr::rename(chr_W = PHI_W, chr_MT = PHI) %>% pivot_longer(!Group) %>% 
  dplyr::rename(FST = value, chr = name) %>% mutate(FST = pmax(0, pmin(1, FST)))
names(f) = c('chr','start','end','snps','FST','Group')
amc = f
f2 = amc %>% mutate(FST = pmax(0, pmin(1, FST))) %>% select(!c(start,end,snps))
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
cors$AvZ = factor(cors$AvZ,levels=c('Autosome','Z','W','mtDNA'))
cors <- cors %>% arrange(AvZ)

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
# AvZ      estimate std.error statistic    df p.value conf.low conf.high   padj signif
# <fct>       <dbl>     <dbl>     <dbl> <dbl>   <dbl>    <dbl>     <dbl>  <dbl> <chr> 
#   1 Autosome    0.448    0.0363     12.3   2.01 0.00640    0.292     0.603 0.0256 *     
#   2 Z           0.391    0.113       3.47  1.44 0.115     -0.328     1.11  0.462  n.s.  
# 3 W           3.76     0.208      18.0   2.05 0.00275    2.88      4.63  0.0110 *     
#   4 mtDNA       1.10     0.339       3.25  1.14 0.165     -2.15      4.35  0.660  n.s. 


# Create the plot
fsize=2
pp1 = ggplot(dfs, aes(x = log(pmax(Distance_km,0.005)), y = log(pmax(mean,0.005)), color = AvZ,shape=AvZ,fill=AvZ)) +
  geom_point(size=1) + #0.25 for main plot 
  geom_smooth(method = "lm", se = TRUE,lwd=0.75,alpha=0.1) +  #0.5 for main plot 
  # geom_text(data = cors[1,], aes(x = Inf, y = -Inf, label = label), size=fsize, vjust = -6, hjust = 1.73, col = cols[1]) +
  # geom_text(data = cors[2,], aes(x = Inf, y = -Inf, label = label), size=fsize, vjust = -4.5, hjust = 1.73, col = cols[2]) +
  # geom_text(data = cors[3,], aes(x = Inf, y = -Inf, label = label), size=fsize, vjust = -3, hjust = 1.1, col = cols[3]) +
  # geom_text(data = cors[4,], aes(x = Inf, y = -Inf, label = label), size=fsize, vjust = -1.5, hjust = 2.02, col = cols[4]) +
  labs(title = "Relationship between Mean and Distance_km",
       x = "log(Geographic Distance)",
       y = "log(FST | ΦST)",
       color = "AvZ") +
  scale_color_brewer(palette='Dark2')+
  scale_fill_brewer(palette='Dark2')+
  facet_wrap(.~AvZ,scales='free',nrow=4,ncol=1)+ #for dual plot sensitivity
  scale_shape_manual(values=c(15,16,17,3))+
  theme_bw(base_size=6) + theme(legend.position = 'none')
pp1

pdf('~/symlinks/host/figures/20250405_Correlations_logFST-logDIST_CO.pdf',height=1.5,width=1.5)
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

write_tsv(cors,file='~/symlinks/host/figures/20250405_Correlations_logFST-logDIST_cortest-results_CO.txt')

pdf('~/symlinks/host/figures/20250405_Correlations-BothMethods.pdf',height=6,width=4)
ggarrange(pp1,pp2,common.legend = TRUE,nrow=1,ncol=2)
dev.off()

```



# Analyses: Egg Associations 

## dbRDA

This analysis uses a distance matrix of covariates as a response variable and egg phenotype as an explanatory variable, as dBRDA uses a matrix as input for response: (e.g. dist(mtDNA) ~ Egg_Phenotype) 

```R
#### dbRDA for e.g. egg type ~ maternal haplogroup associations 
setwd('~/EvoBioWolf/CUCKOO_gentes/gens_associations/')
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
mdf = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt') %>% arrange(ID)
mdf %>% filter(Analysis_FullGensAssociations == 1) %>% count(SpeciesShort,Sex) 

spdat <- NULL
for (sp in c('CC','CO')) {
  
  #only grab samples with known egg (hash out filter for related individuals)
  md_egg = mdf %>%
    filter(Analysis_FullGensAssociations == 1 & SpeciesShort == sp) %>%
    drop_na(Egg) %>%
    #filter(Sex == 'F') %>%  # hash out if you want to do mtDNA! 
    select(ID = ID, Host = HostParentShort, Egg, Habitat, Haplogroup, Ancestry = AncestryA5, Geography = GeographicGroup, Latitude, Longitude, CountryFull, HaplogroupColor,AncestryColor,GeoColor)
  
  ##### Initialize, summary stats on raw data #####
  #ensure they are all factor variables
  md_egg = md_egg %>% mutate_at(c('Host','Egg','Habitat','Haplogroup','Ancestry','Geography'),as.factor)
  
  # If you want to exclude the blue clades W1, W2, W3! 
  #md_egg <-  md_egg %>% filter(!grepl('W1|W2|W3',Haplogroup))
  md_egg = md_egg %>% group_by(Egg) %>% mutate(TotalEgg = n()) %>% ungroup
  minobs=2
  md_egg = md_egg %>% filter(TotalEgg >= minobs)
  md <- md_egg
  md %>% count(Egg)
  
  shp <- mdf %>% filter(SpeciesShort == sp) %>% select(Shape) %>% unique %>% pull(Shape)
  
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
  geo_p = ggplot(geo_mds_input,aes(x=V1,y=V2,fill=Geography)) + 
    geom_point(pch=shp,stroke=0.25) + theme_bw(base_size=6) + xlab('MDS1')+ylab('MDS2')+
    scale_fill_manual(values=md$GeoColor,breaks=md$Geography)+
    guides(fill=guide_legend(override.aes=list(shape=21)))+
    theme(legend.position = "top",
          legend.text = element_text(size = 4),
          legend.title = element_text(size = 4),legend.key.size = unit(0.05, 'cm'))
  geo_p
  
  #load in autosomal distance, calculated: ~/modules/VCF2Dis-1.50/bin/VCF2Dis -InPut autos_canorus_LD.vcf.gz -OutPut autos_canorus_LD.pdist 
  auto = read.table(paste0('~/EvoBioWolf/CUCKOO_gentes/gens_associations/dbrda/autos_',sp,'_LD.pdist'),header=F)
  #add proper names, because VCF2DIS truncates them
  auto = auto %>% mutate(IDNum = gsub('_.*','',V1))
  auto_id = left_join(auto,mdf %>% select(ID) %>% mutate(IDNum = gsub('_.*','',ID)))  
  ord = auto_id$ID
  autos = auto_id %>% select(!c(IDNum,ID,V1))
  names(autos) = ord
  rownames(autos) = ord
  
  #extract columns in same order as geographic distance
  common_names = intersect(rownames(autos), rownames(geo_mat))
  geo_order <- rownames(geo_mat)
  auto_aligned = autos[geo_order, geo_order] #subset 'auto' matrix to match the 'geo_mat' matrix rows and columns
  data.frame(mat1 = rownames(auto_aligned), mat2 = rownames(geo_mat)) %>% filter(mat1 != mat2) #sanity check for alignment 
  auto_mat = as.matrix(auto_aligned)
  
  #visualize autosomal distance matrix in terms of ancestry K = 5
  auto_mds = as.data.frame(cmdscale(auto_mat,2))
  auto_mds$ID = rownames(auto_mds)
  auto_mds_input = left_join(auto_mds,md)
  auto_p = ggplot(auto_mds_input,aes(x=V1,y=V2,fill=Ancestry)) + 
    geom_point(pch=shp,stroke=0.25) + theme_bw(base_size=6) + xlab('MDS1')+ylab('MDS2')+
    scale_fill_manual(values=md$AncestryColor,breaks=md$Ancestry)+
    guides(fill=guide_legend(override.aes=list(shape=21)))+
    theme(legend.position = "top",
          legend.text = element_text(size = 4),
          legend.title = element_text(size = 4),legend.key.size = unit(0.05, 'cm'))
  auto_p
  
  #calculate mtDNA distance (or W, hash out) 
  seq = read.dna(paste0('~/EvoBioWolf/CUCKOO_gentes/gens_associations/dbrda/chr_MT.SNP.DP3-AC1-MQ40-MM1.min4.fasta'),format='fasta')
  #seq = read.dna(paste0('~/EvoBioWolf/CUCKOO_gentes/gens_associations/dbrda/chr_W.SNP.DP3-AC1-MQ40-MM1.min4.fasta'),format='fasta')
  dna_dist = as.data.frame(as.matrix(dist.dna(seq,model='JC69')))
  
  #extract columns in same order as geographic distance
  dna_aligned = dna_dist[geo_order, geo_order] #subset 'auto' matrix to match the 'geo_mat' matrix rows and columns
  data.frame(mat1 = rownames(dna_aligned), mat2 = rownames(geo_mat)) %>% filter(mat1 != mat2) #sanity check for alignment 
  dna_mat = as.matrix(dna_aligned)
  
  #visualize the mtDNA distance matrix in terms of haplogroup 
  dna_mds = as.data.frame(cmdscale(dna_mat,2))
  dna_mds$ID = rownames(dna_mds)
  dna_mds_input = left_join(dna_mds,md)
  dna_p = ggplot(dna_mds_input,aes(x=V1,y=V2,fill=Haplogroup)) + 
    geom_point(pch=shp,stroke=0.25) + theme_bw(base_size=6) + xlab('MDS1')+ylab('MDS2')+
    scale_fill_manual(values=md$HaplogroupColor,breaks=md$Haplogroup)+
    guides(fill=guide_legend(override.aes=list(shape=21)))+
    theme(legend.position = "top",
          legend.text = element_text(size = 4),
          legend.title = element_text(size = 4),legend.key.size = unit(0.05, 'cm'))
  dna_p
  
  ggsave(paste0('~/symlinks/host/figures/20250330_dbRDA-',sp,'-Inputs-MT-Egg.pdf'),
         ggarrange(geo_p,auto_p,dna_p,nrow=1),height=1.5,width=4,dpi=300)
  # ggsave(paste0('~/symlinks/host/figures/20250330_dbRDA-',sp,'-Inputs-W-Egg.pdf'),
  #        ggarrange(geo_p,auto_p,dna_p,nrow=1),height=1.5,width=4,dpi=300)


  geoscat = as.data.frame(geo_mat) %>% mutate(ID = rownames(.)) %>% pivot_longer(!ID,names_to = 'ID2',values_to = 'EDistance') %>% mutate(Explanatory = 'Biogeographical')
  autoscat = as.data.frame(auto_mat) %>% mutate(ID = rownames(.)) %>% pivot_longer(!ID,names_to = 'ID2',values_to = 'EDistance') %>% mutate(Explanatory = 'Autosomal')
  mtscat = as.data.frame(dna_mat) %>% mutate(ID = rownames(.)) %>% pivot_longer(!ID,names_to = 'ID2',values_to = 'EDistance') %>% mutate(Explanatory = 'Matrilineal')
  
  ##### dbRDA: reverse #####
  covars = md %>% select(ID,Egg) %>% mutate_all(as.factor)
  inputs = c('dna_mat','auto_mat','geo_mat')
  dbr = list()
  for (inp in inputs) {
    
    #constrained ordination with step selection 
    null_formula_str = as.formula(paste(inp, "~ 1"))
    m1f = as.formula(paste(inp, "~ Egg"))
    m1 = dbrda(m1f, md, dist="gow",scaling=TRUE) # Model with all explanatory variables

    #by terms
    r1 = as.data.frame(anova(m1, by="terms", permu=10000)) # test for sign. environ. variables
    
    # adjusted R^2
    a1 = round(RsquareAdj(m1)$adj.r.squared,3)
    p1 = round(anova(m1)[1,4],3) # overall test of the significant of the analysis
    
    #save results
    lab = ifelse(inp == 'dna_mat','Haplogroup',ifelse(inp == 'auto_mat','Ancestry','Geography'))
    
    dbrda_results = rbind(r1) %>% drop_na(F) %>% mutate(adjR2 = c(a1), p = c(p1)) %>% dplyr::rename(anova_p = 'Pr(>F)') %>% 
      mutate(Response = rownames(.),Test = lab) 
    
    dbr[[lab]] = dbrda_results
  }
  
  dbrf = rbindlist(dbr) %>% as_tibble %>% mutate(Species = sp)
  spdat <- rbind(spdat,dbrf)
  
}

db_save =  spdat %>% select(-Df,SumOfSqs) %>% mutate(padj = p.adjust(p,method='bonferroni'))
db_save$Species <- factor(db_save$Species,levels=c('CC','CO'))
db_save$Test <- factor(db_save$Test,levels=c('Haplogroup','Ancestry','Geography'))
dbrda_pnt = db_save %>% 
  ggplot(aes(x=Species,y=adjR2,fill=Test))+
  geom_bar(stat='identity',position=position_dodge(width=0.9))+
  scale_fill_manual(values=c(brewer.pal(3,'Set2')))+
  ylab('Adjusted R2')+ylab('')+
  ylim(c(0,1.05))+
  theme_bw(base_size=6)+
  theme(legend.position='top')
dbrda_pnt

ggsave('~/symlinks/host/figures/20250330_dbRDA-MT.pdf',
       dbrda_pnt,height=2,width=1.75,dpi=300)
write.table(db_save,'~/symlinks/host/figures/20250330_dbRDA-MT_Results.txt',quote=F,sep='\t',row.names=F)

# ggsave('~/symlinks/host/figures/20250330_dbRDA-W.pdf',
#        dbrda_pnt,height=2,width=1.75,dpi=300)
# write.table(db_save,'~/symlinks/host/figures/20250330_dbRDA-W_Results.txt',quote=F,sep='\t',row.names=F)


```

### Sensitivity: Excluding Eggs/Haps

```bash
#### dbRDA for e.g. egg type ~ maternal haplogroup associations: sensitivity no E6/E1, no M1/M2/M3
setwd('~/EvoBioWolf/CUCKOO_gentes/gens_associations/')
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
mdf = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt') %>% arrange(ID)
mdf %>% filter(Analysis_FullGensAssociations == 1) %>% count(SpeciesShort,Sex) 

sp <- 'CC'

#only grab samples with known egg (hash out filter for related individuals)
md_egg = mdf %>%
  filter(Analysis_FullGensAssociations == 1 & SpeciesShort == sp) %>%
  drop_na(Egg) %>%
  #filter(Sex == 'F') %>%  # hash out if you want to do mtDNA! 
  #filter(SpeciesShort == 'CC' & !grepl('ECC6|ECC1',Egg)) %>% 
  filter(SpeciesShort == 'CC' & !grepl('MCC1|MCC2|MCC3',Haplogroup)) %>% 
  select(ID = ID, Host = HostParentShort, Egg, Habitat, Haplogroup, Ancestry = AncestryA5, Geography = GeographicGroup, Latitude, Longitude, CountryFull, HaplogroupColor,AncestryColor,GeoColor)

##### Initialize, summary stats on raw data #####
#ensure they are all factor variables
md_egg = md_egg %>% mutate_at(c('Host','Egg','Habitat','Haplogroup','Ancestry','Geography'),as.factor)

# If you want to exclude the blue clades W1, W2, W3! 
#md_egg <-  md_egg %>% filter(!grepl('W1|W2|W3',Haplogroup))
md_egg = md_egg %>% group_by(Egg) %>% mutate(TotalEgg = n()) %>% ungroup
minobs=2
md_egg = md_egg %>% filter(TotalEgg >= minobs)
md <- md_egg
md %>% count(Egg)

shp <- mdf %>% filter(SpeciesShort == sp) %>% select(Shape) %>% unique %>% pull(Shape)

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
geo_p = ggplot(geo_mds_input,aes(x=V1,y=V2,fill=Geography)) + 
  geom_point(pch=shp,stroke=0.25) + theme_bw(base_size=6) + xlab('MDS1')+ylab('MDS2')+
  scale_fill_manual(values=md$GeoColor,breaks=md$Geography)+
  guides(fill=guide_legend(override.aes=list(shape=21)))+
  theme(legend.position = "top",
        legend.text = element_text(size = 4),
        legend.title = element_text(size = 4),legend.key.size = unit(0.05, 'cm'))
geo_p

#load in autosomal distance, calculated: ~/modules/VCF2Dis-1.50/bin/VCF2Dis -InPut autos_canorus_LD.vcf.gz -OutPut autos_canorus_LD.pdist 
auto = read.table(paste0('~/EvoBioWolf/CUCKOO_gentes/gens_associations/dbrda/autos_',sp,'_LD.pdist'),header=F)
#add proper names, because VCF2DIS truncates them
auto = auto %>% mutate(IDNum = gsub('_.*','',V1))
auto_id = left_join(auto,mdf %>% select(ID) %>% mutate(IDNum = gsub('_.*','',ID)))  
ord = auto_id$ID
autos = auto_id %>% select(!c(IDNum,ID,V1))
names(autos) = ord
rownames(autos) = ord

#extract columns in same order as geographic distance
common_names = intersect(rownames(autos), rownames(geo_mat))
geo_order <- rownames(geo_mat)
auto_aligned = autos[geo_order, geo_order] #subset 'auto' matrix to match the 'geo_mat' matrix rows and columns
data.frame(mat1 = rownames(auto_aligned), mat2 = rownames(geo_mat)) %>% filter(mat1 != mat2) #sanity check for alignment 
auto_mat = as.matrix(auto_aligned)

#visualize autosomal distance matrix in terms of ancestry K = 5
auto_mds = as.data.frame(cmdscale(auto_mat,2))
auto_mds$ID = rownames(auto_mds)
auto_mds_input = left_join(auto_mds,md)
auto_p = ggplot(auto_mds_input,aes(x=V1,y=V2,fill=Ancestry)) + 
  geom_point(pch=shp,stroke=0.25) + theme_bw(base_size=6) + xlab('MDS1')+ylab('MDS2')+
  scale_fill_manual(values=md$AncestryColor,breaks=md$Ancestry)+
  guides(fill=guide_legend(override.aes=list(shape=21)))+
  theme(legend.position = "top",
        legend.text = element_text(size = 4),
        legend.title = element_text(size = 4),legend.key.size = unit(0.05, 'cm'))
auto_p

#calculate mtDNA distance (or W, hash out) 
seq = read.dna(paste0('~/EvoBioWolf/CUCKOO_gentes/gens_associations/dbrda/chr_MT.SNP.DP3-AC1-MQ40-MM1.min4.fasta'),format='fasta')
#seq = read.dna(paste0('~/EvoBioWolf/CUCKOO_gentes/gens_associations/dbrda/chr_W.SNP.DP3-AC1-MQ40-MM1.min4.fasta'),format='fasta')
dna_dist = as.data.frame(as.matrix(dist.dna(seq,model='JC69')))

#extract columns in same order as geographic distance
dna_aligned = dna_dist[geo_order, geo_order] #subset 'auto' matrix to match the 'geo_mat' matrix rows and columns
data.frame(mat1 = rownames(dna_aligned), mat2 = rownames(geo_mat)) %>% filter(mat1 != mat2) #sanity check for alignment 
dna_mat = as.matrix(dna_aligned)

#visualize the mtDNA distance matrix in terms of haplogroup 
dna_mds = as.data.frame(cmdscale(dna_mat,2))
dna_mds$ID = rownames(dna_mds)
dna_mds_input = left_join(dna_mds,md)
dna_p = ggplot(dna_mds_input,aes(x=V1,y=V2,fill=Haplogroup)) + 
  geom_point(pch=shp,stroke=0.25) + theme_bw(base_size=6) + xlab('MDS1')+ylab('MDS2')+
  scale_fill_manual(values=md$HaplogroupColor,breaks=md$Haplogroup)+
  guides(fill=guide_legend(override.aes=list(shape=21)))+
  theme(legend.position = "top",
        legend.text = element_text(size = 4),
        legend.title = element_text(size = 4),legend.key.size = unit(0.05, 'cm'))
dna_p

ggsave(paste0('~/symlinks/host/figures/20250330_dbRDA-',sp,'-Inputs-MT-Egg.pdf'),
       ggarrange(geo_p,auto_p,dna_p,nrow=1),height=1.5,width=4,dpi=300)
# ggsave(paste0('~/symlinks/host/figures/20250330_dbRDA-',sp,'-Inputs-W-Egg.pdf'),
#        ggarrange(geo_p,auto_p,dna_p,nrow=1),height=1.5,width=4,dpi=300)


geoscat = as.data.frame(geo_mat) %>% mutate(ID = rownames(.)) %>% pivot_longer(!ID,names_to = 'ID2',values_to = 'EDistance') %>% mutate(Explanatory = 'Biogeographical')
autoscat = as.data.frame(auto_mat) %>% mutate(ID = rownames(.)) %>% pivot_longer(!ID,names_to = 'ID2',values_to = 'EDistance') %>% mutate(Explanatory = 'Autosomal')
mtscat = as.data.frame(dna_mat) %>% mutate(ID = rownames(.)) %>% pivot_longer(!ID,names_to = 'ID2',values_to = 'EDistance') %>% mutate(Explanatory = 'Matrilineal')

##### dbRDA: reverse #####
covars = md %>% select(ID,Egg) %>% mutate_all(as.factor)
inputs = c('dna_mat','auto_mat','geo_mat')
dbr = list()
for (inp in inputs) {
  
  #constrained ordination with step selection 
  null_formula_str = as.formula(paste(inp, "~ 1"))
  m1f = as.formula(paste(inp, "~ Egg"))
  m1 = dbrda(m1f, md, dist="gow",scaling=TRUE) # Model with all explanatory variables
  
  #by terms
  r1 = as.data.frame(anova(m1, by="terms", permu=10000)) # test for sign. environ. variables
  
  # adjusted R^2
  a1 = round(RsquareAdj(m1)$adj.r.squared,3)
  p1 = round(anova(m1)[1,4],3) # overall test of the significant of the analysis
  
  #save results
  lab = ifelse(inp == 'dna_mat','Haplogroup',ifelse(inp == 'auto_mat','Ancestry','Geography'))
  
  dbrda_results = rbind(r1) %>% drop_na(F) %>% mutate(adjR2 = c(a1), p = c(p1)) %>% dplyr::rename(anova_p = 'Pr(>F)') %>% 
    mutate(Response = rownames(.),Test = lab) 
  
  dbr[[lab]] = dbrda_results
}

dbrf = rbindlist(dbr) %>% as_tibble %>% mutate(Filter = 'WithoutECC1_ECC6')
dbrf2 = rbindlist(dbr) %>% as_tibble %>% mutate(Filter = 'WithoutMCC1_MCC2_MCC3')
dats <- rbind(dbrf,dbrf2)

db_save =  dats %>% select(-Df,SumOfSqs) %>% mutate(padj = p.adjust(p,method='bonferroni'))
db_save$Test <- factor(db_save$Test,levels=c('Haplogroup','Ancestry','Geography'))
dbrda_pnt = db_save %>% 
  ggplot(aes(x=Filter,y=adjR2,fill=Test))+
  geom_bar(stat='identity',position=position_dodge(width=0.9))+
  scale_fill_manual(values=c(brewer.pal(3,'Set2')))+
  ylab('Adjusted R2')+ylab('')+
  ylim(c(0,1.05))+
  theme_bw(base_size=6)+
  theme(legend.position='top')
dbrda_pnt

ggsave('~/symlinks/host/figures/20250403_dbRDA-Sensitivity.pdf',
       dbrda_pnt,height=2,width=1.75,dpi=300)
write.table(db_save,'~/symlinks/host/figures/20250403_dbRDA-Sensitivity_Results.txt',quote=F,sep='\t',row.names=F)

```



## MNLR

Create a MNLR with egg or host or habitat as the response variable and Geography + Autosomal K + Haplogroups as the predictors. In short:

* Only retain response variables where there are at least n=2 observations
* Downsample all response classes so that all classes have n=2 observations
* Fit 7 multinomial logistic regression models, each with n=100 bootstraps using all combinations of predictors
* Extract AUC, and use the model to predict response variable on the full dataset again (too small for unseen data prediction)
* Repeat the above procedure 100 times so that different downsampled observations are included 
* Determine which classes are predicted correctly (% correct) from the confusion matrix on real / predicted responses across bootstraps

Run model:

```R
#!/usr/bin/env rscript
args <- commandArgs(trailingOnly = TRUE)

# Egg associations, MNLR, Model Selection c optatus
#### Find associations between MT/W haplotypes and features
setwd('~/EvoBioWolf/CUCKOO_gentes/gens_associations/mnlr/')
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

sp = args[1]

set.seed(123)
md = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt')

#only grab samples with known egg (hash out filter for related individuals)
md_egg = md %>%
  filter(SpeciesShort == sp & Analysis_FullGensAssociations == 1) %>%
  select(ID = ID, Egg, Haplogroup, Geography = GeographicGroup,Ancestry = AncestryA5)
md_egg %>% count(Egg)

md_egg = md_egg %>% group_by(Egg) %>% mutate(TotalEgg = n()) %>% ungroup
minobs=2
md_egg = md_egg %>% filter(TotalEgg >= minobs)
md_egg %>% count(Egg)

##### Initialize, summary stats on raw data #####
#ensure they are all factor variables
md_egg = md_egg %>% mutate_at(c('Egg','Haplogroup','Ancestry','Geography'),as.factor)

#### Count proportions first, count proportions for Egg and habitat and egg
ep = md_egg %>% group_by(Egg) %>% mutate(Total = n()) %>%
  pivot_longer(c(Ancestry,Geography,Haplogroup)) %>% ungroup %>%
  group_by(Egg,name,value) %>%
  summarize(Parts = n(),
            Proportion = Parts / Total,
            Percent = paste0(round( Parts / Total,3)*100,'% (',Total,')')) %>% unique

# Bind them together
ap = rbind(ep %>% ungroup %>% mutate(Response = Egg, variable = 'Egg') %>% select(-Egg))

# Just for ordering the covariates nicely
# Plot proportions
ord <- ap %>% select(name, value) %>%
  distinct() %>%
  mutate(ord = as.numeric(gsub("[^0-9.]", "", value))) %>%
  arrange(name, ord)
ap$value <- factor(ap$value,levels=ord$value)
egglev <- md %>% select(Egg,EggOrder) %>% distinct %>% arrange(EggOrder)
ap$Response <- factor(ap$Response,levels=egglev$Egg)
app = ap %>%
  arrange(value) %>%
  ggplot(aes(y=Response,x=value,fill=Proportion,label=Parts))+
  facet_grid(variable~name,scales='free',space='free')+
  scale_fill_gradient('Proportion Observed',low='yellow',high='red')+
  geom_tile()+
  geom_text(size=1)+
  ylab('')+xlab('')+
  theme_bw(base_size=6)+
  theme(legend.position = 'top')
app

ggsave(paste0('~/symlinks/host/figures/20250330_',sp,'_MNLR_Proportions.pdf'),app,height=3,width=7,dpi=300)

##### Model Selection #####
#assess covariate importance with model selection, using MNLR
vars = 'Egg'

# With bootstrapping
ctrl = trainControl(method = "boot",   # Use bootstrapping
                    number = 100,      # Replicates
                    summaryFunction = multiClassSummary,  #caret attempts to stratify levels
                    classProbs = TRUE, #For classification models, to save class probabilities
                    savePredictions = TRUE) #save predictions for each fold

# Change punctuation e.g. 'A. pal' to A_pal' for host fork
md_cv = md_egg %>% mutate(Egg = gsub('\\. ','_',Egg))

# For determining which predictors improve model:
adat = NULL; conf_dat = NULL; class_dat = NULL; new_preds = NULL; counter = 0

# Filter at the very base level, ensuring that across egg / Egg / habitat we have the same individuals with representative sampling
md_subbed = md_cv %>% mutate_at(c('Egg','Haplogroup','Ancestry','Geography'),as.factor)

for (rep in seq(1,100,1)){  # Create 100 replicate models
  for (var in vars) { counter = counter + 1;

  # Ensure that we have adequate levels, only a sanity since it is already filtered
  retained = md_subbed %>% ungroup %>% count(!!sym(var)) %>% filter(n >= 2) %>% pull(var)
  length(retained)
  # Number of samples to subsample
  subsamp = md_subbed %>% ungroup %>% count(!!sym(var)) %>% filter(n >= 2) %>% summarize(min = min(n)) %>% pull(min)

  cat('Downsampling to n = ',subsamp,', requiring min = ',2,' for variable: ',var,', ','replicate: ',rep,'\n')
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
    dat = data.frame(Model = model, Iteration = rep, Variable = var, Subsampled = subsamp, MinObs = 2,
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
      mutate(Model = model, Iteration = counter, Variable = var, Subsample = subsamp, MinObs = 2, AUC=dat_best$AUC,logloss=dat_best$logLoss,Accuracy = dat_best$Accuracy,AccuracySD=dat_best$AccuracySD)
    new_preds = rbind(new_preds,conf_real)
    rm(conf_real,dat,dat_best)

  } # Exit model loop
  } # Exit response variable loop
} # Exit iteration loop

# Write the no-blue data
write.table(adat,paste0('20250330_Model_Selection_Boot-2Obs-100Reps_',sp,'.txt'),quote=F,sep='\t',row.names=F)
write.table(new_preds,paste0('20250330_ConfusionMatrix_Boot-2Obs-100Reps_',sp,'.txt'),quote=F,sep='\t',row.names=F)


```

### Plot

```bash
# Egg associations, MNLR, Model Selection c optatus
#### Find associations between MT/W haplotypes and features 
setwd('~/EvoBioWolf/CUCKOO_gentes/gens_associations/mnlr/')
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

sp = 'CC'
set.seed(123)
md = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt')

#only grab samples with known egg (hash out filter for related individuals)
md_egg = md %>%
  filter(SpeciesShort == sp & Analysis_FullGensAssociations == 1) %>%
  select(ID = ID, Egg, Haplogroup, Geography = GeographicGroup,Ancestry = AncestryA5)
md_egg %>% count(Egg)

md_egg = md_egg %>% group_by(Egg) %>% mutate(TotalEgg = n()) %>% ungroup
minobs=2
md_egg = md_egg %>% filter(TotalEgg >= minobs)
md_egg %>% count(Egg)

##### Initialize, summary stats on raw data #####
#ensure they are all factor variables
md_egg = md_egg %>% mutate_at(c('Egg','Haplogroup','Ancestry','Geography'),as.factor)

#### Count proportions first, count proportions for Egg and habitat and egg
ep = md_egg %>% group_by(Egg) %>% mutate(Total = n()) %>%
  pivot_longer(c(Ancestry,Geography,Haplogroup)) %>% ungroup %>%
  group_by(Egg,name,value) %>%
  summarize(Parts = n(),
            Proportion = Parts / Total,
            Percent = paste0(round( Parts / Total,3)*100,'% (',Total,')')) %>% unique

# Bind them together
ap = rbind(ep %>% ungroup %>% mutate(Response = Egg, variable = 'Egg') %>% select(-Egg))

# Just for ordering the covariates nicely
# Plot proportions
ord <- ap %>% select(name, value) %>%
  distinct() %>%
  mutate(ord = as.numeric(gsub("[^0-9.]", "", value))) %>%
  arrange(name, ord)
ap$value <- factor(ap$value,levels=ord$value)
egglev <- md %>% select(Egg,EggOrder) %>% distinct %>% arrange(EggOrder)
ap$Response <- factor(ap$Response,levels=egglev$Egg)
app = ap %>%
  arrange(value) %>%
  ggplot(aes(y=Response,x=value,fill=Proportion,label=Parts))+
  facet_grid(variable~name,scales='free',space='free')+
  scale_fill_gradient('Proportion Observed',low='yellow',high='red')+
  geom_tile()+
  geom_text(size=1)+
  ylab('')+xlab('')+
  theme_bw(base_size=6)+
  theme(legend.position = 'top')
app

ggsave(paste0('~/symlinks/host/figures/20250330_',sp,'_MNLR_Proportions.pdf'),app,height=3,width=7,dpi=300)

# Read in saved data 
adat = read_tsv(paste0('20250330_Model_Selection_Boot-2Obs-100Reps_',sp,'.txt'))
conf_mats = read_tsv(paste0('20250330_ConfusionMatrix_Boot-2Obs-100Reps_',sp,'.txt'))

# Plot results 
leg = data.frame(Model = c('m1','m2','m3','m4','m5','m6','m7'),
                 Label = c('A+G+M','G+M','M','A+M','A','A+G','G'))
model_dat = adat %>% left_join(.,leg)

# Plot AUC for each model across all the iterations
model_dat$Label = factor(model_dat$Label, levels=c('A+G+M','A+G','A+M','G+M','A','G','M'))
auc_plot = model_dat %>% ggplot(aes(y=Label,x=AUC,fill=Variable))+
  geom_boxplot()+
  scale_fill_manual(values=brewer.pal(3,'RdBu'))+
  xlab('AUC')+ylab('Model Covariates')+
  theme_bw(base_size=6)+theme(legend.position = 'top',strip.text.y.right = element_text(angle = 0))
auc_plot

ggsave(paste0('~/symlinks/host/figures/20250330_',sp,'_MNLR_ModelSelection_AUC.pdf'),
       auc_plot,height=3,width=7,dpi=300)


# Summarize AUC across the core 3 models 
auc_plot_input <- model_dat %>%
  # %>% 
  group_by(Label,Variable) %>% 
  sum_stats(AUC)

# Plot
cols <- brewer.pal(3,'Set2')[c(1,2,3)]
#model_dat$Label <- factor(model_dat$Label,levels=c('A','G','M'))

auc_summary_plot <- auc_plot_input %>% 
  filter(Label == 'A' | Label == 'G' | Label == 'M') %>% 
  ggplot(aes(y=Label,fill=Label,x=mean,xmin=conf_low,xmax=conf_high))+
  geom_bar(stat='identity',position=position_dodge(width=0.9))+
  scale_fill_manual(values=cols)+
  xlab('')+ylab('AUC (Mean & 95% CI)')+
  geom_errorbar(width=0.25,position=position_dodge(width=0.9))+
  theme_bw(base_size=6)+
  coord_cartesian(xlim=c(0.5,1))+
  theme(legend.position='top')
auc_summary_plot

ggsave(paste0('~/symlinks/host/figures/20250330_',sp,'_MNLR_ModelSelection95CI_AUC.pdf'),
       auc_summary_plot,height=2,width=1.5,dpi=300)

write.table(auc_plot_input,file=paste0('~/symlinks/host/figures/20250330_',sp,'_MNLR_ModelSelection_AUCResults.txt'),quote=F,sep='\t',row.names=F)

#order full, single plot, make sure the 3 variables are in order 
egglev <- md %>% select(Egg,EggOrder) %>% distinct %>% arrange(EggOrder)
conf_mats = conf_mats %>% mutate(Prediction = factor(Prediction,levels=egglev$Egg),
                                 Reference = factor(Reference,levels=egglev$Egg))

### Plot how the addition of haplogroup improves predictions show A+G (m6) vs A+G+M (m1)
auc_vals <- adat %>% group_by(Model) %>% sum_stats(AUC)

### only geography (G; model 7) 
lab <- auc_vals %>% filter(Model == 'm7') %>% mutate(label = paste0('G: ',round(conf_low,2),' - ',round(conf_high,2)))
repredictions_geo = conf_mats %>% 
  filter(Model == 'm7') %>%  # (G ONLY) 
  group_by(Variable,MinObs,Iteration,Reference) %>% 
  mutate(TotalIteration = sum(Freq), #within each of the 100 replicates, calculate the frequency of correct classifications
         Proportion = Freq / TotalIteration) %>% ungroup %>% 
  group_by(Prediction,Reference,Variable) %>% 
  sum_stats(Proportion)
geo_plot = repredictions_geo %>% 
  ggplot(aes(x=Prediction,y=Reference,fill=mean,col=sd))+
  geom_tile(col=NA)+
  ggtitle(lab$label)+
  facet_wrap(lab$label~.,scales='free')+
  scale_fill_continuous(low='white',high='darkblue')+  theme_bw(base_size=6)+
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1),legend.position='top')
geo_plot


### ancestry + geography (A+G; m6)
lab <- auc_vals %>% filter(Model == 'm6') %>% mutate(label = paste0('A+G: ',round(conf_low,2),' - ',round(conf_high,2)))
repredictions_ancgeo = conf_mats %>% 
  filter(Model == 'm6') %>%  # plot model 6 (A+G) 
  group_by(Variable,MinObs,Iteration,Reference) %>% 
  mutate(TotalIteration = sum(Freq), #within each of the 100 replicates, calculate the frequency of correct classifications
         Proportion = Freq / TotalIteration) %>% ungroup %>% 
  group_by(Prediction,Reference,Variable) %>% 
  sum_stats(Proportion)
ancgeo_plot = repredictions_ancgeo %>% 
  ggplot(aes(x=Prediction,y=Reference,fill=mean,col=sd))+
  geom_tile(col=NA)+
  ggtitle(lab$label)+
  facet_wrap(lab$label~.,scales='free')+
  scale_fill_continuous(low='white',high='darkblue')+  theme_bw(base_size=6)+
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1),legend.position='top')
ancgeo_plot

# ancestry + geogrpahy + hap (A+G+M, m1)
lab <- auc_vals %>% filter(Model == 'm1') %>% mutate(label = paste0('A+G+M: ',round(conf_low,2),' - ',round(conf_high,2)))
repredictions_ancgeohap = conf_mats %>% 
  filter(Model == 'm1') %>%  # plot model 1 (A+G+M) 
  group_by(Variable,MinObs,Iteration,Reference) %>% 
  mutate(TotalIteration = sum(Freq), #within each of the 100 replicates, calculate the frequency of correct classifications
         Proportion = Freq / TotalIteration) %>% ungroup %>% 
  group_by(Prediction,Reference,Variable) %>% 
  sum_stats(Proportion)
ancgeohap_plot = repredictions_ancgeohap %>% 
  ggplot(aes(x=Prediction,y=Reference,fill=mean,col=sd))+
  geom_tile(col=NA)+
  ggtitle(lab$label)+
  facet_wrap(lab$label~.,scales='free')+
  scale_fill_continuous(low='white',high='darkblue')+  theme_bw(base_size=6)+
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1),legend.position='top')
ancgeohap_plot

ggarrange(ancgeo_plot,ancgeohap_plot,nrow = 2,common.legend = TRUE)
ggsave(paste0('~/symlinks/host/figures/20250330_MNLR_ConfusionMatrix-Repredictions-',sp,'_M1vsM6.pdf'),
       ggarrange(ancgeo_plot,ancgeohap_plot,nrow = 2,common.legend = TRUE),
       #dpi=300,height=3,width=1) # optatus 
       dpi=300,height=3.5,width=1.5) # canorus

```



### Sensitivity: Egg Exclusions

```bash
# Egg associations, MNLR, Model Selection c optatus
#### Find associations between MT/W haplotypes and features 
setwd('~/EvoBioWolf/CUCKOO_gentes/gens_associations/mnlr/')
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

sp = 'CC'

set.seed(123)
md = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt') 

#only grab samples with known egg (hash out filter for related individuals)
md_egg = md %>%
  filter(SpeciesShort == sp & Analysis_GensAssociations == 1) %>%
  select(ID = ID, Egg, Haplogroup, Geography = GeographicGroup,Ancestry = AncestryA5)
md_egg %>% count(Egg)

md_egg = md_egg %>% group_by(Egg) %>% mutate(TotalEgg = n()) %>% ungroup
minobs=2
md_egg = md_egg %>% filter(TotalEgg >= minobs)
md_egg %>% count(Egg)

# For one CC fork, also drop the most abundant egg types (E1 and E6) to see how the results are impacted
md_egg <- md_egg %>% filter(!Egg %in% c('ECC1','ECC6'))

##### Initialize, summary stats on raw data #####
#ensure they are all factor variables
md_egg = md_egg %>% mutate_at(c('Egg','Haplogroup','Ancestry','Geography'),as.factor)

##### Model Selection #####
#assess covariate importance with model selection, using MNLR

# With bootstrapping
ctrl = trainControl(method = "boot",   # Use bootstrapping
                    number = 100,      # Replicates 
                    summaryFunction = multiClassSummary,  #caret attempts to stratify levels
                    classProbs = TRUE, #For classification models, to save class probabilities
                    savePredictions = TRUE) #save predictions for each fold

# For determining which predictors improve model:
adat = NULL; conf_dat = NULL; class_dat = NULL; new_preds = NULL; counter = 0

# Filter at the very base level, ensuring that across egg / Egg / habitat we have the same individuals with representative sampling
md_subbed = md_egg %>% mutate_at(c('Egg','Haplogroup','Ancestry','Geography'),as.factor)
var = 'Egg'
for (rep in seq(1,100,1)){  # Create 10 replicate models
  counter = counter + 1;
  
  # Ensure that we have adequate levels, only a sanity since it is already filtered
  retained = md_subbed %>% ungroup %>% count(!!sym(var)) %>% filter(n >= 2) %>% pull(var)
  length(retained)
  # Number of samples to subsample
  subsamp = md_subbed %>% ungroup %>% count(!!sym(var)) %>% filter(n >= 2) %>% summarize(min = min(n)) %>% pull(min)
  
  cat('Downsampling to n = ',subsamp,', requiring min = ',2,' for variable: ',var,', ','replicate: ',rep,'\n')
  # Subsampling
  mdi = md_subbed %>%
    filter(!!sym(var) %in% retained) %>%
    group_by(!!sym(var)) %>%
    sample_n(min(n(), subsamp),replace = TRUE) %>%
    ungroup()
  mdi %>% count(!!sym(var))
  
  # Ensure the factor levels are dropped which aren't in the dataset after filtering
  mdi = droplevels(mdi)
  
  # Function to safely train models
  safe_train <- function(formula, data, method, trControl, metric) {
    tryCatch({
      train(formula, data = data, method = method, trControl = trControl, metric = metric, trace = FALSE)
    }, error = function(e) {
      message(paste("Model failed:", as.character(formula), "Error:", e$message))
      return(NULL)
    })
  }
  
  # Training models inside loop
  m1 = safe_train(as.formula(paste(var, "~ Haplogroup + Ancestry + Geography")), mdi, "multinom", ctrl, "AUC")
  m2 = safe_train(as.formula(paste(var, "~ Haplogroup + Geography")), mdi, "multinom", ctrl, "AUC")
  m3 = safe_train(as.formula(paste(var, "~ Haplogroup")), mdi, "multinom", ctrl, "AUC")
  m4 = safe_train(as.formula(paste(var, "~ Haplogroup + Ancestry")), mdi, "multinom", ctrl, "AUC")
  m5 = safe_train(as.formula(paste(var, "~ Ancestry")), mdi, "multinom", ctrl, "AUC")
  m6 = safe_train(as.formula(paste(var, "~ Ancestry + Geography")), mdi, "multinom", ctrl, "AUC")
  m7 = safe_train(as.formula(paste(var, "~ Geography")), mdi, "multinom", ctrl, "AUC")
  
  models = list(m1, m2, m3, m4, m5, m6, m7)
  model_names = c("m1", "m2", "m3", "m4", "m5", "m6", "m7")
  
  for (i in seq_along(models)) {
    mo = models[[i]]
    
    if (is.null(mo)) {
      message(paste("Skipping", model_names[i], "due to failure"))
      next
    }
    
    # Proceed with extracting metrics if model trained successfully
    final_model = mo$finalModel
    AIC = AIC(final_model)
    dat = data.frame(Model = model_names[i], Iteration = rep, Variable = var, Subsampled = subsamp, MinObs = 2,
                     decay = mo$results$decay, logLoss = mo$results$logLoss, logLossSD = mo$results$logLossSD, 
                     AIC = AIC, AUC = mo$results$AUC, Accuracy = mo$results$Accuracy, 
                     AccuracySD = mo$results$AccuracySD, Kappa = mo$results$Kappa)
    
    dat_best = dat %>% slice_min(logLoss)
    adat = rbind(adat, dat_best)
    
    # Confusion matrix
    pred = mo$pred
    pred_best = pred %>% filter(decay == dat_best$decay)
    predicted_classes = predict(mo, newdata = mdi, type = "raw")
    new = mdi %>% select(reference = !!sym(var)) %>% mutate(predicted = predicted_classes)
    conf_new = confusionMatrix(new$predicted, new$reference)
    
    conf_real = as.data.frame(conf_new$table) %>% mutate(
      Prediction = gsub('_','\\. ',Prediction),
      Reference = gsub('_','\\. ',Reference)) %>%
      mutate(Model = model_names[i], Iteration = counter, Variable = var, 
             Subsample = subsamp, MinObs = 2, AUC=dat_best$AUC,
             logloss=dat_best$logLoss, Accuracy = dat_best$Accuracy, 
             AccuracySD=dat_best$AccuracySD)
    
    new_preds = rbind(new_preds, conf_real)
    
  } # Exit model loop
} # Exit iteration loop

# Write the no-blue data
write.table(adat,paste0('20250320_Model_Selection_Boot-2Obs_',sp,'-NoE1E6.txt'),quote=F,sep='\t',row.names=F)
write.table(new_preds,paste0('20250320_ConfusionMatrix_Boot-2Obs_',sp,'-NoE1E6.txt'),quote=F,sep='\t',row.names=F)





######## SECOND
set.seed(123)
md = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt') 

#only grab samples with known egg (hash out filter for related individuals)
md_egg = md %>%
  filter(SpeciesShort == sp & Analysis_GensAssociations == 1) %>%
  select(ID = ID, Egg, Haplogroup, Geography = GeographicGroup,Ancestry = AncestryA5)
md_egg %>% count(Egg)

md_egg = md_egg %>% group_by(Egg) %>% mutate(TotalEgg = n()) %>% ungroup
minobs=2
md_egg = md_egg %>% filter(TotalEgg >= minobs)
md_egg %>% count(Egg)

# For one CC fork, also drop the most abundant egg types (E1 and E6) to see how the results are impacted
md_egg <- md_egg %>% filter(!Haplogroup %in% c('MCC1','MCC2','MCC3'))

##### Initialize, summary stats on raw data #####
#ensure they are all factor variables
md_egg = md_egg %>% mutate_at(c('Egg','Haplogroup','Ancestry','Geography'),as.factor)

##### Model Selection #####
#assess covariate importance with model selection, using MNLR

# With bootstrapping
ctrl = trainControl(method = "boot",   # Use bootstrapping
                    number = 100,      # Replicates 
                    summaryFunction = multiClassSummary,  #caret attempts to stratify levels
                    classProbs = TRUE, #For classification models, to save class probabilities
                    savePredictions = TRUE) #save predictions for each fold

# For determining which predictors improve model:
adat = NULL; conf_dat = NULL; class_dat = NULL; new_preds = NULL; counter = 0

# Filter at the very base level, ensuring that across egg / Egg / habitat we have the same individuals with representative sampling
md_subbed = md_egg %>% mutate_at(c('Egg','Haplogroup','Ancestry','Geography'),as.factor)
var = 'Egg'
for (rep in seq(1,100,1)){  # Create 10 replicate models
  counter = counter + 1;
  
  # Ensure that we have adequate levels, only a sanity since it is already filtered
  retained = md_subbed %>% ungroup %>% count(!!sym(var)) %>% filter(n >= 2) %>% pull(var)
  length(retained)
  # Number of samples to subsample
  subsamp = md_subbed %>% ungroup %>% count(!!sym(var)) %>% filter(n >= 2) %>% summarize(min = min(n)) %>% pull(min)
  
  cat('Downsampling to n = ',subsamp,', requiring min = ',2,' for variable: ',var,', ','replicate: ',rep,'\n')
  # Subsampling
  mdi = md_subbed %>%
    filter(!!sym(var) %in% retained) %>%
    group_by(!!sym(var)) %>%
    sample_n(min(n(), subsamp),replace = TRUE) %>%
    ungroup()
  mdi %>% count(!!sym(var))
  
  # Ensure the factor levels are dropped which aren't in the dataset after filtering
  mdi = droplevels(mdi)
  
  # Function to safely train models
  safe_train <- function(formula, data, method, trControl, metric) {
    tryCatch({
      train(formula, data = data, method = method, trControl = trControl, metric = metric, trace = FALSE)
    }, error = function(e) {
      message(paste("Model failed:", as.character(formula), "Error:", e$message))
      return(NULL)
    })
  }
  
  # Training models inside loop
  m1 = safe_train(as.formula(paste(var, "~ Haplogroup + Ancestry + Geography")), mdi, "multinom", ctrl, "AUC")
  m2 = safe_train(as.formula(paste(var, "~ Haplogroup + Geography")), mdi, "multinom", ctrl, "AUC")
  m3 = safe_train(as.formula(paste(var, "~ Haplogroup")), mdi, "multinom", ctrl, "AUC")
  m4 = safe_train(as.formula(paste(var, "~ Haplogroup + Ancestry")), mdi, "multinom", ctrl, "AUC")
  m5 = safe_train(as.formula(paste(var, "~ Ancestry")), mdi, "multinom", ctrl, "AUC")
  m6 = safe_train(as.formula(paste(var, "~ Ancestry + Geography")), mdi, "multinom", ctrl, "AUC")
  m7 = safe_train(as.formula(paste(var, "~ Geography")), mdi, "multinom", ctrl, "AUC")
  
  models = list(m1, m2, m3, m4, m5, m6, m7)
  model_names = c("m1", "m2", "m3", "m4", "m5", "m6", "m7")
  
  for (i in seq_along(models)) {
    mo = models[[i]]
    
    if (is.null(mo)) {
      message(paste("Skipping", model_names[i], "due to failure"))
      next
    }
    
    # Proceed with extracting metrics if model trained successfully
    final_model = mo$finalModel
    AIC = AIC(final_model)
    dat = data.frame(Model = model_names[i], Iteration = rep, Variable = var, Subsampled = subsamp, MinObs = 2,
                     decay = mo$results$decay, logLoss = mo$results$logLoss, logLossSD = mo$results$logLossSD, 
                     AIC = AIC, AUC = mo$results$AUC, Accuracy = mo$results$Accuracy, 
                     AccuracySD = mo$results$AccuracySD, Kappa = mo$results$Kappa)
    
    dat_best = dat %>% slice_min(logLoss)
    adat = rbind(adat, dat_best)
    
    # Confusion matrix
    pred = mo$pred
    pred_best = pred %>% filter(decay == dat_best$decay)
    predicted_classes = predict(mo, newdata = mdi, type = "raw")
    new = mdi %>% select(reference = !!sym(var)) %>% mutate(predicted = predicted_classes)
    conf_new = confusionMatrix(new$predicted, new$reference)
    
    conf_real = as.data.frame(conf_new$table) %>% mutate(
      Prediction = gsub('_','\\. ',Prediction),
      Reference = gsub('_','\\. ',Reference)) %>%
      mutate(Model = model_names[i], Iteration = counter, Variable = var, 
             Subsample = subsamp, MinObs = 2, AUC=dat_best$AUC,
             logloss=dat_best$logLoss, Accuracy = dat_best$Accuracy, 
             AccuracySD=dat_best$AccuracySD)
    
    new_preds = rbind(new_preds, conf_real)
    
  } # Exit model loop
} # Exit iteration loop

# Write the no-blue data
write.table(adat,paste0('20250320_Model_Selection_Boot-2Obs_',sp,'-NoM1M2M3.txt'),quote=F,sep='\t',row.names=F)
write.table(new_preds,paste0('20250320_ConfusionMatrix_Boot-2Obs_',sp,'-NoM1M2M3.txt'),quote=F,sep='\t',row.names=F)

```

### Plot Sensitivity

```R
# Egg associations, MNLR, Model Selection c optatus
#### Find associations between MT/W haplotypes and features 
setwd('~/EvoBioWolf/CUCKOO_gentes/gens_associations/mnlr/')
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

sp = 'CC'

set.seed(123)
md = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt') 

#only grab samples with known egg (hash out filter for related individuals)
md_egg = md %>%
  filter(SpeciesShort == sp & Analysis_GensAssociations == 1) %>%
  select(ID = ID, Egg, Haplogroup, Geography = GeographicGroup,Ancestry = AncestryA5)
md_egg %>% count(Egg)

md_egg = md_egg %>% group_by(Egg) %>% mutate(TotalEgg = n()) %>% ungroup
minobs=2
md_egg = md_egg %>% filter(TotalEgg >= minobs)
md_egg %>% count(Egg)

# Read in data from sensitivities
# for No M1/M2/M3
adat = read_tsv('20250320_Model_Selection_Boot-2Obs_CC-NoM1M2M3.txt')
conf_mats = read_tsv('20250320_ConfusionMatrix_Boot-2Obs_CC-NoM1M2M3.txt')

# for No E1/E6
adat = read_tsv('20250320_Model_Selection_Boot-2Obs_CC-NoE1E6.txt')
conf_mats = read_tsv('20250320_ConfusionMatrix_Boot-2Obs_CC-NoE1E6.txt')


# Plot results 
leg = data.frame(Model = c('m1','m2','m3','m4','m5','m6','m7'),
                 Label = c('A+G+M','G+M','M','A+M','A','A+G','G'))
model_dat = adat %>% left_join(.,leg)

# Plot AUC for each model across all the iterations
model_dat$Label = factor(model_dat$Label, levels=c('A+G+M','A+G','A+M','G+M','A','G','M'))
auc_plot = model_dat %>% ggplot(aes(y=Label,x=AUC,fill=Variable))+
  geom_boxplot()+
  scale_fill_manual(values=brewer.pal(3,'RdBu'))+
  xlab('AUC')+ylab('Model Covariates')+
  theme_bw(base_size=6)+theme(legend.position = 'top',strip.text.y.right = element_text(angle = 0))
auc_plot

# For no M1/M2/M3
ggsave('~/symlinks/host/figures/20250326_CC_MNLR_ModelSelection_AUC-NoM1M2M3.pdf',
       auc_plot,height=3,width=7,dpi=300)

# For no E1/E6
ggsave('~/symlinks/host/figures/20250326_CC_MNLR_ModelSelection_AUC-NoE1E6.pdf',
       auc_plot,height=3,width=7,dpi=300)


# Summarize AUC across the core 3 models 
auc_plot_input <- model_dat %>%
  # %>% 
  group_by(Label,Variable) %>% 
  sum_stats(AUC)

# Plot
cols <- brewer.pal(3,'Set2')[c(1,2,3)]
#model_dat$Label <- factor(model_dat$Label,levels=c('A','G','M'))

auc_summary_plot <- auc_plot_input %>% 
  filter(Label == 'A' | Label == 'G' | Label == 'M') %>% 
  ggplot(aes(y=Label,fill=Label,x=mean,xmin=conf_low,xmax=conf_high))+
  geom_bar(stat='identity',position=position_dodge(width=0.9))+
  scale_fill_manual(values=cols)+
  xlab('')+ylab('AUC (Mean & 95% CI)')+
  geom_errorbar(width=0.25,position=position_dodge(width=0.9))+
  theme_bw(base_size=6)+
  coord_cartesian(xlim=c(0.5,1))+
  theme(legend.position='top')
auc_summary_plot

# for no M1/M2/M3
ggsave('~/symlinks/host/figures/20250326_CC_MNLR_ModelSelection95CI_AUC-NoM1M2M3.pdf',
       auc_summary_plot,height=2,width=1.5,dpi=300)

# for no E1/E6
ggsave('~/symlinks/host/figures/20250326_CC_MNLR_ModelSelection95CI_AUC-NoE1E6.pdf',
       auc_summary_plot,height=2,width=1.5,dpi=300)

#order full, single plot, make sure the 3 variables are in order 
egglev <- md %>% select(Egg,EggOrder) %>% distinct %>% arrange(EggOrder)
conf_mats = conf_mats %>% mutate(Prediction = factor(Prediction,levels=egglev$Egg),
                                 Reference = factor(Reference,levels=egglev$Egg))

### Plot how the addition of haplogroup improves predictions show A+G (m6) vs A+G+M (m1)
auc_vals <- adat %>% group_by(Model) %>% sum_stats(AUC)


### ancestry + geography (A+G; m6)
lab <- auc_vals %>% filter(Model == 'm6') %>% mutate(label = paste0('A+G: ',round(conf_low,2),' - ',round(conf_high,2)))
repredictions_ancgeo = conf_mats %>% 
  filter(Model == 'm6') %>%  # plot model 6 (A+G) 
  group_by(Variable,MinObs,Iteration,Reference) %>% 
  mutate(TotalIteration = sum(Freq), #within each of the 100 replicates, calculate the frequency of correct classifications
         Proportion = Freq / TotalIteration) %>% ungroup %>% 
  group_by(Prediction,Reference,Variable) %>% 
  sum_stats(Proportion)
ancgeo_plot = repredictions_ancgeo %>% 
  ggplot(aes(x=Prediction,y=Reference,fill=mean,col=sd))+
  geom_tile(col=NA)+
  ggtitle(lab$label)+
  facet_wrap(lab$label~.,scales='free')+
  scale_fill_continuous(low='white',high='darkblue')+  theme_bw(base_size=6)+
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1),legend.position='top')
ancgeo_plot

# ancestry + geogrpahy + hap (A+G+M, m1)
lab <- auc_vals %>% filter(Model == 'm1') %>% mutate(label = paste0('A+G+M: ',round(conf_low,2),' - ',round(conf_high,2)))
repredictions_ancgeohap = conf_mats %>% 
  filter(Model == 'm1') %>%  # plot model 1 (A+G+M) 
  group_by(Variable,MinObs,Iteration,Reference) %>% 
  mutate(TotalIteration = sum(Freq), #within each of the 100 replicates, calculate the frequency of correct classifications
         Proportion = Freq / TotalIteration) %>% ungroup %>% 
  group_by(Prediction,Reference,Variable) %>% 
  sum_stats(Proportion)
ancgeohap_plot = repredictions_ancgeohap %>% 
  ggplot(aes(x=Prediction,y=Reference,fill=mean,col=sd))+
  geom_tile(col=NA)+
  ggtitle(lab$label)+
  facet_wrap(lab$label~.,scales='free')+
  scale_fill_continuous(low='white',high='darkblue')+  theme_bw(base_size=6)+
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1),legend.position='top')
ancgeohap_plot

ggarrange(ancgeo_plot,ancgeohap_plot,nrow = 2,common.legend = TRUE)
# For no M1/M2/M3
ggsave('~/symlinks/host/figures/20250320_MNLR_ConfusionMatrix-Repredictions-CC_M1vsM6-NoM1M2M3.pdf',
       ggarrange(ancgeo_plot,ancgeohap_plot,nrow = 2,common.legend = TRUE),
       dpi=300,height=3.5,width=1.5) # canorus

# For no E1/E6
ggsave('~/symlinks/host/figures/20250320_MNLR_ConfusionMatrix-Repredictions-CC_M1vsM6-NoE1E6.pdf',
       ggarrange(ancgeo_plot,ancgeohap_plot,nrow = 2,common.legend = TRUE),
       dpi=300,height=3.5,width=1.5) # canorus


```





## Discordance Analysis: Egg

Determine the number of shifts for each egg type, binarizing each egg and comparing mtDNA and autosomal trees.

Analyses occurring within 

```bash
/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/phylogenetics/202503_discordancen136
```

Do this separately for CC (n = 86) and CO (n=50):

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --mem=30000mb
#SBATCH --cpus-per-task=10
#SBATCH --time=24:00:00


#mamba activate snps 
SPECIES=$1

mkdir -p ml_trees

echo "WORKING ON ${SPECIES}"

auto=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/merged/snps_only/chr_1.SNPS.vcf.gz
mtdna=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/phylogenetics/202503_fulln136/vcfs/chr_MT.SNP.DP3-AC1-MQ40-MM1.vcf.gz

#### autosomes 
#Subset VCFS
bcftools view --threads 10 --samples-file ${SPECIES}.list -Ou ${auto} | \
        bcftools view --threads 10 --types snps --min-alleles 2 --max-alleles 2 --min-af 0.05 --max-af 0.95 -i "MQ>40 & F_MISSING < 0.1" -Ou | \
        bcftools +prune -m 0.1 --window 5kb -Oz -o autos_LD_${SPECIES}.vcf.gz
bcftools index --threads 10 autos_LD_${SPECIES}.vcf.gz

#Convert to phylip
python ~/modules/vcf2phylip.py -i autos_LD_${SPECIES}.vcf.gz -f --output-folder ml_trees

#create tree, autosomes
iqtree --redo -keep-ident -T 10 -s ml_trees/autos_LD_${SPECIES}.min4.phy --seqtype DNA -m "GTR+ASC" -B 1000
iqtree --redo -keep-ident -T 10 -s ml_trees/autos_LD_${SPECIES}.min4.phy.varsites.phy --seqtype DNA -m "GTR+ASC" -B 1000


#### mtDNA 
#Subset VCFS
bcftools view --force-samples --threads 10 --samples-file ${SPECIES}.list -Ou ${mtdna} | \
        bcftools view --threads 10 --types snps --min-alleles 2 --max-alleles 2 --min-af 0.05 --max-af 0.95 -i "MQ>40 & F_MISSING < 0.1" -Oz -o mtdna_${SPECIES}.vcf.gz
bcftools index --threads 10 mtdna_${SPECIES}.vcf.gz

#Convert to phylip
python ~/modules/vcf2phylip.py -i mtdna_${SPECIES}.vcf.gz -f --output-folder ml_trees

#create tree, autosomes
iqtree --redo -keep-ident -T 10 -s ml_trees/mtdna_${SPECIES}.min4.phy --seqtype DNA -m "GTR+ASC" -B 1000
iqtree --redo -keep-ident -T 10 -s ml_trees/mtdna_${SPECIES}.min4.phy.varsites.phy --seqtype DNA -m "GTR+ASC" -B 1000

```

And then compare egg shifts across both mtDNA and autosomal tree. This can take a while, so submit by egg:

```bash
#!/usr/bin/env rscript
args <- commandArgs(trailingOnly = TRUE)

#### Determine egg shift parsimony using binary classifications
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/phylogenetics/202503_discordancen136/ml_trees')
#.libPaths('~/mambaforge/envs/R/lib/R/library')
library(ape)
library(caper)
library(ggtree)
library(treeio)
library(phytools)
library(tidyverse)
library(viridis)

egg = args[1]

# Read in trees and metadata
md = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt') %>% filter(Analysis_FullGensAssociations == 1)
md$dummy <- 1 # Assign fake dummy var so that all connections are plotted

sp <- ifelse(grepl("^ECO", egg), "CO", "CC")

### mtDNA tree
m = read.iqtree(paste0('mtdna_',sp,'.min4.phy.contree'))
m1 <- midpoint.root(as.phylo(m))

### AUTOSOME tree
a = read.iqtree(paste0('autos_LD_',sp,'.min4.phy.varsites.phy.contree'))
a1 <- midpoint.root(as.phylo(a))

# Store results
spegg <- md %>% filter(SpeciesShort == sp)

results <- list()
for (tree in c('m1','a1')) {

    cat('Working on egg type: ',egg,' for tree: ',tree,'\n')

    # Change egg to binary trait, only target egg is 1 all else is 0
    md_mod <- spegg %>% mutate(Egg = ifelse(Egg == egg,egg,'E0'))

    # Generate base tree
    targ_tree <- get(tree)
    ggt <- ggtree(targ_tree, layout = "circular",branch.length='none') %<+% md_mod

    #grab only egg
    phenos <- as.data.frame(ggt$data %>% filter(isTip == TRUE))
    egg_mat <- as.matrix(phenos %>% select(Egg))
    phenotypes <- setNames(as.vector(egg_mat[, 1]), targ_tree$tip.label)

    #inspect tree
    #ggtree(targ_tree, layout = "rectangular",branch.length='none') %<+% md_mod +
    #  geom_tippoint(aes(fill=Egg),pch=21,size=2)+
    #  scale_fill_brewer(palette='Set2')

    #Plot probabilities
    t2 <- multi2di(targ_tree)

    # Quick ape method
    fitER <- ape::ace(phenotypes,t2,model="ER",type="discrete")

    # Determine ancestral state likelihood number of shifts using MCMC
    sim <- make.simmap(t2, phenotypes, model="ER", nsim=2,Q='mcmc')

    # Output stats on the number of shifts across the 100 replicates
    vec <- data.frame(shifts = countSimmap(sim)$Tr[,1])
    shifts <- vec %>% mutate(Egg = egg, Tree = tree)

    results[[paste0(egg,'_',tree)]] <- shifts

}

full_results <- rbindlist(results)

write.table(full_results,file=paste0('20250331_Results_Binary_EggComparison__',egg,'.txt'),quote=F,sep='\t',row.names=F)

```

Submit:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --time=48:00:00

#mamba activate r
egg=$1

Rscript Transitions.R ${egg}

```

### Plot 



```R
#### Determine egg shift parsimony using binary classifications 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/phylogenetics/202503_discordancen136')
.libPaths('~/mambaforge/envs/R/lib/R/library')
library(ape)
library(caper)
library(ggtree)
library(treeio)
library(phytools)
library(tidyverse)
library(viridis)
library(pheatmap)
library(ggpubr)
library(meRo)
library(data.table)
library(tangler)
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(magick)

# Read in trees and metadata 
md <- read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt') %>% filter(Analysis_PopulationGenetics == 1) %>% drop_na(Egg)
md$dummy <- 1 # Assign fake dummy var so that all connections are plotted 

# Read in 
full_results <- read_tsv('20250331_TransitionsEgg.txt')

# Add egg colors
egglevs = md %>% select(Egg,EggCol,EggShape,EggOrder) %>% unique %>%arrange(EggOrder)
full_results$Egg = factor(full_results$Egg,levels=egglevs$Egg)

# Each egg / tree has 100 bootstrapps, so add an identifier for each
fr <- full_results %>% 
  group_by(Egg, Tree) %>% 
  mutate(rep = row_number()) %>% ungroup

# And for each bootstrap compare mtDNA vs autosomal shifts 
mt_vs_aut <- fr %>% 
  pivot_wider(names_from = Tree,values_from = shifts) %>% 
  mutate(shifts = m1 - a1)

# Calculate summary stats incl. 95% CIs
confs <- mt_vs_aut %>% group_by(Egg) %>% sum_stats(shifts)

plot_mtaut <- mt_vs_aut %>% ungroup %>% 
  ggplot(aes(x=Egg,y=shifts,fill=Egg))+
  geom_hline(yintercept=0,lty=2)+
  geom_hline(yintercept=c(-1,1),lty=3,col='darkgray')+
  geom_errorbar(data=confs,width=0.2,col='black',aes(x=Egg,ymin=conf_low,ymax=conf_high),inherit.aes=FALSE,position=position_nudge(x=-0.3))+
  geom_boxplot(width = .35,outlier.shape = NA, position=position_nudge(x=0.15)) +
  scale_fill_manual(values=egglevs$EggCol,breaks=egglevs$Egg)+
  scale_color_manual(values=egglevs$EggCol,breaks=egglevs$Egg)+
  theme_bw()+
  ylab('')+xlab('')+
  coord_flip()
plot_mtaut

pdf('~/symlinks/host/figures/20250401_mtDNA-vs-Autosomal-shifts-MCMC.pdf',height=4.5,width=4)
plot_mtaut
dev.off()

#### Plot tree discordance ####
# This will also estimate transitions from each egg type to other eggs,  but doesn't seem as reliable as binary above approach

full_tree_results <- list()
full_egg_results <- list()

for (sp in c('CC','CO')) {
  
  ### mtDNA tree
  m = read.iqtree(paste0('ml_trees/mtdna_',sp,'.min4.phy.contree'))
  m1 <- midpoint.root(as.phylo(m))
  
  ### AUTOSOME tree
  a = read.iqtree(paste0('ml_trees/autos_LD_',sp,'.min4.phy.varsites.phy.contree'))
  a1 <- midpoint.root(as.phylo(a))
  
  # Store results
  spegg <- md %>% filter(SpeciesShort == sp)
  
  for (tree in c('m1','a1')) {
    
    cat('Running full ML reconstruction on tree: ',tree,' for species ',sp,'\n')
    
    # Generate base tree 
    targ_tree <- get(tree)
    ggt <- ggtree(targ_tree, layout = "rectangular",branch.length='none') %<+% spegg
    
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
    shifts <- vec %>% mutate(Tree = tree, Species = sp, MLloglik = fitER$loglik, MLest = fitER$rates, MLse = fitER$se) 
    
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
    E_groups <- spegg %>% select(Egg) %>% unique %>% pull(Egg)
    averages <- sapply(E_groups, function(E) calculate_statistics(mat, E))
    
    egg_shifts <- as.data.frame(t(averages)) %>% mutate(Egg = rownames(.), Tree = tree, Species = sp)
    rownames(egg_shifts) <- NULL
    
    full_tree_results[[paste0(sp,'_',tree)]] <- shifts
    full_egg_results[[paste0(sp,'_',tree)]] <- egg_shifts
    
    # Extract nodes and the proportions for pies
    nodes <- data.frame(
      node=1:t2$Nnode+Ntip(t2),
      fitER$lik.anc)
    
    # For stochastic mapping 
    obj <- describe.simmap(simfull,plot=FALSE)
    mcmc_nodes <- as.data.frame(cbind(node=rownames(obj$ace),obj$ace)); rownames(mcmc_nodes) <- NULL
    mcmc_nodes <- mcmc_nodes %>% mutate(across(starts_with('E'), as.numeric))
    nodes_plot <- mcmc_nodes %>% filter(node %in% nodes$node)
    rownames(nodes_plot) <- nodes_plot$node
    nodes_plot$node <- as.integer(nodes_plot$node)
    
    ## cols parameter indicate which columns store stats
    pies <- nodepie(nodes_plot, cols=2:ncol(nodes_plot),outline.color='black',outline.size = 0.1)
    pies <- lapply(pies, function(g) g+scale_fill_manual(values = spegg$EggCol,breaks=spegg$Egg))
    
    t3 <- full_join(t2, data.frame(label = names(phenotypes), stat = phenotypes ), by = 'label')
    tp <- ggtree(t3,layout='rectangular',branch.length = 'none') %<+% md
    tp$data$dummy <- 1
    tp_final <- tp + geom_inset(pies, width = .09, height = .09)
    tp_phenos <- tp_final +
      geom_tippoint(aes(fill=Haplogroup),pch=21,size=1.5)+
      scale_fill_manual(values=md$HaplogroupColor,breaks=md$Haplogroup)
    assign(paste0(sp,'_',tree,'_nodes'),tp_final)
    assign(paste0(sp,'_',tree,'_pies'),tp_phenos)
  }
}

# Grab the model results
full_search_res <- rbindlist(full_tree_results)
full_search_res %>% group_by(Tree,Species,MLloglik,MLest,MLse) %>% sum_stats(shifts)

# Also bind the egg data 
full_search_eggs <- rbindlist(full_egg_results)
full_search_eggs
full_search_eggs$Mean <- unlist(full_search_eggs$Mean)
full_search_eggs$Median <- unlist(full_search_eggs$Median)
full_search_eggs$SD <- unlist(full_search_eggs$SD)

write.table(full_search_res,file='20250401_Full_Search_Results.txt',quote=F,sep='\t',row.names=F)
write.table(full_search_eggs,file='20250401_Full_Search_EggResults.txt',quote=F,sep='\t',row.names=F)

#full_search_res <- read_tsv('20250401_Full_Search_EggResults.txt')

# Plot all connections CC
mt_egg <- CC_m1_nodes +
  geom_tippoint(aes(fill=Egg),pch=22,size=1.5)+
  scale_fill_manual(values=egglevs$EggCol,breaks=egglevs$Egg)
a_egg <- CC_a1_nodes +
  geom_tippoint(aes(fill=Egg),pch=22,size=1.5)+
  scale_fill_manual(values=egglevs$EggCol,breaks=egglevs$Egg)
discord_egg <- simple.tanglegram(tree1=mt_egg,tree2=a_egg,column=dummy,value=1,t2_pad=2,l_color='black',tiplab=F)
discord_egg2 <- simple.tanglegram(tree1=a_egg,tree2=mt_egg,column=dummy,value=1,t2_pad=2,l_color='black',tiplab=F)

pdf('~/symlinks/host/figures/20250401_TreeCompare-mtDNA-Auto-NodePiesML-CC-EGG.pdf',height=5,width=7)
discord_egg
dev.off()

pdf('~/symlinks/host/figures/20250401_TreeCompare-Auto-mtDNA-NodePiesML-CC-EGG.pdf',height=5,width=7)
discord_egg2
dev.off()

### CO 
mt_egg <- CO_m1_nodes +
  geom_tippoint(aes(fill=Egg),pch=22,size=1.5)+
  scale_fill_manual(values=egglevs$EggCol,breaks=egglevs$Egg)
a_egg <- CO_a1_nodes +
  geom_tippoint(aes(fill=Egg),pch=22,size=1.5)+
  scale_fill_manual(values=egglevs$EggCol,breaks=egglevs$Egg)
discord_egg <- simple.tanglegram(tree1=mt_egg,tree2=a_egg,column=dummy,value=1,t2_pad=2,l_color='black',tiplab=F)
discord_egg2 <- simple.tanglegram(tree1=a_egg,tree2=mt_egg,column=dummy,value=1,t2_pad=2,l_color='black',tiplab=F)

pdf('~/symlinks/host/figures/20250401_TreeCompare-mtDNA-Auto-NodePiesML-CO-EGG.pdf',height=5,width=7)
discord_egg
dev.off()

pdf('~/symlinks/host/figures/20250401_TreeCompare-Auto-mtDNA-NodePiesML-CO-EGG.pdf',height=5,width=7)
discord_egg2
dev.off()
```

### Sensitivity: Females & W

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --mem=30000mb
#SBATCH --cpus-per-task=10
#SBATCH --time=24:00:00


#mamba activate snps 
SPECIES=$1

mkdir -p ml_trees

echo "WORKING ON ${SPECIES}"

auto=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/merged/snps_only/chr_1.SNPS.vcf.gz
mtdna=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/phylogenetics/202503_fulln136/vcfs/chr_MT.SNP.DP3-AC1-MQ40-MM1.vcf.gz
w=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/phylogenetics/202503_fulln136/vcfs/chr_W.SNP.DP3-AC1-MQ40-MM1.vcf.gz

#### autosomes 
#Subset VCFS
bcftools view --threads 10 --samples-file ${SPECIES}F.list -Ou autos_LD_${SPECIES}.vcf.gz | \
        bcftools view --threads 10 --types snps --min-alleles 2 --max-alleles 2 --min-af 0.05 --max-af 0.95 -i "MQ>40 & F_MISSING < 0.1" -Oz -o autos_LD_${SPECIES}F.vcf.gz
bcftools index --threads 10 autos_LD_${SPECIES}F.vcf.gz

#Convert to phylip
python ~/modules/vcf2phylip.py -i autos_LD_${SPECIES}F.vcf.gz -f --output-folder ml_trees

#create tree, autosomes
iqtree --redo -keep-ident -T 10 -s ml_trees/autos_LD_${SPECIES}F.min4.phy --seqtype DNA -m "GTR+ASC" -B 1000
iqtree --redo -keep-ident -T 10 -s ml_trees/autos_LD_${SPECIES}F.min4.phy.varsites.phy --seqtype DNA -m "GTR+ASC" -B 1000


#### mtDNA 
#Subset VCFS
bcftools view --force-samples --threads 10 --samples-file ${SPECIES}F.list -Ou ${mtdna} | \
        bcftools view --threads 10 --types snps --min-alleles 2 --max-alleles 2 --min-af 0.05 --max-af 0.95 -i "MQ>40 & F_MISSING < 0.1" -Oz -o mtdna_${SPECIES}F.vcf.gz
bcftools index --threads 10 mtdna_${SPECIES}F.vcf.gz

#Convert to phylip
python ~/modules/vcf2phylip.py -i mtdna_${SPECIES}F.vcf.gz -f --output-folder ml_trees

#create tree, autosomes
iqtree --redo -keep-ident -T 10 -s ml_trees/mtdna_${SPECIES}F.min4.phy --seqtype DNA -m "GTR+ASC" -B 1000
iqtree --redo -keep-ident -T 10 -s ml_trees/mtdna_${SPECIES}F.min4.phy.varsites.phy --seqtype DNA -m "GTR+ASC" -B 1000


#### W 
#Subset VCFS
bcftools view --force-samples --threads 10 --samples-file ${SPECIES}F.list -Ou ${w} | \
        bcftools view --threads 10 --types snps --min-alleles 2 --max-alleles 2 --min-af 0.05 --max-af 0.95 -i "MQ>40 & F_MISSING < 0.1" -Oz -o w_${SPECIES}F.vcf.gz
bcftools index --threads 10 w_${SPECIES}F.vcf.gz

#Convert to phylip
python ~/modules/vcf2phylip.py -i w_${SPECIES}F.vcf.gz -f --output-folder ml_trees

#create tree, autosomes
iqtree --redo -keep-ident -T 10 -s ml_trees/w_${SPECIES}F.min4.phy --seqtype DNA -m "GTR+ASC" -B 1000
iqtree --redo -keep-ident -T 10 -s ml_trees/w_${SPECIES}F.min4.phy.varsites.phy --seqtype DNA -m "GTR+ASC" -B 1000

```

Discordance:

```R
#!/usr/bin/env rscript
args <- commandArgs(trailingOnly = TRUE)

#### Determine egg shift parsimony using binary classifications
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/phylogenetics/202503_discordancen136/ml_trees')
#.libPaths('~/mambaforge/envs/R/lib/R/library')
library(ape)
library(caper)
library(ggtree)
library(treeio)
library(phytools)
library(tidyverse)
library(viridis)
library(data.table)

egg = args[1]

# Read in trees and metadata
md = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt') %>% filter(Analysis_FullGensAssociations == 1)
md$dummy <- 1 # Assign fake dummy var so that all connections are plotted

sp <- ifelse(grepl("^ECO", egg), "CO", "CC")

### mtDNA tree
m = read.iqtree(paste0('mtdna_',sp,'F.min4.phy.contree'))
m1 <- midpoint.root(as.phylo(m))

### W tree
w = read.iqtree(paste0('w_',sp,'F.min4.phy.contree'))
w1 <- midpoint.root(as.phylo(w))

### AUTOSOME tree
a = read.iqtree(paste0('autos_LD_',sp,'F.min4.phy.varsites.phy.contree'))
a1 <- midpoint.root(as.phylo(a))

# Store results
spegg <- md %>% filter(SpeciesShort == sp)

results <- list()
for (tree in c('m1','a1','w1')) {
  
  cat('Working on egg type: ',egg,' for tree: ',tree,'\n')
  
  # Change egg to binary trait, only target egg is 1 all else is 0
  md_mod <- spegg %>% mutate(Egg = ifelse(Egg == egg,egg,'E0'))
  
  # Generate base tree
  targ_tree <- get(tree)
  ggt <- ggtree(targ_tree, layout = "circular",branch.length='none') %<+% md_mod
  
  #grab only egg
  phenos <- as.data.frame(ggt$data %>% filter(isTip == TRUE))
  egg_mat <- as.matrix(phenos %>% dplyr::select(Egg))
  phenotypes <- setNames(as.vector(egg_mat[, 1]), targ_tree$tip.label)
  
  #inspect tree
  #ggtree(targ_tree, layout = "rectangular",branch.length='none') %<+% md_mod +
  #  geom_tippoint(aes(fill=Egg),pch=21,size=2)+
  #  scale_fill_brewer(palette='Set2')
  
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

full_results <- rbindlist(results)

write.table(full_results,file=paste0('20250331_Results_Binary_EggComparison__',egg,'_F.txt'),quote=F,sep='\t',row.names=F)

```

Plot:

```bash
#### Determine egg shift parsimony using binary classifications: FEMALES only, W tree  
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/phylogenetics/202503_discordancen136')
.libPaths('~/mambaforge/envs/R/lib/R/library')
library(ape)
library(caper)
library(ggtree)
library(treeio)
library(phytools)
library(tidyverse)
library(viridis)
library(pheatmap)
library(ggpubr)
library(meRo)
library(data.table)
library(tangler)
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(magick)

# Read in trees and metadata 
md <- read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt') %>% filter(Analysis_PopulationGenetics == 1) %>% drop_na(Egg)
md$dummy <- 1 # Assign fake dummy var so that all connections are plotted 

# Read in 
full_results <- read_tsv('20250331_TransitionsEggF.txt')

# Add egg colors
egglevs = md %>% select(Egg,EggCol,EggShape,EggOrder) %>% unique %>%arrange(EggOrder)
full_results$Egg = factor(full_results$Egg,levels=egglevs$Egg)

# Each egg / tree has 100 bootstrapps, so add an identifier for each
fr <- full_results %>% 
  group_by(Egg, Tree) %>% 
  mutate(rep = row_number()) %>% ungroup

# And for each bootstrap compare mtDNA vs autosomal shifts 
mt_vs_aut <- fr %>% 
  pivot_wider(names_from = Tree,values_from = shifts) %>% 
  mutate(shifts = m1 - a1,
         shiftsw = w1 - a1)

# Calculate summary stats incl. 95% CIs
confs <- mt_vs_aut %>% group_by(Egg) %>% sum_stats(shifts)

plot_mtaut <- mt_vs_aut %>% ungroup %>%
  pivot_longer(c(shifts,shiftsw)) %>%
  ggplot(aes(x=Egg,y=value,fill=Egg))+
  geom_hline(yintercept=0,lty=2)+
  geom_hline(yintercept=c(-1,1),lty=3,col='darkgray')+
  geom_errorbar(data=confs,width=0.2,col='black',aes(x=Egg,ymin=conf_low,ymax=conf_high),inherit.aes=FALSE,position=position_nudge(x=-0.3))+
  geom_boxplot(width = .35,outlier.shape = NA, position=position_nudge(x=0.15)) +
  scale_fill_manual(values=egglevs$EggCol,breaks=egglevs$Egg)+
  scale_color_manual(values=egglevs$EggCol,breaks=egglevs$Egg)+
  facet_grid(.~name,scales='free')+
  theme_bw()+
  ylab('')+xlab('')+
  coord_flip()
plot_mtaut

pdf('~/symlinks/host/figures/20250401_W-vs-mtDNA-vs-Autosomal-shifts-MCMC-FEMALES.pdf',height=4.5,width=6.5)
plot_mtaut
dev.off()

#### Plot tree discordance ####
# This will also estimate transitions from each egg type to other eggs,  but doesn't seem as reliable as binary above approach

full_tree_results <- list()
full_egg_results <- list()

for (sp in c('CC','CO')) {
  
  ### mtDNA tree
  m = read.iqtree(paste0('ml_trees/mtdna_',sp,'F.min4.phy.contree'))
  m1 <- midpoint.root(as.phylo(m))
  
  ### w tree
  w = read.iqtree(paste0('ml_trees/w_',sp,'F.min4.phy.contree'))
  w1 <- midpoint.root(as.phylo(w))
  
  ### AUTOSOME tree
  a = read.iqtree(paste0('ml_trees/autos_LD_',sp,'F.min4.phy.varsites.phy.contree'))
  a1 <- midpoint.root(as.phylo(a))
  
  # Store results
  spegg <- md %>% filter(SpeciesShort == sp)
  
  for (tree in c('m1','a1','w1')) {
    
    cat('Running full ML reconstruction on tree: ',tree,' for species ',sp,'\n')
    
    # Generate base tree 
    targ_tree <- get(tree)
    ggt <- ggtree(targ_tree, layout = "rectangular",branch.length='none') %<+% spegg
    
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
    shifts <- vec %>% mutate(Tree = tree, Species = sp, MLloglik = fitER$loglik, MLest = fitER$rates, MLse = fitER$se) 
    
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
    E_groups <- spegg %>% select(Egg) %>% unique %>% pull(Egg)
    averages <- sapply(E_groups, function(E) calculate_statistics(mat, E))
    
    egg_shifts <- as.data.frame(t(averages)) %>% mutate(Egg = rownames(.), Tree = tree, Species = sp)
    rownames(egg_shifts) <- NULL
    
    full_tree_results[[paste0(sp,'_',tree)]] <- shifts
    full_egg_results[[paste0(sp,'_',tree)]] <- egg_shifts
    
    # Extract nodes and the proportions for pies
    nodes <- data.frame(
      node=1:t2$Nnode+Ntip(t2),
      fitER$lik.anc)
    
    # For stochastic mapping 
    obj <- describe.simmap(simfull,plot=FALSE)
    mcmc_nodes <- as.data.frame(cbind(node=rownames(obj$ace),obj$ace)); rownames(mcmc_nodes) <- NULL
    mcmc_nodes <- mcmc_nodes %>% mutate(across(starts_with('E'), as.numeric))
    nodes_plot <- mcmc_nodes %>% filter(node %in% nodes$node)
    rownames(nodes_plot) <- nodes_plot$node
    nodes_plot$node <- as.integer(nodes_plot$node)
    
    ## cols parameter indicate which columns store stats
    pies <- nodepie(nodes_plot, cols=2:ncol(nodes_plot),outline.color='black',outline.size = 0.1)
    pies <- lapply(pies, function(g) g+scale_fill_manual(values = spegg$EggCol,breaks=spegg$Egg))
    
    t3 <- full_join(t2, data.frame(label = names(phenotypes), stat = phenotypes ), by = 'label')
    tp <- ggtree(t3,layout='rectangular',branch.length = 'none') %<+% md
    tp$data$dummy <- 1
    size = ifelse(sp == 'CO',1,0.09)
    tp_final <- tp + geom_inset(pies, width = size, height = size)
    tp_phenos <- tp_final +
      geom_tippoint(aes(fill=Haplogroup),pch=21,size=1.5)+
      scale_fill_manual(values=md$HaplogroupColor,breaks=md$Haplogroup)
    assign(paste0(sp,'_',tree,'_nodes'),tp_final)
    assign(paste0(sp,'_',tree,'_pies'),tp_phenos)
  }
}

# Grab the model results
full_search_res <- rbindlist(full_tree_results)
full_search_res %>% group_by(Tree,Species,MLloglik,MLest,MLse) %>% sum_stats(shifts)

# Also bind the egg data 
full_search_eggs <- rbindlist(full_egg_results)
full_search_eggs
full_search_eggs$Mean <- unlist(full_search_eggs$Mean)
full_search_eggs$Median <- unlist(full_search_eggs$Median)
full_search_eggs$SD <- unlist(full_search_eggs$SD)

write.table(full_search_res,file='20250401_Full_Search_Results.txt',quote=F,sep='\t',row.names=F)
write.table(full_search_eggs,file='20250401_Full_Search_EggResults.txt',quote=F,sep='\t',row.names=F)

#full_search_res <- read_tsv('20250401_Full_Search_EggResults.txt')

# Plot all connections CC
mt_egg <- CC_m1_nodes +
  geom_tippoint(aes(fill=Egg),pch=22,size=1.5)+
  scale_fill_manual(values=egglevs$EggCol,breaks=egglevs$Egg)
a_egg <- CC_a1_nodes +
  geom_tippoint(aes(fill=Egg),pch=22,size=1.5)+
  scale_fill_manual(values=egglevs$EggCol,breaks=egglevs$Egg)
discord_egg <- simple.tanglegram(tree1=mt_egg,tree2=a_egg,column=dummy,value=1,t2_pad=2,l_color='black',tiplab=F)
discord_egg2 <- simple.tanglegram(tree1=a_egg,tree2=mt_egg,column=dummy,value=1,t2_pad=2,l_color='black',tiplab=F)

pdf('~/symlinks/host/figures/20250401_TreeCompare-mtDNA-Auto-NodePiesML-CC-EGG.pdf',height=5,width=7)
discord_egg
dev.off()

pdf('~/symlinks/host/figures/20250401_TreeCompare-Auto-mtDNA-NodePiesML-CC-EGG.pdf',height=5,width=7)
discord_egg2
dev.off()

### CO 
mt_egg <- CO_m1_nodes +
  geom_tippoint(aes(fill=Egg),pch=22,size=0.5)+
  scale_fill_manual(values=egglevs$EggCol,breaks=egglevs$Egg)
a_egg <- CO_a1_nodes +
  geom_tippoint(aes(fill=Egg),pch=22,size=0.5)+
  scale_fill_manual(values=egglevs$EggCol,breaks=egglevs$Egg)
discord_egg <- simple.tanglegram(tree1=mt_egg,tree2=a_egg,column=dummy,value=1,t2_pad=2,l_color='black',tiplab=F)
discord_egg2 <- simple.tanglegram(tree1=a_egg,tree2=mt_egg,column=dummy,value=1,t2_pad=2,l_color='black',tiplab=F)

pdf('~/symlinks/host/figures/20250401_TreeCompare-mtDNA-Auto-NodePiesML-CO-EGG.pdf',height=5,width=7)
discord_egg
dev.off()

pdf('~/symlinks/host/figures/20250401_TreeCompare-Auto-mtDNA-NodePiesML-CO-EGG.pdf',height=3,width=5)
discord_egg2
dev.off()
```





# Analyses: Egg Hunt

Two primary strategies here:

* Estimate background FST in 10KB windows for each egg type (FSTobs - max(FST[randomized pops n=10])) using a target egg vs all other eggs approach. 
* Estimate base-pair level FST in C. canorus to get blue egg and reverted egg mutations. 

## Preparation

### Sample Selection

```bash
# Egg gene scan, sample selection 
#### Find associations between MT/W haplotypes and features 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/scan_n4')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(viridis)
library(tidyverse)
library(sf)
library(ggspatial)

md = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt') 
md = md %>% drop_na(Egg) %>% filter(Sex == 'F' & Analysis_PopulationGenetics == 1 & Egg != 'NA')
md %>% count(Egg)
keep = md %>% count(Egg) %>% filter(n >= 4)
mds = md %>% filter(Egg %in% keep$Egg)
mds %>% count(Egg)

write.table(mds %>% select(ID,Egg) %>% arrange(Egg),file='Eggs.pop',quote=F,sep='\t',row.names=F,col.names=F)
write.table(mds %>% select(ID,Egg) %>% arrange(Egg) %>% filter(grepl('ECC',Egg)),file='CC.pop',quote=F,sep='\t',row.names=F,col.names=F)
write.table(mds %>% select(ID,Egg) %>% arrange(Egg) %>% filter(grepl('ECO',Egg)),file='CO.pop',quote=F,sep='\t',row.names=F,col.names=F)

# Grab all the pairwise comparisons 
groups = mds %>% select(Egg) %>% unique %>% pull(Egg)
pairwise_combinations <- combn(groups, 2)

# Convert the combinations into a dataframe
pairwise_combinations_df <- data.frame(
  p1 = pairwise_combinations[1,],
  p2 = pairwise_combinations[2,]
)

pairwise_combinations_df <- pairwise_combinations_df %>% filter( (grepl('ECC',p1) & grepl('ECC',p2) | grepl('ECO',p1) & grepl('ECO',p2) ) )
write.table(pairwise_combinations_df,file='Pairwise_Contrasts-WithinSp.list',quote=F,sep='\t',row.names=F,col.names=F)


# Plot the individuals
#jitter points up to 1 lat/long for viewing
mds = mds %>% mutate(LatJit = jitter(Latitude,amount =2),
                     LonJit = jitter(Longitude,amount=2))
sites = st_as_sf(mds, coords = c("LonJit", "LatJit"), 
                 crs = 4326, agr = "constant")

#set up map and convert df to coordinate frame
world = map_data("world")


imm_spat = ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), col='grey90', fill='white') +
  geom_sf(data = sites, 
          aes(fill=Egg,shape=SpeciesShort),
          size=3,show.legend = T) +
  xlab('')+ylab('')+
  coord_sf(xlim = c(min(mds$Longitude)-5, max(mds$Longitude)+5), 
           ylim = c(min(mds$Latitude)-5, max(mds$Latitude)+5), expand = FALSE)+
  scale_fill_manual(values=mds$EggCol,breaks=mds$Egg)+
  scale_shape_manual(values=mds$Shape,breaks=mds$SpeciesShort)+
  theme_classic()+
  facet_grid(SpeciesShort ~.)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"))+
  theme(legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))+
  annotation_scale(line_width=0.5)+
  guides(fill=guide_legend(nrow=6,override.aes=list(shape=21)))
imm_spat

ggsave('~/symlinks/host/figures/20250404_EggHunt_Spatial.pdf',imm_spat,
       dpi=300,height=6,width=8)


```



```bash
  Egg       n
  <chr> <int>
1 ECC1     11
2 ECC10     4
3 ECC6     33
4 ECO1      8
5 ECO3      5
6 ECO4     11
```

### Filter VCF

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --time=48:00:00

# mamba activate snps
# for CHR in $(cat Chromosomes.list); do sbatch -J Filter_${CHR} 2.Refilter_VCF.sh ${CHR}; done
CHR=$1

# Modify WD depending on canorus or optatus
WD=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/gwas_n3/both/females_only
raw_vcfs=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/merged/snps_only/
mask=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/Coverage_Masks/GCA_017976375.1_bCucCan1.pri_genomic.CHR.N75-DoubleCoverage.mask.bed

cd $WD

mkdir -p vcfs

#genotypes BELOW this will be set to missing
MINDP=3

if [[ $CHR = 'chr_W' || $CHR = 'chr_MT' ]]; then
        PLOIDY=1

        echo "FILTERING AND MERGING VARIANT SITES FOR ${CHR}, PLOIDY: ${PLOIDY}"
        bcftools view --threads 5 --samples-file Females.list --force-samples -Ou ${raw_vcfs}/${CHR}.SNPS.vcf.gz | \
                bcftools view --types snps --min-alleles 2 --max-alleles 2 --threads 5 | \
                #remove SNPs in bad coverage regions
                bedtools subtract -header -a - -b ${mask} | \
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
        bcftools view --threads 5 --samples-file Females.list -Ou ${raw_vcfs}/${CHR}.SNPS.vcf.gz | \
                bcftools view --types snps --min-alleles 2 --max-alleles 2 --threads 5 | \
                #set genotypes below MINDP to missing
                bcftools +setGT -- -t q -i "FMT/DP < ${MINDP}" -n "./." | \
                #update AC fields
                bcftools +fill-tags -Ou -- -t AC,AN | \
                bcftools view --min-alleles 2 --max-alleles 2 --types snps --min-ac 1 -i 'MQ > 40 & F_MISSING < 0.1' -Oz -o vcfs/${CHR}.SNP.DP3.vcf.gz
        bcftools index --threads 5 vcfs/${CHR}.SNP.DP3.vcf.gz

fi
```



### Annotate SNPS

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

write.table(dat,file=paste0('raw_dnds/tmp.',chr),quote=F,sep='\t',row.names=F)
dat = read.table(paste0('raw_dnds/tmp.',chr),header=TRUE)
writedat = dat %>%
  as.data.frame %>%
  unique %>%
  group_by(seqnames, start, REF, varAllele, GENEID) %>%
  summarise(CONSEQUENCE = paste(unique(CONSEQUENCE), collapse = ";"), .groups = "drop") %>% ungroup %>%
  arrange(seqnames,start)

write.table(writedat,file=paste0('raw_dnds/Annotated_Variants_',chr,'__20250327.txt'),quote=F,sep='\t',row.names=F)

```

Submit:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_highmem
#SBATCH --cpus-per-task=4
#SBATCH --time=8:00:00

CHR=$1

# Submit: for CHR in $(cat Chromosomes.list); do sbatch -J dnds_${CHR} 3B.Submit_Annotation.sh ${CHR}; done
# mamba activate R
Rscript 3.Annotate_Variants.R ${CHR}

sed '1d' raw_dnds/Annotated_Variants_${CHR}__20250327.txt | awk '{OFS="\t"}{print $1, $2, $2, $5, $6}' | sed 's/gene-//g' | sed 's/ID=//g' > raw_dnds/${CHR}.bed
rm raw_dnds/tmp.${CHR}

```



### Merge autosomal VCFS

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=5
#SBATCH --time=48:00:00

# mamba activate snps
GROUP=$1

mkdir -p sweeps

bcftools view --threads 10 --samples-file pops/${GROUP}.list autosomal_files/autos.vcf.gz | \
        bcftools view --threads 10 --min-alleles 2 --max-alleles 2 --types snps -e 'F_MISSING > 0.1' --min-af 0.05 --max-af 0.95 -Ov -o autosomal_files/AllChromosomes2N.${GROUP}.vcf
bcftools index --threads 10 autosomal_files/AllChromosomes2N.${GROUP}.vcf

cd sweeps

#run raised, -M = missingness, impute per SNP
~/modules/RAiSD/raisd-master/RAiSD -n ${GROUP} -I ../autosomal_files/AllChromosomes2N.${GROUP}.vcf -A 0.995 -M 1 -y 2 -P -f

mkdir -p pop_vcfs 
# And also subset each for xp-EHH vcfs 
for CHR in $(cat Chromosomes.list); do 
        bcftools view --threads 5 --samples-file pops/${GROUP}.list vcfs/${CHR}.SNP.DP3.vcf.gz | \
                bcftools view --threads 5 --min-alleles 2 --max-alleles 2 --types snps -e 'F_MISSING > 0.1' --min-af 0.05 --max-af 0.95 -Oz -o pop_vcfs/${CHR}.${GROUP}.vcf.gz
        bcftools index --threads 5 pop_vcfs/${CHR}.${GROUP}.vcf.gz
        
        java -Xmx16g -jar ~/modules/beagle.28Jun21.220.jar gt=pop_vcfs/${CHR}.${GROUP}.vcf.gz out=pop_vcfs/${CHR}.${GROUP}.phased nthreads=5 impute=false
    bcftools index --threads 5 pop_vcfs/${CHR}.${GROUP}.phased.vcf.gz
    
done 
```

Also run on Z separately:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=5
#SBATCH --time=48:00:00

# mamba activate snps
GROUP=$1

mkdir -p sweeps

bcftools view --threads 10 --samples-file pops/${GROUP}.list vcfs/chr_Z.SNP.DP3.vcf.gz | \
        bcftools view --threads 10 --min-alleles 2 --max-alleles 2 --types snps -e 'F_MISSING > 0.1' --min-af 0.05 --max-af 0.95 -Ov -o autosomal_files/Zchr.${GROUP}.vcf
bcftools index --threads 10 autosomal_files/Zchr.${GROUP}.vcf

```

## bFST: 10-KB Windows [Fig 3]

Estimate background FST in 10-KB windows. 

For each egg type, e.g. ECC1: 

* Identify number of samples for ECC1 and background
* Run 10 replicates: swap the population labels so that there are N=ECC1 samples randomly assigned ECC1
* Estimate FST in 10KB windows for the real comparison (FSTobs) and the n=10 random replicates
* bFST = FSTobs - maximum value observed in the 10 random replicates 

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=3
#SBATCH --time=48:00:00

#mamba activate R
# for egg in $(cat Egg_Types.list); do sbatch -J FST_${egg} FST-Background.sh ${egg}; done

if [ -z "$1" ]; then
    echo "Error: provide population ID"
    exit 1
fi

mkdir -p bgfst/work bgfst/out

TARGET=$1
WIN_SIZE=10000

if [[ $TARGET =~ ^ECO ]]; then
    SPECIES="CO"
else
    SPECIES="CC"
fi

TOTAL_SAMPLES=$(cat ${SPECIES}.pop | wc -l )
TARGET_SAMPLES=$(awk -v e=$TARGET 'awk $2 == e' ${SPECIES}.pop | wc -l)
BACKGROUND_SAMPLES=$((TOTAL_SAMPLES - TARGET_SAMPLES))
awk -v e=$TARGET 'awk $2 == e' ${SPECIES}.pop | awk '{print $1}' > bgfst/work/${TARGET}.TRUE-T.list
awk -v e=$TARGET 'awk $2 != e' ${SPECIES}.pop | awk '{print $1}' > bgfst/work/${TARGET}.TRUE-F.list
echo -e "\e[43m~~~~ FOR ${TARGET} THERE ARE ${TARGET_SAMPLES} TARGETS AND ${BACKGROUND_SAMPLES} BACKGROUND ~~~~\e[0m"

for CHR in $(cat Chromosomes.list); do

    echo -e "\e[42m~~~~ WORKING ON ${TARGET} AND ${CHR} ~~~~\e[0m"

    if [[ $CHR = 'chr_W' || $CHR = 'chr_MT' ]]; then

        echo "WORKING ON HAPLOID: True comparison"
        vcftools --haploid --gzvcf vcfs/${CHR}.SNP.DP3.vcf.gz --out bgfst/work/${CHR}_${TARGET}_TRUE \
            --weir-fst-pop bgfst/work/${TARGET}.TRUE-T.list --weir-fst-pop bgfst/work/${TARGET}.TRUE-F.list --fst-window-size ${WIN_SIZE} --max-missing 0.1 2> bgfst/work/${CHR}_${TARGET}.T.log
        awk -v e=${TARGET} '{OFS="\t"}{print $1, $2, $3, $5, e, "T", "T"}' bgfst/work/${CHR}_${TARGET}_TRUE.windowed.weir.fst | sed '1d' > bgfst/work/${CHR}_${TARGET}_TRUE.out
        echo "WORKING ON HAPLOID: Background comparison"
        for IT in $(seq 1 10); do
            echo "WORKING ON HAPLOID: Background comparison iteration ${IT}"
            cat ${SPECIES}.pop | shuf > bgfst/work/${TARGET}.${IT}.list
            head -n ${TARGET_SAMPLES} bgfst/work/${TARGET}.${IT}.list > bgfst/work/${TARGET}.${IT}-T.list
            tail -n +"$((TARGET_SAMPLES + 1))" bgfst/work/${TARGET}.${IT}.list > bgfst/work/${TARGET}.${IT}-F.list
            vcftools --haploid --gzvcf vcfs/${CHR}.SNP.DP3.vcf.gz --out bgfst/work/${CHR}_${TARGET}_F${IT} \
                --weir-fst-pop bgfst/work/${TARGET}.${IT}-T.list --weir-fst-pop bgfst/work/${TARGET}.${IT}-F.list --fst-window-size ${WIN_SIZE} --max-missing 0.1 2> bgfst/work/${CHR}_${TARGET}.${IT}.log
            awk -v e=${TARGET} -v i=${IT} '{OFS="\t"}{print $1, $2, $3, $5, e, "BG", i}' bgfst/work/${CHR}_${TARGET}_F${IT}.windowed.weir.fst | sed '1d' > bgfst/work/${CHR}_${TARGET}_F${IT}.out

        done

        # Calculate 95% CIs for each SNP and merge with the true comparison
        cat bgfst/work/${CHR}_${TARGET}_*out > bgfst/work/${CHR}_${TARGET}.all.txt
        Rscript FSTBackgroundCalculations.R ${CHR} ${TARGET}

    else

        echo "WORKING ON DIPLOID: True comparison"
        vcftools --gzvcf vcfs/${CHR}.SNP.DP3.vcf.gz --out bgfst/work/${CHR}_${TARGET}_TRUE \
            --weir-fst-pop bgfst/work/${TARGET}.TRUE-T.list --weir-fst-pop bgfst/work/${TARGET}.TRUE-F.list --fst-window-size ${WIN_SIZE} --max-missing 0.1 2> bgfst/work/${CHR}_${TARGET}.T.log
        awk -v e=${TARGET} '{OFS="\t"}{print $1, $2, $3, $5, e, "T", "T"}' bgfst/work/${CHR}_${TARGET}_TRUE.windowed.weir.fst | sed '1d' > bgfst/work/${CHR}_${TARGET}_TRUE.out
        echo "WORKING ON DIPLOID: Background comparison"
        for IT in $(seq 1 10); do
            echo "WORKING ON DIPLOID: Background comparison iteration ${IT}"
            cat ${SPECIES}.pop | shuf > bgfst/work/${TARGET}.${IT}.list
            head -n ${TARGET_SAMPLES} bgfst/work/${TARGET}.${IT}.list > bgfst/work/${TARGET}.${IT}-T.list
            tail -n +"$((TARGET_SAMPLES + 1))" bgfst/work/${TARGET}.${IT}.list > bgfst/work/${TARGET}.${IT}-F.list
            vcftools --gzvcf vcfs/${CHR}.SNP.DP3.vcf.gz --out bgfst/work/${CHR}_${TARGET}_F${IT} \
                --weir-fst-pop bgfst/work/${TARGET}.${IT}-T.list --weir-fst-pop bgfst/work/${TARGET}.${IT}-F.list --fst-window-size ${WIN_SIZE} --max-missing 0.1 2> bgfst/work/${CHR}_${TARGET}.${IT}.log
            awk -v e=${TARGET} -v i=${IT} '{OFS="\t"}{print $1, $2, $3, $5, e, "BG", i}' bgfst/work/${CHR}_${TARGET}_F${IT}.windowed.weir.fst | sed '1d' > bgfst/work/${CHR}_${TARGET}_F${IT}.out

        done

        # Calculate 95% CIs for each SNP and merge with the true comparison
        cat bgfst/work/${CHR}_${TARGET}_*out > bgfst/work/${CHR}_${TARGET}.all.txt
        Rscript FSTBackgroundCalculations.R ${CHR} ${TARGET}
    fi

done

```

Have this script accessible `FSTBackgroundCalculations.R`: 

```bash
#!/usr/bin/env rscript
args <- commandArgs(trailingOnly = TRUE)

setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/scan_n4/bgfst/work')
library(dplyr)
library(meRo)

options(scipen=999)

chr = args[1]
egg = args[2]

d = read.table(paste0(chr,'_',egg,'.all.txt'),header=F)
names(d) = c('chr','start','end','fst','egg','fork','iteration')
d$fst = pmax(0, pmin(1, d$fst))
dt = d %>% filter(fork == 'T')
bg = d %>% filter(fork == 'BG')
cis = bg %>% group_by(chr,start,end,egg,fork) %>% sum_stats(fst) %>% select(chr,start,end,egg,lo=conf_low,hi=conf_high,mean=mean,sd=sd)
cis$lo = pmax(0, pmin(1, cis$lo))
cis$hi = pmax(0, pmin(1, cis$hi))
cis$mean = pmax(0, pmin(1, cis$mean))
maxs = bg %>% group_by(chr,start,end,egg) %>% summarize(max=max(fst))

dtf = left_join(dt,cis) %>% left_join(.,maxs)
dtf$zfst = (dtf$fst - dtf$mean) / (dtf$sd + .Machine$double.eps)
dtf = dtf %>% mutate(across(where(is.numeric), ~ round(.x, 4))) %>% select(chr,start,end,egg,fst,lo,hi,max,mean,sd,zfst)

write.table(dtf,file=paste0('../out/',chr,'_',egg,'.RESULTS.txt'),quote=F,sep='\t',row.names=F,col.names=F)

```

#### Plot bFST: Windows

```bash
#### Plot 10KB FST Gene Scan with Sweeps
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/scan_n4')
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

md = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt')

##### Load and prep FST   ######
tidy_bp = fread('20250325_bFST.txt.gz')
names(tidy_bp) = c('chr','start','end','Egg','FST','lo','hi','max','mean','sd','zfst')
setDT(tidy_bp)
tidy_bp <- tidy_bp %>% select(-lo,-hi,-mean,-zfst,-sd)

#strip chr_, and assign AvZ
tidy_bp <- tidy_bp[
  , chr := gsub('chr_', '', chr)][
    , AvZ := fifelse(chr == 'W', 'W',
                     fifelse(chr == 'Z', 'Z',
                             fifelse(chr == 'MT', 'MT', 'Autosome')))]

#calculate bFST for each window (fst[obs] - max(fst[randomized])) 
tidy_bp <- tidy_bp %>% mutate(FST = pmax(0,pmin(1,FST)),
                              bFST = pmax(0,pmin(1,FST-max)),
                              Species = ifelse(grepl('ECC',Egg,),'CC','CO')) %>% 
  ungroup %>% 
  group_by(chr,start,Species) %>% 
  mutate(second_highest_bFST = nth(sort(bFST, decreasing = TRUE), 2, default = NA_real_),
         dbFST = pmax(0, bFST - second_highest_bFST)) %>% 
  ungroup 

tidy_bp %>% filter(chr == 'MT') %>% arrange(start)

# Formatting
eggcols <- md %>% arrange(EggOrder) %>% select(Egg,EggCol) %>% unique %>% filter(Egg %in% unique(tidy_bp$Egg))
tidy_bp$Egg <- factor(tidy_bp$Egg, levels=eggcols$Egg)
tidy_bp$AvZ <- factor(tidy_bp$AvZ,levels=c('Autosome','Z','MT','W'))

# Summaries
tidy_bp %>% 
  pivot_longer(c(bFST,dbFST)) %>% 
  group_by(Egg,AvZ,name) %>% 
  sum_stats(value) %>% 
  ggplot(aes(y=Egg,x=mean,xmin=conf_low,xmax=conf_high,col=AvZ))+
  geom_point(position=position_dodge(width=0.5))+
  geom_errorbar(width=0.5,position=position_dodge(width=0.5))+
  theme_classic()+
  facet_wrap(name~.,scales='free')

wins <- tidy_bp %>% 
  group_by(Egg) %>% 
  filter(bFST >= quantile(bFST, 0.999,na.rm=TRUE)) %>% 
  count(Egg,AvZ) %>% 
  ungroup() %>% 
  complete(Egg, AvZ, fill = list(n = 0)) %>% 
  ggplot(aes(x = Egg, fill = AvZ, y=n, label = n)) +
  geom_bar(col='black',stat = "identity", aes(y = n), position = position_dodge(width = 0.9)) + 
  geom_text(position = position_dodge(width = 0.9),size=3,vjust=-.25)+
  scale_fill_manual(values = viridis(4,option='mako')) +
  #facet_grid(AvZ~.,scales='free')+
  theme_bw(base_size = 8) +
  coord_cartesian(ylim=c(0,150))+
  ylab('') +
  xlab('') +
  theme(legend.position = "right", 
        legend.text = element_text(size = 6), 
        legend.title = element_text(size = 6)) + 
  guides(fill = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 1)))
wins
ggsave('~/symlinks/host/figures/20250401_bFSTWindowDistributions_AvZ.pdf',
       wins,height=1.25,width=7.5,dpi=300)

##### Plot GWAS #####
library(karyoploteR)
#prep karyoplot
genome = read.table('/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa.bed',header=FALSE) %>% filter(str_detect(V1,'scaff',negate=T)) %>% arrange(desc(V2))
names(genome) = c('chr','start','end'); genome$chr = gsub('chr_','',genome$chr); 
genome = genome %>% mutate(chr = factor(chr, levels = c(as.character(1:39), "Z", "W", "MT"))) %>%  arrange(chr)
G = makeGRangesFromDataFrame(genome,seqnames.field = 'chr',start.field = 'start',end.field = 'end')

#add colors 
chrs = genome %>% select(chr) %>% mutate(Color = if_else(row_number() %% 2 == 0, 'grey60', 'black'))

# Plot huge all SNPs
karyoALL <- tidy_bp %>% filter(!grepl('scaffold',chr))

# Change colors
genome <- genome %>% mutate(color = rep(c("black", "grey70"), length.out = n()))
karyoA <- left_join(karyoALL,genome %>% select(chr,color))
karyoA <- karyoA %>% mutate(color = ifelse(grepl('^ECO',Egg) & color == 'grey70','grey40', color))

pdf('~/symlinks/host/figures/20250327_KARYOPLOT-bFST10KB-MAXlabels.pdf',height=4,width=7.5)
png('~/symlinks/host/figures/20250327_KARYOPLOT-bFST10KB-MAX.png',height=4,width=7.5,units='in',res=300)
png('~/symlinks/host/figures/20250327_KARYOPLOT-bFST10KB-MAX-Wonly.png',height=4,width=7.5,units='in',res=300)
pp = getDefaultPlotParams(plot.type=4)
pp$leftmargin = 0.1
kp = plotKaryotype(plot.type=4, 
                   genome = G,
                   #genome = keepSeqlevels(G, "W", pruning.mode = "coarse"),
                   labels.plotter = NULL,
                   plot.params = pp)
kpAddChromosomeNames(kp, yoffset = -5,cex=0.7)
kpAddBaseNumbers(kp,tick.dist = 10000000,minor.ticks = FALSE,cex=0.4)

#Loop through each track (all species together, CC, and CO), and plot them on their own layer 
counter = 0; tracks = 6
for (egg in rev(eggcols$Egg)) { 
  
  counter = counter + 1; at = autotrack(current.track = counter, total.tracks = tracks);at$r1 <- at$r1-0.02
  
  # Grab target egg 
  cat('Plotting: ',egg,'\n')
  sub <- karyoA %>% filter(Egg == egg) 
  ymin=0;ymax=1
  #ymin=min(sub$bFST,na.rm=TRUE);ymax=max(sub$bFST,na.rm=TRUE)
  kpPoints(kp,chr=sub$chr,x=sub$start+5000,y=sub$bFST,ymin=ymin,ymax=ymax,r0=at$r0,r1=at$r1,col=as.character(sub$color))
  kpAxis(kp,r0=at$r0,r1=at$r1,cex=0.5,numticks = 2,ymin=ymin,ymax=ymax)
  kpAddLabels(kp,cex=0.5,labels = egg,r0=at$r0+.01, r1=at$r1,col="black",srt=0,label.margin = 0.02)
  
}

dev.off()

#### EXAMINE REGIONS OF INTEREST ####
outliers <- tidy_bp %>%
  group_by(Egg,AvZ) %>%  # Only group by Egg, not chr
  filter(!grepl('scaf|MT|W',chr)) %>% 
  slice_max(order_by = bFST, n = 1, with_ties = TRUE) %>%  # Get top 5 per Egg
  ungroup() %>%
  select(chr, start, end, Egg, bFST) %>%
  data.frame() %>% 
  mutate(chr = gsub('chr_scaffold','scaffold',paste0('chr_',chr)))
outliers
write.table(outliers,file='20250331_bFST_outliers.txt',quote=F,sep='\t',row.names=F,col.names=F)
#bedtools sort -i 20250331_bFST_outliers.txt > 20250331_bFST_outliers.sorted.txt


##### ECC1 #####
# ECC1 chrZ : BNC2
deets <- tidy_bp %>% filter(Egg == 'ECC1' & !grepl('scaf',chr)) %>% group_by(chr) %>% slice_max(bFST) %>% arrange(desc(bFST)) %>% head(3) %>% select(chr,start,end)
deets
zchr='Z'
zstart=23520001
pe1 <- tidy_bp %>% filter(chr == zchr & start > zstart-5e5 & end < zstart+5e5) %>%
  ggplot(aes(x=start,y=bFST,col=Egg))+
  geom_text(x=zstart+5000,y=Inf,label='*',col='blue',alpha=1,size=3,vjust=1)+
  scale_color_manual(values=md$EggCol,breaks=md$Egg)+
  ggtitle(paste0(zchr,' ',zstart,'ECC1:BNC2'))+
  geom_line()+scale_x_continuous(breaks = pretty_breaks(n = 3),labels = function(x) paste0(round(x / 1e6, 1), "-Mb")) +scale_y_continuous(breaks=pretty_breaks(n=3))+
  theme_classic(base_size=8)+ylab('')+xlab('')
pe1

##### ECC6 #####
# ECC6 chr2 : BLVRA
deets <- tidy_bp %>% filter(Egg == 'ECC6' & !grepl('scaf',chr)) %>% group_by(chr) %>% slice_max(bFST) %>% arrange(desc(bFST)) %>% head(3) %>% select(chr,start,end)
deets
zchr=2
zstart=104840001
p1 <- tidy_bp %>% filter(chr == zchr & start > zstart-5e5 & end < zstart+5e5) %>%
  ggplot(aes(x=start,y=bFST,col=Egg))+
  geom_rect(aes(xmin = 104902858, xmax = 104937472, ymin = -Inf, ymax = Inf), fill = "gray90", alpha = 0.9,col=NA) + # blvra
  geom_rect(aes(xmin = 104765334, xmax = 104839609, ymin = -Inf, ymax = Inf), fill = "gray90", alpha = 0.9,col=NA) + # vopp1
  geom_text(x=zstart+5000,y=Inf,label='*',col='blue',alpha=1,size=3,vjust=1)+
  scale_color_manual(values=md$EggCol,breaks=md$Egg)+
  ggtitle(paste0(zchr,' ',zstart,'ECC6:BLVRA,VOPP1'))+
  geom_line()+scale_x_continuous(breaks = pretty_breaks(n = 3),labels = function(x) paste0(round(x / 1e6, 1), "-Mb")) +scale_y_continuous(breaks=pretty_breaks(n=3))+
  theme_classic(base_size=8)+ylab('')+xlab('')
p1


# ECC6 chr2 : BLVRA
deets <- tidy_bp %>% filter(Egg == 'ECC6' & !grepl('scaf',chr)) %>% group_by(AvZ) %>% slice_max(bFST) %>% arrange(desc(bFST)) %>% head(3) %>% select(chr,start,end)
deets
zchr='Z'
zstart=18030001
p99 <- tidy_bp %>% filter(chr == zchr & start > zstart-5e5 & end < zstart+5e5) %>%
  ggplot(aes(x=start,y=bFST,col=Egg))+
  geom_rect(aes(xmin = 17935600, xmax = 18066571, ymin = -Inf, ymax = Inf), fill = "gray90", alpha = 0.9,col=NA) + # blvra
  geom_rect(aes(xmin = 18025552, xmax = 18032485, ymin = -Inf, ymax = Inf), fill = "gray90", alpha = 0.9,col=NA) + # vopp1
  geom_text(x=zstart+5000,y=Inf,label='*',col='blue',alpha=1,size=3,vjust=1)+
  scale_color_manual(values=md$EggCol,breaks=md$Egg)+
  ggtitle(paste0(zchr,' ',zstart,'ECC6:IQGAP2'))+
  geom_line()+scale_x_continuous(breaks = pretty_breaks(n = 3),labels = function(x) paste0(round(x / 1e6, 1), "-Mb")) +scale_y_continuous(breaks=pretty_breaks(n=3))+
  theme_classic(base_size=8)+ylab('')+xlab('')
p99

##### ECC10 #####
# ECC10 chrZ : RLN3
deets <- tidy_bp %>% filter(Egg == 'ECC10' & !grepl('scaf',chr)) %>% group_by(chr) %>% slice_max(bFST) %>% arrange(desc(bFST)) %>% head(3) %>% select(chr,start,end)
deets
zchr='Z'
zstart=18980001
p2 <- tidy_bp %>% filter(chr == zchr & start > zstart-2e5 & end < zstart+2.5e5) %>%
  ggplot(aes(x=start,y=bFST,col=Egg))+
  geom_rect(aes(xmin = 18975843, xmax = 18979322, ymin = -Inf, ymax = Inf), fill = "grey90", alpha = 0.4,col=NA) + # rln3
  geom_rect(aes(xmin = 18987657, xmax = 19020048, ymin = -Inf, ymax = Inf), fill = "grey90", alpha = 0.4,col=NA) + # rln3
  geom_text(x=zstart+5000,y=Inf,label='*',col='blue',alpha=1,size=3,vjust=1)+
  scale_color_manual(values=md$EggCol,breaks=md$Egg)+
  ggtitle(paste0(zchr,' ',zstart,'ECC10:RLN3'))+
  geom_line()+scale_x_continuous(breaks = pretty_breaks(n = 3),labels = function(x) paste0(round(x / 1e6, 1), "-Mb")) +scale_y_continuous(breaks=pretty_breaks(n=3))+
  theme_classic(base_size=8)+ylab('')+xlab('')
p2

##### ECO1 #####
# ECO1 chr6 : BOLL
deets <- tidy_bp %>% filter(Egg == 'ECO1' & !grepl('scaf',chr)) %>% group_by(chr) %>% slice_max(bFST) %>% arrange(desc(bFST)) %>% head(3) %>% select(chr,start,end)
deets
zchr=6
zstart=31140001
zstart_second=29970001
p3 <- tidy_bp %>% filter(chr == zchr & start > zstart-1.5e6 & end < zstart+5e5) %>%
  ggplot(aes(x=start,y=bFST,col=Egg))+
  geom_rect(aes(xmin = 31138680, xmax = 31156397, ymin = -Inf, ymax = Inf), fill = "gray90", alpha = 0.9,col=NA) +  #boll
  geom_rect(aes(xmin = 29965430, xmax = 29971141, ymin = -Inf, ymax = Inf), fill = "gray90", alpha = 0.9,col=NA) +  #boll
  geom_text(x=zstart+5000,y=Inf,label='*',col='blue',alpha=1,size=3,vjust=1)+
  geom_text(x=zstart_second+5000,y=Inf,label='*',col='blue',alpha=1,size=3,vjust=1)+
  scale_color_manual(values=md$EggCol,breaks=md$Egg)+
  ggtitle(paste0(zchr,' ',zstart,'ECO1:BOLL'))+
  geom_line()+scale_x_continuous(breaks = pretty_breaks(n = 3),labels = function(x) paste0(round(x / 1e6, 1), "-Mb")) +scale_y_continuous(breaks=pretty_breaks(n=3))+
  theme_classic(base_size=8)+ylab('')+xlab('')
p3

# Second peak?
regiondeets <- tidy_bp %>% filter(Egg == 'ECO1' & chr == zchr & start > zstart-1.5e6 & end < zstart+5e5) %>% slice_max(bFST,n=5) %>% arrange(bFST)
regiondeets

# CPLX1  
zchr='Z'
zstart=45190001
pn <- tidy_bp %>% filter(chr == zchr & start > zstart-5e5 & end < zstart+5e5) %>%
  ggplot(aes(x=start,y=bFST,col=Egg))+
  ggtitle(paste0(zchr,' ',zstart,'ECO1:CPLX1'))+
  geom_text(x=zstart+5000,y=Inf,label='*',col='blue',alpha=1,size=3,vjust=1)+
  scale_color_manual(values=md$EggCol,breaks=md$Egg)+
  geom_line()+scale_x_continuous(breaks = pretty_breaks(n = 3),labels = function(x) paste0(round(x / 1e6, 1), "-Mb")) +scale_y_continuous(breaks=pretty_breaks(n=3))+
  theme_classic(base_size=8)+ylab('')+xlab('')
pn

##### ECO3 #####
# ECO3 chrZ : HINT1
deets <- tidy_bp %>% filter(Egg == 'ECO3' & !grepl('scaf',chr)) %>% group_by(chr) %>% slice_max(bFST) %>% arrange(desc(bFST)) %>% head(5) %>% select(chr,start,end)
deets
zchr='Z'
zstart=34540001
p4 <- tidy_bp %>% filter(chr == zchr & start > zstart-3e5 & end < zstart+3e5) %>%
  ggplot(aes(x=start,y=bFST,col=Egg))+
  geom_rect(aes(xmin = 34549415, xmax = 34556113, ymin = -Inf, ymax = Inf), fill = "gray90", alpha = 0.9,col=NA) +  #boll
  scale_color_manual(values=md$EggCol,breaks=md$Egg)+
  geom_line()+scale_x_continuous(breaks = pretty_breaks(n = 3),labels = function(x) paste0(round(x / 1e6, 1), "-Mb")) +scale_y_continuous(breaks=pretty_breaks(n=3))+
  ggtitle(paste0(zchr,' ',zstart,'ECO3:HINT1'))+
  geom_text(x=zstart+5000,y=Inf,label='*',col='blue',alpha=1,size=3,vjust=1)+
  theme_classic(base_size=8)+ylab('')+xlab('')
p4

##### ECO4 #####
# ECO4 chr 5 TSH
deets <- tidy_bp %>% filter(Egg == 'ECO4' & !grepl('scaf',chr)) %>% group_by(chr) %>% slice_max(bFST) %>% arrange(desc(bFST)) %>% head(5) %>% select(chr,start,end)
deets
zchr=5
zstart=19500001
p5 <- tidy_bp %>% filter(chr == zchr & start > zstart-5e5 & end < zstart+5e5) %>%
  ggplot(aes(x=start,y=bFST,col=Egg))+
  geom_rect(aes(xmin = 19487067, xmax = 19544653, ymin = -Inf, ymax = Inf), fill = "gray90", alpha = 0.9,col=NA) +
  scale_color_manual(values=md$EggCol,breaks=md$Egg)+
  ggtitle(paste0(zchr,' ',zstart,'ECO4:TSHR'))+
  geom_text(x=zstart+5000,y=Inf,label='*',col='blue',alpha=1,size=3,vjust=1)+
  geom_line()+scale_x_continuous(breaks = pretty_breaks(n = 3),labels = function(x) paste0(round(x / 1e6, 1), "-Mb")) +scale_y_continuous(breaks=pretty_breaks(n=3))+
  theme_classic(base_size=8)+ylab('')+xlab('')
p5

# ECO4 chr28 HCN3
deets
zchr=28
zstart=3940001
zstart_3rd=3670001

p6 <- tidy_bp %>% filter(chr == zchr & start > zstart-3.5e6 & end < zstart+8e6) %>%
  ggplot(aes(x=start,y=bFST,col=Egg))+
  geom_rect(aes(xmin = 3941523, xmax = 3958025, ymin = -Inf, ymax = Inf), fill = "gray90", alpha = 0.9,col=NA) + # CLK2
  geom_rect(aes(xmin = 3686372, xmax = 3702332, ymin = -Inf, ymax = Inf), fill = "gray90", alpha = 0.9,col=NA) + # MCL1
  geom_rect(aes(xmin = 3928875, xmax = 3941209, ymin = -Inf, ymax = Inf), fill = "gray90", alpha = 0.9,col=NA) + # HCN3
  #geom_rect(aes(xmin = 2353765, xmax = 2783871, ymin = -Inf, ymax = Inf), fill = "gray90", alpha = 0.9,col=NA) + #weird keratin repeats
  scale_color_manual(values=md$EggCol,breaks=md$Egg)+
  ggtitle(paste0(zchr,' ',zstart,'ECO4:MCL-1'))+
  geom_text(x=zstart+5000,y=Inf,label='*',col='blue',alpha=1,size=3,vjust=1)+
  geom_text(x=zstart_3rd+5000,y=Inf,label='*',col='blue',alpha=1,size=3,vjust=1)+
  geom_line()+scale_x_continuous(breaks = pretty_breaks(n = 3),labels = function(x) paste0(round(x / 1e6, 1), "-Mb")) +scale_y_continuous(breaks=pretty_breaks(n=3))+
  theme_classic(base_size=8)+ylab('')+xlab('')
p6

regiondeets <- tidy_bp %>% filter(Egg == 'ECO4' & chr == zchr & end < 3.9e6) %>% slice_max(bFST,n=5) %>% arrange(bFST) 
regiondeets

ggarrange(p1,p2,p3,p4,p5,p6,nrow=3,ncol=2,common.legend = TRUE)
pdf('~/symlinks/host/figures/20250404_Candidates.pdf',height=3.5,width=6)
ggarrange(p1,p2,p3,p4,p5,p6,nrow=3,ncol=2,common.legend = TRUE)
dev.off()


### Import sweeps
tidy_sweep = fread('20250401_SWEEPS.txt.gz')
names(tidy_sweep) = c('chr','Group','mid','start','end','var','sfs','ld','u')
setDT(tidy_sweep)
tidy_sweep <- tidy_sweep[
  , chr := gsub('chr_', '', chr)][
    , AvZ := fifelse(chr == 'W', 'W',
                     fifelse(chr == 'Z', 'Z',
                             fifelse(chr == 'MT', 'MT', 'Autosome')))]

# Average u in 10kb windows 
window_size <- 1e4
tidy_sweep[, window := floor(start / window_size) * window_size+1]
sweep_avg <- tidy_sweep[, .(mean_u = mean(u, na.rm = TRUE)), 
                        by = .(chr, window, Group)]

# Assign just p1 p2 p3 etc for each population for lty symbology 
sweep_avg <- sweep_avg %>% separate(Group,into=c('Egg','Population'),remove=F) 
pop_id <- sweep_avg %>% select(Egg,Population) %>% distinct %>% group_by(Egg) %>% arrange(Egg) %>%  mutate(PopID = paste0('u',row_number()))
sweeps_pop <- left_join(sweep_avg,pop_id) %>% select(chr,start=window,Group,Egg,Population,PopID,mean_u)

# For visualization, scale mu to bFST. For later on, we also scale it across egg types so that they can be visualized on the same plot (escaled_mean_u)
sweeps <- sweeps_pop %>%
  group_by(Group,Egg,PopID) %>%
  mutate(
    # min/max for mean_u
    min_u = min(mean_u, na.rm = TRUE),
    max_u = max(mean_u, na.rm = TRUE),
    
    # Scale mean_u between 0 and 1 for each Egg
    escaled_mean_u = ifelse(
      min_u == max_u, 
      0.5,  # neutral value if all mean_u are the same
      (mean_u - min_u) / (max_u - min_u))) %>%
  select(-min_u, -max_u) %>%  # don't need this garbage
  ungroup()

# Sanity, check that the scaling works, escaled_mean_u ranges from 0 - 1
sweeps %>%
  group_by(Population,PopID,Egg) %>% 
  summarize(min = min(escaled_mean_u,na.rm=TRUE),max=max(escaled_mean_u,na.rm=TRUE))

# identify sweeps within each egg 
sweepdat <- NULL
for (grp in unique(sweeps$Group)) { 
  wins <- sweeps %>% filter(Group == grp) %>% 
    filter(escaled_mean_u >= quantile(escaled_mean_u, 0.9,na.rm=TRUE)) %>% 
    select(chr,start,Egg)
  sweepdat <- rbind(sweepdat,wins)
}

sweeps %>% mutate(FST = pmax(0,pmin(1,FST)),
                  bFST = pmax(0,pmin(1,FST-max)),
                  Species = ifelse(grepl('ECC',Egg,),'CC','CO')) %>% 
  ungroup %>% 
  group_by(chr,start,Egg) %>% 
  mutate(second_highest_bFST = nth(sort(bFST, decreasing = TRUE), 2, default = NA_real_),
         dbFST = pmax(0, bFST - second_highest_bFST)) %>% 
  ungroup 


pdf('~/symlinks/host/figures/20250327_KARYOPLOT-bFST10KB-MAXlabels.pdf',height=4,width=7.5)
png('~/symlinks/host/figures/20250327_KARYOPLOT-bFST10KB-MAX.png',height=4,width=7.5,units='in',res=300)
png('~/symlinks/host/figures/20250327_KARYOPLOT-bFST10KB-MAX-Wonly.png',height=4,width=7.5,units='in',res=300)
pp = getDefaultPlotParams(plot.type=4)
pp$leftmargin = 0.1
kp = plotKaryotype(plot.type=4, 
                   genome = G,
                   #genome = keepSeqlevels(G, "W", pruning.mode = "coarse"),
                   labels.plotter = NULL,
                   plot.params = pp)
kpAddChromosomeNames(kp, yoffset = -5,cex=0.7)
kpAddBaseNumbers(kp,tick.dist = 10000000,minor.ticks = FALSE,cex=0.4)

#Loop through each track (all species together, CC, and CO), and plot them on their own layer 
counter = 0; tracks = 10
for (grp in unique(sweeps$Group)) { 
  
  counter = counter + 1; at = autotrack(current.track = counter, total.tracks = tracks);at$r1 <- at$r1-0.02
  
  # Grab target group
  cat('Plotting: ',grp,'\n')
  subs <- sweeps %>% filter(Group == grp) 
  ymin=0;ymax=1
  eggcol <- eggcols %>% filter(Egg %in% subs$Egg) %>% pull(EggCol)
  
  # Add sweeps
  kpLines(kp,chr=subs$chr,x=subs$start+5000,y=subs$escaled_mean_u,ymin=ymin,ymax=ymax,r0=at$r0,r1=at$r1,col=eggcol)
  kpAxis(kp,r0=at$r0,r1=at$r1,cex=0.5,numticks = 2,ymin=ymin,ymax=ymax)
  kpAddLabels(kp,cex=0.5,labels = grp,r0=at$r0+.01, r1=at$r1,col="black",srt=0,label.margin = 0.02)
  
}

dev.off()

```

#### Genes



Chr2, ECC6, [chicken oxidative phosphorylation](https://www.nature.com/articles/s41598-021-04077-y) (COA1):  https://www.ncbi.nlm.nih.gov/gene/104065962

Chr2, ECC6 ,  [oxidative response mtDNA dysfunction](https://doi.org/10.1038/labinvest.2011.70) VOPP1): https://www.ncbi.nlm.nih.gov/gene/104065935

Chr2, ECC6, [general biliverdin](https://doi.org/10.1086/694297)  and [BLVRA oxidative response]([10.1016/j.freeradbiomed.2017.11.020](https://doi.org/10.1016/j.freeradbiomed.2017.11.020)) (BLVRA): https://www.ncbi.nlm.nih.gov/gene/104065936

ChrZ, ECC10, [chicken pituitary](https://pmc.ncbi.nlm.nih.gov/articles/PMC9676923/) (RLN3): https://www.ncbi.nlm.nih.gov/gene/104063035

Chr6, ECO1, [germline](https://doi.org/10.1016/j.scr.2017.04.008) (BOLL): https://www.ncbi.nlm.nih.gov/gene/104058472

Chr6, ECO1, psuedogene (CLK1): https://www.ncbi.nlm.nih.gov/gene/128852493

ChrZ, ECO3, [chicken reproduction](https://www.sciencedirect.com/science/article/abs/pii/S0093691X10002712) (HINT1):  

chrZ, ECC6, [chicken follicle](https://www.sciencedirect.com/science/article/pii/S0032579123001244) (IQGAP2): https://www.ncbi.nlm.nih.gov/gene/104056386

Chr28, ECO4, [chicken egg development](https://pubmed.ncbi.nlm.nih.gov/38652954/) (HCN3): https://www.ncbi.nlm.nih.gov/gene/104057086

Chr28, ECO4, [reproduction chickens follicle](https://pmc.ncbi.nlm.nih.gov/articles/PMC4669721/) (MCL-1): https://www.ncbi.nlm.nih.gov/gene/104057095

Chr28, ECO4,[circadian rhythm plants](https://www.pnas.org/doi/10.1073/pnas.96.22.12362?url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org&rfr_dat=cr_pub++0pubmed) and [mammals](https://www.science.org/doi/10.1126/scisignal.2000305?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed) (CLK2): https://www.ncbi.nlm.nih.gov/gene/104057075

Chr28, ECO4, [chicken](mitochondrion) (C28H1orf43): https://www.ncbi.nlm.nih.gov/gene/104061394

Chr5, ECO4, [chickens](https://pmc.ncbi.nlm.nih.gov/articles/PMC8045693) (TSHR): https://www.ncbi.nlm.nih.gov/gene/104059620



Check for excess linkage (e.g. inversions within those regions):

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=5
#SBATCH --time=24:00:00

mkdir -p inversion
for SPECIES in CO CC; do

while read -r CHR START END; do
    fstart=$((START - 2000000))
    fend=$((START + 2000000))

    echo "Working on ${CHR} ${START} ${END} with flank ${fstart} ${fend}"
    # Ensure start is not negative
    if [ "$fstart" -lt 1 ]; then
        fstart=1
    fi

    # Including flanking
    bcftools view --samples-file ${SPECIES}.list vcfs/${CHR}.SNP.DP3.vcf.gz --regions ${CHR}:${fstart}-${fend} -Ou | \
      bcftools view --min-alleles 2 --max-alleles 2 --min-ac 2 -e 'F_MISSING > 0.05' -Oz \
        -o inversion/${CHR}.${SPECIES}.inversion.flank.vcf.gz

    # Only within region
    bcftools view --samples-file ${SPECIES}.list vcfs/${CHR}.SNP.DP3.vcf.gz --regions ${CHR}:${fstart}-${fend} -Ou | \
      bcftools view --min-alleles 2 --max-alleles 2 --min-ac 2 -e 'F_MISSING > 0.05' -Oz \
        -o inversion/${CHR}.${SPECIES}.inversion.vcf.gz

    plink --allow-extra-chr --double-id --vcf inversion/${CHR}.${SPECIES}.inversion.flank.vcf.gz --r2 --out inversion/${CHR}_${SPECIES}_LD --ld-window 999999999 --ld-window-kb 1000000000

    plot_ld --input inversion/${CHR}_${SPECIES}_LD.ld --out inversion/${CHR}_${SPECIES}.png --win_size 10 --highlight ${START}-${END}
    vcf_to_pca --vcf inversion/${CHR}.${SPECIES}.inversion.vcf.gz --metadata ~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt --phenotype Egg --out inversion/${CHR}.${SPECIES}_pca.png

done < regions.bed

done

```

## FST: Base-pair level fixed SNPs [Fig 4]

Hierarchical contrasts, egg:haplotype comparisons:

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



## NDUFAF4 across Cuculiformes

### Fastq Alignment-Based

Species to assay, brood parasitic consortium:

| Species               | Sex    | Tech | Accession   |
| --------------------- | ------ | ---- | ----------- |
| Clamator_glandarius   | Female | HiFi | SRR26807982 |
| Coccyzus_lansbergi    | Female | HiFi | SRR26902471 |
| Dromococcyx_pavoninus | Female | HiFi | SRR26905917 |

All:

| Cuculus  canorus        | Female | Illumina | SRR11394165                          |
| ----------------------- | ------ | -------- | ------------------------------------ |
| Cuculus canorus         | Female | Illumina | SRR11531702                          |
| Cuculus canorus         | Female | Illumina | SRR11531718                          |
| Cuculus canorus         | Male   | Illumina | SRR11531726                          |
| Cuculus canorus         | Male   | Illumina | SRR11531741                          |
| Cuculus micropterus     | Male   | Illumina | SRR14117632                          |
| Cuculus micropterus     | Female | Illumina | SRR14117631                          |
| Cuculus poliocephalus   | Male   | Illumina | SRR14117570                          |
| Cuculus poliocephalus   | Female | Illumina | SRR14117568                          |
| Cuculus canorus         | Female | Illumina | SRR99999129                          |
| Clamator glandarius     | Female | HiFi     | SRR26807982                          |
| Coccyzus lansbergi      | Female | HiFi     | SRR26902471                          |
| Cuculus canorus         | Female | HiFi     | All bCucCan1.pri Reads; PRJNA1008121 |
| Dromococcyx pavoninus   | Female | HiFi     | SRR26905917                          |
| Geococcyx californianus | Female | Illumina | SRR9994302                           |
| Piaya cayana            | Female | Illumina | SRR9947006                           |
| Tapera naevia           | Female | HiFi     | SRR26905805                          |

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
outdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/cuculiformes_NDUFAF4_blasting/Alignment_Based_Approach/Whole_Genome_Alignment/10gb/bams
cuckoo_hifi=/dss/dsslegfs01/pr53da/pr53da-dss-0021/rawdata/External_Data_SRA_GenomeArk/PacBio_GenomeArk/Cuculus_canorus/bCucCan1_pacbio.fastq.gz
subdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/cuculiformes_NDUFAF4_blasting/Alignment_Based_Approach/subset_10gb

# Subset 10gb, and then align, discarding unaligned reads to save space (-F 4)
bbduk.sh in=${cuckoo_hifi} out=${subdir}/${RUN}.10gb.fastq.gz maxbasesout=10000000000
minimap2 -ax map-pb -t 8 ${genome} ${subdir}/${RUN}.10gb.fastq.gz > ${SCRATCH}/${RUN}.sam
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
outdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/cuculiformes_NDUFAF4_blasting/Alignment_Based_Approach/Whole_Genome_Alignment/10gb/bams
subdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/cuculiformes_NDUFAF4_blasting/Alignment_Based_Approach/subset_10gb

# Subset 10gb, and then align, discarding unaligned reads to save space (-F 4)
bbduk.sh in=${qcdata}/${RUN}.trim.fastq.gz out=${subdir}/${RUN}.10gb.fastq.gz maxbasesout=10000000000
bwa mem -M -p -t 10 ${genome} ${subdir}/${RUN}.10gb.fastq.gz | samtools sort -@10 -o ${SCRATCH}/${RUN}.bam -;
samtools view -F 4 -b ${SCRATCH}/${RUN}.bam > ${outdir}/${RUN}.bam
samtools index -b ${outdir}/${RUN}.bam; 
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
outdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/cuculiformes_NDUFAF4_blasting/Alignment_Based_Approach/Whole_Genome_Alignment/10gb/bams
rawdata=/dss/dsslegfs01/pr53da/pr53da-dss-0021/rawdata/External_Data_SRA_GenomeArk/SRA_Cuculiformes_Outgroups
subdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/cuculiformes_NDUFAF4_blasting/Alignment_Based_Approach/subset_10gb

# Subset 10gb, and then align, discarding unaligned reads to save space (-F 4)
bbduk.sh in=${rawdata}/${RUN}_1.fastq.gz out=${subdir}/${RUN}.10gb.fastq.gz maxbasesout=10000000000
minimap2 -ax map-pb -t 8 ${genome} ${subdir}/${RUN}.10gb.fastq.gz > ${SCRATCH}/${RUN}.sam
samtools sort ${SCRATCH}/${RUN}.sam | samtools view -F 4 -b > ${outdir}/${RUN}.bam
samtools index -b ${outdir}/${RUN}.bam; 
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
outdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/cuculiformes_NDUFAF4_blasting/Alignment_Based_Approach/Whole_Genome_Alignment/10gb/bams
subdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/cuculiformes_NDUFAF4_blasting/Alignment_Based_Approach/subset_10gb

# Subset 10gb, and then align, discarding unaligned reads to save space (-F 4)
bbduk.sh in=${rawdata}/${RUN}_1.fastq.gz out=${subdir}/${RUN}.10gb.fastq.gz maxbasesout=10000000000
bwa mem -M -t 10 ${genome} ${subdir}/${RUN}.10gb.fastq.gz | samtools sort -@10 -o ${SCRATCH}/${RUN}.bam -;
samtools view -F 4 -b ${SCRATCH}/${RUN}.bam > ${outdir}/${RUN}.bam
samtools index -b ${outdir}/${RUN}.bam;   
```

Calculate coverage:

```bash
for RUN in $(ls *bam | sed 's/.bam//g'); do 

cds=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/cuculiformes_NDUFAF4_blasting/Alignment_Based_Approach/CDS_Regions_N3_BothParalogs.bed

# Calculate coverage MQ50 
mosdepth --threads 3 --mapq 50 --by ${cds} --fast-mode --no-per-base ../coverage/${RUN}_cds50 ${RUN}.bam

done 

# After, in /coverage/, merge into single file:
for i in $(ls *_cds50.regions.bed | sed 's/_cds50.regions.bed//g'); do awk -v i=${i} '{OFS="\t"}{print $1, $2, $3, $4, i}' ${i}_cds50.regions.bed > ${i}.cds50.cov ; done
cat *.cds50.cov > NDUFAF4_MQ50_Coverage_2024SEPT07.txt
# I then add e.g. species / tech in excel manually since there are not so many fields 
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
```

Ensuring that we only trim and create a tree for the species which have coverage for chrW / chr3: 

```bash
# Align freely
mafft --thread 10 --auto Paralogs_G75.fa > NDUFAF4_freealign.fa

# Trim
trimal -automated1 -in NDUFAF4_freealign.fa -out NDUFAF4_freealign_trimalAuto.fa

# Tree
iqtree --redo -keep-ident -T 20 -s NDUFAF4_freealign_trimalAuto.fa --seqtype DNA -m "MFP" -B 1000

# And then identify gaps 
seqkit stats -a *

# Identify the samples with < 70% gaps, and create a tree:
for i in $(cat Keep.list); do cat ${i} >> NDUFAF4_freealign_trimalAuto.fa; done

# Tree
iqtree --redo -keep-ident -T 20 -s NDUFAF4_freealign_trimalAuto.fa --seqtype DNA -m "GTR" -B 1000
```

Plot Coverage and Tree:

```R
#### Plot NDUFAF4 Coverage 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/cuculiformes_NDUFAF4_blasting/Alignment_Based_Approach/Whole_Genome_Alignment/10gb/coverage/')
.libPaths('~/mambaforge/envs/R/lib/R/library')
library(tidyverse)
library(viridis)

# Load in coverage data 
c <- read_tsv('NDUFAF4_MQ50_CoverageRegions_2024SEPT07.txt')
c$Species <- factor(c$Species,levels=c('Cuculus_canorus','Cuculus_micropterus','Cuculus_poliocephalus','Clamator_glandarius','Coccyzus_lansbergi','Piaya_cayana','Dromococcyx_pavoninus','Geococcyx_californianus'))
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

pdf('../../../../../figures/20240907_NDUFAF4_Coverage.pdf',height=2.75,width=2.25)
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
allcds

pdf('../../../../../figures/NDUFAF4_Coverage-IncludingCDS.pdf',height=6,width=6)
allcds
dev.off()

# Plot Tree afterwards
library(ggtree)
library(treeio)

# Free alignment
iqtree = read.iqtree('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/cuculiformes_NDUFAF4_blasting/Alignment_Based_Approach/Whole_Genome_Alignment/10gb/fastas/results/filtered_outs/NDUFAF4_freealign_trimalAuto.fa.contree')

iqtr = root(as.phylo(iqtree),'Clamator_glandarius_SRR26807982_3')

md <- read_tsv('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/cuculiformes_NDUFAF4_blasting/Alignment_Based_Approach/Whole_Genome_Alignment/bams/Metadata.txt')
g <- ggtree(iqtr,layout='rectangular')
#g <- ggtree(iqtr,layout='ape')
g$data <- g$data %>% mutate(Chromosome = ifelse(grepl('W$',label),'chrW','chr3'),
                            label = gsub('_W$','',label),
                            label = gsub('_3$','',label))
g$data <- left_join(g$data,md %>% dplyr::rename(label = ID))

#plot with outgroups 
gtree <- g +
  geom_tippoint(aes(col=Chromosome,shape=Species),size=1,stroke=1)+
  geom_nodelab(aes(label=label),geom = 'text',size=2,hjust = -3)+ 
  scale_shape_manual(values=c(2,3,1,4,8))+
  xlim(c(0,1))+
  #geom_tiplab(aes(label=Species),size=2,offset = 0.05)+
  geom_tiplab(size=2,offset = 0.02)+
  theme(legend.position='top')
gtree

pdf('../../../../../figures/20240907_NDUFAF4-ParalogTree.pdf',height=2.75,width=2.5)
gtree
dev.off()

```

### Plot Exons

```bash
#### Plot NDUFAF4 genes  
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/')
.libPaths('~/mambaforge/envs/gene_viz/lib/R/library')
library(tidyverse)
library(gggenes)
library(rtracklayer)

# import GFF
gff <- import.gff("/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR-chr_MT.gff")

# filter to features of interest
genes_of_interest <- c("LOC128850245", "NDUFAF4")
gff_df <- as_tibble(gff) %>%
  filter(type == "exon" & gene %in% genes_of_interest) %>%
  mutate(gene = gene,
         strand = as.character(strand),
         start = start,
         end = end,
         length = end - start + 1)

# plot
ggplot(gff_df, aes(xmin = start, xmax = end, y = gene, fill = gene, forward = strand == "+")) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
  theme_classic() +
  labs(x = "Genomic position", y = "") +
  facet_wrap(gene~.,scales='free')+
  theme(legend.position = "none")


# Relative
gff_df_rel <- gff_df %>%
  group_by(gene) %>%
  mutate(rel_start = start - min(start),
         rel_end = end - min(start)) %>%
  ungroup()

# build plot with local coords
gene_plot <- ggplot(gff_df_rel, aes(xmin = rel_start, xmax = rel_end, y = gene, fill = gene, forward = strand == "+")) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
  facet_wrap(~gene, scales = "free", nrow=2) +
  theme_classic() +
  labs(x = "Genomic position (relative)", y = "") +
  theme(legend.position = "none") +
  # add scale bar (e.g., 1kb = 1000 bp) per gene
  geom_segment(data = distinct(gff_df_rel, gene),
               aes(x = 0, xend = 1000, y = -0.5, yend = -0.5), inherit.aes = FALSE) +
  geom_text(data = distinct(gff_df_rel, gene),
            aes(x = 500, y = -1.5, label = "1 kb"), vjust=-5,size = 3, inherit.aes = FALSE)

gene_plot
ggsave('~/symlinks/host/figures/20250425_GeneModelNDUFAF4.pdf',gene_plot,
       height=4,width=5)

```





| Gene ID      | Chromosome | Feature | Start    | End      | Length | Strand | GFF Descriptor                                               |
| ------------ | ---------- | ------- | -------- | -------- | ------ | ------ | ------------------------------------------------------------ |
| NDUFAF4      | chr_3      | gene    | 25511697 | 25515677 | 3980   | -      | NADH:ubiquinone                                              |
| NDUFAF4      | chr_3      | mRNA    | 25511697 | 25515677 | 3980   | -      | NADH:ubiquinone oxidoreductase complex assembly factor 4     |
| NDUFAF4      | chr_3      | exon    | 25515508 | 25515677 | 169    | -      |                                                              |
| NDUFAF4      | chr_3      | exon    | 25514540 | 25514637 | 97     | -      |                                                              |
| NDUFAF4      | chr_3      | exon    | 25511697 | 25512391 | 694    | -      |                                                              |
| NDUFAF4      | chr_3      | CDS     | 25515508 | 25515634 | 126    | -      |                                                              |
| NDUFAF4      | chr_3      | CDS     | 25514540 | 25514637 | 97     | -      |                                                              |
| NDUFAF4      | chr_3      | CDS     | 25512095 | 25512391 | 296    | -      |                                                              |
| LOC128850245 | chr_W      | gene    | 21149910 | 21177468 | 27558  | +      | NADH                                                         |
| LOC128850245 | chr_W      | mRNA    | 21149910 | 21177468 | 27558  | +      | NADH dehydrogenase [ubiquinone] 1 alpha subcomplex assembly factor 4-like |
| LOC128850245 | chr_W      | exon    | 21149910 | 21149957 | 47     | +      |                                                              |
| LOC128850245 | chr_W      | exon    | 21150813 | 21150909 | 96     | +      |                                                              |
| LOC128850245 | chr_W      | exon    | 21177056 | 21177468 | 412    | +      |                                                              |
| LOC128850245 | chr_W      | CDS     | 21149910 | 21149957 | 47     | +      |                                                              |
| LOC128850245 | chr_W      | CDS     | 21150813 | 21150909 | 96     | +      |                                                              |
| LOC128850245 | chr_W      | CDS     | 21177056 | 21177468 | 412    | +      |                                                              |





