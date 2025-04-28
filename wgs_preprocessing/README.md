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

