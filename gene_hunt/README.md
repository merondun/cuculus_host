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



