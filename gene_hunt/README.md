# Analyses: Egg Hunt

Two primary strategies here:

* Estimate background FST in 10KB windows for each egg type with more than 4 females (FSTobs - max(FST[randomized pops n=10])) using a target egg vs all other eggs approach. 
  * Replicate the above, except using a relaxed threshold of at least 3 females per egg for target eggs, and only using nestlings for C. optatus. 
* Estimate base-pair level FST in C. canorus to get blue egg and reverted egg mutations. 

## Preparation

### Sample Selection

```bash
# Egg gene scan, sample selection 
#### Find associations between MT/W haplotypes and features 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/scan_n4')
.libPaths('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/software/mambaforge/envs/r25/lib/R/library')
library(tidyverse)
library(viridis)
library(tidyverse)
library(sf)
library(ggspatial)

md = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt') 
md = md %>% drop_na(Egg) %>% filter(Sex == 'F' & Analysis_PopulationGenetics == 1 & Egg != 'NA')
md %>% count(Egg)
#keep = md %>% count(Egg) %>% filter(n >= 4)
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
mds = md %>% mutate(LatJit = jitter(Latitude,amount =2),
                     LonJit = jitter(Longitude,amount=2))
sites = st_as_sf(mds, coords = c("LonJit", "LatJit"), 
                 crs = 4326, agr = "constant")

#set up map and convert df to coordinate frame
world = map_data("world")

imm_spat = ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), col='grey90', fill='white') +
  geom_sf(data = sites, 
          aes(fill=Egg,shape=Egg),
          size=3,show.legend = T) +
  xlab('')+ylab('')+
  coord_sf(xlim = c(min(mds$Longitude)-5, max(mds$Longitude)+5), 
           ylim = c(min(mds$Latitude)-5, max(mds$Latitude)+5), expand = FALSE)+
  scale_fill_manual(values=mds$EggCol,breaks=mds$Egg)+
  scale_shape_manual(values=mds$EggShape,breaks=mds$Egg)+
  theme_classic()+
  facet_grid(SpeciesShort ~.)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"))+
  theme(legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))+
  annotation_scale(line_width=0.5)
imm_spat

ggsave('~/symlinks/host/figures/20250801_EggHunt_Spatial-Females-AllvONE.pdf',imm_spat,
       dpi=300,height=6,width=8)
```

Targets (using conspecifics n=60 canorus and n=27 optatus as contrast for all FST analsyes).

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
# for CHR in $(cat Chromosomes.list); do sbatch -J Filter_${CHR} 00_Refilter_VCF.sh ${CHR}; done
CHR=$1

# Modify WD depending on canorus or optatus
WD=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/scan_n4/femalesN4_vs_all
raw_vcfs=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/merged/snps_only/

cd $WD

mkdir -p vcfs

#genotypes BELOW this will be set to missing
MINDP=3

if [[ $CHR = 'chr_W' || $CHR = 'chr_MT' ]]; then
        PLOIDY=1

        echo "FILTERING AND MERGING VARIANT SITES FOR ${CHR}, PLOIDY: ${PLOIDY}"
        bcftools view --threads 5 --samples-file AllEggs.list --force-samples -Ou ${raw_vcfs}/${CHR}.SNPS.vcf.gz | \
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
        bcftools view --threads 5 --samples-file AllEggs.list -Ou ${raw_vcfs}/${CHR}.SNPS.vcf.gz | \
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

This will annotate SNPs with synonymous/nonsynonymous using genome GFF. 

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

## bFST: 10-KB One v All [Fig 3]

Estimate background FST in 10-KB windows (100-bp for mtDNA). 

For each egg type, e.g. ECC1: 

* Identify number of samples for ECC1 and background (all canorus females) 
* Run 10 replicates: swap the population labels so that there are N=ECC1 samples randomly assigned ECC1
* Estimate FST in 10KB windows for the real comparison (FSTobs) and the n=10 random replicates
* bFST = FSTobs - maximum value observed in the 10 random replicates 

Subset:

```R
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/scan_n4')
.libPaths('~/symlinks/cuck00/software/mambaforge/envs/r25/lib/R/library')
library(tidyverse)

md = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt') 
md = md %>% drop_na(Egg) %>% filter(Sex == 'F' & Analysis_PopulationGenetics == 1 & Egg != 'NA')
md %>% count(Egg)
write.table(md$ID,file = 'femalesN4_vs_all/AllEggs.list',quote=F,sep='\t',row.names=F,col.names=F)
write.table(md %>% select(ID,Egg),file = 'femalesN4_vs_all/AllEggs.pop',quote=F,sep='\t',row.names=F,col.names=F)
```

#### Estimate FST

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --mem=15000mb
#SBATCH --cpus-per-task=2
#SBATCH --time=1-00:00:00


#mamba activate r25
# cat Egg_Types.list | xargs -I {} sbatch -J FST_{} 01_bFST.sh {}

if [ -z "$1" ]; then
    echo "Error: provide population ID"
    exit 1
fi

mkdir -p bgfst/work bgfst/out

TARGET=$1

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
    
        if [[ $CHR = 'chr_MT' ]]; then
            WIN_SIZE=100
        else
        	WIN_SIZE=10000
        fi

        echo "WORKING ON HAPLOID: True comparison, win size ${WIN_SIZE}"
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

		WIN_SIZE=10000
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

Ensure this script is available which calculates bFST (FSTobs - max(FSTrandomized)):

```R
#!/usr/bin/env rscript
args <- commandArgs(trailingOnly = TRUE)

setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/scan_n4/femalesN4_vs_all/bgfst/work')
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

#### Plot

This will:

* Examine saturation of chromosome class across eggs
* Plot genome scan
* Extract and plot candidate regions 

```R
#### Plot 10KB FST Gene Scan
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/scan_n4/femalesN4_vs_all')
.libPaths('~/r_libs')
library(karyoploteR)
library(ggpubr)
library(scales)
library(tidyverse)
library(viridis)
library(RColorBrewer)
library(zoo)
library(data.table)

md = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt')

##### Load and prep FST   ######
tidy_bp = fread('20250731_bFST.txt.gz')
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
eggcols <- md %>% arrange(EggOrder) %>% select(Egg,EggCol) %>% na.omit %>% unique %>% filter(Egg %in% unique(tidy_bp$Egg))
tidy_bp$Egg <- factor(tidy_bp$Egg, levels=eggcols$Egg)
tidy_bp$AvZ <- factor(tidy_bp$AvZ,levels=c('Autosome','Z','MT','W'))


###### WINDOW SATURATION ######
# calculate top 1% observed... and then see which windows fall above that 
thresh <- tidy_bp %>% 
  group_by(AvZ) %>% 
  summarize(threshold = quantile(bFST, 0.999, na.rm = TRUE))
td <- left_join(tidy_bp,thresh) 

lab <- thresh %>% ungroup  %>% 
  mutate(lab = paste0(AvZ,': ',round(threshold,2))) %>% 
  summarize(comb = paste0(lab,collapse='\n'),AvZ='W',Egg='ECC1')

wins <- td %>% 
  group_by(Egg) %>% 
  filter(bFST >= threshold) %>% 
  count(Egg,AvZ) %>% 
  ungroup() %>% 
  complete(Egg, AvZ, fill = list(n = 0)) %>% 
  ggplot(aes(x = Egg, fill = AvZ, y=n, label = n)) +
  geom_bar(col='black',stat = "identity", aes(y = n), position = position_dodge(width = 0.9)) + 
  geom_text(position = position_dodge(width = 0.9),size=2.75,vjust=-.25)+
  scale_fill_manual(values = viridis(4,option='mako')) +
  theme_bw(base_size = 10) +
  geom_text(data=lab,aes(x=Egg,y=Inf,label=lab$comb),vjust=1,size=2)+
  coord_cartesian(ylim=c(0,450))+
  ylab('') +
  xlab('') +
  theme(legend.position = "right", 
        legend.text = element_text(size = 6), 
        legend.title = element_text(size = 6)) + 
  guides(fill = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 1)))
wins
ggsave('~/symlinks/host/figures/20250731_bFSTWindowDistributions-AvZ-Thresholds999-AllvsONE.pdf',
       wins,height=1.5,width=7.5,dpi=300)

# Alternative: Permutation tests
set.seed(123)  
n_perm <- 1000

### Permutation: within AvZ levels, is there excess of one egg?
permdat <- list()
for (species in c('CC','CO')) { 
  
  cat('Working on: ',species,'\n')
  sp_bp <- tidy_bp %>% filter(Species == species)
  obs <- sp_bp %>%
    group_by(Egg) %>%
    mutate(thresh = quantile(bFST, 0.999, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(is_outlier = bFST >= thresh) %>%
    group_by(AvZ, Egg) %>%
    summarise(observed = sum(is_outlier), .groups = "drop")
  
  ### Across species 
  perm_res <- map_dfr(1:n_perm, function(i) {
    sp_bp %>%
      group_by(Egg) %>%
      mutate(thresh = quantile(bFST, 0.999, na.rm = TRUE)) %>%
      ungroup() %>%
      mutate(is_outlier = bFST >= thresh) %>%
      group_by(AvZ) %>%
      mutate(Egg_perm = sample(Egg)) %>%
      group_by(AvZ, Egg_perm) %>%
      summarise(n_perm = sum(is_outlier), .groups = "drop") %>%
      rename(Egg = Egg_perm) %>%
      mutate(iter = i)
  })
  
  # calculate empirical p-values: is each egg enriched for outliers within each AvZ?
  perm_summary <- perm_res %>%
    group_by(AvZ, Egg) %>%
    summarise(
      mean_perm = mean(n_perm),
      .groups = "drop"
    ) %>%
    left_join(obs, by = c("AvZ", "Egg")) %>%
    rowwise() %>%
    mutate(
      p_empirical = mean(
        abs(perm_res$n_perm[perm_res$AvZ == AvZ & perm_res$Egg == Egg] - mean_perm) >= abs(observed - mean_perm)
      )
    ) %>%
    ungroup()
  
  perm_avz <- perm_summary %>%
    mutate(
      enrichment = observed / mean_perm,
      sig_label = if_else(p_empirical < 0.05, "*", ""),
      Species = species
    )
  
  permdat[[species]] <- perm_avz
}

# Bind 
perms <- rbindlist(permdat) %>% as_tibble

# merge observed and permuted
pe <- perms %>% 
  ggplot(aes(x = Egg, y = enrichment, fill = AvZ)) +
  geom_col(position = position_dodge(width = 0.8), color = "black") +
  geom_text(aes(label = sig_label),position = position_dodge(width = 0.8),vjust = -0.5,size = 4) +
  coord_cartesian(ylim=c(0,max(perm_avz$enrichment)*1.1))+geom_hline(yintercept=1,lty=2)+
  scale_fill_viridis_d(option = "mako") +
  theme_bw(base_size = 10) +
  ylab('')+xlab('') 
pe
ggsave('~/symlinks/host/figures/20250807_bFSTPermutations-ByEggBySpecies-999_AvZ-AllvONE.pdf',
       pe,height=1.5,width=7.5,dpi=300)

# Third Alternative: lowest observed in top 0.1% by species 
outlier_low <- tidy_bp %>% 
  group_by(Species,AvZ) %>% 
  filter(bFST >= quantile(bFST, 0.999, na.rm = TRUE)) %>% 
  summarize(threshold = min(bFST))
tdo <- left_join(tidy_bp,outlier_low) 

lab <- outlier_low %>% ungroup  %>% 
  mutate(lab = paste0(AvZ,': ',round(threshold,2))) %>% 
  group_by(Species) %>% 
  summarize(comb = paste0(lab,collapse='\n'),AvZ='W') %>% 
  mutate(Egg = c('ECC1','ECO1'))
lab
lowin <- tdo %>% 
  group_by(Egg) %>% 
  filter(bFST >= threshold) %>% 
  count(Egg,AvZ) %>% 
  ungroup() %>% 
  complete(Egg, AvZ, fill = list(n = 0)) %>% 
  ggplot(aes(x = Egg, fill = AvZ, y=n, label = n)) +
  geom_bar(col='black',stat = "identity", aes(y = n), position = position_dodge(width = 0.9)) + 
  geom_text(position = position_dodge(width = 0.9),size=3,vjust=-.25)+
  scale_fill_manual(values = viridis(4,option='mako')) +
  geom_text(data=lab,aes(x=Egg,y=Inf,label=lab$comb),vjust=2,size=2)+
  theme_bw(base_size = 10) +
  coord_cartesian(ylim=c(0,325))+
  ylab('') +
  xlab('') +
  theme(legend.position = "right", 
        legend.text = element_text(size = 6), 
        legend.title = element_text(size = 6)) + 
  guides(fill = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 1)))
lowin
ggsave('~/symlinks/host/figures/20250731_bFSTWindowDistributions-OutliersLowestFST-AvZ-Thresholds999-AllvONE.pdf',
       lowin,height=1.5,width=7.5,dpi=300)


###### GWAS #######

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

pdf('~/symlinks/host/figures/20250731_KARYOPLOT-bFST10KB-AllvONE.-labels.pdf',height=3.5,width=7.5)
png('~/symlinks/host/figures/20250731_KARYOPLOT-bFST10KB-AllvONE.png',height=3.5,width=7.5,units='in',res=300)
pp = getDefaultPlotParams(plot.type=4)
pp$leftmargin = 0.1
kp = plotKaryotype(plot.type=4, 
                   genome = G,
                   labels.plotter = NULL,
                   plot.params = pp)
kpAddChromosomeNames(kp, yoffset = -5,cex=0.7)
kpAddBaseNumbers(kp,tick.dist = 25000000,minor.ticks = FALSE,cex=0.4)

#Loop through each track (all species together, CC, and CO), and plot them on their own layer 
counter = 0; tracks = 6
for (egg in rev(eggcols$Egg)) { 
  
  counter = counter + 1; at = autotrack(current.track = counter, total.tracks = tracks);at$r1 <- at$r1-0.02
  
  # Grab target egg 
  cat('Plotting: ',egg,'\n')
  sub <- karyoA %>% filter(Egg == egg) 
  ymin=0;ymax=1
  kpPoints(kp,chr=sub$chr,x=sub$start+5000,y=sub$bFST,ymin=ymin,ymax=ymax,r0=at$r0,r1=at$r1,col=as.character(sub$color))
  kpAxis(kp,r0=at$r0,r1=at$r1,cex=0.5,numticks = 2,ymin=ymin,ymax=ymax)
  kpAddLabels(kp,cex=0.5,labels = egg,r0=at$r0+.01, r1=at$r1,col="black",srt=0,label.margin = 0.02)
  
}

dev.off()

#### EXAMINE REGIONS OF INTEREST ####
leg <- md %>% select(Egg,EggCol) %>% na.omit %>% distinct
genes <- read_tsv('Gene_Lookup.sorted.bed',col_names = F) %>% mutate(X1 = gsub('chr_','',X1))
names(genes) <- c('chr','start','end','strand','gene')
top_hit <- tidy_bp %>%
  group_by(Egg,AvZ) %>%  # Only group by Egg, not chr
  filter(grepl('MT',chr)) %>% 
  slice_max(order_by = bFST, n = 1, with_ties = TRUE) %>%  # Get top hits per egg 
  ungroup() %>%
  select(chr, start, end, Egg, bFST) %>%
  data.frame() %>% filter(grepl('ECC',Egg))
tg <- genes %>% filter(chr == 'MT') %>% mutate(gene = gsub('ID=','',gene))
tg <- tg %>%
  mutate(xstart = if_else(strand == "+", start, end),xend   = if_else(strand == "+", end, start))
ceiling <- tidy_bp %>% filter(chr == 'MT') %>% summarize(ymax=max(bFST)*.99) %>% pull(ymax)
mtp <- tidy_bp %>% filter(chr == 'MT') %>%
  ggplot(aes(x=start,y=bFST,col=Egg))+
  geom_segment(data = tg,aes(x = xstart, xend = xend, y = ceiling*0.95, yend = ceiling*0.95),
               arrow = arrow(length = unit(0.07, "cm"),ends = "last",type = "closed"),inherit.aes = FALSE,linewidth = 0.3,color = "gray40") +
  geom_text(data=tg,aes(x=start+((end-start)/2),y=ceiling,label=gene),angle=30, size=1, inherit.aes=FALSE,hjust = 0,vjust=1)+
  geom_text(data=top_hit,aes(x=start+50),y=Inf,label='*',col='blue',alpha=1,size=3,vjust=1)+
  scale_color_manual(values=leg$EggCol,breaks=leg$Egg)+
  geom_line()+scale_x_continuous(breaks = pretty_breaks(n = 3),labels = function(x) paste0(round(x / 1e3, 1), "-Kb")) +scale_y_continuous(breaks=pretty_breaks(n=3))+
  theme_classic(base_size=8)+ylab('')+xlab('')+coord_cartesian(clip='off')
mtp

# First, inspect mtDNA 
outliers <- tidy_bp %>%
  group_by(Egg,AvZ) %>%  # Only group by Egg, not chr
  filter(!grepl('scaf|MT|W',chr)) %>% 
  slice_max(order_by = bFST, n = 1, with_ties = TRUE) %>%  # Get top hits per egg 
  ungroup() %>%
  select(chr, start, end, Egg, bFST) %>%
  data.frame() %>% 
  mutate(chr = gsub('chr_scaffold','scaffold',paste0('chr_',chr)))
outliers

##### ECC1/ECC6 #####
# ECC6 chr2 : BLVRA
deets <- tidy_bp %>% filter(Egg == 'ECC6' & !grepl('scaf',chr)) %>% group_by(chr) %>% slice_max(bFST) %>% arrange(desc(bFST)) %>% head(3) %>% select(chr,start,end)
deets
zchr=2
zstart=104840001;zs=zstart-5e5;ze=zstart+5e5
tg <- genes %>% filter(chr == zchr & start > zs & end < ze) %>%  mutate(xstart = if_else(strand == "+", start, end),xend   = if_else(strand == "+", end, start))
ceiling <- tidy_bp %>% filter(chr == zchr & start > zs & end < ze) %>% summarize(ymax=max(bFST)*.99) %>% pull(ymax)
p1 <- tidy_bp %>% filter(chr == zchr & start > zs & end < ze) %>%
  ggplot(aes(x=start,y=bFST,col=Egg))+
  geom_segment(data = tg,aes(x = xstart, xend = xend, y = ceiling*0.95, yend = ceiling*0.95),
               arrow = arrow(length = unit(0.07, "cm"),ends = "last",type = "closed"),inherit.aes = FALSE,linewidth = 0.3,color = "gray40") +
  geom_text(data=tg,aes(x=start+((end-start)/2),y=ceiling,label=gene),angle=30, size=1, inherit.aes=FALSE,hjust = 0,vjust=1)+
  geom_text(x=zstart+5000,y=Inf,label='*',col='blue',alpha=1,size=3,vjust=1)+
  scale_color_manual(values=leg$EggCol,breaks=leg$Egg)+
  ggtitle(paste0(zchr,' ',zstart,'ECC6:BLVRA,VOPP1'))+
  geom_line()+scale_x_continuous(breaks = pretty_breaks(n = 3),labels = function(x) paste0(round(x / 1e6, 1), "-Mb")) +scale_y_continuous(breaks=pretty_breaks(n=3))+
  theme_classic(base_size=8)+ylab('')+xlab('')+coord_cartesian(clip='off')
p1

# ECC6/ECO4 chrZ : PTPRD
deets <- tidy_bp %>% filter(Egg == 'ECC6' & !grepl('scaf',chr)) %>% group_by(chr) %>% slice_max(bFST) %>% arrange(desc(bFST)) %>% head(10) 
deets
zchr='11'
zstart=8820001;zs=zstart-1e5;ze=zstart+1e5
tg <- genes %>% filter(chr == zchr & start > zs & end < ze) %>%  mutate(xstart = if_else(strand == "+", start, end),xend   = if_else(strand == "+", end, start))
ceiling <- tidy_bp %>% filter(chr == zchr & start > zs & end < ze) %>% summarize(ymax=max(bFST)*.99) %>% pull(ymax)
ecc6a <- tidy_bp %>% filter(chr == zchr & start > zs & end < ze) %>%
  ggplot(aes(x=start,y=bFST,col=Egg))+
  geom_segment(data = tg,aes(x = xstart, xend = xend, y = ceiling*0.95, yend = ceiling*0.95),
               arrow = arrow(length = unit(0.07, "cm"),ends = "last",type = "closed"),inherit.aes = FALSE,linewidth = 0.3,color = "gray40") +
  geom_text(data=tg,aes(x=start+((end-start)/2),y=ceiling,label=gene),angle=30, size=3, inherit.aes=FALSE,hjust = 0,vjust=1)+
  geom_text(x=zstart+5000,y=Inf,label='*',col='blue',alpha=1,size=3,vjust=1)+
  geom_text(x=21370001+5000,y=Inf,label='*',col='blue',alpha=1,size=3,vjust=1)+ # ECO4
  scale_color_manual(values=leg$EggCol,breaks=leg$Egg)+
  ggtitle(paste0(zchr,' ',zstart,'ECC6/ECC1:MTMR14'))+
  geom_line()+scale_x_continuous(breaks = pretty_breaks(n = 3),labels = function(x) paste0(round(x / 1e6, 1), "-Mb")) +scale_y_continuous(breaks=pretty_breaks(n=3))+
  theme_classic(base_size=8)+ylab('')+xlab('')+coord_cartesian(clip='off')
ecc6a

# alt 
zchr='Z'
zstart=20730001;zs=zstart-8e5;ze=zstart+11e5
tg <- genes %>% filter(chr == zchr & start > zs & end < ze) %>%  mutate(xstart = if_else(strand == "+", start, end),xend   = if_else(strand == "+", end, start))
ceiling <- tidy_bp %>% filter(chr == zchr & start > zs & end < ze) %>% summarize(ymax=max(bFST)*.99) %>% pull(ymax)
p2 <- tidy_bp %>% filter(chr == zchr & start > zs & end < ze) %>%
  ggplot(aes(x=start,y=bFST,col=Egg))+
  geom_segment(data = tg,aes(x = xstart, xend = xend, y = ceiling*0.95, yend = ceiling*0.95),
               arrow = arrow(length = unit(0.07, "cm"),ends = "last",type = "closed"),inherit.aes = FALSE,linewidth = 0.3,color = "gray40") +
  geom_text(data=tg,aes(x=start+((end-start)/2),y=ceiling,label=gene),angle=30, size=1, inherit.aes=FALSE,hjust = 0,vjust=1)+
  geom_text(x=zstart+5000,y=Inf,label='*',col='blue',alpha=1,size=3,vjust=1)+
  geom_text(x=21370001+5000,y=Inf,label='*',col='blue',alpha=1,size=3,vjust=1)+ # ECO4
  scale_color_manual(values=leg$EggCol,breaks=leg$Egg)+
  ggtitle(paste0(zchr,' ',zstart,'ECC6/ECO4:PTPRD'))+
  geom_line()+scale_x_continuous(breaks = pretty_breaks(n = 3),labels = function(x) paste0(round(x / 1e6, 1), "-Mb")) +scale_y_continuous(breaks=pretty_breaks(n=3))+
  theme_classic(base_size=8)+ylab('')+xlab('')+coord_cartesian(clip='off')
p2

##### ECC10 #####
# ECC10 chrZ : RLN3
deets <- tidy_bp %>% filter(Egg == 'ECC10' & !grepl('scaf',chr)) %>% group_by(chr) %>% slice_max(bFST) %>% arrange(desc(bFST)) %>% head(3) %>% select(chr,start,end)
deets
zchr='Z'
zstart=18980001;zs=zstart-2e5;ze=zstart+2.5e5
tg <- genes %>% filter(chr == zchr & start > zs & end < ze) %>%  mutate(xstart = if_else(strand == "+", start, end),xend   = if_else(strand == "+", end, start))
ceiling <- tidy_bp %>% filter(chr == zchr & start > zs & end < ze) %>% summarize(ymax=max(bFST)*.99) %>% pull(ymax)
p3 <- tidy_bp %>% filter(chr == zchr & start > zs & end < ze) %>%
  ggplot(aes(x=start,y=bFST,col=Egg))+
  geom_segment(data = tg,aes(x = xstart, xend = xend, y = ceiling*0.95, yend = ceiling*0.95),
               arrow = arrow(length = unit(0.07, "cm"),ends = "last",type = "closed"),inherit.aes = FALSE,linewidth = 0.3,color = "gray40") +
  geom_text(data=tg,aes(x=start+((end-start)/2),y=ceiling,label=gene),angle=30, size=1, inherit.aes=FALSE,hjust = 0,vjust=1)+
  geom_text(x=zstart+5000,y=ceiling,label='*',col='blue',alpha=1,size=3)+
  scale_color_manual(values=leg$EggCol,breaks=leg$Egg)+
  ggtitle(paste0(zchr,' ',zstart,'ECC10:RLN3'))+
  geom_line()+scale_x_continuous(breaks = pretty_breaks(n = 3),labels = function(x) paste0(round(x / 1e6, 1), "-Mb")) +scale_y_continuous(breaks=pretty_breaks(n=3))+
  theme_classic(base_size=8)+ylab('')+xlab('')+coord_cartesian(clip='off')
p3

# ECC10 inspect other autosomal hits
deets <- tidy_bp %>% filter(Egg == 'ECC10' & !grepl('scaf',chr)) %>% group_by(chr) %>% slice_max(bFST) %>% arrange(desc(bFST)) %>% head(10) 
deets
zchr='30'
zstart=350001;zs=zstart-0.5e5;ze=zstart+0.5e5
tg <- genes %>% filter(chr == zchr & start > zs & end < ze) %>%  mutate(xstart = if_else(strand == "+", start, end),xend   = if_else(strand == "+", end, start))
ceiling <- tidy_bp %>% filter(chr == zchr & start > zs & end < ze) %>% summarize(ymax=max(bFST)*.99) %>% pull(ymax)
ecc10a <- tidy_bp %>% filter(chr == zchr & start > zs & end < ze) %>%
  ggplot(aes(x=start,y=bFST,col=Egg))+
  geom_segment(data = tg,aes(x = xstart, xend = xend, y = ceiling*0.95, yend = ceiling*0.95),
               arrow = arrow(length = unit(0.07, "cm"),ends = "last",type = "closed"),inherit.aes = FALSE,linewidth = 0.3,color = "gray40") +
  geom_text(data=tg,aes(x=start+((end-start)/2),y=ceiling,label=gene),angle=30, size=1, inherit.aes=FALSE,hjust = 0,vjust=1)+
  geom_text(x=zstart+5000,y=ceiling,label='*',col='blue',alpha=1,size=3)+
  scale_color_manual(values=leg$EggCol,breaks=leg$Egg)+
  ggtitle(paste0(zchr,' ',zstart,'ECC10:EIF3G'))+
  geom_line()+scale_x_continuous(breaks = pretty_breaks(n = 3),labels = function(x) paste0(round(x / 1e6, 1), "-Mb")) +scale_y_continuous(breaks=pretty_breaks(n=3))+
  theme_classic(base_size=8)+ylab('')+xlab('')+coord_cartesian(clip='off')
ecc10a 
zchr='31'
zstart=1130001;zs=zstart-1e5;ze=zstart+1e5
tg <- genes %>% filter(chr == zchr & start > zs & end < ze) %>%  mutate(xstart = if_else(strand == "+", start, end),xend   = if_else(strand == "+", end, start))
ceiling <- tidy_bp %>% filter(chr == zchr & start > zs & end < ze) %>% summarize(ymax=max(bFST)*.99) %>% pull(ymax)
ecc10b <- tidy_bp %>% filter(chr == zchr & start > zs & end < ze) %>%
  ggplot(aes(x=start,y=bFST,col=Egg))+
  geom_segment(data = tg,aes(x = xstart, xend = xend, y = ceiling*0.95, yend = ceiling*0.95),
               arrow = arrow(length = unit(0.07, "cm"),ends = "last",type = "closed"),inherit.aes = FALSE,linewidth = 0.3,color = "gray40") +
  geom_text(data=tg,aes(x=start+((end-start)/2),y=ceiling,label=gene),angle=30, size=1, inherit.aes=FALSE,hjust = 0,vjust=1)+
  geom_text(x=zstart+5000,y=ceiling,label='*',col='blue',alpha=1,size=3)+
  scale_color_manual(values=leg$EggCol,breaks=leg$Egg)+
  ggtitle(paste0(zchr,' ',zstart,'ECC10:COL11A2'))+
  geom_line()+scale_x_continuous(breaks = pretty_breaks(n = 3),labels = function(x) paste0(round(x / 1e6, 1), "-Mb")) +scale_y_continuous(breaks=pretty_breaks(n=3))+
  theme_classic(base_size=8)+ylab('')+xlab('')+coord_cartesian(clip='off')
ecc10b
zchr='39'
zstart=760001;zs=zstart-0.5e5;ze=zstart+0.5e5
tg <- genes %>% filter(chr == zchr & start > zs & end < ze) %>%  mutate(xstart = if_else(strand == "+", start, end),xend   = if_else(strand == "+", end, start))
ceiling <- tidy_bp %>% filter(chr == zchr & start > zs & end < ze) %>% summarize(ymax=max(bFST)*.99) %>% pull(ymax)
ecc10c <- tidy_bp %>% filter(chr == zchr & start > zs & end < ze) %>%
  ggplot(aes(x=start,y=bFST,col=Egg))+
  geom_segment(data = tg,aes(x = xstart, xend = xend, y = ceiling*0.95, yend = ceiling*0.95),
               arrow = arrow(length = unit(0.07, "cm"),ends = "last",type = "closed"),inherit.aes = FALSE,linewidth = 0.3,color = "gray40") +
  geom_text(data=tg,aes(x=start+((end-start)/2),y=ceiling,label=gene),angle=30, size=1, inherit.aes=FALSE,hjust = 0,vjust=1)+
  geom_text(x=zstart+5000,y=ceiling,label='*',col='blue',alpha=1,size=3)+
  scale_color_manual(values=leg$EggCol,breaks=leg$Egg)+
  ggtitle(paste0(zchr,' ',zstart,'ECC10:OXA1L,RNF31'))+
  geom_line()+scale_x_continuous(breaks = pretty_breaks(n = 3),labels = function(x) paste0(round(x / 1e6, 1), "-Mb")) +scale_y_continuous(breaks=pretty_breaks(n=3))+
  theme_classic(base_size=8)+ylab('')+xlab('')+coord_cartesian(clip='off')
ecc10c

##### ECO1 #####
# ECO1 chr6 : BOLL
deets <- tidy_bp %>% filter(Egg == 'ECO1' & !grepl('scaf',chr)) %>% group_by(chr) %>% slice_max(bFST) %>% arrange(desc(bFST)) %>% head(10) 
deets
zchr=6
zstart=31140001;zs=zstart-1.5e6;ze=zstart+5e5
tg <- genes %>% filter(chr == zchr & start > zs & end < ze) %>%  mutate(xstart = if_else(strand == "+", start, end),xend   = if_else(strand == "+", end, start))
ceiling <- tidy_bp %>% filter(chr == zchr & start > zs & end < ze) %>% summarize(ymax=max(bFST)*.99) %>% pull(ymax)
p4 <- tidy_bp %>% filter(chr == zchr & start > zs & end < ze) %>%
  ggplot(aes(x=start,y=bFST,col=Egg))+
  geom_segment(data = tg,aes(x = xstart, xend = xend, y = ceiling*0.95, yend = ceiling*0.95),
               arrow = arrow(length = unit(0.07, "cm"),ends = "last",type = "closed"),inherit.aes = FALSE,linewidth = 0.3,color = "gray40") +
  geom_text(data=tg,aes(x=start+((end-start)/2),y=ceiling,label=gene),angle=30, size=1, inherit.aes=FALSE,hjust = 0,vjust=1)+
  geom_text(x=zstart+5000,y=Inf,label='*',col='blue',alpha=1,size=3,vjust=1)+
  scale_color_manual(values=leg$EggCol,breaks=leg$Egg)+
  ggtitle(paste0(zchr,' ',zstart,'ECO1:BOLL'))+
  geom_line()+scale_x_continuous(breaks = pretty_breaks(n = 3),labels = function(x) paste0(round(x / 1e6, 1), "-Mb")) +scale_y_continuous(breaks=pretty_breaks(n=3))+
  theme_classic(base_size=8)+ylab('')+xlab('')+coord_cartesian(clip='off')
p4

#ECo1, other peaks
zchr=30
zstart=920001;zs=zstart-1e5;ze=zstart+1e5
tg <- genes %>% filter(chr == zchr & start > zs & end < ze) %>%  mutate(xstart = if_else(strand == "+", start, end),xend   = if_else(strand == "+", end, start))
ceiling <- tidy_bp %>% filter(chr == zchr & start > zs & end < ze) %>% summarize(ymax=max(bFST)*.99) %>% pull(ymax)
eco1a <- tidy_bp %>% filter(chr == zchr & start > zs & end < ze) %>%
  ggplot(aes(x=start,y=bFST,col=Egg))+
  geom_segment(data = tg,aes(x = xstart, xend = xend, y = ceiling*0.95, yend = ceiling*0.95),
               arrow = arrow(length = unit(0.07, "cm"),ends = "last",type = "closed"),inherit.aes = FALSE,linewidth = 0.3,color = "gray40") +
  geom_text(data=tg,aes(x=start+((end-start)/2),y=ceiling,label=gene),angle=30, size=1, inherit.aes=FALSE,hjust = 0,vjust=1)+
  geom_text(x=zstart+5000,y=Inf,label='*',col='blue',alpha=1,size=3,vjust=1)+
  scale_color_manual(values=leg$EggCol,breaks=leg$Egg)+
  ggtitle(paste0(zchr,' ',zstart,'ECO1:30:FDX2,SLC27A1'))+
  geom_line()+scale_x_continuous(breaks = pretty_breaks(n = 3),labels = function(x) paste0(round(x / 1e6, 1), "-Mb")) +scale_y_continuous(breaks=pretty_breaks(n=3))+
  theme_classic(base_size=8)+ylab('')+xlab('')+coord_cartesian(clip='off')
eco1a
zchr='Z'
zstart=45090001;zs=zstart-2e5;ze=zstart+3e5
tg <- genes %>% filter(chr == zchr & start > zs & end < ze) %>%  mutate(xstart = if_else(strand == "+", start, end),xend   = if_else(strand == "+", end, start))
ceiling <- tidy_bp %>% filter(chr == zchr & start > zs & end < ze) %>% summarize(ymax=max(bFST)*.99) %>% pull(ymax)
eco1b <- tidy_bp %>% filter(chr == zchr & start > zs & end < ze) %>%
  ggplot(aes(x=start,y=bFST,col=Egg))+
  geom_segment(data = tg,aes(x = xstart, xend = xend, y = ceiling*0.95, yend = ceiling*0.95),
               arrow = arrow(length = unit(0.07, "cm"),ends = "last",type = "closed"),inherit.aes = FALSE,linewidth = 0.3,color = "gray40") +
  geom_text(data=tg,aes(x=start+((end-start)/2),y=ceiling,label=gene),angle=30, size=1, inherit.aes=FALSE,hjust = 0,vjust=1)+
  geom_text(x=zstart+5000,y=Inf,label='*',col='blue',alpha=1,size=3,vjust=1)+
  scale_color_manual(values=leg$EggCol,breaks=leg$Egg)+
  ggtitle(paste0(zchr,' ',zstart,'ECO1:CPLX1'))+
  geom_line()+scale_x_continuous(breaks = pretty_breaks(n = 3),labels = function(x) paste0(round(x / 1e6, 1), "-Mb")) +scale_y_continuous(breaks=pretty_breaks(n=3))+
  theme_classic(base_size=8)+ylab('')+xlab('')+coord_cartesian(clip='off')
eco1b

##### ECO3 #####
# ECO3 chrZ : CHSY3 & HINT1 
deets <- tidy_bp %>% filter(Egg == 'ECO3' & !grepl('scaf',chr)) %>% group_by(chr) %>% slice_max(bFST) %>% arrange(desc(bFST)) %>% head(5) %>% select(chr,start,end)
deets
zchr='Z'
zstart=34100001;zs=zstart-2.5e5;ze=zstart+6e5
tg <- genes %>% filter(chr == zchr & start > zs & end < ze) %>%  mutate(xstart = if_else(strand == "+", start, end),xend   = if_else(strand == "+", end, start))
ceiling <- tidy_bp %>% filter(chr == zchr & start > zs & end < ze) %>% summarize(ymax=max(bFST)*.99) %>% pull(ymax)
p5 <- tidy_bp %>% filter(chr == zchr & start > zs & end < ze) %>%
  ggplot(aes(x=start,y=bFST,col=Egg))+
  geom_segment(data = tg,aes(x = xstart, xend = xend, y = ceiling*0.95, yend = ceiling*0.95),
               arrow = arrow(length = unit(0.07, "cm"),ends = "last",type = "closed"),inherit.aes = FALSE,linewidth = 0.3,color = "gray40") +
  geom_text(data=tg,aes(x=start+((end-start)/2),y=ceiling,label=gene),angle=30, size=1, inherit.aes=FALSE,hjust = 0,vjust=1)+
  scale_color_manual(values=leg$EggCol,breaks=leg$Egg)+
  geom_line()+scale_x_continuous(breaks = pretty_breaks(n = 3),labels = function(x) paste0(round(x / 1e6, 1), "-Mb")) +scale_y_continuous(breaks=pretty_breaks(n=3))+
  ggtitle(paste0(zchr,' ',zstart,'ECO3:CHSY3 & HINT1'))+
  geom_text(x=zstart+5000,y=Inf,label='*',col='blue',alpha=1,size=3,vjust=1)+
  theme_classic(base_size=8)+ylab('')+xlab('')+coord_cartesian(clip='off')
p5

##### ECO4 #####
# ECO4 chr 5 TSHR
deets <- tidy_bp %>% filter(Egg == 'ECO4' & !grepl('scaf',chr)) %>% group_by(chr) %>% slice_max(bFST) %>% arrange(desc(bFST)) %>% head(5) 
deets
zchr=5
zstart=19500001;zs=zstart-5e5;ze=zstart+5e5
tg <- genes %>% filter(chr == zchr & start > zs & end < ze) %>%  mutate(xstart = if_else(strand == "+", start, end),xend   = if_else(strand == "+", end, start))
ceiling <- tidy_bp %>% filter(chr == zchr & start > zs & end < ze) %>% summarize(ymax=max(bFST)*.99) %>% pull(ymax)
p6 <- tidy_bp %>% filter(chr == zchr & start > zs & end < ze) %>%
  ggplot(aes(x=start,y=bFST,col=Egg))+
  geom_segment(data = tg,aes(x = xstart, xend = xend, y = ceiling*0.95, yend = ceiling*0.95),
               arrow = arrow(length = unit(0.07, "cm"),ends = "last",type = "closed"),inherit.aes = FALSE,linewidth = 0.3,color = "gray40") +
  geom_text(data=tg,aes(x=start+((end-start)/2),y=ceiling,label=gene),angle=30, size=1, inherit.aes=FALSE,hjust = 0,vjust=1)+
  scale_color_manual(values=leg$EggCol,breaks=leg$Egg)+
  ggtitle(paste0(zchr,' ',zstart,'ECO4:TSHR'))+
  geom_text(x=zstart+5000,y=Inf,label='*',col='blue',alpha=1,size=3,vjust=1)+
  geom_line()+scale_x_continuous(breaks = pretty_breaks(n = 3),labels = function(x) paste0(round(x / 1e6, 1), "-Mb")) +scale_y_continuous(breaks=pretty_breaks(n=3))+
  theme_classic(base_size=8)+ylab('')+xlab('')+coord_cartesian(clip='off')
p6

# ECO4 chr28 HCN3
deets
zchr=28
zstart=3940001;zs=zstart-3.5e6;ze=zstart+8e6
tg <- genes %>% filter(chr == zchr & start > zstart-3.5e6 & end < zstart+8e6) %>%  mutate(xstart = if_else(strand == "+", start, end),xend   = if_else(strand == "+", end, start))
ceiling <- tidy_bp %>% filter(chr == zchr & start > zs & end < ze) %>% summarize(ymax=max(bFST)*.99) %>% pull(ymax)
p7 <- tidy_bp %>% filter(chr == zchr & start > zs & end < ze) %>%
  ggplot(aes(x=start,y=bFST,col=Egg))+
  geom_segment(data = tg,aes(x = xstart, xend = xend, y = ceiling*0.95, yend = ceiling*0.95),
               arrow = arrow(length = unit(0.07, "cm"),ends = "last",type = "closed"),inherit.aes = FALSE,linewidth = 0.3,color = "gray40") +
  geom_text(data=tg,aes(x=start+((end-start)/2),y=ceiling,label=gene),angle=30, size=1, inherit.aes=FALSE,hjust = 0,vjust=1)+
  scale_color_manual(values=leg$EggCol,breaks=leg$Egg)+
  ggtitle(paste0(zchr,' ',zstart,'ECO4:HCN3'))+
  geom_text(x=zstart+5000,y=Inf,label='*',col='blue',alpha=1,size=3,vjust=1)+
  geom_line()+scale_x_continuous(breaks = pretty_breaks(n = 3),labels = function(x) paste0(round(x / 1e6, 1), "-Mb")) +scale_y_continuous(breaks=pretty_breaks(n=3))+
  theme_classic(base_size=8)+ylab('')+xlab('')+coord_cartesian(clip='off')
p7

# alt 
zchr=11
zstart=180001;zs=zstart-1.5e5;ze=zstart+2.5e5
tg <- genes %>% filter(chr == zchr & start > zs & end < ze) %>%  mutate(xstart = if_else(strand == "+", start, end),xend   = if_else(strand == "+", end, start))
ceiling <- tidy_bp %>% filter(chr == zchr & start > zs & end < ze) %>% summarize(ymax=max(bFST)*.99) %>% pull(ymax)
eco4a <- tidy_bp %>% filter(chr == zchr & start > zs & end < ze) %>%
  ggplot(aes(x=start,y=bFST,col=Egg))+
  geom_segment(data = tg,aes(x = xstart, xend = xend, y = ceiling*0.95, yend = ceiling*0.95),
               arrow = arrow(length = unit(0.07, "cm"),ends = "last",type = "closed"),inherit.aes = FALSE,linewidth = 0.3,color = "gray40") +
  geom_text(data=tg,aes(x=start+((end-start)/2),y=ceiling,label=gene),angle=30, size=1, inherit.aes=FALSE,hjust = 0,vjust=1)+
  scale_color_manual(values=leg$EggCol,breaks=leg$Egg)+
  geom_line()+scale_x_continuous(breaks = pretty_breaks(n = 3),labels = function(x) paste0(round(x / 1e6, 1), "-Mb")) +scale_y_continuous(breaks=pretty_breaks(n=3))+
  ggtitle(paste0(zchr,' ',zstart,'ECO4:SLC6A11'))+
  geom_text(x=zstart+5000,y=Inf,label='*',col='blue',alpha=1,size=3,vjust=1)+
  theme_classic(base_size=8)+ylab('')+xlab('')+coord_cartesian(clip='off')
eco4a


pdf('~/symlinks/host/figures/20250801_Candidates-AllvONE.pdf',height=4.5,width=6)
ggarrange(mtp,p1,p2,p3,p4,p5,p6,p7,nrow=4,ncol=2,common.legend = TRUE)
dev.off()

pdf('~/symlinks/host/figures/20250807_Candidates-AllvONE-Alts.pdf',height=4.5,width=6)
ggarrange(ecc6a,ecc10a,ecc10b,ecc10c,eco1a,eco1b,eco4a,nrow=4,ncol=2,common.legend = TRUE)
dev.off()

```

#### Investigate SNP frequencies chr2

Investigate bp-level FST around the BLVRA locus. Create pie charts showing the geographic distribution of alleles for the top SNPs. 

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_highmem
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00

#mamba activate snps 

mkdir -p chr2fst

for TARGET in ECC1 ECC6; do  
SPECIES="CC"

vcftools --gzvcf chr_2.SNP.DP3.vcf.gz --out chr2fst/${TARGET}_chr_2 \
            --weir-fst-pop bgfst/work/${TARGET}.TRUE-T.list --weir-fst-pop bgfst/work/${TARGET}.TRUE-F.list --fst-window-size 1 --max-missing 0.1
        awk -v e=${TARGET} '{OFS="\t"}{print $1, $2, $3, $5, e}' chr2fst/${TARGET}_chr_2.windowed.weir.fst | sed '1d' > chr2fst/${TARGET}_chr_2.out
done 
cat chr2fst/*_chr_2.out > chr2fst_entire_chr.txt

Rscript -e "
df <- read.table('chr2fst_entire_chr.txt');
df\$V4 <- pmax(pmin(df\$V4, 1), 0);
by(df, df\$V5, function(subdf) {
  ci <- quantile(subdf\$V4, c(0.025, 0.975), na.rm = TRUE);
  cat(unique(subdf\$V5), '95% CI:', ci[1], '-', ci[2], '\n')
})
"
```

Plot pies:

```bash
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/scan_n4/femalesN4_vs_all/chr2fst/')
.libPaths('~/r_libs/')
library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(scatterpie)
library(VariantAnnotation)

md = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt')

# Read in top SNPs 
e <- read_tsv('chr2fst_snps.txt',col_names = F)
names(e) <- c('chr','start','end','fst','egg')

genes <- read_tsv('../Gene_Lookup.sorted.bed',col_names = F)
targs <- genes %>% filter(grepl('BLVRA|VOPP1|LOC104065962', X5))
names(targs) <- c('chr','start','end','strand','gene')

snps <- NULL
for (g in targs$gene) {
  cat('Working on ',g,'\n')
  gt <- targs %>% filter(gene == g)
  et <- e %>% filter(start > gt$start & start < gt$end) %>% 
    mutate(fst = pmax(0,pmin(1,fst)), gene = g, start = start - 1) %>% 
    group_by(egg) %>% 
    slice_max(fst,n=1,with_ties = FALSE) %>% dplyr::select(chr,start,end,fst,egg,gene)
  snps <- rbind(snps,et)
  
}
snps <- snps %>% mutate(id = paste0(chr,':',end))
write.table(snps %>% ungroup %>% select(chr,start,end),file='target_snps.bed',quote=F,sep='\t',row.names=F,col.names=F)

# intersect with vcf to get genos 
# bedtools intersect -header -a ../vcfs/chr_2.SNP.DP3.vcf.gz -b target_snps.bed | bcftools view -Oz -o chr2fst_keep.vcf.gz
vcf <- readVcf("chr2fst_keep.vcf.gz", "cuckoo") 
geno_mat <- as.data.frame(geno(vcf)$GT)
geno_mat <- geno_mat %>%
  rownames_to_column("SNP") %>%
  pivot_longer(-SNP, names_to = "ID", values_to = "GT") %>% 
  separate(SNP,into=c('chr','id'),sep = ':', remove=F) %>% 
  separate(id,into=c('end','allele'),sep = '_', remove=T)  
targ_snps <- geno_mat %>% 
  mutate(end = as.numeric(end)) %>% 
  filter(end %in% snps$end) %>% 
  mutate(allele0 = str_count(GT, "0"),
         allele1 = str_count(GT, "1"))
targ_snps %>% dplyr::count(SNP)
merged <- left_join(targ_snps, md, by = "ID") 

agg_data <- merged %>%
  group_by(SNP, GeographicGroup) %>%
  summarise(
    A0 = sum(allele0, na.rm = TRUE),
    A1 = sum(allele1, na.rm = TRUE),
    Latitude = mean(Latitude, na.rm = TRUE),
    Longitude = mean(Longitude, na.rm = TRUE),
    r = sqrt(A0 + A1) * 0.8,  
    .groups = "drop"
  ) %>% filter(grepl('GCC',GeographicGroup))

ag <- agg_data %>% 
  separate(SNP,into=c('chr','start'),sep = ':', remove=F) %>%
  separate(start,into=c('start','allele')) %>% 
  mutate(start = as.numeric(start)-1) %>% 
  left_join(.,snps %>% mutate(lab = paste0(gene,': ',egg)) %>% ungroup %>%  dplyr::select(chr,start,gene,lab))

plots <- ggplot() +
  scatterpie::geom_scatterpie(
    data = ag,
    aes(x = Longitude, y = Latitude, r = r),
    cols = c("A0", "A1"),
    color = NA
  ) +
  scale_fill_manual(values = c("A0" = "#e41a1c", "A1" = "#377eb8")) +
  facet_wrap(~ SNP+lab) +
  coord_fixed() +
  theme_minimal()
plots
ggsave('~/symlinks/host/figures/20250802_Spatial_Distribution_Outlier_SNPS_chr2_candidates-AllvONE.pdf', plots,height=5,width=12)

```

#### Investigate NDUFAF4 canorus

See if any canorus egg types show autosomal differentiation around the BLVRA locus compared to the 500kb region surrounding it. 

```bash
#mamba activate merothon 
WD=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/manyhost_hunt/males/ndufaf4
CD ${WD}
mkdir -p ndufaf4/work ndufaf4/out

grep 'NDUFAF4' Gene_Lookup.bed | \
	awk '{OFS="\t"}{print $1, $2-500000, $3+500000,$4, $5}' > ndufaf4/ndufaf4.bed
bedtools intersect -header -a vcfs/chr_3.SNP.DP3.vcf.gz -b ndufaf4/ndufaf4.bed | bcftools view -Oz -o ndufaf4/ndufaf4.vcf.gz

for TARGET in ECC1 ECC3 ECC6 ECC7 ECC8 ECC10; do

awk -v e=$TARGET 'awk $2 == e' ndufaf4/Eggs.pop | awk '{print $1}' > ndufaf4/work/${TARGET}.T.list
awk -v e=$TARGET 'awk $2 != e' ndufaf4/Eggs.pop | awk '{print $1}' > ndufaf4/work/${TARGET}.F.list

vcftools --gzvcf ndufaf4/ndufaf4.vcf.gz --out ndufaf4/work/${TARGET} \
            --weir-fst-pop ndufaf4/work/${TARGET}.T.list --weir-fst-pop ndufaf4/work/${TARGET}.F.list --fst-window-size 1 --max-missing 0.1
        awk -v e=${TARGET} '{OFS="\t"}{print $1, $2, $3, $5, e}' ndufaf4/work/${TARGET}.windowed.weir.fst | sed '1d' > ndufaf4/out/${TARGET}.out
done
```

in R:

```bash
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/manyhost_hunt/males/ndufaf4')
.libPaths('~/r_libs/')
library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(meRo) #devtools::install_github('merondun/meRo')
library(devtools)

md = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt')
leg <- md %>% select(Egg,EggCol) %>% distinct
n <- read_tsv('NDUFAF4.txt',col_names = F)
names(n) <- c('chr','start','end','fst','Egg')
ceiling <- n %>% summarize(ymax=max(fst)*.99) %>% pull(ymax)
n <- n %>%
  mutate(fst = pmax(0,pmin(1,fst)),
         targ = ifelse(end > 25511697 & end < 25515677, 'NDUFAF4','Background'))
n$Egg <- factor(n$Egg,levels=c('ECC1','ECC3','ECC6','ECC7','ECC8','ECC10'))
np <- n %>% group_by(Egg,targ) %>% 
  sum_stats(fst) %>% 
  ggplot(aes(x=Egg,y=mean,ymin=mean-sd,ymax=mean+sd,col=Egg,shape=targ))+
  geom_point(position=position_dodge(width=0.9),size=5)+
  geom_errorbar(position=position_dodge(width=0.9),width=0.5)+
  ylab('Mean +/- SD FST')+
  scale_color_manual(values=leg$EggCol,breaks=leg$Egg)+
  coord_flip()+
  theme_bw()
np
ggsave('~/symlinks/host/figures/20250802_FST_Near_NDUFAF4-chr3.pdf', np,height=5,width=5)

```

Inspect to see if there are any nonsynonymous mutations:

```R
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/scan_n4/femalesN4_vs_all')
.libPaths('~/r_libs')
library(tidyverse)
library(data.table)
library(VariantAnnotation)
library(GenomicFeatures)

txdb = makeTxDbFromGFF("/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR-chr_MT.gff")
fasta_seq <- readDNAStringSet("/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa")

#load in vcf
vcf = readVcf("vcfs/chr_3.SNP.DP3.vcf.gz", genome = "cuckoo")
chr = 'chr_3'

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

write.table(dat,file='tmp.txt',quote=F,sep='\t',row.names=F)
dat = read.table('tmp.txt',header=TRUE)
writedat = dat %>%
  as.data.frame %>%
  unique %>%
  group_by(seqnames, start, REF, varAllele, GENEID) %>%
  summarise(CONSEQUENCE = paste(unique(CONSEQUENCE), collapse = ";"), .groups = "drop") %>% ungroup %>%
  arrange(seqnames,start)
write.table(writedat,file='ndufaf4/ndufaf4.annotated.txt',quote=F,sep='\t',row.names=F)

pnps <- writedat %>% mutate(gene = gsub('gene-','',GENEID)) %>% 
  group_by(gene) %>% dplyr::count(CONSEQUENCE) %>% 
  filter(CONSEQUENCE == 'nonsynonymous' | CONSEQUENCE == 'synonymous') %>% 
  pivot_wider(names_from = CONSEQUENCE,values_from = n) %>% 
  mutate(pnps = nonsynonymous/synonymous)
t <- pnps %>% filter(grepl('NDUFAF4',gene))
p <- pnps %>% 
  ggplot(aes(x=pnps))+
  geom_histogram()+
  geom_vline(xintercept=t$pnps,lty=2,color='blue')+
  xlab('pN/pS')+ylab('Count (Genes on chr3)')+
  theme_bw()
ggsave('~/symlinks/host/figures/20250805_pNpS_chr3_NDUFAF4.pdf',p,height=2.5,width=4)
```

## bFST: Females Nestlings N=3

Repeat all the above analyses, except only requiring n=3 females for focal eggs. This includes ECC8. This analysis will also only use nestlings, which applies to some optatus egg types (ECO4) which were designated based on a sole egg type existing within a locality. 

```bash
# Egg gene scan, sample selection for females >=3 and nestlings 
#### Find associations between MT/W haplotypes and features 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/scan_n4/females_nestlingsN3_vs_all')
.libPaths('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/software/mambaforge/envs/r25/lib/R/library')
library(tidyverse)
library(viridis)
library(tidyverse)
library(sf)
library(ggspatial)

md = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt') 
md = md %>% drop_na(Egg) %>% filter(Sex == 'F' & Analysis_PopulationGenetics == 1 & Egg != 'NA' & grepl('Young',Age))
md %>% count(Egg)
keep = md %>% count(Egg) %>% filter(n >= 3)
mds = md %>% filter(Egg %in% keep$Egg)
mds %>% count(Egg)

write.table(mds %>% select(ID,Egg) %>% arrange(Egg),file='Eggs.pop',quote=F,sep='\t',row.names=F,col.names=F)
write.table(mds %>% select(ID,Egg) %>% arrange(Egg) %>% filter(grepl('ECC',Egg)),file='CC.pop',quote=F,sep='\t',row.names=F,col.names=F)
write.table(mds %>% select(ID,Egg) %>% arrange(Egg) %>% filter(grepl('ECO',Egg)),file='CO.pop',quote=F,sep='\t',row.names=F,col.names=F)

# Plot the individuals
#jitter points up to 1 lat/long for viewing
mds = md %>% mutate(LatJit = jitter(Latitude,amount =2),
                     LonJit = jitter(Longitude,amount=2))
sites = st_as_sf(mds, coords = c("LonJit", "LatJit"), 
                 crs = 4326, agr = "constant")

#set up map and convert df to coordinate frame
world = map_data("world")

imm_spat = ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), col='grey90', fill='white') +
  geom_sf(data = sites, 
          aes(fill=Egg,shape=Egg),
          size=3,show.legend = T) +
  xlab('')+ylab('')+
  coord_sf(xlim = c(min(mds$Longitude)-5, max(mds$Longitude)+5), 
           ylim = c(min(mds$Latitude)-5, max(mds$Latitude)+5), expand = FALSE)+
  scale_fill_manual(values=mds$EggCol,breaks=mds$Egg)+
  scale_shape_manual(values=mds$EggShape,breaks=mds$Egg)+
  theme_classic()+
  facet_grid(SpeciesShort ~.)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"))+
  theme(legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))+
  annotation_scale(line_width=0.5)
imm_spat

ggsave('~/symlinks/host/figures/20250801_EggHunt_Spatial-FemalesNestlings-AllvONE.pdf',imm_spat,
       dpi=300,height=6,width=8)
```

#### Plot 

Same as above...  use the same filtered VCFs, because vcftools will re-filter based on max-missingness within the cohort:

````R
#### Plot 10KB FST Gene Scan
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/scan_n4/nestlings_females_N3_vs_all/')
.libPaths('~/r_libs')
library(karyoploteR)
library(ggpubr)
library(scales)
library(tidyverse)
library(viridis)
library(RColorBrewer)
library(zoo)
library(data.table)

md = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt')

##### Load and prep FST   ######
tidy_bp = fread('20250802_bFST_Nestling.txt.gz')
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
eggcols <- md %>% arrange(EggOrder) %>% select(Egg,EggCol) %>% na.omit %>% unique %>% filter(Egg %in% unique(tidy_bp$Egg))
tidy_bp$Egg <- factor(tidy_bp$Egg, levels=eggcols$Egg)
tidy_bp$AvZ <- factor(tidy_bp$AvZ,levels=c('Autosome','Z','MT','W'))

###### GWAS #######

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
#karyoA <- karyoA %>% mutate(color = ifelse(grepl('^ECO',Egg) & color == 'grey70','grey40', color))

pdf('~/symlinks/host/figures/20250731_KARYOPLOT-bFST10KB-AllvONE.-labels.pdf',height=3.5,width=7.5)
png('~/symlinks/host/figures/20250802_KARYOPLOT-bFST10KB-AllvONE-NESTLING.png',height=4,width=7.5,units='in',res=300)
pp = getDefaultPlotParams(plot.type=4)
pp$leftmargin = 0.1
kp = plotKaryotype(plot.type=4, 
                   genome = G,
                   #genome = keepSeqlevels(G, "W", pruning.mode = "coarse"),
                   labels.plotter = NULL,
                   plot.params = pp)
kpAddChromosomeNames(kp, yoffset = -5,cex=0.7)
kpAddBaseNumbers(kp,tick.dist = 25000000,minor.ticks = FALSE,cex=0.4)

#Loop through each track (all species together, CC, and CO), and plot them on their own layer 
counter = 0; tracks = 7
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
leg <- md %>% select(Egg,EggCol) %>% na.omit %>% distinct
genes <- read_tsv('Gene_Lookup.sorted.bed',col_names = F) %>% mutate(X1 = gsub('chr_','',X1))
names(genes) <- c('chr','start','end','strand','gene')
top_hit <- tidy_bp %>%
  group_by(Egg,AvZ) %>%  # Only group by Egg, not chr
  filter(grepl('MT',chr)) %>% 
  slice_max(order_by = bFST, n = 1, with_ties = TRUE) %>%  # Get top hits per egg 
  ungroup() %>%
  select(chr, start, end, Egg, bFST) %>%
  data.frame() %>% filter(grepl('ECC',Egg))
tg <- genes %>% filter(chr == 'MT') %>% mutate(gene = gsub('ID=','',gene))
tg <- tg %>%
  #interval_left_join(tg, by = c("start" = "start", "end" = "end")) %>%
  #select(chr=chr.y,start=start.y,end=end.y,gene,Egg,strand) %>% filter(!grepl('ECO',Egg)) %>% 
  mutate(xstart = if_else(strand == "+", start, end),xend   = if_else(strand == "+", end, start))
ceiling <- tidy_bp %>% filter(chr == 'MT') %>% summarize(ymax=max(bFST)*.99) %>% pull(ymax)
mtp <- tidy_bp %>% filter(chr == 'MT') %>%
  ggplot(aes(x=start,y=bFST,col=Egg))+
  geom_segment(data = tg,aes(x = xstart, xend = xend, y = ceiling*0.95, yend = ceiling*0.95),
               arrow = arrow(length = unit(0.07, "cm"),ends = "last",type = "closed"),inherit.aes = FALSE,linewidth = 0.3,color = "gray40") +
  geom_text(data=tg,aes(x=start+((end-start)/2),y=ceiling,label=gene),angle=30, size=1, inherit.aes=FALSE,hjust = 0,vjust=1)+
  geom_text(data=top_hit,aes(x=start+50),y=Inf,label='*',col='blue',alpha=1,size=3,vjust=1)+
  scale_color_manual(values=leg$EggCol,breaks=leg$Egg)+
  geom_line()+scale_x_continuous(breaks = pretty_breaks(n = 3),labels = function(x) paste0(round(x / 1e3, 1), "-Kb")) +scale_y_continuous(breaks=pretty_breaks(n=3))+
  theme_classic(base_size=8)+ylab('')+xlab('')+coord_cartesian(clip='off')
mtp

# First, inspect mtDNA 
outliers <- tidy_bp %>%
  group_by(Egg,AvZ) %>%  # Only group by Egg, not chr
  filter(!grepl('scaf|MT|W',chr) & bFST > 0) %>% 
  slice_max(order_by = bFST, n = 1, with_ties = TRUE) %>%  # Get top hits per egg 
  ungroup() %>%
  select(chr, start, end, Egg, bFST) %>%
  data.frame() 
outliers
#write.table(outliers,file='20250731_bFST_outliers-AllvONE.txt',quote=F,sep='\t',row.names=F,col.names=F)
#bedtools sort -i 20250731_bFST_outliers-AllvONE.txt > 20250731_bFST_outliers-AllvONE.sorted.txt


##### ECC1/ECC6 #####
# ECC6 chr2 : BLVRA
zchr=2
zstart=104840001;zs=zstart-5e5;ze=zstart+5e5
tg <- genes %>% filter(chr == zchr & start > zs & end < ze) %>%  mutate(xstart = if_else(strand == "+", start, end),xend   = if_else(strand == "+", end, start))
ceiling <- tidy_bp %>% filter(chr == zchr & start > zs & end < ze) %>% summarize(ymax=max(bFST)*.99) %>% pull(ymax)
p1 <- tidy_bp %>% filter(chr == zchr & start > zs & end < ze) %>%
  ggplot(aes(x=start,y=bFST,col=Egg))+
  geom_segment(data = tg,aes(x = xstart, xend = xend, y = ceiling*0.95, yend = ceiling*0.95),
               arrow = arrow(length = unit(0.07, "cm"),ends = "last",type = "closed"),inherit.aes = FALSE,linewidth = 0.3,color = "gray40") +
  geom_text(data=tg,aes(x=start+((end-start)/2),y=ceiling,label=gene),angle=30, size=1, inherit.aes=FALSE,hjust = 0,vjust=1)+
  geom_text(x=zstart+5000,y=Inf,label='*',col='blue',alpha=1,size=3,vjust=1)+
  scale_color_manual(values=leg$EggCol,breaks=leg$Egg)+
  ggtitle(paste0(zchr,' ',zstart,'ECC6:BLVRA,VOPP1'))+
  geom_line()+scale_x_continuous(breaks = pretty_breaks(n = 3),labels = function(x) paste0(round(x / 1e6, 1), "-Mb")) +scale_y_continuous(breaks=pretty_breaks(n=3))+
  theme_classic(base_size=8)+ylab('')+xlab('')+coord_cartesian(clip='off')
p1

# ECC6 chrZ : F2R/FR2L/FR2RL1
deets <- tidy_bp %>% filter(Egg == 'ECC6' & !grepl('scaf',chr)) %>% group_by(AvZ) %>% slice_max(bFST) %>% arrange(desc(bFST)) %>% select(chr,start,end)
deets
zchr='Z'
zstart=18030001;zs=zstart-2e5;ze=zstart+2e5
tg <- genes %>% filter(chr == zchr & start > zs & end < ze) %>%  mutate(xstart = if_else(strand == "+", start, end),xend   = if_else(strand == "+", end, start))
ceiling <- tidy_bp %>% filter(chr == zchr & start > zs & end < ze) %>% summarize(ymax=max(bFST)*.99) %>% pull(ymax)
pz1 <- tidy_bp %>% filter(chr == zchr & start > zs & end < ze) %>%
  ggplot(aes(x=start,y=bFST,col=Egg))+
  geom_segment(data = tg,aes(x = xstart, xend = xend, y = ceiling*0.95, yend = ceiling*0.95),
               arrow = arrow(length = unit(0.07, "cm"),ends = "last",type = "closed"),inherit.aes = FALSE,linewidth = 0.3,color = "gray40") +
  geom_text(data=tg,aes(x=start+((end-start)/2),y=ceiling,label=gene),angle=30, size=1, inherit.aes=FALSE,hjust = 0,vjust=1)+
  geom_text(x=zstart+5000,y=Inf,label='*',col='blue',alpha=1,size=3,vjust=1)+
  geom_text(x=21370001+5000,y=Inf,label='*',col='blue',alpha=1,size=3,vjust=1)+ # ECO4
  scale_color_manual(values=leg$EggCol,breaks=leg$Egg)+
  ggtitle(paste0(zchr,' ',zstart,'ECC6/ECO4:PTPRD'))+
  geom_line()+scale_x_continuous(breaks = pretty_breaks(n = 3),labels = function(x) paste0(round(x / 1e6, 1), "-Mb")) +scale_y_continuous(breaks=pretty_breaks(n=3))+
  theme_classic(base_size=8)+ylab('')+xlab('')+coord_cartesian(clip='off')
pz1

# ECC8
deets <- tidy_bp %>% filter(Egg == 'ECC8' & !grepl('scaf|W|MT',chr)) %>% group_by(chr) %>% slice_max(bFST) %>% arrange(desc(bFST)) %>% head(3) %>% select(chr,start,end)
deets
zchr=36
zstart=150001;zs=zstart-1.5e5;ze=zstart+1.5e5
tg <- genes %>% filter(chr == zchr & start > zs & end < ze) %>%  mutate(xstart = if_else(strand == "+", start, end),xend   = if_else(strand == "+", end, start))
ceiling <- tidy_bp %>% filter(chr == zchr & start > zs & end < ze) %>% summarize(ymax=max(bFST)*.99) %>% pull(ymax)
p2 <- tidy_bp %>% filter(chr == zchr & start > zs & end < ze) %>%
  ggplot(aes(x=start,y=bFST,col=Egg))+
  geom_segment(data = tg,aes(x = xstart, xend = xend, y = ceiling*0.95, yend = ceiling*0.95),
               arrow = arrow(length = unit(0.07, "cm"),ends = "last",type = "closed"),inherit.aes = FALSE,linewidth = 0.3,color = "gray40") +
  geom_text(data=tg,aes(x=start+((end-start)/2),y=ceiling,label=gene),angle=30, size=1, inherit.aes=FALSE,hjust = 0,vjust=1)+
  geom_text(x=zstart+5000,y=ceiling,label='*',col='blue',alpha=1,size=3)+
  scale_color_manual(values=leg$EggCol,breaks=leg$Egg)+
  ggtitle(paste0(zchr,' ',zstart,'ECC8:FIS1'))+
  geom_line()+scale_x_continuous(breaks = pretty_breaks(n = 3),labels = function(x) paste0(round(x / 1e6, 1), "-Mb")) +scale_y_continuous(breaks=pretty_breaks(n=3))+
  theme_classic(base_size=8)+ylab('')+xlab('')+coord_cartesian(clip='off')
p2

# ECC10 & ECC8 chrZ : RLN3
deets <- tidy_bp %>% filter(Egg == 'ECC10' & !grepl('scaf|W|MT',chr)) %>% group_by(chr) %>% slice_max(bFST) %>% arrange(desc(bFST)) %>% head(3) %>% select(chr,start,end)
deets
zchr='Z'
zstart=18980001;zs=zstart-2e5;ze=zstart+2e5
tg <- genes %>% filter(chr == zchr & start > zs & end < ze) %>%  mutate(xstart = if_else(strand == "+", start, end),xend   = if_else(strand == "+", end, start))
ceiling <- tidy_bp %>% filter(chr == zchr & start > zs & end < ze) %>% summarize(ymax=max(bFST)*.99) %>% pull(ymax)
p3 <- tidy_bp %>% filter(chr == zchr & start > zs & end < ze) %>%
  ggplot(aes(x=start,y=bFST,col=Egg))+
  geom_segment(data = tg,aes(x = xstart, xend = xend, y = ceiling*0.95, yend = ceiling*0.95),
               arrow = arrow(length = unit(0.07, "cm"),ends = "last",type = "closed"),inherit.aes = FALSE,linewidth = 0.3,color = "gray40") +
  geom_text(data=tg,aes(x=start+((end-start)/2),y=ceiling,label=gene),angle=30, size=1, inherit.aes=FALSE,hjust = 0,vjust=1)+
  geom_text(x=zstart+5000,y=ceiling,label='*',col='blue',alpha=1,size=3)+
  scale_color_manual(values=leg$EggCol,breaks=leg$Egg)+
  ggtitle(paste0(zchr,' ',zstart,'ECC10:RLN3'))+
  geom_line()+scale_x_continuous(breaks = pretty_breaks(n = 3),labels = function(x) paste0(round(x / 1e6, 1), "-Mb")) +scale_y_continuous(breaks=pretty_breaks(n=3))+
  theme_classic(base_size=8)+ylab('')+xlab('')+coord_cartesian(clip='off')
p3

# ECO1 chr35 : FCGBP
zchr=35
zstart=580001;zs=zstart-1e5;ze=zstart+1e5
tg <- genes %>% filter(chr == zchr & start > zs & end < ze) %>%  mutate(xstart = if_else(strand == "+", start, end),xend   = if_else(strand == "+", end, start))
ceiling <- tidy_bp %>% filter(chr == zchr & start > zs & end < ze) %>% summarize(ymax=max(bFST)*.99) %>% pull(ymax)
p4 <- tidy_bp %>% filter(chr == zchr & start > zs & end < ze) %>%
  ggplot(aes(x=start,y=bFST,col=Egg))+
  geom_segment(data = tg,aes(x = xstart, xend = xend, y = ceiling*0.95, yend = ceiling*0.95),
               arrow = arrow(length = unit(0.07, "cm"),ends = "last",type = "closed"),inherit.aes = FALSE,linewidth = 0.3,color = "gray40") +
  geom_text(data=tg,aes(x=start+((end-start)/2),y=ceiling,label=gene),angle=30, size=1, inherit.aes=FALSE,hjust = 0,vjust=1)+
  geom_text(x=zstart+5000,y=Inf,label='*',col='blue',alpha=1,size=3,vjust=1)+
  scale_color_manual(values=leg$EggCol,breaks=leg$Egg)+
  ggtitle(paste0(zchr,' ',zstart,'ECO1:FCGBP'))+
  geom_line()+scale_x_continuous(breaks = pretty_breaks(n = 3),labels = function(x) paste0(round(x / 1e6, 1), "-Mb")) +scale_y_continuous(breaks=pretty_breaks(n=3))+
  theme_classic(base_size=8)+ylab('')+xlab('')+coord_cartesian(clip='off')
p4

##### ECO2 #####
# ECO2 chr6 : SLC40A1
deets <- tidy_bp %>% filter(Egg == 'ECO2' & !grepl('scaf',chr)) %>% group_by(chr) %>% slice_max(bFST) %>% arrange(desc(bFST)) %>% head(3) %>% select(chr,start,end)
deets
zchr=6
zstart=35090001;zs=zstart-2e5;ze=zstart+2e5
tg <- genes %>% filter(chr == zchr & start > zs & end < ze) %>%  mutate(xstart = if_else(strand == "+", start, end),xend   = if_else(strand == "+", end, start))
ceiling <- tidy_bp %>% filter(chr == zchr & start > zs & end < ze) %>% summarize(ymax=max(bFST)*.99) %>% pull(ymax)
p5 <- tidy_bp %>% filter(chr == zchr & start > zs & end < ze) %>%
  ggplot(aes(x=start,y=bFST,col=Egg))+
  geom_segment(data = tg,aes(x = xstart, xend = xend, y = ceiling*0.95, yend = ceiling*0.95),
               arrow = arrow(length = unit(0.07, "cm"),ends = "last",type = "closed"),inherit.aes = FALSE,linewidth = 0.3,color = "gray40") +
  geom_text(data=tg,aes(x=start+((end-start)/2),y=ceiling,label=gene),angle=30, size=1, inherit.aes=FALSE,hjust = 0,vjust=1)+
  geom_text(x=zstart+5000,y=Inf,label='*',col='blue',alpha=1,size=3,vjust=1)+
  scale_color_manual(values=leg$EggCol,breaks=leg$Egg)+
  ggtitle(paste0(zchr,' ',zstart,'ECO2:SLC40A1'))+
  geom_line()+scale_x_continuous(breaks = pretty_breaks(n = 3),labels = function(x) paste0(round(x / 1e6, 1), "-Mb")) +scale_y_continuous(breaks=pretty_breaks(n=3))+
  theme_classic(base_size=8)+ylab('')+xlab('')+coord_cartesian(clip='off')
p5

##### ECO3 #####
# ECO3 chrZ : AIG1
deets <- tidy_bp %>% filter(Egg == 'ECO3' & !grepl('scaf',chr)) %>% group_by(chr) %>% slice_max(bFST) %>% arrange(desc(bFST)) %>% head(5) %>% select(chr,start,end)
deets
zchr=3
zstart=47010001;zs=zstart-1.5e5;ze=zstart+1.5e5
tg <- genes %>% filter(chr == zchr & start > zs & end < ze) %>%  mutate(xstart = if_else(strand == "+", start, end),xend   = if_else(strand == "+", end, start))
ceiling <- tidy_bp %>% filter(chr == zchr & start > zs & end < ze) %>% summarize(ymax=max(bFST)*.99) %>% pull(ymax)
p6 <- tidy_bp %>% filter(chr == zchr & start > zs & end < ze) %>%
  ggplot(aes(x=start,y=bFST,col=Egg))+
  geom_segment(data = tg,aes(x = xstart, xend = xend, y = ceiling*0.95, yend = ceiling*0.95),
               arrow = arrow(length = unit(0.07, "cm"),ends = "last",type = "closed"),inherit.aes = FALSE,linewidth = 0.3,color = "gray40") +
  geom_text(data=tg,aes(x=start+((end-start)/2),y=ceiling,label=gene),angle=30, size=1, inherit.aes=FALSE,hjust = 0,vjust=1)+
  scale_color_manual(values=leg$EggCol,breaks=leg$Egg)+
  geom_line()+scale_x_continuous(breaks = pretty_breaks(n = 3),labels = function(x) paste0(round(x / 1e6, 1), "-Mb")) +scale_y_continuous(breaks=pretty_breaks(n=3))+
  ggtitle(paste0(zchr,' ',zstart,'ECO3:AIG3'))+
  geom_text(x=zstart+5000,y=Inf,label='*',col='blue',alpha=1,size=3,vjust=1)+
  theme_classic(base_size=8)+ylab('')+xlab('')+coord_cartesian(clip='off')
p6


pdf('~/symlinks/host/figures/20250801_Candidates-AllvONE-NESTLINGS.pdf',height=4,width=6)
ggarrange(p1,p2,p3,p4,p5,p6,nrow=3,ncol=2,common.legend = TRUE)
dev.off()

````

## FST: Base-pair level fixed SNPs [Fig 4]

Hierarchical contrasts, egg:haplotype comparisons. Identify samples and name the phylogenetic comparisons: 

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

Base pair resolution 

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

### Plot 

Identify candidate genes with nonsynonymous mutations that are fixed across contrasts (e.g. blue emergence):

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

### Plot NDUFAF4 Individual SNPs

```bash
#### Identify candidates from female-only FST scan 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/manyhost_hunt/males')
.libPaths('~/r_libs/')
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

# within dufaf 4
ndu <- tidy_dat %>% filter(gene == 'NDUFAF4' | gene == 'LOC128850245') %>% filter(grepl('Ancient_Blue|Contemporary_Reversion',Phylo))
ndu_snps <- ndu %>%
  mutate(Phylo = ifelse(Phylo == 'Ancient_Blue', 'Blue Emergence', 'Blue Reversion')) %>%
  ggplot(aes(y = site, x = FST, fill = Phylo)) +
  geom_point(size = 2.5,pch=21,col='black',position = position_jitter(width = 0.05, height = 0.1, seed = 42)) +
  facet_grid(chr ~ effect, scales = 'free', space = 'free') +
  scale_fill_manual(values = brewer.pal(3, 'Set2')[c(3, 1)]) +
  theme_bw(base_size = 10) +
  theme(panel.spacing = unit(0.5, "lines"), clip = "on")
ndu_snps
ggsave('~/symlinks/host/figures/20250805_NDUFAF4_chr3-chrW-SNPs_FST.pdf',ndu_snps,height=4,width=7)

```


