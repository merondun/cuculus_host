# Correlations between genetic and geographic distance

Estimating the strength of the relationship between genetic and geographic distance with separman's Rho and LMMs. 

1: Calculate geographic pairwise distance 

2: Calculate pairwise FST for autosomal and Z data 

3 & 4: AMOVA for distance of mtDNA and W chromosome haplotypes

5: Plot 

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

