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

