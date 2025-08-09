# Analyses: Egg Associations 

## Distance-based redundancy analyses; dbRDA

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

### Sensitivity: Nestlings optatus

```R
#### dbRDA for e.g. egg type ~ maternal haplogroup associations: sensitivity only nestlings optatus 
setwd('~/EvoBioWolf/CUCKOO_gentes/gens_associations/')
.libPaths('~/r_libs')
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
mdf %>% filter(Analysis_GensAssociationsGeneScan == 1) %>% count(SpeciesShort,Sex) 

sp <- 'CO'

#only grab samples with known egg (hash out filter for related individuals)
md_egg = mdf %>%
  filter(Analysis_GensAssociationsGeneScan == 1 & SpeciesShort == sp & Age == 'Young') %>%
  drop_na(Egg) %>%
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

dbrf = rbindlist(dbr) %>% as_tibble %>% mutate(Filter = 'Nestlings Only')

db_save =  dbrf %>% select(-Df,SumOfSqs) %>% mutate(padj = p.adjust(p,method='bonferroni'))
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

ggsave('~/symlinks/host/figures/20250807_dbRDA-Sensitivity-NESTLINGS-optatus.pdf',
       dbrda_pnt,height=2,width=1.75,dpi=300)
```

