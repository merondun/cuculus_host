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


