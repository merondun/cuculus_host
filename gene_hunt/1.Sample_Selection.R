#### Spatial Plotting
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(viridis)
library(RColorBrewer)
library(ggpubr)
library(sf)
library(ggspatial)
library(spThin)
library(factoextra)
library(ggforce)

#Read in metadata: haplotype~egg approach 
md = read_tsv('~/merondun/cuculus_host/Metadata_Host.txt')
md = md %>% drop_na(Egg) %>% 
  filter(Sex == 'F' & Analysis_PopulationGenetics == 1) %>% 
  mutate(Group = paste0(Hap,'__',gsub('. ','',Egg)))
kept = md %>% filter(Analysis_PopulationGenetics == 1) %>% count(Group) %>% filter(n >= 4)
md = md %>% filter(Group %in% kept$Group)
retained = md %>% count(Group)
write.table(md$ID,file='/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/manyhost_hunt/ManyHost.list',quote=F,sep='\t',row.names=F,col.names=F)
write.table(md %>% select(ID,Group),file='/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/manyhost_hunt/ManyHost.pop',quote=F,sep='\t',row.names=F,col.names=F)

#grab males for CNV 
mdm = read_tsv('~/merondun/cuculus_host/Metadata_Host.txt')
mdm = mdm %>% drop_na(Egg) %>% 
  filter(Sex == 'M' & Analysis_PopulationGenetics == 1) %>% 
  mutate(Group = paste0(Hap,'__',gsub('. ','',Egg))) %>% 
  filter(Group %in% retained$Group & Analysis_PopulationGenetics == 1)
keptM = mdm %>% count(Group) %>% filter(n >= 3)
mdm = mdm %>% filter(Group %in% keptM$Group)
femcnv = md %>% filter(Group %in% mdm$Group) 
write.table(c(mdm$ID,femcnv$ID),file='/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/manyhost_hunt/coverage/CNV_Samples.list',quote=F,sep='\t',row.names=F,col.names=F)
cnv = rbind(mdm,md %>% filter(Group%in% mdm$Group)) %>% select(ID,Group)
write.table(cnv,file='/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/manyhost_hunt/coverage/CNV_Samples.pop',quote=F,sep='\t',row.names=F,col.names=F)

#jitter points up to 1 lat/long for viewing
md = md %>% mutate(LatJit = jitter(Latitude,amount =2),
                   LonJit = jitter(Longitude,amount=2))
sites = st_as_sf(md, coords = c("LonJit", "LatJit"), 
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

pdf('figures/GWAS_AllEgg_SpatialDistribution_2024MAR05.pdf',height=3,width=7)
imm_spat
dev.off()

#Read in metadata: Reversion approach 
md = read_tsv('Cuckoo_Full_Metadata_2023OCT3.txt')
md = md %>% filter(!is.na(Reversion))

#jitter points up to 1 lat/long for viewing
md = md %>% mutate(LatJit = jitter(Latitude,amount =2),
                   LonJit = jitter(Longitude,amount=2))
sites = st_as_sf(md, coords = c("LonJit", "LatJit"), 
                 crs = 4326, agr = "constant")
sites = sites %>% filter(!is.na(Reversion))

#set up map and convert df to coordinate frame
world = map_data("world")
cols = brewer.pal(4,'Paired')

egg_spat = ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), col='grey90', fill='white') +
  geom_sf(data = sites, 
          aes(fill=as.factor(Reversion)),
          size=3,alpha=0.8,show.legend = T,pch=21) +
  xlab('')+ylab('')+
  coord_sf(xlim = c(min(md$Longitude)-5, max(md$Longitude)+5), 
           ylim = c(min(md$Latitude)-5, max(md$Latitude)+5), expand = FALSE)+
  scale_fill_manual(values=c(cols[1],cols[2],cols[3],'black'))+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"))+
  theme(legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))+
  annotation_scale(line_width=0.5)

pdf('figures/REV_BNB_SpatialDistribution_2024FEB13.pdf',height=3,width=7)
egg_spat
dev.off()

write.table(md$ID,file='/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/manyhost_hunt/ManyHost.list',quote=F,sep='\t',row.names=F,col.names=F)
write.table(md %>% select(ID,Group),file='/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/manyhost_hunt/ManyHost.pop',quote=F,sep='\t',row.names=F,col.names=F)

