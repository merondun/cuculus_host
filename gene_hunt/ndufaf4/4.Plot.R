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

