#### Plot NDUFAF4 Coverage 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/cuculiformes_NDUFAF4_blasting/Alignment_Based_Approach/Whole_Genome_Alignment/coverage/')
.libPaths('~/mambaforge/envs/R/lib/R/library')
library(tidyverse)
library(viridis)

# Load in coverage data 
c <- read_tsv('NDUFAF4_MQ60_Coverage_2024JULY25_CDS_WGA.txt')
c$Species <- factor(c$Species,levels=c('Cuculus_canorus','Cuculus_micropterus','Cuculus_poliocephalus','Clamator_glandarius','Coccyzus_lansbergi','Piaya_cayana','Dromococcyx_pavoninus','Tapera_naevia','Geococcyx_californianus'))
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

pdf('../../../../figures/20240807_NDUFAF4_Coverage.pdf',height=2.75,width=2.25)
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

png('NDUFAF4_Coverage-IncludingCDS.png',units='in',res=300,height=3.5,width=7)
cp
dev.off()

# Plot Tree afterwards
library(ggtree)
library(treeio)

# Free alignment
iqtree = read.iqtree('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/cuculiformes_NDUFAF4_blasting/Alignment_Based_Approach/Whole_Genome_Alignment/bams/gene_alignment/free_align/NDUFAF4_Gapless50.fa.contree')

iqtr = root(as.phylo(iqtree),'Geococcyx_californianus_SRR9994302_chr_3_25511697_25515677')

md <- read_tsv('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/cuculiformes_NDUFAF4_blasting/Alignment_Based_Approach/Whole_Genome_Alignment/bams/Metadata.txt')
g <- ggtree(iqtr,layout='rectangular')
#g <- ggtree(iqtr,layout='ape')
g$data <- g$data %>% mutate(Chromosome = ifelse(grepl('chr_W',label),'chrW','chr3'),
                            label = gsub('_chr_W.*','',label),
                            label = gsub('_chr_3.*','',label))
g$data <- left_join(g$data,md %>% dplyr::rename(label = ID))

#plot with outgroups 
gtree <- g +
  geom_tippoint(aes(col=Chromosome,shape=Species),size=1,stroke=1)+
  geom_nodelab(aes(label=node),geom = 'text',size=1.5)+ 
  scale_shape_manual(values=c(1,2,3,4,5,6,7,8,9))+
  xlim(c(0,1))+
  #geom_tiplab(aes(label=Species),size=2,offset = 0.05)+
  geom_tiplab(size=2,offset = 0.05)+
  theme(legend.position='none')
gtree

pdf('../../../../figures/20240806_NDUFAF4-ParalogTree.pdf',height=2.75,width=2.5)
gtree
dev.off()
