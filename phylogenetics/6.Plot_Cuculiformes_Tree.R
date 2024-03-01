#### Plot deep cuculiformes tree 
setwd('~/merondun/cuculus_host/phylogenetics/')
.libPaths('~/mambaforge/envs/R/lib/R/library')
library(ggtree)
library(ape)
library(tidyverse)
library(treeio)
library(viridis)
library(ggpubr)
library(RColorBrewer)
library(phytools)
library(ggnewscale)

md = read_tsv('Cuckoo_Full_Metadata_2023OCT3.txt')

###### CIRCULAR W TREE  #######
#plot tree 
iqtree = read.iqtree('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/deep_tree/trees/W.SNP.DP3-AC1-MQ40.min4.phy.contree.NoStriped.contree')
iqtr = midpoint.root(as.phylo(iqtree))

#plot with outgroups 
circ = ggtree(iqtr, layout = "circular",branch.length='none') %<+% md  +
  geom_nodepoint(mapping=aes(subset=(as.numeric(label) >= 95)),col='black',fill='grey90',pch=23,size=1.5,show.legend=F)+
  geom_tippoint(aes(fill = Hap,shape=Species),size=4)+
  scale_fill_manual('W Haplogroup',values=md$HapCol,breaks=md$Hap)+
  scale_shape_manual(values=c(7,21,10,9,25,4,8,24,11,12))+
  guides(fill=guide_legend(nrow=4,override.aes=list(shape=21)),
         shape=guide_legend(nrow=5))+
  theme(legend.text = element_text(size = 8),legend.title = element_text(size = 10),
        legend.key.size = unit(0.2, "cm"),    legend.position = 'top')
circ

apetree = ggtree(iqtr, layout = "ape") %<+% md  +
  geom_nodepoint(mapping=aes(subset=(as.numeric(label) >= 95)),col='black',fill='grey90',pch=23,size=1.5,show.legend=F)+
  geom_tippoint(aes(fill = Hap,shape=Species),size=4)+
  scale_fill_manual('W Haplogroup',values=md$HapCol,breaks=md$Hap)+
  scale_shape_manual(values=c(7,21,10,9,25,4,8,24,11,12))+
  guides(fill=guide_legend(nrow=4,override.aes=list(shape=21)),
         shape=guide_legend(nrow=5))+
  theme(legend.text = element_text(size = 8),legend.title = element_text(size = 10),
        legend.key.size = unit(0.2, "cm"),    legend.position = 'top')+
  geom_treescale(offset=0.001)
apetree

pdf('~/merondun/cuculus_host/phylogenetics/ML_Tree_Deep_Cuculiformes_2024FEB29.pdf',height=6,width=6)
ggarrange(circ,apetree,common.legend = TRUE,nrow=2,heights=c(0.75,0.25))
dev.off()