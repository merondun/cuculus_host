#### Collapse chrW Tree 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host')
.libPaths('~/mambaforge/envs/R/lib/R/library')
library(ggtree)
library(ape)
library(tidyverse)
library(treeio)
library(viridis)
library(adegenet)
library(ggpubr)
library(vcfR)
library(adegenet)
library(RColorBrewer)
library(phytools)
library(ggnewscale)

md = read_tsv('Cuckoo_Full_Metadata_2023OCT3.txt')

##### COLLAPSE W TREE #####
#midpoint root 
iqtree = read.iqtree('beast_dating/variant_only/ml_trees/chr_W_All.SNP.DP3-AC2-MQ40.min4.phy.contree')
iqtr = midpoint.root(as.phylo(iqtree))
p = ggtree(iqtr, layout = "ape") %<+% md 
pdf('figures/W_TREE-COLLAPSE_Inspection_2024FEB27.pdf',height=75,width=55)
p + geom_nodelab(aes(label=node),geom = 'label',size=3) +geom_tippoint(aes(fill = Hap),pch=21,size=4)+scale_fill_manual(values=md$HapCol,breaks=md$Hap)
dev.off()
hapcols = md %>% select(Hap,HapCol) %>% unique %>% arrange(Hap) %>% filter(!grepl('CM|CP',Hap)) %>% na.omit
#collapse
p = ggtree(iqtr, layout = "ape") %<+% md + geom_tippoint(aes(shape=SpeciesShort),size=3) + scale_shape_manual(values=c(21,4,8))
p2 = p %>% collapse(node=201) + geom_point2(aes(subset=(node==201)), shape=21, size=3, fill=hapcols[1,2]); p2
p2 = collapse(p2,node=182) + geom_point2(aes(subset=(node==182)), shape=21, size=3, fill=hapcols[[2,2]]); p2
p2 = collapse(p2,node=180) + geom_point2(aes(subset=(node==180)), shape=21, size=3, fill=hapcols[[3,2]]); p2
p2 = collapse(p2,node=304) + geom_point2(aes(subset=(node==304)), shape=4, size=3); p2
p2 = collapse(p2,node=292) + geom_point2(aes(subset=(node==292)), shape=21, size=3, fill=hapcols[[4,2]]);p2
p2 = collapse(p2,node=275) + geom_point2(aes(subset=(node==275)), shape=21, size=3, fill=hapcols[[5,2]]);p2
p2 = collapse(p2,node=262) + geom_point2(aes(subset=(node==262)), shape=21, size=3, fill=hapcols[[6,2]]);p2
p2 = collapse(p2,node=216) + geom_point2(aes(subset=(node==216)), shape=21, size=3, fill=hapcols[[7,2]]);p2
p2 = p2 + geom_treescale(x=0.05,offset = 0.01) + geom_nodepoint(mapping=aes(subset=(as.numeric(label) >= 95)),col='black',fill='grey90',pch=23,size=1.5,show.legend=F)

pdf('figures/W_TREE-Collapse_Haplogroups_2024FEB27.pdf',height=5,width=6)
p2
dev.off()