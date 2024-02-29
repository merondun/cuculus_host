#### Plot Trees, Assign Haplogroups 
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

md = read_tsv('Cuckoo_Full_Metadata_2023OCT3.txt')
md = md %>% mutate(HostMiss = ifelse(is.na(HostParentShort),'Missing','Present'))

###### MT #######
## Plot Tree 
iqtree = read.iqtree('no_het/chr_MT.IF-GF-MM1-NoHet.min4.phy.varsites.phy.contree')
#iqtr = root(iqtree, outgroup = "387_CP_MBW_RUS_F", resolve.root = TRUE) # rooted..
iqtr = drop.tip(iqtree, c("387_CP_MBW_RUS_F", "386_CP_MBW_RUS_M")) #drop outgroups 

#plot cladogram, first find nodes to collapse 
insp = ggtree(iqtr, layout = "ape") %<+% md  + 
  #geom_nodelab(aes(label=node),geom = 'label',size=3)+  geom_nodepoint(mapping=aes(subset=(SH_aLRT >= 95 & UFboot > 95)),col='black',fill='grey90',pch=23,size=1.5,show.legend=F)+
  geom_tippoint(aes(fill = HostParentShort,shape=Plumage,alpha=HostMiss),size=5)+
  scale_shape_manual(values=c(24,21,25))+
  #scale_fill_manual('W Haplogroup',values=md$Wcol,breaks=md$W)+
  #scale_fill_manual('Distance Group',values=md$KDCol,breaks=md$KDist)+
  #geom_tiplab()+
  guides(fill=guide_legend(nrow=2,override.aes=list(shape=21)))
insp
pdf('figures/MT_TREE-Inspection-NoHet_2024FEB23.pdf',height=25,width=55)
insp
dev.off()


###### W #######
## Plot Tree 
iqtree = read.iqtree('no_het/chr_W.IF-GF-MM1-NoHet.min4.phy.varsites.phy.contree')
#iqtr = root(iqtree, outgroup = "387_CP_MBW_RUS_F", resolve.root = TRUE) # rooted..
iqtr = drop.tip(iqtree, c("387_CP_MBW_RUS_F", "386_CP_MBW_RUS_M")) #drop outgroups 

#plot cladogram, first find nodes to collapse 
insp = ggtree(iqtr, layout = "dendrogram") %<+% md  + 
  #geom_nodelab(aes(label=node),geom = 'label',size=3)+  geom_nodepoint(mapping=aes(subset=(SH_aLRT >= 95 & UFboot > 95)),col='black',fill='grey90',pch=23,size=1.5,show.legend=F)+
  #geom_tippoint(aes(fill = HostParentShort,shape=Plumage,alpha=HostMiss),size=5)+
  geom_tippoint(aes(fill = Hap,shape=Plumage,alpha=HostMiss),size=5)+
  scale_shape_manual(values=c(24,21,25))+
  scale_fill_manual('W Haplogroup',values=md$HapCol,breaks=md$Hap)+
  #scale_fill_manual('Distance Group',values=md$KDCol,breaks=md$KDist)+
  #geom_tiplab()+
  guides(fill=guide_legend(nrow=2,override.aes=list(shape=21)))
insp
pdf('figures/W_TREE-Inspection-NoHet_2024FEB23.pdf',height=25,width=55)
insp
dev.off()


#find the nodes in which all descedents will be assigned a haplogroup 
labs = data.frame(Nodes = c(219,201,165,308,290,277,231))
ds = NULL; counter = 0
for (node in labs$Nodes){ counter = counter + 1 
d = data.frame(ID = get_taxa_name(insp, node = node), W = paste0('W',counter))
ds = rbind(ds,d)
}

#must classify the REDS alone because they are closer to root 
ds = rbind(ds,data.frame(ID = c('023_CC_RED_FIN_F'),
                         W = c('W1')))
#re-inspect 
inspW = ggtree(iqtr, layout = "ape") %<+% ds2  + 
  geom_tippoint(aes(col = W),size=2,pch=16)+
  scale_color_manual('W Haplogroup',values=ds2$Wcol,breaks=ds2$W)+
  geom_treescale(offset = 0.01)
inspW
insp$data %>% filter(is.na(W) & isTip == TRUE)

#add colors
cols = ds %>% select(W) %>% unique() %>% mutate(Wcol = brewer.pal(7,'Paired'))
ds2 = left_join(ds,cols)
write_tsv(ds2,file='figures/W_Haplogroups_VCF_2023DEC19.txt')

#plot for supplementary
circW = ggtree(iqtr, layout = "ape") %<+% md  + 
  geom_nodepoint(mapping=aes(subset=(SH_aLRT >= 95 & UFboot > 95)),col='black',fill='grey90',pch=23,size=1.5,show.legend=F)+
  geom_tippoint(aes(fill = Hap),pch=21,size=2)+
  scale_fill_manual('W Haplogroup',values=md$HapCol,breaks=md$Hap)+
  guides(fill=guide_legend(nrow=2,override.aes=list(shape=21)))+
  theme(legend.position='top')
circW

####### MT #######
## Plot Tree 
iqtree = read.iqtree('trees_full/chr_MT.IF-GF-MM1-AC2.min4.phy.contree')
iqtr = drop.tip(iqtree, c("387_CP_MBW_RUS_F", "386_CP_MBW_RUS_M")) #drop outgroups 

#plot hap, circular
insp = ggtree(iqtr, layout = "rectangular") %<+% ds2  + 
  geom_nodelab(aes(label=node),geom = 'label')+  geom_nodepoint(mapping=aes(subset=(SH_aLRT >= 95 & UFboot > 95)),col='black',fill='grey90',pch=23,size=1.5,show.legend=F)+
  geom_tippoint(aes(fill = W),size=5,pch=21)+
  geom_tiplab()+
  scale_shape_manual(values=c(24,21,25))+  scale_fill_manual('W Haplogroup',values=ds2$Wcol,breaks=ds2$W)+
  guides(fill=guide_legend(nrow=2,override.aes=list(shape=21)))
insp
pdf('figures/MT_TREE-Inspection_2023DEC19.pdf',height=30,width=70)
insp
dev.off()

labs = data.frame(Nodes = c(427,512,539,562,315,380,413),
                  MT_W = c('W7','W6','W4','W5','W3','W2','W1'))
ds = NULL; counter = 0
for (node in labs$Nodes){ counter = counter + 1 
d = data.frame(ID = get_taxa_name(insp, node = node), MT_W = labs$MT_W[counter])
ds = rbind(ds,d)
}

#custom subclades
allsamps = insp$data %>% filter(isTip == TRUE) %>% select(label) %>% dplyr::rename(ID = label)
miss = left_join(allsamps,ds) %>% filter(is.na(MT_W)) 
#they are all W1, BUT samples 023, 086, 149, and 141 are massive chrMT outliers, likely high NUMT 
dsm = left_join(allsamps,ds) %>% replace_na(list(MT_W = 'W1'))

ds3 = left_join(dsm %>% dplyr::rename(Hap = MT_W),cols %>% dplyr::rename(Hap = W, HapCol = Wcol))
#add W for inferring
ds4 = left_join(ds3,ds2) %>% mutate(W_dat = ifelse(is.na(W),'Male','Female'))

#re-inspect 
inspMT = ggtree(iqtr, layout = "ape") %<+% ds4  + 
  geom_tippoint(aes(col = Hap,shape=W_dat),size=2)+
  scale_shape_manual(values=c(16,8))+
  scale_color_manual('W Haplogroup',values=ds3$HapCol,breaks=ds3$Hap)+
  geom_treescale(offset = 0.01)
inspMT
write_tsv(ds3,file='figures/MT-W_Haplogroups_VCF_2023DEC19.txt')

#plot for supplementary
circMT = ggtree(iqtr, layout = "ape") %<+% md  + 
  geom_nodepoint(mapping=aes(subset=(SH_aLRT >= 95 & UFboot > 95)),col='black',fill='grey90',pch=23,size=1.5,show.legend=F)+
  geom_tippoint(aes(fill = Hap),pch=21,size=2)+
  scale_fill_manual('W Haplogroup',values=md$HapCol,breaks=md$Hap)+
  guides(fill=guide_legend(nrow=2,override.aes=list(shape=21)))+
  theme(legend.position='top')
circMT

#plot together
pdf('figures/Haplogroup_Assignments_2023DEC19.pdf',height=4,width=7)
ggarrange(inspMT,inspW,common.legend = TRUE,labels = c('mtDNA','W'))
dev.off()

#plot together
pdf('figures/Haplogroup_Assignments_2024FEB19.pdf',height=8,width=6)
ggarrange(circMT,circW,common.legend = TRUE,labels = c('mtDNA','W'),nrow = 1)
dev.off()
