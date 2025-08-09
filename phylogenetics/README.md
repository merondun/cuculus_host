# Analyses: Matrilineal Phylogenetics

## Estimate Maternal Phylogenies

Subset samples. In general we have 2 sample sets: 

* one with females (n=87) and outgroups (n=6) together for BEAST (n=93) .
* one with all samples for mtDNA (n=136) and females for W (n=87). 

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --mem=30000mb
#SBATCH --cpus-per-task=10
#SBATCH --time=24:00:00

WD=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/phylogenetics/202503_fulln136
raw_vcfs=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/merged/snps_only
#mask with male-biased coverage
mask=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/Coverage_Masks/GCA_017976375.1_bCucCan1.pri_genomic.CHR.N75-DoubleCoverage.mask.bed

cd ${WD}

#mamba activate snps
# sbatch 1.Subset_Samples_Filter.sh
mkdir -p vcfs ml_trees

for CHR in chr_MT chr_W; do

#minimum coverage, LESS than this set to missing
MINDP=3
echo "FILTERING AND MERGING VARIANT SITES FOR ${CHR}"
bcftools view --force-samples --threads 10 --samples-file Samples.list -Ou ${raw_vcfs}/${CHR}.SNPS.vcf.gz | \
        bcftools view --types snps --min-alleles 2 --max-alleles 2 --threads 10 | \
        #remove SNPs in bad coverage regions
        bedtools subtract -header -a - -b ${mask} | \
        #set genotypes below MINDP to missing
        bcftools +setGT -- -t q -i "FMT/DP < ${MINDP}" -n "./." | \
        #set het genotypes to missing based on binomial test
        bcftools +setGT -Ou -- -t "b:AD<1e-5" -n "./." | \
        #set weakly het genotypes to major allele
        bcftools +setGT -Ou -- --target-gt q --new-gt M -i 'GT=="het"' | \
        #set to haploid, can skip this for most purposes
        #bcftools +fixploidy -Ou - -- -f 1 | \
        #update AC fields
        bcftools +fill-tags -Ou -- -t AC,AN | \
        bcftools view --min-alleles 2 --max-alleles 2 --types snps --min-ac 1 -Oz -o vcfs/${CHR}.SNP.DP3.vcf.gz
bcftools index --threads 10 vcfs/${CHR}.SNP.DP3.vcf.gz

#filter, include singletons
bcftools view vcfs/${CHR}.SNP.DP3.vcf.gz --min-alleles 2 --max-alleles 2 --min-ac 1 --max-af 0.999 --types snps -i "MQ>40 & F_MISSING < 0.1" -Oz -o vcfs/${CHR}.SNP.DP3-AC1-MQ40-MM1.vcf.gz
bcftools index --threads 10 vcfs/${CHR}.SNP.DP3-AC1-MQ40-MM1.vcf.gz

#filter, exclude singletons, note that this must be --MIN-AC 3, BECAUSE ITS PSEUDO-HAPLOID
bcftools view vcfs/${CHR}.SNP.DP3.vcf.gz --min-alleles 2 --max-alleles 2 --min-ac 3 --max-af 0.999 --types snps -i "MQ>40 & F_MISSING < 0.1" -Oz -o vcfs/${CHR}.SNP.DP3-AC2-MQ40-MM1.vcf.gz
bcftools index --threads 10 vcfs/${CHR}.SNP.DP3-AC2-MQ40-MM1.vcf.gz

#create tree
python ~/modules/vcf2phylip.py -i vcfs/${CHR}.SNP.DP3-AC1-MQ40-MM1.vcf.gz -f --output-folder ml_trees
iqtree --redo -keep-ident -T 10 -s ml_trees/${CHR}.SNP.DP3-AC1-MQ40-MM1.min4.phy --seqtype DNA -m "MFP+ASC" -alrt 1000 -B 1000
iqtree --redo -keep-ident -T 10 -s ml_trees/${CHR}.SNP.DP3-AC1-MQ40-MM1.min4.phy.varsites.phy --seqtype DNA -m "MFP+ASC" -alrt 1000 -B 1000

#create tree NO SINGLETONS
python ~/modules/vcf2phylip.py -i vcfs/${CHR}.SNP.DP3-AC2-MQ40-MM1.vcf.gz -f --output-folder ml_trees
iqtree --redo -keep-ident -T 10 -s ml_trees/${CHR}.SNP.DP3-AC2-MQ40-MM1.min4.phy --seqtype DNA -m "MFP+ASC" -alrt 1000 -B 1000
iqtree --redo -keep-ident -T 10 -s ml_trees/${CHR}.SNP.DP3-AC2-MQ40-MM1.min4.phy.varsites.phy --seqtype DNA -m "MFP+ASC" -alrt 1000 -B 1000

done

```



## Assign Haplogroups

Using fasta files from VCF, clsuter using tree:

```R
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
library(phytools)
library(ggnewscale)

md = read_tsv('Cuckoo_Full_Metadata_2023OCT3.txt')

###### W #######
## Plot Tree 
iqtree = read.iqtree('beast_dating/variant_only/ml_trees/chr_W_All.SNP.DP3-AC2-MQ40.min4.phy.contree')
iqtr = midpoint.root(as.phylo(iqtree))

#plot cladogram, first find nodes to collapse 
insp = ggtree(iqtr, layout = "dendrogram") %<+% md  + 
  geom_nodelab(aes(label=node),geom = 'label',size=3)+  geom_nodepoint(mapping=aes(subset=(as.numeric(label) >= 95)),col='black',fill='grey90',pch=23,size=1.5,show.legend=F)+
  geom_tippoint(aes(fill = KDist,shape=Plumage),size=5)+
  scale_shape_manual(values=c(24,21,25))+
  #scale_fill_manual('W Haplogroup',values=md$Wcol,breaks=md$W)+
  scale_fill_manual('Distance Group',values=md$KDCol,breaks=md$KDist)+
  geom_tiplab()+
  guides(fill=guide_legend(nrow=2,override.aes=list(shape=21)))
insp
pdf('figures/W_TREE-Inspection-AC2_2024FEB27.pdf',height=25,width=55)
insp
dev.off()

#find the nodes in which all descedents will be assigned a haplogroup 
labs = data.frame(Nodes = c(219,201,165,308,290,277,231))
ds = NULL; counter = 0
for (node in labs$Nodes){ counter = counter + 1 
d = data.frame(ID = get_taxa_name(insp, node = node), W = paste0('W',counter))
ds = rbind(ds,d)
}

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
```

## Plot Phylogenies

### Plot W and mtDNA Trees

Trees were estimated here: 

```bash
/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/phylogenetics/202503_fulln136

#full susbset with n=252
/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/phylogenetics/202503_fulln252
```

Plot many trees:

```R
#### Plot many W and MT trees with different sample subsets 
setwd('~/EvoBioWolf/CUCKOO_gentes/phylogenetics/trees/')
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

md = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt')

###### CIRCULAR W TREE  #######
#plot tree 
trees = list.files('.',pattern='.*AC1.*contree',full.names = TRUE)
counter = 0 
for (tree in trees) { counter = counter + 1;
iqtree = read.iqtree(tree)
iqtr = midpoint.root(as.phylo(iqtree))
lab = gsub('.SNP.*','',gsub('.*chr_','',tree))
cat ('Making tree: ',lab,'\n')

# Haplogroups 
circ = ggtree(iqtr, layout = "circular",branch.length='none') %<+% md  +
  geom_nodepoint(mapping=aes(subset=(as.numeric(label) >= 95)),col='black',fill='grey90',pch=23,size=0.75,show.legend=F)+
  geom_tippoint(aes(fill = Haplogroup,shape=Species),size=1)+
  scale_fill_manual('Haplogroup',values=md$HaplogroupColor,breaks=md$Haplogroup)+
  scale_shape_manual(values=md$Shape,breaks=md$Species)+
  guides(fill=guide_legend(nrow=4,override.aes=list(shape=21)),
         shape=guide_legend(nrow=5))+
  ggtitle(lab)+
  theme(legend.text = element_text(size = 8),legend.title = element_text(size = 10),
        legend.key.size = unit(0.2, "cm"),    legend.position = 'top')
circ

apetree = ggtree(iqtr, layout = "ape") %<+% md  +
  geom_nodepoint(mapping=aes(subset=(as.numeric(label) >= 95)),col='black',fill='grey90',pch=23,size=0.75,show.legend=F)+
  geom_tippoint(aes(fill = Haplogroup,shape=Species),size=1)+
  scale_fill_manual('Haplogroup',values=md$HaplogroupColor,breaks=md$Haplogroup)+
  scale_shape_manual(values=md$Shape,breaks=md$Species)+
  guides(fill=guide_legend(nrow=4,override.aes=list(shape=21)),
         shape=guide_legend(nrow=5))+
  theme(legend.text = element_text(size = 8),legend.title = element_text(size = 10),
        legend.key.size = unit(0.2, "cm"),    legend.position = 'top')
apetree


# Egg
circ_egg = ggtree(iqtr, layout = "circular",branch.length='none') %<+% md  +
  geom_nodepoint(mapping=aes(subset=(as.numeric(label) >= 95)),col='black',fill='grey90',pch=23,size=0.75,show.legend=F)+
  geom_tippoint(aes(fill = Egg,shape=Egg),size=1)+
  scale_fill_manual('Egg',values=md$EggCol,breaks=md$Egg)+
  scale_shape_manual(values=md$EggShape,breaks=md$Egg)+
  ggtitle(lab)+
  theme(legend.text = element_text(size = 8),legend.title = element_text(size = 10),
        legend.key.size = unit(0.2, "cm"),    legend.position = 'top')
circ_egg

apetree_egg = ggtree(iqtr, layout = "ape") %<+% md  +
  geom_nodepoint(mapping=aes(subset=(as.numeric(label) >= 95)),col='black',fill='grey90',pch=23,size=0.75,show.legend=F)+
  geom_tippoint(aes(fill = Egg,shape=Egg),size=1)+
  scale_fill_manual('Egg',values=md$EggCol,breaks=md$Egg)+
  scale_shape_manual(values=md$EggShape,breaks=md$Egg)+
  guides(fill=guide_legend(nrow=4,override.aes=list(shape=21)),
         shape=guide_legend(nrow=5))+
  theme(legend.text = element_text(size = 8),legend.title = element_text(size = 10),
        legend.key.size = unit(0.2, "cm"),    legend.position = 'top')
apetree_egg
ca = ggarrange(circ,apetree,nrow=2,heights=c(0.6,0.4),legend = 'none') 
ca_egg = ggarrange(circ_egg,apetree_egg,nrow=2,heights=c(0.6,0.4),legend = 'none') 
assign(paste0('p',counter),ca)
assign(paste0('e',counter),ca_egg)

}

ggarrange(e1,e2,p1,p2)
ggsave('~/symlinks/host/figures/20250330_MT-W-Trees-All.pdf',
       ggarrange(e1,e2,p1,p2),
       height=9,width=6,dpi=600)

```



### Plot Collapsed Tree [Fig 1A]

Collapse chrW tree on supported nodes for visualization (fig 1)

```R
#### Collapse chrW Tree 
setwd('~/EvoBioWolf/CUCKOO_gentes/phylogenetics')
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

md = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt')

##### COLLAPSE W TREE #####
w = read.iqtree('trees/chr_W.SNP.DP3-AC1-MQ40-MM1.min4.phy.contree')
w <- root(as.phylo(w),outgroup = '387_CP_MBW_RUS_F')
outs <- w$tip.label[grepl('_CM_|_CP_',w$tip.label)]
w_tree <- drop.tip(w,outs)

p = ggtree(w_tree, layout = "dendrogram") %>% ggtree::rotate(166) %>% ggtree::rotate(168) %>% ggtree::rotate(123) %<+% md 
inspect_tree <- p + geom_tippoint(aes(fill = Haplogroup),pch=21,size=4)+scale_fill_manual(values=md$HaplogroupColor,breaks=md$Haplogroup) + 
  geom_nodelab(aes(label=node),geom = 'label',size=2) 
ggsave('~/symlinks/host/figures/20250319_W_TREE-COLLAPSE_Inspection.pdf',
       inspect_tree,height=15,width=15,dpi=300,limitsize = FALSE)

hapcols = md %>% select(Haplogroup,HaplogroupColor) %>% unique %>% arrange(Haplogroup) %>% filter(!grepl('CM|CP',Haplogroup)) %>% na.omit

#collapse
p = ggtree(w_tree, layout = "dendrogram") %>% 
  ggtree::rotate(166) %>% ggtree::rotate(168) %>% ggtree::rotate(123) %<+% md + 
  geom_tippoint(aes(fill=Haplogroup,shape=SpeciesShort),size=3) +
  scale_fill_manual(values=md$HaplogroupColor,breaks=md$Haplogroup) + 
  scale_shape_manual(values=md$Shape,breaks=md$SpeciesShort)+
  geom_nodepoint(mapping=aes(subset=(as.numeric(label) >= 95)),col='black',fill='grey90',pch=23,size=1.5,show.legend=F)
#MCC1
p2 = p %>% collapse(node=117) + geom_point2(aes(subset=(node==117)), shape=21, size=3, fill=hapcols$HaplogroupColor[hapcols$Haplogroup == 'MCC1']); p2
#MCC2
p2 = collapse(p2,node=91) + geom_point2(aes(subset=(node==91)), shape=21, size=3, fill=hapcols$HaplogroupColor[hapcols$Haplogroup == 'MCC2']); p2
#MCC3
p2 = collapse(p2,node=99) + geom_point2(aes(subset=(node==99)), shape=21, size=3, fill=hapcols$HaplogroupColor[hapcols$Haplogroup == 'MCC3']); p2
#BIN MCC4 MCO1 and MCO2 (rufous)
p2 = collapse(p2,node=167) + geom_point2(aes(subset=(node==167)), shape=21, size=3, fill=hapcols$HaplogroupColor[hapcols$Haplogroup == 'MCC4']); p2
#MCC5
p2 = collapse(p2,node=162) + geom_point2(aes(subset=(node==162)), shape=21, size=3, fill=hapcols$HaplogroupColor[hapcols$Haplogroup == 'MCC5']); p2
#MCC6
p2 = collapse(p2,node=160) + geom_point2(aes(subset=(node==160)), shape=21, size=3, fill=hapcols$HaplogroupColor[hapcols$Haplogroup == 'MCC6']); p2
#MCC7
p2 = collapse(p2,node=125) + geom_point2(aes(subset=(node==125)), shape=21, size=3, fill=hapcols$HaplogroupColor[hapcols$Haplogroup == 'MCC7']); p2
#MCO3
p2 = collapse(p2,node=143) + geom_point2(aes(subset=(node==143)), shape=24, size=3, fill=hapcols$HaplogroupColor[hapcols$Haplogroup == 'MCO3']); p2
#MCO4
p2 = collapse(p2,node=145) + geom_point2(aes(subset=(node==145)), shape=24, size=3, fill=hapcols$HaplogroupColor[hapcols$Haplogroup == 'MCO4']); p2
p2 = p2 + geom_treescale(x=0.05,offset = 0.01) 

ggsave('~/symlinks/host/figures/20250319_W_TREE-COLLAPSED.pdf',
       p2,height=1,width=6,dpi=300)

```

### Plot Dendrogram with Phenotype Bars [Fig. 2A]

```bash
#### Plot Trees, Assign Haplogroups 
setwd('~/EvoBioWolf/CUCKOO_gentes/phylogenetics/')
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
library(rcartocolor)

md = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt') 

#only grab samples with known egg (hash out filter for related individuals)
md_egg = md %>%
  filter(Analysis_FullGensAssociations == 1) %>%
  select(ID = ID, Egg, Haplogroup, Geography = GeographicGroup,Ancestry = AncestryA5)
md_egg %>% count(Egg)

md_egg = md_egg %>% group_by(Egg) %>% mutate(TotalEgg = n()) %>% ungroup
minobs=2
md_egg = md_egg %>% filter(TotalEgg >= minobs)
md_egg %>% count(Egg)

## Plot Trees
## mtDNA
mt = read.iqtree('trees/chr_MT.SNP.DP3-AC1-MQ40-MM1.min4.phy.contree')
mt <- midpoint.root(as.phylo(mt))

# voila 
mtp = ggtree(mt, layout = "dendrogram") %<+% md  + 
  #geom_nodelab(aes(label=node),geom = 'label',size=3)+  
  geom_nodepoint(mapping=aes(subset=(as.numeric(label) >= 95)),col='black',fill='grey90',pch=23,size=1.5,show.legend=F)+
  geom_tippoint(aes(fill = Haplogroup,shape=SpeciesShort),size=2)+
  scale_shape_manual(values=c(21,24))+
  scale_fill_manual(values=md$HaplogroupColor,breaks=md$Haplogroup)+
  #scale_fill_manual('Egg',values=md$HostCol,breaks=md$Egg)+
  #geom_tiplab()+
  guides(fill=guide_legend(nrow=2,override.aes=list(shape=21)))+
  theme(legend.position='top')
mtp

## W chromosome
w = read.iqtree('trees/chr_W.SNP.DP3-AC1-MQ40-MM1.min4.phy.contree')
w <- midpoint.root(as.phylo(w))

#plot cladogram, first find nodes to collapse 
wp = ggtree(w, layout = "dendrogram") %<+% md  + 
  geom_nodepoint(mapping=aes(subset=(as.numeric(label) >= 95)),col='black',fill='grey90',pch=23,size=1.5,show.legend=F)+
  geom_tippoint(aes(fill = Haplogroup,shape=SpeciesShort),size=2)+
  scale_shape_manual(values=c(21,24))+
  scale_fill_manual(values=md$HaplogroupColor,breaks=md$Haplogroup)+
  guides(fill=guide_legend(nrow=2,override.aes=list(shape=21)))+
  theme(legend.position='top')
wp

ggsave('~/symlinks/host/figures/20250330_N136-All-MT-W.pdf',
       ggarrange(mtp,wp,nrow=2,common.legend = TRUE),width=6,height=5,dpi=300)


#### W tree with egg: ancestry+geography+haplogroup
md <- md %>% arrange(Haplogroup,HaplogroupColor)
p = ggtree(w, layout = 'dendrogram') %<+% md + 
  geom_nodepoint(mapping=aes(subset=(as.numeric(label) >= 95)),col='black',fill='grey90',pch=23,size=0.75,show.legend=F)+
  geom_tippoint(aes(fill=Haplogroup,shape=SpeciesShort),size=1.25)+
  geom_tiplab(align=TRUE,size=1)+
  scale_shape_manual(values=c(21,24))+
  scale_fill_manual(values=md$HaplogroupColor,breaks=md$Haplogroup)+
  geom_treescale()+
  guides(fill=guide_legend(override.aes=list(shape=21)))
p

# Alt, with egg shape
p = ggtree(w, layout = 'dendrogram',branch.length = 'none') %<+% md + 
  geom_nodepoint(mapping=aes(subset=(as.numeric(label) >= 95)),col='black',fill='grey90',pch=23,size=0.75,show.legend=F)+
  geom_tippoint(aes(fill=Egg,shape=Egg),size=1.25)+
  geom_tiplab(align=TRUE,size=1)+
  scale_shape_manual(values=md$EggShape,breaks=md$Egg)+
  scale_fill_manual(values=md$EggCol,breaks=md$Egg)+
  geom_treescale()
p

phenos = as.data.frame(p$data %>% filter(isTip == TRUE))
rownames(phenos) = phenos$label

#second response variable: egg
p1 = p + new_scale_fill()
egg = phenos %>% select(E = Egg)
md <- md %>% arrange(EggOrder)
p1 = gheatmap(p1, egg,offset=.015,width=0.05) +
  scale_fill_manual('E',values=md$EggCol,breaks=md$Egg,na.value='white')
p1

#third variable: ancestry
p2 = p1 + new_scale_fill()
anc = phenos %>% select(N = AncestryA5)
md <- md %>% arrange(AncestryA5)
p2 = gheatmap(p2, anc,offset=.03,width=0.05) +
  scale_fill_manual('N',values=md$AncestryColor,breaks=md$AncestryA5)
p2

# fourth: geography
p3 = p2 + new_scale_fill()
md <- md %>% arrange(GeographicGroup)
geo = phenos %>% select(G = GeographicGroup)
p3 = gheatmap(p3, geo,offset=.045,width=0.05) +
  scale_fill_manual('G',values=md$GeoColor,breaks=md$GeographicGroup)
p3

# # Save first as haplogroup
# ggsave('~/symlinks/host/figures/20250330_Dendrogram-W-Eggs-Phenotypes-HAPLOGROUP.pdf',
#        p3,
#        height=5,width=8,dpi=300)

# And then save with egg 
ggsave('~/symlinks/host/figures/20250330_Dendrogram-W-Eggs-Phenotypes-EGG.pdf',
       p3,
       height=5,width=8,dpi=300)

```

#### Zoom MCC7

```
#### Plot Trees, Assign Haplogroups 
setwd('~/EvoBioWolf/CUCKOO_gentes/phylogenetics/')
.libPaths('~/symlinks/cuck00/software/mambaforge/envs/r25/lib/R/library')
library(ggtree)
library(ape)
library(tidyverse)
library(treeio)
library(viridis)
#library(adegenet)
library(ggpubr)
#library(vcfR)
library(RColorBrewer)
library(phytools)
library(ggnewscale)
library(rcartocolor)

md = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt') 
leg <- md %>% select(Haplogroup,Shape,HaplogroupColor) %>% na.omit %>% distinct
eggleg <- md %>% select(Egg,Shape,EggCol) %>% na.omit %>% distinct

#only grab samples with known egg (hash out filter for related individuals)
md_egg = md %>%
  filter(Analysis_GensAssociations == 1) %>%
  select(ID = ID, Egg, Haplogroup, Geography = GeographicGroup,Ancestry = AncestryA5)
md_egg %>% count(Egg)

md_egg = md_egg %>% group_by(Egg) %>% mutate(TotalEgg = n()) %>% ungroup
minobs=2
md_egg = md_egg %>% filter(TotalEgg >= minobs)
md_egg %>% count(Egg)

## Plot Trees
## mtDNA
mt = read.iqtree('trees/chr_MT.SNP.DP3-AC1-MQ40-MM1.min4.phy.contree')
mt <- midpoint.root(as.phylo(mt))

# voila 
mtp = ggtree(mt, layout = "dendrogram") %<+% md  + 
  #geom_nodelab(aes(label=node),geom = 'label',size=3)+  
  geom_nodepoint(mapping=aes(subset=(as.numeric(label) >= 95)),col='black',fill='grey90',pch=23,size=1.5,show.legend=F)+
  geom_tippoint(aes(fill = Haplogroup,shape=SpeciesShort),size=2)+
  scale_shape_manual(values=c(21,24))+
  scale_fill_manual(values=leg$HaplogroupColor,breaks=leg$Haplogroup)+
  #scale_fill_manual('Egg',values=md$HostCol,breaks=md$Egg)+
  #geom_tiplab()+
  guides(fill=guide_legend(nrow=2,override.aes=list(shape=21)))+
  theme(legend.position='top')
mtp

# MCC7 tree for mtDNA 
m7m <- ggtree(mt, layout = "dendrogram") %<+% md
ids <- m7m$data %>% filter(isTip == TRUE & grepl('MCC7|MCO3',Haplogroup)) %>% pull(label)
node <- getMRCA(as.phylo(m7m), ids)
m7mt <- tidytree::tree_subset(as.phylo(m7m), node = node, levels_back = 0)
m7mtz <- ggtree(m7mt, layout = "dendrogram") %<+% md +    
  geom_nodepoint(mapping=aes(subset=(as.numeric(label) >= 95)),col='black',fill='grey90',pch=23,size=1.5,show.legend=F)+
  geom_tippoint(aes(fill = Egg,shape=SpeciesShort),size=2)+
  scale_shape_manual(values=c(21,24))+
  scale_fill_manual(values=eggleg$EggCol,breaks=eggleg$Egg)+
  guides(fill=guide_legend(nrow=2,override.aes=list(shape=21)))+
  theme(legend.position='top')
m7mtz

ggsave('~/symlinks/host/figures/20250729_MT-Tree-MCC7-Zoom.pdf',m7mtz,width=4,height=2,dpi=300)

## W chromosome
w = read.iqtree('trees/chr_W.SNP.DP3-AC1-MQ40-MM1.min4.phy.contree')
w <- midpoint.root(as.phylo(w))

#plot cladogram, first find nodes to collapse 
wp = ggtree(w, layout = "dendrogram") %<+% md  + 
  geom_nodepoint(mapping=aes(subset=(as.numeric(label) >= 95)),col='black',fill='grey90',pch=23,size=1.5,show.legend=F)+
  geom_tippoint(aes(fill = Haplogroup,shape=SpeciesShort),size=2)+
  scale_shape_manual(values=c(21,24))+
  scale_fill_manual(values=leg$HaplogroupColor,breaks=leg$Haplogroup)+
  guides(fill=guide_legend(nrow=2,override.aes=list(shape=21)))+
  theme(legend.position='top')
wp

ggsave('~/symlinks/host/figures/20250330_N136-All-MT-W.pdf',
       ggarrange(mtp,wp,nrow=2,common.legend = TRUE),width=6,height=5,dpi=300)

#focus on MCC7 
m7 <- ggtree(w, layout = "dendrogram") %<+% md
ids <- m7$data %>% filter(isTip == TRUE & grepl('MCC7|MCO3',Haplogroup)) %>% pull(label)
node <- getMRCA(as.phylo(m7), ids)
m7t <- tidytree::tree_subset(as.phylo(m7), node = node, levels_back = 0)
m7tz <- ggtree(m7t, layout = "dendrogram") %<+% md +    
  geom_nodepoint(mapping=aes(subset=(as.numeric(label) >= 95)),col='black',fill='grey90',pch=23,size=1.5,show.legend=F)+
  geom_tippoint(aes(fill = Egg,shape=SpeciesShort),size=2)+
  scale_shape_manual(values=c(21,24))+
  scale_fill_manual(values=eggleg$EggCol,breaks=eggleg$Egg)+
  guides(fill=guide_legend(nrow=2,override.aes=list(shape=21)))+
  theme(legend.position='top')
m7tz

ggsave('~/symlinks/host/figures/20250729_W-Tree-MCC7-Zoom.pdf',m7tz,width=4,height=2,dpi=300)


#### W tree with egg: ancestry+geography+haplogroup
md <- md %>% arrange(Haplogroup,HaplogroupColor)
p = ggtree(w, layout = 'dendrogram') %<+% md + 
  geom_nodepoint(mapping=aes(subset=(as.numeric(label) >= 95)),col='black',fill='grey90',pch=23,size=0.75,show.legend=F)+
  geom_tippoint(aes(fill=Haplogroup,shape=SpeciesShort),size=1.25)+
  geom_tiplab(align=TRUE,size=1)+
  scale_shape_manual(values=c(21,24))+
  scale_fill_manual(values=md$HaplogroupColor,breaks=md$Haplogroup)+
  geom_treescale()+
  guides(fill=guide_legend(override.aes=list(shape=21)))
p

# Alt, with egg shape
p = ggtree(w, layout = 'dendrogram',branch.length = 'none') %<+% md + 
  geom_nodepoint(mapping=aes(subset=(as.numeric(label) >= 95)),col='black',fill='grey90',pch=23,size=0.75,show.legend=F)+
  geom_tippoint(aes(fill=Egg,shape=Egg),size=1.25)+
  geom_tiplab(align=TRUE,size=1)+
  scale_shape_manual(values=md$EggShape,breaks=md$Egg)+
  scale_fill_manual(values=md$EggCol,breaks=md$Egg)+
  geom_treescale()
p

phenos = as.data.frame(p$data %>% filter(isTip == TRUE))
rownames(phenos) = phenos$label

#second response variable: egg
p1 = p + new_scale_fill()
egg = phenos %>% select(E = Egg)
md <- md %>% arrange(EggOrder)
p1 = gheatmap(p1, egg,offset=.015,width=0.05) +
  scale_fill_manual('E',values=md$EggCol,breaks=md$Egg,na.value='white')
p1

#third variable: ancestry
p2 = p1 + new_scale_fill()
anc = phenos %>% select(N = AncestryA5)
md <- md %>% arrange(AncestryA5)
p2 = gheatmap(p2, anc,offset=.03,width=0.05) +
  scale_fill_manual('N',values=md$AncestryColor,breaks=md$AncestryA5)
p2

# fourth: geography
p3 = p2 + new_scale_fill()
md <- md %>% arrange(GeographicGroup)
geo = phenos %>% select(G = GeographicGroup)
p3 = gheatmap(p3, geo,offset=.045,width=0.05) +
  scale_fill_manual('G',values=md$GeoColor,breaks=md$GeographicGroup)
p3

# # Save first as haplogroup
# ggsave('~/symlinks/host/figures/20250330_Dendrogram-W-Eggs-Phenotypes-HAPLOGROUP.pdf',
#        p3,
#        height=5,width=8,dpi=300)

# And then save with egg 
ggsave('~/symlinks/host/figures/20250330_Dendrogram-W-Eggs-Phenotypes-EGG.pdf',
       p3,
       height=5,width=8,dpi=300)

```


