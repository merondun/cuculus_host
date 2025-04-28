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



### Plot Collapsed Tree

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

### Plot Dendrogram with Phenotype Bars

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



## TWISST

To go deeper than C. poliocephalus, download Clamator glandarius female data (HiFi) from the interwebs. Align it to the genome, and call SNPs (only on the W for brevity).

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=20
#SBATCH --time=48:00:00

WD=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/phylogenetics/202503_twisst/deeper
READS=/dss/dsslegfs01/pr53da/pr53da-dss-0021/rawdata/External_Data_SRA_GenomeArk/SRA_Cuculiformes_Outgroups/Clamator_glandarius_SRR26807982_1.fastq.gz
GENOME=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa

mkdir -p coverage

# Align reads
pbmm2 align \
    --preset HIFI \
    --sort \
    --num-threads 20 \
    --rg "@RG\tID:1\tSM:999_CG_UNK_UNK_F\tPL:HiFi\tLB:999_CG_UNK_UNK_F" \
    ${GENOME} ${READS} Clamator_glandarius.bam
samtools index Clamator_glandarius.bam

# call SNPs on chrW
bcftools mpileup --threads 20 -Ou -r chr_W -f ${GENOME} Clamator_glandarius.bam | \
	bcftools call --threads 20 -mv -Oz -o Clamator_glandarius.chrW.vcf.gz
bcftools index Clamator_glandarius.chrW.vcf.gz

# coverage
mosdepth --threads 20 --mapq 30 --by 1000000 --fast-mode --no-per-base coverage/Clamator_glandarius Clamator_glandarius.bam
```

Prepare and run:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=10
#SBATCH --time=200:00:00

#mamba activate twisst

WD=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/phylogenetics/202503_twisst
#VCF=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/phylogenetics/202503_snapp/vcfs/chr_W.SNP.DP3-AC1-MQ40-MM1.vcf.gz
VCF=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/phylogenetics/202503_twisst/deeper/vcfs/chr_W.SNP.DP3-AC2-MQ40-MM1.vcf.gz

cd ${WD}

mkdir -p input output

# Force haploid
bcftools view ${VCF} | \
  bcftools +fixploidy -Ou - -- -f 1 | \
  bcftools view -Oz -o input/chr_W.SNP.DP3-AC1-MQ40-MM1-HAPLOID.vcf.gz

# Create geno file
python ~/modules/genomics_general/VCF_processing/parseVCF.py --ploidyMismatchToMissing --ploidy 1 --skipIndels -i input/chr_W.SNP.DP3-AC1-MQ40-MM1-HAPLOID.vcf.gz | \
        bgzip > input/chr_W.SNP.DP3-AC1-MQ40-MM1-HAPLOID.geno.gz

#only retain SNPs with MAC >= 2
python ~/modules/genomics_general/filterGenotypes.py --threads 10 -i input/chr_W.SNP.DP3-AC1-MQ40-MM1-HAPLOID.geno.gz -o input/chr_W.SNP.DP3-AC2-MQ40-MM1-HAPLOID.geno.gz --minVarCount 2 --minAlleles 2 --maxAlleles 2

for win in 50 100 500 1000; do
  miss=$(echo ${win} | awk -v w=${win} '{print w/5}')
  missind=$(echo ${win} | awk -v w=${win} '{print w/10}')

  #calculate trees in windows
  python ~/modules/genomics_general/phylo/phyml_sliding_windows.py -T 10 \
          -g input/chr_W.SNP.DP3-AC2-MQ40-MM1-HAPLOID.geno.gz --prefix output/chr_W.SNP.DP3-AC2-MQ40-MM1-HAPLOID.phyml_bionj.Sw${win}-m${miss}-ind${missind} \
         --windSize ${win} --minSites ${miss} --minPerInd ${missind} --windType sites --model GTR --optimise n

  #run twisst
  ~/modules/twisst/twisst.py -t output/chr_W.SNP.DP3-AC2-MQ40-MM1-HAPLOID.phyml_bionj.Sw${win}-m${miss}-ind${missind}.trees.gz \
    -w output/chr_W.SNP.DP3-AC2-MQ40-MM1-HAPLOID.weights.Sw${win}-m${miss}-ind${missind}.csv.gz --outgroup CG \
    --outputTopos output/chr_W.SNP.DP3-AC2-MQ40-MM1-HAPLOID.topologies.Sw${win}-m${miss}-ind${missind}.trees -g CG -g CP -g CM -g BLUE -g RUFOUS -g OTHER --method complete --groupsFile Females_N94-Simple.pop

done

```

### Plot

```R
#### Plot twisst
setwd('~/EvoBioWolf/CUCKOO_gentes/phylogenetics/twisst/')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(viridis)
library(RColorBrewer)
library(ggpubr)
library(meRo)
library(ggdist)
source('~/modules/twisst/plot_twisst.R')

# Import data 
weights_files = list.files('.',pattern='.*weights.*.csv.gz')
window_data_files = list.files('.',pattern='.*data.tsv')

windat <- list()
for (size in seq(1,4,1)) {
# Read in the data
win <- gsub('-.*','',gsub('.*weights.Sw','',weights_files[size]))
twisst_data <- import.twisst(weights_files=weights_files[size],window_data_files=window_data_files[size])
data <- twisst_data$weights_raw$chr_W
data$win <- as.numeric(win)
data <- data %>% pivot_longer(!win)
windat[[size]] <- data
}
full_dat <- rbindlist(windat) %>% as_tibble

topord <- full_dat %>% select(name) %>% distinct %>% mutate(ord=as.numeric(gsub('topo','',name))) %>% arrange(ord)
full_dat$name <- factor(full_dat$name,levels=topord$name)

# Summarize topos
topo_sums <- full_dat %>% 
  group_by(name,win) %>% 
  sum_stats(value)
top_topos <-  topo_sums %>% ungroup  %>%  filter(win == 500) %>% slice_max(mean,n=10) 
topo_sums$name <- factor(topo_sums$name,levels=top_topos$name)
extract_topos <- top_topos %>% head(n=10) %>% pull(name)
plot <- topo_sums %>% 
  filter(name %in% extract_topos) %>% 
  ggplot(aes(x=name,y=mean,ymin=conf_low,ymax=conf_high,col=name))+
  geom_point()+
  facet_grid(win~.,scales='free')+
  scale_color_manual(values=topo_cols)+
  geom_errorbar()+
  theme_bw(base_size=6)+
  coord_cartesian(ylim=c(0,1))
plot

ggsave('~/symlinks/host/figures/20250325_TWISST_chrW_Sensitivity.pdf',
       plot, height = 5,width=7,dpi=300)

##### Plot topologies using Ape #####
pdf('~/symlinks/host/figures/20250320_Twisst_Plotted_Topologies.pdf',height=1.25,width=15)

par(mfrow = c(1,length(twisst_data$topos)), mar = c(1,1,2,1), xpd=NA)
for (n in 1:length(twisst_data$topos)){
  plot.phylo(twisst_data$topos[[n]], type = "clad", edge.color=topo_cols[n], edge.width=3, label.offset=.2, cex = 1)
  mtext(side=3,text=paste0("topo",n))
}

dev.off()

topo_ids <- as.integer(gsub('topo','',extract_topos))
pdf('~/symlinks/host/figures/20250325_Twisst_Top_Topologies.pdf',height=1.25,width=15)

par(mfrow = c(1,length(topo_ids)), mar = c(1,1,2,1), xpd=NA)
for (i in seq_along(topo_ids)) {
  t <- topo_ids[i]
  plot.phylo(twisst_data$topos[[t]], type = "clad", edge.color=topo_cols[i], edge.width=3, label.offset=.2, cex = 1)
  mtext(side=3, text=paste0("topo", i))
}

dev.off()

```

## BEAST

Following this [tutorial](https://beast2-dev.github.io/beast-docs/beast2/DivergenceDating/DivergenceDatingTutorial.html).

From the host-type fastas:

```bash
seqret -sequence ml_trees/chr_MT_HostOG.SNP.DP3-AC2-MQ40.min4.fasta -outseq nexus/chr_MT_HostOG.SNP.DP3-AC2-MQ40.nex -osformat nexus
seqret -sequence ml_trees/chr_W_HostOG.SNP.DP3-AC2-MQ40.min4.fasta -outseq nexus/chr_W_HostOG.SNP.DP3-AC2-MQ40.nex -osformat nexus
```

For BEAUTI:

```bash
#for beauti
4 gamma categories, HKY, empirical sub rate 
clock rate strict at 5.05E-9
coalescent exponential, priors with lognormal CM / CP ancestor at M = 1.7  S=0.2
# with CG: 3.5 0.25 gives 95% 20.3 - 54.1
100M chains, log every 1k 
```

Run BEAST & Annotate:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=20
#SBATCH --time=200:00:00

RUN=$1
beast -threads 20 -overwrite -java ${RUN}.xml
```

For adding invariant sites: https://groups.google.com/g/beast-users/c/QfBHMOqImFE

Count AGCT:

```bash
grep -v "^>" chr_W.fa | fold -w1 | sort | uniq -c

chr_W:
T	G	C	A
6108268	4562746	4578162	6088177

chr_MT:
T	G	C	A
4749	2471	5775	6703

e.g. 
   <data id='chr_W_HostOG' spec='FilteredAlignment' filter='-' data='@chr_W_HostOGOriginal' constantSiteWeights='6088177 4578162 4562746 6108268'/>
```

Info from the BEAST google groups help forum: 

```bash
You can use a FilteredAlignment to insert constant sites and set the constantSiteWeights attribute. Say, your original alignment is called xyz, so the XML produced by BEAUti contains something like

    <data id="xyz" name="alignment">

It is easiest to rename this to say xyzOriginal,

    <data id="xyzOriginal" name="alignment">

then add another data element, just after the closing </data> element of the alignment would be a good spot, that say

   <data id='xyz' spec='FilteredAlignment' filter='-' data='@xyzOriginal' constantSiteWeights='100 200 300 400'/>

Note id='xyz' and data='@xyzOriginal' should match what you have in the XML.

The constant weights at the end add weights for DNA in order A,C,G,T, so it adds 100 constant sites with all As, 200 with all Cs etc.

In the output to screen, it should report statistics of the xyzOriginal as something like:

6 taxa
768 sites
69 patterns

followed by statistics of the filtered alignment

Filter -
6 taxa
768 sites + 1000 constant sites
69 patterns

where the total number of constant sites added are reported as well.
```

Treeannotator:

```bash
for i in $(ls *trees); do 
/dss/dsshome1/lxc07/di39dux/modules/beast/bin/treeannotator -b 10 -height mean -file ${i} ${i}.ann
done 
```



### Plot

```bash
#### Plot BEAST annotated trees 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/phylogenetics/202503_beast/inputs')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(ggtree)
library(phytools)
library(ape)
library(treeio)
library(viridis)
library(ggpubr)
library(RColorBrewer)
library(tidyverse)

#Read in metadata
md = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt') %>% select(ID,Haplogroup,HaplogroupColor,SpeciesShort,Shape) 
files = list.files('.',paste0('.*ann'))

counter = 0
for (file in files){
  counter = counter +  1 
  iqtree = read.beast(file) 
  gg = ggtree(iqtree,layout='rectangular') %<+% md
  
  #add label for 95% CIs
  lab = gsub('.trees.*','',file)
  heights = gg$data$height_0.95_HPD
  df = as.data.frame(do.call(rbind, heights)) #convert the list to a data frame
  df$node_value = 1:nrow(df) # Add node values as a new column
  colnames(df) = c("value1", "value2", "node")
  df = df[, c("node", "value1", "value2")]
  df = df %>% mutate(lab = paste0(round(value1,1),' - ',round(value2,1))) %>% select(!c(value1,value2))
  
  leg = md %>% select(Haplogroup,HaplogroupColor,SpeciesShort,Shape) %>% unique %>% drop_na(HaplogroupColor)
  gg$data = left_join(gg$data,df)
  ggp = gg  +
    geom_range(range='height_0.95_HPD', color='red', alpha=.6, size=1) +
    geom_tippoint(aes(fill = Haplogroup,shape=SpeciesShort),size=1.5)+
    geom_nodelab(aes(label=lab),size=1.5,vjust=1) +
    ggtitle(lab)+
    #geom_tiplab(size=2)+
    #geom_nodelab(aes(x=branch, label=round(posterior, 2)), vjust=-.5, size=3) +
    scale_color_continuous(low="darkgreen", high="red") +
    scale_fill_manual(values=leg$HaplogroupColor,breaks=leg$Haplogroup)+
    scale_shape_manual(values=leg$Shape,breaks=leg$SpeciesShort)+
    theme(legend.position=c(.1, .8))+
    geom_treescale(x = 5)+
    guides(fill=guide_legend(override.aes=list(shape=21)))+
    theme(legend.position='right')
  ggp
  assign(paste0('p',counter),ggp)
} 

ggarrange(p1,p2,p3,p4,common.legend = TRUE)

pdf('~/symlinks/host/figures/20250424_BEAST_Divergence_Dating_All-Dual.pdf',height=9,width=7)
ggarrange(p1,p2,p3,p4,common.legend = TRUE)
dev.off()

```

