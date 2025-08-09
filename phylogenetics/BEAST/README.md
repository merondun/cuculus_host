## BEAST+Rooting: Using 929-Kb

For BEAST, I will use the same genic subset used in the plumage paper (Merondun 2024), spanning 929-Kb. 

Subset samples:

```bash
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/phylogenetics/202508_Genes10k-75k')
library(tidyverse)

#Read in metadata
md = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt') %>% 
  filter(Sex == 'F' & Analysis_PopulationGenetics == 1 & Egg != 'NA')
ogs <- read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt') %>% filter(grepl('CM|CP', SpeciesShort) & Sex == 'F')
rep1 <- rbind(md,ogs) %>% group_by(Haplogroup) %>% slice_sample(n=1) %>% pull(ID)
rep2 <- rbind(md,ogs) %>% group_by(Haplogroup) %>% slice_sample(n=1) %>% pull(ID)
write.table(rep1,file='Rep1.list',quote=F,sep='\t',row.names=F,col.names=F)
write.table(rep2,file='Rep2.list',quote=F,sep='\t',row.names=F,col.names=F)

```

Grab gene regions:

```bash
grep -wf Retained_Genes_Between10k-75k_BSP.txt genes.gff | sed 's/;.*//g' | sed 's/ID=gene-//g' | awk '{OFS="\t"}{print $1, $4, $5, $9, ".", $7}' > Genes.bed

for SAMPLE in $(cat AllSamples.list); do 
echo "Extracting ${SAMPLE}"
# extract genes
sed '1s/.*/>chr_W/g' consensus_refs/${SAMPLE}.chr_W.All.Mask.fa > consensus_refs/${SAMPLE}.chr_W.All.Mask.tmp.fa
bedtools getfasta -s -fi consensus_refs/${SAMPLE}.chr_W.All.Mask.tmp.fa -fo consensus_refs/${SAMPLE}.chr_W.Genes.Mask.fa -bed Genes.bed -name 
awk '!/^>/ {printf "%s", $0} END {print ""}' consensus_refs/${SAMPLE}.chr_W.Genes.Mask.fa > consensus_refs/${SAMPLE}.chr_W.Genes-Merged.Mask.fa
echo ">${SAMPLE}" | cat - consensus_refs/${SAMPLE}.chr_W.Genes-Merged.Mask.fa > consensus_refs/${SAMPLE}.chr_W.Genes.Mask.fa
sed -i "1s/.*/>${SAMPLE}/" consensus_refs/${SAMPLE}.chr_W.Genes.Mask.fa

done 

# Trim gaps 
samtools faidx chr_W.Genes.Mask.fa $(cat AllEggs.list) > chr_W.Genes.Mask.N87.fa
sed 's/N/-/g; s/n/-/g' chr_W.Genes.Mask.N87.fa > chr_W.Genes.Mask.N87.tmp.fa
trimal -in chr_W.Genes.Mask.N87.tmp.fa -out chr_W.Genes.Mask.N87.tmp2.fa -gt 0.8
sed 's/-/N/g' chr_W.Genes.Mask.N87.tmp2.fa > chr_W.Genes.Mask.TrimalGT8.N87.Ns.fa

# Subset
samtools faidx chr_W.Genes.Mask.TrimalGT8.N87.Ns.fa
samtools faidx chr_W.Genes.Mask.TrimalGT8.N87.Ns.fa $(cat Rep1.list) > rep1.fa
samtools faidx chr_W.Genes.Mask.TrimalGT8.N87.Ns.fa $(cat Rep2.list) > rep2.fa

# prep for beast
seqret -sequence rep1.fa -outseq rep1.nex -osformat nexus
seqret -sequence rep2.fa -outseq rep2.nex -osformat nexus
```

### Root Digger

Root digger on the full tree, and each replicate:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=2-00:00:00

# mamba activate phylos

THREADS=20
FASTA=$1
SEQ=$(basename "$FASTA" .fa)
OUTDIR="roots/${SEQ}"
mkdir -p "$OUTDIR"
grep '>' ${FASTA} | sed 's/>//g' > ${OUTDIR}/outgroups.list

# make tree
MODEL="GTR"

iqtree -s ${FASTA} -m "$MODEL" -pre ${OUTDIR}/${SEQ} \
|| { iqtree -s ${OUTDIR}/${SEQ}.varsites.phy -m "$MODEL" -pre ${OUTDIR}/${SEQ}.varsites.phy && \
mv ${OUTDIR}/${SEQ}.varsites.phy.treefile ${OUTDIR}/${SEQ}.treefile; }

# also use rootdigger
apptainer exec -B "$(pwd)":/work --pwd /work root_digger_v1.7.0.sif \
rd --msa "${FASTA}" --tree "${OUTDIR}/${SEQ}.treefile" --exhaustive
```

Plot: 

```bash
#### Plot root digger trees
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/phylogenetics/202508_Genes10k-75k/')
#.libPaths('~/symlinks/cuck00/software/mambaforge/envs/r/lib/R/library')
.libPaths('~/r_libs/')
library(ggplot2)
library(ggtree)
library(tidyverse)
library(ape)
library(treeio)
library(viridis)
library(ggpubr)
library(RColorBrewer)

md = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt')
leg <- md %>% select(Haplogroup,HaplogroupColor,Species,Shape) %>% na.omit %>% distinct

#plot tree
trees = list.files('trees_root_digger',pattern='.*tree',full.names = TRUE)
counter = 0 
for (tree in trees) { counter = counter + 1;
root_tree = read.tree(tree)
lab = gsub('trees_root_digger/','',gsub('.treefile.rooted.tree','',tree))
cat ('Making tree: ',lab,'\n')

t <- ggtree(root_tree, layout = "rectangular") %<+% md + 
  geom_tippoint(aes(fill = Haplogroup,shape=Species),size=1.5)+
  scale_fill_manual('Haplogroup',values=leg$HaplogroupColor,breaks=leg$Haplogroup)+
  scale_shape_manual(values=leg$Shape,breaks=leg$Species)+
  guides(fill=guide_legend(nrow=4,override.aes=list(shape=21)),
         shape=guide_legend(nrow=5))+
  ggtitle(lab)+
  theme(legend.text = element_text(size = 8),legend.title = element_text(size = 10),
        legend.key.size = unit(0.2, "cm"),    legend.position = 'top')
assign(paste0('p',counter),t)

}
tp <- ggarrange(p2,p3,p1,nrow=3,ncol=1,heights=c(0.25,0.25,0.5),common.legend = TRUE)
tp
ggsave('~/symlinks/host/figures/20250804_Root_digger_Trees.pdf',tp,height=9,width=4.5)

```

### Beast

Prep 929kb alignment with 1 subsampled individual per replicate in beauti:

````bash
#for beauti
4 gamma categories, HKY, empirical sub rate 
clock rate strict at 5.05E-9
coalescent exponential, priors with lognormal CM / CP ancestor at M = 1.4  S=0.15
set prior, NonBlueCP on monophyletic MRCA all except MCC1-MCC3 and CP
set prior, NonCP, monophyletic everything except CP 
10M chains, log every 1k 

# to scale for substitutions/site/Ma: 1.843e-3
#   Rate per Ma= 5.05×10−9/2.74 * 10^6 =1.843*10^−3 substitutions/site/Ma
````

Run:

```
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=10
#SBATCH --time=2-00:00:00

# mamba activate beast

RUN=$1
beast -threads 10 -overwrite -beagle_SSE -seed 777 -java ${RUN}.xml

```

Plot:

```bash
#### Plot BEAST annotated trees 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/phylogenetics/202503_beast/m1.4_s.15/final_trees/')
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

pdf('~/symlinks/host/figures/20250501_BEAST_Divergence_Dating_All-Dual.pdf',height=9,width=7)
ggarrange(p1,p2,p3,p4,common.legend = TRUE)
dev.off()

```



