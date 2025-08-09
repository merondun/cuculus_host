## de novo Assembly mtDNA

Generate circular de novo mtDNA for all samples, first subsample 5gb of illumina seq from each using bbtools: 

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=200:00:00

RUN=$1
reference=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/mtDNA/GCA_017976375.1_bCucCan1.pri_mtDNA.fa
trimdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/trimmed_fastq
subbed=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/trimmed_fastq/5gb
mitodir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/denovo_mito/${RUN}
circdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/denovo_mito/2022_07-Circular

rm -rf $mitodir
mkdir $mitodir

## Generate de novo mitochondrion
zcat $trimdir/${RUN}*gz | bbduk.sh in=stdin.fq out=$subbed/${RUN}.fastq.gz int=t maxbasesout=5000000000
cd $mitodir
~/modules/MITObim/MITObim.pl -start 1 -end 50 -kbait 21 --clean -sample ${RUN} -ref VGP -readpool ${subbed}/${RUN}.fastq.gz --quick ${reference} &> log;

cd $mitodir

## Circularization
mkdir tmp
mkdir raw_mitobim
mkdir circular

#Determine the highest output iteration from MITOBIM
iteration=$(ls iteration* | sort -r --version-sort | grep 'iter' | head -n1 | sed 's/://g')

#Copy the NON-IUPAC fasta to a new directory
cp $iteration/${RUN}-VGP-it*noIUPAC.fasta raw_mitobim/${RUN}.uncirc.fasta

#Determine the longest contiguous circular assembly for each mitogenome
length=$(~/modules/MITObim/misc_scripts/circules.py -f raw_mitobim/${RUN}.uncirc.fasta -k 11 | grep 'suggested' | tail -n 1 | sed "s@.*-l @@g" | sed "s@'.*@@g")
~/modules/MITObim/misc_scripts/circules.py -f raw_mitobim/${RUN}.uncirc.fasta -k 11 -l ${length} -p circular/${RUN}

#Clean up
sed "1 s/^.*$/>${RUN}/" circular/*circular*.fasta > $circdir/${RUN}.fasta
```

Create tree. I will subset, and create a tree with all egg samples (n=90, including outgroups) and with samples with mtDNA/nuclear coverage > 100x (n=39). 

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=2-00:00:00

mkdir -p trees
# mamba activate phylos
mafft --thread 5 --auto Denovo.fa > Denovo_align.fa
sed '/^>/!s/N/-/g' Denovo_align.fa > Denovo_align.tmp.fa
trimal -in Denovo_align.tmp.fa -out Denovo_align.GT95-s08.fa -gt 0.95 -seqoverlap 80 -resoverlap 0.9
sed -i 's/-/N/g' Denovo_align.GT95-s08.fa
iqtree --keep-ident -s Denovo_align.GT95-s08.fa -pre trees/Denovo_align.GT95-s08 -alrt 1000 -B 1000
```

Tree:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=2-00:00:00

THREADS=5
FASTA=$1
SEQ=$(basename "$FASTA" .fa)

iqtree --keep-ident -s ${FASTA} -pre trees/${SEQ} -alrt 1000 -B 1000
```

Plot:

```R
#### Plot de novo alignment
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/phylogenetics/202507_denovo/scan_n4')
.libPaths('~/r_libs/')
library(ggtree)
library(treeio)
library(ape)
library(tidyverse)
library(phytools)
library(Biostrings)
library(ggpubr)

md = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt')
full = read.iqtree('../all_100x/trees/Denovo_align.GT95-s08.contree')
fr = root(as.phylo(full),'387_CP_MBW_RUS_F')
fulltree <- ggtree(fr, layout = "rectangular") %<+% md +
  geom_nodepoint(mapping=aes(subset=(as.numeric(label) >= 95)),col='black',fill='grey90',pch=23,size=1.5,show.legend=F)+
  geom_tippoint(aes(fill = Haplogroup, shape = SpeciesShort), size = 1.5) +
  scale_shape_manual(values = md$Shape, breaks = md$SpeciesShort) +
  scale_fill_manual(values = md$HaplogroupColor, breaks = md$Haplogroup) +
  guides(fill = guide_legend(nrow = 7, override.aes = list(shape = 21))) +
  ggtitle('n = 127, All Samples. Full mtDNA')
fulltree

### ND2: including Payne & Fossoy
gtree = read.tree('../all_100x/trees/Denovo_With_Genbank_align.GT95.contree')
gr = root(as.phylo(gtree),'387_CP_MBW_RUS_F')
gentree <- ggtree(gr, layout = "rectangular") %<+% md 
gentree$data <- gentree$data %>% mutate(SpeciesShort = ifelse(grepl('CM_',label),'CM',
                                                              ifelse(grepl('CP_',label),'CP',
                                                                     ifelse(grepl('_CO_',label),'CO','CC'))),
                                        Haplogroup = ifelse(grepl('Fos',label),'Fossoy2016',
                                                            ifelse(grepl('Payne',label),'Sorenson2005',Haplogroup)),
                                        showlab = ifelse(grepl('Fos|Payne',label),label,NA))
gt <- gentree + 
  geom_nodepoint(mapping=aes(subset=(as.numeric(label) >= 95)),col='black',fill='grey90',pch=23,size=1.5,show.legend=F)+
  geom_tippoint(aes(fill = Haplogroup, shape = SpeciesShort), size = .7) +
  geom_tiplab(aes(label=showlab),size=2)+
  scale_shape_manual(values = md$Shape, breaks = md$SpeciesShort) +
  scale_fill_manual(values = md$HaplogroupColor, breaks = md$Haplogroup) +
  guides(fill = guide_legend(nrow = 7, override.aes = list(shape = 21))) +
  ggtitle('n = 136, ND2. All Samples, Including Genbank (n = 9)')
gt
bp <- ggarrange(fulltree,gt,common.legend = TRUE)
bp
ggsave('~/symlinks/host/figures/20250805_denovo_tree_withGen.pdf',bp,height=8,width=9)


```



