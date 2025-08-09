## Twisst

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

Run twisst:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --mem=30000mb
#SBATCH --cpus-per-task=10
#SBATCH --time=24:00:00

# mamba activate twisst

# First merge samples
for CHR in chr_W chr_MT; do

        bcftools merge --file-list All_${CHR}.vcfs -Oz -o twisst/All_${CHR}.vcf.gz
        bcftools index twisst/All_${CHR}.vcf.gz

done
bcftools concat -Oz -o twisst/All.vcf.gz twisst/All_chr_MT.vcf.gz twisst/All_chr_W.vcf.gz
bcftools index twisst/All.vcf.gz

WD=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/phylogenetics/202507_rooting_tests/twisst
cd ${WD}

mkdir -p input output

# Force haploid
bcftools view All.vcf.gz | \
  bcftools +fixploidy -Ou - -- -f 1 | \
  bcftools view -e 'F_MISSING > 0.1' --min-ac 1 --min-alleles 2 --max-alleles 2 -Oz -o input/All-Haploid.vcf.gz

# Create geno file
python ~/modules/genomics_general/VCF_processing/parseVCF.py --ploidyMismatchToMissing --ploidy 1 --skipIndels -i input/All-Haploid.vcf.gz | \
        bgzip > input/All.geno.gz

for win in 10 25 50 100 250; do
  echo "Working on rep ${REP} for win ${win}"
  miss=$(echo ${win} | awk -v w=${win} '{m = int((w/5)+0.5); print (m < 1 ? 1 : m)}')
  missind=$(echo ${win} | awk -v w=${win} '{m = int((w/10)+0.5); print (m < 1 ? 1 : m)}')

  #calculate trees in windows
  python ~/modules/genomics_general/phylo/phyml_sliding_windows.py -T 10 \
          -g input/All.geno.gz --prefix output/All-T1.Sw${win} \
         --windSize ${win} --minSites ${miss} --minPerInd ${missind} --windType sites --model GTR --optimise n

  #run twisst
  ~/modules/twisst/twisst.py -t output/All-T1.Sw${win}.trees.gz \
    -w output/All-T1.Sw${win}.csv.gz --outgroup CG \
    --outputTopos output/All-T1.topologies.Sw${win}.trees -g CG -g CP -g CM -g BLUE -g GRY -g RUF --method complete --groupsFile all.pop

done

```

Plot twisst:

```bash
#### Plot twisst
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/phylogenetics/202507_rooting_tests/twisst/output')
library(tidyverse)
library(viridis)
library(RColorBrewer)
library(ggpubr)
library(meRo) #devtools::install_github('merondun/meRo')
source('~/modules/twisst/plot_twisst.R')

# Import data 
weights_files = list.files('.',pattern='All-T1.*.csv.gz')
window_data_files = list.files('.',pattern='All-T1.*data.tsv')

coordat <- list()
plots <- list()
for (size in seq(1,4,1)) {
  # Read in the data
  win <- gsub('.csv.gz','',gsub('.*Sw','',weights_files[size]))
  twisst_data <- import.twisst(weights_files=weights_files[size],window_data_files=window_data_files[size])
  
  # Save plot 
  plot.twisst.summary.boxplot(twisst_data, lwd = 3, cex = 0.7, only_best = 5,)
  plots[[paste0("p", size)]] <- recordPlot()
  
  # Also extract coordinate data 
  dat <- as_tibble(rbindlist(twisst_data$window_data))
  names(dat) <- c('chr','start','end','d','sites','lnl')
  dat <- dat %>% select(-d)
  wdat <- cbind(dat,rbindlist(twisst_data$weights))
  wfdat <- wdat %>% 
   pivot_longer(!c(chr,start,end,sites,lnl),names_to = 'topology',values_to='weight') %>% 
   mutate(topology = gsub('chr_W.','',topology),
          win = win)
  
  # Weight 
  coordat[[size]] <- wfdat
}
plots[[1]]
plots[[2]]
plots[[3]]
plots[[4]]

cdat <- rbindlist(coordat) %>% as_tibble
top_topos <- cdat %>% 
  group_by(chr,topology) %>% 
  sum_stats(weight) %>% 
  ungroup %>% 
  group_by(chr) %>% 
  slice_max(mean,n=10,with_ties = FALSE) %>% 
  ungroup %>% 
  group_by(topology) %>% 
  summarize(avg = mean(mean)) %>% 
  arrange(desc(avg)) %>% head(n=10) %>% 
  mutate(col = sample(viridis(10,option='turbo')))
top_topos

# Summarize topos
topo_sums <- cdat %>% 
  group_by(chr,topology,win) %>% 
  sum_stats(weight)
topo_sums$topology <- factor(topo_sums$topology,levels=top_topos$topology)
topo_sums$win <- factor(topo_sums$win,levels=c(25,50,100,250))
plot <- topo_sums %>% 
  filter(topology %in% top_topos$topology) %>% 
  ggplot(aes(x=topology,y=mean,ymin=pmax(conf_low,0),ymax=conf_high,col=topology))+
  geom_point()+
  facet_grid(win~chr,scales='free')+
  scale_color_manual(values=top_topos$col,breaks=top_topos$topology)+
  geom_errorbar()+
  theme_bw(base_size=12)+
  coord_cartesian(ylim=c(0,1))
plot

ggsave('~/symlinks/host/figures/20250731_TWISST_chrMT-W_Sensitivity.pdf',
       plot, height = 5,width=9,dpi=300)

##### Plot topologies using Ape #####
topo_ids <- as.integer(gsub('topo','',top_topos$topology))
pdf('~/symlinks/host/figures/20250731_Twisst_Top_Topologies.pdf',height=1.25,width=15)

par(mfrow = c(1,length(topo_ids)), mar = c(1,1,2,1), xpd=NA)
for (i in seq_along(topo_ids)) {
  t <- topo_ids[i]
  col <- top_topos %>% filter(topology == paste0('topo',t)) %>% pull(col)
  plot.phylo(twisst_data$topos[[t]], type = "clad", edge.color=col, edge.width=3, label.offset=.2, cex = 1)
  mtext(side=3, text=paste0("topo", t))
}

dev.off()

# Find out where the top 4 topos lie: 
cdatc <- left_join(cdat,top_topos %>% select(topology,col))
keep <- c('topo1','topo4','topo13','topo3')
cdatc$win <- factor(cdatc$win,levels=c(25,50,100,250))
leg <- cdatc %>% select(topology,cols) %>% distinct
length_plot <- cdatc %>%
  filter(topology %in% keep) %>%
  ggplot(aes(xmin = start, xmax = end, ymin = 0, ymax = weight, fill = topology)) +
  geom_rect(color = NA) +
  theme_bw() +
  scale_fill_manual(values=cdatc$col,breaks=cdatc$topology)+
  facet_grid(win ~ chr, scales = 'free')
ggsave('~/symlinks/host/figures/20250731_Twisst_TopologyAlongChrs.pdf',height=5,width=10)

```


