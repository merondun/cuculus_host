#### Calculate proportion of Z/W/mt/A fixed SNPs across randomized FST pops 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/manyhost_hunt/males/randomizing_fst')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(viridis)
library(ggpubr)
library(meRo)
library(RColorBrewer)
library(zoo)
library(ggpubr)
library(scales)
library(data.table)

##### Load and prep FST BP data  ######
tidy_bp = read_tsv('~/merondun/cuculus_host/gene_hunt/randomized_fst/Fixed_SNPs_RandomizedEggs.txt',col_names = F)
names(tidy_bp) = c('chr','start','iteration','FST','intgene','effect','gene')

#Assign AvZ
tidy_bp = tidy_bp %>% 
  mutate(
    chr = gsub('chr_','',chr),
    chr = gsub('scaffold.*','scaffolds',chr),
    AvZ = case_when(
      chr == 'W' ~ 'W',
      chr == 'Z' ~ 'Z',
      chr == 'MT' ~ 'MT',
      TRUE ~ 'Autosome'),
    #add floor and ceiling to FST between 0 - 1
    FST = pmax(0, pmin(1, FST)),
    site = paste0(chr,'_',start)
  ) 

##### Proportions of FST=1 SNPs #####
tidy_counts <- tidy_bp %>% ungroup %>% 
  # First count fixed SNPs by chromosome 
  group_by(iteration,AvZ) %>%
  summarize(n = n()) %>% ungroup %>%
  # Ensure all combinations of iteration and AvZ exist
  complete(iteration, AvZ, fill = list(n = 0))

# Import the SNP counts per chromosome
counts <- read.table('~/merondun/cuculus_host/gene_hunt/randomized_fst/Total_SNP_Counts.txt',header=F,sep=' ')
names(counts) <- c('SNPs','chr')
counts <- counts %>% 
  mutate(
    chr = gsub('chr_','',chr),
    chr = gsub('scaffold.*','scaffolds',chr),
    AvZ = case_when(
      chr == 'W' ~ 'W',
      chr == 'Z' ~ 'Z',
      chr == 'MT' ~ 'MT',
      TRUE ~ 'Autosome'))
counts_avz <- counts %>% group_by(AvZ) %>% summarize(SNPs = sum(SNPs))

# Merge counts with totals
tidy_prop <- left_join(tidy_counts,counts_avz) %>% 
  mutate(Proportion = n / SNPs,
         Replicate = ifelse(iteration == 'rf0','True','Randomized'))
tidy_prop$AvZ = factor(tidy_prop$AvZ,levels=c('Autosome','Z','W','MT'))

#plot, divide into top and bottom plots for visualization 
props = tidy_prop %>% ggplot(aes(y=iteration,x=Proportion,fill=Replicate,label=n))+
  geom_bar(stat='identity',position=position_dodge(width=0.9))+
  geom_text(position=position_dodge(width=0.9),hjust=-.25,size=1.5)+
  scale_fill_manual(values=rev(brewer.pal(3,'Set2')))+
  #coord_cartesian(xlim=c(0,1.0))+
  scale_x_continuous(n.breaks = 4)+
  facet_wrap(AvZ~.,scales='free')+
  theme_bw(base_size=6)+xlab('Proportion of FST = 1.0 SNPs Present in At Least 1 Comparison')+ylab('')+
  theme(legend.position = "top", # Moves legend to bottom
        legend.text = element_text(size = 6), # Adjusts text size
        legend.title = element_text(size = 6)) + # Adjusts title size
  guides(fill = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 1)))
props

pdf('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/figures/20240907_FST-BP_Females_Proportions-Randomized.pdf',height=4,width=6)
props
dev.off()
