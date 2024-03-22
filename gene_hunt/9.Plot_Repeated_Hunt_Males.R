#### Plot FST: including males 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/manyhost_hunt/')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(viridis)
library(ggpubr)
library(meRo)
library(RColorBrewer)
library(zoo)
library(ggpubr)
library(scales)

##### Load and prep FST BP data  ######
tidy_bp = read_tsv('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/manyhost_hunt/relatives_hunt/FST1_MirrorRelatives_2024MAR21.txt',col_names=FALSE)
names(tidy_bp) = c('chr','start','p1','p2','FST','genestart','geneend','strand','gene')

groups = tidy_bp %>% select(p1,p2,Group) %>% unique
groups %>% filter(grepl('E1',Group)) %>% data.frame
groups %>% separate(p1,into=c('E1','H1'),remove=F) %>% separate(p2,into=c('E2','H2')) %>% filter(E1 == E2 & H1 != H2)

#BUT, fix misassembly on the Z, assign phylogenetic contrasts 
tidy_bp = tidy_bp  %>% mutate(
  chr = gsub('chr_','',chr),
  chr = ifelse(chr == 'Z' & start > 77e6,'W',chr),
  start = ifelse(chr == 'W' & start > 77e6,start-77e6+22312972,start),
  AvZ = case_when(
    chr == 'W' ~ 'W',
    chr == 'Z' ~ 'Z',
    chr == 'MT' ~ 'MT',
    TRUE ~ 'Autosome'),
  #add floor and ceiling to FST between 0 - 1
  FST = pmax(0, pmin(1, FST)),
  Group = paste0(p1,'__',p2),
  #for each comparison, indicate whether it is ancient comparison (e.g. blue vs noneblue) or intra-clade
  Depth = ifelse(grepl('W1|W2|W3',p1) & !grepl('W1|W2|W3',p2),'Ancient',
                 ifelse(grepl('W1|W2|W3',p2) & !grepl('W1|W2|W3',p1),'Ancient','Contemporary')),
  Comparison = ifelse(grepl('W1|W2',Group) & Depth == 'Ancient','Blue',
                      ifelse(Depth == 'Ancient','Reversion',
                             ifelse(Group == 'E1_W2__E1_W1' | Group == 'E6_W7__E6_W5','Control',
                                    ifelse((grepl('W1|W2|W3',p1) & grepl('W1|W2|W3',p2)) & Depth == 'Contemporary','Reversion','Diversification')))),
  Phylo = paste0(Depth,'_',Comparison),
  site = paste0(chr,'_',start)
) 

tidy_bp %>% count(p1,p2,Phylo) %>% data.frame

#these comparisons do not differ by egg, so remove those sites
tidy_control = tidy_bp %>% filter(Phylo == 'Contemporary_Control')

#remove those
tidy_fixed_targ = tidy_bp %>%
  filter(Group != 'E1_W2__E1_W1' & Group != 'E6_W7__E6_W5') %>% #13 comparisons left
  filter(!site %in% tidy_control$site)

#count 
tidy_fixed_targ %>% ungroup %>% count(p1,p2,Phylo) %>% count(Phylo)
lminput = tidy_fixed_targ %>% select(chr,start,AvZ,gene,Group,p1,p2,Phylo) %>% filter(gene != '.') %>%  separate(p1,into=c('E1','H1'),remove=F) %>% separate(p2,into=c('E2','H2'),remove=F) 

targets = unique(c(lminput$p1,lminput$p2))
target = 'E1_W1'
eggdat = list()

for (targ in targets) {
cat('Finding Overrepresented genes: ',targ,'\n')
t = lminput %>% filter(p1 == targ | p2 == targ) %>% filter(E1 != E2)
phylos = t %>% ungroup %>% select(Group,Phylo) %>% unique %>% count(Phylo) %>% dplyr::rename(num_phylos = n)
tp = t %>% group_by(AvZ,gene,Group,Phylo) %>% 
  summarize(snps = n()) %>% ungroup %>%  #this will reduce so that within each comparison, this gene is now a candidate 
  left_join(.,phylos) %>% 
  group_by(AvZ,Phylo,gene,num_phylos) %>%
  summarize(num_genes = n(),
            proportion_represented = num_genes / num_phylos*100,
            Target = targ)
eggdat[[targ]] = tp
}

library(data.table)
eggs = rbindlist(eggdat) %>% as_tibble %>% unique

#THESE are the genes which in 100% of the phylogenetic contrasts contain an FST=1.0 SNP
eggs100 = eggs %>% filter(proportion_represented == 100)

#so now, let's only examine genes which occur in at least more than 1 contrast across all the different egg:haps
double_obs = eggs100 %>% count(gene) %>% filter(n >= 2)

eggp = eggs100 %>% 
  filter(gene %in% double_obs$gene) %>% ungroup %>% 
  select(AvZ,Phylo,gene,Target) %>% unique %>% 
  group_by(AvZ,Target,gene) %>% 
  mutate(total = n(),gene=toupper(gene)) %>% ungroup %>% 
  mutate(proportion = 1/total) %>% ungroup %>% 
  arrange(desc(total)) %>% 
  mutate(gene = factor(gene, levels = unique(gene))) %>%  
  ggplot(aes(y=gene,fill=Phylo,x=proportion))+
  scale_x_continuous(n.breaks = 2)+
  xlab('Proportion of Phylogenetic Contrasts 100% Observed')+ylab('')+
  facet_grid(AvZ~Target,scales='free',space='free')+
  scale_fill_manual(values=c('cornflowerblue','plum','indianred1','orange'))+
  geom_bar(stat='identity',position='stack')+
  theme_bw(base_size=6)

pdf('../figures/Fixed_SNPs_By_Phylo_RELATIVES-MALES_2024MAR22.pdf',height=4,width=5)
eggp
dev.off()