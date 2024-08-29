#### Identify candidates from female-only FST scan 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/manyhost_hunt/males')
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

md = read_tsv('~/merondun/cuculus_host/Metadata_Host.txt')
eggcols = md %>% ungroup %>% select(Egg) %>% unique %>% mutate(ord = as.numeric(gsub('E','',Egg))) %>% arrange(ord) %>% 
  drop_na(Egg) %>% mutate(col = viridis(12,option='turbo')[1:11])
eggcols = rbind(eggcols,data.frame(Egg='E6W3',ord=5,col='white')) %>% arrange(ord)

##### Load and prep FST BP data  ######
tidy_bp = fread('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/manyhost_hunt/males/FST_Contrasts_2024APR4.txt.gz')
names(tidy_bp) = c('chr','start','p1','p2','FST','intgene','effect','gene')

#Assign AvZ
tidy_bp = tidy_bp %>% 
  mutate(
    chr = gsub('chr_','',chr),
    AvZ = case_when(
      chr == 'W' ~ 'W',
      chr == 'Z' ~ 'Z',
      chr == 'MT' ~ 'MT',
      TRUE ~ 'Autosome'),
    #add floor and ceiling to FST between 0 - 1
    FST = pmax(0, pmin(1, FST)),
    site = paste0(chr,'_',start)
  ) 

phylo = read_tsv('Contrast_Metadata.txt') %>% select(p1,p2,Group,Phylo)
tidy_dat = left_join(tidy_bp,phylo) %>% mutate(gene = toupper(gene),
                                               intgene = toupper(intgene)) 

#grab only the fixed sites!
tidy_fixed = tidy_dat %>% filter(FST == 1) 
write.table(tidy_fixed,file='FST1_TIDY_2024APR4.txt',quote=F,sep='\t',row.names=F)
tidy_fixed = read_tsv('FST1_TIDY_2024APR4.txt')

fwrite(tidy_dat,'FST_TIDY_2024APR4.txt.gz',scipen = getOption('scipen', 999L),sep='\t')
tidy_dat = fread('FST_TIDY_2024APR4.txt.gz')

###### Identify Candidate SNPS ###### 
#these comparisons do not differ by egg, so remove those sites
tidy_control = tidy_fixed %>% filter(Phylo == 'Contemporary_Control')

#remove those
tidy_fixed_targ = tidy_fixed %>%
  filter(Phylo != 'Contemporary_Control') %>% #13 comparisons left
  filter(!site %in% tidy_control$site)

#count the number of contrasts, for now, since we may have duplicate sites (overlapping CDS) assigned, ensure that we only have single
tidy_prop_input = tidy_fixed_targ %>% select(chr,start,AvZ,p1,p2,Group,Phylo,FST) %>% distinct
tidy_prop_input %>% select(Group,Phylo) %>% distinct %>% count(Phylo)
tidy_prop_input %>% count(Phylo,Group,chr,start) %>% filter(n > 1) #this should be zero

#requiring at least one contrast...
tidy_prop_input_filtered = tidy_prop_input %>% group_by(Phylo,AvZ,chr,start) %>% count %>%
  filter(
    ifelse(Phylo == 'Ancient_Blue', n >= 3, 
           ifelse(Phylo == 'Contemporary_Diversification', n >= 1,
                  ifelse(Phylo == 'Contemporary_Reversion', n >= 1, n >= 999))))

#requiring all contrasts for the same SNP (used in primary figure!) 
tidy_prop_input_filtered = tidy_prop_input %>% group_by(Phylo,AvZ,chr,start) %>% count %>%
  filter(
    ifelse(Phylo == 'Ancient_Blue', n >= 6, 
           ifelse(Phylo == 'Contemporary_Diversification', n >= 2,
                  ifelse(Phylo == 'Contemporary_Reversion', n >= 2, n >= 999))))

##### Proportions of FST=1 SNPs #####
#total windows, count the proportion of sites which are FST = 1.0 in ALL contrasts for a certain mutation level
chr_windows = tidy_dat %>% select(chr,start,AvZ,Group,Phylo) %>% distinct %>% count(AvZ,Group,Phylo) %>% ungroup %>% group_by(Phylo,AvZ) %>% slice_max(n) %>% select(Phylo,AvZ,NumberWindows = n) %>% unique %>% filter(!grepl('Control',Phylo))
tidy_fixed_prop = tidy_prop_input_filtered  %>% ungroup %>% 
  count(AvZ,Phylo) %>% distinct %>%  #count number of SNPs by phylo 
  left_join(chr_windows,.) %>% ungroup %>% 
  replace_na(list(n=0)) %>% 
  ungroup %>%
  mutate(Proportion = n/NumberWindows)
tidy_fixed_prop$AvZ = factor(tidy_fixed_prop$AvZ,levels=c('Autosome','Z','W','MT'))
tidy_fixed_prop$Phylo = factor(tidy_fixed_prop$Phylo,levels=c('Ancient_Blue','Contemporary_Reversion','Contemporary_Diversification'))

#check colorblind
# library(colorblindcheck)
# pal <- c('cadetblue3', 'plum', 'indianred1')
# pal_hex <- col2rgb(pal) / 255
# pal_hex <- apply(pal_hex, 2, function(x) rgb(x[1], x[2], x[3]))
# show_col(pal_hex)
# palette_check(pal_hex,plot=TRUE)

#plot, divide into top and bottom plots for visualization 
props = tidy_fixed_prop %>% ggplot(aes(y=AvZ,x=Proportion,fill=Phylo,label=n))+
  geom_bar(stat='identity',position=position_dodge(width=0.9))+
  geom_text(position=position_dodge(width=0.9),hjust=-.25,size=3)+
  scale_fill_manual(values=rev(brewer.pal(3,'Set2')))+
  #coord_cartesian(xlim=c(0,1.0))+
  scale_x_continuous(n.breaks = 4)+
  #facet_grid(AvZ~.,scales='free',nrow=4,ncol=1)+
  theme_bw(base_size=6)+xlab('Proportion of FST = 1.0 SNPs Present in At Least 1 Comparison')+ylab('')+
  theme(legend.position = "top", # Moves legend to bottom
        legend.text = element_text(size = 6), # Adjusts text size
        legend.title = element_text(size = 6)) + # Adjusts title size
  guides(fill = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 1)))
props

pdf('../../figures/20240718_FST-BP_Females_Proportions-Half.pdf',height=2,width=4)
props
dev.off()

#are there regions of MT/W which contain more differentiated windows?
dens = tidy_prop_input_filtered %>% filter((chr == 'W' | chr == 'MT')) %>%
  ggplot(aes(x=start,fill=Phylo))+
  scale_fill_manual(values=c('cadetblue3','goldenrod1','indianred1'))+
  geom_density(alpha=0.8)+facet_wrap(.~chr,scales='free',nrow=2,ncol=1)+theme_bw(base_size=6) +theme(legend.position='top')
dens

pdf('../../figures/FST-BP_Females_Density_2024APR4.pdf',height=4,width=2.5)
dens
dev.off()

##### dN and dS SNPs #####

#these are our strongest candidate SNPs  
candidates = left_join(tidy_prop_input_filtered %>% select(-n),tidy_fixed_targ)

#replace LOC IDs with known homologs
locs = read_tsv('../LOC_Lookup.txt')

#based only on overlapping genes 
cands = left_join(candidates,locs %>% dplyr::rename(intgene=Label)) %>% 
  mutate(gene = ifelse(is.na(ID),intgene,paste0(ID,'-L'))) %>% select(-ID) %>% unique %>% 
  mutate(gene = gsub('NAD','ND',gene), gene = gsub('_A','',gene)) 
cand_genes = cands %>% filter(gene != '.') %>% ungroup %>% count(gene,Phylo,AvZ,chr) 
all_targets = cand_genes %>% ungroup %>% group_by(AvZ,gene) %>% 
  summarize(Group = paste0(Phylo,collapse='__'), total_snps = sum(n),Groups = n()) %>% ungroup 
targets = all_targets %>% filter(Groups >= 1)

#based on nonsynonymous/nonsense (shown in primary figure!) 
cands = left_join(candidates,locs %>% dplyr::rename(gene=Label)) %>% 
  mutate(gene = ifelse(is.na(ID),gene,paste0(ID,'-L'))) %>% distinct %>% 
  mutate(gene = gsub('NAD','ND',gene), gene = gsub('_A','',gene)) 
cand_genes = cands %>% filter(effect != '.' & effect != 'synonymous') %>% ungroup %>% count(gene,Phylo,AvZ,chr) 
cands %>% filter(effect == 'nonsense') 

#plot how many dN SNPs are in each gene, across each phylo category 
cand_genes$Phylo = factor(cand_genes$Phylo,levels=c('Ancient_Blue', 'Contemporary_Reversion','Contemporary_Diversification'))

#Vertical 
p = cand_genes %>% 
  filter(gene %in% targets$gene) %>%   #hash out for dN 
  select(Phylo,gene,AvZ) %>% distinct %>% 
  group_by(gene) %>% 
  mutate(total = n()) %>% ungroup %>% 
  arrange(desc(total)) %>% 
  mutate(gene = factor(gene, levels = unique(gene))) %>%  
  ggplot(aes(y=gene,fill=Phylo))+
  geom_bar(position='stack')+
  facet_grid(AvZ~.,scales='free',space='free')+
  scale_x_continuous(breaks= pretty_breaks(n=3))+
  #scale_shape_manual(values=rev(c(16,15,18,17)))+
  scale_fill_manual(values=rev(brewer.pal(3,'Set2')))+
  theme_bw(base_size=6)+
  theme(legend.position='top')
p

# Horizontal  
cand_genes$AvZ <- factor(cand_genes$AvZ,levels=c('MT','W','Z','Autosome'))
p = cand_genes %>% 
  filter(gene %in% targets$gene) %>%   #hash out for dN 
  select(Phylo,gene,AvZ) %>% distinct %>% 
  group_by(gene) %>% 
  mutate(total = n()) %>% ungroup %>% 
  arrange(desc(total)) %>% 
  mutate(gene = factor(gene, levels = unique(gene))) %>%  
  ggplot(aes(x=gene,fill=Phylo))+
  geom_bar(position='stack')+
  facet_grid(.~AvZ,scales='free',space='free')+
  scale_y_continuous(breaks= pretty_breaks(n=3))+
  scale_fill_manual(values=rev(brewer.pal(3,'Set2')))+
  theme_bw(base_size=6)+xlab('')+ylab('')+
  theme(legend.position='top',
        axis.text.x=element_text(angle=45,hjust=1,vjust=1))
p

pdf('../../figures/20240722_FST-BP_Females_dN_Mutations_Genes_Horizontal.pdf',height=2.5,width=6)
p
dev.off()

#save those genes
genes = targets %>% mutate(gene = gsub('-L','',gene))
write.table(genes$gene,file='../../figures/FST-BP_Females_CandidatesGenesOverlap_2024APR4.txt',quote=F,sep='\t',row.names=F,col.names=F)

#save nonsynonymous genes
dn_genes = cand_genes %>% mutate(gene = gsub('-L','',gene))
write.table(dn_genes$gene,file='../../figures/FST-BP_Females_CandidatesGenesdN_2024APR4.txt',quote=F,sep='\t',row.names=F,col.names=F)

#grab the actual SNPs from the main positions 
inspection = cands %>% filter(effect == 'nonsynonymous' & gene %in% cand_genes$gene) %>% ungroup %>% select(chr,start,gene,Phylo) %>% unique %>% group_by(chr,gene,Phylo) %>% 
  #slice_max(start) %>% 
  data.frame
inspect_bed = inspection %>% filter(grepl('ND2|ND5|ND4|TLDC2|NDUFAF4',gene)) %>% mutate(end = start,chr=paste0('chr_',chr)) %>% select(chr,start,end,Phylo,gene)
write.table(inspect_bed,file='Fixed_SNP_Candidates_PositionsForGenotypePlot_2024JULY18.txt',quote=F,sep='\t',row.names=F,col.names=F)

#in bash:
for chr in $(awk '{print $1}' Fixed_SNP_Candidates_PositionsForGenotypePlot_2024JULY18.txt | sort | uniq); do
echo "WORKING ON CHR: ${chr}"
bedtools intersect -header -a female_vcfs/${chr}.FEM.DP3.vcf.gz -b Fixed_SNP_Candidates_PositionsForGenotypePlot_2024JULY18.txt | bcftools view -Oz -o plot_genos/${chr}.vcf.gz
done
bcftools concat -Oz -o plot_genos/Targets.vcf.gz plot_genos/chr_*.vcf.gz

##### Plot Genotypes #####
library(vcfR)
library(tidyverse)
library(data.table)
library(RColorBrewer)
#plot genotypes 
vcfr = read.vcfR('plot_genos/Targets.vcf.gz')
gt_data <- vcfR2tidy(vcfr, format_fields = 'GT')
gt = gt_data$gt
#remember vcfs are 1-based
names(gt) = c('Key','site','ID','Genotype','Allele')

#add metadata since I want to group by morph
md = read_tsv('~/merondun/cuculus_host/Metadata_Host.txt')
md = md %>% drop_na(Egg)
gtm = left_join(md %>% dplyr::select(ID,Hap,Egg),gt) #ensure that the 'ID' from the vcf matches the 'ID' from your metadata

#now join with the annotated snps
agt = gtm %>% drop_na(site) %>% replace_na(list(Genotype = 'Missing'))
agt %>% count(Key,site)
aga = left_join(agt,inspection %>% dplyr::rename(site=start)) %>% drop_na(gene) %>% mutate(gene = toupper(gsub('ID=','',gene)))
pg = aga %>% arrange(chr,site)
sites = pg %>% select(chr,site) %>% unique
pg$site = factor(pg$site,levels=sites$site)

samples = read_tsv('Female_Hunt.list',col_names=F)
names(samples) = 'ID'

pg = pg %>% filter(ID %in% samples$ID)
pg$Phylo = factor(pg$Phylo,levels=c('Ancient_Blue','Contemporary_Reversion','Contemporary_Diversification'))
pg$gene = factor(pg$gene,levels=c('NDUFAF4-L','ND2','ND4','ND5','TLDC2'))
pg$Egg = factor(pg$Egg,levels=eggcols$Egg)
gtp = pg %>% 
  group_by(chr,gene,site,Genotype,Allele,ID,Hap,Egg) %>% 
  mutate(Genotype = ifelse(Genotype == '0','0/0',
                           ifelse(Genotype == '1','1/1',Genotype))) %>% 
  #summarize(Group = paste0(Phylo,collapse='__')) %>% 
  ggplot(aes(y = ID, x = as.factor(site), fill = Genotype)) +
  geom_tile() +
  ylab("Individuals") +
  facet_grid(Hap+Egg~gene+chr+Phylo,space='free',scales='free')+
  scale_fill_manual(values=c('grey70','grey40','grey10','white'))+
  #scale_fill_manual(values=c('grey70','grey20','red'))+
  theme_minimal(base_size=6) + ylab('')+xlab('') +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.x = element_text(size=4,angle = 90, vjust = 0.5, hjust = 1),
        strip.text.y.right = element_text(angle = 0),legend.position='top',
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
gtp

pdf('../../figures/Genotypes-Targets_2024APR4.pdf',height=3.5,width=1.75)
gtp
dev.off()
