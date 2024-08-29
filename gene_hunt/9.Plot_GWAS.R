#### Identify candidates from GWAS results
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

pops = read_tsv('AllSamples.pop',col_names=F)
names(pops) = c('ID','EggType')
left_join(pops,md) %>% count(Egg,Hap)

##### Load and prep GWAS SNP data   ######
tidy_bp = fread('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/manyhost_hunt/males/GWAS_Results_2024APR2-MultiGene.txt.gz')
names(tidy_bp) = c('chr','start','end','p','Egg','effect','gene','intersect_gene')
setDT(tidy_bp)
genes = unique(tidy_bp, by = c("chr", "start", "p", "Egg", "intersect_gene")) # For sup fig all v one GWAS 
#genes = unique(tidy_bp, by = c("chr", "start", "p", "Egg", "effect", "gene"))

#strip chr_, and assign AvZ
genes[, chr := gsub('chr_', '', chr)][
  , AvZ := fifelse(chr == 'W', 'W',
                   fifelse(chr == 'Z', 'Z',
                           fifelse(chr == 'MT', 'MT', 'Autosome')))]

##### Proportions of Significant SNPs #####
td = genes %>% group_by(Egg) %>% mutate(padj = p.adjust(p,method='bonferroni')) %>% ungroup

candidates = td %>% filter(padj < 0.05)

#total windows, count the proportion of sites which are FST = 1.0 in ALL contrasts for a certain mutation level
chr_windows = td %>% group_by(Egg,AvZ) %>% count() %>% ungroup %>% group_by(Egg,AvZ) %>% slice_max(n) %>% select(Egg,AvZ,NumberWindows = n) %>% distinct
prop = candidates  %>%
  count(AvZ,Egg) %>% unique %>%  #count number of SNPs by phylo 
  left_join(chr_windows,.) %>% ungroup %>% 
  replace_na(list(n=0)) %>% 
  ungroup %>%
  mutate(Proportion = n/NumberWindows)

#plot, divide into top and bottom plots for visualization 
prop$Egg = factor(prop$Egg,levels=eggcols$Egg)
prop$AvZ = factor(prop$AvZ,levels=c('Autosome','Z','W','MT'))
props = prop %>% ggplot(aes(y=Egg,x=Proportion,fill=Egg,label=n))+
  geom_bar(stat='identity',position=position_dodge(width=0.9),col='black')+
  #geom_text(position=position_dodge(width=0.9),hjust=-.25,size=1)+
  scale_fill_manual(values=eggcols$col,breaks=eggcols$Egg)+
  scale_x_continuous(n.breaks = 4)+
  facet_grid(.~AvZ,scales='free')+
  theme_bw(base_size=6)+xlab('Proportion of FST = 1.0 SNPs Present in At Least 1 Comparison')+ylab('')+
  theme(legend.position = "top", # Moves legend to bottom
        legend.text = element_text(size = 6), # Adjusts text size
        legend.title = element_text(size = 6)) + # Adjusts title size
  guides(fill = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 1)))
props

pdf('../../figures/Proportions_GWAS_2024APR3.pdf',height=3,width=6)
props
dev.off()

##### Candidate Genes, here use the nonsynonymnous SNPs after bonferonni correction #####
#filter for bonferonni corrected significant sites 
cands = nonsyn[p <= 0.05/21210834]

#replace LOC IDs with known homologs
locs = read_tsv('../LOC_Lookup.txt')
cands = left_join(cands,locs %>% dplyr::rename(gene=Label)) %>% 
  mutate(gene = ifelse(is.na(ID),gene,paste0(ID,'-L'))) %>% select(-ID) %>% unique

#strip chr_, and assign AvZ
cands[, chr := gsub('chr_', '', chr)][
  , AvZ := fifelse(chr == 'W', 'W',
                   fifelse(chr == 'Z', 'Z',
                           fifelse(chr == 'MT', 'MT', 'Autosome')))]

#out of those, which are overrepresented genes 
length(unique(cands$Egg))
nonsyn_cands = cands %>% filter(gene != '.' & effect != '.' & effect != 'synonymous') 
targets = nonsyn_cands %>% 
  filter(gene != '.' & effect != '.' & effect != 'synonymous') %>% 
  #first count the number of mutations for each gene, for each egg
  group_by(AvZ,Egg,gene) %>% 
  summarize(snps = n()) %>% ungroup %>% 
  #now, find the genes which have multiple eggs with dN 
  group_by(AvZ,gene) %>% 
  summarize(Group = paste0(Egg,collapse='__'), total_snps = sum(snps),Groups = n()) %>% ungroup %>% 
  filter(Groups >= 4) 
targets

#save those genes 
write.table(targets$gene,file='../../figures/Candidate_Genes_GWAS-2Contrast_2024APR2.txt',quote=F,sep='\t',row.names=F,col.names=F)

#plot those 
gene_candidates = nonsyn_cands %>% filter(gene %in% targets$gene) %>% 
  #make all genes uppercase, and replace some of the mtDNA genes 
  mutate(gene = toupper(gsub('ID=','',gene)),
         gene = gsub('COX1','CO1',gene),
         gene = gsub('NAD1','ND1',gene),
         gene = gsub('NAD2','ND2',gene),
         gene = gsub('NAD3','ND3',gene),
         gene = gsub('NAD4','ND4',gene),
         gene = gsub('NAD5','ND5',gene),
         gene = gsub('NAD6','ND6',gene)) %>% 
  #arrange genes according to the highest number of eggs present 
  select(chr,Egg,gene,AvZ) %>% distinct %>% 
  group_by(gene) %>% 
  mutate(total = n()) %>% ungroup %>% 
  arrange(desc(total)) %>% 
  mutate(gene = factor(gene, levels = unique(gene))) %>%  
  ggplot(aes(y=gene,fill=Egg))+
  geom_bar(position='stack',col='black')+
  facet_grid(AvZ~.,scales='free',space='free')+
  scale_x_continuous(breaks= pretty_breaks(n=3))+
  scale_fill_manual(values=eggcols$col,breaks=eggcols$Egg)+
  theme_bw(base_size=6)+
  theme(legend.position='top')
gene_candidates

pdf('../../figures/Fixed_dN-SNPS_GWAS_2024APR2.pdf',height=3.5,width=2)
gene_candidates
dev.off()

#what are those SNPs?
snp_sites = nonsyn_cands %>% filter(gene %in% targets$gene) %>% arrange(chr,start,gene)
write.table(snp_sites,file='../../figures/GWAS_SNP-Sites_2024APR4.txt',quote=F,sep='\t',row.names=F)

##### Plot GWAS #####
library(karyoploteR)
#prep karyoplot
genome = read.table('/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa.bed',header=FALSE) %>% filter(str_detect(V1,'scaff',negate=T)) %>% arrange(desc(V2))
names(genome) = c('chr','start','end'); genome$chr = gsub('chr_','',genome$chr); 
genome = genome %>% mutate(chr = factor(chr, levels = c(as.character(1:39), "Z", "W", "MT"))) %>%  arrange(chr)
G = makeGRangesFromDataFrame(genome,seqnames.field = 'chr',start.field = 'start',end.field = 'end')

#add colors 
chrs = genome %>% select(chr) %>% mutate(Color = if_else(row_number() %% 2 == 0, 'grey60', 'black'))
candsgr = makeGRangesFromDataFrame(cands,keep.extra.columns = TRUE,seqnames.field = 'chr',start.field = 'start',end.field = 'end')

# General karyoplot parameters
pdf('../../figures/KARYOPLOT_GWAS_2024APR3.pdf',height=5,width=7.5)
pp = getDefaultPlotParams(plot.type=4)
pp$leftmargin = 0.1
kp = plotKaryotype(plot.type=4, genome = G,
                   labels.plotter = NULL,
                   plot.params = pp)
kpAddChromosomeNames(kp, yoffset = -5,cex=0.7)
kpAddBaseNumbers(kp,tick.dist = 50000000,minor.ticks = FALSE,cex=0.4)

#Loop through each track (all species together, CC, and CO), and plot them on their own layer 
counter = 0; tracks = 12
for (egg in eggcols$Egg) { 
  col = eggcols %>% filter(Egg == egg) %>% pull(col)
  counter = counter + 1; at = autotrack(current.track = counter, total.tracks = tracks);at$r1 <- at$r1-0.02
  subgr = candsgr[mcols(candsgr)$Egg == egg]
  kpPlotDensity(kp,data=subgr,r0=at$r0,r1=at$r1,col=col)
  kpAxis(kp,r0=at$r0,r1=at$r1,cex=0.3,numticks = 2)
  kpAddLabels(kp,cex=0.25,labels = egg,r0=at$r0+.01, r1=at$r1,col="black",srt=0,label.margin = 0.02)
  
}

dev.off()
