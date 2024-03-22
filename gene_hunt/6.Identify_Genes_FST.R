#### Identify candidates from BP FST, extract CNV, plot genotypes
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
tidy_bp = read_tsv('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/manyhost_hunt/ManyHost_FST_1BP-VCFTOOLS_2024MAR12.txt',col_names=FALSE)
names(tidy_bp) = c('p1','p2','chr','start','FST')

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
                             ifelse(Group == 'W2__Pphoenicurus__W1__Pphoenicurus' | Group == 'W7__Aarundinaceus__W5__Aarundinaceus','Control',
                                    ifelse((grepl('W1|W2|W3',p1) & grepl('W1|W2|W3',p2)) & Depth == 'Contemporary','Reversion','Diversification')))),
  Phylo = paste0(Depth,'_',Comparison),
  site = paste0(chr,'_',start)
) 

#grab only the fixed sites!
tidy_fixed = tidy_bp %>% filter(FST == 1) 
write.table(tidy_fixed,file='FST-BP_ManyHosts_Fixed-VCFTOOLS_FST1_2024MAR12.txt',quote=F,sep='\t',row.names=F)
tidy_fixed = read_tsv('FST-BP_ManyHosts_Fixed-VCFTOOLS_FST1_2024MAR12.txt')

#these comparisons do not differ by egg, so remove those sites
tidy_control = tidy_fixed %>% filter(Comparison == 'Control')

#remove those
tidy_fixed_targ = tidy_fixed %>%
  filter(Group != 'W2__Pphoenicurus__W1__Pphoenicurus' & Group != 'W7__Aarundinaceus__W5__Aarundinaceus') %>% #13 comparisons left
  filter(!site %in% tidy_control$site)

#count 
tidy_bp %>% ungroup %>% count(p1,p2,Phylo) %>% count(Phylo)
#only use SNPs in which ALL of the phylogenetic contrasts are fulfilled 
tidy_fixed_targ_all = tidy_fixed_targ %>% group_by(Phylo,AvZ,chr,start) %>% count %>%
  filter(
    ifelse(Phylo == 'Ancient_Blue', n == 6, 
           ifelse(Phylo == 'Ancient_Reversion', n == 3,
                  ifelse(Phylo == 'Contemporary_Diversification', n == 2,
                         ifelse(Phylo == 'Contemporary_Reversion', n == 2, n >= 999)))))

##### Proportions of FST=1 SNPs #####
#total windows, count the proportion of sites which are FST = 1.0 in ALL contrasts for a certain mutation level
chr_windows = tidy_bp %>% group_by(Phylo,Group,AvZ) %>% count() %>% ungroup %>% group_by(Phylo,AvZ) %>% slice_max(n) %>% select(Phylo,AvZ,NumberWindows = n) %>% unique %>% filter(!grepl('Control',Phylo))
tidy_fixed_prop = tidy_fixed_targ_all %>%
  group_by(Phylo,AvZ,chr,start) %>% count() %>% select(-n) %>% #first, extract only the bp for each group 
  left_join(.,chr_windows) %>% unique %>%
  ungroup %>%
  group_by(Phylo,AvZ,NumberWindows) %>% count %>% #now, count the number of outlier sites in each comparison across each AvZ
  mutate(Proportion = n/NumberWindows)
tidy_fixed_prop = tidy_fixed_prop %>% ungroup %>% complete(Phylo, AvZ, fill = list(n = 0, Proportion = 0))
tidy_fixed_prop$AvZ = factor(tidy_fixed_prop$AvZ,levels=c('Autosome','Z','W','MT'))
tidy_fixed_prop$Phylo = factor(tidy_fixed_prop$Phylo,levels=c('Ancient_Blue','Ancient_Reversion','Contemporary_Reversion','Contemporary_Diversification'))

#plot, divide into top and bottom plots for visualization 
propstop = tidy_fixed_prop %>% ggplot(aes(y=AvZ,x=Proportion,fill=Phylo,label=n))+
  geom_bar(stat='identity',position=position_dodge(width=0.9))+
  geom_text(position=position_dodge(width=0.9),hjust=-.25,size=1)+
  scale_fill_manual(values=c('cornflowerblue','plum','orange','indianred1'))+
  coord_cartesian(xlim=c(0.1,1.0))+
  theme_bw(base_size=6)+xlab('Proportion of FST = 1.0 SNPs Present in At Least 1 Comparison')+ylab('')+
  theme(legend.position = "top", # Moves legend to bottom
        legend.text = element_text(size = 6), # Adjusts text size
        legend.title = element_text(size = 6)) + # Adjusts title size
  guides(fill = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 1)))
propsbottom = tidy_fixed_prop %>% ggplot(aes(y=AvZ,x=Proportion,fill=Phylo,label=n))+
  geom_bar(stat='identity',position=position_dodge(width=0.9))+
  geom_text(position=position_dodge(width=0.9),hjust=-.25,size=1)+
  scale_fill_manual(values=c('cornflowerblue','plum','orange','indianred1'))+
  coord_cartesian(xlim=c(0,0.1))+
  theme_bw(base_size=6)+xlab('Proportion of FST = 1.0 SNPs Present in At Least 1 Comparison')+ylab('')+
  theme(legend.position = "top", # Moves legend to bottom
        legend.text = element_text(size = 6), # Adjusts text size
        legend.title = element_text(size = 6)) + # Adjusts title size
  guides(fill = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 1)))

pdf('../figures/FST-BP_ManyHosts_Proportions-All_2024MAR06.pdf',height=3,width=2.5)
ggarrange(propsbottom,propstop,widths=c(0.7,0.3))
dev.off()

#these are our strongest candidate SNPs  
tidy_fixed_targ_all %>% ungroup %>% count(Phylo,AvZ)

#are there regions of MT/W which contain more differentiated windows?
dens = final_targets %>% filter((chr == 'W' | chr == 'MT')) %>%
  ggplot(aes(x=start,fill=Phylo))+
  scale_fill_manual(values=c('cornflowerblue','plum','indianred1','orange'))+
  geom_density(alpha=0.8)+facet_wrap(.~chr,scales='free',nrow=2,ncol=1)+theme_bw(base_size=6) +theme(legend.position='top')
#dens
pdf('../figures/FST-BP_ManyHosts_Density-All-W-MT_2024MAR19.pdf',height=4,width=2.5)
dens
dev.off()

##### dN and dS SNPs #####
#merge with syn/nonyn annotations from variantannotation 
dnds = read_tsv('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/manyhost_hunt/dnds/DNDS_Annotated_2024MAR07.txt')
names(dnds) = c('chr','start','REF','ALT','effect','gene')
dnds$start = as.integer(dnds$start)
dnds = dnds %>%
  mutate(chr = gsub('chr_','',chr),
         start = ifelse(chr == 'Z' & start > 77e6,start-77e6+22312972,start),
         chr = ifelse(chr == 'Z' & start > 77e6,'W',chr))
fbp = left_join(tidy_fixed_targ_all,dnds)
fbp %>% drop_na(effect) %>% filter(effect != 'synonymous')
fbp %>% drop_na(effect) %>% filter(effect != 'synonymous') %>% ungroup %>% count(Phylo,AvZ,chr,effect,gene) %>% arrange(desc(n)) %>% data.frame
fbp %>% filter(effect=='nonsynonymous') %>% ungroup %>% select(Phylo,AvZ,gene) %>% unique %>% count(AvZ,gene) %>% slice_max(n,with_ties = TRUE)
plot_fixed = fbp %>% drop_na(effect) %>% filter(effect != 'synonymous') %>% ungroup 
#num genes fixed with dN by chr... MANY! 
plot_fixed %>% select(chr,Phylo,gene) %>% unique %>% count(Phylo,chr)

#plot how many dN SNPs are in each gene, across each phylo category 
plot_fixed$Phylo = factor(plot_fixed$Phylo,levels=c('Contemporary_Diversification','Contemporary_Reversion','Ancient_Reversion','Ancient_Blue'))
p = plot_fixed %>% ungroup %>% filter(effect == 'nonsynonymous') %>% count(Phylo,gene,AvZ,effect)  %>% 
  mutate(gene = toupper(gsub('gene-','',gene))) %>% 
  ggplot(aes(x=gene,y=n,fill=Phylo))+
  geom_bar(stat='identity',position='stack')+
  facet_grid(.~AvZ,scales='free',space='free')+
  scale_y_continuous(breaks= pretty_breaks(n=3))+
  #scale_shape_manual(values=rev(c(16,15,18,17)))+
  scale_fill_manual(values=rev(c('cornflowerblue','plum','orange','indianred1')))+
  theme_bw(base_size=4)+
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1))

p

pdf('../figures/dN_Mutations_2024MAR11.pdf',height=1.5,width=4.5)
p
dev.off()

#save
write.table(fbp,file='FST-BP_dNdS_Overlap_MajorityDiverged-All_2024MAR11.txt',quote=F,sep='\t',row.names=F)

#plot pie charts of dN/dS/dNonsense
#count proportion of each by category and by AvZ
conseq_summary = fbp %>% drop_na(effect) %>%
  group_by(Phylo,AvZ,effect) %>% mutate(CountEffect = n()) %>% ungroup %>% select(Phylo,AvZ,effect,CountEffect) %>%
  group_by(Phylo,AvZ) %>%
  mutate(TotalEffect = n()) %>% ungroup %>%
  mutate(prop = CountEffect / TotalEffect) %>% unique %>%
  arrange(Phylo,AvZ)
conseq_summary

#plot pie charts
conseq_sum = ggplot(conseq_summary, aes(x = "", y = prop, fill = effect)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  facet_grid(AvZ~Phylo) +  # Faceting by group
  scale_fill_viridis(discrete=TRUE,option='viridis')+
  geom_label(aes(label = CountEffect),size=1,position = position_stack(vjust = 0.5),fill='white') +  # Add count labels
  theme_void(base_size=6) +  # Remove axis labels and ticks
  theme(legend.position = "bottom")  # Adjust legend position
pdf('../figures/FST-BP_ManyHosts-All_dNdS_Pies_2024MAR06.pdf',height=3,width=2)
conseq_sum
dev.off()

##### Find Top Candidate Genes #####
#write out the SNP positions, overlap those position with genes in bash 
write.table(tidy_fixed_targ_all %>% mutate(chr=paste0('chr_',chr),end=start) %>% ungroup %>% select(chr,start,end,Phylo) %>% arrange(chr,start),file='Fixed_SNP_Candidates_2024MAR20.txt',quote=F,sep='\t',row.names=F,col.names=F)
#intersect with genes
#bedtools intersect -a Fixed_SNP_Candidates_2024MAR20.txt -b Gene_Lookup.bed -wao | awk '{OFS="\t"}{print $1, $2, $3, $4, $6, $7,$8, $9}' > Fixed_SNP_Candidates_Genes_2024MAR20.txt

#load back in 
fixed_genes = read_tsv('Fixed_SNP_Candidates_Genes_2024MAR20.txt',col_names=F)
names(fixed_genes) = c('chr','start','d1','Phylo','GeneStart','GeneEnd','Strand','gene')
fixed_genes = fixed_genes %>% select(-d1) %>% mutate(chr = gsub('chr_','',chr))

#which genes are present in all phylo comparisons? count total SNPs 
overall = fixed_genes %>% ungroup %>% filter(gene != '.') %>% select(Phylo,chr,gene) %>% 
  group_by(chr,Phylo,gene) %>% summarize(snps = n()) %>% ungroup %>% group_by(chr,gene) %>% 
  summarize(Group = paste0(Phylo,collapse='__'), total_snps = sum(snps),Groups = n()) %>% ungroup %>% slice_max(Groups)
#on Z 
ztargs = fixed_genes %>% ungroup %>% filter(gene != '.') %>% select(Phylo,chr,gene) %>% 
  group_by(chr,Phylo,gene) %>% summarize(snps = n()) %>% ungroup %>% group_by(chr,gene) %>% 
  summarize(Group = paste0(Phylo,collapse='__'), total_snps = sum(snps),Groups = n()) %>% ungroup %>% filter(chr == 'Z') %>% slice_max(Groups)
#autosomes
auto_targs = fixed_genes %>% ungroup %>% filter(gene != '.') %>% select(Phylo,chr,gene) %>% 
  group_by(chr,Phylo,gene) %>% summarize(snps = n()) %>% ungroup %>% group_by(chr,gene) %>% 
  summarize(Group = paste0(Phylo,collapse='__'), total_snps = sum(snps),Groups = n()) %>% ungroup %>% filter(!grepl('Z|MT|W',chr)) %>% slice_max(Groups)

#for plot, find which genes have AT LEAST 3 phylogenetic contrasts with a fixed SNP 
fin_targs = fixed_genes %>% ungroup %>% filter(gene != '.') %>% select(Phylo,chr,gene) %>% 
  group_by(chr,Phylo,gene) %>% summarize(snps = n()) %>% ungroup %>% group_by(chr,gene) %>% 
  summarize(Group = paste0(Phylo,collapse='__'), total_snps = sum(snps),Groups = n()) %>% ungroup %>% filter(Groups > 2) 

#plot those 
gene_candidates = fixed_genes %>% filter(gene %in% fin_targs$gene) %>% 
  mutate(gene = toupper(gsub('ID=','',gene)),
         Phylo = factor(Phylo,levels=c('Contemporary_Diversification','Contemporary_Reversion','Ancient_Reversion','Ancient_Blue'))) %>% 
  select(chr,Phylo,gene) %>% unique %>% 
  group_by(gene) %>% 
  mutate(total = n()) %>% ungroup %>% 
  mutate(proportion = 1/total) %>% ungroup %>% 
  arrange(desc(total)) %>% 
  mutate(gene = factor(gene, levels = unique(gene))) %>%  
  ggplot(aes(x=gene,fill=Phylo,y=proportion))+
  geom_bar(stat='identity',position='stack')+
  facet_grid(.~chr,scales='free',space='free')+
  scale_y_continuous(breaks= pretty_breaks(n=3))+
  #scale_shape_manual(values=rev(c(16,15,18,17)))+
  scale_fill_manual(values=rev(c('cornflowerblue','plum','orange','indianred1')))+
  theme_bw(base_size=4)+
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1))

pdf('../figures/Fixed_SNPs_By_Phylo_2024MAR20.pdf',height=1.5,width=5)
gene_candidates
dev.off()

#grab the actual SNPs from the main positions 
inspection = fixed_genes %>% filter(gene %in% c(overall$gene,'NDUFS4')) %>% ungroup %>% select(chr,start,gene,Phylo) %>% unique %>% group_by(chr,gene,Phylo) %>% slice_max(start) %>% data.frame
inspect_bed = inspection %>% mutate(end = start,chr=paste0('chr_',chr)) %>% select(chr,start,end,Phylo,gene)
write.table(inspect_bed,file='Fixed_SNP_Candidates_Genes-PositionsForGenotypePlot_2024MAR20.txt',quote=F,sep='\t',row.names=F,col.names=F)

#in bash:
for chr in $(awk '{print $1}' Fixed_SNP_Candidates_Genes-PositionsForGenotypePlot_2024MAR20.txt | sort | uniq); do
echo "WORKING ON CHR: ${chr}"
bedtools intersect -header -a vcfs/${chr}.SNP.DP3-AC1-MQ40.vcf.gz -b Fixed_SNP_Candidates_Genes-PositionsForGenotypePlot_2024MAR20.txt | bcftools view -Oz -o plot_genos/${chr}.vcf.gz
done
bcftools concat -Oz -o plot_genos/Targets.vcf.gz plot_genos/chr_*vcf.gz

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

gtp = pg %>% 
  group_by(chr,gene,site,Genotype,Allele,ID,Hap,Egg) %>% 
  summarize(Group = paste0(Phylo,collapse='__')) %>% 
  ggplot(aes(y = ID, x = as.factor(site), fill = Genotype)) +
  geom_tile() +
  ylab("Individuals") +
  facet_grid(Hap+Egg~Group+chr+gene,space='free',scales='free')+
  #scale_fill_manual(values=c(brewer.pal(4,'Set2')[c(1,2)],'white'))+
  scale_fill_manual(values=c('grey70','grey20','white'))+
  theme_minimal(base_size=6) + ylab('')+xlab('') +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        #axis.text.x = element_text(size=4,angle = 90, vjust = 0.5, hjust = 1),
        strip.text.y.right = element_text(angle = 0),legend.position='top',
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
gtp

pdf('../figures/Genotypes-Targets_2024MAR20.pdf',height=2.5,width=6)
gtp
dev.off()

###### inspect the NDUFS4 area wher we have blue causative mutations FST=1  #####
tidy_bp %>% group_by(Phylo,Group,AvZ) %>% sum_stats(FST)
z_inspect = tidy_bp %>% filter(AvZ == 'Z' & start >= 7778296-1e4 & start <= 7783688+1e4)

#fixed SNPs?
tidy_fixed %>% filter(AvZ == 'Z' & start >= 7778296-1e4 & start <= 7783688+1e4) %>% count(Phylo)
tidy_fixed %>% filter(AvZ == 'Z' & start >= 7778296-1e4 & start <= 7783688+1e4) %>% count(Phylo)
z_inspect %>% ggplot(aes(x=start,y=FST,col=Phylo,group=Group))+
  geom_point()+
  geom_line()+
  scale_color_manual(values=c('cornflowerblue','plum','indianred1','orange','grey80'))+
  theme_bw()

##### CNV #####
#what's happening with log2FM there?
cnv = read_tsv('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/manyhost_hunt/coverage/CNV_log2_EggSex-ZW_500BP_2024MAR19.txt')
cnv = cnv %>% mutate(chr = gsub('chr_','',chr),
                     chr = ifelse(chr == 'Z' & start > 77e6,'W',chr),
                     start = ifelse(chr == 'W' & start > 77e6,start-77e6+22312972,start)) %>% 
  drop_na(log2)

#assign phylo to each egg and sex comparison 
egg_cnv = cnv %>% filter(comparison == 'Egg') %>% mutate(Group = paste0(p1,'__',p2)) %>% drop_na(log2)
egg_cnv = egg_cnv %>% mutate(Depth = ifelse(grepl('W1|W2|W3',p1) & !grepl('W1|W2|W3',p2),'Ancient',
                                            ifelse(grepl('W1|W2|W3',p2) & !grepl('W1|W2|W3',p1),'Ancient','Contemporary')),
                             Comparison = ifelse(grepl('W1|W2',Group) & Depth == 'Ancient','Blue',
                                                 ifelse(Depth == 'Ancient','Reversion',
                                                        ifelse(Group == 'W2__E1__W1__E1' | Group == 'W7__E6__W5__E6','Control',
                                                               ifelse((grepl('W1|W2|W3',p1) & grepl('W1|W2|W3',p2)) & Depth == 'Contemporary','Reversion','Diversification')))),
                             Phylo = paste0(Depth,'_',Comparison))
sex_cnv = cnv %>% filter(comparison != 'Egg') %>% mutate(Group = comparison) %>% drop_na(log2) %>% 
  mutate(Phylo = Group)

#region elevated in blue:
targ_start = 7778296-5000
targ_end = 7783688+5000

#region between E10 
targ_start = 7739546-5000
targ_end = 7739546+5000

log2 = sex_cnv %>% 
  filter(chr == 'Z'  & start >= targ_start & start <= targ_end) %>% 
  mutate(color = ifelse(comparison == 'Egg',paste0(p1,'_',p2),comparison),
         facet = ifelse(comparison == 'Egg','Egg','Sex')) %>% 
  ggplot(aes(x=start,y=log2,col=color))+
  geom_hline(yintercept=0,lty=3,col='grey50')+
  geom_line(lwd=0.5)+ylab('log2(F/M)')+xlab('Position Z Chromosome')+
  #geom_vline(xintercept=c(7778296,7783688),lty=2)+
  geom_vline(xintercept=c(7739546),lty=2)+
  scale_color_manual(values=brewer.pal(8,'Paired')[c(1:3,7)])+
  theme_bw(base_size=6)+theme(legend.position='top')

bysex = sex_cnv %>% 
  filter(chr == 'Z'  & start >= targ_start & start <= targ_end) %>% 
  mutate(color = ifelse(comparison == 'Egg',paste0(p1,'_',p2),comparison),
         facet = ifelse(comparison == 'Egg','Egg','Sex')) %>% 
  pivot_longer(c(cov1,cov2)) %>% mutate(name = gsub('cov1','Female',name),name=gsub('cov2','Male',name)) %>% 
  ggplot(aes(x=start,y=value,col=name,Group=color))+
  ylab('Normalized Coverage')+xlab('Position Z Chromosome')+
  geom_line(lwd=0.5)+
  #geom_vline(xintercept=c(7778296,7783688),lty=2)+
  geom_vline(xintercept=c(7739546),lty=2)+
  scale_color_manual(values=c('black','grey50'))+
  theme_bw(base_size=6)+theme(legend.position='top')

pdf('../figures/SexCoverage_ZRegion-NDUFS4-E10_2024MAR20.pdf',height=3,width=2)
ggarrange(log2,bysex,nrow=2)
dev.off()

#what's the average of male and female coverage? 
sex_cnv %>% filter(chr == 'Z'  & start >= 7778296 & start <= 7783688) %>% summarize(f = mean(cov1,na.rm=TRUE),m=mean(cov2,na.rm=TRUE))

#for log2FM, what are the biggest outlier windows for sex-differential coverage? 
cnv %>% filter(comparison != 'Egg' & chr == 'Z' & log2 > 0.5) %>% count(chr,start) %>% filter(n == 4)
cnv %>% 
  filter(comparison != 'Egg') %>% 
  filter(chr == 'Z'  & start >= 7778296-1e6 & start <= 7783688+1e6) %>% 
  group_by(comparison) %>% summarize(meanF = mean(cov1),meanM = mean(cov2))

##### Bootstrap FST within +/- 5KB of the NDUFS4 blue egg candidates ###### 
library(data.table)
#Function to perform the comparison for one set. compares value within target bin against background bin: DIFFERENCE! targ - background 
compare_snps <- function(target_snps, background_snps) {
  if (nrow(target_snps) == 0 || nrow(background_snps) == 0) {
    return(NA) #return NA if either group has no SNPs
  }
  
  target_snp <- target_snps[sample(.N, 1)]
  background_snp <- background_snps[sample(.N, 1)]
  return(target_snp$Value - background_snp$Value)
}

#format FST
set.seed(111)
bp_fst = tidy_bp %>% dplyr::rename(Value = FST) %>% drop_na(Value)
setDT(bp_fst)

#specify the region of interest 
targ_start = 7778296-5000
ttarg_end = 7783688+5000
targ_end-targ_start
fst_sub = bp_fst[chr == 'Z' & start > targ_start & start < targ_end]
# Initialize lists to store the results for each group
fst_boots = list()

#for FST 
for (gp in unique(fst_sub$Group)) {
  cat('Working on group: ',gp,'\n')
  target_fst = fst_sub[Group == gp]
  background_fst = bp_fst[chr == 'Z' & !(start > targ_start & start < targ_end) & Group == gp]
  tcomparisons = replicate(10000, compare_snps(target_fst, background_fst)) # perform comparisons
  fst_boots[[gp]] = data.table(Group = gp, Difference = tcomparisons, Variable = 'FST', Phylo = unique(target_fst$Phylo)) #store results for this  group
}

fst_fin = rbindlist(fst_boots) %>% as_tibble

#plot differences
boot_fst_plot = fst_fin %>%
  mutate(Phylo = factor(Phylo, levels = c('Ancient_Blue', 'Ancient_Reversion', 'Contemporary_Reversion', 'Contemporary_Diversification','Contemporary_Control'))) %>%
  group_by(Variable, Group, Phylo) %>%
  sum_stats(Difference) %>%
  ggplot(aes(y=Phylo,x=mean,xmin=conf_low,xmax=conf_high,col=Phylo,Group=Group))+
  geom_point(position=position_dodge(width=0.25),size=0.5)+
  geom_errorbar(position=position_dodge(width=0.25),width=0,lwd=0.25)+ylab('')+xlab('FST Difference Within 5-Kb Z-chr Blue SNPs')+
  scale_color_manual(values=c('cornflowerblue','plum','orange','indianred1','grey80'))+
  geom_vline(xintercept=0,lty=2)+
  theme_bw(base_size=6)+
  theme(legend.position='top')
boot_fst_plot
write.table(fst_fin,file='Bootstrapped_NDFUS4-End_2024MAR20.pdf',quote=F,sep='\t',row.names=F)

pdf('../figures/BootFST_ZRegion_2024MAR20.pdf',height=2,width=2)
boot_fst_plot
dev.off()


##### Deprecated, plot genome-wide ###### 
#prep karyoplot
genome = read.table('/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa.bed',header=FALSE) %>% filter(str_detect(V1,'scaff|MT',negate=T)) %>% arrange(desc(V2))
names(genome) = c('chr','start','end'); genome$chr = gsub('chr_','',genome$chr); genome = genome %>% filter(!grepl('3[0-9]|MT',chr))
genome = genome %>% mutate(chr = factor(chr, levels = c(as.character(1:29), "Z", "W"))) %>%  arrange(chr)
genome = genome %>% mutate(end = ifelse(chr == 'W',end+1326000,end))
G = makeGRangesFromDataFrame(genome,seqnames.field = 'chr',start.field = 'start',end.field = 'end')

#for only plotting Z/W/MT
genome = read.table('/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa.bed',header=FALSE) %>% filter(str_detect(V1,'scaff|MT',negate=T)) %>% arrange(desc(V2))
names(genome) = c('chr','start','end'); genome$chr = gsub('chr_','',genome$chr); 
genomesub = genome %>% filter(grepl('MT|W|Z',chr)) 
genome = genome %>% mutate(end = ifelse(chr == 'W',end+1326000,end))
genomesub = rbind(genomesub,data.frame(chr='MT',start=0,end=19698))
G = makeGRangesFromDataFrame(genomesub,seqnames.field = 'chr',start.field = 'start',end.field = 'end')

#add colors 
chrs = genome %>% select(chr) %>% mutate(Color = if_else(row_number() %% 2 == 0, 'grey60', 'black'))


# General karyoplot parameters
pdf('../figures/KARYOPLOT_correctedFST-All-VCFTOOLS-250KB-Labs_2024MAR06.pdf',height=3,width=7.5)
png('../figures/KARYOPLOT_correctedFST-All-1KB_2024MAR20.png',units='in',res=600,height=3,width=7.5)
pp = getDefaultPlotParams(plot.type=4)
pp$leftmargin = 0.1
kp = plotKaryotype(plot.type=4, genome = G,
                   labels.plotter = NULL,
                   plot.params = pp)
kpAddChromosomeNames(kp, yoffset = -5,cex=0.7)
kpAddBaseNumbers(kp,tick.dist = 50000000,minor.ticks = FALSE,cex=0.4)

#Loop through each track (all species together, CC, and CO), and plot them on their own layer 
counter = 0; tracks = 4

#add fst, for each pop
pops = c('Ancient_Blue','Ancient_Reversion','Contemporary_Reversion','Contemporary_Diversification')
leg=NULL
for (pop in pops) { 
  counter = counter + 1; cat('Current Autotrack: ',counter,' for pop: ',pop,'\n'); at = autotrack(current.track = counter, total.tracks = tracks);at$r1 <- at$r1-0.03
  #subset group
  #d = tidy_all_plot %>% filter(Phylo == pop) %>% ungroup %>% arrange(chr,start,Group)
  
  if(pop == 'Ancient_Blue') {
    col = 'cornflowerblue' } else if(pop =='Ancient_Reversion') {
      col = 'plum' } else if(pop =='Contemporary_Diversification') {
        col = 'indianred1' } else if(pop =='Contemporary_Reversion') {
          col = 'orange' } else { }
  
  ds = tidy_avg %>% filter(Phylo == pop) %>% ungroup %>% arrange(chr,start)
  
  #mean FST
  minv = min(tidy_avg$mean,na.rm=TRUE) ; maxv = max(tidy_avg$mean,na.rm=TRUE)
  kpLines(kp,chr=ds$chr,x=ds$start,y=ds$mean,col=col,ymin=minv,ymax=maxv,r0=at$r0,r1=at$r1-0.06)
  
  # #mean zFST
  # minv = min(tidy_avg$meanZ,na.rm=TRUE) ; maxv = max(tidy_avg$meanZ,na.rm=TRUE)
  # kpLines(kp,chr=ds$chr,x=ds$start,y=ds$meanZ,col=col,ymin=minv,ymax=maxv,r0=at$r0,r1=at$r1-0.06)
  
  #mean zFST_Floor
  # minv = min(tidy_avg$meanZf,na.rm=TRUE) ; maxv = max(tidy_avg$meanZf,na.rm=TRUE) 
  # kpLines(kp,chr=ds$chr,x=ds$start,y=ds$meanZf,col=col,ymin=minv,ymax=maxv,r0=at$r0,r1=at$r1-0.06)
  # 
  
  
  # #also add the FST = 1 targets, colors for rects 
  # if(pop == 'Ancient_Blue') { 
  #   col = 'cornflowerblue' ; d = d %>% mutate(col = brewer.pal(n = n_distinct(Group), name = "Blues")[match(Group, unique(Group))]) } else if(pop =='Ancient_Reversion') { 
  #     col = 'plum' ; d = d %>% mutate(col = brewer.pal(n = n_distinct(Group), name = "Purples")[match(Group, unique(Group))]) } else if(pop =='Contemporary_Diversification') { 
  #       col = 'indianred1' ;  d = d %>% mutate(col = brewer.pal(n = n_distinct(Group), name = "Reds")[match(Group, unique(Group))]) } else if(pop =='Contemporary_Reversion') { 
  #         col = 'orange' ;  d = d %>% mutate(col = brewer.pal(n = n_distinct(Group), name = "Oranges")[match(Group, unique(Group))]) } else { col ='black' ; d = d %>% mutate(col = brewer.pal(n = n_distinct(Group), name = "Greys")[match(Group, unique(Group))])}
  # 
  # leg = rbind(leg, d %>% ungroup %>% select(Phylo,Group,col) %>% unique)
  # #fixedtarg = plot_fixed %>% filter(Phylo == pop)
  # #fixed_names = plot_fixed %>% filter(Phylo == pop) %>% select(chr,start,Phylo,gene) %>% unique %>% mutate(gene = gsub('gene-','',gene))
  # minv = 0 ; maxv = 1
  # 
  # for (group in unique(d$Group)) { 
  #   #plot zfst and axis 
  #   ds = d %>% filter(Group == group)
  #   kpLines(kp,chr=ds$chr,x=ds$start,y=ds$correctedFST,col=ds$col,ymin=minv,ymax=maxv,r0=at$r0,r1=at$r1-0.06)
  #   
  # }
  kpAxis(kp,r0=at$r0,r1=at$r1,cex=0.3,numticks = 2,ymin =minv,ymax=maxv)
  kpAddLabels(kp,cex=0.25,labels = pop,r0=at$r0+.01, r1=at$r1,col="black",srt=0,label.margin = 0.02)
  
  # if(nrow(fixed_names) == 0) next
  #   kpText(kp, chr=fixed_names$chr, x=fixed_names$start, y=1, labels=fixed_names$gene,srt=45,cex=0.2,r0=at$r1-0.15,r1=at$r1-0.05)
  # 
  # #also show FST 1 points
  # if(nrow(fixedtarg) == 0) next
  #   kpPoints(kp,chr=fixedtarg$chr,x=fixedtarg$start,pch=8,y=0.5,cex=0.25,lwd=0.25,col=col,r0=at$r1+0.01,r1=at$r1+0.03)
  # 
}

kpPoints(kp,chr='Z',x=7778296,y=1,r0=0,r1=1,col='black',pch=8,cex=2)

dev.off()
