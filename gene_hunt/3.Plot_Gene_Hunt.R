#### Plot FST, 250KB and BP
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/manyhost_hunt/')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(viridis)
library(ggpubr)
library(meRo)
library(RColorBrewer)
library(karyoploteR)
library(zoo)

#prep karyoplot
genome = read.table('/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa.bed',header=FALSE) %>% filter(str_detect(V1,'scaff|MT',negate=T)) %>% arrange(desc(V2))
names(genome) = c('chr','start','end'); genome$chr = gsub('chr_','',genome$chr); genome = genome %>% filter(!grepl('3[0-9]|MT',chr))
genome = genome %>% mutate(chr = factor(chr, levels = c(as.character(1:29), "Z", "W"))) %>%  arrange(chr)
genome = genome %>% mutate(end = ifelse(chr == 'W',end+1326000,end))
G = makeGRangesFromDataFrame(genome,seqnames.field = 'chr',start.field = 'start',end.field = 'end')

#add colors 
chrs = genome %>% select(chr) %>% mutate(Color = if_else(row_number() %% 2 == 0, 'grey60', 'black'))

#start here 
tidy_all = read_tsv('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/manyhost_hunt/ManyHost_FST_250KB_2024MAR08.txt')
names(tidy_all) = c('p1','p2','chr','start','end','FST','snps')

#BUT, fix misassembly on the Z
tidy_all = tidy_all  %>% 
  drop_na(FST) %>% 
  mutate(
    chr = gsub('chr_','',chr),
    chr = ifelse(chr == 'Z' & start > 77e6,'W',chr),
    start = ifelse(chr == 'W' & start > 77e6,start-77e6+22312972,start),
    end = ifelse(chr == 'W' & end > 77e6,end-77e6+22312972,end),
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
  ) %>% 
  left_join(.,chrs)
tidy_all %>% select(Phylo,Group) %>% unique %>% data.frame
tidy_all %>% count(p1,p2,Phylo)

#okay, so now loop through each group, and identify which regions are elevated compared to the genomic background
#first, calculate zfst 
tidy_all = tidy_all %>% 
  group_by(p1, p2) %>%
  mutate(
    mean_FST = mean(FST, na.rm = TRUE),
    sd_FST = sd(FST, na.rm = TRUE),
    zFST = (FST - mean_FST) / sd_FST
  ) %>%
  ungroup() %>%
  select(-mean_FST, -sd_FST) 

#subtract FST from the mean control FST for plotting 
control20 = tidy_all %>% ungroup %>% filter(Phylo == 'Contemporary_Control') %>% group_by(chr,start) %>% summarize(controlFST = max(FST,na.rm=TRUE))
tidy_all_plot = left_join(tidy_all %>% ungroup %>% filter(Phylo != 'Contemporary_Control') %>% select(chr,start,AvZ,Phylo,Group,FST),control20) %>% mutate(correctedFST = pmax(0,FST - controlFST))
summary(tidy_all_plot$correctedFST)

#Analyze FST base-pair
tidy_bp = read_tsv('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/manyhost_hunt/ManyHost_FST_1BP-VCFTOOLS_2024MAR12.txt',col_names=FALSE)
names(tidy_bp) = c('p1','p2','chr','start','FST')

#BUT, fix misassembly on the Z
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

#grab only the fixed sites
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
#all comparisons fixed
tidy_fixed_targ_all = tidy_fixed_targ %>% group_by(Phylo,AvZ,chr,start) %>% count %>%
  filter(
    ifelse(Phylo == 'Ancient_Blue', n == 6, 
           ifelse(Phylo == 'Ancient_Reversion', n == 3,
                  ifelse(Phylo == 'Contemporary_Diversification', n == 2,
                         ifelse(Phylo == 'Contemporary_Reversion', n == 2, n >= 999)))))

#total windows
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

#underlying blue
tidy_fixed_targ_all %>% ungroup %>% select(Phylo,chr,start,n) %>% pivot_wider(names_from='Phylo',values_from='n') %>% filter(!is.na(Ancient_Blue) & !is.na(Contemporary_Reversion))

#underlying ancient blue 
blu_anc = tidy_fixed_targ_all %>% ungroup %>% select(Phylo,chr,start,n) %>% pivot_wider(names_from='Phylo',values_from='n') %>% filter(!is.na(Ancient_Blue) & !is.na(Ancient_Reversion)) 
blu_anc

#plot
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
#props

pdf('../figures/FST-BP_ManyHosts_Proportions-All_2024MAR06.pdf',height=3,width=2.5)
ggarrange(propsbottom,propstop,widths=c(0.7,0.3))
dev.off()

#make our final targets those which are present in at least half 
final_targets = tidy_fixed_targ_all
tidy_fixed_targ_all %>% ungroup %>% count(Phylo,AvZ)
tidy_fixed_targ_all %>% ungroup %>%
  select(chr,AvZ,start) %>% unique %>% count(chr,AvZ,start) %>%
  filter(n > 2)

#are there regions of MT/W which contain more differentiated windows?
dens = final_targets %>% filter((chr == 'W' | chr == 'MT')) %>%
  ggplot(aes(x=start,fill=Phylo))+
  scale_fill_manual(values=c('cornflowerblue','plum','indianred1','orange'))+
  geom_density(alpha=0.8)+facet_wrap(chr~.,scales='free')+theme_bw(base_size=6) +theme(legend.position='top')
#dens
pdf('../figures/FST-BP_ManyHosts_Density-All--W-MT_2024MAR06.pdf',height=3,width=5)
dens
dev.off()

#merge with syn/nonyn
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
#num genes fixed with dN by chr 
plot_fixed %>% select(chr,Phylo,gene) %>% unique %>% count(Phylo,chr)

#save
write.table(fbp,file='FST-BP_dNdS_Overlap_MajorityDiverged-All_2024MAR11.txt',quote=F,sep='\t',row.names=F)

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

# General karyoplot parameters
pdf('../figures/KARYOPLOT_correctedFST-All-VCFTOOLS-250KB-Labs_2024MAR06.pdf',height=3,width=7.5)
#png('../figures/KARYOPLOT_correctedFST-All-20KB_2024MAR06.png',units='in',res=600,height=3,width=7.5)
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
pops = tidy_all_plot %>% ungroup %>% select(Phylo) %>% unique
pops = pops$Phylo[c(2,3,1,4)]
for (pop in pops) { 
  counter = counter + 1; cat('Current Autotrack: ',counter,' for pop: ',pop,'\n'); at = autotrack(current.track = counter, total.tracks = tracks);at$r1 <- at$r1-0.03
  #subset group
  d = tidy_all_plot %>% filter(Phylo == pop) %>% ungroup %>% arrange(chr,start,Group)
  
  #also add the FST = 1 targets, colors for rects 
  if(pop == 'Ancient_Blue') { 
    col = 'cornflowerblue' ; d = d %>% mutate(col = brewer.pal(n = n_distinct(Group), name = "Blues")[match(Group, unique(Group))]) } else if(pop =='Ancient_Reversion') { 
      col = 'plum' ; d = d %>% mutate(col = brewer.pal(n = n_distinct(Group), name = "Purples")[match(Group, unique(Group))]) } else if(pop =='Contemporary_Diversification') { 
        col = 'indianred1' ;  d = d %>% mutate(col = brewer.pal(n = n_distinct(Group), name = "Reds")[match(Group, unique(Group))]) } else if(pop =='Contemporary_Reversion') { 
          col = 'orange' ;  d = d %>% mutate(col = brewer.pal(n = n_distinct(Group), name = "Oranges")[match(Group, unique(Group))]) } else { col ='black' ; d = d %>% mutate(col = brewer.pal(n = n_distinct(Group), name = "Greys")[match(Group, unique(Group))])}
  
  fixedtarg = plot_fixed %>% filter(Phylo == pop)
  fixed_names = plot_fixed %>% filter(Phylo == pop) %>% select(chr,start,Phylo,gene) %>% unique %>% mutate(gene = gsub('gene-','',gene))
  minv = 0 ; maxv = 1
  #also show FST 1 points
  kpPoints(kp,chr=fixedtarg$chr,x=fixedtarg$start,pch=8,y=0.5,cex=0.25,lwd=0.25,col=col,r0=at$r1+0.01,r1=at$r1+0.03)

  for (group in unique(d$Group)) { 
    #plot zfst and axis 
    ds = d %>% filter(Group == group)
    kpLines(kp,chr=ds$chr,x=ds$start,y=ds$correctedFST,col=ds$col,ymin=minv,ymax=maxv,r0=at$r0,r1=at$r1-0.06)

  }
  kpAxis(kp,r0=at$r0,r1=at$r1,cex=0.3,numticks = 2,ymin =minv,ymax=maxv)
  kpAddLabels(kp,cex=0.25,labels = pop,r0=at$r0+.01, r1=at$r1,col="black",srt=0,label.margin = 0.02)
 
  if(nrow(fixed_names) == 0) next
    kpText(kp, chr=fixed_names$chr, x=fixed_names$start, y=1, labels=fixed_names$gene,srt=45,cex=0.2,r0=at$r1-0.15,r1=at$r1-0.05)
 
}

dev.off()

plot_fixed$Phylo = factor(plot_fixed$Phylo,levels=c('Contemporary_Diversification','Contemporary_Reversion','Ancient_Reversion','Ancient_Blue'))
#plot all targets
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


pdf('../figures/dN_Mutations_2024MAR11.pdf',height=1.5,width=4.5)
p
dev.off()


#genes with the most dN mutations
plot_fixed %>% ungroup %>% count(Phylo,gene) %>% arrange(desc(n)) %>% filter(n > 1) %>% ungroup %>% count(gene) %>% filter(n > 1) %>% data.frame

#underlying blue
plot_fixed %>% ungroup %>% select(Phylo,chr,start,n) %>% pivot_wider(names_from='Phylo',values_from='n') %>% filter(!is.na(Ancient_Blue) & !is.na(Contemporary_Reversion))
blue = plot_fixed %>% ungroup %>% select(Phylo,chr,start,n) %>% pivot_wider(names_from='Phylo',values_from='n') %>% filter(!is.na(Ancient_Blue) & !is.na(Ancient_Reversion))

#any fixed SNPs near the nonsense pseudogene? 
tidy_fixed %>% filter(chr == 'W' & site > 6218150 & site < 6229021)

#lots of SNPs 
tidy_fixed %>% filter(chr == 'W' & start > 6218150 & start < 6229021) %>% count(Phylo)
## A tibble: 5 × 2
#  Phylo                            n
#  <chr>                        <int>
#1 Ancient_Blue                   283
#2 Ancient_Reversion              140
#3 Contemporary_Control             7
#4 Contemporary_Diversification     6
#5 Contemporary_Reversion           1

#any more than background? 
#grep 'chr_W' GCA_017976375.1_bCucCan1.pri_genomic.CHR-chr_MT.gff | awk '$3 == "pseudogene" || $3 == "gene"' | awk '{OFS="\t"}{print $1, $4, $5, $7, $3, $9}' | sed 's/;.*//g' | sed 's/ID=gene-//g' > ~/symlinks/host/manyhost_hunt/W_Genes_Pseudogenes.bed
genes = read_tsv('W_Genes_Pseudogenes.bed',col_names=F)
names(genes) = c('chr','start','end','strand','category','gene')
genes = genes %>% mutate(chr = gsub('chr_','',chr))
fst = tidy_fixed_targ_all %>% filter(chr == 'W')
setDT(fst)

genedat = NULL; alldat = NULL
for (gid in unique(genes$gene)){
coords = wgenes %>% filter(gene == gid)
gsub = fst[start > coords$start & start < coords$end]
props = gsub %>% count(Phylo) %>% mutate(length = coords$end - coords$start, proportion = n/length,gene=gid,category = coords$category)
genedat = rbind(genedat,props)
alldat = rbind(alldat,gsub %>% mutate(gene = gid))
}

#which genes fall inthe top 5% of the mutational distribution in 2 or more phylo comparisons ?
genedat %>% filter(!grepl('Control',Phylo)) %>% group_by(Phylo) %>% mutate(q = quantile(proportion,0.95)) %>% filter(proportion > q) %>% data.frame
phytarg = genedat %>% filter(!grepl('Control',Phylo)) %>% group_by(Phylo) %>% mutate(q = quantile(proportion,0.95)) %>% filter(proportion > q) %>% data.frame %>% group_by(gene) %>% mutate(CountPhylo = n()) %>% ungroup %>% filter(CountPhylo > 1 )
phytarg = genedat %>% filter(!grepl('Control',Phylo)) %>% group_by(Phylo) %>% mutate(q = quantile(proportion,0.95)) %>% filter(proportion > q) %>% data.frame %>% arrange(gene) %>% ungroup %>% count(gene,category) %>% arrange(desc(n)) %>% filter(n > 2)
phytarg
# A tibble: 11 × 8
#   Phylo                   n length proportion gene  category       q CountPhylo
#   <chr>               <int>  <dbl>      <dbl> <chr> <chr>      <dbl>      <int>
# 1 Ancient_Blue           16   4930   0.00325  LOC1… gene     2.83e-3          2
# 2 Ancient_Reversion      17   4930   0.00345  LOC1… gene     2.88e-3          2
# 3 Ancient_Blue           12   2361   0.00508  LOC1… gene     2.83e-3          3
# 4 Ancient_Reversion      13   2361   0.00551  LOC1… gene     2.88e-3          3
# 5 Contemporary_Rever…     1   2361   0.000424 LOC1… gene     2.06e-4          3
# 6 Ancient_Blue           13   3468   0.00375  LOC1… gene     2.83e-3          2
# 7 Ancient_Reversion      13   3468   0.00375  LOC1… gene     2.88e-3          2
# 8 Ancient_Blue           43  10871   0.00396  LOC1… pseudog… 2.83e-3          2
# 9 Ancient_Reversion      43  10871   0.00396  LOC1… pseudog… 2.88e-3          2
#10 Ancient_Blue           51  17889   0.00285  LOC1… pseudog… 2.83e-3          2
#11 Ancient_Reversion      52  17889   0.00291  LOC1… pseudog… 2.88e-3          2

#find the sites from those 
phytarg_genes = genes %>% filter(gene %in% phytarg$gene)

#intersect with fixed data
left_join(plot_fixed %>% select(chr,start,Phylo),alldat %>% filter(gene %in% phytarg$gene))
gene_mutation_targets = alldat %>% filter(gene %in% phytarg$gene) %>% group_by(chr,start,gene) %>% summarize(Group = paste0(Phylo,collapse='__')) %>%
  mutate(chr = paste0('chr_',chr),end = start) %>% select(chr,start,end,gene,Group)

#plot genotypes
plot_genos = plot_fixed %>% filter(grepl('nad2|nad4|LOC128850275',gene)) %>% group_by(chr,start,gene) %>% summarize(Group = paste0(Phylo,collapse='__')) %>% 
  mutate(chr = paste0('chr_',chr),end = start) %>% select(chr,start,end,gene,Group)
write.table(plot_genos,file='plot_genos/Targets.bed',quote=F,sep='\t',row.names=F,col.names=F)

#in bash:
for chr in $(awk '{print $1}' Targets.bed | sort | uniq); do
echo "WORKING ON CHR: ${chr}"
bedtools intersect -header -a ../vcfs/${chr}.SNP.DP3-AC1-MQ40.vcf.gz -b Targets.bed | bcftools view -Oz -o ${chr}.vcf.gz
done
bcftools concat -Oz -o Targets.vcf.gz chr_*vcf.gz

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
agt = agt %>% mutate(chr = ifelse(Key == 2,'W','MT'))

#add gene and group ID
gps = read_tsv('plot_genos/Targets.bed',col_names=F)
names(gps) = c('chr','site','end','gene','Group')
pg  = left_join(agt,gps %>% mutate(chr=gsub('chr_','',chr)))
pg = pg %>% mutate(gene = toupper(gsub('gene-','',gene)))

pg %>% count(chr,gene)
pg = pg %>% arrange(chr,site)
sites = pg %>% select(chr,site) %>% unique
pg$site = factor(pg$site,levels=sites$site)

gtp = pg %>% 
  ggplot(aes(y = ID, x = as.factor(site), fill = Genotype)) +
  geom_tile() +
  ylab("Individuals") +
  facet_grid(Hap+Egg~Group+chr+gene,space='free',scales='free')+
  scale_fill_manual(values=c(brewer.pal(4,'Set2')[c(1,2)],'white'))+
  #scale_fill_manual(values=c(viridis(5,option='turbo')[1:3],'white'))+
  theme_minimal(base_size=6) + ylab('')+xlab('') +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        #axis.text.x = element_text(size=4,angle = 90, vjust = 0.5, hjust = 1),
        strip.text.y.right = element_text(angle = 0),legend.position='top',
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

pdf('../figures/Genotypes-Targets_2024MAR11.pdf',height=2.5,width=5)
gtp
dev.off()

