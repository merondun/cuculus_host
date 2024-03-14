#### calculate log2FM
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/immaculate_hunting/gwas_mtdna/coverage') 
.libPaths('~/mambaforge/envs/R/lib/R/library')
library(pacman)
p_load(tidyverse,ggpubr,RColorBrewer,viridis,zoo,scales,meRo)

#add f/m
fm = read_tsv('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/immaculate_hunting/gwas_mtdna/coverage/FM_500bpMQ20.txt',col_names = F)
names(fm) = c('chr','start','end','coverage','ID')
samples = unique(fm$ID)
fmd = read_tsv('SamplesN5.txt')
fmm = left_join(fm,fmd %>% dplyr::select(ID,Sex,Group))
fmm = fmm %>% mutate(mtG = ifelse(chr == 'chr_MT','MT','nonMT'))

#normalize FM
normalized_fm <- fmm %>%
  group_by(ID,mtG) %>%
  mutate(median_coverage_per_id = median(coverage)) %>% #median coverage by ID
  ungroup() %>%
  group_by(mtG) %>% 
  mutate(overall_median_coverage = median(median_coverage_per_id)) %>% #median coverage overall
  mutate(normalization_factor = overall_median_coverage / median_coverage_per_id) %>% #normalization for each sample
  mutate(normalized_coverage = coverage * normalization_factor) %>% #apply that to each window by sample
  dplyr::select(chr, start, end, normalized_coverage, ID, Sex, Group) %>% ungroup

fmp = normalized_fm %>% group_by(ID,mtG,Group,Sex) %>% summarize(median = median(normalized_coverage)) %>%
  ggplot(aes(y=ID,x=median,col=mtG))+geom_point()+facet_grid(Group+Sex~mtG,scales='free',space='free')+
  theme_bw(base_size=6)+ylab('')+xlab('Median Coverage')+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3))
pdf('../../../figures/FM_Normalization_2024FEB13.pdf',height=5,width=4)
fmp
dev.off()

#calculate f/m by group
fms1 = normalized_fm %>%
  group_by(chr, start, end, Group) %>% 
  #Pre-compute median coverages for F and M
  summarize(median_coverage_F = median(normalized_coverage[Sex == 'F'],na.rm=TRUE),
            median_coverage_M = median(normalized_coverage[Sex == 'M'],na.rm=TRUE))

#replace 0 with 0.1 for division
fms1 = fms1 %>% mutate(across(starts_with('median_'), ~ replace(., . == 0, 0.1)))

#If coverage is 0 in both, make it NA, otherwise calculate log2FM
fms2 = fms1 %>%
  mutate(
    FM = ifelse(median_coverage_F == 0.1 & median_coverage_M == 0.1, NA_real_,
                median_coverage_F / median_coverage_M),
    log2FM = log2(FM),
    log2FM = pmax(pmin(log2FM, 4), -4)  # Apply the limits
  ) %>%
  ungroup()

#plot the output log2fm by chr 
cols = brewer.pal(8,'Paired')
fp = fms2 %>% group_by(chr,Group) %>% 
  filter(!grepl('scaff|chr_3[0-9]',chr)) %>% 
  sum_stats(log2FM) %>%
  mutate(chr = gsub('chr_','',chr)) %>% 
  ggplot(aes(y=chr,x=median,xmin=median-iqr,xmax=median+iqr,col=Group))+
  scale_color_manual(values=c(cols[1],cols[2],cols[3],'black'))+
  geom_point()+ xlab('Median log2(F/M)')+ylab('')+
  theme_bw()
pdf('../../../figures/FM_log2FM-chromosome_2024FEB13.pdf',height=6,width=6)
fp
dev.off()

#summary stats on W 
fms2 %>% group_by(chr,Group) %>% 
  filter(!grepl('scaff|chr_3[0-9]',chr)) %>% 
  sum_stats(log2FM)  %>% filter(grepl('W',chr))

write.table(fms2,'/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/immaculate_hunting/gwas_mtdna/coverage/log2FM_coverageMedian_2024FEB13.txt',quote=F,sep='\t',row.names=F)
