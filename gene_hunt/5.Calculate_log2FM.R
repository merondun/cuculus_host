#### Coverage discrepancies, CNV, calculate log2(Immaculate/Spotted)
#setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/manyhost_hunt')
#.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(viridis)
library(zoo)

cov = read.table('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/manyhost_hunt/coverage/CNV_HappEgg_chrZW_500BP_2024MAR19.txt',header=FALSE)

names(cov) = c('chr','start','end','coverage','ID')
cov = cov %>% mutate(mtG = ifelse(chr == 'chr_MT','MT','nonMT'))

#add grouping
md = read_tsv('~/merondun/cuculus_host/Metadata_Host.txt')
md = md %>% mutate(Group = paste0(Hap,'__',gsub('. ','',Egg)))

covs = left_join(cov,md %>% select(ID,Mapped_Reads_M))

#normalize according to the number of mapped reads 
norm_cov = covs %>%
  mutate(overall_median_coverage = median(Mapped_Reads_M)) %>% #median coverage overall
  mutate(normalization_factor = overall_median_coverage / Mapped_Reads_M) %>% #normalization for each sample
  mutate(normalized_coverage = coverage * normalization_factor) %>% #apply that to each window by sample
  dplyr::select(chr, start, end, normalized_coverage, ID) %>% ungroup

nc = left_join(norm_cov,md %>% select(ID,Group,Sex))
nc = nc %>% mutate(site = paste0(chr,'_',start))

#confirm it worked
fmp = nc %>% group_by(ID,Group,Sex,chr) %>% summarize(median = median(normalized_coverage)) %>%
  ggplot(aes(y=Sex,x=median,fill=Sex))+geom_boxplot()+facet_grid(chr~.,scales='free',space='free')+
  theme_bw(base_size=6)+ylab('')+xlab('Median Coverage')+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3))
fmp
pdf('../figures/CNV_Normalization-500BP_2024MAR19.pdf',height=4,width=3)
fmp
dev.off()

# Generate bootstrap support for each window
bootin = nc %>% dplyr::select(chr,start,site,Group,Sex,normalized_coverage) %>% mutate_at(c('Group','Sex'),as.factor)

#group by chr/start, then summarize to get median coverage by Group and Sex
#comparisons: awk '{print $1, $2}' chr_8.20KB_fst.txt  | tr ' ' '\t' | sort | uniq > ../coverage/Group_Comparisons.list
groups = read_tsv('coverage/Group_Comparisons.list')
groups = groups %>% mutate(pop1 = gsub('Pphoenicurus','E1',pop1),
                           pop1 = gsub('Aarundinaceus','E6',pop1),
                           pop1 = gsub('Phaedon','E10',pop1),
                           pop2 = gsub('Pphoenicurus','E1',pop2),
                           pop2 = gsub('Aarundinaceus','E6',pop2),
                           pop2 = gsub('Phaedon','E10',pop2))
groups = groups %>% filter(!grepl('E10|W5__E6',pop1) & !grepl('E10|W5__E6',pop2))
coverage_summary = bootin %>%
  group_by(chr, start, Group, Sex) %>%
  summarize(median_coverage = median(normalized_coverage, na.rm = TRUE), .groups = 'drop')

#calculate CNV in egg type. Round to the nearest integer, and replace 0 with 0.5 so that comparisons between 1 / 0 do not turn infinite, and so that those comparisons are not huge like 1 / 0.001
# if there is 0 / 0 coverage, then consider that a true NA 
eggcov = NULL
for (i in seq(1,nrow(groups))) {
  p = groups[i,]
  p1 = p$pop1
  p2 = p$pop2
  input = coverage_summary %>% filter((Group == p1 | Group == p2) & Sex == 'F') %>% 
    pivot_wider(names_from = Group,values_from = median_coverage) %>% select(-Sex)
  cat('Working on comparison between: ',p1,' and ',p2,'\n')
  egglog = input %>% mutate(c1 = pmax(round(!!sym(p1),0),0.5),
                            c2 = pmax(round(!!sym(p2),0),0.5),
                            log2 = ifelse(c1 == 0.5 & c2 == 0.5,NA,
                                          log2(c1 / c2))) %>% 
    select(chr,start,log2,cov1=!!sym(p1),cov2=!!sym(p2)) %>% mutate(p1 = p1, p2 = p2, comparison = 'Egg')
  eggcov = rbind(eggcov,egglog)
  
}

#now compare WITHIN each group, for sex 
pops = unique(c(groups$pop1,groups$pop2))
sexcov = NULL
for (p in unique(pops)) {
  input = coverage_summary %>% filter(Group == p)  %>% 
    pivot_wider(names_from = Sex,values_from = median_coverage)
  cat('Working on comparison between sexes for population: ',p,'\n')
  sexlog = input %>% mutate(c1 = pmax(round(F,0),0.5),
                            c2 = pmax(round(M,0),0.5),
                            log2 = ifelse(c1 == 0.5 & c2 == 0.5,NA,
                                          log2(c1 / c2))) %>% 
    select(chr,start,log2,cov1=F,cov2=M) %>% mutate(p1 = 'F', p2 = 'M', comparison = p)
  sexcov = rbind(sexcov,sexlog)
}

#final frame
log2 = rbind(eggcov,sexcov)
log2 %>% count(p1,p2,comparison)

write.table(log2,'/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/manyhost_hunt/coverage/CNV_log2_EggSex-ZW_500BP_2024MAR19.txt',quote=F,sep='\t',row.names=F)
