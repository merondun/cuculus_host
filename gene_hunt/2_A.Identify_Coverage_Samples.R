#Subset samples for coverage calculations: N=3 from each population (based on W-D overlap) 
setwd('~/merondun/cuculus_host/gene_hunt/')
.libPaths('~/mambaforge/envs/R/lib/R/library')
library(tidyverse)
library(RColorBrewer)

set.seed(123)
md = read_tsv('~/merondun/cuculus_host/Metadata_Host.txt') 
md = md %>% mutate(Group = ifelse(Hap == 'W1','W1',
                                  ifelse(Hap == 'W2','W2',
                                         ifelse(Hap == 'W3','W3',
                                                ifelse(grepl('5|6|7',Hap),'Ws','W4')))),
                   GroupDist = paste0(Group,'__',KDist)) 
covsamp = md %>% filter(Analysis_Host == 1) %>% 
  mutate(avg = median(MedianCoverage,na.rm=TRUE)) %>% 
  dplyr::select(ID,Sex,KDist,Hap,Longitude,Latitude,GroupDist,Group,MedianCoverage,avg) %>% 
  mutate(diff = abs(MedianCoverage-avg)) %>% 
  group_by(Sex,GroupDist) %>% 
  slice_min(diff,n=2,with_ties = FALSE) %>% 
  #count(Sex,KDist,Hap) %>% arrange(GroupDist) %>% data.frame %>% filter(n > 1)
  filter(grepl('W1__D5|W2__D3|W3__D14|Ws__D12',GroupDist)) %>% ungroup 

cols=brewer.pal(8,'Paired')
cp = covsamp %>% ggplot(aes(x=Longitude,y=Latitude,col=Group,shape=Sex))+
  geom_point(size=3)+
  scale_color_manual(values=c(cols[1],cols[2],cols[3],'black'))+
  theme_bw()
cp
write.table(covsamp %>% select(ID,Sex,Group),file='~/merondun/cuculus_host/gene_hunt/Coverage_SamplesN2.txt',quote=F,sep='\t',row.names=F)