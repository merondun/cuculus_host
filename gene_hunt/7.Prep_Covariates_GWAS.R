#### Plot the p-distance matrix on the subset GWAS individuals 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/manyhost_hunt/males')
.libPaths('~/mambaforge/envs/R/lib/R/library')
library(tidyverse)
library(viridis)
library(ggpubr)
library(meRo)
library(RColorBrewer)
library(ggpubr)
library(scales)

md = read_tsv('~/merondun/cuculus_host/Metadata_Host.txt')
eggcols = md %>% ungroup %>% select(Egg) %>% unique %>% mutate(ord = as.numeric(gsub('E','',Egg))) %>% arrange(ord) %>% 
	  drop_na(Egg) %>% mutate(col = viridis(11,option='turbo'))
  eggcols = rbind(eggcols,data.frame(Egg='E6W3',ord=12,col='white'))

  #load in autosomal distance, calculated: ~/modules/VCF2Dis-1.50/bin/VCF2Dis -InPut autos_canorus_LD.vcf.gz -OutPut autos_canorus_LD.pdist 
  auto = read.table('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/merged_unrelfull/autos_canorus_LD.pdist',header=F)
  #add proper names, because VCF2DIS truncates them
  auto = auto %>% mutate(IDNum = gsub('_.*','',V1))
  auto_id = left_join(auto,md %>% select(ID) %>% mutate(IDNum = gsub('_.*','',ID)))  
  ord = auto_id$ID
  autos = auto_id %>% select(!c(IDNum,ID,V1))
  names(autos) = ord
  rownames(autos) = ord

  #the order we want, from the bed file 
  ids = read.table('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/manyhost_hunt/males/beds/chr_1.fam')

  #extract columns in same order as geographic distance
  common_names = intersect(rownames(autos),ids$V1)
  auto_aligned = autos[common_names, common_names] #subset 'auto' matrix to match the 'geo_mat' matrix rows and columns
  auto_mat = as.matrix(auto_aligned)

  #visualize autosomal distance matrix in terms of ancestry K = 5
  kcols = data.frame(Kcols = gray.colors(5,start=0.05,end=0.95)[c(1,4,3,5,2)], AncestryK5 = paste0('K',seq(1,5,1)))
  auto_mds = as.data.frame(cmdscale(auto_mat,2))
  auto_mds$ID = rownames(auto_mds)
  auto_mds_input = left_join(auto_mds,md)
  auto_p = ggplot(auto_mds_input,aes(x=-V1,y=V2,fill=AncestryK5,shape=AncestryK5)) + 
	    geom_point(size=1.5) + theme_bw() + xlab('MDS1')+ylab('MDS2')+ggtitle('Autosomal Distance')+
	      scale_fill_manual(name = "Ancestry K5", values = setNames(kcols$Kcols, kcols$Kcluster)) +
	        scale_shape_manual(name = "Ancestry K5", values = setNames(c(21, 22, 23, 24, 25), kcols$Kcluster)) +
		  theme(legend.position = "top", # Moves legend to top
			        legend.text = element_text(size = 6), # Adjusts text size
				        legend.title = element_text(size = 6))
  auto_p

  pdf('../../figures/GWAS_Covariance_Matrix_2024APR3.pdf',height=3,width=3)
  auto_p
  dev.off()

  write.table(auto_mat,file='beds/Covariance_Matrix.txt',quote=F,sep=' ',row.names=F,col.names=F)

  ##### Females
  #the order we want, from the bed file 
  ids = read.table('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/manyhost_hunt/males/beds/chr_W.fam')

  #extract columns in same order as geographic distance
  common_names = intersect(rownames(autos),ids$V1)
  auto_aligned = autos[common_names, common_names] #subset 'auto' matrix to match the 'geo_mat' matrix rows and columns
  auto_mat = as.matrix(auto_aligned)

  #visualize autosomal distance matrix in terms of ancestry K = 5
  kcols = data.frame(Kcols = gray.colors(5,start=0.05,end=0.95)[c(1,4,3,5,2)], Kcluster = paste0('K',seq(1,5,1)))
  auto_mds = as.data.frame(cmdscale(auto_mat,2))
  auto_mds$ID = rownames(auto_mds)
  auto_mds_input = left_join(auto_mds,md)
  auto_p = ggplot(auto_mds_input,aes(x=-V1,y=V2,fill=AncestryK5,shape=AncestryK5)) + 
	    geom_point(size=1.5) + theme_bw() + xlab('MDS1')+ylab('MDS2')+ggtitle('Autosomal Distance')+
	      scale_fill_manual('K',values=kcols$Kcols,breaks=kcols$Kcluster)+
	        scale_shape_manual(values=c(21,22,23,24,25))+
		  guides(fill=guide_legend(override.aes=list(shape=21)))+
		    theme(legend.position = "top", # Moves legend to top
			          legend.text = element_text(size = 6), # Adjusts text size
				          legend.title = element_text(size = 6))
  auto_p

  pdf('../../figures/GWAS_Covariance_Matrix-Females_2024APR3.pdf',height=3,width=3)
  auto_p
  dev.off()

  write.table(auto_mat,file='beds/Covariance_Matrix_N60.txt',quote=F,sep=' ',row.names=F,col.names=F)
