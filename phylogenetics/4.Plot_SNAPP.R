#### Plot SNAPP annotated trees 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/beast_dating/variant_only/snapp/1Mchains/')
.libPaths('~/mambaforge/envs/R/lib/R/library')
library(ggtree)
library(phytools)
library(ape)
library(treeio)
library(viridis)
library(ggpubr)
library(RColorBrewer)
library(tidyverse)

#Constraints give this CI: 3.06 - 8.16 MYA
meanlog <- log(5)  # This is an approximation; see the explanation above
sdlog <- 0.25

# Calculate the 95% CI on the log scale
ci_lower <- qlnorm(0.025, meanlog=meanlog, sdlog=sdlog)
ci_upper <- qlnorm(0.975, meanlog=meanlog, sdlog=sdlog)

ci_lower
ci_upper

#Read in metadata
md = read_tsv('~/merondun/cuculus_host/Metadata_Host.txt')

files = list.files('.',paste0('.*ann'))
file = 'chr_W.trees.ann'
counter = 0
for (file in files){
  counter = counter +  1 
  iqtree = read.beast(file)
  if (counter == 1) { 
    gg = ggtree(iqtree,layout='dendrogram') %<+% md 
    } 
  else {
    #Rotate the nodes for consistency 
    gg = ggtree(iqtree,layout='dendrogram') %<+% md %>% ggtree::rotate(10) %>% ggtree::rotate(15) %>% ggtree::rotate(13)
  }
  
  #for rotating clades
  #ggtree(iqtree,layout='dendrogram') + geom_nodelab(aes(label=node),size=20)
      
  #add label for 95% CIs
  lab = gsub('.trees.*','',file)
  heights = gg$data$height_0.95_HPD
  df = as.data.frame(do.call(rbind, heights)) #convert the list to a data frame
  df$node_value = 1:nrow(df) # Add node values as a new column
  colnames(df) = c("value1", "value2", "node")
  df = df[, c("node", "value1", "value2")]
  df = df %>% mutate(lab = paste0(round(value1,1),' - ',round(value2,1))) %>% select(!c(value1,value2))
  
  leg = md %>% select(Hap,HapCol,SpeciesShort) %>% unique %>% drop_na(Hap) %>% mutate(shape = c(rep(21,6),4,21,8))
  gg$data = left_join(gg$data,df)
  ggp = gg  +
    geom_range(range='height_0.95_HPD', color='red', alpha=.6, size=1) +
    geom_tippoint(aes(fill = label,shape=label),size=1.5)+
    geom_nodelab(aes(label=lab),size=1.5,vjust=1) +
    ggtitle(lab)+
    #geom_tiplab(size=2)+
    #geom_nodelab(aes(x=branch, label=round(posterior, 2)), vjust=-.5, size=3) +
    scale_color_continuous(low="darkgreen", high="red") +
    scale_fill_manual(values=leg$HapCol,breaks=leg$Hap)+
    scale_shape_manual(values=leg$shape,breaks=leg$Hap)+
    theme(legend.position=c(.1, .8))+
    geom_treescale(x = 1)+
    guides(fill=guide_legend(override.aes=list(shape=21)))+
    theme(legend.position='right')
  ggp
  assign(paste0('p',counter),ggp)
} 

ggarrange(p1,p2,common.legend = TRUE)

pdf('../../../../figures/SNAPP_Divergence_Dating_All_2024MAR16.pdf',height=3.5,width=7)
ggarrange(p1,p2,common.legend = TRUE)
dev.off()
