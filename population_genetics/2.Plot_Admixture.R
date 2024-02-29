### Plot ADMIXTURE
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/admixture')
.libPaths('~/mambaforge/envs/R/lib/R/library')
library(tidyverse)
library(viridis)
library(ggpubr)
library(RColorBrewer)
library(meRo) #devtools::install_github('merondun/meRo')

prefix = 'autos_canorus_LD' #1546013 SNPs 
qdir = '.' #directory with Q files

#read in admixture q's 
admix = melt_admixture(prefix = prefix, qdir = qdir)

#read in metadata
md = read_tsv('../Cuckoo_Full_Metadata_2023OCT3.txt')
admixmd = left_join(admix,md)

#visualize corelation of residuals between individuals 0 = good fit, if same population, positively correlated, different histories = negative corrleated
source('~/modules/evalAdmix/visFuns.R')

#this is only if you want to visualize correlation structure of a single run, recommend only for final one
Kval=5
pop <- admixmd %>% select(ID,KDist) %>% unique %>% select(KDist) # N length character vector with each individual population assignment
q <- as.matrix(read.table(paste0(prefix,'.',Kval,'.Q'))) # admixture porpotions q is optional for visualization but if used for ordering plot might look better
r <- as.matrix(read.table(paste0("eval_",Kval)))
plotAdmix(q=q, pop=pop$KDist)
plotCorRes(cor_mat = r, pop = pop$KDist, title = "Admixture evaluation as correlation of residuals", max_z=0.25, min_z=-0.25)

#loop through all admixture runs and extract the average correlation values, we want to MINIMIZE this! (closest to 0)
evaldat = NULL; for (Kval in seq(2,10,1)){
  r <- as.matrix(read.table(paste0("eval_",Kval)))
  mean_value <- mean(r,na.rm=TRUE)
  median_value <- median(r,na.rm=TRUE)
  sd_value <- sd(r,na.rm=TRUE)
  iqr_value <- IQR(r,na.rm=TRUE)
  valdat = data.frame(K = Kval,mean = mean_value,median=median_value,sd=sd_value,iqr=iqr_value)
  evaldat = rbind(valdat,evaldat)
}

#plot, for main figure show the n=3 lowest median
targs = evaldat %>% slice_min(median,n=3)
ep = evaldat %>% 
  ggplot(aes(x=K,y=median,ymin=median-iqr,ymax=median+iqr))+
  geom_rect(data=targs,aes(xmin=K-0.25,xmax=K+0.25,ymin=-Inf,ymax=Inf),fill='darkseagreen3')+
  geom_text(aes(y = 0.014, label = format(signif(median, 2), scientific = TRUE)),size=2) +  ylim(c(-0.015,0.015))+
  geom_hline(yintercept=0,lty=2)+
  geom_point()+ylab('Median +/- IQR Correlation of Residuals') +
  geom_errorbar()+
  theme_bw() + 
  scale_x_continuous(breaks = seq(min(evaldat$K), max(evaldat$K), by = 1)) +
  coord_flip()
ep
pdf('~/merondun/cuculus_host/population_genetics//ADMIXTURE_evalAdmix_MinimizeCorrelation_2024FEB28.pdf',height=5,width=4)
ep
dev.off()

#add a distance vector based on lat/long PC1
admixmd = admixmd %>% 
  mutate(Longitude = as.numeric(Longitude), Latitude = as.numeric(Latitude)) %>%
  { cbind(., pca = prcomp(.[, c("Longitude", "Latitude")], scale. = TRUE)$x[,1]) } %>%
  rename(Distance = pca) 

#if you want to add CV error directly on the label
cv = read.table('~/merondun/cuculus_host/population_genetics/autos_canorus_LD-ADMIXTURE.CVs.txt',header=FALSE)
names(cv) = c('Specified_K','d1','d2','Error')
cv = cv %>% mutate(label = paste0('K',Specified_K,' (',round(Error,2),')')) %>% select(!c(d1,d2,Error))

admixmd = admixmd %>% mutate(ID = fct_reorder(ID,desc(Distance))) 
admixmd = admixmd %>% mutate(Klab = paste0('K',Specified_K))
ad = left_join(admixmd,cv)
#save input
write.table(ad,file='~/merondun/cuculus_host/population_genetics/autos_canorus_LD-ADMIXTURE.Qmatrix.txt',quote=F,sep='\t',row.names=F)
ad = read_tsv('~/merondun/cuculus_host/population_genetics/autos_canorus_LD-ADMIXTURE.Qmatrix.txt')
k_order = ad %>% select(Klab,Specified_K) %>% arrange(Specified_K) %>% unique
ad$Hap = factor(ad$Hap,levels=c('W1','W2','W3','W4','W5','W6','W7'))
ad$Klab = factor(ad$Klab,levels=k_order$Klab)
ad = ad %>% group_by(Hap) %>% mutate(ID = reorder(factor(ID),Distance))
k_order = ad %>% ungroup %>% select(label,Specified_K) %>% arrange(Specified_K) %>% unique
ad$label = factor(ad$label,levels=k_order$label)
ad$KDist = factor(ad$KDist,levels=c('D6','D9','D12','D4','D2','D5','D1','D8','D13','D11','D14','D7','D10','D3'))

#add custom grey color scale 
kcols = data.frame(Kcols = gray.colors(5,start=0.05,end=0.95)[c(1,4,3,5,2)], Kcluster = paste0('K',seq(1,5,1)))

#Now plot
adplot = ad %>% 
  #filter(grepl('^2$|^3|^5', Specified_K)) %>%  # Specify the levels you want 
  #filter(Specified_K < 10) %>%  # Specify the levels you want 
  ggplot(aes(x = ID, y = Q, fill = factor(K))) +
  geom_col(color = NA, width=1) +
  scale_fill_manual(values = viridis(10, option = 'turbo')) +
  #scale_fill_manual(values = kcols$Kcols, breaks = kcols$Kcluster)+
  facet_grid(label ~ KDist, switch = "y", scales = "free", space = "free") +  # Switch facet labels to the left
  #facet_grid(label ~ Hap, switch = "y", scales = "free", space = "free") +  # Switch facet labels to the left
  theme_minimal() + labs(x = "", y = "Autosomal Ancestry Coefficient (Q)") +
  scale_y_continuous(expand = c(0, 0), n.breaks = 3, position = "right") +  # Move y-axis to the right
  scale_x_discrete(expand = expand_scale(add = 1)) +
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 7),
    panel.grid = element_blank(),
    strip.text.y.left = element_text(size = 8, angle = 0),
    strip.placement = "outside",  # Ensure facet labels are outside the plot area
    legend.position = 'bottom'
  ) 

adplot

pdf('~/merondun/cuculus_host/population_genetics/ADMIXTURE_K2K3K5-HAPS_2024FEB28.pdf',height=2.5,width=7.25)
adplot
dev.off()

pdf('~/merondun/cuculus_host/population_genetics/ADMIXTURE_K2-K10-HAPS_2024FEB28.pdf',height=8,width=7.25)
pdf('~/merondun/cuculus_host/population_genetics/ADMIXTURE_K2-K10-DIST_2024FEB28.pdf',height=8,width=7.25)
adplot
dev.off()

#assign individuals a K, and then plot PCA to see how those Ks fit 
#read PCA data
vec = read.table('~/merondun/cuculus_host/population_genetics/autos_canorus_LD-PCA.eigenvectors.txt',header=TRUE,comment.char = '')
vec = vec %>% dplyr::rename(ID = X.IID)
val = read.table('~/merondun/cuculus_host/population_genetics/autos_canorus_LD-PCA.eigenvalues.txt',header=FALSE,comment.char = '')
val = val %>% mutate(VE = paste0(round(V1/sum(V1),2)*100,'%'))

for (k in c(2,3,5)){ 
  #ks = ad %>% filter(Specified_K == k) %>% group_by(ID) %>% 
  #  slice_max(Q) %>% select(ID,K) %>% dplyr::rename(KCluster=K)
  
  ks = ad %>% 
    drop_na(HostParentShort) %>% 
    filter(Specified_K == k) %>%
    group_by(ID) %>%
    arrange(desc(Q)) %>%
    mutate(Rank = row_number()) %>%
    summarise(
      KCluster = case_when(
        max(Q) > 0.55 ~ first(K),
        TRUE ~ paste(K[1], K[2], sep = "_")
      )
    )
  cat('Number of samples intermediate for K',k,': ',ks %>% filter(grepl('_',KCluster)) %>% nrow,'\n')
  
  #merge with metadata, add a $distance vector based on lat/long PC1
  vc = left_join(ks,vec)
  p1 = vc %>% ggplot(aes(x=PC1,y=PC2,fill=KCluster))+
    geom_point(pch=21)+xlab(paste0('PC1 VE: ',val[1,2]))+ylab(paste0('PC2 VE: ',val[2,2]))+
    scale_fill_viridis(discrete=TRUE,option='turbo')+
    theme_bw()
  p2 = vc %>% ggplot(aes(x=PC3,y=PC4,fill=KCluster))+
    geom_point(pch=21)+xlab(paste0('PC3 VE: ',val[3,2]))+ylab(paste0('PC4 VE: ',val[4,2]))+
    scale_fill_viridis(discrete=TRUE,option='turbo')+
    theme_bw()
  ps = ggarrange(p1,p2,common.legend = TRUE,nrow=1,ncol=2)
  assign(paste0('k',k),ps)
}
pdf('~/merondun/cuculus_host/population_genetics/PCA_AssignIndividualsCluster_2024FEB28.pdf',height=7,width=6)
ggarrange(k2,k3,k5,nrow=3,ncol=1)
dev.off()

k5

#simply assign based on K5 
ks = ad %>% filter(Specified_K == 5) %>% group_by(ID) %>% slice_max(Q) %>% dplyr::rename(KCluster=K) %>% select(ID,Hap,KDist,Gens,Habitat,KCluster) %>% drop_na(Gens)
write.table(ks,file='~/merondun/cuculus_host/population_genetics/output_Admixture_AssignedAncestry_K5_Gens_2024FEB28.txt',quote=F,sep='\t',row.names=F)

#assign all samples 
ks = ad %>% filter(Specified_K == 5) %>% group_by(ID) %>% slice_max(Q) %>% dplyr::rename(KCluster=K) %>% select(ID,Hap,KDist,Gens,Habitat,KCluster)
write.table(ks,file='~/merondun/cuculus_host/population_genetics/output_Admixture_AssignedAncestry_K5_AllSamples_2024FEB28.txt',quote=F,sep='\t',row.names=F)