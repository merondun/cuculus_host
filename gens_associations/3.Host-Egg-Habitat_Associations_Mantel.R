#### Calculate pairwise geographic distance between distance groups 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/mantel')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(viridis)
library(geosphere)
library(sf)
library(ggspatial)
library(RColorBrewer)
library(cluster)
library(vegan)
library(ape)
library(ggpubr)

#calculate geographic distance matrix
mdf = read_tsv('~/merondun/cuculus_host/Metadata_Host.txt') 

#filter for unrelated canorus with egg phenotypes
md = mdf %>% filter(Analysis_PopulationGenetics == 1 & SpeciesShort == 'CC') %>% drop_na(Egg) %>% arrange(ID) %>% 
  dplyr::rename(Host=HostParentShort)
kept_order_samples = md$ID

#Calculate pairwise geographic distance 
md_sf = st_as_sf(md, coords = c("Longitude", "Latitude"))
st_crs(md_sf) = 4326 # Set the projection as ESPG 4326 (long_lat)
#Apply st_distance 
geo_m = as.data.frame(st_distance(md_sf))
geo_df = geo_m %>% mutate(across(everything(), ~ as.numeric(gsub(" \\[m\\]$", "", .))/1000))
names(geo_df) = md$ID
rownames(geo_df) = md$ID
geo_mat = as.matrix(geo_df)

#visualize the geographic distance matrix
geo_mds = as.data.frame(cmdscale(geo_mat,2))
geo_mds$ID = rownames(geo_mds)
geo_mds_input = left_join(geo_mds,md)
geo_p = ggplot(geo_mds_input,aes(x=V1,y=V2,color=CountryFull)) + 
  geom_point() + theme_bw() + xlab('MDS1')+ylab('MDS2')+ggtitle('Geographic Distance')+
  scale_color_viridis(discrete=TRUE,option='turbo')+
  theme(legend.position = "top", # Moves legend to top
        legend.text = element_text(size = 6), # Adjusts text size
        legend.title = element_text(size = 6))
geo_p

#load in autosomal distance, calculated: ~/modules/VCF2Dis-1.50/bin/VCF2Dis -InPut autos_canorus_LD.vcf.gz -OutPut autos_canorus_LD.pdist 
auto = read.table('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/merged_unrelfull/autos_canorus_LD.pdist',header=F)
#add proper names, because VCF2DIS truncates them
auto = auto %>% mutate(IDNum = gsub('_.*','',V1))
auto_id = left_join(auto,mdf %>% select(ID) %>% mutate(IDNum = gsub('_.*','',ID)))  
ord = auto_id$ID
autos = auto_id %>% select(!c(IDNum,ID,V1))
names(autos) = ord
rownames(autos) = ord

#extract columns in same order as geographic distance
common_names = intersect(rownames(autos), rownames(geo_mat))
auto_aligned = autos[common_names, common_names] #subset 'auto' matrix to match the 'geo_mat' matrix rows and columns
data.frame(mat1 = rownames(auto_aligned), mat2 = rownames(geo_mat)) %>% filter(mat1 != mat2) #sanity check for alignment 
auto_mat = as.matrix(auto_aligned)

#visualize autosomal distance matrix in terms of ancestry K = 5
kcols = data.frame(Kcols = gray.colors(5,start=0.05,end=0.95)[c(1,4,3,5,2)], Kcluster = paste0('K',seq(1,5,1)))
auto_mds = as.data.frame(cmdscale(auto_mat,2))
auto_mds$ID = rownames(auto_mds)
auto_mds_input = left_join(auto_mds,md)
auto_p = ggplot(auto_mds_input,aes(x=V1,y=V2,fill=AncestryK5,shape=AncestryK5)) + 
  geom_point(size=1.5) + theme_bw() + xlab('MDS1')+ylab('MDS2')+ggtitle('Autosomal Distance')+
  scale_fill_manual('K',values=kcols$Kcols,breaks=kcols$Kcluster)+
  scale_shape_manual(values=c(21,22,23,24,25))+
  guides(fill=guide_legend(override.aes=list(shape=21)))+
  theme(legend.position = "top", # Moves legend to top
        legend.text = element_text(size = 6), # Adjusts text size
        legend.title = element_text(size = 6))
auto_p

#calculate mtDNA distance
seq = read.dna('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/beast_dating/variant_only/ml_trees/chr_MT_All.SNP.DP3-AC1-MQ40.min4.fasta',format='fasta')
dna_dist = as.data.frame(as.matrix(dist.dna(seq,model='JC69')))

#extract columns in same order as geographic distance
common_names2 = intersect(rownames(dna_dist), rownames(geo_mat))
dna_aligned = dna_dist[common_names2, common_names2] #subset 'auto' matrix to match the 'geo_mat' matrix rows and columns
data.frame(mat1 = rownames(dna_aligned), mat2 = rownames(geo_mat)) %>% filter(mat1 != mat2) #sanity check for alignment 
dna_mat = as.matrix(dna_aligned)

#visualize the mtDNA distance matrix in terms of haplogroup 
dna_mds = as.data.frame(cmdscale(dna_mat,2))
dna_mds$ID = rownames(dna_mds)
dna_mds_input = left_join(dna_mds,md)
dna_p = ggplot(dna_mds_input,aes(x=V1,y=V2,color=Hap)) + 
  geom_point() + theme_bw() + xlab('MDS1')+ylab('MDS2')+ggtitle('mtDNA Distance')+
  scale_color_manual(values=md$HapCol,breaks=md$Hap)+
  theme(legend.position = "top", # Moves legend to top
        legend.text = element_text(size = 6), # Adjusts text size
        legend.title = element_text(size = 6))
dna_p

ggarrange(geo_p,auto_p,dna_p)

#loop through response variables
vars = c('Egg','Host','Habitat')
results = NULL; partial_mants = NULL
for (resp in vars) {
  cat('Working on variable: ',resp,'\n')
  #create dummy variables (e.g. E1 1/0/1/1/1/0/0 E2 0/1/0/0/0 wide-form, convert to factor for dissimilarity)
  d1 = md %>% select(ID,!!sym(resp)) %>% 
    mutate(value = 'P') %>%   # Create dummy variables 
    pivot_wider(names_from = !!sym(resp), values_from = value, values_fill = list(value = 'A')) %>% 
    mutate_all(as.factor) %>% 
    as.data.frame
  
  #convert response variables to df, then matrix 
  rownames(d1) = d1$ID
  d1 = d1 %>% select(-ID)
  d1_dist = daisy(d1,metric="gower")
  d1_mat = as.matrix(d1_dist)
  assign(resp,d1_mat)
  
  
  #mantel test
  geo_mantel = mantel(d1_mat, geo_mat, method="spearman", permutations=9999); geo_mantel
  auto_mantel = mantel(d1_mat, auto_mat, method="spearman", permutations=9999); auto_mantel
  dna_mantel = mantel(d1_mat, dna_mat, method="spearman", permutations=9999); dna_mantel
  
  #partial mantels, controlling for pairwise
  mat_names = c("Ancestry", "Geography", "Haplogroup")
  combinations = expand.grid(m1 = mat_names, m2 = mat_names,stringsAsFactors = FALSE)
  combinations = combinations[combinations$m1 != combinations$m2, ] %>% arrange(m1)
  
  #list matrices
  test_matrices = list(Ancestry = auto_mat, Haplogroup = dna_mat, Geography = geo_mat)
  control_matrices = list(Ancestry = auto_mat, Haplogroup = dna_mat, Geography = geo_mat)
  
  # Loop through each combination of test and control matrices
  for (row in seq(1,nrow(combinations),1)) {
    
    #grab the row index comparison
    compare = combinations[row,]
    cat('Working on variable: ',resp,', partial mantels for: ',compare[1,1],' and ', compare[1,2],'\n')
    
    #Perform the partial Mantel test
    result = mantel.partial(d1_mat, test_matrices[[compare$m1]], control_matrices[[compare$m2]],method='spearman', permutations = 9999)
    
    partial_res = data.frame(
      Response = resp,
      Variable = compare$m1,
      Control = compare$m2,
      MantelRho = result$statistic,
      Mantelp = result$signif
    )
    
    partial_mants = rbind(partial_mants,partial_res)
  }
  
  #save the standard mantel tests 
  res = data.frame(
    Response = resp,
    Variable = c('Geography','Ancestry','Haplogroup'),
    Control = 'None',
    MantelRho = c(geo_mantel$statistic, auto_mantel$statistic, dna_mantel$statistic),
    Mantelp = c(geo_mantel$signif, auto_mantel$signif, dna_mantel$signif)
  )
  resall = rbind(res,partial_mants) %>% mutate(MantelRho = abs(MantelRho)) #absolute value of correlation since our response variables are directionless  
  results = rbind(results,resall)
  partial_mants = NULL; res = NULL #ensure we don't carry over any results onto the next loop
  
}

write.table(results,file='Mantel_Results_2024MAR18.txt',quote=F,sep='\t',row.names=F)
results = read.table('Mantel_Results_2024MAR18.txt',header=TRUE)

#correct p-value, take the absolute value of rho 
results = results %>% 
  mutate(padj = p.adjust(Mantelp,method='bonferroni'), # #adjust p 
         signif = ifelse(padj < 0.05,1,0),
         signiflab = ifelse(padj < 0.05,'*','n.s.'))

#plot the partial mantels 
mantel_p = results %>%   
  filter(Control != 'None') %>% 
  mutate(padj = p.adjust(Mantelp,method='bonferroni'), # #adjust p 
         signif = ifelse(padj < 0.05,1,0),
         signiflab = ifelse(padj < 0.05,'*','n.s.')) %>% 
  ggplot(aes(x=Variable,fill=Control,y=MantelRho,label=signiflab))+
  geom_text(vjust=-1,position=position_dodge(width=0.9))+
  geom_bar(col='black',stat='identity',position=position_dodge(width=0.9))+
  scale_fill_manual(values=c(brewer.pal(3,'Set2')))+
  ylim(c(0,0.6))+
  facet_grid(Response~.,scales='free')+
  theme_bw(base_size=6)
mantel_p

#plot the simple mantels 
simple_mantel_p = results %>%   
  filter(Control == 'None') %>% 
  mutate(padj = p.adjust(Mantelp,method='bonferroni'), # #adjust p 
         signif = ifelse(padj < 0.05,1,0),
         signiflab = ifelse(padj < 0.05,'*','n.s.')) %>% 
  ggplot(aes(y=Response,fill=Variable,x=MantelRho,label=signiflab))+
  geom_text(hjust=-1,position=position_dodge(width=0.9))+
  geom_bar(col='black',stat='identity',position=position_dodge(width=0.9))+
  scale_fill_manual(values=c(brewer.pal(3,'Set2')))+
  xlim(c(0,0.6))+
  theme_bw(base_size=6)+theme(legend.position='top')
simple_mantel_p

pdf('../figures/MantelTest_2024MAR18.pdf',height=2,width=1.5)
simple_mantel_p
dev.off()

pdf('../figures/PartialMantelTest_2024MAR18.pdf',height=6,width=4)
mantel_p
dev.off()


#visualize the egg distance matrix for sanity, first ensure that there is 0 variation among eggs 
egg_df = as.data.frame(Egg)
egg_df$ID = rownames(egg_df)
egg_long = egg_df %>% pivot_longer(!ID,names_to='ID2') %>% 
  left_join(.,mdf %>% select(ID,Egg)) %>% 
  left_join(.,mdf %>% select(ID2=ID,Egg2=Egg))
egg_long %>% filter(Egg != Egg2) %>% slice_min(value) #should all be 0.182, meaning all egg clades are diverged similarly
egg_long %>% filter(Egg == Egg2) %>% slice_max(value) #should all be 0, no variation within eggs 

egg_nmds = as.data.frame(metaMDS(Egg, k=2)$points)
egg_nmds$ID = rownames(egg_nmds)
egg_nmds_input = left_join(egg_nmds,md)
egg_p = ggplot(egg_nmds_input,aes(x=MDS1,y=MDS2,col=Egg)) + 
  geom_point(alpha=0.9,size=3) + theme_bw() + xlab('MDS1')+ylab('MDS2')+ggtitle('Egg Distance')+
  scale_color_viridis(discrete=TRUE,option='turbo')+
  theme(legend.position = "top", # Moves legend to bottom
        legend.text = element_text(size = 6), # Adjusts text size
        legend.title = element_text(size = 6))
egg_p


#visualize the host distance matrix for sanity, first ensure that there is 0 variation among hosts
host_df = as.data.frame(Host)
host_df$ID = rownames(host_df)
host_long = host_df %>% pivot_longer(!ID,names_to='ID2') %>% 
  left_join(.,mdf %>% select(ID,Host = HostParentShort)) %>% 
  left_join(.,mdf %>% select(ID2=ID,Host2=HostParentShort))
host_long %>% filter(Host != Host2) %>% slice_min(value) %>% count(value) #should all be 0.0952, meaning all host clades are diverged similarly
host_long %>% filter(Host == Host2) %>% slice_max(value) %>% count(value) #should all be 0, no variation within hosts 

host_nmds = as.data.frame(metaMDS(Host, k=2)$points)
host_nmds$ID = rownames(host_nmds)
host_nmds_input = left_join(host_nmds,md)
host_p = ggplot(host_nmds_input,aes(x=MDS1,y=MDS2,col=Host)) + 
  geom_point(alpha=0.9,size=3) + theme_bw() + xlab('MDS1')+ylab('MDS2')+ggtitle('Host Distance')+
  scale_color_viridis(discrete=TRUE,option='turbo')+
  theme(legend.position = "top", # Moves legend to top
        legend.text = element_text(size = 6), # Adjusts text size
        legend.title = element_text(size = 6))
host_p


#visualize the habitat distance matrix for sanity, first ensure that there is 0 variation among habitat 
habitat_df = as.data.frame(Habitat)
habitat_df$ID = rownames(habitat_df)
habitat_long = habitat_df %>% pivot_longer(!ID,names_to='ID2') %>% 
  left_join(.,mdf %>% select(ID,Habitat = Habitat)) %>% 
  left_join(.,mdf %>% select(ID2=ID,Habitat2=Habitat))
habitat_long %>% filter(Habitat != Habitat2) %>% slice_min(value) %>% count(value) #should all be 0.0952, meaning all habitat clades are diverged similarly
habitat_long %>% filter(Habitat == Habitat2) %>% slice_max(value) %>% count(value) #should all be 0, no variation within habitats 

habitat_nmds = as.data.frame(metaMDS(Habitat, k=2)$points)
habitat_nmds$ID = rownames(habitat_nmds)
habitat_nmds_input = left_join(habitat_nmds,md)
habitat_p = ggplot(habitat_nmds_input,aes(x=MDS1,y=MDS2,col=Habitat)) + 
  geom_point(alpha=0.9,size=3) + theme_bw() + xlab('MDS1')+ylab('MDS2')+ggtitle('Habitat Distance')+
  scale_color_viridis(discrete=TRUE,option='turbo')+
  theme(legend.position = "top", # Moves legend to bottom
        legend.text = element_text(size = 6), # Adjusts text size
        legend.title = element_text(size = 6))
habitat_p

ggarrange(egg_p,host_p,habitat_p,geo_p,auto_p,dna_p,nrow=3,ncol=2)

pdf('../figures/MantelTest_Inputs_2024MAR18.pdf',height=14,width=7)
ggarrange(egg_p,host_p,habitat_p,geo_p,auto_p,dna_p,nrow=3,ncol=2)
dev.off()
