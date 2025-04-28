# Analyses: Population Genetic Differentiation

## Assign Geographic Clades 

Use k-means clustering to assign samples into discrete geographic 'populations' for analysis. 

```R
#### Assign geographic distance groups with k-means
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/relatedness')
.libPaths('~/mambaforge/envs/r/lib/R/library')
#Igraph approach
library(tidyverse)
library(RColorBrewer)
library(geosphere)
library(igraph)
library(spThin)
library(sf)
library(ggspatial)
library(factoextra)
library(ggpubr)

# Read metadata and filter for necessary columns
md = read_tsv('../Cuckoo_Full_Metadata_2023OCT3.txt') 
all = md %>% filter(Retained_Full == 1 & SpeciesShort == 'CC')

ds = all %>% mutate(INDEX = row_number())
thinned <- thin(ds,lat.col = 'Latitude',long.col = 'Longitude',spec.col = 'SpeciesShort',
                thin.par=25,reps=1,out.dir='.',out.base=sp,locs.thinned.list.return=T,write.files = F)
takeindex <- thinned %>% data.frame() %>% row.names() %>% as.vector()
thindat <- ds[grepl(paste0('^',takeindex,'$',collapse='|'),ds$INDEX),]


##### K-Means Clustering #####
sites <- st_as_sf(thindat, coords = c("Longitude", "Latitude"), 
                  crs = 4326, agr = "constant")

dist = st_distance(sites$geometry,  by_element = F)
d1 <- as.data.frame(dist)
d2 <- as.data.frame(lapply(d1, function(y) gsub(" [m]", "", y)))
names(d2) <- thindat$ID
row.names(d2) <- thindat$ID
d3 <- as.dist(as.matrix(d2))
cs1 <- fviz_nbclust(as.matrix(d3),kmeans,method = 'gap_stat',k.max = 25)
opt = cs1$data %>% ggplot(aes(group='Hi',x=clusters,y=gap,ymin=ymin,ymax=ymax))+
  geom_errorbar(col='cadetblue3')+ylab('Gap Statistic')+xlab('Number of Clusters')+
  geom_vline(xintercept=14,lty=2)+
  geom_line(col='cadetblue3')+
  geom_point(pch=16,size=3,col='cadetblue3')+
  theme_bw()
opt
pdf('../figures/OptimalKDistance_Clusters_2023OCT31.pdf',height=5,width=6)
opt
dev.off()

##### Clade Designations #####
#create new clade designation based on species and geographic distances 
world <- map_data("world")
set.seed(9999)

final = kmeans(d2, 14, nstart = 25)
d4 = thindat %>% mutate(Cluster = final$cluster)
dsite = st_as_sf(d4, coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant") 
ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),col='grey50',fill='grey95') +
  geom_sf(data = dsite, 
          aes(fill=as.factor(Cluster),shape=as.factor(Cluster)),
          size=5,show.legend = T) +
  scale_shape_manual(values=rep(c(21,24,25,22),8))+
  scale_fill_manual(values=c(brewer.pal(12,'Paired'),rev(brewer.pal(12,'Paired'))))+
  xlab('Longitude')+ylab('Latitude')+
  coord_sf(xlim = c(min(d4$Longitude)-5, max(d4$Longitude)+5), 
           ylim = c(min(d4$Latitude)-5, max(d4$Latitude)+5), expand = FALSE)+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"),
        legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))+
  annotation_scale(line_width=0.5)

#add names / colors
avgs = d4 %>% group_by(Cluster) %>% summarize(long = mean(Longitude)) %>% arrange(long) %>% 
  mutate(KDist = paste0('D',Cluster), KDCol = rep_len(c(viridis(12, option = 'turbo'),rev(viridis(4, option = 'turbo'))), 14), KDShape = rep_len(c(21, 24, 25, 22), 14)) %>% select(-long)

d5 = left_join(d4,avgs)
dsites = st_as_sf(d5, coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant") 
ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),col='grey50',fill='grey95') +
  geom_sf(data = dsites, 
          aes(fill=KDist,shape=KDist),
          size=5,show.legend = T) +
  xlab('Longitude')+ylab('Latitude')+
  scale_fill_manual(values=d5$KDCol,breaks=d5$KDist)+
  scale_shape_manual(values=d5$KDShape,breaks=d5$KDist)+
  coord_sf(xlim = c(min(d4$Longitude)-5, max(d4$Longitude)+5), 
           ylim = c(min(d4$Latitude)-5, max(d4$Latitude)+5), expand = FALSE)+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"),
        legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))+
  annotation_scale(line_width=0.5)

write.table(d5,'../figures/KDistance_Thinned_2023OCT31.txt',quote=F,sep='\t',row.names=F)

##### Re-assign thinned samples #####
#now, since we thinned the dataset, we need to go back to the full dataset and reassign clades based on the closest point
Bfull <- d5 %>% select(Longitude,Latitude,KDist,KDCol,KDShape,ID) %>% data.frame() 
Bfull <- Bfull %>% mutate(INDEX_B = row_number())
B <- Bfull %>% select(Longitude,Latitude)
ad <- all
Afull <- ad %>% mutate(INDEX_A = row_number())
A <- Afull %>% select(Longitude,Latitude)
for(i in 1:nrow(A)){
  #calucate distance against all of B
  distances<-geosphere::distGeo(A[i,], B)/1000
  #rank the calculated distances
  ranking<-rank(distances, ties.method = "first")
  
  #find the 3 shortest and store the indexes of B back in A
  A$shortest[i]<-which(ranking ==1) #Same as which.min()
  A$shorter[i]<-which(ranking==2)
  A$short[i]<-which(ranking ==3)
  
  #store the distances back in A
  A$shortestD[i]<-distances[A$shortest[i]] #Same as min()
  A$shorterD[i]<-distances[A$shorter[i]]
  A$shortD[i]<-distances[A$short[i]]
}
A <- A %>% mutate(INDEX_A = row_number())
A <- A %>% select(INDEX_A,shortest,Longitude,Latitude)
Aall <- merge(Afull,A,by=Reduce(intersect, list(names(Afull),names(A))))
Aall <- Aall %>% rename(INDEX_B = shortest)
Btake <- Bfull %>% select(INDEX_B,KDist,KDCol,KDShape)
as <- merge(Aall,Btake,by='INDEX_B')

#plot to confirm
as1 = as %>% mutate(LatJit = jitter(Latitude,amount = 1),
                    LonJit = jitter(Longitude,amount = 1)) 
ksites = st_as_sf(as1, coords = c("LonJit", "LatJit"), crs = 4326, agr = "constant") 
allsamp = ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),col='grey50',fill='grey95') +
  geom_sf(data = ksites, 
          aes(fill=KDist,shape=KDist),
          size=5,show.legend = T) +
  xlab('Longitude')+ylab('Latitude')+
  scale_fill_manual(values=d5$KDCol,breaks=d5$KDist)+
  scale_shape_manual(values=d5$KDShape,breaks=d5$KDist)+
  coord_sf(xlim = c(min(d4$Longitude)-5, max(d4$Longitude)+5), 
           ylim = c(min(d4$Latitude)-5, max(d4$Latitude)+5), expand = FALSE)+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"),
        legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))+
  annotation_scale(line_width=0.5)
allsamp

pdf('../figures/GeographicClusters_2023OCT31.pdf',height=8,width=12)
allsamp
dev.off()

write.table(as,'../figures/KDistance_All_2023OCT31.txt',quote=F,sep='\t',row.names=F)

```





## PCA

### C. canorus

Plot PCA from plink on LD-pruned data:

```bash
setwd('~/merondun/cuculus_host/population_genetics/')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(viridis)
library(ggpubr)
library(vegan)
library(ggord)
library(recipes)
library(ggtext)
library(sf)
library(maps)
library(ggspatial)

md = read_tsv('~/merondun/cuculus_host/Metadata_Host.txt')
autos_canorus_LD-PCA.eigenvalues.txt
#read plink
prefix = 'autos_canorus_LD'
vec = read.table('autos_canorus_LD-PCA.eigenvectors.txt',header=F)
names(vec) = c('ID',paste0('PC',seq(1,10,1)))
vec = vec %>% mutate(ID = gsub('_F_.*','_F',ID),ID = gsub('_M_.*','_M',ID))
val = read.table('autos_canorus_LD-PCA.eigenvalues.txt',header=F)
val = val %>% mutate(VE = V1/sum(V1))

#merge with metadata, add a $distance vector based on lat/long PC1
vc = left_join(vec,md)
summary(prcomp(vc[, c("Longitude", "Latitude")], scale. = TRUE)) #66% of variance 
vc = vc %>% 
  mutate(Longitude = as.numeric(Longitude), Latitude = as.numeric(Latitude)) %>%
  { cbind(., pca = prcomp(.[, c("Longitude", "Latitude")], scale. = TRUE)$x[,1]) } %>%
  rename(Distance = pca) 

#superimpose
world <- map_data("world")

# Convert your data to a simple feature object
sites <- st_as_sf(vc, coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant")

#add same custom grey color scale 
kcols = data.frame(Kcols = gray.colors(5,start=0.05,end=0.95)[c(1,4,3,5,2)], AncestryK5 = paste0('K',seq(1,5,1)))

# spatil coordinates
geo_coords <- as.matrix(vc[, c("Longitude", "Latitude")])

# PCA coordinates
vc$PC1_flip = vc$PC1*-1
pca_coords <- as.matrix(vc[, c("PC1", "PC2")])

# Procrustes transformation
proc_transform <- vegan::procrustes(geo_coords, pca_coords,  symmetric = FALSE, scale = TRUE, reflection = 'best')

# Procrustes similarity statistic
procrustes_statistic <- proc_transform$ss / sum(geo_coords^2)

# Extract rotation angle in degrees
rotation_matrix <- proc_transform$rotation
rotation_angle <- atan2(rotation_matrix[2, 1], rotation_matrix[1, 1]) * (180 / pi)

# Permutation test
set.seed(123) # For reproducibility
perm_test <- protest(geo_coords, pca_coords, permutations = 100000)
cat("Procrustes similarity statistic:", procrustes_statistic, "\n")
cat("Rotation angle (degrees counterclockwise):", rotation_angle, "\n")
cat("Permutation test p-value:", perm_test$signif, "\n")

#Procrustes similarity statistic: 0.1368981 
#Rotation angle (degrees counterclockwise): -6.212215 
#Permutation test p-value: 9.9999e-06 

#extract frame 
new_pc = data.frame(nLong = proc_transform$X[,1], nLat = proc_transform$X[,2], nPC1 = proc_transform$Yrot[,1], nPC2 = proc_transform$Yrot[,2], ID = vc$ID, AncestryK5 = vc$AncestryK5, PC1 = vc$PC1, PC2 = vc$PC2, Lat = vc$Latitude, Long = vc$Longitude)

rp = new_pc %>% ggplot(aes(x=PC1,y=-PC2,fill=AncestryK5))+geom_point(pch=21)+scale_fill_manual(values=kcols$Kcols,breaks=kcols$AncestryK5)+
  xlab(paste0('PC1: ',round(val[1,2],3)*100,'%'))+  ylab(paste0('PC2: ',round(val[2,2],3)*100,'%'))+
  theme_bw(base_size=8)+theme(legend.position='top')
rp
new_pc %>% ggplot(aes(x=nPC1,y=nPC2,fill=AncestryK5))+geom_point(pch=21)+scale_fill_manual(values=kcols$Kcols,breaks=kcols$AncestryK5)+theme_bw(base_size=8)+ggtitle('Raw PCA')

#plot 
new_pc %>% ggplot()+
  geom_segment(aes(x=nLong,xend=nPC1,y=nLat,yend=nPC2),alpha=0.1)+
  geom_point(aes(x=nPC1,y=nPC2,fill=AncestryK5),pch=21)+
  scale_fill_manual(values=kcols$Kcols,breaks=kcols$AncestryK5)+
  theme_classic()

# Convert to simple feature object
library(scales)
no_rotation = vc %>% mutate(PC1_scaled = rescale(PC1, to = c(min(Longitude), max(Longitude))),
                            PC2_scaled = rescale(PC2*-1, to = c(min(Latitude), max(Latitude))))
sitesG <- st_as_sf(no_rotation, coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant")
sitesP <- st_as_sf(no_rotation, coords = c("PC1_scaled", "PC2_scaled"), crs = 4326, agr = "constant")

# Get world map data
world <- map_data("world")

# Plotting
pp = no_rotation %>% 
  ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), color = 'grey90', fill = 'white') +
  geom_segment(aes(x=PC1_scaled,xend=Longitude,y=PC2_scaled,yend=Latitude),alpha=0.1,col='darkred')+
  geom_point(aes(x=PC1_scaled,y=PC2_scaled,fill = AncestryK5), size = 1, pch = 21, show.legend = TRUE, color = 'grey20', alpha = 0.9) +
  geom_point(aes(x=Longitude,y=Latitude,fill = AncestryK5), size = 0.5, pch = 21, show.legend = TRUE, color = 'grey20', alpha = 0.2) +
  #geom_sf(data = sitesG, aes(fill = AncestryK5), size = 2, pch = 21, show.legend = TRUE, color = 'grey20', alpha = 0.2) +
  #geom_sf(data = sitesP, aes(fill = AncestryK5), size = 2, pch = 21, show.legend = TRUE, color = 'grey20', alpha = 0.9) +
  #geom_point(data = vc, aes(x = PC1_rot, y = PC2_rot), color = 'red', size = 2) +  # Add rotated PCA points
  scale_fill_manual(values = kcols$Kcols, breaks = kcols$AncestryK5) +  # Use unique Haplotypes for colors
  coord_cartesian(xlim = c(min(no_rotation$Longitude)-5, max(no_rotation$Longitude)+5),
                  ylim = c(min(no_rotation$Latitude)-5, max(no_rotation$Latitude)+5)) +
  theme_classic(base_size=8) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.background = element_rect(fill = "aliceblue")) +
  theme(legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.key.size = unit(0.1, 'cm'),legend.position='none') +
  annotate('text',x=-Inf,y=Inf,label=paste0('Procrustes SS: ',signif(procrustes_statistic,3),' angle = ',signif(rotation_angle,3), ' p = ',signif(perm_test$signif,3)),
           vjust=2,hjust=-0.1,size=2)+
  xlab(paste0('PC1 (',signif(val[1,2],3)*100,'%) / Longitude'))+
  ylab(paste0('PC2 (',signif(val[2,2],3)*100,'%) / Latitude'))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))

pdf('PCA_Superimposed-SpatialDistribution_2024FEB29.pdf',height=3,width=6)
ggarrange(rp,pp,widths=c(0.4,0.6))
dev.off()

#also slightly jitter points by 1degree for visualization 
# Convert your data to a simple feature object
sitesjit <- st_as_sf(vc %>% mutate(Longitude = jitter(Longitude,amount=1),
                                   Latitude = jitter(Latitude,amount=1)), 
                     coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant")

#add same custom grey color scale 
kcols = data.frame(Kcols = gray.colors(5,start=0.05,end=0.95)[c(1,4,3,5,2)], AncestryK5 = paste0('K',seq(1,5,1)))

# Plotting the map with PCA points and lines
kmap = ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),col='grey90',fill='white') + #add basemap
  geom_sf(data = sitesjit , 
          aes(fill=AncestryK5),
          size=2,pch=21,show.legend = T,col='grey20',alpha=0.9) +  #add the geometry for the points, I use pch 21-25 so the points have outlines
  scale_fill_manual(values=kcols$Kcols,breaks=kcols$AncestryK5)+ #custom fill encoded from metadata
  xlab('')+ylab('')+
  coord_sf(xlim = c(min(vc$Longitude)-5, max(vc$Longitude)+5), 
           ylim = c(min(vc$Latitude)-5, max(vc$Latitude)+5), expand = FALSE)+ #expand the map boundaries so that the labels show up, probably an easier way but I just add a gray rectangle in inkscape afterwards in between the labels and the map 
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"))+ #make the ocean blue
  theme(legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))+ 
  annotation_scale(line_width=0.5,location='bl')+  #add scale bar (requires ggspatial)
  guides(fill=guide_legend(override.aes=list(shape=21))) #ensure the legend has a fill-able shape
kmap
hapmap = ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),col='grey90',fill='white') + #add basemap
  geom_sf(data = sitesjit , 
          aes(fill=Hap),
          size=2,pch=21,show.legend = T,col='grey20',alpha=0.9) +  #add the geometry for the points, I use pch 21-25 so the points have outlines
  scale_fill_manual(values=md$HapCol,breaks=md$Hap)+ #custom fill encoded from metadata
  xlab('')+ylab('')+
  coord_sf(xlim = c(min(vc$Longitude)-5, max(vc$Longitude)+5), 
           ylim = c(min(vc$Latitude)-5, max(vc$Latitude)+5), expand = FALSE)+ #expand the map boundaries so that the labels show up, probably an easier way but I just add a gray rectangle in inkscape afterwards in between the labels and the map 
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"))+ #make the ocean blue
  theme(legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))+ 
  annotation_scale(line_width=0.5,location='bl')+  #add scale bar (requires ggspatial)
  guides(fill=guide_legend(override.aes=list(shape=21))) #ensure the legend has a fill-able shape
hapmap

pdf('SpatialDistribution_Haplogroups-AncestryK5_2024FEB29.pdf',height=7,width=7)
ggarrange(kmap,hapmap,nrow=2)
dev.off()


```



### C. optatus

```bash
# Plot PCA for C. optatus 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/manyhost_hunt/optatus')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(viridis)
library(ggpubr)
library(vegan)
library(ggord)
library(recipes)
library(ggtext)
library(sf)
library(maps)
library(ggspatial)

md = read_tsv('metadata_optatus.txt')
#read plink
prefix = 'autos_optatus_LD'
vec = read.table('autosomal_files/autos_optatus_N50_LD.eigenvec',header=F)
names(vec) = c('ID',paste0('PC',seq(1,10,1)))
val = read.table('autosomal_files/autos_optatus_N50_LD.eigenval',header=F)
val = val %>% mutate(VE = V1/sum(V1))

#merge with metadata, add a $distance vector based on lat/long PC1
vc = left_join(vec,md)
summary(prcomp(vc[, c("Longitude", "Latitude")], scale. = TRUE)) #87% of variance 
vc = vc %>% 
  mutate(Longitude = as.numeric(Longitude), Latitude = as.numeric(Latitude)) %>%
  { cbind(., pca = prcomp(.[, c("Longitude", "Latitude")], scale. = TRUE)$x[,1]) } %>%
  rename(Distance = pca) 

#superimpose
world <- map_data("world")

# Convert your data to a simple feature object
sites <- st_as_sf(vc, coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant")

#add same custom grey color scale 
kcols = data.frame(Kcols = gray.colors(5,start=0.05,end=0.95)[c(1,4,3,5,2)], AncestryA5 = paste0('K',seq(1,5,1)))

# spatil coordinates
geo_coords <- as.matrix(vc[, c("Longitude", "Latitude")])

# PCA coordinates
vc$PC1_flip = vc$PC1*-1
pca_coords <- as.matrix(vc[, c("PC1", "PC2")])

# Procrustes transformation
proc_transform <- vegan::procrustes(geo_coords, pca_coords,  symmetric = FALSE, scale = TRUE, reflection = 'best')

# Procrustes similarity statistic
procrustes_statistic <- proc_transform$ss / sum(geo_coords^2)

# Extract rotation angle in degrees
rotation_matrix <- proc_transform$rotation
rotation_angle <- atan2(rotation_matrix[2, 1], rotation_matrix[1, 1]) * (180 / pi)

# Permutation test
set.seed(123) # For reproducibility
perm_test <- protest(geo_coords, pca_coords, permutations = 100000)
cat("Procrustes similarity statistic:", procrustes_statistic, "\n")
cat("Rotation angle (degrees counterclockwise):", rotation_angle, "\n")
cat("Permutation test p-value:", perm_test$signif, "\n")

#Procrustes similarity statistic: 0.1368981 
#Rotation angle (degrees counterclockwise): -6.212215 
#Permutation test p-value: 9.9999e-06 

#extract frame 
new_pc = data.frame(nLong = proc_transform$X[,1], nLat = proc_transform$X[,2], nPC1 = proc_transform$Yrot[,1], nPC2 = proc_transform$Yrot[,2], ID = vc$ID, AncestryA5 = vc$AncestryA5, PC1 = vc$PC1, PC2 = vc$PC2, Lat = vc$Latitude, Long = vc$Longitude)

rp = new_pc %>% ggplot(aes(x=PC1,y=-PC2,fill=AncestryA5))+geom_point(pch=21)+scale_fill_manual(values=kcols$Kcols,breaks=kcols$AncestryA5)+
  xlab(paste0('PC1: ',round(val[1,2],3)*100,'%'))+  ylab(paste0('PC2: ',round(val[2,2],3)*100,'%'))+
  theme_bw(base_size=8)+theme(legend.position='top')
rp
new_pc %>% ggplot(aes(x=nPC1,y=nPC2,fill=AncestryA5))+geom_point(pch=21)+scale_fill_manual(values=kcols$Kcols,breaks=kcols$AncestryA5)+theme_bw(base_size=8)+ggtitle('Raw PCA')

#plot 
new_pc %>% ggplot()+
  geom_segment(aes(x=nLong,xend=nPC1,y=nLat,yend=nPC2),alpha=0.1)+
  geom_point(aes(x=nPC1,y=nPC2,fill=AncestryA5),pch=21)+
  scale_fill_manual(values=kcols$Kcols,breaks=kcols$AncestryA5)+
  theme_classic()

# Convert to simple feature object
library(scales)
no_rotation = vc %>% mutate(PC1_scaled = rescale(PC1*-1, to = c(min(Longitude), max(Longitude))),
                            PC2_scaled = rescale(PC2, to = c(min(Latitude), max(Latitude))))
sitesG <- st_as_sf(no_rotation, coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant")
sitesP <- st_as_sf(no_rotation, coords = c("PC1_scaled", "PC2_scaled"), crs = 4326, agr = "constant")

# Get world map data
world <- map_data("world")

# Plotting
pp = no_rotation %>% 
  ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), color = 'grey90', fill = 'white') +
  geom_segment(aes(x=PC1_scaled,xend=Longitude,y=PC2_scaled,yend=Latitude),alpha=0.1,col='darkred')+
  geom_point(aes(x=PC1_scaled,y=PC2_scaled,fill = AncestryA5), size = 1, pch = 21, show.legend = TRUE, color = 'grey20', alpha = 0.9) +
  geom_point(aes(x=Longitude,y=Latitude,fill = AncestryA5), size = 0.5, pch = 21, show.legend = TRUE, color = 'grey20', alpha = 0.2) +
  #geom_sf(data = sitesG, aes(fill = AncestryA5), size = 2, pch = 21, show.legend = TRUE, color = 'grey20', alpha = 0.2) +
  #geom_sf(data = sitesP, aes(fill = AncestryA5), size = 2, pch = 21, show.legend = TRUE, color = 'grey20', alpha = 0.9) +
  #geom_point(data = vc, aes(x = PC1_rot, y = PC2_rot), color = 'red', size = 2) +  # Add rotated PCA points
  scale_fill_manual(values = kcols$Kcols, breaks = kcols$AncestryA5) +  # Use unique Haplotypes for colors
  coord_cartesian(xlim = c(min(no_rotation$Longitude)-5, max(no_rotation$Longitude)+5),
                  ylim = c(min(no_rotation$Latitude)-5, max(no_rotation$Latitude)+5)) +
  theme_classic(base_size=8) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.background = element_rect(fill = "aliceblue")) +
  theme(legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.key.size = unit(0.1, 'cm'),legend.position='none') +
  annotate('text',x=-Inf,y=Inf,label=paste0('Procrustes SS: ',signif(procrustes_statistic,3),' angle = ',signif(rotation_angle,3), ' p = ',signif(perm_test$signif,3)),
           vjust=2,hjust=-0.1,size=2)+
  xlab(paste0('PC1 (',signif(val[1,2],3)*100,'%) / Longitude'))+
  ylab(paste0('PC2 (',signif(val[2,2],3)*100,'%) / Latitude'))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))
pp

pdf('../../figures/20250127_PCA_Superimposed-SpatialDistribution_Coptatus.pdf',height=3,width=6)
ggarrange(rp,pp,widths=c(0.4,0.6))
dev.off()

#also slightly jitter points by 1degree for visualization 
# Convert your data to a simple feature object
sitesjit <- st_as_sf(vc %>% mutate(Longitude = jitter(Longitude,amount=1),
                                   Latitude = jitter(Latitude,amount=1)), 
                     coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant")

#add same custom grey color scale 
kcols = data.frame(Kcols = gray.colors(5,start=0.05,end=0.95)[c(1,4,3,5,2)], AncestryA5 = paste0('K',seq(1,5,1)))

# Plotting the map with PCA points and lines
kmap = ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),col='grey90',fill='white') + #add basemap
  geom_sf(data = sitesjit , 
          aes(fill=AncestryA5),
          size=2,pch=21,show.legend = T,col='grey20',alpha=0.9) +  #add the geometry for the points, I use pch 21-25 so the points have outlines
  scale_fill_manual(values=kcols$Kcols,breaks=kcols$AncestryA5)+ #custom fill encoded from metadata
  xlab('')+ylab('')+
  coord_sf(xlim = c(min(vc$Longitude)-5, max(vc$Longitude)+5), 
           ylim = c(min(vc$Latitude)-5, max(vc$Latitude)+5), expand = FALSE)+ #expand the map boundaries so that the labels show up, probably an easier way but I just add a gray rectangle in inkscape afterwards in between the labels and the map 
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"))+ #make the ocean blue
  theme(legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))+ 
  annotation_scale(line_width=0.5,location='bl')+  #add scale bar (requires ggspatial)
  guides(fill=guide_legend(override.aes=list(shape=21))) #ensure the legend has a fill-able shape
kmap
hapmap = ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),col='grey90',fill='white') + #add basemap
  geom_sf(data = sitesjit , 
          aes(fill=Haplogroup),
          size=2,pch=21,show.legend = T,col='grey20',alpha=0.9) +  #add the geometry for the points, I use pch 21-25 so the points have outlines
  scale_fill_manual(values=md$HaplogroupColor,breaks=md$Haplogroup)+ #custom fill encoded from metadata
  xlab('')+ylab('')+
  coord_sf(xlim = c(min(vc$Longitude)-5, max(vc$Longitude)+5), 
           ylim = c(min(vc$Latitude)-5, max(vc$Latitude)+5), expand = FALSE)+ #expand the map boundaries so that the labels show up, probably an easier way but I just add a gray rectangle in inkscape afterwards in between the labels and the map 
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"))+ #make the ocean blue
  theme(legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))+ 
  annotation_scale(line_width=0.5,location='bl')+  #add scale bar (requires ggspatial)
  guides(fill=guide_legend(override.aes=list(shape=21))) #ensure the legend has a fill-able shape
hapmap

pdf('../../figures/20250127_SpatialDistribution_Haplogroups-AncestryA4.pdf',height=7,width=7)
ggarrange(kmap,hapmap,nrow=2)
dev.off()


```



## ADMIXTURE

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=10
#SBATCH --time=200:00:00

K=$1

cd admixture

#Run Admixture
admixture -j7 --cv=5 ../merged_unrelfull/autos_canorus_LD.bed ${K} > autos_canorus_LD.log${K}.out

#Evaluate Runs
~/modules/evalAdmix/evalAdmix -plink ../merged_unrelfull/autos_canorus_LD -fname autos_canorus_LD.${K}.P -qname autos_canorus_LD.${K}.Q -P 10 -o eval_${K}
```

```bash
for i in {2..10}; do sbatch -J BAD_BOY_SERGIO_${i} Admixture_Eval.sh ${i}; done
```

CV Error: combine first `grep "CV" autos_canorus_LD*out | sed 's/.*log//g' | sed 's/.out:CV//g' > ~/merondun/cuculus_host/population_genetics/autos_canorus_LD-ADMIXTURE.CVs.txt` 

Also run for C. optatus:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=20
#SBATCH --time=48:00:00

# mamba activate snps

wd=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/manyhost_hunt/optatus
cd $wd

mkdir -p autosomal_files 

#merge autosomes
bcftools concat --file-list Autosomes.list --threads 20 -Oz -o autosomal_files/autos.vcf.gz
bcftools index --threads 20 autosomal_files/autos.vcf.gz

#LD prune
~/modules/plink2 --threads 20 --vcf autosomal_files/autos.vcf.gz --allow-extra-chr --set-missing-var-ids @:# \
        --rm-dup --indep-pairwise 50 5 0.1 --maf 0.05 --hwe 1e-10 --max-alleles 2 --min-alleles 2 --out autosomal_files/autos_optatus_LD
        
#extract, also a vcf and run PCA 
~/modules/plink2 --threads 20 --vcf autosomal_files/autos.vcf.gz --allow-extra-chr --set-missing-var-ids @:# \
        --extract autosomal_files/autos_optatus_LD.prune.in \
        --make-bed --recode vcf bgz --pca --out autosomal_files/autos_optatus_LD
bcftools index --threads 20 autosomal_files/autos_optatus_LD.vcf.gz
sed -i 's/chr_//g' autosomal_files/autos_optatus_LD.bim
```

### Plot, Fig 1B

Plot:

```R
### Plot ADMIXTURE across landscape (tesselation)
setwd('~/EvoBioWolf/CUCKOO_gentes/population_genetics/')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(terra)
library(tidyverse)
library(viridis)
library(stringr)
library(RColorBrewer)
library(meRo) #devtools::install_github('merondun/meRo')
library(LEA) #install with bioconductor, you don't actually need this if you impot your own q-matrix
library(tess3r) #install with github devtools
library(rworldmap) #for ggplot mapping 
library(sf) #for spatial plotting of distance vectors to confirm
library(ggspatial) #to add scale bars onto maps
library(ggpubr)
library(ggnewscale)

# Hash out which species to run, entire script will run afterwards 
admix_run = 'autos_optatus_N50_LD'
sp='CO'
admix_run = 'autos_canorus_LD'
sp='CC'

qdir = 'admixture_q_files' #directory with Q files

admix = melt_admixture(prefix = admix_run, qdir = qdir)

#read in metadata
md = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt')  
admixmd = left_join(admix,md) 

# Reorder individuals baseed on longitude
admixmd = admixmd %>% mutate(ID = fct_reorder(ID,Longitude))
if (sp == "CO") {
  kclust <- 'ACO'
  admixmd <- admixmd %>% mutate(K = gsub('K', kclust, K))
  kcols <- brewer.pal(7,'Purples')[c(1,3,5,6,7)]
  spshape=24
  filt="Cuculus optatus"
} else {
  kclust <- 'ACC'
  admixmd <- admixmd %>% mutate(K = gsub('K', kclust, K))
  kcols <- brewer.pal(6,'Greys')[c(2,3,4,5,6)]
  spshape=21
  filt="Cuculus canorus"
}

# Plot CV error 
cv = read.table(paste0(admix_run,'-ADMIXTURE.CVs.txt'),header=FALSE)
names(cv) <- c('K','d1','d2','Error')
cvs = 
  cv %>% ggplot(aes(x=K,y=Error))+
  geom_line(show.legend = F,col='black')+
  geom_point(show.legend = F,col='black',size=2)+
  scale_color_manual(values=viridis(3))+
  scale_x_continuous(breaks=function(x) pretty(x,n=10))+
  ylab('C-V Error')+
  theme_classic()
cvs
ggsave(paste0('~/symlinks/host/figures/20250318_ADMIXTURE-CV-Error_',sp,'.pdf'),cvs,height=2,width=3,dpi=600)

#if you want to add CV error directly on the label
names(cv) = c('Specified_K','d1','d2','Error')
cv = cv %>% mutate(label = paste0('K',Specified_K,' (',round(Error,2),')')) %>% select(!c(d1,d2,Error)) %>% arrange(Specified_K)
cv

# K = 2-5, by haplogroup  
adplot =
  admixmd %>% filter(Specified_K == 5 | Specified_K == 2) %>%  #specify the levels you want 
  mutate(Specified_K = paste0('K',Specified_K)) %>% 
  ggplot(aes(x = factor(ID), y = Q, fill = factor(K))) +
  geom_col(size = 0.1) +
  facet_grid(Specified_K~Haplogroup, scales = "free", space = "free") +
  theme_minimal(base_size=6) + labs(x = "",y = "") +
  scale_y_continuous(expand = c(0, 0),n.breaks = 3) +
  scale_fill_manual(values=kcols)+
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    #axis.text.x = element_blank(),
    axis.text.x=element_text(angle=90,size=5),
    axis.text.y = element_text(size=3),
    panel.grid = element_blank(),
    legend.position = 'bottom',
    plot.title = element_text(size=6)
  )
adplot

ggsave(paste0('~/symlinks/host/figures/20250330_Admixture-Fig1-K2-5_',sp,'.pdf'),adplot,height=2,width=6,dpi=600)

# K = 2-10, by geographic group 
geo_ord <- md %>% select(GeographicGroup) %>% distinct %>% mutate(geoord = as.numeric(gsub('GCC|GCO','',GeographicGroup))) %>% 
  arrange(geoord) %>% data.frame
admixmd <- admixmd %>% mutate(Specified_K = paste0('K',Specified_K))
k_ord <- admixmd %>%  select(Specified_K,K) %>% distinct %>% mutate(sp_k = as.numeric(gsub('K','',Specified_K))) %>% 
  arrange(sp_k) %>% data.frame
admixmd$GeographicGroup <- factor(admixmd$GeographicGroup,levels=geo_ord$GeographicGroup)
admixmd$Specified_K <- factor(admixmd$Specified_K,levels=unique(k_ord$Specified_K))
admixmd$K <- factor(admixmd$K,levels=unique(k_ord$K))
adplot =
  admixmd %>%  #specify the levels you want 
  ggplot(aes(x = factor(ID), y = Q, fill = factor(K))) +
  geom_col(size = 0.1) +
  facet_grid(Specified_K~GeographicGroup, scales = "free", space = "free") +
  theme_minimal(base_size=6) + labs(x = "",y = "") +
  scale_y_continuous(expand = c(0, 0),n.breaks = 3) +
  #scale_fill_manual(values=viridis(10,option='turbo'))+
  scale_fill_manual(values=brewer.pal(10,'Spectral'))+
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    #axis.text.x = element_blank(),
    axis.text.x=element_text(angle=90,size=5),
    axis.text.y = element_text(size=3),
    panel.grid = element_blank(),
    legend.position = 'bottom',
    plot.title = element_text(size=6)
  )
adplot

ggsave(paste0('~/symlinks/host/figures/20250318_Admixture-GeographicGroupK10_',sp,'.pdf'),adplot,height=8,width=6,dpi=600)


# K = 2-10, by geographic group 
egg_ord <- md %>% select(Egg,EggOrder) %>% distinct %>% arrange(EggOrder)
admixmd$Egg <- factor(admixmd$Egg,levels=egg_ord$Egg)
adplot =
  admixmd %>%  
  drop_na(Egg) %>% #specify the levels you want 
  ggplot(aes(x = factor(ID), y = Q, fill = factor(K))) +
  geom_col(size = 0.1) +
  facet_grid(Specified_K~Egg, scales = "free", space = "free") +
  theme_minimal(base_size=6) + labs(x = "",y = "") +
  scale_y_continuous(expand = c(0, 0),n.breaks = 3) +
  #scale_fill_manual(values=viridis(10,option='turbo'))+
  scale_fill_manual(values=brewer.pal(10,'Spectral'))+
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    #axis.text.x = element_blank(),
    axis.text.x=element_text(angle=90,size=5),
    axis.text.y = element_text(size=3),
    panel.grid = element_blank(),
    legend.position = 'bottom',
    plot.title = element_text(size=6)
  )
adplot

ggsave(paste0('~/symlinks/host/figures/20250318_Admixture-EggK10_',sp,'.pdf'),adplot,height=8,width=4,dpi=600)
                     
# Just K = 5, sorted by ancestry in each geo group
adplot =
  admixmd %>% filter(Specified_K == 5 | Specified_K == 2) %>%  #specify the levels you want 
  mutate(Specified_K = paste0('K',Specified_K)) %>% 
  ggplot(aes(x = factor(ID), y = Q, fill = factor(K))) +
  geom_col(size = 0.1) +
  facet_grid(AncestryA5+GeographicGroup~Specified_K, scales = "free", space = "free") +
  theme_minimal(base_size=6) + labs(x = "",y = "") +
  scale_y_continuous(expand = c(0, 0),n.breaks = 3) +
  scale_fill_manual(values=kcols)+
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    #axis.text.x = element_blank(),
    axis.text.x=element_text(angle=90,size=5),
    axis.text.y = element_text(size=3),
    panel.grid = element_blank(),
    legend.position = 'bottom',
    plot.title = element_text(size=6)
  )+
  coord_flip()
adplot

ggsave(paste0('~/symlinks/host/figures/20250408_Admixture-AncGeoK2-5_',sp,'.pdf'),adplot,height=12,width=6,dpi=600)
```

## Plot Tesselation + mtDNA Pies

```R
### Plot ADMIXTURE across landscape (tesselation)
setwd('~/EvoBioWolf/CUCKOO_gentes/population_genetics/')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(terra)
library(tidyverse)
library(viridis)
library(stringr)
library(RColorBrewer)
library(meRo) #devtools::install_github('merondun/meRo')
library(LEA) #install with bioconductor, you don't actually need this if you impot your own q-matrix
library(tess3r) #install with github devtools
library(rworldmap) #for ggplot mapping 
library(sf) #for spatial plotting of distance vectors to confirm
library(ggspatial) #to add scale bars onto maps
library(ggpubr)
library(ggnewscale)

# Hash out which species to run, entire script will run afterwards 
admix_run = 'autos_optatus_N50_LD'
sp='CO'
admix_run = 'autos_canorus_LD'
sp='CC'

qdir = 'admixture_q_files' #directory with Q files

admix = melt_admixture(prefix = admix_run, qdir = qdir)

#read in metadata
md = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt') %>% 
  select(ID,Haplogroup,Longitude,Latitude,Egg,GeographicGroup,HaplogroupColor) 
admixmd = left_join(admix,md) 

# Reorder individuals baseed on longitude
admixmd = admixmd %>% mutate(ID = fct_reorder(ID,Longitude))
if (sp == "CO") {
  kclust <- 'ACO'
  admixmd <- admixmd %>% mutate(K = gsub('K', kclust, K))
  kcols <- brewer.pal(7,'Purples')[c(1,3,5,6,7)]
  spshape=24
  filt="Cuculus optatus"
} else {
  kclust <- 'ACC'
  admixmd <- admixmd %>% mutate(K = gsub('K', kclust, K))
  kcols <- brewer.pal(6,'Greys')[c(2,3,4,5,6)]
  spshape=21
  filt="Cuculus canorus"
}


# We need to manually re-specify our K orders to reflect geography
admixmd <- admixmd %>% 
  mutate(K = case_when(
    # for optatus, swap 
    Specified_K == 5 & K == "ACO3" ~ "ACO5",
    Specified_K == 5 & K == "ACO5" ~ "ACO3",
    # for canorus, swap 
    Specified_K == 5 & K == "ACC5" ~ "ACC2",
    Specified_K == 5 & K == "ACC4" ~ "ACC5",
    Specified_K == 5 & K == "ACC2" ~ "ACC4",
    TRUE ~ K
  ))


#loop through all admixture runs and extract the average correlation values from evalAdmix, we want to MINIMIZE this! (closest to 0)
evaldat = NULL; for (Kval in seq(2,10,1)){
  r <- as.matrix(read.table(paste0("evalAdmix/",sp,"_eval_",Kval)))
  mean_value <- mean(r,na.rm=TRUE)
  median_value <- median(r,na.rm=TRUE)
  sd_value <- sd(r,na.rm=TRUE)
  iqr_value <- IQR(r,na.rm=TRUE)
  valdat = data.frame(K = Kval,mean = mean_value,median=median_value,sd=sd_value,iqr=iqr_value)
  evaldat = rbind(valdat,evaldat)
}

#plot, for main figure show the n=3 lowest median
targs = evaldat %>% slice_min(abs(median),n=3)
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

ggsave(paste0('~/symlinks/host/figures/20250318_',sp,'_evalAdmix.pdf'),ep,height=5,width=6,dpi=600)

# K = 5 
adplot =
  admixmd %>% filter(Specified_K == 5 ) %>%  #specify the levels you want 
  mutate(Specified_K = paste0('K',Specified_K)) %>% 
  ggplot(aes(x = factor(ID), y = Q, fill = factor(K))) +
  geom_col(size = 0.1) +
  facet_grid(Specified_K~., scales = "free", space = "free") +
  theme_minimal(base_size=6) + labs(x = "",y = "") +
  scale_y_continuous(expand = c(0, 0),n.breaks = 3) +
  scale_fill_manual(values=kcols)+
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    #axis.text.x = element_blank(),
    axis.text.x=element_text(angle=90,size=5),
    axis.text.y = element_text(size=3),
    panel.grid = element_blank(),
    legend.position = 'bottom',
    plot.title = element_text(size=6)
  )
adplot

ggsave(paste0('~/symlinks/host/figures/20250318_Admixture_',sp,'.pdf'),adplot,height=2.5,width=6,dpi=600)

# Assign K1 - K5 
ks = admixmd %>% filter(Specified_K == 5) %>% group_by(ID) %>% slice_max(Q) %>% dplyr::rename(KCluster=K) %>% select(ID,KCluster)

### Tesselation
# Import birdlife shapefiles, breeding == 2, presence ==1 means extant
bg <-  st_read('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/spatial/birdlife_breeding_distributions/SppDataRequest.shp')
filtered_data <- bg[bg$PRESENCE == 1 & bg$SEASONAL == 2 & bg$SCI_NAME == filt, ]

#which K value to plot 
show_k = 5

# Extract the q values and the long/lat 
tesselation_qs <- admixmd %>% filter(Specified_K == show_k) %>% select(ID,K,Q,Longitude,Latitude,Haplogroup) %>% 
  pivot_wider(names_from = K,values_from = Q) %>% select(ID, Longitude, Latitude, Haplogroup, paste0(kclust,seq(1,show_k,1)))
cluster_centroids <- left_join(tesselation_qs,ks) %>% group_by(KCluster) %>% 
  summarize(lat = mean(Latitude),lon=mean(Longitude))

#plot using ggplot
map.polygon <- getMap(resolution = "high")

# Jitter the lat/long and add haplogroup colors 
sitesp = st_as_sf(tesselation_qs %>% mutate(loj = jitter(Longitude,amount=1),laj = jitter(Latitude,amount=1)),remove = F, coords = c("loj", "laj"), crs = 4326, agr = "constant") 
max_col <- 5+show_k-1
pl = ggtess3Q(tesselation_qs[5:max_col], as.matrix(tesselation_qs[2:3]), map.polygon = map.polygon,col.palette = kcols)

# Plot the K1-5, look at color matching, etc 
check_cols = pl +
  geom_path(data = map.polygon, aes(x = long, y = lat, group = group),col='white',lwd=0.2) +
  coord_sf(xlim = c(min(admixmd$Longitude)-5, max(admixmd$Longitude)+5), 
           ylim = c(min(admixmd$Latitude)-5, max(admixmd$Latitude)+5), expand = FALSE)+    
  new_scale_fill()+
  new_scale_color()+
  geom_label(data = cluster_centroids, aes(x = lon, y = lat,label=KCluster,fill=KCluster))+
  scale_fill_manual(values=kcols)+
  xlab("Longitude") + ylab("Latitude") + 
  theme_classic(base_size=8)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"))+
  theme(legend.position = 'top',legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))
check_cols 

# Plot haplogroups 
k5p1 = pl +
  geom_path(data = map.polygon, aes(x = long, y = lat, group = group),col='white',lwd=0.2) +
  coord_sf(xlim = c(min(admixmd$Longitude)-5, max(admixmd$Longitude)+5), 
           ylim = c(min(admixmd$Latitude)-5, max(admixmd$Latitude)+5), expand = FALSE)+    
  new_scale_fill()+
  new_scale_color()+
  geom_point(data = sitesp, aes(x = loj, y = laj,fill=Haplogroup),shape=spshape, size = 2,col='black') +
  scale_fill_manual(values=md$HaplogroupColor,breaks=md$Haplogroup)+
  xlab("Longitude") + ylab("Latitude") + 
  theme_classic(base_size=8)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"))+
  theme(legend.position = 'top',legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))
k5p1

ggsave(paste0('~/symlinks/host/figures/20250318_Tesselation_',sp,'.pdf'),ggarrange(check_cols,k5p1,nrow=2),height=7,width=7,dpi=300)

##### plot all K... this will ANCESTRY Q WITHIN PIES!  #####
maximum_k = 5

for (kval in seq(2,maximum_k,1)) {
  
  cat('Working on K = ',kval,'\n')
  tesselation_qs <- admixmd %>% filter(Specified_K == kval) %>% select(ID,K,Q,Longitude,Latitude,Haplogroup) %>% 
    pivot_wider(names_from = K,values_from = Q) %>% select(ID, Longitude, Latitude, Haplogroup, paste0(kclust,seq(1,kval,1)))
  max_col <- 5 + kval - 1

  # First, grab the individuals and calculate the mean Q values within each cluster. Cluster will be geographic reigon, also calculate mean lat/long for plotting
  kept = admixmd %>% select(ID,Specified_K,K,Q,Latitude,Longitude,Group = GeographicGroup)
  group_summaries = kept %>% 
    filter(Specified_K == kval) %>%
    group_by(Group,K) %>%  #within each group and K, average lat/long/q and count number of individuals 
    summarize(Lat = mean(Latitude),
              Long = mean(Longitude),
              Q = mean(Q),
              N = n_distinct(ID)) %>% 
    ungroup() %>% #
    #Calculate scaling factors for the pies based on num samples
    mutate(MinN = min(N),
           MaxN = max(N)) %>%
    group_by(Group) %>%
    mutate(Scaling_factor = ((N - MinN) / (MaxN - MinN) * 10) + 2) %>%
    select(-MinN, -MaxN) 
  
  ##### Plot pies across the world 
  plot_pie <- function(data) {
    ggplot(data, aes(x = "", y = Q, fill = K,)) +
      geom_bar(col='white',lwd=0.5,width = 1, stat = "identity") +
      coord_polar("y") +
      scale_fill_manual(values=kcols)+
      theme_void() +
      theme(legend.position = "none")
  }
  
  #set up map and make a sf object from the summaries 
  sites = st_as_sf(group_summaries, coords = c("Long", "Lat"), crs = 4326, agr = "constant") 
  
  # Main map plot
  p = 
    ggtess3Q(tesselation_qs[5:max_col], as.matrix(tesselation_qs[2:3]), map.polygon = filtered_data,col.palette = kcols) + 
    geom_path(data = map.polygon, aes(x = long, y = lat, group = group)) +
    geom_sf(data = sites, aes(geometry = geometry), size = 0.1, alpha = 0.1, pch=26) +
    xlab('')+ylab('')+
    coord_sf(xlim = c(min(admixmd$Longitude)-5, max(admixmd$Longitude)+5), 
             ylim = c(min(admixmd$Latitude)-5, max(admixmd$Latitude)+5), expand = FALSE)+    
    theme_classic(base_size = 8)+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"))+
    theme(legend.position = 'top',legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))
  p
  
  # Add pies
  for (i in unique(group_summaries$Group)) {
    subset_data = group_summaries %>% filter(Group == i)
    lon = unique(subset_data$Long)
    lat = unique(subset_data$Lat)
    scale_factor = unique(subset_data$Scaling_factor)
    cat('Scaling factor is : ',scale_factor,' for group : ',i,'\n')
    pie = plot_pie(subset_data)
    p <- p + annotation_custom(ggplotGrob(pie), 
                               xmin = lon - scale_factor, 
                               xmax = lon + scale_factor, 
                               ymin = lat - scale_factor, 
                               ymax = lat + scale_factor)
  }
  p
  assign(paste0('t',kval),p)
  
}

t5
png(paste0('~/symlinks/host/figures/20250318_Tesselation_',sp,'_AllK.png'),res=600,units='in',height=7,width=10)
print(ggarrange(t2,t3,t4,t5,ncol=2,nrow=2))
dev.off()

##### plot all K... this will show HAPLOGROUPS within PIES!  #####
show_hap_k = 5 
tesselation_qs <- admixmd %>% filter(Specified_K == show_hap_k) %>% select(ID,K,Q,Longitude,Latitude,Haplogroup,GeographicGroup) %>% 
  pivot_wider(names_from = K,values_from = Q) %>% select(ID, Longitude, Latitude, Haplogroup, GeographicGroup, paste0(kclust,seq(1,show_hap_k,1)))
max_col <- 6 + show_hap_k - 1

# First, grab the individuals and calculate the mean Q values within each cluster. Cluster will be geographic reigon, also calculate mean lat/long for plotting
#count the proportion of each haplogroup within each distance group 
haps = tesselation_qs %>% count(GeographicGroup,Haplogroup) %>% ungroup %>% group_by(GeographicGroup) %>% 
  mutate(Total = sum(n),
         Proportion = n/Total,
         Percent = paste0(round(n/Total,3)*100,'% (',n,')'))
coords <- tesselation_qs %>% group_by(GeographicGroup) %>% 
  summarize(Latitude=mean(Latitude),
            Longitude=mean(Longitude))
group_summaries <- left_join(haps,coords) %>% left_join(., md %>% select(Haplogroup,HaplogroupColor) %>% unique) %>% 
  arrange(GeographicGroup,Haplogroup)

#we ALSO need to create a scaling factor, based on how many haps are in each $Distance region 
scal = group_summaries %>%
  select(GeographicGroup, Total) %>%
  unique() %>%
  ungroup() %>%
  mutate(Min = min(Total),
         Max = max(Total),
         Scaling_factor = ((Total - Min) / (Max - Min) * 10) + 2)
#make sure the scaling factor is linear
scal %>% ggplot(aes(x=Total,y=Scaling_factor))+geom_point()+theme_bw()
#add that scaling factor back to the haplogroups 
group_summaries = left_join(group_summaries,scal %>% select(GeographicGroup,Scaling_factor))

##### Plot pies across the world 
#set up map and make a sf object from the summaries 
sites = st_as_sf(group_summaries, coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant") 

# Main map plot
p = 
  ggtess3Q(tesselation_qs[6:max_col], as.matrix(tesselation_qs[2:3]), map.polygon = map.polygon,col.palette = kcols) + 
  geom_path(data = map.polygon, aes(x = long, y = lat, group = group)) +
  geom_sf(data = sites, aes(geometry = geometry), size = 0.1, alpha = 0.1, pch=26) +
  xlab('')+ylab('')+
  coord_sf(xlim = c(min(admixmd$Longitude)-5, max(admixmd$Longitude)+5), 
           ylim = c(min(admixmd$Latitude)-5, max(admixmd$Latitude)+5), expand = FALSE)+    
  theme_classic(base_size = 8)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"))+
  theme(legend.position = 'top',legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))
p

# Pie function 
plot_pie <- function(data) {
  ggplot(data, aes(x = "", y = n, fill = Haplogroup,label=paste0(GeographicGroup,': ',Total))) +
    geom_bar(col='black',lwd=0.5,width = 1, stat = "identity") +
    coord_polar("y") +
    scale_fill_manual(values=md$HaplogroupColor,breaks=md$Haplogroup)+
    theme_void() +
    theme(legend.position = "none")
}

# Add pies
for (i in unique(group_summaries$GeographicGroup)) {
  subset_data = group_summaries %>% filter(GeographicGroup == i)
  lon = unique(subset_data$Longitude)
  lat = unique(subset_data$Latitude)
  scale_factor = unique(subset_data$Scaling_factor)
  cat('Scaling factor is : ',scale_factor,' for group : ',i,'\n')
  pie = plot_pie(subset_data)
  p <- p + annotation_custom(ggplotGrob(pie), 
                             xmin = lon - scale_factor, 
                             xmax = lon + scale_factor, 
                             ymin = lat - scale_factor, 
                             ymax = lat + scale_factor)
  }
p

png(paste0('~/symlinks/host/figures/20250318_Tesselation-mtDNA-Pies_',sp,'_K5.png'),height=4,width=7,units='in',res=300)
p
dev.off()

final_cols <- data.frame(KCluster=paste0(kclust,seq(1,5,1)),Kcols=kcols)
final_cols %>% ggplot(aes(x = KCluster, y = 1, fill = Kcols)) +geom_tile(color = "black") + scale_fill_identity() + theme(legend.position = "none")

#simply assign based on K5 
write.table(left_join(ks,final_cols),file=paste0('~/symlinks/host/figures/20250318_',sp,'_assigned_k5.txt'),quote=F,sep='\t',row.names=F)

```

## Plot Spatial: Geography, Eggs, Ancestry

```R
#### Assign geographic distance groups with k-means
setwd('~/EvoBioWolf/CUCKOO_gentes/correlations_geography_genetics')
.libPaths('~/mambaforge/envs/r/lib/R/library')
#Igraph approach
library(tidyverse)
library(RColorBrewer)
library(geosphere)
library(igraph)
library(spThin)
library(sf)
library(ggspatial)
library(factoextra)
library(ggpubr)

# Read metadata and filter for necessary columns
md = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt') %>% drop_na(GeographicGroup) %>% arrange(GeoOrder)
world <- map_data("world")

# Plot Geographic Groups 
as1 = md %>% mutate(LatJit = jitter(Latitude,amount = 1),
                    LonJit = jitter(Longitude,amount = 1)) 
ksites = st_as_sf(as1, coords = c("LonJit", "LatJit"), crs = 4326, agr = "constant") 
allsamp = ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),col='grey50',fill='grey95') +
  geom_sf(data = ksites, 
          aes(fill=GeographicGroup,shape=GeographicGroup),
          size=3,show.legend = T) +
  xlab('Longitude')+ylab('Latitude')+
  scale_fill_manual(values=md$GeoColor,breaks=md$GeographicGroup)+
  scale_shape_manual(values=md$GeoShape,breaks=md$GeographicGroup)+
  coord_sf(xlim = c(min(as1$Longitude)-5, max(as1$Longitude)+5), 
           ylim = c(min(as1$Latitude)-5, max(as1$Latitude)+5), expand = FALSE)+
  theme_classic()+
  facet_grid(SpeciesShort~.)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"),
        legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))+
  annotation_scale(line_width=0.5)+
  theme(legend.position='top')
allsamp

ggsave('~/symlinks/host/figures/20250319_GeographicGroups.pdf',allsamp,dpi=300,height=8,width=9)

# Plot haplogroups 
as1 = md %>% mutate(LatJit = jitter(Latitude,amount = 2),
                    LonJit = jitter(Longitude,amount = 2)) 

hap_colors <- md %>% drop_na(HaplogroupColor) %>% distinct(Haplogroup, HaplogroupColor) %>% deframe()
hap_shapes <- data.frame(names = names(hap_colors),shape = rep_len(c(21,22,23,24),11)) %>% deframe()

# re-plot
ksites = st_as_sf(as1, coords = c("LonJit", "LatJit"), crs = 4326)

plothaps = ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), col = 'grey50', fill = 'grey95') +
  geom_sf(data = ksites, aes(fill = Haplogroup, shape = Haplogroup), size = 2, show.legend = TRUE) +
  scale_fill_manual(values = hap_colors) +
  scale_shape_manual(values = hap_shapes) +
  coord_sf(xlim = c(min(as1$Longitude) - 5, max(as1$Longitude) + 5), 
           ylim = c(min(as1$Latitude) - 5, max(as1$Latitude) + 5), expand = FALSE) +
  facet_grid(SpeciesShort ~ .) +
  labs(x = 'Longitude', y = 'Latitude') +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.background = element_rect(fill = "aliceblue"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.key.size = unit(0.1, 'cm'),
        legend.position = 'top') +
  annotation_scale(line_width = 0.5)
plothaps
ggsave('~/symlinks/host/figures/20250408_HaplogroupLocations.pdf',plothaps,dpi=300,height=8,width=9)



# Plot Egg
md = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt') %>% drop_na(Egg) %>% arrange(EggOrder)
egg = md %>% mutate(LatJit = jitter(Latitude,amount = 2),
                    LonJit = jitter(Longitude,amount = 2)) 
eggsites = st_as_sf(egg, coords = c("LonJit", "LatJit"), crs = 4326, agr = "constant") 
allsamp = ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),col='grey50',fill='grey95') +
  geom_sf(data = eggsites, 
          aes(fill=Egg,shape=Egg),
          size=3,show.legend = T) +
  xlab('Longitude')+ylab('Latitude')+
  scale_fill_manual(values=md$EggCol,breaks=md$Egg)+
  scale_shape_manual(values=md$EggShape,breaks=md$Egg)+
  coord_sf(xlim = c(min(egg$Longitude)-5, max(egg$Longitude)+5), 
           ylim = c(min(egg$Latitude)-5, max(egg$Latitude)+5), expand = FALSE)+
  theme_classic()+
  facet_grid(SpeciesShort~.)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"),
        legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))+
  annotation_scale(line_width=0.5)+
  theme(legend.position='top')
allsamp

ggsave('~/symlinks/host/figures/20250319_EggLocations.pdf',allsamp,dpi=300,height=8,width=9)

# Plot Ancestry
md = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt') %>% drop_na(AncestryA5) %>% arrange(AncestryA5)
anc = md %>% mutate(LatJit = jitter(Latitude,amount = 2),
                    LonJit = jitter(Longitude,amount = 2)) 
ancsites = st_as_sf(anc, coords = c("LonJit", "LatJit"), crs = 4326, agr = "constant") 
plot_anc = ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),col='grey50',fill='grey95') +
  geom_sf(data = ancsites, 
          aes(fill=AncestryA5,shape=SpeciesShort),
          size=3,show.legend = T) +
  xlab('Longitude')+ylab('Latitude')+
  scale_fill_manual(values=md$AncestryColor,breaks=md$AncestryA5)+
  scale_shape_manual(values=md$Shape,breaks=md$SpeciesShort)+
  coord_sf(xlim = c(min(anc$Longitude)-5, max(anc$Longitude)+5), 
           ylim = c(min(anc$Latitude)-5, max(anc$Latitude)+5), expand = FALSE)+
  theme_classic()+
  facet_grid(SpeciesShort~.)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"),
        legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))+
  annotation_scale(line_width=0.5)+
  theme(legend.position='top')+
  guides(fill=guide_legend(override.aes=list(shape=21)))
plot_anc

ggsave('~/symlinks/host/figures/20250319_AncestryLocations.pdf',plot_anc,dpi=300,height=8,width=9)

```

