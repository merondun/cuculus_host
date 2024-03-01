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

