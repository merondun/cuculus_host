#### Calculate pairwise geographic distance between distance groups 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/fst/2024feb_correlations/')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(viridis)
library(geosphere)
library(sf)
library(ggspatial)
library(RColorBrewer)

#Read in metadata
md = read_tsv('../../Cuckoo_Full_Metadata_2023OCT3.txt')

#groups for fst: distance ~ mtDNA
analyze_these = md %>% 
  filter(Analysis_PopulationGenetics == 1 & SpeciesShort == 'CC' & Sex == 'F') %>% select(KDist) %>%
  group_by(KDist) %>% mutate(DistanceHaps = n()) %>% unique %>% filter(DistanceHaps >= 3) %>% pull(KDist)

#calculate geographic distance between the groups
dists = md %>% filter(Analysis_PopulationGenetics == 1 & SpeciesShort == 'CC' & KDist %in% analyze_these) %>% group_by(KDist) %>% summarize(meanLat = mean(Latitude),meanLong = mean(Longitude))

#plot 
world = map_data("world")
sites = st_as_sf(dists, coords = c("meanLong", "meanLat"), 
                 crs = 4326, agr = "constant") 
fst_compars = ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),col='grey90',fill='white') +
  geom_sf(data = sites, 
          aes(fill=KDist,shape=KDist),
          size=4,alpha=0.9,show.legend = T,stroke=0.5) +
  scale_shape_manual(values=md$KDShape,breaks=md$KDist)+
  scale_fill_manual(values=md$KDCol,breaks=md$KDist)+
  xlab('')+ylab('')+
  coord_sf(xlim = c(min(dists$meanLong)-5, max(dists$meanLong)+5), 
           ylim = c(min(dists$meanLat)-5, max(dists$meanLat)+5), expand = FALSE)+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"))+
  theme(legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))+
  annotation_scale(line_width=0.5)

pdf('~/merondun/cuculus_host/correlations_geography_genetics/Geographic_Clusters_FST_Comparisons_2024FEB28.pdf',height=6,width=9)
fst_compars
dev.off()

# initialize an empty data frame to store pairwise distances
pairwise_distances <- data.frame()
sub_data = dists

# calculate pairwise distances for the current 'W' KDist
for(i in 1:(nrow(sub_data) - 1)) {
  for(j in (i + 1):nrow(sub_data)) {
    
    point1 <- c(sub_data$meanLong[i], sub_data$meanLat[i])
    point2 <- c(sub_data$meanLong[j], sub_data$meanLat[j])
    
    distance_km <- distHaversine(point1, point2) / 1000  # convert to km
    
    # append the result to the pairwise_distances data frame
    pairwise_distances <- rbind(pairwise_distances, 
                                data.frame(P1 = sub_data$KDist[i], 
                                           P2 = sub_data$KDist[j],
                                           Distance_km = distance_km))
  }
}

pairwise_distances = pairwise_distances %>% mutate(Group = paste0(P1,'__',P2))

# show the calculated pairwise distances
write_tsv(pairwise_distances,file='~/merondun/cuculus_host/correlations_geography_genetics/Pairwise_GeographicDistance_Km_2024FEB28.txt')

#write out populations
for (pop in unique(analyze_these)) {
  su = md %>% filter(Retained_Full_Unrelated & KDist == pop)
  write.table(su$ID,file=paste0('populations/',pop,'.pop'),quote=F,sep='\t',row.names=F,col.names=F)
}
