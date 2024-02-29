### Plot ADMIXTURE
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/admixture')
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

md = read_tsv('../Cuckoo_Full_Metadata_2023OCT3.txt')
admixmd = read_tsv('~/merondun/cuculus_host/population_genetics/autos_canorus_LD-ADMIXTURE.Qmatrix.txt')

### Tesselation
for (k in c(2,3,5)){
  cat('Tesselating K : ',k,'\n')
  #load in Q matrix 
  show_q = read.table(paste0('~/merondun/cuculus_host/population_genetics/admixture_q_files/autos_canorus_LD.',k,'.Q')) #read in the specific file 
  show_q_mat = as.matrix(show_q) #convert it to a matrix
  class(show_q_mat) = c('tess3Q','matrix','array') #make sure tess3r thinks that it's actually a tess object
  #grab coordinates 
  coords = admixmd %>% select(ID,Longitude,Latitude) %>% unique %>% select(Longitude,Latitude) #convert lat and long 
  coords_mat = as.matrix(coords) #convert coordinates to matrix
  #get map 
  map.polygon <- getMap(resolution = "low")
  pl = ggtess3Q(show_q_mat, coords_mat, map.polygon = map.polygon,col.palette = brewer.pal(show_k,'RdYlBu'))
  
  #plot 
  kp = pl +
    geom_path(data = map.polygon, aes(x = long, y = lat, group = group)) +
    xlim(min(coords$Longitude)-5,max(coords$Longitude)+5) +
    ylim(min(coords$Latitude)-5,max(coords$Latitude)+5) +
    coord_equal() +
    geom_point(data = coords, aes(x = Longitude, y = Latitude), size = 2,pch=21,fill='white') +
    xlab("Longitude") + ylab("Latitude") + ggtitle(paste0('K = ',show_k))+
    theme_classic()
  assign(paste0('k',k),kp)
}

pdf('~/merondun/cuculus_host/population_genetics/TESSELATION_K2-K5_2024FEB28.pdf',height=10,width=6)
ggarrange(k2,k3,k5,nrow=3)
dev.off()

#create  color/shape groups...
md = md %>% 
  filter(Analysis_PopulationGenetics == 1) %>% 
  mutate(HostParentShort = ifelse(Species_Latin == 'C. poliocephalus','Outgroup',HostParentShort)) %>% 
  arrange(Hap,HapCol)

#count the proportion of each haplogroup within each distance group 
mdc = md %>% count(KDist,Hap) %>% ungroup %>% group_by(KDist) %>% 
  mutate(Total = sum(n),
         Proportion = n/Total,
         Percent = paste0(round(n/Total,3)*100,'% (',n,')'))
#join that with lat/long 
mdcc = left_join(mdc,md %>% group_by(KDist) %>% summarize(Latitude=mean(Latitude),Longitude=mean(Longitude))) %>% 
  left_join(. , md %>% select(Hap,HapCol) %>% unique) %>% 
  arrange(KDist,Hap) %>% filter(KDist != 'CP')

#we ALSO need to create a scaling factor, based on how many haps are in each $Distance region 
scal = mdcc %>%
  select(KDist, Total) %>%
  unique() %>%
  ungroup() %>%
  mutate(Min = min(Total),
         Max = max(Total),
         Scaling_factor = ((Total - Min) / (Max - Min) * 10) + 2)
#make sure the scaling factor is linear
scal %>% ggplot(aes(x=Total,y=Scaling_factor))+geom_point()+theme_bw()
#add that scaling factor back to the haplogroups 
mdcc = left_join(mdcc,scal %>% select(KDist,Scaling_factor))

#add custom grey color scale 
kcols = data.frame(Kcols = gray.colors(5,start=0.05,end=0.95)[c(1,4,3,5,2)], Kcluster = paste0('K',seq(1,5,1)))

##### Plot pies across the world 
plot_pie <- function(data) {
  ggplot(data, aes(x = "", y = n, fill = Hap,label=paste0(KDist,': ',Total))) +
    geom_bar(col='black',lwd=0.5,width = 1, stat = "identity") +
    #geom_label(aes(x=Inf,y=Inf),fill='white',vjust=1.5,size=2)+
    #geom_text(size=2,vjust=-1)+
    coord_polar("y") +
    scale_fill_manual(values=md$HapCol,breaks=md$Hap)+
    theme_void() +
    theme(legend.position = "none")
}

#set up map and convert df to coordinate frame
world = map_data("world")
sites = st_as_sf(mdcc, coords = c("Longitude", "Latitude"), 
                 crs = 4326, agr = "constant") 

# Main map plot
#display.brewer.all(5,colorblindFriendly = TRUE)
p = 
  #ggplot()+ #for showing labels only 
  #geom_polygon(data = world, aes(x = long, y = lat, group = group),col='grey90',fill='white') +
  ggtess3Q(show_q_mat, coords_mat, map.polygon = map.polygon,col.palette = kcols$Kcols) + 
  geom_path(data = map.polygon, aes(x = long, y = lat, group = group)) +
  geom_sf(data = sites, aes(geometry = geometry), size = 0.1, alpha = 0.1) +
  coord_sf() +
  xlab('')+ylab('')+
  coord_sf(xlim = c(min(mdcc$Longitude)-5, max(mdcc$Longitude)+5), 
           ylim = c(min(mdcc$Latitude)-5, max(mdcc$Latitude)+5), expand = FALSE)+
  theme_classic(base_size = 8)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"))+
  theme(legend.position = 'top',legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))+
  annotation_scale(line_width=0.5,location='br')
#p

# Add pies
for (i in unique(mdcc$KDist)) {
  subset_data = mdcc %>% filter(KDist == i)
  lon = unique(subset_data$Longitude)
  lat = unique(subset_data$Latitude)
  scale_factor = unique(subset_data$Scaling_factor)
  cat('Scaling factor is : ',scale_factor,'\n')
  pie = plot_pie(subset_data)
  p <- p + annotation_custom(ggplotGrob(pie), 
                             xmin = lon - scale_factor, 
                             xmax = lon + scale_factor, 
                             ymin = lat - scale_factor, 
                             ymax = lat + scale_factor)
}

#p

#pdf('~/merondun/cuculus_host/population_genetics/TESSELATION_mtDNAHapPies_2024FEB28.pdf',height=6,width=9)
png('~/merondun/cuculus_host/population_genetics/TESSELATION_mtDNAHapPies_2024FEB28.png',units='in',res=600,height=3,width=5)
p
dev.off()

#re-run with labels to add the KDistance labels with sample size 
pdf('../figures/ADMIXTURE_TESSELATION_mtDNAHaps_2024FEB19_labs.pdf',height=3,width=5)
p
dev.off()
