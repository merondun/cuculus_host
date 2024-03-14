#Extract land class 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/spatial')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(terra)
library(sf)
library(tidyverse)
library(viridis)

file_list = list.files(path = ".", pattern = "\\.hdf$", full.names = TRUE)

# Initialize an empty list to store SpatRasters
rasters = list()

for (file in file_list) {
  raster_layer = rast(file, lyrs = "LC_Type1")
  
  # Append to the list
  rasters[[length(rasters) + 1]] <- raster_layer
}

# Merge all rasters into one
merged_raster = do.call(merge, rasters)
#plot(merged_raster)

writeRaster(merged_raster, 'LandClass_2010_MODIS.tif',overwrite=TRUE)
merged_raster = rast('LandClass_2010_MODIS.tif')

#Extract points
md = read_tsv('~/merondun/cuculus_host/Metadata_Host.txt')
md = md %>% drop_na(Egg)

# Convert to an sf object
points_sf <- st_as_sf(md, coords = c("Longitude", "Latitude"), crs = 4326)
points_sf_transformed <- st_transform(points_sf, crs = crs(merged_raster))

# Extract the land class values at the given coordinates
land_class_values <- terra::extract(merged_raster, points_sf_transformed)
lc_legend = read_tsv('LC.legend',col_names = F)
names(lc_legend) = c('LC_Type1','Habitat')
lc = left_join(land_class_values,lc_legend)
lcp = lc %>% group_by(Habitat) %>% summarize(total = n()) %>% 
  ggplot(aes(x=Habitat,fill=Habitat,y=total,label=total))+
  geom_bar(stat='identity')+geom_text(vjust=-1)+theme_bw()
pdf('../figures/LandClass_ValuesCounts_SINUproj_2024MAR05.pdf',height=4,width=6)
lcp
dev.off()
lc2 = cbind(md %>% select(-Habitat),lc %>% select(-ID))

#save it 
write.table(lc2,file='../admixture/Admixture_KClusters_K5_LandClass_Input_2024MAR05.txt',quote=F,sep='\t',row.names=F)

#plot to confirm
# Calculate the extent of the points
xmin <- min(st_coordinates(points_sf_transformed)[, "X"]) - 20
xmax <- max(st_coordinates(points_sf_transformed)[, "X"]) + 20
ymin <- min(st_coordinates(points_sf_transformed)[, "Y"]) - 20
ymax <- max(st_coordinates(points_sf_transformed)[, "Y"]) + 20

# Create a SpatExtent object
new_extent <- ext(xmin, xmax, ymin, ymax)

# Clip the raster to this new extent
clipped_raster <- crop(merged_raster, new_extent)

#save the points
st_write(points_sf_transformed, "Cuckoo_Locations_sinu.shp",delete_layer = TRUE)

# Plot to verify
pdf('../figures/LandClass_Values_SINUproj_2024MAR05.pdf',height=5,width=9)
plot(clipped_raster)
plot(points_sf_transformed, add = TRUE, col = "red", pch = 20)
dev.off()

