######################################################################################
#### This code converts geomorphology polygon layer to raster  (N.Hill)            ###
#### Generates distance to canyon metric (C. Gros)                                 ###
#### Modified Nov 2021                                                             ###
######################################################################################

### An updated geomorphology layer was  sourced from Alix Post (GA) and converted to a raster.

# libraries and paths 
library(raster)       ## package for raster manipulation
library(sp)
library(dplyr)
library(rgdal)        ## package for geospatial analysis
library(sf)
library(rasterDT)     ##Faster version of rasterize
library(stars)
library(rgeos)
library(stringr)

#user = "Nicole"
#user = "charley"
user="Jan"

SHOW_PLOT = FALSE

if (user == "Nicole") {
  #VM directory
  VM_path2<-"/perm_storage/shared_space/BioMAS/environmental_data/"
  
  #Nic local directory
  sci.dir <-      "C:/Users/hillna/OneDrive - University of Tasmania/UTAS_work/Projects/Benthic Diversity ARC/"
  env.derived <-  paste0(sci.dir,"data_environmental/derived/")
  tools.dir <-    paste0(sci.dir,"Analysis/Useful_Functions_Tools/")
  ARC_Data.dir <- paste0(sci.dir,"Analysis/ARC_Data/")
  path_geomorphology_gdb <- "C:\\Users\\hillna\\OneDrive - University of Tasmania\\UTAS_work\\Projects\\Benthic Diversity ARC\\data_environmental\\raw\\Geomorphology.gdb"
} else if (user == "charley") {
  sci.dir <- "C:/Users/cgros/code/IMAS/"
  ARC_Data.dir <- paste0(sci.dir,"ARC_Data/")
  env.folder <- "C:/Users/cgros/data/SO_env_layers/"
  env.raw <- paste0(env.folder, "raw/")
  env.derived <- paste0(env.folder, "derived/")
  tools.dir <- paste0(sci.dir,"Useful_Functions_Tools/")
  path_geomorphology_gdb <- paste0(env.raw,"Geomorphology.gdb")
} else if (user == "Jan") {
  sci.dir <-      "C:/Users/jjansen/Desktop/science/"
  env.derived <-  paste0(sci.dir,"data_environmental/derived/")
  tools.dir <-    paste0(sci.dir,"SouthernOceanBiodiversityMapping/Useful_Functions_Tools/")
  ARC_Data.dir <- paste0(sci.dir,"SouthernOceanBiodiversityMapping/ARC_Data/")
  env.raw <- "E:/science/data_environmental/raw/"
  path_geomorphology_gdb <- paste0(env.raw,"Geomorphology.gdb")
} 

path_r_depth <- paste0(env.derived, "Circumpolar_EnvData_500m_shelf_bathy_ibcso2_depth.tif")
path_r_depth2k <- paste0(env.derived, "Circumpolar_EnvData_2km_shelf_bathy_ibcso2_depth.tif")

### 1) set up details for conversion ----
# polar stereographic projection:
stereo <- "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

#naming conventions
string.chr <- "Circumpolar_EnvData_"
string.res <- "500m_"

#template raster
#bathy_shelf<-raster(paste0(VM_path2, "Circumpolar_EnvData_500m_shelf_bathy_gebco_depth"))
bathy_shelf<-raster(path_r_depth)
bathy_shelf2k<-raster(path_r_depth2k)

#coastline
load(paste0(env.derived,"Circumpolar_Coastline.Rdata"))

### 2) read in (as sf object) and rasterise geomporhology ----
#ogrListLayers("/perm_storage/shared_space/BioMAS/environmental_data/Geomorphology.gdb")
ogrListLayers(path_geomorphology_gdb)

#geomorph<-st_read( "/perm_storage/shared_space/BioMAS/environmental_data/Geomorphology.gdb", layer="AntarcticGeomorphology")
geomorph<-st_read(path_geomorphology_gdb, layer="AntarcticGeomorphology")

geomorph_pro<-st_transform(geomorph, CRS(stereo)) #both stereographic, but slightly different projection
geomorph_rast<-fasterizeDT(geomorph_pro, bathy_shelf, field="Feature")
geomorph_rast2k<-fasterizeDT(geomorph_pro, bathy_shelf2k, field="Feature")

if (SHOW_PLOT) {
  plot(geomorph_rast2k)
  plot(coast.proj, add=TRUE) #looks right
}

# names(geomorph_rast) <- "geomorphology"
# names(geomorph_rast2k) <- "geomorphology"
#writeRaster(geomorph_rast, filename = paste0(VM_path2, string.chr, "geomorphology"))
writeRaster(geomorph_rast, filename = paste0(env.derived, string.chr, "geomorphology.tif"), overwrite=TRUE)
writeRaster(geomorph_rast2k, filename = paste0(env.derived, string.chr, "2km_geomorphology.tif"), overwrite=TRUE)


## 3) Distance to canyonheads (for canyons above 3000 m)

# Load coastline
if (SHOW_PLOT) {
  library(RColorBrewer)
  ## Path coastline
  path_coastline <- "C:/Users/cgros/code/IMAS/ARC_Data/prep_environment/Circumpolar_Coastline.Rdata"
  load(path_coastline)
}

# Load bathymetry raster
r_depth <- raster(path_r_depth)

# Load polygons of canyons
p_canyon <- as_Spatial(geomorph_pro[geomorph_pro$Feature == "Canyon", ])

# Check identical CRS
is_same_crs <- identical(sort(unlist(str_split(crs(r_depth), ' '))), 
                         sort(unlist(str_split(crs(p_canyon), ' '))))
stopifnot(is_same_crs)

# Create list to store the canyon heads
lst_heads <- c()

if (SHOW_PLOT) {
  par(mar=c(0,0,0,0))
}
# Iterate through all canyons
for (idx_canyon in 1:nrow(p_canyon)) {
  print(idx_canyon)
  # Select current canyon
  p_canyon_cur <- p_canyon[idx_canyon, ]
  if (SHOW_PLOT) {
    plot(p_canyon_cur, col="blue")
  }
  
  # Get depth values inside canyon
  v_depth_canyon_cur <- raster::extract(r_depth, p_canyon_cur)[[1]]
  v_depth_canyon_cur <- v_depth_canyon_cur[!(is.na(v_depth_canyon_cur))]
  
  # If there is at least a section of the canyon above -3000m
  if (length(v_depth_canyon_cur)) {
    # Define a cut off depth value as the higher decile of the canyon above -3000m
    thr_depth = quantile(v_depth_canyon_cur, probs=c(0.90), na.rm = TRUE)[[1]]
    
    # Rasterize the current canyon
    r_canyon_cur <- rasterize(p_canyon_cur, r_depth)
    
    # Assign NA to raster cells under -3000m (NA)
    values(r_canyon_cur)[is.na(values(r_depth))] <- NA
    
    # Assign NA to raster cells under the cut off value
    values(r_canyon_cur)[values(r_depth) < thr_depth] <- NA
    
    # Convert the non-NA cells to polygons
    # Each polygon represents a canyon head
    p_heads <- sf::as_Spatial(sf::st_as_sf(stars::st_as_stars(r_canyon_cur),
                                           as_points = FALSE, merge = TRUE))
    
    # Get the centroid of each head
    p_centroids_cur <- gCentroid(p_heads, byid=TRUE)
    if (SHOW_PLOT) {
      plot(p_centroids_cur, col="red", add=TRUE)
    }
    
    # Save the result
    lst_heads <- c(lst_heads, p_centroids_cur)
  }
}

# Concatenate all the canyon heads centroids in a SpatialPointsDataFrame
p_centroids <- do.call(rbind, lapply(1:length(lst_heads), function(x) {
  SpatialPointsDataFrame(coords = lst_heads[[x]]@coords, data=data.frame(val=replicate(length(lst_heads[[x]]), 1)))}))
if (SHOW_PLOT) {
  par(mar=c(0,0,0,0))
  plot(p_canyon)
  plot(p_centroids, col="red", add=TRUE)
  plot(coast.proj, col="green", add=TRUE)
}

# Rasterize these points to remove the potential duplicates within the same raster cell
r_centroids <- rasterize(p_centroids, r_depth)
p_centroids_cleaned <- rasterToPoints(r_centroids$val, spatial = TRUE)

# Compute the distance from the canyon heads to all cells of the reference raster
distance_to_canyon <- distanceFromPoints(r_depth, p_centroids_cleaned)
distance_to_canyon[is.na(r_depth)] <- NA

# Save result
names(distance_to_canyon) <- "distance2canyons"
writeRaster(distance_to_canyon, filename = paste0(env.derived, string.chr, "500m_shelf_distance2canyons.tif"), overwrite=TRUE)
#writeRaster(distance_to_canyon, filename = paste0(env.derived, string.chr, "2km_shelf_distance2canyons.tif"), overwrite=TRUE)

if (SHOW_PLOT) {
  par(mar=c(0,0,0,0))
  plot(distance_to_canyon)
  plot(coast.proj, col="green", add=TRUE)
}

