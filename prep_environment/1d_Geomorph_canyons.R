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

#VM directory
VM_path2<-"/perm_storage/shared_space/BioMAS/environmental_data/"

#Nic local directory
sci.dir <-      "C:/Users/hillna/OneDrive - University of Tasmania/UTAS_work/Projects/Benthic Diversity ARC/"
env.derived <-  paste0(sci.dir,"data_environmental/derived/")
tools.dir <-    paste0(sci.dir,"Analysis/Useful_Functions_Tools/")
ARC_Data.dir <- paste0(sci.dir,"Analysis/ARC_Data/")

### 1) set up details for conversion ----
# polar stereographic projection:
stereo <- "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

#naming conventions
string.chr <- "Circumpolar_EnvData_"
string.res <- "500m_"

#template raster
#bathy_shelf<-raster(paste0(VM_path2, "Circumpolar_EnvData_500m_shelf_bathy_gebco_depth"))
bathy_shelf<-raster(paste0(env.derived, "Circumpolar_EnvData_500m_shelf_bathy_gebco_depth"))

#coastline
load(paste0(env.derived,"Circumpolar_Coastline.Rdata"))

### 2) read in (as sf object) and rasterise geomporhology ----
#ogrListLayers("/perm_storage/shared_space/BioMAS/environmental_data/Geomorphology.gdb")
ogrListLayers("C:\\Users\\hillna\\OneDrive - University of Tasmania\\UTAS_work\\Projects\\Benthic Diversity ARC\\data_environmental\\raw\\Geomorphology.gdb")

#geomorph<-st_read( "/perm_storage/shared_space/BioMAS/environmental_data/Geomorphology.gdb", layer="AntarcticGeomorphology")
geomorph<-st_read( "C:\\Users\\hillna\\OneDrive - University of Tasmania\\UTAS_work\\Projects\\Benthic Diversity ARC\\data_environmental\\raw\\Geomorphology.gdb", layer="AntarcticGeomorphology")

geomorph_pro<-st_transform(geomorph, CRS(stereo)) #both stereographic, but slightly different projection
geomorph_rast<-fasterizeDT(geomorph_pro, bathy_shelf, field="Feature")

plot(geomorph_rast)
plot(coast.proj, add=TRUE) #looks right

#writeRaster(geomorph_rast, filename = paste0(VM_path2, string.chr, "geomorphology"))
writeRaster(geomorph_rast, filename = paste0(env.derived, string.chr, "geomorphology"), overwrite=TRUE)

