
## changing .grd files to .tif for more efficient storage (once off, code now at the bottom)
## masking the extent of the environmental layers
## aggregating data from 500m resolution to 2km resolution

##### Setting up----
user = "Jan"
#user = "charley"
#user="nicole"
if (user == "Jan") {
  sci.dir <-      "C:/Users/jjansen/Desktop/science/"
  env.derived <-  paste0(sci.dir,"data_environmental/derived/")
  #bio.dir <-      paste0(sci.dir,"data_biological/")
  ## remote repository (DOESN'T WORK YET):
  # env.dir <- "https://data.imas.utas.edu.au/data_transfer/admin/files/EnvironmentalData/"
  ## common paths (after "sci.dir")
  tools.dir <-    paste0(sci.dir,"SouthernOceanBiodiversityMapping/Useful_Functions_Tools/")
  ARC_Data.dir <- paste0(sci.dir,"SouthernOceanBiodiversityMapping/ARC_Data/")
} 
if (user == "charley") {
  
  sci.dir <- "C:/Users/cgros/code/IMAS/"
  ARC_Data.dir <- paste0(sci.dir,"ARC_Data/")
  env.derived <-  "C:/Users/cgros/data/SO_env_layers/derived/"
  tools.dir <-    paste0(sci.dir,"Useful_Functions_Tools/")
}
if (user == "nicole") {
  
  sci.dir <-    "C:/Users/hillna/OneDrive - University of Tasmania/UTAS_work/Projects/Benthic Diversity ARC/"
  ARC_Data.dir <- paste0(sci.dir,"Analysis/ARC_Data/")
  env.derived <-  paste0(sci.dir,"data_environmental/derived/")
  tools.dir <-    paste0(sci.dir,"Analysis/Useful_Functions_Tools/")
}
##############################################################################################################
##############################################################################################################
library(terra)

## first, load bathy data and draw a mask to include/exclude areas
r2 <- rast(paste0(env.derived,"Circumpolar_EnvData_500m_shelf_bathy_gebco_depth.tif"))
#t <- draw(x="polygon")
#writeVector(t,paste0(env.derived,"Circumpolar_EnvData_mask_shelf.shp"))
t <- vect(paste0(env.derived,"Circumpolar_EnvData_mask_shelf.shp"))

## restrict depth range and mask to remove areas off the shelf
r3 <- r2
r3[r2<(-2500)] <- NA
r4 <- mask(r3,t)
r5 <- aggregate(r4,4)
r6 <- aggregate(r4,8)

## 500m res, 0-3000m depth, raw,    -> 23.6 million cells
length(which(!is.na(r2[])))
## 500m res, 0-2500m depth, raw,    -> 17.6 million cells
length(which(!is.na(r3[])))
## 500m res, 0-2500m depth, masked, -> 16.1 million cells
length(which(!is.na(r4[])))
## 2km res,  0-2500m depth, masked, -> 987k cells
length(which(!is.na(r5[])))
## 4km res,  0-2500m depth, masked, -> 241k cells
length(which(!is.na(r6[])))

### load environmental layers
env_list<-list.files(path = env.derived, pattern="tif$",  full.names=TRUE) 
#subset to  "shelf" files
env_list<-env_list[grep(".500m_shelf", env_list)]
env_list<-env_list[-grep(".500m_shelf_scaled", env_list)]
#for the single rasters layer names are missing. Extract from file name.
env_names<-gsub(".*_|\\..*","",env_list)
#stack all environmental layers and make sure they have appropriate names (currently manual and a bit messy!)
env_stack<-raster::stack(env_list)
names(env_stack) <- env_names



### lower their resolution


### mask both high-res and low-res layers and save the output












## some random code below that I (Jan) used when changing environmental .grd files to .tif 
# env_list<-list.files(path = env.derived, pattern="tif$",  full.names=TRUE) 
# #subset to  "shelf" files
# env_list<-env_list[grep(".500m_shelf", env_list)]
# env_list<-env_list[-grep(".500m_shelf_scaled", env_list)]
# #for the single rasters layer names are missing. Extract from file name.
# env_names<-gsub(".*_|\\..*","",env_list)
# 
# for(i in 1:length(env_names)){
#   ra.loop <- raster(env_list[[i]])
#   writeRaster(ra.loop, sub(".gri*", ".tif", env_list[[i]]))
# }
# 
# 
# nams <- c("SST_mean","SST_sd","SST_sp_mean","SST_sp_sd","SST_su_mean","SST_su_sd")
# 
# ra.loop <- stack(env_list[[1]])
# for(i in 1:nlayers(ra.loop)){
#   writeRaster(ra.loop[[i]], sub("SST.gri*", paste0(nams[i],".tif"), env_list[[1]]))
# }
# 
# #stack all environmental layers and make sure they have appropriate names (currently manual and a bit messy!)
# env_stack<-raster::stack(env_list)
# 
# env_list<-list.files(path = env.derived, pattern="gri$",  full.names=TRUE) 
# #subset to  "shelf" files
# env_list<-env_list[grep(".500m_shelf_scaled", env_list)]
# #for the single rasters layer names are missing. Extract from file name.
# env_names<-gsub(".*_|\\..*","",env_list)
# 
# for(i in 1:length(env_names)){
#   ra.loop <- raster(env_list[[i]])
#   writeRaster(ra.loop, sub(".gri*", ".tif", env_list[[i]]))
# }







