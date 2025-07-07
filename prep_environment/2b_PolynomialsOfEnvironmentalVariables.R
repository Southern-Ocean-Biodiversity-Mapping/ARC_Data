#######################################################################################################
##### This code creates polynomials of environmental rasters
#######################################################################################################

## 1) set up----
library(terra)

user = "Jan"
#user = "charley"
#user="nicole"

if (user == "Jan") {
  
  sci.dir <-      "C:/Users/jjansen/OneDrive - University of Tasmania/science/"
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

####################################

#res = "500m"
res = "2km"

###################################

## 2) get file names of all environmental rasters and bricks and load into one big stack----
env_stack <- rast(paste0(env.derived,"Circumpolar_EnvData_",res,"_shelf_mask_unscaled_variables.tif"))

## creating polynomials of raster layers
depth2 <- rast(env_stack$depth)
sel <- which(!is.na(env_stack$depth[]))
depth2.dat <- poly(c(env_stack$depth[sel])[[1]],2)[,2]
depth2[sel] <- depth2.dat

distance2canyons2 <- rast(env_stack$distance2canyons)
sel <- which(!is.na(env_stack$distance2canyons[]))
distance2canyons2.dat <- poly(c(env_stack$distance2canyons[sel])[[1]],2)[,2]
distance2canyons2[sel] <- distance2canyons2.dat

logslope <- log(env_stack$slope)

poly_stack <- c(depth2, distance2canyons2, logslope)
names(poly_stack) <- c("depth2", "distance2canyons2", "logslope")
writeRaster(poly_stack, filename=paste0(env.derived,"Circumpolar_EnvData_",res,"_shelf_mask_unscaled_polynomials_etc.tif"), overwrite=TRUE)

