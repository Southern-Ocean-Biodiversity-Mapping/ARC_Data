#######################################################################################################
##### This code creates polynomials of environmental rasters
#######################################################################################################

## specify user and setup directory to look up data from
usr <- "VM"
#usr <- "SJ"
#usr <- "JJ
source("0_SourceFile.R")

## set folders
env.derived  <- paste0(usr.main.dir, "data_environmental/derived/")

##############################################################################################################
##############################################################################################################

#res <- "500m"
res <- "2km"

library(terra)

###################################
## get file names of all environmental rasters and load ----
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

