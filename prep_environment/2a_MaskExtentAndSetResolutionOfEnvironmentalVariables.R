
## testing how many cells we have in the rasters for different resolutions (once off, now at the bottom)
## masking the extent of the environmental layers
## aggregating data from 500m resolution to 2km resolution
## saving files into a single tif-file

library(terra)

##### Setting up----
user = "Jan"
#user = "charley"
#user="nicole"
if (user == "Jan") {
  sci.dir <-      "C:/Users/jjansen/OneDrive - University of Tasmania/science/"
  env.derived <-  paste0(sci.dir,"data_environmental/derived/")
  env.raw <- "E:/science/data_environmental/raw/"
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

#res <- "500m"
res <- "2km"

#######################################################
##### Bathymetry
##### - already resampled to 2km in script
##### - limit to 2500m depth
##### - then mask
#######################################################

## first, load bathy data and draw a mask to include/exclude areas
# bathy_list<-list.files(path = env.derived, pattern="tif$",  full.names=TRUE) 
# bathy_list<-bathy_list[grep(paste0(".",res,"_shelf_bathy"), bathy_list)]
# bathy_names <- gsub(".*_|\\..*","",bathy_list)
# r <- rast(bathy_list)
# names(r) <- bathy_names

## load the bedrock bathymetry data
r <- rast(paste0(env.derived,"Circumpolar_EnvData_",res,"_shelf_bathy_ibcso2bed.tif"))

## set cells that are below 2500m depth to NA
r.mask <- r$depth
r.mask[which(r$depth[]<(-2500))] <- NA
## load mask boundaries
t <- vect(paste0(env.derived,"Circumpolar_EnvData_mask_shelf.shp"))
r.mask2 <- mask(r.mask,t)

## mask bathymetry data
r2 <- mask(r,r.mask2)
## save output:
writeRaster(r2, filename=paste0(env.derived,"Circumpolar_EnvData_",res,"_shelf_mask_bathy_ibcso2bed.tif"), overwrite=TRUE)


# ## set cells that are below 2500m depth to NA
# sel.na <- which(r$depth[]<(-2500))
# for(i in 1:nlyr(r)){
#   print(i)
#   r[[i]][sel.na] <- NA
# }
## draw mask around the continental shelf break and save polygon
#t <- draw(x="polygon")
#writeVector(t,paste0(env.derived,"Circumpolar_EnvData_mask_shelf.shp"))
# ## set filenames to save to:
# savestring <- paste0(env.derived,"Circumpolar_EnvData_",res,"_shelf_mask_bathy_ibcso2_")
# 
# ## save bathymetry back to file
# for(i in 1:nlyr(r)){
#   print(i)
#   writeRaster(r2[[i]], filename=paste0(savestring,bathy_names[i],".tif"), overwrite=TRUE)
# }

#######################################################
##### ICE, SSH & SST
##### - resample to 2km
##### - limit to 2500m depth
##### - mask
#######################################################

ice_stack <- rast(paste0(env.raw, "Circumpolar_EnvData_ice.grd"))
sst_stack <- rast(paste0(env.raw, "Circumpolar_EnvData_SST.grd"))
ssh_stack <- rast(paste0(env.raw, "Circumpolar_EnvData_SSH.grd"))

## load bathy file:
bathy_shelf <- rast(paste0(env.derived, "Circumpolar_EnvData_",res,"_shelf_mask_bathy_ibcso2_depth.tif"))
 

if(ext(ice_stack)!=ext(bathy_shelf)){
  ice_stack <- crop(ice_stack, bathy_shelf)
}
if(res=="500m"){
  ice_stack<-terra::project(ice_stack, bathy_shelf)
  ice_stack <- resample(ice_stack, bathy_shelf)
  ice_stack_shelf<-mask(ice_stack, bathy_shelf)
}else{
  ## bring to resolution and 2500m depth range, and mask
ice_stack<-terra::project(ice_stack, bathy_shelf)
ice_stack_shelf<-mask(ice_stack, bathy_shelf)
}

sst_stack<-terra::project(sst_stack, bathy_shelf)
sst_stack_shelf<-mask(sst_stack, bathy_shelf)

ssh_stack<-terra::project(ssh_stack, bathy_shelf)
ssh_stack_shelf<-mask(ssh_stack, bathy_shelf)

## set filenames to save to:
savestring2 <- paste0(env.derived,"Circumpolar_EnvData_",res,"_shelf_mask_")

writeRaster(ice_stack_shelf, filename=paste0(savestring2,"ICE.tif"), overwrite=TRUE)
writeRaster(sst_stack_shelf, filename=paste0(savestring2,"SST.tif"), overwrite=TRUE)
writeRaster(ssh_stack_shelf, filename=paste0(savestring2,"SSH.tif"), overwrite=TRUE)

#######################################################
##### GEOMORPH & DISTANCE TO CANYONS
##### - already resampled in script
##### - limit to 2500m depth
##### - mask
#######################################################
## load bathy file:
bathy_shelf <- rast(paste0(env.derived, "Circumpolar_EnvData_",res,"_shelf_mask_bathy_ibcso2_depth.tif"))

## load files
if(res=="500m"){
  geomorph <- rast(paste0(env.derived, "Circumpolar_EnvData_geomorphology.tif"))
  dist2cany <- rast(paste0(env.derived, "Circumpolar_EnvData_500m_shelf_distance2canyons.tif"))
} else{
  geomorph <- rast(paste0(env.derived, "Circumpolar_EnvData_2km_geomorphology.tif"))
  dist2cany <- rast(paste0(env.derived, "Circumpolar_EnvData_2km_shelf_distance2canyons.tif"))
}

## mask to 2500m depth range and extent
geomorph_shelf<-mask(geomorph, bathy_shelf)
dist2cany_shelf<-mask(dist2cany, bathy_shelf)

## save to:
writeRaster(geomorph_shelf, filename=paste0(savestring2,"geomorphology.tif"), overwrite=TRUE)
writeRaster(dist2cany_shelf, filename=paste0(savestring2,"distance2canyons.tif"), overwrite=TRUE)

#######################################################
##### BGC
##### - project and resample to 2km
##### - limit to 2500m depth
##### - mask
#######################################################

arag <- rast(paste0(env.raw, "bot_arag_BSOSE.grd"))
NO3  <- rast(paste0(env.raw, "bot_NO3_BSOSE.grd"))
O2   <- rast(paste0(env.raw, "bot_O2_BSOSE.grd"))
PO4  <- rast(paste0(env.raw, "bot_PO4_BSOSE.grd"))

bsose_stack <- c(arag,NO3,O2,PO4)

## load bathy file:
bathy_shelf <- rast(paste0(env.derived, "Circumpolar_EnvData_",res,"_shelf_mask_bathy_ibcso2_depth.tif"))

## crop rasters to SO extent, project to polar stereographic, then crop to shelf.
SO.ext2 <- ext(0, 360, -80, -55)
bsose_stack<-crop(bsose_stack, SO.ext2)

#because 0 to 360 degrees- projection doesn't work very well. Clumsy fix
#mean
bsose_temp<-raster::shift(bsose_stack, dx = -360)
origin(bsose_temp)<-origin(bsose_stack)
bsose_stack<-raster::merge(bsose_stack, bsose_temp)

## project and resample
bsose_res <- terra::project(bsose_stack, bathy_shelf)

## mask to 2500m depth range and extent
bsose_shelf<-mask(bsose_res, bathy_shelf)

## save to:
writeRaster(bsose_shelf, filename=paste0(savestring2,"bsose.tif"))

#######################################################
##### NPP and WAOM-OUTPUT
#######################################################

npp <- rast(c(paste0(env.derived,"Circumpolar_EnvData_",res,"_shelf_NPP_Cafe_filled_SummerAverage.tif"),
                 paste0(env.derived,"Circumpolar_EnvData_",res,"_shelf_NPP_Cafe_filled_SummerStandardDeviation.tif")))
names(npp) <- c("npp_mean", "npp_sd")

waom_list<-list.files(path = env.derived, pattern="tif$",  full.names=TRUE) 
waom_list<-waom_list[grep(paste0(".",res,"_shelf_waom"), waom_list)]
waom_names <- gsub(".*waom4k_|\\..*","",waom_list)

waom <- rast(waom_list)
names(waom) <- waom_names

## mask to 2500m depth range and extent
npp_shelf<-mask(npp, bathy_shelf)
waom_shelf<-mask(waom, bathy_shelf)

## save to:
writeRaster(npp_shelf, filename=paste0(savestring2,"NPP.tif"))
writeRaster(waom_shelf, filename=paste0(savestring2,"waom4k.tif"))

#######################################################

##### read in all files and save as one single tif:
env_list<-list.files(path = env.derived, pattern="tif$",  full.names=TRUE) 
#subset to  "shelf" files
env_list<-env_list[grep(paste0(".",res,"_shelf_mask"), env_list)]
# env_list<-env_list[-grep(paste0(".",res,"_shelf_mask_scaled"), env_list)]
# env_list<-env_list[-grep("polynomials", env_list)]

env_stack <- rast(env_list[-c(13)])
writeRaster(env_stack, filename=paste0(savestring2,"unscaled_variables.tif"), overwrite=TRUE)

## remove files, commented to stop accidently deleting things
#file.remove(env_list)

#env_stack2 <- c(subset(test,1:25),subset(env_stack,1:2),subset(test,28:39),subset(env_stack,3:9))

#######################################################
##### FAM and 2km model current speeds (THIS FIXES )
#######################################################
fam.dir <- paste0(sci.dir,"data_environmental/FAM_outputs/")
roms.dir <- paste0(sci.dir,"data_environmental/ROMS_2k_files/")

flux01 <- rast(paste0(fam.dir,"tracking2D_NPP9_200mday_21days_traj_r0001_28days_flux.tif"))
flux02 <- rast(paste0(fam.dir,"tracking2D_NPP9_200mday_21days_traj_r0002_28days_flux.tif"))
flux05 <- rast(paste0(fam.dir,"tracking2D_NPP9_200mday_21days_traj_r0005_28days_flux.tif"))
flux005 <- rast(paste0(fam.dir,"tracking2D_NPP9_200mday_21days_traj_r00005_28days_flux.tif"))
flux.mean <- mean(flux005,flux01,flux02, flux05, na.rm=TRUE)
log.flux.mean <- log(flux.mean)
sed01 <- rast(paste0(fam.dir,"tracking2D_NPP9_200mday_21days_traj_r0001_28days_sed.tif"))
sed02 <- rast(paste0(fam.dir,"tracking2D_NPP9_200mday_21days_traj_r0002_28days_sed.tif"))
sed05 <- rast(paste0(fam.dir,"tracking2D_NPP9_200mday_21days_traj_r0005_28days_sed.tif"))
sed005 <- rast(paste0(fam.dir,"tracking2D_NPP9_200mday_21days_traj_r00005_28days_sed.tif"))
sed.mean <- mean(sed005,sed01,sed02, sed05, na.rm=TRUE)
fam <- c(flux005,flux01,flux02,flux05,flux.mean,log.flux.mean, sed005,sed01,sed02,sed05,sed.mean)
names(fam) <- c("flux00005","flux0001","flux0002","flux0005","flux.mean","log.flux.mean","sed00005","sed0001","sed0002","sed0005","sed.mean")

## these are averages for a single month in summer
uv_absmean <- rast(paste0(roms.dir,"ocean_his_bottom_uv_absmean.tif"))
uv_mean <- rast(paste0(roms.dir,"ocean_his_bottom_uv_mean.tif"))
uv_max <- rast(paste0(roms.dir,"ocean_his_bottom_uv_max.tif"))
w_absmean <- rast(paste0(roms.dir,"ocean_his_bottom_w_absmean.tif"))
w_mean <- rast(paste0(roms.dir,"ocean_his_bottom_w_mean.tif"))
currents <- c(uv_absmean,uv_mean,uv_max,w_absmean,w_mean)
names(currents) <- c("uv_absmean","uv_mean","uv_max","w_absmean","w_mean")

## need to match extents of the rasters first before masking
fam.resampled <- resample(fam, r.mask2, method="near")
currents.resampled <- resample(currents, r.mask2, method="near")

#plot(fam.resampled)
sed.na <- which(is.na(fam.resampled$sed.mean[]))
flux.na <- which(is.na(fam.resampled$flux.mean[]))
fam.resampled$sed.mean[sed.na[-which(sed.na%in%flux.na)]] <- 0

fam2 <- mask(fam.resampled, r.mask2)
currents2 <- mask(currents.resampled, r.mask2)


###############################
unscaled.vars <- rast(paste0(env.derived,"Circumpolar_EnvData_",res,"_shelf_mask_unscaled_variables_old.tif"))
bedbathy <- rast(paste0(env.derived,"Circumpolar_EnvData_",res,"_shelf_mask_bathy_ibcso2bed.tif"))

unscaled.vars[[1:5]] <- bedbathy[[1:5]]

## drop sd and substitute current vals
unscaled.vars <- unscaled.vars[[-44]]
unscaled.vars$seafloorcurrents_max <- currents2$uv_max
unscaled.vars$seafloorcurrents_mean <- currents2$uv_mean
unscaled.vars$seafloorcurrents_absolute <- currents2$uv_absmean
unscaled.vars$seafloorcurrents_residual <- currents2$uv_absmean-currents2$uv_mean

unscaled.vars <- c(unscaled.vars,fam2)

#writeRaster(unscaled.vars, filename=paste0(env.derived,"Circumpolar_EnvData_",res,"_shelf_mask_unscaled_variables.tif"), overwrite=TRUE)
























##########
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

##########
## how many cells per resolution?
# ## restrict depth range and mask to remove areas off the shelf
# r3 <- r2
# r3[r2<(-2500)] <- NA
# r4 <- mask(r3,t)
# r5 <- aggregate(r4,4)
# r6 <- aggregate(r4,8)
# 
# ## 500m res, 0-3000m depth, raw,    -> 23.6 million cells
# length(which(!is.na(r2[])))
# ## 500m res, 0-2500m depth, raw,    -> 17.6 million cells
# length(which(!is.na(r3[])))
# ## 500m res, 0-2500m depth, masked, -> 16.1 million cells
# length(which(!is.na(r4[])))
# ## 2km res,  0-2500m depth, masked, -> 987k cells
# length(which(!is.na(r5[])))
# ## 4km res,  0-2500m depth, masked, -> 241k cells
# length(which(!is.na(r6[])))
