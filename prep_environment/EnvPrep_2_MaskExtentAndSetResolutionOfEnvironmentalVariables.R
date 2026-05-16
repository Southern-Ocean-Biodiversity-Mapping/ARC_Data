
## testing how many cells we have in the rasters for different resolutions (once off, now at the bottom)
## masking the extent of the environmental layers
## aggregating data from 500m resolution to 2km resolution
## saving files into a single tif-file

## specify user and setup directory to look up data from
usr <- "VM"
#usr <- "SJ"
#usr <- "JJ"
source("0_SourceFile.R")

## set folders
env.derived  <- paste0(usr.main.dir, "data_environmental/derived/")
env.raw      <- paste0(usr.main.dir, "data_environmental/raw/")
roms.dir     <- paste0(usr.roms.dir,"data_environmental/derived/ROMS_2k_files/")
roms.dir2    <- paste0(usr.dropbox.dir,"data_environmental/derived/ROMS/")
fam.dir    <- paste0(usr.dropbox.dir,"data_environmental/derived/FAM_outputs/")

##############################################################################################################
##############################################################################################################

#res <- "500m"
res <- "2km"

library(terra)

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
r <- rast(paste0(env.derived,"bathy_outputs/IBCSO_v2_2km_bathymetric_variables.tif"))

## draw mask around the continental shelf break and save polygon
#t <- draw(x="polygon")
#writeVector(t,paste0(env.derived,"Circumpolar_EnvData_mask_shelf.shp"))

## set cells that are below 2500m depth to NA
r.mask <- r$depth
r.mask[which(r$depth[]<(-2500))] <- NA
## load mask boundaries
t <- vect(paste0(env.derived,"bathy_outputs/Circumpolar_EnvData_mask_shelf.shp"))
crs(t) <- crs(r.mask)
r.mask2 <- mask(r.mask,t)

## mask bathymetry data
r2 <- mask(r,r.mask2)
## save output:
writeRaster(r2, filename=paste0(env.derived,"Circumpolar_EnvData_",res,"_shelf_mask_bathy_ibcso2bed.tif"), overwrite=TRUE)

#######################################################
##### ICE, SSH & SST
##### - resample to 2km
##### - limit to 2500m depth
##### - mask
#######################################################

# ice_stack <- rast(paste0(env.raw, "Circumpolar_EnvData_ice.grd"))
# sst_stack <- rast(paste0(env.raw, "Circumpolar_EnvData_SST.grd"))
# ssh_stack <- rast(paste0(env.raw, "Circumpolar_EnvData_SSH.grd"))
# 
# ## load bathy file:
# bathy_shelf <- rast(paste0(env.derived, "Circumpolar_EnvData_",res,"_shelf_mask_bathy_ibcso2_depth.tif"))
#  
# 
# if(ext(ice_stack)!=ext(bathy_shelf)){
#   ice_stack <- crop(ice_stack, bathy_shelf)
# }
# if(res=="500m"){
#   ice_stack<-terra::project(ice_stack, bathy_shelf)
#   ice_stack <- resample(ice_stack, bathy_shelf)
#   ice_stack_shelf<-mask(ice_stack, bathy_shelf)
# }else{
#   ## bring to resolution and 2500m depth range, and mask
# ice_stack<-terra::project(ice_stack, bathy_shelf)
# ice_stack_shelf<-mask(ice_stack, bathy_shelf)
# }
# 
# sst_stack<-terra::project(sst_stack, bathy_shelf)
# sst_stack_shelf<-mask(sst_stack, bathy_shelf)
# 
# ssh_stack<-terra::project(ssh_stack, bathy_shelf)
# ssh_stack_shelf<-mask(ssh_stack, bathy_shelf)
# 
# ## set filenames to save to:
savestring2 <- paste0(env.derived,"Circumpolar_EnvData_",res,"_shelf_mask_")
# 
# writeRaster(ice_stack_shelf, filename=paste0(savestring2,"ICE.tif"), overwrite=TRUE)
# writeRaster(sst_stack_shelf, filename=paste0(savestring2,"SST.tif"), overwrite=TRUE)
# writeRaster(ssh_stack_shelf, filename=paste0(savestring2,"SSH.tif"), overwrite=TRUE)

#######################################################
##### GEOMORPH & DISTANCE TO CANYONS
##### - already resampled in script
##### - limit to 2500m depth
##### - mask
#######################################################
## load bathy file:
bathy_shelf <- rast(paste0(env.derived, "Circumpolar_EnvData_",res,"_shelf_mask_bathy_ibcso2bed.tif"))

## load files
if(res=="500m"){
  # geomorph  <- rast(paste0(env.derived, "Circumpolar_EnvData_geomorphology.tif"))
  # dist2cany <- rast(paste0(env.derived, "Circumpolar_EnvData_500m_shelf_distance2canyons.tif"))
} else{
  geomorph  <- rast(paste0(env.derived, "bathy_outputs/IBCSO_v2_2km_geomorph.tif"))
  dist2cany <- rast(paste0(env.derived, "bathy_outputs/IBCSO_v2_2km_DistanceToCanyons.tif"))
}

## mask to 2500m depth range and extent
geomorph_shelf <-mask(geomorph,  bathy_shelf$depth)
dist2cany_shelf<-mask(dist2cany, bathy_shelf$depth)

names(dist2cany_shelf) <- "distance2canyons"

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
bathy_shelf <- rast(paste0(env.derived, "Circumpolar_EnvData_",res,"_shelf_mask_bathy_ibcso2bed.tif"))

## crop rasters to SO extent, project to polar stereographic, then crop to shelf.
SO.ext2 <- ext(0, 360, -80, -55)
bsose_stack<-crop(bsose_stack, SO.ext2)

#because 0 to 360 degrees- projection doesn't work very well. Clumsy fix
#mean
bsose_temp<-raster::shift(bsose_stack, dx = -360)
origin(bsose_temp)<-origin(bsose_stack)
bsose_stack<-raster::merge(bsose_stack, bsose_temp)

## project and resample
bsose_res <- terra::project(bsose_stack, bathy_shelf$depth)

## mask to 2500m depth range and extent
bsose_shelf<-mask(bsose_res, bathy_shelf$depth)

## save to:
writeRaster(bsose_shelf, filename=paste0(savestring2,"bsose.tif"), overwrite=TRUE)

#######################################################
##### NPP and WAOM-OUTPUT
#######################################################
l <- list.files(paste0(env.derived,"NPP/"), pattern="NPP_climatology_OctMar_2002To2020_filled12boxes", full.names = TRUE)
l.names.1 <- gsub("NPP_climatology_OctMar_2002To2020_filled12boxes_","",basename(l))
l.names <- gsub(".tif","",l.names.1)
npp <- rast(l)
names(npp) <- l.names

## 4k res temp and sal
waom_list<-list.files(path = roms.dir2, pattern="tif$",  full.names=TRUE) 
waom_list<-waom_list[grep(res, waom_list)]
# waom_list<-waom_list[grep(paste0(".",res,"_shelf_waom"), waom_list)]
waom_names <- gsub(".*waom4k_|\\..*","",waom_list)

waom <- rast(waom_list)
names(waom) <- waom_names

## project to SO
npp_proj <-project(npp,  r2$depth)
waom_proj<-project(waom, r2$depth)

## mask to 2500m depth range and extent
npp_shelf <-mask(npp_proj,  r2$depth)
waom_shelf<-mask(waom_proj, r2$depth)

## save to:
writeRaster(npp_shelf, filename=paste0(savestring2,"NPP_climatology_OctMar_2002To2020_filled12boxes.tif"), overwrite=TRUE)
writeRaster(waom_shelf, filename=paste0(savestring2,"waom4k_bottomtempsal.tif"), overwrite=TRUE)

#######################################################
##### FAM and 2km model current speeds
#######################################################
#### current speeds
## these are averages for a single month in summer
uv_absmean <- rast(paste0(roms.dir,"ocean_his_bottom_uv_absmean.tif"))
uv_mean    <- rast(paste0(roms.dir,"ocean_his_bottom_uv_mean.tif"))
uv_max     <- rast(paste0(roms.dir,"ocean_his_bottom_uv_max.tif"))
w_absmean  <- rast(paste0(roms.dir,"ocean_his_bottom_w_absmean.tif"))
w_mean     <- rast(paste0(roms.dir,"ocean_his_bottom_w_mean.tif"))
currents   <- c(uv_absmean,uv_mean,uv_max,w_absmean,w_mean)
names(currents) <- c("seafloorcurrents_absolutemean","seafloorcurrents_mean","seafloorcurrents_maximum","w_absmean","w_mean")
currents$seafloorcurrents_residual <- currents$seafloorcurrents_absolutemean - currents$seafloorcurrents_mean

#### food-availability simulations
all.fam.files <- list.files(fam.dir, full.names=TRUE)
fam.files.sed <- all.fam.files[grep("sed_circumpolar", all.fam.files)]
fam.files.flux<- all.fam.files[grep("flux_circumpolar", all.fam.files)]

flux <- rast(fam.files.flux)
names(flux) <- gsub("_circumpolar.tif","",basename(fam.files.flux))
sed <- rast(fam.files.sed)
names(sed) <- gsub("_circumpolar.tif","",basename(fam.files.sed))

sel.cafe <- grep("cafe", names(flux))
sel.cbpm <- grep("cbpm", names(flux))
sel.eppl <- grep("eppl", names(flux))
sel.vpmg <- grep("vpmg", names(flux))

flux.mean.cafe <- mean(flux[[sel.cafe]], na.rm=TRUE)
flux.mean.cbpm <- mean(flux[[sel.cbpm]], na.rm=TRUE)
flux.mean.eppl <- mean(flux[[sel.eppl]], na.rm=TRUE)
flux.mean.vpmg <- mean(flux[[sel.vpmg]], na.rm=TRUE)
log.flux.mean.cafe <- log(flux.mean.cafe)
log.flux.mean.cbpm <- log(flux.mean.cbpm)
log.flux.mean.eppl <- log(flux.mean.eppl)
log.flux.mean.vpmg <- log(flux.mean.vpmg)
sed.mean.cafe <- mean(sed[[sel.cafe]], na.rm=TRUE)
sed.mean.cbpm <- mean(sed[[sel.cbpm]], na.rm=TRUE)
sed.mean.eppl <- mean(sed[[sel.eppl]], na.rm=TRUE)
sed.mean.vpmg <- mean(sed[[sel.vpmg]], na.rm=TRUE)

fam <- c(flux.mean.cafe,flux.mean.cbpm,flux.mean.eppl,flux.mean.vpmg,
         log.flux.mean.cafe,log.flux.mean.cbpm,log.flux.mean.eppl,log.flux.mean.vpmg,
         sed.mean.cafe,sed.mean.cbpm,sed.mean.eppl,sed.mean.vpmg)
names(fam) <-  c("flux.mean.cafe","flux.mean.cbpm","flux.mean.eppl","flux.mean.vpmg",
                 "log.flux.mean.cafe","log.flux.mean.cbpm","log.flux.mean.eppl","log.flux.mean.vpmg",
                 "sed.mean.cafe","sed.mean.cbpm","sed.mean.eppl","sed.mean.vpmg")

## common projection
fam.proj <- project(fam, r2$depth)
currents.proj <- project(currents, r2$depth)

## replace NAs with zero where not on land and not deeper than 2500m
sel <- which(!is.na(r2$depth[]))
for(i in 1:nlyr(fam.proj)){
  print(i)
  sel.fam <- which(is.na(fam.proj[[i]][sel]))
  fam.proj[[i]][sel[sel.fam]] <- 0
}

fam2 <- mask(fam.proj, r.mask2)
currents2 <- mask(currents.proj, r.mask2)

## save to:
writeRaster(fam2,      filename=paste0(savestring2,"FAM_mean_12boxfilled_NPP9_200mday_21days_r00005tor0005_28days.tif"), overwrite=TRUE)
writeRaster(currents2, filename=paste0(savestring2,"waom2k_bottomcurrents.tif"), overwrite=TRUE)

#######################################################
##### read in all files and save as one single tif:
env_list<-list.files(path = env.derived, pattern="tif$",  full.names=TRUE) 
#subset to  "shelf" files
env_list<-env_list[grep(paste0(".",res,"_shelf_mask"), env_list)]
env_stack <- rast(env_list)
writeRaster(env_stack, filename=paste0(savestring2,"unscaled_variables.tif"), overwrite=TRUE)



