#### reproject NPP to 2km resolution file

#################################################################
#### create and save npp-based seeding file for each region
## Options:
## - 1 particle with individualised value per cell: might get cells with 0 particles -> not good
## - 4 particles with individualised value per cell
## - 9  particles with individualised value per cell
## - make the # of particles dependent on cell value UNFEASIBLE TO RUN, because it needs 50 million particles for the particle distribution to be closely related to the NPP values in the Ross Sea

## specify user and setup directory to look up data from
usr <- "VM"
source("prep_environment/EnvPrep_0_SourceFile.R")

## set input and output folders
env.dir <- paste0(usr.main.dir,"/data_environmental/derived/ROMS_2k_files/")
out.dir <- paste0(usr.main.dir,"/data_environmental/derived/NPP/")
roms.dir <- paste0(usr.roms.dir,"/data_environmental/raw/ROMS_2k_files/")

########################################################
library(spatstat)
library(ppmData) ## for "terra2im"

npp.ca <- rast(paste0(out.dir,"NPP_climatology_OctMar_2002To2020_filled12boxes_cafe_mean.tif"))
npp.ep <- rast(paste0(out.dir,"NPP_climatology_OctMar_2002To2020_filled12boxes_eppl_mean.tif"))
npp.vp <- rast(paste0(out.dir,"NPP_climatology_OctMar_2002To2020_filled12boxes_vpmg_mean.tif"))
npp.cb <- rast(paste0(out.dir,"NPP_climatology_OctMar_2002To2020_filled12boxes_cbpm_mean.tif"))
chla   <- rast(paste0(out.dir,"NPP_climatology_OctMar_2002To2020_chla_mean_strict.tif"))
names(chla) <- "mean"
ra <- rast(paste0(roms.dir,"ocean_his_0001_slices_u_1.tif"), lyrs=1)

model <- "cafe_12boxfilled"
npp2 <- project(npp.ca, ra)
model <- "eppl_12boxfilled"
npp2 <- project(npp.ep, ra)
model <- "vpmg_12boxfilled"
npp2 <- project(npp.vp, ra)
model <- "cbpm_12boxfilled"
npp2 <- project(npp.cb, ra)
model <- "chla"
npp2 <- project(chla, ra)


#### - 9  particles with individualised value per cell
## load ROMS and NPP data
load(paste0(roms.dir,"ocean_his_TrackingSetup_24hcurrents_File01_Region01.Rdata"))
roms.coords.proj <- cbind(c(Rdat$x), c(Rdat$y))
x.range <- round(c(min(roms.coords.proj[,1])-1000,max(roms.coords.proj[,1])+1000),0)
y.range <- round(c(min(roms.coords.proj[,2])-1000,max(roms.coords.proj[,2])+1000),0)
ra.region <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000)
## particles locations
ra.9 <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000/3)
ra.4 <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000/2)
pts9.raw <- data.frame(cbind(crds(ra.9),z=0))
pts4.raw <- data.frame(cbind(crds(ra.4),z=0))
## crop npp to boundaries, and set irrelevant regions to NA
npp2_crop <- crop(npp2, ext(ra.region))
npp2_crop[400:650,1:100] <- NA
npp2_crop[500:650,1:150] <- NA
h.ra <- setValues(ra.region, Rdat$h)
plot(npp2_crop)
contour(h.ra, add=TRUE, levels=c(0,2000,3000))

## particle values
npp.vals9 <- terra::extract(npp2_crop, pts9.raw[,1:2])
pts9.raw$npp <- npp.vals9$mean/9
pts9 <- as.matrix(pts9.raw[-which(is.na(pts9.raw$npp)),])
npp.vals4 <- terra::extract(npp2_crop, pts4.raw[,1:2])
pts4.raw$npp <- npp.vals4$mean/4
pts4 <- as.matrix(pts4.raw[-which(is.na(pts4.raw$npp)),])

save(pts4, file=paste0(out.dir,"ocean_his_TrackingSetup_",model,"_NPP4_Region01.Rdata"))
save(pts9, file=paste0(out.dir,"ocean_his_TrackingSetup_",model,"_NPP9_Region01.Rdata"))

## save input rasters and boundaries for comparison later 
writeRaster(npp2_crop, file=paste0(out.dir,"ocean_his_TrackingSetup_",model,"_NPP_Region01.tif"))

#### SETUP NPP9 FOR ALL REGIONS

#### region 2
load(paste0(roms.dir,"ocean_his_TrackingSetup_24hcurrents_File01_Region02.Rdata"))
roms.coords.proj <- cbind(c(Rdat$x), c(Rdat$y))
x.range <- round(c(min(roms.coords.proj[,1])-1000,max(roms.coords.proj[,1])+1000),0)
y.range <- round(c(min(roms.coords.proj[,2])-1000,max(roms.coords.proj[,2])+1000),0)
ra.region <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000)
## particles locations
ra.9 <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000/3)
pts9.raw <- data.frame(cbind(crds(ra.9),z=0))
npp2_crop <- crop(npp2, ext(ra.region))
## remove unnecessary particles and check distribution
npp2_crop[330:400,c(1:150,320:800)] <- NA
npp2_crop[1:180,1:70] <- NA
h.ra <- setValues(ra.region, Rdat$h)
plot(npp2_crop)
contour(h.ra, add=TRUE, levels=c(0,2000,3000))
## particle values
npp.vals9 <- terra::extract(npp2_crop, pts9.raw[,1:2])
pts9.raw$npp <- npp.vals9$mean/9
pts9 <- as.matrix(pts9.raw[-which(is.na(pts9.raw$npp)),])
save(pts9, file=paste0(out.dir,"ocean_his_TrackingSetup_",model,"_NPP9_Region02.Rdata"))
writeRaster(npp2_crop, file=paste0(out.dir,"ocean_his_TrackingSetup_",model,"_NPP_Region02.tif"))

#### region 3
load(paste0(roms.dir,"ocean_his_TrackingSetup_24hcurrents_File01_Region03.Rdata"))
roms.coords.proj <- cbind(c(Rdat$x), c(Rdat$y))
x.range <- round(c(min(roms.coords.proj[,1])-1000,max(roms.coords.proj[,1])+1000),0)
y.range <- round(c(min(roms.coords.proj[,2])-1000,max(roms.coords.proj[,2])+1000),0)
ra.region <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000)
## particles locations
ra.9 <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000/3)
pts9.raw <- data.frame(cbind(crds(ra.9),z=0))
npp2_crop <- crop(npp2, ext(ra.region))
## remove unnecessary particles and check distribution
npp2_crop[150:950,1:50] <- NA
npp2_crop[530:950,1:170] <- NA
npp2_crop[800:950,1:250] <- NA
h.ra <- setValues(ra.region, Rdat$h)
plot(npp2_crop)
contour(h.ra, add=TRUE, levels=c(0,2000,3000))
## particle values
npp.vals9 <- terra::extract(npp2_crop, pts9.raw[,1:2])
pts9.raw$npp <- npp.vals9$mean/9
pts9 <- as.matrix(pts9.raw[-which(is.na(pts9.raw$npp)),])
save(pts9, file=paste0(out.dir,"ocean_his_TrackingSetup_",model,"_NPP9_Region03.Rdata"))
writeRaster(npp2_crop, file=paste0(out.dir,"ocean_his_TrackingSetup_",model,"_NPP_Region03.tif"))

#### region 4
load(paste0(roms.dir,"ocean_his_TrackingSetup_24hcurrents_File01_Region04.Rdata"))
roms.coords.proj <- cbind(c(Rdat$x), c(Rdat$y))
x.range <- round(c(min(roms.coords.proj[,1])-1000,max(roms.coords.proj[,1])+1000),0)
y.range <- round(c(min(roms.coords.proj[,2])-1000,max(roms.coords.proj[,2])+1000),0)
ra.region <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000)
## particles locations
ra.9 <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000/3)
pts9.raw <- data.frame(cbind(crds(ra.9),z=0))
npp2_crop <- crop(npp2, ext(ra.region))
## remove unnecessary particles and check distribution
npp2_crop[1:900,1:30] <- NA
npp2_crop[c(1:350,700:900),1:70] <- NA
npp2_crop[c(1:250,800:900),1:100] <- NA
npp2_crop[c(1:200,850:900),1:130] <- NA
npp2_crop[1:150,1:200] <- NA
npp2_crop[1:100,1:250] <- NA
npp2_crop[1:50,1:300] <- NA
h.ra <- setValues(ra.region, Rdat$h)
plot(npp2_crop)
contour(h.ra, add=TRUE, levels=c(0,2000,3000))
## particle values
npp.vals9 <- terra::extract(npp2_crop, pts9.raw[,1:2])
pts9.raw$npp <- npp.vals9$mean/9
pts9 <- as.matrix(pts9.raw[-which(is.na(pts9.raw$npp)),])
save(pts9, file=paste0(out.dir,"ocean_his_TrackingSetup_",model,"_NPP9_Region04.Rdata"))
writeRaster(npp2_crop, file=paste0(out.dir,"ocean_his_TrackingSetup_",model,"_NPP_Region04.tif"))

#### region 5
load(paste0(roms.dir,"ocean_his_TrackingSetup_24hcurrents_File01_Region05.Rdata"))
roms.coords.proj <- cbind(c(Rdat$x), c(Rdat$y))
x.range <- round(c(min(roms.coords.proj[,1])-1000,max(roms.coords.proj[,1])+1000),0)
y.range <- round(c(min(roms.coords.proj[,2])-1000,max(roms.coords.proj[,2])+1000),0)
ra.region <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000)
## particles locations
ra.9 <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000/3)
pts9.raw <- data.frame(cbind(crds(ra.9),z=0))
npp2_crop <- crop(npp2, ext(ra.region))
## remove unnecessary particles and check distribution
npp2_crop[400:600,1:100] <- NA
npp2_crop[1:60,50:700] <- NA
npp2_crop[1:100,100:650] <- NA
npp2_crop[1:130,150:550] <- NA
npp2_crop[1:180,200:350] <- NA
h.ra <- setValues(ra.region, Rdat$h)
plot(npp2_crop)
contour(h.ra, add=TRUE, levels=c(0,2000,3000))
## particle values
npp.vals9 <- terra::extract(npp2_crop, pts9.raw[,1:2])
pts9.raw$npp <- npp.vals9$mean/9
pts9 <- as.matrix(pts9.raw[-which(is.na(pts9.raw$npp)),])
save(pts9, file=paste0(out.dir,"ocean_his_TrackingSetup_",model,"_NPP9_Region05.Rdata"))
writeRaster(npp2_crop, file=paste0(out.dir,"ocean_his_TrackingSetup_",model,"_NPP_Region05.tif"))

#### region 6
load(paste0(roms.dir,"ocean_his_TrackingSetup_24hcurrents_File01_Region06.Rdata"))
roms.coords.proj <- cbind(c(Rdat$x), c(Rdat$y))
x.range <- round(c(min(roms.coords.proj[,1])-1000,max(roms.coords.proj[,1])+1000),0)
y.range <- round(c(min(roms.coords.proj[,2])-1000,max(roms.coords.proj[,2])+1000),0)
ra.region <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000)
## particles locations
ra.9 <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000/3)
pts9.raw <- data.frame(cbind(crds(ra.9),z=0))
npp2_crop <- crop(npp2, ext(ra.region))
## remove unnecessary particles and check distribution
npp2_crop[1:100,c(1:550,800:1400)] <- NA
npp2_crop[1:150,c(1:200,800:1100,1250:1400)] <- NA
npp2_crop[1:200,c(1:150,830:1050,1250:1400)] <- NA
npp2_crop[201:250,c(1:100,1200:1400)] <- NA
npp2_crop[251:320,c(1:50,1200:1400)] <- NA
h.ra <- setValues(ra.region, Rdat$h)
plot(npp2_crop)
contour(h.ra, add=TRUE, levels=c(0,2000,3000))
## particle values
npp.vals9 <- terra::extract(npp2_crop, pts9.raw[,1:2])
pts9.raw$npp <- npp.vals9$mean/9
pts9 <- as.matrix(pts9.raw[-which(is.na(pts9.raw$npp)),])
save(pts9, file=paste0(out.dir,"ocean_his_TrackingSetup_",model,"_NPP9_Region06.Rdata"))
writeRaster(npp2_crop, file=paste0(out.dir,"ocean_his_TrackingSetup_",model,"_NPP_Region06.tif"))

#### region 7
load(paste0(roms.dir,"ocean_his_TrackingSetup_24hcurrents_File01_Region07.Rdata"))
roms.coords.proj <- cbind(c(Rdat$x), c(Rdat$y))
x.range <- round(c(min(roms.coords.proj[,1])-1000,max(roms.coords.proj[,1])+1000),0)
y.range <- round(c(min(roms.coords.proj[,2])-1000,max(roms.coords.proj[,2])+1000),0)
ra.region <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000)
## particles locations
ra.9 <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000/3)
pts9.raw <- data.frame(cbind(crds(ra.9),z=0))
npp2_crop <- crop(npp2, ext(ra.region))
## remove unnecessary particles and check distribution
npp2_crop[1:100,250:500] <- NA
npp2_crop[1:300,300:500] <- NA
npp2_crop[1:400,400:500] <- NA
npp2_crop[1:700,450:500] <- NA
h.ra <- setValues(ra.region, Rdat$h)
plot(npp2_crop)
contour(h.ra, add=TRUE, levels=c(0,2000,3000))
## particle values
npp.vals9 <- terra::extract(npp2_crop, pts9.raw[,1:2])
pts9.raw$npp <- npp.vals9$mean/9
pts9 <- as.matrix(pts9.raw[-which(is.na(pts9.raw$npp)),])
save(pts9, file=paste0(out.dir,"ocean_his_TrackingSetup_",model,"_NPP9_Region07.Rdata"))
writeRaster(npp2_crop, file=paste0(out.dir,"ocean_his_TrackingSetup_",model,"_NPP_Region07.tif"))

#### region 8
load(paste0(roms.dir,"ocean_his_TrackingSetup_24hcurrents_File01_Region08.Rdata"))
roms.coords.proj <- cbind(c(Rdat$x), c(Rdat$y))
x.range <- round(c(min(roms.coords.proj[,1])-1000,max(roms.coords.proj[,1])+1000),0)
y.range <- round(c(min(roms.coords.proj[,2])-1000,max(roms.coords.proj[,2])+1000),0)
ra.region <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000)
## particles locations
ra.9 <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000/3)
pts9.raw <- data.frame(cbind(crds(ra.9),z=0))
npp2_crop <- crop(npp2, ext(ra.region))
## remove unnecessary particles and check distribution
npp2_crop[1:70,430:500] <- NA
npp2_crop[500:900,450:500] <- NA
npp2_crop[600:900,380:500] <- NA
npp2_crop[800:900,250:500] <- NA
h.ra <- setValues(ra.region, Rdat$h)
plot(npp2_crop)
contour(h.ra, add=TRUE, levels=c(0,2000,3000))
## particle values
npp.vals9 <- terra::extract(npp2_crop, pts9.raw[,1:2])
pts9.raw$npp <- npp.vals9$mean/9
pts9 <- as.matrix(pts9.raw[-which(is.na(pts9.raw$npp)),])
save(pts9, file=paste0(out.dir,"ocean_his_TrackingSetup_",model,"_NPP9_Region08.Rdata"))
writeRaster(npp2_crop, file=paste0(out.dir,"ocean_his_TrackingSetup_",model,"_NPP_Region08.tif"))

#### region 9
load(paste0(roms.dir,"ocean_his_TrackingSetup_24hcurrents_File01_Region09.Rdata"))
roms.coords.proj <- cbind(c(Rdat$x), c(Rdat$y))
x.range <- round(c(min(roms.coords.proj[,1])-1000,max(roms.coords.proj[,1])+1000),0)
y.range <- round(c(min(roms.coords.proj[,2])-1000,max(roms.coords.proj[,2])+1000),0)
ra.region <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000)
## particles locations
ra.9 <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000/3)
pts9.raw <- data.frame(cbind(crds(ra.9),z=0))
npp2_crop <- crop(npp2, ext(ra.region))
## remove unnecessary particles and check distribution
npp2_crop[500:600,350:800] <- NA
npp2_crop[450:600,500:800] <- NA
npp2_crop[400:600,600:800] <- NA
npp2_crop[350:600,650:800] <- NA
npp2_crop[300:600,700:800] <- NA
npp2_crop[200:600,750:800] <- NA
h.ra <- setValues(ra.region, Rdat$h)
plot(npp2_crop)
contour(h.ra, add=TRUE, levels=c(0,2000,3000))
## particle values
npp.vals9 <- terra::extract(npp2_crop, pts9.raw[,1:2])
pts9.raw$npp <- npp.vals9$mean/9
pts9 <- as.matrix(pts9.raw[-which(is.na(pts9.raw$npp)),])
save(pts9, file=paste0(out.dir,"ocean_his_TrackingSetup_",model,"_NPP9_Region09.Rdata"))
writeRaster(npp2_crop, file=paste0(out.dir,"ocean_his_TrackingSetup_",model,"_NPP_Region09.tif"))

#### region 10
load(paste0(roms.dir,"ocean_his_TrackingSetup_24hcurrents_File01_Region10.Rdata"))
roms.coords.proj <- cbind(c(Rdat$x), c(Rdat$y))
x.range <- round(c(min(roms.coords.proj[,1])-1000,max(roms.coords.proj[,1])+1000),0)
y.range <- round(c(min(roms.coords.proj[,2])-1000,max(roms.coords.proj[,2])+1000),0)
ra.region <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000)
## particles locations
ra.9 <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000/3)
pts9.raw <- data.frame(cbind(crds(ra.9),z=0))
npp2_crop <- crop(npp2, ext(ra.region))
h.ra <- setValues(ra.region, Rdat$h)
plot(npp2_crop)
contour(h.ra, add=TRUE, levels=c(0,2000,3000))
## particle values
npp.vals9 <- terra::extract(npp2_crop, pts9.raw[,1:2])
pts9.raw$npp <- npp.vals9$mean/9
pts9 <- as.matrix(pts9.raw[-which(is.na(pts9.raw$npp)),])
save(pts9, file=paste0(out.dir,"ocean_his_TrackingSetup_",model,"_NPP9_Region10.Rdata"))
writeRaster(npp2_crop, file=paste0(out.dir,"ocean_his_TrackingSetup_",model,"_NPP_Region10.tif"))

# h.ra <- setValues(ra.region, Rdat$h)
# plot(npp2_crop)
# contour(h.ra, add=TRUE, levels=c(0,2000,3000))





