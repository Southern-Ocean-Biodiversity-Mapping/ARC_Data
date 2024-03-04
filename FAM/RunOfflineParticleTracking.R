
#title: "Run Offline Particle-tracking"
#author: "Jan Jansen"
#date: "2024"

#### 1) Set up ----
# load libraries
library(ncdf4)        ## package for netcdf manipulation
library(raadtools)
#library(raster)       ## package for raster manipulation
library(sp)
library(dplyr)
library(blueant)
library(rgdal)        ## package for geospatial analysis
library(ggplot2)      ## package for plotting
library(terra)

library(spatialEco)

# set up directory pointers etc
# Jan's local machine:
env.dir <- "C:/Users/jjansen/Desktop/science/data_environmental/"

# remote repository (DOESN'T WORK YET)
# env.dir <- "https://data.imas.utas.edu.au/data_transfer/admin/files/EnvironmentalData/"

env.raw <- "E:/science/data_environmental/raw/"
env.derived <- paste0(env.dir,"derived/")
AAD_dir <- "E:/science/data_environmental/raw/accessed_through_R"

string.chr <- "Circumpolar_EnvData_"
string.res <- "500m_"

# polar stereographic projection:
stereo <- "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

################################
#### ROMS Currents
## WE NEED A RUN ACROSS SUMMER WITH HIGH-RES HISTORY FILES!!! 8 has one month only, 8_Nov has 10-day intervals...
data.dat100 <- "E:/science/data_environmental/Circumpolar_ROMS/4km_outputs/output_yr10/"
## load lon/lat information from ROMS-grid
grd4k_nc <- nc_open(paste0(env.raw,"waom4extend_grd.nc"))
lon_rho <- ncvar_get(grd4k_nc, varid="lon_rho")
lat_rho <- ncvar_get(grd4k_nc, varid="lat_rho")
## Prepare empty rasters to assign correct projected values to
roms.coords.proj <- rgdal::project(cbind(c(lon_rho), c(lat_rho)), proj=stereo)
x.range <- c(min(roms.coords.proj[,1])-2000,max(roms.coords.proj[,1])+2000)
y.range <- c(min(roms.coords.proj[,2])-2000,max(roms.coords.proj[,2])+2000)
empty.roms.ra <- rast(extent=extent(c(x.range,y.range)), crs=stereo, resolution=4000)
# ## depth
# h <- rast(paste0(data.dat100,"ocean_avg_0001.nc"), subds="h")
# ## seafloor variables
# s <- seq(1,217, by=31)
# 
# # #seafloor currents (seafloor-layer is 1)
# # u.raw1 <- brick(paste0(data.dat100,"ocean_his_0002.nc"), varname="u", level=1)
# # v.raw1 <- brick(paste0(data.dat100,"ocean_his_0002.nc"), varname="v", level=1)
# # u.raw2 <- brick(paste0(data.dat100,"ocean_his_0003.nc"), varname="u", level=1)
# # v.raw2 <- brick(paste0(data.dat100,"ocean_his_0003.nc"), varname="v", level=1)
# # u.raw3 <- brick(paste0(data.dat100,"ocean_his_0004.nc"), varname="u", level=1)
# # v.raw3 <- brick(paste0(data.dat100,"ocean_his_0004.nc"), varname="v", level=1)
# # ## sum up monthly values for a climatology (THIS SHOULD BE DONE ON HIGH RESOLUTION HISTORY FILES)
# # u.sum <- sum(u.raw)
# # v.sum <- sum(v.raw)
# # u.sum.abs <- sum(abs(u.raw))
# # v.sum.abs <- sum(abs(v.raw))
# # u31.sum <- sum(u_31.raw)
# # v31.sum <- sum(v_31.raw)
# 
# #seafloor currents (seafloor-layer is 1)
# u.raw <- c(subset(rast(paste0(data.dat100,"ocean_avg_0001.nc"), subds="u"), subset=s),
#            subset(rast(paste0(data.dat100,"ocean_avg_0002.nc"), subds="u"), subset=s),
#            subset(rast(paste0(data.dat100,"ocean_avg_0003.nc"), subds="u"), subset=s),
#            subset(rast(paste0(data.dat100,"ocean_avg_0004.nc"), subds="u"), subset=s),
#            subset(rast(paste0(data.dat100,"ocean_avg_0005.nc"), subds="u"), subset=s),
#            subset(rast(paste0(data.dat100,"ocean_avg_0006.nc"), subds="u"), subset=s),
#            subset(rast(paste0(data.dat100,"ocean_avg_0007.nc"), subds="u"), subset=s),
#            subset(rast(paste0(data.dat100,"ocean_avg_0008.nc"), subds="u"), subset=s),
#            subset(rast(paste0(data.dat100,"ocean_avg_0009.nc"), subds="u"), subset=s),
#            subset(rast(paste0(data.dat100,"ocean_avg_0010.nc"), subds="u"), subset=s),
#            subset(rast(paste0(data.dat100,"ocean_avg_0011.nc"), subds="u"), subset=s),
#            subset(rast(paste0(data.dat100,"ocean_avg_0012.nc"), subds="u"), subset=s))
# 
# v.raw <- c(subset(rast(paste0(data.dat100,"ocean_avg_0001.nc"), subds="v"), subset=s),
#            subset(rast(paste0(data.dat100,"ocean_avg_0002.nc"), subds="v"), subset=s),
#            subset(rast(paste0(data.dat100,"ocean_avg_0003.nc"), subds="v"), subset=s),
#            subset(rast(paste0(data.dat100,"ocean_avg_0004.nc"), subds="v"), subset=s),
#            subset(rast(paste0(data.dat100,"ocean_avg_0005.nc"), subds="v"), subset=s),
#            subset(rast(paste0(data.dat100,"ocean_avg_0006.nc"), subds="v"), subset=s),
#            subset(rast(paste0(data.dat100,"ocean_avg_0007.nc"), subds="v"), subset=s),
#            subset(rast(paste0(data.dat100,"ocean_avg_0008.nc"), subds="v"), subset=s),
#            subset(rast(paste0(data.dat100,"ocean_avg_0009.nc"), subds="v"), subset=s),
#            subset(rast(paste0(data.dat100,"ocean_avg_0010.nc"), subds="v"), subset=s),
#            subset(rast(paste0(data.dat100,"ocean_avg_0011.nc"), subds="v"), subset=s),
#            subset(rast(paste0(data.dat100,"ocean_avg_0012.nc"), subds="v"), subset=s))
# 
# #seasurface currents (surface-layer is 31)
# s <- seq(31,217, by=31)
# u_31.raw <- c(subset(rast(paste0(data.dat100,"ocean_avg_0001.nc"), subds="u"), subset=s),
#               subset(rast(paste0(data.dat100,"ocean_avg_0002.nc"), subds="u"), subset=s),
#               subset(rast(paste0(data.dat100,"ocean_avg_0003.nc"), subds="u"), subset=s),
#               subset(rast(paste0(data.dat100,"ocean_avg_0004.nc"), subds="u"), subset=s),
#               subset(rast(paste0(data.dat100,"ocean_avg_0005.nc"), subds="u"), subset=s),
#               subset(rast(paste0(data.dat100,"ocean_avg_0006.nc"), subds="u"), subset=s),
#               subset(rast(paste0(data.dat100,"ocean_avg_0007.nc"), subds="u"), subset=s),
#               subset(rast(paste0(data.dat100,"ocean_avg_0008.nc"), subds="u"), subset=s),
#               subset(rast(paste0(data.dat100,"ocean_avg_0009.nc"), subds="u"), subset=s),
#               subset(rast(paste0(data.dat100,"ocean_avg_0010.nc"), subds="u"), subset=s),
#               subset(rast(paste0(data.dat100,"ocean_avg_0011.nc"), subds="u"), subset=s),
#               subset(rast(paste0(data.dat100,"ocean_avg_0012.nc"), subds="u"), subset=s))
# 
# v_31.raw <- c(subset(rast(paste0(data.dat100,"ocean_avg_0001.nc"), subds="v"), subset=s),
#               subset(rast(paste0(data.dat100,"ocean_avg_0002.nc"), subds="v"), subset=s),
#               subset(rast(paste0(data.dat100,"ocean_avg_0003.nc"), subds="v"), subset=s),
#               subset(rast(paste0(data.dat100,"ocean_avg_0004.nc"), subds="v"), subset=s),
#               subset(rast(paste0(data.dat100,"ocean_avg_0005.nc"), subds="v"), subset=s),
#               subset(rast(paste0(data.dat100,"ocean_avg_0006.nc"), subds="v"), subset=s),
#               subset(rast(paste0(data.dat100,"ocean_avg_0007.nc"), subds="v"), subset=s),
#               subset(rast(paste0(data.dat100,"ocean_avg_0008.nc"), subds="v"), subset=s),
#               subset(rast(paste0(data.dat100,"ocean_avg_0009.nc"), subds="v"), subset=s),
#               subset(rast(paste0(data.dat100,"ocean_avg_0010.nc"), subds="v"), subset=s),
#               subset(rast(paste0(data.dat100,"ocean_avg_0011.nc"), subds="v"), subset=s),
#               subset(rast(paste0(data.dat100,"ocean_avg_0012.nc"), subds="v"), subset=s))
# 
# ## interpolate u and v values at the rho points (where h, temp and salt are already defined)
# 
# ## extract current speeds at rho-points where depth is defined
# ## (in ROMS they are all at different locations):
# coord.grd.u <- coord.grd.v <- crds(h)
# ## u has one less column than the grid, so all coordinates need to be moved half a cell to the right for the grids to match up
# coord.grd.u[,1] <- crds(h)[,1]-0.5
# ## v has one less row than the grid, so all coordinates need to be moved half a cell up for the grids to match up
# coord.grd.v[,2] <- crds(h)[,2]-0.5
# ## now extract values at the rho-points and interpolate (because they are 2km away from the nearest original point), and place into projected raster
# u <- v <- u31 <- v31 <- rast(extent=extent(c(x.range,y.range)), crs=stereo, resolution=4000, nlyr=nlyr(u.raw))
# ## need to do this consecutively because of RAM
# u.dat <- extract(u.raw, coord.grd.u, method="bilinear")
# for(i in 1:ncol(u.dat)){
#   print(i)
#   u[[i]][] <- u.dat[,i]
# }
# rm(u.dat, u.raw)
# v.dat <- extract(v.raw, coord.grd.v, method="bilinear")
# for(i in 1:ncol(v.dat)){
#   print(i)
#   v[[i]][] <- v.dat[,i]
# }
# rm(v.dat, v.raw)
# u31.dat <- extract(u_31.raw, coord.grd.u, method="bilinear")
# for(i in 1:ncol(u31.dat)){
#   print(i)
#   u31[[i]][] <- u31.dat[,i]
# }
# rm(u31.dat, u_31.raw)
# v31.dat <- extract(v_31.raw, coord.grd.v, method="bilinear")
# for(i in 1:ncol(v31.dat)){
#   print(i)
#   v31[[i]][] <- v31.dat[,i]
# }
# rm(v31.dat, v_31.raw)
# 
# ## simple current speeds
# uv_31 <- sqrt(u31^2+v31^2)
# uv <- sqrt(u^2+v^2)
# ## and derivatives
# uv.max <- max(uv)
# uv.sd <- stdev(uv)
# 
# ## u and v derivatives
# u.mean <- mean(u)
# u.mean.abs <- mean(abs(u))
# v.mean <- mean(v)
# v.mean.abs <- mean(abs(v))
# ## uv based on mean u and v
# uv.mean <- sqrt(u.mean^2+v.mean^2)
# uv.abs.mean <- sqrt(u.mean.abs^2+v.mean.abs^2) ## this is essentially the same as mean(uv)
# # 
# ## write rasters to file
# writeRaster(u.mean,filename=paste0(env.derived,string.chr,"4km_waom4k_seafloorcurrents_mean_u.tif"))
# writeRaster(v.mean,filename=paste0(env.derived,string.chr,"4km_waom4k_seafloorcurrents_mean_v.tif"))
# writeRaster(uv.mean,filename=paste0(env.derived,string.chr,"4km_waom4k_seafloorcurrents_mean.tif"))
# writeRaster(uv.abs.mean, filename=paste0(env.derived,string.chr,"4km_waom4k_seafloorcurrents_mean_absolute.tif"))
# writeRaster(uv.max, filename=paste0(env.derived,string.chr,"4km_waom4k_seafloorcurrents_max.tif"))
# writeRaster(uv.sd,  filename=paste0(env.derived,string.chr,"4km_waom4k_seafloorcurrents_sd.tif"))

## load current speeds
load(file=paste0(env.derived,string.chr,"4km_waom4k_seafloorcurrents_mean_u.tif"))
load(file=paste0(env.derived,string.chr,"4km_waom4k_seafloorcurrents_mean_v.tif"))
load(file=paste0(env.derived,string.chr,"4km_waom4k_seafloorcurrents_mean.tif"))
load(file=paste0(env.derived,string.chr,"4km_waom4k_seafloorcurrents_mean_absolute.tif"))
load(file=paste0(env.derived,string.chr,"4km_waom4k_seafloorcurrents_max.tif"))
load(file=paste0(env.derived,string.chr,"4km_waom4k_seafloorcurrents_sd.tif"))

###################################
## load NPP
npp.raw <- rast(paste0(env.derived,string.chr,"NPP_Cafe_filled_SummerAverage.tif"))
npp <- project(npp.raw, uv.mean)

###################################
###### test particle tracking
#### load functions
ptrackr_dir <- "C:/Users/jjansen/Desktop/science/Scripts/ptrackR_package_update2024/ptrackr/R/"
## minimum for 3d run
source(paste0(ptrackr_dir,"create_points_pattern.R"))
source(paste0(ptrackr_dir,"setup_knn.R"))
source(paste0(ptrackr_dir,"trackit_3D.R"))
## for 2d run:
source(paste0(ptrackr_dir,"trackit_2D.R"))
source(paste0(ptrackr_dir,"buildparams.R"))

#### toyrun:
## setup ROMS field
load("C:/Users/jjansen/Desktop/science/Scripts/ptrackR_package_update2024/ptrackr/data/toyROMS.rdata")
## setup pattern of points to seed into the model
pts <- cbind(c(toyROMS$lon_u+0.1), c(toyROMS$lat_u+0.1), as.vector(as.numeric(0)))
## run 3D tracking
track <- trackit_3D(pts = pts, romsobject = toyROMS)
## checking results
plot(pts, cex=0.5)
points(track$pnow, col = "red", cex=0.5)
##
track2 <- trackit_2D(pts = track$pnow, romsobject = toyROMS)

#### testrun
## setup ROMS field
nc_4k <- nc_open(paste0(data.dat100,"ocean_avg_0001.nc"))
testROMS <- list()
testROMS$lon_u <- lon_rho
testROMS$lat_u <- lat_rho
testROMS$h <- ncvar_get(nc_4k, varid="h")
#testROMS$hh <- ncvar_get(nc_4k, varid="h")
testROMS$u <- ncvar_get(nc_4k, varid="u") 
testROMS$v <- ncvar_get(nc_4k, varid="v") 
testROMS$w <- ncvar_get(nc_4k, varid="w") 





#### create pattern of points to seed into the model
## ptrackr function doesn't work yet  
#pts <- create_points_pattern(npp, multi=100)
## simplest option is to seed one per cell
pts.all <- cbind(c(lon_rho), c(lat_rho), as.vector(as.numeric(0)))
pts <- pts.all[1:100,]

## run 3D tracking
track <- trackit_3D(pts = pts, romsobject = testROMS)


## checking the results
plot(pts, cex=0.5)
points(track$pnow, col = "red", cex=0.5)







# source(paste0(ptrackr_dir,"create_points_pattern.R"))
# source(paste0(ptrackr_dir,"create_points_pattern.R"))
# source(paste0(ptrackr_dir,"create_points_pattern.R"))
# source(paste0(ptrackr_dir,"create_points_pattern.R"))





# library(rgl)
# plot3d(pts, zlim = c(-1500,1))
# plot3d(track$pnow, col = "red", add = TRUE)

## better:
library(rasterVis)
library(rgdal)

ra <- raster(nrow = 50, ncol = 50, ext = extent(surface_chl))
r_roms <- rasterize(x = cbind(as.vector(toyROMS$lon_u), as.vector(toyROMS$lat_u)), y = ra, field = as.vector(-toyROMS$h))
pr <- projectRaster(r_roms, crs = "+proj=laea +lon_0=137 +lat_0=-66")  #get the right projection (through the centre)

plot3D(pr, adjust = FALSE, zfac = 50)                    # plot bathymetry with 50x exaggerated depth
points <- matrix(NA, ncol=3, nrow=dim(track$ptrack)[1])  # get Tracking-points
for(i in seq_len(dim(track$ptrack)[1])){
  points[i,] <- track$ptrack[i,,track$stopindex[i]] 
}
pointsxy <- project(as.matrix(points[,1:2]), projection(pr))  #projection on Tracking-points
points3d(pointsxy[,1], pointsxy[,2], points[,3]*50)

ptsxy <- project(as.matrix(pts[,1:2]), projection(pr))  #projection on Tracking-points
points3d(ptsxy[,1], ptsxy[,2], pts[,3]*50, col = "red")



##################################
## assign NPP values to each particle



#### 2) Source AAD datafiles----
#Code to download data from other remote repositories to your local machine
#https://ropensci.org/blog/2018/11/13/antarctic/

#This is a very large file
## seaice data needs login...see separate code file for sea-ice
# src <- bind_rows(sources("NSIDC SMMR-SSM/I Nasateam sea ice concentration", hemisphere = "south", time_resolutions = "day", years = c(2013)),
#                  sources("Southern Ocean summer chlorophyll-a climatology (Johnson)"),
#                  sources("IBCSO bathymetry"),
#                  bb_modify_source(sources("Oceandata MODIS Aqua Level-3 mapped daily 4km chl-a"), method = list(search = "A2014*L3m_DAY_CHL_chlor_a_4km.nc")))
src <- bind_rows(sources("Southern Ocean summer chlorophyll-a climatology (Johnson)"),
                 sources("IBCSO bathymetry"),
                 )
result <- bb_get(src, local_file_root = AAD_dir, clobber = 0, verbose = TRUE, confirm = NULL)                                                                 


#But usually we only need to set the directory where the AAD functions look for the data
set_data_roots(AAD_dir)


#### 3) Extract bathymetry and generate TPI----
# choose which data-layer to use, which resolution and projection.

## Choose dataset
#data.name <- "ibcso"  ## ibcso can be found using Ben's AAd dataset-wrapper
data.name <- "ibcso2"  ## DOESN'T EXIST YET, DOESN'T FUNCTION IN THE CODE YET
#data.name <- "gebco"  ## gebco is on file

## Choose resolution and projection (only relevant for gebco bathymetry)
data.format <- "project.500" ## project polar stereographic to 500m resolution
#data.format <- "project.raw" ## project polar stereographic but keep resolution as it is
#data.format <- "unprojected" ## keep raw format

## Choose range of depth
data.depth <- "shelf"
#data.depth <- "full"

if(data.depth=="shelf"){
  depth.range <- c(0,-3000)
}else depth.range <- c(0,-Inf)


#NOTE: gebco, takes 30-40mins to run!
## read in bathymetry data and project/reproject if needed
ri <- readtopo("ibcso") ## already projected and at 500m resolution
if(data.name=="gebco"){
  ## load raw (unprojected) data
  # g <- raster(paste0(env.raw,"GEBCO_2020/gebco_2020_SO.tif"))
  # g <- raster(paste0(env.raw,"gebco_2021_sub_ice_topo/GEBCO_2021_sub_ice_topo_CF.nc"))
  # crs(g) <- "+proj=longlat +datum=WGS84 +no_defs"
  # g <- crop(g,extent(c(-180,180,-90,-60)))
  g <- raster(paste0(env.raw,"GEBCO_2021/gebco_2021_SO.tif"))
  if(data.format=="project.500"|data.format=="project.raw"){
    ## project raster to polar stereographic
    g <- projectRaster(g,crs=stereo)
  }
  if(data.format=="project.500"){
    ## change resolution to 500m grid cells
    g <- projectRaster(g,r)
  }
  r <- g  
}

## OR, load the latest IBCSO version here:
r.raw <- raster(paste0(env.raw,"IBCSO_2022/IBCSO_v2_ice-surface.tif"))
r.raw.c <- crop(r.raw, ri)
r <- projectRaster(r.raw.c, ri)

## create layers for depth, slope and topographic position index (TPI) at different scales
r.depth <- r
#r.depth[r.depth>=200] <- NA
#r.depth[r.depth<=-4000] <- NA

r.slope <- terrain(r.depth)
r.tpi <- tpi(r.depth)
r.tpi5 <- tpi(r.depth, scale=5)
r.tpi11 <- tpi(r.depth, scale=11)
# r.tpi21 <- tpi(r, scale=21)
# r.tpi31 <- tpi(r, scale=31)

## set anything on land to NA, and optionally set abyssal zone to NA
r.depth[r>0] <- NA
r.depth[r<=depth.range[2]] <- NA
r.slope[is.na(r.depth[])] <- NA
r.tpi[is.na(r.depth[])] <- NA
r.tpi5[is.na(r.depth[])] <- NA
r.tpi11[is.na(r.depth[])] <- NA
# r.tpi21[is.na(r.depth[])] <- NA
# r.tpi31[is.na(r.depth[])] <- NA


###

r <- raster(paste0(env.derived,string.chr,"500m_bathy_gebco.grd"))
r[r>=0] <- NA
r[r<=depth.range[2]] <- NA

r.depth <- raster(paste0(env.derived,string.chr,"500m_shelf_bathy_ibcso2_depth.tif"))







