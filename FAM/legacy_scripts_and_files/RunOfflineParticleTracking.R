#### 1) Set up ----
# load libraries
library(ncdf4)        ## package for netcdf manipulation
library(terra)
library(reproj) ## reproject coordinates
library(dplyr)
library(reshape) ## prep data for plotting
library(ggplot2) ## plotting

# library(raadtools)
# library(sp)
# library(blueant)
# library(ggplot2)      ## package for plotting
# library(spatialEco)

## old packages to replace
#library(rgdal)        ## package for geospatial analysis
#library(raster)       ## package for raster manipulation

# set up directory pointers etc
# Jan's local machine:
sci.dir <- "C:/Users/jjansen/OneDrive - University of Tasmania/science/"
env.dir <- paste0(sci.dir,"data_environmental/")
ARC_dir <- paste0(sci.dir,"SouthernOceanBiodiversityMapping/ARC_Data/")

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
nc_close(grd4k_nc)
## Prepare empty rasters to assign correct projected values to
#roms.coords.proj <- rgdal::project(cbind(c(lon_rho), c(lat_rho)), proj=stereo)
roms.coords.proj <- reproj(cbind(c(lon_rho), c(lat_rho)), target=stereo)
roms.coords.proj.lon <- matrix(roms.coords.proj[,1], ncol=1400)
roms.coords.proj.lat <- matrix(roms.coords.proj[,2], ncol=1400)
x.range <- c(min(roms.coords.proj[,1])-2000,max(roms.coords.proj[,1])+2000)
y.range <- c(min(roms.coords.proj[,2])-2000,max(roms.coords.proj[,2])+2000)
#empty.roms.ra <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=4000)

## depth
h <- rast(paste0(data.dat100,"ocean_avg_0001.nc"), subds="h")
## seafloor variables
s <- seq(1,217, by=31)

# #seafloor currents (seafloor-layer is 1)
# u.raw1 <- brick(paste0(data.dat100,"ocean_his_0002.nc"), varname="u", level=1)
# v.raw1 <- brick(paste0(data.dat100,"ocean_his_0002.nc"), varname="v", level=1)
# u.raw2 <- brick(paste0(data.dat100,"ocean_his_0003.nc"), varname="u", level=1)
# v.raw2 <- brick(paste0(data.dat100,"ocean_his_0003.nc"), varname="v", level=1)
# u.raw3 <- brick(paste0(data.dat100,"ocean_his_0004.nc"), varname="u", level=1)
# v.raw3 <- brick(paste0(data.dat100,"ocean_his_0004.nc"), varname="v", level=1)
# ## sum up monthly values for a climatology (THIS SHOULD BE DONE ON HIGH RESOLUTION HISTORY FILES)
# u.sum <- sum(u.raw)
# v.sum <- sum(v.raw)
# u.sum.abs <- sum(abs(u.raw))
# v.sum.abs <- sum(abs(v.raw))
# u31.sum <- sum(u_31.raw)
# v31.sum <- sum(v_31.raw)

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
# u <- v <- u31 <- v31 <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=4000, nlyr=nlyr(u.raw))
# ## need to do this consecutively because of RAM
# u.dat <- extract(u.raw, coord.grd.u, method="bilinear")
# for(i in 1:ncol(u.dat)){
#   print(i)
#   u[[i]][] <- u.dat[,i]
# }
# rm(u.dat, u.raw)
# u31.dat <- extract(u_31.raw, coord.grd.u, method="bilinear")
# for(i in 1:ncol(u31.dat)){
#   print(i)
#   u31[[i]][] <- u31.dat[,i]
# }
# rm(u31.dat, u_31.raw, coord.grd.u)
# v.dat <- extract(v.raw, coord.grd.v, method="bilinear")
# for(i in 1:ncol(v.dat)){
#   print(i)
#   v[[i]][] <- v.dat[,i]
# }
# rm(v.dat, v.raw)
# v31.dat <- extract(v_31.raw, coord.grd.v, method="bilinear")
# for(i in 1:ncol(v31.dat)){
#   print(i)
#   v31[[i]][] <- v31.dat[,i]
# }
# rm(v31.dat, v_31.raw, coord.grd.v)
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
# writeRaster(u.mean,filename=paste0(env.derived,string.chr,"4km_waom4k_seafloorcurrents_mean_u.tif"), overwrite=TRUE)
# writeRaster(v.mean,filename=paste0(env.derived,string.chr,"4km_waom4k_seafloorcurrents_mean_v.tif"), overwrite=TRUE)
# writeRaster(uv.mean,filename=paste0(env.derived,string.chr,"4km_waom4k_seafloorcurrents_mean.tif"), overwrite=TRUE)
# writeRaster(uv.abs.mean, filename=paste0(env.derived,string.chr,"4km_waom4k_seafloorcurrents_mean_absolute.tif"), overwrite=TRUE)
# writeRaster(uv.max, filename=paste0(env.derived,string.chr,"4km_waom4k_seafloorcurrents_max.tif"), overwrite=TRUE)
# writeRaster(uv.sd,  filename=paste0(env.derived,string.chr,"4km_waom4k_seafloorcurrents_sd.tif"), overwrite=TRUE)

## load current speeds
u.mean <- rast(paste0(env.derived,string.chr,"4km_waom4k_seafloorcurrents_mean_u.tif"))
v.mean <- rast(paste0(env.derived,string.chr,"4km_waom4k_seafloorcurrents_mean_v.tif"))
uv.mean <- rast(paste0(env.derived,string.chr,"4km_waom4k_seafloorcurrents_mean.tif"))
uv.abs.mean <- rast(paste0(env.derived,string.chr,"4km_waom4k_seafloorcurrents_mean_absolute.tif"))
uv.max <- rast(paste0(env.derived,string.chr,"4km_waom4k_seafloorcurrents_max.tif"))
uv.sd <- rast(paste0(env.derived,string.chr,"4km_waom4k_seafloorcurrents_sd.tif"))

###################################
## load NPP
npp.raw <- rast(paste0(env.derived,string.chr,"NPP_Cafe_filled_SummerAverage.tif"))
npp <- project(npp.raw, uv.mean)

###################################
## load depth field
#### run python code:
# library(reticulate)
# source_python(paste0(ARC_dir,"FAM/RunOfflineParticleTracking_PythonPart.py"))
# 
# z_rho

###################################
###### test particle tracking
#### load functions
ptrackr_dir <- paste0(sci.dir,"Scripts/ptrackR_package_update2024/ptrackr/R/")
## minimum for 3d run
source(paste0(ptrackr_dir,"create_points_pattern.R"))
source(paste0(ptrackr_dir,"setup_knn.R"))
source(paste0(ptrackr_dir,"trackit_3D.R"))
## for 2d run:
source(paste0(ptrackr_dir,"trackit_2D.R"))
source(paste0(ptrackr_dir,"buildparams.R"))

#### toyrun:
## setup ROMS field
load(paste0(sci.dir,"Scripts/ptrackR_package_update2024/ptrackr/data/toyROMS.rdata"))
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
# testROMS$lon_u <- lon_rho
# testROMS$lat_u <- lat_rho
testROMS$lon_u <- roms.coords.proj.lon[-1575,-1400]
testROMS$lat_u <- roms.coords.proj.lat[-1575,-1400]
testROMS$h <- ncvar_get(nc_4k, varid="h")[-1575,-1400]
#testROMS$hh <- ncvar_get(nc_4k, varid="h")
testROMS$i_u <- ncvar_get(nc_4k, varid="u")[,-1400,,]
testROMS$i_v <- ncvar_get(nc_4k, varid="v")[-1575,,,]
testROMS$i_w <- ncvar_get(nc_4k, varid="w")[-1575,-1400,-32,]

s_rho <- ncvar_get(nc_4k, varid="s_rho") ## 31 vertical levels varying between 0 to 1
Cs_r <- ncvar_get(nc_4k, varid="Cs_r")   ## 31 "S-coordinate stretching curves at RHO-points
hc <- ncvar_get(nc_4k, varid="hc")       ## value of 250
h <- ncvar_get(nc_4k, varid="h")         ## 1575 1400 model bathymetry
zeta <- ncvar_get(nc_4k, varid="zeta")   ## 1575 1400 7 free surface elevation
zice <- ncvar_get(nc_4k, varid="zice")   ## 1575 1400 ice draft

nc_close(nc_4k)

library(abind)
## simplest calculation:
hh <- array(data=h*(Cs_r[1]), dim=c(1575,1400,1))
for(i in 2:31){
  hh <- abind(hh, h*(Cs_r[i]))
}
testROMS$hh <- hh[-1575,-1400,]

rm(s_rho, Cs_r, hc, zeta, grd4k_nc, nc_4k)

# ## "proper" calculation, but gives wrong values?
# ## calculate the parameters for the terrain-following vertical coordinates
# z0_rho <- array(data=((hc*s_rho)[1]+h*(Cs_r[1]))/(hc+h), dim=c(1575,1400,1))
# for(i in 2:31){
#   z0_rho <- abind(z0_rho, ((hc*s_rho)[i]+h*(Cs_r[i]))/(hc+h))
# }
# ## z0_rho has 1575 1400 31
# ##
# z_rho <- array(data=zeta[,,1]+(zeta[,,1]+h)*z0_rho[,,1]*zice, dim=c(1575,1400,1))
# for(i in 2:31){
#   z_rho <- abind(z_rho, zeta[,,1]+(zeta[,,1]+h)*z0_rho[,,i]*zice)
# }

#### create pattern of points to seed into the model
## which points are not on land and not deeper than 3000m?
hz <- h[-1575,-1400]+zice[-1575,-1400] ## need to combine bathymetry with ice draft
hz[which(is.na(testROMS$i_w[,,1,1]))] <- NA
hz[hz>4000] <- NA
hz[1:1574,1:150] <- NA
hz[1500:1574,1:1399] <- NA
hz[1:600,1:300] <- NA
h.df <- melt(hz, varnames = c("row", "col"), value.name = "value")
ggplot(h.df, aes(x = row, y = col, fill = value)) + geom_tile()

## simplest option is to seed one per cell
not.na <- which(!is.na(hz))
#pts.all <- cbind(c(lon_rho), c(lat_rho), as.vector(as.numeric(0)))
pts.all <- cbind(c(roms.coords.proj.lon[-1575,-1400]), c(roms.coords.proj.lat[-1575,-1400]), as.vector(as.numeric(0)))
pts.sel <- sample(not.na, 200000)
pts <- pts.all[pts.sel,]

## run 3D tracking
start.time <- Sys.time()
track <- trackit_3D(pts = pts, romsobject = testROMS, time=1)
end.time <- Sys.time()
end.time-start.time

## the below times are trackit_3D only
## 51s for 1k particles and 1 day
## 48s for 10k particles and 1 day
## 3min for 100k particles and 1 day (printed steps starting at 1s per step, then slowing down to 4s per step by day 1)
## 3.5min for 200k particles and 1 day (printed steps starting at 3s per step, then slowing down to 6s per step by day 1)

## checking the results
plot(pts, cex=0.5)
points(track$pnow, col = "red", cex=0.5)



## NOTES:
## - STILL AN UGLY VERSION, BECAUSE I DELETE ROWS & COLUMS FROM ROMS DATA TO MAKE THEM COMPARABLE
## - DELETE MORE POINTS FROM AREAS FAR OFF THE SHELF
## - DELETE POINTS FROM UNDER ICE-SHELVES










###################################
## fix NA issue (may be kdtree related)

## TRACKIT_3D

pts = pts
romsobject = testROMS2
w_sink=100
time=0.5
romsparams=NULL
loop_trackit=FALSE
time_steps_in_s=1800

## We need an id for each particle to follow individual tracks
id_vec <- seq_len(nrow(pts))
  
## build a kdtree
if(loop_trackit==TRUE){
  kdtree <- romsobject$kdtree
  kdxy <- romsobject$kdxy
}else{
  sknn <- with(romsobject, setup_knn(lon_u, lat_u, hh))             # (lon_roms=lon_u, lat_roms=lat_u, depth_roms=hh)
  kdtree <- sknn$kdtree
  kdxy <- sknn$kdxy
}
  
## assign current speeds and depth for each ROMS-cell (lat/lon position), and boundaries of the region
if(!is.null(romsparams)){
  i_u <- romsparams$i_u
  i_v <- romsparams$i_v
  i_w <- romsparams$i_w
  h <- romsparams$h
  #roms_ext <- romsparams$roms_ext
}else{
  i_u <- romsobject$i_u[,,1]
  i_v <- romsobject$i_v[,,1]
  i_w <- romsobject$i_w[,,1]
  h <- romsobject$h
  #roms_ext <- c(min(romsobject$lon_u), max(romsobject$lon_u), min(romsobject$lat_u), max(romsobject$lat_u))
}
  
## w_sink is m/days, time is days
w_sink <- -w_sink/(60*60*24)                               ## sinking speed transformation into m/sec
ntime <- time*24*2                                         ## days transformation into 0.5h-intervals

## empty objects for the loop
ptrack <- array(0, c(length(as.vector(pts))/3, 3, ntime))  ## create an empty array to store particle-tracks
stopped <- rep(FALSE, length(as.vector(pts))/3)            ## create a stopping-vector
stopindex <- rep(0, length(as.vector(pts))/3)              ## a vector to store indices of when particles stopped
indices <- vector("list", ntime)                           ## a list of indices to store which 3D-cell a particle is in
indices_2D <- vector("list", ntime)                        ## a list of indices to store which 2D-cell a particle is in                       

pnow <- plast <- pts                ## copies of the starting points for updating in the loop

for (itime in seq_len(ntime)) {
  
  if(loop_trackit==FALSE){
    if(itime==1) message(paste0("starting # of particles: ",dim(pts)[1]))
    print(itime)
  }
  
  ## index 1st nearest neighbour of trace points to grid points
  dmap <- kdtree$query(plast, k = 1, eps = 0, radius=0)           ## one kdtree
  head(dmap$nn.idx)
  head(dmap$nn.dists)
  
  ## and to 2D space
  two_dim_pos <- kdxy$query(plast[,1:2], k = 1, eps = 0, radius=0)
  head(two_dim_pos$nn.idx)
  head(two_dim_pos$nn.dists)
  
  ## store indices for tracing particle positions
  indices[[itime]] <- dmap$nn.idx
  indices_2D[[itime]] <- two_dim_pos$nn.idx
  
  ## different to 2D in this line:
  idx_for_roms <- dmap$nn.idx
  
  ## extract component values from the vars
  thisu <- i_u[idx_for_roms]                             ## u-component of ROMS
  thisv <- i_v[idx_for_roms]                             ## v-component of ROMS
  thisw <- i_w[idx_for_roms]                             ## w-component of ROMS
  
  cbind(thisu[1:20], thisv[1:20], thisw[1:20])
  
  ## update this time step longitude, latitude, depth
  pnow[,1] <- plast[,1] + (thisu * time_steps_in_s) / (1.852 * 60 * 1000 * cos(plast[,2] * pi/180))
  pnow[,2] <- plast[,2] + (thisv * time_steps_in_s) / (1.852 * 60 * 1000)
  pnow[,3] <- pmin(0, plast[,3])  + ((thisw + w_sink)* time_steps_in_s )
  
  ## different to 2D here
  #     stopped <- (pnow[,1] < roms_ext[1] | pnow[,1] > roms_ext[2] |
  #                   pnow[,2] < roms_ext[3] | pnow[,2] > roms_ext[4]
  #                 )
  #     
  ##########---- only in trackit_3D:
  ## stopping conditions (hit the bottom)
  stopped <- pnow[,3] <= -h[two_dim_pos$nn.idx]
  stopindex[stopindex == 0 & stopped] <- itime
  ##########----
  
  ## assign stopping location of points to ptrack
  ptrack[,,itime] <- pnow
  plast <- pnow
  if (all(stopped)) {
    message("exiting, all stopped")
    break;
  }
}
ptrack <- ptrack[,,seq(itime),drop=FALSE]
list(ptrack = ptrack, pnow = pnow, plast = plast, stopindex = stopindex, indices = indices, indices_2D = indices_2D)
























