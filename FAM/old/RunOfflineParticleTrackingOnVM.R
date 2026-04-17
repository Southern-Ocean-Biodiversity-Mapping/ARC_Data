#### 1) Set up ----
# load libraries
library(ncdf4)        ## package for netcdf manipulation
library(terra)
library(reproj) ## reproject coordinates
library(dplyr)
library(reshape) ## prep data for plotting
library(ggplot2) ## plotting

# set up directory pointers etc
env.dir <- "/pvol/data_environmental/"

string.chr <- "Circumpolar_EnvData_"
string.res <- "500m_"

# polar stereographic projection:
stereo <- "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

################################
#### ROMS Currents
## WE NEED A RUN ACROSS SUMMER WITH HIGH-RES HISTORY FILES!!! 8 has one month only, 8_Nov has 10-day intervals...
## load lon/lat information from ROMS-grid
grd4k_nc <- nc_open(paste0(env.dir,"waom4extend_grd.nc"))
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
h <- rast(paste0(env.dir,"ocean_avg_0001.nc"), subds="h")

###################################
## load NPP
# npp.raw <- rast(paste0(env.derived,string.chr,"NPP_Cafe_filled_SummerAverage.tif"))
# npp <- project(npp.raw, uv.mean)

###################################
###### test particle tracking
#### load functions
ptrackr_dir <- "/pvol/Scripts/ptrackR_package_update2024/ptrackr/R/"
## minimum for 3d run
source(paste0(ptrackr_dir,"create_points_pattern.R"))
source(paste0(ptrackr_dir,"setup_knn.R"))
source(paste0(ptrackr_dir,"trackit_3D.R"))
## for 2d run:
source(paste0(ptrackr_dir,"trackit_2D.R"))
source(paste0(ptrackr_dir,"buildparams.r"))

#### toyrun:
## setup ROMS field
load("/pvol/Scripts/ptrackR_package_update2024/ptrackr/data/toyROMS.rdata")
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
nc_4k <- nc_open(paste0(env.dir,"ocean_avg_0001.nc"))
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
#pts.sel <- sample(not.na, 1000000)
pts.sel <- not.na
pts <- pts.all[pts.sel,]

## run 3D tracking
start.time <- Sys.time()
track <- trackit_3D(pts = pts, romsobject = testROMS, time=15)
end.time <- Sys.time()
end.time-start.time

save(track, file="/pvol/3_model_analysis/testrun_500kparticles15days.Rdata")

## the below times are trackit_3D only
## 1min for 1k particles and 1 day
## 1.8min for 10k particles and 1 day
## 5.8min for 100k particles and 1 day (printed steps starting at 1s per step, then slowing down to 4s per step by day 1)
## 7.3hrs for 550k particles and 15 days (= 12GB file output)

## checking the results
sel.pl <- sample(1:nrow(pts), 1000)
pts.sub <- pts[sel.pl,]
pnow.sub <- track$pnow[sel.pl,]
plot(pts.sub[,1:2], cex=0.5)
points(pnow.sub[,1:2], col = "red", cex=0.5)

## checking the results
sel.pl <- which(pts[,1]>1050000&pts[,1]<1120000&pts[,2]<(-2080000)&pts[,2]>(-2100000))
pts.sub <- pts[sel.pl,]
pnow.sub <- track$pnow[sel.pl,]
plot(pts.sub[,1:2], cex=0.5)
points(pnow.sub[,1:2], col = "red", cex=0.5)


plot(pts, cex=0.5)
points(track$pnow, col = "red", cex=0.5)



## NOTES:
## - STILL AN UGLY VERSION, BECAUSE I DELETE ROWS & COLUMS FROM ROMS DATA TO MAKE THEM COMPARABLE
## - DELETE MORE POINTS FROM AREAS FAR OFF THE SHELF
## - DELETE POINTS FROM UNDER ICE-SHELVES
