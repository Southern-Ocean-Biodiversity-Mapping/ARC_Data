#### 1) Set up ----
# load libraries
library(ncdf4)        ## package for netcdf manipulation
library(terra)
library(reproj) ## reproject coordinates
library(dplyr)
library(reshape) ## prep data for plotting
library(ggplot2) ## plotting
library(abind) ## combine arrays

# set up directory pointers etc
env.dir <- "/pvol/data_environmental/ROMS_4k_files/"

string.chr <- "Circumpolar_EnvData_"
string.res <- "500m_"

# polar stereographic projection:
stereo <- "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

################################
#### ROMS Currents
## WE NEED A RUN ACROSS SUMMER WITH HIGH-RES HISTORY FILES!!! 8 has one month only, 8_Nov has 10-day intervals...
## load lon/lat information from ROMS-grid
grd4k_nc <- nc_open(paste0(env.dir,"waom4extend_grd.nc"))
lon_rho <- ncvar_get(grd4k_nc, varid="lon_rho")[-1575,-1400]
lat_rho <- ncvar_get(grd4k_nc, varid="lat_rho")[-1575,-1400]
nc_close(grd4k_nc)
## Prepare empty rasters to assign correct projected values to
#roms.coords.proj <- rgdal::project(cbind(c(lon_rho), c(lat_rho)), proj=stereo)
roms.coords.proj <- reproj(cbind(c(lon_rho), c(lat_rho)), target=stereo)
roms.coords.proj.lon <- matrix(roms.coords.proj[,1], ncol=1399)
roms.coords.proj.lat <- matrix(roms.coords.proj[,2], ncol=1399)
x.range <- c(min(roms.coords.proj[,1])-2000,max(roms.coords.proj[,1])+2000)
y.range <- c(min(roms.coords.proj[,2])-2000,max(roms.coords.proj[,2])+2000)
empty.roms.ra <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=4000)
uv.ra <- empty.roms.ra
uv.surf.ra <- empty.roms.ra
h.ra <- empty.roms.ra

## depth
#h.ra <- rast(paste0(env.dir,"ocean_avg_0001.nc"), subds="h")

###################################
## load NPP
# npp.raw <- rast(paste0(env.derived,string.chr,"NPP_Cafe_filled_SummerAverage.tif"))
# npp <- project(npp.raw, uv.mean)

###################################
#### setup ROMS field
nc_4k <- nc_open(paste0(env.dir,"ocean_avg_0001.nc"))
s_rho <- ncvar_get(nc_4k, varid="s_rho") ## 31 vertical levels varying between 0 to 1
Cs_r <- ncvar_get(nc_4k, varid="Cs_r")   ## 31 "S-coordinate stretching curves at RHO-points
hc <- ncvar_get(nc_4k, varid="hc")       ## value of 250
h <- ncvar_get(nc_4k, varid="h")         ## 1575 1400 model bathymetry
zeta <- ncvar_get(nc_4k, varid="zeta")   ## 1575 1400 7 free surface elevation
zice <- ncvar_get(nc_4k, varid="zice")   ## 1575 1400 ice draft

testROMS <- list()
# testROMS$lon_u <- lon_rho
# testROMS$lat_u <- lat_rho
#testROMS$hh <- ncvar_get(nc_4k, varid="h")
testROMS$x <- roms.coords.proj.lon[-1575,-1400]
testROMS$y <- roms.coords.proj.lat[-1575,-1400]
testROMS$h <- h[-1575,-1400]
testROMS$i_u <- ncvar_get(nc_4k, varid="u")[,-1400,,]
testROMS$i_v <- ncvar_get(nc_4k, varid="v")[-1575,,,]
testROMS$i_w <- ncvar_get(nc_4k, varid="w")[-1575,-1400,-32,]
nc_close(nc_4k)

### create uv and h raster for plotting
uv.ra[] <- c(sqrt(testROMS$i_u[,1399:1,1,1]^2+c(testROMS$i_v[,1399:1,1,1])^2))
uv.surf.ra[] <- c(sqrt(testROMS$i_u[,1399:1,31,1]^2+c(testROMS$i_v[,1399:1,31,1])^2))
h.ra[] <- c(testROMS$h[,1399:1])

### calculate depths of each vertical ROMS cell
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

hz <- h[-1575,-1400]+zice[-1575,-1400] ## need to combine bathymetry with ice draft
hz[which(is.na(testROMS$i_w[,,1,1]))] <- NA
hz[hz>4000] <- NA
hz[1:1574,1:150] <- NA
hz[1500:1574,1:1399] <- NA
hz[1:600,1:300] <- NA
h.df <- melt(hz, varnames = c("row", "col"), value.name = "value")
#ggplot(h.df, aes(x = row, y = col, fill = value)) + geom_tile()

###################################
#### visually check model setup
## 3D current field:
plot(uv.ra)
## 3D slice through Ross Sea (x=c(-450000,-300000), y=-1460000)
plot(uv.ra)


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
## to loop:
source(paste0(ptrackr_dir,"loopit_2D3D.R"))

###################################

#### testrun

#### create pattern of points to seed into the model
## which points are not on land and not deeper than 3000m?
## CUT OFF BY DISTANCE TO COAST???

## simplest option is to seed one per cell
# not.na <- which(!is.na(hz))
# #pts.all <- cbind(c(lon_rho), c(lat_rho), as.vector(as.numeric(0)))
# pts.all <- cbind(c(roms.coords.proj.lon[-1575,-1400]), c(roms.coords.proj.lat[-1575,-1400]), as.vector(as.numeric(0)))
# pts.sel <- sample(not.na, 100)
# #pts.sel <- not.na
# pts <- pts.all[pts.sel,]

## pts seeded only near the shelf break in the Ross Sea
not.na <- which(!is.na(hz))
pts.all <- cbind(c(roms.coords.proj.lon[-1575,-1400]), c(roms.coords.proj.lat[-1575,-1400]), as.vector(as.numeric(0)))
pts.all2 <- pts.all[not.na,]
RS.sel <- which(pts.all2[,1]>-500000 & pts.all2[,1]<0 & pts.all2[,2]<=-1400000 & pts.all2[,2]>=-1700000)
set.seed(1)
pts <- pts.all2[sample(RS.sel,1000),]

### check distribution of points
# plot(uv.ra)
# points(pts)
plot(uv.ra, xlim=c(-550000,50000), ylim=c(-1750000,-1350000))
points(pts)

####
### run 3D tracking
## the below times are trackit_3D only
## 1min for 1k particles and 1 day
## 1.8min for 10k particles and 1 day
## 5.8min for 100k particles and 1 day (printed steps starting at 1s per step, then slowing down to 4s per step by day 1)
## 7.3hrs for 550k particles and 15 days (= 12GB file output)
start.time <- Sys.time()
track <- trackit_3D(pts = pts, romsobject = testROMS, time=30, w_sink=100, projected=TRUE)
end.time <- Sys.time()
end.time-start.time


### checking 3D results
plot(uv.ra, xlim=c(-550000,50000), ylim=c(-1750000,-1350000))
## seeding points
points(pts, cex=0.5, pch=16)
## all points that haven't stopped have the pnow coordinates
pstop <- matrix(NA, nrow = 1000, ncol = 3)
not.stopped <- which(track$stopindex==0)
pstop[not.stopped,] <- track$pnow[not.stopped,]
## where the code says they stopped:
for (i in c(1:1000)[-not.stopped]) {
  pstop[i, ] <- track$ptrack[i, 1:3, track$stopindex[i]]
}
points(pstop, col="blue", cex=0.5)

### tracking some individuals
ind <- sample(1:1000,30)#c(31,87,425,613,957)
## all points
plot(uv.ra, xlim=c(-550000,50000), ylim=c(-1750000,-1350000))
# plot(pts, cex=0.5)
points(pts[ind,1:2], col="red")
## the full track of those points
for(i in 1:length(ind)){
  points(t(track$ptrack[ind[i],1:2,]), col="blue")
}
## only the track until they stopped
for(i in 1:length(ind)){
  points(t(track$ptrack[ind[i],1:2,1:track$stopindex[ind[i]]]), col="red")
}
## the end locations
points(track$pnow[ind,1:2], col = "green", pch=16)
points(pstop[ind,1:2], col="orange", pch=16)

####
### 2D run (seeding with the correct locations)
track2 <- trackit_2D(pts = pstop, romsobject = testROMS, projected=TRUE, sedimentation=TRUE)

### checking 2D results
plot(uv.ra, xlim=c(-550000,50000), ylim=c(-1750000,-1350000))
#plot(pts, cex=0.5)
points(pts[ind,1:2], col="red")
points(pstop[ind,1:2], col="orange", pch=16)
## all points that haven't stopped have the pnow coordinates
pstop2 <- matrix(NA, nrow = 1000, ncol = 3)
not.stopped <- which(track2$stopindex==0)
pstop2[not.stopped,] <- track2$pnow[not.stopped,]
## where the code says they stopped:
for (i in c(1:1000)) {
  pstop2[i, ] <- track2$ptrack[i, 1:3, track2$stopindex[i]]
}
points(pstop2, col="blue", cex=0.7)

### tracking some individuals
## only the track until they stopped
for(i in 1:length(ind)){
  points(t(track$ptrack[ind[i],1:2,1:track$stopindex[ind[i]]]), col="red", cex=0.1)
}
## the end locations
points(pstop[ind,1:2], col="orange", pch=16)
## the full track of those points
for(i in 1:length(ind)){
  points(t(track2$ptrack[ind[i],1:2,]), col="purple", cex=0.1)
}
## only the track until they stopped
for(i in 1:length(ind)){
  points(t(track2$ptrack[ind[i],1:2,1:track2$stopindex[ind[i]]]), col="gold", pch=16, cex=0.1)
}

# ############################
# ### tracking all individuals
# plot(uv.surf.ra, xlim=c(-550000,50000), ylim=c(-1750000,-1350000))
# ## only the track until they stopped
# for(i in 1:1000){
#   points(t(track$ptrack[i,1:2,1:track$stopindex[i]]), col="red", cex=0.1)
# }
# ## the end locations
# points(pstop[ind,1:2], col="orange", pch=16)
# 
# plot(uv.ra, xlim=c(-550000,50000), ylim=c(-1750000,-1350000))
# ## the full track of those points
# for(i in 1:1000){
#   points(t(track2$ptrack[i,1:2,]), col="purple", cex=0.1)
# }
# ## only the track until they stopped
# for(i in 1:1000){
#   points(t(track2$ptrack[i,1:2,1:track2$stopindex[i]]), col="gold", cex=0.1)
# }




###################################

#### SECOND testrun: smaller region, more particles seeded only near the shelf break in the Ross Sea

xlim=c(-150000,-50000)
ylim=c(-1700000,-1600000)
p.xlim=c(-200000,0)
p.ylim=c(-1750000,-1550000)

RS.sel <- which(pts.all2[,1]>xlim[1] & pts.all2[,1]<xlim[2] & pts.all2[,2]>=ylim[1] & pts.all2[,2]<=ylim[2])
set.seed(1)
pts <- pts.all2[RS.sel,]

### check distribution of points
plot(uv.ra, xlim=p.xlim, ylim=p.ylim)
points(pts, cex=0.1)

####
### run 3D tracking
start.time <- Sys.time()
track <- trackit_3D(pts = pts, romsobject = testROMS, time=30, w_sink=100, projected=TRUE)
end.time <- Sys.time()
end.time-start.time


### checking 3D results
plot(uv.ra, xlim=p.xlim, ylim=p.ylim)
## seeding points
points(pts, cex=0.5, pch=16)
## all points that haven't stopped have the pnow coordinates
pstop <- matrix(NA, nrow = nrow(track$ptrack), ncol = 3)
# not.stopped <- which(track$stopindex==0)
# pstop[not.stopped,] <- track$pnow[not.stopped,]
## where the code says they stopped:
for (i in 1:nrow(track$ptrack)) {
  pstop[i, ] <- track$ptrack[i, 1:3, track$stopindex[i]]
}
points(pstop, col="blue", cex=0.5)

### tracking some individuals
ind <- sample(1:nrow(track$ptrack),30)#c(31,87,425,613,957)
## all points
plot(uv.ra, xlim=p.xlim, ylim=p.ylim)
# plot(pts, cex=0.5)
points(pts[ind,1:2], col="red")
## the full track of those points
for(i in 1:length(ind)){
  points(t(track$ptrack[ind[i],1:2,]), col="blue")
}
## only the track until they stopped
for(i in 1:length(ind)){
  points(t(track$ptrack[ind[i],1:2,1:track$stopindex[ind[i]]]), col="red")
}
## the end locations
points(track$pnow[ind,1:2], col = "green", pch=16)
points(pstop[ind,1:2], col="orange", pch=16)

####
### 2D run (seeding with the correct locations)
track2 <- trackit_2D(pts = pstop, romsobject = testROMS, projected=TRUE, sedimentation=TRUE)

### checking 2D results
plot(uv.ra, xlim=p.xlim, ylim=p.ylim)
#plot(pts, cex=0.5)
points(pts[ind,1:2], col="red")
points(pstop[ind,1:2], col="orange", pch=16)
## all points that haven't stopped have the pnow coordinates
pstop2 <- matrix(NA, nrow = nrow(track$ptrack), ncol = 3)
# not.stopped <- which(track2$stopindex==0)
# pstop2[not.stopped,] <- track2$pnow[not.stopped,]
## where the code says they stopped:
for (i in 1:nrow(track$ptrack)) {
  pstop2[i, ] <- track2$ptrack[i, 1:3, track2$stopindex[i]]
}
points(pstop2, col="blue", cex=0.7)

### tracking some individuals
## only the track until they stopped
for(i in 1:length(ind)){
  points(t(track$ptrack[ind[i],1:2,1:track$stopindex[ind[i]]]), col="red", cex=0.1)
}
## the end locations
points(pstop[ind,1:2], col="orange", pch=16)
## the full track of those points
for(i in 1:length(ind)){
  points(t(track2$ptrack[ind[i],1:2,]), col="purple", cex=0.1)
}
## only the track until they stopped
for(i in 1:length(ind)){
  points(t(track2$ptrack[ind[i],1:2,1:625]), col="gold", pch=16, cex=0.1)
}







###################################

#### THIRD testrun: small region, one particle per skm (so 16 per 4km cell)

xlim=c(-150000,-50000)
ylim=c(-1700000,-1600000)
p.xlim=c(-200000,0)
p.ylim=c(-1750000,-1550000)

## create points in that region
x.seq <- seq(xlim[1],xlim[2], by=1000)
y.seq <- seq(ylim[1],ylim[2], by=1000)
x.dat <- rep(x.seq, length(y.seq))
y.dat <- rep(y.seq[1], length(x.seq))
for(i in 1:length(y.seq)){
  y.dat <- c(y.dat, rep(y.seq[i], length(x.seq)))
}
pts <- cbind(x.dat, y.dat, 0)

### check distribution of points
plot(uv.ra, xlim=p.xlim, ylim=p.ylim)
points(pts, cex=0.1)

####
### run 3D tracking
start.time <- Sys.time()
track <- trackit_3D(pts = pts, romsobject = testROMS, time=30, w_sink=100, projected=TRUE)
end.time <- Sys.time()
end.time-start.time


### checking 3D results
plot(uv.ra, xlim=p.xlim, ylim=p.ylim)
## seeding points
points(pts, cex=0.5, pch=16)
## all points that haven't stopped have the pnow coordinates
pstop <- matrix(NA, nrow = nrow(track$ptrack), ncol = 3)
not.stopped <- which(track$stopindex==0)
pstop[not.stopped,] <- track$pnow[not.stopped,]
## where the code says they stopped:
for (i in c(1:nrow(track$ptrack))[-not.stopped]) {
  pstop[i, ] <- track$ptrack[i, 1:3, track$stopindex[i]]
}
points(pstop, col="blue", cex=0.5)

### tracking some individuals
ind <- sample(1:nrow(track$ptrack),30)#c(31,87,425,613,957)
## all points
plot(uv.ra, xlim=p.xlim, ylim=p.ylim)
# plot(pts, cex=0.5)
points(pts[ind,1:2], col="red")
## the full track of those points
for(i in 1:length(ind)){
  points(t(track$ptrack[ind[i],1:2,]), col="blue")
}
## only the track until they stopped
for(i in 1:length(ind)){
  points(t(track$ptrack[ind[i],1:2,1:track$stopindex[ind[i]]]), col="red")
}
## the end locations
points(track$pnow[ind,1:2], col = "green", pch=16)
points(pstop[ind,1:2], col="orange", pch=16)

####
### 2D run (seeding with the correct locations)
track2 <- trackit_2D(pts = pstop, romsobject = testROMS, projected=TRUE, sedimentation=TRUE)

### checking 2D results
plot(uv.ra, xlim=p.xlim, ylim=p.ylim)
#plot(pts, cex=0.5)
points(pts[ind,1:2], col="red")
points(pstop[ind,1:2], col="orange", pch=16)
## all points that haven't stopped have the pnow coordinates
pstop2 <- matrix(NA, nrow = nrow(track$ptrack), ncol = 3)
# not.stopped <- which(track2$stopindex==0)
# pstop2[not.stopped,] <- track2$pnow[not.stopped,]
## where the code says they stopped:
for (i in 1:nrow(track$ptrack)) {
  pstop2[i, ] <- track2$ptrack[i, 1:3, track2$stopindex[i]]
}
points(pstop2, col="blue", cex=0.7)

### tracking some individuals
## the 3D track until they stopped
for(i in 1:length(ind)){
  points(t(track$ptrack[ind[i],1:2,1:track$stopindex[ind[i]]]), col="red", cex=0.1)
}
## the 3D end locations
points(pstop[ind,1:2], col="red", pch=16)

## the 2D track from those points
for(i in 1:length(ind)){
  points(t(track2$ptrack[ind[i],1:2,]), col="purple", cex=0.1)
}

## the 2D end locations
points(pstop2[ind,1:2], col="blue", pch=16)

## the 2D track until they actually stopped
for(i in 1:length(ind)){
  points(t(track2$ptrack[ind[i],1:2,1:track2$stopindex[ind[i]]]), col="gold", cex=0.1)
}






###################################

#### FOURTH testrun: using loopit_2D3D

####
### 2D run (seeding with the correct locations)
track.l <-loopit_2D3D(pts = pstop, romsobject = testROMS, projected=TRUE, sedimentation=TRUE, speed=100, domain="2D", runtime=200)#, trajectories=TRUE)


## the 2D end locations
plot(uv.ra, xlim=p.xlim, ylim=p.ylim)
points(pstop2[,1:2], col="blue", pch=16, cex=0.7)
points(track.l$pend, col="grey", pch=16, cex=0.5)





###################################

#### FIFTH testrun: loop through ROMS files

####
### 2D run (seeding with the correct locations)
track.l7 <-loopit_2D3D(pts = pstop, romsobject = testROMS, projected=TRUE, sedimentation=TRUE, speed=100, domain="2D", runtime=200,
                      roms_slices = 7)


## the 2D end locations
plot(uv.ra, xlim=p.xlim, ylim=p.ylim)
points(pstop2[,1:2], col="blue", pch=16, cex=0.7)
points(track.l$pend, col="grey", pch=16, cex=0.5)
points(track.l7$pend, col="purple", pch=16, cex=0.5)


uv.ra.all <- uv.ra
uv.ra.all$lyr.2 <- c(sqrt(testROMS$i_u[,1399:1,1,2]^2+c(testROMS$i_v[,1399:1,1,2])^2))
uv.ra.all$lyr.3 <- c(sqrt(testROMS$i_u[,1399:1,1,3]^2+c(testROMS$i_v[,1399:1,1,3])^2))
uv.ra.all$lyr.4 <- c(sqrt(testROMS$i_u[,1399:1,1,4]^2+c(testROMS$i_v[,1399:1,1,4])^2))
uv.ra.all$lyr.5 <- c(sqrt(testROMS$i_u[,1399:1,1,5]^2+c(testROMS$i_v[,1399:1,1,5])^2))
uv.ra.all$lyr.6 <- c(sqrt(testROMS$i_u[,1399:1,1,6]^2+c(testROMS$i_v[,1399:1,1,6])^2))
uv.ra.all$lyr.7 <- c(sqrt(testROMS$i_u[,1399:1,1,7]^2+c(testROMS$i_v[,1399:1,1,7])^2))
plot(uv.ra.all, xlim=p.xlim, ylim=p.ylim)


#save(track, file="/pvol/3_model_analysis/testrun_500kparticles15days.Rdata")

## NOTES:
## - STILL AN UGLY VERSION, BECAUSE I DELETE ROWS & COLUMS FROM ROMS DATA TO MAKE THEM COMPARABLE
## - DELETE MORE POINTS FROM AREAS FAR OFF THE SHELF
## - DELETE POINTS FROM UNDER ICE-SHELVES
