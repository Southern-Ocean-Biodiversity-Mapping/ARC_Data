#### 1) Set up ----
# load libraries
library(ncdf4)        ## package for netcdf manipulation
library(terra)
library(reproj) ## reproject coordinates
library(dplyr)
library(reshape) ## prep data for plotting
library(ggplot2) ## plotting
library(spatstat) ## poisson point process
library(ppmData) ## poisson point process

## set up directory pointers etc
env.dir <- "/pvol/data_environmental/ROMS_2k_files/"

string.chr <- "Circumpolar_EnvData_"
string.res <- "500m_"

## polar stereographic projection:
stereo <- "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

################################
###### particle tracking code
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

################################
#### ROMS Currents
load(paste0(env.dir,"ocean_his_0001_2slices_TrackingSetup_transposed.Rdata"))
## reduce the ROMS by excluding areas that are not relevant for the shelf, such as North of the Ross Sea
Rdat <- ROMS_2k_t
# 
## reduce the array to speed up runtime
Rdat <- list()
Rdat$x    <- ROMS_2k_t$x[1:2500,1:3000]
Rdat$y    <- ROMS_2k_t$y[1:2500,1:3000]
Rdat$h    <- ROMS_2k_t$h[1:2500,1:3000]
Rdat$hh   <- ROMS_2k_t$hh[1:2500,1:3000,]
Rdat$zice <- ROMS_2k_t$zice[1:2500,1:3000]
Rdat$i_u  <- ROMS_2k_t$i_u[1:2500,1:3000,,]
Rdat$i_v  <- ROMS_2k_t$i_v[1:2500,1:3000,,]
Rdat$i_w  <- ROMS_2k_t$i_w[1:2500,1:3000,,]
rm(ROMS_2k_t)

################################
###### create uv and h raster for plotting
roms.coords.proj <- cbind(c(Rdat$x), c(Rdat$y))
x.range <- c(min(roms.coords.proj[,1])-1000,max(roms.coords.proj[,1])+1000)
y.range <- c(min(roms.coords.proj[,2])-1000,max(roms.coords.proj[,2])+1000)
empty.roms.ra <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000)

uv.ra      <- setValues(empty.roms.ra, sqrt(Rdat$i_u[,,1,1]^2 +c(Rdat$i_v[,,1,1])^2))
uv.surf.ra <- setValues(empty.roms.ra, sqrt(Rdat$i_u[,,31,1]^2+c(Rdat$i_v[,,31,1])^2))
h.ra       <- setValues(empty.roms.ra, Rdat$h)

################################
###### coastline for plotting
# coast.unproj <- vect("/pvol/data_environmental/antarctic_coastline_2023/add_coastline_high_res_polygon_v7_8.shp")
# coast.proj <- project(coast.unproj, uv.ra)

###################################
#### visually check model setup
## 3D current field:
plot(uv.ra)
# plot(coast.proj, add=TRUE)

################################
## scroll down from here to find individual test runs and continue there
################################

###################################
## pts seeded only near the shelf break in the Ross Sea
hz <- Rdat$h+Rdat$zice ## need to combine bathymetry with ice draft
hz[which(is.na(Rdat$i_w[,,1,1]))] <- NA
hz[hz>4000] <- NA
hz[2500:2800,] <- NA
hz[,3000:3150] <- NA
hz[2200:2800,1:1200] <- NA
hz[1400:1600,2000:2300] <- NA
h.df <- melt(hz, varnames = c("row", "col"), value.name = "value")
##
not.na <- which(!is.na(hz))
pts.all <- cbind(c(roms.coords.proj[,1]), c(roms.coords.proj[,2]), as.vector(as.numeric(0)))
pts.all2 <- pts.all[not.na,]
RS.sel <- which(pts.all2[,1]>-500000 & pts.all2[,1]<0 & pts.all2[,2]<=-1400000 & pts.all2[,2]>=-1700000)
set.seed(1)
pts <- pts.all2[sample(RS.sel,1000),]

### check distribution of points
# plot(uv.ra)
# points(pts)
plot(uv.ra, xlim=c(-550000,50000), ylim=c(-1750000,-1350000))
points(pts)



#####################################
## FIRST TESTRUN

####
### run 3D tracking
## the below times are trackit_3D only
## 1min for 1k particles and 1 day
## 1.8min for 10k particles and 1 day
## 5.8min for 100k particles and 1 day (printed steps starting at 1s per step, then slowing down to 4s per step by day 1)
## 7.3hrs for 550k particles and 15 days (= 12GB file output)
start.time <- Sys.time()
track <- trackit_3D(pts = pts, romsobject = Rdat, time=30, w_sink=100, projected=TRUE)
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
track2 <- trackit_2D(pts = pstop, romsobject = Rdat, projected=TRUE, sedimentation=TRUE)

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
for (i in c(1:1000)[-not.stopped]) {
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


#####################################
## SECOND TESTRUN: small region, one particle per skm (so 4 per 2km cell)

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
### run 3D tracking (5+1min for 30 days and 10kparticles)
start.time <- Sys.time()
track <- trackit_3D(pts = pts, romsobject = Rdat, time=30, w_sink=100, projected=TRUE)
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
  points(track$ptrack[ind[i],1:2,], col="blue")
}
## only the track until they stopped
for(i in 1:length(ind)){
  points(t(track$ptrack[ind[i],1:2,1:track$stopindex[ind[i]]]), col="red")
}
## the end locations
points(track$pnow[ind,1:2], col = "green", pch=16)
points(pstop[ind,1:2], col="orange", pch=16)

### tracking some individuals AFTER 1 DAY
ind <- sample(1:nrow(track$ptrack),30)#c(31,87,425,613,957)
## all points
plot(uv.ra, xlim=p.xlim, ylim=p.ylim)
# plot(pts, cex=0.5)
points(pts[ind,1:2], col="red")
## the track of those points until after 24h
for(i in 1:length(ind)){
  points(t(track$ptrack[ind[i],1:2,1:48]), col="blue")
}
## the track of those points until after 48h
for(i in 1:length(ind)){
  points(t(track$ptrack[ind[i],1:2,49:96]), col="purple")
}
## the depth of those points until after 48h
plot(track$ptrack[ind,3,1:96])

####
### 2D run 5+2min for 50 days ()
track2 <- trackit_2D(pts = pstop, romsobject = Rdat, projected=TRUE, sedimentation=TRUE)

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
#### THIRD testrun: loop through ROMS files

####
### 2D run = 41min for 10k particles with sedimentation, finished after 130 loops
start.time <- Sys.time()
track.l2 <-loopit_2D3D(pts = pstop, romsobject = Rdat, projected=TRUE, sedimentation=TRUE, speed=100, domain="2D", runtime=200,
                       roms_slices = 2)
end.time <- Sys.time()
end.time-start.time

## the 2D end locations
plot(uv.ra, xlim=p.xlim, ylim=p.ylim)
points(pstop2[,1:2], col="blue", pch=16, cex=0.7)
points(track.l$pend, col="grey", pch=16, cex=0.5)
points(track.l2$pend, col="purple", pch=16, cex=0.5)

#save(track, track2, track.l2, file="/pvol/3_model_analysis/testrun_10kparticles_sedimentation.Rdata")


###################################
#### Fourth testrun: loop through ROMS files in one region, but deleting areas of the shelf, or other regions entirely

#### ROMS Currents
load(paste0(env.dir,"ocean_his_0001_2slices_TrackingSetup_transposed.Rdata"))
## reduce the ROMS by excluding areas that are not relevant for the shelf, such as North of the Ross Sea
Rdat <- list()
Rdat$x    <- ROMS_2k_t$x[1:2500,1:3000]
Rdat$y    <- ROMS_2k_t$y[1:2500,1:3000]
Rdat$h    <- ROMS_2k_t$h[1:2500,1:3000]
Rdat$hh   <- ROMS_2k_t$hh[1:2500,1:3000,]
Rdat$zice <- ROMS_2k_t$zice[1:2500,1:3000]
Rdat$i_u  <- ROMS_2k_t$i_u[1:2500,1:3000,,]
Rdat$i_v  <- ROMS_2k_t$i_v[1:2500,1:3000,,]
Rdat$i_w  <- ROMS_2k_t$i_w[1:2500,1:3000,,]

rm(ROMS_2k_t)

###### create uv and h raster for plotting
roms.coords.proj <- cbind(c(Rdat$x), c(Rdat$y))
x.range <- c(min(roms.coords.proj[,1])-1000,max(roms.coords.proj[,1])+1000)
y.range <- c(min(roms.coords.proj[,2])-1000,max(roms.coords.proj[,2])+1000)
empty.roms.ra <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000)
uv.ra      <- setValues(empty.roms.ra, sqrt(Rdat$i_u[,,1,1]^2 +c(Rdat$i_v[,,1,1])^2))
uv.surf.ra <- setValues(empty.roms.ra, sqrt(Rdat$i_u[,,31,1]^2+c(Rdat$i_v[,,31,1])^2))
h.ra       <- setValues(empty.roms.ra, Rdat$h)
plot(uv.ra)

## pts seeded only near the shelf break in the Ross Sea
hz <- Rdat$h+Rdat$zice ## need to combine bathymetry with ice draft
hz[which(is.na(Rdat$i_w[,,1,1]))] <- NA
hz[hz>4000] <- NA
hz[2200:2500,1:1200] <- NA
hz[1400:1600,2000:2300] <- NA
h.df <- melt(hz, varnames = c("row", "col"), value.name = "value")
##
not.na <- which(!is.na(hz))
pts.all <- cbind(c(roms.coords.proj[,1]), c(roms.coords.proj[,2]), as.vector(as.numeric(0)))
pts.all2 <- pts.all[not.na,]
RS.sel <- which(pts.all2[,1]>-500000 & pts.all2[,1]<0 & pts.all2[,2]<=-1400000 & pts.all2[,2]>=-1700000)
set.seed(1)
pts <- pts.all2[sample(RS.sel,1000),]
plot(uv.ra)
points(pts.all2[sample(1:nrow(pts.all2),10000),])


####
### 2D run
start.time <- Sys.time()
track.l2 <-loopit_2D3D(pts = pstop, romsobject = Rdat, projected=TRUE, sedimentation=TRUE, speed=100, domain="2D", runtime=200,
                       roms_slices = 2)
end.time <- Sys.time()
end.time-start.time

## the 2D end locations
plot(uv.ra, xlim=p.xlim, ylim=p.ylim)
points(pstop2[,1:2], col="blue", pch=16, cex=0.7)
points(track.l$pend, col="grey", pch=16, cex=0.5)
points(track.l2$pend, col="purple", pch=16, cex=0.5)




###################################
###################################
###################################
#### FIFTH testrun: ###############
## loop through ROMS files in one region, but deleting other regions entirely

#### ROMS Currents
load(paste0(env.dir,"ocean_his_0001_2slices_TrackingSetup_transposed.Rdata"))
## reduce the ROMS by excluding entire areas
## e.g. Ross Sea only
row.lim <- 1501:2500
col.lim <- 1251:1750
Rdat <- list()
Rdat$x    <- ROMS_2k_t$x[row.lim, col.lim]
Rdat$y    <- ROMS_2k_t$y[row.lim, col.lim]
Rdat$h    <- ROMS_2k_t$h[row.lim, col.lim]
Rdat$hh   <- ROMS_2k_t$hh[row.lim, col.lim,]
Rdat$zice <- ROMS_2k_t$zice[row.lim, col.lim]
Rdat$i_u  <- ROMS_2k_t$i_u[row.lim, col.lim,,]
Rdat$i_v  <- ROMS_2k_t$i_v[row.lim, col.lim,,]
Rdat$i_w  <- ROMS_2k_t$i_w[row.lim, col.lim,,]

rm(ROMS_2k_t)

###### create uv and h raster for plotting
roms.coords.proj <- cbind(c(Rdat$x), c(Rdat$y))
x.range <- c(min(roms.coords.proj[,1])-1000,max(roms.coords.proj[,1])+1000)
y.range <- c(min(roms.coords.proj[,2])-1000,max(roms.coords.proj[,2])+1000)
empty.roms.ra <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000)
uv.ra      <- setValues(empty.roms.ra, sqrt(Rdat$i_u[,,1,1]^2 +c(Rdat$i_v[,,1,1])^2))
uv.surf.ra <- setValues(empty.roms.ra, sqrt(Rdat$i_u[,,31,1]^2+c(Rdat$i_v[,,31,1])^2))
h.ra       <- setValues(empty.roms.ra, Rdat$h)
plot(uv.ra)

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
### run 3D tracking (1min for 30 days and 10kparticles)
start.time <- Sys.time()
track <- trackit_3D(pts = pts, romsobject = Rdat, time=30, w_sink=100, projected=TRUE)
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

####
### 2D run = 2.5min for 10k particles with sedimentation, finished after 130 loops
start.time <- Sys.time()
track.l2 <-loopit_2D3D(pts = pstop, romsobject = Rdat, projected=TRUE, sedimentation=TRUE, speed=100, domain="2D", runtime=200,
                       roms_slices = 2)
end.time <- Sys.time()
end.time-start.time

## the 2D end locations
plot(uv.ra, xlim=p.xlim, ylim=p.ylim)
points(pstop2[,1:2], col="blue", pch=16, cex=0.7)
points(track.l2$pend, col="purple", pch=16, cex=0.5)


###################################
#### SIXTH testrun: ###############
## loop through ROMS files in one region, more particles in a larger region

xlim=c(-500000,300000)
ylim=c(-2100000,-1350000)

## create points in that region
x.seq <- seq(xlim[1],xlim[2], by=4000)
y.seq <- seq(ylim[1],ylim[2], by=4000)
x.dat <- rep(x.seq, length(y.seq))
y.dat <- rep(y.seq[1], length(x.seq))
for(i in 1:length(y.seq)){
  y.dat <- c(y.dat, rep(y.seq[i], length(x.seq)))
}
pts <- cbind(x.dat, y.dat, 0)

### check distribution of points
plot(uv.ra)
points(pts, cex=0.1)

#### 3D tracking ####
### 13min for 30 days and 38k particles
start.time <- Sys.time()
track <- trackit_3D(pts = pts, romsobject = Rdat, time=30, w_sink=100, projected=TRUE)
end.time <- Sys.time()
end.time-start.time

### checking 3D results
plot(uv.ra, ylim=c(-2200000,-1300000))
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
points(pstop, col="blue", cex=0.1)

#### 2D run
### 7.5min for 38k particles with sedimentation, finished after 400 loops
start.time <- Sys.time()
track.l <-loopit_2D3D(pts = pstop, romsobject = Rdat, projected=TRUE, sedimentation=TRUE, speed=100, domain="2D", runtime=200, roms_slices = 2)
end.time <- Sys.time()
end.time-start.time

## the 2D end locations
plot(uv.ra)
points(pstop[,1:2], col="blue", pch=16, cex=0.7)
points(track.l$pend, col="purple", pch=16, cex=0.5)


###################################
#### SEVENTH testrun: ###############
## reduce ROMS dimensions for 2D run, by removing upper water column

Rdat.l <- list()
Rdat.l$x    <- Rdat$x
Rdat.l$y    <- Rdat$y
Rdat.l$h    <- Rdat$h
Rdat.l$hh   <- Rdat$hh[,,1:10]
Rdat.l$zice <- Rdat$zice
Rdat.l$i_u  <- Rdat$i_u[,,1:10,]
Rdat.l$i_v  <- Rdat$i_v[,,1:10,]
Rdat.l$i_w  <- Rdat$i_w[,,1:10,]

#### 2D run
### 3 min for 38k particles with sedimentation, finished after 400 loops
start.time <- Sys.time()
track.l.l <-loopit_2D3D(pts = pstop, romsobject = Rdat.l, projected=TRUE, sedimentation=TRUE, speed=100, domain="2D", runtime=200, roms_slices = 2)
end.time <- Sys.time()
end.time-start.time

## the 2D end locations
plot(uv.ra)
points(pstop[,1:2], col="blue", pch=16, cex=0.7)
points(track.l.l$pend, col="purple", pch=16, cex=0.5)

#### 2D run - larger particles
### 3 min for 38k particles with sedimentation, finished after 400 loops
start.time <- Sys.time()
track.l.l32 <-loopit_2D3D(pts = pstop, romsobject = Rdat.l, projected=TRUE, sedimentation=TRUE, speed=100, domain="2D", runtime=200, roms_slices = 2, particle_radius = 0.00032)
end.time <- Sys.time()
end.time-start.time

## the 2D end locations
points(track.l.l32$pend, col="gold", pch=16, cex=0.5)

###################################
###################################
###################################

#### SUMMARY ####
## 1. reduce domain to smaller region for faster processing
## 2. remove upper layers of ROMS file for 2D run for faster processing

###################################
###################################
###################################

#### FINAL testrun ###############

## reduce the ROMS by excluding entire areas
# load(paste0(env.dir,"ocean_his_0001_2slices_TrackingSetup_transposed.Rdata"))
# row.lim <- 1501:2500
# col.lim <- 1251:1750
# Rdat <- list()
# Rdat$x    <- ROMS_2k_t$x[row.lim, col.lim]
# Rdat$y    <- ROMS_2k_t$y[row.lim, col.lim]
# Rdat$h    <- ROMS_2k_t$h[row.lim, col.lim]
# Rdat$hh   <- ROMS_2k_t$hh[row.lim, col.lim,]
# Rdat$zice <- ROMS_2k_t$zice[row.lim, col.lim]
# Rdat$i_u  <- ROMS_2k_t$i_u[row.lim, col.lim,,]
# Rdat$i_v  <- ROMS_2k_t$i_v[row.lim, col.lim,,]
# Rdat$i_w  <- ROMS_2k_t$i_w[row.lim, col.lim,,]
# rm(ROMS_2k_t)
# save(Rdat, file=paste0(env.dir,"ocean_his_0001_2slices_TrackingSetup_transposed_RossSea.Rdata"))
load(paste0(env.dir,"ocean_his_0001_2slices_TrackingSetup_transposed_RossSea.Rdata"))
          
## create uv and h raster for plotting
roms.coords.proj <- cbind(c(Rdat$x), c(Rdat$y))
x.range <- c(min(roms.coords.proj[,1])-1000,max(roms.coords.proj[,1])+1000)
y.range <- c(min(roms.coords.proj[,2])-1000,max(roms.coords.proj[,2])+1000)
empty.roms.ra <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000)
uv.ra      <- setValues(empty.roms.ra, sqrt(Rdat$i_u[,,1,1]^2 +c(Rdat$i_v[,,1,1])^2))
uv.surf.ra <- setValues(empty.roms.ra, sqrt(Rdat$i_u[,,31,1]^2+c(Rdat$i_v[,,31,1])^2))
h.ra       <- setValues(empty.roms.ra, Rdat$h)
hh.ra      <- setValues(empty.roms.ra, Rdat$hh[,,1])
plot(uv.ra)

## seeded with particles based on NPP values
npp <- rast("/pvol/data_environmental/Circumpolar_EnvData_2km_NPP_Cafe_filled_SummerAverage.tif")
npp2 <- project(npp, h.ra)
plot(npp2)
Z <- terra2im(npp2/100000000) ## library(ppmData)
pp <- rpoispp(Z) ## library(spatstat)
pp
length(which(!is.na(npp2[])))
npp.pts <- cbind(pp$x, pp$y,0)

### check distribution of points
# plot(uv.ra)
# points(npp.pts, cex=0.1)

#### 3D tracking ####
### 20 min for 5 days and 46k particles
### 2.5h for 5 days and 541k particles, 300MB file
### 1.06 days for 5 days and 5.4m particles (~20 per cell), 3GB file
start.time <- Sys.time()
track.t <- loopit_2D3D(pts = npp.pts, romsobject = Rdat, projected=TRUE, speed=100, domain="3D", runtime=5, roms_slices = 2, looping_time = 0.5)
end.time <- Sys.time()
end.time-start.time
### 18 min for 5 days and 46k particles
### 3h for 5 days and 541k particles, 4.6GB file for Ross Sea data
start.time2 <- Sys.time()
track2 <- trackit_3D(pts = npp.pts, romsobject = Rdat, time=5, w_sink=100, projected=TRUE)
end.time2 <- Sys.time()
end.time2-start.time2
#save(track2, file="/pvol/3_model_analysis/testrun_540kparticles_3Dtracking_5days.Rdata")
save(track.t, file="/pvol/3_model_analysis/testruns_5404kparticles_3Dtracking_5days.Rdata")



#### 3D tracking ####
### 20 min for 5 days and 46k particles
start.time <- Sys.time()
track <- loopit_2D3D(pts = npp.pts, romsobject = Rdat, projected=TRUE, speed=100, domain="3D", runtime=5, roms_slices = 2, looping_time = 0.5)
end.time <- Sys.time()
end.time-start.time
### 10 min for 5 days and 46k particles
start.time2 <- Sys.time()
track2 <- trackit_3D(pts = npp.pts, romsobject = Rdat, time=5, w_sink=100, projected=TRUE)
end.time2 <- Sys.time()
end.time2-start.time2

### checking 3D results
plot(uv.ra, ylim=c(-2200000,-1300000))
## seeding points
points(npp.pts, cex=0.1, pch=16)
## all points that haven't stopped have the pnow coordinates
pstop <- matrix(NA, nrow = nrow(track$ptrack), ncol = 3)
not.stopped <- which(track$stopindex==0)
pstop[not.stopped,] <- track$pnow[not.stopped,]
## where the code says they stopped:
for (i in c(1:nrow(track$ptrack))[-not.stopped]) {
  pstop[i, ] <- track$ptrack[i, 1:3, track$stopindex[i]]
}
points(pstop, col="blue", cex=0.1)


pstop <- matrix(NA, nrow = nrow(track.t$ptrack), ncol = 3)
not.stopped <- which(track.t$stopindex==0)
pstop[not.stopped,] <- track.t$pnow[not.stopped,]
## where the code says they stopped:
for (i in c(1:nrow(track.t$ptrack))[-not.stopped]) {
  pstop[i, ] <- track.t$ptrack[i, 1:3, track.t$stopindex[i]]
}

#### 2D tracking ####
Rdat.l <- Rdat
Rdat.l$hh <- Rdat$hh[,,1:8]
Rdat.l$i_u <- Rdat$i_u[,,1:8,]
Rdat.l$i_v <- Rdat$i_v[,,1:8,]
Rdat.l$i_w <- Rdat$i_w[,,1:8,]

### ?? min for 540k particles with sedimentation, finished after ?? loops (400 loops = 50 days)
### ?? min for 540k particles with sedimentation, finished after ?? loops (40 loops = 5 days)
start.time <- Sys.time()
track.t.2D <-loopit_2D3D(pts = pstop, romsobject = Rdat.l, projected=TRUE, sedimentation=TRUE, speed=100, domain="2D", runtime=5, looping_time = 0.5, roms_slices = 2)
end.time <- Sys.time()
end.time-start.time

## the 2D end locations
plot(uv.ra)
points(pstop[,1:2], col="blue", pch=16, cex=0.1)
points(track.l.l$pend, col="purple", pch=16, cex=0.1)

# #### 2D run - larger particles
# ### ? min for 264k particles with sedimentation, finished after ?? loops
# start.time <- Sys.time()
# track.l.l32 <-loopit_2D3D(pts = pstop, romsobject = Rdat.l, projected=TRUE, sedimentation=TRUE, speed=100, domain="2D", runtime=200, roms_slices = 2, particle_radius = 0.00032)
# end.time <- Sys.time()
# end.time-start.time
# ## the 2D end locations
# points(track.l.l32$pend, col="gold", pch=16, cex=0.5)



## !!! TEST LOOPIT_2D3D FOR THE 3D PART!! IS IT FASTER???
