#### 1) Set up ----
# load libraries
library(ncdf4)## package for netcdf manipulation
library(tidync)## package for netcdf manipulation
library(terra)
library(reproj) ## reproject coordinates
library(dplyr)
library(reshape) ## prep data for plotting
library(ggplot2) ## plotting
library(abind) 

## set up directory pointers etc
env.dir <- "/pvol/data_environmental/ROMS_2k_files/"
env.dir3 <- "/pvol3TB/data_environmental/ROMS_2k_files/"

string.chr <- "Circumpolar_EnvData_"
string.res <- "500m_"

## polar stereographic projection:
stereo <- "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

################################
#### Reduce ROMS to fewer time-slices:
u.tifs <- list.files(env.dir, pattern="_slices_u", full.names=TRUE)
ra.u <- rast(u.tifs[1], lyrs=1)
for(i in 2:length(u.tifs)){
  ra.u <- c(ra.u, rast(u.tifs[i], lyrs=1))
}
ra.u

v.tifs <- list.files(env.dir, pattern="_slices_v", full.names=TRUE)
ra.v <- rast(v.tifs[1], lyrs=1)
for(i in 2:length(v.tifs)){
  ra.v <- c(ra.v, rast(v.tifs[i], lyrs=1))
}
ra.v

xlim=c(-500000,400000)
ylim=c(-2100000,-1200000)
range=c(-1,1)
par(mfrow=c(4,4))
plot(ra.u, xlim=xlim, ylim=ylim, range=range)
plot(subset(ra.u,c(1,9,17,24,25,26,32,33,34)), xlim=xlim, ylim=ylim, range=range)

plot(ra1, xlim=xlim, ylim=ylim, range=range)
plot(ra2, xlim=xlim, ylim=ylim, range=range)
plot(ra3, xlim=xlim, ylim=ylim, range=range)
plot(ra4, xlim=xlim, ylim=ylim, range=range)
plot(ra5, xlim=xlim, ylim=ylim, range=range)
plot(ra6, xlim=xlim, ylim=ylim, range=range)
plot(ra7, xlim=xlim, ylim=ylim, range=range)
plot(ra8, xlim=xlim, ylim=ylim, range=range)
plot(ra9, xlim=xlim, ylim=ylim, range=range)
plot(ra10, xlim=xlim, ylim=ylim, range=range)
plot(ra11, xlim=xlim, ylim=ylim, range=range)
plot(ra12, xlim=xlim, ylim=ylim, range=range)
plot(ra13, xlim=xlim, ylim=ylim, range=range)
plot(ra14, xlim=xlim, ylim=ylim, range=range)
plot(ra15, xlim=xlim, ylim=ylim, range=range)
plot(ra16, xlim=xlim, ylim=ylim, range=range)

plot(ra17, xlim=xlim, ylim=ylim, range=range)

plot(ra20, xlim=xlim, ylim=ylim, range=range)
plot(ra21, xlim=xlim, ylim=ylim, range=range)
plot(ra22, xlim=xlim, ylim=ylim, range=range)
plot(ra23, xlim=xlim, ylim=ylim, range=range)
plot(ra24, xlim=xlim, ylim=ylim, range=range)
plot(ra25, xlim=xlim, ylim=ylim, range=range)
plot(ra26, xlim=xlim, ylim=ylim, range=range)
plot(ra27, xlim=xlim, ylim=ylim, range=range)


## Which slices are most similar, and can be used to close a loop?
sum(Rdat$i_u[,,,1] - Rdat$i_u[,,,2], na.rm=TRUE)
sum(Rdat$i_u[,,,1] - Rdat$i_u[,,,3], na.rm=TRUE)
sum(Rdat$i_u[,,,1] - Rdat$i_u[,,,4], na.rm=TRUE)
sum(Rdat$i_u[,,,1] - Rdat$i_u[,,,5], na.rm=TRUE)
sum(Rdat$i_u[,,,1] - Rdat$i_u[,,,6], na.rm=TRUE)
sum(Rdat$i_u[,,,1] - Rdat$i_u[,,,7], na.rm=TRUE)
sum(Rdat$i_u[,,,1] - Rdat$i_u[,,,8], na.rm=TRUE)
sum(Rdat$i_u[,,,1] - Rdat$i_u[,,,9], na.rm=TRUE)







