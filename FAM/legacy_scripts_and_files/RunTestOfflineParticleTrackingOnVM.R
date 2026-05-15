#### 1) Set up ----
# load libraries
library(ncdf4)        ## package for netcdf manipulation
library(terra)
library(reproj) ## reproject coordinates
library(dplyr)
library(reshape) ## prep data for plotting
library(ggplot2) ## plotting

# set up directory pointers etc
env.dir <- "/pvol/data_environmental/ROMS_4k_files/"

string.chr <- "Circumpolar_EnvData_"
string.res <- "500m_"

# polar stereographic projection:
stereo <- "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

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

#### load data
load("/pvol/Scripts/ptrackR_package_update2024/ptrackr/data/toyROMS.rdata")

names(toyROMS)[1:2] <- c("x", "y")

# #### prep currents data for plotting NOT WORKING
# roms.coords.proj <- reproj(cbind(c(toyROMS$lon_u), c(toyROMS$lat_u)), target=stereo)
# x.range <- c(min(roms.coords.proj[,1])-200,max(roms.coords.proj[,1])+200)
# y.range <- c(min(roms.coords.proj[,2])-200,max(roms.coords.proj[,2])+200)
# empty.roms.ra <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=400)

#### toyrun:

## setup ROMS field
## setup pattern of points to seed into the model
pts <- cbind(c(toyROMS$x+0.1), c(toyROMS$y+0.1), as.vector(as.numeric(0)))

plot(pts)

## run 3D tracking
track <- trackit_3D(pts = pts, romsobject = toyROMS, projected=FALSE)
## checking results
plot(pts, cex=0.5)
points(track$pnow, col = "red", cex=0.5)
##
track2 <- trackit_2D(pts = track$pnow, romsobject = toyROMS)

###################################
#### analyse toyrun
dim(track$ptrack)
## compare depth levels where particles have stopped
track$ptrack[2,3,track$stopindex[2]]
track$stopindex[2]
track$pnow[2,3]

## pnow shows particles as if they had never encountered the bottom
plot(track$ptrack[1:10,1:2,track$stopindex])

points(track$pnow, col = "red", cex=0.5)


###################################
###### full practice run with toyROMS

####
### 3D run
track <- trackit_3D(pts = pts, romsobject = toyROMS)

### checking 3D results
## seeding points
plot(pts, cex=0.5)
## now
points(track$pnow, col = "red", cex=0.5)
## where the code says they stopped:
pstop <- matrix(NA, nrow = 2000, ncol = 3)
for (i in 1:2000) {
  pstop[i, ] <- track$ptrack[i, 1:3, track$stopindex[i]]
}
points(pstop, col="blue", cex=0.5)

### tracking some individuals
ind <- c(31,87,613,1425,1957)
## all points
plot(pts, cex=0.5)
points(pts[ind,1:2], col="red")
## the full track of those points
points(t(track$ptrack[ind[1],1:2,]), col="blue")
points(t(track$ptrack[ind[2],1:2,]), col="blue")
points(t(track$ptrack[ind[3],1:2,]), col="blue")
points(t(track$ptrack[ind[4],1:2,]), col="blue")
points(t(track$ptrack[ind[5],1:2,]), col="blue")
## only the track until they stopped
points(t(track$ptrack[ind[1],1:2,1:track$stopindex[ind[1]]]), col="red")
points(t(track$ptrack[ind[2],1:2,1:track$stopindex[ind[2]]]), col="red")
points(t(track$ptrack[ind[3],1:2,1:track$stopindex[ind[3]]]), col="red")
points(t(track$ptrack[ind[4],1:2,1:track$stopindex[ind[4]]]), col="red")
points(t(track$ptrack[ind[5],1:2,1:track$stopindex[ind[5]]]), col="red")
## the end locations
points(track$pnow[ind,1:2], col = "green", pch=16)
points(pstop[ind,1:2], col="orange", pch=16)

####
### 2D run (seeding with the correct locations)
track2 <- trackit_2D(pts = pstop, romsobject = toyROMS, projected=FALSE)

### checking 2D results
plot(pts, cex=0.5)
points(pts[ind,1:2], col="red")
points(pstop[ind,1:2], col="orange", pch=16)
## all points that haven't stopped have the pnow coordinates
pstop2 <- matrix(NA, nrow = 2000, ncol = 3)
not.stopped <- which(track2$stopindex==0)
pstop2[not.stopped,] <- track2$pnow[not.stopped,]
## where the code says they stopped:
for (i in c(1:2000)[-not.stopped]) {
  pstop2[i, ] <- track2$ptrack[i, 1:3, track2$stopindex[i]]
}
points(pstop2, col="blue", cex=0.7)

### tracking some individuals
## only the track until they stopped
for(i in 1:length(ind)){
points(t(track$ptrack[ind[i],1:2,1:track$stopindex[ind[i]]]), col="red")
}
## the end locations
points(pstop[ind,1:2], col="orange", pch=16)
## the full track of those points
for(i in 1:length(ind)){
points(t(track2$ptrack[ind[i],1:2,]), col="purple")
}
## only the track until they stopped
for(i in 1:length(ind)){
points(t(track2$ptrack[ind[i],1:2,1:track2$stopindex[ind[i]]]), col="gold")
}

####
### run and analyse toyrun with loopit_2D3D
## each step in runtime represents 6h, so e.g. 60 means 15 days
track.l <- loopit_2D3D(pts = pstop, romsobject = toyROMS, speed=100, domain="2D", runtime=200, trajectories=TRUE)

## pend has the location of all stopped particles after 200 steps except those in track.l$id_list[200]
## where to find the stopping location of e.g. particle 31? Manually checked and it's step 55
#which(sapply(track.l$id_list, FUN=function(X) 31 %in% X))

## plot the end position over the plot created earlier
points(track.l$pend, col="grey", pch=16, cex=0.5)

#### GREAT, STOPPING LOCATIONS MATCH!!!


