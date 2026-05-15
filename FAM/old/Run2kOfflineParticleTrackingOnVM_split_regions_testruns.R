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
#source(paste0(ptrackr_dir,"loopit_2D3D.R"))
source(paste0(ptrackr_dir,"loopit_2D3D_traj.R"))

################################
#### ROMS Currents
# load(paste0(env.dir,"ocean_his_0001_2slices_TrackingSetup_transposed.Rdata"))
## reduce the ROMS by excluding areas that are not relevant for the shelf, such as North of the Ross Sea
# Rdat <- ROMS_2k_t
# # 
# ## reduce the array to speed up runtime
# Rdat <- list()
# Rdat$x    <- ROMS_2k_t$x[1:2500,1:3000]
# Rdat$y    <- ROMS_2k_t$y[1:2500,1:3000]
# Rdat$h    <- ROMS_2k_t$h[1:2500,1:3000]
# Rdat$hh   <- ROMS_2k_t$hh[1:2500,1:3000,]
# Rdat$zice <- ROMS_2k_t$zice[1:2500,1:3000]
# Rdat$i_u  <- ROMS_2k_t$i_u[1:2500,1:3000,,]
# Rdat$i_v  <- ROMS_2k_t$i_v[1:2500,1:3000,,]
# Rdat$i_w  <- ROMS_2k_t$i_w[1:2500,1:3000,,]
# rm(ROMS_2k_t)

## reduce the ROMS by excluding entire areas
# load(paste0(env.dir,"ocean_his_0001_2slices_TrackingSetup_transposed.Rdata"))
# row.lim <- 1501:2500
# col.lim <- 1201:1760
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

###############################
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

#### One particle per cell with npp value in it ###############
## create points in that region
# x.seq <- seq(min(roms.coords.proj[,1]),max(roms.coords.proj[,1]), by=1000)
# y.seq <- seq(min(roms.coords.proj[,2]),max(roms.coords.proj[,2]), by=1000)
# x.dat <- rep(x.seq, length(y.seq))
# y.dat <- rep(y.seq[1], length(x.seq))
# for(i in 1:length(y.seq)){
#   y.dat <- c(y.dat, rep(y.seq[i], length(x.seq)))
# }
# pts <- cbind(x.dat, y.dat, 0)
## seeded with particles based on NPP values
npp <- rast("/pvol/data_environmental/Circumpolar_EnvData_2km_NPP_Cafe_filled_SummerAverage.tif")
npp2 <- project(npp, h.ra)
npp2[npp2[]>0] <- 1
npp2[600:1000,1:100] <- NA
npp2[900:1000,400:560] <- NA
## thin out to one particle every 4 cells
npp2[seq(1,1000, by=2),] <- NA
npp2[,seq(1,560, by=2)] <- NA
##
cell.pts <- cbind(crds(as.points(npp2)),0)


###############################
#### 3D tracking ##############
###############################

start.time <- Sys.time() #1.5h
track.3D <- loopit_2D3D(pts = cell.pts, romsobject = Rdat, projected=TRUE, speed=100, domain="3D", runtime=30, roms_slices = 2, looping_time = 0.5, trajectories=TRUE)
end.time <- Sys.time()
end.time-start.time
start.time <- Sys.time() #2.8h
track.3D.s <- trackit_3D(pts = cell.pts, romsobject = Rdat, projected=TRUE, w_sink=100, time=30)
end.time <- Sys.time()
end.time-start.time
#save(track.3D, file="/pvol/3_model_analysis/split_regions_runs_3Dtracking_30days.Rdata")
#save(track.3D.s, file="/pvol/3_model_analysis/split_regions_runs_3Dtracking2_30days.Rdata")
# load("/pvol/3_model_analysis/split_regions_runs_3Dtracking_30days.Rdata")
# load("/pvol/3_model_analysis/split_regions_runs_3Dtracking2_30days.Rdata")

#### 3D: extract results ####

### from LOOPIT_2D3D
## all points that haven't stopped have stopindex 0
## "pts" and "pnow" have some that have settled in the last step (calculated below as "pfloat")
not.stopped <- which(track.3D$stopindex==0)
pfloat <- matrix(track.3D$ptrack[not.stopped,,dim(track.3D$ptrack)[3]],ncol=3)
## the ones that stopped are "pend"
pstop <- track.3D$pend

### from TRACKIT_3D
## floating
not.stopped.s <- which(track.3D.s$stopindex==0)
pfloat.s <- matrix(track.3D.s$ptrack[not.stopped.s,,dim(track.3D.s$ptrack)[3]],ncol=3)
## stopped
pstop.s <- matrix(NA, nrow = nrow(track.3D.s$ptrack[-not.stopped.s,,1]), ncol = 3)
stop.sel <- c(1:nrow(track.3D.s$ptrack))[-not.stopped.s]
for (i in 1:nrow(pstop.s)) {
  sel <- stop.sel[i]
  pstop.s[i, ] <- track.3D.s$ptrack[sel, 1:3, track.3D.s$stopindex[sel]]
}

#### plot 3D results ####

xlim.p <- c(-500000,-100000)
ylim.p <- c(-2000000,-1400000)
par(mfrow=c(1,2))

plot(uv.ra, xlim=xlim.p, ylim=ylim.p, main="loopit_2D3D results")
points(cell.pts, cex=0.1, pch=16)
#points(track.3D$pend[,1:2], col="blue", cex=0.1)
points(pstop, col="blue", cex=0.1)
points(pfloat, col="purple", cex=0.1)

plot(uv.ra, xlim=xlim.p, ylim=ylim.p, main="trackit_3D results")
points(cell.pts, cex=0.1, pch=16)
points(pstop.s, col="blue", cex=0.1)
points(pfloat.s, col="purple", cex=0.1)

#### 3D particle tracks ####

ind <- sample(seq(1,nrow(cell.pts)), 10000)

### from LOOPIT_2D3D
## x_list and y_list create "pend"
## x_list - the particles that stopped in each loop - list length is the number of loops - list contents are x coordinates of points
## y_list - as above but for y coordinates of points
## idx_list_2D - cell-indices of each moving pts from each time-slice - list length is the number of loops - contains lists with 24? half-hour outputs
## idx_list_2D - these cells have been visited by particles... relevant for flux calculations!!!

par(mfrow=c(1,2))
rainbow.cols <- rainbow(60)
plot(uv.ra, xlim=xlim.p, ylim=ylim.p, main="loopit_2D3D results")
for(i in 1:60){
  ## identify which row the positions of the selected particles are stored (because particles get deleted every time-step)
  pts.sel <- which(track.3D$id_list[[i]]%in%ind)
  ## simplify the list of cells into a matrix
  idx.m <- matrix(unlist(track.3D$idx_list_2D[[i]]), ncol=length(track.3D$idx_list_2D[[i]]))
  ## the cell indices for all selected particles at that time-step
  c.idx <- c(idx.m[pts.sel])
  ## plot
  points(Rdat$x[c.idx], Rdat$y[c.idx], cex=0.1, col=rainbow.cols[i])
}

### from TRACKIT_3D
plot(uv.ra, xlim=xlim.p, ylim=ylim.p, main="trackit_3D results")
## only the track until they stopped
for(i in 1:length(ind)){
  points(t(track.3D.s$ptrack[ind[i],1:2,1:track.3D.s$stopindex[ind[i]]]), col=rainbow(dim(track.3D.s$ptrack)[3]), cex=0.1)
}

###############################
#### 2D tracking ##############
###############################

## NO SEDIMENTATION BECAUSE WE WANT POTENTIAL TRACKS ONLY TO DECIDE DOMAIN OVERLAP

pts.3D <- rbind(pstop, pfloat) ## use the loopit output
#pts.3D.s <- rbind(pstop.s, pfloat.s)

###############################
#### 2D tracking #### 
Rdat.l <- Rdat
Rdat.l$hh <- Rdat$hh[,,1:8]
Rdat.l$i_u <- Rdat$i_u[,,1:8,]
Rdat.l$i_v <- Rdat$i_v[,,1:8,]
Rdat.l$i_w <- Rdat$i_w[,,1:8,]

start.time <- Sys.time() ## 33min, (3.6min with sedimentation=TRUE)
track.2D <-loopit_2D3D(pts = pts.3D, romsobject = Rdat.l, projected=TRUE, speed=100, domain="2D", runtime=30, looping_time = 0.5, roms_slices = 2, trajectories = TRUE)
end.time <- Sys.time()
end.time-start.time
#save(track.2D, file="/pvol/3_model_analysis/split_regions_runs_2Dtracking_30days.Rdata")

start.time <- Sys.time() ## 33min
track.2D.traj <-loopit_2D3D(pts = pts.3D, romsobject = Rdat.l, projected=TRUE, speed=100, domain="2D", runtime=30, looping_time = 0.5, roms_slices = 2, trajectories = TRUE, detailed_trajectories = TRUE)
end.time <- Sys.time()
end.time-start.time
#save(track.2D.traj, file="/pvol/3_model_analysis/split_regions_runs_2Dtracking_30days_detailed.Rdata")

# start.time <- Sys.time() ## 38.5min
# track.2D.s <-trackit_2D(pts = pts.3D, romsobject = Rdat.l, projected=TRUE, w_sink=100, time=30)
# end.time <- Sys.time()
# end.time-start.time
# #save(track.2D.s, file="/pvol/3_model_analysis/split_regions_runs_2Dtracking2_30days.Rdata")

## WHICH PARTICLES HAVE LEFT THE DOMAIN?
## MAX DISTANCE TRAVELLED
## MAX DISTANCE FROM LEFT/RIGHT SIDE
## ROSS SEA BOUNDARIES

### from LOOPIT_2D3D
## all points that haven't stopped have stopindex 0
## "pts" and "pnow" have some that have settled in the last step (calculated below as "pfloat")
not.stopped <- which(track.2D$stopindex==0)
pfloat <- matrix(track.2D$ptrack[not.stopped,,dim(track.2D$ptrack)[3]],ncol=3)
## the ones that stopped are "pend"
pstop <- track.2D$pend

# ### from TRACKIT_2D !!! DOESN'T WORK !!!
# ## floating
# not.stopped.s <- which(track.2D.s$stopindex==0)
# pfloat.s <- matrix(track.2D.s$ptrack[not.stopped.s,,dim(track.2D.s$ptrack)[3]],ncol=3)
# ## stopped
# pstop.s <- matrix(NA, nrow = nrow(track.2D.s$ptrack[-not.stopped.s,,1]), ncol = 3)
# stop.sel <- c(1:nrow(track.2D.s$ptrack))[-not.stopped.s]
# for (i in 1:nrow(pstop.s)) {
#   sel <- stop.sel[i]
#   pstop.s[i, ] <- track.2D.s$ptrack[sel, 1:3, track.2D.s$stopindex[sel]]
# }

#### plot 2D results ####
xlim.p <- c(-500000,-100000)
ylim.p <- c(-2000000,-1400000)
par(mfrow=c(1,2))

plot(uv.ra, xlim=xlim.p, ylim=ylim.p, main="loopit_2D3D results")
points(cell.pts, cex=0.1, pch=16)
points(pstop, col="blue", cex=0.1)
points(pfloat, col="purple", cex=0.1)

#### 2D particle tracks ####

ind <- sample(seq(1,nrow(cell.pts)), 500)

rainbow.cols <- rainbow(60)
plot(uv.ra)#, xlim=xlim.p, ylim=ylim.p, main="loopit_2D3D results")
for(i in 1:60){
  ## identify which row the positions of the selected particles are stored (because particles get deleted every time-step)
  pts.sel <- which(track.2D$id_list[[i]]%in%ind)
  ## simplify the list of cells into a matrix
  idx.m <- matrix(unlist(track.2D$idx_list_2D[[i]]), ncol=length(track.2D$idx_list_2D[[i]]))
  ## the cell indices for all selected particles at that time-step
  c.idx <- c(idx.m[pts.sel])
  ## plot
  points(Rdat$x[c.idx], Rdat$y[c.idx], cex=0.1, col=rainbow.cols[i])
}

rainbow.cols <- rainbow(60)
plot(uv.ra)#, xlim=xlim.p, ylim=ylim.p, main="loopit_2D3D results")
for(i in 1:60){
  ## identify which row the positions of the selected particles are stored (because particles get deleted every time-step)
  pts.sel <- which(track.2D.traj$id_list[[i]]%in%ind)
  ## simplify the list of cells into a matrix
  xy.m <- cbind(track.2D.traj$ptrack_x_list[[i]], track.2D.traj$ptrack_y_list[[i]])
  ## plot
  points(xy.m[pts.sel,], cex=0.1, col=rainbow.cols[i])
}

###################################
#### 3D and 2D particle tracks ####
###################################

# load("/pvol/3_model_analysis/split_regions_runs_3Dtracking_30days.Rdata")
# load("/pvol/3_model_analysis/split_regions_runs_2Dtracking_30days.Rdata")

ind <- sample(seq(1,nrow(cell.pts)), 1000)

### from LOOPIT_2D3D
par(mfrow=c(1,2))
rainbow.cols <- rainbow(60)
plot(uv.ra, main="3D results")
for(i in 1:60){
  ## identify which row the positions of the selected particles are stored (because particles get deleted every time-step)
  pts.sel <- which(track.3D$id_list[[i]]%in%ind)
  ## simplify the list of cells into a matrix
  idx.m <- matrix(unlist(track.3D$idx_list_2D[[i]]), ncol=length(track.3D$idx_list_2D[[i]]))
  ## the cell indices for all selected particles at that time-step
  c.idx <- c(idx.m[pts.sel])
  ## plot
  points(Rdat$x[c.idx], Rdat$y[c.idx], cex=0.1, col=rainbow.cols[i])
}
plot(uv.ra, main="2D results")
for(i in 1:60){
  ## identify which row the positions of the selected particles are stored (because particles get deleted every time-step)
  pts.sel <- which(track.2D$id_list[[i]]%in%ind)
  ## simplify the list of cells into a matrix
  idx.m <- matrix(unlist(track.2D$idx_list_2D[[i]]), ncol=length(track.2D$idx_list_2D[[i]]))
  ## the cell indices for all selected particles at that time-step
  c.idx <- c(idx.m[pts.sel])
  ## plot
  points(Rdat$x[c.idx], Rdat$y[c.idx], cex=0.1, col=rainbow.cols[i])
}


###################################
#### for each cell, calculate the max distance travelled to get there
###################################
## need particle ID in each cell
## particle location at start
## straight line euclidian distance from start to finish

#### DETAILED TRAJECTORIES FOR 3D
# start.time <- Sys.time() #2.8h
# track.3D <- loopit_2D3D(pts = cell.pts, romsobject = Rdat, projected=TRUE, speed=100, domain="3D", runtime=30, roms_slices = 2, looping_time = 0.5, trajectories=TRUE, detailed_trajectories = TRUE)
# end.time <- Sys.time()
# end.time-start.time
#save(track.3D, file="/pvol/3_model_analysis/split_regions_runs_3Dtracking_30days_detailed.Rdata")
load("/pvol/3_model_analysis/split_regions_runs_3Dtracking_30days_detailed.Rdata")
Rdat$dist3D <- matrix(0, ncol=ncol(Rdat$x),nrow=nrow(Rdat$x))
Rdat$dist3Ds <- matrix(0, ncol=ncol(Rdat$x),nrow=nrow(Rdat$x))

## identify which row the positions of the particles are stored (particles get deleted at every time-step)
pts.sel.start <- track.3D$id_list[[1]]
## ptrack x and y positions
c.crds.start.xy <- cbind(track.3D$ptrack_x_list[[1]],track.3D$ptrack_y_list[[1]])
## starting cell indices
c.idx.start <- matrix(unlist(track.3D$idx_list_2D[[1]]), ncol=length(track.3D$idx_list_2D[[1]]))[,1]
## update values for cells if larger than previous value
for(i in 2:60){
  print(i)
  ## where are particles now
  pts.sel <- track.3D$id_list[[i]]
  ##
  c.crds.now.xy <- cbind(track.3D$ptrack_x_list[[i]],track.3D$ptrack_y_list[[i]])
  ## distance to where they started
  dist.dat <- rep(0, nrow(c.crds.now.xy))
  for(k in 1:nrow(c.crds.now.xy)){
    dist.dat[k] <- sqrt(sum((c.crds.now.xy[k, ] - c.crds.start.xy[pts.sel,][k, ])^2))
  }
  ## if larger than values already stored, then replace
  ## matrix of cell indices for currently observed particles
  idx.m <- matrix(unlist(track.3D$idx_list_2D[[i]]), ncol=length(track.3D$idx_list_2D[[i]]))
  c.idx <- idx.m[,24]
  repl.sel <- which(Rdat$dist3D[c.idx]<=dist.dat)
  Rdat$dist3D[c.idx][repl.sel] <- dist.dat[repl.sel]
 
  ## distance moved away from starting cell
  c.idx.start.update <- c.idx.start[pts.sel]
  repl.sel.start <- which(Rdat$dist3Ds[c.idx.start.update]<=dist.dat)
  Rdat$dist3Ds[c.idx.start.update][repl.sel.start] <- dist.dat[repl.sel.start]
}
dist3D.ra  <- setValues(empty.roms.ra, Rdat$dist3D)
dist3Ds.ra <- setValues(empty.roms.ra, Rdat$dist3Ds)

par(mfrow=c(1,2))
ylim=c(-2400000,-1000000)
plot(dist3D.ra, ylim=ylim, main="3D - distance travelled TO each cell")
contour(h.ra, add=TRUE, levels=c(0,1000,2000,3000))
plot(dist3Ds.ra, ylim=ylim, main="3D - distance travelled FROM each cell", maxcell=2000000)
contour(h.ra, add=TRUE, levels=c(0,1000,2000,3000))

#### DETAILED TRAJECTORIES FOR 2D
load("/pvol/3_model_analysis/split_regions_runs_2Dtracking_30days_detailed.Rdata")
Rdat$dist2D <- matrix(0, ncol=ncol(Rdat$x),nrow=nrow(Rdat$x))
Rdat$dist2Ds <- matrix(0, ncol=ncol(Rdat$x),nrow=nrow(Rdat$x))

## identify which row the positions of the particles are stored (particles get deleted at every time-step)
pts.sel.start <- track.2D.traj$id_list[[1]]
## ptrack x and y positions
c.crds.start.xy <- cbind(track.2D.traj$ptrack_x_list[[1]],track.2D.traj$ptrack_y_list[[1]])
## starting cell indices
c.idx.start <- matrix(unlist(track.2D.traj$idx_list_2D[[1]]), ncol=length(track.2D.traj$idx_list_2D[[1]]))[,1]
## update values for cells if larger than previous value
for(i in 2:60){
  print(i)
  ## where are particles now
  pts.sel <- track.2D.traj$id_list[[i]]
  ##
  c.crds.now.xy <- cbind(track.2D.traj$ptrack_x_list[[i]],track.2D.traj$ptrack_y_list[[i]])
  ## distance to where they started
  dist.dat <- rep(0, nrow(c.crds.now.xy))
  for(k in 1:nrow(c.crds.now.xy)){
    dist.dat[k] <- sqrt(sum((c.crds.now.xy[k, ] - c.crds.start.xy[pts.sel,][k, ])^2))
  }
  ## if larger than values already stored, then replace
  idx.m <- matrix(unlist(track.2D.traj$idx_list_2D[[i]]), ncol=length(track.2D.traj$idx_list_2D[[i]]))
  c.idx <- c(idx.m[pts.sel])
  repl.sel <- which(Rdat$dist2D[c.idx]<=dist.dat)
  Rdat$dist2D[c.idx][repl.sel] <- dist.dat[repl.sel]
  ## distance moved away from starting cell
  c.idx.start.update <- c.idx.start[pts.sel]
  repl.sel.start <- which(Rdat$dist2Ds[c.idx.start.update]<=dist.dat)
  Rdat$dist2Ds[c.idx.start.update][repl.sel.start] <- dist.dat[repl.sel.start]
  
}
dist2D.ra  <- setValues(empty.roms.ra, Rdat$dist2Ds)
dist2Ds.ra <- setValues(empty.roms.ra, Rdat$dist2Ds)

par(mfrow=c(1,2))
plot(dist2D.ra)
contour(h.ra, add=TRUE, levels=c(0,1000,2000,3000))
plot(dist2Ds.ra)
contour(h.ra, add=TRUE, levels=c(0,1000,2000,3000))























##############################################################################

# load("/pvol/3_model_analysis/split_regions_runs_3Dtracking_30days.Rdata")
# load("/pvol/3_model_analysis/split_regions_runs_2Dtracking_30days.Rdata")

'%!in%' <- function(x,y)!('%in%'(x,y))
library(geosphere)
library(sf)

### assign 0 to each cell in the domain

Rdat$dist <- matrix(0, ncol=ncol(Rdat$x),nrow=nrow(Rdat$x))

## identify which row the positions of the particles are stored (particles get deleted at every time-step)
pts.sel.start <- track.3D$id_list[[1]]
## simplify the list of cells into a matrix
idx.m.start <- matrix(unlist(track.3D$idx_list_2D[[1]]), ncol=length(track.3D$idx_list_2D[[1]]))
## the cell indices for all particles at that time-step
c.idx.start <- c(idx.m.start[,1])
## coordinates for those cells
c.crds.start.xy <- cbind(Rdat$x[c.idx.start], Rdat$y[c.idx.start])

## update values for cells if larger than previous value
for(i in 2){
  print(i)
  ## where are particles now
  pts.sel <- track.3D$id_list[[i]]
  #idx.m <- matrix(unlist(track.3D$idx_list_2D[[i]]), ncol=length(track.3D$idx_list_2D[[i]]))
  idx.m <- matrix(unlist(track.3D$idx_list[[i]]), ncol=length(track.3D$idx_list[[i]]))
  c.idx <- c(idx.m[,24])
  c.crds.now.xy <- cbind(Rdat$x[c.idx], Rdat$y[c.idx])
  ## distance to where they started
  dist.dat <- rep(0, nrow(c.crds.now.xy))
  for(k in 1:nrow(c.crds.now.xy)){
    dist.dat[k] <- sqrt(sum((c.crds.now.xy[k, ] - c.crds.start.xy[pts.sel,][k, ])^2))
  }
  ## if larger than values already stored, then replace
  repl.sel <- which(Rdat$dist[c.idx]<=dist.dat)
  Rdat$dist[c.idx][repl.sel] <- dist.dat[repl.sel]
}
dist.ra       <- setValues(empty.roms.ra, Rdat$dist)
#### SOMETHING IS WEIRD, DISTANCES ARE VEEEEERYYYYY LOOOONG

ind <- c(74,100,1145)
ind <- sample(seq(1,nrow(cell.pts)), 10000)
rainbow.cols <- rainbow(60)
plot(uv.ra, main="3D results")
for(i in 1:60){
  ## identify which row the positions of the selected particles are stored (because particles get deleted every time-step)
  pts.sel <- which(track.3D$id_list[[i]]%in%ind)
  ## simplify the list of cells into a matrix
  idx.m <- matrix(unlist(track.3D$idx_list_2D[[i]]), ncol=length(track.3D$idx_list_2D[[i]]))
  ## the cell indices for all selected particles at that time-step
  c.idx <- c(idx.m[pts.sel])
  ## plot
  points(Rdat$x[c.idx], Rdat$y[c.idx], cex=0.1, col=rainbow.cols[i])
}


#### TRACKIT_3D ONLY WORKS WELL
### 
Rdat$dist <- matrix(0, ncol=ncol(Rdat$x),nrow=nrow(Rdat$x))
dist.dat <- rep(0, nrow(track.3D.s$ptrack))
for(k in 1:nrow(track.3D.s$ptrack)){
  dist.dat[k] <- sqrt(sum((track.3D.s$ptrack[k,1:2,1] - track.3D.s$ptrack[k,1:2,track.3D.s$stopindex[k]])^2))
}
plot(uv.surf.ra)
sel <- which(dist.dat>100000)
for(i in sel){
  points(t(track.3D.s$ptrack[i,1:2,1:track.3D.s$stopindex[i]]), col="blue", cex=0.1)
}
points(track.3D.s$ptrack[sel,1:2,1], cex=0.1)
sel <- which(dist.dat>200000)
for(i in sel){
  points(t(track.3D.s$ptrack[i,1:2,1:track.3D.s$stopindex[i]]), col="red", cex=0.1)
}
points(track.3D.s$ptrack[sel,1:2,1], cex=0.1)


























### assign 0 to each cell in the domain
dist.ra <- uv.ra
dist.ra[!is.na(uv.ra[])] <- 0

## identify which row the positions of the particles are stored (particles get deleted at every time-step)
pts.sel.start <- track.3D$id_list[[1]]
## simplify the list of cells into a matrix
idx.m.start <- matrix(unlist(track.3D$idx_list_2D[[1]]), ncol=length(track.3D$idx_list_2D[[1]]))
## the cell indices for all particles at that time-step
c.idx.start <- c(idx.m.start[,1])
## coordinates for those cells
c.crds.start.xy <- cbind(Rdat$x[c.idx.start], Rdat$y[c.idx.start])
# c.crds.start.xy <- xyFromCell(uv.ra, c.idx.start)
# c.crds.start.xy <- st_as_sf(data.frame(xyFromCell(uv.ra, c.idx.start)), coords=c("x","y"), crs=crs(dist.ra))
# c.crds.start.lonlat <- st_transform(c.crds.start.xy, crs="+proj=longlat")

## update values for cells if larger than previous value
for(i in 2:60){
  pts.sel <- track.3D$id_list[[i]]
  idx.m <- matrix(unlist(track.3D$idx_list_2D[[i]]), ncol=length(track.3D$idx_list_2D[[i]]))
  c.idx <- c(idx.m[pts.sel])
  ## find which c.idx relates to which c.idx.start
  #pts.sel.start[which(pts.sel.start%in%pts.sel)]
  ## distance to where they started
  c.crds.now.xy <- cbind(Rdat$x[c.idx], Rdat$y[c.idx])
  # c.crds.now.xy <- st_as_sf(data.frame(xyFromCell(uv.ra, c.idx)), coords=c("x","y"), crs=crs(uv.ra))
  # c.crds.now.lonlat <- st_transform(c.crds.now.xy, crs="+proj=longlat")
  dist.dat <- NA
  for(k in 1:nrow(c.crds.now.xy)){
    dist.dat[k] <- sqrt(sum((c.crds.now.xy[k, ] - c.crds.start.xy[pts.sel,][k, ])^2))
  }
  ## if larger than values already stored, then replace
  dist.ra[]
  
}





euclidean_distance <- function(row) {
  sqrt(sum((m.now[row, ] - m.start[row, ])^2))
}
# Apply the function over the rows
distances <- apply(c.crds.now.xy, 1, euclidean_distance)




c.crds.start.xy <- st_as_sf(data.frame(xyFromCell(uv.ra, c.idx.start)), coords=c("x","y"), crs=crs(uv.ra))




xys <- terra::vect(c.crds.start, geom=c('Lon','Lat'),crs="+proj=longlat")


lcc <- terra::project(latlons,"+proj=lcc +lat_0=38.5 +lon_0=262.5 +lat_1=38.5 +lat_2=38.5 +x_0=0 +y_0=0 +R=6371229 +units=m +no_defs")





ind <- sample(seq(1,nrow(cell.pts)), 1000)

### from LOOPIT_2D3D
par(mfrow=c(1,2))
rainbow.cols <- rainbow(60)
plot(uv.ra, main="3D results")
for(i in 1:60){
  ## identify which row the positions of the selected particles are stored (because particles get deleted every time-step)
  pts.sel <- which(track.3D$id_list[[i]]%in%ind)
  ## simplify the list of cells into a matrix
  idx.m <- matrix(unlist(track.3D$idx_list_2D[[i]]), ncol=length(track.3D$idx_list_2D[[i]]))
  ## the cell indices for all selected particles at that time-step
  c.idx <- c(idx.m[pts.sel])
  ## plot
  points(Rdat$x[c.idx], Rdat$y[c.idx], cex=0.1, col=rainbow.cols[i])
}
plot(uv.ra, main="2D results")
for(i in 1:60){
  ## identify which row the positions of the selected particles are stored (because particles get deleted every time-step)
  pts.sel <- which(track.2D$id_list[[i]]%in%ind)
  ## simplify the list of cells into a matrix
  idx.m <- matrix(unlist(track.2D$idx_list_2D[[i]]), ncol=length(track.2D$idx_list_2D[[i]]))
  ## the cell indices for all selected particles at that time-step
  c.idx <- c(idx.m[pts.sel])
  ## plot
  points(Rdat$x[c.idx], Rdat$y[c.idx], cex=0.1, col=rainbow.cols[i])
}


















#####################################################################################
######################################################################################
#### 2D particle tracks on 3D bathymetry ####
library(rayshader)

prepare_overlays <- function(height.dat, topo.dat){
  topo.ext <- raster::extent(as.vector(ext(topo.dat)))
  
  ## change data to matrix
  topo_mat=raster_to_matrix(topo.dat)
  height_mat=raster_to_matrix(height.dat)
  
  ## calculate sphereshades and lineoverlays
  sphereshades  = sphere_shade( topo_mat, texture="desert", colorintensity=3)
  rayshades     = ray_shade(    topo_mat, sunaltitude=6, zscale=zscale, lambert=FALSE)
  lambshades    = lamb_shade(   topo_mat, sunaltitude=6, zscale=zscale)
  ambientshades = ambient_shade(topo_mat)
  textureshades = texture_shade(topo_mat, detail=6/10, contrast=10, brightness = 11)
  
  heightshades  = height_shade(height_mat, texture=rev(heat.colors(256)))
#  coastoverlay = generate_line_overlay(sf_coast, extent=topo.ext, heightmap=topo_mat)
#  lineoverlay  = generate_line_overlay(sf_lines, extent=topo.ext, heightmap=topo_mat, color="black", offset=c(1500,-1500))
  
#  SF_overlay_shadow = generate_point_overlay(sf_points, extent=topo.ext, heightmap=topo_mat, color="black", offset=c(1500,-1500), size=5, pch=21)
#  SF_overlay = generate_point_overlay(sf_points, extent=topo.ext, heightmap=topo_mat, color=SF.dat$coloursSF, size=5, pch=21)
  
  return(list(topo_mat = topo_mat,
              sphereshades = sphereshades,
              rayshades = rayshades,
              lambshades = lambshades,
              ambientshades = ambientshades,
              textureshades = textureshades,
              heightshades = heightshades#,
              #coastoverlay = coastoverlay,
              #lineoverlay = lineoverlay,
              #SF_overlay_shadow = SF_overlay_shadow,
              #SF_overlay = SF_overlay
  ))
}

## make raster suitable for plotting
topo.dat <- aggregate(h.ra,2)
curr.dat <- resample(uv.ra, h.ra)
#rich.dat <- resample(rich.dat.ross, topo.dat.ross)

## parameters for 3d-plotting
zscale <- 50 ## the higher the less exaggerated the height values are
phi <- 40    ## vertical angle of view
theta <- 150 ## horizontal angle of view
zoom <- 0.7  ## zoom between 0 (close) and 1 (far)
windowsize <- c(1000,1000)
rs <- prepare_overlays(curr.dat, topo.dat)

rs$topo_mat %>%
  sphere_shade(texture = "bw") %>%
  add_overlay(rs$heightshades, alphalayer = 0.6) %>%
  #add_overlay(rs$coastoverlay) %>% 
  #add_overlay(rs$lineoverlay) %>%   
  #add_overlay(rs$SF_overlay_shadow) %>% 
  #add_overlay(rs$SF_overlay) %>% 
  plot_3d(rs$topo_mat, zscale=zscale, phi=phi, theta=theta, zoom=zoom, windowsize=windowsize)
render_snapshot()












## 3D output plot

#' ra <- raster(nrow=50,ncol=50,ext=extent(surface_chl))
#' r_roms <- rasterize(x = cbind(as.vector(toyROMS$lon_u), as.vector(toyROMS$lat_u)), y= ra, field = as.vector(-toyROMS$h))
#' pr <- projectRaster(r_roms, crs = "+proj=laea +lon_0=137 +lat_0=-66")  #get the right projection (through the centre)
#' 
#' plot3D(pr, adjust = FALSE, zfac = 50)                    # plot bathymetry with 50x exaggerated depth
#' pointsxy <- project(as.matrix(run$pend[,1:2]), projection(pr))  #projection on Tracking-points
#' points3d(pointsxy[,1],pointsxy[,2],run$pend[,3]*50)#,xlim=xlim,ylim=ylim)



## 2D output plot
#' ########## 2D-tracking with storing trajectories:
#' pts_seeded <- create_points_pattern(surface_chl, multi=100)
#' run <- loopit_2D3D(pts_seeded = pts_seeded, roms_slices = 4, particle_radius = 0.00001, romsobject = toyROMS, speed = 100, runtime = 50, sedimentation = TRUE, trajectories = TRUE)
#' 
#' plot(pts_seeded)
#' points(run$pend, col="red", cex=0.6)
#' points(run$pts , col="blue", cex=0.6)
#' 
#' ## looking at the horizontal flux: this should be another function to handle the output
#' ra <- raster(nrow=50,ncol=50,ext=extent(surface_chl))
#' mat_list <- list()
#' for(islices in 1:length(run$idx_list_2D)){
#'   mat_list[[islices]] <- matrix(unlist(run$idx_list_2D[[islices]]),ncol=12)
#' }
#' testmatrix <- do.call(rbind, mat_list)
#' testid <- unlist(run$id_list)
#' flux_list <- split(testmatrix,testid)
#' for(k in 1:length(flux_list)){
#'   ## cells visited by a particle ("presence-only")
#'   flux_list[[k]] <- unique(flux_list[[k]])
#'   ## drop first and last value (input and setting cell)
#'   flux_list[[k]] <- flux_list[[k]][-c(1,length(flux_list[[k]]))]
#' } 
#' flux <- as.vector(unlist(flux_list))
#' 
#' xlim <- c(xmin(ra),xmax(ra))
#' ylim <- c(ymin(ra),ymax(ra))
#' df <- data.frame(cbind(toyROMS$lon_u[flux],toyROMS$lat_u[flux]))
#' df$cell <- cellFromXY(ra, df)
#' ra[] <- tabulate(df$cell, ncell(ra))
#' plot(ra)
#' 
#' ## looking at the current-slices
#' roms_list <- list()
#' par(mfrow=c(2,2))
#' for(i in 1:4){
#'   roms_list[[i]] <- rasterize(x = cbind(as.vector(toyROMS$lon_u), as.vector(toyROMS$lat_u)), y= ra, field = as.vector(sqrt((toyROMS$i_u[,,,i]^2)+(toyROMS$i_v[,,,i])^2)))
#'   plot(roms_list[[i]])
#' }
#' par(mfrow=c(1,1))













#########################################
########################################
#########################################
# ## reduce the ROMS by excluding entire areas
# load(paste0(env.dir,"ocean_his_0001_2slices_TrackingSetup_transposed.Rdata"))
# row.lim <- 1701:2400
# col.lim <- 1001:1400
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
# save(Rdat, file=paste0(env.dir,"ocean_his_0001_2slices_TrackingSetup_transposed_RossSeaE.Rdata"))
load(paste0(env.dir,"ocean_his_0001_2slices_TrackingSetup_transposed_RossSeaE.Rdata"))

###############################
###### create uv and h raster for plotting
roms.coords.proj <- cbind(c(Rdat$x), c(Rdat$y))
x.range <- c(min(roms.coords.proj[,1])-1000,max(roms.coords.proj[,1])+1000)
y.range <- c(min(roms.coords.proj[,2])-1000,max(roms.coords.proj[,2])+1000)
empty.roms.ra <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000)

uv.ra      <- setValues(empty.roms.ra, sqrt(Rdat$i_u[,,1,1]^2 +c(Rdat$i_v[,,1,1])^2))
uv.surf.ra <- setValues(empty.roms.ra, sqrt(Rdat$i_u[,,31,1]^2+c(Rdat$i_v[,,31,1])^2))
h.ra       <- setValues(empty.roms.ra, Rdat$h)

###################################
#### visually check model setup
## 3D current field:
plot(uv.ra)
# plot(coast.proj, add=TRUE)

#### One particle per cell with npp value in it ###############
## seeded with particles based on NPP values
npp <- rast("/pvol/data_environmental/Circumpolar_EnvData_2km_NPP_Cafe_filled_SummerAverage.tif")
npp2 <- project(npp, h.ra)
npp2[npp2[]>0] <- 1
npp2[400:700,1:200] <- NA
## thin out to one particle every 4 cells
npp2[seq(1,700, by=2),] <- NA
npp2[,seq(1,400, by=2)] <- NA
##
cell.pts <- cbind(crds(as.points(npp2)),0)


###############################
#### 3D tracking ##############
###############################

start.time <- Sys.time() #1.5h
track.3D <- loopit_2D3D(pts = cell.pts, romsobject = Rdat, projected=TRUE, speed=100, domain="3D", runtime=30, roms_slices = 2, looping_time = 0.5, trajectories=TRUE)
end.time <- Sys.time()
end.time-start.time
#save(track.3D, file="/pvol/3_model_analysis/split_regions2_runs_3Dtracking_30days.Rdata")
# load("/pvol/3_model_analysis/split_regions2_runs_3Dtracking_30days.Rdata")

#### 3D: extract results ####

### from LOOPIT_2D3D
## all points that haven't stopped have stopindex 0
## "pts" and "pnow" have some that have settled in the last step (calculated below as "pfloat")
not.stopped <- which(track.3D$stopindex==0)
pfloat <- matrix(track.3D$ptrack[not.stopped,,dim(track.3D$ptrack)[3]],ncol=3)
## the ones that stopped are "pend"
pstop <- track.3D$pend

#### plot 3D results ####

par(mfrow=c(1,2))

plot(uv.ra, main="loopit_2D3D results")
points(cell.pts, cex=0.1, pch=16)
#points(track.3D$pend[,1:2], col="blue", cex=0.1)
points(pstop, col="blue", cex=0.1)
points(pfloat, col="purple", cex=0.1)

#### 3D particle tracks ####

ind <- sample(seq(1,nrow(cell.pts)), 1000)

### from LOOPIT_2D3D

rainbow.cols <- rainbow(60)
plot(uv.ra, main="loopit_2D3D results")
for(i in 1:60){
  ## identify which row the positions of the selected particles are stored (because particles get deleted every time-step)
  pts.sel <- which(track.3D$id_list[[i]]%in%ind)
  ## simplify the list of cells into a matrix
  idx.m <- matrix(unlist(track.3D$idx_list_2D[[i]]), ncol=length(track.3D$idx_list_2D[[i]]))
  ## the cell indices for all selected particles at that time-step
  c.idx <- c(idx.m[pts.sel])
  ## plot
  points(Rdat$x[c.idx], Rdat$y[c.idx], cex=0.1, col=rainbow.cols[i])
}

###############################
#### 2D tracking ##############
###############################

pts.3D <- rbind(pstop, pfloat) ## use the loopit output

###############################
#### 2D tracking #### 
Rdat.l <- Rdat
Rdat.l$hh <- Rdat$hh[,,1:8]
Rdat.l$i_u <- Rdat$i_u[,,1:8,]
Rdat.l$i_v <- Rdat$i_v[,,1:8,]
Rdat.l$i_w <- Rdat$i_w[,,1:8,]

start.time <- Sys.time() ## 3.4min
track.2D <-loopit_2D3D(pts = pts.3D, romsobject = Rdat.l, projected=TRUE, speed=100, domain="2D", runtime=30, looping_time = 0.5, roms_slices = 2, trajectories = TRUE)
end.time <- Sys.time()
end.time-start.time
#save(track.2D, file="/pvol/3_model_analysis/split_regions2_runs_2Dtracking_30days.Rdata")

### from LOOPIT_2D3D
## all points that haven't stopped have stopindex 0
## "pts" and "pnow" have some that have settled in the last step (calculated below as "pfloat")
not.stopped <- which(track.2D$stopindex==0)
pfloat <- matrix(track.2D$ptrack[not.stopped,,dim(track.2D$ptrack)[3]],ncol=3)
## the ones that stopped are "pend"
pstop <- track.2D$pend

#### plot 2D results ####
par(mfrow=c(1,2))

plot(uv.ra, main="loopit_2D3D results")
points(cell.pts, cex=0.1, pch=16)
points(pstop, col="blue", cex=0.1)
points(pfloat, col="purple", cex=0.1)

#### 2D particle tracks ####

ind <- sample(seq(1,nrow(cell.pts)), 500)

rainbow.cols <- rainbow(60)
plot(uv.ra)#, xlim=xlim.p, ylim=ylim.p, main="loopit_2D3D results")
for(i in 1:60){
  ## identify which row the positions of the selected particles are stored (because particles get deleted every time-step)
  pts.sel <- which(track.2D$id_list[[i]]%in%ind)
  ## simplify the list of cells into a matrix
  idx.m <- matrix(unlist(track.2D$idx_list_2D[[i]]), ncol=length(track.2D$idx_list_2D[[i]]))
  ## the cell indices for all selected particles at that time-step
  c.idx <- c(idx.m[pts.sel])
  ## plot
  points(Rdat$x[c.idx], Rdat$y[c.idx], cex=0.1, col=rainbow.cols[i])
}
