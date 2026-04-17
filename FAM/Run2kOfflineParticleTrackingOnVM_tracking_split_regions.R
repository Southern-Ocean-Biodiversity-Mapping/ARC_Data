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
env.dir3 <- "/pvol3TB/data_environmental/ROMS_2k_files/"

string.chr <- "Circumpolar_EnvData_"
string.res <- "500m_"

## polar stereographic projection:
stereo <- "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

'%!in%' <- function(x,y)!('%in%'(x,y))

################################
###### particle tracking code
#### load functions
#ptrackr_dir <- "/pvol/Scripts/ptrackR_package_update2024/ptrackr/R/"
ptrackr_dir <- "/pvol3TB/Scripts/ptrackR_package_update2024/ptrackr/R/"

## minimum for 3d run
source(paste0(ptrackr_dir,"create_points_pattern.R"))
source(paste0(ptrackr_dir,"setup_knn.R"))
source(paste0(ptrackr_dir,"trackit_3D.R"))
## for 2d run:
source(paste0(ptrackr_dir,"trackit_2D.R"))
source(paste0(ptrackr_dir,"buildparams.r"))
## to loop:
#source(paste0(ptrackr_dir,"loopit_2D3D.R"))
source(paste0(ptrackr_dir,"loopit_2D3D_traj_ROMSreduced.R"))

################################
#### 3D tracking looping through regions 1:10
speed=200


############################
## QUICK TESTING CODE:
i=1
message(i)
## load ROMS file
load(paste0(env.dir3,"ocean_his_TrackingSetup_6hourlycurrents_28days_Region",sprintf("%02d",i),".Rdata"))
str(Rdat6h)
## load NPP file
load(paste0(env.dir3,"ocean_his_TrackingSetup_eppl_12boxfilled_NPP9_Region",sprintf("%02d",i),".Rdata"))
str(pts9)
pts1 <- pts9[seq(1,nrow(pts9), by=9),]
pts001 <- pts1[seq(1,nrow(pts1), by=100),]


#### run 3D tracking
start.time <- Sys.time()
track.3D.red <- loopit_2D3D(pts=as.matrix(pts1[,1:3]), romsobject=Rdat6h, projected=TRUE, speed=speed, domain="3D", runtime=1, roms_slices=30, looping_time=0.25, ROMSreduced = TRUE)
print(Sys.time()-start.time)
start.time <- Sys.time()
track.3D.notred <- loopit_2D3D(pts=as.matrix(pts1[,1:3]), romsobject=Rdat6h, projected=TRUE, speed=speed, domain="3D", runtime=1, roms_slices=30, looping_time=0.25, ROMSreduced = FALSE)
# track.3D <- loopit_2D3D(pts=as.matrix(pts001[,1:3]), romsobject=Rdat6h, projected=TRUE, speed=speed, domain="3D", runtime=1, roms_slices=113, looping_time=0.25, trajectories=TRUE, detailed_trajectories=TRUE)
print(Sys.time()-start.time)

###########################
#### now back to the full run
## this code follows particles as they sink (with a given speed) through the water column and get advected by currents
npp.model <- "eppl_12boxfilled"
for(i in 1:10){
  message(i)
  ## load ocean currents file (here currents for every 6h)
  load(paste0(env.dir3,"ocean_his_TrackingSetup_6hourlycurrents_28days_Region",sprintf("%02d",i),".Rdata"))
  ## load particle positions (here 9 particles in every cell, each will later be given a value representing 1/9 of the total NPP from that cell)
  load(paste0(env.dir3,"ocean_his_TrackingSetup_",npp.model,"_NPP9_Region",sprintf("%02d",i),".Rdata"))
  
  #### run 3D tracking
  start.time <- Sys.time()
  #track.3D <- loopit_2D3D(pts=as.matrix(pts4[,1:3]), romsobject=Rdat, projected=TRUE, speed=speed, domain="3D", runtime=1, roms_slices=113, looping_time=0.25, trajectories=TRUE, detailed_trajectories=TRUE)
  #track.3D <- loopit_2D3D(pts=as.matrix(pts9[,1:3]), romsobject=Rdat6h, projected=TRUE, speed=speed, domain="3D", runtime=1, roms_slices=113, looping_time=0.25, trajectories=TRUE, detailed_trajectories=TRUE)
  track.3D <- loopit_2D3D(pts=as.matrix(pts9[,1:3]), romsobject=Rdat6h, projected=TRUE, speed=speed, domain="3D", runtime=1, roms_slices=84, looping_time=0.25, ROMSreduced=TRUE)
  print(Sys.time()-start.time)
  #save(track.3D, file=paste0("/pvol3TB/FAM_outputs/tracking3D_Region",sprintf("%02d",i),"_NPP9_",speed,"mday_28days_trajectoriesdetailed.Rdata"))
  save(track.3D, file=paste0("/pvol3TB/FAM_outputs/tracking3D_Region",sprintf("%02d",i),"_",npp.model,"_NPP9_",speed,"mday_21days_red.Rdata"))
}

#### 3D translating looping through regions 1:10
## this is to 
for(i in 1:10){
  message(i)
  ## load input NPP9 particle locations
  load(paste0(env.dir3,"ocean_his_TrackingSetup_NPP9_Region",sprintf("%02d",i),".Rdata"))
  ## load original NPP raster from that region, to be used as reference for creating a new raster
  npp_crop <- rast(paste0(env.dir3,"ocean_his_TrackingSetup_NPP_Region",sprintf("%02d",i),".tif"))
  # ## end particle locations need to be the same format as input locations
  # xy.end <- pts9
  # ## set the x & y positions to NA, then replace with end locations
  # xy.end[,1:2] <- NA
  ## load tracking run:
  load(paste0("/pvol3TB/FAM_outputs/tracking3D_Region",sprintf("%02d",i),"_NPP9_",speed,"mday_21days_red.Rdata"))
  ## add NPP to xyz_end
  xyz.end <- track.3D$xyz_end
  xyz.end$npp <- pts9[,4] ##$npp
  # ## fill x and y end positions from the back (stopped ones get removed, so we want to have their last record before removal)
  # pts.sel <- track.3D$id_list[[113]]
  # xy.end$x[pts.sel] <- track.3D$ptrack_x_list[[113]]
  # xy.end$y[pts.sel] <- track.3D$ptrack_y_list[[113]]
  # for(k in 112:1){
  #   print(k)
  #   ## the original row indices of the points at each time-step
  #   pts.sel <- track.3D$id_list[[k]]
  #   ## we only want to replace those values that are NA (haven't already received an x and y position)
  #   na.sel <- which(is.na(xy.end$x[pts.sel]))
  #   ## now replace NAs
  #   xy.end$x[pts.sel[na.sel]] <- track.3D$ptrack_x_list[[k]][na.sel]
  #   xy.end$y[pts.sel[na.sel]] <- track.3D$ptrack_y_list[[k]][na.sel]
  # }
  ## add cellID to each point
  xyz.end$cellID <- extract(npp_crop, xyz.end[,1:2], cells=TRUE)$cell
  ## sum up npp for each cellID
  cell.npp <- aggregate(xyz.end$npp, by=list(cellID=xyz.end$cellID), FUN=sum)
  ## add these values to the raster
  ra.3d <- rast(npp_crop)
  ra.3d[cell.npp$cellID] <- cell.npp$x
  ## save output
  save(xyz.end, cell.npp, file=paste0("/pvol3TB/FAM_outputs/tracking3D_Region",sprintf("%02d",i),"_NPP9_",speed,"mday_21days_xyzbottomnpp.Rdata"))
  writeRaster(ra.3d, file=paste0("/pvol3TB/FAM_outputs/tracking3D_Region",sprintf("%02d",i),"_NPP9_",speed,"mday_21days_bottomnpp.tif"))
}

#### 2D looping through regions 1:10
speed <- 200
pr <- 0.0001
r <- "0001"
traj <- TRUE
if(traj==TRUE){t <- "traj_"} else {t <- ""}

for(i in 10){
  message(i)
  ## why this?
  # load("/pvol3TB/FAM_outputs/tracking3D_Region",sprintf("%02d",i),"_NPP9_",speed,"mday_21days_xyzbottomnpp.Rdata")

  ## load bottom npp distribution
  ra9 <- rast(paste0("/pvol3TB/FAM_outputs/tracking3D_Region",sprintf("%02d",i),"_NPP9_",speed,"mday_21days_bottomnpp.tif"))
  
  ## translate bottom npp distribution to particle distribution
  Z <- terra2im(ra9/10000000) ## need to reduce by at least 10m
  #Z <- terra2im(ra9/500000000) ## TEST
  
  pp <- rpoispp(Z) ## library(spatstat)
  npp.pts <- cbind(pp$x, pp$y,0)
  ## transform pts back to raster, to check how well the points represent the original data
  ra.pts <- rasterize(npp.pts[,1:2], ra9, fun=sum)
  #plot(ra.pts[],ra9[], cex=0.1)

  ## load Rdat file
  load(paste0("/pvol3TB/data_environmental/ROMS_2k_files/ocean_his_TrackingSetup_6hourlycurrents_28days_Region",sprintf("%02d",i),".Rdata"))
  
  ## set z-values to depth of cells
  roms.coords.proj <- cbind(c(Rdat6h$x), c(Rdat6h$y))
  x.range <- c(min(roms.coords.proj[,1])-1000,max(roms.coords.proj[,1])+1000)
  y.range <- c(min(roms.coords.proj[,2])-1000,max(roms.coords.proj[,2])+1000)
  empty.roms.ra <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000)
  h.ra       <- setValues(empty.roms.ra, Rdat6h$h)
  npp.pts.h <- extract(h.ra, npp.pts[,1:2])
  npp.pts[,3] <- npp.pts.h$lyr.1
  
  ## reduce size of Rdat
  Rdat6h$hh <- Rdat6h$hh [,,-c(4:31)]
  Rdat6h$i_u<- Rdat6h$i_u[,,-c(4:31),]
  Rdat6h$i_v<- Rdat6h$i_v[,,-c(4:31),]
  Rdat6h$i_w<- Rdat6h$i_w[,,-c(4:31),]
  
  start.time <- Sys.time()
  track.2D <-loopit_2D3D(pts=npp.pts, romsobject=Rdat6h, projected=TRUE, sedimentation=TRUE, speed=speed, particle_radius=pr, domain="2D", runtime=1, looping_time=0.25, roms_slices=113, trajectories=traj)
  end.time <- Sys.time()
  end.time-start.time
  #save(track.2D, file=paste0("/pvol3TB/FAM_outputs/tracking2D_Region_TEST.Rdata"))
  ## no run before 11/2024 has npp.pts saved...
  save(track.2D, npp.pts, file=paste0("/pvol3TB/FAM_outputs/tracking2D_Region",sprintf("%02d",i),"_NPP9_",speed,"mday_21days_",t,"r",r,"_28days_Region",sprintf("%02d",i),".Rdata"))
}

#### 2D detailed trajectories for 1 particle per cell
speed <- 200
pr <- 0.0001
r <- "0001"
traj <- TRUE
if(traj==TRUE){t <- "traj_"} else {t <- ""}
for(i in 1:4){
  message(i)
  ## load npp distribution
  npp_crop <- rast(paste0(env.dir3,"ocean_his_TrackingSetup_NPP_Region",sprintf("%02d",i),".tif"))
  ## translate bottom npp distribution to 1 particle per cell
  npp_crop[npp_crop>0] <- 1
  ## and extract the positions:
  npp.pts <- cbind(crds(npp_crop),0)
  ## load Rdat file
  load(paste0("/pvol3TB/data_environmental/ROMS_2k_files/ocean_his_TrackingSetup_6hourlycurrents_28days_Region",sprintf("%02d",i),".Rdata"))
  # ## set z-values to depth of cells
  # roms.coords.proj <- cbind(c(Rdat6h$x), c(Rdat6h$y))
  # x.range <- c(min(roms.coords.proj[,1])-1000,max(roms.coords.proj[,1])+1000)
  # y.range <- c(min(roms.coords.proj[,2])-1000,max(roms.coords.proj[,2])+1000)
  # empty.roms.ra <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000)
  # h.ra       <- setValues(empty.roms.ra, Rdat6h$h)
  # npp.pts.h <- extract(h.ra, npp.pts[,1:2])
  # npp.pts[,3] <- npp.pts.h$lyr.1
  ## reduce size of Rdat
  Rdat6h$hh <- Rdat6h$hh [,,-c(4:31)]
  Rdat6h$i_u<- Rdat6h$i_u[,,-c(4:31),]
  Rdat6h$i_v<- Rdat6h$i_v[,,-c(4:31),]
  Rdat6h$i_w<- Rdat6h$i_w[,,-c(4:31),]
  
  start.time <- Sys.time()
  track.2D <-loopit_2D3D(pts = npp.pts, romsobject = Rdat6h, projected=TRUE, sedimentation=FALSE, speed=speed, particle_radius=pr, domain="2D", runtime=1, looping_time=0.25, roms_slices = 113, trajectories=TRUE, detailed_trajectories=TRUE)
  end.time <- Sys.time()
  end.time-start.time
  save(track.2D, file=paste0("/pvol3TB/FAM_outputs/tracking2D_Region",sprintf("%02d",i),"_NPP1_",speed,"mday_21days_",t,"28days.Rdata"))
}


# pts_seeded=npp.pts[sample(1:nrow(npp.pts),1000000),]
# romsobject=Rdat6h
# projected=TRUE
# sedimentation=TRUE
# speed=speed
# particle_radius=pr
# domain="2D"
# runtime=1
# looping_time=0.25
# roms_slices=113
# trajectories=traj
# start_slice=1
# time_steps_in_s=1800
# uphill_restricted=NULL
# sed_at_max_speed=FALSE
# mean_move=FALSE
# projected=TRUE
# detailed_trajectories = FALSE
# ROMSreduced=FALSE
# irun=1


#### 2D analysis of NPP9, creating tif files of sedimentation
speed <- 200
r <- "0001"
for(i in c(1,5:10)){
  ## load particle run file
  load(paste0("/pvol3TB/FAM_outputs/tracking2D_Region",sprintf("%02d",i),"_NPP9_",speed,"mday_21days_",t,"r",r,"_28days_Region",sprintf("%02d",i),".Rdata"))
  ## load raster specs
  npp_crop <- rast(paste0(env.dir3,"ocean_his_TrackingSetup_NPP_Region",sprintf("%02d",i),".tif"))
  ## add cellID to each point
  xy.end.2D <- data.frame(track.2D$pend)
  cellIDs.df <- extract(npp_crop, xy.end.2D[,1:2], cells=TRUE)
  xy.end.2D$cellID <-cellIDs.df$cell
  ## sum up points for each cellID
  library(dplyr)
  dat <- xy.end.2D %>%
    group_by(cellID) %>%
    summarise(freq = n())
  ## add these values to the raster
  ra.2d <- rast(npp_crop)
  sel.na <- which(is.na(dat$cellID))
  ra.2d[dat$cellID[-sel.na]] <- dat$freq[-sel.na]
  # ra.2d.adj <- ra.2d
  # ra.2d.adj[ra.2d[]>2000] <- 2000
  ## save sedimentation output
  writeRaster(ra.2d, file=paste0("/pvol3TB/FAM_outputs/tracking2D_Region",sprintf("%02d",i),"_NPP9_",speed,"mday_21days_",t,"r",r,"_28days_Region",sprintf("%02d",i),"_sed.tif"))
}

#### 2D analysis of NPP9, creating tif files of flux
## cells visited at each time-step (6hourly)
for(i in c(5)){
  ## load particle tracking run (takes 5minutes)
  load(paste0("/pvol3TB/FAM_outputs/tracking2D_Region",sprintf("%02d",i),"_NPP9_",speed,"mday_21days_",t,"r",r,"_28days_Region",sprintf("%02d",i),".Rdata"))
  ## summing frequencies up in one go breaks the machine. For example, a tibble can handle 2 billion counts, but Ross Sea has 5 billion
  #knnIDs.v <- unlist(track.2D$idx_list_2D)
  #knnIDs <- tibble(knnID=knnIDs.v)
  ## looping into smaller bits instead
  knnIDs_tibbles <- list()
  for(k in 1:length(track.2D$idx_list_2D)){
    print(k)
    #tib <- data.frame("knnID"=unlist(track.2D$idx_list_2D[[k]]))
    tib <- tibble(knnID=unlist(track.2D$idx_list_2D[[k]]))
    knnIDs_tibbles[[k]] <- tib %>%
      count(knnID, name = "freq")
  }
  dat.flux <- bind_rows(knnIDs_tibbles) %>%
    group_by(knnID) %>%
    summarise(freq = sum(freq))
  #### need to translate knn IDs to raster cell IDs:
  rm(track.2D)
  ## load Rdat file (takes 5 minutes) and add these values to the correct raster
  load(paste0("/pvol3TB/data_environmental/ROMS_2k_files/ocean_his_TrackingSetup_6hourlycurrents_28days_Region",sprintf("%02d",i),".Rdata"))
  npp_crop <- rast(paste0(env.dir3,"ocean_his_TrackingSetup_NPP_Region",sprintf("%02d",i),".tif"))
  xy.pos <- cbind(Rdat6h$x[dat.flux$knnID], Rdat6h$y[dat.flux$knnID])
  dat.flux$cellID <- cellFromXY(npp_crop, xy.pos)
  rm(Rdat6h)
  ra.2d.flux <- rast(npp_crop)
  #sel.na <- which(is.na(dat.flux$cellID))
  # ra.2d.flux[dat.flux$cellID[-sel.na]] <- dat.flux$freq[-sel.na]
  ra.2d.flux[dat.flux$cellID] <- dat.flux$freq
  ## lastly, save tif file
  writeRaster(ra.2d.flux, file=paste0("/pvol3TB/FAM_outputs/tracking2D_Region",sprintf("%02d",i),"_NPP9_",speed,"mday_21days_",t,"r",r,"_28days_flux.tif"))
}
 

###################################################################################################
####### OPTIONAL: SHOW PARTICLE TRACKS UNTIL SEDIMENTATION (TAKES ~1h per region for 1% of tracks)
###################################################################################################
#### 2D analysis of NPP9, particle trajectories (subset of 100,000)
for(k in 1:10){
## load particle tracking run (takes 5-10minutes)
load(paste0("/pvol3TB/FAM_outputs/tracking2D_Region",sprintf("%02d",k),"_NPP9_",speed,"mday_21days_",t,"r",r,"_28days_Region",sprintf("%02d",k),".Rdata"))
## matrix that stores 1% of particle tracks
scaling <- 100
tracks <- matrix(NA, ncol=1356, nrow=ceiling(length(track.2D$id_list[[1]])/scaling))
## read first and second 6h, and store the ones that settle into list, keep others for next step
track.b <- as.matrix(do.call(cbind, track.2D$idx_list_2D[[1]])) ## all cells visited at that time step
track.c <- track.b                                                  ## just to make it work with the below loop
stopped <- which(track.2D$id_list[[1]]%!in%track.2D$id_list[[2]])   ## particles that stopped after first loop
stopped.reduced <- ceiling(length(stopped)/scaling)                     ## reduce to manageable size
tracks.ind <- track.c[sample(stopped,stopped.reduced),]             ## 
tracks[1:nrow(tracks.ind),1:ncol(tracks.ind)] <- tracks.ind         ## fill in matrix
row.num <- 0
## now loop across the rest and save output
for(i in 2:112){
  print(i)
  prev.stopped <- stopped
  prev.c <- track.c
  row.num <- row.num + nrow(tracks.ind)
  track.b <- as.matrix(do.call(cbind, track.2D$idx_list_2D[[i]]))
  track.c  <- cbind(track.c[-prev.stopped,],track.b)
  stopped <- which(track.2D$id_list[[i]]%!in%track.2D$id_list[[i+1]])
  stopped.reduced <- ceiling(length(stopped)/scaling)                          
  tracks.ind <- track.c[sample(stopped,stopped.reduced),]     
  tracks[(row.num+1):(row.num+nrow(tracks.ind)),1:ncol(tracks.ind)] <- tracks.ind
}
prev.stopped <- stopped
prev.c <- track.c
row.num <- row.num + nrow(tracks.ind)
track.b <- as.matrix(do.call(cbind, track.2D$idx_list_2D[[113]]))
track.c <- cbind(track.c[-prev.stopped,],track.b)
stopped.reduced <- ceiling(length(track.2D$id_list[[113]])/scaling)                        
tracks.ind <- track.c[sample(1:length(track.2D$id_list[[113]]),stopped.reduced),]
na.left <- which(is.na(tracks[,1]))
tracks[na.left,1:ncol(tracks.ind)] <- tracks.ind[1:length(na.left),]


###################
save(tracks, file=paste0("/pvol3TB/FAM_outputs/tracking2D_Region",sprintf("%02d",k),"_NPP9_",speed,"mday_21days_",t,"r",r,"_28days_Region",sprintf("%02d",i),"_1in100tracks.Rdata"))
rm(track.b, track.c, tracks.ind, prev.c, prev.stopped, stopped, track.2D)

#### need to translate matrix of IDs to raster cell IDs:
## load Rdat file (takes 5 minutes) and add these values to the correct raster
load(paste0("/pvol3TB/data_environmental/ROMS_2k_files/ocean_his_TrackingSetup_6hourlycurrents_28days_Region",sprintf("%02d",k),".Rdata"))
npp_crop <- rast(paste0(env.dir3,"ocean_his_TrackingSetup_NPP_Region",sprintf("%02d",k),".tif"))

tracks.x <- matrix(NA, ncol=ncol(tracks), nrow=nrow(tracks))
tracks.y <- matrix(NA, ncol=ncol(tracks), nrow=nrow(tracks))
tracks.x[] <- Rdat6h$x[tracks]
tracks.y[] <- Rdat6h$y[tracks]

## convert tracks to lines
sp.lines <- list()
for(i in 1:nrow(tracks)){
  xy.raw <- cbind(tracks.x[i,],tracks.y[i,])
  xy.raw2 <- xy.raw[!duplicated(xy.raw),]
  xy <- matrix(xy.raw2[complete.cases(xy.raw2),], ncol=2)
  pts <- sp::SpatialPoints(xy)
  sp.lines[[i]] <- as(pts,"SpatialLines")
}
merged.lines <- do.call(rbind, sp.lines)
save(merged.lines, file=paste0("/pvol3TB/FAM_outputs/tracking2D_Region",sprintf("%02d",k),"_NPP9_",speed,"mday_21days_",t,"r",r,"_28days_1in100tracksLines.Rdata"))

k=1
load(paste0("/pvol3TB/FAM_outputs/tracking2D_Region",sprintf("%02d",k),"_NPP9_",speed,"mday_21days_",t,"r",r,"_28days_1in100tracksLines.Rdata"))
plot(merged.lines)



###################################################################################################
####### OPTIONAL: SHOW 28-day DISTANCE TRAVELLED
###################################################################################################
#### for each cell, calculate the max distance travelled to get there

## No starting seafloor npp pts saved, so cannot get exact particle movements from positions.
## But we do know cell IDs where particles were seeded...

## need particle ID in each cell
## particle location at start
## straight line euclidian distance from start to finish

#### DETAILED TRAJECTORIES FOR 3D DON'T EXIST!!! FOR CODE CHECK "Run2kOffline...testruns.R" script
##load("/pvol3TB/FAM_outputs/tracking3D_Region01_NPP9_200mday_21days_xyzbottomnpp.Rdata")
speed <- 200
pr <- 0.0001
r <- "0001"
traj <- TRUE
if(traj==TRUE){t <- "traj_"} else {t <- ""}

## 3 AND 4 MISSING, DO ON PTRACKVM

for(j in 5:10){
message(j)
#### DETAILED TRAJECTORIES FOR 2D, USING CELL INDICES
load(paste0("/pvol3TB/FAM_outputs/tracking2D_Region",sprintf("%02d",j),"_NPP9_",speed,"mday_21days_",t,"r",r,"_28days.Rdata"))

## add starting cell-indices and create a column for stopping indices:
c.idx.mat <- matrix(unlist(track.2D$idx_list_2D[[1]][[1]]))
c.idx.mat <- cbind(c.idx.mat, NA)

## which particles stop at the end of the first loop
pts.not.stopped <- which(track.2D$id_list[[1]]%in%track.2D$id_list[[2]])
## and fill the indices of those that have stopped
c.idx.mat[-pts.not.stopped,2] <- unlist(track.2D$idx_list_2D[[1]][[12]])[-pts.not.stopped]
## update values for cells if larger than previous value
for(i in 2:112){
  #print(i)
  ## which particles stop at the end of the tracking-loop
  pts.not.stopped <- which(track.2D$id_list[[i]]%in%track.2D$id_list[[i+1]])
  ## translate these to overall row-indices
  pts.sel <- track.2D$id_list[[i]][-pts.not.stopped]
  ## and fill the indices of those that have stopped
  c.idx.mat[pts.sel,2] <- unlist(track.2D$idx_list_2D[[i]][[12]])[-pts.not.stopped]
}  
## and for the last loop, score the indices of all remaining particles
c.idx.mat[track.2D$id_list[[113]],2] <- unlist(track.2D$idx_list_2D[[113]][[12]])
print("particle IDs done")
## extract x-y-positions of cell indices, then calculate distances between start and end
load(paste0("/pvol3TB/data_environmental/ROMS_2k_files/ocean_his_TrackingSetup_6hourlycurrents_28days_Region",sprintf("%02d",j),".Rdata"))
xy.pos.start <- cbind(Rdat6h$x[c.idx.mat[,1]], Rdat6h$y[c.idx.mat[,1]])
xy.pos.end   <- cbind(Rdat6h$x[c.idx.mat[,2]], Rdat6h$y[c.idx.mat[,2]])
rm(Rdat6h)
npp_crop <- rast(paste0(env.dir3,"ocean_his_TrackingSetup_NPP_Region",sprintf("%02d",j),".tif"))
## triangle distance
xy.dist <- sqrt((xy.pos.end[,1]-xy.pos.start[,1])^2 + (xy.pos.start[,2]-xy.pos.end[,2])^2)
## now add into raster for plotting
distance.dat <- data.frame("knncell.idx.start" = c.idx.mat[,1],
                           "knncell.idx.end"   = c.idx.mat[,2],
                           "ra.idx.start"= cellFromXY(npp_crop, xy.pos.start),
                           "ra.idx.end"  = cellFromXY(npp_crop, xy.pos.end),
                           "x.crds.start" = xy.pos.start[,1],
                           "y.crds.start" = xy.pos.start[,2],
                           "x.crds.end" = xy.pos.end[,1],
                           "y.crds.end" = xy.pos.end[,2],
                           "xy.dist" = xy.dist)
rm(xy.dist)
save(distance.dat,         file=paste0("/pvol3TB/FAM_outputs/tracking2D_Region",sprintf("%02d",j),"_NPP9_",speed,"mday_21days_",t,"r",r,"_28days_distdat.Rdata"))

## need to take into account that some cells have multiple values (from multiple particles)
## we take the maximum distance and average
## here's for the distance to source
dist.to.source.max <- rast(npp_crop)
dist.to.source.avg <- rast(npp_crop)
## Aggregate the max and avg distances for each cell
max.dist <- aggregate(distance.dat$xy.dist, by = list(distance.dat$ra.idx.end), FUN = max)
avg.dist <- aggregate(distance.dat$xy.dist, by = list(distance.dat$ra.idx.end), FUN = mean)
## Fill the raster with the max and avg distances
dist.to.source.max[max.dist$Group.1] <- max.dist$x
dist.to.source.avg[avg.dist$Group.1] <- avg.dist$x

## and here for distance to sink
dist.to.sink.max <- rast(npp_crop)
dist.to.sink.avg <- rast(npp_crop)
## Aggregate the max and avg distances for each cell
max.dist.s <- aggregate(distance.dat$xy.dist, by = list(distance.dat$ra.idx.start), FUN = max)
avg.dist.s <- aggregate(distance.dat$xy.dist, by = list(distance.dat$ra.idx.start), FUN = mean)
## Fill the raster with the max and avg distances
dist.to.sink.max[max.dist.s$Group.1] <- max.dist.s$x
dist.to.sink.avg[avg.dist.s$Group.1] <- avg.dist.s$x

writeRaster(dist.to.source.max,file=paste0("/pvol3TB/FAM_outputs/tracking2D_Region",sprintf("%02d",j),"_NPP9_",speed,"mday_21days_",t,"r",r,"_28days_disttosourcemax.tif"))
writeRaster(dist.to.source.avg,file=paste0("/pvol3TB/FAM_outputs/tracking2D_Region",sprintf("%02d",j),"_NPP9_",speed,"mday_21days_",t,"r",r,"_28days_disttosourceavg.tif"))
writeRaster(dist.to.sink.max,  file=paste0("/pvol3TB/FAM_outputs/tracking2D_Region",sprintf("%02d",j),"_NPP9_",speed,"mday_21days_",t,"r",r,"_28days_disttosinkmax.tif"))
writeRaster(dist.to.sink.avg,  file=paste0("/pvol3TB/FAM_outputs/tracking2D_Region",sprintf("%02d",j),"_NPP9_",speed,"mday_21days_",t,"r",r,"_28days_disttosinkavg.tif"))
}



########################################################################
##### Compile regional FAM raster files into a single circumpolar file
########################################################################

#### start with the current speed files which are circumpolar and extract the coordinates for the values
uv.max <- rast(paste0(env.dir,"ocean_his_bottom_uv_max.tif"))
circant.coords <- crds(uv.max)

#### regional FAM files
fam.dir <- "/pvol/FAM_outputs/"
fam.dir3 <- "/pvol3TB/FAM_outputs/"
fam.sed.list <- list()
fam.flux.list <- list()
for(i in 1:10){
  fam.file.base <- paste0("tracking2D_Region",sprintf("%02d",i),"_NPP9_200mday_21days_traj_r0005_28days_")
  ## sed files
  fam.file <- paste0(fam.dir3, fam.file.base, "sed.tif")
  if(file.exists(fam.file)){
    fam.sed.list[[i]] <- rast(fam.file)
  }else{
    fam.sed.list[[i]] <- rast(paste0(fam.dir, fam.file.base, "sed.tif"))
  }
  ## flux files
  fam.file <- paste0(fam.dir3, fam.file.base, "flux.tif")
  if(file.exists(fam.file)){
    fam.flux.list[[i]] <- rast(fam.file)
  }else{
    fam.flux.list[[i]] <- rast(paste0(fam.dir, fam.file.base, "flux.tif"))
  }
}

#### replace values in the current speed files with FAM values
flux.ra <- rast(uv.max)
not.na <- which(!is.na(uv.max[]))
for(i in 1:10){
  print(i)
  ## the values from that region
  flux.vals.added <- terra::extract(fam.flux.list[[i]], circant.coords)[,1]
  ## all values that are na in from that extract:
  not.na.in.added <- which(!is.na(flux.vals.added))
  ## remove any of these values that are smaller than those already added earlier
  flux.vals.existing <- terra::extract(flux.ra, circant.coords)[,1]
  if(any(flux.vals.added[not.na.in.added]<flux.vals.existing[not.na.in.added], na.rm=TRUE)){
    print("small vals detected")
    smaller.vals <- which(flux.vals.added[not.na.in.added]<flux.vals.existing[not.na.in.added])
    flux.ra[not.na][not.na.in.added,][-smaller.vals] <- flux.vals.added[not.na.in.added][-smaller.vals]
  }else{
  ## add the rest
  flux.ra[not.na][not.na.in.added,] <- flux.vals.added[not.na.in.added]
  }
}
plot(flux.ra)

sed.ra <- rast(uv.max)
not.na <- which(!is.na(uv.max[]))
for(i in 1:10){
  print(i)
  ## the values from that region
  sed.vals.added <- terra::extract(fam.sed.list[[i]], circant.coords)[,1]
  ## all values that are na in from that extract:
  not.na.in.added <- which(!is.na(sed.vals.added))
  ## remove any of these values that are smaller than those already added earlier
  sed.vals.existing <- terra::extract(sed.ra, circant.coords)[,1]
  if(any(sed.vals.added[not.na.in.added]<sed.vals.existing[not.na.in.added], na.rm=TRUE)){
    print("small vals detected")
    smaller.vals <- which(sed.vals.added[not.na.in.added]<sed.vals.existing[not.na.in.added])
    sed.ra[not.na][not.na.in.added,][-smaller.vals] <- sed.vals.added[not.na.in.added][-smaller.vals]
  }else{
    ## add the rest
    sed.ra[not.na][not.na.in.added,] <- sed.vals.added[not.na.in.added]
  }
}
plot(sed.ra)

writeRaster(flux.ra, file=paste0("/pvol/FAM_outputs/tracking2D_NPP9_200mday_21days_traj_r0005_28days_flux.tif"))
writeRaster(sed.ra, file=paste0("/pvol/FAM_outputs/tracking2D_NPP9_200mday_21days_traj_r0005_28days_sed.tif"))






##

