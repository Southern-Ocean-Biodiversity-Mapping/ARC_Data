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
#### setting up preliminary boundaries for each region, and reduce ROMS domain accordingly

load(paste0(env.dir,"ocean_his_0001_2slices_TrackingSetup_transposed.Rdata"))
## reduce the ROMS by excluding areas that are not relevant for the shelf, such as North of the Ross Sea
Rdat <- ROMS_2k_t

roms.coords.proj <- cbind(c(Rdat$x), c(Rdat$y))
x.range <- c(min(roms.coords.proj[,1])-1000,max(roms.coords.proj[,1])+1000)
y.range <- c(min(roms.coords.proj[,2])-1000,max(roms.coords.proj[,2])+1000)
empty.roms.ra <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000)

uv.ra      <- setValues(empty.roms.ra, sqrt(Rdat$i_u[,,1,1]^2 +c(Rdat$i_v[,,1,1])^2))
uv.surf.ra <- setValues(empty.roms.ra, sqrt(Rdat$i_u[,,31,1]^2+c(Rdat$i_v[,,31,1])^2))
h.ra       <- setValues(empty.roms.ra, Rdat$h)
h.ra[is.na(uv.ra[])] <- NA

#writeRaster(h.ra, filename=paste0(env.dir,"waom2_h.tif"))

lim <- matrix(NA, nrow=4, ncol=10)
lim[,1] <- c(1601,2400,1201,1800) ## Ross Sea
lim[,2] <- c(1701,2100,401,1300)  ## Amundsen Sea
lim[,3] <- c(851,1800,201,760)    ## Bellinghausen Sea
lim[,4] <- c(1,900,51,600)        ## Peninsula
lim[,5] <- c(201,1000,501,1300)   ## Weddell
lim[,6] <- c(1,550,1201,2500)   ## Lazarev
lim[,7] <- c(401,1300,2401,2900)   ## Amery
lim[,8] <- c(1201,2100,2501,3000)   ## Mawson
lim[,9] <- c(2001,2600,1901,2700)   ## D'Urville
lim[,10] <- c(2301,2650,1301,2000)   ## Seamounts
ext1 <- ext(Rdat$x[lim[1,1],lim[3,1]], Rdat$x[lim[1,1],lim[4,1]],
            Rdat$y[lim[2,1],lim[3,1]], Rdat$y[lim[1,1],lim[3,1]])
ext2 <- ext(Rdat$x[lim[1,2],lim[3,2]], Rdat$x[lim[1,2],lim[4,2]],
            Rdat$y[lim[2,2],lim[3,2]], Rdat$y[lim[1,2],lim[3,2]])
ext3 <- ext(Rdat$x[lim[1,3],lim[3,3]], Rdat$x[lim[1,3],lim[4,3]],
            Rdat$y[lim[2,3],lim[3,3]], Rdat$y[lim[1,3],lim[3,3]])
ext4 <- ext(Rdat$x[lim[1,4],lim[3,4]], Rdat$x[lim[1,4],lim[4,4]],
            Rdat$y[lim[2,4],lim[3,4]], Rdat$y[lim[1,4],lim[3,4]])
ext5 <- ext(Rdat$x[lim[1,5],lim[3,5]], Rdat$x[lim[1,5],lim[4,5]],
            Rdat$y[lim[2,5],lim[3,5]], Rdat$y[lim[1,5],lim[3,5]])
ext6 <- ext(Rdat$x[lim[1,6],lim[3,6]], Rdat$x[lim[1,6],lim[4,6]],
            Rdat$y[lim[2,6],lim[3,6]], Rdat$y[lim[1,6],lim[3,6]])
ext7 <- ext(Rdat$x[lim[1,7],lim[3,7]], Rdat$x[lim[1,7],lim[4,7]],
            Rdat$y[lim[2,7],lim[3,7]], Rdat$y[lim[1,7],lim[3,7]])
ext8 <- ext(Rdat$x[lim[1,8],lim[3,8]], Rdat$x[lim[1,8],lim[4,8]],
            Rdat$y[lim[2,8],lim[3,8]], Rdat$y[lim[1,8],lim[3,8]])
ext9 <- ext(Rdat$x[lim[1,9],lim[3,9]], Rdat$x[lim[1,9],lim[4,9]],
            Rdat$y[lim[2,9],lim[3,9]], Rdat$y[lim[1,9],lim[3,9]])
ext10 <- ext(Rdat$x[lim[1,10],lim[3,10]], Rdat$x[lim[1,10],lim[4,10]],
             Rdat$y[lim[2,10],lim[3,10]], Rdat$y[lim[1,10],lim[3,10]])
plot(h.ra)
lines(ext1)
lines(ext2)
lines(ext3)
lines(ext4)
lines(ext5)
lines(ext6)
lines(ext7)
lines(ext8)
lines(ext9)
lines(ext10)

## setup ROMS and npp data and save for quick loading
# na.sel <- matrix(NA, nrow=4, ncol=10)
# na.sel[,1] <- c(601,800,1,100) ## Ross Sea 55k points
# na.sel[,2] <- c(351,400,301,800)  ## Amundsen Sea 56k points
# na.sel[,3] <- c()  ## Bellinghausen Sea
# na.sel[,4] <- c()  ## Peninsula
# na.sel[,5] <- c()   ## Weddell
# na.sel[,6] <- c()   ## Lazarev
# na.sel[,7] <- c()   ## Amery
# na.sel[,8] <- c()   ## Mawson
# na.sel[,9] <- c()   ## D'Urville
# na.sel[,10]<- c()   ## Seamounts

npp <- rast("/pvol/data_environmental/Circumpolar_EnvData_2km_NPP_Cafe_filled_SummerAverage.tif")
for(i in 1:10){
  print(i)
  row.lim <- lim[1,i]:lim[2,i]
  col.lim <- lim[3,i]:lim[4,i]
  ## ROMS
  Rdat <- list()
  Rdat$x    <- ROMS_2k_t$x[row.lim, col.lim]
  Rdat$y    <- ROMS_2k_t$y[row.lim, col.lim]
  Rdat$h    <- ROMS_2k_t$h[row.lim, col.lim]
  Rdat$hh   <- ROMS_2k_t$hh[row.lim, col.lim,]
  Rdat$zice <- ROMS_2k_t$zice[row.lim, col.lim]
  Rdat$i_u  <- ROMS_2k_t$i_u[row.lim, col.lim,,]
  Rdat$i_v  <- ROMS_2k_t$i_v[row.lim, col.lim,,]
  Rdat$i_w  <- ROMS_2k_t$i_w[row.lim, col.lim,,]
  save(Rdat, file=paste0(env.dir,"ocean_his_0001_2slices_TrackingSetup_transposed_",i,".Rdata"))
  ## npp
  roms.coords.proj <- cbind(c(Rdat$x), c(Rdat$y))
  x.range <- c(min(roms.coords.proj[,1])-1000,max(roms.coords.proj[,1])+1000)
  y.range <- c(min(roms.coords.proj[,2])-1000,max(roms.coords.proj[,2])+1000)
  empty.roms.ra <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000)
  npp2 <- project(npp, empty.roms.ra)
  npp2[npp2[]>0] <- 1
#  npp2[na.sel[1,i]:na.sel[2,i],na.sel[3,i]:na.sel[4,i]] <- NA
  ## thin out to one particle every 4 cells
  npp2[seq(1,nrow(npp2), by=2),] <- NA
  npp2[,seq(1,ncol(npp2), by=2)] <- NA
  ##
  cell.pts <- cbind(crds(as.points(npp2)),0)
  save(cell.pts, file=paste0(env.dir,"distancetravelled_seeded_cellpts_",i,".Rdata"))
}

uv.ra      <- setValues(empty.roms.ra, sqrt(Rdat$i_u[,,1,1]^2 +c(Rdat$i_v[,,1,1])^2))
h.ra <- setValues(empty.roms.ra, Rdat$h)
h.ra[is.na(uv.ra[])] <- NA

par(mfrow=c(1,2))
plot(h.ra)
plot(npp2)


###############################
########## LOOPS ##########
#### PARTICLE TRACKING
for(i in 1:10){
  message(i)
  #### load the ROMS data
  load(paste0(env.dir,"ocean_his_0001_2slices_TrackingSetup_transposed_",i,".Rdata"))
  
  ## load points to seed into 3D tracking
  load(paste0(env.dir,"distancetravelled_seeded_cellpts_",i,".Rdata"))
  print(dim(cell.pts))
  ## run 3D tracking
  start.time <- Sys.time()
  track.3D <- loopit_2D3D(pts=cell.pts, romsobject=Rdat, projected=TRUE, speed=100, domain="3D", runtime=30, roms_slices=2, looping_time=0.5, trajectories=TRUE, detailed_trajectories=TRUE)
  save(track.3D, file=paste0("/pvol/3_model_analysis/split_regions_runs_3Dtracking_30days_detailed_",i,".Rdata"))
  print(Sys.time()-start.time)
  
  ## points to seed into 2D tracking
  not.stopped <- which(track.3D$stopindex==0)
  pfloat <- matrix(track.3D$ptrack[not.stopped,,dim(track.3D$ptrack)[3]],ncol=3)
  pts.3D <- rbind(track.3D$pend, pfloat) ## stopped and floating points for 2D tracking
  ## run 2D tracking
  start.time <- Sys.time()
  track.2D <-loopit_2D3D(pts=pts.3D, romsobject=Rdat, projected=TRUE, speed=100, domain="2D", runtime=30, roms_slices=2, looping_time=0.5, trajectories=TRUE, detailed_trajectories=TRUE)
  save(track.2D, file=paste0("/pvol/3_model_analysis/split_regions_runs_2Dtracking_30days_detailed_",i,".Rdata"))
  print(Sys.time()-start.time)
}
  
#### DISTANCE CALCULATION
for(k in 1:10){
  message(k)
  
  #### load particle tracking output
  load(paste0("/pvol/3_model_analysis/split_regions_runs_3Dtracking_30days_detailed_",k,".Rdata"))
  load(paste0("/pvol/3_model_analysis/split_regions_runs_2Dtracking_30days_detailed_",k,".Rdata"))
  
  #### load the ROMS data
  load(paste0(env.dir,"ocean_his_0001_2slices_TrackingSetup_transposed_",k,".Rdata"))

  #### new distance data to be added
  Rdat$dist3D <- matrix(0, ncol=ncol(Rdat$x),nrow=nrow(Rdat$x))  ## distance to start
  Rdat$dist3Ds <- matrix(0, ncol=ncol(Rdat$x),nrow=nrow(Rdat$x)) ## distance from start
  Rdat$dist2D <- matrix(0, ncol=ncol(Rdat$x),nrow=nrow(Rdat$x))  ## distance to start
  Rdat$dist2Ds <- matrix(0, ncol=ncol(Rdat$x),nrow=nrow(Rdat$x)) ## distance from start
  
  #### 3D distances
  ## identify which row the positions of the particles are stored (particles get deleted at every time-step)
  pts.sel.start <- track.3D$id_list[[1]]
  ## ptrack x and y positions
  c.crds.start.xy <- cbind(track.3D$ptrack_x_list[[1]],track.3D$ptrack_y_list[[1]])
  ## starting cell indices
  c.idx.start <- matrix(unlist(track.3D$idx_list_2D[[1]]), ncol=length(track.3D$idx_list_2D[[1]]))[,1]
  ## update values for cells if larger than previous value
  for(i in 2:60){
    ## where are particles now
    pts.sel <- track.3D$id_list[[i]]
    ##
    c.crds.now.xy <- cbind(track.3D$ptrack_x_list[[i]],track.3D$ptrack_y_list[[i]])
    ## distance to where they started
    dist.dat <- rep(0, nrow(c.crds.now.xy))
    for(j in 1:nrow(c.crds.now.xy)){
      dist.dat[j] <- sqrt(sum((c.crds.now.xy[j, ] - c.crds.start.xy[pts.sel,][j, ])^2))
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
  
  #### 2D distances
  ## identify which row the positions of the particles are stored (particles get deleted at every time-step)
  pts.sel.start <- track.2D$id_list[[1]]
  ## ptrack x and y positions
  c.crds.start.xy <- cbind(track.2D$ptrack_x_list[[1]],track.2D$ptrack_y_list[[1]])
  ## starting cell indices
  c.idx.start <- matrix(unlist(track.2D$idx_list_2D[[1]]), ncol=length(track.2D$idx_list_2D[[1]]))[,1]
  ## update values for cells if larger than previous value
  for(i in 2:60){
    print(i)
    ## where are particles now
    pts.sel <- track.2D$id_list[[i]]
    ##
    c.crds.now.xy <- cbind(track.2D$ptrack_x_list[[i]],track.2D$ptrack_y_list[[i]])
    ## distance to where they started
    dist.dat <- rep(0, nrow(c.crds.now.xy))
    for(j in 1:nrow(c.crds.now.xy)){
      dist.dat[j] <- sqrt(sum((c.crds.now.xy[j, ] - c.crds.start.xy[pts.sel,][j, ])^2))
    }
    ## if larger than values already stored, then replace
    idx.m <- matrix(unlist(track.2D$idx_list_2D[[i]]), ncol=length(track.2D$idx_list_2D[[i]]))
    c.idx <- c(idx.m[pts.sel])
    repl.sel <- which(Rdat$dist2D[c.idx]<=dist.dat)
    Rdat$dist2D[c.idx][repl.sel] <- dist.dat[repl.sel]
    
    ## distance moved away from starting cell
    c.idx.start.update <- c.idx.start[pts.sel]
    repl.sel.start <- which(Rdat$dist2Ds[c.idx.start.update]<=dist.dat)
    Rdat$dist2Ds[c.idx.start.update][repl.sel.start] <- dist.dat[repl.sel.start]
    #### save output
  }
 save(Rdat, file=paste0(env.dir,"ocean_his_0001_2slices_TrackingSetup_transposed_distancerunscompleted_",k,".Rdata"))
}


#############################################################
#############################################################
## Analyse outputs:
#############################################################
stereo <- "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

ra_extents <- function(){
  
}

ra.extents <- list()
for(i in 1:10){
  print(i)
  load(paste0("/pvol/data_environmental/ROMS_2k_files/ocean_his_0001_2slices_TrackingSetup_transposed_distancerunscompleted_",i,".Rdata"))
  x.range <- c(min(Rdat$x)-1000,max(Rdat$x)+1000)
  y.range <- c(min(Rdat$y)-1000,max(Rdat$y)+1000)
  ra.extents[[i]] <- ext(c(x.range, y.range))
}

ra_setup_and_plot <- function(k){
  load(paste0("/pvol/data_environmental/ROMS_2k_files/ocean_his_0001_2slices_TrackingSetup_transposed_distancerunscompleted_",k,".Rdata"))
  dist.2D  <- rast(extent=ra.extents[[k]], crs=stereo, resolution=2000, vals=Rdat$dist2D)
  dist.3D  <- rast(extent=ra.extents[[k]], crs=stereo, resolution=2000, vals=Rdat$dist3D)
  dist.2Ds <- rast(extent=ra.extents[[k]], crs=stereo, resolution=2000, vals=Rdat$dist2Ds)
  dist.3Ds <- rast(extent=ra.extents[[k]], crs=stereo, resolution=2000, vals=Rdat$dist3Ds)
  h.ra     <- rast(extent=ra.extents[[k]], crs=stereo, resolution=2000, vals=Rdat$h)

  if(k==1){
    sel <- c(10,2)
  } else if(k==10){
    sel <- c(9,1)
  } else sel <- c(k-1, k+1)
  par(mfrow=c(2,2))
  plot(dist.2D,  main="2D final pos, dist to start", range=c(0,500000))
  contour(h.ra, levels=1500, add=TRUE, lty=2)
  lines(ra.extents[[sel[1]]], lty=2)
  lines(ra.extents[[sel[2]]], lty=2)
  plot(dist.2Ds, main="2D start pos, dist to final", range=c(0,500000))
  contour(h.ra, levels=1500, add=TRUE, lty=2)
  lines(ra.extents[[sel[1]]], lty=2)
  lines(ra.extents[[sel[2]]], lty=2)
  plot(dist.3D,  main="3D final pos, dist to start", range=c(0,600000))
  contour(h.ra, levels=1500, add=TRUE, lty=2)
  lines(ra.extents[[sel[1]]], lty=2)
  lines(ra.extents[[sel[2]]], lty=2)
  plot(dist.3Ds, main="3D start pos, dist to final", range=c(0,600000))
  contour(h.ra, levels=1500, add=TRUE, lty=2)
  lines(ra.extents[[sel[1]]], lty=2)
  lines(ra.extents[[sel[2]]], lty=2)
  # par(mfrow=c(1,1))
  # plot(dist.3Ds, main="3D start pos, dist to final")
  # contour(h.ra, levels=1500, add=TRUE, lty=2)
}

## 1
ra_setup_and_plot(1)
ra_setup_and_plot(2)
ra_setup_and_plot(3)
ra_setup_and_plot(4) ## AP
ra_setup_and_plot(5) ## Weddell
ra_setup_and_plot(6)
ra_setup_and_plot(7)
ra_setup_and_plot(8)
ra_setup_and_plot(9)
ra_setup_and_plot(10)







