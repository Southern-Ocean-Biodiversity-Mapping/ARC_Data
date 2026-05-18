####
## This script:
## - reads in the original ROMS history files produced from the model run ("OcVel_his_...")
## - shift u- and v-points to rho-points and extract current velocities, save output temporarily (huge storage requirements)
## - setup base ROMS object
## - combine current velocity files to 24h (daily) files
## - split 24h files into regions
## - reduce 24h files to 6hourly outputs
## - calculate bottom temperature and salinity from 4k model
####

#### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #####
#### DON'T JUST RUN THE SCRIPT, IT'S HOURS OF COMPUTING TIME!!! #####
#### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #####


## specify user and setup directory to look up data from
usr <- "VM"
source("0_SourceFile.R")

## set folders
env.dir <- paste0(usr.main.dir,"data_environmental/derived/ROMS_2k_files/")
roms.dir <- paste0(usr.roms.dir,"data_environmental/raw/ROMS_2k_files/")
roms.dir2 <- paste0(usr.dropbox.dir,"data_environmental/raw/")
out.dir <- paste0(usr.roms.dir,"data_environmental/derived/ROMS_2k_files/")

############################

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

string.chr <- "Circumpolar_EnvData_"
string.res <- "500m_"

################################
#### ROMS Currents
### WE NEED A RUN ACROSS SUMMER WITH HIGH-RES HISTORY FILES!!! 8 has one month only, 8_Nov has 10-day intervals...
### load lon/lat information from ROMS-grid

## tidync code:
grd2k_nc <- tidync(paste0(roms.dir,"waom2extend_grd.nc"))
lon_rho.raw <- hyper_tibble(grd2k_nc, select_var="lon_rho") ## 3150 2800 - longitude
lat_rho.raw <- hyper_tibble(grd2k_nc, select_var="lat_rho") ## 3150 2800 - latitude
roms.coords.proj <- reproj(cbind(lon_rho.raw$lon_rho, lat_rho.raw$lat_rho), target=stereo)
rm(grd2k_nc,lon_rho.raw,lat_rho.raw)

## ncdf4 code
# grd2k_nc <- nc_open(paste0(env.dir,"waom2extend_grd.nc"))
# lon_rho.raw <- ncvar_get(grd2k_nc, varid="lon_rho") ## 3150 2800 - longitude
# lat_rho.raw <- ncvar_get(grd2k_nc, varid="lat_rho") ## 3150 2800 - latitude
# nc_close(grd2k_nc)
# rm(grd2k_nc)
# roms.coords.proj <- reproj(cbind(c(lon_rho.raw), c(lat_rho.raw)), target=stereo)

### create rasters
roms.coords.proj.lon <- matrix(roms.coords.proj[,1], ncol=2800)
roms.coords.proj.lat <- matrix(roms.coords.proj[,2], ncol=2800)
x.range <- c(min(roms.coords.proj[,1])-1000,max(roms.coords.proj[,1])+1000)
y.range <- c(min(roms.coords.proj[,2])-1000,max(roms.coords.proj[,2])+1000)
empty.roms.ra <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000)
empty.roms.ra2 <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000, nlyr=2)
empty.roms.ra31 <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000, nlyr=31)

################################
#### TO DO ONCE for each history file:
### shifted u and v current velocities to match rho-points, then extract interpolated current velocities at those rho-points
### for w, we simply remove the surface layer (32)
## 1. load current velocities
## 2. create a filled raster each for the original locations of the u, v, w points
## 3. extract interpolated u and v values at the rho points

x.u.range <- c(min(roms.coords.proj[,1]),max(roms.coords.proj[,1]))
empty.u.ra <- rast(extent=ext(c(x.u.range,y.range)), crs=stereo, resolution=2000, nlyr=31)
empty.u.ra.topbottom <- rast(extent=ext(c(x.u.range,y.range)), crs=stereo, resolution=2000, nlyr=2)
y.v.range <- c(min(roms.coords.proj[,2]),max(roms.coords.proj[,2]))
empty.v.ra <- rast(extent=ext(c(x.range,y.v.range)), crs=stereo, resolution=2000, nlyr=31)
empty.v.ra.topbottom <- rast(extent=ext(c(x.range,y.v.range)), crs=stereo, resolution=2000, nlyr=2)
ra.crds <- crds(empty.roms.ra)

### u and v together
## below code repeats for each time-slice
surface_and_seafloor_only <- TRUE

for(f in 1:31){ ## 31 days
  message(paste0("file ",f))
  #sel <- sel.list[[f]]
  dgts <- formatC(f, width=2, format="d", flag="0")
  nc_2k <- tidync(paste0(env.dir,"OcVel_his_00",dgts,".nc"))
  if(f==1) {sel <- 1:9
  }else sel <- 1:8
  for(k in 1:length(sel)){
    message(paste0("slice ",f,"-",k))
    idx <- sel[k]
    ####
    ### u first
    ## extract the values from the netcdf file
    u.raw.loop <- nc_2k %>% activate("u") %>% hyper_filter(ocean_time = index==idx)  %>% hyper_array(force=TRUE)
    if(surface_and_seafloor_only){
      ## empty raster to fill values into
      u.raw.loop.ra <- empty.u.ra.topbottom
      ##
      print("u top and bottom layer only")
      u.raw.loop.ra[[1]][] <- c(u.raw.loop$u[,2800:1,1]) #flip rows upside down, so the representation is correct
      u.raw.loop.ra[[2]][] <- c(u.raw.loop$u[,2800:1,31]) #flip rows upside down, so the representation is correct
      ## extract values from the raster at rho-coordinates: 1-2min per slice
      u.vals.loop <- extract(u.raw.loop.ra, ra.crds, method="bilinear")
      ## assign these values to a raster so we get the right format of rows and columns
      u.ra.loop <- setValues(empty.roms.ra2, u.vals.loop)  
    }else{
      ## empty raster to fill values into
      u.raw.loop.ra <- empty.u.ra
      ## fill values into the raw raster: ~ 3 min per slice
      for(i in 1:31){ ## 31 depth layers
        print(i)
        u.raw.loop.ra[[i]][] <- c(u.raw.loop$u[,2800:1,i]) #flip rows upside down, so the representation is correct
      }
      ## extract values from the raster at rho-coordinates: 1-2min per slice
      u.vals.loop <- extract(u.raw.loop.ra, ra.crds, method="bilinear")
      ## assign these values to a raster so we get the right format of rows and columns
      u.ra.loop <- setValues(empty.roms.ra31, u.vals.loop)  
    }
    writeRaster(u.ra.loop, file=paste0(roms.dir,"ocean_his_00",dgts,"_slices_u_",k,".tif"), overwrite=TRUE)
    ## save output data for each time-slice in an array so we can use them in the particle tracking
    # u.vals.array <- as.array(u.ra.loop)
    # save(u.vals.array, file=paste0(roms.dir,"ocean_his_00",dgts,"_slices_u_",k,".Rdata"))
    ##
    rm(u.vals.loop, u.raw.loop, u.raw.loop.ra, u.ra.loop)
    
    ####
    ### v second
    ## extract the values from the netcdf file
    v.raw.loop <- nc_2k %>% activate("v") %>% hyper_filter(ocean_time = index==idx)  %>% hyper_array(force=TRUE)
    if(surface_and_seafloor_only){
      ## empty raster to fill values into
      v.raw.loop.ra <- empty.v.ra.topbottom
      ##
      print("v top and bottom layer only")
      v.raw.loop.ra[[1]][] <- c(v.raw.loop$v[,2799:1,1]) #flip rows upside down, so the representation is correct
      v.raw.loop.ra[[2]][] <- c(v.raw.loop$v[,2799:1,31]) #flip rows upside down, so the representation is correct
      ## extract values from the raster at rho-coordinates: 1-2min per slice
      v.vals.loop <- extract(v.raw.loop.ra, ra.crds, method="bilinear")
      ## assign these values to a raster so we get the right format of rows and columns
      v.ra.loop <- setValues(empty.roms.ra2, v.vals.loop)  
    }else{
      ## empty raster to fill values into
      v.raw.loop.ra <- empty.v.ra
      ## fill values into the raw raster: ~ 3 min per slice
      for(i in 1:31){ ## 31 depth layers
        print(i)
        v.raw.loop.ra[[i]][] <- c(v.raw.loop$v[,2799:1,i]) #flip rows upside down, so the representation is correct
      }
      ## extract values from the raster at rho-coordinates: 1-2min per slice
      v.vals.loop <- extract(v.raw.loop.ra, ra.crds, method="bilinear")
      ## assign these values to a raster so we get the right format of rows and columns
      v.ra.loop <- setValues(empty.roms.ra31, v.vals.loop)  
    }
    writeRaster(v.ra.loop, file=paste0(roms.dir,"ocean_his_00",dgts,"_slices_v_",k,".tif"), overwrite=TRUE)
    ## save output data for each time-slice in an array so we can use them in the particle tracking
    # v.vals.array <- as.array(v.ra.loop)
    # save(v.vals.array, file=paste0(roms.dir,"ocean_his_00",dgts,"_slices_v_",k,".Rdata"))
    ##
    rm(v.vals.loop, v.raw.loop, v.raw.loop.ra, v.ra.loop)
    # v.ra.t <- setValues(empty.roms.ra31, v.vals.array)
    # plot(v.ra.t$lyr.1)
    
    ####
    ### w third
    ## extract the values from the netcdf file
    w.raw.loop <- nc_2k %>% activate("w") %>% hyper_filter(ocean_time = index==idx)  %>% hyper_array(force=TRUE)
    if(surface_and_seafloor_only){
      ## empty raster to fill values into
      w.ra.loop <- empty.roms.ra2
      ##
      print("w top and bottom layer only")
      w.ra.loop[[1]][] <- c(w.raw.loop$w[,2800:1,1]) #flip rows upside down, so the representation is correct
      w.ra.loop[[2]][] <- c(w.raw.loop$w[,2800:1,31]) #flip rows upside down, so the representation is correct
    }else{
      ## empty raster to fill values into
      w.ra.loop <- empty.roms.ra31
      ## fill values into the raw raster: ~ 3 min per slice
      for(i in 1:31){ ## 31 depth layers
        print(i)
        w.ra.loop[[i]][] <- c(w.raw.loop$w[,2800:1,i]) #flip rows upside down, so the representation is correct
      }
    }
    writeRaster(w.ra.loop, file=paste0(roms.dir,"ocean_his_00",dgts,"_slices_w_",k,".tif"), overwrite=TRUE)
    ## save output data for each time-slice in an array so we can use them in the particle tracking
    # w.vals.array <- as.array(w.ra.loop)
    # save(w.vals.array, file=paste0(roms.dir,"ocean_his_00",dgts,"_slices_w_",k,".Rdata"))
    # w.ra.t <- setValues(empty.roms.ra31, w.vals.array)
    # plot(w.ra.t$lyr.1)
  }}


################################
#### setup base ROMS object
# nc_2k <- nc_open(paste0(roms.dir,"ocean_his_0003.nc"))
nc_2k <- nc_open(paste0(roms.dir,"OcVel_his_0001.nc"))
s_rho <- ncvar_get(nc_2k, varid="s_rho") ## 31             - vertical levels varying between 0 to 1
Cs_r  <- ncvar_get(nc_2k, varid="Cs_r")  ## 31             - "S-coordinate stretching curves at RHO-points
hc    <- ncvar_get(nc_2k, varid="hc")    ##                - value of 250
h     <- ncvar_get(nc_2k, varid="h")     ## 3150 2800      - model bathymetry
zeta  <- ncvar_get(nc_2k, varid="zeta")  ## 3150 2800 24    - free surface elevation
zice  <- ncvar_get(nc_2k, varid="zice")  ## 3150 2800      - ice draft
nc_close(nc_2k)

# nc_2k <- tidync(paste0(roms.dir,"OcVel_his_0001.nc"))
# h.test <- nc_2k %>% activate("h") %>% hyper_array(force=TRUE)
# h.test.t <- t(h.test$h[,2800:1])
# h.test.ra      <- setValues(empty.roms.ra, h.test.t)
# plot(h.test.ra)
# 
# u.test <- nc_2k %>% activate("u") %>% hyper_filter(ocean_time = index==1)  %>% hyper_array(force=TRUE)
# u.test.t <- aperm(u.test$u[,2800:1,], c(2,1,3))
# u.test.ra      <- setValues(empty.roms.ra, u.test.t)
# plot(h.test.ra)

### calculate depths of each vertical ROMS cell
## simplest calculation:
hh <- array(data=h*(Cs_r[1]), dim=c(3150,2800,1))
for(i in 2:31){
  hh <- abind(hh, h*(Cs_r[i]))
}

### common properties
## need to transpose the array so that the continent appears in the most commonly seen orientation
ROMS_2k_base <- list()
ROMS_2k_base$x <- t(roms.coords.proj.lon[,2800:1])
ROMS_2k_base$y <- t(roms.coords.proj.lat[,2800:1])
ROMS_2k_base$h <- t(h[,2800:1])
ROMS_2k_base$hh <- aperm(hh[,2800:1,], c(2,1,3))
ROMS_2k_base$zice <- t(zice[,2800:1])
# ## check if all correct
# h.ra      <- setValues(empty.roms.ra, ROMS_2k_base$h)
# plot(h.ra)
# hh.ra      <- setValues(empty.roms.ra, ROMS_2k_base$hh[,,31])
# plot(hh.ra)
save(ROMS_2k_base, file=paste0(out.dir,"ocean_his_TrackingSetup_base.Rdata"))
rm(ROMS_2k_base)

###################################
#### Current speeds for each day:
### compile into bottom and surface speed tifs:
curr.files.u <- list.files(roms.dir, pattern="slices_u", full.names=TRUE)
u.ra.bottom  <- rast(curr.files.u)[[seq(1,488, by=2)]]
u.ra.top     <- rast(curr.files.u)[[seq(2,488, by=2)]]
curr.files.v <- list.files(roms.dir, pattern="slices_v", full.names=TRUE)
v.ra.bottom  <- rast(curr.files.v)[[seq(1,488, by=2)]]
v.ra.top     <- rast(curr.files.v)[[seq(2,488, by=2)]]
curr.files.w <- list.files(roms.dir, pattern="slices_w", full.names=TRUE)
w.ra.bottom  <- rast(curr.files.w)[[seq(1,488, by=2)]]
w.ra.top     <- rast(curr.files.w)[[seq(2,488, by=2)]]
### calculate mean, max, and abs.mean current speeds
## first creating pisitive values to extract the absolute mean (so plus and minus don't cancel out)
u.ra.bottom.abs1 <- abs(u.ra.bottom[[1:61]])## higher values here will be due to tidal components
u.ra.bottom.abs2 <- abs(u.ra.bottom[[62:122]])
u.ra.bottom.abs3 <- abs(u.ra.bottom[[123:183]])
u.ra.bottom.abs4 <- abs(u.ra.bottom[[184:244]])
u.ra.bottom.abs <- c(u.ra.bottom.abs1,u.ra.bottom.abs2,u.ra.bottom.abs3,u.ra.bottom.abs4)
u.ra.bottom.abs.mean <- mean(u.ra.bottom.abs)
writeRaster(u.ra.bottom.abs.mean, file=paste0(out.dir,"ocean_his_bottom_u_absmean.tif"), overwrite=TRUE)
u.ra.bottom.mean <- mean(u.ra.bottom)
writeRaster(u.ra.bottom.mean, file=paste0(out.dir,"ocean_his_bottom_u_mean.tif"), overwrite=TRUE)
u.ra.bottom.max <- max(u.ra.bottom)
writeRaster(u.ra.bottom.max, file=paste0(out.dir,"ocean_his_bottom_u_max.tif"), overwrite=TRUE)
rm(u.ra.bottom.abs1,u.ra.bottom.abs2,u.ra.bottom.abs3,u.ra.bottom.abs4,v.ra.bottom.abs,u.ra.bottom.abs.mean,u.ra.bottom.mean,u.ra.bottom.max)
##
v.ra.bottom.abs1 <- abs(v.ra.bottom[[1:61]])
v.ra.bottom.abs2 <- abs(v.ra.bottom[[62:122]])
v.ra.bottom.abs3 <- abs(v.ra.bottom[[123:183]])
v.ra.bottom.abs4 <- abs(v.ra.bottom[[184:244]])
v.ra.bottom.abs <- c(v.ra.bottom.abs1,v.ra.bottom.abs2,v.ra.bottom.abs3,v.ra.bottom.abs4)
v.ra.bottom.abs.mean <- mean(v.ra.bottom.abs)
writeRaster(v.ra.bottom.abs.mean, file=paste0(out.dir,"ocean_his_bottom_v_absmean.tif"), overwrite=TRUE)
v.ra.bottom.mean <- mean(v.ra.bottom)
writeRaster(v.ra.bottom.mean, file=paste0(out.dir,"ocean_his_bottom_v_mean.tif"), overwrite=TRUE)
v.ra.bottom.max <- max(v.ra.bottom)
writeRaster(v.ra.bottom.max, file=paste0(out.dir,"ocean_his_bottom_v_max.tif"), overwrite=TRUE)
rm(v.ra.bottom.abs1,v.ra.bottom.abs2,v.ra.bottom.abs3,v.ra.bottom.abs4,v.ra.bottom.abs,v.ra.bottom.abs.mean,v.ra.bottom.mean)
##
w.ra.bottom.abs1 <- abs(w.ra.bottom[[1:61]])
w.ra.bottom.abs2 <- abs(w.ra.bottom[[62:122]])
w.ra.bottom.abs3 <- abs(w.ra.bottom[[123:183]])
w.ra.bottom.abs4 <- abs(w.ra.bottom[[184:244]])
w.ra.bottom.abs <- c(w.ra.bottom.abs1,w.ra.bottom.abs2,w.ra.bottom.abs3,w.ra.bottom.abs4)
w.ra.bottom.abs.mean <- mean(w.ra.bottom.abs)
writeRaster(w.ra.bottom.abs.mean, file=paste0(out.dir,"ocean_his_bottom_w_absmean.tif"), overwrite=TRUE)
w.ra.bottom.mean <- mean(w.ra.bottom)
writeRaster(w.ra.bottom.mean, file=paste0(out.dir,"ocean_his_bottom_w_mean.tif"), overwrite=TRUE)
w.ra.bottom.max <- max(w.ra.bottom)
writeRaster(w.ra.bottom.max, file=paste0(out.dir,"ocean_his_bottom_w_max.tif"), overwrite=TRUE)
rm(w.ra.bottom.abs1,w.ra.bottom.abs2,w.ra.bottom.abs3,w.ra.bottom.abs4,w.ra.bottom.abs.mean,w.ra.bottom.mean)
##
u.ra.bottom.mean <- rast(paste0(out.dir,"ocean_his_bottom_u_mean.tif"))
v.ra.bottom.mean <- rast(paste0(out.dir,"ocean_his_bottom_v_mean.tif"))
uv.ra.bottom.mean <- sqrt(u.ra.bottom.mean^2+v.ra.bottom.mean^2)
writeRaster(uv.ra.bottom.mean, file=paste0(out.dir,"ocean_his_bottom_uv_mean.tif"), overwrite=TRUE)
##
u.ra.bottom.abs.mean <- rast(paste0(out.dir,"ocean_his_bottom_u_absmean.tif"))
v.ra.bottom.abs.mean <- rast(paste0(out.dir,"ocean_his_bottom_v_absmean.tif"))
uv.ra.bottom.abs.mean <- sqrt(u.ra.bottom.abs.mean^2+v.ra.bottom.abs.mean^2)
writeRaster(uv.ra.bottom.abs.mean, file=paste0(out.dir,"ocean_his_bottom_uv_absmean.tif"), overwrite=TRUE)
v.ra.bottom.max <- rast(paste0(out.dir,"ocean_his_bottom_v_max.tif"))
u.ra.bottom.max <- rast(paste0(out.dir,"ocean_his_bottom_u_max.tif"))
uv.ra.bottom.max <- sqrt(u.ra.bottom.max^2+v.ra.bottom.max^2)
writeRaster(uv.ra.bottom.max, file=paste0(out.dir,"ocean_his_bottom_uv_max.tif"), overwrite=TRUE)
##

### current speeds: 8 slices at 3h give a 24 compilation
ROMS_2k <- list()
ROMS_2k$i_u <- array(NA, c(2800,3150,31,9))
ROMS_2k$i_v <- array(NA, c(2800,3150,31,9))
ROMS_2k$i_w <- array(NA, c(2800,3150,31,9))
for(j in 1:9){
  print(j)
  ROMS_2k$i_u[,,,j] <- as.array(rast(paste0(roms.dir,"ocean_his_0001_slices_u_",j,".tif")))
  ROMS_2k$i_v[,,,j] <- as.array(rast(paste0(roms.dir,"ocean_his_0001_slices_v_",j,".tif")))
  ROMS_2k$i_w[,,,j] <- as.array(rast(paste0(roms.dir,"ocean_his_0001_slices_w_",j,".tif")))
}
save(ROMS_2k, file=paste0(roms.dir,"ocean_his_TrackingSetup_24hcurrents_File01.Rdata"))

for(i in 2:30){
  print(i)  
  ROMS_2k <- list()
  ROMS_2k$i_u <- array(NA, c(2800,3150,31,8))
  ROMS_2k$i_v <- array(NA, c(2800,3150,31,8))
  ROMS_2k$i_w <- array(NA, c(2800,3150,31,8))
  dgts <- formatC(i, width=2, format="d", flag="0")
  for(j in 1:8){
    #  for(j in 1:3){
    ROMS_2k$i_u[,,,j] <- as.array(rast(paste0(env.dir,"ocean_his_00",dgts,"_slices_u_",j,".tif")))
    ROMS_2k$i_v[,,,j] <- as.array(rast(paste0(env.dir,"ocean_his_00",dgts,"_slices_v_",j,".tif")))
    ROMS_2k$i_w[,,,j] <- as.array(rast(paste0(env.dir,"ocean_his_00",dgts,"_slices_w_",j,".tif")))
  }
  save(ROMS_2k, file=paste0(roms.dir,"ocean_his_TrackingSetup_24hcurrents_File",dgts,".Rdata"))
}

roms.coords.proj <- cbind(c(ROMS_2k_base$x), c(ROMS_2k_base$y))
x.range <- c(min(roms.coords.proj[,1])-1000,max(roms.coords.proj[,1])+1000)
y.range <- c(min(roms.coords.proj[,2])-1000,max(roms.coords.proj[,2])+1000)
empty.roms.ra <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000)
h.ra      <- setValues(empty.roms.ra, ROMS_2k_base$h)
u.ra      <- setValues(empty.roms.ra, ROMS_2k$i_u[,,1,1])


#################################
#### split into regional datasets
load(paste0(roms.dir,"ocean_his_TrackingSetup_24hcurrents_File01.Rdata"))
load(paste0(out.dir,"ocean_his_TrackingSetup_base.Rdata"))

lim <- matrix(NA, nrow=4, ncol=10)
lim[,1] <- c(1751,2400,1201,1800) ## Ross Sea           650 x 600 = 390k
lim[,2] <- c(1701,2100,401,1300)  ## Amundsen Sea       400 x 900 = 360k
lim[,3] <- c(851,1800,201,760)    ## Bellinghausen Sea  950 x 560 = 532k
lim[,4] <- c(1,900,51,600)        ## Peninsula          900 x 550 = 495k
lim[,5] <- c(401,1000,401,1300)   ## Weddell            800 x 800 = 640k
lim[,6] <- c(1,550,1101,2500)     ## Lazarev            550 x 1300= 715k
lim[,7] <- c(351,1300,2401,2900)   ## Amery             900 x 500 = 450k
lim[,8] <- c(1201,2100,2501,3000)   ## Mawson           900 x 500 = 450k
lim[,9] <- c(2001,2600,1901,2700)   ## D'Urville        600 x 800 = 480k
lim[,10] <- c(2301,2650,1301,2000)   ## Seamounts       350 x 700 = 245k

for(i in 1:10){
  print(i)
  row.lim <- lim[1,i]:lim[2,i]
  col.lim <- lim[3,i]:lim[4,i]
  ## ROMS
  Rdat <- list()
  Rdat$x    <- ROMS_2k_base$x[row.lim, col.lim]
  Rdat$y    <- ROMS_2k_base$y[row.lim, col.lim]
  Rdat$h    <- ROMS_2k_base$h[row.lim, col.lim]
  Rdat$hh   <- ROMS_2k_base$hh[row.lim, col.lim,]
  Rdat$zice <- ROMS_2k_base$zice[row.lim, col.lim]
  Rdat$i_u  <- ROMS_2k$i_u[row.lim, col.lim,,]
  Rdat$i_v  <- ROMS_2k$i_v[row.lim, col.lim,,]
  Rdat$i_w  <- ROMS_2k$i_w[row.lim, col.lim,,]
  save(Rdat, file=paste0(roms.dir,"ocean_his_TrackingSetup_24hcurrents_File01_Region",sprintf("%02d",i),".Rdata"))
}

#### split into regional datasets
for(k in 1:30){
  message(k)
  load(paste0(roms.dir,"ocean_his_TrackingSetup_24hcurrents_File",sprintf("%02d",k),".Rdata"))
  for(i in 1:10){
    print(i)
    row.lim <- lim[1,i]:lim[2,i]
    col.lim <- lim[3,i]:lim[4,i]
    ## ROMS
    Rdat <- list()
    Rdat$x    <- ROMS_2k_base$x[row.lim, col.lim]
    Rdat$y    <- ROMS_2k_base$y[row.lim, col.lim]
    Rdat$h    <- ROMS_2k_base$h[row.lim, col.lim]
    Rdat$hh   <- ROMS_2k_base$hh[row.lim, col.lim,]
    Rdat$zice <- ROMS_2k_base$zice[row.lim, col.lim]
    Rdat$i_u  <- ROMS_2k$i_u[row.lim, col.lim,,]
    Rdat$i_v  <- ROMS_2k$i_v[row.lim, col.lim,,]
    Rdat$i_w  <- ROMS_2k$i_w[row.lim, col.lim,,]
    save(Rdat, file=paste0(roms.dir,"ocean_his_TrackingSetup_24hcurrents_File",sprintf("%02d",k),"_Region",sprintf("%02d",i),".Rdata"))
  }
}


#### combine current files into one for tracking
### the full 28 days crashes the VM
# load(paste0(roms.dir,"ocean_his_TrackingSetup_24hcurrents_Region01.Rdata"))
# Rdat.raw <- Rdat
# load(paste0(roms.dir,"ocean_his_TrackingSetup_24hcurrents_Region01_b.Rdata"))
# #load(paste0(roms.dir,"ocean_his_TrackingSetup_24hcurrents_Region01_c.Rdata"))
# Rdat.raw$i_u <- abind(Rdat.raw$i_u, Rdat$i_u)
# Rdat.raw$i_v <- abind(Rdat.raw$i_v, Rdat$i_v)
# Rdat.raw$i_w <- abind(Rdat.raw$i_w, Rdat$i_w)
# Rdat <- Rdat.raw
# rm(Rdat.raw)
# save(Rdat, file=paste0(roms.dir,"ocean_his_TrackingSetup_24hcurrents_Region01a.Rdata"))
### reducing to 6h outputs for 28 days for the 3D tracking
### 3h outputs for the 2D tracking should work!!!
library(abind)
for(k in 1:10){
  message(k)
  load(paste0(roms.dir,"ocean_his_TrackingSetup_24hcurrents_File01_Region",sprintf("%02d",k),".Rdata"))
  Rdat6h <- Rdat
  sel <- seq(1, 9, by=2)
  Rdat6h$i_u <- Rdat6h$i_u[,,,sel]
  Rdat6h$i_v <- Rdat6h$i_v[,,,sel]
  Rdat6h$i_w <- Rdat6h$i_w[,,,sel]
  sel <- seq(2, 8, by=2)
  for(i in 2:28){
    print(i)
    load(paste0(roms.dir,"ocean_his_TrackingSetup_24hcurrents_File",sprintf("%02d",i),"_Region",sprintf("%02d",k),".Rdata"))
    Rdat6h$i_u <- abind(Rdat6h$i_u, Rdat$i_u[,,,sel])
    Rdat6h$i_v <- abind(Rdat6h$i_v, Rdat$i_v[,,,sel])
    Rdat6h$i_w <- abind(Rdat6h$i_w, Rdat$i_w[,,,sel])
  }
  save(Rdat6h, file=paste0(roms.dir,"ocean_his_TrackingSetup_6hourlycurrents_28days_Region",sprintf("%02d",k),".Rdata"))
}

#### region 6 is too big to do it all at once:
load(paste0(roms.dir,"ocean_his_TrackingSetup_24hcurrents_File01_Region06.Rdata"))
## a
Rdat6h <- Rdat
sel <- seq(1, 9, by=2)
Rdat6h$i_u <- Rdat6h$i_u[,,,sel]
Rdat6h$i_v <- Rdat6h$i_v[,,,sel]
Rdat6h$i_w <- Rdat6h$i_w[,,,sel]
sel <- seq(2, 8, by=2)
for(i in 2:10){
  print(i)
  load(paste0(roms.dir,"ocean_his_TrackingSetup_24hcurrents_File",sprintf("%02d",i),"_Region06.Rdata"))
  Rdat6h$i_u <- abind(Rdat6h$i_u, Rdat$i_u[,,,sel])
  Rdat6h$i_v <- abind(Rdat6h$i_v, Rdat$i_v[,,,sel])
  Rdat6h$i_w <- abind(Rdat6h$i_w, Rdat$i_w[,,,sel])
}
save(Rdat6h, file=paste0(roms.dir,"ocean_his_TrackingSetup_6hourlycurrents_28days_Region06a.Rdata"))
## b
load(paste0(roms.dir,"ocean_his_TrackingSetup_24hcurrents_File11_Region06.Rdata"))
Rdat6h <- Rdat
sel <- seq(2, 8, by=2)
Rdat6h$i_u <- Rdat6h$i_u[,,,sel]
Rdat6h$i_v <- Rdat6h$i_v[,,,sel]
Rdat6h$i_w <- Rdat6h$i_w[,,,sel]
for(i in 12:20){
  print(i)
  load(paste0(roms.dir,"ocean_his_TrackingSetup_24hcurrents_File",sprintf("%02d",i),"_Region06.Rdata"))
  Rdat6h$i_u <- abind(Rdat6h$i_u, Rdat$i_u[,,,sel])
  Rdat6h$i_v <- abind(Rdat6h$i_v, Rdat$i_v[,,,sel])
  Rdat6h$i_w <- abind(Rdat6h$i_w, Rdat$i_w[,,,sel])
}
save(Rdat6h, file=paste0(roms.dir,"ocean_his_TrackingSetup_6hourlycurrents_28days_Region06b.Rdata"))
## c
load(paste0(roms.dir,"ocean_his_TrackingSetup_24hcurrents_File21_Region06.Rdata"))
Rdat6h <- Rdat
sel <- seq(2, 8, by=2)
Rdat6h$i_u <- Rdat6h$i_u[,,,sel]
Rdat6h$i_v <- Rdat6h$i_v[,,,sel]
Rdat6h$i_w <- Rdat6h$i_w[,,,sel]
for(i in 22:28){
  print(i)
  load(paste0(roms.dir,"ocean_his_TrackingSetup_24hcurrents_File",sprintf("%02d",i),"_Region06.Rdata"))
  Rdat6h$i_u <- abind(Rdat6h$i_u, Rdat$i_u[,,,sel])
  Rdat6h$i_v <- abind(Rdat6h$i_v, Rdat$i_v[,,,sel])
  Rdat6h$i_w <- abind(Rdat6h$i_w, Rdat$i_w[,,,sel])
}
save(Rdat6h, file=paste0(roms.dir,"ocean_his_TrackingSetup_6hourlycurrents_28days_Region06c.Rdata"))
## bind together
load(paste0(roms.dir,"ocean_his_TrackingSetup_6hourlycurrents_28days_Region06c.Rdata"))
c <- Rdat6h
load(paste0(roms.dir,"ocean_his_TrackingSetup_6hourlycurrents_28days_Region06b.Rdata"))
Rdat6h$i_u <- abind(Rdat6h$i_u, c$i_u)
Rdat6h$i_v <- abind(Rdat6h$i_v, c$i_v)
Rdat6h$i_w <- abind(Rdat6h$i_w, c$i_w)
save(Rdat6h, file=paste0(roms.dir,"ocean_his_TrackingSetup_6hourlycurrents_28days_Region06bc.Rdata"))
## step by ...
load(paste0(roms.dir,"ocean_his_TrackingSetup_6hourlycurrents_28days_Region06bc.Rdata"))
b <- Rdat6h
load(paste0(roms.dir,"ocean_his_TrackingSetup_6hourlycurrents_28days_Region06a.Rdata"))
Rdat6h$i_u <- abind(Rdat6h$i_u, b$i_u)
Rdat6h$i_v <- abind(Rdat6h$i_v, b$i_v)
save(Rdat6h, file=paste0(roms.dir,"ocean_his_TrackingSetup_6hourlycurrents_28days_Region06abcuv.Rdata"))
## ... step
load(paste0(roms.dir,"ocean_his_TrackingSetup_6hourlycurrents_28days_Region06bc.Rdata"))
b <- Rdat6h
rm(Rdat6h)
b$i_u <- NA
b$i_v <- NA
b$hh <- NA
load(paste0(roms.dir,"ocean_his_TrackingSetup_6hourlycurrents_28days_Region06abcuv.Rdata"))
Rdat6h$i_w <- abind(Rdat6h$i_w, b$i_w)
save(Rdat6h, file=paste0(roms.dir,"ocean_his_TrackingSetup_6hourlycurrents_28days_Region06.Rdata"))
##


################################################
library(ncdf4)
library(terra)
#### ROMS Bottom Temperature & Salinity
## 4k models for now
data.dat100 <- paste0(roms.dir2, "ROMS_4k_files/output_yr10/")
#### load lon/lat information from ROMS-grid
grd4k_nc <- nc_open(paste0(roms.dir2,"waom4extend_grd.nc"))
lon_rho <- ncvar_get(grd4k_nc, varid="lon_rho")
lat_rho <- ncvar_get(grd4k_nc, varid="lat_rho")
#### Prepare empty rasters to assign correct projected values to
roms.coords.proj <- rgdal::project(cbind(c(lon_rho), c(lat_rho)), proj=stereo)
x.range <- c(min(roms.coords.proj[,1])-2000,max(roms.coords.proj[,1])+2000)
y.range <- c(min(roms.coords.proj[,2])-2000,max(roms.coords.proj[,2])+2000)
empty.roms.ra <- rast(ext=ext(c(x.range,y.range)), crs=stereo, resolution=4000)

#depth
h <- rast(paste0(data.dat100,"ocean_avg_0001.nc"), subds="h")
## seafloor variables
s <- seq(1,217, by=31)
#salinity
salt.raw  <- c(subset(rast(paste0(data.dat100,"ocean_avg_0001.nc"), subds="salt"), subset=s),
               subset(rast(paste0(data.dat100,"ocean_avg_0002.nc"), subds="salt"), subset=s),
               subset(rast(paste0(data.dat100,"ocean_avg_0003.nc"), subds="salt"), subset=s),
               subset(rast(paste0(data.dat100,"ocean_avg_0004.nc"), subds="salt"), subset=s),
               subset(rast(paste0(data.dat100,"ocean_avg_0005.nc"), subds="salt"), subset=s),
               subset(rast(paste0(data.dat100,"ocean_avg_0006.nc"), subds="salt"), subset=s),
               subset(rast(paste0(data.dat100,"ocean_avg_0007.nc"), subds="salt"), subset=s),
               subset(rast(paste0(data.dat100,"ocean_avg_0008.nc"), subds="salt"), subset=s),
               subset(rast(paste0(data.dat100,"ocean_avg_0009.nc"), subds="salt"), subset=s),
               subset(rast(paste0(data.dat100,"ocean_avg_0010.nc"), subds="salt"), subset=s),
               subset(rast(paste0(data.dat100,"ocean_avg_0011.nc"), subds="salt"), subset=s),
               subset(rast(paste0(data.dat100,"ocean_avg_0012.nc"), subds="salt"), subset=s))
salt <- mean(salt.raw)
#temperature
temp.raw <- c(subset(rast(paste0(data.dat100,"ocean_avg_0001.nc"), subds="temp"), subset=s),
              subset(rast(paste0(data.dat100,"ocean_avg_0002.nc"), subds="temp"), subset=s),
              subset(rast(paste0(data.dat100,"ocean_avg_0003.nc"), subds="temp"), subset=s),
              subset(rast(paste0(data.dat100,"ocean_avg_0004.nc"), subds="temp"), subset=s),
              subset(rast(paste0(data.dat100,"ocean_avg_0005.nc"), subds="temp"), subset=s),
              subset(rast(paste0(data.dat100,"ocean_avg_0006.nc"), subds="temp"), subset=s),
              subset(rast(paste0(data.dat100,"ocean_avg_0007.nc"), subds="temp"), subset=s),
              subset(rast(paste0(data.dat100,"ocean_avg_0008.nc"), subds="temp"), subset=s),
              subset(rast(paste0(data.dat100,"ocean_avg_0009.nc"), subds="temp"), subset=s),
              subset(rast(paste0(data.dat100,"ocean_avg_0010.nc"), subds="temp"), subset=s),
              subset(rast(paste0(data.dat100,"ocean_avg_0011.nc"), subds="temp"), subset=s),
              subset(rast(paste0(data.dat100,"ocean_avg_0012.nc"), subds="temp"), subset=s))
temp <- mean(temp.raw)
## 
sa <- te <- empty.roms.ra
salt.dat <- extract(salt, coord.grd.u, method="bilinear")
temp.dat <- extract(temp, coord.grd.u, method="bilinear")
sa[] <- salt.dat[,1]
te[] <- temp.dat[,1]
# se[] <- extract(settle_08, coordinates(h), method="bilinear")
# su[] <- extract(susp_08, coordinates(h), method="bilinear")
# fl[] <- extract(flux_08, coordinates(h), method="bilinear")

## resample to standard 500m resolution of other environmental variables
t_500 <- resample(te,r)
s_500 <- resample(sa,r)

## shelf only
t_500_shelf <- t_500
s_500_shelf <- s_500

t_500_shelf[is.na(r)] <- NA
s_500_shelf[is.na(r)] <- NA

## write rasters to file
roms4k.dir <- paste0(usr.dropbox.dir,"data_environmental/derived/ROMS/")
writeRaster(te,               overwrite=TRUE, filename=paste0(roms4k.dir,string.chr,"waom4k_seafloortemperature.tif"))
writeRaster(t_500,            overwrite=TRUE, filename=paste0(roms4k.dir,string.chr,"500m_waom4k_seafloortemperature.tif"))
writeRaster(t_500_shelf,      overwrite=TRUE, filename=paste0(roms4k.dir,string.chr,"500m_shelf_waom4k_seafloortemperature.tif"))
writeRaster(sa,               overwrite=TRUE, filename=paste0(roms4k.dir,string.chr,"waom4k_seafloorsalinity.tif"))
writeRaster(s_500,            overwrite=TRUE, filename=paste0(roms4k.dir,string.chr,"500m_waom4k_seafloorsalinity.tif"))
writeRaster(s_500_shelf,      overwrite=TRUE, filename=paste0(roms4k.dir,string.chr,"500m_shelf_waom4k_seafloorsalinity.tif"))

## resample to 2km resolution of other environmental variables
t_2k <-       resample(te,r2k.depth)
s_2k <-       resample(sa,r2k.depth)

## shelf only
t_2k_shelf <- t_2k
s_2k_shelf <- s_2k
t_2k_shelf[is.na(r2k.depth)] <- NA
s_2k_shelf[is.na(r2k.depth)] <- NA

## write rasters to file
writeRaster(t_2k,            overwrite=TRUE, filename=paste0(roms4k.dir,string.chr,"2km_waom4k_seafloortemperature.tif"))
writeRaster(t_2k_shelf,      overwrite=TRUE, filename=paste0(roms4k.dir,string.chr,"2km_shelf_waom4k_seafloortemperature.tif"))
writeRaster(s_2k,            overwrite=TRUE, filename=paste0(roms4k.dir,string.chr,"2km_waom4k_seafloorsalinity.tif"))
writeRaster(s_2k_shelf,      overwrite=TRUE, filename=paste0(roms4k.dir,string.chr,"2km_shelf_waom4k_seafloorsalinity.tif"))




