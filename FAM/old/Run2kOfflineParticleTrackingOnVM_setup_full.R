####
## This script:
## - reads in the original ROMS history files ("OcVel_his_...")
## - shift u- and v-points to rho-points and extract current velocities, save output temporarily (huge storage requirements)
## - setup base ROMS object
## - combine current files to 24h (daily) files, surface and seafloor layers only
## - split 24h files into regions
## - reduce 24h files to 6hourly outputs
####

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
#### ROMS Currents
### WE NEED A RUN ACROSS SUMMER WITH HIGH-RES HISTORY FILES!!! 8 has one month only, 8_Nov has 10-day intervals...
### load lon/lat information from ROMS-grid

## tidync code:
grd2k_nc <- tidync(paste0(env.dir3,"waom2extend_grd.nc"))
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
    writeRaster(u.ra.loop, file=paste0(env.dir3,"ocean_his_00",dgts,"_slices_u_",k,".tif"), overwrite=TRUE)
    ## save output data for each time-slice in an array so we can use them in the particle tracking
    # u.vals.array <- as.array(u.ra.loop)
    # save(u.vals.array, file=paste0(env.dir3,"ocean_his_00",dgts,"_slices_u_",k,".Rdata"))
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
    writeRaster(v.ra.loop, file=paste0(env.dir3,"ocean_his_00",dgts,"_slices_v_",k,".tif"), overwrite=TRUE)
    ## save output data for each time-slice in an array so we can use them in the particle tracking
    # v.vals.array <- as.array(v.ra.loop)
    # save(v.vals.array, file=paste0(env.dir3,"ocean_his_00",dgts,"_slices_v_",k,".Rdata"))
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
    writeRaster(w.ra.loop, file=paste0(env.dir3,"ocean_his_00",dgts,"_slices_w_",k,".tif"), overwrite=TRUE)
    ## save output data for each time-slice in an array so we can use them in the particle tracking
    # w.vals.array <- as.array(w.ra.loop)
    # save(w.vals.array, file=paste0(env.dir3,"ocean_his_00",dgts,"_slices_w_",k,".Rdata"))
    # w.ra.t <- setValues(empty.roms.ra31, w.vals.array)
    # plot(w.ra.t$lyr.1)
}}


################################
#### setup base ROMS object
# nc_2k <- nc_open(paste0(env.dir3,"ocean_his_0003.nc"))
nc_2k <- nc_open(paste0(env.dir3,"OcVel_his_0001.nc"))
s_rho <- ncvar_get(nc_2k, varid="s_rho") ## 31             - vertical levels varying between 0 to 1
Cs_r  <- ncvar_get(nc_2k, varid="Cs_r")  ## 31             - "S-coordinate stretching curves at RHO-points
hc    <- ncvar_get(nc_2k, varid="hc")    ##                - value of 250
h     <- ncvar_get(nc_2k, varid="h")     ## 3150 2800      - model bathymetry
zeta  <- ncvar_get(nc_2k, varid="zeta")  ## 3150 2800 24    - free surface elevation
zice  <- ncvar_get(nc_2k, varid="zice")  ## 3150 2800      - ice draft
nc_close(nc_2k)

# nc_2k <- tidync(paste0(env.dir3,"OcVel_his_0001.nc"))
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
save(ROMS_2k_base, file=paste0(env.dir3,"ocean_his_TrackingSetup_base.Rdata"))
rm(ROMS_2k_base)

###################################
#### Current speeds for each day:
### compile into bottom and surface speed tifs:

# curr.files.u <- list.files(env.dir3, pattern="slices_u", full.names=TRUE)
# u.ra.bottom  <- rast(curr.files.u)[[seq(1,488, by=2)]]
# u.ra.top     <- rast(curr.files.u)[[seq(2,488, by=2)]]
# curr.files.v <- list.files(env.dir3, pattern="slices_v", full.names=TRUE)
# v.ra.bottom  <- rast(curr.files.v)[[seq(1,488, by=2)]]
# v.ra.top     <- rast(curr.files.v)[[seq(2,488, by=2)]]
# curr.files.w <- list.files(env.dir3, pattern="slices_w", full.names=TRUE)
# w.ra.bottom  <- rast(curr.files.w)[[seq(1,488, by=2)]]
# w.ra.top     <- rast(curr.files.w)[[seq(2,488, by=2)]]
# ### calculate mean, max, and abs.mean current speeds
# ##
# u.ra.bottom.abs1 <- abs(u.ra.bottom[[1:61]])## higher values here will be due to tidal components
# u.ra.bottom.abs2 <- abs(u.ra.bottom[[62:122]])
# u.ra.bottom.abs3 <- abs(u.ra.bottom[[123:183]])
# u.ra.bottom.abs4 <- abs(u.ra.bottom[[184:244]])
# u.ra.bottom.abs <- c(u.ra.bottom.abs1,u.ra.bottom.abs2,u.ra.bottom.abs3,u.ra.bottom.abs4)
# u.ra.bottom.abs.mean <- mean(u.ra.bottom.abs)
# writeRaster(u.ra.bottom.abs.mean, file=paste0(env.dir3,"ocean_his_bottom_u_absmean.tif"), overwrite=TRUE)
# u.ra.bottom.mean <- mean(u.ra.bottom)
# writeRaster(u.ra.bottom.mean, file=paste0(env.dir3,"ocean_his_bottom_u_mean.tif"), overwrite=TRUE)
# u.ra.bottom.max <- max(u.ra.bottom)
# writeRaster(u.ra.bottom.max, file=paste0(env.dir3,"ocean_his_bottom_u_max.tif"), overwrite=TRUE)
# rm(u.ra.bottom.abs1,u.ra.bottom.abs2,u.ra.bottom.abs3,u.ra.bottom.abs4,v.ra.bottom.abs,u.ra.bottom.abs.mean,u.ra.bottom.mean,u.ra.bottom.max)
# ##
# v.ra.bottom.abs1 <- abs(v.ra.bottom[[1:61]])
# v.ra.bottom.abs2 <- abs(v.ra.bottom[[62:122]])
# v.ra.bottom.abs3 <- abs(v.ra.bottom[[123:183]])
# v.ra.bottom.abs4 <- abs(v.ra.bottom[[184:244]])
# v.ra.bottom.abs <- c(v.ra.bottom.abs1,v.ra.bottom.abs2,v.ra.bottom.abs3,v.ra.bottom.abs4)
# v.ra.bottom.abs.mean <- mean(v.ra.bottom.abs)
# writeRaster(v.ra.bottom.abs.mean, file=paste0(env.dir3,"ocean_his_bottom_v_absmean.tif"), overwrite=TRUE)
# v.ra.bottom.mean <- mean(v.ra.bottom)
# writeRaster(v.ra.bottom.mean, file=paste0(env.dir3,"ocean_his_bottom_v_mean.tif"), overwrite=TRUE)
# v.ra.bottom.max <- max(v.ra.bottom)
# writeRaster(v.ra.bottom.max, file=paste0(env.dir3,"ocean_his_bottom_v_max.tif"), overwrite=TRUE)
# rm(v.ra.bottom.abs1,v.ra.bottom.abs2,v.ra.bottom.abs3,v.ra.bottom.abs4,v.ra.bottom.abs,v.ra.bottom.abs.mean,v.ra.bottom.mean)
# ##
# w.ra.bottom.abs1 <- abs(w.ra.bottom[[1:61]])
# w.ra.bottom.abs2 <- abs(w.ra.bottom[[62:122]])
# w.ra.bottom.abs3 <- abs(w.ra.bottom[[123:183]])
# w.ra.bottom.abs4 <- abs(w.ra.bottom[[184:244]])
# w.ra.bottom.abs <- c(w.ra.bottom.abs1,w.ra.bottom.abs2,w.ra.bottom.abs3,w.ra.bottom.abs4)
# w.ra.bottom.abs.mean <- mean(w.ra.bottom.abs)
# writeRaster(w.ra.bottom.abs.mean, file=paste0(env.dir3,"ocean_his_bottom_w_absmean.tif"), overwrite=TRUE)
# w.ra.bottom.mean <- mean(w.ra.bottom)
# writeRaster(w.ra.bottom.mean, file=paste0(env.dir3,"ocean_his_bottom_w_mean.tif"), overwrite=TRUE)
# w.ra.bottom.max <- max(w.ra.bottom)
# writeRaster(w.ra.bottom.max, file=paste0(env.dir3,"ocean_his_bottom_w_max.tif"), overwrite=TRUE)
# rm(w.ra.bottom.abs1,w.ra.bottom.abs2,w.ra.bottom.abs3,w.ra.bottom.abs4,w.ra.bottom.abs.mean,w.ra.bottom.mean)
# ##
# u.ra.bottom.mean <- rast(paste0(env.dir3,"ocean_his_bottom_u_mean.tif"))
# v.ra.bottom.mean <- rast(paste0(env.dir3,"ocean_his_bottom_v_mean.tif"))
# uv.ra.bottom.mean <- sqrt(u.ra.bottom.mean^2+v.ra.bottom.mean^2)
# writeRaster(uv.ra.bottom.mean, file=paste0(env.dir3,"ocean_his_bottom_uv_mean.tif"), overwrite=TRUE)
# ##
# u.ra.bottom.abs.mean <- rast(paste0(env.dir3,"ocean_his_bottom_u_absmean.tif"))
# v.ra.bottom.abs.mean <- rast(paste0(env.dir3,"ocean_his_bottom_v_absmean.tif"))
# uv.ra.bottom.abs.mean <- sqrt(u.ra.bottom.abs.mean^2+v.ra.bottom.abs.mean^2)
# writeRaster(uv.ra.bottom.abs.mean, file=paste0(env.dir3,"ocean_his_bottom_uv_absmean.tif"), overwrite=TRUE)
# v.ra.bottom.max <- rast(paste0(env.dir3,"ocean_his_bottom_v_max.tif"))
# u.ra.bottom.max <- rast(paste0(env.dir3,"ocean_his_bottom_u_max.tif"))
# uv.ra.bottom.max <- sqrt(u.ra.bottom.max^2+v.ra.bottom.max^2)
# writeRaster(uv.ra.bottom.max, file=paste0(env.dir3,"ocean_his_bottom_uv_max.tif"), overwrite=TRUE)
##




### current speeds: 8 slices at 3h give a 24 compilation
ROMS_2k <- list()
ROMS_2k$i_u <- array(NA, c(2800,3150,31,9))
ROMS_2k$i_v <- array(NA, c(2800,3150,31,9))
ROMS_2k$i_w <- array(NA, c(2800,3150,31,9))
for(j in 1:9){
  print(j)
  ROMS_2k$i_u[,,,j] <- as.array(rast(paste0(env.dir3,"ocean_his_0001_slices_u_",j,".tif")))
  ROMS_2k$i_v[,,,j] <- as.array(rast(paste0(env.dir3,"ocean_his_0001_slices_v_",j,".tif")))
  ROMS_2k$i_w[,,,j] <- as.array(rast(paste0(env.dir3,"ocean_his_0001_slices_w_",j,".tif")))
}
save(ROMS_2k, file=paste0(env.dir3,"ocean_his_TrackingSetup_24hcurrents_File01.Rdata"))

### current speeds: 8 slices at 3h give a 24 compilation
#for(i in 2:30){
for(i in 29){
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
  save(ROMS_2k, file=paste0(env.dir3,"ocean_his_TrackingSetup_24hcurrents_File",dgts,".Rdata"))
}


roms.coords.proj <- cbind(c(ROMS_2k_base$x), c(ROMS_2k_base$y))
x.range <- c(min(roms.coords.proj[,1])-1000,max(roms.coords.proj[,1])+1000)
y.range <- c(min(roms.coords.proj[,2])-1000,max(roms.coords.proj[,2])+1000)
empty.roms.ra <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000)
h.ra      <- setValues(empty.roms.ra, ROMS_2k_base$h)
u.ra      <- setValues(empty.roms.ra, ROMS_2k$i_u[,,1,1])



#################################
#### split into regional datasets
load(paste0(env.dir3,"ocean_his_TrackingSetup_24hcurrents_File01.Rdata"))
load(paste0(env.dir3,"ocean_his_TrackingSetup_base.Rdata"))

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
  save(Rdat, file=paste0(env.dir3,"ocean_his_TrackingSetup_24hcurrents_File01_Region",sprintf("%02d",i),".Rdata"))
}

#### split into regional datasets
for(k in 29:30){
  message(k)
  load(paste0(env.dir3,"ocean_his_TrackingSetup_24hcurrents_File",sprintf("%02d",k),".Rdata"))
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
    save(Rdat, file=paste0(env.dir3,"ocean_his_TrackingSetup_24hcurrents_File",sprintf("%02d",k),"_Region",sprintf("%02d",i),".Rdata"))
  }
}


#### combine current files into one for tracking
### the full 28 days crashes the VM
# load(paste0(env.dir3,"ocean_his_TrackingSetup_24hcurrents_Region01.Rdata"))
# Rdat.raw <- Rdat
# load(paste0(env.dir3,"ocean_his_TrackingSetup_24hcurrents_Region01_b.Rdata"))
# #load(paste0(env.dir3,"ocean_his_TrackingSetup_24hcurrents_Region01_c.Rdata"))
# Rdat.raw$i_u <- abind(Rdat.raw$i_u, Rdat$i_u)
# Rdat.raw$i_v <- abind(Rdat.raw$i_v, Rdat$i_v)
# Rdat.raw$i_w <- abind(Rdat.raw$i_w, Rdat$i_w)
# Rdat <- Rdat.raw
# rm(Rdat.raw)
# save(Rdat, file=paste0(env.dir3,"ocean_his_TrackingSetup_24hcurrents_Region01a.Rdata"))
### reducing to 6h outputs for 28 days for the 3D tracking
### 3h outputs for the 2D tracking should work!!!
library(abind)
for(k in 7:10){
  message(k)
  load(paste0(env.dir3,"ocean_his_TrackingSetup_24hcurrents_File01_Region",sprintf("%02d",k),".Rdata"))
  Rdat6h <- Rdat
  sel <- seq(1, 9, by=2)
  Rdat6h$i_u <- Rdat6h$i_u[,,,sel]
  Rdat6h$i_v <- Rdat6h$i_v[,,,sel]
  Rdat6h$i_w <- Rdat6h$i_w[,,,sel]
  sel <- seq(2, 8, by=2)
  for(i in 2:28){
    print(i)
    load(paste0(env.dir3,"ocean_his_TrackingSetup_24hcurrents_File",sprintf("%02d",i),"_Region",sprintf("%02d",k),".Rdata"))
    Rdat6h$i_u <- abind(Rdat6h$i_u, Rdat$i_u[,,,sel])
    Rdat6h$i_v <- abind(Rdat6h$i_v, Rdat$i_v[,,,sel])
    Rdat6h$i_w <- abind(Rdat6h$i_w, Rdat$i_w[,,,sel])
  }
  save(Rdat6h, file=paste0(env.dir3,"ocean_his_TrackingSetup_6hourlycurrents_28days_Region",sprintf("%02d",k),".Rdata"))
}

#### region 6 is too big to do it all at once:
load(paste0(env.dir3,"ocean_his_TrackingSetup_24hcurrents_File01_Region06.Rdata"))
## a
Rdat6h <- Rdat
sel <- seq(1, 9, by=2)
Rdat6h$i_u <- Rdat6h$i_u[,,,sel]
Rdat6h$i_v <- Rdat6h$i_v[,,,sel]
Rdat6h$i_w <- Rdat6h$i_w[,,,sel]
sel <- seq(2, 8, by=2)
for(i in 2:10){
    print(i)
    load(paste0(env.dir3,"ocean_his_TrackingSetup_24hcurrents_File",sprintf("%02d",i),"_Region06.Rdata"))
    Rdat6h$i_u <- abind(Rdat6h$i_u, Rdat$i_u[,,,sel])
    Rdat6h$i_v <- abind(Rdat6h$i_v, Rdat$i_v[,,,sel])
    Rdat6h$i_w <- abind(Rdat6h$i_w, Rdat$i_w[,,,sel])
}
save(Rdat6h, file=paste0(env.dir3,"ocean_his_TrackingSetup_6hourlycurrents_28days_Region06a.Rdata"))
## b
load(paste0(env.dir3,"ocean_his_TrackingSetup_24hcurrents_File11_Region06.Rdata"))
Rdat6h <- Rdat
sel <- seq(2, 8, by=2)
Rdat6h$i_u <- Rdat6h$i_u[,,,sel]
Rdat6h$i_v <- Rdat6h$i_v[,,,sel]
Rdat6h$i_w <- Rdat6h$i_w[,,,sel]
for(i in 12:20){
  print(i)
  load(paste0(env.dir3,"ocean_his_TrackingSetup_24hcurrents_File",sprintf("%02d",i),"_Region06.Rdata"))
  Rdat6h$i_u <- abind(Rdat6h$i_u, Rdat$i_u[,,,sel])
  Rdat6h$i_v <- abind(Rdat6h$i_v, Rdat$i_v[,,,sel])
  Rdat6h$i_w <- abind(Rdat6h$i_w, Rdat$i_w[,,,sel])
}
save(Rdat6h, file=paste0(env.dir3,"ocean_his_TrackingSetup_6hourlycurrents_28days_Region06b.Rdata"))
## c
load(paste0(env.dir3,"ocean_his_TrackingSetup_24hcurrents_File21_Region06.Rdata"))
Rdat6h <- Rdat
sel <- seq(2, 8, by=2)
Rdat6h$i_u <- Rdat6h$i_u[,,,sel]
Rdat6h$i_v <- Rdat6h$i_v[,,,sel]
Rdat6h$i_w <- Rdat6h$i_w[,,,sel]
for(i in 22:28){
  print(i)
  load(paste0(env.dir3,"ocean_his_TrackingSetup_24hcurrents_File",sprintf("%02d",i),"_Region06.Rdata"))
  Rdat6h$i_u <- abind(Rdat6h$i_u, Rdat$i_u[,,,sel])
  Rdat6h$i_v <- abind(Rdat6h$i_v, Rdat$i_v[,,,sel])
  Rdat6h$i_w <- abind(Rdat6h$i_w, Rdat$i_w[,,,sel])
}
save(Rdat6h, file=paste0(env.dir3,"ocean_his_TrackingSetup_6hourlycurrents_28days_Region06c.Rdata"))
## bind together
load(paste0(env.dir3,"ocean_his_TrackingSetup_6hourlycurrents_28days_Region06c.Rdata"))
c <- Rdat6h
load(paste0(env.dir3,"ocean_his_TrackingSetup_6hourlycurrents_28days_Region06b.Rdata"))
Rdat6h$i_u <- abind(Rdat6h$i_u, c$i_u)
Rdat6h$i_v <- abind(Rdat6h$i_v, c$i_v)
Rdat6h$i_w <- abind(Rdat6h$i_w, c$i_w)
save(Rdat6h, file=paste0(env.dir3,"ocean_his_TrackingSetup_6hourlycurrents_28days_Region06bc.Rdata"))
## step by ...
load(paste0(env.dir3,"ocean_his_TrackingSetup_6hourlycurrents_28days_Region06bc.Rdata"))
b <- Rdat6h
load(paste0(env.dir3,"ocean_his_TrackingSetup_6hourlycurrents_28days_Region06a.Rdata"))
Rdat6h$i_u <- abind(Rdat6h$i_u, b$i_u)
Rdat6h$i_v <- abind(Rdat6h$i_v, b$i_v)
save(Rdat6h, file=paste0(env.dir3,"ocean_his_TrackingSetup_6hourlycurrents_28days_Region06abcuv.Rdata"))
## ... step
load(paste0(env.dir3,"ocean_his_TrackingSetup_6hourlycurrents_28days_Region06bc.Rdata"))
b <- Rdat6h
rm(Rdat6h)
b$i_u <- NA
b$i_v <- NA
b$hh <- NA
load(paste0(env.dir3,"ocean_his_TrackingSetup_6hourlycurrents_28days_Region06abcuv.Rdata"))
Rdat6h$i_w <- abind(Rdat6h$i_w, b$i_w)
save(Rdat6h, file=paste0(env.dir3,"ocean_his_TrackingSetup_6hourlycurrents_28days_Region06.Rdata"))
##




#### we only need the bottom few layers for the 2D tracking

#### Calculate average current speeds
## THE BELOW CODE RESULTS IN TIF-FILES WITH VARYING RESOLUTION! NOT GOOD...

# ### Full files already deleted, so rather than producing them again, taking the regional files an dmerging the data
# ### only the bottom current layer
# base.str <- "/pvol3TB/data_environmental/ROMS_2k_files/ocean_his_TrackingSetup_24hcurrents_"
# ## Function to load data files and extract relevant objects
# load_data <- function(region, index) {
#   file_name <- sprintf(paste0(base.str,"File%02d_Region%02d.Rdata"), index, region)
#   load(file_name)
#   list(x = Rdat$x, y = Rdat$y, i_u = Rdat$i_u[,,1,], i_v = Rdat$i_v[,,1,], i_w = Rdat$i_w[,,1,])
# }
# ## Function to calculate mean, max, and absolute mean across the third dimension
# calculate_stats <- function(combined_i_u, combined_i_v, combined_i_w) {
#   mean_i_u <- apply(combined_i_u, c(1, 2), mean, na.rm = TRUE)
#   max_i_u <- apply(combined_i_u, c(1, 2), function(x) ifelse(all(is.na(x)), NA, max(x, na.rm = TRUE)))
#   abs_mean_i_u <- apply(abs(combined_i_u), c(1, 2), mean, na.rm = TRUE)
#   
#   mean_i_v <- apply(combined_i_v, c(1, 2), mean, na.rm = TRUE)
#   max_i_v <- apply(combined_i_v, c(1, 2), function(x) ifelse(all(is.na(x)), NA, max(x, na.rm = TRUE)))
#   abs_mean_i_v <- apply(abs(combined_i_v), c(1, 2), mean, na.rm = TRUE)
#   
#   mean_i_w <- apply(combined_i_w, c(1, 2), mean, na.rm = TRUE)
#   max_i_w <- apply(combined_i_w, c(1, 2), function(x) ifelse(all(is.na(x)), NA, max(x, na.rm = TRUE)))
#   abs_mean_i_w <- apply(abs(combined_i_w), c(1, 2), mean, na.rm = TRUE)
#   
#   list(
#     mean_i_u = mean_i_u, max_i_u = max_i_u, abs_mean_i_u = abs_mean_i_u,
#     mean_i_v = mean_i_v, max_i_v = max_i_v, abs_mean_i_v = abs_mean_i_v,
#     mean_i_w = mean_i_w, max_i_w = max_i_w, abs_mean_i_w = abs_mean_i_w
#   )
# }
# ## Process each region individually and save results
# for (region in 1:10) {
#   # Load one file to get coordinates
#   data <- load_data(region, 1)
#   x_coords <- data$x
#   y_coords <- data$y
#   
#   combined_i_u <- list()
#   combined_i_v <- list()
#   combined_i_w <- list()
# 
#   for (index in 1:30) {
#     print(index)
#     data <- load_data(region, index)
#     combined_i_u <- c(combined_i_u, list(data$i_u))
#     combined_i_v <- c(combined_i_v, list(data$i_v))
#     combined_i_w <- c(combined_i_w, list(data$i_w))
#   }
#   print("data loaded for the region")
#   ## Combine the lists into single arrays
#   combined_i_u <- do.call(abind, c(combined_i_u, along = 3))
#   combined_i_v <- do.call(abind, c(combined_i_v, along = 3))
#   combined_i_w <- do.call(abind, c(combined_i_w, along = 3))
#   ## Calculate statistics
#   stats <- calculate_stats(combined_i_u, combined_i_v, combined_i_w)
#   print("stats calculated")
#   save(stats, file = sprintf(paste0(base.str,"Region%02d_bottomspeedstats.Rdata"), region))
#   cat("Processed and saved results for Region", region, "\n")
#   
#   ## Calculate speed from mean, max, and absolute mean values
#   mean_speed <- sqrt(stats$mean_i_u^2 + stats$mean_i_v^2)
#   max_speed  <- sqrt(stats$max_i_u^2 + stats$max_i_v^2)
#   abs_mean_speed <- sqrt(stats$abs_mean_i_u^2 + stats$abs_mean_i_v^2)
#   
#   ## Create SpatRaster objects for mean, max, and absolute mean values
#   mean_raster <- rast(nrows = nrow(mean_speed), ncols = ncol(mean_speed), 
#                       xmin = min(x_coords), xmax = max(x_coords), 
#                       ymin = min(y_coords), ymax = max(y_coords), crs = stereo)
#   values(mean_raster) <- mean_speed
#   max_raster <- rast(nrows = nrow(max_speed), ncols = ncol(max_speed), 
#                      xmin = min(x_coords), xmax = max(x_coords), 
#                      ymin = min(y_coords), ymax = max(y_coords), crs = stereo)
#   values(max_raster) <- max_speed
#   abs_mean_raster <- rast(nrows = nrow(abs_mean_speed), ncols = ncol(abs_mean_speed), 
#                           xmin = min(x_coords), xmax = max(x_coords), 
#                           ymin = min(y_coords), ymax = max(y_coords), crs = stereo)
#   values(abs_mean_raster) <- abs_mean_speed
# 
#   ## Save the rasters as GeoTIFF files
#   writeRaster(mean_raster,     sprintf(paste0(base.str,"Region%02d_mean.tif"), region), overwrite = TRUE)
#   writeRaster(max_raster,      sprintf(paste0(base.str,"Region%02d_max.tif"), region), overwrite = TRUE)
#   writeRaster(abs_mean_raster, sprintf(paste0(base.str,"Region%02d_abs_mean.tif"), region), overwrite = TRUE)
#   
#   cat("Processed and saved GeoTIFF files for Region", region, "\n")
# }
# ## then merge the 10 regions
# region=1
# load(sprintf(paste0(base.str,"Region%02d_bottomspeedstats.Rdata"), region))
# 
# ra1 <- rast(sprintf(paste0(base.str,"Region%02d_mean.tif"), 1))
# ra2 <- rast(sprintf(paste0(base.str,"Region%02d_mean.tif"), 2))
# ra3 <- rast(sprintf(paste0(base.str,"Region%02d_mean.tif"), 3))
# ra4 <- rast(sprintf(paste0(base.str,"Region%02d_mean.tif"), 4))
# ra5 <- rast(sprintf(paste0(base.str,"Region%02d_mean.tif"), 5))
# ra6 <- rast(sprintf(paste0(base.str,"Region%02d_mean.tif"), 6))
# ra7 <- rast(sprintf(paste0(base.str,"Region%02d_mean.tif"), 7))
# ra8 <- rast(sprintf(paste0(base.str,"Region%02d_mean.tif"), 8))
# ra9 <- rast(sprintf(paste0(base.str,"Region%02d_mean.tif"), 9))
# ra10 <- rast(sprintf(paste0(base.str,"Region%02d_mean.tif"), 10))
# 
# merge_antarctic_rasters <- function(base.str, regions = 1:10) {
#   # Create an empty list to store rasters
#   raster_list <- list()
#   
#   # Read all raster files
#   for (region in regions) {
#     filename <- sprintf(paste0(base.str, "Region%02d_mean.tif"), region)
#     rast_obj <- rast(filename)
#     
#     # Ensure projection is correct
#     stereo <- "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
#     if (crs(rast_obj) != stereo) {
#       rast_obj <- project(rast_obj, stereo)
#     }
#     
#     raster_list[[region]] <- rast_obj
#   }
#   
#   # Merge rasters with mean for overlapping areas
#   merged_raster <- mosaic(sprc(raster_list), fun = "mean")
#   
#   return(merged_raster)
# }
# 
# merged_currents <- merge_antarctic_rasters(base.str)





#################################################################
#### create and save npp-based seeding file for each region
## Options:
## - 1 particle with individualised value per cell: might get cells with 0 particles -> not good
## - 4 particles with individualised value per cell
## - 9  particles with individualised value per cell
## - make the # of particles dependent on cell value UNFEASIBLE TO RUN, because it needs 50 million particles for the particle distribution to be closely related to the NPP values in the Ross Sea
library(spatstat)
library(ppmData) ## for "terra2im"
npp <- rast("/pvol/data_environmental/Circumpolar_EnvData_2km_NPP_Cafe_filled_SummerAverage.tif")
ra <- rast(paste0(env.dir,"ocean_his_0001_slices_u_1.tif"), lyrs=31)
npp2 <- project(npp, ra)

#### - 9  particles with individualised value per cell
## load ROMS and NPP data
load(paste0(env.dir3,"ocean_his_TrackingSetup_24hcurrents_File01_Region01.Rdata"))
roms.coords.proj <- cbind(c(Rdat$x), c(Rdat$y))
x.range <- round(c(min(roms.coords.proj[,1])-1000,max(roms.coords.proj[,1])+1000),0)
y.range <- round(c(min(roms.coords.proj[,2])-1000,max(roms.coords.proj[,2])+1000),0)
ra.region <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000)
## particles locations
ra.9 <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000/3)
ra.4 <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000/2)
pts9.raw <- data.frame(cbind(crds(ra.9),z=0))
pts4.raw <- data.frame(cbind(crds(ra.4),z=0))
## crop npp to boundaries, and set irrelevant regions to NA
npp2_crop <- crop(npp2, ext(ra.region))
npp2_crop[400:650,1:100] <- NA
npp2_crop[500:650,1:150] <- NA
h.ra <- setValues(ra.region, Rdat$h)
plot(npp2_crop)
contour(h.ra, add=TRUE, levels=c(0,2000,3000))

## particle values
npp.vals9 <- extract(npp2_crop, pts9.raw[,1:2])
pts9.raw$npp <- npp.vals9$Circumpolar_EnvData_2km_NPP_Cafe_filled_SummerAverage/9
pts9 <- as.matrix(pts9.raw[-which(is.na(pts9.raw$npp)),])
npp.vals4 <- extract(npp2_crop, pts4.raw[,1:2])
pts4.raw$npp <- npp.vals4$Circumpolar_EnvData_2km_NPP_Cafe_filled_SummerAverage/4
pts4 <- as.matrix(pts4.raw[-which(is.na(pts4.raw$npp)),])

save(pts4, file=paste0(env.dir3,"ocean_his_TrackingSetup_NPP4_Region01.Rdata"))
save(pts9, file=paste0(env.dir3,"ocean_his_TrackingSetup_NPP9_Region01.Rdata"))

## save input rasters and boundaries for comparison later 
writeRaster(npp2_crop, file=paste0(env.dir3,"ocean_his_TrackingSetup_NPP_Region01.tif"))

#### SETUP NPP9 FOR ALL REGIONS

npp <- rast("/pvol/data_environmental/Circumpolar_EnvData_2km_NPP_Cafe_filled_SummerAverage.tif")
ra <- rast(paste0(env.dir,"ocean_his_0001_slices_u_1.tif"), lyrs=31)
npp2 <- project(npp, ra)

#### region 2
load(paste0(env.dir3,"ocean_his_TrackingSetup_24hcurrents_File01_Region02.Rdata"))
roms.coords.proj <- cbind(c(Rdat$x), c(Rdat$y))
x.range <- round(c(min(roms.coords.proj[,1])-1000,max(roms.coords.proj[,1])+1000),0)
y.range <- round(c(min(roms.coords.proj[,2])-1000,max(roms.coords.proj[,2])+1000),0)
ra.region <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000)
## particles locations
ra.9 <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000/3)
pts9.raw <- data.frame(cbind(crds(ra.9),z=0))
npp2_crop <- crop(npp2, ext(ra.region))
## remove unnecessary particles and check distribution
npp2_crop[330:400,c(1:150,320:800)] <- NA
npp2_crop[1:180,1:70] <- NA
## particle values
npp.vals9 <- extract(npp2_crop, pts9.raw[,1:2])
pts9.raw$npp <- npp.vals9$Circumpolar_EnvData_2km_NPP_Cafe_filled_SummerAverage/9
pts9 <- as.matrix(pts9.raw[-which(is.na(pts9.raw$npp)),])
save(pts9, file=paste0(env.dir3,"ocean_his_TrackingSetup_NPP9_Region02.Rdata"))
writeRaster(npp2_crop, file=paste0(env.dir3,"ocean_his_TrackingSetup_NPP_Region02.tif"))

#### region 3
load(paste0(env.dir3,"ocean_his_TrackingSetup_24hcurrents_File01_Region03.Rdata"))
roms.coords.proj <- cbind(c(Rdat$x), c(Rdat$y))
x.range <- round(c(min(roms.coords.proj[,1])-1000,max(roms.coords.proj[,1])+1000),0)
y.range <- round(c(min(roms.coords.proj[,2])-1000,max(roms.coords.proj[,2])+1000),0)
ra.region <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000)
## particles locations
ra.9 <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000/3)
pts9.raw <- data.frame(cbind(crds(ra.9),z=0))
npp2_crop <- crop(npp2, ext(ra.region))
## remove unnecessary particles and check distribution
npp2_crop[150:950,1:50] <- NA
npp2_crop[530:950,1:170] <- NA
npp2_crop[800:950,1:250] <- NA
## particle values
npp.vals9 <- extract(npp2_crop, pts9.raw[,1:2])
pts9.raw$npp <- npp.vals9$Circumpolar_EnvData_2km_NPP_Cafe_filled_SummerAverage/9
pts9 <- as.matrix(pts9.raw[-which(is.na(pts9.raw$npp)),])
save(pts9, file=paste0(env.dir3,"ocean_his_TrackingSetup_NPP9_Region03.Rdata"))
writeRaster(npp2_crop, file=paste0(env.dir3,"ocean_his_TrackingSetup_NPP_Region03.tif"))

#### region 4
load(paste0(env.dir3,"ocean_his_TrackingSetup_24hcurrents_File01_Region04.Rdata"))
roms.coords.proj <- cbind(c(Rdat$x), c(Rdat$y))
x.range <- round(c(min(roms.coords.proj[,1])-1000,max(roms.coords.proj[,1])+1000),0)
y.range <- round(c(min(roms.coords.proj[,2])-1000,max(roms.coords.proj[,2])+1000),0)
ra.region <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000)
## particles locations
ra.9 <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000/3)
pts9.raw <- data.frame(cbind(crds(ra.9),z=0))
npp2_crop <- crop(npp2, ext(ra.region))
## remove unnecessary particles and check distribution
npp2_crop[1:900,1:30] <- NA
npp2_crop[c(1:350,700:900),1:70] <- NA
npp2_crop[c(1:250,800:900),1:100] <- NA
npp2_crop[c(1:200,850:900),1:130] <- NA
npp2_crop[1:150,1:200] <- NA
npp2_crop[1:100,1:250] <- NA
npp2_crop[1:50,1:300] <- NA
## particle values
npp.vals9 <- extract(npp2_crop, pts9.raw[,1:2])
pts9.raw$npp <- npp.vals9$Circumpolar_EnvData_2km_NPP_Cafe_filled_SummerAverage/9
pts9 <- as.matrix(pts9.raw[-which(is.na(pts9.raw$npp)),])
save(pts9, file=paste0(env.dir3,"ocean_his_TrackingSetup_NPP9_Region04.Rdata"))
writeRaster(npp2_crop, file=paste0(env.dir3,"ocean_his_TrackingSetup_NPP_Region04.tif"))

#### region 5
load(paste0(env.dir3,"ocean_his_TrackingSetup_24hcurrents_File01_Region05.Rdata"))
roms.coords.proj <- cbind(c(Rdat$x), c(Rdat$y))
x.range <- round(c(min(roms.coords.proj[,1])-1000,max(roms.coords.proj[,1])+1000),0)
y.range <- round(c(min(roms.coords.proj[,2])-1000,max(roms.coords.proj[,2])+1000),0)
ra.region <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000)
## particles locations
ra.9 <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000/3)
pts9.raw <- data.frame(cbind(crds(ra.9),z=0))
npp2_crop <- crop(npp2, ext(ra.region))
## remove unnecessary particles and check distribution
npp2_crop[400:600,1:100] <- NA
npp2_crop[1:60,50:700] <- NA
npp2_crop[1:100,100:650] <- NA
npp2_crop[1:130,150:550] <- NA
npp2_crop[1:180,200:350] <- NA
## particle values
npp.vals9 <- extract(npp2_crop, pts9.raw[,1:2])
pts9.raw$npp <- npp.vals9$Circumpolar_EnvData_2km_NPP_Cafe_filled_SummerAverage/9
pts9 <- as.matrix(pts9.raw[-which(is.na(pts9.raw$npp)),])
save(pts9, file=paste0(env.dir3,"ocean_his_TrackingSetup_NPP9_Region05.Rdata"))
writeRaster(npp2_crop, file=paste0(env.dir3,"ocean_his_TrackingSetup_NPP_Region05.tif"))

#### region 6
load(paste0(env.dir3,"ocean_his_TrackingSetup_24hcurrents_File01_Region06.Rdata"))
roms.coords.proj <- cbind(c(Rdat$x), c(Rdat$y))
x.range <- round(c(min(roms.coords.proj[,1])-1000,max(roms.coords.proj[,1])+1000),0)
y.range <- round(c(min(roms.coords.proj[,2])-1000,max(roms.coords.proj[,2])+1000),0)
ra.region <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000)
## particles locations
ra.9 <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000/3)
pts9.raw <- data.frame(cbind(crds(ra.9),z=0))
npp2_crop <- crop(npp2, ext(ra.region))
## remove unnecessary particles and check distribution
npp2_crop[1:100,c(1:550,800:1400)] <- NA
npp2_crop[1:150,c(1:200,800:1100,1250:1400)] <- NA
npp2_crop[1:200,c(1:150,830:1050,1250:1400)] <- NA
npp2_crop[201:250,c(1:100,1200:1400)] <- NA
npp2_crop[251:320,c(1:50,1200:1400)] <- NA
## particle values
npp.vals9 <- extract(npp2_crop, pts9.raw[,1:2])
pts9.raw$npp <- npp.vals9$Circumpolar_EnvData_2km_NPP_Cafe_filled_SummerAverage/9
pts9 <- as.matrix(pts9.raw[-which(is.na(pts9.raw$npp)),])
save(pts9, file=paste0(env.dir3,"ocean_his_TrackingSetup_NPP9_Region06.Rdata"))
writeRaster(npp2_crop, file=paste0(env.dir3,"ocean_his_TrackingSetup_NPP_Region06.tif"))

#### region 7
load(paste0(env.dir3,"ocean_his_TrackingSetup_24hcurrents_File01_Region07.Rdata"))
roms.coords.proj <- cbind(c(Rdat$x), c(Rdat$y))
x.range <- round(c(min(roms.coords.proj[,1])-1000,max(roms.coords.proj[,1])+1000),0)
y.range <- round(c(min(roms.coords.proj[,2])-1000,max(roms.coords.proj[,2])+1000),0)
ra.region <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000)
## particles locations
ra.9 <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000/3)
pts9.raw <- data.frame(cbind(crds(ra.9),z=0))
npp2_crop <- crop(npp2, ext(ra.region))
## remove unnecessary particles and check distribution
npp2_crop[1:100,250:500] <- NA
npp2_crop[1:300,300:500] <- NA
npp2_crop[1:400,400:500] <- NA
npp2_crop[1:700,450:500] <- NA
## particle values
npp.vals9 <- extract(npp2_crop, pts9.raw[,1:2])
pts9.raw$npp <- npp.vals9$Circumpolar_EnvData_2km_NPP_Cafe_filled_SummerAverage/9
pts9 <- as.matrix(pts9.raw[-which(is.na(pts9.raw$npp)),])
save(pts9, file=paste0(env.dir3,"ocean_his_TrackingSetup_NPP9_Region07.Rdata"))
writeRaster(npp2_crop, file=paste0(env.dir3,"ocean_his_TrackingSetup_NPP_Region07.tif"))

#### region 8
load(paste0(env.dir3,"ocean_his_TrackingSetup_24hcurrents_File01_Region08.Rdata"))
roms.coords.proj <- cbind(c(Rdat$x), c(Rdat$y))
x.range <- round(c(min(roms.coords.proj[,1])-1000,max(roms.coords.proj[,1])+1000),0)
y.range <- round(c(min(roms.coords.proj[,2])-1000,max(roms.coords.proj[,2])+1000),0)
ra.region <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000)
## particles locations
ra.9 <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000/3)
pts9.raw <- data.frame(cbind(crds(ra.9),z=0))
npp2_crop <- crop(npp2, ext(ra.region))
## remove unnecessary particles and check distribution
npp2_crop[1:70,430:500] <- NA
npp2_crop[500:900,450:500] <- NA
npp2_crop[600:900,380:500] <- NA
npp2_crop[800:900,250:500] <- NA
## particle values
npp.vals9 <- extract(npp2_crop, pts9.raw[,1:2])
pts9.raw$npp <- npp.vals9$Circumpolar_EnvData_2km_NPP_Cafe_filled_SummerAverage/9
pts9 <- as.matrix(pts9.raw[-which(is.na(pts9.raw$npp)),])
save(pts9, file=paste0(env.dir3,"ocean_his_TrackingSetup_NPP9_Region08.Rdata"))
writeRaster(npp2_crop, file=paste0(env.dir3,"ocean_his_TrackingSetup_NPP_Region08.tif"))

#### region 9
load(paste0(env.dir3,"ocean_his_TrackingSetup_24hcurrents_File01_Region09.Rdata"))
roms.coords.proj <- cbind(c(Rdat$x), c(Rdat$y))
x.range <- round(c(min(roms.coords.proj[,1])-1000,max(roms.coords.proj[,1])+1000),0)
y.range <- round(c(min(roms.coords.proj[,2])-1000,max(roms.coords.proj[,2])+1000),0)
ra.region <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000)
## particles locations
ra.9 <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000/3)
pts9.raw <- data.frame(cbind(crds(ra.9),z=0))
npp2_crop <- crop(npp2, ext(ra.region))
## remove unnecessary particles and check distribution
npp2_crop[500:600,350:800] <- NA
npp2_crop[450:600,500:800] <- NA
npp2_crop[400:600,600:800] <- NA
npp2_crop[350:600,650:800] <- NA
npp2_crop[300:600,700:800] <- NA
npp2_crop[200:600,750:800] <- NA
## particle values
npp.vals9 <- extract(npp2_crop, pts9.raw[,1:2])
pts9.raw$npp <- npp.vals9$Circumpolar_EnvData_2km_NPP_Cafe_filled_SummerAverage/9
pts9 <- as.matrix(pts9.raw[-which(is.na(pts9.raw$npp)),])
save(pts9, file=paste0(env.dir3,"ocean_his_TrackingSetup_NPP9_Region09.Rdata"))
writeRaster(npp2_crop, file=paste0(env.dir3,"ocean_his_TrackingSetup_NPP_Region09.tif"))

#### region 10
load(paste0(env.dir3,"ocean_his_TrackingSetup_24hcurrents_File01_Region10.Rdata"))
roms.coords.proj <- cbind(c(Rdat$x), c(Rdat$y))
x.range <- round(c(min(roms.coords.proj[,1])-1000,max(roms.coords.proj[,1])+1000),0)
y.range <- round(c(min(roms.coords.proj[,2])-1000,max(roms.coords.proj[,2])+1000),0)
ra.region <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000)
## particles locations
ra.9 <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000/3)
pts9.raw <- data.frame(cbind(crds(ra.9),z=0))
npp2_crop <- crop(npp2, ext(ra.region))
## particle values
npp.vals9 <- extract(npp2_crop, pts9.raw[,1:2])
pts9.raw$npp <- npp.vals9$Circumpolar_EnvData_2km_NPP_Cafe_filled_SummerAverage/9
pts9 <- as.matrix(pts9.raw[-which(is.na(pts9.raw$npp)),])
save(pts9, file=paste0(env.dir3,"ocean_his_TrackingSetup_NPP9_Region10.Rdata"))
writeRaster(npp2_crop, file=paste0(env.dir3,"ocean_his_TrackingSetup_NPP_Region10.tif"))

# h.ra <- setValues(ra.region, Rdat$h)
# plot(npp2_crop)
# contour(h.ra, add=TRUE, levels=c(0,2000,3000))















#### UNFEASIBLE: particle numbers correspond to value
# npp2_crop <- crop(npp2, ext(-500000,400000,-2000000, -1200000))
# Z <- terra2im(npp2_crop/10000000) ## need to reduce by at least 10m
# pp <- rpoispp(Z) ## library(spatstat)
# npp.pts <- cbind(pp$x, pp$y,0)
# ## transform pts back to raster, to check how well the points represent the original data
# ra.pts <- rasterize(npp.pts[,1:2], npp2_crop, fun=sum)
# plot(ra.pts[],npp2_crop[], cex=0.1)
# 
# ## dilute points to make it more manageable
# diluted10 <- seq(1,nrow(npp.pts),by=2)
# diluted10 <- sample(1:nrow(npp.pts),4000000)
# diluted10 <- 1:4000000
# npp.pts.red <- npp.pts[diluted10,] ## sampling randomly results in the same pattern as seeding fewer with rpoispp
# ra.pts.red <- rasterize(npp.pts.red[,1:2], npp2_crop, fun=sum)
# plot(ra.pts.red[],ra.pts[], cex=0.1)
# plot(ra.pts.red[],npp2_crop[], cex=0.1)
# for(i in 1:10){
#   ## setup boundaries
#   load(paste0(env.dir3,"ocean_his_TrackingSetup_24hcurrents_File01_Region",sprintf("%02d",i),".Rdata"))
#   roms.coords.proj <- cbind(c(Rdat$x), c(Rdat$y))
#   x.range <- c(min(roms.coords.proj[,1])-1000,max(roms.coords.proj[,1])+1000)
#   y.range <- c(min(roms.coords.proj[,2])-1000,max(roms.coords.proj[,2])+1000)
#   ra.region <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000)
#   ## crop npp to boundaries
#   npp2_crop <- crop(npp2, ext(ra.region))
#   ## convert to image and create random points
#   Z <- terra2im(npp2_crop/10000000) ## library(ppmData)
#   pp <- rpoispp(Z) ## library(spatstat)
#   npp.pts <- cbind(pp$x, pp$y,0)
#   save(npp.pts, file=paste0(env.dir3,"ocean_his_TrackingSetup_NPP_Region",sprintf("%02d",i),".Rdata"))
# }

####


















roms.coords.proj <- cbind(c(Rdat$x), c(Rdat$y))
x.range <- c(min(roms.coords.proj[,1])-1000,max(roms.coords.proj[,1])+1000)
y.range <- c(min(roms.coords.proj[,2])-1000,max(roms.coords.proj[,2])+1000)
empty.roms.ra <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000)
empty.roms.ra31 <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000, nlyr=31)
empty.roms.ra8 <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000, nlyr=8)
h.ra      <- setValues(empty.roms.ra, Rdat$h)
u.ra      <- setValues(empty.roms.ra, Rdat$i_u[,,1,1])
uv.ra      <- setValues(empty.roms.ra, sqrt(Rdat$i_u[,,1,1]^2 +c(Rdat$i_v[,,1,1])^2))

u31.ra      <- setValues(empty.roms.ra31, Rdat$i_u[,,,1])

#### Particle seeding
# npp <- rast("/pvol/data_environmental/Circumpolar_EnvData_2km_NPP_Cafe_filled_SummerAverage.tif")
# for(i in 1:10){
#   print(i)
#   row.lim <- lim[1,i]:lim[2,i]
#   col.lim <- lim[3,i]:lim[4,i]
#   
# ## npp
# roms.coords.proj <- cbind(c(Rdat$x), c(Rdat$y))
# x.range <- c(min(roms.coords.proj[,1])-1000,max(roms.coords.proj[,1])+1000)
# y.range <- c(min(roms.coords.proj[,2])-1000,max(roms.coords.proj[,2])+1000)
# empty.roms.ra <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000)
# npp2 <- project(npp, empty.roms.ra)
# npp2[npp2[]>0] <- 1
# #  npp2[na.sel[1,i]:na.sel[2,i],na.sel[3,i]:na.sel[4,i]] <- NA
# ## thin out to one particle every 4 cells
# npp2[seq(1,nrow(npp2), by=2),] <- NA
# npp2[,seq(1,ncol(npp2), by=2)] <- NA
# ##
# cell.pts <- cbind(crds(as.points(npp2)),0)
# save(cell.pts, file=paste0(env.dir3,"NPP_seeded_particles_Region",sprintf("%02d",i),".Rdata"))
# }














