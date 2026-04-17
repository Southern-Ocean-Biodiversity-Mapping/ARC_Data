#### 1) Set up ----
# load libraries
library(ncdf4)        ## package for netcdf manipulation
library(terra)
library(reproj) ## reproject coordinates
library(dplyr)
library(reshape) ## prep data for plotting
library(ggplot2) ## plotting
library(abind) 

## set up directory pointers etc
env.dir <- "/pvol/data_environmental/ROMS_2k_files/"

string.chr <- "Circumpolar_EnvData_"
string.res <- "500m_"

## polar stereographic projection:
stereo <- "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

################################
#### ROMS Currents
### WE NEED A RUN ACROSS SUMMER WITH HIGH-RES HISTORY FILES!!! 8 has one month only, 8_Nov has 10-day intervals...

### load lon/lat information from ROMS-grid
grd2k_nc <- nc_open(paste0(env.dir,"waom2extend_grd.nc"))
lon_rho.raw <- ncvar_get(grd2k_nc, varid="lon_rho") ## 3150 2800 - longitude
lat_rho.raw <- ncvar_get(grd2k_nc, varid="lat_rho") ## 3150 2800 - latitude
nc_close(grd2k_nc)

### create rasters
roms.coords.proj <- reproj(cbind(c(lon_rho.raw), c(lat_rho.raw)), target=stereo)
roms.coords.proj.lon <- matrix(roms.coords.proj[,1], ncol=2800)
roms.coords.proj.lat <- matrix(roms.coords.proj[,2], ncol=2800)
x.range <- c(min(roms.coords.proj[,1])-1000,max(roms.coords.proj[,1])+1000)
y.range <- c(min(roms.coords.proj[,2])-1000,max(roms.coords.proj[,2])+1000)
empty.roms.ra <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000)
empty.roms.ra31 <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000, nlyr=31)

################################
#### TO DO ONLY ONCE:

### shifted u and v current velocities to match rho-points, then extract interpolated current velocities at those rho-points
### for w, we simply removee the surface layer (32)
## 1. load current velocities
## 2. create a filled raster each for the original locations of the u, v, w points
## 3. extract interpolated u and v values at the rho points

### load u data from ROMS file
nc_2k <- nc_open(paste0(env.dir,"ocean_his_0001.nc"))
u.raw <- ncvar_get(nc_2k, varid="u")     ## 3149 2800 31 7 - u speed
nc_close(nc_2k)

# #plot(rast(u.raw[,,1,1]))
# ## a small test if the code is doing the right thing
# empty.roms.ra2 <- rast(extent=ext(c(x.range,y.range)), crs=stereo, resolution=2000, nlyr=2)
# x.u.range <- c(min(roms.coords.proj[,1]),max(roms.coords.proj[,1]))
# empty.u.ra <- rast(extent=ext(c(x.u.range,y.range)), crs=stereo, resolution=2000, nlyr=2)
# u.raw.ra.list <- list() ## to fill with rasters with original u positions
# u.ra.list <- list()     ## to fill with rasters with one more column and interpolated u positions
# ra.crds <- crds(empty.roms.ra)
# u.vals.array <- array(numeric(),c(3150,2800,2,3))
# ## this
# start.time <- Sys.time()
# for(k in 1:3){
#   print(k)
#   u.raw.ra.list[[k]] <- empty.u.ra
#   for(i in 1:2){
#     u.raw.ra.list[[k]][[i]][] <- c(u.raw[,,i,k])
#   }
#   u.vals.loop <- extract(u.raw.ra.list[[k]], ra.crds, method="bilinear")
#   u.ra.list[[k]] <- setValues(empty.roms.ra2, u.vals.loop)
#   u.vals.array[,,,k] <- u.ra.list[[k]][]
# }
# Sys.time()-start.time
# ## vs
# start.time <- Sys.time()
# for(k in 1:3){
#   print(k)
#   u.raw.ra.loop <- empty.u.ra
#   for(i in 1:2){
#     u.raw.ra.loop[[i]][] <- c(u.raw[,,i,k])
#   }
#   u.vals.loop <- extract(u.raw.ra.loop, ra.crds, method="bilinear")
#   u.ra.loop <- setValues(empty.roms.ra2, u.vals.loop)
#   u.vals.array[,,,k] <- u.ra.loop[]
# }
# Sys.time()-start.time

## for u-velocity, the x-dimension is reduced by one, so we simply don't extend the raster so that the boundary stops exactly where the rho-point is. But not for the y-dimension
## setup u raster
x.u.range <- c(min(roms.coords.proj[,1]),max(roms.coords.proj[,1]))
empty.u.ra <- rast(extent=ext(c(x.u.range,y.range)), crs=stereo, resolution=2000, nlyr=31)
#u.raw.ra.list <- list() ## to fill with rasters with original u positions
#u.ra.list <- list()     ## to fill with rasters with one more column and interpolated u positions
ra.crds <- crds(empty.roms.ra)
u.vals.array <- array(numeric(),c(3150,2800,31,2))
for(k in 1:2){
  message(k)
  # u.raw.ra.list[[k]] <- empty.u.ra
  u.raw.ra.loop <- empty.u.ra
  for(i in 1:31){
    print(i)
    # u.raw.ra.list[[k]][[i]][] <- c(u.raw[,2800:1,i,k])
    u.raw.ra.loop[[i]][] <- c(u.raw[,,i,k])
  }
  # u.vals.loop <- extract(u.raw.ra.list[[k]], ra.crds, method="bilinear")
  u.vals.loop <- extract(u.raw.ra.loop, ra.crds, method="bilinear")
  ## assign these values to a raster
  # u.ra.list[[k]] <- setValues(empty.roms.ra31, u.vals.loop)
  u.ra.loop <- setValues(empty.roms.ra31, u.vals.loop)
  ## and save them in an array so we can use them in the particle tracking
  # u.vals.array[,,,k] <- u.ra.list[[k]][]
  u.vals.array[,,,k] <- u.ra.loop[]
}
## save output
save(u.vals.array, file=paste0(env.dir,"ocean_his_0001_array_2slices_u.Rdata"))
rm(u.raw, u.vals.array, u.ra.loop, u.raw.ra.loop)

### same procedure for v data:
## load v data from ROMS file
nc_2k <- nc_open(paste0(env.dir,"ocean_his_0001.nc"))
v.raw <- ncvar_get(nc_2k, varid="v")     ## 3150 2799 31 7 - v speed
nc_close(nc_2k)
## for v-velocity, the y-dimension is reduced by one
## setup u raster
y.v.range <- c(min(roms.coords.proj[,2]),max(roms.coords.proj[,2]))
empty.v.ra <- rast(extent=ext(c(x.range,y.v.range)), crs=stereo, resolution=2000, nlyr=31)
v.vals.array <- array(numeric(),c(3150,2800,31,2))
for(k in 1:2){
  message(k)
   v.raw.ra.loop <- empty.v.ra
  for(i in 1:31){
    print(i)
     v.raw.ra.loop[[i]][] <- c(v.raw[,,i,k])
  }
  v.vals.loop <- extract(v.raw.ra.loop, ra.crds, method="bilinear")
  ## assign these values to a raster
  v.ra.loop <- setValues(empty.roms.ra31, v.vals.loop)
  ## and save them in an array so we can use them in the particle tracking
  v.vals.array[,,,k] <- v.ra.loop[]
}
## save output
save(v.vals.array, file=paste0(env.dir,"ocean_his_0001_array_2slices_v.Rdata"))
rm(v.raw, v.vals.array, v.ra.loop, v.raw.ra.loop)


################################
### load ROMS data including w velocities
nc_2k <- nc_open(paste0(env.dir,"ocean_his_0001.nc"))
s_rho <- ncvar_get(nc_2k, varid="s_rho") ## 31             - vertical levels varying between 0 to 1
Cs_r  <- ncvar_get(nc_2k, varid="Cs_r")  ## 31             - "S-coordinate stretching curves at RHO-points
hc    <- ncvar_get(nc_2k, varid="hc")    ##                - value of 250
h     <- ncvar_get(nc_2k, varid="h")     ## 3150 2800      - model bathymetry
zeta  <- ncvar_get(nc_2k, varid="zeta")  ## 3150 2800 7    - free surface elevation
zice  <- ncvar_get(nc_2k, varid="zice")  ## 3150 2800      - ice draft
w.raw <- ncvar_get(nc_2k, varid="w")     ## 3150 2800 32 7 - w speed
nc_close(nc_2k)

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

### setup ROMS object
ROMS_2k <- list()
ROMS_2k$x <- roms.coords.proj.lon#[,2800:1]
ROMS_2k$y <- roms.coords.proj.lat#[,2800:1]
ROMS_2k$h <- h#[,2800:1]
ROMS_2k$hh <- hh#[,,31:1] ## for some reason this was upside down.. so need to invert
ROMS_2k$zice <- zice
load(paste0(env.dir,"ocean_his_0001_array_2slices_u.Rdata"))
ROMS_2k$i_u <- u.vals.array
rm(u.vals.array)
load(paste0(env.dir,"ocean_his_0001_array_2slices_v.Rdata"))
ROMS_2k$i_v <- v.vals.array
rm(v.vals.array)
ROMS_2k$i_w <- w.raw[,,-32,1:2]
rm(w.raw)

## check if it all makes sense
u.ra      <- setValues(empty.roms.ra, t(ROMS_2k$i_u[,,1,1]))
u.surf.ra <- setValues(empty.roms.ra, ROMS_2k$i_u[,,31,1])
h.ra       <- setValues(empty.roms.ra, c(ROMS_2k$h))

## need to transpose the array so that the continent appears in the most commonly seen orientation
ROMS_2k_t <- ROMS_2k
ROMS_2k_t$x <- t(ROMS_2k$x[,2800:1])
ROMS_2k_t$y <- t(ROMS_2k$y[,2800:1])
ROMS_2k_t$h <- t(ROMS_2k$h[,2800:1])
ROMS_2k_t$hh <- aperm(ROMS_2k$hh[,2800:1,], c(2,1,3))
ROMS_2k_t$zice <- t(ROMS_2k$zice[,2800:1])
ROMS_2k_t$i_u <- aperm(ROMS_2k$i_u[,2800:1,,], c(2,1,3,4))
ROMS_2k_t$i_v <- aperm(ROMS_2k$i_v[,2800:1,,], c(2,1,3,4))
ROMS_2k_t$i_w <- aperm(ROMS_2k$i_w[,2800:1,,], c(2,1,3,4))

u.ra      <- setValues(empty.roms.ra, ROMS_2k_t$i_u[,,1,1])
u.surf.ra <- setValues(empty.roms.ra, ROMS_2k$i_u[,,31,1])
h.ra       <- setValues(empty.roms.ra, c(ROMS_2k$h))

save(ROMS_2k_t, file=paste0(env.dir,"ocean_his_0001_2slices_TrackingSetup_transposed.Rdata"))
####################################################
