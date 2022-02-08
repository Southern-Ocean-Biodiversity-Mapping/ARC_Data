
#title: "Circumpolar Environmental Data"
#author: "Jan Jansen"
#date: "15/07/2021"

### This code extracts circumpolar bathymetry data and generate the topographic position index at varying scales
### as well as generates Net Primary Productivity climatologies from previous files extracted for FAM modelling
### and seafloor variables from the circumpolar Regional Oceanographic Modelling System (ROMS) developed by Ben Galton-Fenzi and Fabio ***


## Things you need to do to run this script

# 1. specify the directory where to find the data

# - specify your local directory or point to the owncloud directories (THIS DOESN'T WORK YET)

# - Download the data from the link below and change "env.dir" in the code snippet below to match your local machine.

#Owncloud repository folder is "EnvironmentalData". All data can be found here:
#https://owncloud.imas-data-service.cloud.edu.au/index.php/apps/files/?dir=/EnvironmentalData&fileid=186349

#### 1) Set up ----
# load libraries
library(ncdf4)        ## package for netcdf manipulation
library(raadtools)
library(raster)       ## package for raster manipulation
library(sp)
library(dplyr)
library(blueant)
library(rgdal)        ## package for geospatial analysis
library(ggplot2)      ## package for plotting

library(spatialEco)

# set up directory pointers etc
# Jan's local machine:
env.dir <- "C:/Users/jjansen/Desktop/science/data_environmental/"

# remote repository (DOESN'T WORK YET)
# env.dir <- "https://data.imas.utas.edu.au/data_transfer/admin/files/EnvironmentalData/"

env.raw <- paste0(env.dir,"raw/")
env.derived <- paste0(env.dir,"derived/")
AAD_dir <- paste0(env.dir,"raw/accessed_through_R")

string.chr <- "Circumpolar_EnvData_"
string.res <- "500m_"

# polar stereographic projection:
stereo <- "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

                                                                

#### 2) Source AAD datafiles----
#Code to download data from other remote repositories to your local machine
#https://ropensci.org/blog/2018/11/13/antarctic/

#This is a very large file
## seaice data needs login...see separate code file for sea-ice
src <- bind_rows(sources("NSIDC SMMR-SSM/I Nasateam sea ice concentration", hemisphere = "south", time_resolutions = "day", years = c(2013)),
                 sources("Southern Ocean summer chlorophyll-a climatology (Johnson)"),
                 sources("IBCSO bathymetry"),
                 bb_modify_source(sources("Oceandata MODIS Aqua Level-3 mapped daily 4km chl-a"), method = list(search = "A2014*L3m_DAY_CHL_chlor_a_4km.nc")))
result <- bb_get(src, local_file_root = AAD_dir, clobber = 0, verbose = TRUE, confirm = NULL)                                                                 


#But usually we only need to set the directory where the AAD functions look for the data
set_data_roots(AAD_dir)


#### 3) Extract bathymetry and generate TPI----
# choose which data-layer to use, which resolution and projection.

## Choose dataset
#data.name <- "ibcso"  ## ibcso can be found using Ben's AAd dataset-wrapper
#data.name <- "ibcso2"  ## DOESN'T EXIST YET, DOESN'T FUNCTION IN THE CODE YET
data.name <- "gebco"  ## gebco is on file

## Choose resolution and projection (only relevant for gebco bathymetry)
data.format <- "project.500" ## project polar stereographic to 500m resolution
#data.format <- "project.raw" ## project polar stereographic but keep resolution as it is
#data.format <- "unprojected" ## keep raw format

## Choose range of depth
data.depth <- "shelf"
#data.depth <- "full"

if(data.depth=="shelf"){
  depth.range <- c(0,-3000)
}else depth.range <- c(0,-Inf)


#NOTE: gebco, takes 30-40mins to run!

## read in bathymetry data and project/reproject if needed
r <- readtopo("ibcso") ## already projected and at 500m resolution
if(data.name=="gebco"){
  ## load raw (unprojected) data
  # g <- raster(paste0(env.raw,"GEBCO_2020/gebco_2020_SO.tif"))
  # g <- raster(paste0(env.raw,"gebco_2021_sub_ice_topo/GEBCO_2021_sub_ice_topo_CF.nc"))
  # crs(g) <- "+proj=longlat +datum=WGS84 +no_defs"
  # g <- crop(g,extent(c(-180,180,-90,-60)))
  g <- raster(paste0(env.raw,"GEBCO_2021/gebco_2021_SO.tif"))
  if(data.format=="project.500"|data.format=="project.raw"){
    ## project raster to polar stereographic
    g <- projectRaster(g,crs=stereo)
  }
  if(data.format=="project.500"){
    ## change resolution to 500m grid cells
    g <- projectRaster(g,r)
  }
  r <- g  
}

## create layers for depth, slope and topographic position index (TPI) at different scales
r.depth <- r
r.slope <- terrain(r)
r.tpi <- tpi(r)
r.tpi5 <- tpi(r, scale=5)
r.tpi11 <- tpi(r, scale=11)
# r.tpi21 <- tpi(r, scale=21)
# r.tpi31 <- tpi(r, scale=31)

## set anything on land to NA, and optionally set abyssal zone to NA
r.depth[r>=0] <- NA
r.depth[r<=depth.range[2]] <- NA
r.slope[is.na(r.depth[])] <- NA
r.tpi[is.na(r.depth[])] <- NA
r.tpi5[is.na(r.depth[])] <- NA
r.tpi11[is.na(r.depth[])] <- NA
# r.tpi21[is.na(r.depth[])] <- NA
# r.tpi31[is.na(r.depth[])] <- NA

r.bathy <- stack(r.depth, r.slope, r.tpi, r.tpi5, r.tpi11)#, r.tpi21, r.tpi31)
names(r.bathy) <- c("depth","slope","tpi","tpi5","tpi11")#"tpi21","tpi31")

# # xlim=c(128,148)
# # ylim=c(-67.5,-64)
# par(mfrow=c(2,2))
# plot(r.depth) #, xlim=xlim, ylim=ylim
# plot(r.slope) #, xlim=xlim, ylim=ylim
# plot(r.tpi5)  #, xlim=xlim, ylim=ylim
# plot(r.tpi31) #, xlim=xlim, ylim=ylim

# Save data

if(data.name=="ibcso"){
  res.string <- string.res
  ra.string <- ""
  if(data.depth=="shelf"){
      res.string <- string.res
      ra.string <- "shelf_"
  }
}
if(data.name=="gebco"){
  if(data.depth=="full"){
    if(data.format=="project.raw"){ ## original resolution (119m x 460m)
      res.string <- ""
      ra.string <- ""
    }
    if(data.format=="project.500"){ ## 500m resolution
      res.string <- string.res
      ra.string <- ""      
    }
  }
  if(data.depth=="shelf"){
    if(data.format=="project.raw"){ ## original resolution (119m x 460m)
      res.string <- ""
      ra.string <- "shelf_"
    }
    if(data.format=="project.500"){ ## 500m resolution
      res.string <- string.res
      ra.string <- "shelf_"      
    }
  }
}
save.string <- paste0(env.derived, string.chr, res.string, ra.string, "bathy_", data.name)

# writeRaster(g, filename=paste0("C:/Users/jjansen/Desktop/science/data_environmental/derived/Circumpolar_EnvData_500m_bathy_gebco.Rdata"), overwrite=TRUE)
writeRaster(r.depth, filename=paste0(save.string,"_depth.grd"), overwrite=TRUE)
writeRaster(r.slope, filename=paste0(save.string,"_slope.grd"), overwrite=TRUE)
writeRaster(r.tpi,   filename=paste0(save.string,"_tpi.grd"), overwrite=TRUE)
writeRaster(r.tpi5,  filename=paste0(save.string,"_tpi5.grd"), overwrite=TRUE)
writeRaster(r.tpi11, filename=paste0(save.string,"_tpi11.grd"), overwrite=TRUE)


#### 4) Net Primary Production ----

#Raw NPP is being read in and manipulated using a different script as part of preparing the input for the circumpolar ocean model ("ReadIn_Circumpolar_Environmental_Data_ROMS_NPP.Rmd").    
#Here, we're using the standardised files found in the derived data folder to produce NPP climatologies.  
#Files that are "filled": missing npp observations where chla observation were present, were replaced with predicted values from the npp-chla relationship in the region.  

# Read data from NetCDF files

years <- as.character(2002:2020)

## npp list
npp.list <- list()

## load data:
npp.list$npp_2002 <- get(load(paste0(env.derived,string.chr,"NPP_Cafe_filled_2002.Rdata")))
npp.list$npp_2003 <- get(load(paste0(env.derived,string.chr,"NPP_Cafe_filled_2003.Rdata")))
npp.list$npp_2004 <- get(load(paste0(env.derived,string.chr,"NPP_Cafe_filled_2004.Rdata")))
npp.list$npp_2005 <- get(load(paste0(env.derived,string.chr,"NPP_Cafe_filled_2005.Rdata")))
npp.list$npp_2006 <- get(load(paste0(env.derived,string.chr,"NPP_Cafe_filled_2006.Rdata")))
npp.list$npp_2007 <- get(load(paste0(env.derived,string.chr,"NPP_Cafe_filled_2007.Rdata")))
npp.list$npp_2008 <- get(load(paste0(env.derived,string.chr,"NPP_Cafe_filled_2008.Rdata")))
npp.list$npp_2009 <- get(load(paste0(env.derived,string.chr,"NPP_Cafe_filled_2009.Rdata")))
npp.list$npp_2010 <- get(load(paste0(env.derived,string.chr,"NPP_Cafe_filled_2010.Rdata")))
npp.list$npp_2011 <- get(load(paste0(env.derived,string.chr,"NPP_Cafe_filled_2011.Rdata")))
npp.list$npp_2012 <- get(load(paste0(env.derived,string.chr,"NPP_Cafe_filled_2012.Rdata")))
npp.list$npp_2013 <- get(load(paste0(env.derived,string.chr,"NPP_Cafe_filled_2013.Rdata")))
npp.list$npp_2014 <- get(load(paste0(env.derived,string.chr,"NPP_Cafe_filled_2014.Rdata")))
npp.list$npp_2015 <- get(load(paste0(env.derived,string.chr,"NPP_Cafe_filled_2015.Rdata")))
npp.list$npp_2016 <- get(load(paste0(env.derived,string.chr,"NPP_Cafe_filled_2016.Rdata")))
npp.list$npp_2017 <- get(load(paste0(env.derived,string.chr,"NPP_Cafe_filled_2017.Rdata")))
npp.list$npp_2018 <- get(load(paste0(env.derived,string.chr,"NPP_Cafe_filled_2018.Rdata")))
npp.list$npp_2019 <- get(load(paste0(env.derived,string.chr,"NPP_Cafe_filled_2019.Rdata")))
npp.list$npp_2020 <- get(load(paste0(env.derived,string.chr,"NPP_Cafe_filled_2020.Rdata")))

## stack each month, keeping in mind:
## the year 2002 only has data for the second half of the year
## 2020 is missing September

list_01 <- list()
list_02 <- list()
list_03 <- list()
list_10 <- list()
list_11 <- list()
list_12 <- list()

for(i in 1:length(npp.list)){
  if(years[i]=="2002"){
    list_01[[i]] <- NA
    list_02[[i]] <- NA
    list_03[[i]] <- NA
    list_10[[i]] <- npp.list[[i]][[2]]
    list_11[[i]] <- npp.list[[i]][[3]]
    list_12[[i]] <- npp.list[[i]][[4]]
  }else if(years[i]=="2020"){
    list_01[[i]] <- npp.list[[i]][[1]]
    list_02[[i]] <- npp.list[[i]][[2]]
    list_03[[i]] <- npp.list[[i]][[3]]
    list_10[[i]] <- npp.list[[i]][[5]]
    list_11[[i]] <- npp.list[[i]][[6]]
    list_12[[i]] <- npp.list[[i]][[7]]
  }else{
    list_01[[i]] <- npp.list[[i]][[1]]
    list_02[[i]] <- npp.list[[i]][[2]]
    list_03[[i]] <- npp.list[[i]][[3]]
    list_10[[i]] <- npp.list[[i]][[6]]
    list_11[[i]] <- npp.list[[i]][[7]]
    list_12[[i]] <- npp.list[[i]][[8]]
  }
  names(list_01)[i] <- names(list_02)[i] <- names(list_03)[i] <- names(list_10)[i] <- names(list_11)[i] <- names(list_12)[i] <- names(npp.list)[i]
}

Jan_ave <- mean(brick(list_01[-1]), na.rm=TRUE)
Feb_ave <- mean(brick(list_02[-1]), na.rm=TRUE)
Mar_ave <- mean(brick(list_03[-1]), na.rm=TRUE)
Oct_ave <- mean(brick(list_10), na.rm=TRUE)
Nov_ave <- mean(brick(list_11), na.rm=TRUE)
Dec_ave <- mean(brick(list_12), na.rm=TRUE)

Su_ave <- mean(Jan_ave,Feb_ave,Mar_ave,Oct_ave,Nov_ave,Dec_ave, na.rm=TRUE)
writeRaster(Su_ave, filename=paste0(env.derived,string.chr,"NPP_Cafe_filled_SummerAverage.Rdata"), overwrite=TRUE)

npp_su <- projectRaster(Su_ave,r)
writeRaster(npp_su, filename=paste0(env.derived,string.chr,"500m_NPP_Cafe_filled_SummerAverage.Rdata"), overwrite=TRUE)

npp_su_shelf <- npp_su
npp_su_shelf[is.na(r.depth)] <- NA
writeRaster(npp_su_shelf, filename=paste0(env.derived,string.chr,"500m_shelf_NPP_Cafe_filled_SummerAverage.Rdata"), overwrite=TRUE)


#### 5) ROMS Currents & Temperature & FAM ----

## 4k models for now, 2k to follow, and proper FAM to follow!!
data.dat100 <- paste0(env.dir,"Circumpolar_ROMS/4km_outputs/output_sed_test1/")
#### load lon/lat information from ROMS-grid
grd4k_nc <- nc_open(paste0(env.raw,"waom4extend_grd.nc"))
lon_rho <- ncvar_get(grd4k_nc, varid="lon_rho")
lat_rho <- ncvar_get(grd4k_nc, varid="lat_rho")
#### Prepare empty rasters to assign correct projected values to
roms.coords.proj <- project(cbind(c(lon_rho), c(lat_rho)), proj=stereo)
x.range <- c(min(roms.coords.proj[,1])-2000,max(roms.coords.proj[,1])+2000)
y.range <- c(min(roms.coords.proj[,2])-2000,max(roms.coords.proj[,2])+2000)
empty.roms.ra <- raster(extent(c(x.range,y.range)), crs=stereo, resolution=4000)

#depth
h <- raster(paste0(data.dat100,"ocean_avg_0001.nc"), varname="h", level=1)
#salinity
salt <- raster(paste0(data.dat100,"ocean_avg_0001.nc"), varname="salt", level=1)
#temperature
temp <- raster(paste0(data.dat100,"ocean_avg_0001.nc"), varname="temp", level=1)
#settling FAM
settle6 <- raster(paste0(data.dat100,"ocean_avg_0001.nc"), varname="sand_06", level=1)
#seafloor currents (seafloor-layer is 1)
u.raw <- brick(paste0(data.dat100,"ocean_his_0001.nc"), varname="u", level=1)
v.raw <- brick(paste0(data.dat100,"ocean_his_0001.nc"), varname="v", level=1)
#seasurface currents (surface-layer is 31)
u_31.raw <- brick(paste0(data.dat100,"ocean_avg_0001.nc"), varname="u", level=31)
v_31.raw <- brick(paste0(data.dat100,"ocean_avg_0001.nc"), varname="v", level=31)
## sum up monthly values for a climatology (THIS SHOULD BE DONE ON HIGH RESOLUTION HISTORY FILES)
u.sum <- sum(u.raw)
v.sum <- sum(v.raw)
u.sum.abs <- sum(abs(u.raw))
v.sum.abs <- sum(abs(v.raw))
u31.sum <- sum(u_31.raw)
v31.sum <- sum(v_31.raw)
## extract current speeds at rho-points where depth is defined
## (in ROMS they are all at different locations):
## u has one less column than the grid, so all coordinates need to be moved half a cell to the right for the grids to match up
coord.grd.u <- coordinates(h)
coord.grd.u[,1] <- coordinates(h)[,1]-0.5
## v has one less row than the grid, so all coordinates need to be moved half a cell up for the grids to match up
coord.grd.v <- coordinates(h)
coord.grd.v[,2] <- coordinates(h)[,2]-0.5
## now extract values at the rho-points and interpolate (because they are 2km away from the nearest original point), and place into projected raster
u <- v <- u.abs <- v.abs <- u31 <- v31 <- empty.roms.ra
u[] <- extract(u.sum, coord.grd.u, method="bilinear")
v[] <- extract(v.sum, coord.grd.v, method="bilinear")
u.abs[] <- extract(u.sum.abs, coord.grd.u, method="bilinear")
v.abs[] <- extract(v.sum.abs, coord.grd.v, method="bilinear")
u31[] <- extract(u31.sum, coord.grd.u, method="bilinear")
v31[] <- extract(v31.sum, coord.grd.v, method="bilinear")
#seasurface current speeds
uv_31 <- sqrt(u31^2+v31^2)
#temporal mean seafloor current speed
mean.uv <- sqrt(u^2+v^2)
#absolute mean seafloor current speed
abs.uv <- sqrt(u.abs^2+v.abs^2)
#residual seafloor current speed
res.uv <- abs.uv-mean.uv
## remove inland values for depth
h2 <- empty.roms.ra
h2[] <- h[]
h2[is.na(mean.uv)] <- NA

## resample to standard 500m resolution of other environmental variables
mean.uv_500 <- resample(mean.uv,r)
abs.uv_500 <- resample(abs.uv,r)
res.uv_500 <- resample(res.uv,r)
t_500 <- resample(temp,r)
s_500 <- resample(salt,r)
settle6_500 <- resample(settle6,r)

## shelf only
mean.uv_500_shelf <- mean.uv_500
abs.uv_500_shelf <- abs.uv_500
res.uv_500_shelf <- res.uv_500
t_500_shelf <- t_500
s_500_shelf <- s_500
settle6_500_shelf <- settle6_500

mean.uv_500_shelf[is.na(r.depth)] <- NA
abs.uv_500_shelf[is.na(r.depth)] <- NA
res.uv_500_shelf[is.na(r.depth)] <- NA
t_500_shelf[is.na(r.depth)] <- NA
s_500_shelf[is.na(r.depth)] <- NA
settle6_500_shelf[is.na(r.depth)] <- NA

## write rasters to file
writeRaster(mean.uv_500,       filename=paste0(env.derived,string.chr,"500m_waom4k_seafloorcurrents_mean.Rdata"))
writeRaster(mean.uv_500_shelf, filename=paste0(env.derived,string.chr,"500m_shelf_waom4k_seafloorcurrents_mean.Rdata"))
writeRaster(abs.uv_500,       filename=paste0(env.derived,string.chr,"500m_waom4k_seafloorcurrents_absolute.Rdata"))
writeRaster(abs.uv_500_shelf, filename=paste0(env.derived,string.chr,"500m_shelf_waom4k_seafloorcurrents_absolute.Rdata"))
writeRaster(res.uv_500,       filename=paste0(env.derived,string.chr,"500m_waom4k_seafloorcurrents_residual.Rdata"))
writeRaster(res.uv_500_shelf, filename=paste0(env.derived,string.chr,"500m_shelf_waom4k_seafloorcurrents_residual.Rdata"))
writeRaster(temp,            filename=paste0(env.derived,string.chr,"waom4k_seafloortemperature.Rdata"))
writeRaster(t_500,        filename=paste0(env.derived,string.chr,"500m_waom4k_seafloortemperature.Rdata"))
writeRaster(t_500_shelf,  filename=paste0(env.derived,string.chr,"500m_shelf_waom4k_seafloortemperature.Rdata"))
writeRaster(salt,            filename=paste0(env.derived,string.chr,"waom4k_seafloorsalinity.Rdata"))
writeRaster(s_500,        filename=paste0(env.derived,string.chr,"500m_waom4k_seafloorsalinity.Rdata"))
writeRaster(s_500_shelf,  filename=paste0(env.derived,string.chr,"500m_shelf_waom4k_seafloorsalinity.Rdata"))
writeRaster(settle6,            filename=paste0(env.derived,string.chr,"waom4k_settle6test.Rdata"))
writeRaster(settle6_500,        filename=paste0(env.derived,string.chr,"500m_waom4k_settle6test.Rdata"))
writeRaster(settle6_500_shelf,  filename=paste0(env.derived,string.chr,"500m_shelf_waom4k_settle6test.Rdata"))


# #NOTE: 2k res (UPDATE ONCE FAM HAS RUN)  
# #NetCDF-files are downloaded from GADI
# 
# ## file paths
# f.grd <- paste0(env.raw,"waom2_grd.nc")
# f.u <- paste0(env.raw,"ocean_avg_0538-0610_u_avg.nc")
# f.v <- paste0(env.raw,"ocean_avg_0538-0610_v_avg.nc")
# f.t <- paste0(env.raw,"ocean_avg_0538-0610_temp_avg.nc")
# 
# ## read as raster
# lon <- raster(f.grd, varname="lon_rho")
# lat <- raster(f.grd, varname="lat_rho")
# u <- raster(f.u, lvar=3, level=1, varname="u")
# v <- raster(f.v, lvar=3, level=1, varname="v")
# t <- raster(f.t, lvar=3, level=1, varname="temp")
# 
# ## bring to same extent (address one of the quirks of ROMS)
# ext <- extent(1,3149,1,2649)
# lon <- crop(lon,y=ext,snap="out")
# lat <- crop(lat,y=ext,snap="out")
# u <- crop(u,y=ext,snap="out")
# v <- crop(v,y=ext,snap="out")
# t <- crop(t,y=ext,snap="out")
# #w <- crop(w,y=ext,snap="out")
# 
# ## calculate a single seafloor current speed value
# uv <- sqrt(abs(u)^2+abs(v)^2)
# ## projection and extent for the raster (netcdf files were already polar-projected with true south at -71S)
# crs <- stereo ##"+proj=stere +lat_ts=-71 +lat_0=-90 +datum=WGS84"
# pts <- rgdal::project(cbind(values(lon), values(lat)), crs)
# ex <- extent(pts)
# uv <- setExtent(uv, ex)
# t <- setExtent(t, ex)
# projection(uv) <- crs
# projection(t) <- crs
# 
# ## resample to standard 500m resolution of other environmental variables
# uv_500 <- resample(uv,r)
# t_500 <- resample(t,r)
# 
# ## shelf only
# uv_500_shelf <- uv_500
# uv_500_shelf[is.na(r.depth)] <- NA
# t_500_shelf <- t_500
# t_500_shelf[is.na(r.depth)] <- NA
# 
# ## write rasters to file
# writeRaster(uv,           filename=paste0(env.derived,string.chr,"waom2k_seafloorcurrents.Rdata"))
# writeRaster(uv_500,       filename=paste0(env.derived,string.chr,"500m_waom2k_seafloorcurrents.Rdata"))
# writeRaster(uv_500_shelf, filename=paste0(env.derived,string.chr,"500m_shelf_waom2k_seafloorcurrents.Rdata"))
# writeRaster(t,            filename=paste0(env.derived,string.chr,"waom2k_seafloortemperature.Rdata"))
# writeRaster(t_500,        filename=paste0(env.derived,string.chr,"500m_waom2k_seafloortemperature.Rdata"))
# writeRaster(t_500_shelf,  filename=paste0(env.derived,string.chr,"500m_shelf_waom2k_seafloortemperature.Rdata"))

