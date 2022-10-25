
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
#library(terra)

library(spatialEco)

# set up directory pointers etc
# Jan's local machine:
env.dir <- "C:/Users/jjansen/Desktop/science/data_environmental/"

# remote repository (DOESN'T WORK YET)
# env.dir <- "https://data.imas.utas.edu.au/data_transfer/admin/files/EnvironmentalData/"

env.raw <- "E:/science/data_environmental/raw/"
env.derived <- paste0(env.dir,"derived/")
AAD_dir <- "E:/science/data_environmental/raw/accessed_through_R"

string.chr <- "Circumpolar_EnvData_"
string.res <- "500m_"

# polar stereographic projection:
stereo <- "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

                                                                

#### 2) Source AAD datafiles----
#Code to download data from other remote repositories to your local machine
#https://ropensci.org/blog/2018/11/13/antarctic/

#This is a very large file
## seaice data needs login...see separate code file for sea-ice
# src <- bind_rows(sources("NSIDC SMMR-SSM/I Nasateam sea ice concentration", hemisphere = "south", time_resolutions = "day", years = c(2013)),
#                  sources("Southern Ocean summer chlorophyll-a climatology (Johnson)"),
#                  sources("IBCSO bathymetry"),
#                  bb_modify_source(sources("Oceandata MODIS Aqua Level-3 mapped daily 4km chl-a"), method = list(search = "A2014*L3m_DAY_CHL_chlor_a_4km.nc")))
src <- bind_rows(sources("Southern Ocean summer chlorophyll-a climatology (Johnson)"),
                 sources("IBCSO bathymetry"),
                 )
result <- bb_get(src, local_file_root = AAD_dir, clobber = 0, verbose = TRUE, confirm = NULL)                                                                 


#But usually we only need to set the directory where the AAD functions look for the data
set_data_roots(AAD_dir)


#### 3) Extract bathymetry and generate TPI----
# choose which data-layer to use, which resolution and projection.

## Choose dataset
#data.name <- "ibcso"  ## ibcso can be found using Ben's AAd dataset-wrapper
data.name <- "ibcso2"  ## DOESN'T EXIST YET, DOESN'T FUNCTION IN THE CODE YET
#data.name <- "gebco"  ## gebco is on file

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
ri <- readtopo("ibcso") ## already projected and at 500m resolution
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

## OR, load the latest IBCSO version here:
r.raw <- raster(paste0(env.raw,"IBCSO_2022/IBCSO_v2_ice-surface.tif"))
r.raw.c <- crop(r.raw, ri)
r <- projectRaster(r.raw.c, ri)

## create layers for depth, slope and topographic position index (TPI) at different scales
r.depth <- r
#r.depth[r.depth>=200] <- NA
#r.depth[r.depth<=-4000] <- NA

r.slope <- terrain(r.depth)
r.tpi <- tpi(r.depth)
r.tpi5 <- tpi(r.depth, scale=5)
r.tpi11 <- tpi(r.depth, scale=11)
# r.tpi21 <- tpi(r, scale=21)
# r.tpi31 <- tpi(r, scale=31)

## set anything on land to NA, and optionally set abyssal zone to NA
r.depth[r>0] <- NA
r.depth[r<=depth.range[2]] <- NA
r.slope[is.na(r.depth[])] <- NA
r.tpi[is.na(r.depth[])] <- NA
r.tpi5[is.na(r.depth[])] <- NA
r.tpi11[is.na(r.depth[])] <- NA
# r.tpi21[is.na(r.depth[])] <- NA
# r.tpi31[is.na(r.depth[])] <- NA

# r.bathy <- stack(r.depth, r.slope, r.tpi, r.tpi5, r.tpi11)#, r.tpi21, r.tpi31)
# names(r.bathy) <- c("depth","slope","tpi","tpi5","tpi11")#"tpi21","tpi31")

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
# writeRaster(r.depth, filename=paste0(save.string,"_depth.grd"), overwrite=TRUE)
# writeRaster(r.slope, filename=paste0(save.string,"_slope.grd"), overwrite=TRUE)
# writeRaster(r.tpi,   filename=paste0(save.string,"_tpi.grd"), overwrite=TRUE)
# writeRaster(r.tpi5,  filename=paste0(save.string,"_tpi5.grd"), overwrite=TRUE)
# writeRaster(r.tpi11, filename=paste0(save.string,"_tpi11.grd"), overwrite=TRUE)
writeRaster(r.depth, filename=paste0(save.string,"_depth.tif"), overwrite=TRUE)
writeRaster(r.slope, filename=paste0(save.string,"_slope.tif"), overwrite=TRUE)
writeRaster(r.tpi,   filename=paste0(save.string,"_tpi.tif"), overwrite=TRUE)
writeRaster(r.tpi5,  filename=paste0(save.string,"_tpi5.tif"), overwrite=TRUE)
writeRaster(r.tpi11, filename=paste0(save.string,"_tpi11.tif"), overwrite=TRUE)

##### redo the same for 2km resolution data:
rm(r.depth, r.slope, r.tpi, r.tpi5, r.tpi11)
r2k.depth <- aggregate(r, 4)
# r2k.depth[r2k.depth>=200] <- NA
# r2k.depth[r2k.depth<=-3500] <- NA

r2k.slope <- terrain(r2k.depth)
r2k.tpi <- tpi(r2k.depth)
r2k.tpi5 <- tpi(r2k.depth, scale=5)
r2k.tpi11 <- tpi(r2k.depth, scale=11)
# r.tpi21 <- tpi(r, scale=21)
# r.tpi31 <- tpi(r, scale=31)

## set anything on land to NA, and optionally set abyssal zone to NA
r2k.depth[r2k.depth>=0] <- NA
r2k.depth[r2k.depth<=depth.range[2]] <- NA
r2k.slope[is.na(r2k.depth[])] <- NA
r2k.tpi[is.na(r2k.depth[])] <- NA
r2k.tpi5[is.na(r2k.depth[])] <- NA
r2k.tpi11[is.na(r2k.depth[])] <- NA

save.string <- paste0(env.derived, string.chr, "2km_", ra.string, "bathy_", data.name)
writeRaster(r2k.depth, filename=paste0(save.string,"_depth.tif"), overwrite=TRUE)
writeRaster(r2k.slope, filename=paste0(save.string,"_slope.tif"), overwrite=TRUE)
writeRaster(r2k.tpi,   filename=paste0(save.string,"_tpi.tif"), overwrite=TRUE)
writeRaster(r2k.tpi5,  filename=paste0(save.string,"_tpi5.tif"), overwrite=TRUE)
writeRaster(r2k.tpi11, filename=paste0(save.string,"_tpi11.tif"), overwrite=TRUE)


###

r <- raster(paste0(env.derived,string.chr,"500m_bathy_gebco.grd"))
r[r>=0] <- NA
r[r<=depth.range[2]] <- NA

r.depth <- raster(paste0(env.derived,string.chr,"500m_shelf_bathy_ibcso2_depth.tif"))

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
Jan_sd <- calc(brick(list_01[-1]), fun=sd, na.rm=TRUE)
Feb_sd <- calc(brick(list_02[-1]), fun=sd, na.rm=TRUE)
Mar_sd <- calc(brick(list_03[-1]), fun=sd, na.rm=TRUE)
Oct_sd <- calc(brick(list_10), fun=sd, na.rm=TRUE)
Nov_sd <- calc(brick(list_11), fun=sd, na.rm=TRUE)
Dec_sd <- calc(brick(list_12), fun=sd, na.rm=TRUE)

Su_ave <- calc(stack(Jan_ave,Feb_ave,Mar_ave,Oct_ave,Nov_ave,Dec_ave), fun=mean, na.rm=TRUE)
Su_sd <- calc(stack(brick(list_01[-1]),brick(list_02[-1]),brick(list_03[-1]),
                    brick(list_10),brick(list_11), brick(list_12)), fun=sd, na.rm=TRUE)
writeRaster(Su_ave, filename=paste0(env.derived,string.chr,"NPP_Cafe_filled_SummerAverage.tif"), overwrite=TRUE)
writeRaster(Su_sd, filename=paste0(env.derived,string.chr,"NPP_Cafe_filled_SummerStandardDeviation.tif"), overwrite=TRUE)

npp_su <- projectRaster(Su_ave,r)
npp_su_sd <- projectRaster(Su_sd,r)
writeRaster(npp_su, filename=paste0(env.derived,string.chr,"500m_NPP_Cafe_filled_SummerAverage.tif"), overwrite=TRUE)
writeRaster(npp_su_sd, filename=paste0(env.derived,string.chr,"500m_NPP_Cafe_filled_SummerStandardDeviation.tif"), overwrite=TRUE)

npp2k_su <- projectRaster(Su_ave,r2k.depth)
npp2k_su_sd <- projectRaster(Su_sd,r2k.depth)
writeRaster(npp2k_su, filename=paste0(env.derived,string.chr,"2km_NPP_Cafe_filled_SummerAverage.tif"), overwrite=TRUE)
writeRaster(npp2k_su_sd, filename=paste0(env.derived,string.chr,"2km_NPP_Cafe_filled_SummerStandardDeviation.tif"), overwrite=TRUE)

r.depth <- raster(paste0(env.derived,string.chr,"500m_shelf_bathy_ibcso2_depth.tif"))
npp_su_shelf <- npp_su
npp_su_shelf[is.na(r.depth)] <- NA
npp_su_sd_shelf <- npp_su_sd
npp_su_sd_shelf[is.na(r.depth)] <- NA
writeRaster(npp_su_shelf, filename=paste0(env.derived,string.chr,"500m_shelf_NPP_Cafe_filled_SummerAverage.tif"), overwrite=TRUE)
writeRaster(npp_su_sd_shelf, filename=paste0(env.derived,string.chr,"500m_shelf_NPP_Cafe_filled_SummerStandardDeviation.tif"), overwrite=TRUE)

npp2k_su_shelf <- npp2k_su
npp2k_su_shelf[is.na(r2k.depth)] <- NA
npp2k_su_sd_shelf <- npp2k_su_sd
npp2k_su_sd_shelf[is.na(r2k.depth)] <- NA
writeRaster(npp2k_su_shelf, filename=paste0(env.derived,string.chr,"2km_shelf_NPP_Cafe_filled_SummerAverage.tif"), overwrite=TRUE)
writeRaster(npp2k_su_sd_shelf, filename=paste0(env.derived,string.chr,"2km_shelf_NPP_Cafe_filled_SummerStandardDeviation.tif"), overwrite=TRUE)


#### 5) ROMS Currents & Temperature & FAM ----

## 4k models for now, 2k to follow, and proper FAM to follow!!
data.dat100 <- "E:/science/data_environmental/Circumpolar_ROMS/4km_outputs/output_sed_float_test8/"
#### load lon/lat information from ROMS-grid
grd4k_nc <- nc_open(paste0(env.raw,"waom4extend_grd.nc"))
lon_rho <- ncvar_get(grd4k_nc, varid="lon_rho")
lat_rho <- ncvar_get(grd4k_nc, varid="lat_rho")
#### Prepare empty rasters to assign correct projected values to
roms.coords.proj <- rgdal::project(cbind(c(lon_rho), c(lat_rho)), proj=stereo)
x.range <- c(min(roms.coords.proj[,1])-2000,max(roms.coords.proj[,1])+2000)
y.range <- c(min(roms.coords.proj[,2])-2000,max(roms.coords.proj[,2])+2000)
empty.roms.ra <- raster(extent(c(x.range,y.range)), crs=stereo, resolution=4000)

#depth
h <- raster(paste0(data.dat100,"ocean_avg_0001.nc"), varname="h", level=1)
#salinity
salt <- raster(paste0(data.dat100,"ocean_avg_0001.nc"), varname="salt", level=1)
#temperature
temp <- raster(paste0(data.dat100,"ocean_avg_0001.nc"), varname="temp", level=1)
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
u <- v <- u.abs <- v.abs <- u31 <- v31 <- sa <- te <- se <- su <- fl <- empty.roms.ra
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

#settling FAM
susp_08_full <- brick(paste0(data.dat100,"ocean_his_0001.nc"), varname="sand_08", level=1)*86400
settle_08_full <- brick(paste0(data.dat100,"ocean_his_0001.nc"), varname="sandfrac_08", level=1)*86400
susp_08 <- susp_08_full[[31]]
settle_08 <- settle_08_full[[31]]-settle_08_full[[1]]
flux_08 <- susp_08-settle_08

# susp_his_08 <-  brick(empty.roms.ra,nl=nlayers(u.raw))
# settle_08_full <- brick(empty.roms.ra,nl=nlayers(u.raw))
# #settle6 <- raster(paste0(data.dat100,"ocean_avg_0001.nc"), varname="sand_06", level=1)
# susp_his_08[] <- brick(paste0(data.dat100,"ocean_his_0001.nc"), varname="sand_08", level=1)[]*86400
# settle_08_full[] <- brick(paste0(data.dat100,"ocean_his_0001.nc"), varname="sandfrac_08", level=1)[]*86400
# susp_08 <- susp_his_08[[31]]
# settle_08 <- settle_08_full[[31]]-settle_08_full[[1]]
# flux_08 <- susp_08-settle_08

## remove inland values for depth
h2 <- empty.roms.ra
h2[] <- h[]
h2[is.na(mean.uv)] <- NA
sa[] <- extract(salt, coordinates(h), method="bilinear")
te[] <- extract(temp, coordinates(h), method="bilinear")
se[] <- extract(settle_08, coordinates(h), method="bilinear")
su[] <- extract(susp_08, coordinates(h), method="bilinear")
fl[] <- extract(flux_08, coordinates(h), method="bilinear")

## resample to standard 500m resolution of other environmental variables
mean.uv_500 <- resample(mean.uv,r)
abs.uv_500 <- resample(abs.uv,r)
res.uv_500 <- resample(res.uv,r)
t_500 <- resample(te,r)
s_500 <- resample(sa,r)
settle_08_500 <- resample(se,r)
susp_08_500 <- resample(su,r)
flux_08_500 <- resample(fl,r)

## shelf only
mean.uv_500_shelf <- mean.uv_500
abs.uv_500_shelf <- abs.uv_500
res.uv_500_shelf <- res.uv_500
t_500_shelf <- t_500
s_500_shelf <- s_500
settle_08_500_shelf <- settle_08_500
susp_08_500_shelf <- susp_08_500
flux_08_500_shelf <- flux_08_500

mean.uv_500_shelf[is.na(r)] <- NA
abs.uv_500_shelf[is.na(r)] <- NA
res.uv_500_shelf[is.na(r)] <- NA
t_500_shelf[is.na(r)] <- NA
s_500_shelf[is.na(r)] <- NA
settle_08_500_shelf[is.na(r)] <- NA
susp_08_500_shelf[is.na(r)] <- NA
flux_08_500_shelf[is.na(r)] <- NA

## write rasters to file
writeRaster(mean.uv_500,      overwrite=TRUE, filename=paste0(env.derived,string.chr,"500m_waom4k_seafloorcurrents_mean.tif"))
writeRaster(mean.uv_500_shelf,overwrite=TRUE, filename=paste0(env.derived,string.chr,"500m_shelf_waom4k_seafloorcurrents_mean.tif"))
writeRaster(abs.uv_500,       overwrite=TRUE, filename=paste0(env.derived,string.chr,"500m_waom4k_seafloorcurrents_absolute.tif"))
writeRaster(abs.uv_500_shelf, overwrite=TRUE, filename=paste0(env.derived,string.chr,"500m_shelf_waom4k_seafloorcurrents_absolute.tif"))
writeRaster(res.uv_500,       overwrite=TRUE, filename=paste0(env.derived,string.chr,"500m_waom4k_seafloorcurrents_residual.tif"))
writeRaster(res.uv_500_shelf, overwrite=TRUE, filename=paste0(env.derived,string.chr,"500m_shelf_waom4k_seafloorcurrents_residual.tif"))
writeRaster(temp,             overwrite=TRUE, filename=paste0(env.derived,string.chr,"waom4k_seafloortemperature.tif"))
writeRaster(t_500,            overwrite=TRUE, filename=paste0(env.derived,string.chr,"500m_waom4k_seafloortemperature.tif"))
writeRaster(t_500_shelf,      overwrite=TRUE, filename=paste0(env.derived,string.chr,"500m_shelf_waom4k_seafloortemperature.tif"))
writeRaster(salt,             overwrite=TRUE, filename=paste0(env.derived,string.chr,"waom4k_seafloorsalinity.tif"))
writeRaster(s_500,            overwrite=TRUE, filename=paste0(env.derived,string.chr,"500m_waom4k_seafloorsalinity.tif"))
writeRaster(s_500_shelf,      overwrite=TRUE, filename=paste0(env.derived,string.chr,"500m_shelf_waom4k_seafloorsalinity.tif"))
writeRaster(settle_08,          overwrite=TRUE, filename=paste0(env.derived,string.chr,"waom4k_test_settle08.tif"))
writeRaster(settle_08_500,      overwrite=TRUE, filename=paste0(env.derived,string.chr,"500m_waom4k_test_settle08.tif"))
writeRaster(settle_08_500_shelf,overwrite=TRUE, filename=paste0(env.derived,string.chr,"500m_shelf_waom4k_test_settle08.tif"))
writeRaster(susp_08,            overwrite=TRUE, filename=paste0(env.derived,string.chr,"waom4k_test_susp08.tif"))
writeRaster(susp_08_500,        overwrite=TRUE, filename=paste0(env.derived,string.chr,"500m_waom4k_test_susp08.tif"))
writeRaster(susp_08_500_shelf,  overwrite=TRUE, filename=paste0(env.derived,string.chr,"500m_shelf_waom4k_test_susp08.tif"))
writeRaster(flux_08,            overwrite=TRUE, filename=paste0(env.derived,string.chr,"waom4k_test_flux08.tif"))
writeRaster(flux_08_500,        overwrite=TRUE, filename=paste0(env.derived,string.chr,"500m_waom4k_test_flux08.tif"))
writeRaster(flux_08_500_shelf,  overwrite=TRUE, filename=paste0(env.derived,string.chr,"500m_shelf_waom4k_test_flux08.tif"))


## resample to standard 500m resolution of other environmental variables
mean.uv_2k <- resample(mean.uv,r2k.depth)
abs.uv_2k <-  resample(abs.uv,r2k.depth)
res.uv_2k <-  resample(res.uv,r2k.depth)
t_2k <-       resample(te,r2k.depth)
s_2k <-       resample(sa,r2k.depth)
settle_08_2k<-resample(se,r2k.depth)
susp_08_2k <- resample(su,r2k.depth)
flux_08_2k <- resample(fl,r2k.depth)

## shelf only
mean.uv_2k_shelf <- mean.uv_2k
abs.uv_2k_shelf <- abs.uv_2k
res.uv_2k_shelf <- res.uv_2k
t_2k_shelf <- t_2k
s_2k_shelf <- s_2k
settle_08_2k_shelf <- settle_08_2k
susp_08_2k_shelf <- susp_08_2k
flux_08_2k_shelf <- flux_08_2k

mean.uv_2k_shelf[is.na(r2k.depth)] <- NA
abs.uv_2k_shelf[is.na(r2k.depth)] <- NA
res.uv_2k_shelf[is.na(r2k.depth)] <- NA
t_2k_shelf[is.na(r2k.depth)] <- NA
s_2k_shelf[is.na(r2k.depth)] <- NA
settle_08_2k_shelf[is.na(r2k.depth)] <- NA
susp_08_2k_shelf[is.na(r2k.depth)] <- NA
flux_08_2k_shelf[is.na(r2k.depth)] <- NA

## write rasters to file
writeRaster(mean.uv_2k,      overwrite=TRUE, filename=paste0(env.derived,string.chr,"2km_waom4k_seafloorcurrents_mean.tif"))
writeRaster(mean.uv_2k_shelf,overwrite=TRUE, filename=paste0(env.derived,string.chr,"2km_shelf_waom4k_seafloorcurrents_mean.tif"))
writeRaster(abs.uv_2k,       overwrite=TRUE, filename=paste0(env.derived,string.chr,"2km_waom4k_seafloorcurrents_absolute.tif"))
writeRaster(abs.uv_2k_shelf, overwrite=TRUE, filename=paste0(env.derived,string.chr,"2km_shelf_waom4k_seafloorcurrents_absolute.tif"))
writeRaster(res.uv_2k,       overwrite=TRUE, filename=paste0(env.derived,string.chr,"2km_waom4k_seafloorcurrents_residual.tif"))
writeRaster(res.uv_2k_shelf, overwrite=TRUE, filename=paste0(env.derived,string.chr,"2km_shelf_waom4k_seafloorcurrents_residual.tif"))
writeRaster(t_2k,            overwrite=TRUE, filename=paste0(env.derived,string.chr,"2km_waom4k_seafloortemperature.tif"))
writeRaster(t_2k_shelf,      overwrite=TRUE, filename=paste0(env.derived,string.chr,"2km_shelf_waom4k_seafloortemperature.tif"))
writeRaster(s_2k,            overwrite=TRUE, filename=paste0(env.derived,string.chr,"2km_waom4k_seafloorsalinity.tif"))
writeRaster(s_2k_shelf,      overwrite=TRUE, filename=paste0(env.derived,string.chr,"2km_shelf_waom4k_seafloorsalinity.tif"))
writeRaster(settle_08_2k,      overwrite=TRUE, filename=paste0(env.derived,string.chr,"2km_waom4k_test_settle08.tif"))
writeRaster(settle_08_2k_shelf,overwrite=TRUE, filename=paste0(env.derived,string.chr,"2km_shelf_waom4k_test_settle08.tif"))
writeRaster(susp_08_2k,        overwrite=TRUE, filename=paste0(env.derived,string.chr,"2km_waom4k_test_susp08.tif"))
writeRaster(susp_08_2k_shelf,  overwrite=TRUE, filename=paste0(env.derived,string.chr,"2km_shelf_waom4k_test_susp08.tif"))
writeRaster(flux_08_2k,        overwrite=TRUE, filename=paste0(env.derived,string.chr,"2km_waom4k_test_flux08.tif"))
writeRaster(flux_08_2k_shelf,  overwrite=TRUE, filename=paste0(env.derived,string.chr,"2km_shelf_waom4k_test_flux08.tif"))


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

