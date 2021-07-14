##### read in circumpolar environmental data

library(SOmap)
library(ncdf4)
library(raadtools)
library(spatialEco)
library(raster)
library(sp)
library(dplyr)
library(blueant)
library(raadtools)
library(ggplot2)

##### specify directories

env.raw <- "C:/Users/jjansen/Desktop/science/data_environmental/raw/"
env.derived <- "C:/Users/jjansen/Desktop/science/data_environmental/derived/"
my_data_dir <- "C:/Users/jjansen/Desktop/science/data_environmental/raw/accessed_through_R"

###################################################
####### bathymetry and seascapes
###################################################

##### to download data from repositories
# ## FROM https://ropensci.org/blog/2018/11/13/antarctic/
# src <- bind_rows(sources("NSIDC SMMR-SSM/I Nasateam sea ice concentration", hemisphere = "south", time_resolutions = "day", years = 2013),
#                 sources("Southern Ocean summer chlorophyll-a climatology (Johnson)"),
#                 sources("IBCSO bathymetry"),
#                 bb_modify_source(sources("Oceandata MODIS Aqua Level-3 mapped daily 4km chl-a"), method = list(search = "A2014*L3m_DAY_CHL_chlor_a_4km.nc")))
# result <- bb_get(src, local_file_root = my_data_dir, clobber = 0, verbose = TRUE, confirm = NULL)

##### load standard IBCSO bathymetry (2013 version)
set_data_roots(my_data_dir)

## bathy
r <- readtopo("ibcso")

r.depth <- r
r.slope <- terrain(r)
r.tpi <- tpi(r)
r.tpi5 <- tpi(r, scale=5)
r.tpi21 <- tpi(r, scale=21)
r.tpi31 <- tpi(r, scale=31)

r.depth[r>0] <- NA
r.depth[r<(-3000)] <- NA ## comment this line to get a raster across the full depth
r.slope[is.na(r.depth[])] <- NA
r.tpi[is.na(r.depth[])] <- NA
r.tpi5[is.na(r.depth[])] <- NA
r.tpi21[is.na(r.depth[])] <- NA
r.tpi31[is.na(r.depth[])] <- NA

# xlim=c(128,148)
# ylim=c(-67.5,-64)
par(mfrow=c(2,2))
plot(r.depth) #, xlim=xlim, ylim=ylim
plot(r.slope) #, xlim=xlim, ylim=ylim
plot(r.tpi5)  #, xlim=xlim, ylim=ylim
plot(r.tpi31) #, xlim=xlim, ylim=ylim

r.bathy <- stack(r.depth, r.slope, r.tpi, r.tpi5, r.tpi21, r.tpi31)
names(r.bathy) <- c("depth","slope","tpi","tpi5","tpi21","tpi31")
writeRaster(r.bathy, filename=paste0(env.derived,"Circumpolar_EnvData_bathy.Rdata"))
#save(r.bathy, file=paste0(env.derived,"Circumpolar_EnvData_bathy.Rdata"))

## as mentioned above, for a full depth raster comment the line "r.depth[r<(-3000)] <- NA"
#writeRaster(r.bathy, filename=paste0(env.derived,"Circumpolar_EnvData_bathy_fulldepth.Rdata"))
#save(r.bathy, file=paste0(env.derived,"Circumpolar_EnvData_bathy_fulldepth.Rdata"))

## creating a smaller file for the Weddell Sea analysis only
r.bathy_smaller <- crop(r.bathy, y=extent(-2800000,0,400000,2250000))
writeRaster(r.bathy_smaller, filename=paste0(env.derived,"Circumpolar_EnvData_bathy_Weddell.Rdata"))
# r.bathy_smaller <- raster(paste0(env.derived,"Circumpolar_EnvData_bathy_Weddell.grd"))
plot(r.bathy_smaller$depth,xlim=c(-2800000,-500000),ylim=c(500000,2000000))

##### load newer bathymetry data from GEBCO and reproject gebco2020 grid to polar stereographic
## original (unprojected) data
g <- raster(paste0(env.raw,"GEBCO_2020/gebco_2020_SO.tif"))

## data projected with original resolution
#g2 <- projectRaster(g,crs="+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")

## data projected and resampled to 500m grid cells
g2_500 <- projectRaster(g,r)

#r <- g ## unprojected data
#r <- g2 ## original resolution
r <- g2_500 ## 500m resolution
r.depth <- r
r.slope <- terrain(r)
r.tpi <- tpi(r)
r.tpi5 <- tpi(r, scale=5)
r.tpi11 <- tpi(r, scale=11)
#r.tpi21 <- tpi(r, scale=21)
#r.tpi31 <- tpi(r, scale=31)

r.depth[r>0] <- NA
r.depth[r<(-2500)] <- NA ## comment this line to get a raster across the full depth
r.slope[is.na(r.depth[])] <- NA
r.tpi[is.na(r.depth[])] <- NA
r.tpi5[is.na(r.depth[])] <- NA
r.tpi11[is.na(r.depth[])] <- NA
#r.tpi21[is.na(r.depth[])] <- NA
#r.tpi31[is.na(r.depth[])] <- NA

# xlim=c(128,148)
# ylim=c(-67.5,-64)
par(mfrow=c(2,2))
plot(r.depth) #, xlim=xlim, ylim=ylim
plot(r.slope) #, xlim=xlim, ylim=ylim
plot(r.tpi5)  #, xlim=xlim, ylim=ylim
plot(r.tpi31) #, xlim=xlim, ylim=ylim

# r.bathy <- stack(r.depth, r.slope, r.tpi, r.tpi5, r.tpi11)#, r.tpi31)
# names(r.bathy) <- c("depth","slope","tpi","tpi5","tpi11")#,"tpi31")

# writeRaster(r.bathy, filename=paste0(env.derived,"Circumpolar_EnvData_bathy500m_gebco2020.Rdata"))
# save(r.bathy, file=paste0(env.derived,"Circumpolar_EnvData_bathy500m_gebco2020.Rdata"))
# r <- raster(paste0(env.derived,"Circumpolar_EnvData_bathy500m_gebco2020.grd"))

#writeRaster(g2, filename=paste0(env.derived,"Circumpolar_EnvData_gebco2020_fulldepth.Rdata"), overwrite=TRUE)
#save(g2, file=paste0(env.derived,"Circumpolar_EnvData_gebco2020_fulldepth.Rdata"))
load(paste0(env.derived,"Circumpolar_EnvData_gebco2020_fulldepth.Rdata"))

## original projected resolution (119m x 460m)
writeRaster(r.depth, filename=paste0(env.dir,"Circumpolar_EnvData_bathy_gebco2020_depth.Rdata"))
writeRaster(r.slope, filename=paste0(env.dir,"Circumpolar_EnvData_bathy_gebco2020_slope.Rdata"))
writeRaster(r.tpi, filename=paste0(env.dir,"Circumpolar_EnvData_bathy_gebco2020_tpi.Rdata"))
writeRaster(r.tpi5, filename=paste0(env.dir,"Circumpolar_EnvData_bathy_gebco2020_tpi5.Rdata"))
writeRaster(r.tpi11, filename=paste0(env.dir,"Circumpolar_EnvData_bathy_gebco2020_tpi11.Rdata"))

## 500m resolution
writeRaster(r.depth, filename=paste0(env.dir,"Circumpolar_EnvData_bathy500m_gebco2020_depth.Rdata"), overwrite=TRUE)
writeRaster(r.slope, filename=paste0(env.dir,"Circumpolar_EnvData_bathy500m_gebco2020_slope.Rdata"), overwrite=TRUE)
writeRaster(r.tpi, filename=paste0(env.dir,"Circumpolar_EnvData_bathy500m_gebco2020_tpi.Rdata"), overwrite=TRUE)
writeRaster(r.tpi5, filename=paste0(env.dir,"Circumpolar_EnvData_bathy500m_gebco2020_tpi5.Rdata"), overwrite=TRUE)
writeRaster(r.tpi11, filename=paste0(env.dir,"Circumpolar_EnvData_bathy500m_gebco2020_tpi11.Rdata"), overwrite=TRUE)

## 500m resolution, restricted to the shelf only
writeRaster(r.depth, filename=paste0(env.dir,"Circumpolar_EnvData_bathy500m_shelf_gebco2020_depth.Rdata"), overwrite=TRUE)
writeRaster(r.slope, filename=paste0(env.dir,"Circumpolar_EnvData_bathy500m_shelf_gebco2020_slope.Rdata"), overwrite=TRUE)
writeRaster(r.tpi, filename=paste0(env.dir,"Circumpolar_EnvData_bathy500m_shelf_gebco2020_tpi.Rdata"), overwrite=TRUE)
writeRaster(r.tpi5, filename=paste0(env.dir,"Circumpolar_EnvData_bathy500m_shelf_gebco2020_tpi5.Rdata"), overwrite=TRUE)
writeRaster(r.tpi11, filename=paste0(env.dir,"Circumpolar_EnvData_bathy500m_shelf_gebco2020_tpi11.Rdata"), overwrite=TRUE)


###################################################
####### primary production
###################################################
npp_2002 <- get(load(paste0(env.dir,"Circumpolar_NPP_Cafe_filled_2002.Rdata")))
npp_2003 <- get(load(paste0(env.dir,"Circumpolar_NPP_Cafe_filled_2003.Rdata")))
npp_2004 <- get(load(paste0(env.dir,"Circumpolar_NPP_Cafe_filled_2004.Rdata")))
npp_2005 <- get(load(paste0(env.dir,"Circumpolar_NPP_Cafe_filled_2005.Rdata")))
npp_2006 <- get(load(paste0(env.dir,"Circumpolar_NPP_Cafe_filled_2006.Rdata")))
npp_2007 <- get(load(paste0(env.dir,"Circumpolar_NPP_Cafe_filled_2007.Rdata")))
npp_2008 <- get(load(paste0(env.dir,"Circumpolar_NPP_Cafe_filled_2008.Rdata")))
npp_2009 <- get(load(paste0(env.dir,"Circumpolar_NPP_Cafe_filled_2009.Rdata")))
npp_2010 <- get(load(paste0(env.dir,"Circumpolar_NPP_Cafe_filled_2010.Rdata")))

npp_su_stack <- stack(npp_2002$cafe.2002335,npp_2003$cafe.2003001,npp_2003$cafe.2003032,npp_2003$cafe.2003060,
                      npp_2003$cafe.2003335,npp_2004$cafe.2004001,npp_2004$cafe.2004032,npp_2004$cafe.2004061,
                      npp_2004$cafe.2004336,npp_2005$cafe.2005001,npp_2005$cafe.2005032,npp_2005$cafe.2005060,
                      npp_2005$cafe.2005335,npp_2006$cafe.2006001,npp_2006$cafe.2006032,npp_2006$cafe.2006060,
                      npp_2006$cafe.2006335,npp_2007$cafe.2007001,npp_2007$cafe.2007032,npp_2007$cafe.2007060,
                      npp_2007$cafe.2007335,npp_2008$cafe.2008001,npp_2008$cafe.2008032,npp_2008$cafe.2008061,
                      npp_2008$cafe.2008336,npp_2009$cafe.2009001,npp_2009$cafe.2009032,npp_2009$cafe.2009060,
                      npp_2009$cafe.2009335,npp_2010$cafe.2010001,npp_2010$cafe.2010032,npp_2010$cafe.2010060)

npp_su_ave <- mean(npp_su_stack, na.rm=TRUE)
npp_su <- projectRaster(npp_su_ave,depth)
writeRaster(npp_su, filename=paste0(env.dir,"Circumpolar_EnvData_NPP_SummerAverage.Rdata"))
# npp[is.na(r.depth)] <- NA
# writeRaster(npp, filename=paste0(env.dir,"Circumpolar_EnvData_NPP_SummerAverage_shelf.Rdata"))


###################################################
####### ocean currents & temperature 2k res (UPDATE ONCE FAM HAS RUN)
###################################################
## file paths
f.grd <- paste0(env.dir,"Circumpolar_ROMS/waom2_grd.nc")
f.u <- paste0(env.dir,"Circumpolar_ROMS/ocean_avg_0538-0610_u_avg.nc")
f.v <- paste0(env.dir,"Circumpolar_ROMS/ocean_avg_0538-0610_v_avg.nc")
f.t <- paste0(env.dir,"Circumpolar_ROMS/ocean_avg_0538-0610_temp_avg.nc")
## read as raster
lon <- raster(f.grd, varname="lon_rho")
lat <- raster(f.grd, varname="lat_rho")
u <- raster(f.u, lvar=3, level=1, varname="u")
v <- raster(f.v, lvar=3, level=1, varname="v")
t <- raster(f.t, lvar=3, level=1, varname="temp")
## bring to same extent (address one of the quirks of ROMS)
ext <- extent(1,3149,1,2649)
lon <- crop(lon,y=ext,snap="out")
lat <- crop(lat,y=ext,snap="out")
u <- crop(u,y=ext,snap="out")
v <- crop(v,y=ext,snap="out")
t <- crop(t,y=ext,snap="out")
#w <- crop(w,y=ext,snap="out")
## calculate a single seafloor current speed value
uv <- sqrt(abs(u)^2+abs(v)^2)
## projection and extent for the raster (netcdf files were already polar-projected with true south at -71S)
crs <- "+proj=stere +lat_ts=-71 +lat_0=-90 +datum=WGS84"
pts <- rgdal::project(cbind(values(lon), values(lat)), crs)
ex <- extent(pts)
uv <- setExtent(uv, ex)
t <- setExtent(t, ex)
projection(uv) <- crs
projection(t) <- crs
## resample to standard 500m resolution of other environmental variables
uv_500 <- resample(uv,r)
t_500 <- resample(t,r)
## shelf only
uv_500_shelf <- uv_500
uv_500_shelf[is.na(depth)] <- NA
t_500_shelf <- t_500
t_500_shelf[is.na(depth)] <- NA
## write rasters to file
writeRaster(uv, filename=paste0(env.dir,"Circumpolar_EnvData_waom2k_seafloorcurrents.Rdata"))
writeRaster(uv_500, filename=paste0(env.dir,"Circumpolar_EnvData_waom2k_seafloorcurrents_500mInterpolation.Rdata"))
writeRaster(uv_500_shelf, filename=paste0(env.dir,"Circumpolar_EnvData_waom2k_seafloorcurrents_500mInterpolation_shelf.Rdata"))
writeRaster(t, filename=paste0(env.dir,"Circumpolar_EnvData_waom2k_seafloortemperature.Rdata"))
writeRaster(t_500, filename=paste0(env.dir,"Circumpolar_EnvData_waom2k_seafloortemperature_500mInterpolation.Rdata"))
writeRaster(t_500_shelf, filename=paste0(env.dir,"Circumpolar_EnvData_waom2k_seafloortemperature_500mInterpolation_shelf.Rdata"))
# plot(uv,xlim=c(0,1000000),ylim=c(-2300000,-1500000))
# plot(coast.proj, add=TRUE)
# contour(r, add=TRUE)


###################################################
####### seaice
###################################################
SO.ext <- extent(-180, 180, -80, -55)

## create a raster with the native extent of the file
timespan <- as.Date(c("2013-01-01","2013-01-02","2013-01-03"))
full.extent.raster <- readice(product="nsidc", date=timespan)
full.extent.projection <- projection(full.extent.raster)

####### SEASONAL
## dates
sp_dates <- list()
su_dates <- list()
au_dates <- list()
sp_dates[[1]] <- seq(as.Date("2002-09-23"), as.Date("2002-12-20"), by = "1 day")
sp_dates[[2]] <- seq(as.Date("2003-09-23"), as.Date("2003-12-20"), by = "1 day")
sp_dates[[3]] <- seq(as.Date("2004-09-23"), as.Date("2004-12-20"), by = "1 day")
sp_dates[[4]] <- seq(as.Date("2005-09-23"), as.Date("2005-12-20"), by = "1 day")
sp_dates[[5]] <- seq(as.Date("2006-09-23"), as.Date("2006-12-20"), by = "1 day")
sp_dates[[6]] <- seq(as.Date("2007-09-23"), as.Date("2007-12-20"), by = "1 day")
sp_dates[[7]] <- seq(as.Date("2008-09-23"), as.Date("2008-12-20"), by = "1 day")
sp_dates[[8]] <- seq(as.Date("2009-09-23"), as.Date("2009-12-20"), by = "1 day")
sp_dates[[9]] <- seq(as.Date("2010-09-23"), as.Date("2010-12-20"), by = "1 day")
sp_dates[[10]] <- seq(as.Date("2011-09-23"), as.Date("2011-12-20"), by = "1 day")
sp_dates[[11]] <- seq(as.Date("2012-09-23"), as.Date("2012-12-20"), by = "1 day")
sp_dates[[12]] <- seq(as.Date("2013-09-23"), as.Date("2013-12-20"), by = "1 day")
sp_dates[[13]] <- seq(as.Date("2014-09-23"), as.Date("2014-12-20"), by = "1 day")
sp_dates[[14]] <- seq(as.Date("2015-09-23"), as.Date("2015-12-20"), by = "1 day")
sp_dates[[15]] <- seq(as.Date("2016-09-23"), as.Date("2016-12-20"), by = "1 day")
sp_dates[[16]] <- seq(as.Date("2017-09-23"), as.Date("2017-12-20"), by = "1 day")
sp_dates[[17]] <- seq(as.Date("2018-09-23"), as.Date("2018-12-20"), by = "1 day")
sp_dates[[18]] <- seq(as.Date("2019-09-23"), as.Date("2019-12-20"), by = "1 day")

su_dates[[1]] <- seq(as.Date("2002-12-21"), as.Date("2003-03-20"), by = "1 day")
su_dates[[2]] <- seq(as.Date("2003-12-21"), as.Date("2004-03-20"), by = "1 day")
su_dates[[3]] <- seq(as.Date("2004-12-21"), as.Date("2005-03-20"), by = "1 day")
su_dates[[4]] <- seq(as.Date("2005-12-21"), as.Date("2006-03-20"), by = "1 day")
su_dates[[5]] <- seq(as.Date("2006-12-21"), as.Date("2007-03-20"), by = "1 day")
su_dates[[6]] <- seq(as.Date("2007-12-21"), as.Date("2008-03-20"), by = "1 day")
su_dates[[7]] <- seq(as.Date("2008-12-21"), as.Date("2009-03-20"), by = "1 day")
su_dates[[8]] <- seq(as.Date("2009-12-21"), as.Date("2010-03-20"), by = "1 day")
su_dates[[9]] <- seq(as.Date("2010-12-21"), as.Date("2011-03-20"), by = "1 day")
su_dates[[10]] <- seq(as.Date("2011-12-21"), as.Date("2012-03-20"), by = "1 day")
su_dates[[11]] <- seq(as.Date("2012-12-21"), as.Date("2013-03-20"), by = "1 day")
su_dates[[12]] <- seq(as.Date("2013-12-21"), as.Date("2014-03-20"), by = "1 day")
su_dates[[13]] <- seq(as.Date("2014-12-21"), as.Date("2015-03-20"), by = "1 day")
su_dates[[14]] <- seq(as.Date("2015-12-21"), as.Date("2016-03-20"), by = "1 day")
su_dates[[15]] <- seq(as.Date("2016-12-21"), as.Date("2017-03-20"), by = "1 day")
su_dates[[16]] <- seq(as.Date("2017-12-21"), as.Date("2018-03-20"), by = "1 day")
su_dates[[17]] <- seq(as.Date("2018-12-21"), as.Date("2019-03-20"), by = "1 day")

au_dates[[1]] <- seq(as.Date("2003-03-21"), as.Date("2003-06-21"), by = "1 day")
au_dates[[2]] <- seq(as.Date("2004-03-21"), as.Date("2004-06-21"), by = "1 day")
au_dates[[3]] <- seq(as.Date("2005-03-21"), as.Date("2005-06-21"), by = "1 day")
au_dates[[4]] <- seq(as.Date("2006-03-21"), as.Date("2006-06-21"), by = "1 day")
au_dates[[5]] <- seq(as.Date("2007-03-21"), as.Date("2007-06-21"), by = "1 day")
au_dates[[6]] <- seq(as.Date("2008-03-21"), as.Date("2008-06-21"), by = "1 day")
au_dates[[7]] <- seq(as.Date("2009-03-21"), as.Date("2009-06-21"), by = "1 day")
au_dates[[8]] <- seq(as.Date("2010-03-21"), as.Date("2010-06-21"), by = "1 day")
au_dates[[9]] <- seq(as.Date("2011-03-21"), as.Date("2011-06-21"), by = "1 day")
au_dates[[10]] <- seq(as.Date("2012-03-21"), as.Date("2012-06-21"), by = "1 day")
au_dates[[11]] <- seq(as.Date("2013-03-21"), as.Date("2013-06-21"), by = "1 day")
au_dates[[12]] <- seq(as.Date("2014-03-21"), as.Date("2014-06-21"), by = "1 day")
au_dates[[13]] <- seq(as.Date("2015-03-21"), as.Date("2015-06-21"), by = "1 day")
au_dates[[14]] <- seq(as.Date("2016-03-21"), as.Date("2016-06-21"), by = "1 day")
au_dates[[15]] <- seq(as.Date("2017-03-21"), as.Date("2017-06-21"), by = "1 day")
au_dates[[16]] <- seq(as.Date("2018-03-21"), as.Date("2018-06-21"), by = "1 day")
au_dates[[17]] <- seq(as.Date("2019-03-21"), as.Date("2019-06-21"), by = "1 day")

## spring-raster
ice_spring <- stack(readice(sp_dates[[1]], grid = raster(ext=SO.ext, crs=full.extent.projection, res=1/24)))
names(chla_spring)[1] <- paste0("chla_",su_years[1])
for(i in 2:18){
  message(i)
  chla_spring[[i]] <- readchla(sp_dates[[i]], algorithm="johnson", product="MODISA", grid = raster(ext=SO.ext, crs=full.extent.projection, res=1/24))
  names(chla_spring)[i] <- paste0("chla_",years[i])
}
save(chla_spring, file="Circumpolar_EnvData_ChlA_spring.Rdata")
rm(chla_spring)

## summer-raster
chla_summer <- stack(readchla(su_dates[[1]], algorithm="johnson", product="MODISA", grid = raster(ext=SO.ext, crs=full.extent.projection, res=1/24)))
names(chla_summer)[1] <- paste0("chla_",su_years[1])
for(i in 2:17){
  message(i)
  chla_summer[[i]] <- readchla(su_dates[[i]], algorithm="johnson", product="MODISA", grid = raster(ext=SO.ext, crs=full.extent.projection, res=1/24))
  names(chla_summer)[i] <- paste0("chla_",su_years[i])
}
save(chla_summer, file="Circumpolar_EnvData_ChlA_summer.Rdata")
rm(chla_summer)

## autumn-raster
chla_autumn <- stack(readchla(au_dates[[1]], algorithm="johnson", product="MODISA", grid = raster(ext=SO.ext, crs=full.extent.projection, res=1/24)))
names(chla_autumn)[1] <- paste0("chla_",years[2])
for(i in 2:17){
  message(i)
  chla_autumn[[i]] <- readchla(au_dates[[i]], algorithm="johnson", product="MODISA", grid = raster(ext=SO.ext, crs=full.extent.projection, res=1/24))
  names(chla_autumn)[i] <- paste0("chla_",years[i+1])
}
save(chla_autumn, file="Circumpolar_EnvData_ChlA_autumn.Rdata")
rm(chla_autumn)


















#### seascapes
# tpi_small <- r.tpi5
# tpi_large <- r.tpi31
# ## calculate mean:
# m.s <- cellStats(tpi_small, 'mean')
# m.l <- cellStats(tpi_large, 'mean')
# ## calculate sd:
# sd.s <- cellStats(tpi_small, 'sd')
# sd.l <- cellStats(tpi_large, 'sd')
# ## Standardize the TPI grids using the formula:
# ## tpi<sf>_stdi = int((((tpi<sf> – mean) / stdv) * 100) + 0.5)
# t.s_stdi <- calc(tpi_small, function(x) (((x - m.s) / sd.s) * 100) + 0.5)
# t.l_stdi <- calc(tpi_large, function(x) (((x - m.l) / sd.l) * 100) + 0.5)
# ## set low and high values to classify the landforms
# v_low <- -100
# v_high <- 100
# ## specify the function, using ifelse to decide between cases: (x is t31..., y is t201..., z is slope)
# my_fun <- function(x,y,z){
#   ifelse (x > v_low & x < v_high & y > v_low & y < v_high & z <= 5, 5, 
#           ifelse(x > v_low & x < v_high & y > v_low & y < v_high & z >= 5, 6,   #Jeroen:>5, vorher>=6
#                  ifelse(x > v_low & x < v_high & y >= v_high, 7,
#                         ifelse(x > v_low & x < v_high & y <= v_low, 4,
#                                ifelse(x <= v_low & y > v_low & y < v_high, 2,
#                                       ifelse(x >= v_high & y > v_low & y < v_high, 9,
#                                              ifelse(x <= v_low & y >= v_high, 3,
#                                                     ifelse(x <= v_low & y <= v_low, 1,
#                                                            ifelse(x >= v_high & y >= v_high, 10,
#                                                                   ifelse(x >= v_high & y <= v_low, 8, NA)
#                                                            )))))))))
# }
# ## do the calculation (we run my_fun on the three raster layers)
# land_s_l <- overlay(stack(t.s_stdi,t.l_stdi,r.slope), fun=my_fun)
# plot(land_s_l, xlim=c(138,148), ylim=c(-68,-65))


