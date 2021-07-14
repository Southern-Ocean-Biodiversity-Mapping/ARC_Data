library(magick)
library(geosphere)
library(ggplot2)
library(MBHdesign)
library(raster)
library(lubridate)
library(raadtools)
library(readxl)
'%!in%' <- function(x,y)!('%in%'(x,y))

## environmental data for plotting
env.dir <- "C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_environmental/"

depth <- raster(paste0(env.dir,"Circumpolar_EnvData_bathy500m_shelf_gebco2020_depth.grd"))
load(paste0(env.dir,"Circumpolar_Coastline_W.Rdata"))

bio.path <- "C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/Stills/ANA08D/"

image.dir <- "D:/ARC_DP_data/a_RawData_DirectFromContributors/ANA08D/"

dat.raw <- read.csv(paste0(bio.path,"ANA08D_yoyo_elog.csv"))
dat <- dat.raw[,c(3,7,8)]
names(dat)[2:3] <- c("Latitude","Longitude")

spatial.dat <- data.frame(dat[,c(3,2)])
coordinates(spatial.dat) <- c("Longitude","Latitude")
proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
polar.dat.ANA08D.raw <- spTransform(spatial.dat, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
save(polar.dat.ANA08D.raw, file="ANA08D_locations.Rdata")



















































image.dir_2 <- "D:/ARC_DP_data/a_RawData_DirectFromContributors/JR15005/"
image.dir_3 <- "D:/ARC_DP_data/a_RawData_DirectFromContributors/JR17001/"
image.dir_4 <- "D:/ARC_DP_data/a_RawData_DirectFromContributors/JR17003/"

path.bad.images_1 <- "D:/ARC_DP_data/adjusted_SUCS/JR262/"
path.bad.images_2 <- "D:/ARC_DP_data/adjusted_SUCS/JR15005/"
path.bad.images_3 <- "D:/ARC_DP_data/adjusted_SUCS/JR17001/"
path.bad.images_4 <- "D:/ARC_DP_data/adjusted_SUCS/JR17003/"

dat.JR262.raw <- read_xlsx(paste0(image.dir_1,"JR262_SUCS_SouthOrkneys.xlsx"))
dat.JR262.raw <- dat.JR262.raw[-c(2,4,6),]
dat.JR262.raw[,c(6,5)] <- dat.JR262.raw[,c(6,5)]*-1
dat.JR15005.raw <- read.table(paste0(image.dir_2,"JR15005_SUCS.txt"),sep="\t", header=TRUE)
names(dat.JR15005.raw)[3:4] <- c("Latitude","Longitude")
dat.JR17001.raw <- read_xlsx(paste0(image.dir_3,"JR17001_SUCS.xlsx"))
dat.JR17003.raw <- read.csv(paste0(image.dir_4,"JR17003_SUCS.csv"))
dat.JR17003.raw$Event.No..Built.In...Integer.[which(dat.JR17003.raw$Event.No..Built.In...Integer.=="O11")] <- 11
names(dat.JR17003.raw)[9:10] <- c("Latitude","Longitude")
lonlat <- rbind((dat.JR262.raw[,c(6,5)]),dat.JR15005.raw[,c(4,3)],dat.JR17001.raw[,c(6,5)],dat.JR17003.raw[,c(10,9)])

# spatial.dat <- data.frame(lonlat)
# spatial.dat <- spatial.dat[-which(is.na(rowSums(spatial.dat))),]
# coordinates(spatial.dat) <- c("Longitude","Latitude")
# proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
# polar.dat <- spTransform(spatial.dat, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# 
# plot(r2, xlim=c(-3000000,-1500000), ylim=c(700000,2500000), col=terrain.colors(99))
# plot(coast.proj, add=TRUE)
# points(polar.dat, pch=16, col="blue")
# 
# plot(r2, xlim=c(-2450000,-2050000), ylim=c(2100000,2450000), col=terrain.colors(99))
# plot(coast.proj, add=TRUE)
# points(polar.dat, pch=16, col="blue")
# 
# plot(r2, xlim=c(-2415000,-2370000), ylim=c(2245000,2270000), col=terrain.colors(99))
# points(polar.dat, pch=16, col="blue")

## JR267: lonlat per site, take 1 random image per site
## JR15005:
## JR17001: lonlat per image, but coded differently, 3 images per transect?
## JR17003: lonlat per image, but coded differently


# ## JR262
# spatial.dat <- data.frame(dat.JR262.raw[,c(6,5)])
# coordinates(spatial.dat) <- c("Longitude","Latitude")
# proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
# polar.dat.JR262.raw <- spTransform(spatial.dat, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# #save(polar.dat.JR262.raw, file="JR262_locations.Rdata")
# 
# plot(r2, xlim=c(-2500000,-2000000), ylim=c(2000000,2500000), col=terrain.colors(99))
# plot(coast.proj, add=TRUE)
# points(polar.dat.JR262.raw, pch=16, col="blue")
# 
# ## JR15005
# spatial.dat <- data.frame(dat.JR15005.raw[,c(4,3)])
# coordinates(spatial.dat) <- c("Longitude","Latitude")
# proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
# polar.dat.JR15005.raw <- spTransform(spatial.dat, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# #save(polar.dat.JR15005.raw, file="JR15005_locations.Rdata")
# 
# ## JR17001 & JR17003
# spatial.dat <- data.frame(dat.JR17001.raw[,c(6,5)])
# spatial.dat <- spatial.dat[-which(is.na(rowSums(spatial.dat))),]
# coordinates(spatial.dat) <- c("Longitude","Latitude")
# proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
# polar.dat.JR17001.raw <- spTransform(spatial.dat, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# 
# spatial.dat <- data.frame(dat.JR17003.raw[,c(10,9)])
# spatial.dat <- spatial.dat[-which(is.na(rowSums(spatial.dat))),]
# coordinates(spatial.dat) <- c("Longitude","Latitude")
# proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
# polar.dat.JR17003.raw <- spTransform(spatial.dat, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# 
# plot(r2, xlim=c(-2700000,-2300000), ylim=c(700000,1700000), col=terrain.colors(99))
# plot(coast.proj, add=TRUE)
# points(polar.dat.JR17001.raw, pch=16, col="blue")
# points(polar.dat.JR17003.raw, pch=16, col="cyan")
# 
# #save(polar.dat.JR17001.raw, file="JR17001_locations.Rdata")
# #save(polar.dat.JR17003.raw, file="JR17003_locations.Rdata")
# 
# plot(r2, xlim=c(-2635000,-2630000), ylim=c(1592000,1596000), col=terrain.colors(99))
# plot(coast.proj, add=TRUE)
# points(polar.dat.JR17001.raw, pch=16, col="blue")
# points(polar.dat.JR17003.raw, pch=16, col="cyan")


## dataframe with Filename, lon, lat, transectID, surveyID



