## Evaluate ROMS outputs
library(raster)
library(colorspace)
library(viridis)
library(ncdf4)
library(raadtools)

env.dir <- "C:/Users/jjansen/Desktop/science/data_environmental/"
env.raw <- paste0(env.dir,"raw/")
env.derived <- paste0(env.dir,"derived/")
#data.dat <- paste0("Circumpolar_ROMS/10km_outputs/sed_test1/")
data.dat <- paste0(env.dir,"Circumpolar_ROMS/10km_outputs/sed_test5/")
AAD_dir <- paste0(env.dir,"raw/accessed_through_R")


xlim=c(250,340)
ylim=c(100,180)

## load projected coastline for plotting
load(paste0(env.derived,"Circumpolar_Coastline.Rdata"))


#### load lon/lat information from ROMS-grid
grd10k_nc <- nc_open(paste0(env.raw,"waom10extend_grd.nc"))
lon_rho <- ncvar_get(grd10k_nc, varid="lon_rho")
lat_rho <- ncvar_get(grd10k_nc, varid="lat_rho")


##### EVALUATE FAM - FLUX OF NPP #####

#### seafloor-layer is 1, while surface-layer is 31

## depth
h <- raster(paste0(data.dat,"ocean_avg_0003.nc"), varname="h", level=1)

## currents
u <- brick(paste0(data.dat,"ocean_avg_0003.nc"), varname="u", level=1)
v <- brick(paste0(data.dat,"ocean_avg_0003.nc"), varname="v", level=1)
uv <- sqrt(u^2+v^2)
# plot(uv[[1]])
# u_31 <- brick(paste0(data.dat,"ocean_avg_0003.nc"), varname="u", level=31)
# v_31 <- brick(paste0(data.dat,"ocean_avg_0003.nc"), varname="v", level=31)
# uv_31 <- sqrt(u_31^2+v_31^2)
# plot(uv_31[[1]])

## surface production
surf_01 <- brick(paste0(data.dat,"ocean_avg_0003.nc"), varname="sand_01", level=31)
surf_02 <- brick(paste0(data.dat,"ocean_avg_0003.nc"), varname="sand_02", level=31)
surf_03 <- brick(paste0(data.dat,"ocean_avg_0003.nc"), varname="sand_03", level=31)
surf_04 <- brick(paste0(data.dat,"ocean_avg_0003.nc"), varname="sand_04", level=31)
surf_05 <- brick(paste0(data.dat,"ocean_avg_0003.nc"), varname="sand_05", level=31)
surf_06 <- brick(paste0(data.dat,"ocean_avg_0003.nc"), varname="sand_06", level=31)
surf_07 <- brick(paste0(data.dat,"ocean_avg_0003.nc"), varname="sand_07", level=31)
surf_08 <- brick(paste0(data.dat,"ocean_avg_0003.nc"), varname="sand_08", level=31)

## bottom ocean layer
susp_01 <- brick(paste0(data.dat,"ocean_avg_0003.nc"), varname="sand_01", level=1)
susp_02 <- brick(paste0(data.dat,"ocean_avg_0003.nc"), varname="sand_02", level=1)
susp_03 <- brick(paste0(data.dat,"ocean_avg_0003.nc"), varname="sand_03", level=1)
susp_04 <- brick(paste0(data.dat,"ocean_avg_0003.nc"), varname="sand_04", level=1)
susp_05 <- brick(paste0(data.dat,"ocean_avg_0003.nc"), varname="sand_05", level=1)
susp_06 <- brick(paste0(data.dat,"ocean_avg_0003.nc"), varname="sand_06", level=1)
susp_07 <- brick(paste0(data.dat,"ocean_avg_0003.nc"), varname="sand_07", level=1)
susp_08 <- brick(paste0(data.dat,"ocean_avg_0003.nc"), varname="sand_08", level=1)

## settled particles
settle_01 <- brick(paste0(data.dat,"ocean_his_0003.nc"), varname="sandfrac_01", level=1)
settle_02 <- brick(paste0(data.dat,"ocean_his_0003.nc"), varname="sandfrac_02", level=1)
settle_03 <- brick(paste0(data.dat,"ocean_his_0003.nc"), varname="sandfrac_03", level=1)
settle_04 <- brick(paste0(data.dat,"ocean_his_0003.nc"), varname="sandfrac_04", level=1)
settle_05 <- brick(paste0(data.dat,"ocean_his_0003.nc"), varname="sandfrac_05", level=1)
settle_06 <- brick(paste0(data.dat,"ocean_his_0003.nc"), varname="sandfrac_06", level=1)
settle_07 <- brick(paste0(data.dat,"ocean_his_0003.nc"), varname="sandfrac_07", level=1)
settle_08 <- brick(paste0(data.dat,"ocean_his_0003.nc"), varname="sandfrac_08", level=1)

##### Overview for Ross Sea
breaks <- seq(-0.000025,0.000025,length.out=99)
breaks.susp <- seq(0,0.000025, length.out=256)
breaks.settle <- seq(0.1245,0.1255, length.out=256)
#ticks.susp <- c(0,0.00001,0.00002)

par(mfrow=c(2,2))
plot(surf_01[[6]], xlim=xlim, ylim=ylim, main="NPP-avg", breaks=breaks.susp)
contour(h, add=TRUE)
plot(susp_01[[6]], xlim=xlim, ylim=ylim, main="FAM-avg-suspended", breaks=breaks.susp)
contour(h, add=TRUE)
plot(settle_01[[6]], xlim=xlim, ylim=ylim, main="FAM-his-settled", breaks=breaks.settle)
contour(h, add=TRUE)
plot(uv[[6]], xlim=xlim, ylim=ylim, main="currents_avg")
contour(h, add=TRUE)
# 
# par(mfrow=c(2,2))
# plot(surf_01[[6]], xlim=xlim, ylim=ylim, main="NPP", breaks=breaks.susp)
# contour(h, add=TRUE)
# plot(susp_01[[6]], xlim=xlim, ylim=ylim, main="FAM-suspended", breaks=breaks.susp)
# contour(h, add=TRUE)
# plot(settle_01[[6]], xlim=xlim, ylim=ylim, main="FAM-settled", breaks=breaks.settle)
# contour(h, add=TRUE)
# plot(uv2[[6]], xlim=xlim, ylim=ylim, main="currents_his")
# contour(h, add=TRUE)

par(mfrow=c(2,2))
plot(h, xlim=xlim, ylim=ylim, main="depth")
plot(uv[[6]], xlim=xlim, ylim=ylim, main="currents")
plot(surf_01[[6]], xlim=xlim, ylim=ylim, main="NPP")
plot(susp_01[[6]], xlim=xlim, ylim=ylim, main="FAM-suspended")

par(mfrow=c(2,2))
plot(surf_05[[6]], xlim=xlim, ylim=ylim, main="NPP")
contour(h, add=TRUE)
plot(sand_05[[6]], xlim=xlim, ylim=ylim, main="FAM")
contour(h, add=TRUE)
plot(uv[[6]], xlim=xlim, ylim=ylim, main="currents")
contour(h, add=TRUE)
plot(sand_05[[6]]-surf_05[[6]], xlim=xlim, ylim=ylim, col=diverge_hsv(99), main="FAM-NPP")
contour(h, add=TRUE)

par(mfrow=c(2,4))
plot(susp_01[[6]]-surf_01[[6]], xlim=xlim, ylim=ylim, breaks=breaks, col=diverge_hsv(99), main="FAM-NPP")
plot(susp_02[[6]]-surf_02[[6]], xlim=xlim, ylim=ylim, breaks=breaks, col=diverge_hsv(99), main="FAM-NPP")
plot(susp_03[[6]]-surf_03[[6]], xlim=xlim, ylim=ylim, breaks=breaks, col=diverge_hsv(99), main="FAM-NPP")
plot(susp_04[[6]]-surf_04[[6]], xlim=xlim, ylim=ylim, breaks=breaks, col=diverge_hsv(99), main="FAM-NPP")
plot(susp_05[[6]]-surf_05[[6]], xlim=xlim, ylim=ylim, breaks=breaks, col=diverge_hsv(99), main="FAM-NPP")
plot(susp_06[[6]]-surf_06[[6]], xlim=xlim, ylim=ylim, breaks=breaks, col=diverge_hsv(99), main="FAM-NPP")
plot(susp_07[[6]]-surf_07[[6]], xlim=xlim, ylim=ylim, breaks=breaks, col=diverge_hsv(99), main="FAM-NPP")
plot(susp_08[[6]]-surf_08[[6]], xlim=xlim, ylim=ylim, breaks=breaks, col=diverge_hsv(99), main="FAM-NPP")

par(mfrow=c(2,2))
plot(surf_01[[6]], xlim=xlim, ylim=ylim, main="NPP")
contour(h, add=TRUE)
plot(sand5_01[[6]], xlim=xlim, ylim=ylim, main="FAM")
contour(h, add=TRUE)
plot(uv[[6]], xlim=xlim, ylim=ylim, main="currents")
contour(h, add=TRUE)
plot(sand5_01[[6]]-surf_01[[6]], xlim=xlim, ylim=ylim, breaks=breaks, col=diverge_hsv(99), main="FAM-NPP")
contour(h, add=TRUE)

##### Evaluate sedimentation patterns:
breaks=seq(-0.0001,0.0001,length.out = 255)
col=diverge_hcl(256)

set8.6 <- settle_08[[6]]-settle_08[[1]]
set6.6 <- settle_06[[6]]-settle_06[[1]]
set4.6 <- settle_04[[6]]-settle_04[[1]]
set2.6 <- settle_02[[6]]-settle_02[[1]]

plot(set8.6, breaks=breaks, col=col, main="sand_08")
plot(set6.6, breaks=breaks, col=col, main="sand_06")
plot(set4.6, breaks=breaks, col=col, main="sand_04")
plot(set2.6, breaks=breaks, col=col, main="sand_02")


##### EVALUATE FAM - PARTICLE TRAJECTORIES #####

##### Evaluate floats
### it's being stored in a rather weird way...
### if you read in as a raster, the rows are the time-steps
### everything between 153-183 is NA for some reason

#data.flts <- "C:/Users/jjansen/Desktop/science/data_environmental/Circumpolar_ROMS/10km_outputs/"
data.flts <- "C:/Users/jjansen/Desktop/science/data_environmental/Circumpolar_ROMS/10km_outputs/output_sed_float_test1/"
data.flts <- "C:/Users/jjansen/Desktop/science/data_environmental/Circumpolar_ROMS/10km_outputs/output_sed_float_test3/"
data.flts <- "C:/Users/jjansen/Desktop/science/data_environmental/Circumpolar_ROMS/10km_outputs/output_sed_float_test4/"
data.flts <- "C:/Users/jjansen/Desktop/science/data_environmental/Circumpolar_ROMS/10km_outputs/output_sed_float_test5/"
grd.x <- brick(paste0(data.flts,"ocean_flt.nc"), varname="Xgrid", level=1)
grd.y <- brick(paste0(data.flts,"ocean_flt.nc"), varname="Ygrid", level=1)
grd.z <- brick(paste0(data.flts,"ocean_flt.nc"), varname="Zgrid", level=1)
flts.x <- brick(paste0(data.flts,"ocean_flt.nc"), varname="x", level=1)
flts.y <- brick(paste0(data.flts,"ocean_flt.nc"), varname="y", level=1)
flts.z <- brick(paste0(data.flts,"ocean_flt.nc"), varname="depth", level=1)
flts.s1 <- brick(paste0(data.flts,"ocean_flt.nc"), varname="sand_01", level=1)
flts.s2 <- brick(paste0(data.flts,"ocean_flt.nc"), varname="sand_02", level=1)
flts.s3 <- brick(paste0(data.flts,"ocean_flt.nc"), varname="sand_03", level=1)
flts.s4 <- brick(paste0(data.flts,"ocean_flt.nc"), varname="sand_04", level=1)
flts.s5 <- brick(paste0(data.flts,"ocean_flt.nc"), varname="sand_05", level=1)
flts.s6 <- brick(paste0(data.flts,"ocean_flt.nc"), varname="sand_06", level=1)
flts.s7 <- brick(paste0(data.flts,"ocean_flt.nc"), varname="sand_07", level=1)
flts.s8 <- brick(paste0(data.flts,"ocean_flt.nc"), varname="sand_08", level=1)

#### start with a simple 2D plot
## only use a subset of 500 traces
set.seed(1)
s <- sample(1:length(which(!is.na(grd.x[151,]))), 500)

## plot depth vs index
par(mfrow=c(1,1))
plot(flts.z[1:152,s])

## trajectories:
plot(uv$X662904000, xlim=c(200,350), ylim=c(50,200), col=rev(magma(99)))
points(cbind(grd.x[152,s], grd.y[152,s]), pch=16, col="blue")
points(cbind(grd.x[1:151,s], grd.y[1:151,s]), pch=16, col="blue", cex=0.2)
## end-points
plot(uv$X662904000, xlim=c(200,350), ylim=c(50,200), col=rev(magma(99)))
points(cbind(grd.x[1,], grd.y[1,]), pch=16, col="blue", cex=0.2)

## Mertz Glacier end-points
# plot(uv$X662904000, xlim=c(430,480), ylim=c(70,110), col=rev(magma(99)))
plot(h$bathymetry.at.RHO.points, xlim=c(430,480), ylim=c(70,110), col=rev(magma(99)))
points(cbind(grd.x[1,], grd.y[1,]), pch=16, col="blue", cex=0.2)

s2 <- sample(1:ncol(flts.z),9)
## plot a single float through time:
par(mfrow=c(3,3))
plot(flts.z[1:152,s2[1]])
plot(flts.z[1:152,s2[2]])
plot(flts.z[1:152,s2[3]])
plot(flts.z[1:152,s2[4]])
plot(flts.z[1:152,s2[5]])
plot(flts.z[1:152,s2[6]])
plot(flts.z[1:152,s2[7]])
plot(flts.z[1:152,s2[8]])
plot(flts.z[1:152,s2[9]])




#### Need to translate grid points into polar stereographic projection values to allow proper plotting:
r2 <- raster(paste0(env.derived,"Circumpolar_EnvData_500m_shelf_bathy_gebco_depth.grd"))
## ROMS lon/lats
lon.ra <- raster(paste0(env.raw,"waom10extend_grd.nc"), varname="lon_rho")
lat.ra <- raster(paste0(env.raw,"waom10extend_grd.nc"), varname="lat_rho")
## lons/lats of the floats at each time-interval
lons.1 <- extract(lon.ra, cbind(grd.x[1,],grd.y[1,]))
lons.not.nas <- which(!is.na(lons.1))
pts.list <- list()
for(i in 1:152){
  print(i)
  lons <- extract(lon.ra, cbind(grd.x[i,],grd.y[i,]))
  lats <- extract(lat.ra, cbind(grd.x[i,],grd.y[i,]))
  ## spatialpoints:
  sp <- SpatialPoints(coords=cbind(lons[lons.not.nas],lats[lons.not.nas]),
                      proj4string=sp::CRS("+proj=longlat +datum=WGS84"))
  # spdf <- SpatialPointsDataFrame(sp, 
  pts.list[[i]] <- spTransform(sp, crs(r2))
}


### start with plotting in 2D:
##
plot(pts.list[[1]][s], pch=16, cex=0.5, col=viridis(500))
points(list.rbind(lapply(pts.list,"[",s)), pch=16, cex=0.2, col=viridis(500))
plot(coast.proj, add=TRUE)
##
plot(pts.list[[1]][s], pch=16, cex=0.5, col=rep(viridis(500)[1],500))
points(list.rbind(lapply(pts.list,"[",s)), pch=16, cex=0.2, col=viridis(75500))
plot(coast.proj, add=TRUE)
# points(pts.list[1:152][s])
# points(cbind(grd.x[1:151,s], grd.y[1:151,s]), pch=16, col=rep(viridis(100),10), cex=0.5)
# points(cbind(grd.x[152,s], grd.y[152,s]), pch=16, col=rep(viridis(100),10))

plot(uv$X662904000, xlim=c())
plot(pts.list[[1]][s], pch=16, cex=0.5, col=viridis(500))
points(list.rbind(lapply(pts.list,"[",s)), pch=16, cex=0.2, col=viridis(500))
plot(coast.proj, add=TRUE)



pts.subset.list <- list()
for(i in 1:152){
  pts.subset.list[[i]] <- pts.list[[i]]@coords[s,]
}
## traces of individual particles
pts.subset.df <- list.rbind(pts.subset.list)
traces.list <- list()
for(i in 1:500){
  traces.list[[i]] <- pts.subset.df[seq(i,nrow(pts.subset.df),by=500),]
}
traces.df <- list.rbind(traces.list)



# plot(cbind(grd.x[152,], grd.y[152,]), pch=16, cex=0.2, xlim=c(220,370), ylim=c(90,190))
# points(cbind(grd.x[1:151,109875], grd.y[1:151,109875]), pch=16, col="green3", cex=0.5)
# points(cbind(grd.x[1:151,109900], grd.y[1:151,109900]), pch=16, col="salmon3", cex=0.5)
# points(cbind(grd.x[1:151,109925], grd.y[1:151,109925]), pch=16, col="cyan3", cex=0.5)
# points(cbind(grd.x[1:151,109950], grd.y[1:151,109950]), pch=16, col="plum3", cex=0.5)
# points(cbind(grd.x[1:151,109975], grd.y[1:151,109975]), pch=16, col="red3", cex=0.5)
# points(cbind(grd.x[1:151,109995], grd.y[1:151,109995]), pch=16, col="purple3", cex=0.5)
# points(cbind(grd.x[1:151,112660], grd.y[1:151,112660]), pch=16, col="green3", cex=0.5)
# points(cbind(grd.x[1:151,112685], grd.y[1:151,112785]), pch=16, col="salmon3", cex=0.5)
# points(cbind(grd.x[1:151,112710], grd.y[1:151,112710]), pch=16, col="cyan3", cex=0.5)
# points(cbind(grd.x[1:151,112735], grd.y[1:151,112735]), pch=16, col="plum3", cex=0.5)
# points(cbind(grd.x[1:151,112760], grd.y[1:151,112760]), pch=16, col="red3", cex=0.5)
# points(cbind(grd.x[1:151,112785], grd.y[1:151,112785]), pch=16, col="purple3", cex=0.5)

##### animate 2D points
library(SOmap)
library(gganimate) ## remotes::install_github("thomasp85/gganimate")
require(transformr) ## remotes::install_github("thomasp85/transformr")
require(sf)
require(rgdal)
require(rgeos)

# ## use a rotated polar projection
# polar_proj_rot <- "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=90 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
# ele[c("x", "y")] <- rgdal::project(as.matrix(ele[c("lon", "lat")]), polar_proj_rot)
# 
# ## use elevation data from ETOPO2 to show the land masses
# ## read it and reproject from long-lat to our polar projection
# cx <- projectRaster(readtopo("etopo2", xylim = c(-180, 180, -90, -40)),
#                     raster(extent(c(-4e6, 4e6, -5e6, 5e6)), nrows = 400, ncols = 400, crs = polar_proj_rot))
# cx <- as.data.frame(cx, xy = TRUE)
# cx$z[cx$z <= 0] <- NA_real_

## create a lagged version of the track that we can use to show a trailing "worm"
tail_length <- 8
pts_lagged <- ele %>% mutate(lag_n = 0)
for (li in seq_len(tail_length)) ele_lagged <- rbind(ele_lagged, ele %>% mutate(date = lead(date, li), lag_n = -li))
ele_lagged <- ele_lagged %>% dplyr::filter(!is.na(date))

g <- ggplot() + 
  geom_raster(data = cx, aes(x, y, fill = z)) +
  scale_fill_distiller(palette = "Greys", guide = FALSE, na.value = "#FFFFFF00") +
  geom_path(data = ele_lagged, aes(x, y, alpha = lag_n), colour = "orange", size = 1) +
  scale_alpha_continuous(guide = FALSE) +
  geom_path(data = ele, aes(x, y), colour = "orange", size = 2) +
  theme_void() +
  # xlim(c(-4e6, 4e6)) +
  # ylim(c(-4.75e6, 4.75e6)) +
  transition_time(date)
animate(g)




























s1.152 <- cbind(grd.x[152,], grd.y[152,], flts.s1[152,])
s8.152 <- cbind(grd.x[152,], grd.y[152,], flts.s8[152,])

##
sink.1 <- cbind(grd.x[1,], grd.y[1,], flts.z[1,])
sink.31 <- cbind(grd.x[31,], grd.y[31,], flts.z[31,])
sink.51 <- cbind(grd.x[51,], grd.y[51,], flts.z[51,])
sink.71 <- cbind(grd.x[71,], grd.y[71,], flts.z[71,])
sink.91 <- cbind(grd.x[91,], grd.y[91,], flts.z[91,])
sink.111 <- cbind(grd.x[111,], grd.y[111,], flts.z[111,])
sink.121 <- cbind(grd.x[121,], grd.y[121,], flts.z[121,])
sink.131 <- cbind(grd.x[131,], grd.y[131,], flts.z[131,])
sink.100 <- cbind(grd.x[100,], grd.y[100,], flts.z[100,])
sink.152 <- cbind(grd.x[152,], grd.y[152,], flts.z[152,])
sink.182 <- cbind(grd.x[182,], grd.y[182,], flts.z[182,])

par(mfrow=c(3,3))
plot(sink.1[,c(1,3)], cex=0.5, ylim=c(-900,0))
plot(sink.31[,c(1,3)], cex=0.5, ylim=c(-900,0))
plot(sink.51[,c(1,3)], cex=0.5, ylim=c(-900,0))
plot(sink.71[,c(1,3)], cex=0.5, ylim=c(-900,0))
plot(sink.91[,c(1,3)], cex=0.5, ylim=c(-900,0))
plot(sink.111[,c(1,3)], cex=0.5, ylim=c(-900,0))
plot(sink.131[,c(1,3)], cex=0.5, ylim=c(-900,0))
plot(sink.152[,c(1,3)], cex=0.5, ylim=c(-900,0))
plot(sink.182[,c(1,3)], cex=0.5, ylim=c(-900,0))

## plot a single float through time:
par(mfrow=c(3,3))
plot(flts.z[1:152,100001])
plot(flts.z[1:152,100010])
plot(flts.z[1:152,100011])
plot(flts.z[1:152,100100])
plot(flts.z[1:152,100101])
plot(flts.z[1:152,100110])
plot(flts.z[1:152,100111])
plot(flts.z[1:152,101000])


## What else could we do??? 
nc_data <- nc_open(paste0(data.flts,"ocean_flt.nc"))
nc.depth <- ncvar_get(nc_data, "depth")





## tried previously:

pts.1 <- cbind(grd.x[1,], grd.y[1,], grd.z[1,])
pts.31 <- cbind(grd.x[31,], grd.y[31,], grd.z[31,])
pts.51 <- cbind(grd.x[51,], grd.y[51,], grd.z[51,])
pts.71 <- cbind(grd.x[71,], grd.y[71,], grd.z[71,])
pts.91 <- cbind(grd.x[91,], grd.y[91,], grd.z[91,])
pts.111 <- cbind(grd.x[111,], grd.y[111,], grd.z[111,])
pts.121 <- cbind(grd.x[121,], grd.y[121,], grd.z[121,])
pts.131 <- cbind(grd.x[131,], grd.y[131,], grd.z[131,])
pts.100 <- cbind(grd.x[100,], grd.y[100,], grd.z[100,])
pts.152 <- cbind(grd.x[152,], grd.y[152,], grd.z[152,])


plot(pts.1[,1:2], cex=0.5)
points(pts.31[,1:2], col="red", cex=0.5)
points(pts.100[,1:2], col="blue", cex=0.3)
points(pts.152[,1:2], col="green", cex=0.3)

library(rgl)

plot3d(pts.1[10000:20000,])
points3d(pts.31[10000:20000,], col="red")
points3d(pts.152[10000:20000,], col="red")

flts.1   <- cbind(flts.x[1,],   flts.y[1,],   flts.z[1,])
flts.31  <- cbind(flts.x[31,],  flts.y[31,],  flts.z[31,])
flts.51  <- cbind(flts.x[51,],  flts.y[51,],  flts.z[51,])
flts.71  <- cbind(flts.x[71,],  flts.y[71,],  flts.z[71,])
flts.90  <- cbind(flts.x[90,],  flts.y[90,],  flts.z[90,])
flts.91  <- cbind(flts.x[91,],  flts.y[91,],  flts.z[91,])
flts.92  <- cbind(flts.x[92,],  flts.y[92,],  flts.z[92,])
flts.111 <- cbind(flts.x[111,], flts.y[111,], flts.z[111,])
flts.121 <- cbind(flts.x[121,], flts.y[121,], flts.z[121,])
flts.131 <- cbind(flts.x[131,], flts.y[131,], flts.z[131,])
flts.100 <- cbind(flts.x[100,], flts.y[100,], flts.z[100,])
flts.152 <- cbind(flts.x[152,], flts.y[152,], flts.z[152,])


plot(flts.1[25000:35000,3])

plot(flts.91[25000:35000,3])
points(flts.90[25000:35000,3],col="blue")
points(flts.92[25000:35000,3],col="red")


plot(flts.91[25385,3], ylim=c(-100,0))
points(flts.90[25385,3], col="red")
points(flts.92[25385,3], col="blue")
points(flts.111[25385,3])
points(flts.131[25385,3])
points(flts.152[25385,3])
points(flts.71[25385,3])
points(flts.51[25385,3])
points(flts.31[25385,3])


plot(pts.1)
