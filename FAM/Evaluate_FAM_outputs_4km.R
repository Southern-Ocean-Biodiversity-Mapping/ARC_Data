###########################################
########## Evaluate ROMS outputs ##########
###########################################

########## Setup environment and prepare data ##########

#### setup environment
library(raster)
library(colorspace)
library(viridis)
library(ncdf4)
library(raadtools)
## data-directories
env.dir <- "C:/Users/jjansen/Desktop/science/data_environmental/"
env.raw <- "E:/science/data_environmental/raw/"
env.derived <- paste0(env.dir,"derived/")
AAD_dir <- paste0(env.dir,"raw/accessed_through_R")
ARC_data_dir <- "C:/Users/jjansen/Desktop/science/SouthernOceanBiodiversityMapping/ARC_Data/"
ant.proj <- "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

#### load non-ROMS data
## load projected coastline for plotting
load(paste0(env.derived,"Circumpolar_Coastline.Rdata"))
r2 <- raster(paste0(env.derived,"Circumpolar_EnvData_500m_shelf_bathy_gebco_depth.grd"))
## load circumpolar Diatom data
load(paste0(ARC_data_dir,"FAM/Circumpolar_Diatom_metadata.Rdata"))
## load Mertz Diatom data
mertz.diatom.dat <- read.csv("C:/Users/jjansen/Desktop/science/PhD/1 bentho-pelagic-coupling/DiatomAbundance.csv")
mertz.diatom.locs <- project(as.matrix(mertz.diatom.dat[,c(4,3)]),proj=ant.proj)

#### load ROMS-runs:
data.dat100 <- paste0(env.dir,"Circumpolar_ROMS/4km_outputs/output_sed_test1/")
# data.dat100 <- paste0(env.dir,"Circumpolar_ROMS/10km_outputs/output_sed_float_test3/")
# data.dat200 <- paste0(env.dir,"Circumpolar_ROMS/10km_outputs/output_sed_float_test4/")
# data.dat400 <- paste0(env.dir,"Circumpolar_ROMS/10km_outputs/output_sed_float_test5/")
# data.dat <- paste0("Circumpolar_ROMS/10km_outputs/sed_test1/")
# data.dat <- paste0(env.dir,"Circumpolar_ROMS/10km_outputs/sed_test5/")

#### load lon/lat information from ROMS-grid
grd4k_nc <- nc_open(paste0(env.raw,"waom4extend_grd.nc"))
lon_rho <- ncvar_get(grd4k_nc, varid="lon_rho")
lat_rho <- ncvar_get(grd4k_nc, varid="lat_rho")

#### Prepare empty rasters to assign correct projected values to
roms.coords.proj <- project(cbind(c(lon_rho), c(lat_rho)), proj=ant.proj)
x.range <- c(min(roms.coords.proj[,1])-2000,max(roms.coords.proj[,1])+2000)
y.range <- c(min(roms.coords.proj[,2])-2000,max(roms.coords.proj[,2])+2000)
empty.roms.ra <- raster(extent(c(x.range,y.range)), crs=ant.proj, resolution=4000)

#### load other ROMS data
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

# #temporal mean seafloor current speed
# mean.u <- sum(u)
# mean.v <- sum(v)
# #absolute mean seafloor current speed
# abs.u <- sum(abs(u))
# abs.v <- sum(abs(v))
# #residual seafloor current speed
# res.u <- abs.u-mean.u
# res.v <- abs.v-mean.v
# mean.uv2 <- sqrt(mean.u^2+mean.v^2)
# abs.uv2 <- sqrt(abs.u^2+abs.v^2)
# res.uv2 <- abs.uv2-mean.uv2

########## Overview plots - Ocean currents ##########
## ocean current data
#cols <- rev(terrain.colors(96))
cols <- rev(magma(96)) #12, 24, 48, 96
breaks <- seq(0,6,length.out=97)
par(mfrow=c(2,2))
plot(mean.uv, breaks=breaks, col=cols, main="seafloor currents - mean")
plot(abs.uv, breaks=breaks, col=cols, main="seafloor currents - absolute speed")
plot(res.uv, breaks=breaks, col=cols, main="seafloor currents - residual")
plot(uv_31, breaks=breaks, col=cols, main="seasurface currents")
## ocean current data in the Mertz region
xlim <- c(1300000,1900000)
ylim <- c(-2300000,-1900000)
par(mfrow=c(2,2))
plot(mean.uv, breaks=breaks, col=cols, main="seafloor currents - mean", xlim=xlim, ylim=ylim)
#plot(coast.proj, add=TRUE)
plot(abs.uv, breaks=breaks, col=cols, main="seafloor currents - absolute speed", xlim=xlim, ylim=ylim)
#plot(coast.proj, add=TRUE)
plot(res.uv, breaks=breaks, col=cols, main="seafloor currents - residual", xlim=xlim, ylim=ylim)
#plot(coast.proj, add=TRUE)
plot(h2, col=cols, main="depth", xlim=xlim, ylim=ylim)
#plot(coast.proj, add=TRUE)

########## EVALUATE FAM ##########
## values originally in kg/m2/s -> change to kg/m2/day (*86400)
### NOTE THAT CURRENTLY DIFFERENT MODEL RUNS HAVE ALL THE SAME OUTPUT APART FROM THE "ocean_flt.nc" file
surf_01 <- surf_02 <- surf_03 <- surf_04 <- surf_05 <- surf_06 <- surf_07 <- surf_08 <- brick(empty.roms.ra,nl=6)
susp_01 <- susp_02 <- susp_03 <- susp_04 <- susp_05 <- susp_06 <- susp_07 <- susp_08 <- brick(empty.roms.ra,nl=6)
settle_01 <- settle_02 <- settle_03 <- settle_04 <- settle_05 <- settle_06 <- settle_07 <- settle_08 <- brick(empty.roms.ra,nl=7)
## surface production
surf_01[] <- brick(paste0(data.dat100,"ocean_avg_0001.nc"), varname="sand_01", level=31)[]*86400
surf_02[] <- brick(paste0(data.dat100,"ocean_avg_0001.nc"), varname="sand_02", level=31)[]*86400
surf_03[] <- brick(paste0(data.dat100,"ocean_avg_0001.nc"), varname="sand_03", level=31)[]*86400
surf_04[] <- brick(paste0(data.dat100,"ocean_avg_0001.nc"), varname="sand_04", level=31)[]*86400
surf_05[] <- brick(paste0(data.dat100,"ocean_avg_0001.nc"), varname="sand_05", level=31)[]*86400
surf_06[] <- brick(paste0(data.dat100,"ocean_avg_0001.nc"), varname="sand_06", level=31)[]*86400
surf_07[] <- brick(paste0(data.dat100,"ocean_avg_0001.nc"), varname="sand_07", level=31)[]*86400
surf_08[] <- brick(paste0(data.dat100,"ocean_avg_0001.nc"), varname="sand_08", level=31)[]*86400
## bottom ocean layer
susp_01[] <- brick(paste0(data.dat100,"ocean_avg_0001.nc"), varname="sand_01", level=1)[]*86400
susp_02[] <- brick(paste0(data.dat100,"ocean_avg_0001.nc"), varname="sand_02", level=1)[]*86400
susp_03[] <- brick(paste0(data.dat100,"ocean_avg_0001.nc"), varname="sand_03", level=1)[]*86400
susp_04[] <- brick(paste0(data.dat100,"ocean_avg_0001.nc"), varname="sand_04", level=1)[]*86400
susp_05[] <- brick(paste0(data.dat100,"ocean_avg_0001.nc"), varname="sand_05", level=1)[]*86400
susp_06[] <- brick(paste0(data.dat100,"ocean_avg_0001.nc"), varname="sand_06", level=1)[]*86400
susp_07[] <- brick(paste0(data.dat100,"ocean_avg_0001.nc"), varname="sand_07", level=1)[]*86400
susp_08[] <- brick(paste0(data.dat100,"ocean_avg_0001.nc"), varname="sand_08", level=1)[]*86400
## settled particles
settle_01[] <- brick(paste0(data.dat100,"ocean_his_0001.nc"), varname="sandfrac_01", level=1)[]*86400
settle_02[] <- brick(paste0(data.dat100,"ocean_his_0001.nc"), varname="sandfrac_02", level=1)[]*86400
settle_03[] <- brick(paste0(data.dat100,"ocean_his_0001.nc"), varname="sandfrac_03", level=1)[]*86400
settle_04[] <- brick(paste0(data.dat100,"ocean_his_0001.nc"), varname="sandfrac_04", level=1)[]*86400
settle_05[] <- brick(paste0(data.dat100,"ocean_his_0001.nc"), varname="sandfrac_05", level=1)[]*86400
settle_06[] <- brick(paste0(data.dat100,"ocean_his_0001.nc"), varname="sandfrac_06", level=1)[]*86400
settle_07[] <- brick(paste0(data.dat100,"ocean_his_0001.nc"), varname="sandfrac_07", level=1)[]*86400
settle_08[] <- brick(paste0(data.dat100,"ocean_his_0001.nc"), varname="sandfrac_08", level=1)[]*86400

########## Overview plots - FAM ##########
## surface productivity 
# plot(surf_08*86400, breaks=seq(0,0.24,length.out=25), col=rev(heat.colors(25)))
# #plot(susp_08*86400, breaks=seq(0,0.24,length.out=25), col=rev(heat.colors(25)))
# plot(susp_08*86400, breaks=round(seq(-0.1,0.24,length.out=35),3), col=c(rev(blues9),rev(heat.colors(25))))
# plot(settle_08*86400, breaks=c(10780,seq(10795,10810,length.out=25),10935), col=c("green",rev(heat.colors(24)),"black"))

#### comparisons
breaks1 <- seq(0,0.25,length.out=51)
breaks2 <- c(-0.5,round(seq(0,0.25,length.out=51),3))
#breaks3 <- c(10780,seq(10795,10810,length.out=25),10935)
breaks3 <- c(10795,seq(10800,10805,length.out=51))
cols1 <- rev(viridis(50))#rev(magma(48))
cols2 <- c("grey",rev(viridis(50)))
cols3 <- c("grey",rev(viridis(50)))
par(mfrow=c(3,3), mar=c(3,4,2,1))
plot(surf_08[[2]], breaks=breaks1, col=cols1, main=paste0("surf_08 ",names(surf_08)[2]))
plot(surf_08[[3]], breaks=breaks1, col=cols1, main=paste0("surf_08 ",names(surf_08)[3]))
plot(surf_08[[4]], breaks=breaks1, col=cols1, main=paste0("surf_08 ",names(surf_08)[4]))
plot(susp_08[[2]], breaks=breaks2, col=cols2, main=paste0("susp_08 ",names(susp_08)[2]))
plot(susp_08[[3]], breaks=breaks2, col=cols2, main=paste0("susp_08 ",names(susp_08)[3]))
plot(susp_08[[4]], breaks=breaks2, col=cols2, main=paste0("susp_08 ",names(susp_08)[4]))
plot(settle_08[[2]], breaks=breaks3, col=cols3, main=paste0("settle_08 ",names(settle_08)[2]))
plot(settle_08[[3]], breaks=breaks3, col=cols3, main=paste0("settle_08 ",names(settle_08)[3]))
plot(settle_08[[4]], breaks=breaks3, col=cols3, main=paste0("settle_08 ",names(settle_08)[4]))

#### comparisons Mertz
par(mfrow=c(3,3))
plot(surf_08[[2]], xlim=xlim, ylim=ylim, breaks=breaks1, col=cols1, main=paste0("surf_08 ",names(surf_08)[2]))
plot(surf_08[[3]], xlim=xlim, ylim=ylim, breaks=breaks1, col=cols1, main=paste0("surf_08 ",names(surf_08)[3]))
plot(surf_08[[4]], xlim=xlim, ylim=ylim, breaks=breaks1, col=cols1, main=paste0("surf_08 ",names(surf_08)[4]))
plot(susp_08[[2]], xlim=xlim, ylim=ylim, breaks=breaks2, col=cols2, main=paste0("susp_08 ",names(susp_08)[2]))
plot(susp_08[[3]], xlim=xlim, ylim=ylim, breaks=breaks2, col=cols2, main=paste0("susp_08 ",names(susp_08)[3]))
plot(susp_08[[4]], xlim=xlim, ylim=ylim, breaks=breaks2, col=cols2, main=paste0("susp_08 ",names(susp_08)[4]))
plot(settle_08[[2]], xlim=xlim, ylim=ylim, breaks=breaks3, col=cols3, main=paste0("settle_08 ",names(settle_08)[2]))
plot(settle_08[[3]], xlim=xlim, ylim=ylim, breaks=breaks3, col=cols3, main=paste0("settle_08 ",names(settle_08)[3]))
plot(settle_08[[4]], xlim=xlim, ylim=ylim, breaks=breaks3, col=cols3, main=paste0("settle_08 ",names(settle_08)[4]))

#### comparisons surf / sink / settle for Mertz
par(mfrow=c(3,3))
plot(settle_04[[2]], xlim=xlim, ylim=ylim, breaks=breaks3, col=cols1, main=paste0("settle_04 ",names(settle_04)[2]))
plot(settle_04[[4]], xlim=xlim, ylim=ylim, breaks=breaks3, col=cols1, main=paste0("settle_04 ",names(settle_04)[4]))
plot(settle_04[[6]], xlim=xlim, ylim=ylim, breaks=breaks3, col=cols1, main=paste0("settle_04 ",names(settle_04)[6]))
plot(settle_06[[2]], xlim=xlim, ylim=ylim, breaks=breaks3, col=cols1, main=paste0("settle_06 ",names(settle_06)[2]))
plot(settle_06[[4]], xlim=xlim, ylim=ylim, breaks=breaks3, col=cols1, main=paste0("settle_06 ",names(settle_06)[4]))
plot(settle_06[[6]], xlim=xlim, ylim=ylim, breaks=breaks3, col=cols1, main=paste0("settle_06 ",names(settle_06)[6]))
plot(settle_08[[2]], xlim=xlim, ylim=ylim, breaks=breaks3, col=cols1, main=paste0("settle_08 ",names(settle_08)[2]))
plot(settle_08[[4]], xlim=xlim, ylim=ylim, breaks=breaks3, col=cols1, main=paste0("settle_08 ",names(settle_08)[4]))
plot(settle_08[[6]], xlim=xlim, ylim=ylim, breaks=breaks3, col=cols1, main=paste0("settle_08 ",names(settle_08)[6]))

########## Comparison to Diatom values ##########
#### Extract flux values at Diatom locations
## circumpolar
surf.vals <- extract(surf_08, polar.dat.diatom)
susp.vals <- extract(susp_08, polar.dat.diatom)
settle.vals <- extract(settle_08, polar.dat.diatom)
diatom.vals <- diatom.dat$Absolute_diatom_abundance_in_million_valves_per_gram_dry_sediment
## Mertz region
surf.vals.m <- extract(surf_08,mertz.diatom.locs)
susp.vals.m <- extract(susp_08,mertz.diatom.locs)
settle.vals.m <- extract(settle_08,mertz.diatom.locs)
diatom.vals.m <- mertz.diatom.dat$Abundance_NO_benthics__valves_g_dry_weight_

#### circumpolar relationship
par(mfrow=c(2,2))
plot(surf.vals[,5], log(diatom.vals))
plot(susp.vals[,5], log(diatom.vals))
plot(settle.vals[,5], log(diatom.vals))
#### compare to Mertz data:
par(mfrow=c(2,2))
plot(surf.vals.m[,5], log(diatom.vals.m))
plot(susp.vals.m[,5], log(diatom.vals.m))
plot(settle.vals.m[,5], log(diatom.vals.m))

# library(MASS)
# library(modEvA)
# summary(glm.nb(diatom.vals.m~surf.vals.m))
# summary(glm.nb(diatom.vals.m~susp.vals.m))
# summary(glm.nb(diatom.vals.m~settle.vals.m))
# Dsquared(glm.nb(diatom.vals~surf.vals))
# Dsquared(glm.nb(diatom.vals~susp.vals))
# Dsquared(glm.nb(diatom.vals~settle.vals))

## sort diatom cruises by region
dc.ap.w <- c("DF82","DF86","LMG1311","NBP1001","PD88-III","PD90-7")
dc.ap.e <- c("LMG0502","NBP0003","NBP0603","NBP1203")
dc.ap.both <- c("NBP0107")
dc.eant <- c("NBP1402", "NBP0101")
dc.ross <- c("DF83 III", "NBP9501", "DF83 II", "NBP9401", "PD92")
dc.amund <- "NBP0702"

dat <- data.frame("surf.vals"=surf.vals[,5], "diatom.vals"=diatom.vals)
dat$cruise <- diatom.dat$Cruise
dat$cruise.short <- strtrim(diatom.dat$Cruise,3)
dat$region <- "empty"
dat[diatom.dat$Cruise%in%dc.ap.w,5] <- "dc.ap.w"
dat[diatom.dat$Cruise%in%dc.ap.e,5] <- "dc.ap.e"
dat[diatom.dat$Cruise%in%dc.ap.both,5] <- "dc.ap.both"
dat[diatom.dat$Cruise%in%dc.eant,5] <- "dc.eant"
dat[diatom.dat$Cruise%in%dc.ross,5] <- "dc.ross"
dat[diatom.dat$Cruise%in%dc.amund,5] <- "dc.amund"

library(ggplot2)
# ggplot(dat, aes(x=surf.vals, y=log(diatom.vals), color=cruise)) +
#   geom_point() 
# ggplot(dat, aes(x=surf.vals, y=log(diatom.vals), color=cruise.short)) +
#   geom_point() 
# ggplot(dat, aes(x=surf.vals, y=log(diatom.vals), color=region)) +
#   geom_point() 
ggplot(dat, aes(x=surf.vals, y=log(diatom.vals), color=region)) +
  geom_point() +
  facet_wrap(~region)
ggplot(dat, aes(x=surf.vals, y=log(diatom.vals), color=region)) +
  geom_point() +
  facet_wrap(~cruise)


####################################################################



## particle-tracks
#it's being stored in a rather weird way...
#if you read in as a raster, the rows are the time-steps
#everything between 153-183 is NA for some reason
flts.z100 <- brick(paste0(data.dat100,"ocean_flt.nc"), varname="depth", level=1)
grd.x100 <- brick(paste0(data.dat100,"ocean_flt.nc"), varname="Xgrid", level=1)
grd.y100 <- brick(paste0(data.dat100,"ocean_flt.nc"), varname="Ygrid", level=1)
# flts.z200 <- brick(paste0(data.dat200,"ocean_flt.nc"), varname="depth", level=1)
# grd.x200 <- brick(paste0(data.dat200,"ocean_flt.nc"), varname="Xgrid", level=1)
# grd.y200 <- brick(paste0(data.dat200,"ocean_flt.nc"), varname="Ygrid", level=1)
# flts.z400 <- brick(paste0(data.dat400,"ocean_flt.nc"), varname="depth", level=1)
# grd.x400 <- brick(paste0(data.dat400,"ocean_flt.nc"), varname="Xgrid", level=1)
# grd.y400 <- brick(paste0(data.dat400,"ocean_flt.nc"), varname="Ygrid", level=1)

##### EVALUATE FAM - Plot results #####
#### Overview for Ross Sea
xlim=c(250,340)
ylim=c(100,180)
breaks <- seq(-0.000025,0.000025,length.out=99)
breaks.susp <- seq(0,0.000025, length.out=256)
breaks.settle <- seq(0.1245,0.1255, length.out=256)
#ticks.susp <- c(0,0.00001,0.00002)

### NPP, currents and FAM
par(mfrow=c(2,2))
plot(surf_01[[1]], xlim=xlim, ylim=ylim, main="NPP-avg", breaks=breaks.susp)
contour(h, add=TRUE)
plot(susp_01[[6]], xlim=xlim, ylim=ylim, main="FAM-avg-suspended", breaks=breaks.susp)
contour(h, add=TRUE)
plot(settle_01[[6]], xlim=xlim, ylim=ylim, main="FAM-his-settled", breaks=breaks.settle)
contour(h, add=TRUE)
plot(uv[[6]], xlim=xlim, ylim=ylim, main="currents_avg")
contour(h, add=TRUE)

### particle trajectories (only use a subset of 500 traces)
set.seed(1)
s <- sample(1:length(which(!is.na(grd.x100[151,]))), 500)

## trajectories of 500 floats
par(mfrow=c(2,3), mar=c(3,0,4,0), oma=c(0,2,0,3))
plot(uv[[6]], xlim=c(200,350), ylim=c(50,200), col=rev(magma(99)), main="trajectories - 100m/day")
points(cbind(grd.x100[152,s], grd.y100[152,s]), pch=16, col="blue")
points(cbind(grd.x100[1:151,s], grd.y100[1:151,s]), pch=16, col="blue", cex=0.2)
plot(uv[[6]], xlim=c(200,350), ylim=c(50,200), col=rev(magma(99)), main="trajectories - 200m/day")
points(cbind(grd.x200[152,s], grd.y200[152,s]), pch=16, col="blue")
points(cbind(grd.x200[1:151,s], grd.y200[1:151,s]), pch=16, col="blue", cex=0.2)
plot(uv[[6]], xlim=c(200,350), ylim=c(50,200), col=rev(magma(99)), main="trajectories - 400m/day")
points(cbind(grd.x400[152,s], grd.y400[152,s]), pch=16, col="blue")
points(cbind(grd.x400[1:151,s], grd.y400[1:151,s]), pch=16, col="blue", cex=0.2)

## end-points of all floats after 2 months
plot(uv[[6]], xlim=c(200,350), ylim=c(50,200), col=rev(magma(99)), main="end-points - 100m/day")
points(cbind(grd.x100[121,], grd.y100[121,]), pch=16, col="blue", cex=0.2)
plot(uv[[6]], xlim=c(200,350), ylim=c(50,200), col=rev(magma(99)), main="end-points - 200m/day")
points(cbind(grd.x200[121,], grd.y200[121,]), pch=16, col="blue", cex=0.2)
plot(uv[[6]], xlim=c(200,350), ylim=c(50,200), col=rev(magma(99)), main="end-points - 400m/day")
points(cbind(grd.x400[121,], grd.y400[121,]), pch=16, col="blue", cex=0.2)

# ## end-points of all floats after 6 months
# plot(uv[[6]], xlim=c(200,350), ylim=c(50,200), col=rev(magma(99)), main="end-points - 100m/day")
# points(cbind(grd.x100[1,], grd.y100[1,]), pch=16, col="blue", cex=0.2)
# plot(uv[[6]], xlim=c(200,350), ylim=c(50,200), col=rev(magma(99)), main="end-points - 200m/day")
# points(cbind(grd.x200[1,], grd.y200[1,]), pch=16, col="blue", cex=0.2)
# plot(uv[[6]], xlim=c(200,350), ylim=c(50,200), col=rev(magma(99)), main="end-points - 400m/day")
# points(cbind(grd.x400[1,], grd.y400[1,]), pch=16, col="blue", cex=0.2)
# 
# ## end-points of all floats after 1 month (essentially the starting point)
# plot(uv[[6]], xlim=c(200,350), ylim=c(50,200), col=rev(magma(99)), main="end-points - 100m/day")
# points(cbind(grd.x100[152,], grd.y100[152,]), pch=16, col="blue", cex=0.2)
# plot(uv[[6]], xlim=c(200,350), ylim=c(50,200), col=rev(magma(99)), main="end-points - 200m/day")
# points(cbind(grd.x200[152,], grd.y200[152,]), pch=16, col="blue", cex=0.2)
# plot(uv[[6]], xlim=c(200,350), ylim=c(50,200), col=rev(magma(99)), main="end-points - 400m/day")
# points(cbind(grd.x400[152,], grd.y400[152,]), pch=16, col="blue", cex=0.2)


#### Particles move indefinitely... Stop when they reach the seafloor and have low current speeds
## speed of particles

## plot a single float through time:
par(mfrow=c(3,3))
plot(flts.z100[1:152,s[1]])
plot(flts.z100[1:152,s[2]])
plot(flts.z100[1:152,s[3]])
plot(flts.z100[1:152,s[4]])
plot(flts.z100[1:152,s[5]])
plot(flts.z100[1:152,s[6]])
plot(flts.z100[1:152,s[7]])
plot(flts.z100[1:152,s[8]])
plot(flts.z100[1:152,s[9]])

## change particle locations to projected values

## calculate movement speed using changes in m in grid.x and grid.y per day

## plot depth against movement speed


#### 3D plot of particle subset












#### start with a simple 2D plot

## plot depth vs index
par(mfrow=c(1,1))
plot(flts.z[1:152,s])



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
