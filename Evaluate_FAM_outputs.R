## Evaluate ROMS outputs
library(raster)
library(colorspace)
data.dat <- "C:/Users/jjansen/Desktop/science/data_environmental/Circumpolar_ROMS/10km_outputs/sed_test1/"
data.dat <- "C:/Users/jjansen/Desktop/science/data_environmental/Circumpolar_ROMS/10km_outputs/sed_test5/"

xlim=c(250,340)
ylim=c(100,180)

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


##### Evaluate something?
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

par(mfrow=c(2,2))
plot(surf_01[[6]], xlim=xlim, ylim=ylim, main="NPP", breaks=breaks.susp)
contour(h, add=TRUE)
plot(susp_01[[6]], xlim=xlim, ylim=ylim, main="FAM-suspended", breaks=breaks.susp)
contour(h, add=TRUE)
plot(settle_01[[6]], xlim=xlim, ylim=ylim, main="FAM-settled", breaks=breaks.settle)
contour(h, add=TRUE)
plot(uv2[[6]], xlim=xlim, ylim=ylim, main="currents_his")
contour(h, add=TRUE)

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
plot(sand_01[[6]]-surf_01[[6]], xlim=xlim, ylim=ylim, breaks=breaks, col=diverge_hsv(99), main="FAM-NPP")
# contour(h, add=TRUE)
plot(sand_02[[6]]-surf_02[[6]], xlim=xlim, ylim=ylim, breaks=breaks, col=diverge_hsv(99), main="FAM-NPP")
# contour(h, add=TRUE)
plot(sand_03[[6]]-surf_03[[6]], xlim=xlim, ylim=ylim, breaks=breaks, col=diverge_hsv(99), main="FAM-NPP")
# contour(h, add=TRUE)
plot(sand_04[[6]]-surf_04[[6]], xlim=xlim, ylim=ylim, breaks=breaks, col=diverge_hsv(99), main="FAM-NPP")
# contour(h, add=TRUE)
plot(sand_05[[6]]-surf_05[[6]], xlim=xlim, ylim=ylim, breaks=breaks, col=diverge_hsv(99), main="FAM-NPP")
# contour(h, add=TRUE)
plot(sand_06[[6]]-surf_06[[6]], xlim=xlim, ylim=ylim, breaks=breaks, col=diverge_hsv(99), main="FAM-NPP")
# contour(h, add=TRUE)
plot(sand_07[[6]]-surf_07[[6]], xlim=xlim, ylim=ylim, breaks=breaks, col=diverge_hsv(99), main="FAM-NPP")
# contour(h, add=TRUE)
plot(sand_08[[6]]-surf_08[[6]], xlim=xlim, ylim=ylim, breaks=breaks, col=diverge_hsv(99), main="FAM-NPP")
# contour(h, add=TRUE)



data.dat5 <- "C:/Users/jjansen/Desktop/science/data_environmental/Circumpolar_ROMS/10km_outputs/sed_test5/"
## sediments
sand5_01 <- brick(paste0(data.dat5,"ocean_avg_0003.nc"), varname="sand_01", level=1)
sand5_02 <- brick(paste0(data.dat5,"ocean_avg_0003.nc"), varname="sand_02", level=1)
sand5_03 <- brick(paste0(data.dat5,"ocean_avg_0003.nc"), varname="sand_03", level=1)
sand5_04 <- brick(paste0(data.dat5,"ocean_avg_0003.nc"), varname="sand_04", level=1)
sand5_05 <- brick(paste0(data.dat5,"ocean_avg_0003.nc"), varname="sand_05", level=1)
sand5_06 <- brick(paste0(data.dat5,"ocean_avg_0003.nc"), varname="sand_06", level=1)
sand5_07 <- brick(paste0(data.dat5,"ocean_avg_0003.nc"), varname="sand_07", level=1)
sand5_08 <- brick(paste0(data.dat5,"ocean_avg_0003.nc"), varname="sand_08", level=1)

par(mfrow=c(2,4))
plot(sand5_01[[6]]-surf_01[[6]], xlim=xlim, ylim=ylim, breaks=breaks, col=diverge_hsv(99), main="FAM-NPP")
plot(sand5_02[[6]]-surf_02[[6]], xlim=xlim, ylim=ylim, breaks=breaks, col=diverge_hsv(99), main="FAM-NPP")
plot(sand5_03[[6]]-surf_03[[6]], xlim=xlim, ylim=ylim, breaks=breaks, col=diverge_hsv(99), main="FAM-NPP")
plot(sand5_04[[6]]-surf_04[[6]], xlim=xlim, ylim=ylim, breaks=breaks, col=diverge_hsv(99), main="FAM-NPP")
plot(sand5_05[[6]]-surf_05[[6]], xlim=xlim, ylim=ylim, breaks=breaks, col=diverge_hsv(99), main="FAM-NPP")
plot(sand5_06[[6]]-surf_06[[6]], xlim=xlim, ylim=ylim, breaks=breaks, col=diverge_hsv(99), main="FAM-NPP")
plot(sand5_07[[6]]-surf_07[[6]], xlim=xlim, ylim=ylim, breaks=breaks, col=diverge_hsv(99), main="FAM-NPP")
plot(sand5_08[[6]]-surf_08[[6]], xlim=xlim, ylim=ylim, breaks=breaks, col=diverge_hsv(99), main="FAM-NPP")


par(mfrow=c(1,1))
plot(sand_02[[6]]-sand5_02[[6]], xlim=xlim, ylim=ylim, breaks=breaks, col=diverge_hsv(99), main="FAM-NPP")


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
settle_01 <- brick(paste0(data.dat,"ocean_his_0003.nc"), varname="sandfrac_01", level=1)
plot(settle_01[[1]])

breaks=seq(-0.0001,0.0001,length.out = 255)
col=diverge_hcl(256)

set8.6 <- settle_08[[6]]-settle_08[[1]]
set6.6 <- settle_06[[6]]-settle_06[[1]]
set4.6 <- settle_04[[6]]-settle_04[[1]]
set2.6 <- settle_02[[6]]-settle_02[[1]]

plot(set8.6, breaks=breaks, col=col)
plot(set6.6, breaks=breaks, col=col)
plot(set4.6, breaks=breaks, col=col)
plot(set2.6, breaks=breaks, col=col)


##### Evaluate floats
### it's being stored in a rather weird way...
### if you read in as a raster, the rows are the time-steps
data.flts <- "C:/Users/jjansen/Desktop/science/data_environmental/Circumpolar_ROMS/10km_outputs/"
grd.x <- brick(paste0(data.flts,"ocean_flt.nc"), varname="Xgrid", level=1)
grd.y <- brick(paste0(data.flts,"ocean_flt.nc"), varname="Ygrid", level=1)
grd.z <- brick(paste0(data.flts,"ocean_flt.nc"), varname="Zgrid", level=1)
flts.x <- brick(paste0(data.flts,"ocean_flt.nc"), varname="x", level=1)
flts.y <- brick(paste0(data.flts,"ocean_flt.nc"), varname="y", level=1)
flts.z <- brick(paste0(data.flts,"ocean_flt.nc"), varname="depth", level=1)

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
