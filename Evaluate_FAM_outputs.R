## Evaluate ROMS outputs
library(raster)
library(colorspace)
data.dat <- "C:/Users/jjansen/Desktop/science/data_environmental/Circumpolar_ROMS/10km_outputs/sed_test1/"

xlim=c(250,350)
ylim=c(100,200)

#### seafloor-layer is 1, while surface-layer is 31
## depth
h <- raster(paste0(data.dat,"ocean_avg_0003.nc"), varname="h", level=1)
## currents
u <- brick(paste0(data.dat,"ocean_avg_0003.nc"), varname="u", level=1)
v <- brick(paste0(data.dat,"ocean_avg_0003.nc"), varname="v", level=1)
uv <- sqrt(u^2+v^2)
plot(uv[[1]])
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
## sediments
sand_01 <- brick(paste0(data.dat,"ocean_avg_0003.nc"), varname="sand_01", level=1)
sand_02 <- brick(paste0(data.dat,"ocean_avg_0003.nc"), varname="sand_02", level=1)
sand_03 <- brick(paste0(data.dat,"ocean_avg_0003.nc"), varname="sand_03", level=1)
sand_04 <- brick(paste0(data.dat,"ocean_avg_0003.nc"), varname="sand_04", level=1)
sand_05 <- brick(paste0(data.dat,"ocean_avg_0003.nc"), varname="sand_05", level=1)
sand_06 <- brick(paste0(data.dat,"ocean_avg_0003.nc"), varname="sand_06", level=1)
sand_07 <- brick(paste0(data.dat,"ocean_avg_0003.nc"), varname="sand_07", level=1)
sand_08 <- brick(paste0(data.dat,"ocean_avg_0003.nc"), varname="sand_08", level=1)

breaks <- seq(-0.000025,0.000025,length.out=99)

par(mfrow=c(2,2))
plot(h, xlim=xlim, ylim=ylim, main="depth")
plot(uv[[6]], xlim=xlim, ylim=ylim, main="currents")
plot(surf_01[[6]], xlim=xlim, ylim=ylim, main="NPP")
plot(sand_01[[6]], xlim=xlim, ylim=ylim, main="FAM")

par(mfrow=c(2,2))
plot(surf_01[[6]], xlim=xlim, ylim=ylim, main="NPP")
contour(h, add=TRUE)
plot(sand_01[[6]], xlim=xlim, ylim=ylim, main="FAM")
contour(h, add=TRUE)
plot(uv[[6]], xlim=xlim, ylim=ylim, main="currents")
contour(h, add=TRUE)
plot(sand_01[[6]]-surf_01[[6]], xlim=xlim, ylim=ylim, breaks=breaks, col=diverge_hsv(99), main="FAM-NPP")
contour(h, add=TRUE)

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
