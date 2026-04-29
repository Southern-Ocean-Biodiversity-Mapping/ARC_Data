## specify user and setup directory to look up data from
usr <- "VM"
source("0_SourceFile.R")

## set input and output folders
env.dir <- paste0(usr.main.dir,"data_environmental/derived/NPP")

###############################################
library(terra)

## load monthly data
chl   <- rast(file.path(env.dir, "NPP_monthly_OctMar_chla.tif"))
cafe.filled   <- rast(file.path(env.dir, "NPP_filled_monthly_OctMar_12boxes_cafe.tif"))
cafe.unfilled <- rast(file.path(env.dir, "NPP_unfilled_monthly_OctMar_cafe.tif"))
eppl.filled   <- rast(file.path(env.dir, "NPP_filled_monthly_OctMar_12boxes_eppley.tif"))
eppl.unfilled <- rast(file.path(env.dir, "NPP_unfilled_monthly_OctMar_eppley.tif"))
vpmg.filled   <- rast(file.path(env.dir, "NPP_filled_monthly_OctMar_12boxes_vpmg.tif"))
vpmg.unfilled <- rast(file.path(env.dir, "NPP_unfilled_monthly_OctMar_vpmg.tif"))
cbpm.filled   <- rast(file.path(env.dir, "NPP_filled_monthly_OctMar_12boxes_cbpm.tif"))
cbpm.unfilled <- rast(file.path(env.dir, "NPP_unfilled_monthly_OctMar_cbpm.tif"))

## long-term climatology
clim.cafe.filled.sum   <- sum(cafe.filled[[1:108]], na.rm=TRUE)
clim.cafe.filled.mean <- clim.cafe.filled.sum/108
clim.eppl.filled.sum   <- sum(eppl.filled[[1:108]], na.rm=TRUE)
clim.eppl.filled.mean <- clim.eppl.filled.sum/108
clim.vpmg.filled.sum   <- sum(vpmg.filled[[1:108]], na.rm=TRUE)
clim.vpmg.filled.mean <- clim.vpmg.filled.sum/108
clim.cbpm.filled.sum   <- sum(cbpm.filled[[1:108]], na.rm=TRUE)
clim.cbpm.filled.mean <- clim.cbpm.filled.sum/108
clim.cafe.filled.sd   <- stdev(cafe.filled[[1:108]], na.rm=TRUE)
clim.eppl.filled.sd   <- stdev(eppl.filled[[1:108]], na.rm=TRUE)
clim.vpmg.filled.sd   <- stdev(vpmg.filled[[1:108]], na.rm=TRUE)
clim.cbpm.filled.sd   <- stdev(cbpm.filled[[1:108]], na.rm=TRUE)

names(clim.cafe.filled.mean) <- names(clim.eppl.filled.mean) <- names(clim.vpmg.filled.mean) <- names(clim.cbpm.filled.mean) <- "mean"
names(clim.cafe.filled.sd) <- names(clim.eppl.filled.sd) <- names(clim.vpmg.filled.sd) <- names(clim.cbpm.filled.sd) <- "sd"

writeRaster(clim.cafe.filled.sum,  filename=file.path(env.dir,"NPP_climatology_OctMar_2002To2020_filled12boxes_cafe_sum.tif"))
writeRaster(clim.cafe.filled.mean, filename=file.path(env.dir,"NPP_climatology_OctMar_2002To2020_filled12boxes_cafe_mean.tif"))
writeRaster(clim.eppl.filled.sum,  filename=file.path(env.dir,"NPP_climatology_OctMar_2002To2020_filled12boxes_eppl_sum.tif"))
writeRaster(clim.eppl.filled.mean, filename=file.path(env.dir,"NPP_climatology_OctMar_2002To2020_filled12boxes_eppl_mean.tif"))
writeRaster(clim.vpmg.filled.sum,  filename=file.path(env.dir,"NPP_climatology_OctMar_2002To2020_filled12boxes_vpmg_sum.tif"))
writeRaster(clim.vpmg.filled.mean, filename=file.path(env.dir,"NPP_climatology_OctMar_2002To2020_filled12boxes_vpmg_mean.tif"))
writeRaster(clim.cbpm.filled.sum,  filename=file.path(env.dir,"NPP_climatology_OctMar_2002To2020_filled12boxes_cbpm_sum.tif"))
writeRaster(clim.cbpm.filled.mean, filename=file.path(env.dir,"NPP_climatology_OctMar_2002To2020_filled12boxes_cbpm_mean.tif"))
writeRaster(clim.cafe.filled.sd,  filename=file.path(env.dir,"NPP_climatology_OctMar_2002To2020_filled12boxes_cafe_sd.tif"))
writeRaster(clim.eppl.filled.sd,  filename=file.path(env.dir,"NPP_climatology_OctMar_2002To2020_filled12boxes_eppl_sd.tif"))
writeRaster(clim.vpmg.filled.sd,  filename=file.path(env.dir,"NPP_climatology_OctMar_2002To2020_filled12boxes_vpmg_sd.tif"))
writeRaster(clim.cbpm.filled.sd,  filename=file.path(env.dir,"NPP_climatology_OctMar_2002To2020_filled12boxes_cbpm_sd.tif"))
