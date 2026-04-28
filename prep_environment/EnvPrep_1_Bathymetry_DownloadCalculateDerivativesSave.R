###############################################################################
###### Environmental data preparation - Bathymetry & distance to canyons ######
###############################################################################

## specify user and setup directory to look up data from
usr <- "VM"
#usr <- "SJ"
#usr <- "JJ
source("prep_environment/EnvPrep_0_SourceFile.R")

## set input and output folders
env.dir <- paste0(usr.main.dir,"/data_environmental/raw/")
out.dir <- paste0(usr.main.dir,"/data_environmental/derived/bathy_outputs/")
local.out <- "~/temp"

########################################
## download IBCSO V2 bed and ice layers
url <- "https://download.pangaea.de/dataset/937574/files/IBCSO_v2_bed.tif"
url2<- "https://download.pangaea.de/dataset/937574/files/IBCSO_v2_ice-surface.tif"
# Define the destination path
destfile <- paste0(env.dir,"IBCSO_v2_bed.tif")
destfile2<- paste0(env.dir,"IBCSO_v2_ice-surface.tif")
# Download the file
download.file(url, destfile)
download.file(url2, destfile2)

library(terra)
## read ibcso files
ibcso_bed <- rast(paste0(env.dir,"IBCSO_v2_bed.tif"))
ibcso_ice <- rast(paste0(env.dir,"IBCSO_v2_ice-surface.tif"))

## we use the bed layer and calculate slope and topographic position index at different scales
r.slope <- terrain(ibcso_bed)
r.tpi   <- terrain(ibcso_bed, v="TPI")
r.tpi5  <- terrain(ibcso_bed, v="TPI", scale=5)
r.tpi11 <- terrain(ibcso_bed, v="TPI", scale=11)

r <- c(ibcso_bed, ibcso_ice, r.slope, r.tpi, r.tpi5, r.tpi11)
names(r) <- c("depth", "ibcso_ice", "slope", "tpi", "tpi5", "tpi11")

## at 2km resolution:
r.2km <- terra::aggregate(r, 4)
r.2km$depth[r.2km$ibcso_ice>0] <- NA
r.2km$depth[r.2km$depth>0] <- NA
r.2km$slope[is.na(r.2km$depth[])] <- NA
r.2km$tpi  [is.na(r.2km$depth[])] <- NA
r.2km$tpi5 [is.na(r.2km$depth[])] <- NA
r.2km$tpi11[is.na(r.2km$depth[])] <- NA
writeRaster(r.2km[[c(1,3:6)]], filename=paste0(local.out, "IBCSO_v2_2km_bathymetric_variables.tif"), overwrite=TRUE)
localr.2km <- paste0(local.out, "IBCSO_v2_2km_bathymetric_variables.tif")
# command line code to transfer output from temp to dropbox
system2(
  "rclone",
  args = c(
    "copy",
    localr.2km,
    "dropbox:Data/data_environmental/derived/bathy_outputs/",
    "--progress"
  )
)
## at 500m resolution:
r.500m <- r
r.500m$depth[r.500m$ibcso_ice>0] <- NA
r.500m$depth[r.500m$depth>0] <- NA
r.500m$slope[is.na(r.500m$depth[])] <- NA
r.500m$tpi  [is.na(r.500m$depth[])] <- NA
r.500m$tpi5 [is.na(r.500m$depth[])] <- NA
r.500m$tpi11[is.na(r.500m$depth[])] <- NA
writeRaster(r.500m[[c(1,3:6)]], filename=paste0(local.out, "IBCSO_v2_500m_bathymetric_variables.tif"), overwrite=TRUE)
localr.500m <- paste0(local.out, "IBCSO_v2_500m_bathymetric_variables.tif")
# command line code to transfer output from temp to dropbox
system2(
  "rclone",
  args = c(
    "copy",
    localr.500m,
    "dropbox:Data/data_environmental/derived/bathy_outputs/",
    "--progress"
  )
)
##################################################
#### Distance to underwater canyons identified from Arosio & Amblas 2025
r.500m <- rast(paste0(env.dir, "IBCSO_v2_500m_bathymetric_variables.tif"))
canyons <- vect("/perm_storage/shared_space/BioMAS/environmental_data/ArosioAmblasAntarcticCanyons/Antarctica_drainage_2025_complete.shp")
## filter to canyons that reach shallower than 3000m
canyons_shallow <- subset(canyons, canyons$Z_Max > -3000)
## Create ocean mask from which to calculate distances: water=1, land/ice=NA
ocean <- classify(r.500m$depth, rbind(c(-Inf, -1e-6, 1), c(0, Inf, NA)))
## Rasterize canyons on the grid and restrict to ocean
canyon_r <- rasterize(canyons_shallow, r.500m$depth, field=1, touches=TRUE)
## Identify target cells (canyon-on-water cells)
ra.comb <- sum(ocean,canyon_r, na.rm=TRUE)
## calculate distance to canyons (which have value = 2)
dist_water_m <- gridDist(ra.comb, target=2)

writeRaster(dist_water_m, filename=paste0(out.dir, "IBCSO_v2_500m_DistanceToCanyons.tif"), overwrite=TRUE)
dist_water_m_2km <- aggregate(dist_water_m, 4)
writeRaster(dist_water_m, filename=paste0(out.dir, "IBCSO_v2_2km_DistanceToCanyons.tif"), overwrite=TRUE)


###############################################
#### Geomorphology, based on https://github.com/jacquomo/Geomorphmetry_SeamapAus/blob/main/Geomorphometry_for_Australian_Marine_Parks.Rmd
library(whitebox)

# 2km res
# need to write 2km depth raster to disk as whitebox wont read the in memory raster
writeRaster(r.2km$depth, paste0(env.dir, "IBCSO_v2_2km_depth.tif"), overwrite = TRUE)
dem <- rast(paste0(env.dir, "IBCSO_v2_2km_depth.tif"))
wbt_geomorphons(
  dem = dem,
  output = paste0("IBCSO_v2_2km_geomorph", ".tif"),
  search = 20,
  threshold = 5,
  fdist = 20,
  wd = out.dir
)

# 500m res
# need to write 2km depth raster to disk as whitebox wont read the in memory raster
writeRaster(r.2km$depth, paste0(env.dir, "IBCSO_v2_2km_depth.tif"), overwrite = TRUE)
dem <- rast(paste0(env.dir, "IBCSO_v2_2km_depth.tif"))
wbt_geomorphons(
  dem = dem,
  output = paste0("IBCSO_v2_2km_geomorph", ".tif"),
  search = 20,
  threshold = 5,
  fdist = 20,
  wd = out.dir
)
