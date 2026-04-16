####################################################################################################################
### Code to extract and format bottom layer of B-SOSE BGC variables (netcdf files provided by Matthew Marzloff)  ###
### Code to format CSIRO Atlas of Regional Seas (CARS) BGC bottom variables                                      ###
####################################################################################################################


##########################################################
## 1) Set-up ---
###########################################################
library(raster)
library(ncdf4)
library(tidyr)
library(dplyr)
#library(ncdfgeom)

path<-"C:\\Users\\hillna\\Downloads\\"
#raw netcdf files can be found on cloudstor folder:
#https://owncloud.imas-data-service.cloud.edu.au/index.php/s/ORxSWb6xbJRWfNI


sci.dir <-      "C:/Users/hillna/OneDrive - University of Tasmania/UTAS_work/Projects/Benthic Diversity ARC/"
env.dir <-  paste0(sci.dir,"data_environmental/")
tools.dir <-    paste0(sci.dir,"Analysis/Useful_Functions_Tools/")
ARC_Data.dir <- paste0(sci.dir,"Analysis/ARC_Data/")

path2<-"C:\\Users\\hillna\\OneDrive - University of Tasmania\\UTAS_work\\Projects\\Benthic Diversity ARC\\data_environmental\\"

SO.ext <- extent(-180, 180, -80, -55)
SO.ext2 <- extent(0, 360, -80, -55)
# bathy data as template for projected rasters
bathy_shelf<-raster(paste0(path2, "derived\\Circumpolar_EnvData_500m_shelf_bathy_gebco_depth"))

#naming conventions
string.chr <- "Circumpolar_EnvData_"
string.res <- "500m_"

#################################################################
## 2) B-SOSE variables----
#################################################################
#each of these files is 21 GB. Open one at at time, run through code, save outputs and close netcdf; and repeat for next variables

#arag_nc<-nc_open(paste0(path, "bsose_I135_2013to2019_monthly_OmegaArag.nc"))
#no3_nc<-nc_open(paste0(path, "bsose_I135_2013to2019_monthly_NO3.nc"))
#o2_nc<-nc_open(paste0(path, "bsose_I135_2013to2019_monthly_O2.nc"))
po4_nc<-nc_open(paste0(path, "bsose_I135_2013to2019_monthly_PO4.nc"))


#get space and time variables. Same for all net cdfs
lon.nc<-ncvar_get(po4_nc, "XC")
lat.nc<-ncvar_get(po4_nc, "YC")
depth.nc<-ncvar_get(po4_nc, "Z")
time.nc<-ncvar_get(po4_nc, "time")


#extract BGC values for every time slice where value !=0 (which is land/seafloor)

#set up array for values for each time slice
Kbtm = array(dim=c(length(lon.nc), length(lat.nc), 84))

# for each time slice evaluate whether value for each lat/lon is zero (land/seafloor) or a valid value- leave NAs in array where is land/seafloor
for(t in 1:84){
 
 #temp_nc<-ncvar_get(nc=arag_nc, varid="BLGOMAR", start=c(1,1,1,t), count=c(-1,-1,-1,1))
  #temp_nc<-ncvar_get(nc=no3_nc, varid="TRAC04", start=c(1,1,1,t), count=c(-1,-1,-1,1))
  #temp_nc<-ncvar_get(nc=o2_nc, varid="TRAC03", start=c(1,1,1,t), count=c(-1,-1,-1,1))
  temp_nc<-ncvar_get(nc=po4_nc, varid="TRAC05", start=c(1,1,1,t), count=c(-1,-1,-1,1))
  
print(t)

for (k in 1:dim(temp_nc)[3]) 
  {
for (j in 1:dim(temp_nc)[2])
    {
for (i in 1:dim(temp_nc)[1]) 
      {

  if (temp_nc[i,j,k] != 0) {
  Kbtm[i,j,t] = temp_nc[i,j,k]} #gives actual variable values
  #Kbtm[i,j] = k}               #gives index of seafloor depth layer
        }
    }
  }
}

#save(Kbtm, file=paste0(env.dir, "raw/arag_bot_vals.RData"))
#save(Kbtm, file=paste0(env.dir, "raw/no3_bot_vals.RData"))
#save(Kbtm, file=paste0(env.dir, "raw/o2_bot_vals.RData"))
save(Kbtm, file=paste0(env.dir, "raw/po4_bot_vals.RData"))

#nc_close(arag_nc)
#nc_close(no3_nc)
#nc_close(o2_nc)
nc_close(po4_nc)

#average and sd of values and turn into spatial points dataframe
bot_mean<-apply(Kbtm, c(1,2), mean)
bot_sd<-apply(Kbtm, c(1,2), sd)

##turn into long format 
rownames(bot_mean)<-rownames(bot_sd)<-lon.nc
colnames(bot_mean)<-colnames(bot_sd)<-lat.nc

#mean
bot_mean_long<-tibble::rownames_to_column(as.data.frame(bot_mean), "lon") %>%
    pivot_longer(., cols= -lon, names_to="lat", values_to="BGC")

bot_mean_long$lon<-as.numeric(bot_mean_long$lon)
bot_mean_long$lat<-as.numeric(bot_mean_long$lat)

#sd
bot_sd_long<-tibble::rownames_to_column(as.data.frame(bot_sd), "lon") %>%
  pivot_longer(., cols= -lon, names_to="lat", values_to="BGC")

bot_sd_long$lon<-as.numeric(bot_sd_long$lon)
bot_sd_long$lat<-as.numeric(bot_sd_long$lat)

## then spatial points data frame then rasterize
#mean
bot_mean_sp<-SpatialPointsDataFrame(coords=bot_mean_long[,1:2], data=bot_mean_long[,3])
temp_rast<-raster(ext=extent(bot_mean_sp), resolution= 0.1667)
bot_mean_rast<-rasterize(bot_mean_sp, temp_rast, fun=mean, field="BGC")

#sd
bot_sd_sp<-SpatialPointsDataFrame(coords=bot_sd_long[,1:2], data=bot_sd_long[,3])
bot_sd_rast<-rasterize(bot_sd_sp, temp_rast, fun=mean, field="BGC")


## stack, then crop rasters to SO extent, project to polar stereographic, then crop to shelf.
bot_stack<-stack(bot_mean_rast, bot_sd_rast)
bot_stack<-crop(bot_stack,SO.ext2)

#because 0 to 360 degrees- projection doesn't work very well. Clumsy fix
#mean
bot_temp<-raster::shift(bot_stack, dx = -360)
origin(bot_temp)<-origin(bot_stack)
bot_stack<-raster::merge(bot_stack, bot_temp)


bot_500<-projectRaster(bot_stack, bathy_shelf)

bot_500_shelf<-mask(bot_500, bathy_shelf)

#names(bot_500_shelf)<-c("arag_mean", "arag_sd")
#names(bot_500_shelf)<-c("no3_mean", "no3_sd")
#names(bot_500_shelf)<-c("o2_mean", "o2_sd")
names(bot_500_shelf)<-c("po4_mean", "po4_sd")


## Save final rasters
#arag
#writeRaster(bot_500_shelf, paste0(env.dir, "derived/", string.chr, string.res, "shelf_bot_arag_BSOSE"))
#no3
#writeRaster(bot_500_shelf, paste0(env.dir, "derived/", string.chr, string.res, "shelf_bot_NO3_BSOSE"))
#o2
#writeRaster(bot_500_shelf, paste0(env.dir, "derived/", string.chr, string.res, "shelf_bot_O2_BSOSE"))
#po4
writeRaster(bot_500_shelf, paste0(env.dir, "derived/", string.chr, string.res, "shelf_bot_PO4_BSOSE"))

###########################################################
## 3) CARS variables ----
###########################################################

#WGS84 lat/lon, 0.5 *0.5 deg resolution

NO3<-stack(brick(paste0(env.dir, "raw\\nitrate_cars2009_bot.nc"),varname="mean"),
           brick(paste0(env.dir, "raw\\nitrate_cars2009_bot.nc"),varname="seas_range"),
           brick(paste0(env.dir, "raw\\nitrate_cars2009_bot.nc"),varname="std_dev"))

PO4<-stack(brick(paste0(env.dir, "raw\\phosphate_cars2009_bot.nc"),varname="mean"),
           brick(paste0(env.dir, "raw\\phosphate_cars2009_bot.nc"),varname="seas_range"),
           brick(paste0(env.dir, "raw\\phosphate_cars2009_bot.nc"),varname="std_dev"))

O2<-stack(brick(paste0(env.dir, "raw\\oxygen_cars2009_bot.nc"),varname="mean"),
           brick(paste0(env.dir, "raw\\oxygen_cars2009_bot.nc"),varname="seas_range"),
           brick(paste0(env.dir, "raw\\oxygen_cars2009_bot.nc"),varname="std_dev"))

names(NO3)<-names(PO4)<-names(O2)<-c("mean", "seas_range", "std_dev")

#crop rasters to SO extent, project to polar stereographic, then crop to shelf.
NO3<-crop(NO3, SO.ext2)
PO4<-crop(PO4, SO.ext2)
O2<-crop(O2, SO.ext2)

#raster seems to have issues with 0-360 deg!!??
NO3<-raster::merge(NO3, raster::shift(NO3, dx = -360))
NO3_500<-projectRaster(NO3, bathy_shelf)

PO4<-raster::merge(PO4, raster::shift(PO4, dx = -360))
PO4_500<-projectRaster(PO4, bathy_shelf)

O2<-raster::merge(O2, raster::shift(O2, dx = -360))
O2_500<-projectRaster(O2, bathy_shelf)

names(NO3_500)<- names(PO4_500)<- names(O2_500)<- c("mean", "seas_range", "std_dev")

NO3_500_shelf<-mask(NO3_500, bathy_shelf)
PO4_500_shelf<-mask(PO4_500, bathy_shelf)
O2_500_shelf<-mask(O2_500, bathy_shelf)

#save rasters
writeRaster(NO3_500_shelf, filename=paste0(env.dir, "derived\\", string.chr, string.res, "shelf_CARS_NO3"), overwrite=TRUE)
writeRaster(PO4_500_shelf, filename=paste0(env.dir, "derived\\", string.chr, string.res, "shelf_CARS_PO4"), overwrite=TRUE)
writeRaster(O2_500_shelf, filename=paste0(env.dir, "derived\\", string.chr, string.res, "shelf_CARS_O2"), overwrite=TRUE)


