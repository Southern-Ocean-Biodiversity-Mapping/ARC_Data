
#title: "Circumpolar Environmental Data- part 2"
#author: "Nicole Hill"
#date: "15/09/2021"

### This code extracts satellite- derived circumpolar SST, SSH, and sea-ice using the package 'raadtools' and was run on the project's AntPro VM:
### https://ant.antarctic-biodiversity-modelling.cloud.edu.au/
### Daily files were extracted where possibel and used to calculate climatologies
### An update geomorphology layer was also sourced from Alix Post (GA) and converted to a raster.
### extracting these files and generating climatologies takes a LONG time- even on the VM. if no parallel processing then a couple of hours.

# libraries and paths 
library(raadtools)
library(raster)       ## package for raster manipulation
library(sp)
library(dplyr)
library(rgdal)        ## package for geospatial analysis
library(foreach)      ## parallel processing on VM
library(doParallel)   ## parallel processing on VM


VM_path2<-"/perm_storage/shared_space/BioMAS/environmental_data/"


### 1) set up details for extraction ----
# polar stereographic projection:
stereo <- "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# extent
SO.ext <- extent(-180, 180, -80, -55)
SO.ext2 <- extent(0, 360, -80, -55)
SO.stereo<-extent(-3333500, 3334000, -3337000, 3333500) 

#dates
dts<-seq(from= as.Date("2002-01-01"), to= as.Date("2020-12-31"), by= "1 day")

#Create indices and function to generate seasonal stats
# Based on austral summer 
years<-2002:2020

#Spring: Sep, Oct, Nov
spring_ind<-c()
for(i in 1:length(years)){
  spring_ind<-append(spring_ind, 
                     seq(as.Date(paste0(years[i], "-09-01")), as.Date(paste0(years[i] , "-11-30")), by = "1 day"))
}

#Summer: Dec, Jan, Feb
summer_ind<-seq(as.Date("2002-12-01"), as.Date("2002-12-31"), by = "1 day")
for(i in 2:length(years)){
  summer_ind<-append(summer_ind, 
                     seq(as.Date(paste0(years[i], "-01-01")), as.Date(paste0(years[i], "-03-31")), by = "1 day"))
}


# Function to generate seasonal stats

seas_stats<-function(stack,                   # name of raster stack or brick
                     season_dates,            # sequence of dates to extract ( identified in Zvalues of raster layers)
                     stat)                    # statsitic to calculate (mean, sd, min, max or anything currently recognised by fun)
  {
  get_rasts<-subset(stack, which(getZ(stack) %in% season_dates))
  stat<-calc(get_rasts, fun= stat, na.rm=TRUE)                  
  
}


#naming conventions
string.chr <- "Circumpolar_EnvData_"
string.res <- "500m_"

# bathy data as template for projected rasters
#bathy<-raster(VM_path2, "Circumpolar_EnvData_500m_bathy_gebco_depth")
bathy_shelf<-raster(paste0(VM_path2, "Circumpolar_EnvData_500m_shelf_bathy_gebco_depth"))
#bathy_shelf<-raster("C:\\Users\\hillna\\OneDrive - University of Tasmania\\UTAS_work\\Projects\\Benthic Diversity ARC\\data_environmental\\derived\\Circumpolar_EnvData_500m_shelf_bathy_gebco_depth")


## 2) Extract Optimally Interpolated SST (OISST) ----
#native resolution =0.25deg
#crs=  +proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 
# No need to mask with ice cover as interpolated product where is ice is  ~ -1.7 deg

## Parallelised to use multiple cores and reduce processing time

#Extract daily data
UseCores=12
c1<-makeCluster(UseCores)
registerDoParallel(c1)

  sst<- foreach(i =1: length(dts)) %dopar% {
  library(raadtools)
  library(raster)
  temp_sst<-readsst(dts[i], time.resolution = "daily", xylim = SO.ext, varname="sst", setNA=TRUE)
}

daily_sst<-stack(sst)
daily_sst<-setZ(daily_sst, dts) #this allows to pull out select dates to calculate seasonal average etc
names(daily_sst)<-dts

stopCluster(c1)

writeRaster(daily_sst, filename=paste0(VM_path2, "daily_sst")) 

# generate yearly, spring, summer, climatologies (mean and sd)
sst_mean<-calc(daily_sst, mean)
sst_sd<-calc(daily_sst, sd)

sst_sp_mean<-seas_stats(stack= daily_sst, season_dates = spring_ind, stat=mean )
sst_sp_sd<-seas_stats(stack= daily_sst, season_dates = spring_ind, stat=sd )

sst_su_mean<-seas_stats(stack= daily_sst, season_dates = summer_ind, stat=mean )
sst_su_sd<-seas_stats(stack= daily_sst, season_dates = summer_ind, stat= sd )

sst_stack<-stack(sst_mean, sst_sd, sst_sp_mean,sst_sp_sd, sst_su_mean, sst_su_sd)
names(sst_stack)<-c ("sst_mean", "sst_sd", "sst_sp_mean", "sst_sp_sd", "sst_su_mean", "sst_su_sd")

#project raster and interpolate to 500m
sst_stack_500<-projectRaster(sst_stack, bathy_shelf)

#crop rasters to shelf
sst_stack_500m_shelf<-mask(sst_stack_500, bathy_shelf)

#save rasters
writeRaster(sst_stack, filename=paste0(VM_path2, string.chr, "SST"), overwrite=TRUE)
writeRaster(sst_stack_500, filename=paste0(VM_path2, string.chr, string.res, "SST"), overwrite=TRUE)
writeRaster(sst_stack_500m_shelf, filename=paste0(VM_path2, string.chr, string.res, "shelf_SST"), overwrite=TRUE)

#remove intermediate rasters
rm(daily_sst, sst_mean, sst_sd, sst_sp_mean,sst_sp_sd, sst_su_mean, sst_su_sd)


## 3) Sea surface height----
#resolution 0.25 x 0.25 deg
#crs= +proj=longlat +rf=298.257 +a=6378136.3 

#issue with raster extents in first attempt
ssh_names<-as.Date(sshfiles()$date)
date_range <- seq(min(ssh_names), max(ssh_names), by = 1) 
date_range[!date_range %in% ssh_names] 
#no dates missing!

test<-unlist(sapply(ssh, function(x) crs(x), simplify=TRUE))
test_ext<-sapply(ssh, function(x) extent(x), simplify=TRUE)

allSame <- function(x) length(unique(x)) == 1

allSame(test_ext)
unique(test_ext)# one set of rasters has extent of -180, 180. Now to find it!

get_ind<-lapply(ssh, function(x) if (xmin(extent(x)) == -180) print (getZ(x)))
# from "2020-06-04" extent is -180 to 180, before that it is 0 to 360 and so stack has been cut-off at 180....
ssh_dts<-seq(as.Date("2002-01-01"), as.Date("2020-06-03"), by = 1) 


#Issues with memory when extracting and stacking daily rasters. So try extracting daily raster for one month, generate aveage value and keeping only monthly rasters
ssh_dts<-as.data.frame(ssh_dts)
ssh_dts$MonthYear<-format(ssh_dts$ssh_dts, "%Y-%m")

ssh_monthly<-stack()
for(i in 1:length(unique(ssh_dts$MonthYear))){
  print(unique(ssh_dts$MonthYear)[i])
mth_dts<-ssh_dts$ssh_dts[ssh_dts$MonthYear %in% unique(ssh_dts$MonthYear)[i]]
daily_ssh<-stack()
for (j in 1:length(mth_dts)){
temp_ssh<-readssh(mth_dts[j], time.resolution = "daily", xylim = SO.ext2, setNA=TRUE)
daily_ssh<-stack(daily_ssh, temp_ssh)
}
ssh_monthly<-stack(ssh_monthly, calc(daily_ssh, mean))
}
names(ssh_monthly)<-unique(ssh_dts$MonthYear)
ssh_monthly<-setZ(ssh_monthly,unique(ssh_dts$MonthYear))

#writeRaster(ssh_monthly, filename = paste0(VM_path2, "monthly_ssh"))

## generate yearly and seasonal climatology stats
# monthly indices
spring_mth_ind<-c()
for(i in 1:length(years)){
  spring_mth_ind<-append(spring_mth_ind, 
                     format(seq(as.Date(paste0(years[i], "-09-01")), as.Date(paste0(years[i] , "-11-01")), by = "1 month"), "%Y-%m"))
}

#Summer: Dec, Jan, Feb
summer_mth_ind<-format(as.Date("2002-12-01"), "%Y-%m")
for(i in 2:length(years)){
  summer_mth_ind<-append(summer_mth_ind, 
                         format(seq(as.Date(paste0(years[i], "-01-01")), as.Date(paste0(years[i] , "-02-01")), by = "1 month"), "%Y-%m"))
}

#raster layers
ssh_mean<-calc(monthly_ssh, mean, na.rm=TRUE)
ssh_sd<-calc(monthly_ssh, sd, na.rm=TRUE)

ssh_sp_mean<-seas_stats(stack= monthly_ssh, season_dates = spring_mth_ind, stat=mean )
ssh_sp_sd<-seas_stats(stack= monthly_ssh, season_dates = spring_mth_ind, stat=sd )

ssh_su_mean<-seas_stats(stack= monthly_ssh, season_dates = summer_mth_ind, stat=mean )
ssh_su_sd<-seas_stats(stack= monthly_ssh, season_dates = summer_mth_ind, stat= sd )

ssh_stack<-stack(ssh_mean, ssh_sd, ssh_sp_mean,ssh_sp_sd, ssh_su_mean, ssh_su_sd)
names(ssh_stack)<-c ("ssh_mean", "ssh_sd", "ssh_sp_mean", "ssh_sp_sd", "ssh_su_mean", "ssh_su_sd")

#project raster and interpolate to 500m
#issues projecting because data is 0 to 360 deg; 'fix' from Ben R below
ssh_stack_temp<-raster::merge(ssh_stack, raster::shift(ssh_stack, dx = -360))
ssh_stack_500<-projectRaster(ssh_stack_temp, bathy_shelf)

#crop rasters to shelf
ssh_stack_500m_shelf<-mask(ssh_stack_500, bathy_shelf)

#save rasters
writeRaster(ssh_stack, filename=paste0(VM_path2, string.chr, "SSH"))
writeRaster(ssh_stack_500, filename=paste0(VM_path2, string.chr, string.res, "SSH")) 
writeRaster(ssh_stack_500m_shelf, filename=paste0(VM_path2, string.chr, string.res, "shelf_SSH"), overwrite=TRUE)

#remove intermediate rasters
rm(monthly_ssh, ssh_mean, ssh_sd, ssh_sp_mean,ssh_sp_sd, ssh_su_mean, ssh_su_sd)


## 4) Sea ice concentration: amsr product ----
#native resolution: 6250, 6250 m
#CRS arguments: +proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 +k=1 +x_0=0 +y_0=0 +a=6378273 +b=6356889.449 +units=m +no_defs 
#sequence starts 2002- 6-19
#early part of time series includes ice on land, latter has it masked out.
#memory issues so calculate monthly average and then stack

#check no gaps in time series of ice files
ice_dts<-as.Date(raadfiles::amsr_daily_files()$date)
date_range <- seq(min(ice_dts), max(ice_dts), by = 1) 
date_range[!date_range %in% ice_dts] # ~300 days missing : period between Oct 2011 and June 2012 (change of satellite?) and some other random dates

#restrict to Dec 2020
ice_dts<-ice_dts[ice_dts < as.Date("2021-01-01")]
ice_dts<-as.data.frame(ice_dts)
ice_dts$MonthYear<-format(ice_dts$ice_dts, "%Y-%m")

#extract daily file and generate monthly summaries (memory issues when tryign to extract entire series at daily resolution)
monthly_ice<-stack()
for(i in 1:length(unique(ice_dts$MonthYear))){
  print(unique(ice_dts$MonthYear)[i])
  mth_dts<-ice_dts$ice_dts[ice_dts$MonthYear %in% unique(ice_dts$MonthYear)[i]]
  daily_ice<-stack()
  for (j in 1:length(mth_dts)){
    temp_ice<-read_amsr_ice(mth_dts[j], time.resolution = "daily", xylim = SO.stereo, setNA=TRUE)
    daily_ice<-stack(daily_ice, temp_ice)
  }
  monthly_ice<-stack(monthly_ice, calc(daily_ice, mean))
}
names(monthly_ice)<-unique(ice_dts$MonthYear)
monthly_ice<-setZ(monthly_ice,unique(ice_dts$MonthYear))

writeRaster(monthly_ice, filename = paste0(VM_path2, "monthly_ice"))


# generate yearly, spring, summer climatologies (mean and sd)
ice_mean<- calc(monthly_ice, mean, na.rm=TRUE)
ice_max<- calc(monthly_ice, max, na.rm=TRUE)
ice_sd<-calc(monthly_ice, sd, na.rm=TRUE)

#remove missing months from spring and summer indices
spring_mth_ind<-spring_mth_ind[! spring_mth_ind %in% c("2011-10", "2011-11") ]
summer_mth_ind<-summer_mth_ind[ ! summer_mth_ind %in% c("2011-21", "2021-01", "2012-02")]

ice_sp_mean<-seas_stats(stack= monthly_ice, season_dates = spring_mth_ind, stat=mean )
ice_sp_max<-seas_stats(stack= monthly_ice, season_dates = spring_mth_ind, stat=max )
ice_sp_sd<-seas_stats(stack= monthly_ice, season_dates = spring_mth_ind, stat=sd )

ice_su_mean<-seas_stats(stack= monthly_ice, season_dates = summer_mth_ind, stat=mean )
ice_su_max<-seas_stats(stack= monthly_ice, season_dates = summer_mth_ind, stat=max )
ice_su_sd<-seas_stats(stack= monthly_ice, season_dates = summer_mth_ind, stat= sd )


#Proportion months ice >85% 
ice_cov<-calc(monthly_ice, function(x) ifelse(x>85, 1, 0))
ice_prop<-sum(ice_cov, na.rm=TRUE)/215

ice_stack<-stack(ice_prop, ice_mean, ice_max, ice_sd, ice_sp_mean, ice_sp_max, ice_sp_sd, ice_su_mean,  ice_su_max, ice_su_sd)
names(ice_stack)<-c("ice_prop", "ice_mean", "ice_max", "ice_sd", "ice_sp_mean", "ice_sp_max", "ice_sp_sd", "ice_su_mean", "ice_su_max", "ice_su_sd")

#mask out land- later layers have land masked
ice_stack<-mask(ice_stack, subset(monthly_ice, 215) )

# slight reprojection
ice_stack_500<-projectRaster(ice_stack, bathy_shelf)

#crop rasters to shelf
ice_stack_500m_shelf<-mask(ice_stack_500, bathy_shelf)

#save rasters
writeRaster(ice_stack, filename=paste0(VM_path2, string.chr, "ice"), overwrite=TRUE)
writeRaster(ice_stack_500, filename=paste0(VM_path2, string.chr, string.res, "ice"), overwrite=TRUE)
writeRaster(ice_stack_500m_shelf, filename=paste0(VM_path2, string.chr, string.res, "shelf_ice"), overwrite=TRUE)

#remove intermediate rasters
rm(monthly_ice, ice_mean, ice_max, ice_sd, ice_sp_mean, ice_sp_max, ice_sp_sd, ice_su_mean, ice_su_max, ice_su_sd)

######################################################
#### 4) 2012 geomorphology from Alix Post ---
#library(sf)
#library(fasterize)
library(rasterDT)
library(rgdal)
#ogrListLayers("/perm_storage/shared_space/BioMAS/environmental_data/Geomorphology.gdb")
ogrListLayers("C:\\Users\\hillna\\OneDrive - University of Tasmania\\UTAS_work\\Projects\\Benthic Diversity ARC\\data_environmental\\raw\\Geomorphology.gdb")
#geomorph<-readOGR("/perm_storage/shared_space/BioMAS/environmental_data/Geomorphology.gdb", layer="AntarcticGeomorphology")
geomorph<-readOGR("C:\\Users\\hillna\\OneDrive - University of Tasmania\\UTAS_work\\Projects\\Benthic Diversity ARC\\data_environmental\\raw\\Geomorphology.gdb", layer="AntarcticGeomorphology")

#geomorph<-st_read( "/perm_storage/shared_space/BioMAS/environmental_data/Geomorphology.gdb", layer="AntarcticGeomorphology")
#geomorph_rast<-fasterizeDT(geomorph, bathy_shelf, field="Feature")
geomorph_rast<-rasterize(geomorph, bathy_shelf, field=3)
plot(geomorph_rast)
writeRaster(geomorph_rast, filename = paste0(VM_path2, string.chr, "geomorphology"))





