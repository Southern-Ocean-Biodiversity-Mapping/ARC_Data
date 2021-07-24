#########################################################################################################################
## Ross Sea images: Assigning lat/lon to each image, removing bad images and then subsetting the full dataset
#########################################################################################################################

## TAN1802 ##

## comments from Dave Bowden:
#As noted earlier, the txt-files are raw logs. For matching images to position coordinates, using the time-stamp in the jpeg files matched to times in the *_posi files is likely to be the more reliable way of working. For reference/orientation, the camera vehicle seabed position is given by ‘SUB1_Lon’ and ‘SUB1_Lat’ and the three types of text file contain these data:
# *_prot: summary metadata for deployment, all manually recorded observations plus camera instruction messages (including ‘photo taken’), selected nav data [This file designed for people to read as summary of deployment].
# *_obser: all manually recorded observations plus camera instruction messages, plus USBL beacon messages [This file designed for replotting transect and observations in OFOP or GIS]
# *_posi: all navigational data from before camera vehicle goes in the water to when it returns to deck.

## FOR SOME TRANSECTS THE CAMERA TIMESTAMP IS TOTALLY WRONG
## 092 094 096 097 098 102 105 106 107 108 121 122 125 126 127
weird_transect_dates_a <- c("092","094","096","097","102","105","121","122","125")
weird_transect_dates_b <- c("098","106","107","108","126","127")
## dat$CreateDate plus 5 years, 8 months, 1 day, 1 hour, 23 minutes, 32sec TO PC_time
## plus minus a few seconds

## OTHER TRANSECTS ARE OFF by PC_TIME plus 1-5 SEC

## USBL data for transect 122 is bad, using ship coordinates instead
## first 53 images of transect 170 have no valid GPS data
## images images 162 to 223 on transect 191 have no USBL data, using interpolation to fill the gap

###############################################################################
#### libraries and functions
library(exifr)
library(dplyr)
library(leaflet)
library(lubridate)
library(SOmap)
library(raadtools)
library(spatialEco)
library(raster)
library(sp)
library(blueant)
library(geosphere)
library(MBHdesign)
'%!in%' <- function(x,y)!('%in%'(x,y))

#### specify paths
img.path <- "D:/ARC_DP_data/a_RawData_DirectFromContributors/TAN1802/"
txt.path <- "D:/ARC_DP_data/a_RawData_DirectFromContributors/TAN1802/"
gps.txt.path <- "D:/ARC_DP_data/a_RawData_DirectFromContributors/TAN1802/"

my_data_dir <- "C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_environmental/accessed_through_R"

###############################################################################
#### read gps data from posi.txt files ####
txt.files <- list.files(txt.path,pattern = "*posi.txt", full.names = TRUE, recursive=TRUE)
txt.files.short <- list.files(txt.path,pattern = "*posi.txt", recursive=TRUE)
txt.files.list <- list()
for(i in 1:length(txt.files)){
  txt.files.list[[i]] <- read.csv(txt.files[i], sep="\t", row.names=NULL)
  names(txt.files.list[[i]]) <- c(names(txt.files.list[[i]])[-1],NA)
  names(txt.files.list[[i]])[1] <- substr(names(txt.files.list[[i]])[1],3,6)
  message(i)
  print(names(txt.files.list[[i]]))
  txt.files.list[[i]]$SUB1_Lat.interpolated <- NA
  txt.files.list[[i]]$SUB1_Lon.interpolated <- NA
  ## some images have no proper USBL GPS. Interpolate between images that have USBL gps-points
  ## use ship gps to infer image locations on transect 191 (images 162 to 223)
  if(i==24){
    txt.files.list[[i]]$SUB1_Lat[3549:4536] <- 0
    txt.files.list[[i]]$SUB1_Lon[3549:4536] <- 0
  }
  if(any(txt.files.list[[i]]$SUB1_Lat==0|txt.files.list[[i]]$SUB1_Lat==-999.9|txt.files.list[[i]]$SUB1_Lon==0|txt.files.list[[i]]$SUB1_Lon==-999.9)){
    message("interpolating missing GPS values")
    no.gps <- which(txt.files.list[[i]]$SUB1_Lat==0|txt.files.list[[i]]$SUB1_Lat==-999.9|txt.files.list[[i]]$SUB1_Lon==0|txt.files.list[[i]]$SUB1_Lon==-999.9)
    #### generate vector of numbers to inform which values need replacing and which values to use for the interpolation
    ### all these gps-points are missing in a sequence:
    seq.miss.gps <- sort(unique(c(which(diff(no.gps)==1),which(diff(no.gps)==1)+1)))
    ## fist missing point in each sequence:
    seq.start.pt <- c(1,which(diff(no.gps[seq.miss.gps])!=1)+1)
    ## last missing point in each sequence:
    seq.end.pt <- c(which(diff(no.gps[seq.miss.gps])!=1),length(seq.miss.gps))
    ## fill values by sequence, store also separate for future reference
    if(length(seq.miss.gps)>0){
      for(k in 1:length(seq.start.pt)){
        seq.v <- no.gps[seq.miss.gps[seq.start.pt[k]:seq.end.pt[k]]]
        a <- seq.v[1]-1
        b <- seq.v[length(seq.v)]+1
        if(a==0|b>nrow(txt.files.list[[i]])) next
        l <- length(seq.v)
        txt.files.list[[i]]$SUB1_Lat[seq.v] <- seq(txt.files.list[[i]]$SUB1_Lat[a],txt.files.list[[i]]$SUB1_Lat[b],length.out=l)
        txt.files.list[[i]]$SUB1_Lon[seq.v] <- seq(txt.files.list[[i]]$SUB1_Lon[a],txt.files.list[[i]]$SUB1_Lon[b],length.out=l)
        txt.files.list[[i]]$SUB1_Lat.interpolated[seq.v] <- seq(txt.files.list[[i]]$SUB1_Lat[a],txt.files.list[[i]]$SUB1_Lat[b],length.out=l)
        txt.files.list[[i]]$SUB1_Lon.interpolated[seq.v] <- seq(txt.files.list[[i]]$SUB1_Lon[a],txt.files.list[[i]]$SUB1_Lon[b],length.out=l)
      }
    }
    message("sequentially missing GPS points interpolated")
  }
}

#### read image-metadata and extract time-stamp ####
files <- list.files(img.path,pattern = "*.jpg", full.names = TRUE, recursive=TRUE)
dat.full <- read_exif(files)
dat <- dat.full
dat$CreateDate
dat$CreateDate.date <- ymd(substr(dat$CreateDate,1,10))
dat$CreateDate.time <- substr(dat$CreateDate,12,19)
dat$CreateDate.datetime <- paste0(dat$CreateDate.date," ",substr(dat$CreateDate.time,1,5))
dat$FileName_uniform <- dat$FileName
dat$FileName_uniform <- gsub("_94_","_094_",dat$FileName_uniform)
dat$FileName_uniform <- gsub("_Stn_","_",dat$FileName_uniform)
dat$transectID <- factor(substr(dat$FileName_uniform,9,11))

#### extract GPS coordinates from posi.txt lining up image-timetamps with UTC-time + 20 seconds (as calculated above) ####
dat$GPS_lat <- NA
dat$GPS_lon <- NA
dat$GPS_lat.interpolated <- NA
dat$GPS_lon.interpolated <- NA
dat$GPS_ship_lat <- NA
dat$GPS_ship_lon <- NA
dat$SUB1_Altitude <- NA
dat$SUB1_USBL_Depth <- NA
for(i in 1:length(txt.files)){
  loop.transect.ID <- substr(txt.files.short[i],9,11)
  loop.dat.select <- which(dat$transectID==loop.transect.ID)
  ## remove rows showing the same time-stamps as the row above, and reformat txt-file time-stamps for comparison
  txt.loop.full <- txt.files.list[[i]]
  txt.loop.full$lookup.time <- as.character(hms(substr(txt.loop.full$PC_Time,11,19)))
  if(any(duplicated(txt.loop.full$lookup.time))){
    txt.loop <- txt.loop.full[-which(duplicated(txt.loop.full$lookup.time)),]
  } else txt.loop <- txt.loop.full
  txt.loop$SUB1_Lat.interpolated <- NA
  txt.loop$SUB1_Lon.interpolated <- NA
  message(paste0("loop ",i))
  message(paste0("transect ",loop.transect.ID))
  if(loop.transect.ID=="113") next
  print(txt.loop.full$lookup.time[seq(1,length(txt.loop.full$lookup.time), by=100)])
  print(hms(dat$CreateDate.time[loop.dat.select][seq(1,length(dat$CreateDate.time[loop.dat.select]), by=10)]))
  ## adjusting time-stamp difference by adjusting the timestamp of the image-metadata, and reformat image data time-stamp for comparison
  if(loop.transect.ID%in%weird_transect_dates_a){
    td <- seconds_to_period(as.numeric(hms(dat$CreateDate.time[loop.dat.select])+hms("13:23:32")))
  }else if(loop.transect.ID%in%weird_transect_dates_b){
    td <- seconds_to_period(as.numeric(hms(dat$CreateDate.time[loop.dat.select])+hms("01:23:32")))
  }else if(loop.transect.ID%in%c("170","179","180","183","184","191","193","195","196","213")){
    td <- seconds_to_period(as.numeric(hms(dat$CreateDate.time[loop.dat.select])+hms("12:00:01")))
  }else{
    td <- seconds_to_period(as.numeric(hms(dat$CreateDate.time[loop.dat.select])+hms("00:00:01")))
  }
  
  ## make everything comparable
  loop.dat.time <- as.character(hms(sprintf('%02d:%02d:%02d', td@hour, minute(td), second(td))))
  tt <- hms(txt.loop$lookup.time)
  txt.loop.time <- as.character(hms(sprintf('%02d:%02d:%02d', tt@hour, minute(tt), second(tt))))
  
  ## sometimes we have a mixup around midnight when 0H is just being dropped
  if(loop.transect.ID%in%c("106","207","213")){
    txt.loop.time <- gsub("12H ","",txt.loop.time)
  }
  
  ## create a lookup object to find gps points in txt-file
  sel.txt <- which(txt.loop.time%in%loop.dat.time)
  ## create the lookup object to find corresponding images (which images actually have a proper ship-log for time and gps)
  sel.dat <- which(loop.dat.time%in%txt.loop.time)
  if(any(duplicated(loop.dat.time[sel.dat]))){
    sel.dat <- sel.dat[duplicated(loop.dat.time[sel.dat])==FALSE]
  }
  ## some images don't have a gps point at that exact second, so add a second and see if there a gps coordinate available
  sel.dat.unsuccessful <- which(loop.dat.time%!in%txt.loop$lookup.time)
  if(length(sel.dat.unsuccessful>0)){
    loop.dat.time2 <- as.character(hms(hms(loop.dat.time[sel.dat.unsuccessful])+1))
    added.txt <- which(txt.loop$lookup.time%in%loop.dat.time2)
    sel.txt <- c(sel.txt,added.txt)
    sel.txt <- sel.txt[order(sel.txt)]
    sel.dat.unsuccessful.fixed <- sel.dat.unsuccessful[which(loop.dat.time2%in%txt.loop$lookup.time)]
    sel.dat <- c(sel.dat,sel.dat.unsuccessful.fixed)
    sel.dat <- sel.dat[order(sel.dat)]
  }
  ## some images have no proper USBL GPS. Interpolate between images that have USBL gps-points
  message(paste0("images in folder: ",length(loop.dat.time)))
  message(paste0("# of gps coordinates: ",length(sel.txt),";",length(sel.dat)," images"))
  if(length(sel.txt)!=length(sel.dat)){
    message("something went wrong !!!!!!")
    break
  }
  message("                      ")
  message("                      ")  
  ## assign gps coordinates
  dat$GPS_lat[loop.dat.select][sel.dat] <- txt.loop$SUB1_Lat[sel.txt]
  dat$GPS_lon[loop.dat.select][sel.dat] <- txt.loop$SUB1_Lon[sel.txt]
  dat$GPS_lat.interpolated[loop.dat.select][sel.dat] <- txt.loop$SUB1_Lat.interpolated[sel.txt]
  dat$GPS_lon.interpolated[loop.dat.select][sel.dat] <- txt.loop$SUB1_Lon.interpolated[sel.txt]
  dat$GPS_ship_lat[loop.dat.select][sel.dat] <- txt.loop$SHIP_Lat[sel.txt]
  dat$GPS_ship_lon[loop.dat.select][sel.dat] <- txt.loop$SHIP_Lon[sel.txt]
  ##
  dat$SUB1_Altitude[loop.dat.select][sel.dat] <- txt.loop$SUB1_Altitude[sel.txt]
  dat$SUB1_USBL_Depth[loop.dat.select][sel.dat] <- txt.loop$SUB1_USBL_Depth[sel.txt]
}

#### plot all positions of the actual images, and underlay all tracked gps locations to check if anything is off ####
# set_data_roots(my_data_dir)
# r <- readtopo("ibcso")
# r2 <- r
# r2[r2>0] <- NA

dat$GPS_lat[dat$GPS_lat==-999.9] <- NA
dat$GPS_lon[dat$GPS_lon==-999.9] <- NA
dat$GPS_lon[dat$GPS_lon==0] <- NA

## correct coordinates for transect 122:
dat$GPS_lat[dat$transectID=="122"] <- dat$GPS_ship_lat[dat$transectID=="122"]
dat$GPS_lon[dat$transectID=="122"] <- dat$GPS_ship_lon[dat$transectID=="122"]
## remove first 53 images of transect 170:
dat$GPS_lat[dat$transectID=="170"][1:53] <- NA
dat$GPS_lon[dat$transectID=="170"][1:53] <- NA

## remove images without gps data
dat2 <- dat[-which(is.na(dat$GPS_lat+dat$GPS_lon)),]

# ## overview
# spatial.dat <- data.frame(cbind(dat2$GPS_lon, dat2$GPS_lat))
# names(spatial.dat) <- c("Longitude","Latitude")
# coordinates(spatial.dat) <- c("Longitude","Latitude")
# proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
# polar.dat <- spTransform(spatial.dat, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# 
# spatial.dat.ship <- data.frame(cbind(dat2$GPS_ship_lon, dat2$GPS_ship_lat))[!is.na(dat2$GPS_ship_lat+dat2$GPS_ship_lon),]
# names(spatial.dat.ship) <- c("Longitude","Latitude")
# coordinates(spatial.dat.ship) <- c("Longitude","Latitude")
# proj4string(spatial.dat.ship) <- CRS("+proj=longlat +datum=WGS84")
# polar.dat.ship <- spTransform(spatial.dat.ship, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# 
# plot(r2, xlim=c(-600000,500000), ylim=c(-2600000,-1200000))
# points(polar.dat)
# #points(polar.dat.ship, col="red")
# text(polar.dat,labels=substr(dat2$FileName,9,11), adj=c(2,2))
# 
# ## each transect individually
# for(i in 1:length(txt.files)){
#   loop.transect.ID <- substr(txt.files.short[i],9,11)
#   loop.dat.select <- which(dat2$transectID==loop.transect.ID)
#   loop.dat <- dat2[loop.dat.select,]
#   message(i)
#   print(paste0("transect ",loop.transect.ID))
#   #print(loop.dat$GPS_lon)
#   if(nrow(loop.dat)==0) next
#   ## transform coordinates
#   spatial.dat <- data.frame(cbind(loop.dat$GPS_lon, loop.dat$GPS_lat))
#   names(spatial.dat) <- c("Longitude","Latitude")
#   coordinates(spatial.dat) <- c("Longitude","Latitude")
#   proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
#   polar.dat <- spTransform(spatial.dat, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
#   
#   spatial.dat.ship <- data.frame(cbind(loop.dat$GPS_ship_lon, loop.dat$GPS_ship_lat))[!is.na(loop.dat$GPS_ship_lat+loop.dat$GPS_ship_lon),]
#   names(spatial.dat.ship) <- c("Longitude","Latitude")
#   coordinates(spatial.dat.ship) <- c("Longitude","Latitude")
#   proj4string(spatial.dat.ship) <- CRS("+proj=longlat +datum=WGS84")
#   polar.dat.ship <- spTransform(spatial.dat.ship, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
#   
#   ## plot
#   plot(polar.dat, main=paste0("transect ",loop.transect.ID))
#   points(polar.dat.ship, col="blue")
#   scalebar(200, type="bar")
# }

dat3 <- dat2
## list of all filenames with gps data
#write.table(dat3$FileName, file="TAN1802_filenames_allfileswithgps.txt", eol=",", col.names=FALSE, row.names=FALSE)

# ## check quality of images that have interpolated gps points
# interpolated_images <- t(read.table("TAN1901_filenames_allfileswithgps.txt", sep=","))
# noninterpolated_images <- t(read.table("TAN1901_filenames_allfileswithgps_nointerpolation.txt", sep=","))
# new.images <- interpolated_images[which(interpolated_images%!in%noninterpolated_images)]
# file.copy(dat3$SourceFile[which(dat3$FileName%in%new.images)], "D:/ARC_DP_data/Temp_Folder_imageswithinterpolatedGPS/")


#####################################################################
##### 2. SUBSET IMAGES FROM IMAGE LOCATIONS (& STORE FILENAMES) #####
#####################################################################

## dat3 already has all images without gps removed

# ## now remove images that are taken too close or too ar off the seafloor and copy them into a different folder to visually check
# plot(dat3$SUB1_Altitude,ylim=c(0,50), col=dat3$transectID)
# length(dat3$FileName[dat3$SUB1_Altitude>35|dat3$SUB1_Altitude<10])
# path.bad.images <- "D:/ARC_DP_data/adjusted_TAN1802/bad_quality/"
#bad_images <- list.files(path.bad.images)
# save(bad_images, file="TAN1802_filenames_badimages.RData")
load("C:/Users/jjansen/Desktop/science/DP190101858_MappingAntarcticSeafloorBiodiversity/TAN1802_filenames_badimages.RData")
# file.copy(dat3$SourceFile[dat3$SUB1_Altitude>35|dat3$SUB1_Altitude<10], path.bad.images)
# file.copy(dat3$SourceFile[dat3$SUB1_Altitude<=35&dat3$SUB1_Altitude>=10], path.good.images)
# dat3$SUB1_Altitude[dat3$FileName=="tan1901_164_225.jpg"]
#file.copy(dat3$SourceFile[dat3$SUB1_Altitude>35|dat3$SUB1_Altitude<10], "D:/ARC_DP_data/Temp_Folder_imageswithoutdistance/")
## remove bad images first
usable.images.raw <- dat3[which(dat3$FileName%!in%bad_images),]
## and also those without gps data
usable.images <- usable.images.raw[!is.na(usable.images.raw$GPS_lat+usable.images.raw$GPS_lon),]

## all transects are straight, but some gps variation: use distance to start point to calculate subset
dat <- data.frame(usable.images[,which(names(usable.images)%in%c("FileName","GPS_lon", "GPS_lat","GPS_ship_lon","GPS_ship_lat","SUB1_Altitude","SUB1_USBL_Depth","CreateDate.datetime"))])[,c(4,3,1,6,5,7,8,2)]
dat$FileName_uniform <- dat$FileName
dat$FileName_uniform <- gsub("_94_","_094_",dat$FileName_uniform)
dat$FileName_uniform <- gsub("_Stn_","_",dat$FileName_uniform)
dat$transectID <- factor(substr(dat$FileName_uniform,9,11))


################################################
## INSERT AFTER I'VE RUN IT:
dat.new <- dat
load("C:/Users/jjansen/Desktop/science/SouthernOceanBiodiversityMapping/ARC_Data/ReadIn_Circumpolar_DownwardImage_Data_And_Subset_Surveys/TAN1802_dat.Rdata")
dim(dat.new)
dim(dat)
all(dat$FileName%in%dat.new$FileName)

dat$depth <- dat.new$SUB1_USBL_Depth
dat$time <- dat.new$CreateDate.datetime
save(dat,total.t.length.v, t.l, file="C:/Users/jjansen/Desktop/science/data_biological/TAN1802_dat.Rdata")
################################################



## calculate transect length, define how many images to select and give images a random number
dat$image.select <- NA
total.t.length.v <- NA
t.l <- NA
dat$dist.from.start <- NA
dat$dist.to.start <- NA
dat$dist.to.end <- NA
dat$mean.dist.to.10.nearest.images <- NA
dat$dist.from.center <- NA
samp.v.list <- list()
center.idx.v <- NA
for(i in 1:length(levels(dat$transectID))){
  print(levels(dat$transectID)[i])
  ## all data except for badly illuminated images
  subset.v <- which(dat$transectID==levels(dat$transectID)[i])
  dat.subset.v <- subset.v
  
  t.counts <- length(dat.subset.v)
  ## select first and last coordinate to calculate transect distance
  v.first <- dat.subset.v[1]
  v.last <- tail(dat.subset.v, n=1)
  t.l[i] <- max(distm(dat[subset.v,1:2], fun=distVincentyEllipsoid))
  
  # total.t.length.v[i] <- t.l
  ## calculate distance of each image to the start of the transect along the transect line
  dat$dist.from.start[v.first] <- 0
  dat$dist.to.start[v.first] <- 0
  dat$dist.to.end[v.first] <- t.l[i]
  for(k in 2:t.counts){
    dist.to.previous.image <- distVincentyEllipsoid(dat[dat.subset.v[k-1],1:2], dat[dat.subset.v[k],1:2])
    dat$dist.from.start[dat.subset.v[k]] <- dat$dist.from.start[dat.subset.v[k-1]]+dist.to.previous.image
    dat$dist.to.start[dat.subset.v[k]] <- distVincentyEllipsoid(dat[v.first,1:2], dat[dat.subset.v[k],1:2])
    dat$dist.to.end[dat.subset.v[k]] <- distVincentyEllipsoid(dat[v.last,1:2], dat[dat.subset.v[k],1:2])
    }
  total.t.length.v[i] <- dat$dist.from.start[dat.subset.v[t.counts]]
  ## choose a spatial random subset of images in each transect using quasiSamp on the distance from the start of the transect
  set.seed(42)
  ## use distance to start:
  samp <- quasiSamp(length(dat.subset.v)*5,dimension=1,potential.sites=dat$dist.to.start[dat.subset.v])
  # if(levels(dat$transectID)[i]%in%transects.dist.to.start){
  #  samp <- quasiSamp(length(dat.subset.v)*5,dimension=1,potential.sites=dat$dist.to.start[dat.subset.v])
  #  plot(dat$dist.from.start[dat.subset.v],dat$dist.to.start[dat.subset.v], main=paste0(levels(dat$transectID)[i]," from start vs to start"))
  # }else if(levels(dat$transectID)[i]%in%transects.dist.to.end){
  #  samp <- quasiSamp(length(dat.subset.v)*5,dimension=1,potential.sites=dat$dist.to.end[dat.subset.v])
  #  plot(dat$dist.from.start[dat.subset.v],dat$dist.to.end[dat.subset.v], main=paste0(levels(dat$transectID)[i]," from start vs to end"))
  # } else samp <- quasiSamp(length(dat.subset.v)*5,dimension=1,potential.sites=dat$dist.from.start[dat.subset.v])

  samp.v <- samp$ID[which(duplicated(samp)==FALSE)]
  print(paste0("first double selection for image # ", which(duplicated(samp)==TRUE)[1]))
  print(paste0("# of images in first batch: ", ceiling(total.t.length.v[i]/50)))
  #samp.v[(length(samp.v)+1):length(subset.v)] <- 9998
  dat$image.select[dat.subset.v[samp.v]] <-  1:length(samp.v)
  dat$image.select[dat.subset.v[-samp.v]] <-  999
  # ## choose a random subset of images in each transect by giving them random numbers starting from 1
  #dat$image.select[dat.subset.v] <- sample(1:length(dat.subset.v),length(dat.subset.v))
  samp.v.list[[i]] <- samp.v
}
## 
dat$image.select[is.na(dat$image.select)] <- 9999

## SAVE OUTPUT FOR FUTURE REFENCE (i.e. start here to add more images to the analysis)
#save(dat,total.t.length.v, t.l, file="C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/TAN1802_dat.Rdata")
#load(file="C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/TAN1802_dat.Rdata")

t.images <- ceiling(total.t.length.v/100)
t.images.full <- ceiling(t.l/100)
## 500 images in total
#t.images <- ceiling((total.t.length.v/total.t.length)*500)
## 

spatial.dat <- data.frame(dat[,1:2])#[!is.na(rowSums(dat[,4:5])),]
coordinates(spatial.dat) <- c("GPS_lon","GPS_lat")
proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
polar.dat <- spTransform(spatial.dat, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

ship.subset <- which(!is.na(rowSums(dat[,4:5])))
spatial.dat.ship <- data.frame(dat[,4:5])[ship.subset,]
names(spatial.dat.ship) <- c("Longitude","Latitude")
coordinates(spatial.dat.ship) <- c("Longitude","Latitude")
proj4string(spatial.dat.ship) <- CRS("+proj=longlat +datum=WGS84")
polar.dat.ship <- spTransform(spatial.dat.ship, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

#par(mfrow=c(4,4), mar=c(5,4,2,1))
for(i in 1:length(levels(dat$transectID))){
  subset.v <- which(dat$transectID==levels(dat$transectID)[i])
  ship.subset.v <- which(dat$transectID[ship.subset]==levels(dat$transectID)[i])
  ## select the correct number of images for each transect
  # if(levels(dat$transectID)[i]%in%transects.dist.to.start){
  #   sel <- order(dat$image.select[subset.v])[1:t.images.full[i]]
  # }else if(levels(dat$transectID)[i]%in%transects.dist.to.end){
  #   sel <- order(dat$image.select[subset.v])[1:t.images.full[i]]
  # }else 
    sel <- order(dat$image.select[subset.v])[1:t.images[i]]
  ## create vector of selected images
  if(i==1){dat.sel <<-  subset.v[sel]
  }else dat.sel <-  c(dat.sel,subset.v[sel])
  png(filename=paste0("TAN1901_subset_transect_",unique(dat$transectID[subset.v]),".png"),width=1000, height=1000)
  plot(polar.dat[subset.v], main=levels(dat$transectID)[i])
  text(polar.dat[subset.v][1], labels=subset.v[1], adj=1)
  points(polar.dat[subset.v][1], col="green")
  points(polar.dat[tail(subset.v)[6]], col="deeppink")
  #points(polar.dat.ship[ship.subset.v], col="blue", pch=15)
  a <- round(total.t.length.v[i])
  b <- round(t.l[i])
  mtext(paste0(round(a/b,2),"  ---  ",a,"m along transect; ",b,"m direct from start to finish"))
  points(polar.dat[subset.v[sel]], col="red", pch=16,cex=2)
  points(polar.dat[subset.v[sel[1:ceiling(length(sel)/2)]]], col="deepskyblue", pch=16)
  legend("topright",legend=c("%-cover","counts"),col=c("red","deepskyblue"),pch=16,pt.cex=c(2,1),bty="n")
  #points(polar.dat[subset.v[center.idx.v[i]],1:2], col="green", pch=16)
  scalebar(100, type="bar")
  dev.off()
  # plot(dat[subset.v,1:2])
  # points(dat[subset.v[sel],1:2], col="red", pch=15)
  # points(dat[subset.v[center.idx.v[i]],1:2], col="green", pch=16)
}

dat$FileName[dat.sel]
## filenames
selected.filenames <- dat$FileName[dat.sel]
## create new filenames and identify folders where to find them (might be different for each survey)
selected.filenames.renamed <- paste0("TAN1802_",dat$transectID[dat.sel],"_",sprintf("%04d",dat$image.select[dat.sel]),"__",selected.filenames)

## copy files into Annotation folder
img.path.origin <- paste0(img.path,selected.filenames)
img.path.destin <- "C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/Stills/Annotation_images/TAN1802/"
file.copy(img.path.origin,img.path.destin)

## write list of filenames into cropped Annotation folder for upload
img.path.destin_crop <- "C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/Stills/Annotation_images_cropped/TAN1802/"
write.table(selected.filenames.renamed, paste0(img.path.destin_crop,"filenames.txt"), eol=",", col.names=FALSE, row.names=FALSE)

filenames.in.folder_original <- paste0(img.path.destin,selected.filenames)
filenames.in.folder_changed <- paste0(img.path.destin,selected.filenames.renamed)
file.rename(filenames.in.folder_original,filenames.in.folder_changed)


# a <- list.files(img.path,pattern = "*.jpg")
# img.path.origin.good <- paste0(img.path,a[a%!in%bad_images])
# img.path.destin.good <- paste0(img.path,"TAN0802_good_images/")
# file.copy(img.path.origin.good,img.path.destin.good)
##

########################################################
##### 3. CROP SELECTED IMAGES AND STORE SEPARATELY #####
########################################################

## TAN 1901: DON'T CROP BUT JUST COPY

## @@@ !!! @@@ 
## crop images, do this in the terminal using the following code (after navigating into the "Stills" folder):
# cd C:
# magick mogrify -path Annotation_images_cropped\TAN1901\ -gravity Center -crop 90% +repage Annotation_images\TAN1901\*jpg
## @@@ !!! @@@ 
















