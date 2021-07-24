## Ross Sea images: Assign lat/lon to each image:
## comments from Dave Bowden:
# 1.	‘\DTIS_splined track_Jan09’. This contains a navigation log file for each transect, with positions derived after smoothing to correct for noise and splining to yield records at 1 s intervals. SUB1_lon and SUB1_lat are the USBL camera coordinates.
# 2.	‘\tan0802_OFOP_QA_prot_files’. This contains ‘protocol’ log files, which record locations and times of all observations and camera commands during a transect, including the key points of ‘AT THE BOTTOM’, ‘Start video recording’, ‘OFF THE BOTTOM’, ‘Stop video recording’, and ‘photo’. The SUB1 coordinates in these files are the raw data i.e., before the smoothing and splining steps above.
# Using these two sets of files, you can extract position coordinates for individual images. This can be done either by:
# a)	Lining up the first seabed image (the first one that isn’t a deck test shot) with the first ‘photo’ message after the ‘start video recording’ message, then working through in sequence, matching ‘photo’ messages to image files. The system takes a still image every 15 s, so once the first image is matched, you can just assign by elapsed time.
# b)	Hoovering each set of files into a single dataframe and using ‘photo’ times in the protocol files to extract corrected, splined, coordinates from the splined tracks: more accurate positions but more work and perhaps more accuracy than is needed for this application. Also note that date-times in our OFOP recording system are notoriously unreliable … .

#### FIRST: LINE UP IMAGE-INFO WITH TXT-FILES USING TIME-STAMPS
#### PC-TIME VARIES BETWEEN TRANSECTS!!! -> USE UTC-TIME
#### CAN'T USE PROT.TXT FILES BECAUSE NOT ALL IMAGES ARE LOGGED
## jpg-image-time is between 10-30 20 seconds in front of UTC-time in this survey (but 12 hours behind, so use 11:59:40)

# transect 25: min 0s, max 12s 
# transect 40: min 7s, max 
# transect 55: 11s
# transect 65: 12s
# ... later its 30s
## go for 20s as middle ground

###############
############### STATION 274 is a SPECIAL CASE!!! Calculted GPS by interpolation between start and finish
###############

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

sci.dir <-      "C:/Users/jjansen/Desktop/science/"
env.derived <-  paste0(sci.dir,"data_environmental/derived/")

r.path <- "R:/IMAS/Antarctic_Seafloor/Clean_Data_For_Permanent_Storage/TAN0802/TAN0802_1_raw_images_and_metadata/"

ps.path <- paste0(r.path,"images_colourcorrected/")
path.bad.images <- paste0(ps.path,"bad_quality/")

#### specify paths
img.path <- paste0(r.path,"images_original/")
#txt.path <- "D:/ARC_DP_data/TAN0802/tan0802_OFOP_QA_prot_files/"
#gps.txt.path <- "D:/ARC_DP_data/TAN0802/DTIS_splined track_Jan09/"
txt.path <- paste0(r.path,"metadata/DTIS_splined track_Jan09/")

my_data_dir <- "C:/Users/jjansen/Desktop/science/data_environmental/accessed_through_R"

################################################
## INSERT AFTER I'VE RUN IT
load("C:/Users/jjansen/Desktop/science/SouthernOceanBiodiversityMapping/ARC_Data/ReadIn_Circumpolar_DownwardImage_Data_And_Subset_Surveys/TAN0802_dat.Rdata")
dat.old <- dat
t.old <- total.t.length.v
## TBC further below: THIS IS A MESS!!!
################################################

#### read in metadata associated with the images
files <- list.files(img.path,pattern = "*.jpg", full.names = TRUE)
dat.full <- read_exif(files)
dat <- dat.full
dat$CreateDate
dat$CreateDate.date <- ymd(substr(dat$CreateDate,1,10))
dat$CreateDate.time <- substr(dat$CreateDate,12,19)
dat$CreateDate.datetime <- paste0(dat$CreateDate.date," ",substr(dat$CreateDate.time,1,5))

#### read gps data from posi.txt files ####
txt.files <- list.files(txt.path,pattern = "*posi.txt", full.names = TRUE)
txt.files.short <- list.files(txt.path,pattern = "*posi.txt")
txt.files.list <- list()
for(i in 1:length(txt.files)){
  txt.files.list[[i]] <- read.csv(txt.files[i], sep="\t", row.names=NULL)
  #names(txt.files.list[[i]]) <- c(names(txt.files.list[[i]])[-1],NA)
  names(txt.files.list[[i]])[1] <- gsub("X..","",names(txt.files.list[[i]])[1])
  message(i)
  message(txt.files.short[i])
  print(names(txt.files.list[[i]]))
}

#### extract GPS coordinates from posi.txt lining up image-timetamps with UTC-time + 20 seconds ####
dat$transectID <- factor(substr(dat$FileName,9,11))
dat$GPS_lat <- NA
dat$GPS_lon <- NA
dat$GPS_ship_lat <- NA
dat$GPS_ship_lon <- NA
dat$depth <- NA
for(i in 1:length(txt.files)){
  loop.transect.ID <- substr(txt.files.short[i],9,11)
  loop.dat.select <- which(dat$transectID==loop.transect.ID)
  ## remove rows showing the same time-stamps as the row above, and reformat txt-file time-stamps for comparison
  txt.loop.full <- txt.files.list[[i]]
  txt.loop.full$lookup.time <- as.character(hms(txt.loop.full$Time))
  if(any(duplicated(txt.loop.full$lookup.time))){
    txt.loop <- txt.loop.full[-which(duplicated(txt.loop.full$lookup.time)),]
  } else txt.loop <- txt.loop.full
  ## adjusting time-stamp difference by adjusting the timestamp of the image-metadata, and reformat image data time-stamp for comparison
  td <- seconds_to_period(as.numeric(hms(dat$CreateDate.time[loop.dat.select])+hms("11:59:40")))
  loop.dat.time <- as.character(hms(sprintf('%02d:%02d:%02d', td@hour, minute(td), second(td))))
  ## create a lookup object to find gps points in txt-file
  sel.txt <- which(txt.loop$lookup.time%in%loop.dat.time)
  ## create the lookup object to find corresponding images (which images actually have a proper ship-log for time and gps)
  sel.dat <- which(loop.dat.time%in%txt.loop$lookup.time)
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
  print(i)
  message(paste0("transect ",loop.transect.ID))
  message(paste0("# of gps coordinates: ",length(sel.txt),";",length(sel.dat)," images"))
  if(length(sel.txt)!=length(sel.dat)){
    message("something went wrong")
    break
  }
  ## assign gps coordinates
  dat$GPS_lat[loop.dat.select][sel.dat] <- txt.loop$SUB1_Lat[sel.txt]
  dat$GPS_lon[loop.dat.select][sel.dat] <- txt.loop$SUB1_Lon[sel.txt]
  dat$GPS_ship_lat[loop.dat.select][sel.dat] <- txt.loop$SHIP_Lat[sel.txt]
  dat$GPS_ship_lon[loop.dat.select][sel.dat] <- txt.loop$SHIP_Lon[sel.txt]
  dat$depth[loop.dat.select][sel.dat] <- txt.loop$Water_Depth[sel.txt]
}

#### plot all positions of the actual images, and underlay all tracked gps locations to check if anything is off ####
# set_data_roots(my_data_dir)
# r <- readtopo("ibcso")
# r2 <- r
# r2[r2>0] <- NA

dat$GPS_lat[dat$GPS_lat==-999.9] <- NA
dat$GPS_lon[dat$GPS_lon==-999.9] <- NA
dat$GPS_lon[dat$GPS_lon==0] <- NA

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

#####################################################################
##### 2. SUBSET IMAGES FROM IMAGE LOCATIONS (& STORE FILENAMES) #####
#####################################################################

## some transect must have had a drunken captain... we don't want to oversample because the transects run in loops so use distance to start rather than distance along the transect line:
transects.dist.to.start <- c("076","166","169","246","248","285","294")
transects.dist.to.end <- c("141","219","244")

## badly illuminated images have been copied into the folder "TAN0802_bad_images"
#bad_images.path <- "D:/ARC_DP_data/adjusted_TAN0802/bad_quality/"
#bad_images <- list.files(bad_images.path)
# save(bad_images, file="TAN0802_filenames_badimages.RData")
load("C:/Users/jjansen/Desktop/science/DP190101858_MappingAntarcticSeafloorBiodiversity/TAN0802_filenames_badimages.RData")
# ## good images have been copied into the folder "TAN0802_good_images"
# good_images.path <- "D:/ARC_DP_data/TAN0802/TAN0802_good_images/"
# good_images <- list.files(good_images.path)

## remove bad images first
usable.images.raw <- dat[which(dat$FileName%!in%bad_images),]
## and also those without gps data
usable.images <- usable.images.raw[!is.na(usable.images.raw$GPS_lat+usable.images.raw$GPS_lon),]

dat <- data.frame(usable.images[,which(names(usable.images)%in%c("FileName","GPS_lon", "GPS_lat","GPS_ship_lon","GPS_ship_lat","depth","CreateDate.datetime"))])[,c(4,3,1,6,5,7,2)]
dat$transectID <- as.factor(substr(dat$FileName,9,11))



################################################
## SECOND INSERT AFTER I'VE RUN IT:
dat.old$depth <- dat$depth
dat.old$time <- dat$CreateDate.datetime
dat <- dat.old
total.t.length.v <- t.old
save(dat,total.t.length.v, t.l, file="C:/Users/jjansen/Desktop/science/data_biological/TAN0802_dat.Rdata")
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
  if(levels(dat$transectID)[i]%in%transects.dist.to.start){
    samp <- quasiSamp(length(dat.subset.v)*5,dimension=1,potential.sites=dat$dist.to.start[dat.subset.v])
    plot(dat$dist.from.start[dat.subset.v],dat$dist.to.start[dat.subset.v], main=paste0(levels(dat$transectID)[i]," from start vs to start"))
  }else if(levels(dat$transectID)[i]%in%transects.dist.to.end){
    samp <- quasiSamp(length(dat.subset.v)*5,dimension=1,potential.sites=dat$dist.to.end[dat.subset.v])
    plot(dat$dist.from.start[dat.subset.v],dat$dist.to.end[dat.subset.v], main=paste0(levels(dat$transectID)[i]," from start vs to end"))
  } else samp <- quasiSamp(length(dat.subset.v)*5,dimension=1,potential.sites=dat$dist.from.start[dat.subset.v])
  
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
#save(dat,total.t.length.v, t.l, file="C:/Users/jjansen/Desktop/science/data_biological/TAN0802_dat.Rdata")
#load(file="C:/Users/jjansen/Desktop/science/data_biological/TAN0802_dat.Rdata")

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
  if(levels(dat$transectID)[i]%in%transects.dist.to.start){
    sel <- order(dat$image.select[subset.v])[1:t.images.full[i]]
  }else if(levels(dat$transectID)[i]%in%transects.dist.to.end){
    sel <- order(dat$image.select[subset.v])[1:t.images.full[i]]
  }else sel <- order(dat$image.select[subset.v])[1:t.images[i]]
  ## create vector of selected images
  if(i==1){dat.sel <<-  subset.v[sel]
  }else dat.sel <-  c(dat.sel,subset.v[sel])
  png(filename=paste0("TAN0802_subset_transect_",unique(dat$transectID[subset.v]),".png"),width=1000, height=1000)
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
selected.filenames.renamed <- paste0(substr(dat$FileName[dat.sel],1,12),sprintf("%04d",dat$image.select[dat.sel]),"__",selected.filenames)
selected.filenames.folders <- substr(selected.filenames,1,7)

## copy files into Annotation folder
img.path.origin <- paste0(img.path,selected.filenames)
img.path.destin <- "C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/Stills/Annotation_images/TAN0802/"
file.copy(img.path.origin,img.path.destin)

## write list of filenames into cropped Annotation folder for upload
img.path.destin_crop <- "C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/Stills/Annotation_images_cropped/TAN0802/"
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

#### IN TAN0802 NO CROPPING DONE, EDGES LOOK FINE

# # # @@@ !!! @@@
# # # crop images, do this in the terminal using the following code (after navigating into the "Stills" folder):
# # cd C:
# # magick mogrify -path Annotation_images_cropped\TAN0802\ -gravity Center -crop 90% +repage Annotation_images\TAN0802\*jpg
# # # @@@ !!! @@@
