#########################################################################################################################
## Sabrine Coast images: Assigning lat/lon to each image, removing bad images and then subsetting the full dataset
#########################################################################################################################

## NBP1402 ##

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
img.path <- "E:/ARC_DP_data/a_RawData_DirectFromContributors/NBP14_02_Post2016/"
txt.path <- "E:/ARC_DP_data/a_RawData_DirectFromContributors/NBP14_02_Post2016/"
gps.txt.path <- "E:/ARC_DP_data/a_RawData_DirectFromContributors/NBP14_02_Post2016/Totten_all/"

r.path <- "R:/IMAS/Antarctic_Seafloor/Clean_Data_For_Permanent_Storage/"
nbp.path <- paste0(r.path,"NBP1402/NBP1402_1_raw_images_and_metadata/images_colourcorrected/")

my_data_dir <- "C:/Users/jjansen/Desktop/science/data_environmental/raw/accessed_through_R"

###############################################################################
##
dat <- read.csv(paste0(txt.path,"NBP14-02_camera_results_metadata_Post et al 2017.csv"))
names(dat) <- c("FileName","Survey","Date","Time","transectID","Depth","GPS_lat","GPS_lon")

dat$FileName <- as.character(dat$FileName)
add.yoyo <- which(grepl("yoyo",dat$FileName)==FALSE)
dat$FileName[add.yoyo] <- paste0("yoyo",as.character(dat$FileName[add.yoyo]))
dat$FileName <- paste0(dat$FileName,".jpg")
## fixing missing timestamps in transect yoyo36
dat$Time[365:367] <- c("16:55","16:56","16:57")
dat$Time[369:372] <- c("16:57","16:57","16:58","16:58")
dat$Time[374:377] <- c("16:58","16:59","16:59","17:00")
dat$Time[379:382] <- c("17:00","17:01","17:01","17:02")
dat$Time[384:387] <- c("17:02","17:03","17:03","17:04")
dat$Time[389:392] <- c("17:04","17:05","17:05","17:06")
dat$Time[394:397] <- c("17:06","17:07","17:07","17:08")
dat$Time[399:402] <- c("17:08","17:09","17:09","17:10")
dat$Time[404:407] <- c("17:10","17:11","17:11","17:12")
dat$Time[409:412] <- c("17:12","17:12","17:13","17:13")
dat$Time[414:417] <- c("17:13","17:14","17:14","17:15")
dat$Time[419:422] <- c("17:15","17:16","17:16","17:17")
dat$Time[424:427] <- c("17:17","17:17","17:18","17:18")
dat$Time[429:432] <- c("17:18","17:18","17:19","17:19")
dat$Time[434:437] <- c("17:19","17:20","17:20","17:21")
dat$Time[439:446] <- c("17:21","17:22","17:22","17:22","17:23","17:23","17:23","17:24")
dat$Time[448:449] <- c("17:24","17:25")
dat$time <- dmy_hm(paste0(dat$Date," ",dat$Time))

#### plot all positions of the actual images, and underlay all tracked gps locations to check if anything is off ####
# set_data_roots(my_data_dir)
# r <- readtopo("ibcso")
# r2 <- r
# r2[r2>0] <- NA

## remove images without gps data
dat2 <- dat[-which(is.na(dat$GPS_lat+dat$GPS_lon)),]
dat2$transectID <- factor(dat2$transectID)

#### some images clearly have wrong GPS data
## transect 35: 133:146 (interpolate position using images on transect)
bad.gps.select <- which(dat2$transectID==35)[133:146]
last.good.gps <- which(dat2$transectID==35)[132]
first.good.gps <- which(dat2$transectID==35)[147]
dat2$GPS_lat[bad.gps.select] <- seq(dat2$GPS_lat[last.good.gps],dat2$GPS_lat[first.good.gps],length.out=length(bad.gps.select)+2)[-c(1,length(bad.gps.select)+2)]
dat2$GPS_lon[bad.gps.select] <- seq(dat2$GPS_lon[last.good.gps],dat2$GPS_lon[first.good.gps],length.out=length(bad.gps.select)+2)[-c(1,length(bad.gps.select)+2)]
## transect 36: 76:86 (remove images (nothing to interpolate position from))
bad.gps.select <- which(dat2$transectID==36)[76:86]
dat2$GPS_lat[bad.gps.select] <- NA
dat2$GPS_lon[bad.gps.select] <- NA
## transect 40: 2 (interpolate position)
bad.gps.select <- which(dat2$transectID==40)[2]
last.good.gps <- which(dat2$transectID==40)[1]
first.good.gps <- which(dat2$transectID==40)[3]
dat2$GPS_lat[bad.gps.select] <- seq(dat2$GPS_lat[last.good.gps],dat2$GPS_lat[first.good.gps],length.out=length(bad.gps.select)+2)[-c(1,length(bad.gps.select)+2)]
dat2$GPS_lon[bad.gps.select] <- seq(dat2$GPS_lon[last.good.gps],dat2$GPS_lon[first.good.gps],length.out=length(bad.gps.select)+2)[-c(1,length(bad.gps.select)+2)]
## transect 49: 5 (interpolate position)
bad.gps.select <- which(dat2$transectID==49)[5]
last.good.gps <- which(dat2$transectID==49)[4]
first.good.gps <- which(dat2$transectID==49)[6]
dat2$GPS_lat[bad.gps.select] <- seq(dat2$GPS_lat[last.good.gps],dat2$GPS_lat[first.good.gps],length.out=length(bad.gps.select)+2)[-c(1,length(bad.gps.select)+2)]
dat2$GPS_lon[bad.gps.select] <- seq(dat2$GPS_lon[last.good.gps],dat2$GPS_lon[first.good.gps],length.out=length(bad.gps.select)+2)[-c(1,length(bad.gps.select)+2)]

dat2 <- dat2[-which(is.na(dat2$GPS_lat+dat2$GPS_lon)),]
dat2$transectID <- factor(dat2$transectID)

#### load ship track
ship.track <- read.table(paste0(gps.txt.path,"Totten_all.txt"),header=TRUE)
## transform am/pm to 24h format
ship.track.datetimeAMPM <- paste0(ship.track[,1]," ",ship.track[,2]," ",ship.track[,3])
ship.track.datetime24h <- strptime(ship.track.datetimeAMPM, "%e/%m/%Y %I:%M:%S %p", tz="UTC")
ship.track.date <- strptime(ship.track.datetimeAMPM, "%e/%m/%Y", tz="UTC")
transect.date <- strptime(dat2$Date, "%e-%m-%y", tz="UTC")
transect.time <- format(strptime(dat2$Time, "%H:%M", tz="UTC"), format="%H:%M:%S")
#transect.datetime24h <- strptime(test, "%e-%m-%y %H:%M", tz="UTC")
idx <- which(ship.track.date%in%transect.date)

spatial.dat <- data.frame(cbind(ship.track$longitude[idx], ship.track$latitude[idx]))
names(spatial.dat) <- c("Longitude","Latitude")
coordinates(spatial.dat) <- c("Longitude","Latitude")
proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
polar.dat.ship <- spTransform(spatial.dat, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

## overview
spatial.dat <- data.frame(cbind(dat2$GPS_lon, dat2$GPS_lat))
names(spatial.dat) <- c("Longitude","Latitude")
coordinates(spatial.dat) <- c("Longitude","Latitude")
proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
polar.dat <- spTransform(spatial.dat, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

# plot(r2, xlim=c(2170000,2320000), ylim=c(-1370000,-1220000))
# points(polar.dat)
# points(polar.dat.ship, col="red")
# text(polar.dat,labels=dat2$transectID, adj=c(2,2))

#par(mfrow=c(2,3),mar=c(5,4,4,1))
## each transect individually
for(i in 1:length(levels(dat2$transectID))){
#for(i in c(2,3,4,5,6,8)){
  loop.transect.ID <- levels(dat2$transectID)[i]
  loop.dat.select <- which(dat2$transectID==loop.transect.ID)
  loop.dat <- dat2[loop.dat.select,]
  message(i)
  print(paste0("transect ",loop.transect.ID))
  #print(loop.dat$GPS_lon)
  if(nrow(loop.dat)==0) next
  ## transform coordinates
  spatial.dat <- data.frame(cbind(loop.dat$GPS_lon, loop.dat$GPS_lat))
  names(spatial.dat) <- c("Longitude","Latitude")
  coordinates(spatial.dat) <- c("Longitude","Latitude")
  proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
  polar.dat <- spTransform(spatial.dat, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

  ## plot
  plot(polar.dat, main=paste0("transect ",loop.transect.ID))
  #points(polar.dat.ship, col="red", cex=0.5)
  text(polar.dat, labels=1:nrow(spatial.dat@coords), cex=0.8,col="blue",adj=-1)
  #text(polar.dat, labels=loop.dat$FileName, adj=-1,cex=0.5,col="blue")
  scalebar(100, type="bar",below="m")
}

dat3 <- dat2
## list of all filenames with gps data
#write.table(dat3$FileName, file="NBP1402_filenames_allfileswithgps.txt", eol=",", col.names=FALSE, row.names=FALSE)

#####################################################################
##### 2. SUBSET IMAGES FROM IMAGE LOCATIONS (& STORE FILENAMES) #####
#####################################################################

## dat3 already has all images without gps removed

## now remove images that are taken too close or too ar off the seafloor and copy them into a different folder to visually check
path.bad.images <- paste0(nbp.path, "bad_quality/")
#bad_images <- list.files(path.bad.images)
# save(bad_images, file="NBP1402_filenames_badimages.RData")
load("C:/Users/jjansen/Desktop/science/DP190101858_MappingAntarcticSeafloorBiodiversity/NBP1402_filenames_badimages.RData")
# file.copy(dat3$SourceFile[dat3$SUB1_Altitude>35|dat3$SUB1_Altitude<10], path.bad.images)
# file.copy(dat3$SourceFile[dat3$SUB1_Altitude<=35&dat3$SUB1_Altitude>=10], path.good.images)
# dat3$SUB1_Altitude[dat3$FileName=="tan1901_164_225.jpg"]
#file.copy(dat3$SourceFile[dat3$SUB1_Altitude>35|dat3$SUB1_Altitude<10], "D:/ARC_DP_data/Temp_Folder_imageswithoutdistance/")
## remove bad images first
usable.images.raw <- dat3[which(dat3$FileName%!in%bad_images),]
## and also those without gps data
usable.images.raw2 <- usable.images.raw[!is.na(usable.images.raw$GPS_lat+usable.images.raw$GPS_lon),]
## and those that seem to not extist in file
usable.images <- usable.images.raw2[-which(usable.images.raw2$FileName%!in%list.files(nbp.path)),]

## all transects are straight, but some gps variation: use distance to start point to calculate subset
dat <- data.frame(usable.images[,which(names(usable.images)%in%c("FileName","GPS_lon", "GPS_lat","Depth","time"))])[,c(4,3,1,5,2)]
dat$FileName_uniform <- substr(dat$FileName,1,6)
dat$FileName_uniform <- gsub("yoyo","",dat$FileName_uniform)
dat$FileName_uniform <- gsub("_DSC","",dat$FileName_uniform)
dat$FileName_uniform <- paste0("Transect",dat$FileName_uniform)
dat$transectID <- factor(dat$FileName_uniform)
## calculate transect length, define how many images to select and give images a random number
## SUBSET NOT REPRODUCIBLE BECAUSE set.seed() DIDN'T WORK
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
#save(dat,total.t.length.v, t.l, file="C:/Users/jjansen/Desktop/science/data_biological/NBP1402_dat.Rdata")
#load(file="C:/Users/jjansen/Desktop/science/data_biological/NBP1402_dat.Rdata")

t.images <- ceiling(total.t.length.v/100)
t.images.full <- ceiling(t.l/100)
## 500 images in total
#t.images <- ceiling((total.t.length.v/total.t.length)*500)
## 

spatial.dat <- data.frame(dat[,1:2])#[!is.na(rowSums(dat[,4:5])),]
coordinates(spatial.dat) <- c("GPS_lon","GPS_lat")
proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
polar.dat <- spTransform(spatial.dat, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

#par(mfrow=c(4,4), mar=c(5,4,2,1))
for(i in 1:length(levels(dat$transectID))){
  subset.v <- which(dat$transectID==levels(dat$transectID)[i])
  ## select the correct number of images for each transect
  # if(levels(dat$transectID)[i]%in%transects.dist.to.start){
  #   sel <- order(dat$image.select[subset.v])[1:t.images.full[i]]
  # }else if(levels(dat$transectID)[i]%in%transects.dist.to.end){
  #   sel <- order(dat$image.select[subset.v])[1:t.images.full[i]]
  # }else 
    sel <- order(dat$image.select[subset.v])[1:t.images[i]]
    if(any(is.na(sel))) sel <- sel[-which(is.na(sel))]
  ## create vector of selected images
  if(i==1){dat.sel <<-  subset.v[sel]
  }else dat.sel <-  c(dat.sel,subset.v[sel])
  png(filename=paste0("NBP1402_subset_transect_",unique(dat$transectID[subset.v]),".png"),width=1000, height=1000)
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
selected.filenames.renamed <- paste0("NBP1402_",dat$transectID[dat.sel],"_",sprintf("%04d",dat$image.select[dat.sel]),"__",selected.filenames)

## copy files into Annotation folder
img.path.origin <- paste0(img.path,selected.filenames)
img.path.destin <- "C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/Stills/Annotation_images/NBP1402/"
file.copy(img.path.origin,img.path.destin)

## write list of filenames into cropped Annotation folder for upload
img.path.destin_crop <- "C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/Stills/Annotation_images_cropped/NBP1402/"
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

## NBP1402

## @@@ !!! @@@ 
## crop images, do this in the terminal using the following code (after navigating into the "Stills" folder):
# cd C:
# cd C:\Users\jjansen\OneDrive - University of Tasmania\Desktop\science\data_biological\Stills\
# magick mogrify -path Annotation_images_cropped\NBP1402\ -gravity Center -crop 90% +repage Annotation_images\NBP1402\*jpg
## @@@ !!! @@@ 
















