library(magick)
library(geosphere)
library(ggplot2)
library(MBHdesign)
library(raster)
library(lubridate)
library(raadtools)
library(readxl)
'%!in%' <- function(x,y)!('%in%'(x,y))
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

env.dir <- "C:/Users/jjansen/Desktop/science/data_environmental/derived/"
image.dir <- "E:/ARC_DP_data/a_RawData_DirectFromContributors/"

## environmental data for plotting
# my_data_dir <- "C:/Users/jjansen/Desktop/science/data_environmental/accessed_through_R"
# set_data_roots(my_data_dir)
# r <- readtopo("ibcso")
# r2 <- r
# r2[r2>0] <- NA
# r2[r2<(-2500)] <- NA
r2 <- raster(paste0(env.dir,"Circumpolar_EnvData_500m_shelf_bathy_gebco_depth.grd"))
load(paste0(env.dir,"Circumpolar_Coastline.Rdata"))
stereo <- crs(coast.proj)#"+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

## 4 different surveys, so a number of metadata records and image folders need to be loaded:
bio.path <- "C:/Users/jjansen/Desktop/science/data_biological/"

image.dir_1 <- "E:/ARC_DP_data/a_RawData_DirectFromContributors/JR262/"
image.dir_2 <- "E:/ARC_DP_data/a_RawData_DirectFromContributors/JR15005/"
image.dir_3 <- "E:/ARC_DP_data/a_RawData_DirectFromContributors/JR17001/"
image.dir_4 <- "E:/ARC_DP_data/a_RawData_DirectFromContributors/JR17003/"

path.bad.images_1 <- "E:/ARC_DP_data/adjusted_SUCS/JR262/"
path.bad.images_2 <- "E:/ARC_DP_data/adjusted_SUCS/JR15005/"
path.bad.images_3 <- "E:/ARC_DP_data/adjusted_SUCS/JR17001/"
path.bad.images_4 <- "E:/ARC_DP_data/adjusted_SUCS/JR17003/"

dat.JR262.raw <- read_xlsx(paste0(image.dir_1,"JR262_SUCS_SouthOrkneys.xlsx"))
dat.JR262.raw <- dat.JR262.raw[-c(2,4,6),]
dat.JR262.raw[,c(6,5)] <- dat.JR262.raw[,c(6,5)]*-1
dat.JR15005.raw <- read.table(paste0(image.dir_2,"JR15005_SUCS.txt"),sep="\t", header=TRUE)
names(dat.JR15005.raw)[3:4] <- c("Latitude","Longitude")
dat.JR15005.raw$Time <- dmy_hm(dat.JR15005.raw$Time)
dat.JR17001.raw <- read_xlsx(paste0(image.dir_3,"JR17001_SUCS.xlsx"))
dat.JR17001.raw$Time <- ymd_hms(dat.JR17001.raw$Time)
dat.JR17003.raw <- read.csv(paste0(image.dir_4,"JR17003_SUCS.csv"))
dat.JR17003.raw$Event.No..Built.In...Integer.[which(dat.JR17003.raw$Event.No..Built.In...Integer.=="O11")] <- 11
names(dat.JR17003.raw)[9:10] <- c("Latitude","Longitude")
dat.JR17003.raw$Time <- dmy_hm(dat.JR17003.raw$Time)
lonlat <- rbind((dat.JR262.raw[,c(6,5)]),dat.JR15005.raw[,c(4,3)],dat.JR17001.raw[,c(6,5)],dat.JR17003.raw[,c(10,9)])

# spatial.dat <- data.frame(lonlat)
# spatial.dat <- spatial.dat[-which(is.na(rowSums(spatial.dat))),]
# coordinates(spatial.dat) <- c("Longitude","Latitude")
# proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
# polar.dat <- spTransform(spatial.dat, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# 
# plot(r2, xlim=c(-3000000,-1500000), ylim=c(700000,2500000), col=terrain.colors(99))
# plot(coast.proj, add=TRUE)
# points(polar.dat, pch=16, col="blue")
# 
# plot(r2, xlim=c(-2450000,-2050000), ylim=c(2100000,2450000), col=terrain.colors(99))
# plot(coast.proj, add=TRUE)
# points(polar.dat, pch=16, col="blue")
# 
# plot(r2, xlim=c(-2415000,-2370000), ylim=c(2245000,2270000), col=terrain.colors(99))
# points(polar.dat, pch=16, col="blue")

## JR267: lonlat per site, take 1 random image per site
## JR15005:
## JR17001: lonlat per image, but coded differently, 3 images per transect?
## JR17003: lonlat per image, but coded differently

# ## JR262
# spatial.dat <- data.frame(dat.JR262.raw[,c(6,5)])
# coordinates(spatial.dat) <- c("Longitude","Latitude")
# proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
# polar.dat.JR262.raw <- spTransform(spatial.dat, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# #save(polar.dat.JR262.raw, file="JR262_locations.Rdata")
# 
# plot(r2, xlim=c(-2500000,-2000000), ylim=c(2000000,2500000), col=terrain.colors(99))
# plot(coast.proj, add=TRUE)
# points(polar.dat.JR262.raw, pch=16, col="blue")
# 
# ## JR15005
# spatial.dat <- data.frame(dat.JR15005.raw[,c(4,3)])
# coordinates(spatial.dat) <- c("Longitude","Latitude")
# proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
# polar.dat.JR15005.raw <- spTransform(spatial.dat, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# #save(polar.dat.JR15005.raw, file="JR15005_locations.Rdata")
# 
# ## JR17001 & JR17003
# spatial.dat <- data.frame(dat.JR17001.raw[,c(6,5)])
# spatial.dat <- spatial.dat[-which(is.na(rowSums(spatial.dat))),]
# coordinates(spatial.dat) <- c("Longitude","Latitude")
# proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
# polar.dat.JR17001.raw <- spTransform(spatial.dat, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# 
# spatial.dat <- data.frame(dat.JR17003.raw[,c(10,9)])
# spatial.dat <- spatial.dat[-which(is.na(rowSums(spatial.dat))),]
# coordinates(spatial.dat) <- c("Longitude","Latitude")
# proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
# polar.dat.JR17003.raw <- spTransform(spatial.dat, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# 
# plot(r2, xlim=c(-2700000,-2300000), ylim=c(700000,1700000), col=terrain.colors(99))
# plot(coast.proj, add=TRUE)
# points(polar.dat.JR17001.raw, pch=16, col="blue")
# points(polar.dat.JR17003.raw, pch=16, col="cyan")
# 
# #save(polar.dat.JR17001.raw, file="JR17001_locations.Rdata")
# #save(polar.dat.JR17003.raw, file="JR17003_locations.Rdata")
# 
# plot(r2, xlim=c(-2635000,-2630000), ylim=c(1592000,1596000), col=terrain.colors(99))
# plot(coast.proj, add=TRUE)
# points(polar.dat.JR17001.raw, pch=16, col="blue")
# points(polar.dat.JR17003.raw, pch=16, col="cyan")

## image names
dir.files_1 <- list.files(image.dir_1,pattern=".", recursive=TRUE)
dir.files_1 <- dir.files_1[which(grepl(".png|.jpg", dir.files_1))]
dir.files_2 <- list.files(image.dir_2,pattern=".", recursive=TRUE)
dir.files_2 <- dir.files_2[which(grepl(".png|.jpg", dir.files_2))]
dir.files_3 <- list.files(image.dir_3,pattern=".", recursive=TRUE)
dir.files_3 <- dir.files_3[which(grepl(".png|.jpg", dir.files_3))]
dir.files_4 <- list.files(image.dir_4,pattern=".", recursive=TRUE)
dir.files_4 <- dir.files_4[which(grepl(".png|.jpg", dir.files_4))]

## bad images vs good images
bad.dir.files_1 <- list.files(path.bad.images_1,pattern=".", recursive=TRUE)
bad.dir.files_1 <- bad.dir.files_1[which(grepl(".png|.jpg", bad.dir.files_1))]
good.images_1 <- bad.dir.files_1[-grep("bad_quality",bad.dir.files_1)]
bad.dir.files_2 <- list.files(path.bad.images_2,pattern=".", recursive=TRUE)
bad.dir.files_2 <- bad.dir.files_2[which(grepl(".png|.jpg", bad.dir.files_2))]
good.images_2 <- bad.dir.files_2[-grep("bad_quality",bad.dir.files_2)]
bad.dir.files_3 <- list.files(path.bad.images_3,pattern=".", recursive=TRUE)
bad.dir.files_3 <- bad.dir.files_3[which(grepl(".png|.jpg", bad.dir.files_3))]
good.images_3 <- bad.dir.files_3[-grep("bad_quality",bad.dir.files_3)]
bad.dir.files_4 <- list.files(path.bad.images_4,pattern=".", recursive=TRUE)
bad.dir.files_4 <- bad.dir.files_4[which(grepl(".png|.jpg", bad.dir.files_4))]
good.images_4 <- bad.dir.files_4[-grep("bad_quality",bad.dir.files_4)]

dat.JR262 <- data.frame(cbind("JR262",dir.files_1,sub("/.*","",dir.files_1),sub(".*/","",dir.files_1)),NA,NA,NA,NA)
dat.JR15005 <- data.frame(cbind("JR15005",dir.files_2,sub("/.*","",dir.files_2),sub(".*/","",dir.files_2)),NA,NA,NA,NA)
dat.JR17001 <- data.frame(cbind("JR17001",dir.files_3,sub("/.*","",dir.files_3),sub(".*/","",dir.files_3)),NA,NA,NA,NA)
dat.JR17003 <- data.frame(cbind("JR17003",dir.files_4,sub("/.*","",dir.files_4),sub(".*/","",dir.files_4)),NA,NA,NA,NA)
names(dat.JR262) <- names(dat.JR15005) <- names(dat.JR17001) <- names(dat.JR17003) <- 
  c("surveyID", "filename_in_folder","transectID.raw", "filename","lon","lat","time","depth")
dat.JR15005$time <- as_datetime(dat.JR15005$time)
dat.JR17001$time <- as_datetime(dat.JR17001$time)
dat.JR17003$time <- as_datetime(dat.JR17003$time)

t <- sub("] stills","",dat.JR262$transectID.raw)
t <- sub(") stills","",t)
t <- sub("]","",t)
t <- substrRight(t,3)
t <- sub("t","",t)
t <- sub(" ","",t)
t <- as.factor(sprintf("%03d",as.numeric(t)))
dat.JR262$transectID <- t

t <- sub("Event_","",dat.JR15005$transectID.raw)
t <- as.factor(sprintf("%03d",as.numeric(t)))
dat.JR15005$transectID <- t

t <- sub("Anvers Island ","AI",dat.JR17001$transectID.raw)
t <- sub("King George Island ","KGI",t)
t <- gsub("\\([^\\)]+\\)","",t)
t <- as.factor(t)
dat.JR17001$transectID <- t

t <- sub("Event_","",dat.JR17003$transectID.raw)
t <- as.factor(t)
dat.JR17003$transectID <- t

##### assign lon-lat to each image
## JR262 (single value per station)
for(i in 1:4){
  sel <- which(grepl(dat.JR262.raw$Location[i], dat.JR262$transectID.raw))
  dat.JR262[sel,5:6] <- dat.JR262.raw[i,c(6,5)]
  dat.JR262$time[sel] <- dat.JR262.raw$Time[i]
  dat.JR262$depth[sel] <- dat.JR262.raw$Depth[i]
}
## JR15005 (interpolation)
dat.JR15005$image.number <- as.numeric(substr(sub(".*P","",dat.JR15005$filename),1,2))
dat.JR15005.raw$transectID <- as.factor(paste0("Event_",sprintf("%02d",dat.JR15005.raw$Event)))
dat.JR15005.raw$t.length <- NA
for(i in 1:nrow(dat.JR15005.raw)){
  dat.JR15005.raw$t.length[i] <- max(dat.JR15005$image.number[dat.JR15005$transectID.raw==dat.JR15005.raw$transectID[i]])
}
for(i in 1:length(levels(dat.JR15005.raw$transectID))){
  t.ID <- levels(dat.JR15005.raw$transectID)[i]
  message(t.ID)
  r.sel <- which(dat.JR15005.raw$transectID==t.ID)
  sel <- which(dat.JR15005$transectID.raw==t.ID)
  l <- dat.JR15005.raw$t.length[r.sel[1]]
  lats <- seq(dat.JR15005.raw$Latitude[r.sel[1]], dat.JR15005.raw$Latitude[r.sel[2]],length.out=l)
  lons <- seq(dat.JR15005.raw$Longitude[r.sel[1]], dat.JR15005.raw$Longitude[r.sel[2]],length.out=l)
  print(lats)
  time.seq <- seq(dat.JR15005.raw$Time[r.sel[1]], dat.JR15005.raw$Time[r.sel[2]],length.out=l)
  for(k in 1:l){
    n.sel <- which(!is.na(match(dat.JR15005$image.number[sel],k)))
    dat.JR15005$lat[sel][n.sel] <- lats[k]
    dat.JR15005$lon[sel][n.sel] <- lons[k]
    dat.JR15005$time[sel][n.sel] <- time.seq[k]
  }
}
## JR17001 (lat-lons for each image)
m.sel <- match(substr(dat.JR17001$filename,1,8),dat.JR17001.raw$`Image ID`)
dat.JR17001$lat <- dat.JR17001.raw$Latitude[m.sel]
dat.JR17001$lon <- dat.JR17001.raw$Longitude[m.sel]
dat.JR17001$time <- dat.JR17001.raw$Time[m.sel]
dat.JR17001$depth <- dat.JR17001.raw$Depth[m.sel]
dat.JR17001$temperature <- dat.JR17001.raw$temp[m.sel]
dat.JR17001$salinity <- dat.JR17001.raw$salinity[m.sel]
dat.JR17001$oxygen <- dat.JR17001.raw$oxygen[m.sel]
dat.JR17001$chla <- dat.JR17001.raw$chla[m.sel]

## JR17003 (lat-lons for each image)
dat.JR17003$image.code <- paste0(substr(dat.JR17003$filename,4,5),"_",substr(sub(".*_","",dat.JR17003$filename),1,4))
dat.JR17003.raw$image.code <- paste0(dat.JR17003.raw$Event.No..Built.In...Integer.,"_",dat.JR17003.raw$Photo.ID.Built.In...String.)
m.sel <- match(dat.JR17003$image.code,dat.JR17003.raw$image.code)
dat.JR17003$lat <- dat.JR17003.raw$Latitude.seatex.gga...seatex.gga.lat.[m.sel]
dat.JR17003$lon <- dat.JR17003.raw$Longitude.seatex.gga...seatex.gga.lon[m.sel]
dat.JR17003$depth <- dat.JR17003.raw$Water.depth.ea600...ea600.depth.[m.sel]
dat.JR17003$time <- dat.JR17003.raw$Time[m.sel]

##### bad images
good.dat.JR262 <- dat.JR262[which(dat.JR262$filename_in_folder%in%good.images_1),]
good.dat.JR15005 <- dat.JR15005[which(dat.JR15005$filename_in_folder%in%good.images_2),]
good.dat.JR17001 <- dat.JR17001[which(dat.JR17001$filename_in_folder%in%good.images_3),]
good.dat.JR17003 <- dat.JR17003[which(dat.JR17003$filename_in_folder%in%good.images_4),]

#####################################################################
##### 2. SUBSET IMAGES FROM IMAGE LOCATIONS (& STORE FILENAMES) #####
#####################################################################
## image subsetting and final dat file:
## JR262: 5 random images
## JR15002: spatially balanced random, using random subset of the duplicated image numbers
## JR17001: spatially balanced random
## JR17003: 

## metadata ???
#dat.sucs <- get(load(file=paste0(bio.path,"Stills/SUCS/SUCS_metadata.Rdata")))


####################################################################################################################
################################################ JR262 #############################################################
####################################################################################################################
###### RANDOM SUBSET
dat <- good.dat.JR262
dat$transectID <- factor(dat$transectID)
  
## calculate transect length, define how many images to select and give images a random number
dat$image.select <- NA
for(i in 1:length(levels(dat$transectID))){
  t.select <- levels(dat$transectID)[i]
  print(t.select)
  ##
  dat.subset.v <- which(dat$transectID==levels(dat$transectID)[i])
  t.counts <- length(dat.subset.v)
  #### subset images
  ## randomly select 5 images
    set.seed(43)
    dat$image.select[dat.subset.v][1:t.counts] <- sample(1:t.counts)
}

## add time from cruise report:
dat$time <- c(rep(ymd("2011-11-18"),8),rep(ymd("2011-11-19"),23))
## SAVE OUTPUT FOR FUTURE REFENCE (i.e. start here to add more images to the analysis)
# save(dat, file="C:/Users/jjansen/Desktop/science/data_biological/JR262_dat.Rdata")
#load(file="C:/Users/jjansen/Desktop/science/data_biological/JR262_dat.Rdata")

spatial.dat <- data.frame(dat[,5:6])
coordinates(spatial.dat) <- c("lon","lat")
proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
polar.dat <- spTransform(spatial.dat, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

col.breaks <- seq(-100,-1400,length.out=101)
cols <- terrain.colors(99)

par(mfrow=c(2,2), mar=c(5,4,2,1))
#par(mfrow=c(1,1), mar=c(5,4,2,1))
for(i in 1:length(levels(dat$transectID))){
  subset.v <- which(dat$transectID==levels(dat$transectID)[i])
  ## select the correct number of images for each transect
  t.counts <- length(subset.v)
  if(t.counts >=5){
    sel <- order(dat$image.select[subset.v])[1:5]
  }else sel <- order(dat$image.select[subset.v])[1:t.counts]
  print(sel)
  ## create vector of selected images
  if(i==1){dat.sel <<-  subset.v[sel]
  }else dat.sel <-  c(dat.sel,subset.v[sel])
  xlim <- NA
  ylim <- NA
  xlim[1] <- extent(polar.dat[subset.v])[1]-1000
  xlim[2] <- extent(polar.dat[subset.v])[2]+1000
  ylim[1] <- extent(polar.dat[subset.v])[3]-1000
  ylim[2] <- extent(polar.dat[subset.v])[4]+1000
  plot(r2, xlim=xlim,ylim=ylim, main=paste0("transect ",levels(dat$transectID)[i]),col=cols, breaks=col.breaks)
  plot(polar.dat[subset.v],add=TRUE)
  points(polar.dat[subset.v[sel]], col="red", pch=15)
  plot(coast.proj, add=TRUE)
  scalebar(100, type="bar")
}

dat$filename[dat.sel]
## filenames
selected.filenames <- dat$filename[dat.sel]
## create new filenames and identify folders where to find them (might be different for each survey)
selected.filenames.renamed <- paste0(dat$surveyID[dat.sel],"_",dat$transectID[dat.sel],"_",sprintf("%04d",dat$image.select[dat.sel]),"__",selected.filenames)
selected.filenames.folders <- paste0(gsub('/.*', '', dat$filename_in_folder)[dat.sel],"/")

## copy files into Annotation folder
img.path.origin <- paste0(image.dir,"JR262/",selected.filenames.folders,selected.filenames)
img.path.destin <- "C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/Stills/Annotation_images/JR262/"
file.copy(img.path.origin,img.path.destin)

## write list of filenames into cropped Annotation folder for upload
img.path.destin_crop <- "C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/Stills/Annotation_images_cropped/JR262/"
write.table(selected.filenames.renamed, paste0(img.path.destin_crop,"filenames.txt"), eol=",", col.names=FALSE, row.names=FALSE)

filenames.in.folder_original <- paste0(img.path.destin,selected.filenames)
filenames.in.folder_changed <- paste0(img.path.destin,selected.filenames.renamed)
file.rename(filenames.in.folder_original,filenames.in.folder_changed)


################################################################################################################
################################################ JR15005 #######################################################
################################################################################################################
dat <- good.dat.JR15005
dat$transectID <- factor(dat$transectID)

## calculate transect length, define how many images to select and give images a random number
dat$image.select <- NA
total.t.length.v <- NA
dat$dist.from.start <- NA
samp.v.list <- list()
center.idx.v <- NA
for(i in 1:length(levels(dat$transectID))){
  t.select <- levels(dat$transectID)[i]
  print(t.select)
  ##
  dat.subset.v <- which(dat$transectID==levels(dat$transectID)[i])
  t.counts <- length(dat.subset.v)
  ## select first and last coordinate to calculate transect distance
  v.first <- dat.subset.v[1]
  v.last <- tail(dat.subset.v, n=1)
  ## calculate distance of each image to the start of the transect along the transect line
  dat$dist.from.start[v.first] <- 0
  for(k in 2:t.counts){
    dist.to.previous.image <- distVincentyEllipsoid(dat[dat.subset.v[k-1],5:6], dat[dat.subset.v[k],5:6])
    dat$dist.from.start[dat.subset.v[k]] <- dat$dist.from.start[dat.subset.v[k-1]]+dist.to.previous.image
  }
  total.t.length.v[i] <- dat$dist.from.start[dat.subset.v[t.counts]]
  #### subset images
    ## choose a spatial random subset of images in each transect using quasiSamp on the distance from the start of the transect
    set.seed(2)
    samp <- quasiSamp(length(dat.subset.v)*5,dimension=1,potential.sites=dat$dist.from.start[dat.subset.v])
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
#dat$image.select[is.na(dat$image.select)] <- 9999

## SAVE OUTPUT FOR FUTURE REFENCE (i.e. start here to add more images to the analysis)
# save(dat,total.t.length.v, file="C:/Users/jjansen/Desktop/science/data_biological/JR15005_dat.Rdata")
#load(file="C:/Users/jjansen/Desktop/science/data_biological/JR15005_dat.Rdata")

barplot(round(total.t.length.v), names.arg=as.character(levels(dat$transectID)), las=2, main="JR15005", xlab="TransectID", ylab="length in m")

## total transect length across all the survey
#total.t.length <- sum(total.t.length.v)
## how many images from each transect, if we select on average 1 every 50m?
t.images <- ceiling(total.t.length.v/100)
t.images[t.images==0] <- 5
## 500 images in total
#t.images <- ceiling((total.t.length.v/total.t.length)*500)
## 

spatial.dat <- data.frame(dat[,5:6])
coordinates(spatial.dat) <- c("lon","lat")
proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
polar.dat <- spTransform(spatial.dat, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

col.breaks <- seq(-100,-1400,length.out=101)
cols <- terrain.colors(99)

par(mfrow=c(3,3), mar=c(4,3,2,1))
#par(mfrow=c(1,1), mar=c(5,4,2,1))
for(i in 1:length(levels(dat$transectID))){
  subset.v <- which(dat$transectID==levels(dat$transectID)[i])
  ## select the correct number of images for each transect
  sel <- order(dat$image.select[subset.v])[1:t.images[i]]
  print(sel)
  ## create vector of selected images
  if(i==1){dat.sel <<-  subset.v[sel]
  }else dat.sel <-  c(dat.sel,subset.v[sel])
  xlim <- NA
  ylim <- NA
  xlim[1] <- extent(polar.dat[subset.v])[1]-500
  xlim[2] <- extent(polar.dat[subset.v])[2]+500
  ylim[1] <- extent(polar.dat[subset.v])[3]-500
  ylim[2] <- extent(polar.dat[subset.v])[4]+500
  plot(r2, xlim=xlim,ylim=ylim, main=paste0("transect ",levels(dat$transectID)[i]),col=cols, breaks=col.breaks)
  plot(polar.dat[subset.v],add=TRUE)
  points(polar.dat[subset.v[sel]], col="red", pch=15)
  plot(coast.proj, add=TRUE)
  scalebar(100, type="bar")
}

dat$filename[dat.sel]
## filenames
selected.filenames <- dat$filename[dat.sel]
## create new filenames and identify folders where to find them (might be different for each survey)
selected.filenames.renamed <- paste0(dat$surveyID[dat.sel],"_",dat$transectID[dat.sel],"_",sprintf("%04d",dat$image.select[dat.sel]),"__",selected.filenames)
selected.filenames.folders <- paste0(gsub('/.*', '', dat$filename_in_folder)[dat.sel],"/")

## copy files into Annotation folder
img.path.origin <- paste0(image.dir,"JR15005/",selected.filenames.folders,selected.filenames)
img.path.destin <- "C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/Stills/Annotation_images/JR15005/"
file.copy(img.path.origin,img.path.destin)

## write list of filenames into cropped Annotation folder for upload
img.path.destin_crop <- "C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/Stills/Annotation_images_cropped/JR15005/"
write.table(selected.filenames.renamed, paste0(img.path.destin_crop,"filenames.txt"), eol=",", col.names=FALSE, row.names=FALSE)

filenames.in.folder_original <- paste0(img.path.destin,selected.filenames)
filenames.in.folder_changed <- paste0(img.path.destin,selected.filenames.renamed)
file.rename(filenames.in.folder_original,filenames.in.folder_changed)


################################################################################################################
################################################ JR17001 #######################################################
################################################################################################################
dat <- good.dat.JR17001[-which(is.na(good.dat.JR17001$lon)),]
dat$transectID <- factor(dat$transectID)

## calculate transect length, define how many images to select and give images a random number
## IMAGE SELECT IS NOT REPRODUCIBLE BECAUSE THE SEED WAS NOT SET WHEN I RAN IT!!!
dat$image.select <- NA
total.t.length.v <- NA
dat$dist.from.start <- NA
samp.v.list <- list()
center.idx.v <- NA
for(i in 1:length(levels(dat$transectID))){
  t.select <- levels(dat$transectID)[i]
  print(t.select)
  ##
  dat.subset.v <- which(dat$transectID==levels(dat$transectID)[i])
  t.counts <- length(dat.subset.v)
  ## select first and last coordinate to calculate transect distance
  v.first <- dat.subset.v[1]
  v.last <- tail(dat.subset.v, n=1)
  ## calculate distance of each image to the start of the transect along the transect line
  dat$dist.from.start[v.first] <- 0
  for(k in 2:t.counts){
    dist.to.previous.image <- distVincentyEllipsoid(dat[dat.subset.v[k-1],5:6], dat[dat.subset.v[k],5:6])
    dat$dist.from.start[dat.subset.v[k]] <- dat$dist.from.start[dat.subset.v[k-1]]+dist.to.previous.image
  }
  total.t.length.v[i] <- dat$dist.from.start[dat.subset.v[t.counts]]
  #### subset images
  ## choose a spatial random subset of images in each transect using quasiSamp on the distance from the start of the transect
  #set.seed(4)
  samp <- quasiSamp(length(dat.subset.v)*5,dimension=1,potential.sites=dat$dist.from.start[dat.subset.v])
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
#dat$image.select[is.na(dat$image.select)] <- 9999

## SAVE OUTPUT FOR FUTURE REFENCE (i.e. start here to add more images to the analysis)
# save(dat,total.t.length.v, file="C:/Users/jjansen/Desktop/science/data_biological/JR17001_dat.Rdata")
#load(file="C:/Users/jjansen/Desktop/science/data_biological/JR17001_dat.Rdata")

barplot(round(total.t.length.v), names.arg=as.character(levels(dat$transectID)), las=2, main="JR17001", xlab="TransectID", ylab="length in m")

## total transect length across all the survey
#total.t.length <- sum(total.t.length.v)
## how many images from each transect, if we select on average 1 every 50m?
t.images <- ceiling(total.t.length.v/100)
t.images[t.images==0] <- 5
## 500 images in total
#t.images <- ceiling((total.t.length.v/total.t.length)*500)
## 

spatial.dat <- data.frame(dat[,5:6])
coordinates(spatial.dat) <- c("lon","lat")
proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
polar.dat <- spTransform(spatial.dat, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

col.breaks <- seq(-100,-1400,length.out=101)
cols <- terrain.colors(99)

par(mfrow=c(2,3), mar=c(5,4,2,1))
#par(mfrow=c(1,1), mar=c(5,4,2,1))
for(i in 1:length(levels(dat$transectID))){
  subset.v <- which(dat$transectID==levels(dat$transectID)[i])
  ## select the correct number of images for each transect
  sel <- order(dat$image.select[subset.v])[1:t.images[i]]
  print(sel)
  ## create vector of selected images
  if(i==1){dat.sel <<-  subset.v[sel]
  }else dat.sel <-  c(dat.sel,subset.v[sel])
  xlim <- NA
  ylim <- NA
  xlim[1] <- extent(polar.dat[subset.v])[1]-500
  xlim[2] <- extent(polar.dat[subset.v])[2]+500
  ylim[1] <- extent(polar.dat[subset.v])[3]-500
  ylim[2] <- extent(polar.dat[subset.v])[4]+500
  plot(r2, xlim=xlim,ylim=ylim, main=paste0("transect ",levels(dat$transectID)[i]),col=cols, breaks=col.breaks)
  plot(polar.dat[subset.v],add=TRUE)
  points(polar.dat[subset.v[sel]], col="red", pch=15)
  plot(coast.proj, add=TRUE)
  scalebar(100, type="bar")
}

dat$filename[dat.sel]
## filenames
selected.filenames <- dat$filename[dat.sel]
## create new filenames and identify folders where to find them (might be different for each survey)
selected.filenames.renamed <- paste0(dat$surveyID[dat.sel],"_",dat$transectID[dat.sel],"_",sprintf("%04d",dat$image.select[dat.sel]),"__",selected.filenames)
selected.filenames.folders <- paste0(gsub('/.*', '', dat$filename_in_folder)[dat.sel],"/")

## copy files into Annotation folder
img.path.origin <- paste0(image.dir,"JR17001/",selected.filenames.folders,selected.filenames)
img.path.destin <- "C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/Stills/Annotation_images/JR17001/"
file.copy(img.path.origin,img.path.destin)

## write list of filenames into cropped Annotation folder for upload
img.path.destin_crop <- "C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/Stills/Annotation_images_cropped/JR17001/"
write.table(selected.filenames.renamed, paste0(img.path.destin_crop,"filenames.txt"), eol=",", col.names=FALSE, row.names=FALSE)

filenames.in.folder_original <- paste0(img.path.destin,selected.filenames)
filenames.in.folder_changed <- paste0(img.path.destin,selected.filenames.renamed)
file.rename(filenames.in.folder_original,filenames.in.folder_changed)

################################################################################################################
################################################ JR17003 #######################################################
################################################################################################################
dat <- good.dat.JR17003[-which(is.na(good.dat.JR17003$lon)),]
dat$transectID <- factor(dat$transectID)

## calculate transect length, define how many images to select and give images a random number
## IMAGE SELECT IS NOT REPRODUCIBLE BECAUSE THE SEED WAS NOT SET WHEN I RAN IT!!!
dat$image.select <- NA
total.t.length.v <- NA
dat$dist.from.start <- NA
samp.v.list <- list()
center.idx.v <- NA
for(i in 1:length(levels(dat$transectID))){
  t.select <- levels(dat$transectID)[i]
  print(t.select)
  ##
  dat.subset.v <- which(dat$transectID==levels(dat$transectID)[i])
  t.counts <- length(dat.subset.v)
  ## select first and last coordinate to calculate transect distance
  v.first <- dat.subset.v[1]
  v.last <- tail(dat.subset.v, n=1)
  ## calculate distance of each image to the start of the transect along the transect line
  dat$dist.from.start[v.first] <- 0
  for(k in 2:t.counts){
    dist.to.previous.image <- distVincentyEllipsoid(dat[dat.subset.v[k-1],5:6], dat[dat.subset.v[k],5:6])
    dat$dist.from.start[dat.subset.v[k]] <- dat$dist.from.start[dat.subset.v[k-1]]+dist.to.previous.image
  }
  total.t.length.v[i] <- dat$dist.from.start[dat.subset.v[t.counts]]
  #### subset images
  ## choose a spatial random subset of images in each transect using quasiSamp on the distance from the start of the transect
  #set.seed(2)
  samp <- quasiSamp(length(dat.subset.v)*5,dimension=1,potential.sites=dat$dist.from.start[dat.subset.v])
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
#dat$image.select[is.na(dat$image.select)] <- 9999

## SAVE OUTPUT FOR FUTURE REFENCE (i.e. start here to add more images to the analysis)
# save(dat,total.t.length.v, file="C:/Users/jjansen/Desktop/science/data_biological/JR17003_dat.Rdata")
#load(file="C:/Users/jjansen/Desktop/science/data_biological/JR17003_dat.Rdata")

barplot(round(total.t.length.v), names.arg=as.character(levels(dat$transectID)), las=2, main="JR17003", xlab="TransectID", ylab="length in m")

## total transect length across all the survey
#total.t.length <- sum(total.t.length.v)
## how many images from each transect, if we select on average 1 every 50m?
t.images <- ceiling(total.t.length.v/100)
t.images[t.images==0] <- 5
## 500 images in total
#t.images <- ceiling((total.t.length.v/total.t.length)*500)
## 

## transect 011 should have less images
## transect length without the 9th image:
total.t.length.v2 <- distVincentyEllipsoid(dat[1,5:6], dat[8,5:6])
t.images2 <- ceiling(total.t.length.v2/100)

spatial.dat <- data.frame(dat[,5:6])
coordinates(spatial.dat) <- c("lon","lat")
proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
polar.dat <- spTransform(spatial.dat, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

col.breaks <- seq(-100,-1400,length.out=101)
cols <- terrain.colors(99)

par(mfrow=c(3,4), mar=c(5,4,2,1))
#par(mfrow=c(1,1), mar=c(5,4,2,1))
for(i in 1:length(levels(dat$transectID))){
  subset.v <- which(dat$transectID==levels(dat$transectID)[i])
  ## select the correct number of images for each transect
  if(i==1){
    sel <- order(dat$image.select[subset.v])[2]
  }else sel <- order(dat$image.select[subset.v])[1:t.images[i]]
  print(sel)
  ## create vector of selected images
  if(i==1){dat.sel <<-  c(9,subset.v[sel])
  }else dat.sel <-  c(dat.sel,subset.v[sel])
  ## create vector of selected images
  xlim <- NA
  ylim <- NA
  xlim[1] <- extent(polar.dat[subset.v])[1]-500
  xlim[2] <- extent(polar.dat[subset.v])[2]+500
  ylim[1] <- extent(polar.dat[subset.v])[3]-500
  ylim[2] <- extent(polar.dat[subset.v])[4]+500
  plot(r2, xlim=xlim,ylim=ylim, main=paste0("transect ",levels(dat$transectID)[i]),col=cols, breaks=col.breaks)
  plot(polar.dat[subset.v],add=TRUE)
  if(i==1){points(polar.dat[c(9,subset.v[sel])], col="red", pch=15)
  }else points(polar.dat[subset.v[sel]], col="red", pch=15)
   plot(coast.proj, add=TRUE)
  scalebar(100, type="bar")
}

dat$filename[dat.sel]
## filenames
selected.filenames <- dat$filename[dat.sel]
## create new filenames and identify folders where to find them (might be different for each survey)
selected.filenames.renamed <- paste0(dat$surveyID[dat.sel],"_",dat$transectID[dat.sel],"_",sprintf("%04d",dat$image.select[dat.sel]),"__",selected.filenames)
selected.filenames.folders <- paste0(gsub('/.*', '', dat$filename_in_folder)[dat.sel],"/")

## copy files into Annotation folder
img.path.origin <- paste0(image.dir,"JR17003/",selected.filenames.folders,selected.filenames)
img.path.destin <- "C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/Stills/Annotation_images/JR17003/"
file.copy(img.path.origin,img.path.destin)

## write list of filenames into cropped Annotation folder for upload
img.path.destin_crop <- "C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/Stills/Annotation_images_cropped/JR17003/"
write.table(selected.filenames.renamed, paste0(img.path.destin_crop,"filenames.txt"), eol=",", col.names=FALSE, row.names=FALSE)

filenames.in.folder_original <- paste0(img.path.destin,selected.filenames)
filenames.in.folder_changed <- paste0(img.path.destin,selected.filenames.renamed)
file.rename(filenames.in.folder_original,filenames.in.folder_changed)






########################################################
##### 3. CROP SELECTED IMAGES AND STORE SEPARATELY #####
########################################################

## NO CROPPING

































##################################
######### PS18
##################################
## sort by descending transectID
dat.PS18 <- dat.PS18[order(substring(dat.PS18[,1],1,3)),]

dat.PS18$transect_duration <- as.numeric(ymd_hms(dat.PS18$`Date/Time_End`)-ymd_hms(dat.PS18$`Date/Time_Start`))
dat.PS18$transect_distance <- distHaversine(p1=dat.PS18[,c(11,8)], p2=dat.PS18[,c(10,7)])
dat.PS18$transect_speed <- 0.54*((dat.PS18$transect_distance/1000)/(dat.PS18$transect_duration/60))
plot(dat.PS18$transect_speed, ylab="speed in knots", xlab="transect-index")
text(dat.PS18$transect_speed,adj=c(0.5,-0.5))
abline(1,0, col="red",lty=2)

t.toofast <- which(dat.PS18$transect_speed>1)
## TRANSECTS 12 and 13 ARE VERY SUSPICIOUSLY FAST

#### load image names
files.dat <- data.frame(matrix(ncol=3,nrow=length(dir.files)))
files.dat[,1] <- dir.files
files.dat[,2] <- substring(dir.files,6,10)
files.dat[,2] <- as.factor(files.dat[,2])
files.dat[,3] <- substring(dir.files,6,20)
files.dat[,3] <- sub("^[^_]*", "",files.dat[,3])
files.dat[,3] <- substring(files.dat[,3],2,20)
files.dat[,3] <- sub(".jpg", "", files.dat[,3])
files.dat[,4] <- as.numeric(substring(files.dat[,3],1,4))
files.dat[,5] <- NA
transect.length <- NA
ntransects <- length(levels(files.dat[,2]))
for(i in 1:ntransects){
  t.sel <- which(files.dat[,2]==levels(files.dat[,2])[i])
  a <- files.dat[t.sel,4]
  b <- a-head(a,1)+1
  files.dat[t.sel,5] <- b
  transect.length[i] <- tail(b,1)
  message(i)
  message(levels(files.dat[,2])[i])
  print(a)
  print(b)
}
names(files.dat) <- c("Filename","transectID","imageID","imageNumber","imageSequence")

## generate lat-lon positions for each potential image using start/end points and the total number of images taken in that transect
new.pos.lat <- list()
new.pos.lon <- list()
for(i in c(1:23,25)){
  new.pos.lat[[i]] <- seq(dat.PS18$Latitude_Start[i], dat.PS18$Latitude_End[i],length.out=transect.length[i])
  new.pos.lon[[i]] <- seq(dat.PS18$Longitude_Start[i], dat.PS18$Longitude_End[i],length.out=transect.length[i])
}
## assing lon/lat for the actual images taken
## but replace lon/lat positions for all images on transects that were unrealistically fast with starting lon-lat only
files.dat[,6] <- NA
files.dat[,7] <- NA
names(files.dat)[6:7] <- c("lon","lat")
for(i in c(1:23,25)){
  t.sel <- which(files.dat[,2]==levels(files.dat[,2])[i])
  img.sel <- files.dat[t.sel,5]
  if(i%in%t.toofast){
    files.dat[t.sel,6] <- dat.PS18$Longitude_Start[i]
    files.dat[t.sel,7] <- dat.PS18$Latitude_Start[i]
  }else{
    files.dat[t.sel,6] <- new.pos.lon[[i]][img.sel]
    files.dat[t.sel,7] <- new.pos.lat[[i]][img.sel]
  }
}
PS18_image_metadata <- files.dat
## transect 24 has no start/end, so replace with median points
PS18_image_metadata[1558:1624,6] <- dat.PS18$Median_Longitude[24]
PS18_image_metadata[1558:1624,7] <- dat.PS18$Median_Latitude[24]
#save(PS18_image_metadata, file=paste0(bio.path,"Stills/PS18/PS18_image_metadata.Rdata"))




