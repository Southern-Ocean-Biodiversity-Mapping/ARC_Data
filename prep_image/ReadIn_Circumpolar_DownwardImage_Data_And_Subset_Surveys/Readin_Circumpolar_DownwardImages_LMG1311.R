library(magick)
library(geosphere)
library(ggplot2)
library(MBHdesign)
library(raster)
library(lubridate)
library(raadtools)
library(readxl)
library(exifr)
'%!in%' <- function(x,y)!('%in%'(x,y))
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

env.dir <- "C:/Users/jjansen/Desktop/science/data_environmental/derived/"
image.dir <- "D:/ARC_DP_data/a_RawData_DirectFromContributors/"

## environmental data for plotting
# my_data_dir <- "C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_environmental/accessed_through_R"
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

image.dir <- "D:/ARC_DP_data/a_RawData_DirectFromContributors/LMG1311/"
path.bad.images <- "D:/ARC_DP_data/adjusted_LMG1311/"

dat.LMG1311.raw <- read.csv(paste0(image.dir,"LMG13-11 Station Log -DOM-IMA-8XNTXT2.csv"))[-c(3:4),-11]
dat.LMG1311.raw$surveyID <- "LMG1311"
dat.LMG1311.raw$lon_start <- as.numeric(measurements::conv_unit(gsub('°',' ',dat.LMG1311.raw$Longitude.W.start), from='deg_dec_min',to='dec_deg'))*-1
dat.LMG1311.raw$lat_start <- as.numeric(measurements::conv_unit(gsub('°',' ',dat.LMG1311.raw$Latitude.S.start), from='deg_dec_min',to='dec_deg'))*-1
dat.LMG1311.raw$lon_end <- as.numeric(measurements::conv_unit(gsub('°',' ',dat.LMG1311.raw$Longitude.W.end), from='deg_dec_min',to='dec_deg'))*-1
dat.LMG1311.raw$lat_end <- as.numeric(measurements::conv_unit(gsub('°',' ',dat.LMG1311.raw$Latitude.S.end), from='deg_dec_min',to='dec_deg'))*-1
lons <- as.numeric(c(dat.LMG1311.raw$lon_start, dat.LMG1311.raw$lon_end))
lats <- as.numeric(c(dat.LMG1311.raw$lat_start, dat.LMG1311.raw$lat_end))
dat.LMG1311.raw$depth.start <- c(355,446)
dat.LMG1311.raw$depth.end <- c(403,675.8)

spatial.dat <- data.frame(cbind(lons,lats))
coordinates(spatial.dat) <- c("lons","lats")
proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
polar.dat <- spTransform(spatial.dat, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
plot(r2,xlim=c(-2510000,-2460000),ylim=c(1250000,1300000))
points(polar.dat)

## image names
dir.files <- list.files(image.dir,pattern=".JPG", recursive=TRUE)

## bad images vs good images
bad.dir.files <- list.files(path.bad.images,pattern=".JPG", recursive=TRUE)
good.images <- sub(".*/","",bad.dir.files[-grep("bad_quality",bad.dir.files)])

dat.LMG1311 <- data.frame(cbind("LMG1311",NA,dir.files,sub(".*/","",dir.files)),NA,NA,NA,NA,NA)
names(dat.LMG1311) <- c("surveyID", "transectID","filename_in_folder","filename","lon","lat","DateTime", "depth.start", "depth.end")
dat.LMG1311$DateTime <- read_exif(paste0(image.dir,dir.files))$CreateDate
dat.LMG1311$transectID[1:36] <- "2"
dat.LMG1311$transectID[37:263] <- "3"
dat.LMG1311$depth.start[1:36] <- dat.LMG1311.raw$depth.start[1]
dat.LMG1311$depth.start[37:263] <- dat.LMG1311.raw$depth.start[2]
dat.LMG1311$depth.end[1:36] <- dat.LMG1311.raw$depth.end[1]
dat.LMG1311$depth.end[37:263] <- dat.LMG1311.raw$depth.end[2]

transect2_start_image <- 5
transect2_end_image <- 36
transect3_start_image <- 52
transect3_end_image <- 261
dat.LMG1311$lon[transect2_start_image] <- dat.LMG1311.raw$lon_start[1]
dat.LMG1311$lon[transect3_start_image] <- dat.LMG1311.raw$lon_start[2]
dat.LMG1311$lon[transect2_end_image] <- dat.LMG1311.raw$lon_end[1]
dat.LMG1311$lon[transect3_end_image] <- dat.LMG1311.raw$lon_end[2]
dat.LMG1311$lat[transect2_start_image] <- dat.LMG1311.raw$lat_start[1]
dat.LMG1311$lat[transect3_start_image] <- dat.LMG1311.raw$lat_start[2]
dat.LMG1311$lat[transect2_end_image] <- dat.LMG1311.raw$lat_end[1]
dat.LMG1311$lat[transect3_end_image] <- dat.LMG1311.raw$lat_end[2]

## time intervall between images
dat.LMG1311$Timer <- NA
times <- ymd_hms(dat.LMG1311$DateTime)
s <- seconds(times)
dat.LMG1311$Timer[5:36] <- s[5:36]-s[5]
dat.LMG1311$Timer[52:261] <- s[52:261]-s[52]

## interpolate coordinates according to time-intervals
lon_seq1 <- seq(dat.LMG1311.raw$lon_start[1],dat.LMG1311.raw$lon_end[1], length.out=dat.LMG1311$Timer[36])
lat_seq1 <- seq(dat.LMG1311.raw$lat_start[1],dat.LMG1311.raw$lat_end[1], length.out=dat.LMG1311$Timer[36])
lon_seq2 <- seq(dat.LMG1311.raw$lon_start[2],dat.LMG1311.raw$lon_end[2], length.out=dat.LMG1311$Timer[261])
lat_seq2 <- seq(dat.LMG1311.raw$lat_start[2],dat.LMG1311.raw$lat_end[2], length.out=dat.LMG1311$Timer[261])
dat.LMG1311$lon[5:36] <- lon_seq1[c(1,dat.LMG1311$Timer[6:36])]
dat.LMG1311$lat[5:36] <- lat_seq1[c(1,dat.LMG1311$Timer[6:36])]
dat.LMG1311$lon[52:261] <- lon_seq2[c(1,dat.LMG1311$Timer[53:261])]
dat.LMG1311$lat[52:261] <- lat_seq2[c(1,dat.LMG1311$Timer[53:261])]

spatial.dat <- data.frame(dat.LMG1311[which(!is.na(dat.LMG1311$lon)),5:6])
coordinates(spatial.dat) <- c("lon","lat")
proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
polar.dat <- spTransform(spatial.dat, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
plot(r2,xlim=c(-2500000,-2480000),ylim=c(1265000,1280000))
points(polar.dat)

##### bad images
good.dat.LMG1311 <- dat.LMG1311[which(dat.LMG1311$filename%in%good.images),]

#####################################################################
##### 2. SUBSET IMAGES FROM IMAGE LOCATIONS (& STORE FILENAMES) #####
#####################################################################
dat <- good.dat.LMG1311
dat$transectID <- factor(dat$transectID)
## calculate transect length, define how many images to select and give images a random number

## transect 2 has 3 images separated from the rest
## select these images and then base the subset on the remaining transect 
## THIS IS NOT REPRODUCIBLE, SOMETHING WENT WRONG WHEN SETTING THE SEED
dat$image.select <- NA
#dat$image.select[1:3] <- 1:3
total.t.length.v <- NA
dat$dist.from.start <- NA
samp.v.list <- list()
center.idx.v <- NA
for(i in 1:length(levels(dat$transectID))){
  t.select <- levels(dat$transectID)[i]
  print(t.select)
  ##
  # if(i==1){
  #   dat.subset.v <- which(dat$transectID==levels(dat$transectID)[i])[-c(1:3)]
  # }else
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
  #set.seed(42)
  samp <- quasiSamp(length(dat.subset.v)*5,dimension=1,potential.sites=dat$dist.from.start[dat.subset.v])
  samp.v <- samp$ID[which(duplicated(samp)==FALSE)]
  print(paste0("first double selection for image # ", which(duplicated(samp)==TRUE)[1]))
  print(paste0("# of images in first batch: ", ceiling(total.t.length.v[i]/50)))
  #samp.v[(length(samp.v)+1):length(subset.v)] <- 9998
  # if(i==1){
  #   dat$image.select[dat.subset.v[samp.v]] <-  4:(length(samp.v)+3)
  # } else
  dat$image.select[dat.subset.v[samp.v]] <-  1:length(samp.v)
  dat$image.select[dat.subset.v[-samp.v]] <-  999
  # ## choose a random subset of images in each transect by giving them random numbers starting from 1
  #dat$image.select[dat.subset.v] <- sample(1:length(dat.subset.v),length(dat.subset.v))
  samp.v.list[[i]] <- samp.v
}
## 
#dat$image.select[is.na(dat$image.select)] <- 9999

## SAVE OUTPUT FOR FUTURE REFENCE (i.e. start here to add more images to the analysis)
# save(dat,total.t.length.v, file="C:/Users/jjansen/Desktop/science/data_biological/LMG1311_dat.Rdata")
#load(file="C:/Users/jjansen/Desktop/science/data_biological/LMG1311_dat.Rdata")

barplot(round(total.t.length.v), names.arg=as.character(levels(dat$transectID)), las=2, main="LMG1311", xlab="TransectID", ylab="length in m")

## total transect length across all the survey
#total.t.length <- sum(total.t.length.v)
## how many images from each transect, if we select on average 1 every 50m?
t.images <- ceiling(total.t.length.v/100)

## transect 2 should have less images
## transect length without the first three images:
total.t.length.v2 <- distVincentyEllipsoid(dat[4,5:6], dat[31,5:6])
t.images2 <- ceiling(total.t.length.v2/100)

spatial.dat <- data.frame(dat[,5:6])
coordinates(spatial.dat) <- c("lon","lat")
proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
polar.dat <- spTransform(spatial.dat, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

col.breaks <- seq(-100,-1400,length.out=101)
cols <- terrain.colors(99)

par(mfrow=c(1,1), mar=c(5,4,2,1))
#par(mfrow=c(1,1), mar=c(5,4,2,1))
for(i in 1:length(levels(dat$transectID))){
  if(i==1){
    subset.v <- which(dat$transectID==levels(dat$transectID)[i])[-c(1:3)]
  }else subset.v <- which(dat$transectID==levels(dat$transectID)[i])
  
  ## select the correct number of images for each transect
  if(i==1){
    sel <- order(dat$image.select[subset.v])[1:(t.images2)]
  }else sel <- order(dat$image.select[subset.v])[1:t.images[i]]
  print(sel)
  ## create vector of selected images
  if(i==1){dat.sel <<-  c(1:3,subset.v[sel])
  }else dat.sel <-  c(dat.sel,subset.v[sel])
  xlim <- NA
  ylim <- NA
  xlim[1] <- extent(polar.dat[subset.v])[1]-1000
  xlim[2] <- extent(polar.dat[subset.v])[2]+1000
  ylim[1] <- extent(polar.dat[subset.v])[3]-1000
  ylim[2] <- extent(polar.dat[subset.v])[4]+1000
  plot(r2, xlim=xlim,ylim=ylim, main=paste0("transect ",levels(dat$transectID)[i]),col=cols, breaks=col.breaks)
  plot(polar.dat[subset.v],add=TRUE)
  if(i==1){points(polar.dat[c(1:3,subset.v[sel])], col="red", pch=15)
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
img.path.origin <- paste0(image.dir,selected.filenames.folders,selected.filenames)
img.path.destin <- "C:/Users/jjansen/Desktop/science/data_biological/Stills/Annotation_images/LMG1311/"
file.copy(img.path.origin,img.path.destin)

## write list of filenames into cropped Annotation folder for upload
img.path.destin_crop <- "C:/Users/jjansen/Desktop/science/data_biological/Stills/Annotation_images_cropped/LMG1311/"
write.table(selected.filenames.renamed, paste0(img.path.destin_crop,"LMG1311_filenames.txt"), eol=",", col.names=FALSE, row.names=FALSE)

filenames.in.folder_original <- paste0(img.path.destin,selected.filenames)
filenames.in.folder_changed <- paste0(img.path.destin,selected.filenames.renamed)
file.rename(filenames.in.folder_original,filenames.in.folder_changed)

########################################################
##### 3. CROP SELECTED IMAGES AND STORE SEPARATELY #####
########################################################

## NO CROPPING






























