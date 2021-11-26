#########################################################################################################################
## Ross Sea images: Assigning lat/lon to each image, removing bad images and then subsetting the full dataset
#########################################################################################################################

## NBP1502 ##

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
library(measurements)
'%!in%' <- function(x,y)!('%in%'(x,y))
env.dir <- "C:/Users/jjansen/Desktop/science/data_environmental/derived/"

#### specify paths
img.path2 <- "D:/ARC_DP_data/a_RawData_DirectFromContributors/NBP1502/"
img.path <- "D:/ARC_DP_data/adjusted_NBP1502/"
txt.path <- "D:/ARC_DP_data/a_RawData_DirectFromContributors/NBP1502/"
gps.txt.path <- "D:/ARC_DP_data/a_RawData_DirectFromContributors/NBP1502/"

my_data_dir <- "C:/Users/jjansen/Desktop/science/data_environmental/accessed_through_R"

r2 <- raster(paste0(env.dir,"Circumpolar_EnvData_500m_shelf_bathy_gebco_depth.grd"))
load(paste0(env.dir,"Circumpolar_Coastline.Rdata"))
stereo <- crs(coast.proj)#"+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

###############################################################################
##
dat.raw <- read.table(paste0(txt.path,"NBP1502_locations.txt"), sep="\t", header=TRUE)
names(dat.raw)[3:4] <- c("Lat_deg", "Lon_deg")
dat.raw$Lat <- measurements::conv_unit(dat.raw$Lat_deg, from = 'deg_dec_min', to = 'dec_deg')
dat.raw$Lon <- measurements::conv_unit(dat.raw$Lon_deg, from = 'deg_dec_min', to = 'dec_deg')

NBP1502_folders <- c("NBP1502_YoYo_photos_2015-01-27_BC1/",
                     "NBP1502_YoYo_photos_2015-02-03A_GZW4_KC3/",
                     "NBP1502_YoYo_photos_2015-02-03B_KC5",
                     "NBP1502_YoYo_photos_2015-02-10_KC6",
                     "NBP1502_YoYo_photos_2015-02-17_KC12")

transect_images <- list.files(paste0(img.path,NBP1502_folders[1]),pattern=".jpg")
image_numbers <- as.numeric(gsub(".jpg"," ",gsub("DSC","",transect_images)))
image_positions <- image_numbers-image_numbers[1]+1
transect_length <- tail(image_numbers,1)-image_numbers[1]+1
transect_lats <- seq(dat.raw$Lat[1], dat.raw$Lat[2], length.out=transect_length)
transect_lons <- seq(dat.raw$Lon[1], dat.raw$Lon[2], length.out=transect_length)
dat.1 <- data.frame("FileName"=transect_images,"transectID"="BC1","GPS_lat"=NA,"GPS_lon"=NA, foldername=NBP1502_folders[1])
dat.1$GPS_lat <- transect_lats[image_positions]
dat.1$GPS_lon <- transect_lons[image_positions]

transect_images <- list.files(paste0(img.path,NBP1502_folders[2]),pattern=".jpg")
image_numbers <- as.numeric(gsub(".jpg"," ",gsub("DSC","",transect_images)))
image_positions <- image_numbers-image_numbers[1]+1
transect_length <- tail(image_numbers,1)-image_numbers[1]+1
transect_lats <- seq(dat.raw$Lat[3], dat.raw$Lat[4], length.out=transect_length)
transect_lons <- seq(dat.raw$Lon[3], dat.raw$Lon[4], length.out=transect_length)
dat.2 <- data.frame("FileName"=transect_images,"transectID"="KC3","GPS_lat"=NA,"GPS_lon"=NA, foldername=NBP1502_folders[2])
dat.2$GPS_lat <- transect_lats[image_positions]
dat.2$GPS_lon <- transect_lons[image_positions]

transect_images <- list.files(paste0(img.path,NBP1502_folders[3]),pattern=".jpg")
image_numbers <- as.numeric(gsub(".jpg"," ",gsub("DSC","",transect_images)))
image_positions <- image_numbers-image_numbers[1]+1
transect_length <- tail(image_numbers,1)-image_numbers[1]+1
transect_lats <- seq(dat.raw$Lat[5], dat.raw$Lat[6], length.out=transect_length)
transect_lons <- seq(dat.raw$Lon[5], dat.raw$Lon[6], length.out=transect_length)
dat.3 <- data.frame("FileName"=transect_images,"transectID"="KC5","GPS_lat"=NA,"GPS_lon"=NA, foldername=NBP1502_folders[3])
dat.3$GPS_lat <- transect_lats[image_positions]
dat.3$GPS_lon <- transect_lons[image_positions]

transect_images <- list.files(paste0(img.path,NBP1502_folders[4]),pattern=".jpg")
image_numbers <- as.numeric(gsub(".jpg"," ",gsub("DSC","",transect_images)))
image_positions <- image_numbers-image_numbers[1]+1
transect_length <- tail(image_numbers,1)-image_numbers[1]+1
transect_lats <- seq(dat.raw$Lat[7], dat.raw$Lat[8], length.out=transect_length)
transect_lons <- seq(dat.raw$Lon[7], dat.raw$Lon[8], length.out=transect_length)
dat.4 <- data.frame("FileName"=transect_images,"transectID"="KC6","GPS_lat"=NA,"GPS_lon"=NA, foldername=NBP1502_folders[4])
dat.4$GPS_lat <- transect_lats[image_positions]
dat.4$GPS_lon <- transect_lons[image_positions]

transect_images <- list.files(paste0(img.path,NBP1502_folders[5]),pattern=".jpg")
image_numbers <- as.numeric(gsub(".jpg"," ",gsub("DSC","",transect_images)))
image_positions <- image_numbers-image_numbers[1]+1
transect_length <- tail(image_numbers,1)-image_numbers[1]+1
transect_lats <- seq(dat.raw$Lat[9], dat.raw$Lat[10], length.out=transect_length)
transect_lons <- seq(dat.raw$Lon[9], dat.raw$Lon[10], length.out=transect_length)
dat.5 <- data.frame("FileName"=transect_images,"transectID"="KC12","GPS_lat"=NA,"GPS_lon"=NA, foldername=NBP1502_folders[5])
dat.5$GPS_lat <- transect_lats[image_positions]
dat.5$GPS_lon <- transect_lons[image_positions]

dat <- rbind(dat.1,dat.2,dat.3,dat.4,dat.5)
dat$Survey <- "NBP1502"
dat$transectID <- as.factor(dat$transectID)

#names(dat) <- c("FileName","Survey","Date","Time","transectID","Depth","GPS_lat","GPS_lon")

## overview
spatial.dat <- data.frame(cbind(dat$GPS_lon, dat$GPS_lat))
names(spatial.dat) <- c("Longitude","Latitude")
coordinates(spatial.dat) <- c("Longitude","Latitude")
proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
polar.dat <- spTransform(spatial.dat, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))


#####################################################################
##### 2. SUBSET IMAGES FROM IMAGE LOCATIONS (& STORE FILENAMES) #####
#####################################################################

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
  print("_________________________")
  print(levels(dat$transectID)[i])
  ## all data except for badly illuminated images
  subset.v <- which(dat$transectID==levels(dat$transectID)[i])
  dat.subset.v <- subset.v
  
  t.counts <- length(dat.subset.v)
  ## select first and last coordinate to calculate transect distance
  v.first <- dat.subset.v[1]
  v.last <- tail(dat.subset.v, n=1)
  t.l[i] <- max(distm(dat[subset.v,c(4,3)], fun=distVincentyEllipsoid))
  print(paste0("t.l: ",round(t.l[i],0),"m"))
  # total.t.length.v[i] <- t.l
  ## calculate distance of each image to the start of the transect along the transect line
  dat$dist.from.start[v.first] <- 0
  dat$dist.to.start[v.first] <- 0
  dat$dist.to.end[v.first] <- t.l[i]
  for(k in 2:t.counts){
    dist.to.previous.image <- distVincentyEllipsoid(dat[dat.subset.v[k-1],c(4,3)], dat[dat.subset.v[k],c(4,3)])
    dat$dist.from.start[dat.subset.v[k]] <- dat$dist.from.start[dat.subset.v[k-1]]+dist.to.previous.image
    dat$dist.to.start[dat.subset.v[k]] <- distVincentyEllipsoid(dat[v.first,c(4,3)], dat[dat.subset.v[k],c(4,3)])
    dat$dist.to.end[dat.subset.v[k]] <- distVincentyEllipsoid(dat[v.last,c(4,3)], dat[dat.subset.v[k],c(4,3)])
    }
  total.t.length.v[i] <- dat$dist.from.start[dat.subset.v[t.counts]]
  print(paste0("total.t.length.v: ",round(total.t.length.v[i],0),"m"))
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
  print(paste0("# of images in first batch: ", ceiling(total.t.length.v[i]/100)))
  #samp.v[(length(samp.v)+1):length(subset.v)] <- 9998
  dat$image.select[dat.subset.v[samp.v]] <-  1:length(samp.v)
  dat$image.select[dat.subset.v[-samp.v]] <-  999
  # ## choose a random subset of images in each transect by giving them random numbers starting from 1
  #dat$image.select[dat.subset.v] <- sample(1:length(dat.subset.v),length(dat.subset.v))
  samp.v.list[[i]] <- samp.v
}

## add dates:
library(lubridate)
dat$Date <- as.Date(NA)
unique(dat$transectID)
dat$Date[dat$transectID=="BC1"] <- ymd("2015:01:27")
dat$Date[dat$transectID=="KC3"] <- ymd("2015:02:03")
dat$Date[dat$transectID=="KC5"] <- ymd("2015:02:03")
dat$Date[dat$transectID=="KC6"] <- ymd("2015:02:10")
dat$Date[dat$transectID=="KC12"] <- ymd("2015:02:17")

## SAVE OUTPUT FOR FUTURE REFENCE (i.e. start here to add more images to the analysis)
#save(dat,total.t.length.v, t.l, file="C:/Users/jjansen/Desktop/science/SouthernOceanBiodiversityMapping/ARC_Data/prep_image/ReadIn_Circumpolar_DownwardImage_Data_And_Subset_Surveys/NBP1502_dat.Rdata")
#load(file="C:/Users/jjansen//Desktop/science/SouthernOceanBiodiversityMapping/ARC_Data/prep_image/ReadIn_Circumpolar_DownwardImage_Data_And_Subset_Surveys/NBP1502_dat.Rdata")

t.images <- ceiling(total.t.length.v/100)
t.images.full <- ceiling(t.l/100)
## 500 images in total
#t.images <- ceiling((total.t.length.v/total.t.length)*500)
## 

spatial.dat <- data.frame(dat[,c(4,3)])#[!is.na(rowSums(dat[,4:5])),]
coordinates(spatial.dat) <- c("GPS_lon","GPS_lat")
proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
polar.dat <- spTransform(spatial.dat, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

col.breaks <- seq(-100,-1400,length.out=101)
cols <- terrain.colors(99)

par(mfrow=c(1,1), mar=c(5,4,2,1))
for(i in 1:length(levels(dat$transectID))){
  subset.v <- which(dat$transectID==levels(dat$transectID)[i])
  sel <- order(dat$image.select[subset.v])[1:t.images[i]]
  if(any(is.na(sel))) sel <- sel[-which(is.na(sel))]
  ## create vector of selected images
  if(i==1){dat.sel <<-  subset.v[sel]
  }else {dat.sel <-  c(dat.sel,subset.v[sel])}
  
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

dat$FileName[dat.sel]
## filenames
selected.filenames <- dat$FileName[dat.sel]
selected.folders <- dat$foldername[dat.sel]
## create new filenames and identify folders where to find them (might be different for each survey)
selected.filenames.renamed <- paste0("NBP1502_",dat$transectID[dat.sel],"_",sprintf("%04d",dat$image.select[dat.sel]),"__",selected.filenames)

## copy files into Annotation folder
img.path.origin <- paste0(img.path2,selected.folders,"/",selected.filenames)
img.path.destin <- paste0("C:/Users/jjansen/Desktop/science/data_biological/Stills/Annotation_images/NBP1502/",selected.filenames.renamed)
file.copy(img.path.origin,img.path.destin)

## write list of filenames into cropped Annotation folder for upload
img.path.destin_crop <- "C:/Users/jjansen/Desktop/science/data_biological/Stills/Annotation_images_cropped/NBP1502/"
write.table(selected.filenames.renamed, paste0(img.path.destin_crop,"filenames.txt"), eol=",", col.names=FALSE, row.names=FALSE)

# filenames.in.folder_original <- paste0(img.path.destin,selected.filenames)
# filenames.in.folder_changed <- paste0(img.path.destin,selected.filenames.renamed)
# file.rename(filenames.in.folder_original,filenames.in.folder_changed)


# a <- list.files(img.path,pattern = "*.jpg")
# img.path.origin.good <- paste0(img.path,a[a%!in%bad_images])
# img.path.destin.good <- paste0(img.path,"TAN0802_good_images/")
# file.copy(img.path.origin.good,img.path.destin.good)
##

########################################################
##### 3. CROP SELECTED IMAGES AND STORE SEPARATELY #####
########################################################

## NBP1502

## @@@ !!! @@@ 
## crop images, do this in the terminal using the following code (after navigating into the "Stills" folder):
# cd C:
# cd C:\Users\jjansen\Desktop\science\data_biological\Stills\
# magick mogrify -path Annotation_images_cropped\NBP1502\ -gravity Center -crop 90% +repage Annotation_images\NBP1502\*jpg
## @@@ !!! @@@ 
















