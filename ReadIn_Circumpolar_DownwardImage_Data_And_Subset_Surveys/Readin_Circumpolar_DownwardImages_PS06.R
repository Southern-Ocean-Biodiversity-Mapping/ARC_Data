#################################################
########## Prepare images for analysis ##########
#################################################

## !!! THIS SCRIPT REQUIRES THE USE OF BOTH R AND THE TERMINAL !!!
## command line code is highlighted with "@@@ !!! @@@"

## NOTE SOME GPS COORDINATES ARE OFF
## check using ship-speed and time to calculate distance from starting point

library(magick)
library(geosphere)
library(ggplot2)
library(MBHdesign)
library(raster)
library(lubridate)
'%!in%' <- function(x,y)!('%in%'(x,y))

# path.bad.images <- paste0(bio.path,"Stills/PS06/bad_quality/")
# bio.path <- "C:/Users/jjansen/Desktop/science/data_biological/"
# image.dir <- ...
# dat.PS06 <- get(load(file=paste0(bio.path,"Stills/PS06/PS06_metadata.Rdata")))

r.path <- "R:/IMAS/Antarctic_Seafloor/Clean_Data_For_Permanent_Storage/PS06/"
path.bad.images <- paste0(r.path,"PS06_1_raw_images_and_metadata/images_colourcorrected/bad_quality/")
image.dir <- paste0(r.path,"PS06_1_raw_images_and_metadata/images_original/")
dat.PS06 <- get(load(file=paste0(r.path,"PS06_1_raw_images_and_metadata/metadata/PS06_metadata.Rdata")))

## sort by descending transectID
dat.PS06 <- dat.PS06[order(substring(dat.PS06[,1],1,3)),]

##################################
######### PS06
##################################
#### load image names
#dir.files <- list.files(paste0(bio.path,"Stills/PS06/originals/"),pattern=".jpg")
dir.files <- list.files(image.dir,pattern=".jpg")
files.dat <- data.frame(matrix(ncol=3,nrow=length(dir.files)))
files.dat[,1] <- dir.files
files.dat[,2] <- substring(dir.files,6,8)
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

## assign median lon/lat and time to each image
files.dat[,6] <- NA
files.dat[,7] <- NA
files.dat[,8] <- as_datetime(NA)
files.dat[,9] <- as_datetime(NA)
files.dat[,10] <- NA
names(files.dat)[6:10] <- c("lon","lat","time_start","time_end","depth")
for(i in 1:ntransects){
  t.sel <- which(files.dat[,2]==levels(files.dat[,2])[i])
  img.sel <- files.dat[t.sel,5]
  files.dat[t.sel,6] <- dat.PS06$Median_Longitude[i]
  files.dat[t.sel,7] <- dat.PS06$Median_Latitude[i]
  files.dat[t.sel,8] <- ymd_hms(dat.PS06$`Date/Time_Start`[i])
  files.dat[t.sel,9] <- ymd_hms(dat.PS06$`Date/Time_End`[i])
  files.dat[t.sel,10] <- dat.PS06$Elevation[i]
}
PS06_image_metadata <- files.dat
#save(PS06_image_metadata, file=paste0(r.path,"PS06_1_raw_images_and_metadata/metadata/PS06_image_metadata.Rdata"))

#####################################################################
##### 2. SUBSET IMAGES FROM IMAGE LOCATIONS (& STORE FILENAMES) #####
#####################################################################
## dataframe containing: coordinates, filename, transectID

## now remove bad quality images
bad_images <- list.files(path.bad.images, pattern=".jpg")
dat.raw <- PS06_image_metadata[which(PS06_image_metadata$Filename%!in%bad_images),]
## remove transect 312 from the south atlantic
## also remove transect 288 because there is no environmental data for it
dat <- dat.raw[-which(dat.raw$transectID=="288"),]
dat$transectID <- factor(dat$transectID)

## remove transects with images that are too blurry
t.quality <- data.frame(read_excel(paste0(r.path,"PS06_1_raw_images_and_metadata/metadata/PS06_transect_quality.xlsx")))
bad_transects <- which(t.quality$quality=="blurry")
bad_transect.sel <- which(dat$transectID%in%t.quality$transectID[bad_transects])
dat2 <- dat[-bad_transect.sel,]
dat2$transectID <- factor(dat2$transectID)
dat <- dat2

## sort out badly illumianted images, calculate transect length, define how many images to select and give images a random number
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
  print(t.counts)
  #### subset images
  ## randomly select 5 images from each transect
  # set.seed(43)
  # dat$image.select[dat.subset.v][1:t.counts] <- sample(1:t.counts)
  ## choose a spatial random subset of images in each transect using quasiSamp on the distance from the start of the transect
  set.seed(42)
  samp <- quasiSamp(t.counts*5,dimension=1,potential.sites=1:t.counts)
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
# save(dat, file="C:/Users/jjansen/Desktop/science/data_biological/PS06_dat.Rdata")
#load(file="C:/Users/jjansen/Desktop/science/data_biological/PS06_dat.Rdata")

spatial.dat <- data.frame(dat[,6:7])
coordinates(spatial.dat) <- c("lon","lat")
proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
polar.dat <- spTransform(spatial.dat, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

col.breaks <- seq(-100,-1400,length.out=101)
cols <- terrain.colors(99)

# #### load depth
library(raadtools)
my_data_dir <- "C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_environmental/accessed_through_R"
set_data_roots(my_data_dir)
r <- readtopo("ibcso")
r2 <- r
r2[r2>0] <- NA
r2[r2<(-2000)] <- NA
# 
#### load coastline
stereo <- "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
coast.lonlat <- shapefile("C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_environmental/antarctic coastline/addv5_sc_coast_ln_gg.shp")
coast.proj <- spTransform(coast.lonlat, CRS(stereo))

par(mfrow=c(3,4), mar=c(5,4,2,1))
#par(mfrow=c(1,1), mar=c(5,4,2,1))
for(i in 1:length(levels(dat$transectID))){
  subset.v <- which(dat$transectID==levels(dat$transectID)[i])
  ## select the correct number of images for each transect
  sel <- order(dat$image.select[subset.v])[1:5]
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
  #plot(coast.proj, add=TRUE)
  scalebar(100, type="bar")
}


dat$Filename[dat.sel]
## filenames
selected.filenames <- dat$Filename[dat.sel]
## create new filenames and identify folders where to find them (might be different for each survey)
selected.filenames.renamed <- paste0(substr(dat$Filename[dat.sel],1,8),"_",sprintf("%04d",dat$image.select[dat.sel]),"__",selected.filenames)
selected.filenames.folders <- substr(selected.filenames,1,8)

## copy files into Annotation folder
img.path.origin <- paste0(image.dir,selected.filenames)
img.path.destin <- "C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/Stills/Annotation_images/PS06/"
file.copy(img.path.origin,img.path.destin)

## write list of filenames into cropped Annotation folder for upload
img.path.destin_crop <- "C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/Stills/Annotation_images_cropped/PS06/"
write.table(selected.filenames.renamed, paste0(img.path.destin_crop,"filenames.txt"), eol=",", col.names=FALSE, row.names=FALSE)

filenames.in.folder_original <- paste0(img.path.destin,selected.filenames)
filenames.in.folder_changed <- paste0(img.path.destin,selected.filenames.renamed)
file.rename(filenames.in.folder_original,filenames.in.folder_changed)

########################################################
##### 3. CROP SELECTED IMAGES AND STORE SEPARATELY #####
########################################################

## NO CROPPING

