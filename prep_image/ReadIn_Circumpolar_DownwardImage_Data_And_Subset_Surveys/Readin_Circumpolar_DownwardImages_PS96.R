#################################################
########## Prepare images for analysis ##########
#################################################

## !!! THIS SCRIPT REQUIRES THE USE OF BOTH R AND THE TERMINAL !!!
## command line code is highlighted with "@@@ !!! @@@"

library(magick)
library(geosphere)
library(ggplot2)
library(MBHdesign)
library(raster)

########## DO THIS SEPARATELY FOR EACH SURVEY ##########

########## PS96 ##########

## HOW ANNOYING, each dataset start with the data in a different row...
data.start <- c(26,26,26,25,25,24,25,25,25,26,25,26,25)

image.dir <- "D:/ARC_DP_data/PS96_Piepenburg2016/"
link <- "_links-to-photographs.tab"
pre.chars <- "PS96_"
transect.names <- c(
  "001-4","007-1","008-2","010-3","026-3","027-2","037-3","048-2","057-3","061-1","072-4","090-4","106-2"
)
dat.header <- names(read.table(paste0(image.dir,pre.chars,transect.names[1],link), skip=24, header=TRUE, sep="\t", stringsAsFactors=FALSE, fill = TRUE)[,c(1:8)])
dat <- read.table(paste0(image.dir,pre.chars,transect.names[1],link), skip=data.start[1]-1, header=FALSE, sep="\t", stringsAsFactors=FALSE, fill = TRUE)[,c(1:8)]
dat[,9] <- as.character(transect.names[1])
for(i in 2:length(transect.names)){
  dat.temp <- read.table(paste0(image.dir,pre.chars,transect.names[i],link), skip=data.start[i]-1, header=FALSE, sep="\t",stringsAsFactors=FALSE, fill = TRUE)[,c(1:8)]
  dat.temp[,9] <- as.character(transect.names[i])
  dat <- rbind(dat, dat.temp)
}
names(dat) <- c(dat.header[1:5], "Area","Type","Filename","transectID")
PS96.dat <- cbind(dat[,c(3,2,4,6,7)],NA,dat[,c(9,8,1,5)])
names(PS96.dat)[6] <- "SurveyID"
PS96.dat[6] <- "PS96"
PS96.dat$SurveyID <- as.factor(PS96.dat$SurveyID)
PS96.dat$transectID <- as.factor(PS96.dat$transectID)
PS96.dat$Type <- as.factor(PS96.dat$Type)

########################################################
##### 1. ADD INFO ABOUT ILLUMINATION OF THE IMAGES #####
########################################################

## @@@ !!! @@@ 
###run the following commands in the cmd line to store files containing brightness measures:
#cd C:\Users\jjansen\"OneDrive - University of Tasmania"\Desktop\science/data_biological/Stills
#magick convert PS96_Piepenburg2016\PS96_001\*.jpg -colorspace Gray -format "%[fx:quantumrange*image.mean]," info:- >PS96_Piepenburg2016\PS96_001_brightness.txt
#magick convert PS96_Piepenburg2016\PS96_007\*.jpg -colorspace Gray -format "%[fx:quantumrange*image.mean]," info:- >PS96_Piepenburg2016\PS96_007_brightness.txt
#magick convert PS96_Piepenburg2016\PS96_008\*.jpg -colorspace Gray -format "%[fx:quantumrange*image.mean]," info:- >PS96_Piepenburg2016\PS96_008_brightness.txt
#...
## @@@ !!! @@@
##THEN manually delete the last comma and add a new line

## read the file in and store the information
t.names <- basename(list.dirs(image.dir))[-1]
PS96.dat$brightness <- NA
for(i in 1:length(t.names)){
  message(i)
  message(length(list.files(paste0(image.dir,t.names[i]))))
  bright <- read.table(paste0(image.dir,t.names[i],"_brightness.txt"),sep=",")
  print(paste0(length(bright)," .txt-files"))
  sel <- which(PS96.dat$transectID==transect.names[i])
  print(paste0(length(sel)," dataframe"))
  PS96.dat$brightness[sel] <- unlist(bright)
  rm(bright)
}

## plot variability in brighness in each transect
ggplot(data = PS96.dat[,c(7,11)], aes(x=transectID, y=brightness)) + geom_boxplot(aes())
## define which images to exclude because of bad illumination
#bad.illumination <- which(PS96.dat$brightness<20000|PS96.dat$brightness>48000)

plot(PS96.dat$Height..m.,PS96.dat$brightness)

ggplot(data = PS96.dat[,c(7,10)], aes(x=transectID, y=Height..m.)) + geom_boxplot(aes())+ theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1))

#####################################################################
##### 2. SUBSET IMAGES FROM IMAGE LOCATIONS (& STORE FILENAMES) #####
#####################################################################

## images 1 and 2 in the first transect have the wrong geolocations, and so does image 2466, other imges don't show the seafloor
dat <- PS96.dat[-c(1,2,702,2311:2315,2468),]
dat[781,1:2] <- cbind(mean(dat[780,1],dat[782,1]),mean(dat[780,2],dat[782,2]))
## sort out badly illumianted images, calculate transect length, define how many images to select and give images a random number
dat$image.select <- NA
total.t.length.v <- NA
dat$dist.from.start <- NA
dat$mean.dist.to.10.nearest.images <- NA
dat$dist.from.center <- NA
samp.v.list <- list()
center.idx.v <- NA
for(i in 1:length(levels(dat$transectID))){
  print(levels(dat$t.ID_temp)[i])
  ## all data except for badly illuminated images
  subset.v <- which(dat$transectID==levels(dat$transectID)[i])
  # ## define which images to exclude because of bad illumination
  #bad.light <- which(dat$brightness[subset.v]<20000|dat$brightness[subset.v]>48000)
  bad.light <- which(dat$Height..m.[subset.v]<1|dat$Height..m.[subset.v]>2)
  if (length(bad.light>0)){
    dat.subset.v <- subset.v[-bad.light]
  } else dat.subset.v <- subset.v
  # dat.subset.v <- subset.v
  
  t.counts <- length(dat.subset.v)
  ## select first and last coordinate to calculate transect distance
  v.first <- dat.subset.v[1]
  v.last <- tail(dat.subset.v, n=1)
  # t.l <- distVincentyEllipsoid(dat[v.first,1:2], dat[v.last,1:2])
  # total.t.length.v[i] <- t.l
  ## calculate distance of each image to the start of the transect along the transect line
  dat$dist.from.start[v.first] <- 0
  for(k in 2:t.counts){
    dist.to.previous.image <- distVincentyEllipsoid(dat[dat.subset.v[k-1],1:2], dat[dat.subset.v[k],1:2])
    dat$dist.from.start[dat.subset.v[k]] <- dat$dist.from.start[dat.subset.v[k-1]]+dist.to.previous.image
  }
  total.t.length.v[i] <- dat$dist.from.start[dat.subset.v[t.counts]]

  ## choose a spatial random subset of images in each transect using quasiSamp on the distance from the start of the transect
  set.seed(42)
  samp <- quasiSamp(length(dat.subset.v)*5,dimension=1,potential.sites=dat$dist.from.start[dat.subset.v])
  samp.v <- samp$ID[which(duplicated(samp)==FALSE)]
  print(paste0("first double selection for image # ", which(duplicated(samp)==TRUE)[1]))
  print(paste0("# of images in first batch: ", ceiling(total.t.length.v[i]/50)))
  #samp.v[(length(samp.v)+1):length(subset.v)] <- 9998
  dat$image.select[dat.subset.v[samp.v]] <-  1:length(samp.v)
  dat$image.select[dat.subset.v[-samp.v]] <-  999
  dat$image.select[subset.v][bad.light] <- 999
  # ## choose a random subset of images in each transect by giving them random numbers starting from 1
  #dat$image.select[dat.subset.v] <- sample(1:length(dat.subset.v),length(dat.subset.v))
  samp.v.list[[i]] <- samp.v
}
## 
dat$image.select[is.na(dat$image.select)] <- 9999

## SAVE OUTPUT FOR FUTURE REFENCE (i.e. start here to add more images to the analysis)
# save(dat,total.t.length.v, file="C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/PS96_dat_2ndrun.Rdata")
#load(file="C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/PS96_dat.Rdata")

## total transect length across all the survey
#total.t.length <- sum(total.t.length.v)
## how many images from each transect, if we select on average 1 every 50m?
t.images <- ceiling(total.t.length.v/100)
## 500 images in total
#t.images <- ceiling((total.t.length.v/total.t.length)*500)
## 

spatial.dat <- data.frame(dat[,1:2])
coordinates(spatial.dat) <- c("Longitude","Latitude")
proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
polar.dat <- spTransform(spatial.dat, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

#par(mfrow=c(4,4), mar=c(5,4,2,1))
for(i in 1:length(levels(dat$transectID))){
  subset.v <- which(dat$transectID==levels(dat$transectID)[i])
  subset.v.good <- which(dat$Height..m.[subset.v]>=1&dat$Height..m.[subset.v]<=2)
  ## select the correct number of images for each transect
  sel <- order(dat$image.select[subset.v])[1:t.images[i]]
  ## create vector of selected images
  if(i==1){dat.sel <<-  subset.v[sel]
  }else dat.sel <-  c(dat.sel,subset.v[sel])
  plot(polar.dat[subset.v[subset.v.good]])
  points(polar.dat[subset.v[sel]], col="red", pch=15)
  #points(polar.dat[subset.v[center.idx.v[i]],1:2], col="green", pch=16)
  scalebar(500, type="bar")
  # plot(dat[subset.v,1:2])
  # points(dat[subset.v[sel],1:2], col="red", pch=15)
  # points(dat[subset.v[center.idx.v[i]],1:2], col="green", pch=16)
}

dat$Filename[dat.sel]
## filenames
selected.filenames <- paste0(dat$Filename[dat.sel],".jpg")
## create new filenames and identify folders where to find them (might be different for each survey)
selected.filenames.renamed <- paste0(substr(dat$Filename[dat.sel],1,8),"_",sprintf("%04d",dat$image.select[dat.sel]),"__",selected.filenames)
selected.filenames.folders <- substr(selected.filenames,1,8)

## copy files into Annotation folder
img.path.origin <- paste0(image.dir,selected.filenames.folders,"/",selected.filenames)
img.path.destin <- "C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/Stills/Annotation_images/PS96/"
file.copy(img.path.origin,img.path.destin)

## write list of filenames into cropped Annotation folder for upload
img.path.destin_crop <- "C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/Stills/Annotation_images_cropped/PS96/"
write.table(selected.filenames.renamed, paste0(img.path.destin_crop,"filenames.txt"), eol=",", col.names=FALSE, row.names=FALSE)

filenames.in.folder_original <- paste0(img.path.destin,selected.filenames)
filenames.in.folder_changed <- paste0(img.path.destin,selected.filenames.renamed)
file.rename(filenames.in.folder_original,filenames.in.folder_changed)

########################################################
##### 3. CROP SELECTED IMAGES AND STORE SEPARATELY #####
########################################################

## @@@ !!! @@@ 
## crop images, do this in the terminal using the following code (after navigating into the "Stills" folder):
# C:
# cd C:\Users\jjansen\"OneDrive - University of Tasmania"\Desktop\science/data_biological/Stills
#magick mogrify -path Annotation_images_cropped\PS96\ -gravity Center -crop 80% +repage Annotation_images\PS96\*jpg
## @@@ !!! @@@ 
