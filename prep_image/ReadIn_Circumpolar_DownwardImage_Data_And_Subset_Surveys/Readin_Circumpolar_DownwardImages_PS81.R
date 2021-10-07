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

#### PS81
## transect 159 super low vis
## need to split transect 222


#########################################
## INSERT AFTER I'VE DONE EVERYTHING: THIS IS TO ADD PS-AREA TO "PS81_dat_FINAL.Rdata"
met.dat <- read_excel("R:/IMAS/Antarctic_Seafloor/Clean_Data_For_Permanent_Storage/PS81/PS81_1_raw_images_and_metadata/metadata/140520_PS81_OFOS_allStnsalleFotos_MetaDaten_GIS.xls",skip=9)
met.dat$`Photo name` <- gsub("_at_","T",met.dat$`Photo name`)
met.dat$`Photo name` <- gsub("2013_01_","2013-01-",met.dat$`Photo name`)
met.dat$`Photo name` <- gsub("2013_02_","2013-02-",met.dat$`Photo name`)
met.dat$`Photo name` <- gsub("2013_03_","2013-03-",met.dat$`Photo name`)

load(file="C:/Users/jjansen/Desktop/science/SouthernOceanBiodiversityMapping/ARC_Data/ReadIn_Circumpolar_DownwardImage_Data_And_Subset_Surveys/PS81_dat_FINAL.Rdata")
#PS81.dat$Filename%in%met.dat$`Photo name`
chrs <- charmatch(PS81.dat$Filename,met.dat$`Photo name`)
PS81.dat$Area_from_metadata <- met.dat$`area (m^2)`[chrs]
#save(PS81.dat,total.t.length.v, file="C:/Users/jjansen/Desktop/science/data_biological/PS81_dat_FINAL.Rdata")

load(file="C:/Users/jjansen/Desktop/science/SouthernOceanBiodiversityMapping/ARC_Data/ReadIn_Circumpolar_DownwardImage_Data_And_Subset_Surveys/PS81_dat_shallow.Rdata")
chrs <- charmatch(dat$Filename,met.dat$`Photo name`)
dat$Area_from_metadata <- met.dat$`area (m^2)`[chrs]
#save(dat,total.t.length.v, file="C:/Users/jjansen/Desktop/science/data_biological/PS81_dat_shallow.Rdata")




#########################################




########## PS81 ##########
image.dir <- "D:/ARC_DP_data/a_RawData_DirectFromContributors/PS81_Piepenburg2013"
image.dir.adj <- "D:/ARC_DP_data/adjusted_PS81/bad_quality/"
dat.files.path <- paste0(image.dir,"/",list.files(image.dir, pattern=".tab"))
dat.imagenames <- list.files(image.dir, pattern=".jpg")

## 
data.start <- 26

dat.header <- names(read.table(dat.files.path[1], skip=24, header=TRUE, sep="\t", stringsAsFactors=FALSE, fill = TRUE)[,c(1:8)])
dat <- read.table(dat.files.path[1], skip=data.start-1, header=FALSE, sep="\t", stringsAsFactors=FALSE, fill = TRUE)[,c(1:8)]
dat[,9] <- substr(dat[,8],1,8)
for(i in 2:length(dat.files.path)){
  dat.temp <- read.table(dat.files.path[i], skip=data.start-1, header=FALSE, sep="\t",stringsAsFactors=FALSE, fill = TRUE)[,c(1:8)]
  dat.temp[,9] <- substr(dat.temp[,8],1,8)
  dat <- rbind(dat, dat.temp)
}
names(dat) <- c(dat.header[1:5], "Area","Type","Filename","transectID")
dat$Area <- NA
dat$Type <- NA
PS81.dat <- cbind(dat[,c(3,2,4,6,7)],NA,dat[,c(9,8,1,5)])
names(PS81.dat)[6] <- "SurveyID"
PS81.dat[6] <- "PS81"
PS81.dat$SurveyID <- as.factor(PS81.dat$SurveyID)
PS81.dat$transectID <- as.factor(PS81.dat$transectID)
PS81.dat$Type <- as.factor(PS81.dat$Type)

########################################################
##### 1. ADD INFO ABOUT ILLUMINATION OF THE IMAGES #####
########################################################

## @@@ !!! @@@ 
###run the following command in the cmd line to store files containing brightness measures:
#cd D:\ARC_DP_data\PS81_Piepenburg2013
#for %I in (D:\ARC_DP_data\PS81_Piepenburg2013\*.jpg) do magick convert %I -colorspace Gray -format "%[fx:quantumrange*image.mean]\n" info:- > %~nI.txt
#magick convert *.jpg -colorspace Gray -format "%[fx:quantumrange*image.mean]," info:- >PS81_brightness.txt
## @@@ !!! @@@
##
# txt.filenames <- list.files(image.dir, pattern=".jpg.txt")
# txt.filenames.renamed <- paste0(substr(txt.filenames,1,35),".txt")
# filenames.in.folder_original <- paste0(image.dir,"/",txt.filenames)
# filenames.in.folder_changed <- paste0(image.dir,"/",txt.filenames.renamed)
# file.rename(filenames.in.folder_original,filenames.in.folder_changed)

## move files into separate folder, then:
## read the file in and store the information
# txt.dir <- "D:/ARC_DP_data/PS81_Piepenburg2013/brightness_files"
# txt.imagenames <- list.files(txt.dir, pattern=".txt")
# txt.imagenames.path <- paste0(txt.dir,"/",txt.imagenames)
# bright <- NA
# for(i in 1:length(txt.imagenames)){
#   bright[i] <- unlist(read.table(txt.imagenames.path[i]))
# }
# PS81.dat$brightness <- bright

# ## plot variability in brighness in each transect
# ggplot(data = PS81.dat[,c(7,11)], aes(x=transectID, y=brightness)) + geom_boxplot(aes())+ theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1))
# ## define which images to exclude because of bad illumination
# #bad.illumination <- which(PS96.dat$brightness<20000|PS96.dat$brightness>48000)
# 
# plot(PS81.dat$Dist.subs..m.,PS81.dat$brightness)

ggplot(data = PS81.dat[,c(7,10)], aes(x=transectID, y=Dist.subs..m.)) + geom_boxplot(aes())+ theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1))

#####################################################################
##### 2. SUBSET IMAGES FROM IMAGE LOCATIONS (& STORE FILENAMES) #####
#####################################################################
## remove images that are shallower than 100m
dat <- PS81.dat[-which(PS81.dat$Bathy.depth..m.<=100),]
dat$t.ID_temp <- as.character(dat$transectID)
dat$t.ID_temp[which(dat$transectID=="PS81_222")[(1:323)]] <- "PS81_222a"
dat$t.ID_temp[which(dat$transectID=="PS81_222")[(324:872)]] <- "PS81_222b"
dat$t.ID_temp <- as.factor(dat$t.ID_temp)

## sort out badly illumianted images, calculate transect length, define how many images to select and give images a random number
dat$image.select <- NA
total.t.length.v <- NA
dat$dist.from.start <- NA
dat$mean.dist.to.10.nearest.images <- NA
dat$dist.from.center <- NA
samp.v.list <- list()
center.idx.v <- NA
ends <- 1
for(i in 1:length(levels(dat$t.ID_temp))){
  print(levels(dat$t.ID_temp)[i])
  ## all data except for badly illuminated images
  subset.v <- which(dat$t.ID_temp==levels(dat$t.ID_temp)[i])
  ## define which images to exclude because of bad illumination
  #bad.light <- which(dat$brightness[subset.v]<7000|dat$brightness[subset.v]>48000)
  bad.light <- which(dat$Dist.subs..m.[subset.v]<1|dat$Dist.subs..m.[subset.v]>2)
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
  # ## calculate distance of each image to the start of the transect along the transect line
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

#save(dat,total.t.length.v,samp.v.list, file="C:/Users/jjansen/Desktop/science/data_biological/PS81_dat_2ndrun.Rdata")
#load(file="C:/Users/jjansen/Desktop/science/data_biological/PS81_dat_FINAL.Rdata")

## total transect length across all the survey
#total.t.length <- sum(total.t.length.v)
## how many images from each transect, if we select on average 1 every 50m?
t.images <- ceiling(total.t.length.v/100)
## 500 images in total
#t.images <- ceiling((total.t.length.v/total.t.length)*500)
##

par(mfrow=c(1,1))
# for(i in 1:length(levels(dat$t.ID_temp))){
#   plot(dat$mean.dist.to.10.nearest.images[dat$t.ID_temp==levels(dat$t.ID_temp)[i]])
# }

spatial.dat <- data.frame(dat[,1:2])
coordinates(spatial.dat) <- c("Longitude","Latitude")
proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
polar.dat <- spTransform(spatial.dat, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

#par(mfrow=c(4,4), mar=c(5,4,2,1))
for(i in 1:length(levels(dat$t.ID_temp))){
  subset.v <- which(dat$t.ID_temp==levels(dat$t.ID_temp)[i])
  subset.v.good <- which(dat$Dist.subs..m.[subset.v]>=1&dat$Dist.subs..m.[subset.v]<=2)
  ## select the correct number of images for each transect
  sel <- order(dat$image.select[subset.v])[1:t.images[i]]
  ## create vector of selected images
  if(i==1){dat.sel <<-  subset.v[sel]
  }else dat.sel <-  c(dat.sel,subset.v[sel])
  plot(polar.dat[subset.v[subset.v.good]], main=levels(dat$t.ID_temp)[i])
  points(polar.dat[subset.v[sel]], col="red", pch=15)
  #points(polar.dat[subset.v[center.idx.v[i]],1:2], col="green", pch=16)
  scalebar(500, type="bar")
  # plot(dat[subset.v,1:2])
  # points(dat[subset.v[sel],1:2], col="red", pch=15)
  # points(dat[subset.v[center.idx.v[i]],1:2], col="green", pch=16)
}

dat$Filename[dat.sel]
## filenames
selected.filenames <- paste0(dat$Filename[dat.sel])
## create new filenames and identify folders where to find them (might be different for each survey)
selected.filenames.renamed <- paste0(substr(dat$Filename[dat.sel],1,8),"_",sprintf("%04d",dat$image.select[dat.sel]),"__",selected.filenames)
selected.filenames.folders <- substr(selected.filenames,1,8)

## copy files into Annotation folder
img.path.origin <- paste0(image.dir,"/",selected.filenames)
img.path.destin <- "C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/Stills/Annotation_images/PS81/"
file.copy(img.path.origin,img.path.destin)

## write list of filenames into cropped Annotation folder for upload
img.path.destin_crop <- "C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/Stills/Annotation_images_cropped/PS81/"
write.table(selected.filenames.renamed, paste0(img.path.destin_crop,"PS81_filenames1.txt"), eol=",", col.names=FALSE, row.names=FALSE)

filenames.in.folder_original <- paste0(img.path.destin,selected.filenames)
filenames.in.folder_changed <- paste0(img.path.destin,selected.filenames.renamed)
file.rename(filenames.in.folder_original,filenames.in.folder_changed)

########################################################
##### 3. CROP SELECTED IMAGES AND STORE SEPARATELY #####
########################################################

## @@@ !!! @@@ 
## crop images, do this in the terminal using the following code (after navigating into the "Stills" folder):
#cd C:\Users\jjansen\"OneDrive - University of Tasmania"\Desktop\science/data_biological/Stills
#magick mogrify -path Annotation_images_cropped\PS81\ -gravity Center -crop 80% +repage Annotation_images\PS81\*jpg
## @@@ !!! @@@ 





################################################################################################################
##### IMAGES SHALLOWER THAN 100m DEPTH #########################################################################
################################################################################################################
dat <- PS81.dat[which(PS81.dat$Bathy.depth..m.<=100),]
dat.NA <- PS81.dat[which(is.na(PS81.dat$Bathy.depth..m.)),]
dat$t.ID_temp <- as.character(dat$transectID)
dat$t.ID_temp <- as.factor(dat$t.ID_temp)

## remove bad images. Hadn't done this before for the deeper dataset
bad_images <- list.files(image.dir.adj)
filenames <- dat$Filename
idx <- which(filenames%in%bad_images)
dat <- dat[-idx,]

l <- levels(dat$t.ID_temp)
par(mfrow=c(2,2))
for(i in 1:length(l)){
  plot(PS81.dat[PS81.dat$transectID==l[i],1:2], main=l[i])
  points(dat[dat$transectID==l[i],1:2], col="red")
  points(dat.NA[dat.NA$transectID==l[i],1:2], col="green")
}
par(mfrow=c(1,1))

## sort out badly illumianted images, calculate transect length, define how many images to select and give images a random number
dat$image.select <- NA
total.t.length.v <- NA
dat$dist.from.start <- NA
dat$mean.dist.to.10.nearest.images <- NA
dat$dist.from.center <- NA
samp.v.list <- list()
center.idx.v <- NA
ends <- 1
for(i in 1:length(levels(dat$t.ID_temp))){
  print(levels(dat$t.ID_temp)[i])
  ## all data except for badly illuminated images
  subset.v <- which(dat$t.ID_temp==levels(dat$t.ID_temp)[i])
  ## define which images to exclude because of bad illumination
  #bad.light <- which(dat$brightness[subset.v]<7000|dat$brightness[subset.v]>48000)
  bad.light <- which(dat$Dist.subs..m.[subset.v]<1|dat$Dist.subs..m.[subset.v]>2)
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
  # ## calculate distance of each image to the start of the transect along the transect line
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
dat$transectID <- factor(dat$transectID)
#save(dat,total.t.length.v,samp.v.list, file="C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/PS81_dat_shallow.Rdata")
#load(file="C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/PS81_dat_shallow.Rdata")

## total transect length across all the survey
#total.t.length <- sum(total.t.length.v)
## how many images from each transect, if we select on average 1 every 50m?
t.images <- ceiling(total.t.length.v/100)
## 500 images in total
#t.images <- ceiling((total.t.length.v/total.t.length)*500)
##

par(mfrow=c(1,1))
# for(i in 1:length(levels(dat$t.ID_temp))){
#   plot(dat$mean.dist.to.10.nearest.images[dat$t.ID_temp==levels(dat$t.ID_temp)[i]])
# }

spatial.dat <- data.frame(dat[,1:2])
coordinates(spatial.dat) <- c("Longitude","Latitude")
proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
polar.dat <- spTransform(spatial.dat, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

#par(mfrow=c(4,4), mar=c(5,4,2,1))
for(i in 1:length(levels(dat$t.ID_temp))){
  subset.v <- which(dat$t.ID_temp==levels(dat$t.ID_temp)[i])
  subset.v.good <- which(dat$Dist.subs..m.[subset.v]>=1&dat$Dist.subs..m.[subset.v]<=2)
  ## select the correct number of images for each transect
  sel <- order(dat$image.select[subset.v])[1:t.images[i]]
  ## create vector of selected images
  if(i==1){dat.sel <<-  subset.v[sel]
  }else dat.sel <-  c(dat.sel,subset.v[sel])
  plot(polar.dat[subset.v[subset.v.good]], main=levels(dat$t.ID_temp)[i])
  points(polar.dat[subset.v[sel]], col="red", pch=15)
  #points(polar.dat[subset.v[center.idx.v[i]],1:2], col="green", pch=16)
  scalebar(100, type="bar")
  # plot(dat[subset.v,1:2])
  # points(dat[subset.v[sel],1:2], col="red", pch=15)
  # points(dat[subset.v[center.idx.v[i]],1:2], col="green", pch=16)
}

dat$Filename[dat.sel]
## filenames
selected.filenames <- paste0(dat$Filename[dat.sel])
## create new filenames and identify folders where to find them (might be different for each survey)
selected.filenames.renamed <- paste0(substr(dat$Filename[dat.sel],1,8),"_",sprintf("%04d",dat$image.select[dat.sel]),"__",selected.filenames)
selected.filenames.folders <- substr(selected.filenames,1,8)

## copy files into Annotation folder
img.path.origin <- paste0(image.dir,"/",selected.filenames)
img.path.destin <- "C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/Stills/Annotation_images/PS81_shallow/"
file.copy(img.path.origin,img.path.destin)

## write list of filenames into cropped Annotation folder for upload
img.path.destin_crop <- "C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/Stills/Annotation_images_cropped/PS81_shallow/"
write.table(selected.filenames.renamed, paste0(img.path.destin_crop,"PS81_shallow_filenames.txt"), eol=",", col.names=FALSE, row.names=FALSE)

filenames.in.folder_original <- paste0(img.path.destin,selected.filenames)
filenames.in.folder_changed <- paste0(img.path.destin,selected.filenames.renamed)
file.rename(filenames.in.folder_original,filenames.in.folder_changed)

## @@@ !!! @@@ 
## crop images, do this in the terminal using the following code (after navigating into the "Stills" folder):
#cd C:\Users\jjansen\"OneDrive - University of Tasmania"\Desktop\science/data_biological/Stills
#magick mogrify -path Annotation_images_cropped\PS81_shallow\ -gravity Center -crop 80% +repage Annotation_images\PS81_shallow\*jpg
## @@@ !!! @@@ 




