## ARC-DP:
## Selecting additional images for annotation after initial subsetting
## 1. to replace shitty images
## 2. to increase the number of images scored

library(tidyverse)
library(raster)
'%!in%' <- function(x,y)!('%in%'(x,y))

##
load(file="C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/PS96_dat_FINAL.Rdata")
PS96.t.length <- total.t.length.v
load(file="C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/PS81_dat_FINAL.Rdata")
PS81.t.length <- total.t.length.v
rm(total.t.length.v)

## specify the path to the folder with the annotated images
img.path.PS96 <- "C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/Stills/Annotation_images_cropped/PS96/"
img.path.PS81 <- "C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/Stills/Annotation_images_cropped/PS81/"


#############################################################################
#############################################################################
#### PS96
#########
no.of.selected.images.CN <- ceiling(PS96.t.length/100)
no.of.selected.images.CN.planned <- ceiling(PS96.t.length/90)
no.of.selected.images.B <- ceiling(PS96.t.length/300)
no.of.selected.images.B.planned <- ceiling(PS96.t.length/200)

spatial.dat <- data.frame(PS96.dat[,1:2])
coordinates(spatial.dat) <- c("Longitude","Latitude")
proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
polar.dat <- spTransform(spatial.dat, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

for(i in 1:length(levels(PS96.dat$transectID))){
  t <- levels(PS96.dat$transectID)[i]
  ## all images in transect
  t.subset <- which(PS96.dat$transectID==t)
  dat.subset.full <- PS96.dat[t.subset,]
  ## annotated CN images within that transect
  ann.subset <- order(PS96.dat$image.select[PS96.dat$transectID==t])[1:no.of.selected.images.CN[i]]
  dat.subset <- dat.subset.full[ann.subset,]
  ## planned CN annotations within that transect
  ann.subset.planned <- order(PS96.dat$image.select[PS96.dat$transectID==t])[1:no.of.selected.images.CN.planned[i]]
  dat.subset.planned <- dat.subset.full[ann.subset.planned,]
  ## annotated BIIGLE images within that transect
  ann.subset.B <- order(PS96.dat$image.select[PS96.dat$transectID==t])[1:no.of.selected.images.B[i]]
  dat.subset.B <- dat.subset.full[ann.subset.B,]
  ## planned BIIGLE annotations within that transect
  ann.subset.B.planned <- order(PS96.dat$image.select[PS96.dat$transectID==t])[1:no.of.selected.images.B.planned[i]]
  dat.subset.B.planned <- dat.subset.full[ann.subset.B.planned,]
  ##
  message(paste0("TRANSECT ",t))
  names.list.in.dat <- dat.subset$Filename
  names.list.in.folder.raw <- list.files(img.path.PS96)[substr(list.files(img.path.PS96),21,25)==t]
  names.list.in.folder <- substr(names.list.in.folder.raw,16,41)
  message("images that still need to be added and annotated:")
  add.v <- names.list.in.dat[names.list.in.dat%!in%names.list.in.folder]
  print(add.v)
  message("delete:")
  delete.v <- names.list.in.folder[names.list.in.folder%!in%names.list.in.dat]
  print(delete.v)
  message("")
  message("list of all images that are marked for annotation in this transect")
  print(names.list.in.dat)
  message("")
  message("")
  ## store list of images to delete and replace for the survey
  if(i==1){
    add.list <<- add.v
    delete.list <<- delete.v
  }else {
    add.list <- c(add.list,add.v)
    delete.list <- c(delete.list,delete.v)
  }
  ## THE PLOTS ARE TO CHECK IF THE SUBSET OF BIIGLE ANNOTATIONS ARE SPATIALLY BALANCED
  plot(polar.dat[t.subset], pch=16) ## all images
  points(polar.dat[t.subset[ann.subset.planned]], col="orange", pch=16, cex=1.5) ## next CN images
  points(polar.dat[t.subset[ann.subset]], col="red", pch=16, cex=1.5) ## annotated CN images
  text(polar.dat[t.subset[ann.subset]], dat.subset$Filename, adj=c(1,-2))
  text(polar.dat[t.subset[ann.subset]], dat.subset$image.select, adj=c(-1,2))
  points(polar.dat[t.subset[ann.subset.B.planned]], col="skyblue", pch=16, cex=1.5) ## next images to annotate in Biigle
  points(polar.dat[t.subset[ann.subset.B]], col="royalblue", pch=16, cex=1.5) ## Biigle images
}

# ##########################################
# ####### ONCE ONLY
# ##### copy images and crop
# selected.files.idx <- which(PS96.dat$Filename%in%add.list)
# ## filenames
# selected.filenames <- paste0(add.list,".jpg")
# ## create new filenames and identify folders where to find them (might be different for each survey)
# selected.filenames.renamed <- paste0(substr(selected.filenames,1,8),"_",sprintf("%04d",PS96.dat$image.select[selected.files.idx]),"__",selected.filenames)
# selected.filenames.folders <- substr(selected.filenames,1,8)
# 
# ## copy files into Annotation folder
# image.dir <- "D:/ARC_DP_data/PS96_Piepenburg2016/"
# img.path.origin <- paste0(image.dir,selected.filenames.folders,"/",selected.filenames)
# img.path.destin <- "C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/Stills/Annotation_images/PS96ReplacementImagesOnceOnly/"
# file.copy(img.path.origin,img.path.destin)
# 
# ## write list of filenames into cropped Annotation folder for upload
# img.path.destin_crop <- "C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/Stills/Annotation_images_cropped/PS96ReplacementImagesOnceOnly/"
# write.table(selected.filenames.renamed, paste0(img.path.destin_crop,"additional_files_filenames.txt"), eol=",", col.names=FALSE, row.names=FALSE)
# 
# filenames.in.folder_original <- paste0(img.path.destin,selected.filenames)
# filenames.in.folder_changed <- paste0(img.path.destin,selected.filenames.renamed)
# file.rename(filenames.in.folder_original,filenames.in.folder_changed)
# 
# ## @@@ !!! @@@ 
# ## crop images, do this in the terminal using the following code (after navigating into the "Stills" folder):
# #cd C:\Users\jjansen\"OneDrive - University of Tasmania"\Desktop\science/data_biological/Stills
# #magick mogrify -path Annotation_images_cropped\PS96ReplacementImagesOnceOnly\ -gravity Center -crop 80% +repage Annotation_images\PS96ReplacementImagesOnceOnly\*jpg
# ## @@@ !!! @@@ 
# ##########################################

## BIIGLE-ANNOTATIONS LIST EDITS FOR HAVING approx. 1 image per 200m IN A BALANCED WAY
##
## 001:
## already: 1,2,4,5,6,7,8,
## next: 9:11
## SWAP 9 with 31
## SO: 10,11,31
##
## 007:
## already: 1:9
## next: 10:13
## SWAP 11,13 with 17,26
## SO: 10,17,26
##
## 008:
## already: 1:6
## next: 7:8
## SWAP 8 with 12
## SO: 7,12
##
## 010:
## already: 1:3,6:9
## next: 10:12
## SWAP 10,12 with 16,23
## SO: 11,16,23
##
## 026:
## already: 1:3,6:9
## next: 10:12
## SWAP 10 with 16
## SO: 11,12,16
##
## 027:
## already: 1:3
## next: 6:7
## SWAP NOTHING
## SO: 
##
## 037:
## already: 2,3,6,7
## next: 8
## SWAP NOTHING
## SO:
##
## 048:
## already: 1:3,6:10
## next: 11,12,13,16
## SWAP 11,16
## SO: 12,13,22,26
##
## 057:
## already: 1,2,6,7
## next: 8,9
## SWAP 9 with 15
## SO: 8,15
##
## 061:
## already: 1:3,6:8
## next: 9:11
## SWAP 9 with 20
## SO: 10,11,20
##
## 072:
## already: 1:3,6,7
## next: 8,9
## SWAP 9 with 16
## SO: 8,16
##
## 090:
## already: 1:3,6:8,10,11
## next: 12,13,15,16
## SWAP NOTHING
## SO:
##
## 106:
## already: 1:3,6
## next: 7
## SWAP 7 with 8
## SO: 8

#############################################################################
#############################################################################

#############################################################################
#############################################################################
#### PS81
#########
no.of.selected.images.CN <- ceiling(PS81.t.length/100)
no.of.selected.images.CN.planned <- ceiling(PS81.t.length/100)
no.of.selected.images.B <- ceiling(PS81.t.length/300)
no.of.selected.images.B.planned <- ceiling(PS81.t.length/200)

spatial.dat <- data.frame(PS81.dat[,1:2])
coordinates(spatial.dat) <- c("Longitude","Latitude")
proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
polar.dat <- spTransform(spatial.dat, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

for(i in 1:length(levels(PS81.dat$t.ID_temp))){
  t <- levels(PS81.dat$t.ID_temp)[i]
  ## all images in transect
  t.subset <- which(PS81.dat$t.ID_temp==t)
  dat.subset.full <- PS81.dat[t.subset,]
  ## annotated CN images within that transect
  ann.subset <- order(PS81.dat$image.select[PS81.dat$t.ID_temp==t])[1:no.of.selected.images.CN[i]]
  dat.subset <- dat.subset.full[ann.subset,]
  ## planned CN annotations within that transect
  ann.subset.planned <- order(PS81.dat$image.select[PS81.dat$t.ID_temp==t])[1:no.of.selected.images.CN.planned[i]]
  dat.subset.planned <- dat.subset.full[ann.subset.planned,]
  ## annotated BIIGLE images within that transect
  ann.subset.B <- order(PS81.dat$image.select[PS81.dat$t.ID_temp==t])[1:no.of.selected.images.B[i]]
  dat.subset.B <- dat.subset.full[ann.subset.B,]
  ## planned BIIGLE annotations within that transect
  ann.subset.B.planned <- order(PS81.dat$image.select[PS81.dat$t.ID_temp==t])[1:no.of.selected.images.B.planned[i]]
  dat.subset.B.planned <- dat.subset.full[ann.subset.B.planned,]
  ##
  message(paste0("TRANSECT ",t))
  names.list.in.dat <- dat.subset$Filename
  names.list.in.folder.raw <- list.files(img.path.PS81)[paste0("PS81_",substr(list.files(img.path.PS81),21,23))==t]
  if(t=="PS81_222a") {
    t.2="PS81_222-2a"
    names.list.in.folder.raw <- list.files(img.path.PS81)[paste0("PS81_",substr(list.files(img.path.PS81),21,26))==t.2]
  }
  if(t=="PS81_222b"){
    t.2="PS81_222-2b"
    names.list.in.folder.raw <- list.files(img.path.PS81)[paste0("PS81_",substr(list.files(img.path.PS81),21,26))==t.2]
  }
  names.list.in.folder <- substr(names.list.in.folder.raw,16,nchar(names.list.in.folder.raw))
  message("images that still need to be added and annotated:")
  add.v <- names.list.in.dat[names.list.in.dat%!in%names.list.in.folder]
  print(add.v)
  message("delete:")
  delete.v <- names.list.in.folder[names.list.in.folder%!in%names.list.in.dat]
  print(delete.v)
  message("")
  message("list of all images that are marked for annotation in this transect")
  #print(names.list.in.dat)
  message("")
  message("")
  ## store list of images to delete and replace for the survey
  if(i==1){
    add.list <<- add.v
    delete.list <<- delete.v
  }else {
    add.list <- c(add.list,add.v)
    delete.list <- c(delete.list,delete.v)
  }
  ## THE PLOTS ARE TO CHECK IF THE SUBSET OF BIIGLE ANNOTATIONS ARE SPATIALLY BALANCED
  plot(polar.dat[t.subset], pch=16) ## all images
  points(polar.dat[t.subset[ann.subset.planned]], col="orange", pch=16, cex=1.5) ## next CN images
  points(polar.dat[t.subset[ann.subset]], col="red", pch=16, cex=1.5) ## annotated CN images
  text(polar.dat[t.subset[ann.subset]], dat.subset$Filename, adj=c(1,-2))
  text(polar.dat[t.subset[ann.subset]], dat.subset$image.select, adj=c(-1,2))
  points(polar.dat[t.subset[ann.subset.B.planned]], col="skyblue", pch=16, cex=1.5) ## next images to annotate in Biigle
  points(polar.dat[t.subset[ann.subset.B]], col="royalblue", pch=16, cex=1.5) ## Biigle images
  scalebar(500)
}

##########################################
####### ONCE ONLY
##### copy images and crop
selected.files.idx <- which(PS81.dat$Filename%in%add.list)
## filenames
selected.filenames <- add.list
## create new filenames and identify folders where to find them (might be different for each survey)
selected.filenames.renamed <- paste0(substr(selected.filenames,1,8),"_",sprintf("%04d",PS81.dat$image.select[selected.files.idx]),"__",selected.filenames)

## copy files into Annotation folder
image.dir <- "D:/ARC_DP_data/PS81_Piepenburg2013/"
img.path.origin <- paste0(image.dir,selected.filenames)
img.path.destin <- "C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/Stills/Annotation_images/PS81ReplacementImagesOnceOnly/"
file.copy(img.path.origin,img.path.destin)

## write list of filenames into cropped Annotation folder for upload
img.path.destin_crop <- "C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/Stills/Annotation_images_cropped/PS81ReplacementImagesOnceOnly/"
write.table(selected.filenames.renamed, paste0(img.path.destin_crop,"additional_files_filenames.txt"), eol=",", col.names=FALSE, row.names=FALSE)

filenames.in.folder_original <- paste0(img.path.destin,selected.filenames)
filenames.in.folder_changed <- paste0(img.path.destin,selected.filenames.renamed)
file.rename(filenames.in.folder_original,filenames.in.folder_changed)

## @@@ !!! @@@
## crop images, do this in the terminal using the following code (after navigating into the "Stills" folder):
#cd C:\Users\jjansen\"OneDrive - University of Tasmania"\Desktop\science/data_biological/Stills
#magick mogrify -path Annotation_images_cropped\PS81ReplacementImagesOnceOnly\ -gravity Center -crop 80% +repage Annotation_images\PS81ReplacementImagesOnceOnly\*jpg
## @@@ !!! @@@
##########################################

## BIIGLE-ANNOTATIONS LIST EDITS FOR HAVING approx. 1 image per 200m IN A BALANCED WAY
##
## :116
## already: 1:8, 10:12,14,15
## next: 16:21
## SWAP 16,17,19,20
## SO: 18,21,34,39,44,65
##
## :118
## already: 2:5,7:15
## next: 16:22
## SWAP 16,17,21
## SO: 18:20,24,30,33
##
## :159
## already: 2:14
## next: 15,17:21
## SWAP 15,18:21
## SO: 17,24,25,31,33
##
## :160
## already: 1:8
## next: 9:12
## SWAP 10,11,12
## SO: 9,16,23,24
##
## :161
## already: 1:11
## next: 12:15,17,18
## SWAP 15,18
## SO:  12:14,17,22,23
##
## :163
## already: 1:13
## next: 14:19
## SWAP 16,18,19
## SO: 14,15,17,22,26,33
##
## :164
## already: 1:14
## next: 15:21
## SWAP: 17:21
## SO: 15,16,23,34,36,41,43

























## BIIGLE-ANNOTATIONS PRIORITY LIST FOR ACHIEVING approx. 1 image per 180m





## look-up image-names from annotation folder to check which images have been annotated
ceiling(PS96.t.length/100)


## which ones are the already annotated shitty images that need substitution?
PS96.dat[order(PS96.dat$image.select[PS96.dat$transectID==levels(PS96.dat$transectID)[1]])[1:19],c(8,11,15)]



## check total number of images required for each transect (increase this number if more images should be annotated)

## fill missing ones with new images from the list

## crop and colour correct new images

## add to list of images to annotate


