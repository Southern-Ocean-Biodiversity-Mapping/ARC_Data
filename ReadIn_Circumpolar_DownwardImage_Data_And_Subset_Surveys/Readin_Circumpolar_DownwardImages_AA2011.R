#########################################################################################################################
## Mertz images: Finding lat/lon of each image, removing bad images and then subsetting to 1 image per station
#########################################################################################################################

## AA2011 ##

###############################################################################
#### libraries and functions
# library(dplyr)
# library(leaflet)
# library(lubridate)
# library(SOmap)
# #library(raadtools)
# library(spatialEco)
library(raster)
library(readxl)
# library(sp)
# #library(blueant)
# library(geosphere)
# library(MBHdesign)
'%!in%' <- function(x,y)!('%in%'(x,y))

#### specify paths
img.path <- txt.path <- "C:/Users/jjansen/Desktop/science/data_biological/Stills/AA2011/"

# img.path <- "D:/ARC_DP_data/a_RawData_DirectFromContributors/AA2011/"
# txt.path <- "D:/ARC_DP_data/a_RawData_DirectFromContributors/AA2011/"
# gps.txt.path <- "D:/ARC_DP_data/a_RawData_DirectFromContributors/AA2011/"

###############################################################################
## load metadata
dat <- read.csv(paste0(txt.path,"ctd_images.csv"))
names(dat)[2:7] <- c("transectID","transectID_folder","Filename","DateTime","lat","lon")
dat$FilePath <- paste0(dat$transectID_folder,"/",dat$Filename)

## create list with names (and paths) of good images only
img.names.infile <- list.files(paste0(img.path,"images_corrected"), recursive=TRUE, pattern=".jpg", ignore.case=TRUE)
good.images <- img.names.infile[-which(grepl("bad",img.names.infile,fixed = TRUE))]

## find gps coordinates for each image
dat2 <- dat[dat$FilePath%in%good.images,]
dat2$transectID <- as.factor(dat2$transectID)

## remove images that are too blurry
t.quality <- data.frame(read_excel(paste0(img.path,"Transect_quality.xlsx")))
names(t.quality)[1] <- "transectID_folder"
bad_transects <- which(t.quality$image.quality=="blurry")
dat3 <- dat2[-which(dat2$transectID_folder%in%t.quality$transectID_folder[bad_transects]),]
dat3$transectID <- as.factor(dat3$transectID)

## subset to 1 image per station
dat3$image.select <- NA
for(i in 1:length(levels(dat3$transectID))){
  transect <- levels(dat3$transectID)[i]
  sel <- which(dat3$transectID==transect)
  set.seed(42)
  dat3$image.select[sel] <- sample(1:length(sel),length(sel))
}

## save output
dat <- dat3
dat$transectID <- factor(dat$transectID)

#save(dat, file="C:/Users/jjansen/Desktop/science/data_biological/AA2011_dat.Rdata")
#load(file="C:/Users/jjansen/Desktop/science/data_biological/AA2011_dat.Rdata")

spatial.dat <- data.frame(dat[which(dat$image.select==1),c(7,6)])
coordinates(spatial.dat) <- c("lon","lat")
proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
polar.dat <- spTransform(spatial.dat, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

plot(polar.dat)
scalebar(100000, type="bar", label=c(0,50,100), below="km")

dat.sel <- which(dat$image.select==1)

## filenames
selected.filenames <- dat$Filename[dat.sel]
selected.filenames.path <- dat$FilePath[dat.sel]
selected.filenames.seq <- paste0(dat$transectID_folder[dat.sel],"__",dat$Filename[dat.sel])

## create new filenames and identify folders where to find them (might be different for each survey)
selected.filenames.renamed <- paste0("AA2011_",dat$transectID_folder[dat.sel],"_",sprintf("%04d",dat$image.select[dat.sel]),"__",selected.filenames)

## copy files into Annotation folder
img.path.origin <- paste0(img.path,"camera/",selected.filenames.path)
img.path.destin <- paste0("C:/Users/jjansen/Desktop/science/data_biological/Stills/Annotation_images/AA2011/",selected.filenames.seq)
file.copy(img.path.origin,img.path.destin)

## write list of filenames into cropped Annotation folder for upload
img.path.destin_crop <- "C:/Users/jjansen/Desktop/science/data_biological/Stills/Annotation_images_cropped/AA2011/"
write.table(selected.filenames.renamed, paste0(img.path.destin_crop,"filenames.txt"), eol=",", col.names=FALSE, row.names=FALSE)

filenames.in.folder_original <- img.path.destin
filenames.in.folder_changed <- paste0("C:/Users/jjansen/Desktop/science/data_biological/Stills/Annotation_images/AA2011/",selected.filenames.renamed)
file.rename(filenames.in.folder_original,filenames.in.folder_changed)

########################################################
##### 3. CROP SELECTED IMAGES AND STORE SEPARATELY #####
########################################################

## AA2011

## @@@ !!! @@@ 
## crop images, do this in the terminal using the following code (after navigating into the "Stills" folder):
C:
cd C:\Users\jjansen\Desktop\science\data_biological\Stills\
magick mogrify -path Annotation_images_cropped\AA2011\ -gravity Center -crop 90% +repage Annotation_images\AA2011\*jpg
## @@@ !!! @@@ 
















