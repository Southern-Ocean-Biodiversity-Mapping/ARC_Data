library(magick)
library(geosphere)
library(ggplot2)
library(MBHdesign)
library(raster)
library(lubridate)
library(raadtools)
'%!in%' <- function(x,y)!('%in%'(x,y))

sci.dir <-      "C:/Users/jjansen/Desktop/science/"
env.derived <-  paste0(sci.dir,"data_environmental/derived/")

r.path <- "R:/IMAS/Antarctic_Seafloor/Clean_Data_For_Permanent_Storage/PS18/PS18_1_raw_images_and_metadata/"

ps.path <- paste0(r.path,"images_colourcorrected/")
path.bad.images <- paste0(ps.path,"bad_quality/")

ps.path.img <- paste0(r.path,"images_original/")
ps.path.met <- paste0(r.path,"metadata/")

## metadata
dat.PS18 <- get(load(file=paste0(ps.path.met,"PS18_metadata.Rdata")))
## image names
#dir.files <- list.files(paste0(bio.path,"Stills/PS18/originals/"),pattern=".jpg")
dir.files <- list.files(ps.path.img,pattern=".jpg")

# ## environmental data for plotting
# my_data_dir <- "C:/Users/jjansen/Desktop/science/data_environmental/accessed_through_R"
# set_data_roots(my_data_dir)
# r <- readtopo("ibcso")
# r2 <- r
# r2[r2>0] <- NA
# r2[r2<(-2000)] <- NA
# ## coastline
# stereo <- "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
# coast.lonlat <- shapefile("C:/Users/jjansen/Desktop/science/data_environmental/antarctic coastline/addv5_sc_coast_ln_gg.shp")
# coast.proj <- spTransform(coast.lonlat, CRS(stereo))

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
files.dat[,8] <- as_datetime(NA)
files.dat[,9] <- as_datetime(NA)
files.dat[,10] <- NA
files.dat[,11] <- NA
names(files.dat)[6:11] <- c("lon","lat","time_start", "time_end","depth_start","depth_end")
for(i in c(1:23,25)){
  t.sel <- which(files.dat[,2]==levels(files.dat[,2])[i])
  img.sel <- files.dat[t.sel,5]
  files.dat[t.sel,8] <- ymd_hms(dat.PS18$`Date/Time_Start`[i])
  files.dat[t.sel,9] <- ymd_hms(dat.PS18$`Date/Time_End`[i])
  files.dat[t.sel,10] <- dat.PS18$Elevation_Start[i]
  files.dat[t.sel,11] <- dat.PS18$Elevation_End[i]
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
# save(PS18_image_metadata, file=paste0(ps.path.met,"PS18_image_metadata.Rdata"))

#####################################################################
##### 2. SUBSET IMAGES FROM IMAGE LOCATIONS (& STORE FILENAMES) #####
#####################################################################
## dataframe containing: coordinates, filename, transectID

## RANDOM NUMBERS IN THE SUBSETTING NOT FULLY REPRODUCIBLE!!! EVERYTHING ELSE IS...

## now remove bad quality images
bad_images <- list.files(path.bad.images)
dat.raw <- PS18_image_metadata[which(PS18_image_metadata$Filename%!in%bad_images),]
## remove transects 162, 173 & 207
# near ice-shelf without data: 162-1 (9)
# under ice-shelf without data: 173-2 (13)
# under ice-shelf without data: 207-1 (21)
dat <- dat.raw[-which(dat.raw$transectID=="162-1"|dat.raw$transectID=="173-2"|dat.raw$transectID=="207-1"),]
dat$transectID <- factor(dat$transectID)

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
  ## select first and last coordinate to calculate transect distance
  v.first <- dat.subset.v[1]
  v.last <- tail(dat.subset.v, n=1)
  # t.l <- distVincentyEllipsoid(dat[v.first,6:7], dat[v.last,6:7])
  # total.t.length.v[i] <- t.l
  ## calculate distance of each image to the start of the transect along the transect line
  dat$dist.from.start[v.first] <- 0
  for(k in 2:t.counts){
    dist.to.previous.image <- distVincentyEllipsoid(dat[dat.subset.v[k-1],6:7], dat[dat.subset.v[k],6:7])
    dat$dist.from.start[dat.subset.v[k]] <- dat$dist.from.start[dat.subset.v[k-1]]+dist.to.previous.image
  }
  total.t.length.v[i] <- dat$dist.from.start[dat.subset.v[t.counts]]
  #### subset images
    ## randomly select 5 images from the first 25 in transects that don't have good start/end points
  if(t.select%in%dat.PS18$TransectID[t.toofast]){
    set.seed(43)
    dat$image.select[dat.subset.v][1:25] <- sample(1:25)
    dat$image.select[dat.subset.v[26:t.counts]] <-  999
  }else{
    ## choose a spatial random subset of images in each transect using quasiSamp on the distance from the start of the transect
    set.seed(42)
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
}
## 
#dat$image.select[is.na(dat$image.select)] <- 9999

## SAVE OUTPUT FOR FUTURE REFENCE (i.e. start here to add more images to the analysis)
# save(dat,total.t.length.v, file="C:/Users/jjansen/Desktop/science/data_biological/PS18_dat.Rdata")
#load(file="C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/PS18_dat.Rdata")

barplot(round(total.t.length.v), names.arg=as.character(levels(dat$transectID)), las=2, main="PS18", xlab="TransectID", ylab="length in m")

## total transect length across all the survey
#total.t.length <- sum(total.t.length.v)
## how many images from each transect, if we select on average 1 every 50m?
t.images <- ceiling(total.t.length.v/100)
t.images[t.images==0] <- 5
## 500 images in total
#t.images <- ceiling((total.t.length.v/total.t.length)*500)
## 

spatial.dat <- data.frame(dat[,6:7])
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
  sel <- order(dat$image.select[subset.v])[1:t.images[i]]
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
  #plot(r2, xlim=xlim,ylim=ylim, main=paste0("transect ",levels(dat$transectID)[i]),col=cols, breaks=col.breaks)
  #plot(polar.dat[subset.v],add=TRUE)
  #points(polar.dat[subset.v[sel]], col="red", pch=15)
  #plot(coast.proj, add=TRUE)
  #scalebar(100, type="bar")
}

dat$Filename[dat.sel]
## filenames
selected.filenames <- dat$Filename[dat.sel]
## create new filenames and identify folders where to find them (might be different for each survey)
selected.filenames.renamed <- paste0(substr(dat$Filename[dat.sel],1,8),"_",sprintf("%04d",dat$image.select[dat.sel]),"__",selected.filenames)
selected.filenames.folders <- substr(selected.filenames,1,8)

## copy files into Annotation folder
img.path.origin <- paste0(image.dir,selected.filenames)
img.path.destin <- "C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/Stills/Annotation_images/PS18/"
file.copy(img.path.origin,img.path.destin)

## write list of filenames into cropped Annotation folder for upload
img.path.destin_crop <- "C:/Users/jjansen/Desktop/science/data_biological/Stills/Annotation_images_cropped/PS18/"
write.table(selected.filenames.renamed, paste0(img.path.destin_crop,"filenames.txt"), eol=",", col.names=FALSE, row.names=FALSE)

filenames.in.folder_original <- paste0(img.path.destin,selected.filenames)
filenames.in.folder_changed <- paste0(img.path.destin,selected.filenames.renamed)
file.rename(filenames.in.folder_original,filenames.in.folder_changed)

########################################################
##### 3. CROP SELECTED IMAGES AND STORE SEPARATELY #####
########################################################

## NO CROPPING


