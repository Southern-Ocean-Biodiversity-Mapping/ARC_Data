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
'%!in%' <- function(x,y)!('%in%'(x,y))

image.dir <- "D:/ARC_DP_data/a_RawData_DirectFromContributors/PS118/"
d <- "R:/IMAS/Antarctic_Seafloor/Clean_Data_For_Permanent_Storage/PS118/PS118_1_raw_images_and_metadata/metadata/"
dest <- "R:/IMAS/Antarctic_Seafloor/Clean_Data_For_Permanent_Storage/PS118/PS118_1_raw_images_and_metadata/images_original/"
env.dir <- "C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_environmental/"

# #### load depth
r2 <- raster(paste0(env.dir,"Circumpolar_EnvData_bathy500m_shelf_gebco2020_depth.grd"))

#### load coastline
stereo <- "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
# coast.lonlat <- shapefile("C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_environmental/antarctic coastline/addv5_sc_coast_ln_gg.shp")
# coast.proj <- spTransform(coast.lonlat, CRS(stereo))
load(paste0(env.dir,"Circumpolar_Coastline.Rdata"))

########## DO THIS SEPARATELY FOR EACH SURVEY ##########

########## PS118 ##########
t.ID <- c("6-9","7-1","8-1","9-1","11-2","12-12","38-13","39-1","69-1","77-1","81-1")
tab1 <- read.table(paste0(d,"PS118_6-9_links-to-images.tab"),skip=21, fill=TRUE, header=TRUE,sep="\t")
tab2 <- read.table(paste0(d,"PS118_7-1_links-to-images.tab"),skip=21, fill=TRUE, header=TRUE,sep="\t")
tab3 <- read.table(paste0(d,"PS118_8-1_links-to-images.tab"),skip=22, fill=TRUE, header=TRUE,sep="\t")
tab4 <- read.table(paste0(d,"PS118_9-1_links-to-images.tab"),skip=22, fill=TRUE, header=TRUE,sep="\t")
tab5 <- read.table(paste0(d,"PS118_11-2_links-to-images.tab"),skip=22, fill=TRUE, header=TRUE,sep="\t")
tab6 <- read.table(paste0(d,"PS118_12-12_links-to-images.tab"),skip=22, fill=TRUE, header=TRUE,sep="\t")
tab7 <- read.table(paste0(d,"PS118_38-13_links-to-images.tab"),skip=22, fill=TRUE, header=TRUE,sep="\t")
tab8 <- read.table(paste0(d,"PS118_39-1_links-to-images.tab"),skip=22, fill=TRUE, header=TRUE,sep="\t")
tab9 <- read.table(paste0(d,"PS118_69-1_links-to-images.tab"),skip=22, fill=TRUE, header=TRUE,sep="\t")
tab10 <- read.table(paste0(d,"PS118_77-1_links-to-images.tab"),skip=22, fill=TRUE, header=TRUE,sep="\t")
tab11 <- read.table(paste0(d,"PS118_81-1_links-to-images.tab"),skip=22, fill=TRUE, header=TRUE,sep="\t")
tab10[,2:3] <- tab10[,2:3]*-1

tab1$transectID <- t.ID[1]
tab2$transectID <- t.ID[2]
tab3$transectID <- t.ID[3]
tab4$transectID <- t.ID[4]
tab5$transectID <- t.ID[5]
tab6$transectID <- t.ID[6]
tab7$transectID <- t.ID[7]
tab8$transectID <- t.ID[8]
tab9$transectID <- t.ID[9]
tab10$transectID <- t.ID[10]
tab11$transectID <- t.ID[11]

tab <- rbind(tab1,tab2,tab3,tab4,tab5,tab6,tab7,tab8,tab9,tab10,tab11)
tab$transectID <- as.factor(tab$transectID)
tab$Area <- NA
tab$SurveyID <- as.factor("PS118")
tab$File.name <- as.character(tab$File.name)
tab$Date.Time <- as.character(tab$Date.Time)
names(tab)[6] <- "Filename"
tab.reordered <- tab[,c(3,2,4,11,5,12,10,6,1)]

PS118.dat <- tab.reordered

#####################################################################
##### 2. SUBSET IMAGES FROM IMAGE LOCATIONS (& STORE FILENAMES) #####
#####################################################################
## identify bad quality images
#pth <- "R:/IMAS/Antarctic_Seafloor/Clean_Data_For_Permanent_Storage/PS118/PS118_1_raw_images_and_metadata/images_colourcorrected/"
pth <- "D:/ARC_DP_data/adjusted_PS118/"
bad.images1 <- list.files(paste0(pth,"PS118_06/bad_quality/"))
bad.images2 <- list.files(paste0(pth,"PS118_07/bad_quality/"))
bad.images3 <- list.files(paste0(pth,"PS118_08/bad_quality/"))
bad.images4 <- list.files(paste0(pth,"PS118_09/bad_quality/"))
bad.images5 <- list.files(paste0(pth,"PS118_11/bad_quality/"))
bad.images6 <- list.files(paste0(pth,"PS118_12/bad_quality/"))
bad.images7 <- list.files(paste0(pth,"PS118_38/bad_quality/"))
bad.images8 <- list.files(paste0(pth,"PS118_39/bad_quality/"))
bad.images9 <- list.files(paste0(pth,"PS118_69/bad_quality/"))
bad.images10 <- list.files(paste0(pth,"PS118_77/bad_quality/"))
bad.images11 <- list.files(paste0(pth,"PS118_81/bad_quality/"))
bad_images <- c(bad.images1,bad.images2,bad.images3,bad.images4,bad.images5,bad.images6,bad.images7,bad.images8,bad.images9,bad.images10,bad.images11)
## identify images that are shallower than 100m
## KEEP THEM!!!
#bad_depth <- paste0(PS118.dat$Filename[which(PS118.dat$Depth.water..m.<=100)],".jpg")
## identify images taken using the hotkey rather than the timer
evil_hotkey <- paste0(PS118.dat$Filename[which(PS118.dat$Type=="HOTKEY")],".jpg")

## remove all of these bad images
filenames <- paste0(PS118.dat$Filename,".jpg")
#idx <- which(filenames%in%bad_images|filenames%in%bad_depth|filenames%in%evil_hotkey)
idx <- which(filenames%in%bad_images|filenames%in%evil_hotkey)
dat <- PS118.dat[-idx,]

## remove transect 38
dat <- dat[-which(dat$transectID=="38-13"),]

## some images without lon/lat:
dat <- dat[-which(is.na(rowSums(dat[,1:2]))),]

dat$transectID <- factor(dat$transectID)

## plot all transects
spatial.dat <- data.frame(dat[,1:2])
coordinates(spatial.dat) <- c("Longitude","Latitude")
proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
polar.dat <- spTransform(spatial.dat, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
par(mfrow=c(4,3), mar=c(3,2,2,1),oma=c(0,0,0,2))
for(i in 1:length(levels(dat$transectID))){
  subset.v <- which(dat$transectID==levels(dat$transectID)[i])
  xlim <- extent(polar.dat[subset.v])[1:2]
  xlim[1] <- xlim[1]-500
  xlim[2] <- xlim[2]+500
  ylim <- extent(polar.dat[subset.v])[3:4]
  ylim[1] <- ylim[1]-500
  ylim[2] <- ylim[2]+500
  plot(r2, xlim=xlim, ylim=ylim,main=levels(dat$transectID)[i])
  points(polar.dat[subset.v])
  #text(polar.dat[subset.v])
  scalebar(500, type="bar")
}

#### replace bad coordinates either with NAs or by interpolation
## shit image locations for interpolation
interpl.and.repl <- function(a,b,t.ID){
  subset.v <- which(dat$transectID==t.ID)
  start.idx <- a-1
  end.idx <- b+1
  N <- length(a:b)
  dat[subset.v,1][a:b] <<- seq(dat[subset.v,1][start.idx],dat[subset.v,1][end.idx],length.out=N+2)[seq(2,N+1)]
  dat[subset.v,2][a:b] <<- seq(dat[subset.v,2][start.idx],dat[subset.v,2][end.idx],length.out=N+2)[seq(2,N+1)]
}

#### transect 39 causes problems!!!
## which images are further away than the sum of distance between the 5 previous images
dat.subset.v <- which(dat$transectID=="39-1")
t.counts <- length(dat.subset.v)
## select first and last coordinate to calculate transect distance
v.first <- dat.subset.v[1]
v.last <- tail(dat.subset.v, n=1)
## calculate distance of each image to the start of the transect along the transect line
dist.from.start <- NA
dist.from.start[v.first] <- 0
for(k in 2:t.counts){
  dist.to.previous.image <- distVincentyEllipsoid(dat[dat.subset.v[k-1],1:2], dat[dat.subset.v[k],1:2])
  dist.from.start[dat.subset.v[k]] <- dist.from.start[dat.subset.v[k-1]]+dist.to.previous.image
}
dist.from.start[dat.subset.v]
dist.to.previous <- dist.from.start[dat.subset.v]-dist.from.start[dat.subset.v-1]
red.points <- which(dist.to.previous>100)

## plot transect 39:
spatial.dat <- data.frame(dat[,1:2])
coordinates(spatial.dat) <- c("Longitude","Latitude")
proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
polar.dat <- spTransform(spatial.dat, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
par(mfrow=c(1,1), mar=c(3,2,2,1),oma=c(0,0,0,2))
for(i in 3){
  subset.v <- which(dat$transectID==levels(dat$transectID)[i])
  xlim <- extent(polar.dat[subset.v])[1:2]
  xlim[1] <- xlim[1]+3000
  xlim[2] <- xlim[2]-4500
  ylim <- extent(polar.dat[subset.v])[3:4]
  ylim[1] <- ylim[1]
  ylim[2] <- ylim[2]-3000
  plot(r2, xlim=xlim, ylim=ylim,main=levels(dat$transectID)[i])
  #points(polar.dat[subset.v])
  text(polar.dat[subset.v],cex=0.5,col="grey")
  text(polar.dat[subset.v][red.points],labels=red.points,cex=0.8, col="red")
  scalebar(500, type="bar")
}

interpl.and.repl(1554,1555,"39-1") #1832 - 244
interpl.and.repl(1563,1563,"39-1")#2073,2076
interpl.and.repl(1598,1603,"39-1")
interpl.and.repl(1624,1626,"39-1")
interpl.and.repl(1636,1642,"39-1")
interpl.and.repl(1663,1676,"39-1")
interpl.and.repl(1692,1696,"39-1")
interpl.and.repl(1703,1703,"39-1")
interpl.and.repl(2088,2095,"39-1")
interpl.and.repl(1428,1429,"39-1")

#### transect 6-9 has werid image locations as well
dat.subset.v <- which(dat$transectID=="6-9")
plot(dat[dat.subset.v,1:2])
text(dat[dat.subset.v,1:2],adj=c(0,-0.5),cex=0.5)
# plot(dat[dat.subset.v,1:2],xlim=c(-57.81,-57.807), ylim=c(-64.933,-64.931))
# text(dat[dat.subset.v,1:2],adj=c(0,-0.5),cex=0.5)
# plot(dat[dat.subset.v,1:2],xlim=c(-57.801,-57.796), ylim=c(-64.937,-64.935))
# text(dat[dat.subset.v,1:2],adj=c(0,-0.5),cex=0.5)
##
interpl.and.repl(59,60,"6-9")
interpl.and.repl(104,106,"6-9")
dat[dat.subset.v[168:172],1:2] <- NA## remove
interpl.and.repl(206,211,"6-9")
interpl.and.repl(298,309,"6-9")
interpl.and.repl(323,324,"6-9")
interpl.and.repl(490,492,"6-9")
interpl.and.repl(643,647,"6-9")
interpl.and.repl(658,662,"6-9")
dat.subset.v <- which(dat$transectID=="6-9")
plot(dat[dat.subset.v,1:2])
text(dat[dat.subset.v,1:2],adj=c(0,-0.5),cex=0.5)

#### transect 7-1 has weird image locations as well
dat.subset.v <- which(dat$transectID=="7-1")
plot(dat[dat.subset.v,1:2])
text(dat[dat.subset.v,1:2],adj=c(0,-0.5),cex=0.5)
##
interpl.and.repl(39,39,"7-1")
dat <- dat[-which(is.na(dat[,1])),]
dat.subset.v <- which(dat$transectID=="7-1")
plot(dat[dat.subset.v,1:2])
text(dat[dat.subset.v,1:2],adj=c(0,-0.5),cex=0.5)

#### transect 8-1 has weird image locations as well
dat.subset.v <- which(dat$transectID=="8-1")
plot(dat[dat.subset.v,1:2])
text(dat[dat.subset.v,1:2],adj=c(0,-0.5),cex=0.5)
# spatial.dat <- data.frame(dat[dat.subset.v,1:2])
# coordinates(spatial.dat) <- c("Longitude","Latitude")
# proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
# polar.dat <- spTransform(spatial.dat, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# plot(polar.dat)
#text(polar.dat,adj=c(0,-0.5),cex=0.5)
# plot(polar.dat,xlim=c(-2378144+2000,-2374168),ylim=c(1607012+500, 1609632-2000))
# plot(polar.dat,xlim=c(-2378144+2000,-2374168),ylim=c(1607012+700, 1609632-1600))
# plot(polar.dat,xlim=c(-2378144+2000,-2374168-1300),ylim=c(1607012+800, 1609632-1600))
# text(polar.dat,adj=c(0,-0.5),cex=0.5)
##
dat[dat.subset.v[1],1:2] <- NA## remove
interpl.and.repl(165,165,"8-1")
interpl.and.repl(185,185,"8-1")
interpl.and.repl(317,317,"8-1")
interpl.and.repl(330,330,"8-1")
interpl.and.repl(334,340,"8-1")
interpl.and.repl(517,517,"8-1")
interpl.and.repl(594,594,"8-1")
dat.subset.v <- which(dat$transectID=="8-1")
plot(dat[dat.subset.v,1:2])
text(dat[dat.subset.v,1:2],adj=c(0,-0.5),cex=0.5)

#### transect 69
dat.subset.v <- which(dat$transectID=="69-1")
plot(dat[dat.subset.v,1:2])
text(dat[dat.subset.v,1:2],adj=c(0,-0.5),cex=0.5)
# spatial.dat <- data.frame(dat[dat.subset.v,1:2])
# coordinates(spatial.dat) <- c("Longitude","Latitude")
# proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
# polar.dat <- spTransform(spatial.dat, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# plot(polar.dat)
# text(polar.dat,adj=c(0,-0.5),cex=0.5)
#plot(polar.dat,xlim=c(-2485902+1000, -2482558-1500),ylim=c(2006812, 2010837-3300))
##
dat[dat.subset.v[30],1:2] <- NA## remove
#dat[dat.subset.v[180:181],1:2] <- NA## remove
interpl.and.repl(177,178,"69-1")
dat.subset.v <- which(dat$transectID=="69-1")
plot(dat[dat.subset.v,1:2])
text(dat[dat.subset.v,1:2],adj=c(0,-0.5),cex=0.5)

#### transect 81
dat.subset.v <- which(dat$transectID=="81-1")
plot(dat[dat.subset.v,1:2])
text(dat[dat.subset.v,1:2],adj=c(0,-0.5),cex=0.5)
# spatial.dat <- data.frame(dat[dat.subset.v,1:2])
# coordinates(spatial.dat) <- c("Longitude","Latitude")
# proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
# polar.dat <- spTransform(spatial.dat, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# plot(polar.dat)
# text(polar.dat,adj=c(0,-0.5),cex=0.5)
##
dat[dat.subset.v[118:122],1:2] <- NA## remove
dat.subset.v <- which(dat$transectID=="81-1")
plot(dat[dat.subset.v,1:2])
text(dat[dat.subset.v,1:2],adj=c(0,-0.5),cex=0.5)

dat <- dat[-which(is.na(dat[,1])),]

## Use distance to start for transects 9,11,12,81

####
## calculate transect length, define how many images to select and give images a random number
dat$image.select <- NA
total.t.length.v <- NA
dat$dist.from.start <- NA
dat$mean.dist.to.10.nearest.images <- NA
dat$dist.from.center <- NA
samp.v.list <- list()
center.idx.v <- NA
ends <- 1
for(i in 1:length(levels(dat$transectID))){
  print(levels(dat$transectID)[i])
  ## all data except for badly illuminated images
  dat.subset.v <- which(dat$transectID==levels(dat$transectID)[i])
  t.counts <- length(dat.subset.v)
  ## select first and last coordinate to calculate transect distance
  v.first <- dat.subset.v[1]
  v.last <- tail(dat.subset.v, n=1)
  # t.l <- distVincentyEllipsoid(dat[v.first,1:2], dat[v.last,1:2])
  # total.t.length.v[i] <- t.l
  # ## calculate distance of each image to the start of the transect along the transect line
  dat$dist.from.start[v.first] <- 0
  # if(levels(dat$transectID)[i]=="8-1"){
  #   dat$dist.from.start[v.last] <- 0
  #   for(k in (t.counts-1):1){
  #     dist.to.start <- distVincentyEllipsoid(dat[dat.subset.v[t.counts],1:2], dat[dat.subset.v[k],1:2])
  #     dat$dist.from.start[dat.subset.v[k]] <- dist.to.start
  #   }
#  }else 
    if(levels(dat$transectID)[i]%in%c("11-2","6-1","7-1","8-1","9-1","12-12","81-1")){
    for(k in 2:t.counts){
      dist.to.start <- distVincentyEllipsoid(dat[dat.subset.v[1],1:2], dat[dat.subset.v[k],1:2])
      dat$dist.from.start[dat.subset.v[k]] <- dist.to.start
    }
  }else{
    for(k in 2:t.counts){
      dist.to.previous.image <- distVincentyEllipsoid(dat[dat.subset.v[k-1],1:2], dat[dat.subset.v[k],1:2])
      dat$dist.from.start[dat.subset.v[k]] <- dat$dist.from.start[dat.subset.v[k-1]]+dist.to.previous.image
    }
  }
  total.t.length.v[i] <- dat$dist.from.start[dat.subset.v[t.counts]]
  print(dat$dist.from.start[dat.subset.v])

  ## choose a spatial random subset of images in each transect using quasiSamp on the distance from the start of the transect
  set.seed(48)
  samp <- quasiSamp(length(dat.subset.v),dimension=1,potential.sites=dat$dist.from.start[dat.subset.v])
  samp.v <- samp$ID[which(duplicated(samp)==FALSE)]
  print(paste0("first double selection for image # ", which(duplicated(samp)==TRUE)[1]))
  print(paste0("# of images in first batch: ", ceiling(total.t.length.v[i]/50)))
  #samp.v[(length(samp.v)+1):length(subset.v)] <- 9998
  dat$image.select[dat.subset.v[samp.v]] <-  1:length(samp.v)
  dat$image.select[dat.subset.v[-samp.v]] <-  9995
  # ## choose a random subset of images in each transect by giving them random numbers starting from 1
  #dat$image.select[dat.subset.v] <- sample(1:length(dat.subset.v),length(dat.subset.v))
  samp.v.list[[i]] <- samp.v
}
##
dat$image.select[is.na(dat$image.select)] <- 9999
dat$transectID.forfile <- substr(dat$transectID,1,2)
dat$transectID.forfile <- gsub("-","",dat$transectID.forfile)
dat$transectID.forfile <- sprintf("%02d",as.numeric(as.character(dat$transectID.forfile)))
dat$transectID.forfile <- as.factor(dat$transectID.forfile)

## total transect length across all the survey
#total.t.length <- sum(total.t.length.v)
## how many images from each transect, if we select on average 1 every 50m?
t.images <- ceiling(total.t.length.v/100)
## 500 images in total
#t.images <- ceiling((total.t.length.v/total.t.length)*500)

## crop transects 39 & 69 to become straight (replace numbers in image select with higher numbers)
t.images.straight <- t.images
crop.39 <- c(1:800,1780:2266)
crop.69 <- c(1:29,254:303)
nrow(dat)-length(c(crop.39,crop.69))

dat.subset.v <- which(dat$transectID=="39-1")
sel.raw <- order(dat$image.select[dat.subset.v])[1:t.images[3]]
sel <- sel.raw[which(sel.raw%!in%crop.39)]
plot(dat[dat.subset.v,1:2])
points(dat[dat.subset.v[sel],1:2],col="red")
dat$image.select[dat.subset.v[-sel]] <- dat$image.select[dat.subset.v[-sel]]+1000
t.images.straight[3] <- length(sel)

dat.subset.v <- which(dat$transectID=="69-1")
sel.raw <- order(dat$image.select[dat.subset.v])[1:t.images[5]]
sel <- sel.raw[which(sel.raw%!in%crop.69)]
plot(dat[dat.subset.v,1:2])
points(dat[dat.subset.v[sel],1:2],col="red")
dat$image.select[dat.subset.v[-sel]] <- dat$image.select[dat.subset.v[-sel]]+1000
t.images.straight[5] <- length(sel)

## correct image.select numbers that are larger than 4 digits
correct <- dat$image.select[dat$image.select>10000]-1000
dat$image.select[dat$image.select>10000] <- correct

#save(dat,total.t.length.v,t.images.straight,samp.v.list,crop.39,crop.69, file="C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/PS118_dat.Rdata")
#load(file="C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/PS118_dat.Rdata")


## USING STRAIGHT TRANSECT SUBSET ONLY:
t.images <- t.images.straight

par(mfrow=c(1,1))
# for(i in 1:length(levels(dat$transectID))){
#   plot(dat$mean.dist.to.10.nearest.images[dat$transectID==levels(dat$transectID)[i]])
# }

spatial.dat <- data.frame(dat[,1:2])
coordinates(spatial.dat) <- c("Longitude","Latitude")
proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
polar.dat <- spTransform(spatial.dat, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

#par(mfrow=c(4,4), mar=c(5,4,2,1))
for(i in 1:length(levels(dat$transectID))){
  subset.v <- which(dat$transectID==levels(dat$transectID)[i])
  ## select the correct number of images for each transect
  sel <- order(dat$image.select[subset.v])[1:t.images[i]]
  ## create vector of selected images
  if(i==1){dat.sel <<-  subset.v[sel]
  }else dat.sel <-  c(dat.sel,subset.v[sel])
  plot(polar.dat[subset.v], main=levels(dat$transectID)[i])
  points(polar.dat[subset.v[sel]], col="red", pch=15)
  #points(polar.dat[subset.v[center.idx.v[i]],1:2], col="green", pch=16)
  scalebar(200, type="bar")
  # plot(dat[subset.v,1:2])
  # points(dat[subset.v[sel],1:2], col="red", pch=15)
  # points(dat[subset.v[center.idx.v[i]],1:2], col="green", pch=16)
}

dat$Filename[dat.sel]
## filenames
selected.filenames <- paste0(dat$Filename[dat.sel],".jpg")
## create new filenames and identify folders where to find them (might be different for each survey)
selected.filenames.renamed <- paste0("PS118_",dat$transectID.forfile[dat.sel],"_",sprintf("%04d",dat$image.select[dat.sel]),"__",selected.filenames)
selected.filenames.folders <- paste0("PS118_",dat$transectID.forfile[dat.sel])

## copy files into Annotation folder
img.path.origin <- paste0(image.dir,selected.filenames.folders,"/",selected.filenames)
img.path.destin <- "C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/Stills/Annotation_images/PS118/"
file.copy(img.path.origin,img.path.destin)

## write list of filenames into cropped Annotation folder for upload
img.path.destin_crop <- "C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/Stills/Annotation_images_cropped/PS118/"
write.table(selected.filenames.renamed, paste0(img.path.destin_crop,"PS118_filenames.txt"), eol=",", col.names=FALSE, row.names=FALSE)

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
#magick mogrify -path Annotation_images_cropped\PS118\ -gravity Center -crop 80% +repage Annotation_images\PS118\*jpg
## @@@ !!! @@@ 































################################################################################################################
##### IMAGES SHALLOWER THAN 100m DEPTH #########################################################################
################################################################################################################
## identify images that are shallower than 100m
bad_depth <- paste0(PS118.dat$Filename[which(PS118.dat$Depth.water..m.>100)],".jpg")
## identify images taken using the hotkey rather than the timer
evil_hotkey <- paste0(PS118.dat$Filename[which(PS118.dat$Type=="HOTKEY")],".jpg")

## remove all of these bad images
filenames <- paste0(PS118.dat$Filename,".jpg")
idx <- which(filenames%in%bad_images|filenames%in%bad_depth|filenames%in%evil_hotkey)
dat <- PS118.dat[-idx,]

## only transect 11 truly has images shallower than 100m, the rest are errors in the data
dat <- dat[dat$transectID=="11-2",]
dat$transectID <- factor(dat$transectID)

## plot all transects
spatial.dat <- data.frame(dat[,1:2])
coordinates(spatial.dat) <- c("Longitude","Latitude")
proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
polar.dat <- spTransform(spatial.dat, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
par(mfrow=c(1,1), mar=c(3,2,2,1),oma=c(0,0,0,2))
for(i in 1:length(levels(dat$transectID))){
  subset.v <- which(dat$transectID==levels(dat$transectID)[i])
  xlim <- extent(polar.dat[subset.v])[1:2]
  xlim[1] <- xlim[1]-2000
  xlim[2] <- xlim[2]+2000
  ylim <- extent(polar.dat[subset.v])[3:4]
  ylim[1] <- ylim[1]-2000
  ylim[2] <- ylim[2]+2000
  plot(r2, xlim=xlim, ylim=ylim,main=levels(dat$transectID)[i])
  points(polar.dat[subset.v])
  #text(polar.dat[subset.v])
  scalebar(500, type="bar")
}

####
## calculate transect length, define how many images to select and give images a random number
dat$image.select <- NA
total.t.length.v <- NA
dat$dist.from.start <- NA
dat$mean.dist.to.10.nearest.images <- NA
dat$dist.from.center <- NA
samp.v.list <- list()
center.idx.v <- NA
ends <- 1
for(i in 1:length(levels(dat$transectID))){
  print(levels(dat$transectID)[i])
  ## all data except for badly illuminated images
  dat.subset.v <- which(dat$transectID==levels(dat$transectID)[i])
  t.counts <- length(dat.subset.v)
  ## select first and last coordinate to calculate transect distance
  v.first <- dat.subset.v[1]
  v.last <- tail(dat.subset.v, n=1)
  # t.l <- distVincentyEllipsoid(dat[v.first,1:2], dat[v.last,1:2])
  # total.t.length.v[i] <- t.l
  # ## calculate distance of each image to the start of the transect NOT along the transect line
  dat$dist.from.start[v.first] <- 0
  for(k in 2:t.counts){
    dist.to.start <- distVincentyEllipsoid(dat[dat.subset.v[1],1:2], dat[dat.subset.v[k],1:2])
    dat$dist.from.start[dat.subset.v[k]] <- dist.to.start
  }
  total.t.length.v[i] <- dat$dist.from.start[dat.subset.v[t.counts]]
  print(dat$dist.from.start[dat.subset.v])
  ## choose a spatial random subset of images in each transect using quasiSamp on the distance from the start of the transect
  set.seed(42)
  samp <- quasiSamp(length(dat.subset.v),dimension=1,potential.sites=dat$dist.from.start[dat.subset.v])
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
dat$transectID.forfile <- substr(dat$transectID,1,2)
dat$transectID.forfile <- gsub("-","",dat$transectID.forfile)
dat$transectID.forfile <- sprintf("%02d",as.numeric(as.character(dat$transectID.forfile)))
dat$transectID.forfile <- as.factor(dat$transectID.forfile)

#save(dat,total.t.length.v,samp.v.list, file="C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/PS118_dat_shallow.Rdata")
#load(file="C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/PS118_dat_shallow.Rdata")

## total transect length across all the survey
#total.t.length <- sum(total.t.length.v)
## how many images from each transect, if we select on average 1 every 50m?
t.images <- ceiling(total.t.length.v/100)
## 500 images in total
#t.images <- ceiling((total.t.length.v/total.t.length)*500)
##

par(mfrow=c(1,1))
# for(i in 1:length(levels(dat$transectID))){
#   plot(dat$mean.dist.to.10.nearest.images[dat$transectID==levels(dat$transectID)[i]])
# }

spatial.dat <- data.frame(dat[,1:2])
coordinates(spatial.dat) <- c("Longitude","Latitude")
proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
polar.dat <- spTransform(spatial.dat, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

#par(mfrow=c(4,4), mar=c(5,4,2,1))
for(i in 1:length(levels(dat$transectID))){
  subset.v <- which(dat$transectID==levels(dat$transectID)[i])
  ## select the correct number of images for each transect
  sel <- order(dat$image.select[subset.v])[1:t.images[i]]
  ## create vector of selected images
  if(i==1){dat.sel <<-  subset.v[sel]
  }else dat.sel <-  c(dat.sel,subset.v[sel])
  plot(polar.dat[subset.v], main=levels(dat$transectID)[i])
  points(polar.dat[subset.v[sel]], col="red", pch=15)
  #points(polar.dat[subset.v[center.idx.v[i]],1:2], col="green", pch=16)
  scalebar(100, type="bar")
  # plot(dat[subset.v,1:2])
  # points(dat[subset.v[sel],1:2], col="red", pch=15)
  # points(dat[subset.v[center.idx.v[i]],1:2], col="green", pch=16)
}

dat$Filename[dat.sel]
## filenames
selected.filenames <- paste0(dat$Filename[dat.sel],".jpg")
## create new filenames and identify folders where to find them (might be different for each survey)
selected.filenames.renamed <- paste0("PS118_",dat$transectID.forfile[dat.sel],"_",sprintf("%04d",dat$image.select[dat.sel]),"__",selected.filenames)
selected.filenames.folders <- paste0("PS118_",dat$transectID.forfile[dat.sel])

## copy files into Annotation folder
img.path.origin <- paste0(image.dir,selected.filenames.folders,"/",selected.filenames)
img.path.destin <- "C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/Stills/Annotation_images/PS118_shallow/"
file.copy(img.path.origin,img.path.destin)

## write list of filenames into cropped Annotation folder for upload
img.path.destin_crop <- "C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/Stills/Annotation_images_cropped/PS118_shallow/"
write.table(selected.filenames.renamed, paste0(img.path.destin_crop,"PS118_filenames.txt"), eol=",", col.names=FALSE, row.names=FALSE)

filenames.in.folder_original <- paste0(img.path.destin,selected.filenames)
filenames.in.folder_changed <- paste0(img.path.destin,selected.filenames.renamed)
file.rename(filenames.in.folder_original,filenames.in.folder_changed)
