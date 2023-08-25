'%!in%' <- function(x,y)!('%in%'(x,y))
library(geosphere)
library(ggplot2)
library(MBHdesign)
library(raster)
library(lubridate)
library(openxlsx)
library(measurements)


# #### load depth
library(raadtools)
my_data_dir <- "C:/Users/jjansen/Desktop/science/data_environmental/raw/accessed_through_R"
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

##################################
######### CRS (WAP data)
##################################

## INTERPOALTE BETWEEN START AND END POINTS OF TRANSECTS TO GET GPS COORDINATES FOR EACH IMAGE
## USE FIRST AND LAST IMAGES FROM CORRECTED FOLDERS BUT CHECK BECAUSE SOMETIMES A NUMBER OF IMAGES AT THE START/END ARE MISSING
## 1069: 1,2 missing
## 1297: use images between 4:111 (check original folder)

bio.path <- "R:/IMAS/Antarctic_Seafloor/Clean_Data_For_Permanent_Storage/CRS/CRS_1_raw_images_and_metadata/"
dat.WAP <- read.xlsx(paste0(bio.path,"metadata/Site data.xlsx"),sheet=2)
dat.stations <- read.xlsx(paste0(bio.path,"metadata/Site data.xlsx"),sheet=1)
names(dat.WAP)[3:7] <- c("transectID","Latitude_Start","Longitude_Start","Latitude_End","Longitude_End")

## correct some random values in dms to decimal degrees
for(i in 4:7){
  y <- dat.WAP[30:34,i]
  dat.WAP[30:34,i] <- conv_unit(y, from='deg_dec_min',to='dec_deg')
}
## fill missing stations etc
dat.WAP[2:8,1] <- dat.WAP[1,1]
dat.WAP[10:14,1] <- dat.WAP[9,1]
dat.WAP[16:17,1] <- dat.WAP[15,1]
dat.WAP[19:21,1] <- dat.WAP[18,1]
dat.WAP[23:25,1] <- dat.WAP[22,1]
dat.WAP[27:29,1] <- dat.WAP[26,1]
dat.WAP[32:34,1] <- dat.WAP[31,1]
dat.WAP[2,2] <- dat.WAP[1,2]
dat.WAP[4,2] <- dat.WAP[3,2]
dat.WAP[6,2] <- dat.WAP[5,2]
dat.WAP[8,2] <- dat.WAP[7,2]
dat.WAP[10,2] <- dat.WAP[9,2]
dat.WAP[12,2] <- dat.WAP[11,2]
dat.WAP[14,2] <- dat.WAP[13,2]
dat.WAP[17,2] <- dat.WAP[16,2]
dat.WAP[19,2] <- dat.WAP[18,2]
dat.WAP[21,2] <- dat.WAP[20,2]
dat.WAP[23,2] <- dat.WAP[22,2]
dat.WAP[25,2] <- dat.WAP[24,2]
dat.WAP[27,2] <- dat.WAP[26,2]
dat.WAP[29,2] <- dat.WAP[28,2]
for(i in 1:3){dat.WAP[,i] <- as.factor(dat.WAP[,i])}
for(i in 4:7){dat.WAP[,i] <- as.numeric(dat.WAP[,i])}

## create dataframe with images, transectID, individual coordinates and whether images are good or bad
WAP_metadata <- data.frame(matrix(NA, ncol=6,nrow=1))
names(WAP_metadata) <- c("Filename","transectID","lon","lat","imageQuality", "depth")
#### load raw image names
dir.folders <- list.files(paste0(bio.path,"images_original/"))
# dir.files.raw <- list.files(paste0(bio.path,"images_original/"), recursive=T)
# dir.files.cor <- list.files(paste0(bio.path,"images_colourcorrected/"), recursive=T)
## run first transect (1069) separately as images 1+2 are missing
loop.files.raw <- list.files(paste0(bio.path,"images_original/",dir.folders[1],"/"))
loop.files.cor <- list.files(paste0(bio.path,"images_colourcorrected/",dir.folders[1],"/"), pattern="S")
total.n <- 49
transect.sel <- which(dat.WAP$transectID==substr(dir.folders[1],5,8))
lon.seq <- seq(dat.WAP[transect.sel,5],dat.WAP[transect.sel,7], length.out=49)
lat.seq <- seq(dat.WAP[transect.sel,4],dat.WAP[transect.sel,6], length.out=49)
WAP_metadata[1:49,1] <- loop.files.raw
WAP_metadata[1:49,2] <- substr(dir.folders[1],5,8)
WAP_metadata[1:49,3] <- lon.seq
WAP_metadata[1:49,4] <- lat.seq
WAP_metadata[1:49,6] <- dat.WAP$`Mean.depth.(m)`[transect.sel]
bad <- which(loop.files.raw%!in%loop.files.cor)
WAP_metadata[c(1:49)[bad],5] <- "bad"
## now for all other transects
for(i in 2:length(dir.folders)){
  t.ID <- substr(dir.folders[i],5,8)
  print(t.ID)
  ## find image names (it's a bit annoying that it tries to include the folder "bad_quality". Can't simply use .jpg cause some images are tif)
  loop.files.raw <- list.files(paste0(bio.path,"images_original/",dir.folders[i],"/"))
  loop.files.cor <- list.files(paste0(bio.path,"images_colourcorrected/",dir.folders[i],"/"), pattern="S")
  ## now look how many images in the raw folder between transect start and end, where to find them and where to put them in the new dataframe
  raw.start <- which(loop.files.raw==loop.files.cor[1])
  raw.end <- which(loop.files.raw==tail(loop.files.cor,1))
  folder.sel <- raw.start:raw.end
  #for transect 1297: use images between 4:111 (check original folder)
  if(t.ID=="1297"){folder.sel <- 4:111}
  ##
  n <- length(folder.sel)
  dat.sel <- (total.n+1):(total.n+n)
  ## find transect and interpolate coordinates
  transect.sel <- which(dat.WAP$transectID==t.ID)
  lon.seq <- seq(dat.WAP[transect.sel,5],dat.WAP[transect.sel,7], length.out=n)
  lat.seq <- seq(dat.WAP[transect.sel,4],dat.WAP[transect.sel,6], length.out=n)
  ## assign to dataframe
  WAP_metadata[dat.sel,1] <- loop.files.raw[folder.sel]
  WAP_metadata[dat.sel,2] <- t.ID
  WAP_metadata[dat.sel,3] <- lon.seq
  WAP_metadata[dat.sel,4] <- lat.seq
  WAP_metadata[dat.sel,6] <- dat.WAP$`Mean.depth.(m)`[transect.sel]
  bad <- which(loop.files.raw[folder.sel]%!in%loop.files.cor)
  WAP_metadata[dat.sel[bad],5] <- "bad"
  ## update the running total
  total.n <- total.n+n
}
WAP_metadata$transectID <- as.factor(WAP_metadata$transectID)
WAP_metadata$imageQuality[is.na(WAP_metadata$imageQuality)] <- "good"
WAP_metadata$imageQuality <- as.factor(WAP_metadata$imageQuality)


##### plots
# dat <- WAP_metadata
# 
# spatial.dat <- data.frame(dat[,3:4])
# coordinates(spatial.dat) <- c("lon","lat")
# proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
# polar.dat <- spTransform(spatial.dat, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# 
# spatial.dat <- data.frame(dat.WAP[,c(5,4)])
# coordinates(spatial.dat) <- c("Longitude_Start","Latitude_Start")
# proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
# polar.dat.start <- spTransform(spatial.dat, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# 
# spatial.dat <- data.frame(dat.WAP[,c(7,6)])
# coordinates(spatial.dat) <- c("Longitude_End","Latitude_End")
# proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
# polar.dat.end <- spTransform(spatial.dat, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# 
# spatial.dat <- data.frame(dat.stations[,c(3,2)])
# coordinates(spatial.dat) <- c("Longitude","Latitude")
# proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
# polar.dat.station <- spTransform(spatial.dat, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# 
# col.breaks <- seq(-300,-700,length.out=101)
# cols <- terrain.colors(99)
# 
# par(mfrow=c(3,4), mar=c(5,4,2,1))
# #par(mfrow=c(1,1), mar=c(5,4,2,1))
# for(i in 1:length(levels(dat$transectID))){
#   subset.v <- which(dat$transectID==levels(dat$transectID)[i])
#   ## create vector of selected images
#   if(i==1){dat.sel <<-  subset.v
#   }else dat.sel <-  c(dat.sel,subset.v)
#   xlim <- NA
#   ylim <- NA
#   xlim[1] <- extent(polar.dat[subset.v])[1]-1000
#   xlim[2] <- extent(polar.dat[subset.v])[2]+1000
#   ylim[1] <- extent(polar.dat[subset.v])[3]-1000
#   ylim[2] <- extent(polar.dat[subset.v])[4]+1000
#   plot(r2, xlim=xlim,ylim=ylim, main=paste0("transect ",levels(dat$transectID)[i]),col=cols, breaks=col.breaks)
#   plot(polar.dat[subset.v],add=TRUE)
#   points(polar.dat.start,col="red")
#   points(polar.dat.end, col="blue")
#   points(polar.dat.station, col="maroon",cex=2,pch=15)
#   #plot(coast.proj, add=TRUE)
#   scalebar(100, type="bar")
# }
# par(mfrow=c(1,1))
# plot(r2, xlim=c(-2700000,-2200000),ylim=c(800000,1550000),main="Overview")#, col=cols, breaks=col.breaks)
# plot(polar.dat,add=TRUE)
# points(polar.dat.start,col="red")
# points(polar.dat.end, col="blue")
# points(polar.dat.station, col="maroon",pch=15)
# #plot(coast.proj, add=TRUE)
# scalebar(200000, type="bar", label=c(0,100,200),below="km")
# 
# ## weddell site only
# plot(r2, xlim=c(-2400000,-2370000),ylim=c(1490000,1510000), main="Lockyer")#, col=cols, breaks=col.breaks)
# plot(polar.dat,add=TRUE,cex=0.5)
# points(polar.dat.start,col="red")
# points(polar.dat.end, col="blue")
# points(polar.dat.station, col="maroon",pch=15)
# #plot(coast.proj, add=TRUE)
# scalebar(5000, type="bar", below="m")
# 
# # ## fjords (3 northern ones)
# # plot(r2, xlim=c(-2510000,-2430000),ylim=c(1230000,1390000))#, col=cols, breaks=col.breaks)
# # plot(polar.dat,add=TRUE)
# # points(polar.dat.start,col="red")
# # points(polar.dat.end, col="blue")
# # points(polar.dat.station, col="maroon",pch=15)
# # #plot(coast.proj, add=TRUE)
# # scalebar(10000, type="bar")
# 
# ## fjords top N
# plot(r2, xlim=c(-2510000,-2480000),ylim=c(1360000,1380000), main="Hughes")#, col=cols, breaks=col.breaks)
# plot(polar.dat,add=TRUE,cex=0.5)
# points(polar.dat.start,col="red")
# points(polar.dat.end, col="blue")
# #points(polar.dat.station, col="maroon",pch=15)
# text(polar.dat.end,labels=dat.WAP$transectID, adj=c(1,-0.5))
# #plot(coast.proj, add=TRUE)
# scalebar(5000, type="bar",below="m")
# 
# # ## fjords top middle + bottom middle
# # plot(r2, xlim=c(-2490000,-2440000),ylim=c(1230000,1285000))#, col=cols, breaks=col.breaks)
# # plot(polar.dat,add=TRUE)
# # points(polar.dat.start,col="red")
# # points(polar.dat.end, col="blue")
# # points(polar.dat.station, col="maroon",pch=15)
# # #plot(coast.proj, add=TRUE)
# # scalebar(10000, type="bar")
# 
# ## fjords top middle
# plot(r2, xlim=c(-2485000,-2455000),ylim=c(1263000,1283000),main="Anvord Bay")#, col=cols, breaks=col.breaks)
# plot(polar.dat,add=TRUE, cex=0.5)
# points(polar.dat.start,col="red")
# points(polar.dat.end, col="blue")
# points(polar.dat.station, col="maroon",pch=15)
# #plot(coast.proj, add=TRUE)
# scalebar(5000, type="bar", below="m")
# 
# ## fjords bottom middle
# plot(r2, xlim=c(-2473000,-2443000),ylim=c(1230000,1250000),main="Flandres Bay")#, col=cols, breaks=col.breaks)
# plot(polar.dat,add=TRUE,cex=0.5)
# points(polar.dat.start,col="red")
# points(polar.dat.end, col="blue")
# points(polar.dat.station, col="maroon",pch=15)
# #plot(coast.proj, add=TRUE)
# scalebar(5000, type="bar", below="m")
# 
# ## fjords south
# plot(r2, xlim=c(-2420000,-2390000),ylim=c(1125000,1145000),main="Barilari Bay")#, col=cols, breaks=col.breaks)
# plot(polar.dat,add=TRUE,cex=0.5)
# points(polar.dat.start,col="red")
# points(polar.dat.end, col="blue")
# points(polar.dat.station, col="maroon",pch=15)
# #plot(coast.proj, add=TRUE)
# scalebar(5000, type="bar", below="m")
# 
# # ## shelf sites
# # plot(r2, xlim=c(-2550000,-2350000),ylim=c(850000,1170000))#, col=cols, breaks=col.breaks)
# # plot(polar.dat,add=TRUE,cex=0.5)
# # points(polar.dat.start,col="red")
# # points(polar.dat.end, col="blue")
# # points(polar.dat.station, col="maroon",pch=15)
# # #plot(coast.proj, add=TRUE)
# # scalebar(5000, type="bar", below="m")
# 
# ## shelf sites N
# plot(r2, xlim=c(-2540000,-2510000),ylim=c(1150000,1170000),main="Station B")#, col=cols, breaks=col.breaks)
# plot(polar.dat,add=TRUE,cex=0.5)
# points(polar.dat.start,col="red")
# points(polar.dat.end, col="blue")
# points(polar.dat.station, col="maroon",pch=15)
# #plot(coast.proj, add=TRUE)
# scalebar(5000, type="bar", below="m")
# 
# ## shelf sites Middle
# plot(r2, xlim=c(-2460000,-2430000),ylim=c(1010000,1030000),main="Station E")#, col=cols, breaks=col.breaks)
# plot(polar.dat,add=TRUE,cex=0.5)
# points(polar.dat.start,col="red")
# points(polar.dat.end, col="blue")
# points(polar.dat.station, col="maroon",pch=15)
# #plot(coast.proj, add=TRUE)
# scalebar(5000, type="bar", below="m")
# 
# ## shelf sites S
# plot(r2, xlim=c(-2390000,-2360000),ylim=c(870000,890000), main="Station F")#, col=cols, breaks=col.breaks)
# plot(polar.dat,add=TRUE,cex=0.5)
# points(polar.dat.start,col="red")
# points(polar.dat.end, col="blue")
# points(polar.dat.station, col="maroon",pch=15)
# #plot(coast.proj, add=TRUE)
# scalebar(5000, type="bar", below="m")


#####################################################################
##### 2. SUBSET IMAGES FROM IMAGE LOCATIONS (& STORE FILENAMES) #####
#####################################################################
## dataframe containing: coordinates, filename, transectID

## now remove bad quality images
dat <- WAP_metadata[which(WAP_metadata$imageQuality=="good"),c(3,4,1,2,5,6)]
dat$transectID <- factor(dat$transectID)

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
  for(k in 2:t.counts){
    dist.to.previous.image <- distVincentyEllipsoid(dat[dat.subset.v[k-1],1:2], dat[dat.subset.v[k],1:2])
    dat$dist.from.start[dat.subset.v[k]] <- dat$dist.from.start[dat.subset.v[k-1]]+dist.to.previous.image
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

#### CRS survey information 
## Anvord, Flanders and Barilla in 2010 during NBP10-01 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5214585/)
library(lubridate)
dat$Date <- as.Date(NA)
unique(dat$transectID)
dat$Date[dat$transectID=="1289"] <- ymd("2010:01:20")
dat$Date[dat$transectID=="1290"] <- ymd("2010:01:20")
dat$Date[dat$transectID=="1283"] <- ymd("2010:01:18")
dat$Date[dat$transectID=="1284"] <- ymd("2010:01:19")
dat$Date[dat$transectID=="1337"] <- ymd("2010:02:20")
dat$Date[dat$transectID=="1338"] <- ymd("2010:02:20")
dat$Date[dat$transectID=="1285"] <- ymd("2010:01:19")
dat$Date[dat$transectID=="1286"] <- ymd("2010:01:19")
dat$Date[dat$transectID=="1281"] <- ymd("2010:01:18")
dat$Date[dat$transectID=="1282"] <- ymd("2010:01:18")
dat$Date[dat$transectID=="1279"] <- ymd("2010:01:17")
dat$Date[dat$transectID=="1280"] <- ymd("2010:01:17")
dat$Date[dat$transectID=="1276"] <- ymd("2010:01:14")
dat$Date[dat$transectID=="1278"] <- ymd("2010:01:17")
dat$Date[dat$transectID=="1300"] <- ymd("2010:01:27")
dat$Date[dat$transectID=="1295"] <- ymd("2010:01:24")
dat$Date[dat$transectID=="1297"] <- ymd("2010:01:26")
## LMG0902: https://www.marine-geo.org/tools/entry/LMG0902
dat$Date[dat$transectID=="1207"] <- ymd("2009:02:27")
dat$Date[dat$transectID=="1208"] <- ymd("2009:02:27")
dat$Date[dat$transectID=="1217"] <- ymd("2009:02:28")
dat$Date[dat$transectID=="1219"] <- ymd("2009:02:28")
dat$Date[dat$transectID=="1255"] <- ymd("2009:03:05")
dat$Date[dat$transectID=="1267"] <- ymd("2009:03:06")

## added later (August 2023)
dat$Date[dat$transectID=="1069"] <- ymd("2008:07:19")
dat$Date[dat$transectID=="1072"] <- ymd("2008:07:19")
dat$Date[dat$transectID=="1091"] <- ymd("2008:07:21")
dat$Date[dat$transectID=="1103"] <- ymd("2008:07:22")
dat$Date[dat$transectID=="1130"] <- ymd("2008:07:26")
dat$Date[dat$transectID=="1132"] <- ymd("2008:07:26")
dat$Date[dat$transectID=="1315"] <- ymd("2010:02:08")
dat$Date[dat$transectID=="1320"] <- ymd("2010:02:15")
dat$Date[dat$transectID=="1323"] <- ymd("2010:02:15")
dat$Date[dat$transectID=="1324"] <- ymd("2010:02:15")
dat$Date[dat$transectID=="1325"] <- ymd("2010:02:16")


## SAVE OUTPUT FOR FUTURE REFENCE (i.e. start here to add more images to the analysis)
# save(dat,total.t.length.v, file="C:/Users/jjansen/Desktop/science/SouthernOceanBiodiversityMapping/ARC_Data/prep_image/ReadIn_Circumpolar_DownwardImage_Data_And_Subset_Surveys/CRS_dat.Rdata")
#load(file="C:/Users/jjansen/Desktop/science/SouthernOceanBiodiversityMapping/ARC_Data/prep_image/ReadIn_Circumpolar_DownwardImage_Data_And_Subset_Surveys/CRS_dat.Rdata")

## how many images from each transect, if we select on average 1 every 50m?
t.images <- ceiling(total.t.length.v/100)
#t.images <- ceiling((total.t.length.v/total.t.length)*500)
## 

spatial.dat <- data.frame(dat[,1:2])
coordinates(spatial.dat) <- c("lon","lat")
proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
polar.dat <- spTransform(spatial.dat, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

col.breaks <- seq(-100,-1400,length.out=101)
cols <- terrain.colors(99)

par(mfrow=c(1,1), mar=c(5,4,2,1))
for(i in 1:length(levels(dat$transectID))){
  subset.v <- which(dat$transectID==levels(dat$transectID)[i])
  ## select the correct number of images for each transect
  sel <- order(dat$image.select[subset.v])[1:t.images[i]]
  if(t.images[i]>length(dat$image.select[subset.v])){
    sel <- order(dat$image.select[subset.v])
  }
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
selected.filenames.renamed <- paste0("CRS_",dat$transectID[dat.sel],"_",sprintf("%04d",dat$image.select[dat.sel]),"__",selected.filenames)
## list folders, match with transectIDs to look up each image from the correct folder
folder.names_full <- list.dirs(paste0(bio.path,"images_original/"),recursive=FALSE,full.names=FALSE)
selected.filenames.folders <- folder.names_full[match(dat$transectID[dat.sel],levels(dat$transectID))]

# ## rename original files because they are not unique between transects
# original.files.before <- paste0(bio.path,"images_original/",dat$Filename,"/",selected.filenames) #bio.path,"images_original/"
# original.files.renamed <- paste0("CRS_",dat$transectID[dat.sel],"_",sprintf("%04d",dat$image.select[dat.sel]),"__",selected.filenames)
# 
# file.rename(filenames.in.folder_original,filenames.in.folder_changed)


## copy files into Annotation folder and rename immediately because filenames are not unique between transects
path_D <- "D:/ARC_DP_data/adjusted_CRS/"
img.path.origin <- paste0(path_D,selected.filenames.folders,"/",selected.filenames) #bio.path,"images_original/"
img.path.destin <- "C:/Users/jjansen/Desktop/science/data_biological/Stills/Annotation_images/CRS/"
#file.copy(img.path.origin,img.path.destin)
for(i in 1:length(img.path.origin)){
  file.copy(img.path.origin[i],img.path.destin)
  org <- paste0(img.path.destin,selected.filenames[i])
  ren <- paste0(img.path.destin,selected.filenames.renamed[i])
  file.rename(org,ren)
}
# ## rename files
# filenames.in.folder_original <- paste0(img.path.destin,selected.filenames)
# filenames.in.folder_changed <- paste0(img.path.destin,selected.filenames.renamed)
# file.rename(filenames.in.folder_original,filenames.in.folder_changed)

## write list of filenames into cropped Annotation folder for upload
img.path.destin_crop <- "C:/Users/jjansen/Desktop/science/data_biological/Stills/Annotation_images_cropped/CRS/"
write.table(selected.filenames.renamed, paste0(img.path.destin_crop,"filenames.txt"), eol=",", col.names=FALSE, row.names=FALSE)



########################################################
##### 3. CROP SELECTED IMAGES AND STORE SEPARATELY #####
########################################################

## NO CROPPING












