########
## This script is to fix the image selection for the first runs in PS96 and PS81
## output will be stored as PS81_dat_FINAL.Rdata and PS96_dat_FINAL.Rdata

library(raster)
'%!in%' <- function(x,y)!('%in%'(x,y))

load(file="C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/PS96_dat_1strun.Rdata")
dat[1:699,] <- first_three_transects[,-13]
PS96.dat <- dat
PS96_total.t.length.v <- total.t.length.v

load(file="C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/PS96_dat_2ndrun.Rdata")
PS96b.dat <- dat
PS96b_total.t.length.v <- total.t.length.v
PS96.dat$Height..m. <- PS96b.dat$Height..m.

load(file="C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/PS81_dat_1strun.Rdata")
PS81a.dat <- dat
PS81a_total.t.length.v <- total.t.length.v
PS81a_samp.v.list <- samp.v.list
PS81a.dat$image.select[PS81a.dat$image.select==9998] <- 999

load(file="C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/PS81_dat_2ndrun.Rdata")
PS81b.dat <- dat
PS81b_total.t.length.v <- total.t.length.v
PS81b_samp.v.list <- samp.v.list
PS81b.dat$image.select[PS81b.dat$image.select==9998] <- 999

PS81.dat <- PS81a.dat[1:3523,]
PS81.dat[3524:14177,] <- PS81b.dat[3524:14177,c(1:9,11:16)]
PS81.dat$Dist.subs..m. <- PS81b.dat$Dist.subs..m.

PS81_total.t.length.v <- c(PS81a_total.t.length.v[1:7],PS81b_total.t.length.v[8:31])

rm(dat,total.t.length.v,samp.v.list,first_three_transects,PS81a.dat)

##### ######################################## #####
##### ###### @@@@@@@@@@ PS81 @@@@@@@@@@ ###### #####
##### ######################################## #####
## check existing images and the ones to replace
## choosing legacy sites doesn't work properly, so I went manual

dat <- PS81.dat
total.t.length.v <- PS81_total.t.length.v
## store original image select separately
org.img.select <- PS81.dat$image.select
## wipe image select column
dat$image.select <- NA
## fill with PS81b data for the b-transects
dat$image.select[3524:14177] <- PS81b.dat$image.select[3524:14177]


## first plot all transects
t.images <- ceiling(total.t.length.v/100)
spatial.dat <- data.frame(dat[,1:2])
coordinates(spatial.dat) <- c("Longitude","Latitude")
proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
polar.dat <- spTransform(spatial.dat, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
#par(mfrow=c(4,4), mar=c(5,4,2,1))
for(i in 1:length(levels(dat$t.ID_temp))){
  subset.v <- which(dat$t.ID_temp==levels(dat$t.ID_temp)[i])
  subset.v.good <- which(dat$Dist.subs..m.[subset.v]>=1&dat$Dist.subs..m.[subset.v]<=2)
  ## select the correct number of images for each transect
  sel <- order(org.img.select[subset.v])[1:t.images[i]]
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
store.dat.sel <- dat.sel

## check which images from the first few transects are further away than 2m range
not.in.range <- which(dat$Dist.subs..m.[store.dat.sel]>2)
dat[store.dat.sel,][not.in.range,]
rmv.idx.full <- as.numeric(rownames(dat[store.dat.sel,][not.in.range,]))

plot(dat$Dist.subs..m.[store.dat.sel])
text(dat$Dist.subs..m.[store.dat.sel],labels=paste0(dat$transectID[store.dat.sel],"_00",org.img.select[store.dat.sel]))

dat$Filename[rmv.idx.full]

##################
#### TRANSECT 116:
rmv.idx <- rmv.idx.full[1:8]
i=1
subset.v <- which(dat$t.ID_temp==levels(dat$t.ID_temp)[i])
subset.v.good <- which(dat$Dist.subs..m.[subset.v]>=1&dat$Dist.subs..m.[subset.v]<=2)
sel <- order(org.img.select[subset.v])[1:t.images[i]]
keep.idx <- which(subset.v[sel]%!in%rmv.idx)

## plot transect 116 and determine which images should be scored instead of the badly illuminated ones
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[rmv.idx], col="green", pch=17)
#text(polar.dat[rmv.idx,1:2], col="green", adj=-0.5, labels=paste0(dat$transectID[rmv.idx],"_00",org.img.select[rmv.idx]))#, pch=15)
scalebar(500, type="bar")

#### order of all next images in line according to PS81b image selection:
already.scored <- order(org.img.select[subset.v])[1:t.images[i]]
bad.dist <- which(PS81b.dat$Dist.subs..m.[subset.v]<1|PS81b.dat$Dist.subs..m.[subset.v]>2)
scored.and.bad <- bad.dist[which(bad.dist%in%already.scored)]
scored.or.bad <- union(already.scored, bad.dist)
sel.add <- subset.v[-scored.or.bad]
## sort the images into a spatially balance sequence
img.seq <- order(PS81b.dat$image.select[sel.add])
## check where the next images are positioned
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[1:50]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[1:50]]],adj=2)
## select which images are too close to already scored images
dontUse <- c(1,3,4,5,6,7,9,11,12,14,15,16,18,23,24,25,33,35,36,38,42,47,56,79)
idx.range <- 1:80
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[idx.range[dontUse]]]], col="yellow", pch=16)
text(polar.dat[sel.add[img.seq[idx.range[dontUse]]]],adj=c(-1,-2), labels=idx.range[dontUse])
points(polar.dat[sel.add[img.seq[idx.range[-dontUse]]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[idx.range[-dontUse]]]],adj=c(1,2), labels=idx.range[-dontUse])

## plot next 8 images to pick only:
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[c(1:21)[-dontUse]]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[c(1:21)[-dontUse]]]],adj=2, labels=c(1:21)[-dontUse])

#### edit image selection column
## last image of the so far scored images
no.in.last.img <- tail(org.img.select[subset.v[sel]],1)
## keep already scored good-distance images at the top of the list
sel.clean <- subset.v[sel[sel%!in%scored.and.bad]]
dat$image.select[sel.clean] <- org.img.select[sel.clean]
## put the selected images next in line
dat$image.select[sel.add[img.seq[-dontUse]]] <- PS81b.dat$image.select[sel.add[img.seq[-dontUse]]]+no.in.last.img
## move dontUse images from PS81b to the back of the valid list (between the list and the NAs)
dat$image.select[sel.add[img.seq[dontUse]]] <- PS81b.dat$image.select[sel.add[img.seq[dontUse]]]+2000
## move bad.dist all the way to the back of the list
dat$image.select[subset.v[bad.dist]] <- 6666

#### check output by plotting all images in their sequence
sel.final <- order(dat$image.select[subset.v])[1:100]
plot(polar.dat[subset.v])
points(polar.dat[subset.v[sel.final]], col="blue", pch=15)
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[c(1:100)[-dontUse]]]], col="deepskyblue", pch=16)


##################
#### TRANSECT 118:
## all images in range here
i=2
subset.v <- which(dat$t.ID_temp==levels(dat$t.ID_temp)[i])
subset.v.good <- which(dat$Dist.subs..m.[subset.v]>=1&dat$Dist.subs..m.[subset.v]<=2)
sel <- order(org.img.select[subset.v])[1:t.images[i]]
keep.idx <- which(subset.v[sel]%!in%rmv.idx)

#### order of all next images in line according to PS81b image selection:
already.scored <- order(org.img.select[subset.v])[1:t.images[i]]
bad.dist <- which(PS81b.dat$Dist.subs..m.[subset.v]<1|PS81b.dat$Dist.subs..m.[subset.v]>2)
scored.and.bad <- bad.dist[which(bad.dist%in%already.scored)]
scored.or.bad <- union(already.scored, bad.dist)
sel.add <- subset.v[-scored.or.bad]
## sort the images into a spatially balance sequence
img.seq <- order(PS81b.dat$image.select[sel.add])
## check where the next images are positioned
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[1:50]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[1:50]]],adj=c(1,2))
## select which images are too close to already scored images
dontUse <- c(1,2,4,5,6,7,8,13,14,16,18,19,20,21,22,24,26,28,33,40,42,46,47,55,56,64,68,78)#,9,11,12,14,15,16,18,23,24,25,33,35,36,38,42,47,56,79)
idx.range <- 1:80
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[idx.range[dontUse]]]], col="yellow", pch=16)
text(polar.dat[sel.add[img.seq[idx.range[dontUse]]]],adj=c(-1,-2), labels=idx.range[dontUse])
points(polar.dat[sel.add[img.seq[idx.range[-dontUse]]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[idx.range[-dontUse]]]],adj=c(1,2), labels=idx.range[-dontUse])

## plot next images to pick only:
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[c(1:31)[-dontUse]]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[c(1:31)[-dontUse]]]],adj=2, labels=c(1:31)[-dontUse])

#### edit image selection column
## last image of the so far scored images
no.in.last.img <- tail(org.img.select[subset.v[sel]],1)
## keep already scored good-distance images at the top of the list
sel.clean <- subset.v[sel[sel%!in%scored.and.bad]]
dat$image.select[sel.clean] <- org.img.select[sel.clean]
## put the selected images next in line
dat$image.select[sel.add[img.seq[-dontUse]]] <- PS81b.dat$image.select[sel.add[img.seq[-dontUse]]]+no.in.last.img
## move dontUse images from PS81b to the back of the valid list (between the list and the NAs)
dat$image.select[sel.add[img.seq[dontUse]]] <- PS81b.dat$image.select[sel.add[img.seq[dontUse]]]+2000
## move bad.dist all the way to the back of the list
dat$image.select[subset.v[bad.dist]] <- 6666

#### check output by plotting all images in their sequence
sel.final <- order(dat$image.select[subset.v])[1:100]
plot(polar.dat[subset.v])
points(polar.dat[subset.v[sel.final]], col="blue", pch=15)
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[c(1:100)[-dontUse]]]], col="deepskyblue", pch=16)


##################
#### TRANSECT 159:
## all images in range here
i=3
subset.v <- which(dat$t.ID_temp==levels(dat$t.ID_temp)[i])
subset.v.good <- which(dat$Dist.subs..m.[subset.v]>=1&dat$Dist.subs..m.[subset.v]<=2)
sel <- order(org.img.select[subset.v])[1:t.images[i]]
keep.idx <- which(subset.v[sel]%!in%rmv.idx)

#### order of all next images in line according to PS81b image selection:
already.scored <- order(org.img.select[subset.v])[1:t.images[i]]
bad.dist <- which(PS81b.dat$Dist.subs..m.[subset.v]<1|PS81b.dat$Dist.subs..m.[subset.v]>2)
scored.and.bad <- bad.dist[which(bad.dist%in%already.scored)]
scored.or.bad <- union(already.scored, bad.dist)
sel.add <- subset.v[-scored.or.bad]
## sort the images into a spatially balance sequence
img.seq <- order(PS81b.dat$image.select[sel.add])
## check where the next images are positioned
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[1:50]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[1:50]]],adj=c(1,2))
## select which images are too close to already scored images
dontUse <- c(1,2,3,5,6,11,12,14,15,16,17,19,21,23,24,25,32,37,39,49,66,68,71,72,74,75,76)
idx.range <- 1:80
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[idx.range[dontUse]]]], col="yellow", pch=16)
text(polar.dat[sel.add[img.seq[idx.range[dontUse]]]],adj=c(0.5,-2.5), labels=idx.range[dontUse])
points(polar.dat[sel.add[img.seq[idx.range[-dontUse]]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[idx.range[-dontUse]]]],adj=c(0.5,2.5), labels=idx.range[-dontUse])

## plot next images to pick only:
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[c(1:31)[-dontUse]]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[c(1:31)[-dontUse]]]],adj=2, labels=c(1:31)[-dontUse])

#### edit image selection column
## last image of the so far scored images
no.in.last.img <- tail(org.img.select[subset.v[sel]],1)
## keep already scored good-distance images at the top of the list
sel.clean <- subset.v[sel[sel%!in%scored.and.bad]]
dat$image.select[sel.clean] <- org.img.select[sel.clean]
## put the selected images next in line
dat$image.select[sel.add[img.seq[-dontUse]]] <- PS81b.dat$image.select[sel.add[img.seq[-dontUse]]]+no.in.last.img
## move dontUse images from PS81b to the back of the valid list (between the list and the NAs)
dat$image.select[sel.add[img.seq[dontUse]]] <- PS81b.dat$image.select[sel.add[img.seq[dontUse]]]+2000
## move bad.dist all the way to the back of the list
dat$image.select[subset.v[bad.dist]] <- 6666

#### check output by plotting all images in their sequence
sel.final <- order(dat$image.select[subset.v])[1:100]
plot(polar.dat[subset.v])
points(polar.dat[subset.v[sel.final]], col="blue", pch=15)
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[c(1:100)[-dontUse]]]], col="deepskyblue", pch=16)


##################
#### TRANSECT 160:
rmv.idx <- rmv.idx.full[9]
i=4
subset.v <- which(dat$t.ID_temp==levels(dat$t.ID_temp)[i])
subset.v.good <- which(dat$Dist.subs..m.[subset.v]>=1&dat$Dist.subs..m.[subset.v]<=2)
sel <- order(org.img.select[subset.v])[1:t.images[i]]
keep.idx <- which(subset.v[sel]%!in%rmv.idx)

## plot transect 160 and determine which images should be scored instead of the badly illuminated ones
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[rmv.idx], col="green", pch=17)
#text(polar.dat[rmv.idx,1:2], col="green", adj=-0.5, labels=paste0(dat$transectID[rmv.idx],"_00",org.img.select[rmv.idx]))#, pch=15)
scalebar(500, type="bar")

#### order of all next images in line according to PS81b image selection:
already.scored <- order(org.img.select[subset.v])[1:t.images[i]]
bad.dist <- which(PS81b.dat$Dist.subs..m.[subset.v]<1|PS81b.dat$Dist.subs..m.[subset.v]>2)
scored.and.bad <- bad.dist[which(bad.dist%in%already.scored)]
scored.or.bad <- union(already.scored, bad.dist)
sel.add <- subset.v[-scored.or.bad]
## sort the images into a spatially balance sequence
img.seq <- order(PS81b.dat$image.select[sel.add])
## check where the next images are positioned
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[1:50]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[1:50]]],adj=2)
## select which images are too close to already scored images
dontUse <- c(2,5,6,7,8,11,12,13,14,15,17,18,20,22,23,24,26,27,29,37,40,57)
idx.range <- 1:60
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[idx.range[dontUse]]]], col="yellow", pch=16)
text(polar.dat[sel.add[img.seq[idx.range[dontUse]]]],adj=c(0.5,-3), labels=idx.range[dontUse])
points(polar.dat[sel.add[img.seq[idx.range[-dontUse]]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[idx.range[-dontUse]]]],adj=c(0.5,3), labels=idx.range[-dontUse])

## plot next 8 images to pick only:
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[c(1:21)[-dontUse]]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[c(1:21)[-dontUse]]]],adj=c(0.5,3), labels=c(1:21)[-dontUse])

#### edit image selection column
## last image of the so far scored images
no.in.last.img <- tail(org.img.select[subset.v[sel]],1)
## keep already scored good-distance images at the top of the list
sel.clean <- subset.v[sel[sel%!in%scored.and.bad]]
dat$image.select[sel.clean] <- org.img.select[sel.clean]
## put the selected images next in line
dat$image.select[sel.add[img.seq[-dontUse]]] <- PS81b.dat$image.select[sel.add[img.seq[-dontUse]]]+no.in.last.img
## move dontUse images from PS81b to the back of the valid list (between the list and the NAs)
dat$image.select[sel.add[img.seq[dontUse]]] <- PS81b.dat$image.select[sel.add[img.seq[dontUse]]]+2000
## move bad.dist all the way to the back of the list
dat$image.select[subset.v[bad.dist]] <- 6666

#### check output by plotting all images in their sequence
sel.final <- order(dat$image.select[subset.v])[1:100]
plot(polar.dat[subset.v])
points(polar.dat[subset.v[sel.final]], col="blue", pch=15)
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[c(1:100)[-dontUse]]]], col="deepskyblue", pch=16)


##################
#### TRANSECT 161:
rmv.idx <- rmv.idx.full[10]
i=5
subset.v <- which(dat$t.ID_temp==levels(dat$t.ID_temp)[i])
subset.v.good <- which(dat$Dist.subs..m.[subset.v]>=1&dat$Dist.subs..m.[subset.v]<=2)
sel <- order(org.img.select[subset.v])[1:t.images[i]]
keep.idx <- which(subset.v[sel]%!in%rmv.idx)

## plot transect 160 and determine which images should be scored instead of the badly illuminated ones
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[rmv.idx], col="green", pch=17)
#text(polar.dat[rmv.idx,1:2], col="green", adj=-0.5, labels=paste0(dat$transectID[rmv.idx],"_00",org.img.select[rmv.idx]))#, pch=15)
scalebar(500, type="bar")

#### order of all next images in line according to PS81b image selection:
already.scored <- order(org.img.select[subset.v])[1:t.images[i]]
bad.dist <- which(PS81b.dat$Dist.subs..m.[subset.v]<1|PS81b.dat$Dist.subs..m.[subset.v]>2)
scored.and.bad <- bad.dist[which(bad.dist%in%already.scored)]
scored.or.bad <- union(already.scored, bad.dist)
sel.add <- subset.v[-scored.or.bad]
## sort the images into a spatially balance sequence
img.seq <- order(PS81b.dat$image.select[sel.add])
## check where the next images are positioned
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[1:50]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[1:50]]],adj=2)
## select which images are too close to already scored images
dontUse <- c(1,3,4,5,7,9,10,12,14,15,16,17,19,23,26,27,29,31,34,40)
idx.range <- 1:60
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[idx.range[dontUse]]]], col="yellow", pch=16)
text(polar.dat[sel.add[img.seq[idx.range[dontUse]]]],adj=-3, labels=idx.range[dontUse])
points(polar.dat[sel.add[img.seq[idx.range[-dontUse]]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[idx.range[-dontUse]]]],adj=3, labels=idx.range[-dontUse])

## plot next 8 images to pick only:
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[c(1:21)[-dontUse]]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[c(1:21)[-dontUse]]]],adj=3, labels=c(1:21)[-dontUse])

#### edit image selection column
## last image of the so far scored images
no.in.last.img <- tail(org.img.select[subset.v[sel]],1)
## keep already scored good-distance images at the top of the list
sel.clean <- subset.v[sel[sel%!in%scored.and.bad]]
dat$image.select[sel.clean] <- org.img.select[sel.clean]
## put the selected images next in line
dat$image.select[sel.add[img.seq[-dontUse]]] <- PS81b.dat$image.select[sel.add[img.seq[-dontUse]]]+no.in.last.img
## move dontUse images from PS81b to the back of the valid list (between the list and the NAs)
dat$image.select[sel.add[img.seq[dontUse]]] <- PS81b.dat$image.select[sel.add[img.seq[dontUse]]]+2000
## move bad.dist all the way to the back of the list
dat$image.select[subset.v[bad.dist]] <- 6666

#### check output by plotting all images in their sequence
sel.final <- order(dat$image.select[subset.v])[1:100]
plot(polar.dat[subset.v])
points(polar.dat[subset.v[sel.final]], col="blue", pch=15)
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[c(1:100)[-dontUse]]]], col="deepskyblue", pch=16)


##################
#### TRANSECT 163:
## all images in range here
i=6
subset.v <- which(dat$t.ID_temp==levels(dat$t.ID_temp)[i])
subset.v.good <- which(dat$Dist.subs..m.[subset.v]>=1&dat$Dist.subs..m.[subset.v]<=2)
sel <- order(org.img.select[subset.v])[1:t.images[i]]
keep.idx <- which(subset.v[sel]%!in%rmv.idx)

#### order of all next images in line according to PS81b image selection:
already.scored <- order(org.img.select[subset.v])[1:t.images[i]]
bad.dist <- which(PS81b.dat$Dist.subs..m.[subset.v]<1|PS81b.dat$Dist.subs..m.[subset.v]>2)
scored.and.bad <- bad.dist[which(bad.dist%in%already.scored)]
scored.or.bad <- union(already.scored, bad.dist)
sel.add <- subset.v[-scored.or.bad]
## sort the images into a spatially balance sequence
img.seq <- order(PS81b.dat$image.select[sel.add])
## check where the next images are positioned
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[1:50]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[1:50]]],adj=c(1,2))
## select which images are too close to already scored images
dontUse <- c(1,3,4,5,6,7,10,11,12,13,14,20,21,23,25,26,29,33,34,35,36,39,45,49,61,62,64,67,77)
idx.range <- 1:80
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[idx.range[dontUse]]]], col="yellow", pch=16)
text(polar.dat[sel.add[img.seq[idx.range[dontUse]]]],adj=c(0.5,-2.5), labels=idx.range[dontUse])
points(polar.dat[sel.add[img.seq[idx.range[-dontUse]]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[idx.range[-dontUse]]]],adj=c(0.5,2.5), labels=idx.range[-dontUse])

## plot next images to pick only:
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[c(1:21)[-dontUse]]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[c(1:21)[-dontUse]]]],adj=2, labels=c(1:21)[-dontUse])

#### edit image selection column
## last image of the so far scored images
no.in.last.img <- tail(org.img.select[subset.v[sel]],1)
## keep already scored good-distance images at the top of the list
sel.clean <- subset.v[sel[sel%!in%scored.and.bad]]
dat$image.select[sel.clean] <- org.img.select[sel.clean]
## put the selected images next in line
dat$image.select[sel.add[img.seq[-dontUse]]] <- PS81b.dat$image.select[sel.add[img.seq[-dontUse]]]+no.in.last.img
## move dontUse images from PS81b to the back of the valid list (between the list and the NAs)
dat$image.select[sel.add[img.seq[dontUse]]] <- PS81b.dat$image.select[sel.add[img.seq[dontUse]]]+2000
## move bad.dist all the way to the back of the list
dat$image.select[subset.v[bad.dist]] <- 6666

#### check output by plotting all images in their sequence
sel.final <- order(dat$image.select[subset.v])[1:100]
plot(polar.dat[subset.v])
points(polar.dat[subset.v[sel.final]], col="blue", pch=15)
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[c(1:100)[-dontUse]]]], col="deepskyblue", pch=16)


##################
#### TRANSECT 164:
## all images in range here
i=7
subset.v <- which(dat$t.ID_temp==levels(dat$t.ID_temp)[i])
subset.v.good <- which(dat$Dist.subs..m.[subset.v]>=1&dat$Dist.subs..m.[subset.v]<=2)
sel <- order(org.img.select[subset.v])[1:t.images[i]]
keep.idx <- which(subset.v[sel]%!in%rmv.idx)

#### order of all next images in line according to PS81b image selection:
already.scored <- order(org.img.select[subset.v])[1:t.images[i]]
bad.dist <- which(PS81b.dat$Dist.subs..m.[subset.v]<1|PS81b.dat$Dist.subs..m.[subset.v]>2)
scored.and.bad <- bad.dist[which(bad.dist%in%already.scored)]
scored.or.bad <- union(already.scored, bad.dist)
sel.add <- subset.v[-scored.or.bad]
## sort the images into a spatially balance sequence
img.seq <- order(PS81b.dat$image.select[sel.add])
## check where the next images are positioned
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[1:50]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[1:50]]],adj=2)
## select which images are too close to already scored images
dontUse <- c(1,3,5,6,7,9,10,11,12,16,17,18,19,20,21,23,24,25,27,28,36,42,47,51,57,59,61,62,63,68,70,78)
idx.range <- 1:80
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[idx.range[dontUse]]]], col="yellow", pch=16)
text(polar.dat[sel.add[img.seq[idx.range[dontUse]]]],adj=-2, labels=idx.range[dontUse])
points(polar.dat[sel.add[img.seq[idx.range[-dontUse]]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[idx.range[-dontUse]]]],adj=2, labels=idx.range[-dontUse])

## plot next images to pick only:
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[c(1:21)[-dontUse]]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[c(1:21)[-dontUse]]]],adj=2, labels=c(1:21)[-dontUse])

#### edit image selection column
## last image of the so far scored images
no.in.last.img <- tail(org.img.select[subset.v[sel]],1)
## keep already scored good-distance images at the top of the list
sel.clean <- subset.v[sel[sel%!in%scored.and.bad]]
dat$image.select[sel.clean] <- org.img.select[sel.clean]
## put the selected images next in line
dat$image.select[sel.add[img.seq[-dontUse]]] <- PS81b.dat$image.select[sel.add[img.seq[-dontUse]]]+no.in.last.img
## move dontUse images from PS81b to the back of the valid list (between the list and the NAs)
dat$image.select[sel.add[img.seq[dontUse]]] <- PS81b.dat$image.select[sel.add[img.seq[dontUse]]]+2000
## move bad.dist all the way to the back of the list
dat$image.select[subset.v[bad.dist]] <- 6666

#### check output by plotting all images in their sequence
sel.final <- order(dat$image.select[subset.v])[1:100]
plot(polar.dat[subset.v])
points(polar.dat[subset.v[sel.final]], col="blue", pch=15)
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[c(1:100)[-dontUse]]]], col="deepskyblue", pch=16)


PS81.dat <- dat
save(PS81.dat,total.t.length.v, file="C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/PS81_dat_FINAL.Rdata")




#############################################################################################################
#############################################################################################################

##### ######################################## #####
##### ###### @@@@@@@@@@ PS96 @@@@@@@@@@ ###### #####
##### ######################################## #####

dat <- PS96.dat
total.t.length.v <- PS96_total.t.length.v
## store original image select separately
org.img.select <- PS96.dat$image.select
## wipe image select column
dat$image.select <- NA
## fill with PS96 data for the b-transects
#dat$image.select[3524:14177] <- PS81b.dat$image.select[3524:14177]


## first plot all transects
t.images <- ceiling(total.t.length.v/100)
spatial.dat <- data.frame(dat[,1:2])
coordinates(spatial.dat) <- c("Longitude","Latitude")
proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
polar.dat <- spTransform(spatial.dat, CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
#par(mfrow=c(4,4), mar=c(5,4,2,1))
for(i in 1:length(levels(dat$transectID))){
  subset.v <- which(dat$transectID==levels(dat$transectID)[i])
  subset.v.good <- which(dat$Height..m.[subset.v]>=1&dat$Height..m.[subset.v]<=2)
  ## select the correct number of images for each transect
  sel <- order(org.img.select[subset.v])[1:t.images[i]]
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
store.dat.sel <- dat.sel

## check which images from the transects are further away than 2m range
not.in.range <- which(dat$Height..m.[store.dat.sel]>2)
dat[store.dat.sel,][not.in.range,]
rmv.idx.full <- as.numeric(rownames(dat[store.dat.sel,][not.in.range,]))
dat[rmv.idx.full,]

plot(dat$Height..m.[store.dat.sel])
text(dat$Height..m.[store.dat.sel],labels=paste0(dat$transectID[store.dat.sel],"_00",org.img.select[store.dat.sel]))


##################
#### TRANSECT 001:
rmv.idx <- rmv.idx.full[1:2]
i=1
subset.v <- which(dat$transectID==levels(dat$transectID)[i])
subset.v.good <- which(dat$Height..m.[subset.v]>=1&dat$Height..m.[subset.v]<=2)
sel <- order(org.img.select[subset.v])[1:t.images[i]]
keep.idx <- which(subset.v[sel]%!in%rmv.idx)

## plot transect 001 and determine which images should be scored instead of the badly illuminated ones
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[rmv.idx], col="green", pch=17)
#text(polar.dat[rmv.idx,1:2], col="green", adj=-0.5, labels=paste0(dat$transectID[rmv.idx],"_00",org.img.select[rmv.idx]))#, pch=15)
scalebar(500, type="bar")

#### order of all next images in line according to PS81b image selection:
already.scored <- order(org.img.select[subset.v])[1:t.images[i]]
bad.dist <- which(PS96b.dat$Height..m.[subset.v]<1|PS96b.dat$Height..m.[subset.v]>2)
scored.and.bad <- bad.dist[which(bad.dist%in%already.scored)]
scored.or.bad <- union(already.scored, bad.dist)
sel.add <- subset.v[-scored.or.bad]
## sort the images into a spatially balance sequence
img.seq <- order(PS96b.dat$image.select[sel.add])
#img.seq <- order(PS96.dat$image.select[sel.add])
## check where the next images are positioned
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[1:40]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[1:40]]],adj=2)
## select which images are too close to already scored images
dontUse <- c(1:8,10,13,19,28,29)
idx.range <- 1:50
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[idx.range[dontUse]]]], col="yellow", pch=16)
text(polar.dat[sel.add[img.seq[idx.range[dontUse]]]],adj=-2, labels=idx.range[dontUse])
points(polar.dat[sel.add[img.seq[idx.range[-dontUse]]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[idx.range[-dontUse]]]],adj=2, labels=idx.range[-dontUse])

## plot next 8 images to pick only:
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[c(1:21)[-dontUse]]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[c(1:21)[-dontUse]]]],adj=2, labels=c(1:21)[-dontUse])

#### edit image selection column
## last image of the so far scored images
no.in.last.img <- tail(org.img.select[subset.v[sel]],1)
## keep already scored good-distance images at the top of the list
sel.clean <- subset.v[sel[sel%!in%scored.and.bad]]
dat$image.select[sel.clean] <- org.img.select[sel.clean]
## put the selected images next in line
dat$image.select[sel.add[img.seq[-dontUse]]] <- PS96b.dat$image.select[sel.add[img.seq[-dontUse]]]+no.in.last.img
## move dontUse images from PS81b to the back of the valid list (between the list and the NAs)
dat$image.select[sel.add[img.seq[dontUse]]] <- PS96b.dat$image.select[sel.add[img.seq[dontUse]]]+2000
## move bad.dist all the way to the back of the list
dat$image.select[subset.v[bad.dist]] <- 6666

#### check output by plotting all images in their sequence
sel.final <- order(dat$image.select[subset.v])[1:100]
plot(polar.dat[subset.v])
points(polar.dat[subset.v[sel.final]], col="blue", pch=15)
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[c(1:100)[-dontUse]]]], col="deepskyblue", pch=16)


##################
#### TRANSECT 007:
i=2
subset.v <- which(dat$transectID==levels(dat$transectID)[i])
subset.v.good <- which(dat$Height..m.[subset.v]>=1&dat$Height..m.[subset.v]<=2)
sel <- order(org.img.select[subset.v])[1:t.images[i]]
keep.idx <- which(subset.v[sel]%!in%rmv.idx)

#### order of all next images in line according to PS81b image selection:
already.scored <- order(org.img.select[subset.v])[1:t.images[i]]
bad.dist <- which(PS96b.dat$Height..m.[subset.v]<1|PS96b.dat$Height..m.[subset.v]>2)
scored.and.bad <- bad.dist[which(bad.dist%in%already.scored)]
scored.or.bad <- union(already.scored, bad.dist)
sel.add <- subset.v[-scored.or.bad]
## sort the images into a spatially balance sequence
img.seq <- order(PS96b.dat$image.select[sel.add])
#img.seq <- order(PS96.dat$image.select[sel.add])
## check where the next images are positioned
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[1:40]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[1:40]]],adj=2)
## select which images are too close to already scored images
dontUse <- c(3,6,8,9,10,11,16,20,21,23,24,25,26,30,34,37,56)
idx.range <- 1:60
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[idx.range[dontUse]]]], col="yellow", pch=16)
text(polar.dat[sel.add[img.seq[idx.range[dontUse]]]],adj=-2, labels=idx.range[dontUse])
points(polar.dat[sel.add[img.seq[idx.range[-dontUse]]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[idx.range[-dontUse]]]],adj=2, labels=idx.range[-dontUse])

## plot next 8 images to pick only:
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[c(1:21)[-dontUse]]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[c(1:21)[-dontUse]]]],adj=2, labels=c(1:21)[-dontUse])

#### edit image selection column
## last image of the so far scored images
no.in.last.img <- tail(org.img.select[subset.v[sel]],1)
## keep already scored good-distance images at the top of the list
sel.clean <- subset.v[sel[sel%!in%scored.and.bad]]
dat$image.select[sel.clean] <- org.img.select[sel.clean]
## put the selected images next in line
dat$image.select[sel.add[img.seq[-dontUse]]] <- PS96b.dat$image.select[sel.add[img.seq[-dontUse]]]+no.in.last.img
## move dontUse images from PS81b to the back of the valid list (between the list and the NAs)
dat$image.select[sel.add[img.seq[dontUse]]] <- PS96b.dat$image.select[sel.add[img.seq[dontUse]]]+2000
## move bad.dist all the way to the back of the list
dat$image.select[subset.v[bad.dist]] <- 6666

#### check output by plotting all images in their sequence
sel.final <- order(dat$image.select[subset.v])[1:100]
plot(polar.dat[subset.v])
points(polar.dat[subset.v[sel.final]], col="blue", pch=15)
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[c(1:100)[-dontUse]]]], col="deepskyblue", pch=16)


##################
#### TRANSECT 008:
i=3
subset.v <- which(dat$transectID==levels(dat$transectID)[i])
subset.v.good <- which(dat$Height..m.[subset.v]>=1&dat$Height..m.[subset.v]<=2)
sel <- order(org.img.select[subset.v])[1:t.images[i]]
keep.idx <- which(subset.v[sel]%!in%rmv.idx)

#### order of all next images in line according to PS81b image selection:
already.scored <- order(org.img.select[subset.v])[1:t.images[i]]
bad.dist <- which(PS96b.dat$Height..m.[subset.v]<1|PS96b.dat$Height..m.[subset.v]>2)
scored.and.bad <- bad.dist[which(bad.dist%in%already.scored)]
scored.or.bad <- union(already.scored, bad.dist)
sel.add <- subset.v[-scored.or.bad]
## sort the images into a spatially balance sequence
img.seq <- order(PS96b.dat$image.select[sel.add])
#img.seq <- order(PS96.dat$image.select[sel.add])
## check where the next images are positioned
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[1:30]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[1:30]]],adj=2)
## select which images are too close to already scored images
dontUse <- c(1:5,7,9,10,11,16,17,20,21,22)
idx.range <- 1:40
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[idx.range[dontUse]]]], col="yellow", pch=16)
text(polar.dat[sel.add[img.seq[idx.range[dontUse]]]],adj=-2, labels=idx.range[dontUse])
points(polar.dat[sel.add[img.seq[idx.range[-dontUse]]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[idx.range[-dontUse]]]],adj=2, labels=idx.range[-dontUse])

## plot next 8 images to pick only:
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[c(1:21)[-dontUse]]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[c(1:21)[-dontUse]]]],adj=2, labels=c(1:21)[-dontUse])

#### edit image selection column
## last image of the so far scored images
no.in.last.img <- tail(org.img.select[subset.v[sel]],1)
## keep already scored good-distance images at the top of the list
sel.clean <- subset.v[sel[sel%!in%scored.and.bad]]
dat$image.select[sel.clean] <- org.img.select[sel.clean]
## put the selected images next in line
dat$image.select[sel.add[img.seq[-dontUse]]] <- PS96b.dat$image.select[sel.add[img.seq[-dontUse]]]+no.in.last.img
## move dontUse images from PS81b to the back of the valid list (between the list and the NAs)
dat$image.select[sel.add[img.seq[dontUse]]] <- PS96b.dat$image.select[sel.add[img.seq[dontUse]]]+2000
## move bad.dist all the way to the back of the list
dat$image.select[subset.v[bad.dist]] <- 6666

#### check output by plotting all images in their sequence
sel.final <- order(dat$image.select[subset.v])[1:100]
plot(polar.dat[subset.v])
points(polar.dat[subset.v[sel.final]], col="blue", pch=15)
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[c(1:100)[-dontUse]]]], col="deepskyblue", pch=16)


##################
#### TRANSECT 010:
i=4
subset.v <- which(dat$transectID==levels(dat$transectID)[i])
subset.v.good <- which(dat$Height..m.[subset.v]>=1&dat$Height..m.[subset.v]<=2)
sel <- order(org.img.select[subset.v])[1:t.images[i]]
keep.idx <- which(subset.v[sel]%!in%rmv.idx)

#### order of all next images in line according to PS81b image selection:
already.scored <- order(org.img.select[subset.v])[1:t.images[i]]
bad.dist <- which(PS96b.dat$Height..m.[subset.v]<1|PS96b.dat$Height..m.[subset.v]>2)
scored.and.bad <- bad.dist[which(bad.dist%in%already.scored)]
scored.or.bad <- union(already.scored, bad.dist)
sel.add <- subset.v[-scored.or.bad]
## sort the images into a spatially balance sequence
img.seq <- order(PS96b.dat$image.select[sel.add])
#img.seq <- order(PS96.dat$image.select[sel.add])
## check where the next images are positioned
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[1:30]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[1:30]]],adj=2)
## select which images are too close to already scored images
dontUse <- c(2,3,5,7:10,12,14,16,17,19,20,26,28,31,36,37)
idx.range <- 1:40
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[idx.range[dontUse]]]], col="yellow", pch=16)
text(polar.dat[sel.add[img.seq[idx.range[dontUse]]]],adj=-2, labels=idx.range[dontUse])
points(polar.dat[sel.add[img.seq[idx.range[-dontUse]]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[idx.range[-dontUse]]]],adj=2, labels=idx.range[-dontUse])

## plot next 8 images to pick only:
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[c(1:21)[-dontUse]]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[c(1:21)[-dontUse]]]],adj=2, labels=c(1:21)[-dontUse])

#### edit image selection column
## last image of the so far scored images
no.in.last.img <- tail(org.img.select[subset.v[sel]],1)
## keep already scored good-distance images at the top of the list
sel.clean <- subset.v[sel[sel%!in%scored.and.bad]]
dat$image.select[sel.clean] <- org.img.select[sel.clean]
## put the selected images next in line
dat$image.select[sel.add[img.seq[-dontUse]]] <- PS96b.dat$image.select[sel.add[img.seq[-dontUse]]]+no.in.last.img
## move dontUse images from PS81b to the back of the valid list (between the list and the NAs)
dat$image.select[sel.add[img.seq[dontUse]]] <- PS96b.dat$image.select[sel.add[img.seq[dontUse]]]+2000
## move bad.dist all the way to the back of the list
dat$image.select[subset.v[bad.dist]] <- 6666

#### check output by plotting all images in their sequence
sel.final <- order(dat$image.select[subset.v])[1:100]
plot(polar.dat[subset.v])
points(polar.dat[subset.v[sel.final]], col="blue", pch=15)
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[c(1:100)[-dontUse]]]], col="deepskyblue", pch=16)


##################
#### TRANSECT 026:
i=5
subset.v <- which(dat$transectID==levels(dat$transectID)[i])
subset.v.good <- which(dat$Height..m.[subset.v]>=1&dat$Height..m.[subset.v]<=2)
sel <- order(org.img.select[subset.v])[1:t.images[i]]
keep.idx <- which(subset.v[sel]%!in%rmv.idx)

#### order of all next images in line according to PS81b image selection:
already.scored <- order(org.img.select[subset.v])[1:t.images[i]]
bad.dist <- which(PS96b.dat$Height..m.[subset.v]<1|PS96b.dat$Height..m.[subset.v]>2)
scored.and.bad <- bad.dist[which(bad.dist%in%already.scored)]
scored.or.bad <- union(already.scored, bad.dist)
sel.add <- subset.v[-scored.or.bad]
## sort the images into a spatially balance sequence
img.seq <- order(PS96b.dat$image.select[sel.add])
#img.seq <- order(PS96.dat$image.select[sel.add])
## check where the next images are positioned
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[1:30]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[1:30]]],adj=2)
## select which images are too close to already scored images
dontUse <- c(2:6,7,17,21,29,32,40,44,46,47)
idx.range <- 1:50
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[idx.range[dontUse]]]], col="yellow", pch=16)
text(polar.dat[sel.add[img.seq[idx.range[dontUse]]]],adj=-2, labels=idx.range[dontUse])
points(polar.dat[sel.add[img.seq[idx.range[-dontUse]]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[idx.range[-dontUse]]]],adj=2, labels=idx.range[-dontUse])

## plot next 8 images to pick only:
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[c(1:21)[-dontUse]]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[c(1:21)[-dontUse]]]],adj=2, labels=c(1:21)[-dontUse])

#### edit image selection column
## last image of the so far scored images
no.in.last.img <- tail(org.img.select[subset.v[sel]],1)
## keep already scored good-distance images at the top of the list
sel.clean <- subset.v[sel[sel%!in%scored.and.bad]]
dat$image.select[sel.clean] <- org.img.select[sel.clean]
## put the selected images next in line
dat$image.select[sel.add[img.seq[-dontUse]]] <- PS96b.dat$image.select[sel.add[img.seq[-dontUse]]]+no.in.last.img
## move dontUse images from PS81b to the back of the valid list (between the list and the NAs)
dat$image.select[sel.add[img.seq[dontUse]]] <- PS96b.dat$image.select[sel.add[img.seq[dontUse]]]+2000
## move bad.dist all the way to the back of the list
dat$image.select[subset.v[bad.dist]] <- 6666

#### check output by plotting all images in their sequence
sel.final <- order(dat$image.select[subset.v])[1:100]
plot(polar.dat[subset.v])
points(polar.dat[subset.v[sel.final]], col="blue", pch=15)
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[c(1:100)[-dontUse]]]], col="deepskyblue", pch=16)


##################
#### TRANSECT 027:
i=6
subset.v <- which(dat$transectID==levels(dat$transectID)[i])
subset.v.good <- which(dat$Height..m.[subset.v]>=1&dat$Height..m.[subset.v]<=2)
sel <- order(org.img.select[subset.v])[1:t.images[i]]
keep.idx <- which(subset.v[sel]%!in%rmv.idx)

#### order of all next images in line according to PS81b image selection:
already.scored <- order(org.img.select[subset.v])[1:t.images[i]]
bad.dist <- which(PS96b.dat$Height..m.[subset.v]<1|PS96b.dat$Height..m.[subset.v]>2)
scored.and.bad <- bad.dist[which(bad.dist%in%already.scored)]
scored.or.bad <- union(already.scored, bad.dist)
sel.add <- subset.v[-scored.or.bad]
## sort the images into a spatially balance sequence
img.seq <- order(PS96b.dat$image.select[sel.add])
#img.seq <- order(PS96.dat$image.select[sel.add])
## check where the next images are positioned
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[1:30]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[1:30]]],adj=2)
## select which images are too close to already scored images
dontUse <- c(3,6,8,22)
idx.range <- 1:40
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[idx.range[dontUse]]]], col="yellow", pch=16)
text(polar.dat[sel.add[img.seq[idx.range[dontUse]]]],adj=-2, labels=idx.range[dontUse])
points(polar.dat[sel.add[img.seq[idx.range[-dontUse]]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[idx.range[-dontUse]]]],adj=2, labels=idx.range[-dontUse])

## plot next 8 images to pick only:
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[c(1:21)[-dontUse]]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[c(1:21)[-dontUse]]]],adj=2, labels=c(1:21)[-dontUse])

#### edit image selection column
## last image of the so far scored images
no.in.last.img <- tail(org.img.select[subset.v[sel]],1)
## keep already scored good-distance images at the top of the list
sel.clean <- subset.v[sel[sel%!in%scored.and.bad]]
dat$image.select[sel.clean] <- org.img.select[sel.clean]
## put the selected images next in line
dat$image.select[sel.add[img.seq[-dontUse]]] <- PS96b.dat$image.select[sel.add[img.seq[-dontUse]]]+no.in.last.img
## move dontUse images from PS81b to the back of the valid list (between the list and the NAs)
dat$image.select[sel.add[img.seq[dontUse]]] <- PS96b.dat$image.select[sel.add[img.seq[dontUse]]]+2000
## move bad.dist all the way to the back of the list
dat$image.select[subset.v[bad.dist]] <- 6666

#### check output by plotting all images in their sequence
sel.final <- order(dat$image.select[subset.v])[1:100]
plot(polar.dat[subset.v])
points(polar.dat[subset.v[sel.final]], col="blue", pch=15)
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[c(1:100)[-dontUse]]]], col="deepskyblue", pch=16)



##################
#### TRANSECT 037:
i=7
subset.v <- which(dat$transectID==levels(dat$transectID)[i])
subset.v.good <- which(dat$Height..m.[subset.v]>=1&dat$Height..m.[subset.v]<=2)
sel <- order(org.img.select[subset.v])[1:t.images[i]]
keep.idx <- which(subset.v[sel]%!in%rmv.idx)

#### order of all next images in line according to PS81b image selection:
already.scored <- order(org.img.select[subset.v])[1:t.images[i]]
bad.dist <- which(PS96b.dat$Height..m.[subset.v]<1|PS96b.dat$Height..m.[subset.v]>2)
scored.and.bad <- bad.dist[which(bad.dist%in%already.scored)]
scored.or.bad <- union(already.scored, bad.dist)
sel.add <- subset.v[-scored.or.bad]
## sort the images into a spatially balance sequence
img.seq <- order(PS96b.dat$image.select[sel.add])
#img.seq <- order(PS96.dat$image.select[sel.add])
## check where the next images are positioned
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[1:30]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[1:30]]],adj=2)
## select which images are too close to already scored images
dontUse <- c(6,10,20,24,29,33)
idx.range <- 1:40
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[idx.range[dontUse]]]], col="yellow", pch=16)
text(polar.dat[sel.add[img.seq[idx.range[dontUse]]]],adj=-2, labels=idx.range[dontUse])
points(polar.dat[sel.add[img.seq[idx.range[-dontUse]]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[idx.range[-dontUse]]]],adj=2, labels=idx.range[-dontUse])

## plot next 8 images to pick only:
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[c(1:21)[-dontUse]]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[c(1:21)[-dontUse]]]],adj=2, labels=c(1:21)[-dontUse])

#### edit image selection column
## last image of the so far scored images
no.in.last.img <- tail(org.img.select[subset.v[sel]],1)
## keep already scored good-distance images at the top of the list
sel.clean <- subset.v[sel[sel%!in%scored.and.bad]]
dat$image.select[sel.clean] <- org.img.select[sel.clean]
## put the selected images next in line
dat$image.select[sel.add[img.seq[-dontUse]]] <- PS96b.dat$image.select[sel.add[img.seq[-dontUse]]]+no.in.last.img
## move dontUse images from PS81b to the back of the valid list (between the list and the NAs)
dat$image.select[sel.add[img.seq[dontUse]]] <- PS96b.dat$image.select[sel.add[img.seq[dontUse]]]+2000
## move bad.dist all the way to the back of the list
dat$image.select[subset.v[bad.dist]] <- 6666

#### check output by plotting all images in their sequence
sel.final <- order(dat$image.select[subset.v])[1:100]
plot(polar.dat[subset.v])
points(polar.dat[subset.v[sel.final]], col="blue", pch=15)
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[c(1:100)[-dontUse]]]], col="deepskyblue", pch=16)


##################
#### TRANSECT 048:
rmv.idx <- rmv.idx.full[3]
i=8
subset.v <- which(dat$transectID==levels(dat$transectID)[i])
subset.v.good <- which(dat$Height..m.[subset.v]>=1&dat$Height..m.[subset.v]<=2)
sel <- order(org.img.select[subset.v])[1:t.images[i]]
keep.idx <- which(subset.v[sel]%!in%rmv.idx)

## plot transect 048 and determine which images should be scored instead of the badly illuminated ones
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[rmv.idx], col="green", pch=17)
#text(polar.dat[rmv.idx,1:2], col="green", adj=-0.5, labels=paste0(dat$transectID[rmv.idx],"_00",org.img.select[rmv.idx]))#, pch=15)
scalebar(500, type="bar")

#### order of all next images in line according to PS81b image selection:
already.scored <- order(org.img.select[subset.v])[1:t.images[i]]
bad.dist <- which(PS96b.dat$Height..m.[subset.v]<1|PS96b.dat$Height..m.[subset.v]>2)
scored.and.bad <- bad.dist[which(bad.dist%in%already.scored)]
scored.or.bad <- union(already.scored, bad.dist)
sel.add <- subset.v[-scored.or.bad]
## sort the images into a spatially balance sequence
img.seq <- order(PS96b.dat$image.select[sel.add])
#img.seq <- order(PS96.dat$image.select[sel.add])
## check where the next images are positioned
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[1:30]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[1:30]]],adj=c(0.5,-2))
## select which images are too close to already scored images
dontUse <- c(4,5,6,7,8,9,10,14,17,24,26,28)
idx.range <- 1:40
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[idx.range[dontUse]]]], col="yellow", pch=16)
text(polar.dat[sel.add[img.seq[idx.range[dontUse]]]],adj=c(0.5,-2), labels=idx.range[dontUse])
points(polar.dat[sel.add[img.seq[idx.range[-dontUse]]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[idx.range[-dontUse]]]],adj=c(0.5,2), labels=idx.range[-dontUse])

## plot next 8 images to pick only:
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[c(1:21)[-dontUse]]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[c(1:21)[-dontUse]]]],adj=c(0.5,2), labels=c(1:21)[-dontUse])

#### edit image selection column
## last image of the so far scored images
no.in.last.img <- tail(org.img.select[subset.v[sel]],1)
## keep already scored good-distance images at the top of the list
sel.clean <- subset.v[sel[sel%!in%scored.and.bad]]
dat$image.select[sel.clean] <- org.img.select[sel.clean]
## put the selected images next in line
dat$image.select[sel.add[img.seq[-dontUse]]] <- PS96b.dat$image.select[sel.add[img.seq[-dontUse]]]+no.in.last.img
## move dontUse images from PS81b to the back of the valid list (between the list and the NAs)
dat$image.select[sel.add[img.seq[dontUse]]] <- PS96b.dat$image.select[sel.add[img.seq[dontUse]]]+2000
## move bad.dist all the way to the back of the list
dat$image.select[subset.v[bad.dist]] <- 6666

#### check output by plotting all images in their sequence
sel.final <- order(dat$image.select[subset.v])[1:100]
plot(polar.dat[subset.v])
points(polar.dat[subset.v[sel.final]], col="blue", pch=15)
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[c(1:100)[-dontUse]]]], col="deepskyblue", pch=16)


##################
#### TRANSECT 057:
rmv.idx <- rmv.idx.full[4]#:10]
i=9
subset.v <- which(dat$transectID==levels(dat$transectID)[i])
subset.v.good <- which(dat$Height..m.[subset.v]>=1&dat$Height..m.[subset.v]<=2)
sel <- order(org.img.select[subset.v])[1:t.images[i]]
keep.idx <- which(subset.v[sel]%!in%rmv.idx)

## plot transect 057 and determine which images should be scored instead of the badly illuminated ones
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[rmv.idx], col="green", pch=17)
#text(polar.dat[rmv.idx,1:2], col="green", adj=-0.5, labels=paste0(dat$transectID[rmv.idx],"_00",org.img.select[rmv.idx]))#, pch=15)
scalebar(500, type="bar")

#### order of all next images in line according to PS81b image selection:
already.scored <- order(org.img.select[subset.v])[1:t.images[i]]
bad.dist <- which(PS96b.dat$Height..m.[subset.v]<1|PS96b.dat$Height..m.[subset.v]>2)
scored.and.bad <- bad.dist[which(bad.dist%in%already.scored)]
scored.or.bad <- union(already.scored, bad.dist)
sel.add <- subset.v[-scored.or.bad]
## sort the images into a spatially balance sequence
img.seq <- order(PS96b.dat$image.select[sel.add])
#img.seq <- order(PS96.dat$image.select[sel.add])
## check where the next images are positioned
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[1:20]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[1:20]]],adj=c(0.5,-2))
## select which images are too close to already scored images
dontUse <- c(2,4,5,6,7,8,10)
idx.range <- 1:30
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[idx.range[dontUse]]]], col="yellow", pch=16)
text(polar.dat[sel.add[img.seq[idx.range[dontUse]]]],adj=c(0.5,-2), labels=idx.range[dontUse])
points(polar.dat[sel.add[img.seq[idx.range[-dontUse]]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[idx.range[-dontUse]]]],adj=c(0.5,2), labels=idx.range[-dontUse])

## plot next 8 images to pick only:
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[c(1:21)[-dontUse]]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[c(1:21)[-dontUse]]]],adj=c(0.5,2), labels=c(1:21)[-dontUse])

#### edit image selection column
## last image of the so far scored images
no.in.last.img <- tail(org.img.select[subset.v[sel]],1)
## keep already scored good-distance images at the top of the list. ONLY HERE: REMOVE FIRST IMAGE ONLY, NOT ALL 7
#sel.clean <- subset.v[sel[sel%!in%scored.and.bad]]
sel.clean <- subset.v[sel[sel%!in%scored.and.bad[1]]]
dat$image.select[sel.clean] <- org.img.select[sel.clean]
## put the selected images next in line
dat$image.select[sel.add[img.seq[-dontUse]]] <- PS96b.dat$image.select[sel.add[img.seq[-dontUse]]]+no.in.last.img
## move dontUse images from PS81b to the back of the valid list (between the list and the NAs)
dat$image.select[sel.add[img.seq[dontUse]]] <- PS96b.dat$image.select[sel.add[img.seq[dontUse]]]+2000
## move bad.dist all the way to the back of the list
#dat$image.select[subset.v[bad.dist]] <- 6666
dat$image.select[subset.v[which(is.na(dat$image.select[subset.v]))]] <- 6666

#### check output by plotting all images in their sequence
sel.final <- order(dat$image.select[subset.v])[1:100]
plot(polar.dat[subset.v])
points(polar.dat[subset.v[sel.final]], col="blue", pch=15)
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[c(1:58)[-dontUse]]]], col="deepskyblue", pch=16)


##################
#### TRANSECT 061:
i=10
subset.v <- which(dat$transectID==levels(dat$transectID)[i])
subset.v.good <- which(dat$Height..m.[subset.v]>=1&dat$Height..m.[subset.v]<=2)
sel <- order(org.img.select[subset.v])[1:t.images[i]]
keep.idx <- which(subset.v[sel]%!in%rmv.idx)

#### order of all next images in line according to PS81b image selection:
already.scored <- order(org.img.select[subset.v])[1:t.images[i]]
bad.dist <- which(PS96b.dat$Height..m.[subset.v]<1|PS96b.dat$Height..m.[subset.v]>2)
scored.and.bad <- bad.dist[which(bad.dist%in%already.scored)]
scored.or.bad <- union(already.scored, bad.dist)
sel.add <- subset.v[-scored.or.bad]
## sort the images into a spatially balance sequence
img.seq <- order(PS96b.dat$image.select[sel.add])
#img.seq <- order(PS96.dat$image.select[sel.add])
## check where the next images are positioned
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[1:20]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[1:20]]],adj=c(0.5,2))
## select which images are too close to already scored images
dontUse <- c(2,3,6,7,8,10,12,14,18,20,23,24,29,33,37)
idx.range <- 1:40
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[idx.range[dontUse]]]], col="yellow", pch=16)
text(polar.dat[sel.add[img.seq[idx.range[dontUse]]]],adj=c(0.5,-2), labels=idx.range[dontUse])
points(polar.dat[sel.add[img.seq[idx.range[-dontUse]]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[idx.range[-dontUse]]]],adj=c(0.5,2), labels=idx.range[-dontUse])

## plot next 8 images to pick only:
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[c(1:21)[-dontUse]]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[c(1:21)[-dontUse]]]],adj=c(0.5,2), labels=c(1:21)[-dontUse])

#### edit image selection column
## last image of the so far scored images
no.in.last.img <- tail(org.img.select[subset.v[sel]],1)
## keep already scored good-distance images at the top of the list
sel.clean <- subset.v[sel[sel%!in%scored.and.bad]]
dat$image.select[sel.clean] <- org.img.select[sel.clean]
## put the selected images next in line
dat$image.select[sel.add[img.seq[-dontUse]]] <- PS96b.dat$image.select[sel.add[img.seq[-dontUse]]]+no.in.last.img
## move dontUse images from PS81b to the back of the valid list (between the list and the NAs)
dat$image.select[sel.add[img.seq[dontUse]]] <- PS96b.dat$image.select[sel.add[img.seq[dontUse]]]+2000
## move bad.dist all the way to the back of the list
dat$image.select[subset.v[bad.dist]] <- 6666

#### check output by plotting all images in their sequence
sel.final <- order(dat$image.select[subset.v])[1:100]
plot(polar.dat[subset.v])
points(polar.dat[subset.v[sel.final]], col="blue", pch=15)
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[c(1:100)[-dontUse]]]], col="deepskyblue", pch=16)


##################
#### TRANSECT 072:
i=11
subset.v <- which(dat$transectID==levels(dat$transectID)[i])
subset.v.good <- which(dat$Height..m.[subset.v]>=1&dat$Height..m.[subset.v]<=2)
sel <- order(org.img.select[subset.v])[1:t.images[i]]
keep.idx <- which(subset.v[sel]%!in%rmv.idx)

#### order of all next images in line according to PS81b image selection:
already.scored <- order(org.img.select[subset.v])[1:t.images[i]]
bad.dist <- which(PS96b.dat$Height..m.[subset.v]<1|PS96b.dat$Height..m.[subset.v]>2)
scored.and.bad <- bad.dist[which(bad.dist%in%already.scored)]
scored.or.bad <- union(already.scored, bad.dist)
sel.add <- subset.v[-scored.or.bad]
## sort the images into a spatially balance sequence
img.seq <- order(PS96b.dat$image.select[sel.add])
#img.seq <- order(PS96.dat$image.select[sel.add])
## check where the next images are positioned
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[1:20]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[1:20]]],adj=c(2,0))
## select which images are too close to already scored images
dontUse <- c(2,5:9,18:20,30,34)
idx.range <- 1:40
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[idx.range[dontUse]]]], col="yellow", pch=16)
text(polar.dat[sel.add[img.seq[idx.range[dontUse]]]],adj=c(-2,0), labels=idx.range[dontUse])
points(polar.dat[sel.add[img.seq[idx.range[-dontUse]]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[idx.range[-dontUse]]]],adj=c(2,0), labels=idx.range[-dontUse])

## plot next 8 images to pick only:
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[c(1:21)[-dontUse]]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[c(1:21)[-dontUse]]]],adj=c(0.5,2), labels=c(1:21)[-dontUse])

#### edit image selection column
## last image of the so far scored images
no.in.last.img <- tail(org.img.select[subset.v[sel]],1)
## keep already scored good-distance images at the top of the list
sel.clean <- subset.v[sel[sel%!in%scored.and.bad]]
dat$image.select[sel.clean] <- org.img.select[sel.clean]
## put the selected images next in line
dat$image.select[sel.add[img.seq[-dontUse]]] <- PS96b.dat$image.select[sel.add[img.seq[-dontUse]]]+no.in.last.img
## move dontUse images from PS81b to the back of the valid list (between the list and the NAs)
dat$image.select[sel.add[img.seq[dontUse]]] <- PS96b.dat$image.select[sel.add[img.seq[dontUse]]]+2000
## move bad.dist all the way to the back of the list
dat$image.select[subset.v[bad.dist]] <- 6666

#### check output by plotting all images in their sequence
sel.final <- order(dat$image.select[subset.v])[1:100]
plot(polar.dat[subset.v])
points(polar.dat[subset.v[sel.final]], col="blue", pch=15)
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[c(1:100)[-dontUse]]]], col="deepskyblue", pch=16)


##################
#### TRANSECT 090:
i=12
subset.v <- which(dat$transectID==levels(dat$transectID)[i])
subset.v.good <- which(dat$Height..m.[subset.v]>=1&dat$Height..m.[subset.v]<=2)
sel <- order(org.img.select[subset.v])[1:t.images[i]]
keep.idx <- which(subset.v[sel]%!in%rmv.idx)

#### order of all next images in line according to PS81b image selection:
already.scored <- order(org.img.select[subset.v])[1:t.images[i]]
bad.dist <- which(PS96b.dat$Height..m.[subset.v]<1|PS96b.dat$Height..m.[subset.v]>2)
scored.and.bad <- bad.dist[which(bad.dist%in%already.scored)]
scored.or.bad <- union(already.scored, bad.dist)
sel.add <- subset.v[-scored.or.bad]
## sort the images into a spatially balance sequence
img.seq <- order(PS96b.dat$image.select[sel.add])
#img.seq <- order(PS96.dat$image.select[sel.add])
## check where the next images are positioned
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[1:20]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[1:20]]],adj=c(2,0))
## select which images are too close to already scored images
dontUse <- c(2:10,14,15,20,25,31,32)
idx.range <- 1:40
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[idx.range[dontUse]]]], col="yellow", pch=16)
text(polar.dat[sel.add[img.seq[idx.range[dontUse]]]],adj=c(-2,0), labels=idx.range[dontUse])
points(polar.dat[sel.add[img.seq[idx.range[-dontUse]]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[idx.range[-dontUse]]]],adj=c(2,0), labels=idx.range[-dontUse])

## plot next 8 images to pick only:
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[c(1:21)[-dontUse]]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[c(1:21)[-dontUse]]]],adj=c(0.5,2), labels=c(1:21)[-dontUse])

#### edit image selection column
## last image of the so far scored images
no.in.last.img <- tail(org.img.select[subset.v[sel]],1)
## keep already scored good-distance images at the top of the list
sel.clean <- subset.v[sel[sel%!in%scored.and.bad]]
dat$image.select[sel.clean] <- org.img.select[sel.clean]
## put the selected images next in line
dat$image.select[sel.add[img.seq[-dontUse]]] <- PS96b.dat$image.select[sel.add[img.seq[-dontUse]]]+no.in.last.img
## move dontUse images from PS81b to the back of the valid list (between the list and the NAs)
dat$image.select[sel.add[img.seq[dontUse]]] <- PS96b.dat$image.select[sel.add[img.seq[dontUse]]]+2000
## move bad.dist all the way to the back of the list
dat$image.select[subset.v[bad.dist]] <- 6666

#### check output by plotting all images in their sequence
sel.final <- order(dat$image.select[subset.v])[1:100]
plot(polar.dat[subset.v])
points(polar.dat[subset.v[sel.final]], col="blue", pch=15)
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[c(1:100)[-dontUse]]]], col="deepskyblue", pch=16)


##################
#### TRANSECT 106:
i=13
subset.v <- which(dat$transectID==levels(dat$transectID)[i])
subset.v.good <- which(dat$Height..m.[subset.v]>=1&dat$Height..m.[subset.v]<=2)
sel <- order(org.img.select[subset.v])[1:t.images[i]]
keep.idx <- which(subset.v[sel]%!in%rmv.idx)

#### order of all next images in line according to PS81b image selection:
already.scored <- order(org.img.select[subset.v])[1:t.images[i]]
bad.dist <- which(PS96b.dat$Height..m.[subset.v]<1|PS96b.dat$Height..m.[subset.v]>2)
scored.and.bad <- bad.dist[which(bad.dist%in%already.scored)]
scored.or.bad <- union(already.scored, bad.dist)
sel.add <- subset.v[-scored.or.bad]
## sort the images into a spatially balance sequence
img.seq <- order(PS96b.dat$image.select[sel.add])
#img.seq <- order(PS96.dat$image.select[sel.add])
## check where the next images are positioned
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[1:20]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[1:20]]],adj=c(0.5,2))
## select which images are too close to already scored images
dontUse <- c(2,6,8,18,22,23,25)
idx.range <- 1:30
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[idx.range[dontUse]]]], col="yellow", pch=16)
text(polar.dat[sel.add[img.seq[idx.range[dontUse]]]],adj=c(0.5,-2), labels=idx.range[dontUse])
points(polar.dat[sel.add[img.seq[idx.range[-dontUse]]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[idx.range[-dontUse]]]],adj=c(0.5,2), labels=idx.range[-dontUse])

## plot next 8 images to pick only:
plot(polar.dat[subset.v[subset.v.good]])
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[c(1:21)[-dontUse]]]], col="deepskyblue", pch=16)
text(polar.dat[sel.add[img.seq[c(1:21)[-dontUse]]]],adj=c(0.5,2), labels=c(1:21)[-dontUse])

#### edit image selection column
## last image of the so far scored images
no.in.last.img <- tail(org.img.select[subset.v[sel]],1)
## keep already scored good-distance images at the top of the list
sel.clean <- subset.v[sel[sel%!in%scored.and.bad]]
dat$image.select[sel.clean] <- org.img.select[sel.clean]
## put the selected images next in line
dat$image.select[sel.add[img.seq[-dontUse]]] <- PS96b.dat$image.select[sel.add[img.seq[-dontUse]]]+no.in.last.img
## move dontUse images from PS81b to the back of the valid list (between the list and the NAs)
dat$image.select[sel.add[img.seq[dontUse]]] <- PS96b.dat$image.select[sel.add[img.seq[dontUse]]]+2000
## move bad.dist all the way to the back of the list
dat$image.select[subset.v[bad.dist]] <- 6666

#### check output by plotting all images in their sequence
sel.final <- order(dat$image.select[subset.v])[1:100]
plot(polar.dat[subset.v])
points(polar.dat[subset.v[sel.final]], col="blue", pch=15)
points(polar.dat[subset.v[sel][keep.idx]], col="red", pch=15)
points(polar.dat[sel.add[img.seq[c(1:100)[-dontUse]]]], col="deepskyblue", pch=16)


PS96.dat <- dat
PS96.dat$dist.from.start <- PS96b.dat$dist.from.start
save(PS96.dat,total.t.length.v, file="C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/PS96_dat_FINAL.Rdata")





































