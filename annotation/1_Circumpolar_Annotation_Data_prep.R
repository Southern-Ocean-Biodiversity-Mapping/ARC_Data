### 1) ### Setting up----
#library(raster)
library(readxl)
library(readr)
library(dplyr)
library(data.table)
library(proj4)
library(stringr)
library(terra)

user = "Jan"
#user = "charley"
#user="nicole"

if (user == "Jan") {
  
  sci.dir <-      "C:/Users/jjansen/OneDrive - University of Tasmania/science/"
  env.derived <-  paste0(sci.dir,"data_environmental/derived/")
  #bio.dir <-      paste0(sci.dir,"data_biological/")
  
  ## remote repository (DOESN'T WORK YET):
  # env.dir <- "https://data.imas.utas.edu.au/data_transfer/admin/files/EnvironmentalData/"
  
  ## common paths (after "sci.dir")
  tools.dir <-    paste0(sci.dir,"SouthernOceanBiodiversityMapping/Useful_Functions_Tools/")
  ARC_Data.dir <- paste0(sci.dir,"SouthernOceanBiodiversityMapping/ARC_Data/")
  
} 
if (user == "charley") {
  
  sci.dir <- "C:/Users/cgros/code/IMAS/"
  ARC_Data.dir <- paste0(sci.dir,"ARC_Data/")
  env.derived <-  "C:/Users/cgros/data/SO_env_layers/derived/"
  tools.dir <-    paste0(sci.dir,"Useful_Functions_Tools/")
  
}
if (user == "nicole") {
  
  sci.dir <-    "C:/Users/hillna/OneDrive - University of Tasmania/UTAS_work/Projects/Benthic Diversity ARC/"
  ARC_Data.dir <- paste0(sci.dir,"Analysis/ARC_Data/")
  env.derived <-  paste0(sci.dir,"data_environmental/derived/")
  tools.dir <-    paste0(sci.dir,"Analysis/Useful_Functions_Tools/")
  
}

'%!in%' <- function(x,y)!('%in%'(x,y))

###############
## choose resolution of environmental variables:
res <- "500m"
res <- "2km"
###############

## R-drive paths
RS.dir <- "R:/IMAS/Antarctic_Seafloor/Clean_Data_For_Permanent_Storage/"
ann.dir <- paste0(RS.dir,"AnnotationLibrary_AllFinishedSurveys/")

# Image quality path
image.quality.path <- "R:/IMAS/Antarctic_Seafloor/image_quality_analysis/image_quality_score.csv"

##### load still and diatom sample locations and bathymetry:

## from "Readin_Circumpolar_DownwardImage_Data.Rmd"
load(paste0(ARC_Data.dir,"prep_image/Circumpolar_DownwardImages_metadata_2024.Rdata"))

## from "ReadIn_Circumpolar_Environmental_Data.Rmd"
r2 <- rast(paste0(env.derived,"Circumpolar_EnvData_",res,"_shelf_bathy_ibcso2_depth.tif"))

##### load coastline
#stereo <- as.character(r2@crs)
#load(paste0(env.derived,"Circumpolar_Coastline.Rdata"))

source(paste0(tools.dir,"bubbleplot.R"))





### 2) Data Preparation

## load area estimates for each survey (?Why this comment: "no separate area-files for PS96 & PS18")
area_xls_files <- list.files(RS.dir,full.names=TRUE,pattern=".xlsx")
dat_area <- sapply(area_xls_files, readxl::read_excel, simplify=FALSE, skip=1) %>% bind_rows(.id="id")
names(dat_area)[8] <- "image area in m2" ## label for area contains superscript which is annoying
dat_area_subset <- dat_area[match(unique(dat_area$`image filename`),dat_area$`image filename`),]

## load Biigle image annotations
Biigle_csv_files <- list.files(paste0(ann.dir,"Counts_annotation_files_individual_surveys/"),full.names=TRUE,recursive=TRUE,pattern="254")
dat_counts <- sapply(Biigle_csv_files, readr::read_csv, simplify=FALSE) %>% bind_rows(.id="id")
dat_counts$Label <- gsub(" ","",dat_counts$label_hierarchy) # fixing label names
dat_counts$Label <- gsub(">","__",dat_counts$Label)
dat_counts$Label <- gsub("-","_",dat_counts$Label)

## load %-cover image annotations
dat_cover.raw <- read_csv(paste0(ann.dir,"Cover/Circumpolar_DownwardImages_PointScore_Annotations_20220307.csv"))
#remove uncorrected labels
dat_cover <- dat_cover.raw[,-c(4,8,11)]
names(dat_cover)[8] <- "Label"

## metadata for each image, including cellIDs (things like filename, full_path, surveyID, transectID, cellID, lon, lat, x, y, area, CoralNet, Biigle)
## dat.list.clean contains only the images on the shelf, but if we remove the seamount images here we can't test against them later!!!
image_metadata <- dplyr::select(dat.list[[1]],Filename.standardised,lon,lat,transectID)
for(i in 2:length(dat.list)){
  message(i)
  print(names(dat.list[[i]]))
  if("lon" %in% names(dat.list[[i]])){ dat.temp <- dplyr::select(dat.list[[i]],Filename.standardised,lon,lat,transectID) }
  if("Lon" %in% names(dat.list[[i]])){ dat.temp <- dplyr::select(dat.list[[i]],Filename.standardised,Lon,Lat,transectID) }
  if("Longitude" %in% names(dat.list[[i]])){ dat.temp <- dplyr::select(dat.list[[i]],Filename.standardised,Longitude,Latitude,transectID) }
  if("GPS_lon" %in% names(dat.list[[i]])){ dat.temp <- dplyr::select(dat.list[[i]],Filename.standardised,GPS_lon,GPS_lat,transectID) }  
  names(dat.temp) <- c("Filename.standardised","lon","lat","transectID")
  image_metadata <- rbind(image_metadata,dat.temp)
}
# image_metadata$Filename.standardised <- gsub(".tif",".jpg",image_metadata$Filename.standardised)
image_metadata$proj_coord_x <- proj4::project(image_metadata[,2:3], proj=crs(r2))$x
image_metadata$proj_coord_y <- proj4::project(image_metadata[,2:3], proj=crs(r2))$y
image_metadata$cellID <- extract(r2, image_metadata[,5:6], cellnumbers=TRUE)[,1]

name_split <- str_split(image_metadata$Filename.standardised,"_", simplify=T)
image_metadata$survey <- name_split[,1]
image_metadata$transectID_full <- paste0(name_split[,1],"_",name_split[,2])
image_metadata$transectID_full <- as.factor(image_metadata$transectID_full)
# image_metadata$survey[which(image_metadata$survey=="tan0802")] <- "TAN0802"
# image_metadata$survey[which(image_metadata$survey=="tan1901")] <- "TAN1901"

ids <- unique(image_metadata$cellID)

# ## checking if raster cells match up with transects
# r3 <- r2
# r3[setdiff(seq_len(ncell(r3)), ids)] <- NA
# r3[!is.na(r3)] <- 1
# cell_outlines <- rasterToPolygons(r3, dissolve=TRUE)
# 
# plot(r2,xlim=c(-2480000,-2460000),ylim=c(1260000,1280000))
# plot(cell_outlines, add=TRUE, border='red', lwd=2)
# points(image_metadata[which(image_metadata$survey=="CRS"),5:6])
# points(image_metadata[which(image_metadata$survey=="CRS"&image_metadata$cover=="yes"),5:6],col="blue")
# 
# plot(r2,xlim=c(-2475000,-2473000),ylim=c(1273000,1276000))
# plot(cell_outlines, add=TRUE, border='red', lwd=2)
# points(image_metadata[which(image_metadata$survey=="CRS"),5:6])
# points(image_metadata[which(image_metadata$survey=="CRS"&image_metadata$cover=="yes"),5:6],col="blue", pch=15)
# 
# plot(r2,xlim=c(-2478000,-2476000),ylim=c(1267500,1270500))
# plot(cell_outlines, add=TRUE, border='red', lwd=2)
# points(image_metadata[which(image_metadata$survey=="CRS"),5:6])
# points(image_metadata[which(image_metadata$survey=="CRS"&image_metadata$cover=="yes"),5:6],col="blue", pch=15)
# text(image_metadata[which(image_metadata$survey=="CRS"),5:6], labels=image_metadata$transectID[which(image_metadata$survey=="CRS")])

# Init area and area_source
image_metadata$area <- NA
image_metadata$area_source <- NA
# Pulling area info from metadata
for (survey_current in names(dat.list)) {
  if ("Area" %in% names(dat.list[[survey_current]])) {
    for (fname_current in dat.list[[survey_current]]$Filename.standardised) {
      area_current <- dat.list[[survey_current]][dat.list[[survey_current]]$Filename.standardised == fname_current, ]$Area
      if (!is.na(area_current)) {
        image_metadata[image_metadata$Filename.standardised == fname_current, "area"] <- area_current
        image_metadata[image_metadata$Filename.standardised == fname_current, "area_source"] <- "metadata"
      }
    }
  }
  if ("Area_from_metadata" %in% names(dat.list[[survey_current]])) {
    for (fname_current in dat.list[[survey_current]]$Filename.standardised) {
      area_current <- dat.list[[survey_current]][dat.list[[survey_current]]$Filename.standardised == fname_current, ]$Area_from_metadata
      if (!is.na(area_current)) {
        image_metadata[image_metadata$Filename.standardised == fname_current, "area"] <- area_current
        image_metadata[image_metadata$Filename.standardised == fname_current, "area_source"] <- "metadata"
      }
    }
  }
  if ("area" %in% names(dat.list[[survey_current]])) {
    for (fname_current in dat.list[[survey_current]]$Filename.standardised) {
      area_current <- dat.list[[survey_current]][dat.list[[survey_current]]$Filename.standardised == fname_current, ]$area
      if (!is.na(area_current)) {
        image_metadata[image_metadata$Filename.standardised == fname_current, "area"] <- area_current
        image_metadata[image_metadata$Filename.standardised == fname_current, "area_source"] <- "metadata"
      }
    }
  }
}
# Pulling area info from laser points
idx <- match(image_metadata$Filename.standardised,dat_area_subset$`image filename`)
fill.idx <- which(!is.na(idx))
search.idx <- idx[fill.idx]
image_metadata$area[fill.idx] <- dat_area_subset$`image area in m2`[search.idx] # replace area values where filenames match
image_metadata$area <- as.numeric(image_metadata$area)
image_metadata$area_source[fill.idx] <- "laser_points"
# When area data is not available, assign transect (or survey average?)
# Get indexes of images where area data not available
idx_area_not_available <- which(is.na(image_metadata$area))
image_metadata$area_source[idx_area_not_available] <- NA
# Iterate through surveys
for (survey_current in unique(image_metadata$survey)) {
  # Check if missing area values
  if (sum(is.na(image_metadata[image_metadata$survey == survey_current, ]$area))) {
    # Iterate through transects
    lst_transect <- unique(image_metadata[image_metadata$survey == survey_current, ]$transectID)
    for (transect_current in lst_transect) {
      area_transect_cur <- image_metadata[(image_metadata$survey == survey_current) 
                                          & (image_metadata$transectID == transect_current), ]$area
      if (sum(is.na(area_transect_cur)) & (sum(is.na(area_transect_cur)) != length(area_transect_cur))) {
        # Get indexes
        idx_NA_values <- which(is.na(image_metadata[(image_metadata$survey == survey_current) 
                                                    & (image_metadata$transectID == transect_current), ]$area))
        idx_nonNA_values <- which(!is.na(image_metadata[(image_metadata$survey == survey_current) 
                                                        & (image_metadata$transectID == transect_current), ]$area))
        # Compute average value of non NA images
        avg_value <- mean(image_metadata[(image_metadata$survey == survey_current) 
                                         & (image_metadata$transectID == transect_current), ]$area[idx_nonNA_values])
        ## FOR PS81,96,118 and AA2011 THIS SHOULD HAVE BEEN ONLY AVERAGES ACROSS METADATA IMAGES BECAUSE OF CROPPING (FIXED IN 5_SetupFinalARCDPDataset...)
        # Assign this value to images with NA area
        image_metadata[(image_metadata$survey == survey_current) 
                       & (image_metadata$transectID == transect_current), ]$area[idx_NA_values] <- avg_value
        image_metadata[(image_metadata$survey == survey_current) 
                       & (image_metadata$transectID == transect_current), ]$area_source[idx_NA_values] <- "averaged_across_transect"
      }
    }
    area_survey_cur <- image_metadata[(image_metadata$survey == survey_current), ]$area
    if (sum(is.na(area_survey_cur)) & (sum(is.na(area_survey_cur)) != length(area_survey_cur))) {
      # Get indexes
      idx_NA_values <- which(is.na(image_metadata[(image_metadata$survey == survey_current), ]$area))
      idx_nonNA_values <- which(!is.na(image_metadata[(image_metadata$survey == survey_current), ]$area))
      # Compute average value of non NA images
      avg_value <- mean(image_metadata[(image_metadata$survey == survey_current), ]$area[idx_nonNA_values])
      # Assign this value to images with NA area
      image_metadata[(image_metadata$survey == survey_current), ]$area[idx_NA_values] <- avg_value
      image_metadata[(image_metadata$survey == survey_current), ]$area_source[idx_NA_values] <- "averaged_across_survey"
    }
  }
}

#
# which images have been annotated by which platform?
image_metadata$cover <- "no"
image_metadata$counts <- "no"

# find matching filenames for cover (simple: all images in coralnet)
idx <- match(image_metadata$Filename.standardised,unique(dat_cover$Name))
fill.idx <- which(!is.na(idx))
search.idx <- idx[fill.idx]
image_metadata$cover[fill.idx] <- "yes" # replace area values where filenames match

## check if we selected all CoralNet images:
not.sel <- which(unique(dat_cover$Name)%!in%image_metadata$Filename.standardised)
unique(dat_cover$Name)[not.sel] ## some transects show up that have been disregarded but still uploaeded to coralnet

# find matching filenames for counts (not simple, images might have been annotated but no animals found, therefore no annotations...!)
# idx <- match(image_metadata$Filename.standardised,unique(dat_counts$filename))
# fill.idx <- which(!is.na(idx))
# search.idx <- idx[fill.idx]
# image_metadata$counts[fill.idx] <- "yes" # replace area values where filenames match
t.lvl <- levels(image_metadata$transectID_full)
for(i in 1:length(t.lvl)){
  print(i)
  sel <- which(image_metadata$transectID_full==t.lvl[i]) # select the transect
  sel.cov <- which(image_metadata$cover[sel]=="yes") # only images annotated in cover
  cts.img.sel <- order(image_metadata$Filename.standardised[sel[sel.cov]]) # order them by filename
  cts.length <- 1:ceiling(length(sel.cov)/2) # calculate how many images have been annotated in Biigle
  cts.files <- image_metadata$Filename.standardised[sel[sel.cov]][cts.img.sel][cts.length] # check filenames
  #print(cts.files)
  image_metadata$counts[sel[sel.cov]][cts.img.sel][cts.length] <- "yes"
}
## check if we selected all Biigle images:
not.sel <- which(unique(dat_counts$filename)%!in%image_metadata$Filename.standardised)
unique(dat_counts$filename)[not.sel]  ## some transects show up that have been disregarded but still uploaeded to coralnet

# Add image quality score
# Read csv file
df_image_quality_score <- read.csv(image.quality.path)
# Add column to image_metadata by matching filename columns
image_metadata$image_quality_score <- df_image_quality_score$image_quality_score[match(image_metadata$Filename.standardised, df_image_quality_score$filename)]

### add survey/year/gear information based on filename
image_metadata$year <- NA
image_metadata$year[image_metadata$survey=="PS06"] <- 1984
image_metadata$year[image_metadata$survey=="PS14"] <- 1989
image_metadata$year[image_metadata$survey=="PS18"] <- 1990
image_metadata$year[image_metadata$survey=="PS61"] <- 2002
image_metadata$year[image_metadata$survey=="PS81"] <- 2013
image_metadata$year[image_metadata$survey=="PS96"] <- 2015
image_metadata$year[image_metadata$survey=="PS118"] <- 2019
image_metadata$year[image_metadata$survey=="tan0802"] <- 2008
image_metadata$year[image_metadata$survey=="TAN1802"] <- 2018
image_metadata$year[image_metadata$survey=="tan1901"] <- 2019
image_metadata$year[image_metadata$survey=="AA2011"] <- 2011
image_metadata$year[image_metadata$survey=="CRS"] <- NA
image_metadata$year[image_metadata$survey=="NBP1402"] <- 2014
image_metadata$year[image_metadata$survey=="NBP1502"] <- 2015
image_metadata$year[image_metadata$survey=="LMG1311"] <- 2013
image_metadata$year[image_metadata$survey=="JR262"] <- 2011
image_metadata$year[image_metadata$survey=="JR15005"] <- 2015
image_metadata$year[image_metadata$survey=="JR17001"] <- 2017
image_metadata$year[image_metadata$survey=="JR17003"] <- 2018
image_metadata$gear <- NA
image_metadata$gear[image_metadata$survey=="PS06"] <- "FTS"
image_metadata$gear[image_metadata$survey=="PS14"] <- "FTS"
image_metadata$gear[image_metadata$survey=="PS18"] <- "FTS"
image_metadata$gear[image_metadata$survey=="PS61"] <- "FTS"
image_metadata$gear[image_metadata$survey=="PS81"] <- "OFOS"
image_metadata$gear[image_metadata$survey=="PS96"] <- "OFOS"
image_metadata$gear[image_metadata$survey=="PS118"] <- "OFOBS"
image_metadata$gear[image_metadata$survey=="tan0802"] <- "DTIS"
image_metadata$gear[image_metadata$survey=="TAN1802"] <- "DTIS"
image_metadata$gear[image_metadata$survey=="tan1901"] <- "DTIS"
image_metadata$gear[image_metadata$survey=="AA2011"] <- "CTD"
image_metadata$gear[image_metadata$survey=="CRS"] <- "YOYO"
image_metadata$gear[image_metadata$survey=="NBP1402"] <- "YOYO"
image_metadata$gear[image_metadata$survey=="NBP1502"] <- "YOYO"
image_metadata$gear[image_metadata$survey=="LMG1311"] <- "YOYO"
image_metadata$gear[image_metadata$survey=="JR262"] <- "SUCS"
image_metadata$gear[image_metadata$survey=="JR15005"] <- "SUCS"
image_metadata$gear[image_metadata$survey=="JR17001"] <- "SUCS"
image_metadata$gear[image_metadata$survey=="JR17003"] <- "SUCS"

# ## CRS surveys are from multiple years
dates <- unique(dat.list$WAP$Date)
dates
# sel.dates <- which(dat.list$WAP$Date==dates[14])
# dat.list$WAP$transectID[sel.dates]
  
image_metadata$year[image_metadata$survey=="CRS" & image_metadata$transectID=="1069"] <- 2008
image_metadata$year[image_metadata$survey=="CRS" & image_metadata$transectID=="1072"] <- 2008
image_metadata$year[image_metadata$survey=="CRS" & image_metadata$transectID=="1091"] <- 2008
image_metadata$year[image_metadata$survey=="CRS" & image_metadata$transectID=="1103"] <- 2008
image_metadata$year[image_metadata$survey=="CRS" & image_metadata$transectID=="1130"] <- 2008
image_metadata$year[image_metadata$survey=="CRS" & image_metadata$transectID=="1132"] <- 2008
image_metadata$year[image_metadata$survey=="CRS" & image_metadata$transectID=="1207"] <- 2009
image_metadata$year[image_metadata$survey=="CRS" & image_metadata$transectID=="1208"] <- 2009
image_metadata$year[image_metadata$survey=="CRS" & image_metadata$transectID=="1217"] <- 2009
image_metadata$year[image_metadata$survey=="CRS" & image_metadata$transectID=="1219"] <- 2009
image_metadata$year[image_metadata$survey=="CRS" & image_metadata$transectID=="1255"] <- 2009
image_metadata$year[image_metadata$survey=="CRS" & image_metadata$transectID=="1267"] <- 2009
image_metadata$year[image_metadata$survey=="CRS" & image_metadata$transectID=="1276"] <- 2010
image_metadata$year[image_metadata$survey=="CRS" & image_metadata$transectID=="1278"] <- 2010
image_metadata$year[image_metadata$survey=="CRS" & image_metadata$transectID=="1279"] <- 2010
image_metadata$year[image_metadata$survey=="CRS" & image_metadata$transectID=="1280"] <- 2010
image_metadata$year[image_metadata$survey=="CRS" & image_metadata$transectID=="1281"] <- 2010
image_metadata$year[image_metadata$survey=="CRS" & image_metadata$transectID=="1282"] <- 2010
image_metadata$year[image_metadata$survey=="CRS" & image_metadata$transectID=="1283"] <- 2010
image_metadata$year[image_metadata$survey=="CRS" & image_metadata$transectID=="1284"] <- 2010
image_metadata$year[image_metadata$survey=="CRS" & image_metadata$transectID=="1285"] <- 2010
image_metadata$year[image_metadata$survey=="CRS" & image_metadata$transectID=="1286"] <- 2010
image_metadata$year[image_metadata$survey=="CRS" & image_metadata$transectID=="1289"] <- 2010
image_metadata$year[image_metadata$survey=="CRS" & image_metadata$transectID=="1290"] <- 2010
image_metadata$year[image_metadata$survey=="CRS" & image_metadata$transectID=="1295"] <- 2010
image_metadata$year[image_metadata$survey=="CRS" & image_metadata$transectID=="1297"] <- 2010
image_metadata$year[image_metadata$survey=="CRS" & image_metadata$transectID=="1300"] <- 2010
image_metadata$year[image_metadata$survey=="CRS" & image_metadata$transectID=="1315"] <- 2010
image_metadata$year[image_metadata$survey=="CRS" & image_metadata$transectID=="1320"] <- 2010
image_metadata$year[image_metadata$survey=="CRS" & image_metadata$transectID=="1323"] <- 2010
image_metadata$year[image_metadata$survey=="CRS" & image_metadata$transectID=="1324"] <- 2010
image_metadata$year[image_metadata$survey=="CRS" & image_metadata$transectID=="1325"] <- 2010
image_metadata$year[image_metadata$survey=="CRS" & image_metadata$transectID=="1337"] <- 2010
image_metadata$year[image_metadata$survey=="CRS" & image_metadata$transectID=="1338"] <- 2010


## add transect location names where available
image_metadata$transect_location <- NA
image_metadata$transect_location[image_metadata$survey=="CRS" & image_metadata$transectID==1130] <- "StationB2"
image_metadata$transect_location[image_metadata$survey=="CRS" & image_metadata$transectID==1132] <- "StationB2"
image_metadata$transect_location[image_metadata$survey=="CRS" & image_metadata$transectID==1255] <- "StationB3"
image_metadata$transect_location[image_metadata$survey=="CRS" & image_metadata$transectID==1267] <- "StationB3"
image_metadata$transect_location[image_metadata$survey=="CRS" & image_metadata$transectID==1091] <- "StationE2"
image_metadata$transect_location[image_metadata$survey=="CRS" & image_metadata$transectID==1103] <- "StationE2"
image_metadata$transect_location[image_metadata$survey=="CRS" & image_metadata$transectID==1217] <- "StationE3"
image_metadata$transect_location[image_metadata$survey=="CRS" & image_metadata$transectID==1219] <- "StationE3"
image_metadata$transect_location[image_metadata$survey=="CRS" & image_metadata$transectID==1069] <- "StationF2"
image_metadata$transect_location[image_metadata$survey=="CRS" & image_metadata$transectID==1072] <- "StationF2"
image_metadata$transect_location[image_metadata$survey=="CRS" & image_metadata$transectID==1207] <- "StationF3"
image_metadata$transect_location[image_metadata$survey=="CRS" & image_metadata$transectID==1208] <- "StationF3"
image_metadata$transect_location[image_metadata$survey=="CRS" & image_metadata$transectID==1315] <- "Lockyer"
image_metadata$transect_location[image_metadata$survey=="CRS" & image_metadata$transectID==1320] <- "Hughes"
image_metadata$transect_location[image_metadata$survey=="CRS" & image_metadata$transectID==1323] <- "Hughes"
image_metadata$transect_location[image_metadata$survey=="CRS" & image_metadata$transectID==1324] <- "Hughes"
image_metadata$transect_location[image_metadata$survey=="CRS" & image_metadata$transectID==1325] <- "Hughes"
image_metadata$transect_location[image_metadata$survey=="CRS" & image_metadata$transectID==1285] <- "Andvord_Inner"
image_metadata$transect_location[image_metadata$survey=="CRS" & image_metadata$transectID==1286] <- "Andvord_Inner"
image_metadata$transect_location[image_metadata$survey=="CRS" & image_metadata$transectID==1337] <- "Andvord_Middle"
image_metadata$transect_location[image_metadata$survey=="CRS" & image_metadata$transectID==1338] <- "Andvord_Middle"
image_metadata$transect_location[image_metadata$survey=="CRS" & image_metadata$transectID==1289] <- "Andvord_Mouth"
image_metadata$transect_location[image_metadata$survey=="CRS" & image_metadata$transectID==1290] <- "Andvord_Mouth"
image_metadata$transect_location[image_metadata$survey=="CRS" & image_metadata$transectID==1283] <- "Andvord_Outer"
image_metadata$transect_location[image_metadata$survey=="CRS" & image_metadata$transectID==1284] <- "Andvord_Outer"
image_metadata$transect_location[image_metadata$survey=="CRS" & image_metadata$transectID==1295] <- "Barilari_Inner"
image_metadata$transect_location[image_metadata$survey=="CRS" & image_metadata$transectID==1297] <- "Barilari_Inner"
image_metadata$transect_location[image_metadata$survey=="CRS" & image_metadata$transectID==1300] <- "Barilari_Outer"
image_metadata$transect_location[image_metadata$survey=="CRS" & image_metadata$transectID==1279] <- "Flandres_InnerA"
image_metadata$transect_location[image_metadata$survey=="CRS" & image_metadata$transectID==1280] <- "Flandres_InnerA"
image_metadata$transect_location[image_metadata$survey=="CRS" & image_metadata$transectID==1276] <- "Flandres_InnerB"
image_metadata$transect_location[image_metadata$survey=="CRS" & image_metadata$transectID==1278] <- "Flandres_InnerB"
image_metadata$transect_location[image_metadata$survey=="CRS" & image_metadata$transectID==1281] <- "Flandres_Mouth"
image_metadata$transect_location[image_metadata$survey=="CRS" & image_metadata$transectID==1282] <- "Flandres_Mouth"

image_metadata$transect_location[image_metadata$survey=="JR17001" & image_metadata$transectID=="AI"] <- "AnversIsland"
image_metadata$transect_location[image_metadata$survey=="JR17001" & image_metadata$transectID=="KGI"] <- "KingGeorgeIsland"
image_metadata$transect_location[image_metadata$survey=="JR17001" & image_metadata$transectID=="MT1"] <- "MargueriteBayTrough1"
image_metadata$transect_location[image_metadata$survey=="JR17001" & image_metadata$transectID=="MT2"] <- "MargueriteBayTrough2"
image_metadata$transect_location[image_metadata$survey=="JR17001" & image_metadata$transectID=="MT3"] <- "MargueriteBayTrough3"
image_metadata$transect_location[image_metadata$survey=="JR17001" & image_metadata$transectID=="MT4"] <- "MargueriteBayTrough4"
image_metadata$transect_location[image_metadata$survey=="JR17001" & image_metadata$transectID=="MT5"] <- "MargueriteBayTrough5"

image_metadata$transect_location[image_metadata$survey=="LMG1311" & image_metadata$transectID=="2"] <- "GerlacheStrait"
image_metadata$transect_location[image_metadata$survey=="LMG1311" & image_metadata$transectID=="3"] <- "GerlacheStrait"

## change survey names for CRS cruises to the actual cruise-names where available
NBP0808_transects <- c("1069","1072","1091","1103","1130","1132")
NBP1001_transects <- c("1276","1278","1279","1280","1281","1282","1283","1284","1285","1286","1289","1290","1295","1297","1300","1315","1320","1323","1324","1325","1337","1338")
image_metadata$survey[image_metadata$survey=="CRS" & image_metadata$transectID%in%NBP0808_transects] <- "NBP0808"
image_metadata$survey[image_metadata$survey=="CRS" & image_metadata$transectID%in%NBP1001_transects] <- "NBP1001"

# image_metadata$year[image_metadata$survey=="CRS"] <- NA

for(i in c(1,8,15,17)){
  image_metadata[,i] <- as.factor(image_metadata[,i])
}


### 3) Generate Counts data  ----

### 3a) site(image)-by-species matrix

labs.counts <- unique(dat_counts$Label)
nams.counts <- unique(dat_counts$filename)
dat_counts_image_by_species <- data.frame(matrix(NA,nrow=length(nams.counts),ncol=length(labs.counts)))
rownames(dat_counts_image_by_species) <- nams.counts
colnames(dat_counts_image_by_species) <- labs.counts
for(i in 1:length(nams.counts)){
  sel.r <- which(dat_counts$filename==nams.counts[i])  
  for(j in 1:length(labs.counts)){
    sel.c <- which(dat_counts$Label[sel.r]==labs.counts[j])  
    dat_counts_image_by_species[i,j] <- length(sel.c)
  }
}

### 3b) site(cell)-by-species matrix

## create new dataframe using image locations
dat_counts_cell_by_species <- data.frame(matrix(NA,nrow=length(ids),ncol=length(labs.counts)))
rownames(dat_counts_cell_by_species) <- ids
colnames(dat_counts_cell_by_species) <- labs.counts
counts_N <- rep(NA, length(ids))
counts_area <- rep(NA, length(ids))
counts_cells_survey1 <- rep(NA, length(ids))
counts_cells_survey2 <- rep(NA, length(ids))
counts_cells_transect1 <- as.character(rep(NA, length(ids)))
counts_cells_transect2 <- rep(NA, length(ids))
counts_cells_transect3 <- rep(NA, length(ids))
counts_cells_transect4 <- rep(NA, length(ids))
counts_cells_mean_image_quality_score <- rep(NA, length(ids))

for(i in 1:length(ids)){
  sel.r <- which(image_metadata$cellID==ids[i]&image_metadata$counts=="yes")  # find images that are part of that cell and annotated in biigle
  if(length(sel.r)==0){ # if none of the images are annotated, skip to the next iteration
    #print("0")
    next
  }
  sel.names <- image_metadata$Filename.standardised[sel.r] # find the names of these images
  dat.temp <- dat_counts_image_by_species[which(nams.counts%in%sel.names),] # find annotations from these images from the image-dataset
  counts_N[i] <-  nrow(dat.temp) # how many images are there
  counts_area[i] <- sum(image_metadata$area[sel.r]) # total area across these images
  counts_cells_mean_image_quality_score[i] <- mean(image_metadata$image_quality_score[sel.r]) # ave quality across these images
  dat_counts_cell_by_species[i,] <- colSums(dat.temp) # how many individuals per species
  ## check if all images are from the same survey
  counts_cells_survey1[i] <- as.character(unique(image_metadata$survey[sel.r])[1])
  print(counts_cells_survey1[i])
  if(length(unique(image_metadata$survey[sel.r]))>1){
    message("TWO surveys!")
    counts_cells_survey2[i] <- unique(image_metadata$survey[sel.r])[2]
  }
  ## ... and same transect
  counts_cells_transect1[i] <- unique(as.character(image_metadata$transectID_full)[sel.r])[1]
  message(counts_cells_transect1[i])
  # print(length(unique(as.character(image_metadata$transectID_full)[sel.r])))
  if(length(unique(image_metadata$transectID_full[sel.r]))>1){
    counts_cells_transect2[i] <- unique(as.character(image_metadata$transectID_full)[sel.r])[2]
    if(length(unique(image_metadata$transectID_full[sel.r]))>2){
      counts_cells_transect3[i] <- unique(as.character(image_metadata$transectID_full)[sel.r])[3]
      if(length(unique(image_metadata$transectID_full[sel.r]))>3){
        message("FOUR transects!")
        counts_cells_transect4[i] <- unique(as.character(image_metadata$transectID_full)[sel.r])[4]
      }
    }}
}

### 4) %-cover ----

### 4a) site(image)-by-species matrix

labs.cov <- unique(dat_cover$Label)
nams.cov <- unique(dat_cover$Name)
dat_cover_image_by_species <- data.frame(matrix(NA,nrow=length(nams.cov),ncol=length(labs.cov)))
rownames(dat_cover_image_by_species) <- nams.cov
colnames(dat_cover_image_by_species) <- labs.cov
for(i in 1:length(nams.cov)){
  sel.r <- which(dat_cover$Name==nams.cov[i])  
  for(j in 1:length(labs.cov)){
    sel.c <- which(dat_cover$Label[sel.r]==labs.cov[j])  
    dat_cover_image_by_species[i,j] <- length(sel.c)
  }
}

## sites on the seamounts are still incuded here
sel.seamounts <- which(rownames(dat_cover_image_by_species)%!in%image_metadata$Filename.standardised[image_metadata$cover=="yes"])
rownames(dat_cover_image_by_species)[sel.seamounts]


### 4b) site(cell)-by-species matrix

## create new dataframe using image locations
dat_cover_cell_by_species <- data.frame(matrix(NA,nrow=length(ids),ncol=length(labs.cov)))
rownames(dat_cover_cell_by_species) <- ids
colnames(dat_cover_cell_by_species) <- labs.cov
cover_N <- rep(NA, length(ids))
cover_area <- rep(NA, length(ids))
cover_cells_survey1 <- rep(NA, length(ids))
cover_cells_survey2 <- rep(NA, length(ids))
cover_cells_transect1 <- rep(NA, length(ids))
cover_cells_transect2 <- rep(NA, length(ids))
cover_cells_transect3 <- rep(NA, length(ids))
cover_cells_transect4 <- rep(NA, length(ids))
cover_cells_mean_image_quality_score <- rep(NA, length(ids))
for(i in 1:length(ids)){
  #print(i)
  sel.r <- which(image_metadata$cellID==ids[i]&image_metadata$cover=="yes") # find images that are part of that cell and annotated in coralnet
  if(length(sel.r)==0) next # if none of the images are annotated, skip to the next iteration
  sel.names <- image_metadata$Filename.standardised[sel.r] # find the names of these images
  dat.temp <- dat_cover_image_by_species[which(nams.cov%in%sel.names),] # find annotations from these images from the image-dataset
  cover_N[i] <-  nrow(dat.temp)
  cover_area[i] <- sum(image_metadata$area[sel.r])
  cover_cells_mean_image_quality_score[i] <- mean(image_metadata$image_quality_score[sel.r]) # ave quality across these images
  dat_cover_cell_by_species[i,] <- colSums(dat.temp)
  ## check if all images are from the same survey
  cover_cells_survey1[i] <- as.character(unique(image_metadata$survey[sel.r])[1])
  if(length(unique(image_metadata$survey[sel.r]))>1){
    message("Two surveys!")
    cover_cells_survey2[i] <- unique(image_metadata$survey[sel.r])[2]
  }
  ## ... and same transect
  cell.transects <- unique(as.character(image_metadata$transectID_full)[sel.r])
  cover_cells_transect1[i] <- cell.transects[1]
  t.per.cell <- length(cell.transects)
  if(t.per.cell>1){
    cover_cells_transect2[i] <- cell.transects[2]
    if(t.per.cell>2){
      message("Three transects!")
      message(t.per.cell)
      cover_cells_transect3[i] <- cell.transects[3]
      if(t.per.cell>3){
        message("FOUR transects!")
        message(t.per.cell)
        cover_cells_transect4[i] <- cell.transects[4]
      }
    }}
}

### 4b) cell metadata

cell.coords <- xyFromCell(r2, ids)
cell.lonlat <- project(cell.coords, proj=crs(r2), inverse=TRUE) 

cell_metadata <- data.frame(cbind(ids,cell.lonlat,cell.coords,cover_N, counts_N, cover_area, counts_area, 
                                  cover_cells_survey1, cover_cells_transect1, cover_cells_transect2, cover_cells_transect3, cover_cells_transect4, 
                                  counts_cells_survey1, counts_cells_transect1, counts_cells_transect2, counts_cells_transect3,counts_cells_transect4,cover_cells_mean_image_quality_score))
names(cell_metadata) <- c("cellID", "lon", "lat", "proj_coord_x", "proj_coord_y", "cover_N", "counts_N", "cover_area", "counts_area",
                          "cover_cells_survey", "cover_cells_transect1", "cover_cells_transect2", "cover_cells_transect3","cover_cells_transect4", 
                          "counts_cells_survey", "counts_cells_transect1", "counts_cells_transect2", "counts_cells_transect3","counts_cells_transect4","image_quality_score")

### add survey/year/gear information to cell metadata
cell_metadata$year <- NA
cell_metadata$year[cell_metadata$cover_cells_survey=="PS06"] <- 1984
cell_metadata$year[cell_metadata$cover_cells_survey=="PS14"] <- 1989
cell_metadata$year[cell_metadata$cover_cells_survey=="PS18"] <- 1990
cell_metadata$year[cell_metadata$cover_cells_survey=="PS61"] <- 2002
cell_metadata$year[cell_metadata$cover_cells_survey=="PS81"] <- 2013
cell_metadata$year[cell_metadata$cover_cells_survey=="PS96"] <- 2015
cell_metadata$year[cell_metadata$cover_cells_survey=="PS118"] <- 2019
cell_metadata$year[cell_metadata$cover_cells_survey=="tan0802"] <- 2008
cell_metadata$year[cell_metadata$cover_cells_survey=="TAN1802"] <- 2018
cell_metadata$year[cell_metadata$cover_cells_survey=="tan1901"] <- 2019
cell_metadata$year[cell_metadata$cover_cells_survey=="AA2011"] <- 2011
cell_metadata$year[cell_metadata$cover_cells_survey=="NBP0808"] <- 2008
cell_metadata$year[cell_metadata$cover_cells_survey=="CRS"] <- 2009
cell_metadata$year[cell_metadata$cover_cells_survey=="NBP1001"] <- 2010
cell_metadata$year[cell_metadata$cover_cells_survey=="NBP1402"] <- 2014
cell_metadata$year[cell_metadata$cover_cells_survey=="NBP1502"] <- 2015
cell_metadata$year[cell_metadata$cover_cells_survey=="LMG1311"] <- 2013
cell_metadata$year[cell_metadata$cover_cells_survey=="JR262"] <- 2011
cell_metadata$year[cell_metadata$cover_cells_survey=="JR15005"] <- 2015
cell_metadata$year[cell_metadata$cover_cells_survey=="JR17001"] <- 2017
cell_metadata$year[cell_metadata$cover_cells_survey=="JR17003"] <- 2018
cell_metadata$gear <- NA
cell_metadata$gear[cell_metadata$cover_cells_survey=="PS06"] <- "FTS"
cell_metadata$gear[cell_metadata$cover_cells_survey=="PS14"] <- "FTS"
cell_metadata$gear[cell_metadata$cover_cells_survey=="PS18"] <- "FTS"
cell_metadata$gear[cell_metadata$cover_cells_survey=="PS61"] <- "FTS"
cell_metadata$gear[cell_metadata$cover_cells_survey=="PS81"] <- "OFOS"
cell_metadata$gear[cell_metadata$cover_cells_survey=="PS96"] <- "OFOS"
cell_metadata$gear[cell_metadata$cover_cells_survey=="PS118"] <- "OFOBS"
cell_metadata$gear[cell_metadata$cover_cells_survey=="tan0802"] <- "DTIS"
cell_metadata$gear[cell_metadata$cover_cells_survey=="TAN1802"] <- "DTIS"
cell_metadata$gear[cell_metadata$cover_cells_survey=="tan1901"] <- "DTIS"
cell_metadata$gear[cell_metadata$cover_cells_survey=="AA2011"] <- "CTD"
cell_metadata$gear[cell_metadata$cover_cells_survey=="NBP0808"] <- "YOYO"
cell_metadata$gear[cell_metadata$cover_cells_survey=="CRS"] <- "YOYO"
cell_metadata$gear[cell_metadata$cover_cells_survey=="NBP1001"] <- "YOYO"
cell_metadata$gear[cell_metadata$cover_cells_survey=="NBP1402"] <- "YOYO"
cell_metadata$gear[cell_metadata$cover_cells_survey=="NBP1502"] <- "YOYO"
cell_metadata$gear[cell_metadata$cover_cells_survey=="LMG1311"] <- "YOYO"
cell_metadata$gear[cell_metadata$cover_cells_survey=="JR262"] <- "SUCS"
cell_metadata$gear[cell_metadata$cover_cells_survey=="JR15005"] <- "SUCS"
cell_metadata$gear[cell_metadata$cover_cells_survey=="JR17001"] <- "SUCS"
cell_metadata$gear[cell_metadata$cover_cells_survey=="JR17003"] <- "SUCS"

## add transect location names where available
cell_metadata$transect_location <- ""
cell_metadata$transect_location[cell_metadata$cover_cells_transect1=="CRS_1130"] <- "StationB2"
cell_metadata$transect_location[cell_metadata$cover_cells_transect1=="CRS_1132"] <- "StationB2"
cell_metadata$transect_location[cell_metadata$cover_cells_transect1=="CRS_1255"] <- "StationB3"
cell_metadata$transect_location[cell_metadata$cover_cells_transect1=="CRS_1267"] <- "StationB3"
cell_metadata$transect_location[cell_metadata$cover_cells_transect1=="CRS_1091"] <- "StationE2"
cell_metadata$transect_location[cell_metadata$cover_cells_transect1=="CRS_1103"] <- "StationE2"
cell_metadata$transect_location[cell_metadata$cover_cells_transect1=="CRS_1217"] <- "StationE3"
cell_metadata$transect_location[cell_metadata$cover_cells_transect1=="CRS_1219"] <- "StationE3"
cell_metadata$transect_location[cell_metadata$cover_cells_transect1=="CRS_1069"] <- "StationF2"
cell_metadata$transect_location[cell_metadata$cover_cells_transect1=="CRS_1072"] <- "StationF2"
cell_metadata$transect_location[cell_metadata$cover_cells_transect1=="CRS_1207"] <- "StationF3"
cell_metadata$transect_location[cell_metadata$cover_cells_transect1=="CRS_1208"] <- "StationF3"
cell_metadata$transect_location[cell_metadata$cover_cells_transect1=="CRS_1315"] <- "Lockyer"
cell_metadata$transect_location[cell_metadata$cover_cells_transect1=="CRS_1320"] <- "Hughes"
cell_metadata$transect_location[cell_metadata$cover_cells_transect1=="CRS_1323"] <- "Hughes"
cell_metadata$transect_location[cell_metadata$cover_cells_transect1=="CRS_1324"] <- "Hughes"
cell_metadata$transect_location[cell_metadata$cover_cells_transect1=="CRS_1325"] <- "Hughes"
cell_metadata$transect_location[cell_metadata$cover_cells_transect1=="CRS_1285"] <- "Andvord_Inner"
cell_metadata$transect_location[cell_metadata$cover_cells_transect1=="CRS_1286"] <- "Andvord_Inner"
cell_metadata$transect_location[cell_metadata$cover_cells_transect1=="CRS_1337"] <- "Andvord_Middle"
cell_metadata$transect_location[cell_metadata$cover_cells_transect1=="CRS_1338"] <- "Andvord_Middle"
cell_metadata$transect_location[cell_metadata$cover_cells_transect1=="CRS_1289"] <- "Andvord_Mouth"
cell_metadata$transect_location[cell_metadata$cover_cells_transect1=="CRS_1290"] <- "Andvord_Mouth"
cell_metadata$transect_location[cell_metadata$cover_cells_transect1=="CRS_1283"] <- "Andvord_Outer"
cell_metadata$transect_location[cell_metadata$cover_cells_transect1=="CRS_1284"] <- "Andvord_Outer"
cell_metadata$transect_location[cell_metadata$cover_cells_transect1=="CRS_1295"] <- "Barilari_Inner"
cell_metadata$transect_location[cell_metadata$cover_cells_transect1=="CRS_1297"] <- "Barilari_Inner"
cell_metadata$transect_location[cell_metadata$cover_cells_transect1=="CRS_1300"] <- "Barilari_Outer"
cell_metadata$transect_location[cell_metadata$cover_cells_transect1=="CRS_1279"] <- "Flandres_InnerA"
cell_metadata$transect_location[cell_metadata$cover_cells_transect1=="CRS_1280"] <- "Flandres_InnerA"
cell_metadata$transect_location[cell_metadata$cover_cells_transect1=="CRS_1276"] <- "Flandres_InnerB"
cell_metadata$transect_location[cell_metadata$cover_cells_transect1=="CRS_1278"] <- "Flandres_InnerB"
cell_metadata$transect_location[cell_metadata$cover_cells_transect1=="CRS_1281"] <- "Flandres_Mouth"
cell_metadata$transect_location[cell_metadata$cover_cells_transect1=="CRS_1282"] <- "Flandres_Mouth"

cell_metadata$transect_location[cell_metadata$cover_cells_survey=="JR17001" & cell_metadata$cover_cells_transect1=="AI"] <- "AnversIsland"
cell_metadata$transect_location[cell_metadata$cover_cells_survey=="JR17001" & cell_metadata$cover_cells_transect1=="KGI"] <- "KingGeorgeIsland"
cell_metadata$transect_location[cell_metadata$cover_cells_survey=="JR17001" & cell_metadata$cover_cells_transect1=="MT1"] <- "MargueriteBayTrough1"
cell_metadata$transect_location[cell_metadata$cover_cells_survey=="JR17001" & cell_metadata$cover_cells_transect1=="MT2"] <- "MargueriteBayTrough2"
cell_metadata$transect_location[cell_metadata$cover_cells_survey=="JR17001" & cell_metadata$cover_cells_transect1=="MT3"] <- "MargueriteBayTrough3"
cell_metadata$transect_location[cell_metadata$cover_cells_survey=="JR17001" & cell_metadata$cover_cells_transect1=="MT4"] <- "MargueriteBayTrough4"
cell_metadata$transect_location[cell_metadata$cover_cells_survey=="JR17001" & cell_metadata$cover_cells_transect1=="MT5"] <- "MargueriteBayTrough5"

cell_metadata$transect_location[cell_metadata$cover_cells_survey=="LMG1311" & cell_metadata$cover_cells_transect1=="3"] <- "GerlacheStrait"
cell_metadata$transect_location[cell_metadata$cover_cells_survey=="LMG1311" & cell_metadata$cover_cells_transect1=="2"] <- "GerlacheStrait"

## change character data to factors
for(i in c(1,10:19,22,23)){
  cell_metadata[,i] <- as.factor(cell_metadata[,i])
}
for(i in c(1,8,15)){
  image_metadata[,i] <- as.factor(image_metadata[,i])
}
## change character data to numeric
for(i in c(2:9,20)){
  cell_metadata[,i] <- as.numeric(cell_metadata[,i])
}

### 5) save output----
cover_cells <- dat_cover_cell_by_species[-which(is.na(rowSums(dat_cover_cell_by_species))),]
count_cells <- dat_counts_cell_by_species[-which(is.na(rowSums(dat_counts_cell_by_species))),]
cover_images <- dat_cover_image_by_species
count_images <- dat_counts_image_by_species

head(image_metadata)
head(cell_metadata)

save(image_metadata, cell_metadata, cover_cells, cover_images, count_cells, count_images,
     file=paste0(ARC_Data.dir,"annotation/Circumpolar_Annotation_Data_",res,"_202312.Rdata"))


# ### 6) Inspect resulting dataframes ----
# 
# ### 6a) locations across a region
# 
# ## check if things are correct:
# plot(r2, xlim=c(0,300000),ylim=c(-2100000,-1800000))
# points(image_metadata[,5:6])
# 
# 
# ### 6b) a single transect
# 
# plot(r2, xlim=c(282500,283500),ylim=c(-2011500,-2010000))
# points(image_metadata[,5:6])
# points(image_metadata[which(image_metadata$cover=="yes"),5:6],col="red",cex=2)
# points(image_metadata[which(image_metadata$counts=="yes"),5:6],col="blue",cex=3)
# text(image_metadata[which(image_metadata$cover=="yes"),5:6], labels=image_metadata[which(image_metadata$cover=="yes"),1], adj=-0.2)
# 
# ### 6c) transects that were faulty before: (e.g. CRS 1207, 1217)
# sel <- which(image_metadata$transectID_full=="CRS_1207")
# sel.cov <- which(image_metadata$cover[sel]=="yes")
# sel.cts <- which(image_metadata$counts[sel]=="yes")
# plot(r2, xlim=c(-2378000,-2374000), ylim=c(876000,880000))
# points(image_metadata[sel,5:6])
# points(image_metadata[sel[sel.cov],5:6],col="red",cex=2)
# points(image_metadata[sel[sel.cts],5:6],col="blue",cex=3)
# #text(image_metadata[sel[sel.cov],5:6], labels=image_metadata[sel[sel.cov],1], adj=-0.2)
# 
# sel <- which(image_metadata$transectID_full=="CRS_1217")
# sel.cov <- which(image_metadata$cover[sel]=="yes")
# sel.cts <- which(image_metadata$counts[sel]=="yes")
# plot(r2, xlim=c(-2443000,-2441000), ylim=c(1021000,1023000))
# points(image_metadata[sel,5:6])
# points(image_metadata[sel,5:6])
# points(image_metadata[sel[sel.cov],5:6],col="red",cex=2)
# points(image_metadata[sel[sel.cts],5:6],col="blue",cex=3)
# ## add other transects
# points(image_metadata[,5:6])
# points(image_metadata[which(image_metadata$cover=="yes"),5:6],col="red",cex=2)
# points(image_metadata[which(image_metadata$counts=="yes"),5:6],col="blue",cex=3)
# #text(image_metadata[which(image_metadata$cover=="yes"),5:6], labels=image_metadata[which(image_metadata$cover=="yes"),1], adj=-0.2)
# 
# sel <- which(image_metadata$transectID_full=="LMG1311_3")
# sel.cov <- which(image_metadata$cover[sel]=="yes")
# sel.cts <- which(image_metadata$counts[sel]=="yes")
# plot(r2, xlim=c(-2492000,-2487000), ylim=c(1268000,1273000))
# points(image_metadata[sel,5:6])
# points(image_metadata[sel[sel.cov],5:6],col="red",cex=2)
# points(image_metadata[sel[sel.cts],5:6],col="blue",cex=3)
# ## add other transects
# points(image_metadata[,5:6])
# points(image_metadata[which(image_metadata$cover=="yes"),5:6],col="red",cex=2)
# points(image_metadata[which(image_metadata$counts=="yes"),5:6],col="blue",cex=3)
# #text(image_metadata[which(image_metadata$cover=="yes"),5:6], labels=image_metadata[which(image_metadata$cover=="yes"),1], adj=-0.2)
# 
# #### aggregated cover of some living things on a single transect
# 
# ## projected coordinates for plotting
# img.coord <- SpatialPoints(coords=image_metadata[,5:6], proj4string=crs(stereo))
# 
# #### create subset of data to plot
# ## all images in that transect
# sel <- which(str_detect(image_metadata$Filename.standardised,"tan1901_065"))
# ## annotated images in that transect
# sel2 <- which(str_detect(image_metadata$Filename.standardised,"tan1901_065")&image_metadata$cover=="yes")
# # unique(image_metadata$cellID[sel])
# # unique(image_metadata$cellID[sel2])
# 
# #### just one transect
# sel.images <- grep("tan1901_065", rownames(cover_images))
# val <- rowSums(cover_images[sel.images,-c(1,3,4,6,7,10)])
# 
# 
# plot(r2, xlim=c(282500,283500),ylim=c(-2011500,-2010000))
# points(image_metadata[,5:6])
# points(image_metadata[which(image_metadata$cover=="yes"),5:6],col="red",cex=0.5)
# points(img.coord[sel2], cex=log(val), col="blue")
# 

##########################################################
library(lubridate)
load(paste0(ARC_Data.dir,"prep_image/Circumpolar_DownwardImages_metadata_2024.Rdata"))
load(file=paste0(ARC_Data.dir,"annotation/Circumpolar_Annotation_Data_500m_202312.Rdata"))

head(image_metadata)

### standardise image metadata for csv-files

# We need:
# * SurveyID  
# * TransectID  
# * Filename  
# * Filename-standardised  
# * Longitude  
# * Latitude  
# * Date taken 
# * Measured depth
# * Link to original source  
# * License  


##### first, read in all metadata and store as a list. Standardise depth to be a positive value
csv.list <- list()
csv.names <- c("Filename","Filename.standardised","Survey.ID","Transect.ID","Longitude","Latitude",
               "Date","Date_start","Date_end","Source.link","License",
               "Depth","Depth_start","Depth_end","HeightAboveSeafloor","Image_area","Image_area_source")

csv.list$CRS <- cbind(dat.list$WAP[,c('Filename','Filename.standardised')], 
                      "CRS", dat.list$WAP[,c('transectID','lon','lat')],
                      ymd(dat.list$WAP$Date), NA, NA, 
                      NA, NA,
                      dat.list$WAP[,c('depth.mean')], NA, NA, 2.5, 3, "metadata")

csv.list$NBP0808 <- csv.list$CRS[-which(csv.list$CRS$transectID%!in%NBP0808_transects),]
csv.list$NBP0808[,3] <- "NBP0808"
csv.list$NBP1001 <- csv.list$CRS[-which(csv.list$CRS$transectID%!in%NBP1001_transects),]
csv.list$NBP1001[,3] <- "NBP1001"
csv.list$CRS <- csv.list$CRS[-which(csv.list$CRS$transectID%in%c(NBP0808_transects,NBP1001_transects)),]

csv.list$PS06 <- cbind(dat.list$PS06[,c('Filename','Filename.standardised')], 
                       "PS06", dat.list$PS06[,c('transectID','lon','lat')],
                       NA, ymd_hms(dat.list$PS06$time_start), ymd_hms(dat.list$PS06$time_end), 
                       NA, "CC-BY-3.0",
                       dat.list$PS06[,c('depth')]*-1, NA, NA, NA, NA, NA)

csv.list$PS14 <- cbind(dat.list$PS14[,c('Filename','Filename.standardised')], 
                       "PS14",dat.list$PS14[,c('transectID','lon','lat')],
                       NA, ymd_hms(dat.list$PS14$time_start, truncated=3), ymd_hms(dat.list$PS14$time_end), 
                       NA, "CC-BY-3.0",
                       NA, dat.list$PS14[,c('depth_start','depth_end')]*-1, NA, 0.56, "metadata")

csv.list$PS18 <- cbind(dat.list$PS18[,c('Filename','Filename.standardised')],
                       "PS18",dat.list$PS18[,c('transectID','lon','lat')],
                       NA, ymd_hms(dat.list$PS18$time_start), ymd_hms(dat.list$PS18$time_end),
                       NA,"CC-BY-3.0",
                       NA,dat.list$PS18[,c('depth_start','depth_end')]*-1, NA,0.9, "metadata")

csv.list$PS61 <- cbind(dat.list$PS61[,c('Filename','Filename.standardised')],
                       "PS61",dat.list$PS61[,c('transectID','lon','lat')],
                       NA, ymd_hms(dat.list$PS61$time_start), ymd_hms(dat.list$PS61$time_end),
                       NA,"CC-BY-3.0",
                       NA,as.numeric(substr(dat.list$PS61$depth_start,1,6))*-1,
                       as.numeric(substr(dat.list$PS61$depth_end,1,6))*-1, NA,1, "metadata")

temp.PS81 <- cbind(dat.list$PS81[,c('Filename','Filename.standardised')],
                   "PS81", dat.list$PS81[,c('transectID','Longitude','Latitude')],
                   ymd_hms(dat.list$PS81$Date.Time), NA, NA,
                   "doi.pangaea.de/10.1594/PANGAEA.872719", "CC-BY-3.0",
                   dat.list$PS81[,c('Bathy.depth..m.')], NA, NA, NA, dat.list$PS81$Area_from_metadata, "metadata")
temp.PS81_s <- cbind(dat.list$PS81_shallow[,c('Filename','Filename.standardised')],
                     "PS81", dat.list$PS81_shallow[,c('transectID','Longitude','Latitude')],
                     ymd_hms(dat.list$PS81_shallow$Date.Time), NA, NA,
                     "doi.pangaea.de/10.1594/PANGAEA.872719", "CC-BY-3.0",
                     dat.list$PS81_shallow[,c('Bathy.depth..m.')], NA, NA, NA, dat.list$PS81_shallow$Area_from_metadata, "metadata")
names(temp.PS81) <- names(temp.PS81_s) <- csv.names
csv.list$PS81 <- rbind(temp.PS81,temp.PS81_s)

csv.list$PS96 <- cbind(dat.list$PS96[,c('Filename','Filename.standardised')],
                       "PS96", dat.list$PS96[,c('transectID','Longitude','Latitude')],
                       ymd_hms(dat.list$PS96$Date.Time), NA, NA,
                       "doi.pangaea.de/10.1594/PANGAEA.862097", "CC-BY-3.0",
                       dat.list$PS96[,c('Depth.water..m.')], NA, NA, dat.list$PS96[,c('Height..m.','Area')], "laserpoints")

csv.list$PS118 <- cbind(dat.list$PS118[,c('Filename','Filename.standardised')],
                        "PS118", dat.list$PS118[,c('transectID','Longitude','Latitude')], 
                        ymd_hms(dat.list$PS118$Date.Time, truncated=1), NA, NA,
                        "doi.pangaea.de/10.1594/PANGAEA.911904", "CC-BY-4.0",
                        dat.list$PS118[,c('Depth.water..m.')], NA, NA, NA, dat.list$PS118[,c('Area')], "metadata")
csv.list$PS118[1875:1972,7] <- ymd_hm(substr(dat.list$PS118$Date.Time[1875:1972],1,19))

csv.list$TAN0802 <- cbind(dat.list$TAN0802[,c('FileName','Filename.standardised')],
                          "TAN0802", dat.list$TAN0802[,c('transectID','GPS_lon','GPS_lat')],
                          ymd_hm(dat.list$TAN0802$time), NA, NA,
                          NA ,NA,
                          dat.list$TAN0802[,c('depth')], NA, NA, NA, NA, NA)

csv.list$TAN1802 <- cbind(dat.list$TAN1802[,c('FileName','Filename.standardised')],
                          "TAN1802", dat.list$TAN1802[,c('transectID','GPS_lon','GPS_lat')],
                          ymd_hm(dat.list$TAN1802$time), NA, NA,
                          NA, NA,
                          dat.list$TAN1802[,c('depth')], NA, NA, NA, NA, NA)

csv.list$TAN1901 <- cbind(dat.list$TAN1901[,c('FileName','Filename.standardised')],
                          "TAN1901",dat.list$TAN1901[,c('transectID','GPS_lon','GPS_lat')],
                          ymd_hm(dat.list$TAN1901$time), NA, NA,
                          NA, NA,
                          dat.list$TAN1901[,c('depth')], NA, NA, NA, NA, NA)

csv.list$NBP1402 <- cbind(dat.list$NBP1402[,c('FileName','Filename.standardised')],
                          "NBP1402", dat.list$NBP1402[,c('transectID','GPS_lon','GPS_lat')],
                          ymd_hms(dat.list$NBP1402$time), NA, NA,
                          "www.usap-dc.org/view/dataset/601310", "CC-BY-NC 4.0",
                          dat.list$NBP1402[,c('Depth')], NA, NA, 2.5, 4.8, "metadata")

csv.list$NBP1502 <- cbind(dat.list$NBP1502[,c('FileName','Filename.standardised')],
                          "NBP1502", dat.list$NBP1502[,c('transectID','GPS_lon','GPS_lat')],
                          ymd(dat.list$NBP1502$Date), NA, NA,
                          "doi.org/10.15784/601182", "CC-BY-NC 4.0",
                          NA, NA, NA, NA, NA, NA)

csv.list$AA2011 <- cbind(dat.list$AA2011[,c('Filename','Filename.standardised')],
                         "AA2011",dat.list$AA2011[,c('transectID','lon','lat')],
                         dmy_hms(dat.list$AA2011$DateTime),
                         NA, NA,
                         "doi.org/doi:10.4225/15/59acda196ccfb", "CC-BY-4.0",
                         dat.list$AA2011[,c('pressure')], NA, NA, dat.list$AA2011[,c('altimeter')], NA, NA)

csv.list$JR262 <- cbind(dat.list$JR262[,c('filename','Filename.standardised')],
                        "JR262",dat.list$JR262[,c('transectID','lon','lat')],
                        ymd(dat.list$JR262$time), NA, NA,
                        NA,NA,
                        dat.list$JR262[,c('depth')], NA, NA, 1.2, 0.51, "metadata")

csv.list$JR15005 <- cbind(dat.list$JR15005[,c('filename','Filename.standardised')],
                          "JR15005", dat.list$JR15005[,c('transectID','lon','lat')],
                          dat.list$JR15005$time, NA, NA,
                          NA,NA,
                          dat.list$JR15005[,c('depth')],NA,  NA, 1.2, 0.51, "metadata")

csv.list$JR17001 <- cbind(dat.list$JR17001[,c('filename','Filename.standardised')],
                          "JR17001", dat.list$JR17001[,c('transectID','lon','lat')],
                          dat.list$JR17001$time, NA, NA,
                          NA, NA,
                          dat.list$JR17001[,c('depth')], NA, NA, 1.2, 0.51, "metadata")

csv.list$JR17003 <- cbind(dat.list$JR17003[,c('filename','Filename.standardised')],
                          "JR17003", dat.list$JR17003[,c('transectID','lon','lat')],
                          dat.list$JR17003$time, NA, NA,
                          "doi.org/10.5285/48dcef16-6719-45e5-a335-3a97f099e451", "http://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/",
                          dat.list$JR17003[,c('depth')], NA, NA, 1.2, 0.51, "metadata")

csv.list$LMG1311 <- cbind(dat.list$LMG1311[,c('filename','Filename.standardised')],
                          "LMG1311",dat.list$LMG1311[,c('transectID','lon','lat')],
                          ymd_hms(dat.list$LMG1311$DateTime), NA, NA,
                          "doi.org/10.15784/601311", "CC-BY-NC 4.0",
                          NA, dat.list$LMG1311[,c('depth.start','depth.end')], 2.5, NA, NA)

for(i in 1:length(csv.list)){
  names(csv.list[[i]]) <- csv.names
}

##### some surveys need special attention
## PS06, PS14, PS18, PS61 all have a separate doi/link for each transect
unique(csv.list$PS06$Transect.ID)
csv.list$PS06$Source.link[csv.list$PS06$Transect.ID=="289"] <- "doi.pangaea.de/10.1594/PANGAEA.713324"
csv.list$PS06$Source.link[csv.list$PS06$Transect.ID=="292"] <- "doi.pangaea.de/10.1594/PANGAEA.713325"
unique(csv.list$PS14$Transect.ID)
csv.list$PS14$Source.link[csv.list$PS14$Transect.ID=="245"] <- "doi.pangaea.de/10.1594/PANGAEA.691559"
csv.list$PS14$Source.link[csv.list$PS14$Transect.ID=="246"] <- "doi.pangaea.de/10.1594/PANGAEA.691560"
csv.list$PS14$Source.link[csv.list$PS14$Transect.ID=="250"] <- "doi.pangaea.de/10.1594/PANGAEA.691561"
csv.list$PS14$Source.link[csv.list$PS14$Transect.ID=="256"] <- "doi.pangaea.de/10.1594/PANGAEA.691562"
csv.list$PS14$Source.link[csv.list$PS14$Transect.ID=="259"] <- "doi.pangaea.de/10.1594/PANGAEA.691563"
csv.list$PS14$Source.link[csv.list$PS14$Transect.ID=="260"] <- "doi.pangaea.de/10.1594/PANGAEA.691564"
csv.list$PS14$Source.link[csv.list$PS14$Transect.ID=="261"] <- "doi.pangaea.de/10.1594/PANGAEA.691565"
csv.list$PS14$Source.link[csv.list$PS14$Transect.ID=="270"] <- "doi.pangaea.de/10.1594/PANGAEA.691566"
csv.list$PS14$Source.link[csv.list$PS14$Transect.ID=="275"] <- "doi.pangaea.de/10.1594/PANGAEA.691568"
csv.list$PS14$Source.link[csv.list$PS14$Transect.ID=="276"] <- "doi.pangaea.de/10.1594/PANGAEA.691569"
csv.list$PS14$Source.link[csv.list$PS14$Transect.ID=="277"] <- "doi.pangaea.de/10.1594/PANGAEA.691570"
csv.list$PS14$Source.link[csv.list$PS14$Transect.ID=="278-2"] <- "doi.pangaea.de/10.1594/PANGAEA.691571"
csv.list$PS14$Source.link[csv.list$PS14$Transect.ID=="278-3"] <- "doi.pangaea.de/10.1594/PANGAEA.691572"
csv.list$PS14$Source.link[csv.list$PS14$Transect.ID=="280"] <- "doi.pangaea.de/10.1594/PANGAEA.691573"
csv.list$PS14$Source.link[csv.list$PS14$Transect.ID=="285"] <- "doi.pangaea.de/10.1594/PANGAEA.691574"
csv.list$PS14$Source.link[csv.list$PS14$Transect.ID=="293"] <- "doi.pangaea.de/10.1594/PANGAEA.691575"
csv.list$PS14$Source.link[csv.list$PS14$Transect.ID=="294"] <- "doi.pangaea.de/10.1594/PANGAEA.691576"
csv.list$PS14$Source.link[csv.list$PS14$Transect.ID=="304"] <- "doi.pangaea.de/10.1594/PANGAEA.691577"
csv.list$PS14$Source.link[csv.list$PS14$Transect.ID=="305"] <- "doi.pangaea.de/10.1594/PANGAEA.691578"
csv.list$PS14$Source.link[csv.list$PS14$Transect.ID=="306"] <- "doi.pangaea.de/10.1594/PANGAEA.691579"
csv.list$PS14$Source.link[csv.list$PS14$Transect.ID=="307"] <- "doi.pangaea.de/10.1594/PANGAEA.691580"
unique(csv.list$PS18$Transect.ID)
csv.list$PS18$Source.link[csv.list$PS18$Transect.ID=="126-2"] <- "doi.pangaea.de/10.1594/PANGAEA.667008"
csv.list$PS18$Source.link[csv.list$PS18$Transect.ID=="126-6"] <- "doi.pangaea.de/10.1594/PANGAEA.667009"
csv.list$PS18$Source.link[csv.list$PS18$Transect.ID=="129-1"] <- "doi.pangaea.de/10.1594/PANGAEA.667010"
csv.list$PS18$Source.link[csv.list$PS18$Transect.ID=="131-1"] <- "doi.pangaea.de/10.1594/PANGAEA.667011"
csv.list$PS18$Source.link[csv.list$PS18$Transect.ID=="134-1"] <- "doi.pangaea.de/10.1594/PANGAEA.667012"
csv.list$PS18$Source.link[csv.list$PS18$Transect.ID=="135-1"] <- "doi.pangaea.de/10.1594/PANGAEA.667013"
csv.list$PS18$Source.link[csv.list$PS18$Transect.ID=="136-4"] <- "doi.pangaea.de/10.1594/PANGAEA.667014"
csv.list$PS18$Source.link[csv.list$PS18$Transect.ID=="160-2"] <- "doi.pangaea.de/10.1594/PANGAEA.667015"
csv.list$PS18$Source.link[csv.list$PS18$Transect.ID=="165-1"] <- "doi.pangaea.de/10.1594/PANGAEA.667017"
csv.list$PS18$Source.link[csv.list$PS18$Transect.ID=="169-2"] <- "doi.pangaea.de/10.1594/PANGAEA.667018"
csv.list$PS18$Source.link[csv.list$PS18$Transect.ID=="171-1"] <- "doi.pangaea.de/10.1594/PANGAEA.667019"
csv.list$PS18$Source.link[csv.list$PS18$Transect.ID=="174-1"] <- "doi.pangaea.de/10.1594/PANGAEA.667021"
csv.list$PS18$Source.link[csv.list$PS18$Transect.ID=="175-2"] <- "doi.pangaea.de/10.1594/PANGAEA.667022"
csv.list$PS18$Source.link[csv.list$PS18$Transect.ID=="179-1"] <- "doi.pangaea.de/10.1594/PANGAEA.667023"
csv.list$PS18$Source.link[csv.list$PS18$Transect.ID=="180-2"] <- "doi.pangaea.de/10.1594/PANGAEA.667024"
csv.list$PS18$Source.link[csv.list$PS18$Transect.ID=="182-1"] <- "doi.pangaea.de/10.1594/PANGAEA.667025"
csv.list$PS18$Source.link[csv.list$PS18$Transect.ID=="189-2"] <- "doi.pangaea.de/10.1594/PANGAEA.667026"
csv.list$PS18$Source.link[csv.list$PS18$Transect.ID=="206-2"] <- "doi.pangaea.de/10.1594/PANGAEA.667027"
csv.list$PS18$Source.link[csv.list$PS18$Transect.ID=="211-1"] <- "doi.pangaea.de/10.1594/PANGAEA.667029"
csv.list$PS18$Source.link[csv.list$PS18$Transect.ID=="212-2"] <- "doi.pangaea.de/10.1594/PANGAEA.667030"
csv.list$PS18$Source.link[csv.list$PS18$Transect.ID=="220-1"] <- "doi.pangaea.de/10.1594/PANGAEA.667031"
csv.list$PS18$Source.link[csv.list$PS18$Transect.ID=="222-1"] <- "doi.pangaea.de/10.1594/PANGAEA.667032"
unique(csv.list$PS61$Transect.ID)
csv.list$PS61$Source.link[csv.list$PS61$Transect.ID=="235-1"] <- "doi.org/10.1594/PANGAEA.220746"
csv.list$PS61$Source.link[csv.list$PS61$Transect.ID=="249-1"] <- "doi.org/10.1594/PANGAEA.220747"


for(i in 1:21){
  message(names(csv.list)[i])
  print(head(csv.list[[i]]))
}


### replace area with annotated area where available and record this
for(i in 1:21){
  message(names(csv.list)[i])
  ## merge csv_list and image_metadata by filename and replace area values in csv_list
  dat.merged <- merge(csv.list[[i]][,c(2,16,17)], image_metadata[,c(1,10,11)], by="Filename.standardised", all.x=TRUE)
  
  ## Match the merged data back to the original dataframe order
  matched_data <- dat.merged[match(csv.list[[i]]$Filename.standardised, dat.merged$Filename.standardised),]
  csv.list[[i]]['Image_area'] <- ifelse(!is.na(matched_data$area),matched_data$area, matched_data$Image_area)
  csv.list[[i]]['Image_area_source'] <- ifelse(!is.na(matched_data$area), matched_data$area_source, matched_data$Image_area_source)
  # csv.list[[i]][c('Image_area','Image_area_source')] <- merge(csv.list[[i]][,c(2,16,17)], image_metadata[,c(1,10,11)], by="Filename.standardised", all.x=TRUE)[,3:4]
  print(head(csv.list[[i]]))
}

## check if all dim are the same (tan0802 and tan1901 are in caps vs lowercase...)
for(i in 1:21){
  survID <- names(csv.list)[i]
  message(survID)
  metadat.sel <- which(image_metadata$survey==survID)
  print(dim(image_metadata[metadat.sel,]))
  print(dim(csv.list[[i]]))
}

### limit csv-files to annotated images only:
sel.names <- names(csv.list)
sel.names[c(11,13)] <- c("tan0802","tan1901")
ann.csv.list <- list()
for(i in 1:21){
  survID <- sel.names[i]
  message(survID)
  metadat.sel <- which(image_metadata$survey==survID&image_metadata$cover=="yes")
  ann.sel <- image_metadata$Filename.standardised[metadat.sel]
  ann.csv.list[[i]] <- csv.list[[i]][csv.list[[i]]$Filename.standardised%in%ann.sel,]
}
names(ann.csv.list) <- names(csv.list)

## check if correct
for(i in 1:21){
  message(names(csv.list)[i])
  print(dim(ann.csv.list[[i]]))
  print(dim(csv.list[[i]]))
}

### dat.subset.list has a few images marked as annotated that are not:
## PS81 (1041 vs 1045)
dat.subset.list$PS81$Filename.standardised[which(dat.subset.list$PS81$Filename.standardised%!in%ann.csv.list$PS81$Filename.standardised)]
## TAN1901 (363 vs 365)
dat.subset.list$TAN1901$Filename.standardised[which(dat.subset.list$TAN1901$Filename.standardised%!in%ann.csv.list$TAN1901$Filename.standardised)]

### save individual csv files
r.dir <- "R:/IMAS/Antarctic_Seafloor/Clean_Data_For_Permanent_Storage/"
for(i in 1:21){
  survID <- sel.names[i]
  message(survID)
  write.csv(csv.list[[i]], file=paste0(r.dir,"metadata_full_dataset_",survID,".csv"), row.names=FALSE)
  write.csv(ann.csv.list[[i]], file=paste0(r.dir,"metadata_annotated_dataset_",survID,".csv"), row.names=FALSE)
}

# ### add height above seafloor to PS96 data
# csv.dat <- read.csv(paste0(r.dir,"metadata_full_dataset_PS96.csv"))
# ## using code from Readin_Circumpolar_DownwardImages_PS96.R in the prep_image folder
# data.start <- c(26,26,26,25,25,24,25,25,25,26,25,26,25)
# image.dir <- "E:/ARC_DP_data/a_RawData_DirectFromContributors/PS96_Piepenburg2016/"
# link <- "_links-to-photographs.tab"
# pre.chars <- "PS96_"
# transect.names <- c(
#   "001-4","007-1","008-2","010-3","026-3","027-2","037-3","048-2","057-3","061-1","072-4","090-4","106-2"
# )
# dat.header <- names(read.table(paste0(image.dir,pre.chars,transect.names[1],link), skip=24, header=TRUE, sep="\t", stringsAsFactors=FALSE, fill = TRUE)[,c(1:8)])
# dat <- read.table(paste0(image.dir,pre.chars,transect.names[1],link), skip=data.start[1]-1, header=FALSE, sep="\t", stringsAsFactors=FALSE, fill = TRUE)[,c(1:8)]
# dat[,9] <- as.character(transect.names[1])
# for(i in 2:length(transect.names)){
#   dat.temp <- read.table(paste0(image.dir,pre.chars,transect.names[i],link), skip=data.start[i]-1, header=FALSE, sep="\t",stringsAsFactors=FALSE, fill = TRUE)[,c(1:8)]
#   dat.temp[,9] <- as.character(transect.names[i])
#   dat <- rbind(dat, dat.temp)
# }
# names(dat) <- c(dat.header[1:5], "Area","Type","Filename","transectID")
# PS96.dat <- cbind(dat[,c(3,2,4,6,7)],NA,dat[,c(9,8,1,5)])
# names(PS96.dat)[6] <- "SurveyID"
# PS96.dat[6] <- "PS96"
# PS96.dat$SurveyID <- as.factor(PS96.dat$SurveyID)
# PS96.dat$transectID <- as.factor(PS96.dat$transectID)
# PS96.dat$Type <- as.factor(PS96.dat$Type)
# 
# ## now add height to csv data
# csv.dat$HeigthAboveSeafloor <- PS96.dat$Height..m.[match(csv.dat$Filename, PS96.dat$Filename)]
# head(csv.dat)
# write.csv(csv.dat, file=paste0(r.dir,"metadata_full_dataset_PS96.csv"), row.names=FALSE)
