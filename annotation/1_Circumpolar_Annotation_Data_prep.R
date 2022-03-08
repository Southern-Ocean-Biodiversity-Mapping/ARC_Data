### 1) ### Setting up----
library(raster)
library(readxl)
library(readr)
library(dplyr)
library(data.table)
library(proj4)
library(stringr)

user = "Jan"
#user = "charley"
#user="nicole"

if (user == "Jan") {
  
  sci.dir <-      "C:/Users/jjansen/Desktop/science/"
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

## R-drive paths
RS.dir <- "R:/IMAS/Antarctic_Seafloor/Clean_Data_For_Permanent_Storage/"
ann.dir <- paste0(RS.dir,"AnnotationLibrary_AllFinishedSurveys/")

# Image quality path
image.quality.path <- "R:/IMAS/Antarctic_Seafloor/image_quality_analysis/image_quality_score.csv"

##### load still and diatom sample locations and bathymetry:

## from "Readin_Circumpolar_DownwardImage_Data.Rmd"
load(paste0(ARC_Data.dir,"prep_image/Circumpolar_DownwardImages_metadata.Rdata"))

## from "ReadIn_Circumpolar_Environmental_Data.Rmd"
r2 <- raster(paste0(env.derived,"Circumpolar_EnvData_500m_shelf_bathy_gebco_depth.grd"))

##### load coastline
#stereo <- as.character(r2@crs)
#load(paste0(env.derived,"Circumpolar_Coastline.Rdata"))

source(paste0(tools.dir,"bubbleplot.R"))



### 2) Data Preparation - Cover is a corrected file, but counts is still BASED ON UNCORRECTED INDIVIDUAL SURVEY FILES ----

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
image_metadata <- select(dat.list.clean[[1]],Filename.standardised,lon,lat,transectID)
for(i in 2:length(dat.list.clean)){
  message(i)
  print(names(dat.list.clean[[i]]))
  if("lon" %in% names(dat.list.clean[[i]])){ dat.temp <- select(dat.list.clean[[i]],Filename.standardised,lon,lat,transectID) }
  if("Lon" %in% names(dat.list.clean[[i]])){ dat.temp <- select(dat.list.clean[[i]],Filename.standardised,Lon,Lat,transectID) }
  if("Longitude" %in% names(dat.list.clean[[i]])){ dat.temp <- select(dat.list.clean[[i]],Filename.standardised,Longitude,Latitude,transectID) }
  if("GPS_lon" %in% names(dat.list.clean[[i]])){ dat.temp <- select(dat.list.clean[[i]],Filename.standardised,GPS_lon,GPS_lat,transectID) }  
  names(dat.temp) <- c("Filename.standardised","lon","lat","transectID")
  image_metadata <- rbind(image_metadata,dat.temp)
}
image_metadata$Filename.standardised <- gsub(".tif",".jpg",image_metadata$Filename.standardised)
image_metadata$proj_coord_x <- project(image_metadata[,2:3], proj=crs(r2))$x
image_metadata$proj_coord_y <- project(image_metadata[,2:3], proj=crs(r2))$y
image_metadata$cellID <- extract(r2, image_metadata[,5:6], cellnumbers=TRUE)[,1]

image_metadata$survey <- str_split(image_metadata$Filename.standardised,"_", simplify=T)[,1]

ids <- unique(image_metadata$cellID)

# Init area and area_source
image_metadata$area <- NA
image_metadata$area_source <- NA
# Pulling area info from metadata
for (survey_current in names(dat.list.clean)) {
  if ("Area" %in% names(dat.list.clean[[survey_current]])) {
    for (fname_current in dat.list.clean[[survey_current]]$Filename.standardised) {
      area_current <- dat.list.clean[[survey_current]][dat.list.clean[[survey_current]]$Filename.standardised == fname_current, ]$Area
      if (!is.na(area_current)) {
        image_metadata[image_metadata$Filename.standardised == fname_current, "area"] <- area_current
        image_metadata[image_metadata$Filename.standardised == fname_current, "area_source"] <- "metadata"
      }
    }
  }
  if ("Area_from_metadata" %in% names(dat.list.clean[[survey_current]])) {
    for (fname_current in dat.list.clean[[survey_current]]$Filename.standardised) {
      area_current <- dat.list.clean[[survey_current]][dat.list.clean[[survey_current]]$Filename.standardised == fname_current, ]$Area_from_metadata
      if (!is.na(area_current)) {
        image_metadata[image_metadata$Filename.standardised == fname_current, "area"] <- area_current
        image_metadata[image_metadata$Filename.standardised == fname_current, "area_source"] <- "metadata"
      }
    }
  }
  if ("area" %in% names(dat.list.clean[[survey_current]])) {
    for (fname_current in dat.list.clean[[survey_current]]$Filename.standardised) {
      area_current <- dat.list.clean[[survey_current]][dat.list.clean[[survey_current]]$Filename.standardised == fname_current, ]$area
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
image_metadata$cover <- "no"
image_metadata$counts <- "no"
# find matching filenames for cover
idx <- match(image_metadata$Filename.standardised,unique(dat_cover$Name))
fill.idx <- which(!is.na(idx))
search.idx <- idx[fill.idx]
image_metadata$cover[fill.idx] <- "yes" # replace area values where filenames match
# find matching filenames for counts
idx <- match(image_metadata$Filename.standardised,unique(dat_counts$filename))
fill.idx <- which(!is.na(idx))
search.idx <- idx[fill.idx]
image_metadata$counts[fill.idx] <- "yes" # replace area values where filenames match
# check:
image_metadata[which(image_metadata$cover=="yes"),]

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

for(i in c(1,8,15)){
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
counts_cells_survey <- rep(NA, length(ids))
counts_cells_transect1 <- rep(NA, length(ids))
counts_cells_transect2 <- rep(NA, length(ids))
counts_cells_transect3 <- rep(NA, length(ids))

for(i in 1:length(ids)){
  #print(i)
  sel.r <- which(image_metadata$cellID==ids[i]&image_metadata$counts=="yes") # find images that are part of that cell and annotated in biigle
  if(length(sel.r)==0) next # if none of the images are annotated, skip to the next iteration
  sel.names <- image_metadata$Filename.standardised[sel.r] # find the names of these images
  dat.temp <- dat_counts_image_by_species[which(nams.counts%in%sel.names),] # find annotations from these images from the image-dataset
  counts_N[i] <-  nrow(dat.temp)
  counts_area[i] <- sum(image_metadata$area[sel.r])
  dat_counts_cell_by_species[i,] <- colSums(dat.temp)
  ## check if all images are from the same survey and same transect
  counts_cells_survey[i] <- unique(image_metadata$survey[sel.r])
  counts_cells_transect1[i] <- unique(image_metadata$transectID[sel.r])[1]
  print(length(unique(image_metadata$transectID[sel.r])))
  if(length(unique(image_metadata$transectID[sel.r]))>1){
    counts_cells_transect2[i] <- unique(image_metadata$transectID[sel.r])[2]
    if(length(unique(image_metadata$transectID[sel.r]))>2){
      counts_cells_transect3[i] <- unique(image_metadata$transectID[sel.r])[3]
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


### 4b) site(cell)-by-species matrix

## create new dataframe using image locations
dat_cover_cell_by_species <- data.frame(matrix(NA,nrow=length(ids),ncol=length(labs.cov)))
rownames(dat_cover_cell_by_species) <- ids
colnames(dat_cover_cell_by_species) <- labs.cov
cover_N <- rep(NA, length(ids))
cover_area <- rep(NA, length(ids))
cover_cells_survey <- rep(NA, length(ids))
cover_cells_transect1 <- rep(NA, length(ids))
cover_cells_transect2 <- rep(NA, length(ids))
cover_cells_transect3 <- rep(NA, length(ids))
for(i in 1:length(ids)){
  #print(i)
  sel.r <- which(image_metadata$cellID==ids[i]&image_metadata$cover=="yes") # find images that are part of that cell and annotated in coralnet
  if(length(sel.r)==0) next # if none of the images are annotated, skip to the next iteration
  sel.names <- image_metadata$Filename.standardised[sel.r] # find the names of these images
  dat.temp <- dat_cover_image_by_species[which(nams.cov%in%sel.names),] # find annotations from these images from the image-dataset
  cover_N[i] <-  nrow(dat.temp)
  cover_area[i] <- sum(image_metadata$area[sel.r])
  dat_cover_cell_by_species[i,] <- colSums(dat.temp)
  ## check if all images are from the same survey and same transect
  cover_cells_survey[i] <- unique(image_metadata$survey[sel.r])
  cover_cells_transect1[i] <- unique(image_metadata$transectID[sel.r])[1]
  if(length(unique(image_metadata$transectID[sel.r]))>1){
    cover_cells_transect2[i] <- unique(image_metadata$transectID[sel.r])[2]
    if(length(unique(image_metadata$transectID[sel.r]))>2){
      cover_cells_transect3[i] <- unique(image_metadata$transectID[sel.r])[3]
    }}
}

### 4b) cell metadata

cell.coords <- xyFromCell(r2, ids)
cell.lonlat <- project(cell.coords, proj=crs(r2), inverse=TRUE) 

cell_metadata <- data.frame(cbind(ids,cell.lonlat,cell.coords,cover_N, counts_N, cover_area, counts_area, 
                                  cover_cells_survey, cover_cells_transect1, cover_cells_transect2, cover_cells_transect3, 
                                  counts_cells_survey, counts_cells_transect1, counts_cells_transect2, counts_cells_transect3))
names(cell_metadata) <- c("cellID", "lon", "lat", "proj_coord_x", "proj_coord_y", "cover_N", "counts_N", "cover_area", "counts_area",
                          "cover_cells_survey", "cover_cells_transect1", "cover_cells_transect2", "cover_cells_transect3", 
                          "counts_cells_survey", "counts_cells_transect1", "counts_cells_transect2", "counts_cells_transect3")

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
cell_metadata$year[cell_metadata$cover_cells_survey=="CRS"] <- 2010
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
cell_metadata$gear[cell_metadata$cover_cells_survey=="CRS"] <- "YOYO"
cell_metadata$gear[cell_metadata$cover_cells_survey=="NBP1402"] <- "YOYO"
cell_metadata$gear[cell_metadata$cover_cells_survey=="NBP1502"] <- "YOYO"
cell_metadata$gear[cell_metadata$cover_cells_survey=="LMG1311"] <- "YOYO"
cell_metadata$gear[cell_metadata$cover_cells_survey=="JR262"] <- "SUCS"
cell_metadata$gear[cell_metadata$cover_cells_survey=="JR15005"] <- "SUCS"
cell_metadata$gear[cell_metadata$cover_cells_survey=="JR17001"] <- "SUCS"
cell_metadata$gear[cell_metadata$cover_cells_survey=="JR17003"] <- "SUCS"

for(i in c(1,10:17,19)){
  cell_metadata[,i] <- as.factor(cell_metadata[,i])
}

### 5) save output----
cover_cells <- dat_cover_cell_by_species[-which(is.na(rowSums(dat_cover_cell_by_species))),]
count_cells <- dat_counts_cell_by_species[-which(is.na(rowSums(dat_counts_cell_by_species))),]
cover_images <- dat_cover_image_by_species
count_images <- dat_counts_image_by_species

head(image_metadata)
head(cell_metadata)

save(image_metadata, cell_metadata, cover_cells, cover_images, count_cells, count_images,
     file=paste0(ARC_Data.dir,"annotation/Circumpolar_Annotation_Data.Rdata"))



### 6) Inspect resulting dataframes ----

### 6a) locations across a region

## check if things are correct:
plot(r2, xlim=c(0,300000),ylim=c(-2100000,-1800000))
points(image_metadata[,5:6])


### 6b) a single transect

plot(r2, xlim=c(282500,283500),ylim=c(-2011500,-2010000))
points(image_metadata[,5:6])
points(image_metadata[which(image_metadata$cover=="yes"),5:6],col="red",cex=2)
points(image_metadata[which(image_metadata$counts=="yes"),5:6],col="blue",cex=3)
text(image_metadata[which(image_metadata$cover=="yes"),5:6], labels=image_metadata[which(image_metadata$cover=="yes"),1], adj=-0.2)


#### aggregated cover of some living things on a single transect

## projected coordinates for plotting
img.coord <- SpatialPoints(coords=image_metadata[,5:6], proj4string=crs(stereo))

#### create subset of data to plot
## all images in that transect
sel <- which(str_detect(image_metadata$Filename.standardised,"tan1901_065"))
## annotated images in that transect
sel2 <- which(str_detect(image_metadata$Filename.standardised,"tan1901_065")&image_metadata$cover=="yes")
# unique(image_metadata$cellID[sel])
# unique(image_metadata$cellID[sel2])

#### just one transect
sel.images <- grep("tan1901_065", rownames(cover_images))
val <- rowSums(cover_images[sel.images,-c(1,3,4,6,7,10)])


plot(r2, xlim=c(282500,283500),ylim=c(-2011500,-2010000))
points(image_metadata[,5:6])
points(image_metadata[which(image_metadata$cover=="yes"),5:6],col="red",cex=0.5)
points(img.coord[sel2], cex=log(val), col="blue")






