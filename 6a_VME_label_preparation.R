#### Reading in and preparing VME annotations for upload to Squidle
'%!in%' <- function(x,y)!('%in%'(x,y))
so.dir <- "C:/Users/jjansen/OneDrive - University of Tasmania/science/SouthernOceanBiodiversityMapping/"
vme.dir <- paste0(so.dir,"ARC_Benthic_Mapping/vulnerable_marine_ecosystems/")

## load image annotations
vme.dat.raw <- read.csv(paste0(vme.dir,"data/biodata/biodata_step1.csv"))
unique(vme.dat.raw$label_names)

## load label translation to AMC scheme
label.translation <- read.csv(paste0(so.dir,"ARC_Data/annotation/VME_labels_vs_AMC_2025_07.csv"))

## remove mobile annotations to create a vme-only subset
lab.subset <- unique(vme.dat.raw$label_names)[-c(1:49,72:75, 79:102, 105:113,117,119,120,126:132)]
vme.dat <- vme.dat.raw[vme.dat.raw$label_names %in% lab.subset,]

## translate labels
vme.dat$label_amc <- NA
vme.dat$label_amc_id <- NA
for(i in 1:nrow(label.translation)){
  sel <- which(vme.dat$label_names==label.translation$name[i])
  vme.dat$label_amc[sel] <- label.translation$AMC[i]
  vme.dat$label_amc_id[sel] <- label.translation$AMC_ID[i]
}

any(is.na(vme.dat$label_amc))

####################
library(dplyr)
library(tidyr)
## PIXEL CENTRE COORDINATES MISSING FROM "biodata_step1.csv"
## Also, image width and height missing
## read in csv files from Biigle with image annotation reports to get the pixel coordinates of the annotation centres
csv.files <- list.files(paste0(vme.dir,"data/rawdata/292_csv_image_annotation_report/"), recursive=TRUE, full.names=TRUE)
csv.dat.list <- list()
for(i in 1:19){
  csv.dat.list[[i]] <- read.csv(csv.files[i])
}
csv.dat <- do.call(rbind, csv.dat.list)
## csv.dat has x,y,radius stored under $points in the format "[x,y,radius]"
## extract x,y,radius into separate columns
csv.dat <- csv.dat %>%
  separate(points, into=c("points.x", "points.y", "points.radius"), sep=",", remove=FALSE) %>%
  mutate(points.x = as.numeric(gsub("\\[","",points.x)),
         points.y = as.numeric(gsub("]","",points.y)),
         points.radius = as.numeric(gsub("]","",gsub(" ","",points.radius))))
## rectangles are defined by 8 coordinate values.
## When "shape_name is "rectangle", set calculate from coordinates the centre points.x, points.y, but set points.radius to NA
## Then, create points.x1, , points.y1, points.x2, points.y2, etc
rect.sel <- which(csv.dat$shape_name=="Rectangle")
csv.dat$points.x[rect.sel] <- NA
csv.dat$points.y[rect.sel] <- NA
csv.dat$points.radius[rect.sel] <- NA
for(i in rect.sel){
  coords <- as.numeric(unlist(strsplit(gsub("\\[|\\]","",csv.dat$points[i]),",")))
  csv.dat$points.x1[i] <- coords[1]
  csv.dat$points.y1[i] <- coords[2]
  csv.dat$points.x2[i] <- coords[3]
  csv.dat$points.y2[i] <- coords[4]
  csv.dat$points.x3[i] <- coords[5]
  csv.dat$points.y3[i] <- coords[6]
  csv.dat$points.x4[i] <- coords[7]
  csv.dat$points.y4[i] <- coords[8]
  csv.dat$points.x[i] <- mean(c(coords[1], coords[3], coords[5], coords[7]))
  csv.dat$points.y[i] <- mean(c(coords[2], coords[4], coords[6], coords[8]))
}
## there are 3 "LineString" annotations
## we translate them to points and radius
## take the center point of these coordinates (center of max in min along each of x and y) and calculate a radius as half the distance of the maximum extent in either x or y direction
line.sel <- which(csv.dat$shape_name=="LineString")
for(i in line.sel){
  coords <- as.numeric(unlist(strsplit(gsub("\\[|\\]","",csv.dat$points[i]),",")))
  x.coords <- coords[seq(1,length(coords),by=2)]
  y.coords <- coords[seq(2,length(coords),by=2)]
  ## take the mean of the maximum extents as center point
  csv.dat$points.x[i] <- mean(c(max(x.coords), min(x.coords)))
  csv.dat$points.y[i] <- mean(c(max(y.coords), min(y.coords)))
  x.range <- max(x.coords)-min(x.coords)
  y.range <- max(y.coords)-min(y.coords)
  csv.dat$points.radius[i] <- max(x.range,y.range)/2
}
##why are points.x1 etc filled for "LineString"?
csv.dat$shape_name[line.sel] <- "Circle"
csv.dat[line.sel,21:28] <- NA

#### now that all points are prepared, we can:
#### 1st: start matching all columns to the vme data
#### 2nd: prepare vme data to the format for upload to Squidle 
## match up with vme.dat based on annotation_id and add csv.dat$points to vme.dat
vme.dat$point.pixels.x <- NA
vme.dat$point.pixels.y <- NA
vme.dat$point.pixels.radius <- NA
for(i in 1:nrow(vme.dat)){
  sel <- which(csv.dat$annotation_id==vme.dat$annotation_id[i])
  ## if sel is not length 1, annotations in biostep1 have since been deleted and we skip to next loop
  if(length(sel)!=1) next
  vme.dat$point.pixels.x[i] <- csv.dat$points.x[sel]
  vme.dat$point.pixels.y[i] <- csv.dat$points.y[sel]
  vme.dat$point.pixels.radius[i] <- csv.dat$points.radius[sel]
  vme.dat$user.username[i] <- csv.dat$lastname[sel]
  ## point.media.deployment.key is the part of the filename up to the second underscore
  filename.parts <- unlist(strsplit(csv.dat$filename[i],"_"))[1:2]
  vme.dat$point.media.deployment.key[i] = paste0(filename.parts, collapse="_")
}
## remove 12442832 (Bryozoan LineString) from vme.dat
vme.dat <- vme.dat[-which(vme.dat$annotation_id==12442832),]
## change "LineString" to "Circle" in vme.dat$shape_name
vme.dat$shape_name[vme.dat$shape_name=="LineString"] <- "Circle"
## for the rectangles, we add coordinates as point.pixels.rectangle with the format "[[x1,y1],[x2,y2],[x3,y3],[x4,y4]]"
rect.sel.vme <- which(vme.dat$shape_name=="Rectangle")
vme.dat$point.pixels.rectangle <- NA
for(i in rect.sel.vme){
  sel <- which(csv.dat$annotation_id==vme.dat$annotation_id[i])
  coords <- c(csv.dat$points.x1[sel], csv.dat$points.y1[sel],
              csv.dat$points.x2[sel], csv.dat$points.y2[sel],
              csv.dat$points.x3[sel], csv.dat$points.y3[sel],
              csv.dat$points.x4[sel], csv.dat$points.y4[sel])
  poly.string <- paste0("[[",coords[1],",",coords[2],"],[",
                        coords[3],",",coords[4],"],[",
                        coords[5],",",coords[6],"],[",
                        coords[7],",",coords[8],"]]")
  vme.dat$point.pixels.rectangle[i] <- poly.string
}
head(vme.dat)

## change image_filename to point.media.key
vme.dat <- vme.dat %>%
  rename(point.media.key = image_filename)


#################################################################
#### now copying the code from the mobile species to create polygons for each circle annotation
#### "5c_CreatePolygonsForCircleAnnotations.R"
## Load necessary libraries
library(dplyr)
library(tidyr)
library(readr)
## Read the mobiles CSV file, it contains information about relevant image dimensions
asaid.dir <- "R:/IMAS/Antarctic_Seafloor/ASAID_data/AnnotationLibrary_AllFinishedSurveys/"
#annotations <- read.csv(paste0(asaid.dir,"Counts/Circumpolar_DownwardImages_ExhaustiveSearch_Annotations_ForSquidle_202507.csv"))
annotations <- read.csv(paste0(asaid.dir,"Cover/Circumpolar_DownwardImages_PointScore_Annotations_ForSquidle_202507.csv"))

## adding point.pixels.width and point.pixels.height to vme.dat
## we need to look up for each unique "point.media.key" in annotations the widths and heights, and then add these to vme.dat by "image_filename"
unique_images <- unique(vme.dat$point.media.key)
image_dims <- data.frame(point.media.key = character(),
                         point.pixels.width = numeric(),
                         point.pixels.height = numeric(),
                         stringsAsFactors = FALSE)
vme.dat$point.pixels.width <- NA
vme.dat$point.pixels.height <- NA
for(i in 1:length(unique_images)){
  img_key <- unique_images[i]
  sel <- which(annotations$Name == img_key)
  img_width <- annotations$Width[sel][1]
  img_height <- annotations$Height[sel][1]
  image_dims <- rbind(image_dims,
                      data.frame(point.media.key = img_key,
                                 point.pixels.width = img_width,
                                 point.pixels.height = img_height,
                                 stringsAsFactors = FALSE))
  ## add these dimensions to vme.dat
  sel_dim <- which(vme.dat$point.media.key == img_key)
  sel_img <- which(image_dims$point.media.key == img_key)
  vme.dat$point.pixels.width[sel_dim] <- image_dims$point.pixels.width[sel_img]
  vme.dat$point.pixels.height[sel_dim] <- image_dims$point.pixels.height[sel_img]
}
## add the following ones because they are missing:
# PS81_159_0005__PS81_159-1_2013-02-08T04_44_24_5900.jpg #Width=5616, Height=3744
# PS81_159_0017__PS81_159-1_2013-02-08T05_10_20_5961.jpg #Width=5616, Height=3744
# PS81_163_0003__PS81_163-2_2013-02-10T20_35_00_7468.jpg #Width=5616, Height=3744
# PS81_163_0017__PS81_163-2_2013-02-10T20_23_28_7443.jpg #Width=5616, Height=3744
# CRS_1103_0007__CRS_1103_104.jpg #Width=3872, Height=2592
# CRS_1103_0001__CRS_1103_060.jpg #Width=3872, Height=2592
PS18.sel <- which(vme.dat$point.media.key%in%c(
  "PS81_159_0005__PS81_159-1_2013-02-08T04_44_24_5900.jpg",
  "PS81_159_0017__PS81_159-1_2013-02-08T05_10_20_5961.jpg",
  "PS81_163_0003__PS81_163-2_2013-02-10T20_35_00_7468.jpg",
  "PS81_163_0017__PS81_163-2_2013-02-10T20_23_28_7443.jpg"))
vme.dat$point.pixels.width[PS18.sel] <- 5616
vme.dat$point.pixels.height[PS18.sel] <- 3744
CRS.sel <- which(vme.dat$point.media.key%in%c(
  "CRS_1103_0007__CRS_1103_104.jpg",
  "CRS_1103_0001__CRS_1103_060.jpg"))
vme.dat$point.pixels.width[CRS.sel] <- 3872
vme.dat$point.pixels.height[CRS.sel] <- 2592

missing_dims <- vme.dat %>%
  filter(is.na(point.pixels.width) | is.na(point.pixels.height)) %>%
  distinct(point.media.key)

print(missing_dims)

#################################################
#### AA2011 is special, because of rotated images
## filter vme.dat for AA2011 images only
vme.dat.AA2011 <- vme.dat %>%
  filter(grepl("AA2011",point.media.key))
rotated_images <- c("AA2011_CTD54_CAM02_0001__IMG_0233.JPG", "AA2011_CTD56_CAM05_0001__IMG_0358.JPG")

## adjust for cropped part and re-calculate positions for rotated images
annotations_uncropped <- vme.dat.AA2011 %>%
  mutate(
    point.pixels.x = point.pixels.x + 0.05 * point.pixels.width,
    point.pixels.y = point.pixels.y + 0.05 * point.pixels.height,
    point.x = if_else(point.media.key %in% rotated_images,
                      point.pixels.y / point.pixels.height,
                      point.pixels.x / point.pixels.width),
    point.y = if_else(point.media.key %in% rotated_images,
                      1 - (point.pixels.x / point.pixels.width),
                      point.pixels.y / point.pixels.height)
  )

## Function to generate relative polygon coordinates
generate_relative_polygon <- function(radius, width, height, media_key, rotated_images, n_points = 16) {
  # Swap width and height if image is rotated
  if (media_key %in% rotated_images) {
    temp <- width
    width <- height
    height <- temp
  }
  r_x <- radius / width
  r_y <- radius / height
  angles <- seq(0, 2 * pi, length.out = n_points + 1)[-1]
  x_offsets <- round(r_x * cos(angles), 10)
  y_offsets <- round(r_y * sin(angles), 10)
  coords <- Map(function(x, y) c(x, y), x_offsets, y_offsets)
  list(point.polygon = coords)
}

## Apply the function and unpack the list into a new column
polygon_data <- annotations_uncropped %>%
  rowwise() %>%
  mutate(point.polygon = paste0("[", paste(
    sapply(generate_relative_polygon(point.pixels.radius, point.pixels.width, point.pixels.height, point.media.key, rotated_images)$point.polygon,
           function(pt) paste0("[", pt[1], ",", pt[2], "]")),
    collapse = ", "), "]")) %>%
  ungroup()

## Write to CSV
vme.save.dir <- "R:/IMAS/Antarctic_Seafloor/ASAID_data/AnnotationLibrary_AllFinishedSurveys/VME/"
write_csv(polygon_data, file=paste0(vme.save.dir,"Circumpolar_DownwardImages_VME_Annotations_ForSquidle_202511_AA2011.csv"))

#################################################
## Function to generate relative polygon coordinates
generate_relative_polygon <- function(radius, width, height, n_points = 16) {
  r_x <- radius / width
  r_y <- radius / height
  angles <- seq(0, 2 * pi, length.out = n_points + 1)[-1]
  x_offsets <- round(r_x * cos(angles), 10)
  y_offsets <- round(r_y * sin(angles), 10)
  coords <- Map(function(x, y) c(x, y), x_offsets, y_offsets)
  list(point.polygon = coords)
}

#### SURVEYS WITH 20% CROP (PS81, PS96, PS118)
annotations_uncropped <- vme.dat
annotations_uncropped$point.pixels.x <- vme.dat$point.pixels.x + 0.1*vme.dat$point.pixels.width
annotations_uncropped$point.pixels.y <- vme.dat$point.pixels.y + 0.1*vme.dat$point.pixels.height
annotations_uncropped$point.x <- annotations_uncropped$point.pixels.x/vme.dat$point.pixels.width
annotations_uncropped$point.y <- annotations_uncropped$point.pixels.y/vme.dat$point.pixels.height
## Apply the function and unpack the list into a new column
polygon_data <- annotations_uncropped %>%
  rowwise() %>%
  mutate(point.polygon = paste0("[", paste(
    sapply(generate_relative_polygon(point.pixels.radius, point.pixels.width, point.pixels.height)$point.polygon,
           function(pt) paste0("[", pt[1], ",", pt[2], "]")),
    collapse = ", "), "]")) %>%
  ungroup()
## subset to PS81, PS96 and PS118 surveys
polygon_data_PS81 <- polygon_data %>%  filter(grepl("PS81", point.media.key))
polygon_data_PS96 <- polygon_data %>%  filter(grepl("PS96", point.media.key))
polygon_data_PS118 <- polygon_data %>%  filter(grepl("PS118", point.media.key))
## subset and write to csv
write_csv(polygon_data_PS81, paste0(asaid.dir, "VME/Circumpolar_DownwardImages_VME_Annotations_ForSquidle_202511_PS81.csv"))
write_csv(polygon_data_PS96, paste0(asaid.dir, "VME/Circumpolar_DownwardImages_VME_Annotations_ForSquidle_202511_PS96.csv"))
write_csv(polygon_data_PS118, paste0(asaid.dir, "VME/Circumpolar_DownwardImages_VME_Annotations_ForSquidle_202511_PS118.csv"))

##########################################
#### now process all other surveys
annotations_uncropped <- vme.dat
annotations_uncropped$point.x <- annotations_uncropped$point.pixels.x/vme.dat$point.pixels.width
annotations_uncropped$point.y <- annotations_uncropped$point.pixels.y/vme.dat$point.pixels.height
## Apply the function and unpack the list into a new column
polygon_data <- annotations_uncropped %>%
  rowwise() %>%
  mutate(point.polygon = paste0("[", paste(
    sapply(generate_relative_polygon(point.pixels.radius, point.pixels.width, point.pixels.height)$point.polygon,
           function(pt) paste0("[", pt[1], ",", pt[2], "]")),
    collapse = ", "), "]")) %>%
  ungroup()
## subset to all other surveys (PS06, PS14, PS18, PS61, NBP1402, NBP1502, LMG1311, TAN0802, TAN1802, TAN1901, CRS, JR262, JR15005, JR17001, JR17003)

####

polygon_dat <- polygon_data %>%  filter(grepl("PS06", point.media.key))
write_csv(polygon_dat, paste0(asaid.dir, "VME/Circumpolar_DownwardImages_VME_Annotations_ForSquidle_202511_PS06.csv"))
polygon_dat <- polygon_data %>%  filter(grepl("PS14", point.media.key))
write_csv(polygon_dat, paste0(asaid.dir, "VME/Circumpolar_DownwardImages_VME_Annotations_ForSquidle_202511_PS14.csv"))
polygon_dat <- polygon_data %>%  filter(grepl("PS18", point.media.key))
write_csv(polygon_dat, paste0(asaid.dir, "VME/Circumpolar_DownwardImages_VME_Annotations_ForSquidle_202511_PS18.csv"))
polygon_dat <- polygon_data %>%  filter(grepl("PS61", point.media.key))
write_csv(polygon_dat, paste0(asaid.dir, "VME/Circumpolar_DownwardImages_VME_Annotations_ForSquidle_202511_PS61.csv"))

polygon_dat <- polygon_data %>%  filter(grepl("tan0802", point.media.key))
write_csv(polygon_dat, paste0(asaid.dir, "VME/Circumpolar_DownwardImages_VME_Annotations_ForSquidle_202511_TAN0802.csv"))
polygon_dat <- polygon_data %>%  filter(grepl("TAN1802", point.media.key))
write_csv(polygon_dat, paste0(asaid.dir, "VME/Circumpolar_DownwardImages_VME_Annotations_ForSquidle_202511_TAN1802.csv"))
polygon_dat <- polygon_data %>%  filter(grepl("tan1901", point.media.key))
write_csv(polygon_dat, paste0(asaid.dir, "VME/Circumpolar_DownwardImages_VME_Annotations_ForSquidle_202511_TAN1901.csv"))

polygon_dat <- polygon_data %>% filter(grepl("JR262", point.media.key))
write_csv(polygon_dat, paste0(asaid.dir, "VME/Circumpolar_DownwardImages_VME_Annotations_ForSquidle_202511_JR262.csv"))
polygon_dat <- polygon_data %>% filter(grepl("JR15005", point.media.key))
write_csv(polygon_dat, paste0(asaid.dir, "VME/Circumpolar_DownwardImages_VME_Annotations_ForSquidle_202511_JR15005.csv"))
polygon_dat <- polygon_data %>% filter(grepl("JR17001", point.media.key))
write_csv(polygon_dat, paste0(asaid.dir, "VME/Circumpolar_DownwardImages_VME_Annotations_ForSquidle_202511_JR17001.csv"))
polygon_dat <- polygon_data %>% filter(grepl("JR17003", point.media.key))
write_csv(polygon_dat, paste0(asaid.dir, "VME/Circumpolar_DownwardImages_VME_Annotations_ForSquidle_202511_JR17003.csv"))

polygon_dat <- polygon_data %>% filter(grepl("NBP1402", point.media.key))
write_csv(polygon_dat, paste0(asaid.dir, "VME/Circumpolar_DownwardImages_VME_Annotations_ForSquidle_202511_NBP1402.csv"))
polygon_dat <- polygon_data %>% filter(grepl("NBP1502", point.media.key))
write_csv(polygon_dat, paste0(asaid.dir, "VME/Circumpolar_DownwardImages_VME_Annotations_ForSquidle_202511_NBP1502.csv"))
polygon_dat <- polygon_data %>% filter(grepl("LMG1311", point.media.key))
write_csv(polygon_dat, paste0(asaid.dir, "VME/Circumpolar_DownwardImages_VME_Annotations_ForSquidle_202511_LMG1311.csv"))

## subset to CRS
polygon_dat <- polygon_data %>% filter(grepl("CRS", point.media.key))
## create transect ID column based on the number between the first two underscores in point.media.key 
polygon_dat$transectID <- sapply(polygon_dat$point.media.key, function(x) {
  parts <- unlist(strsplit(x, "_"))
  return(parts[2])
})
##
sel.lmg0902 <- which(polygon_dat$transectID%in%c(1207,1208,1217,1219,1255,1267))
sel.nbp0808 <- which(polygon_dat$transectID%in%c(1069,1072,1091,1130,1132))
sel.nbp1001 <- (1:nrow(polygon_dat))[-c(sel.lmg0902, sel.nbp0808)]

polygon_dat.nbp0808 <- polygon_dat[sel.nbp0808,]
polygon_dat.lmg0902 <- polygon_dat[sel.lmg0902,]
polygon_dat.nbp1001 <- polygon_dat[sel.nbp1001,]

write_csv(polygon_dat.nbp0808, paste0(asaid.dir, "VME/Circumpolar_DownwardImages_VME_Annotations_ForSquidle_202511_NBP0808.csv"))
write_csv(polygon_dat.lmg0902, paste0(asaid.dir, "VME/Circumpolar_DownwardImages_VME_Annotations_ForSquidle_202511_LMG0902.csv"))
write_csv(polygon_dat.nbp1001, paste0(asaid.dir, "VME/Circumpolar_DownwardImages_VME_Annotations_ForSquidle_202511_NBP1001.csv"))












## adding point.pixels.width and point.pixels.height to vme.dat.AA2011
## we need to look up for each unique "point.media.key" in annotations the widths and heights, and then add these to vme.dat.AA2011 by "image_filename"
unique_images <- unique(vme.dat.AA2011$point.media.key)
image_dims <- data.frame(point.media.key = character(),
                         point.pixels.width = numeric(),
                         point.pixels.height = numeric(),
                         stringsAsFactors = FALSE)
vme.dat.AA2011$point.pixels.width <- NA
vme.dat.AA2011$point.pixels.height <- NA
for(i in 1:length(unique_images)){
  img_key <- unique_images[i]
  sel <- which(annotations$point.media.key == img_key)
  img_width <- annotations$point.pixels.width[sel][1]
  img_height <- annotations$point.pixels.height[sel][1]
  image_dims <- rbind(image_dims,
                      data.frame(point.media.key = img_key,
                                 point.pixels.width = img_width,
                                 point.pixels.height = img_height,
                                 stringsAsFactors = FALSE))
  ## add these dimensions to vme.dat.AA2011
  sel_dim <- which(vme.dat.AA2011$point.media.key == img_key)
  sel_img <- which(image_dims$point.media.key == img_key)
  vme.dat.AA2011$point.pixels.width[sel_dim] <- image_dims$point.pixels.width[sel_img]
  vme.dat.AA2011$point.pixels.height[sel_dim] <- image_dims$point.pixels.height[sel_img]
}


