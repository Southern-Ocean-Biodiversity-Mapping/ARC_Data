##### CREATE ANNOTATION FILES FOR UPLOAD TO SQUIDLE+ ######
## Load necessary libraries
library(dplyr)
library(tidyr)
library(readr)
counts.dir <- "R:/IMAS/Antarctic_Seafloor/ASAID_data/AnnotationLibrary_AllFinishedSurveys/Counts/"

#################################################
#### AA2011 is special, because of rotated images
## Read the CSV file
annotations <- read_csv(paste0(counts.dir,"Circumpolar_DownwardImages_ExhaustiveSearch_Annotations_ForSquidle_202502_AA2011.csv"))
rotated_images <- c("AA2011_CTD54_CAM02_0001__IMG_0233.JPG", "AA2011_CTD56_CAM05_0001__IMG_0358.JPG")
## adjust for cropped part and re-calculate positions for rotated images
annotations_uncropped <- annotations %>%
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
    sapply(generate_relative_polygon(radius, point.pixels.width, point.pixels.height, point.media.key, rotated_images)$point.polygon,
           function(pt) paste0("[", pt[1], ",", pt[2], "]")),
    collapse = ", "), "]")) %>%
  ungroup()

## Write to CSV
write_csv(polygon_data, paste0(counts.dir, "Circumpolar_DownwardImages_ExhaustiveSearch_Annotations_ForSquidle_202507_AA2011.csv"))


#################################################################################
#### SURVEYS WITH 20% CROP (PS81, PS96, PS118)
## Read the CSV file
annotations <- read_csv(paste0(counts.dir,"Circumpolar_DownwardImages_ExhaustiveSearch_Annotations_ForSquidle_202507.csv"))
annotations <- annotations %>%
  rename(point.media.key = Name,
         label.id = AMC_ID,
         user.username = Annotator,
         point.media.deployment.key = TransectID)

annotations_uncropped <- annotations
# annotations_uncropped$point.pixels.x <- annotations$point.pixels.x + 0.1*annotations$point.pixels.width
# annotations_uncropped$point.pixels.y <- annotations$point.pixels.y + 0.1*annotations$point.pixels.height
# annotations_uncropped$point.x <- annotations_uncropped$point.pixels.x/annotations$point.pixels.width
# annotations_uncropped$point.y <- annotations_uncropped$point.pixels.y/annotations$point.pixels.height

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

## Apply the function and unpack the list into a new column
polygon_data <- annotations_uncropped %>%
  rowwise() %>%
  mutate(point.polygon = paste0("[", paste(
    sapply(generate_relative_polygon(radius, Width, Height)$point.polygon,
           function(pt) paste0("[", pt[1], ",", pt[2], "]")),
    collapse = ", "), "]")) %>%
  ungroup()

## subset and write to csv
polygon_dat <- polygon_data[which(polygon_data$SurveyID=="PS06"),]
write_csv(polygon_dat, paste0(counts.dir, "Circumpolar_DownwardImages_ExhaustiveSearch_Annotations_ForSquidle_202507_PS06.csv"))
polygon_dat <- polygon_data[which(polygon_data$SurveyID=="PS14"),]
write_csv(polygon_dat, paste0(counts.dir, "Circumpolar_DownwardImages_ExhaustiveSearch_Annotations_ForSquidle_202507_PS14.csv"))
polygon_dat <- polygon_data[which(polygon_data$SurveyID=="PS18"),]
write_csv(polygon_dat, paste0(counts.dir, "Circumpolar_DownwardImages_ExhaustiveSearch_Annotations_ForSquidle_202507_PS18.csv"))
polygon_dat <- polygon_data[which(polygon_data$SurveyID=="PS61"),]
write_csv(polygon_dat, paste0(counts.dir, "Circumpolar_DownwardImages_ExhaustiveSearch_Annotations_ForSquidle_202507_PS61.csv"))
polygon_data_PS81 <- polygon_data[which(polygon_data$SurveyID=="PS81"),]
write_csv(polygon_data_PS81, paste0(counts.dir, "Circumpolar_DownwardImages_ExhaustiveSearch_Annotations_ForSquidle_202507_PS81.csv"))
polygon_data_PS96 <- polygon_data[which(polygon_data$SurveyID=="PS96"),]
write_csv(polygon_data_PS96, paste0(counts.dir, "Circumpolar_DownwardImages_ExhaustiveSearch_Annotations_ForSquidle_202507_PS96.csv"))
polygon_dat <- polygon_data[which(polygon_data$SurveyID=="PS118"),]
write_csv(polygon_dat, paste0(counts.dir, "Circumpolar_DownwardImages_ExhaustiveSearch_Annotations_ForSquidle_202507_PS118.csv"))

polygon_dat <- polygon_data[which(polygon_data$SurveyID=="TAN0802"),]
write_csv(polygon_dat, paste0(counts.dir, "Circumpolar_DownwardImages_ExhaustiveSearch_Annotations_ForSquidle_202507_TAN0802.csv"))
polygon_dat <- polygon_data[which(polygon_data$SurveyID=="TAN1802"),]
write_csv(polygon_dat, paste0(counts.dir, "Circumpolar_DownwardImages_ExhaustiveSearch_Annotations_ForSquidle_202507_TAN1802.csv"))
polygon_dat <- polygon_data[which(polygon_data$SurveyID=="TAN1901"),]
write_csv(polygon_dat, paste0(counts.dir, "Circumpolar_DownwardImages_ExhaustiveSearch_Annotations_ForSquidle_202507_TAN1901.csv"))

polygon_dat <- polygon_data[which(polygon_data$SurveyID=="JR262"),]
write_csv(polygon_dat, paste0(counts.dir, "Circumpolar_DownwardImages_ExhaustiveSearch_Annotations_ForSquidle_202507_JR262.csv"))
polygon_dat <- polygon_data[which(polygon_data$SurveyID=="JR15005"),]
write_csv(polygon_dat, paste0(counts.dir, "Circumpolar_DownwardImages_ExhaustiveSearch_Annotations_ForSquidle_202507_JR15005.csv"))
polygon_dat <- polygon_data[which(polygon_data$SurveyID=="JR17001"),]
write_csv(polygon_dat, paste0(counts.dir, "Circumpolar_DownwardImages_ExhaustiveSearch_Annotations_ForSquidle_202507_JR17001.csv"))
polygon_dat <- polygon_data[which(polygon_data$SurveyID=="JR17003"),]
write_csv(polygon_dat, paste0(counts.dir, "Circumpolar_DownwardImages_ExhaustiveSearch_Annotations_ForSquidle_202507_JR17003.csv"))

polygon_dat <- polygon_data[which(polygon_data$SurveyID=="NBP0808"),]
write_csv(polygon_dat, paste0(counts.dir, "Circumpolar_DownwardImages_ExhaustiveSearch_Annotations_ForSquidle_202507_NBP0808.csv"))
polygon_dat <- polygon_data[which(polygon_data$SurveyID=="NBP1001"),]
write_csv(polygon_dat, paste0(counts.dir, "Circumpolar_DownwardImages_ExhaustiveSearch_Annotations_ForSquidle_202507_NBP1001.csv"))
polygon_dat <- polygon_data[which(polygon_data$SurveyID=="NBP1402"),]
write_csv(polygon_dat, paste0(counts.dir, "Circumpolar_DownwardImages_ExhaustiveSearch_Annotations_ForSquidle_202507_NBP1402.csv"))
polygon_dat <- polygon_data[which(polygon_data$SurveyID=="NBP1502"),]
write_csv(polygon_dat, paste0(counts.dir, "Circumpolar_DownwardImages_ExhaustiveSearch_Annotations_ForSquidle_202507_NBP1502.csv"))

polygon_dat <- polygon_data[which(polygon_data$SurveyID=="CRS"),]
write_csv(polygon_dat, paste0(counts.dir, "Circumpolar_DownwardImages_ExhaustiveSearch_Annotations_ForSquidle_202507_CRS.csv"))
polygon_dat <- polygon_data[which(polygon_data$SurveyID=="LMG1311"),]
write_csv(polygon_dat, paste0(counts.dir, "Circumpolar_DownwardImages_ExhaustiveSearch_Annotations_ForSquidle_202507_LMG1311.csv"))























#######################################
### OLD CODE THAT CREATE POLYGONS BASED ON PIXELS
# # Function to generate x and y coordinate strings
# generate_polygon_xy_strings <- function(x, y, r, n_points = 16) {
#   angles <- seq(0, 2 * pi, length.out = n_points + 1)[-1]
#   x_coords <- round(x + r * cos(angles), 2)
#   y_coords <- round(y + r * sin(angles), 2)
#   list(
#     point.polygon.x = paste(x_coords, collapse = ";"),
#     point.polygon.y = paste(y_coords, collapse = ";")
#   )
# }
# 
# # Apply the function and unpack the list into two columns
# polygon_data <- annotations %>%
#   rowwise() %>%
#   mutate(coords = list(generate_polygon_xy_strings(point.pixels.x, point.pixels.y, radius))) %>%
#   unnest_wider(coords) %>%
#   ungroup()
# 
# # Write to CSV without row names
# write.csv(polygon_data, file = paste0(counts.dir, "Circumpolar_DownwardImages_ExhaustiveSearch_Annotations_ForSquidle_202507_AA2011.csv"), row.names = FALSE)
