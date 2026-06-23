## Some issues with the colourcorrected images in AA2011.
## Navdata is correct, but sometimes the wrong image selected
## See email from Emma Flukes on 02.06.2025

library(readr)
library(stringr)

## directories
r.dir   <- "R:/IMAS/Antarctic_Seafloor/"
sq.dir  <- paste0(r.dir,"SQUIDLE_dataset_202411/")
dat.dir <- paste0(r.dir,"Clean_Data_For_Permanent_Storage/")
## new folder for fixed AA2011 data
sq.fix.dir  <- paste0(r.dir, "SQUIDLE_dataset_fix_AA2011/")
## folder where original images are stored
img.src.dir <- paste0(dat.dir, "AA2011/AA2011_1_raw_images_and_metadata/images_original/")

'%!in%' <- function(x,y)!('%in%'(x,y))

## list navdata files
nav_files <- list.files(
  path = paste0(sq.dir, "Antarctica_East_2011_AA2011/CTD/"),
  recursive = TRUE,
  pattern = "navdata.csv",
  full.names = TRUE
)

## create subfolders, copy navdata files and original images into these folders
## Process each navdata.csv
for (file_path in nav_files) {
  print(file_path)
  folder_name <- dirname(file_path)
  folder_base <- basename(folder_name)
  dest_folder <- file.path(sq.fix.dir, folder_base)
  
  ## Ensure destination folders exist
  dir.create(dest_folder, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dest_folder, "images"), showWarnings = FALSE)
  dir.create(file.path(dest_folder, "thumbnails"), showWarnings = FALSE)
  
  ## Copy navdata.csv
  file.copy(file_path, file.path(dest_folder, "navdata.csv"), overwrite = TRUE)
  
  ## Read navdata.csv
  navdata <- read_csv(file_path, show_col_types = FALSE)
  
  ## Extract image names
  if ("key" %in% names(navdata)) {
    image_keys <- navdata$key
    
    for (img_key in image_keys) {
      ## Extract subfolder name from image key (e.g., CTD54_CAM02)
      subfolder <- str_extract(img_key, "CTD\\d+_CAM\\d+")
      image_name <- str_extract(img_key, "IMG_\\d+\\.JPG")
      
      if (!is.na(subfolder) && !is.na(image_name)) {
        src_img_path <- file.path(img.src.dir, subfolder, image_name)
        dest_img_path <- file.path(dest_folder, "images", img_key)  # Rename to full key
        
        if (file.exists(src_img_path)) {
          file.copy(src_img_path, dest_img_path, overwrite = TRUE)
        } else {
          warning(paste("Image not found:", src_img_path))
        }
      }
    }
  } else {
    warning(paste("No 'key' column in:", file_path))
  }
}

## colour correct images IN PHOTOSHOP

## create thumbnails
library(magick)
## Function to resize images for thumbnails
resize_image <- function(image_path, output_path, size = "150x150") {
  image <- image_read(image_path)
  thumbnail <- image_resize(image, size)
  image_destroy(image)  # Free memory
  image_write(thumbnail, output_path)
}

## Loop through each destination folder and create thumbnails
for (file_path in nav_files) {
  print(file_path)
  folder_base <- basename(dirname(file_path))
  t.dir <- file.path(sq.fix.dir, folder_base)
  
  ## Just to prevent errors we only include .jpg or .jpeg files
  input_images <- list.files(
    file.path(t.dir, "images"),
    pattern = "\\.(jpg|jpeg)$", ignore.case = TRUE,
    full.names = TRUE
  )
  
  output_dir <- file.path(t.dir, "thumbnails")
  
  for (image_path in input_images) {
    output_path <- file.path(output_dir, basename(image_path))
    tryCatch({
      resize_image(image_path, output_path)
    }, error = function(e) {
      message("Failed to process: ", image_path)
      message("Reason: ", e$message)
    })
  }
}






