###### SETUP FOLDER STRUCTURE FOR SQUIDLE+ ######
library(stringr)
library(lubridate)

## directories
r.dir   <- "R:/IMAS/Antarctic_Seafloor/"
sq.dir  <- paste0(r.dir,"SQUIDLE_dataset_202411/")
dat.dir <- paste0(r.dir,"ASAID_data/")

'%!in%' <- function(x,y)!('%in%'(x,y))

## name each campaign (copy folder names from existing setup)
## and add the gear into each campaign
campaign.names <- list.dirs(paste0(r.dir,"SQUIDLE_dataset_202411"), recursive=FALSE, full.names=FALSE)
gear.names <- c("CTD","YOYO","FTS",rep("YOYO",3),"SUCS","YOYO","OFOS",rep("SUCS",3),"OFOBS","DTIS","YOYO",rep("DTIS",2),rep("FTS",3),"OFOS")
# for(i in 1:length(campaign.names)){
#   dir.create(paste0(sq.dir, campaign.names[i]))
#   dir.create(paste0(sq.dir, campaign.names[i],"/",gear.names[i]))
# }

####
## read in annotation files (ann.count.dat & ann.cover.dat)
#ann.pts <- read.csv(paste0(dat.dir,"AnnotationLibrary_AllFinishedSurveys/Cover/Circumpolar_DownwardImages_PointScore_Annotations_202312_rawforpublication.csv"))
# load(paste0(dat.dir,"AnnotationLibrary_AllFinishedSurveys/Cover/Circumpolar_DownwardImages_PointScore_Annotations_202312.Rdata"))
# load(paste0(dat.dir,"AnnotationLibrary_AllFinishedSurveys/Counts/Circumpolar_DownwardImages_ExhaustiveSearch_Annotations_202312.Rdata"))
# ann.vme <- "..."

## read in image metadata (annotated and not annotated)
full.metadata.list <- list()
full.metadata.paths <- list.files(dat.dir, pattern="metadata_full", full.names=TRUE)
for(i in 1:length(full.metadata.paths)){
  full.metadata.list[[i]] <- read.csv(full.metadata.paths[i])
}
ann.metadata.list <- list()
ann.metadata.paths <- list.files(dat.dir, pattern="metadata_annotated", full.names=TRUE)
for(i in 1:length(ann.metadata.paths)){
  ann.metadata.list[[i]] <- read.csv(ann.metadata.paths[i])
}
## we need to fix a few images that are in the metadata, but don't exist:
full.metadata.list[[17]] <- full.metadata.list[[17]][-which(full.metadata.list[[17]]$Filename=="PS81_188-1_2013-02-20T09_59_19_0076.jpg"),]
## PS118_39_9995__TIMER_2019_03_23_at_23_53_37_IMG_0980.jpg
## AND THESE:
## PS81_160-1_2013-02-08T17_26_39_6493.jpg
## PS81_160_0152__PS81_160-1_2013-02-08T17_27_09_6494.jpg
## PS81_160_0193__PS81_160-1_2013-02-08T20_06_25_6842.jpg
## PS81_160_0349__PS81_160-1_2013-02-08T20_06_55_6843.jpg
## PS81_215-2_2013-03-01T23_33_38_4610.jpg
## PS81_215-2_2013-03-01T23_34_24_4611.jpg

#### campaign/gear/transect/
## deployment csv file in each gear, showing transect start lon/lat (MAYBE NOT NEEDED)

#### Setup transect IDs for folders (first part of filename) AND
#### Fix some dates AND
#### Need to fix some image area estimates, because they were wrongly based on the cropped annotated images:
## create a function to automatically read transect IDs
#Function to subset string until just before the second underscore
subset_string <- function(text) {
  #Find the positions of all underscores
  underscore_positions <- str_locate_all(text, "_")[[1]]
    #Check if there are at least two underscores
  if (nrow(underscore_positions) >= 2) {
    second_underscore_pos <- underscore_positions[2, 1]
    #Subset the string from the beginning to just before the second underscore
    return(str_sub(text, 1, second_underscore_pos - 1))
  } else {
    #If there are less than two underscores, return the original string
    return(text)
  }
  #  return(str_sub(text, 1, underscore_positions[2,1] - 1))
}
for(i in 1:length(campaign.names)){
  message(i)
  dat <- full.metadata.list[[i]]
  ## standardised transect ID names:
  dat$Transect.ID.standardised <- unlist(lapply(dat$Filename.standardised, subset_string))
  t.IDs <- unique(dat$Transect.ID.standardised)
  #print(t.IDs)
  s.ID <- unique(dat$Survey.ID)
  
  #### fix date format
  if(s.ID %in% c("TAN0802", "TAN1802","TAN1901")){
    dat$timestamp <- format(ymd_hms(dat$Date, truncated=3)-hours(12), "%Y-%m-%d %H:%M:%S UTC")
    print(head(dat$timestamp))
  }else if(i %in% c(2,6,8,9,11)){ ## surveys with Date only
    dat$timestamp <- format(ymd(dat$Date), "%Y-%m-%d %H:%M:%S UTC")
    print(head(dat$timestamp))
  }else if(i %in% c(12,14:16)){ ## surveys with start/end of transects only
    dat$timestamp <- NA
    for(k in 1:length(t.IDs)){
      t.rows <- which(dat$Transect.ID.standardised==t.IDs[k])
      if(t.IDs[k]=="PS18_220"){
        int.start <- ymd_hms("1991-03-13 00:00:00")
        int.end   <- ymd_hms("1991-03-13 15:20:00")
      }else{
        int.start <- ymd_hms(dat$Date_start[t.rows], truncated=3)[1]
        int.end   <- ymd_hms(dat$Date_end[t.rows], truncated=3)[1]
      }
      int.values <- approx(x=1:nrow(dat[t.rows,]), y=seq(int.start, int.end, length.out=nrow(dat[t.rows,])), xout=1:nrow(dat[t.rows,]))$y
      int.times <- as.POSIXct(int.values, tz="UTC")
      int.times.date <- format(ymd_hms(int.times), "%Y-%m-%d %H:%M:%S UTC")
      if(is.na(int.times.date[1])){
        int.times.date[1] <- format(ymd(int.start), "%Y-%m-%d %H:%M:%S UTC")
        message("NA in date fixed")
        }
      dat$timestamp[t.rows] <- int.times.date
    }
    print(head(dat$timestamp))
  }else{
    dat$timestamp <- format(ymd_hms(dat$Date, truncated=3), "%Y-%m-%d %H:%M:%S UTC")
    print(head(dat$timestamp))
  }
  #### fix area estimates in PS81, PS96, PS118 and AA2011 (averages had been calculated on the cropped images, so actual area is bigger)
  ## PS81 and 96 areas taken from metadata:
  if(s.ID %in% c("PS81", "PS96", "PS118", "AA2011")){
    ## which observations need to be re-calculated
    ave.t.sel <- which(dat$Image_area_source=="averaged_across_transect")
    ## which transects do they belong to
    t.sel <- dat$Transect.ID[ave.t.sel]
    ## how many transects are there
    t.length <- length(unique(t.sel))
    ## PS81 and PS96 have existing area information from metadata to fill the blanks
    if(s.ID %in% c("PS81", "PS96")){
      if(t.length>1){
        for(k in 1:t.length){
          ## the transect
          t.ID.now <- unique(t.sel)[k]
          ## replace those observations with the mean from the metadata at those transects
          dat$Image_area[ave.t.sel][t.sel==t.ID.now] <- mean(dat$Image_area[dat$Image_area_source=="metadata" & dat$Transect.ID==t.ID.now])
        }
      }else{
        t.ID.now <- unique(t.sel)
        dat$Image_area[ave.t.sel][t.sel==unique(t.sel)] <- mean(dat$Image_area[dat$Image_area_source=="metadata" & dat$Transect.ID==t.ID.now])
      }
    } ## PS118 only has laser points to fill the blanks (need to add the cropped area back in)
    if(s.ID == "PS118"){
      if(t.length>1){
        for(k in 1:t.length){
          t.ID.now <- unique(t.sel)[k]
          ## replace those observations with the laser points but add back the cropped area
          dat$Image_area[ave.t.sel][t.sel==t.ID.now] <- mean(dat$Image_area[dat$Image_area_source=="laser_points" & dat$Transect.ID==t.ID.now])/0.64
        }
      }else{
        t.ID.now <- unique(t.sel)
        dat$Image_area[ave.t.sel][t.sel==unique(t.sel)] <- mean(dat$Image_area[dat$Image_area_source=="laser_points" & dat$Transect.ID==t.ID.now])/0.64
      }
    } ## AA2011 has most images without data from transects
    if(s.ID == "AA2011"){
      if(t.length>1){
        for(k in 1:t.length){
          t.ID.now <- unique(t.sel)[k]
          ## replace those observations with the laser points but add back the cropped area (unless no laser points!!!, then use survey average)
          lasers <- any(dat$Image_area_source=="laser_points" & dat$Transect.ID==t.ID.now)
          if(lasers){
            dat$Image_area[ave.t.sel][t.sel==t.ID.now] <- mean(dat$Image_area[dat$Image_area_source=="laser_points" & dat$Transect.ID==t.ID.now])/0.81
          }else {
            dat$Image_area[ave.t.sel][t.sel==t.ID.now] <- mean(dat$Image_area[dat$Image_area_source=="laser_points"])/0.81
            dat$Image_area_source[ave.t.sel][t.sel==t.ID.now] <- "averaged_across_survey"
          }
        }
      }else{
        t.ID.now <- unique(t.sel)
        dat$Image_area[ave.t.sel][t.sel==unique(t.sel)] <- mean(dat$Image_area[dat$Image_area_source=="laser_points" & dat$Transect.ID==t.ID.now])/0.81
      }
      ## 
      ave.s.sel <- which(dat$Image_area_source=="averaged_across_survey")
      dat$Image_area[ave.s.sel] <- mean(dat$Image_area[dat$Image_area_source=="laser_points"])/0.81
    }
  }
  ## fix filenames in PS96 and PS118
  if(s.ID %in% c("PS118", "PS96")){
    dat$Filename <- paste0(dat$Filename,".jpg")
  }
  full.metadata.list[[i]] <- dat
}


#### to fix TAN1802 and TAN1901 depths for some images (they are -999.9)
for(i in 16:17){
  message(i)
  dat <- full.metadata.list[[i]]
  ## standardised transect ID names:
  dat$Transect.ID.standardised <- unlist(lapply(dat$Filename.standardised, subset_string))
  t.IDs <- unique(dat$Transect.ID.standardised)
  #print(t.IDs)
  s.ID <- unique(dat$Survey.ID)
  #### fix date format
  if(s.ID %in% c("TAN0802", "TAN1802","TAN1901")){
    dat$timestamp <- format(ymd_hms(dat$Date, truncated=3)-hours(12), "%Y-%m-%d %H:%M:%S UTC")
    print(head(dat$timestamp))
  }
  full.metadata.list[[i]] <- dat
}




########
library(magick)
# Function to resize images for thumbnails
resize_image <- function(image_path, output_path, size = "150x150") {
  image <- image_read(image_path)
  thumbnail <- image_resize(image, size)
  image_destroy(image) ## this is needed because otherwise it runs out of temporary memory when processing many images
  image_write(thumbnail, output_path)
}

## now add folders, images and navdata file
for(i in 1:length(campaign.names)){ #
  message(i)
  dat <- full.metadata.list[[i]]
  s.ID <- unique(dat$Survey.ID) 
  
  ## we need to reduce the navdata files to annotated images only for "CRS","NBP0808" and "NBP1001"
  if(s.ID %in% c("CRS", "NBP0808", "NBP1001")) {
    # Get the corresponding annotated metadata for this campaign
    ann_dat <- ann.metadata.list[[i]]
    ## Filter dat to only keep rows where Filename.standardised exists in ann_dat
    dat <- dat[dat$Filename.standardised %in% ann_dat$Filename.standardised, ]
  }
  
  #### adding main folders
  t.IDs <- unique(dat$Transect.ID.standardised)
  ## need to match the survey IDs between the two dataframes first
  sel <- which(str_extract(campaign.names, "(?<=_)[^_]*$")==unique(dat$Survey.ID))
  
  #### now add images, thumbnails and navdata file to each folder
  survey.path <- paste0(dat.dir,unique(dat$Survey.ID))
  
  if(s.ID %in% c("AA2011","JR15005","JR17001","JR17003","JR262","LMG1311","NBP1402","NBP1502","PS06","PS118","PS14","PS18","PS61","PS96","TAN0802","TAN1802","TAN1901")){
    ## where to look for annotated images
    folder1 <- list.files(survey.path, full.names=TRUE)[3]
    ## where to look for other images
    folder2 <- paste0(list.files(survey.path, full.names=TRUE)[1],"/images_colourcorrected")
    oth.files.full <- list.files(folder2, recursive=TRUE, full.names=TRUE)
  }
  
  ## loop across transects
  for(k in 1:length(t.IDs)){
    t.ID.now <- t.IDs[k]
    print(paste0(k," of ",length(t.IDs),": ",t.ID.now))
    t.dir <- paste0(sq.dir, campaign.names[sel],"/",gear.names[sel],"/",t.ID.now)
    # dir.create(t.dir)
    # dir.create(paste0(t.dir, "/images"))
    # dir.create(paste0(t.dir, "/thumbnails"))
    
    ##subset data to transects only
    dat.t <- dat[dat$Transect.ID.standardised==t.ID.now,]
    # ## navdata file:
    navdata <- data.frame(key       = dat.t$Filename.standardised,
                          pose.lon  = dat.t$Longitude,
                          pose.lat  = dat.t$Latitude,
                          pose.dep  = dat.t$Depth,
                          pose.alt = dat.t$HeightAboveSeafloor,
                          timestamp_start = dat.t$timestamp,
                          #data.survey_ID         = dat.t$Survey.ID,
                          data.transect_ID       = dat.t$Transect.ID.standardised,
                          data.image_area        = dat.t$Image_area,
                          data.image_area_source = dat.t$Image_area_source,
                          data.license     = dat.t$License,
                          data.source_link = dat.t$Source.link)
    ## Exclude columns that are entirely NA
    navdata <- Filter(function(x) !all(is.na(x)), navdata)
    #if(!is.null(dat.t$timestamp.end)) navdata$timestamp.end <- dat.t$timestamp.end
    write.csv(navdata, file=paste0(t.dir, "/navdata.csv"), row.names=FALSE)
    
    ## fill images, first use the ones that are colour-corrected and cropped, then fill with all other ones
    #for the CRS images we only upload the annotated ones
    #for PS81, we have two folders of colourcorrected images (PS81 and PS81_shallow) for transects 185,186,189
    if(unique(dat$Survey.ID) == "PS81" & t.ID.now %in% c("PS81_185","PS81_186","PS81_189")){
      ## where to look for images
      folder1 <- paste0(dat.dir,"PS81/PS81_3_cropped_and_colourcorrected_images_for_annotation/")
      folder2 <- paste0(dat.dir,"PS81_shallow/PS81_shallow_3_cropped_and_colourcorrected_images_for_annotation/")
      folder3 <- paste0(dat.dir,"PS81/PS81_1_raw_images_and_metadata/images_original/")
      ## which images to copy
      f1.sel <- which(dat.t$Filename.standardised %in% list.files(folder1))
      f2.sel <- which(dat.t$Filename.standardised %in% list.files(folder2))
      ## define names and directories
      ann.images1 <- paste0(folder1, dat.t$Filename.standardised[f1.sel])
      ann.images2 <- paste0(folder2, dat.t$Filename.standardised[f2.sel])
      ann.images3 <- paste0(folder3, dat.t$Filename[-c(f1.sel,f2.sel)])
      ann.images3.target <- paste0(t.dir, "/images/", dat.t$Filename.standardised[-c(f1.sel,f2.sel)])
      ## copy files
      file.copy(ann.images2, paste0(t.dir, "/images"))
      file.copy(ann.images3, ann.images3.target)
    }else if(unique(dat$Survey.ID) == "PS81" & t.ID.now %!in% c("PS81_185","PS81_186","PS81_189")){
      ## where to look for images
      folder1 <- paste0(dat.dir,"PS81/PS81_3_cropped_and_colourcorrected_images_for_annotation/")
      folder3 <- paste0(dat.dir,"PS81/PS81_1_raw_images_and_metadata/images_original/")
      ## which images to copy
      f1.sel <- which(dat.t$Filename.standardised %in% list.files(folder1))
      ## define names and directories
      ann.images1 <- paste0(folder1, dat.t$Filename.standardised[f1.sel])
      ann.images3 <- paste0(folder3, dat.t$Filename[-f1.sel])
      ann.images3.target <- paste0(t.dir, "/images/", dat.t$Filename.standardised[-f1.sel])
      ## copy files
      file.copy(ann.images3, ann.images3.target)
    }else if(unique(dat$Survey.ID) %in% c("CRS","NBP0808", "NBP1001")){
      ## where to look for images
      folder1 <- paste0(dat.dir,"CRS/CRS_3_colourcorrected_images_for_annotation/")
      ## which images to copy
      f1.sel <- which(dat.t$Filename.standardised %in% list.files(folder1))
      ## define names and directories
      ann.images1 <- paste0(folder1, dat.t$Filename.standardised[f1.sel])
    }else{
      # ## where to look for annotated images
      # folder1 <- list.files(survey.path, full.names=TRUE)[3]
      # ## where to look for other images
      # folder2 <- paste0(list.files(survey.path, full.names=TRUE)[1],"/images_colourcorrected")
      # oth.files.full <- list.files(folder2, recursive=TRUE, full.names=TRUE)
      if(length(oth.files.full)==0){
        folder2 <- paste0(list.files(survey.path, full.names=TRUE)[1],"/images_original")
        oth.files.full <- list.files(folder2, recursive=TRUE, full.names=TRUE)
      }
      ## which images to copy
      f1.sel <- which(dat.t$Filename.standardised %in% list.files(folder1))
      oth.files.names <- basename(oth.files.full) ## remove directories from names
      f2.sel <- which(dat.t$Filename %in% oth.files.names)
      f2.sel.rev <- which(oth.files.names %in% dat.t$Filename)
      ## define names and directories
      ann.images1 <- paste0(folder1, "/", dat.t$Filename.standardised[f1.sel])
      ann.images2 <- oth.files.full[f2.sel.rev]
      img.match <- match(sub(".*/", "", ann.images2),dat.t$Filename)
      ann.images2.target <- paste0(t.dir, "/images/", dat.t$Filename.standardised[img.match])
      ## copy files
      file.copy(ann.images2, ann.images2.target)
    }
    file.copy(ann.images1, paste0(t.dir, "/images"))
    
    ## create thumbnails from images
    input <- list.files(paste0(t.dir, "/images"), full.names=TRUE)
    output_dir <- paste0(t.dir, "/thumbnails")
    for(image_path in input) {
      output_path <- file.path(output_dir, basename(image_path))
      resize_image(image_path, output_path)
    }
  }
  
}

## NEED TO MANUALLY REMOVE THIS ENTRY FROM PS118 navdata.csv:
## PS118_39_9995__TIMER_2019_03_23_at_23_53_37_IMG_0980.jpg

## AND THESE:
## PS81_160-1_2013-02-08T17_26_39_6493.jpg
## PS81_160_0152__PS81_160-1_2013-02-08T17_27_09_6494.jpg
## PS81_160_0193__PS81_160-1_2013-02-08T20_06_25_6842.jpg
## PS81_160_0349__PS81_160-1_2013-02-08T20_06_55_6843.jpg
## PS81_215-2_2013-03-01T23_33_38_4610.jpg
## PS81_215-2_2013-03-01T23_34_24_4611.jpg


#### Edit: change TAN1802 and TAN1901 depths for some images (they are -999.9)












## change names to match what Squidle expects




##
for(i in 19){
  message(i)
  dat <- full.metadata.list[[i]]
  ##
  s.ID <- unique(dat$Survey.ID)
  t.IDs <- unique(dat$Transect.ID.standardised)
  print(s.ID)
  ## match the survey IDs between the two dataframes
  sel <- which(str_extract(campaign.names, "(?<=_)[^_]*$")==s.ID)
  ##
  for(k in c(19,20,22:26,30:34)){
    t.ID.now <- t.IDs[k]
    print(paste0(k," of ",length(t.IDs),": ",t.ID.now))
    ##subset data to transects only
    dat.t <- dat[dat$Transect.ID.standardised==t.ID.now,]
    print(any(dat.t$Longitude>180))
  ##
  }}

## transect names
tr.names <- c("tan0802_183","tan0802_186","tan0802_200","tan0802_202","tan0802_205","tan0802_207",
              "tan0802_214","tan0802_228","tan0802_239","tan0802_244","tan0802_246","tan0802_248")
full.metadata.list[[19]]$Longitude[full.metadata.list[[19]]$Transect.ID.standardised%in%tr.names]



