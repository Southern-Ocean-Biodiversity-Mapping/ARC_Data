###################################################################################################
### this script updates the labels in the csv annotation files to publishable labels            ###
###################################################################################################

# 1) libraries and paths 
library(tidyverse)
library(writexl)
library(readxl)
library(plyr)
'%!in%' <- function(x,y)!('%in%'(x,y))

user = "Jan"
#user = "charley"
#user="nicole"

if (user == "Jan") {
  sci.dir <-      "C:/Users/jjansen/OneDrive - University of Tasmania/science/"
  env.derived <-  paste0(sci.dir,"data_environmental/derived/")
  tools.dir <-    paste0(sci.dir,"SouthernOceanBiodiversityMapping/Useful_Functions_Tools/")
  ARC_Data.dir <- paste0(sci.dir,"SouthernOceanBiodiversityMapping/ARC_Data/")
  r.path <- "R:/IMAS/Antarctic_Seafloor/Clean_Data_For_Permanent_Storage/"
  cover.path <- paste0(r.path,"AnnotationLibrary_AllFinishedSurveys/Cover/")
  count.path <- paste0(r.path,"AnnotationLibrary_AllFinishedSurveys/Counts/")
  counts.surveys.path <- paste0(r.path,"AnnotationLibrary_AllFinishedSurveys/Counts_annotation_files_individual_surveys/")
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

###################################

## read in excel sheets
cover_list<-read_xlsx(path=paste0(ARC_Data.dir, "Annotation/Species_list_2023_01.xlsx"),sheet=1)
count_list<-read_xlsx(path=paste0(ARC_Data.dir, "Annotation/Species_list_2023_01.xlsx"),sheet=2)
## read in CATAMI reference sheets
cover_list2<-read_xlsx(path=paste0(ARC_Data.dir, "Annotation/Species_list_vs_CATAMI_2025_01.xlsx"),sheet=1)
count_list2<-read_xlsx(path=paste0(ARC_Data.dir, "Annotation/Species_list_vs_CATAMI_2025_01.xlsx"),sheet=2)

## copied from "AnnotationQualityControl_CoralNet_AllSurveys_03_final_library.R":
## specify survey and folder names
survey.IDs <- c("AA2011","CRS","JR262","JR15005","JR17001","JR17003","LMG1311","NBP1402","NBP1502",
                "PS06","PS14","PS18","PS61","PS81","PS81_shallow","PS96","PS118",
                "TAN0802","TAN1802","TAN1901")
folder.IDs <- c("AA2011_3_colourcorrected_images_for_annotation"
                ,"CRS_3_colourcorrected_images_for_annotation"
                ,"JR262_3_cropped_and_colourcorrected_images_for_annotation"
                ,"JR15005_3_cropped_and_colourcorrected_images_for_annotation"
                ,"JR17001_3_cropped_and_colourcorrected_images_for_annotation"
                ,"JR17003_3_cropped_and_colourcorrected_images_for_annotation"
                ,"LMG1311_3_colourcorrected_images_for_annotation"
                ,"NBP1402_3_colourcorrected_images_for_annotation" 
                ,"NBP1502_3_colourcorrected_images_for_annotation" 
                ,"PS06_3_colourcorrected_images_for_annotation"
                ,"PS14_3_colourcorrected_images_for_annotation"
                ,"PS18_3_colourcorrected_images_for_annotation"
                ,"PS61_3_colourcorrected_images_for_annotation"
                ,"PS81_3_cropped_and_colourcorrected_images_for_annotation"
                ,"PS81_shallow_3_cropped_and_colourcorrected_images_for_annotation"
                ,"PS96_3_cropped_images_for_annotation_NoColourCorrectionNeeded"
                ,"PS118_3_cropped_and_colourcorrected_images_for_annotation"
                ,"TAN0802_3_colourcorrected_images_for_annotation"
                ,"TAN1802_3_colourcorrected_images_for_annotation"
                ,"TAN1901_3_colourcorrected_images_for_annotation")
## path to full images
img.path.raw <- paste0(r.path,survey.IDs,"/",folder.IDs)
load(paste0(ARC_Data.dir,"Image_level_bio_202501.Rdata"))

##############################

##### Point-score data #####

## read in annotation file
ann.cover.dat.raw <- read.csv(paste0(cover.path,"Circumpolar_DownwardImages_PointScore_Annotations_20220307.csv"))
ann.cover.dat <- ann.cover.dat.raw[,c(1:3,5:7,9,4,10)]
ann.cover.dat$Label_after_review_clean <- NA
  
## clean up labels
ann.cover.dat.temp <- ann.cover.dat %>%
  left_join(cover_list[,c("Label", "Merge_With")], by=c("Label_after_review"="Label"))
ann.cover.dat.temp$Label_after_review_clean <- ifelse(!is.na(ann.cover.dat.temp$Merge_With), ann.cover.dat.temp$Merge_With, ann.cover.dat.temp$Label_after_review)
ann.cover.dat2 <- ann.cover.dat.temp[,1:10]

## add publishable labels, CATAMI and AMC labels
ann.cover.dat.temp <- ann.cover.dat2 %>%
  left_join(cover_list2[,c("Label", "Name_to_publish", "AMC", "AMC_ID", "CATAMI_broad", "CATAMI")], by=c("Label_after_review_clean"="Label"))
ann.cover.dat.temp$CATAMI_broad[ann.cover.dat.temp$CATAMI_broad=="NA"] <- ""
ann.cover.dat.temp$CATAMI[is.na(ann.cover.dat.temp$CATAMI)] <- ""
# ann.cover.dat.temp$Name_to_publish <- ifelse(!is.na(ann.cover.dat.temp$Name_to_publish), ann.cover.dat.temp$Name_to_publish, ann.cover.dat.temp$Label_after_review_clean)
# ann.cover.dat.temp$CATAMI_broad <- ifelse(!is.na(ann.cover.dat.temp$CATAMI_broad), ann.cover.dat.temp$CATAMI_broad, ann.cover.dat.temp$Label_after_review_clean)
# ann.cover.dat.temp$CATAMI_fine <- ifelse(!is.na(ann.cover.dat.temp$CATAMI), ann.cover.dat.temp$CATAMI, ann.cover.dat.temp$Label_after_review_clean)
ann.cover.dat.temp$CATAMI <- paste0(ann.cover.dat.temp$CATAMI_broad," ",ann.cover.dat.temp$CATAMI)
#ann.cover.dat3 <- ann.cover.dat.temp#[,c(1:11,13,14)]

## add survey names:
ann.cover.dat3 <- ann.cover.dat.temp %>% left_join(img.metadata[,c("Filename.standardised", "survey")], by=c("Name"="Filename.standardised"))
names(ann.cover.dat3)[16] <- "SurveyID"
ann.cover.dat3$SurveyID <- sub("tan","TAN",ann.cover.dat3$SurveyID)

write.csv(ann.cover.dat3,file=paste0(cover.path,"Circumpolar_DownwardImages_PointScore_Annotations_202501_rawforpublication.csv"))

## for building the classification catalog, add the folder path to each image
ann.cover.dat4 <- ann.cover.dat3
## create lookup-table to find full images in the folders
ann.cover.dat4$TransectID <- sub("(^[^_]+[_][^_]+)(.+$)","\\1",ann.cover.dat4$Name)
survey_lookup <- sub("_.*","",ann.cover.dat4$Name)
survey_lookup <- sub("tan","TAN",survey_lookup)

ann.cover.dat4$Folder.path <- NA
for(i in 1:length(survey.IDs)){
  sel.temp <- which(survey_lookup==survey.IDs[i])
  ann.cover.dat4$Folder.path[sel.temp] <- paste0(r.path,survey.IDs[i],"/",folder.IDs[i],"/")
}
## PS81_shallow is annoying...:
files.PS81shallow <- list.files(paste0(r.path,survey.IDs[15],"/",folder.IDs[15],"/"), pattern=".jpg")
ann.cover.dat4$Folder.path[which(ann.cover.dat4$Name%in%files.PS81shallow)] <- paste0(r.path,survey.IDs[15],"/",folder.IDs[15],"/")
#head(ann.cover.dat4)

## add CAAB data and make character vectors a factor
# options(scipen = 999) ## we don't want scientific abbreviations in the CAAB code
# ann.cover.dat4$CAAB <- paste0("CAAB ", ann.cover.dat4$CAAB)
# factor.cols <- c("Annotator","Label","Label_after_review","Label_after_review_clean","Name_to_publish","AMC","AMC_ID","CATAMI","SurveyID")
# ann.cover.dat4[,factor.cols] <- lapply(ann.cover.dat4[,factor.cols], factor)
ann.cover.dat <- ann.cover.dat4

save(ann.cover.dat, file=paste0(cover.path,"Circumpolar_DownwardImages_PointScore_Annotations_202501.Rdata"))


##### Exhaustive-search data #####

## load the data
ann.files.csv <- list.files(counts.surveys.path,pattern="species.csv", recursive=TRUE, full.names=TRUE)
ann.count.dat.raw <- ldply(ann.files.csv, read_csv)

## lookup table to find images to generate crops from
ann.count.dat.raw$TransectID <- as.factor(sub("(^[^_]+[_][^_]+)(.+$)","\\1",ann.count.dat.raw$filename))
ann.count.dat.raw2 <- ann.count.dat.raw %>% left_join(img.metadata[,c("Filename.standardised", "survey")], by=c("filename"="Filename.standardised"))
names(ann.count.dat.raw2)[19] <- "SurveyID"
ann.count.dat.raw2$SurveyID <- sub("tan","TAN",ann.count.dat.raw2$SurveyID)
## somehow survey doesn't get copied for transect 163 in P81
sel163 <- which(is.na(ann.count.dat.raw2$SurveyID))
ann.count.dat.raw2$TransectID[sel163]
ann.count.dat.raw2$SurveyID[sel163] <- "PS81"
##
Survey_lookup2 <- sub("_.*","",ann.count.dat.raw$filename)
Survey_lookup2 <- sub("tan","TAN",Survey_lookup2)
#unique(ann.count.dat.raw$SurveyID)
ann.count.dat.raw2$Folder.path <- NA
for(i in 1:length(survey.IDs)){
  sel.temp <- which(Survey_lookup2==survey.IDs[i])
  ann.count.dat.raw2$Folder.path[sel.temp] <- paste0(r.path,survey.IDs[i],"/",folder.IDs[i],"/")
}
## PS81_shallow is annoying...:
files.PS81shallow <- list.files(paste0(r.path,survey.IDs[15],"/",folder.IDs[15],"/"), pattern=".jpg")
ann.count.dat.raw2$Folder.path[which(ann.count.dat.raw2$filename%in%files.PS81shallow)] <- paste0(r.path,survey.IDs[15],"/",folder.IDs[15],"/")

## change character-string of points into a readable format
ann.count.dat.raw2$points  <- as.character(ann.count.dat.raw2$points)
ann.count.dat.raw2$points <- gsub("\\[","",ann.count.dat.raw2$points)
ann.count.dat.raw2$points <- gsub("\\]","",ann.count.dat.raw2$points)
pointvector <- as.numeric(unlist(strsplit(ann.count.dat.raw2$points,",")))
ann.count.dat.raw2$x      <- pointvector[seq(1,length(pointvector),3)]
ann.count.dat.raw2$y      <- pointvector[seq(2,length(pointvector),3)]
ann.count.dat.raw2$radius <- pointvector[seq(3,length(pointvector),3)]

ann.count.dat <- ann.count.dat.raw2[names(ann.count.dat.raw2)%in%c("label_hierarchy","filename","image_longitude","image_latitude","shape_name","x","y","radius","Folder.path", "lastname","TransectID","SurveyID")]
names(ann.count.dat)[1] <- "Label"
names(ann.count.dat)[2] <- "Annotator"
names(ann.count.dat)[3] <- "Name"
ann.count.dat$Label <- gsub(" ","",  ann.count.dat$Label)
ann.count.dat$Label <- gsub(">","__",ann.count.dat$Label)
ann.count.dat$Label <- gsub("-","_",ann.count.dat$Label)
ann.count.dat$Label         <- factor(ann.count.dat$Label)
ann.count.dat$Name           <- factor(ann.count.dat$Name)
ann.count.dat$shape_name     <- factor(ann.count.dat$shape_name)

## clean up labels
ann.count.dat.temp <- ann.count.dat %>%
  left_join(count_list[,c("Label", "Merge_With")], by=c("Label"="Label"))
ann.count.dat.temp$Label_clean <- ifelse(!is.na(ann.count.dat.temp$Merge_With), ann.count.dat.temp$Merge_With, ann.count.dat.temp$Label)
ann.count.dat2 <- ann.count.dat.temp[,c(3,10:12,6,2,1,14,9,7,8)]

## add publishable labels, CATAMI labels and CAAB code
ann.count.dat.temp <- ann.count.dat2 %>%
  left_join(count_list2[,c("Label", "Name_to_publish","AMC","AMC_ID", "CATAMI_broad", "CATAMI")], by=c("Label_clean"="Label"))
ann.count.dat.temp$CATAMI_broad[ann.count.dat.temp$CATAMI_broad=="NA"] <- ""
ann.count.dat.temp$CATAMI[is.na(ann.count.dat.temp$CATAMI)] <- ""
# ann.count.dat.temp$Name_to_publish <- ifelse(!is.na(ann.count.dat.temp$Name_to_publish), ann.count.dat.temp$Name_to_publish, ann.count.dat.temp$Label_after_review_clean)
# ann.count.dat.temp$CATAMI_broad <- ifelse(!is.na(ann.count.dat.temp$CATAMI_broad), ann.count.dat.temp$CATAMI_broad, ann.count.dat.temp$Label_after_review_clean)
# ann.count.dat.temp$CATAMI_fine <- ifelse(!is.na(ann.count.dat.temp$CATAMI), ann.count.dat.temp$CATAMI, ann.count.dat.temp$Label_after_review_clean)
ann.count.dat.temp$CATAMI <- paste0(ann.count.dat.temp$CATAMI_broad," ",ann.count.dat.temp$CATAMI)
ann.count.dat3 <- ann.count.dat.temp#[,c(1:12,14,15)]

ann.count.dat <- ann.count.dat3
save(ann.count.dat, file=paste0(count.path,"Circumpolar_DownwardImages_ExhaustiveSearch_Annotations_202501.Rdata"))

write.csv(ann.count.dat[,-9],file=paste0(count.path,"Circumpolar_DownwardImages_ExhaustiveSearch_Annotations_202501_rawforpublication.csv"))


####################################################################################

## to check annotation numbers for the data paper:

## load annotation file
load(paste0(cover.path,"Circumpolar_DownwardImages_PointScore_Annotations_202312.Rdata"))
load(paste0(count.path,"Circumpolar_DownwardImages_ExhaustiveSearch_Annotations_202412.Rdata"))


seamount_transects <- c("TAN1802_160","TAN1802_170","TAN1802_179","TAN1802_180","TAN1802_183",
                        "TAN1802_184","TAN1802_185","TAN1802_191","TAN1802_193","TAN1802_195",
                        "TAN1802_196","TAN1802_197","TAN1802_207","TAN1802_208","TAN1802_209",
                        "TAN1802_213","tan1901_209")

seamounts.remove <- which(ann.cover.dat$TransectID%in%seamount_transects)
table(ann.cover.dat$SurveyID[-seamounts.remove])
sum(table(ann.cover.dat$SurveyID[-seamounts.remove]))
length(table(ann.cover.dat$Name[-seamounts.remove]))

seamounts.remove2 <- which(ann.count.dat$TransectID%in%seamount_transects)
table(ann.count.dat$SurveyID[-seamounts.remove2])
sum(table(ann.count.dat$SurveyID[-seamounts.remove2]))


########################################################################################
#### Prep annotation data for upload
#### Add image width and height to each dataframe
#### Add transectIDs to cover data

count.dat.raw <- read.csv(paste0(count.path,"Circumpolar_DownwardImages_ExhaustiveSearch_Annotations_202501_rawforpublication.csv"))
cover.dat.raw <- read.csv(paste0(cover.path,"Circumpolar_DownwardImages_PointScore_Annotations_202501_rawforpublication.csv"))
count.dat.raw$TransectID <- as.factor(count.dat.raw$TransectID)
library(stringr)
cover.dat.raw$TransectID <- str_extract(cover.dat.raw$Name, "^[^_]+_[^_]+")
cover.dat.raw$TransectID <- as.factor(cover.dat.raw$TransectID)
head(count.dat.raw)
head(cover.dat.raw)

library(magick)
library(dplyr)
library(progress)
# Function to get image dimensions
get_image_dimensions <- function(image_path) {
  image_info <- image_info(image_read(image_path))
  return(data.frame(Filename = basename(image_path), Width = image_info$width, Height = image_info$height))
}
# Directories containing the images
# image_dirs <- c(
#   "R:/IMAS/Antarctic_Seafloor/SQUIDLE_dataset_202411/",
#   "R:/IMAS/Antarctic_Seafloor/Clean_Data_For_Permanent_Storage/"
# )
image_dirs <- "R:/IMAS/Antarctic_Seafloor/SQUIDLE_dataset_202411/"
# Combine filenames from both dataframes
all_filenames <- unique(c(count.dat.raw$Name, cover.dat.raw$Name))
# List to store image dimensions
image_dimensions <- list()
# Initialize progress bar
total_images <- length(all_filenames)
pb <- progress_bar$new(
  format = "  Processing [:bar] :percent in :elapsed",
  total = total_images, clear = FALSE, width = 60
)
# Traverse through directories and get dimensions of each image, ignoring "thumbnails" folders
for (image_dir in image_dirs) {
  image_files <- list.files(image_dir, pattern = "\\.(JPG|jpg|jpeg|png)$", full.names = TRUE, recursive = TRUE)
  image_files <- image_files[!grepl("/thumbnails/", image_files)]
  for (image_file in image_files) {
    if (basename(image_file) %in% all_filenames) {
      image_dimensions[[length(image_dimensions) + 1]] <- get_image_dimensions(image_file)
      pb$tick()
    }
  }
}
# Combine all dimensions into a single dataframe
image_dimensions_df <- bind_rows(image_dimensions)
# Check for unmatched filenames
unmatched_filenames <- setdiff(all_filenames, image_dimensions_df$Filename)
if (length(unmatched_filenames) > 0) {
  warning("The following filenames did not find a matching image: ", paste(unmatched_filenames, collapse = ", "))
}
# Merge image dimensions with count.dat.raw
count.dat.raw <- count.dat.raw %>%
  left_join(image_dimensions_df, by = c("Name" = "Filename"))
# Merge image dimensions with cover.dat.raw
cover.dat.raw <- cover.dat.raw %>%
  left_join(image_dimensions_df, by = c("Name" = "Filename"))

## rename "Column" to "x"
names(cover.dat.raw)[4] <- "x"
## "Row" is y, but it's starting from the top. We need to measure from the bottom, so needs to be: y= image.height - Row
cover.dat.raw$y <- cover.dat.raw$Height-cover.dat.raw$Row

## remove NAs, which are placed into the AMC where we identified empty tubes
## OR CHANGE TO RELEVANT CATEGORY?
sel <- which(is.na(count.dat.raw$AMC))
count.dat.raw$Label_clean[sel]
count.dat.raw$AMC[sel] <- "TubesEmpty TemporaryLabelDoNotUse"
count.dat.raw$AMC_ID[sel] <- 28898

## subset to relevant data
cover.dat <- cover.dat.raw[,c(2,4,21,5,13,14,18:20)]
count.dat <- count.dat.raw[,c(2:7,13,14,10,17,18)]

# head(count.dat)
# head(cover.dat)

# Save the updated dataframes to CSV files
write.csv(count.dat, paste0(count.path,"Circumpolar_DownwardImages_ExhaustiveSearch_Annotations_ForSquidle_202412.csv"), row.names = FALSE)
write.csv(cover.dat, paste0(cover.path,"Circumpolar_DownwardImages_PointScore_Annotations_ForSquidle_202412.csv"), row.names = FALSE)

########################################################
## fixing x/y because I've calculated them wrong earlier (now fixed (Feb 2025))
count.dims <- read.csv(paste0(count.path,"Circumpolar_DownwardImages_ExhaustiveSearch_Annotations_ForSquidle_202412.csv"))
cover.dims <- read.csv(paste0(cover.path,"Circumpolar_DownwardImages_PointScore_Annotations_ForSquidle_202412.csv"))

## Width/Height mixed up in: PS118_09 PS118_12 PS118_81
sel <- which(cover.dims$TransectID%in%c("PS118_09","PS118_12","PS118_81"))
cover.dims$Width[sel] <- 5760
cover.dims$Height[sel] <- 3840

sel <- which(count.dims$TransectID%in%c("PS118_09","PS118_12","PS118_81"))
count.dims$Width[sel] <- 5760
count.dims$Height[sel] <- 3840

## AA2011 CTD54 and 56 were annotated portrait mode
sel <- which(cover.dims$TransectID%in%c("AA2011_CTD54","AA2011_CTD56"))
cover.dims$Width[sel] <- 2336
cover.dims$Height[sel] <- 3504

sel <- which(count.dims$TransectID%in%c("AA2011_CTD54","AA2011_CTD56"))
count.dims$Width[sel] <- 2336
count.dims$Height[sel] <- 3504

cover.dat.raw$Width <- cover.dims$Width
cover.dat.raw$Height <- cover.dims$Height
count.dat.raw$Width <- count.dims$Width
count.dat.raw$Height <- count.dims$Height

## rename "Row" to "y"
names(cover.dat.raw)[3] <- "point.pixels.row"
## rename "Column" to "x"
names(cover.dat.raw)[4] <- "point.pixels.col"

## pixels and points need to be defined starting from top left
# ## "Row" is y, but it's starting from the top. We need to measure from the bottom, so needs to be: y= image.height - Row
# cover.dat.raw$point.pixel.y <- cover.dat.raw$Height-cover.dat.raw$Row
cover.dat.raw$point.x <- cover.dat.raw$point.pixels.col/cover.dat.raw$Width
cover.dat.raw$point.y <- cover.dat.raw$point.pixels.row/cover.dat.raw$Height

## count data measured y from top down
names(count.dat.raw)[3:4] <- c("point.pixels.x","point.pixels.y")
count.dat.raw$point.x <- count.dat.raw$point.pixels.x/count.dat.raw$Width
count.dat.raw$point.y <- count.dat.raw$point.pixels.y/count.dat.raw$Height

#### AA2011: 10% crop, so 5% from top and 5% from left need to be adjusted for
## in AA 2011, the images on Squidle are already the cropped ones
sel <- which(cover.dat.raw$SurveyID=="AA2011")
width.adjustment <- round(0.1*cover.dat.raw$Width[sel])
height.adjustment <- round(0.1*cover.dat.raw$Height[sel])
cover.dat.raw$point.x[sel] <- cover.dat.raw$point.pixels.col[sel]/(cover.dat.raw$Width[sel]-width.adjustment)
cover.dat.raw$point.y[sel] <- cover.dat.raw$point.pixels.row[sel]/(cover.dat.raw$Height[sel]-height.adjustment)

## for the portrait annotations, y becomes x, and inverted x becomes y
sel2 <- which(cover.dat.raw$TransectID%in%c("AA2011_CTD54","AA2011_CTD56"))
width.adjustment <- round(0.1*cover.dat.raw$Width[sel2])
height.adjustment <- round(0.1*cover.dat.raw$Height[sel2])
cover.dat.raw$point.x[sel2] <- cover.dat.raw$point.pixels.row[sel2]/(cover.dat.raw$Height[sel2]-height.adjustment)
new.y <- rev(cover.dat.raw$point.pixels.col[sel2])
cover.dat.raw$point.y[sel2] <- new.y/(cover.dat.raw$Width[sel2]-width.adjustment)

##
sel <- which(count.dat.raw$SurveyID=="AA2011")
width.adjustment <-  round(0.1*count.dat.raw$Width[sel])
height.adjustment <- round(0.1*count.dat.raw$Height[sel])
count.dat.raw$point.x[sel] <- count.dat.raw$point.pixels.x[sel]/(count.dat.raw$Width[sel]-width.adjustment)
count.dat.raw$point.y[sel] <- count.dat.raw$point.pixels.y[sel]/(count.dat.raw$Height[sel]-height.adjustment)

## for the portrait annotations, y becomes x, and inverted x becomes y
sel2 <- which(count.dat.raw$TransectID%in%c("AA2011_CTD54","AA2011_CTD56"))
width.adjustment <-  round(0.1*count.dat.raw$Width[sel2])
height.adjustment <- round(0.1*count.dat.raw$Height[sel2])
count.dat.raw$point.x[sel2] <- count.dat.raw$point.pixels.y[sel2]/(count.dat.raw$Height[sel2]-height.adjustment)
new.y <- rev(count.dat.raw$point.pixels.x[sel2])
count.dat.raw$point.y[sel2] <- new.y/(count.dat.raw$Width[sel2]-width.adjustment)

#### PS96 & PS118: 20% crop, so 10% from top and 10% from left need to be adjusted for
sel <- which(cover.dat.raw$SurveyID%in%c("PS96","PS118"))
width.adjustment <- round(0.1*cover.dat.raw$Width[sel])
height.adjustment <- round(0.1*cover.dat.raw$Height[sel])
cover.dat.raw$point.x[sel] <- (cover.dat.raw$point.pixels.col[sel]+width.adjustment)/cover.dat.raw$Width[sel]
cover.dat.raw$point.y[sel] <- (cover.dat.raw$point.pixels.row[sel]+height.adjustment)/cover.dat.raw$Height[sel]

sel <- which(count.dat.raw$SurveyID%in%c("PS96","PS118"))
width.adjustment <-  round(0.1*count.dat.raw$Width[sel])
height.adjustment <- round(0.1*count.dat.raw$Height[sel])
count.dat.raw$point.x[sel] <- (count.dat.raw$point.pixels.x[sel]+width.adjustment)/ count.dat.raw$Width[sel]
count.dat.raw$point.y[sel] <- (count.dat.raw$point.pixels.y[sel]+height.adjustment)/count.dat.raw$Height[sel]

#### PS81: 20% crop, so 10% from top and 10% from left need to be adjusted for
## But, not in transects 185, 186, 188, 189, 197
sel <- which(cover.dat.raw$SurveyID=="PS81" & cover.dat.raw$TransectID%!in%c("PS81_185","PS81_186","PS81_188","PS81_189","PS81_197"))
width.adjustment <- round(0.1*cover.dat.raw$Width[sel])
height.adjustment <- round(0.1*cover.dat.raw$Height[sel])
cover.dat.raw$point.x[sel] <- (cover.dat.raw$point.pixels.col[sel]+width.adjustment)/cover.dat.raw$Width[sel]
cover.dat.raw$point.y[sel] <- (cover.dat.raw$point.pixels.row[sel]+height.adjustment)/cover.dat.raw$Height[sel]

sel <- which(count.dat.raw$SurveyID=="PS81" & count.dat.raw$TransectID%!in%c("PS81_185","PS81_186","PS81_188","PS81_189","PS81_197"))
width.adjustment <-  round(0.1*count.dat.raw$Width[sel])
height.adjustment <- round(0.1*count.dat.raw$Height[sel])
count.dat.raw$point.x[sel] <- (count.dat.raw$point.pixels.x[sel]+width.adjustment)/ count.dat.raw$Width[sel]
count.dat.raw$point.y[sel] <- (count.dat.raw$point.pixels.y[sel]+height.adjustment)/count.dat.raw$Height[sel]

## remove NAs, which are placed into the AMC where we identified empty tubes
## OR CHANGE TO RELEVANT CATEGORY?
sel <- which(is.na(count.dat.raw$AMC))
count.dat.raw$Label_clean[sel]
count.dat.raw$AMC[sel] <- "TubesEmpty TemporaryLabelDoNotUse"
count.dat.raw$AMC_ID[sel] <- 28898

## subset to relevant data
cover.dat <- cover.dat.raw[,c(2:4,5,13,14,18:22,17)]
count.dat <- count.dat.raw[,c(2:7,13,14,10,17:20,11)]

# Save the updated dataframes to CSV files
write.csv(count.dat, paste0(count.path,"Circumpolar_DownwardImages_ExhaustiveSearch_Annotations_ForSquidle_202502.csv"), row.names = FALSE)
write.csv(cover.dat, paste0(cover.path,"Circumpolar_DownwardImages_PointScore_Annotations_ForSquidle_202502.csv"), row.names = FALSE)

#####################################################################
## Another update to fix AA2011 10% bleed. This is already on the dimensions of the full image, so the calculations are different
## Split into separate campaigns, and rename to match squidle names for faster upload
count.dat <- read.csv(paste0(count.path,"Circumpolar_DownwardImages_ExhaustiveSearch_Annotations_ForSquidle_202502.csv"))
cover.dat <- read.csv(paste0(cover.path,"Circumpolar_DownwardImages_PointScore_Annotations_ForSquidle_202502.csv"))

#### AA2011: 10% crop, so 5% from top and 5% from left need to be adjusted for
## in AA 2011, the images on Squidle are already the cropped ones
sel <- which(cover.dat$SurveyID=="AA2011")
width.adjustment <- round(0.05*cover.dat$Width[sel])
height.adjustment <- round(0.05*cover.dat$Height[sel])
cover.dat$point.x[sel] <- (cover.dat$point.pixels.col[sel]+width.adjustment)/cover.dat$Width[sel]
cover.dat$point.y[sel] <- (cover.dat$point.pixels.row[sel]+height.adjustment)/cover.dat$Height[sel]

## for the portrait annotations, y becomes x, and inverted x becomes y
sel2 <- which(cover.dat$TransectID%in%c("AA2011_CTD54","AA2011_CTD56"))
width.adjustment <- round(0.05*cover.dat$Width[sel2])
height.adjustment <- round(0.05*cover.dat$Height[sel2])
cover.dat$point.x[sel2] <- (cover.dat$point.pixels.row[sel2]+height.adjustment)/cover.dat$Height[sel2]
new.y <- rev(cover.dat$point.pixels.col[sel2])
cover.dat$point.y[sel2] <- (new.y+width.adjustment)/cover.dat$Width[sel2]

##
sel <- which(count.dat$SurveyID=="AA2011")
width.adjustment <-  round(0.05*count.dat$Width[sel])
height.adjustment <- round(0.05*count.dat$Height[sel])
count.dat$point.x[sel] <- (count.dat$point.pixels.x[sel]+width.adjustment)/count.dat$Width[sel]
count.dat$point.y[sel] <- (count.dat$point.pixels.y[sel]+height.adjustment)/count.dat$Height[sel]

## for the portrait annotations, y becomes x, and inverted x becomes y
sel2 <- which(count.dat$TransectID%in%c("AA2011_CTD54","AA2011_CTD56"))
width.adjustment <-  round(0.05*count.dat$Width[sel2])
height.adjustment <- round(0.05*count.dat$Height[sel2])
count.dat$point.x[sel2] <- (count.dat$point.pixels.y[sel2]+height.adjustment)/count.dat$Height[sel2]
new.y <- rev(count.dat$point.pixels.x[sel2])
count.dat$point.y[sel2] <- (new.y+width.adjustment)/count.dat$Width[sel2]

# Save the updated dataframes to CSV files
write.csv(count.dat, paste0(count.path,"Circumpolar_DownwardImages_ExhaustiveSearch_Annotations_ForSquidle_202507.csv"), row.names = FALSE)
write.csv(cover.dat, paste0(cover.path,"Circumpolar_DownwardImages_PointScore_Annotations_ForSquidle_202507.csv"), row.names = FALSE)


## Split into separate campaigns, and rename to match squidle names for faster upload
count.dat <- read.csv(paste0(count.path,"Circumpolar_DownwardImages_ExhaustiveSearch_Annotations_ForSquidle_202507.csv"))
cover.dat <- read.csv(paste0(cover.path,"Circumpolar_DownwardImages_PointScore_Annotations_ForSquidle_202507.csv"))

head(count.dat)
names(count.dat)[1] <- "point.media.key"
names(count.dat)[8] <- "label.id"
names(count.dat)[9] <- "point.media.deployment.key"
names(count.dat)[10:11] <- c("point.pixels.width","point.pixels.height")

head(cover.dat)
names(cover.dat)[1] <- "point.media.key"
names(cover.dat)[6] <- "label.id"
names(cover.dat)[7] <- "point.media.deployment.key"
names(cover.dat)[8:9] <- c("point.pixels.width","point.pixels.height")

## subset by survey
count.dat$SurveyID <- as.factor(count.dat$SurveyID)
cover.dat$SurveyID <- as.factor(cover.dat$SurveyID)

for(i in c(1)){#,17,18,13)){ #1:21){
  print(i)
  s.ID <- levels(cover.dat$SurveyID)[i]
  cover.sel <- which(cover.dat$SurveyID==s.ID)
  cover.dat.subset <- cover.dat[cover.sel,]
  write.csv(cover.dat.subset, paste0(cover.path,"Circumpolar_DownwardImages_PointScore_Annotations_ForSquidle_202507_",s.ID,".csv"), row.names = FALSE)
  count.sel <- which(count.dat$SurveyID==s.ID)
  count.dat.subset <- count.dat[count.sel,]
  write.csv(count.dat.subset, paste0(count.path,"Circumpolar_DownwardImages_ExhaustiveSearch_Annotations_ForSquidle_202507_",s.ID,".csv"), row.names = FALSE)
}



