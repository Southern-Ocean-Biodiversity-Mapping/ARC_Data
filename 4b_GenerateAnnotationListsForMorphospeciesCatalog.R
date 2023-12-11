###################################################################################################
### this script updates the labels in the csv annotation files to publishable labels            ###
###################################################################################################

# 1) libraries and paths 
library(tidyverse)
library(writexl)
library(readxl)
library(plyr)

user = "Jan"
#user = "charley"
#user="nicole"

if (user == "Jan") {
  sci.dir <-      "C:/Users/jjansen/Desktop/science/"
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

## read in excel sheets
cover_list<-read_xlsx(path=paste0(ARC_Data.dir, "Annotation/Species_list_2023_01.xlsx"),sheet=1)
count_list<-read_xlsx(path=paste0(ARC_Data.dir, "Annotation/Species_list_2023_01.xlsx"),sheet=2)
## read in CATAMI reference sheets
cover_list2<-read_xlsx(path=paste0(ARC_Data.dir, "Annotation/Species_list_vs_CATAMI_2023_10.xlsx"),sheet=1)
count_list2<-read_xlsx(path=paste0(ARC_Data.dir, "Annotation/Species_list_vs_CATAMI_2023_10.xlsx"),sheet=2)

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
load(paste0(ARC_Data.dir,"Image_level_bio_202312.Rdata"))
ann.cover.dat3 <- ann.cover.dat.temp %>% left_join(img.metadata[,c("Filename.standardised", "survey")], by=c("Name"="Filename.standardised"))
names(ann.cover.dat3)[16] <- "SurveyID"
ann.cover.dat3$SurveyID <- sub("tan","TAN",ann.cover.dat3$SurveyID)

write.csv(ann.cover.dat3,file=paste0(cover.path,"Circumpolar_DownwardImages_PointScore_Annotations_202312_rawforpublication.csv"))

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

save(ann.cover.dat, file=paste0(cover.path,"Circumpolar_DownwardImages_PointScore_Annotations_202312.Rdata"))


##### Exhaustive-search data #####

## load the data
ann.files.csv <- list.files(counts.surveys.path,pattern="species.csv", recursive=TRUE, full.names=TRUE)
ann.count.dat.raw <- ldply(ann.files.csv, read_csv)

## lookup table to find images to generate crops from
ann.count.dat.raw$TransectID <- sub("(^[^_]+[_][^_]+)(.+$)","\\1",ann.count.dat.raw$filename)
ann.count.dat.raw2 <- ann.count.dat.raw %>% left_join(img.metadata[,c("Filename.standardised", "survey")], by=c("filename"="Filename.standardised"))
names(ann.count.dat.raw2)[16] <- "SurveyID"
ann.count.dat.raw2$SurveyID <- sub("tan","TAN",ann.count.dat.raw2$SurveyID)
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
save(ann.count.dat, file=paste0(count.path,"Circumpolar_DownwardImages_ExhaustiveSearch_Annotations_202312.Rdata"))

####################################################################################

## to check annotation numbers for the data paper:

## load annotation file
load(paste0(cover.path,"Circumpolar_DownwardImages_PointScore_Annotations_202312.Rdata"))
load(paste0(count.path,"Circumpolar_DownwardImages_ExhaustiveSearch_Annotations_202312.Rdata"))


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




