###################################################################################################
### This code converts removes 'unscorable' points and calculates prevalence on overall cover   ###
### This forms initial species list to make decisions about which species to keep/merge/exclude ###
###                                                                                             ###
### file with comments is read back in and exclusions/ aggregations made to form final datafile ###
### save richness & abundances of functional groups before species are removed from the data    ###
### N.Hill- Modified Jan 2022                                                                   ###
### J.Jansen: edited the script to create species lists for:                                    ###
###     the annotation library, the image data, 500m res and 2km res data                       ###
###################################################################################################

# 1) libraries and paths 
library(tidyverse)
library(writexl)
library(readxl)

user = "Jan"
#user = "charley"
#user="nicole"

if (user == "Jan") {
  sci.dir <-      "C:/Users/jjansen/OneDrive - University of Tasmania/science/"
  env.derived <-  paste0(sci.dir,"data_environmental/derived/")
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

output <- "AMC"
#output <- "ASAID"

##### 1: load image data and create table of count/prevalence per morphotype

##### 2: save excel-file and group relevant labels for final publishable annotation library

##### 3: read in excel-sheet, group relevant labels, calculate prevalence again and save into excel again for final publishable annotation library

##### 4: read back in the species excel-sheet and make groupings for modelling of COVER data



################################################################################
##### 1: load image data and create table of count/prevalence per morphotype

## load image data (saved within the annotation data files for the 500m and 2km data)
load(paste0(ARC_Data.dir, "annotation/Circumpolar_Annotation_Data_500m_202312.Rdata"))
 cover_cells_500m <- cover_cells
count_cells_500m <- count_cells
cell_metadata_500m <- cell_metadata
image_metadata_500m <- image_metadata
load(paste0(ARC_Data.dir, "annotation/Circumpolar_Annotation_Data_2km_202312.Rdata"))
cover_cells_2km <- cover_cells
count_cells_2km <- count_cells
cell_metadata_2km <- cell_metadata
image_metadata_2km <- image_metadata
rm(cell_metadata, image_metadata, count_cells, cover_cells)

## COVER - Prevalence
cover_prev<-data.frame(count=colSums(cover_images>0)) %>%
  mutate(., prev= round(count/nrow(cover_images), 3 )) %>%
  rownames_to_column(., var="Label") %>%
  arrange(., desc(count))%>%
  add_column(Exclude= "",
             Merge_With= "")

## COVER - Percent cover
# # remove unscorable points
# cover_images2<-cover_images %>%
#   dplyr::select(., -Unscorable) 
# 
# # percent cover of taxa at the image level:
# img.cover_pc<-prop.table(as.matrix(cover_images), margin=1)*100  
# 
# # proportion of all scored points assigned to each taxa
# img.cover_pc_overall<- data.frame(perc_overall=colSums(prop.table(as.matrix(cover_images))*100)) %>%
#   rownames_to_column(., var="Taxa") %>%
#   arrange(., desc(perc_overall))%>%
#   mutate(perc_overall= round(perc_overall, 3)) 
# 
# ## average percent cover at the image level:
# img.cover_pc_mean<- data.frame(perc_mean= apply(img.cover_pc, 2, mean))%>%
#   rownames_to_column(., var="Taxa") %>%
#   arrange(., desc(perc_mean))%>%
#   mutate(perc_mean= round(perc_mean, 3)) 

## COUNTS - Prevalence
count_prev<-data.frame(count=colSums(count_images>0)) %>%
  mutate(., prev= round(count/nrow(count_images), 3 )) %>%
  rownames_to_column(., var="Label") %>%
  arrange(., desc(count)) %>%
  #add extra columns for notes
  add_column(Exclude= "",
             Merge_With= "")

## COUNTS - Overall abundance
# img.count_ab_overall<-data.frame(total_count=colSums(count_images))%>%
#   rownames_to_column(., var="Taxa") %>%
#   arrange(., desc(total_count))



################################################################################
##### 2: save excel-file and edit/delete relevant labels for image annotation library

## write out excel-file to comment on
# write_xlsx(x= list(COVER_prevalence = cover_prev,
#           COUNT_prevalence = count_prev),
#           path=paste0(ARC_Data.dir, "Annotation/Species_list_2023_01.xlsx"))

################################################################################
##### 3: read in excel-sheet, group relevant labels for each resolution, calculate prevalence again, and save into excel again for final publishable annotation library

## read in excel sheets
cover_list<-read_xlsx(path=paste0(ARC_Data.dir, "Annotation/Species_list_2023_01.xlsx"),sheet=1)
count_list<-read_xlsx(path=paste0(ARC_Data.dir, "Annotation/Species_list_2023_01.xlsx"),sheet=2)

##### image data #####

### reformat data to long, merge, change names and convert back to wide
## Cover data
cover_images_long <- cover_images %>%
  mutate(cellID=rownames(cover_images))  %>% #add cellID
  pivot_longer(cols=`1Sub_Fine`:Echinoderms_Crinoids_Stalked, 
               names_to ="Label", values_to = "count") %>%           #long format and merge names to change
  left_join(cover_list[,c("Label", "Merge_With")])
cover_images_long$new <- ifelse(!is.na(cover_images_long$Merge_With), cover_images_long$Merge_With, cover_images_long$Label)
cover_images_renamed <- pivot_wider(cover_images_long, id_cols=cellID, names_from = new, values_from = count, values_fn=sum,values_fill = 0)
cover_mod <- cover_images_renamed
cover_mod$cellID <- as.factor(cover_mod$cellID)

## Count data
count_images_long <- count_images %>%
  mutate(cellID=rownames(count_images))  %>% #add cellID
  pivot_longer(cols=Echinoderms__Crinoid_unstalked:Molluscs__Gastropods_shell__LimpetLike, 
               names_to ="Label", values_to = "count") %>%           #long format and merge names to change
  left_join(count_list[,c("Label", "Merge_With")])
count_images_long$new <- ifelse(!is.na(count_images_long$Merge_With), count_images_long$Merge_With, count_images_long$Label)
count_images_renamed <- pivot_wider(count_images_long, id_cols=cellID, names_from = new, values_from = count, values_fn=sum,values_fill = 0)
#remove species to exclude
count_mod <- count_images_renamed %>%
  dplyr::select( - count_list$Label[which(count_list$Exclude =='x')])
count_mod$cellID <- as.factor(count_mod$cellID)

### prevalences and species names at image level for comparison with CATAMI
## cover at the image level:
cover_prev2 <- data.frame(count_img=colSums(cover_mod[,-1]>0)) %>%
  mutate(., prev_img= round(count_img/nrow(cover_mod), 3 )) %>%
  rownames_to_column(., var="Label") %>%
  arrange(., desc(count_img))%>%
  add_column(Name_to_publish= "", CATAMI_broad= "", CATAMI= "", CAAB= "",
             count_500m= "", prev_500m= "", Merge_With_2pc_500m= "",
             count_2km= "", prev_2km= "", Merge_With_2pc_2km= "")
cover.reordered <- order(cover_prev2$Label)
cover_prev.img <- cover_prev2[cover.reordered,]

## counts at the image level:
count_prev2<-data.frame(count_img=colSums(count_mod[,-1]>0)) %>%
  mutate(., prev_img= round(count_img/nrow(count_mod), 3 )) %>%
  rownames_to_column(., var="Label") %>%
  arrange(., desc(count_img)) %>%
  #add extra columns for notes
  add_column(Name_to_publish= "", CATAMI_broad= "", CATAMI= "", CAAB= "",
             count_500m= "", prev_500m= "", Exclude_500m="", Merge_With_2pc_500m= "",
             count_2km= "", prev_2km= "", Exclude_2km="", Merge_With_2pc_2km= "")
counts.reordered <- order(count_prev2$Label)
count_prev.img <- count_prev2[counts.reordered,]

##### 500m and 2km data #####

### reformat data to long, merge, change names and convert back to wide
## Cover data
cover_images_long <- cover_cells_500m %>%
  mutate(cellID=rownames(cover_cells_500m))  %>% #add cellID
  pivot_longer(cols=`1Sub_Fine`:Echinoderms_Crinoids_Stalked, 
               names_to ="Label", values_to = "count") %>%           #long format and merge names to change
  left_join(cover_list[,c("Label", "Merge_With")])
cover_images_long$new <- ifelse(!is.na(cover_images_long$Merge_With), cover_images_long$Merge_With, cover_images_long$Label)
cover_images_renamed <- pivot_wider(cover_images_long, id_cols=cellID, names_from = new, values_from = count, values_fn=sum,values_fill = 0)
cover_mod.500m <- cover_images_renamed
cover_mod.500m$cellID <- as.factor(cover_mod.500m$cellID)

cover_images_long <- cover_cells_2km %>%
  mutate(cellID=rownames(cover_cells_2km))  %>% #add cellID
  pivot_longer(cols=`1Sub_Fine`:Echinoderms_Crinoids_Stalked, 
               names_to ="Label", values_to = "count") %>%           #long format and merge names to change
  left_join(cover_list[,c("Label", "Merge_With")])
cover_images_long$new <- ifelse(!is.na(cover_images_long$Merge_With), cover_images_long$Merge_With, cover_images_long$Label)
cover_images_renamed <- pivot_wider(cover_images_long, id_cols=cellID, names_from = new, values_from = count, values_fn=sum,values_fill = 0)
cover_mod.2km <- cover_images_renamed
cover_mod.2km$cellID <- as.factor(cover_mod.2km$cellID)

## Count data
count_images_long <- count_cells_500m %>%
  mutate(cellID=rownames(count_cells_500m))  %>% #add cellID
  pivot_longer(cols=Echinoderms__Crinoid_unstalked:Molluscs__Gastropods_shell__LimpetLike, 
               names_to ="Label", values_to = "count") %>%           #long format and merge names to change
  left_join(count_list[,c("Label", "Merge_With")])
count_images_long$new <- ifelse(!is.na(count_images_long$Merge_With), count_images_long$Merge_With, count_images_long$Label)
count_images_renamed <- pivot_wider(count_images_long, id_cols=cellID, names_from = new, values_from = count, values_fn=sum,values_fill = 0)
#remove species to exclude
count_mod.500m <- count_images_renamed %>%
  dplyr::select( - count_list$Label[which(count_list$Exclude =='x')])
count_mod.500m$cellID <- as.factor(count_mod.500m$cellID)

count_images_long <- count_cells_2km %>%
  mutate(cellID=rownames(count_cells_2km))  %>% #add cellID
  pivot_longer(cols=Echinoderms__Crinoid_unstalked:Molluscs__Gastropods_shell__LimpetLike, 
               names_to ="Label", values_to = "count") %>%           #long format and merge names to change
  left_join(count_list[,c("Label", "Merge_With")])
count_images_long$new <- ifelse(!is.na(count_images_long$Merge_With), count_images_long$Merge_With, count_images_long$Label)
count_images_renamed <- pivot_wider(count_images_long, id_cols=cellID, names_from = new, values_from = count, values_fn=sum,values_fill = 0)
#remove species to exclude
count_mod.2km <- count_images_renamed %>%
  dplyr::select( - count_list$Label[which(count_list$Exclude =='x')])
count_mod.2km$cellID <- as.factor(count_mod.2km$cellID)


### prevalence and species names at 500m level
## cover:
cover_prev2 <- data.frame(count_500m=colSums(cover_mod.500m[,-1]>0)) %>%
  mutate(., prev_500m= round(count_500m/nrow(cover_mod.500m), 3 )) %>%
  rownames_to_column(., var="Label") %>%
  arrange(., desc(count_500m))
cover.reordered <- order(cover_prev2$Label)
cover_prev.500m <- cover_prev2[cover.reordered,]

## counts:
count_prev2<-data.frame(count_500m=colSums(count_mod.500m[,-1]>0)) %>%
  mutate(., prev_500m= round(count_500m/nrow(count_mod.500m), 3 )) %>%
  rownames_to_column(., var="Label") %>%
  arrange(., desc(count_500m))
counts.reordered <- order(count_prev2$Label)
count_prev.500m <- count_prev2[counts.reordered,]

### prevalences and species names at 2km level
## cover:
cover_prev2 <- data.frame(count_2km=colSums(cover_mod.2km[,-1]>0)) %>%
  mutate(., prev_2km= round(count_2km/nrow(cover_mod.2km), 3 )) %>%
  rownames_to_column(., var="Label") %>%
  arrange(., desc(count_2km))
cover.reordered <- order(cover_prev2$Label)
cover_prev.2km <- cover_prev2[cover.reordered,]

## counts:
count_prev2<-data.frame(count_2km=colSums(count_mod.2km[,-1]>0)) %>%
  mutate(., prev_2km= round(count_2km/nrow(count_mod.2km), 3 )) %>%
  rownames_to_column(., var="Label") %>%
  arrange(., desc(count_2km))
counts.reordered <- order(count_prev2$Label)
count_prev.2km <- count_prev2[counts.reordered,]

### merge dataframes into one:
cover_prev3 <- cover_prev.img
cover_prev3$count_500m <- cover_prev.500m$count_500m
cover_prev3$prev_500m <- cover_prev.500m$prev_500m
cover_prev3$count_2km <- cover_prev.2km$count_2km
cover_prev3$prev_2km <- cover_prev.2km$prev_2km

count_prev3 <- count_prev.img
count_prev3$count_500m <- count_prev.500m$count_500m
count_prev3$prev_500m  <- count_prev.500m$prev_500m
count_prev3$count_2km  <- count_prev.2km$count_2km
count_prev3$prev_2km   <- count_prev.2km$prev_2km

### save excel sheet
# write_xlsx(x=list(COVER_naming = cover_prev3, COUNT_naming = count_prev3),
#           path=paste0(ARC_Data.dir, "Annotation/Species_list_vs_CATAMI.xlsx"))

################################################################################
##### species to group or exclude being marked on the excel sheet, then save with date


##### 4: read back in the species excel-sheet and make groupings for modelling of COVER data
## read in excel sheets
cover_list<-read_xlsx(path=paste0(ARC_Data.dir, "Annotation/Species_list_vs_CATAMI_2023_10.xlsx"),sheet=1)
count_list<-read_xlsx(path=paste0(ARC_Data.dir, "Annotation/Species_list_vs_CATAMI_2023_10.xlsx"),sheet=2)

## read in cell-metadata
load(paste0(ARC_Data.dir,"Cell_level_env_500m_202312.Rdata"))
meta_env_500m <- cell_metadata_env
load(paste0(ARC_Data.dir,"Cell_level_env_2km_202412.Rdata"))
meta_env_2km <- cell_metadata_env
rm(cell_metadata_env, cell_metadata_env_scaled, transect.xy)

## number of scorable points
sc.img <- 108-cover_mod$Unscorable
sc.500m <- meta_env_500m$cover_points_scorable
sc.2km  <- meta_env_2km$cover_points_scorable

#### group morphospecies in img, 500m and 2km data as decided in the excel file
if(output=="AMC"){
  label_str <- "AMC"
}else label_str <- "Final_labels_img"
  
### reformat data to long, merge, change names und update labels to publishable names and convert back to wide
## Cover data
cover_images_long <- cover_mod[,-1] %>%
  mutate(cellID=cover_mod$cellID)  %>% #add cellID
  pivot_longer(cols=`Sub_Fine`:Echinoderms_Crinoids_Stalked,names_to ="Label", values_to = "count") %>% #long format and merge names to change
  left_join(cover_list[,c("Label", "AMC", "Final_labels_img", "Merge_With_1pc_img")])
cover_images_long$new <- ifelse(!is.na(cover_images_long$Merge_With_1pc_img), cover_images_long$Merge_With_1pc_img, cover_images_long$AMC)
cover_images_renamed <- pivot_wider(cover_images_long, id_cols=cellID, names_from = new, values_from = count, values_fn=sum, values_fill = 0)
#cover_images_renamed <- pivot_wider(cover_images_long, id_cols=cellID, names_from=Final_labels_img, values_from=count, values_fn=sum, values_fill=0)
#cover_images_renamed <- pivot_wider(cover_images_long, id_cols=cellID, names_from="AMC", values_from=count, values_fn=sum, values_fill=0)
cover_mod <- cover_images_renamed
cover_mod$cellID <- as.factor(cover_mod$cellID)

cover_images_long <- cover_mod.500m[,-1] %>%
  mutate(cellID=cover_mod.500m$cellID)  %>% #add cellID
  pivot_longer(cols=`Sub_Fine`:Echinoderms_Crinoids_Stalked, names_to ="Label", values_to = "count") %>% #long format and merge names to change
  left_join(cover_list[,c("Label", "AMC","Merge_With_2pc_500m")])
cover_images_long$new <- ifelse(!is.na(cover_images_long$Merge_With_2pc_500m), cover_images_long$Merge_With_2pc_500m, cover_images_long$AMC)
cover_images_renamed <- pivot_wider(cover_images_long, id_cols=cellID, names_from = new, values_from = count, values_fn=sum,values_fill = 0)
# cover_images_renamed <- pivot_wider(cover_images_long, id_cols=cellID, names_from="AMC", values_from=count, values_fn=sum, values_fill=0)
cover_mod.500m <- cover_images_renamed
cover_mod.500m$cellID <- as.factor(cover_mod.500m$cellID)

cover_images_long <- cover_mod.2km[,-1] %>%
  mutate(cellID=cover_mod.2km$cellID)  %>% #add cellID
  pivot_longer(cols=`Sub_Fine`:Echinoderms_Crinoids_Stalked, names_to ="Label", values_to = "count") %>% #long format and merge names to change
  left_join(cover_list[,c("Label", "AMC","Merge_With_2pc_2km")])
cover_images_long$new <- ifelse(!is.na(cover_images_long$Merge_With_2pc_2km), cover_images_long$Merge_With_2pc_2km, cover_images_long$AMC)
cover_images_renamed <- pivot_wider(cover_images_long, id_cols=cellID, names_from = new, values_from = count, values_fn=sum,values_fill = 0)
# cover_images_renamed <- pivot_wider(cover_images_long, id_cols=cellID, names_from="AMC", values_from=count, values_fn=sum, values_fill=0)
cover_mod.2km <- cover_images_renamed
cover_mod.2km$cellID <- as.factor(cover_mod.2km$cellID)

## Count data
count_images_long <- count_mod[,-1] %>%
  mutate(cellID=count_mod$cellID)  %>% #add cellID
  pivot_longer(cols=Echinoderms__Crinoid_unstalked:Molluscs__Gastropods_shell__LimpetLike, names_to ="Label", values_to = "count") %>% #long format and merge names to change
  left_join(count_list[,c("Label", "AMC","Final_labels_img", "Merge_With_1pc_img")])
count_images_long$new <- ifelse(!is.na(count_images_long$Merge_With_1pc_img), count_images_long$Merge_With_1pc_img, count_images_long$AMC)
count_images_renamed <- pivot_wider(count_images_long, id_cols=cellID, names_from = new, values_from = count, values_fn=sum,values_fill = 0)
# count_images_renamed <- pivot_wider(count_images_long, id_cols=cellID, names_from = Final_labels_img, values_from = count, values_fn=sum, values_fill = 0)
# count_images_renamed <- pivot_wider(count_images_long, id_cols=cellID, names_from = AMC, values_from = count, values_fn=sum, values_fill = 0)
#remove species to exclude
count_mod <- count_images_renamed %>%
  dplyr::select( - count_list$Label[which(count_list$Exclude_img =='x')])
count_mod$cellID <- as.factor(count_mod$cellID)

count_images_long <- count_mod.500m[,-1] %>%
  mutate(cellID=count_mod.500m$cellID)  %>% #add cellID
  pivot_longer(cols=Echinoderms__Crinoid_unstalked:Molluscs__Gastropods_shell, names_to ="Label", values_to = "count") %>% #long format and merge names to change
  left_join(count_list[,c("Label", "AMC","Merge_With_2pc_500m")])
count_images_long$new <- ifelse(!is.na(count_images_long$Merge_With_2pc_500m), count_images_long$Merge_With_2pc_500m, count_images_long$AMC)
count_images_renamed <- pivot_wider(count_images_long, id_cols=cellID, names_from = new, values_from = count, values_fn=sum,values_fill = 0)
# count_images_renamed <- pivot_wider(count_images_long, id_cols=cellID, names_from =AMC, values_from = count, values_fn=sum,values_fill = 0)
#remove species to exclude
count_mod.500m <- count_images_renamed %>%
  dplyr::select( - count_list$Label[which(count_list$Exclude_500m =='x')])
count_mod.500m$cellID <- as.factor(count_mod.500m$cellID)

count_images_long <- count_mod.2km[,-1] %>%
  mutate(cellID=count_mod.2km$cellID)  %>% #add cellID
  pivot_longer(cols=Echinoderms__Crinoid_unstalked:Molluscs__Gastropods_shell, names_to ="Label", values_to = "count") %>% #long format and merge names to change
  left_join(count_list[,c("Label", "AMC","Merge_With_2pc_2km")])
count_images_long$new <- ifelse(!is.na(count_images_long$Merge_With_2pc_2km), count_images_long$Merge_With_2pc_2km, count_images_long$Label)
count_images_renamed <- pivot_wider(count_images_long, id_cols=cellID, names_from = new, values_from = count, values_fn=sum,values_fill = 0)
#remove species to exclude
count_mod.2km <- count_images_renamed %>%
  dplyr::select( - count_list$Label[which(count_list$Exclude_2km =='x')])
count_mod.2km$cellID <- as.factor(count_mod.2km$cellID)

###########################################################################################
######## THE BELOW CODE IS REPETITIVE BUT WORKS...!!! COULD CLEAN UP BUT NO TIME... #######

#### 4a) IMAGES:
## Calculate abundances of functional groups and richness on the raw data (before species are excluded)
## for the count data
dat.count <- count_images
## for the cover data
## need to remove the cellID from the columns to properly calculate rowSums etc
dat <- data.frame(cover_mod[,-1])
rownames(dat) <- cover_mod$cellID
scorable.pts <- sc.img

### cover data
## names of faunal groups for cover_cells_renamed (minus cellID):
dataset.names <- names(dat)
dataset.names[order(dataset.names)]
## selector for each faunal class
sel_S <- grep("Sp",substr(dataset.names,1,2))
sel_O <- grep("Oc",substr(dataset.names,1,2))
sel_B <- grep("Bry",substr(dataset.names,1,3))
sel_M <- grep("Mo",substr(dataset.names,1,2))
sel_E <- grep("Ech",substr(dataset.names,1,3))
sel_Asc <- grep("Asc",substr(dataset.names,1,3))
sel_TW <- grep("Worms...Polychaetes...T",dataset.names)
sel_Hy <- grep("Hyd",dataset.names)
## selector for functional group
sel_SF <- c(sel_S,sel_O,sel_B,sel_Asc,sel_TW,sel_Hy)
## selector for sediment class etc
sel_sed_soft <- c(grep("Sand",dataset.names),grep("Pebble",dataset.names))
sel_sed_loose <- c(grep("Cobble",dataset.names),grep("Biologenic",dataset.names)) 
sel_sed_hard <- c(grep("Bould",dataset.names),grep("Rock",dataset.names))
sel_sed <- c(sel_sed_soft,sel_sed_loose,sel_sed_hard)

sel_noid.cov <- grep("Unidentif",dataset.names)
sel_unsc.cov <- grep("Unscorable",dataset.names)

## calculating species richness etc, keeping in mind that the first column is the cellID
cover_pa <- dat
cover_pa[cover_pa>0] <- 1
cover_SF.prop <- rowSums(dat[,sel_SF])/scorable.pts
cover_SF      <- rowSums(dat[,sel_SF])
cover_SF_pa   <- cover_SF
cover_SF_pa[cover_SF>0] <- 1
richness      <- rowSums(cover_pa[,-sel_sed])
richness.l    <- rowSums(cover_pa[,-sel_sed])/log(scorable.pts)
cover_all.prop<- rowSums(dat[,-sel_sed])/scorable.pts
cover_all     <- rowSums(dat[,-sel_sed])
## Bryozoans
cover_B.prop <- rowSums(dat[,sel_B])/scorable.pts
cover_B      <- rowSums(dat[,sel_B])
cover_B_pa   <- cover_B
cover_B_pa[cover_B>0] <- 1
## Sponges
cover_S.prop <- rowSums(dat[,sel_S])/scorable.pts
cover_S      <- rowSums(dat[,sel_S])
cover_S_pa   <- cover_S
cover_S_pa[cover_S>0] <- 1
## large species groupings - cover
cover_groupings <- data.frame(cbind(cover_all, cover_SF, cover_SF_pa, richness, richness.l, cover_B, cover_B_pa, cover_S, cover_S_pa))

##### setup biological data - counts
## names of faunal groups for dat.count:
count.names <- names(dat.count)

## selector for each faunal class
sel_noid <- grep("NoID",count.names)
sel_echino <- grep("Echinoderms",count.names)
sel_crust <- grep("Crustacea",count.names)

count_mobile <- rowSums(dat.count[,-sel_noid]) ## remove NoIDs
count_echino <- rowSums(dat.count[,sel_echino])
count_crust  <- rowSums(dat.count[,sel_crust])

count_pa <- dat.count
count_pa[count_pa>0] <- 1
count_richness <- rowSums(count_pa[,-sel_sed])

## large species groupings - counts
count_groupings <- data.frame(cbind(count_mobile, count_echino, count_crust, count_richness))

#########################
# save outputs, but first add metadata of images
sel.metadat <- match(rownames(cover_images),image_metadata_500m$Filename.standardised)
img.metadata <- image_metadata_500m[sel.metadat,1:7]
names(img.metadata)[7] <- "cellID_500m"
img.metadata$cellID_2km <- image_metadata_2km$cellID[sel.metadat]
img.metadata[,9:(ncol(image_metadata_500m)+1)] <- image_metadata_500m[sel.metadat,8:ncol(image_metadata_500m)]

save(cover_mod, count_mod, cover_groupings, count_groupings, img.metadata, file=paste0(ARC_Data.dir,"Image_level_bio_202412.Rdata"))
###

#################################
#### 4b) 500m:
## Calculate abundances of functional groups and richness on the raw data (before species are excluded)
## for the count data
dat.count <- count_cells_500m
## for the cover data
## need to remove the cellID from the columns to properly calculate rowSums etc
dat <- data.frame(cover_mod.500m[,-1])
rownames(dat) <- cover_mod.500m$cellID
scorable.pts <- sc.500m

### cover data
## names of faunal groups for cover_cells_renamed (minus cellID):
dataset.names <- names(dat)
dataset.names[order(dataset.names)]
## selector for each faunal class
sel_S <- grep("Sp",substr(dataset.names,1,2))
sel_O <- grep("Oc",substr(dataset.names,1,2))
sel_B <- grep("Bry",substr(dataset.names,1,3))
sel_M <- grep("Mo",substr(dataset.names,1,2))
sel_E <- grep("Ech",substr(dataset.names,1,3))
sel_Asc <- grep("Asc",substr(dataset.names,1,3))
sel_TW <- grep("Worms_Polychaetes_T",dataset.names)
sel_Hy <- grep("Hyd",dataset.names)
## selector for functional group
sel_SF <- c(sel_S,sel_O,sel_B,sel_Asc,sel_TW,sel_Hy)
## selector for sediment class etc
sel_sed_soft <- c(grep("Fine",dataset.names),grep("PbGrv",dataset.names))
sel_sed_loose <- c(grep("Cbble",dataset.names),grep("BioRu",dataset.names),grep("BioShl",dataset.names),grep("BioOth",dataset.names)) 
sel_sed_hard <- c(grep("Bould",dataset.names),grep("Rock",dataset.names))
sel_sed <- grep("Sub_",dataset.names)
sel_noid.cov <- grep("NoID",dataset.names)
sel_unsc.cov <- grep("Unscorable",dataset.names)

## calculating species richness etc, keeping in mind that the first column is the cellID
cover_pa <- dat
cover_pa[cover_pa>0] <- 1
cover_SF.prop <- rowSums(dat[,sel_SF])/scorable.pts
cover_SF      <- rowSums(dat[,sel_SF])
cover_SF_pa   <- cover_SF
cover_SF_pa[cover_SF>0] <- 1
richness      <- rowSums(cover_pa[,-sel_sed])
richness.l    <- rowSums(cover_pa[,-sel_sed])/log(scorable.pts)
cover_all.prop<- rowSums(dat[,-sel_sed])/scorable.pts
cover_all     <- rowSums(dat[,-sel_sed])
## Bryozoans
cover_B.prop <- rowSums(dat[,sel_B])/scorable.pts
cover_B      <- rowSums(dat[,sel_B])
cover_B_pa   <- cover_B
cover_B_pa[cover_B>0] <- 1
## Sponges
cover_S.prop <- rowSums(dat[,sel_S])/scorable.pts
cover_S      <- rowSums(dat[,sel_S])
cover_S_pa   <- cover_S
cover_S_pa[cover_S>0] <- 1
## large species groupings - cover
cover_groupings <- data.frame(cbind(cover_all, cover_SF, cover_SF_pa, richness, richness.l, cover_B, cover_B_pa, cover_S, cover_S_pa))

##### setup biological data - counts
## names of faunal groups for dat.count:
count.names <- names(dat.count)

## selector for each faunal class
sel_noid <- grep("NoID",count.names)
sel_echino <- grep("Echinoderms",count.names)
sel_crust <- grep("Crustacea",count.names)

count_mobile <- rowSums(dat.count[,-sel_noid]) ## remove NoIDs
count_echino <- rowSums(dat.count[,sel_echino])
count_crust  <- rowSums(dat.count[,sel_crust])

count_pa <- dat.count
count_pa[count_pa>0] <- 1
count_richness <- rowSums(count_pa[,-sel_sed])
count_richness.l <- rowSums(count_pa[,-sel_sed])/log(meta_env_500m$counts_N[!is.na(meta_env_500m$counts_N)])

## large species groupings - counts
count_groupings <- data.frame(cbind(count_mobile, count_echino, count_crust, count_richness, count_richness.l))

#############
# load(paste0(ARC_Data.dir, "Cell_level_env_",res,".Rdata"))
# #join cover data back to cell metadata and environmental data
# cover_mod_env<-left_join(cell_metadata_env, cover_mod, by="cellID")
# #join count data back to cell metadata and environmental data
# count_mod_env<-left_join(cell_metadata_env, count_mod, by="cellID")
# save outputs
save(cover_mod.500m, count_mod.500m, cover_groupings, count_groupings, file=paste0(ARC_Data.dir,"Cell_level_bio_2pc_500m_202312.Rdata"))


#################################
#### 4c) 2km:
## Calculate abundances of functional groups and richness on the raw data (before species are excluded)
## for the count data
dat.count <- count_cells_2km
## for the cover data
## need to remove the cellID from the columns to properly calculate rowSums etc
dat <- data.frame(cover_mod.2km[,-1])
rownames(dat) <- cover_mod.2km$cellID
scorable.pts <- sc.2km

### cover data
## names of faunal groups for cover_cells_renamed (minus cellID):
dataset.names <- names(dat)
dataset.names[order(dataset.names)]
## selector for each faunal class
sel_S <- grep("Sp",substr(dataset.names,1,2))
sel_O <- grep("Oc",substr(dataset.names,1,2))
sel_B <- grep("Bry",substr(dataset.names,1,3))
sel_M <- grep("Mo",substr(dataset.names,1,2))
sel_E <- grep("Ech",substr(dataset.names,1,3))
sel_Asc <- grep("Asc",substr(dataset.names,1,3))
sel_TW <- grep("Worms_Polychaetes_T",dataset.names)
sel_Hy <- grep("Hyd",dataset.names)
## selector for functional group
sel_SF <- c(sel_S,sel_O,sel_B,sel_Asc,sel_TW,sel_Hy)
## selector for sediment class etc
sel_sed_soft <- c(grep("Fine",dataset.names),grep("PbGrv",dataset.names))
sel_sed_loose <- c(grep("Cbble",dataset.names),grep("BioRu",dataset.names),grep("BioShl",dataset.names),grep("BioOth",dataset.names)) 
sel_sed_hard <- c(grep("Bould",dataset.names),grep("Rock",dataset.names))
sel_sed <- grep("Sub_",dataset.names)
sel_noid.cov <- grep("NoID",dataset.names)
sel_unsc.cov <- grep("Unscorable",dataset.names)

## calculating species richness etc, keeping in mind that the first column is the cellID
cover_pa <- dat
cover_pa[cover_pa>0] <- 1
cover_SF.prop <- rowSums(dat[,sel_SF])/scorable.pts
cover_SF      <- rowSums(dat[,sel_SF])
cover_SF_pa   <- cover_SF
cover_SF_pa[cover_SF>0] <- 1
richness      <- rowSums(cover_pa[,-sel_sed])
richness.l    <- rowSums(cover_pa[,-sel_sed])/log(scorable.pts)
cover_all.prop<- rowSums(dat[,-sel_sed])/scorable.pts
cover_all     <- rowSums(dat[,-sel_sed])
## Bryozoans
cover_B.prop <- rowSums(dat[,sel_B])/scorable.pts
cover_B      <- rowSums(dat[,sel_B])
cover_B_pa   <- cover_B
cover_B_pa[cover_B>0] <- 1
## Sponges
cover_S.prop <- rowSums(dat[,sel_S])/scorable.pts
cover_S      <- rowSums(dat[,sel_S])
cover_S_pa   <- cover_S
cover_S_pa[cover_S>0] <- 1
## large species groupings - cover
cover_groupings <- data.frame(cbind(cover_all, cover_SF, cover_SF_pa, richness, richness.l, cover_B, cover_B_pa, cover_S, cover_S_pa))

##### setup biological data - counts
## names of faunal groups for dat.count:
count.names <- names(dat.count)

## selector for each faunal class
sel_noid <- grep("NoID",count.names)
sel_echino <- grep("Echinoderms",count.names)
sel_crust <- grep("Crustacea",count.names)

count_mobile <- rowSums(dat.count[,-sel_noid]) ## remove NoIDs
count_echino <- rowSums(dat.count[,sel_echino])
count_crust  <- rowSums(dat.count[,sel_crust])

count_pa <- dat.count
count_pa[count_pa>0] <- 1
count_richness <- rowSums(count_pa[,-sel_sed])

count_richness.l <- rowSums(count_pa[,-sel_sed])/log(meta_env_2km$counts_N[!is.na(meta_env_2km$counts_N)])
## large species groupings - counts
count_groupings <- data.frame(cbind(count_mobile, count_echino, count_crust, count_richness, count_richness.l))

#############
# load(paste0(ARC_Data.dir, "Cell_level_env_",res,".Rdata"))
# #join cover data back to cell metadata and environmental data
# cover_mod_env<-left_join(cell_metadata_env, cover_mod, by="cellID")
# #join count data back to cell metadata and environmental data
# count_mod_env<-left_join(cell_metadata_env, count_mod, by="cellID")
# save outputs
save(cover_mod.2km, count_mod.2km, cover_groupings, count_groupings, file=paste0(ARC_Data.dir,"Cell_level_bio_2pc_2km_202412.Rdata"))




