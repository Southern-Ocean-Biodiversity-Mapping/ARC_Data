###################################################################################################
#### This code converts removes 'unscorable' points and calculates prevalence on overall cover  ###
### This forms initial species list to make decisions about which species to keep/merge/exclude ###
###                                                                                             ###
### file with comments is read back in and exclusions/ aggregations made to form final datafile ###
### save richness & abundances of functional groups before species are removed from the data    ###
### N.Hill- Modified Jan 2022                                                                   ###
###################################################################################################


# 1) libraries and paths 
library(tidyverse)
library(writexl)
library(readxl)

user = "Jan"
#user = "charley"
#user="nicole"

if (user == "Jan") {
  sci.dir <-      "C:/Users/jjansen/Desktop/science/"
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

######
#res <- "500m"
res <- "2km"
######


load(paste0(ARC_Data.dir, "annotation/Circumpolar_Annotation_Data_",res,".Rdata"))


## 2) COVER ----
## 2a) Prevalence
cover_prev<-data.frame(count=colSums(cover_cells>0)) %>%
  mutate(., prev= round(count/nrow(cover_cells), 3 )) %>%
  rownames_to_column(., var="Label") %>%
  arrange(., desc(count))%>%
  add_column(Exclude= "",
             Merge_With= "",
             Other_Notes="")
#... & at the image level:
img.cover_prev<-data.frame(count=colSums(cover_images>0)) %>%
  mutate(., prev= round(count/nrow(cover_images), 3 )) %>%
  rownames_to_column(., var="Label") %>%
  arrange(., desc(count))%>%
  add_column(Exclude= "",
             Merge_With= "",
             Other_Notes="")

##2b) Percent cover
# remove unscorable points
cover_cells2<-cover_cells %>%
  dplyr::select(., -Unscorable) 
#... & at the image level:
cover_images2<-cover_images %>%
  dplyr::select(., -Unscorable) 

# percent cover of taxa in each cell
cover_pc<-prop.table(as.matrix(cover_cells), margin=1)*100 
#... & at the image level:
img.cover_pc<-prop.table(as.matrix(cover_images), margin=1)*100  

# proportion of all scored points assigned to each taxa
cover_pc_overall<- data.frame(perc_overall=colSums(prop.table(as.matrix(cover_cells))*100)) %>%
  rownames_to_column(., var="Taxa") %>%
  arrange(., desc(perc_overall))%>%
  mutate(perc_overall= round(perc_overall, 3)) 
#... & at the image level:
img.cover_pc_overall<- data.frame(perc_overall=colSums(prop.table(as.matrix(cover_images))*100)) %>%
  rownames_to_column(., var="Taxa") %>%
  arrange(., desc(perc_overall))%>%
  mutate(perc_overall= round(perc_overall, 3)) 

## average percent cover across all cells
cover_pc_mean<- data.frame(perc_mean= apply(cover_pc, 2, mean))%>%
rownames_to_column(., var="Taxa") %>%
  arrange(., desc(perc_mean))%>%
  mutate(perc_mean= round(perc_mean, 3)) 
#very similar to overall percent cover
#... & at the image level:
img.cover_pc_mean<- data.frame(perc_mean= apply(img.cover_pc, 2, mean))%>%
  rownames_to_column(., var="Taxa") %>%
  arrange(., desc(perc_mean))%>%
  mutate(perc_mean= round(perc_mean, 3)) 


## 3) COUNTS ----
## 3a) Prevalence
count_prev<-data.frame(count=colSums(count_cells>0)) %>%
  mutate(., prev= round(count/nrow(count_cells), 3 )) %>%
  rownames_to_column(., var="Label") %>%
  arrange(., desc(count)) %>%
#add extra columns for notes
  add_column(Exclude= "",
             Merge_With= "",
             Other_Notes="")
#... & at the image level:
img.count_prev<-data.frame(count=colSums(count_images>0)) %>%
  mutate(., prev= round(count/nrow(count_images), 3 )) %>%
  rownames_to_column(., var="Label") %>%
  arrange(., desc(count)) %>%
  #add extra columns for notes
  add_column(Exclude= "",
             Merge_With= "",
             Other_Notes="")

## 3b) Overall abundance
count_ab_overall<-data.frame(total_count=colSums(count_cells))%>%
  rownames_to_column(., var="Taxa") %>%
  arrange(., desc(total_count))
#... & at the image level:
img.count_ab_overall<-data.frame(total_count=colSums(count_images))%>%
  rownames_to_column(., var="Taxa") %>%
  arrange(., desc(total_count))


## 4) Save to excel to comments and notes on aggregation etc.

# write_xlsx(x= list(COVER_prevalence=cover_prev,
#           COVER_overall=cover_pc_overall,
#           COUNT_prevalence= count_prev,
#           COUNT_totAb= count_ab_overall),
#           path=paste0(ARC_Data.dir, "Annotation/Species_list_",res,"_2022_10.xlsx"))

# write_xlsx(x= list(COVER_prevalence=img.cover_prev,
#           COVER_overall=img.cover_pc_overall,
#           COUNT_prevalence= img.count_prev,
#           COUNT_totAb= img.count_ab_overall),
#           path=paste0(ARC_Data.dir, "Annotation/Species_list_images_2022_12.xlsx"))


#######################################################################

## 5) Read commented file back in and make exclusions and aggregations
# blanks read in as NAs

## 5a) cover data
mod_cover_list<-read_xlsx(path=paste0(ARC_Data.dir, "Annotation/Species_list_",res,"_2022_12.xlsx"),sheet=1)
names(mod_cover_list)[5:6]<- c("Merge_1pc", "Merge_2pc")

img.mod_cover_list<-read_xlsx(path=paste0(ARC_Data.dir, "Annotation/Species_list_images_2022_12.xlsx"),sheet=1)

#reformat to long, merge, change names and convert back to wide
cover_cells_long<-cover_cells %>%
  mutate(cellID=as.numeric(rownames(cover_cells)))  %>% #add cellID
  pivot_longer( cols=`1Sub_Fine`:Echinoderms_Crinoids_Stalked, 
                names_to ="Label", values_to = "count") %>%           #long format and merge names to change
  left_join( mod_cover_list[,c("Label", "Merge_2pc")])
cover_cells_long$new<-  ifelse(!is.na(cover_cells_long$Merge_2pc), cover_cells_long$Merge_2pc, cover_cells_long$Label)
cover_cells_renamed<-pivot_wider(cover_cells_long, id_cols=cellID, names_from = new, values_from = count, values_fn=sum,values_fill = 0)

#also for images
cover_images_long<-cover_images %>%
  mutate(cellID=rownames(cover_images))  %>% #add cellID
  pivot_longer( cols=`1Sub_Fine`:Echinoderms_Crinoids_Stalked, 
                names_to ="Label", values_to = "count") %>%           #long format and merge names to change
  left_join(img.mod_cover_list[,c("Label", "Merge_2pc")])
cover_images_long$new<-  ifelse(!is.na(cover_images_long$Merge_2pc), cover_images_long$Merge_2pc, cover_images_long$Label)
cover_images_renamed<-pivot_wider(cover_images_long, id_cols=cellID, names_from = new, values_from = count, values_fn=sum,values_fill = 0)

#also for list of annotation labels for library
cover_annimages_long<-cover_images %>%
  mutate(cellID=rownames(cover_images))  %>% #add cellID
  pivot_longer( cols=`1Sub_Fine`:Echinoderms_Crinoids_Stalked, 
                names_to ="Label", values_to = "count") %>%           #long format and merge names to change
  left_join(img.mod_cover_list[,c("Label", "Merge_For_Annotation_Library")])
cover_annimages_long$new<-  ifelse(!is.na(cover_annimages_long$Merge_For_Annotation_Library),cover_annimages_long$Merge_For_Annotation_Library, cover_annimages_long$Label)
cover_annimages_renamed<-pivot_wider(cover_annimages_long, id_cols=cellID, names_from = new, values_from = count, values_fn=sum,values_fill = 0)


#remove species to exclude (DON'T REMOVE BECAUSE WE MESS UP TOTAL COVER DATA?)
cover_mod<-cover_cells_renamed
# cover_mod<-cover_cells_renamed %>%
#   select( - mod_cover_list$Label[which(mod_cover_list$Exclude =='x')])
cover_mod$cellID<-as.factor(cover_mod$cellID)

img.cover_mod<-cover_images_renamed
img.cover_mod$cellID<-as.factor(img.cover_mod$cellID)

img.ann.cover_mod<-cover_annimages_renamed
img.ann.cover_mod$cellID<-as.factor(img.ann.cover_mod$cellID)

## 5b) count data
mod_count_list<-read_xlsx(path=paste0(ARC_Data.dir, "Annotation/Species_list_",res,"_2022_12.xlsx"),sheet=3)
names(mod_count_list)[5:6]<- c("Merge_1pc", "Merge_2pc")

img.mod_count_list<-read_xlsx(path=paste0(ARC_Data.dir, "Annotation/Species_list_images_2022_12.xlsx"),sheet=3)

#reformat to long, merge, change names and convert back to wide
count_cells_long<-count_cells %>%
  mutate(cellID=as.numeric(rownames(count_cells)))  %>% #add cellID
  pivot_longer( cols=Echinoderms__Crinoid_unstalked:Molluscs__Gastropods_shell__LimpetLike, 
                        names_to ="Label", values_to = "count") %>%           #long format and merge names to change
  left_join( mod_count_list[,c("Label", "Merge_2pc")])
count_cells_long$new<-  ifelse(!is.na(count_cells_long$Merge_2pc), count_cells_long$Merge_2pc, count_cells_long$Label)
count_cells_renamed<-pivot_wider(count_cells_long, id_cols=cellID, names_from = new, values_from = count, values_fn=sum,values_fill = 0)

#also for images
count_images_long<-count_images %>%
  mutate(cellID=rownames(count_images))  %>% #add cellID
  pivot_longer( cols=Echinoderms__Crinoid_unstalked:Molluscs__Gastropods_shell__LimpetLike, 
                names_to ="Label", values_to = "count") %>%           #long format and merge names to change
  left_join( img.mod_count_list[,c("Label", "Merge_2pc")])
count_images_long$new<-  ifelse(!is.na(count_images_long$Merge_2pc), count_images_long$Merge_2pc, count_images_long$Label)
count_images_renamed<-pivot_wider(count_images_long, id_cols=cellID, names_from = new, values_from = count, values_fn=sum,values_fill = 0)

#also list of annotation labels for library
count_annimages_long<-count_images %>%
  mutate(cellID=rownames(count_images))  %>% #add cellID
  pivot_longer( cols=Echinoderms__Crinoid_unstalked:Molluscs__Gastropods_shell__LimpetLike, 
                names_to ="Label", values_to = "count") %>%           #long format and merge names to change
  left_join( img.mod_count_list[,c("Label", "Merge_For_Annotation_Library")])
count_annimages_long$new<-  ifelse(!is.na(count_annimages_long$Merge_For_Annotation_Library), count_annimages_long$Merge_For_Annotation_Library, count_annimages_long$Label)
count_annimages_renamed<-pivot_wider(count_annimages_long, id_cols=cellID, names_from = new, values_from = count, values_fn=sum,values_fill = 0)

#remove species to exclude
count_mod<-count_cells_renamed %>%
 dplyr::select( - mod_count_list$Label[which(mod_count_list$Exclude =='x')])
img.count_mod<-count_images_renamed %>%
  dplyr::select( - img.mod_count_list$Label[which(img.mod_count_list$Exclude =='x')])
img.ann.count_mod<-count_annimages_renamed %>%
  dplyr::select( - img.mod_count_list$Label[which(img.mod_count_list$Exclude =='x')])
# #some of these categories no longer exist. Only need to exclude 'Tube"
# count_mod <- count_cells_renamed %>%
#   select(!Tube)

count_mod$cellID<-as.factor(count_mod$cellID)
img.count_mod$cellID<-as.factor(img.count_mod$cellID)
img.ann.count_mod$cellID<-as.factor(img.ann.count_mod$cellID)

##### 

## 6) Calculate abundances of functional groups and richness on the raw data (before species are excluded)
## first read in cell-metadata
load(paste0(ARC_Data.dir,"Cell_level_env_",res,".Rdata"))

## need to remove the cellID from the columns to properly calculate rowSums etc
cover_cells_renamed2 <- data.frame(cover_cells_renamed[,-1])
rownames(cover_cells_renamed2) <- cover_cells_renamed$cellID
  
## 6a) cover data
## names of faunal groups for cover_cells_renamed (minus cellID):
dataset.names <- names(cover_cells_renamed2)
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
scorable.pts <- cell_metadata_env$cover_points_scorable

cover_cells_pa <- cover_cells_renamed2
cover_cells_pa[cover_cells_pa>0] <- 1

cover_SF.prop <- rowSums(cover_cells_renamed2[,sel_SF])/scorable.pts
cover_SF <- rowSums(cover_cells_renamed2[,sel_SF])
cover_SF_pa <- cover_SF
cover_SF_pa[cover_SF>0] <- 1
richness <- rowSums(cover_cells_pa[,-sel_sed])
richness.l <- rowSums(cover_cells_pa[,-sel_sed])/log(scorable.pts)#(cell_metadata_env$cover_points_total)
cover_all.prop <- rowSums(cover_cells_renamed2[,-sel_sed])/scorable.pts
cover_all <- rowSums(cover_cells_renamed2[,-sel_sed])

cover_B.prop <- rowSums(cover_cells_renamed2[,sel_B])/scorable.pts
cover_B <- rowSums(cover_cells_renamed2[,sel_B])
cover_B_pa <- cover_B
cover_B_pa[cover_B>0] <- 1

cover_S.prop <- rowSums(cover_cells_renamed2[,sel_S])/scorable.pts
cover_S <- rowSums(cover_cells_renamed2[,sel_S])
cover_S_pa <- cover_S
cover_S_pa[cover_S>0] <- 1

## large species groupings - cover
cover_groupings <- data.frame(cbind(cover_all, cover_SF, cover_SF_pa, richness, richness.l, cover_B, cover_B_pa, cover_S, cover_S_pa))

##### setup biological data - counts
## names of faunal groups for count_cells:
count.names <- names(count_cells)

## selector for each faunal class
sel_noid <- grep("NoID",count.names)
sel_echino <- grep("Echinoderms",count.names)
sel_crust <- grep("Crustacea",count.names)

count_mobile <- rowSums(count_cells[,-sel_noid]) ## remove NoIDs
count_echino <- rowSums(count_cells[,sel_echino])
count_crust <- rowSums(count_cells[,sel_crust])

count_cells_pa <- count_cells
count_cells_pa[count_cells_pa>0] <- 1
count_richness <- rowSums(count_cells_pa[,-sel_sed])
count_richness.l <- rowSums(count_cells_pa[,-sel_sed])/log(cell_metadata_env$counts_N[!is.na(cell_metadata_env$counts_N)])

## large species groupings - counts
count_groupings <- data.frame(cbind(count_mobile, count_echino, count_crust, count_richness, count_richness.l))

#############
# load(paste0(ARC_Data.dir, "Cell_level_env_",res,".Rdata"))
# #join cover data back to cell metadata and environmental data
# cover_mod_env<-left_join(cell_metadata_env, cover_mod, by="cellID")
# #join count data back to cell metadata and environmental data
# count_mod_env<-left_join(cell_metadata_env, count_mod, by="cellID")
# save outputs
save(cover_mod, count_mod, cover_groupings, count_groupings, file=paste0(ARC_Data.dir,"Cell_level_bio_2pc_",res,".Rdata"))


##### FOR IMAGES INSTEAD OF CELLS:
## 6) Calculate abundances of functional groups and richness on the raw data (before species are excluded)
# ## first read in cell-metadata
load(paste0(ARC_Data.dir,"Cell_level_env_500m.Rdata"))

## need to remove the cellID from the columns to properly calculate rowSums etc
cover_images_renamed2 <- data.frame(cover_images_renamed[,-1])
rownames(cover_images_renamed2) <- cover_images_renamed$cellID

## 6a) cover data
## names of faunal groups for cover_cells_renamed (minus cellID):
dataset.names <- names(cover_images_renamed2)
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
scorable.pts <- 108-cover_images_renamed2$Unscorable

cover_images_pa <- cover_images_renamed2
cover_images_pa[cover_images_pa>0] <- 1

img.cover_SF.prop <- rowSums(cover_images_renamed2[,sel_SF])/scorable.pts
img.cover_SF <- rowSums(cover_images_renamed2[,sel_SF])
img.cover_SF_pa <- img.cover_SF
img.cover_SF_pa[img.cover_SF>0] <- 1
img.richness <- rowSums(cover_images_pa[,-sel_sed])
img.richness.l <- rowSums(cover_images_pa[,-sel_sed])/log(scorable.pts)
img.cover_all.prop <- rowSums(cover_images_renamed2[,-sel_sed])/scorable.pts
img.cover_all <- rowSums(cover_images_renamed2[,-sel_sed])

img.cover_B.prop <- rowSums(cover_images_renamed2[,sel_B])/scorable.pts
img.cover_B <- rowSums(cover_images_renamed2[,sel_B])
img.cover_B_pa <- img.cover_B
img.cover_B_pa[img.cover_B>0] <- 1

img.cover_S.prop <- rowSums(cover_images_renamed2[,sel_S])/scorable.pts
img.cover_S <- rowSums(cover_images_renamed2[,sel_S])
img.cover_S_pa <- img.cover_S
img.cover_S_pa[img.cover_S>0] <- 1

## large species groupings - cover
img.cover_groupings <- data.frame(cbind(img.cover_all, img.cover_SF, img.cover_SF_pa, img.richness, img.richness.l, img.cover_B, img.cover_B_pa, img.cover_S, img.cover_S_pa))

##### setup biological data - counts
## names of faunal groups for count_images:
count.names <- names(count_images)

## selector for each faunal class
sel_noid <- grep("NoID",count.names)
sel_echino <- grep("Echinoderms",count.names)
sel_crust <- grep("Crustacea",count.names)

img.count_mobile <- rowSums(count_images[,-sel_noid]) ## remove NoIDs
img.count_echino <- rowSums(count_images[,sel_echino])
img.count_crust <- rowSums(count_images[,sel_crust])

count_images_pa <- count_images
count_images_pa[count_images_pa>0] <- 1
img.count_richness <- rowSums(count_images_pa[,-sel_sed])

## large species groupings - counts
img.count_groupings <- data.frame(cbind(img.count_mobile, img.count_echino, img.count_crust, img.count_richness))

# save outputs, but first add metadata of images
sel.metadat <- match(rownames(cover_images),image_metadata$Filename.standardised)
img.metadata <- image_metadata[sel.metadat,]

save(img.cover_mod, img.count_mod, img.cover_groupings, img.count_groupings, img.metadata, file=paste0(ARC_Data.dir,"Image_level_bio.Rdata"))
###


####################################################
### species names for comparison with CATAMI
## cover at the image level:
img.ann.cover_prev2 <- data.frame(count=colSums(img.ann.cover_mod[,-1]>0)) %>%
  mutate(., prev= round(count/nrow(img.ann.cover_mod), 3 )) %>%
  rownames_to_column(., var="Label") %>%
  arrange(., desc(count))%>%
  add_column(CATAMI= "",
             CAAB= "")
cover.reordered <- order(img.ann.cover_prev2$Label)
img.ann.cover_prev3 <- img.ann.cover_prev2[cover.reordered,]

## counts at the image level:
img.ann.count_prev2<-data.frame(count=colSums(img.ann.count_mod[,-1]>0)) %>%
  mutate(., prev= round(count/nrow(img.ann.count_mod), 3 )) %>%
  rownames_to_column(., var="Label") %>%
  arrange(., desc(count)) %>%
  #add extra columns for notes
  add_column(CATAMI= "",
             CAAB= "")
counts.reordered <- order(img.ann.count_prev2$Label)
img.ann.count_prev3 <- img.ann.count_prev2[counts.reordered,]

# write_xlsx(x=list(COVER_naming=img.cover_prev3, COUNT_naming= img.count_prev3),
#           path=paste0(ARC_Data.dir, "Annotation/Species_list_images_vs_CATAMI_.xlsx"))





##




















































# #### old code previously used to calculate abundances of faunal groups
# 
# ## names of faunal groups for cover_cells_clean:
# dataset.names <- names(cover_cells_clean)
# ## selector for each faunal class
# sel_S <- grep("Sp",substr(dataset.names,1,2))
# sel_O <- grep("Oc",substr(dataset.names,1,2))
# sel_B <- grep("Bry",substr(dataset.names,1,3))
# sel_M <- grep("Mo",substr(dataset.names,1,2))
# sel_E <- grep("Ech",substr(dataset.names,1,3))
# sel_Asc <- grep("Asc",substr(dataset.names,1,3))
# sel_TW <- grep("Worms_Polychaetes_T",dataset.names)
# sel_Hy <- grep("Hyd",dataset.names)
# 
# sel_SF <- c(sel_S,sel_O,sel_B,sel_Asc,sel_TW,sel_Hy)
# 
# sel_sed_soft <- c(grep("Fine",dataset.names),grep("PbGrv",dataset.names))
# sel_sed_loose <- c(grep("Cbble",dataset.names),grep("BioRu",dataset.names),grep("BioShl",dataset.names),grep("BioOth",dataset.names)) 
# sel_sed_hard <- c(grep("Bould",dataset.names),grep("Rock",dataset.names))
# sel_sed <- grep("Sub_",dataset.names)
# sel_noid.cov <- grep("NoID",dataset.names)
# sel_unsc.cov <- grep("Unscorable",dataset.names)
# 
# cover_cells_clean_pa <- cover_cells_clean
# cover_cells_clean_pa[cover_cells_clean_pa>0] <- 1
# 
# cover_SF.prop <- rowSums(cover_cells_clean[,sel_SF])/cell_metadata_env_clean$cover_points_scorable
# cover_SF <- rowSums(cover_cells_clean[,sel_SF])
# cover_SF_pa <- cover_SF
# cover_SF_pa[cover_SF>0] <- 1
# richness <- rowSums(cover_cells_clean_pa[,-sel_sed])
# richness.l <- rowSums(cover_cells_clean_pa[,-sel_sed])/log(cell_metadata_env_clean$cover_points_total)
# cover_all.prop <- rowSums(cover_cells_clean[,-sel_sed])/cell_metadata_env_clean$cover_points_scorable
# cover_all <- rowSums(cover_cells_clean[,-sel_sed])
# 
# cover_B.prop <- rowSums(cover_cells_clean[,sel_B])/cell_metadata_env_clean$cover_points_scorable
# cover_B <- rowSums(cover_cells_clean[,sel_B])
# cover_B_pa <- cover_B
# cover_B_pa[cover_B>0] <- 1
# 
# cover_S.prop <- rowSums(cover_cells_clean[,sel_S])/cell_metadata_env_clean$cover_points_scorable
# cover_S <- rowSums(cover_cells_clean[,sel_S])
# cover_S_pa <- cover_S
# cover_S_pa[cover_S>0] <- 1
# 
# ##### setup biological data - counts
# 
# ## names of faunal groups for count_cells:
# count.names <- names(count_cells)
# 
# ## selector for each faunal class
# sel_noid <- grep("NoID",count.names)
# sel_echino <- grep("Echinoderms",count.names)
# sel_crust <- grep("Crustacea",count.names)
# 
# count_mobile <- rowSums(count_cells[,-sel_noid]) ## remove NoIDs
# count_echino <- rowSums(count_cells[,sel_echino])
# count_crust <- rowSums(count_cells[,sel_crust])
# 
# count_cells_pa <- count_cells
# count_cells_pa[count_cells_pa>0] <- 1
# count_richness <- rowSums(count_cells_pa[,-sel_sed])
# count_richness.l <- rowSums(count_cells_pa[,-sel_sed])/log(cell_metadata_env$counts_N[!is.na(cell_metadata_env$counts_N)])
# 
# ######################################################################################################
# ##### cover data
# ## individual species
# dat_cov_species <- cover_cells_clean[,-c(sel_sed, sel_noid.cov, sel_unsc.cov)]
# ## large species groupings - cover
# dat_cov_sum <- data.frame(cbind(cover_all, cover_SF, cover_SF_pa, richness, richness.l, cover_B, cover_B_pa, cover_S, cover_S_pa))
# ## presence-absence data
# dat_cov_pa <- cover_cells_clean_pa[,-c(sel_sed,sel_noid.cov,sel_unsc.cov)]
# 
# ##### count data
# ## individual species
# dat_count_species <- count_cells[,-sel_noid]
# ## large species groupings - counts
# dat_count_sum <- data.frame(cbind(count_mobile, count_echino, count_crust, count_richness, count_richness.l))
# ## presence-absence data
# dat_count_pa <- count_cells_pa[,-sel_noid]
# 
# ######################################################################################################
# 
# save(dat_cov_species, dat_cov_sum, dat_cov_pa,
#      dat_count_species, dat_count_sum, dat_count_pa,
#      file=paste0(biodiv.dir,"biodiversity_bio_dat.Rdata"))