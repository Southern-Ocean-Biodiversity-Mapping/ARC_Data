###################################################################################################
#### This code converts removes 'unscorable' points and calculates prevalence on overall cover  ###
### This forms initial species list to make decisions about which species to keep/merge/exclude ###
###                                                                                             ###
### file with comments is read back in and exclusions/ aggregations made to form final datafile ###
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

load(paste0(ARC_Data.dir, "annotation/Circumpolar_Annotation_Data.RData"))


## 2) COVER ----
## 2a) Prevalence
cover_prev<-data.frame(count=colSums(cover_cells>0)) %>%
  mutate(., prev= round(count/961, 3 )) %>%
  rownames_to_column(., var="Label") %>%
  arrange(., desc(count))%>%
  add_column(Exclude= "",
             Merge_With= "",
             Other_Notes="")


##2b) Percent cover
# remove unscorable points
cover_cells2<-cover_cells %>%
  dplyr::select(., -Unscorable) 

# percent cover of taxa in each cell
cover_pc<-prop.table(as.matrix(cover_cells ), margin=1)*100 
  
# proportion of all scored points assigned to each taxa
cover_pc_overall<- data.frame(perc_overall=colSums(prop.table(as.matrix(cover_cells))*100)) %>%
  rownames_to_column(., var="Taxa") %>%
  arrange(., desc(perc_overall))%>%
  mutate(perc_overall= round(perc_overall, 3)) 

## average percent cover across all cells
cover_pc_mean<- data.frame(perc_mean= apply(cover_pc, 2, mean))%>%
rownames_to_column(., var="Taxa") %>%
  arrange(., desc(perc_mean))%>%
  mutate(perc_mean= round(perc_mean, 3)) 
#very similar to overall percent cover

## 3) COUNTS ----
## 3a) Prevalence
  
count_prev<-data.frame(count=colSums(count_cells>0)) %>%
  mutate(., prev= round(count/897, 3 )) %>%
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


## 4) Save to excel to comments and notes on aggregation etc.

# write_xlsx(x= list(COVER_prevalence=cover_prev,
#           COVER_overall=cover_pc_overall,
#           COUNT_prevalence= count_prev,
#           COUNT_totAb= count_ab_overall),
#           path=paste0(ARC_Data.dir, "Annotation/Species_list_2022_10_06.xlsx"))


#######################################################################

## 5) Read commented file back in and make exclusions and aggregations
# blanks read in as NAs
load(paste0(ARC_Data.dir, "annotation/Circumpolar_Annotation_Env_Data.RData"))


## 5a) cover data
mod_cover_list<-read_xlsx(path=paste0(ARC_Data.dir, "Annotation/Species_list_2022_10_07.xlsx"),
sheet=1)
names(mod_cover_list)[5:6]<- c("Merge_1pc", "Merge_2pc")


#reformat to long, merge, change names and convert back to wide
cover_cells_long<-cover_cells %>%
  mutate(cellID=as.numeric(rownames(cover_cells)))  %>% #add cellID
  pivot_longer( cols=`1Sub_Fine`:Echinoderms_Crinoids_Stalked, 
                names_to ="Label", values_to = "count") %>%           #long format and merge names to change
  left_join( mod_cover_list[,c("Label", "Merge_2pc")])

cover_cells_long$new<-  ifelse(!is.na(cover_cells_long$Merge_2pc), cover_cells_long$Merge_2pc, cover_cells_long$Label)

cover_cells_renamed<-pivot_wider(cover_cells_long, id_cols=cellID, names_from = new, values_from = count, values_fn=sum,values_fill = 0)

#remove species to exclude (DON'T REMOVE BECAUSE WE MESS UP TOTAL COVER DATA?)
cover_mod<-cover_cells_renamed
# cover_mod<-cover_cells_renamed %>%
#   select( - mod_cover_list$Label[which(mod_cover_list$Exclude =='x')])
cover_mod$cellID<-as.factor(cover_mod$cellID)

#join back to cell metadata and environmental data

cover_mod_env<-left_join(cell_metadata_env, cover_mod, by="cellID")


## 5b) count data
mod_count_list<-read_xlsx(path=paste0(ARC_Data.dir, "Annotation/Species_list_2022_10_07.xlsx"),
                          sheet=3)
names(mod_count_list)[5:6]<- c("Merge_1pc", "Merge_2pc")

#reformat to long, merge, change names and convert back to wide
count_cells_long<-count_cells %>%
  mutate(cellID=as.numeric(rownames(count_cells)))  %>% #add cellID
  pivot_longer( cols=Echinoderms__Crinoid_unstalked:Molluscs__Gastropods_shell__LimpetLike, 
                        names_to ="Label", values_to = "count") %>%           #long format and merge names to change
  left_join( mod_count_list[,c("Label", "Merge_2pc")])

count_cells_long$new<-  ifelse(!is.na(count_cells_long$Merge_2pc), count_cells_long$Merge_2pc, count_cells_long$Label)

count_cells_renamed<-pivot_wider(count_cells_long, id_cols=cellID, names_from = new, values_from = count, values_fn=sum,values_fill = 0)

#remove species to exclude
count_mod<-count_cells_renamed %>%
 select( - mod_count_list$Label[which(mod_count_list$Exclude =='x')])
# #some of these categories no longer exist. Only need to exclude 'Tube"
# count_mod <- count_cells_renamed %>%
#   select(!Tube)

count_mod$cellID<-as.factor(count_mod$cellID)

#join back to cell metadata and environmental data
count_mod_env<-left_join(cell_metadata_env, count_mod, by="cellID")



# save outputs
save(cover_mod, cover_cells_env, count_mod, count_cells_env, file=paste0(ARC_Data.dir,"Cell_level_bioenv_2pc.RData"))
