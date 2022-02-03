###################################################################################################
#### This code converts removes 'unscorable' points and calculates prevalence on overall cover  ###
### This forms initial species list to make decisions about which species to keep/merge/exclude ###
###                                                                                             ###
### N.Hill- Modified Jan 2022                                                                   ###
###################################################################################################



# 1) libraries and paths 
library(tidyverse)
library(writexl)

sci.dir <-      "C:/Users/hillna/OneDrive - University of Tasmania/UTAS_work/Projects/Benthic Diversity ARC/"
env.derived <-  paste0(sci.dir,"data_environmental/derived/")
tools.dir <-    paste0(sci.dir,"Analysis/Useful_Functions_Tools/")
ARC_Data.dir <- paste0(sci.dir,"Analysis/ARC_Data/")


load(paste0(ARC_Data.dir, "annotation/Circumpolar_Annotation_Data.RData"))
#load(paste0(ARC_Data.dir, "annotation/Circumpolar_Annotation_Env_Data.RData"))


## 2) COVER ----
## 2a) Prevalence
cover_prev<-data.frame(count=colSums(cover_cells>0)) %>%
  mutate(., prev= round(count/950, 3 )) %>%
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
  mutate(., prev= round(count/891, 3 )) %>%
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

write_xlsx(x= list(COVER_prevalence=cover_prev,
          COVER_overall=cover_pc_overall,
          COUNT_prevalence= count_prev,
          COUNT_totAb= count_ab_overall),
          path=paste0(ARC_Data.dir, "Annotation/Species_list.xlsx"))
