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
  mutate(., prev= round(count/nrow(count_cells), 3 )) %>%
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
#           path=paste0(ARC_Data.dir, "Annotation/Species_list_",res,"_2022_10.xlsx"))


#######################################################################

## 5) Read commented file back in and make exclusions and aggregations
# blanks read in as NAs

## 5a) cover data
mod_cover_list<-read_xlsx(path=paste0(ARC_Data.dir, "Annotation/Species_list_",res,"_2022_10.xlsx"),
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


## 5b) count data
mod_count_list<-read_xlsx(path=paste0(ARC_Data.dir, "Annotation/Species_list_",res,"_2022_10.xlsx"),
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
 dplyr::select( - mod_count_list$Label[which(mod_count_list$Exclude =='x')])
# #some of these categories no longer exist. Only need to exclude 'Tube"
# count_mod <- count_cells_renamed %>%
#   select(!Tube)

count_mod$cellID<-as.factor(count_mod$cellID)

##### 

## 6) Calculate abundances of functional groups and richness on the raw data (before species are excluded)
## first read in cell-metadata
load(paste0(ARC_Data.dir,"Cell_level_env_",res,".Rdata"))

## 6a) cover data
## names of faunal groups for cover_cells_renamed:
dataset.names <- names(cover_cells_renamed)
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

## calculating species richness
scorable.pts <- cell_metadata_env$cover_points_scorable

cover_cells_pa <- cover_cells_renamed
cover_cells_pa[cover_cells_pa>0] <- 1

cover_SF.prop <- rowSums(cover_cells_renamed[,sel_SF])/scorable.pts
cover_SF <- rowSums(cover_cells_renamed[,sel_SF])
cover_SF_pa <- cover_SF
cover_SF_pa[cover_SF>0] <- 1
richness <- rowSums(cover_cells_pa[,-sel_sed])
richness.l <- rowSums(cover_cells_pa[,-sel_sed])/log(cell_metadata_env$cover_points_total)
cover_all.prop <- rowSums(cover_cells_renamed[,-sel_sed])/scorable.pts
cover_all <- rowSums(cover_cells_renamed[,-sel_sed])

cover_B.prop <- rowSums(cover_cells_renamed[,sel_B])/scorable.pts
cover_B <- rowSums(cover_cells_renamed[,sel_B])
cover_B_pa <- cover_B
cover_B_pa[cover_B>0] <- 1

cover_S.prop <- rowSums(cover_cells_renamed[,sel_S])/scorable.pts
cover_S <- rowSums(cover_cells_renamed[,sel_S])
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