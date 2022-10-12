#######################################################################################################
##### This code extracts the derived environmental rasters at the location of annotated images.    ####
##### Currently matched at the 500m cell ID level.                                                 ####
#### Then merges the annoations with the cell level environmental data                              ####
#### Author Nicole Hill October 2021                                                               ####
#######################################################################################################


## 1) set up----
library(tidyverse)
library(raster)
library(rasterVis)
library(stringr)

#user = "Jan"
#user = "charley"
user="nicole"

if (user == "Jan") {
  
  sci.dir <-      "C:/Users/jjansen/Desktop/science/"
  env.derived <-  paste0(sci.dir,"data_environmental/derived/")
  #bio.dir <-      paste0(sci.dir,"data_biological/")
  
  ## remote repository (DOESN'T WORK YET):
  # env.dir <- "https://data.imas.utas.edu.au/data_transfer/admin/files/EnvironmentalData/"
  
  ## common paths (after "sci.dir")
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


## 2) get file names of all environmental rasters and bricks and load into one big stack----
#all files with "gri" extension
env_list<-list.files(path = env.derived, pattern="gri$",  full.names=TRUE) 
#subset to  "shelf" files
env_list<-env_list[grep(".500m_shelf", env_list)]

#for the single rasters layer names are missing. Extract from file name.
env_names<-gsub(".*_|\\..*","",env_list)


#stack all environmental layers and make sure they have appropriate names (currently manual and a bit messy!)
env_stack<-stack(env_list)
names(env_stack)
names(env_stack)[1:5]<-env_names[1:5]
names(env_stack)[14:22]<-paste(rep(c("CARS_NO3", "CARS_O2", "CARS_PO4"),each=3),c("mean", "seas_range", "std_dev"), sep="_")
names(env_stack)[23] <-"distance2canyons"
names(env_stack)[34]<-"NPP_su_mean"
names(env_stack)[35:40]<-c("ssh_mean","ssh_sd","ssh_sp_mean","ssh_sp_sd","ssh_su_mean","ssh_su_sd")
names(env_stack)[41:46]<-c("sst_mean","sst_sd","sst_sp_mean","sst_sp_sd","sst_su_mean","sst_su_sd")
names(env_stack)[47:56]<-c("waom2k_seafloorcurrents", "waom2k_seafloortemperature", "waom4k_seafloorcurrents_absolute", "waom4k_seafloorcurrents_mean", 
                           "waom4k_seafloorcurrents_residual", "waom4k_seafloorsalinity", "waom4k_seafloortemperature",
                           "waom4k_test_flux08","waom4k_test_settle08","waom4k_test_susp08")


#add environmental data with non-conformant names- 
#### remember to update column index if changes!!!
env_stack<-stack( env_stack,
                  raster(paste0(env.derived, "Circumpolar_EnvData_geomorphology")))
names(env_stack)[57]<-"geomorph"

geomorph_cat<-levels(env_stack[[57]])[[1]]


## 3) Match environmental data to image data (at cell level) ----
#can run a image level too if needed
load(paste0(ARC_Data.dir, "annotation/Circumpolar_Annotation_Data.RData"))


 #for some reason nearly every column of cell_metadata are now character strings
 #cell_metadata<- cell_metadata %>% 
#   mutate(across(cellID:counts_area, as.numeric))%>%
#   mutate(across(cover_cells_transect1:cover_cells_transect3, as.numeric))%>%
#   mutate(across(counts_cells_transect1:counts_cells_transect3, as.numeric))


# subset to only cells that have scored images
cell_metadata_env<- cell_metadata %>%
  filter(! is.na (cover_N))

#extract environmental data
cell_metadata_env<-cbind(cell_metadata_env, 
                         raster::extract(env_stack, cell_metadata_env[,c("proj_coord_x", "proj_coord_y")]))
#add geomorph name
cell_metadata_env<-cell_metadata_env %>%
  left_join(., geomorph_cat, by=c("geomorph"= "ID"))
 
cell_metadata_env<- rename(cell_metadata_env, geomorph_cat=VALUE)


## 4) Merge annotation save combined Annotation and environmental data as RData file
cover_cells_env<- left_join(cell_metadata_env,
                            cover_cells %>%
                              mutate(cellID=as.factor(rownames(cover_cells))),  
                              #mutate(cellID=as.numeric(rownames(cover_cells))),
                            by='cellID')

count_cells_env<- right_join(cell_metadata_env,
                            count_cells %>%
                            mutate(cellID=as.factor(rownames(count_cells))),
                            #mutate(cellID=as.numeric(rownames(count_cells))),
                            by='cellID')


save(cell_metadata_env, cover_cells_env, count_cells_env,
           file =paste0(ARC_Data.dir, "annotation/Circumpolar_Annotation_Env_Data.RData"))


##########################################
