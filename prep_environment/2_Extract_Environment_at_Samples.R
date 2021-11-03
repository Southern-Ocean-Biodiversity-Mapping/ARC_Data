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


sci.dir <-      "C:/Users/hillna/OneDrive - University of Tasmania/UTAS_work/Projects/Benthic Diversity ARC/"
env.derived <-  paste0(sci.dir,"data_environmental/derived/")
tools.dir <-    paste0(sci.dir,"Analysis/Useful_Functions_Tools/")
ARC_Data.dir <- paste0(sci.dir,"Analysis/ARC_Data/")


## 2) get file names of all environmental rasters and bricks and load into one big stack----
#all files with "gri" extension
env_list<-list.files(path = env.derived, pattern="gri$",  full.names=TRUE) 
#subset to  "shelf" files
env_list<-env_list[grep(".500m_shelf", env_list)]

#for the single raters layer names are missing. Extract from file name.
env_names<-gsub(".*_|\\..*","",env_list)


#stack all environmental layers and make sure they have appropriate names (currently manual and a bit messy!)
env_stack<-stack(env_list)
names(env_stack)
names(env_stack)[1:5]<-env_names[1:5]
names(env_stack)[6:14]<-paste(rep(c("CARS_NO3", "CARS_O2", "CARS_PO4"),each=3),c("mean", "seas_range", "std_dev"), sep="_")
names(env_stack)[25]<-"NPP_su_mean"
names(env_stack)[38:39]<-env_names[9:10]


#add environmental data with non-conformant names
env_stack<-stack( env_stack,
                  raster(paste0(env.derived, "Circumpolar_EnvData_geomorphology")))
names(env_stack)[40]<-"geomorph"

geomorph_cat<-levels(env_stack[[40]])[[1]]


## 3) Match environmental data to image data (at cell level) ----
#can run a image level too if needed
load(paste0(ARC_Data.dir, "annotation/Circumpolar_Annotation_Data.RData"))

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
                                      mutate(cellID=as.numeric(rownames(cover_cells))),  
                            by='cellID')

count_cells_env<- right_join(cell_metadata_env,
                            count_cells %>%
                                    mutate(cellID=as.numeric(rownames(count_cells))),  
                            by='cellID')


save(cell_metadata_env, cover_cells_env, count_cells_env,
           file =paste0(ARC_Data.dir, "annotation/Circumpolar_Annotation_Env_Data.RData"))


##########################################
