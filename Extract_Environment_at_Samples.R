##### This code takes 

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


#stack all environmental layers and make sure they have appropriate names (bit messy!)
env_stack<-stack(env_list)
names(env_stack)[1:5]<-env_names[1:5]
names(env_stack)[15]<-"NPP_su_mean"
names(env_stack)[22:23]<-env_names[9:10]

#add environmental data with non-conformant names
env_stack<-stack( env_stack,
                  raster(paste0(env.derived, "Circumpolar_EnvData_geomorphology")))
names(env_stack)[24]<-"geomorph"



## 3) Match environmental data to image data (at cell level) ----
load(paste0(ARC_Data.dir, "Circumpolar_Annotation_Data.RData"))

# subset to only cells that have scored images
cell_metadata_env<- cell_metadata %>%
  filter(! is.na (cover_N))

#extract environmental data
cell_metadata_env<-cbind(cell_metadata_env, 
                         raster::extract(env_stack, cell_metadata_env[,c("proj_coord_x", "proj_coord_y")]))

## 4) Add file to Annotation data RData
save(cell_metadata, cell_metadata_env, count_cells, count_images, cover_cells, cover_images, 
           file =paste0(ARC_Data.dir, "Circumpolar_Annotation_Data.RData"))


cover_im_md<- cover_images %>% 
  rownames_to_column( var="Filename.standardised") %>%
  left_join( y=image_metadata) %>%
  write.csv(., file = paste0(path, "cover_image_md.csv" ))

#first import all files in a single folder as a list 
rastlist <- list.files(path = "/path/to/wd", pattern='.TIF$', all.files=TRUE, full.names=FALSE)

library(raster)
allrasters <- stack(rastlist)