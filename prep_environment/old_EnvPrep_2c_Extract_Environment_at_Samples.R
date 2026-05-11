#######################################################################################################
##### This code extracts the derived environmental rasters at the location of annotated images.    ####
##### matched at the 2km cell ID level.                                                            ####
##### Then merges the annotations with the cell level environmental data                            ####
#######################################################################################################

## specify user and setup directory to look up data from
usr <- "VM"
#usr <- "SJ"
#usr <- "JJ
source("0_SourceFile.R")

## set folders
env.derived  <- paste0(usr.main.dir, "data_environmental/derived/")

## 1) set up----
library(tidyverse)
library(terra)
library(rasterVis)
library(stringr)

################
#res <- "500m"
res <- "2km"
##############

## 2) get file names of all environmental rasters and bricks and load into one big stack----
env_stack <- rast(c(paste0(env.derived,"Circumpolar_EnvData_",res,"_shelf_mask_unscaled_variables.tif"),
                       paste0(env.derived,"Circumpolar_EnvData_",res,"_shelf_mask_unscaled_polynomials_etc.tif")))

## 3) Match environmental data to image data (at cell level) ----
#can run a image level too if needed
load(paste0(ARC_Data.dir, "annotation/Circumpolar_Annotation_Data_",res,"_202312.RData"))


 #for some reason nearly every column of cell_metadata are now character strings
 #cell_metadata<- cell_metadata %>% 
#   mutate(across(cellID:counts_area, as.numeric))%>%
#   mutate(across(cover_cells_transect1:cover_cells_transect3, as.numeric))%>%
#   mutate(across(counts_cells_transect1:counts_cells_transect3, as.numeric))


# subset to only cells that have scored images
cell_metadata_env<- cell_metadata %>%
  filter(! is.na (cover_N))

#extract environmental data
# cell_metadata_env<-cbind(cell_metadata_env, 
#                          raster::extract(env_stack, cell_metadata_env[,c("proj_coord_x", "proj_coord_y")]))
cell_metadata_env<-cbind(cell_metadata_env, 
                         terra::extract(env_stack, cell_metadata_env[,c("proj_coord_x", "proj_coord_y")])[,-1])

# #add geomorph name
# cell_metadata_env<-cell_metadata_env %>%
#   left_join(., geomorph_cat, by=c("geomorph"= "ID"))
#  
# cell_metadata_env<- rename(cell_metadata_env, geomorph_cat=VALUE)


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
           file =paste0(ARC_Data.dir, "annotation/Circumpolar_Annotation_Env_Data_",res,"_202412.RData"))


##########################################
## check where NAs are

#### ice:
plot(env_stack$ice_mean)
points(cell_metadata_env[,4:5])
points(cell_metadata_env[which(is.na(cell_metadata_env$ice_mean)),4:5], col="blue", pch=16)

## Ross Sea: cutoff from mask/depth-range, seamounts are gone
plot(env_stack$ice_mean, xlim=c(-500000,500000), ylim=c(-2500000,-2000000))
points(cell_metadata_env[,4:5])
points(cell_metadata_env[which(is.na(cell_metadata_env$ice_mean)),4:5], col="blue")

## Mertz: Satellite-land-mask?
plot(env_stack$ice_mean, xlim=c(1000000,2000000), ylim=c(-2500000,-2000000))
points(cell_metadata_env[,4:5])
points(cell_metadata_env[which(is.na(cell_metadata_env$ice_mean)),4:5], col="blue", pch=16)

## AP: Satellite-land-mask?
plot(env_stack$ice_mean, xlim=c(-2700000,-2200000), ylim=c(1100000,1900000))
points(cell_metadata_env[,4:5])
points(cell_metadata_env[which(is.na(cell_metadata_env$ice_mean)),4:5], col="blue", pch=16)

## Lazarev: Satellite-land-mask?
plot(env_stack$ice_mean, xlim=c(-700000,500000), ylim=c(1500000,2300000))
points(cell_metadata_env[,4:5])
points(cell_metadata_env[which(is.na(cell_metadata_env$ice_mean)),4:5], col="blue", pch=16)



