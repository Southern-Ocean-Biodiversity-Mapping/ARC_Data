#############################################################
## Presence- Absence image-level data                      ##
## restricted to 0.5% prevalent species                    ##
## matched to @ 2km resolution  scaled environmental data  ##
## For input into HMSC and RCP analyses                    ##
#############################################################


####################
#### 1) set up----
####################
library(dplyr)
library(tidyr)
#library(raster)
library(terra)
#devtools::install_github('skiptoniam/ecomix')
library(ecomix)
library(ggplot2)
library(forcats)
library(ggpubr)
library(RColorBrewer)


# load same files as Jan and Charley using for diversity and VME analysis
#comp= "nicole"
comp= "vm"

if(comp=="nicole"){
  path<- "C:\\Users\\hillna\\OneDrive - University of Tasmania\\UTAS_work\\Projects\\Benthic Diversity ARC\\Analysis\\ARC_Data\\"
  
  # biological data at 2km cells and 2pc prevalence
  #load(paste0(path, "Cell_level_bio_2pc_2km.RData"))
  #load(paste0(path, "Cell_level_env_2km.Rdata"))
  # add path for environmental data
}

if(comp=="vm"){
  path="/perm_storage/shared_space/BioMAS/"
  
  # biological data at 2km cells and 2pc prevalence
  #load(paste0(path, "ARC_Data/Cell_level_bio_2pc_2km.Rdata"))
  load(paste0(path, "ARC_Data/Image_level_bio.Rdata"))
  load(paste0(path, "ARC_Data/Cell_level_env_2km.Rdata"))
  load(paste0(path,"environmental_data/Circumpolar_Coastline.Rdata"))
  source( paste0(path, "ARC_Benthic_Mapping/bioregions/RCP_Helper_funcs_CircPolar.R"))
}

######################################################
#### 2) Format biological and environmental data ----
######################################################

## 2a) remove NA's and seamount transects (rows)
#note: matching IDs at image level are:
# img.metadata$Filename.standardised & cover_mod$CellID

#merge image metadata and cover info
img_cov<-left_join(img.metadata, cover_mod, by=join_by(Filename.standardised==cellID))

## remove seamount transects (tan1802 & tan1901)
seamount.transects <- c("TAN1802_160","TAN1802_170","TAN1802_179","TAN1802_180","TAN1802_184","TAN1802_185","TAN1802_191","TAN1802_193",
                        "TAN1802_195","TAN1802_196","TAN1802_197","TAN1802_207","TAN1802_208","TAN1802_209","TAN1802_213","tan1901_209")


## remove seamount transects
img_cov <- img_cov %>% 
  filter(! transectID_full %in% seamount.transects)

#join to scaled environmental data via img_cov$CellID (which was CellID originally in cover_mod) assume this is 2km resolution?


## 2b)Biological data: remove and combine some taxa, convert to PA, remove rare species (cols)

#names(img_cov)
#remove substrate categories, No ID and unscorable
rem<-c("Sand / Mud", "Pebble / Gravel" ,"Boulders" , "Cobbles" ,"Biologenic Rubble","Rock" , 
       "Unidentifiable","Unscorable")

cover_ims_PA<- img_cov [,18:115]%>%
  dplyr::select(., ! all_of (rem)) 

## combine UBS_B with Bryozoan_Hard_Branching_Antler
cover_ims_PA$`Bryozoa - Hard - Branching - Morphotype 1 - Antler` <- cover_ims_PA$`Bryozoa - Hard - Branching - Morphotype 1 - Antler`+cover_ims_PA$`Unidentified Biological Matrix - Bryozoan associated`
cover_ims_PA<- cover_ims_PA %>%
  dplyr::select( ., - `Unidentified Biological Matrix - Bryozoan associated`)

#Combine Unidentified biological matrix to Hydroid Matrix
cover_ims_PA$`Hydroid Matrix`<- cover_ims_PA$`Hydroid Matrix` + cover_ims_PA$`Unidentified Biological Matrix`
cover_ims_PA<- cover_ims_PA %>%
  dplyr::select( ., -`Unidentified Biological Matrix` ) 

#combine 2 taxa that are same but slightly different name
cover_ims_PA$`Sponge - Massive forms - Simple - Other`<-cover_ims_PA$`Sponge - Massive forms - Simple - Other` + 
          cover_ims_PA$`Sponge - Massive forms - Simple -Other`
cover_ims_PA<- cover_ims_PA %>%
  dplyr::select( ., -`Sponge - Massive forms - Simple -Other` ) 

# convert bio to presence-absence data
cover_ims_PA[cover_ims_PA >0] <- 1

#check rare species- #3427 images
prev<-colSums(cover_ims_PA)
summary(prev)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00   37.75   76.00  126.39  150.25  888.00 

#remove species with few occurrences (17 or less than 0.5%)
prev<-names(prev[prev>17])

cover_ims_PA<-cover_ims_PA %>%
  dplyr::select(., all_of(prev))

#add back in 2km CellID and then use to match to environmental data
cover_ims_PA<-bind_cols(img_cov[,1:17], cover_ims_PA)


#vars <- c("cellID","cover_cells_survey", "cover_N",
#          "depth","depth2","logslope","tpi","distance2canyons","distance2canyons2",
#          "seafloortemperature","seafloorcurrents_mean","seafloorsalinity","npp_mean")

#vars2 <- c("cellID","cover_cells_survey", "cover_N",
#          "depth","depth2","logslope","tpi","distance2canyons","distance2canyons2",
#          "seafloortemperature","seafloorcurrents_mean","seafloorsalinity","npp_mean", "sst_mean", 
#          "sst_sd", "sst_su_mean", "sst_su_sd", "ice_su_mean", "ice_su_sd", "o2_mean", "o2_sd")

#keep the following non-highly correlated vars
vars3<-c("cellID", "depth","depth2","logslope","tpi","distance2canyons","distance2canyons2",
         "seafloortemperature","seafloorcurrents_mean","seafloorsalinity","npp_mean", 
         "sst_su_mean", "ice_su_mean")

cover_ims_PA$cellID_2km <-as.factor(cover_ims_PA$cellID_2km)

bioenv<-left_join(cover_ims_PA, cell_metadata_env_scaled[,vars3], by=join_by(cellID_2km==cellID))
#note there will be some NAs in the environmental data due to the inclusion of sst and ice, which I am investigating in the RCP analyses.
save(bioenv, file = paste0(path, "ARC_Data/BioEnv_Data/Bioenv_imgPA_005pc_2km.RData"))
