##########
## WHAT THIS SCRIPT DOES:
## - loading and cleaning biological and environmental data
## - removing correlated environmental variables
## - scaling environmental data (at the sampling locations)
## - setting up derivatives of biological data, such as cover, richness etc, and saving for analysis later
## - scaling and adding polynomials to environmental rasters (circumpolar)
## - create a environmental dataframe with one row per cell (circumpolar), one column per variable
##########

##### Setting up----
library(PerformanceAnalytics) ## plotting correlations
library(terra)

user = "Jan"
#user = "charley"
#user="nicole"
if (user == "Jan") {
  sci.dir <-      "C:/Users/jjansen/OneDrive - University of Tasmania/science/"
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

biodiv.dir <- paste0(sci.dir,"SouthernOceanBiodiversityMapping/ARC_Benthic_Mapping/biodiversity/")

##############################################################################################################
##############################################################################################################
## Running this script for both data at 500m res and at 2km res

#res <- "500m"
res <- "2km"

######################################
##### load biological and environmental data
load(paste0(ARC_Data.dir,"annotation/Circumpolar_Annotation_Data_",res,"_202312.Rdata"))
## cell_metadata, count_cells, cover_cells
## image_metadata, count_images, cover_images

load(paste0(ARC_Data.dir,"annotation/Circumpolar_Annotation_Env_Data_",res,"_202412.Rdata"))
## cell_metadata_env, count_cells_env, cover_cells_env
## image_metadata_env

##############################################################################################################
## environmental data
##############################################################################################################

##### add information about number of unscorable points per cell
cell_metadata_env$cover_points_total <- rowSums(cover_cells)
cell_metadata_env$cover_points_scorable <- rowSums(cover_cells)-cover_cells$Unscorable

######################################################################################################
##### check correlations (essentially the same between 500m and 2km resolution)
## for 2km res:
## first batch:
chart.Correlation(cell_metadata_env[,24:37])
# remove tpi5, arag_mean, no3_mean & sd, po4_mean & sd
chart.Correlation(dplyr::select(cell_metadata_env[,24:37],-c("tpi5", "tpi11", "arag_mean", "no3_mean", "no3_sd", "po4_mean", "po4_sd")))

## second batch:
chart.Correlation(cell_metadata_env[,39:50])
#remove ice prop & mean & max, ice summer sd
chart.Correlation(dplyr::select(cell_metadata_env[,39:50],-c("ice_prop", "ice_mean", "ice_max", "ice_su_sd","npp_sd")))

## third batch:
chart.Correlation(cell_metadata_env[,51:62])
#remove ssh summer mean & sd, ssh spring mean, sst seasonal means, yearly sd and spring sd,  flux0001, flux0002
chart.Correlation(dplyr::select(cell_metadata_env[,51:62],-c("ssh_su_mean","ssh_su_sd","ssh_sp_mean","sst_sd","sst_sp_mean","sst_sp_sd","sst_su_mean")))#,"flux0001","flux0002")))

## fourth batch:
chart.Correlation(cell_metadata_env[,63:82])
#remove seafloorcurrents_absolute, individual flux and sed runs
chart.Correlation(dplyr::select(cell_metadata_env[,63:82],-c("seafloorcurrents_absolute","flux0001","flux0002","flux0005","flux00005","sed0001","sed0002","sed0005","sed00005")))

# ## together: note that npp_mean-susp08 are correlated, and tpi-tpi11
env.remove <- c("tpi5", "tpi11", "arag_mean", "no3_mean", "no3_sd", "po4_mean", "po4_sd",
                "ice_prop", "ice_mean", "ice_max", "ice_su_sd","npp_sd",
                "ssh_su_mean","ssh_su_sd","ssh_sp_mean","sst_sd","sst_sp_mean","sst_sp_sd","sst_su_mean",
                "seafloorcurrents_absolute","flux00005","flux0001","flux0002","flux0005","sed00005","sed0001","sed0002","sed0005")
env.sel.remove <- which(names(cell_metadata_env)%in%env.remove)
env.sel.remove.metadata <- c(1:23,83,84)
# chart.Correlation(cell_metadata_env[,-c(env.sel.remove.metadata,env.sel.remove,37)][,1:15])
# chart.Correlation(cell_metadata_env[,-c(env.sel.remove.metadata,env.sel.remove,37)][,16:29])

'%!in%' <- function(x,y)!('%in%'(x,y))
env.sel <- which(names(cell_metadata_env)%!in%env.remove)
sel.not.correlated <- (1:length(names(cell_metadata_env)))[-c(env.sel.remove.metadata,env.sel.remove)]

######################################################################################################
##### scale environmental data
cell_metadata_env_scaled <- cell_metadata_env
scale.means <- NA
scale.sd <- NA
sel.not.to.be.scaled <- c(env.sel.remove.metadata,38)
for(i in (1:ncol(cell_metadata_env_scaled))[-sel.not.to.be.scaled]){
  scale.means[i] <- mean(cell_metadata_env_scaled[,i], na.rm=TRUE)
  scale.sd[i] <- sd(cell_metadata_env_scaled[,i], na.rm=TRUE)
  cell_metadata_env_scaled[,i] <- (cell_metadata_env_scaled[,i]-scale.means[i])/scale.sd[i]
}

######################################################################################################
## to specify a spatial latent factor we need coordinates for each transect, calculated here:
transect.xy <- aggregate(image_metadata$proj_coord_x~image_metadata$transectID_full, FUN=mean)
transect.xy[,3] <- aggregate(image_metadata$proj_coord_y~image_metadata$transectID_full, FUN=mean)[,2]
names(transect.xy) <- c("transectID_full", "proj_coord_x", "proj_coord_y")

######################################################################################################
save(cell_metadata_env, transect.xy, sel.not.correlated,
     cell_metadata_env_scaled, scale.means, scale.sd,
     file=paste0(ARC_Data.dir,"Cell_level_env_",res,"_202412.Rdata"))


######################################################################################################
##### scaling rasters and preparing circumpolar cell data for predictions
######################################################################################################
## we only need res, env.derived and cell_metadata_env:
rm(list=setdiff(ls(), c("res","scale.means","scale.sd","env.derived","cell_metadata_env_scaled","sel.not.correlated","env.sel.remove.metadata")))

## get file names of all environmental rasters and bricks and load into one big stack----
pred_stack <- rast(c(paste0(env.derived,"Circumpolar_EnvData_",res,"_shelf_mask_unscaled_variables.tif"),
                    paste0(env.derived,"Circumpolar_EnvData_",res,"_shelf_mask_unscaled_polynomials_etc.tif")))

## matching up metadata with pred-layers
sel.vars <- names(cell_metadata_env_scaled)[sel.not.correlated]
sel.ra <- which(names(pred_stack)%in%sel.vars)
#names(pred_stack)[which(names(pred_stack)%in%sel.vars)]
pred_stack.nc <- subset(pred_stack,sel.ra)
names(pred_stack.nc)

## split scaling into runs of 4:
seq_split <- split(1:nlyr(pred_stack.nc), ceiling(seq_along(1:nlyr(pred_stack.nc))/4))
for(j in 1:length(seq_split)){
  ## creating an empty stack of 10 layers
  pred_stack_scaled <- rast(pred_stack.nc[[seq_split[[j]]]])
  for(i in 1:4){
    l <- seq_split[[j]][i]
    print(l)
    if(is.na(l)) break
    k <- names(pred_stack.nc)[l]
    ## we don't want to scale geomorphology
    if(k=="geomorphology"){
      pred_stack_scaled[[i]] <- pred_stack.nc[[l]]
    }else{
      ## select which layers match between the stack and the cell_metadata
      c.sel <- which(names(cell_metadata_env_scaled)==k)
      s.sel <- which(names(pred_stack.nc)==k)
      # pred_stack_scaled[[i]] <- rast(pred_stack.nc$depth)
      ## scale the raster
      pred_stack_scaled[[i]] <- (pred_stack.nc[[s.sel]]-scale.means[c.sel])/scale.sd[c.sel]
    }
  }
  writeRaster(pred_stack_scaled, filename=paste0(env.derived,"Circumpolar_EnvData_",res,"_shelf_mask_scaled_temporaryfile",j,".tif"), overwrite=TRUE)
}

## combine the temporary files into one
file_list<-list.files(path = env.derived, pattern="tif$",  full.names=TRUE) 
#subset to  "shelf" files
file_list<-file_list[grep(paste0(".",res,"_shelf_mask_scaled_temporary"), file_list)]
# ## we need to reorder the list so that #10 comes after 9...
#file_list <- file_list[c(1,3:10,2)]
## read in and save as one file
all_temporary_files <- rast(file_list)
writeRaster(all_temporary_files, filename=paste0(env.derived,"Circumpolar_EnvData_",res,"_shelf_mask_scaled.tif"), overwrite=TRUE)
#file.remove(file_list)

##############################################
## creating a dataframe with scaled values, one row per environmental cell
env_stack_scaled <- rast(paste0(env.derived,"Circumpolar_EnvData_",res,"_shelf_mask_scaled.tif"))
sel <- which(!is.na(env_stack_scaled$depth[]))

## res=="2km" takes about 2min
## res=="500m" takes about 1-2hrs
## split into groups of 3:
seq_split <- split(1:nlyr(env_stack_scaled), ceiling(seq_along(1:nlyr(env_stack_scaled))/3))
tempdat1 <- cbind(env_stack_scaled[[seq_split[[1]]]][sel])
tempdat2 <- cbind(env_stack_scaled[[seq_split[[2]]]][sel])
tempdat3 <- cbind(env_stack_scaled[[seq_split[[3]]]][sel])
tempdat4 <- cbind(env_stack_scaled[[seq_split[[4]]]][sel])
tempdat5 <- cbind(env_stack_scaled[[seq_split[[5]]]][sel])
tempdat6 <- cbind(env_stack_scaled[[seq_split[[6]]]][sel])
tempdat7 <- cbind(env_stack_scaled[[seq_split[[7]]]][sel])
tempdat8 <- cbind(env_stack_scaled[[seq_split[[8]]]][sel])
tempdat9 <- cbind(env_stack_scaled[[seq_split[[9]]]][sel])
tempdat10 <- cbind(env_stack_scaled[[seq_split[[10]]]][sel])
tempdat11 <- cbind(env_stack_scaled[[seq_split[[11]]]][sel])
# tempdat12 <- cbind(env_stack_scaled[[seq_split[[12]]]][sel])
# tempdat13 <- cbind(env_stack_scaled[[seq_split[[13]]]][sel])
# tempdat14 <- cbind(env_stack_scaled[[seq_split[[14]]]][sel])
# tempdat15 <- cbind(env_stack_scaled[[seq_split[[15]]]][sel])
# tempdat16 <- cbind(env_stack_scaled[[seq_split[[16]]]][sel])
# tempdat17 <- cbind(env_stack_scaled[[seq_split[[17]]]][sel])
pred_stack.dat <- data.frame(cbind(tempdat1,tempdat2,tempdat3,tempdat4,tempdat5,tempdat6,tempdat7,tempdat8,tempdat9,tempdat10,tempdat11))#,tempdat12,tempdat13,tempdat14,tempdat15,tempdat16,tempdat17))
pred_stack.dat$gear <- "OFOS"
pred_stack.dat$cover_cells_survey <- "PS96"
pred_stack.dat$cover_cells_transect1 <- "PS96_001"
save(pred_stack.dat, file=paste0(env.derived,"Circumpolar_EnvData_",res,"_shelf_mask_scaled_dataframe.Rdata"))


  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

# pred_stack.dat1 <- cbind(env_stack_scaled$depth[sel], env_stack_scaled$depth2[sel],
#                          env_stack_scaled$distance2canyons[sel], env_stack_scaled$distance2canyons2[sel])
# pred_stack.dat2 <- cbind(env_stack_scaled$logslope[sel], env_stack_scaled$slope[sel],
#                          env_stack_scaled$tpi[sel], env_stack_scaled$tpi11[sel])
# pred_stack.dat3 <- cbind(env_stack_scaled$tpi5[sel], env_stack_scaled$waom4k_test_flux08[sel],
#                          env_stack_scaled$waom4k_test_settle08[sel], env_stack_scaled$waom4k_test_susp08[sel])
# pred_stack.dat4 <- cbind(env_stack_scaled$waom4k_seafloorcurrents_absolute[sel], env_stack_scaled$waom4k_seafloorcurrents_mean[sel],
#                          env_stack_scaled$waom4k_seafloorcurrents_residual[sel])
# pred_stack.dat5 <- cbind(env_stack_scaled$waom4k_seafloorsalinity[sel], env_stack_scaled$waom4k_seafloortemperature[sel])
# 
# pred_stack.dat <- data.frame(cbind(pred_stack.dat1, pred_stack.dat2, pred_stack.dat3, pred_stack.dat4, pred_stack.dat5, 10))
# names(pred_stack.dat) <- c("depth", "depth2", "distance2canyons","distance2canyons2",
#                            "logslope","slope","tpi","tpi11",
#                            "tpi5", "waom4k_test_flux08","waom4k_test_settle08","waom4k_test_susp08",
#                            "waom4k_seafloorcurrents_absolute","waom4k_seafloorcurrents_mean","waom4k_seafloorcurrents_residual",
#                            "waom4k_seafloorsalinity","waom4k_seafloortemperature",
#                            "annotated_area")


















































# ## get file names of all environmental rasters and bricks and load into one big stack----
# #all files with "gri" extension
# env_list<-list.files(path = env.derived, pattern="gri$",  full.names=TRUE) 
# #subset to  "shelf" files
# env_list<-env_list[grep(".500m_shelf", env_list)]
# #for the single rasters layer names are missing. Extract from file name.
# env_names<-gsub(".*_|\\..*","",env_list)
# #stack all environmental layers and make sure they have appropriate names (currently manual and a bit messy!)
# env_stack<-stack(env_list)
# names(env_stack)
# names(env_stack)[1:6]<-env_names[1:6]
# names(env_stack)[15:23]<-paste(rep(c("CARS_NO3", "CARS_O2", "CARS_PO4"),each=3),c("mean", "seas_range", "std_dev"), sep="_")
# names(env_stack)[24] <-"distance2canyons"
# names(env_stack)[35]<-"NPP_su_mean"
# names(env_stack)[36:41]<-c("ssh_mean","ssh_sd","ssh_sp_mean","ssh_sp_sd","ssh_su_mean","ssh_su_sd")
# names(env_stack)[42:47]<-c("sst_mean","sst_sd","sst_sp_mean","sst_sp_sd","sst_su_mean","sst_su_sd")
# names(env_stack)[48:57]<-c("waom2k_seafloorcurrents", "waom2k_seafloortemperature", "waom4k_seafloorcurrents_absolute", "waom4k_seafloorcurrents_mean", 
#                            "waom4k_seafloorcurrents_residual", "waom4k_seafloorsalinity", "waom4k_seafloortemperature",
#                            "waom4k_test_flux08","waom4k_test_settle08","waom4k_test_susp08")
# 
# 
# ## only select rasters we actually need
# pred_stack <- raster::subset(env_stack, c(1,2,4,23,50,52,53,55))
# #logslope
# names(pred_stack)[2] <- "logslope"
# pred_stack$logslope <- log(env_stack$slope)
# #depth2
# pred_stack$depth2 <- pred_stack$depth
# pred_stack$depth2 <- raster(pred_stack$depth)
# sel <- which(!is.na(pred_stack$depth[]))
# depth2.dat <- poly(pred_stack$depth[sel],2)[,2] ## takes 10min or so
# pred_stack$depth2[sel] <- depth2.dat
# #dist2cany2
# pred_stack$distance2canyons2 <- pred_stack$distance2canyons
# pred_stack$distance2canyons2 <- raster(pred_stack$distance2canyons)
# sel <- which(!is.na(pred_stack$distance2canyons[]))
# distance2canyons2.dat <- poly(pred_stack$distance2canyons[sel],2)[,2] ## takes 10min or so
# pred_stack$distance2canyons2[sel] <- distance2canyons2.dat
# 
# 
# plot(pred_stack)
# 
# pred_stack_scaled <- pred_stack
# for(i in 1:nlayers(pred_stack_scaled)){
#   k <- names(pred_stack_scaled)[i]
#   c.sel <- which(names(cell_metadata_env_scaled)==k)
#   pred_stack_scaled[[i]] <- raster(pred_stack$depth)
#   pred_stack_scaled[[i]] <- (pred_stack[[i]]-scale.means[c.sel])/scale.sd[c.sel]
# }
# 
# plot(pred_stack_scaled)
# 
