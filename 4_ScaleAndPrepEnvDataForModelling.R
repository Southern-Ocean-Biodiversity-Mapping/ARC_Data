##########
## WHAT THIS SCRIPT DOES:
## - loading and cleaning biological and environmental data
## - scaling environmental data (at the sampling locations)
## - setting up derivatives of biological data, such as cover, richness etc, and saving for analysis later
## - scaling and adding polynomials to environmental rasters (circumpolar)
## - create a environmental dataframe with one row per cell (circumpolar), one column per variable
##########

##### Setting up----
library(PerformanceAnalytics) ## plotting correlations
user = "Jan"
#user = "charley"
#user="nicole"
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
##############################################################################################################
##############################################################################################################

res <- "500m"
# res <- "2km"


######################################

## Running this scipt for both data at 500m res and at 2km res

##### load biological and environmental data
load(paste0(ARC_Data.dir,"annotation/Circumpolar_Annotation_Data.Rdata"))
## cell_metadata, count_cells, cover_cells
## image_metadata, count_images, cover_images

load(paste0(ARC_Data.dir,"annotation/Circumpolar_Annotation_Env_Data_",res,".RData"))
## cell_metadata_env, count_cells_env, cover_cells_env
## image_metadata_env

##############################################################################################################
##### CONSIDERING ONLY CELL DATA: the only dataframes we really need are:
## cell_metadata_env
## cover_cells, count_cells
##############################################################################################################

## find NAs in waom data
waom.na.sel <- which(!is.na(cell_metadata_env$seafloorcurrents_absolute))

cell_metadata_env_clean <- cell_metadata_env[waom.na.sel,]
cover_cells_clean <- cover_cells[waom.na.sel,]

##### add information about number of unscorable points per cell
cell_metadata_env_clean$cover_points_total <- rowSums(cover_cells_clean)
cell_metadata_env_clean$cover_points_scorable <- rowSums(cover_cells_clean)-cover_cells_clean$Unscorable

##### scale environmental data
cell_metadata_env_clean_scaled <- cell_metadata_env_clean
scale.means <- NA
scale.sd <- NA
for(i in (1:ncol(cell_metadata_env_clean_scaled))[-c(1:21,68,69,78,79)]){
  scale.means[i] <- mean(cell_metadata_env_clean_scaled[,i], na.rm=TRUE)
  scale.sd[i] <- sd(cell_metadata_env_clean_scaled[,i], na.rm=TRUE)
  cell_metadata_env_clean_scaled[,i] <- (cell_metadata_env_clean_scaled[,i]-scale.means[i])/scale.sd[i]
}


















##### setup biological data - cover

## images per cell
cell_metadata_env$cover_N

## names of faunal groups for cover_cells_clean:
dataset.names <- names(cover_cells_clean)
## selector for each faunal class
sel_S <- grep("Sp",substr(dataset.names,1,2))
sel_O <- grep("Oc",substr(dataset.names,1,2))
sel_B <- grep("Bry",substr(dataset.names,1,3))
sel_M <- grep("Mo",substr(dataset.names,1,2))
sel_E <- grep("Ech",substr(dataset.names,1,3))
sel_Asc <- grep("Asc",substr(dataset.names,1,3))
sel_TW <- grep("Worms_Polychaetes_T",dataset.names)
sel_Hy <- grep("Hyd",dataset.names)

sel_SF <- c(sel_S,sel_O,sel_B,sel_Asc,sel_TW,sel_Hy)

sel_sed_soft <- c(grep("Fine",dataset.names),grep("PbGrv",dataset.names))
sel_sed_loose <- c(grep("Cbble",dataset.names),grep("BioRu",dataset.names),grep("BioShl",dataset.names),grep("BioOth",dataset.names)) 
sel_sed_hard <- c(grep("Bould",dataset.names),grep("Rock",dataset.names))
sel_sed <- grep("Sub_",dataset.names)
sel_noid.cov <- grep("NoID",dataset.names)
sel_unsc.cov <- grep("Unscorable",dataset.names)

cover_cells_clean_pa <- cover_cells_clean
cover_cells_clean_pa[cover_cells_clean_pa>0] <- 1

cover_SF.prop <- rowSums(cover_cells_clean[,sel_SF])/cell_metadata_env_clean$cover_points_scorable
cover_SF <- rowSums(cover_cells_clean[,sel_SF])
cover_SF_pa <- cover_SF
cover_SF_pa[cover_SF>0] <- 1
richness <- rowSums(cover_cells_clean_pa[,-sel_sed])
richness.l <- rowSums(cover_cells_clean_pa[,-sel_sed])/log(cell_metadata_env_clean$cover_points_total)
cover_all.prop <- rowSums(cover_cells_clean[,-sel_sed])/cell_metadata_env_clean$cover_points_scorable
cover_all <- rowSums(cover_cells_clean[,-sel_sed])

cover_B.prop <- rowSums(cover_cells_clean[,sel_B])/cell_metadata_env_clean$cover_points_scorable
cover_B <- rowSums(cover_cells_clean[,sel_B])
cover_B_pa <- cover_B
cover_B_pa[cover_B>0] <- 1

cover_S.prop <- rowSums(cover_cells_clean[,sel_S])/cell_metadata_env_clean$cover_points_scorable
cover_S <- rowSums(cover_cells_clean[,sel_S])
cover_S_pa <- cover_S
cover_S_pa[cover_S>0] <- 1

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

######################################################################################################
##### cover data
## individual species
dat_cov_species <- cover_cells_clean[,-c(sel_sed, sel_noid.cov, sel_unsc.cov)]
## large species groupings - cover
dat_cov_sum <- data.frame(cbind(cover_all, cover_SF, cover_SF_pa, richness, richness.l, cover_B, cover_B_pa, cover_S, cover_S_pa))
## presence-absence data
dat_cov_pa <- cover_cells_clean_pa[,-c(sel_sed,sel_noid.cov,sel_unsc.cov)]

##### count data
## individual species
dat_count_species <- count_cells[,-sel_noid]
## large species groupings - counts
dat_count_sum <- data.frame(cbind(count_mobile, count_echino, count_crust, count_richness, count_richness.l))
## presence-absence data
dat_count_pa <- count_cells_pa[,-sel_noid]

######################################################################################################
## to specify a spatial latent factor we need coordinates for each transect, calculated here:
transect.xy <- aggregate(image_metadata$proj_coord_x~image_metadata$transectID_full, FUN=mean)
transect.xy[,3] <- aggregate(image_metadata$proj_coord_y~image_metadata$transectID_full, FUN=mean)[,2]
names(transect.xy) <- c("transectID_full", "proj_coord_x", "proj_coord_y")

######################################################################################################
## check correlations

# chart.Correlation(cell_metadata_env_clean[,21:57])
# #first remove sst & seasonal ssh
# chart.Correlation(cell_metadata_env_clean[,21:47])
# #remove ice mean, ice spring mean & max, ice summer mean & sd, npp mean
# chart.Correlation(cell_metadata_env_clean[,c(21:35,37,38,41,43,46,47)])
# #remove tpi5, arag_mean, no3_mean & sd, po4_mean & sd
# chart.Correlation(cell_metadata_env_clean[,c(21:24,27,30,31,34,35,37,38,41,43,46,47)])
# #remove flux, 2k-temperature & 2k-currents
# chart.Correlation(cell_metadata_env_clean[,c(60:64,66,67,70:77)])
# 
# chart.Correlation(cell_metadata_env_clean[,c(21:24,27,30,31,34,35,37,38,41,43,46,47,60:64,66,67,70:77)])

sel.not.correlated <- c(21:24,27,30,31,34,35,37,38,41,43,46,47,60:64,66,67,70:77)
######################################################################################################
biodiv.dir <- paste0(sci.dir,"SouthernOceanBiodiversityMapping/ARC_Benthic_Mapping/biodiversity/")

save(dat_cov_species, dat_cov_sum, dat_cov_pa,
     dat_count_species, dat_count_sum, dat_count_pa,
     file=paste0(biodiv.dir,"biodiversity_bio_dat.Rdata"))

save(cell_metadata_env_clean, transect.xy, sel.not.correlated,
     cell_metadata_env_clean_scaled, scale.means, scale.sd,
     file=paste0(biodiv.dir,"biodiversity_env_dat.Rdata"))

######################################################################################################
##### scaling rasters and preparing circumpolar cell data for predictions
######################################################################################################
load(file="biodiversity_bio_dat.Rdata")
load(file="biodiversity_env_dat.Rdata")

## get file names of all environmental rasters and bricks and load into one big stack----
#all files with "gri" extension
env_list<-list.files(path = env.derived, pattern="gri$",  full.names=TRUE) 
#subset to  "shelf" files
env_list<-env_list[grep(".500m_shelf", env_list)]
env_list<-env_list[-grep(".500m_shelf_scaled", env_list)]
#for the single rasters layer names are missing. Extract from file name.
env_names<-gsub(".*_|\\..*","",env_list)
#stack all environmental layers and make sure they have appropriate names (currently manual and a bit messy!)
env_stack<-stack(env_list)
names(env_stack) <- env_names
names(env_stack)[10:17]<-c("waom4k_seafloorcurrents_absolute", "waom4k_seafloorcurrents_mean", 
                           "waom4k_seafloorcurrents_residual", "waom4k_seafloorsalinity", "waom4k_seafloortemperature",
                           "waom4k_test_flux08","waom4k_test_settle08","waom4k_test_susp08")



## only select rasters we actually need
pred_stack <- env_stack
plot(pred_stack)

# pred_stack_scaled1 <- subset(pred_stack, 1:8)
# for(i in 1:nlayers(pred_stack_scaled1)){
#   print(i)
#   k <- names(pred_stack_scaled1)[i]
#   c.sel <- which(names(cell_metadata_env_clean_scaled)==k)
#   s.sel <- which(names(pred_stack)==k)
#   pred_stack_scaled1[[i]] <- raster(pred_stack$depth)
#   pred_stack_scaled1[[i]] <- (pred_stack[[s.sel]]-scale.means[c.sel])/scale.sd[c.sel]
# }
# writeRaster(pred_stack_scaled1[[1]], filename=paste0(env.derived,"Circumpolar_EnvData_500m_shelf_scaled_bathy_gebco_depth.Rdata"), overwrite=TRUE)
# writeRaster(pred_stack_scaled1[[2]], filename=paste0(env.derived,"Circumpolar_EnvData_500m_shelf_scaled_bathy_gebco_depth2.Rdata"), overwrite=TRUE)
# writeRaster(pred_stack_scaled1[[3]], filename=paste0(env.derived,"Circumpolar_EnvData_500m_shelf_scaled_bathy_gebco_logslope.Rdata"), overwrite=TRUE)
# writeRaster(pred_stack_scaled1[[4]], filename=paste0(env.derived,"Circumpolar_EnvData_500m_shelf_scaled_bathy_gebco_slope.Rdata"), overwrite=TRUE)
# writeRaster(pred_stack_scaled1[[5]], filename=paste0(env.derived,"Circumpolar_EnvData_500m_shelf_scaled_bathy_gebco_tpi.Rdata"), overwrite=TRUE)
# writeRaster(pred_stack_scaled1[[6]], filename=paste0(env.derived,"Circumpolar_EnvData_500m_shelf_scaled_bathy_gebco_tpi11.Rdata"), overwrite=TRUE)
# writeRaster(pred_stack_scaled1[[7]], filename=paste0(env.derived,"Circumpolar_EnvData_500m_shelf_scaled_bathy_gebco_tpi5.Rdata"), overwrite=TRUE)
# writeRaster(pred_stack_scaled1[[8]], filename=paste0(env.derived,"Circumpolar_EnvData_500m_shelf_scaled_bathy_gebco_distance2canyons.Rdata"), overwrite=TRUE)

# pred_stack_scaled2 <- subset(pred_stack, 9:nlayers(pred_stack))
# for(i in 1:nlayers(pred_stack_scaled2)){
#   print(i)
#   k <- names(pred_stack_scaled2)[i]
#   c.sel <- which(names(cell_metadata_env_clean_scaled)==k)
#   s.sel <- which(names(pred_stack)==k)
#   pred_stack_scaled2[[i]] <- raster(pred_stack$depth)
#   pred_stack_scaled2[[i]] <- (pred_stack[[s.sel]]-scale.means[c.sel])/scale.sd[c.sel]
# }
# writeRaster(pred_stack_scaled2[[1]], filename=paste0(env.derived,"Circumpolar_EnvData_500m_shelf_scaled_bathy_gebco_distance2canyons2.Rdata"), overwrite=TRUE)
# writeRaster(pred_stack_scaled2[[2]], filename=paste0(env.derived,"Circumpolar_EnvData_500m_shelf_scaled_waom4k_seafloorcurrents_absolute.Rdata"), overwrite=TRUE)
# writeRaster(pred_stack_scaled2[[3]], filename=paste0(env.derived,"Circumpolar_EnvData_500m_shelf_scaled_waom4k_seafloorcurrents_mean.Rdata"), overwrite=TRUE)
# writeRaster(pred_stack_scaled2[[4]], filename=paste0(env.derived,"Circumpolar_EnvData_500m_shelf_scaled_waom4k_seafloorcurrents_residual.Rdata"), overwrite=TRUE)
# writeRaster(pred_stack_scaled2[[5]], filename=paste0(env.derived,"Circumpolar_EnvData_500m_shelf_scaled_waom4k_seafloorsalinity.Rdata"), overwrite=TRUE)
# writeRaster(pred_stack_scaled2[[6]], filename=paste0(env.derived,"Circumpolar_EnvData_500m_shelf_scaled_waom4k_seafloortemperature.Rdata"), overwrite=TRUE)
# writeRaster(pred_stack_scaled2[[7]], filename=paste0(env.derived,"Circumpolar_EnvData_500m_shelf_scaled_waom4k_seafloor_test_flux08.Rdata"), overwrite=TRUE)
# writeRaster(pred_stack_scaled2[[8]], filename=paste0(env.derived,"Circumpolar_EnvData_500m_shelf_scaled_waom4k_seafloor_test_settle08.Rdata"), overwrite=TRUE)
# writeRaster(pred_stack_scaled2[[9]], filename=paste0(env.derived,"Circumpolar_EnvData_500m_shelf_scaled_waom4k_seafloor_test_susp08.Rdata"), overwrite=TRUE)

## get file names of all environmental rasters and bricks and load into one big stack----
#all files with "gri" extension
env_list<-list.files(path = env.derived, pattern="gri$",  full.names=TRUE) 
#subset to  "shelf" files
env_list<-env_list[grep(".500m_shelf_scaled", env_list)]
#for the single rasters layer names are missing. Extract from file name.
env_names<-gsub(".*_|\\..*","",env_list)
#stack all environmental layers and make sure they have appropriate names (currently manual and a bit messy!)
env_stack_scaled<-stack(env_list)
names(env_stack_scaled) <- env_names
names(env_stack_scaled)[10:17]<-c("waom4k_seafloorcurrents_absolute", "waom4k_seafloorcurrents_mean", 
                                  "waom4k_seafloorcurrents_residual", "waom4k_seafloorsalinity", "waom4k_seafloortemperature",
                                  "waom4k_test_flux08","waom4k_test_settle08","waom4k_test_susp08")


##############################################
## 1-2hrs:
sel <- which(!is.na(env_stack_scaled$depth[]))

pred_stack.dat1 <- cbind(env_stack_scaled$depth[sel], env_stack_scaled$depth2[sel],
                         env_stack_scaled$distance2canyons[sel], env_stack_scaled$distance2canyons2[sel])
pred_stack.dat2 <- cbind(env_stack_scaled$logslope[sel], env_stack_scaled$slope[sel],
                         env_stack_scaled$tpi[sel], env_stack_scaled$tpi11[sel])
pred_stack.dat3 <- cbind(env_stack_scaled$tpi5[sel], env_stack_scaled$waom4k_test_flux08[sel],
                         env_stack_scaled$waom4k_test_settle08[sel], env_stack_scaled$waom4k_test_susp08[sel])
pred_stack.dat4 <- cbind(env_stack_scaled$waom4k_seafloorcurrents_absolute[sel], env_stack_scaled$waom4k_seafloorcurrents_mean[sel],
                         env_stack_scaled$waom4k_seafloorcurrents_residual[sel])
pred_stack.dat5 <- cbind(env_stack_scaled$waom4k_seafloorsalinity[sel], env_stack_scaled$waom4k_seafloortemperature[sel])

pred_stack.dat <- data.frame(cbind(pred_stack.dat1, pred_stack.dat2, pred_stack.dat3, pred_stack.dat4, pred_stack.dat5, 10))
names(pred_stack.dat) <- c("depth", "depth2", "distance2canyons","distance2canyons2",
                           "logslope","slope","tpi","tpi11",
                           "tpi5", "waom4k_test_flux08","waom4k_test_settle08","waom4k_test_susp08",
                           "waom4k_seafloorcurrents_absolute","waom4k_seafloorcurrents_mean","waom4k_seafloorcurrents_residual",
                           "waom4k_seafloorsalinity","waom4k_seafloortemperature",
                           "annotated_area")
pred_stack.dat$gear <- "OFOS"
pred_stack.dat$cover_cells_survey <- "PS96"
pred_stack.dat$cover_cells_transect1 <- "PS96_001"
save(pred_stack.dat, file=paste0(biodiv.dir,"biodiversity_pred_stack_scaled_dat.Rdata"))



###################################################################
## mask for environmental space

















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
#   c.sel <- which(names(cell_metadata_env_clean_scaled)==k)
#   pred_stack_scaled[[i]] <- raster(pred_stack$depth)
#   pred_stack_scaled[[i]] <- (pred_stack[[i]]-scale.means[c.sel])/scale.sd[c.sel]
# }
# 
# plot(pred_stack_scaled)
# 
