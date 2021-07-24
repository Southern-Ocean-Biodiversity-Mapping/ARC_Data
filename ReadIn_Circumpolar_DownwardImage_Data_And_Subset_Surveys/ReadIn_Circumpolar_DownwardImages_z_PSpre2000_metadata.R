

### Survey PS06:
image_dir <- "C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/Stills/PS06/"
dat.raw <- read.table(paste0(image_dir,"PS06_metadata.txt"), fill=TRUE, skip=6)
n <- length(levels(dat.raw[,1]))
dat.full <- data.frame(matrix(ncol=n,nrow=10))
names(dat.full) <- levels(dat.raw[,1])
for(i in 1:n){
  dat.full[,i] <- factor(dat.raw[dat.raw[,1]==names(dat.full)[i],2])
}
names(dat.full)[1] <- "TransectID" 
names(dat.full)[5] <- "Median_Latitude" 
names(dat.full)[7] <- "Median_Longitude" 
dat.full$SurveyID <- as.factor("PS06")
dat.full[,1] <- as.factor(gsub("PS06/","",as.character(dat.full[,1])))
dat.full[,4] <- as.numeric(gsub("_m","",as.character(dat.full[,4])))
dat.full[,8] <- as.numeric(gsub("_m","",as.character(dat.full[,8])))
dat.full[,10] <- as.numeric(gsub("_m","",as.character(dat.full[,10])))
for(i in c(5,7)){
  dat.full[,i] <- as.numeric(as.character(dat.full[,i]))
}
PS06_metadata <- dat.full
save(PS06_metadata, file=paste0(image_dir,"PS06_metadata.Rdata"))


### Survey PS14:
image_dir <- "C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/Stills/PS14/"
dat.raw <- read.table(paste0(image_dir,"PS14_metadata.txt"), fill=TRUE, skip=7)
n <- length(levels(dat.raw[,1]))
dat.full <- data.frame(matrix(ncol=n,nrow=24))
names(dat.full) <- levels(dat.raw[,1])
for(i in 1:n){
  dat.full[,i] <- factor(dat.raw[dat.raw[,1]==names(dat.full)[i],2])
}
names(dat.full)[1] <- "TransectID" 
dat.full$SurveyID <- as.factor("PS14")
dat.full[,1] <- as.factor(gsub("PS14/","",as.character(dat.full[,1])))
dat.full[,5] <- as.numeric(gsub("_m","",as.character(dat.full[,5])))
dat.full[,6] <- as.numeric(gsub("_m","",as.character(dat.full[,6])))
dat.full[,12] <- as.numeric(gsub("_m","",as.character(dat.full[,12])))
dat.full[,16] <- as.numeric(gsub("_m","",as.character(dat.full[,16])))
for(i in c(4:8,10:14,16:19)){
  dat.full[,i] <- as.numeric(as.character(dat.full[,i]))
}
PS14_metadata <- dat.full
save(PS14_metadata, file=paste0(image_dir,"PS14_metadata.Rdata"))


### Survey PS18:
image_dir <- "C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/Stills/PS18/"
dat.raw <- read.table(paste0(image_dir,"PS18_metadata.txt"), fill=TRUE, skip=7)
n <- length(levels(dat.raw[,1]))
dat.full <- data.frame(matrix(ncol=n, nrow=25))
names(dat.full) <- levels(dat.raw[,1])
for(i in 1:n){
  dat.full[,i] <- factor(dat.raw[dat.raw[,1]==names(dat.full)[i],2])
}
names(dat.full)[1] <- "TransectID" 
dat.full$SurveyID <- as.factor("PS18")
dat.full[,1] <- as.factor(gsub("PS18/","",as.character(dat.full[,1])))
dat.full[,5] <- as.numeric(gsub("_m","",as.character(dat.full[,5])))
dat.full[,6] <- as.numeric(gsub("_m","",as.character(dat.full[,6])))
dat.full[,12] <- as.numeric(gsub("_m","",as.character(dat.full[,12])))
dat.full[,16] <- as.numeric(gsub("_m","",as.character(dat.full[,16])))
for(i in c(4:8,10:14,16:19)){
  dat.full[,i] <- as.numeric(as.character(dat.full[,i]))
}
PS18_metadata <- dat.full
save(PS18_metadata, file=paste0(image_dir,"PS18_metadata.Rdata"))


### Survey PS61:
image_dir <- "C:/Users/jjansen/OneDrive - University of Tasmania/Desktop/science/data_biological/Stills/PS61/"
dat.raw <- read.table(paste0(image_dir,"PS61_metadata.txt"), fill=TRUE, skip=6)
dat.full.t <- data.frame(t(dat.raw))
names(dat.full.t) <- dat.raw[,1]
dat.full <- dat.full.t[-1,]
rownames(dat.full) <- ""
names(dat.full)[1] <- "TransectID" 
dat.full[,1] <- as.factor(gsub("PS61/","",as.character(dat.full[,1])))
dat.full[,8] <- as.numeric(gsub("_m","",as.character(dat.full[,8])))
dat.full[,9] <- as.numeric(gsub("_m","",as.character(dat.full[,9])))
dat.full[,16] <- as.numeric(gsub("_m","",as.character(dat.full[,16])))
dat.full[,17] <- as.numeric(gsub("_m","",as.character(dat.full[,17])))
for(i in 2:13){
  dat.full[,i] <- as.numeric(as.character(dat.full[,i]))
}
for(i in c(14,15,18,19)){
  dat.full[,i] <- factor(dat.full[,i])
}
dat.full$SurveyID <- as.factor("PS61")
PS61_metadata <- dat.full
save(PS61_metadata, file=paste0(image_dir,"PS61_metadata.Rdata"))
