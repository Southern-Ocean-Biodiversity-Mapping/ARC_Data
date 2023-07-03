## Adding thumbnails to the image folders that we prepare for squidle


library(magick)
##
r.dir <- "R:/IMAS/Antarctic_Seafloor/SQUIDLE_dataset_20230703/"
## create a vector with all image names and locations
f1 <- c(list.files(r.dir, recursive=TRUE, pattern=".jpg"),
        list.files(r.dir, recursive=TRUE, pattern=".png"),
        list.files(r.dir, recursive=TRUE, pattern=".JPG"))
target.dirs <- gsub("images", "thumbnails", f1)

for(i in 1:length(f1)){
  print(i)
  img <- image_read(paste0(r.dir,f1[i]))
  img_thumb <- image_scale(img, "500x300")
  image_write(img_thumb, paste0(r.dir,target.dirs[i]))
}

