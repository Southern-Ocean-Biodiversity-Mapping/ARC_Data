# setup_knn.R (updated)
# Setup kNN lookup. Can build 2D-only or both 2D and 3D trees.

setup_knn <- function(x_roms = x, y_roms = y, depth_roms = hh, build_3d = TRUE) {
  require(nabor)
  
  # 2D tree uses horizontal coordinates only
  xy <- cbind(as.vector(x_roms), as.vector(y_roms))
  colnames(xy) <- c("x","y")
  kdxy <- WKNND(xy)
  
  if (!build_3d) {
    return(list(kdtree = NULL, kdxy = kdxy))
  }
  
  # 3D tree uses full (x,y,z) set
  allxyz <- cbind(as.vector(x_roms), as.vector(y_roms), as.vector(depth_roms))
  colnames(allxyz) <- c("x","y","z")
  kdtree <- WKNND(allxyz)
  
  list(kdtree = kdtree, kdxy = kdxy)
}