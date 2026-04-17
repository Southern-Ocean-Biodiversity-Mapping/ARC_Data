# analysis_helpers_3D_backtransform.R

suppressWarnings({
  if (!requireNamespace("terra", quietly = TRUE)) stop("Package 'terra' is required.")
})

#' Backtransform 3D particle endpoints to a bottom NPP raster.
#'
#' @param xyz_end Matrix/data.frame with columns x,y,z for all original particles.
#' @param npp_values Numeric vector of NPP values for each particle (same length as rows in xyz_end).
#' @param npp_template terra SpatRaster used as template for output.
#'
#' @return list with:
#' - xyz_end_npp: data.frame including npp and cellID
#' - cell_npp: data.frame with summed npp per cellID
#' - bottom_npp_raster: SpatRaster with summed npp per cell
backtransform_3D_to_bottom_raster <- function(xyz_end, npp_values, npp_template) {
  xyz_end <- as.data.frame(xyz_end)
  if (!all(c("x", "y") %in% names(xyz_end))) {
    stop("xyz_end must have columns named 'x' and 'y' (and optionally 'z').")
  }
  if (length(npp_values) != nrow(xyz_end)) {
    stop("npp_values must match nrow(xyz_end).")
  }
  
  xyz_end$npp <- npp_values
  
  # Raster cell IDs for each final particle position
  cell_id <- terra::extract(npp_template, xyz_end[, c("x", "y")], cells = TRUE)$cell
  xyz_end$cellID <- cell_id
  
  ok <- !is.na(xyz_end$cellID) & !is.na(xyz_end$npp)
  xyz_ok <- xyz_end[ok, ]
  
  # Sum NPP per raster cell
  cell_npp <- aggregate(xyz_ok$npp, by = list(cellID = xyz_ok$cellID), FUN = sum)
  names(cell_npp)[2] <- "npp_sum"
  
  out <- terra::rast(npp_template)
  out[] <- NA_real_
  out[cell_npp$cellID] <- cell_npp$npp_sum
  
  list(
    xyz_end_npp = xyz_end,
    cell_npp = cell_npp,
    bottom_npp_raster = out
  )
}