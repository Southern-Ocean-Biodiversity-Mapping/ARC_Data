# analysis_helpers_flux_sed.R

# This script provides memory-safe post-processing utilities:
# - sedimentation raster from particle end positions (pend)
# - flux raster from streaming counts (kNN index space)
# - legacy flux raster from stored indices (fallback, chunked)
#
# Conventions:
# - kNN indices refer to ROMS horizontal grid point indices (linear indexing into Rdat6h$x / Rdat6h$y).
# - Raster cell IDs refer to terra raster cell numbers (1..ncell(r)).
# - The mapping kNN -> raster cell is computed using the kNN cell coordinate (x,y) and cellFromXY().

suppressWarnings({
  if (!requireNamespace("terra", quietly = TRUE)) stop("Package 'terra' is required.")
})

#' Build a mapping from ROMS kNN cell indices to raster cell IDs.
#'
#' This computes raster cell IDs for every ROMS horizontal grid point, using its x/y coordinates.
#' The result is a vector where map_knn_to_cell[knn_id] = raster_cell_id (or NA if outside raster).
#'
#' @param npp_raster terra SpatRaster defining target grid/projection.
#' @param roms_x 2D matrix of ROMS x coordinates (same shape as roms_y).
#' @param roms_y 2D matrix of ROMS y coordinates.
#' @param chunk_size How many kNN points to map per chunk (controls peak RAM).
#'
#' @return Integer vector of length prod(dim(roms_x)).
build_knn_to_raster_cell_map <- function(npp_raster, roms_x, roms_y, chunk_size = 200000L) {
  stopifnot(inherits(npp_raster, "SpatRaster"))
  n_knn <- length(roms_x)
  map_knn_to_cell <- rep(NA_integer_, n_knn)
  
  # Flatten coordinates once (these are vectors, not copies of big arrays)
  x <- as.vector(roms_x)
  y <- as.vector(roms_y)
  
  idx <- seq_len(n_knn)
  chunk_size <- as.integer(chunk_size)
  
  for (i1 in seq(1L, n_knn, by = chunk_size)) {
    i2 <- min(i1 + chunk_size - 1L, n_knn)
    ii <- idx[i1:i2]
    xy <- cbind(x[ii], y[ii])
    
    # cellFromXY is faster and lighter than extract(..., cells=TRUE) for pure indexing
    map_knn_to_cell[ii] <- terra::cellFromXY(npp_raster, xy)
  }
  
  map_knn_to_cell
}

#' Create a sedimentation raster from particle end positions (pend).
#'
#' Counts how many particles settle in each raster cell and writes those counts into a raster.
#' Uses streaming-friendly tabulate aggregation.
#'
#' @param npp_raster terra SpatRaster template defining output grid.
#' @param pend Matrix/data.frame with columns x,y (and optionally z).
#' @param chunk_size Chunk size for coordinate-to-cell mapping.
#'
#' @return SpatRaster with per-cell sedimentation counts (integer).
sedimentation_raster_from_pend <- function(npp_raster, pend, chunk_size = 200000L) {
  stopifnot(inherits(npp_raster, "SpatRaster"))
  pend <- as.matrix(pend)
  if (ncol(pend) < 2) stop("pend must have at least 2 columns (x,y)")
  
  n <- nrow(pend)
  cell_id <- rep(NA_integer_, n)
  
  chunk_size <- as.integer(chunk_size)
  for (i1 in seq(1L, n, by = chunk_size)) {
    i2 <- min(i1 + chunk_size - 1L, n)
    ii <- i1:i2
    cell_id[ii] <- terra::cellFromXY(npp_raster, pend[ii, 1:2, drop = FALSE])
  }
  
  # Drop NA cells (outside raster)
  ok <- !is.na(cell_id)
  counts <- tabulate(cell_id[ok], nbins = terra::ncell(npp_raster))
  
  out <- terra::rast(npp_raster)
  out[] <- NA_real_
  out[counts > 0] <- counts[counts > 0]
  out
}

#' Create a flux raster from streaming flux counts in kNN-index space.
#'
#' This matches your previous approach:
#'  1) flux_counts is indexed by ROMS kNN cell IDs
#'  2) translate kNN IDs -> raster cell IDs (precomputed map)
#'  3) aggregate to raster cells and write the raster
#'
#' @param npp_raster terra SpatRaster template defining output grid.
#' @param flux_counts Integer/numeric vector where flux_counts[knn_id] is the visit count.
#' @param map_knn_to_cell Integer vector mapping knn_id -> raster cell ID (from build_knn_to_raster_cell_map()).
#'
#' @return SpatRaster with per-cell flux values.
flux_raster_from_knn_counts <- function(npp_raster, flux_counts, map_knn_to_cell) {
  stopifnot(inherits(npp_raster, "SpatRaster"))
  
  if (is.null(flux_counts)) stop("flux_counts is NULL. Run tracking with flux_mode='entry'/'total'/'presence' and return flux_counts.")
  if (length(map_knn_to_cell) < length(flux_counts)) {
    stop("map_knn_to_cell must be at least as long as flux_counts (should be length prod(dim(roms_x))).")
  }
  
  # Which kNN cells have nonzero flux?
  knn_idx <- which(flux_counts > 0)
  if (length(knn_idx) == 0L) {
    out <- terra::rast(npp_raster)
    out[] <- NA_real_
    return(out)
  }
  
  # Translate kNN -> raster cell
  cell_id <- map_knn_to_cell[knn_idx]
  ok <- !is.na(cell_id)
  if (!any(ok)) {
    out <- terra::rast(npp_raster)
    out[] <- NA_real_
    return(out)
  }
  
  cell_id <- cell_id[ok]
  vals <- flux_counts[knn_idx[ok]]
  
  # Aggregate to raster cells (some kNN points can map to same raster cell)
  # Use rowsum for speed and low memory.
  agg <- rowsum(vals, group = cell_id, reorder = FALSE)
  cell_unique <- as.integer(rownames(agg))
  val_unique <- as.numeric(agg[, 1])
  
  out <- terra::rast(npp_raster)
  out[] <- NA_real_
  out[cell_unique] <- val_unique
  out
}

#' Legacy fallback: derive flux_counts from stored idx_list_2D (chunked, no tibble/dplyr).
#'
#' idx_list_2D is expected to be a list over wrapper loops, each containing a list of timesteps.
#' Each timestep element is an integer vector of kNN IDs for active particles at that step.
#'
#' This produces total or entry flux (entry cannot be reconstructed exactly without previous-cell state
#' unless the stored indices are per-particle and aligned; here we implement "total" visits).
#'
#' @param idx_list_2D Nested list structure from older runs (track.2D$idx_list_2D).
#' @param n_knn Total number of kNN cells (typically prod(dim(roms_x))).
#' @param verbose Print progress.
#'
#' @return Integer vector flux_counts (length n_knn).
legacy_flux_counts_from_idx_list <- function(idx_list_2D, n_knn, verbose = TRUE) {
  n_knn <- as.integer(n_knn)
  flux_counts <- integer(n_knn)
  
  for (irun in seq_along(idx_list_2D)) {
    if (verbose) message("Legacy flux: processing wrapper loop ", irun, "/", length(idx_list_2D))
    steps <- idx_list_2D[[irun]]
    if (is.null(steps) || length(steps) == 0L) next
    
    for (itime in seq_along(steps)) {
      idx <- steps[[itime]]
      if (is.null(idx) || length(idx) == 0L) next
      # Update counts for this timestep
      flux_counts <- flux_counts + tabulate(idx, nbins = n_knn)
    }
  }
  
  flux_counts
}