# analysis_helpers_io_cache.R
# Cohesive I/O + caching utilities for the ptrackR workflow.

suppressWarnings({
  if (!requireNamespace("terra", quietly = TRUE)) stop("Package 'terra' is required.")
})

#' Create directories if needed.
ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

#' Standardised filename builder for tracking outputs.
#'
#' @param out_dir Output directory.
#' @param stage "2D" or "3D"
#' @param region Integer region index.
#' @param speed Sinking speed (m/day)
#' @param days_label Character label for days in filename (e.g., "21days")
#' @param npp_label Character label e.g. "NPP9"
#' @param extra Additional string inserted before extension (e.g. "traj_", "r0001_", "red").
#' @param ext Extension including dot, e.g. ".Rdata" or ".tif"
#'
#' @return Full filepath.
fname_tracking <- function(out_dir,
                           stage = c("2D", "3D"),
                           region,
                           speed,
                           days_label = "21days",
                           npp_label = "NPP9",
                           extra = "",
                           ext = ".Rdata") {
  stage <- match.arg(stage)
  ensure_dir(out_dir)
  paste0(out_dir,
         "/tracking", stage,
         "_Region", sprintf("%02d", region),
         "_", npp_label, "_", speed, "mday_", days_label,
         if (nzchar(extra)) paste0("_", extra) else "",
         ext)
}

#' Standardised filename builder for regional ROMS object.
fname_roms_region <- function(env_dir3, region) {
  paste0(env_dir3,
         "ocean_his_TrackingSetup_6hourlycurrents_28days_Region",
         sprintf("%02d", region), ".Rdata")
}

#' Standardised filename builder for NPP raster template per region.
fname_npp_raster_region <- function(env_dir3, region, base = "ocean_his_TrackingSetup_NPP_Region") {
  paste0(env_dir3, base, sprintf("%02d", region), ".tif")
}

#' Standardised filename builder for bottom NPP raster produced by 3D backtransform.
fname_bottom_npp_region <- function(out_dir, region, speed, days_label = "21days", npp_label = "NPP9") {
  paste0(out_dir,
         "/tracking3D_Region", sprintf("%02d", region),
         "_", npp_label, "_", speed, "mday_", days_label,
         "_bottomnpp.tif")
}

#' Cache filename for kNN->raster-cell mapping.
fname_knn_map_cache <- function(cache_dir, region, speed, tag = "nppcrop") {
  ensure_dir(cache_dir)
  paste0(cache_dir, "/knn_to_cell_map_Region", sprintf("%02d", region),
         "_", tag, "_", speed, "mday.rds")
}

#' Load ROMS object (expects object name Rdat6h inside .Rdata).
load_roms_region <- function(env_dir3, region) {
  f <- fname_roms_region(env_dir3, region)
  e <- new.env(parent = emptyenv())
  load(f, envir = e)
  if (!exists("Rdat6h", envir = e)) stop("Expected object 'Rdat6h' in: ", f)
  get("Rdat6h", envir = e)
}

#' Load raster template for region.
load_npp_raster_template <- function(env_dir3, region) {
  terra::rast(fname_npp_raster_region(env_dir3, region))
}

#' Load or build kNN->raster-cell map, with caching.
#'
#' Requires build_knn_to_raster_cell_map() from analysis_helpers_flux_sed.R
#'
#' @param cache_dir Directory to store map RDS files.
#' @param region Region index.
#' @param speed Speed used only for cache naming (helps avoid collisions across experiments).
#' @param npp_raster terra SpatRaster template.
#' @param roms_x,roms_y ROMS x/y matrices.
#' @param chunk_size Chunk size for mapping.
#'
#' @return Integer vector: map_knn_to_cell[knn_id] = raster cell ID (or NA).
load_or_build_knn_map <- function(cache_dir,
                                  region,
                                  speed,
                                  npp_raster,
                                  roms_x,
                                  roms_y,
                                  chunk_size = 200000L) {
  cache_file <- fname_knn_map_cache(cache_dir, region, speed)
  if (file.exists(cache_file)) {
    return(readRDS(cache_file))
  }
  map_knn_to_cell <- build_knn_to_raster_cell_map(
    npp_raster = npp_raster,
    roms_x = roms_x,
    roms_y = roms_y,
    chunk_size = chunk_size
  )
  saveRDS(map_knn_to_cell, cache_file)
  map_knn_to_cell
}