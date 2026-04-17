# analysis_helpers_compile.R

suppressWarnings({
  if (!requireNamespace("terra", quietly = TRUE)) stop("Package 'terra' is required.")
})

#' Compile regional rasters into a circumpolar raster.
#'
#' Fast path:
#' - if all regional rasters share the same geometry as the template, use mosaic with fun="max"
#'   (or another function), which is much faster than repeated extract().
#'
#' Fallback:
#' - extract regional values onto template coordinates and fill output with a max-like rule.
#'   This is slower but works when grids differ.
#'
#' @param template terra SpatRaster used as the circumpolar grid.
#' @param region_rasters list of SpatRaster objects (some may be NULL).
#' @param fun function name for mosaic aggregation (e.g., "max", "sum").
#' @param use_fast_if_aligned if TRUE, attempt fast mosaic when geometries match.
#' @param fallback_extract if TRUE, use extract-based fallback when not aligned.
#' @param verbose print progress
#'
#' @return SpatRaster on template geometry.
compile_regions_to_template <- function(template,
                                        region_rasters,
                                        fun = "max",
                                        use_fast_if_aligned = TRUE,
                                        fallback_extract = TRUE,
                                        verbose = TRUE) {
  stopifnot(inherits(template, "SpatRaster"))
  region_rasters <- Filter(Negate(is.null), region_rasters)
  if (length(region_rasters) == 0L) {
    out <- terra::rast(template)
    out[] <- NA_real_
    return(out)
  }
  
  # Check alignment (geometry)
  aligned <- TRUE
  if (use_fast_if_aligned) {
    for (r in region_rasters) {
      if (!terra::compareGeom(template, r, stopOnError = FALSE)) {
        aligned <- FALSE
        break
      }
    }
  } else {
    aligned <- FALSE
  }
  
  if (aligned) {
    if (verbose) message("Compiling using fast mosaic (aligned grids).")
    # Use terra::mosaic; fun controls overlap rule.
    out <- do.call(terra::mosaic, c(region_rasters, list(fun = fun)))
    # Ensure template geometry (mosaic should already match, but enforce if needed)
    out <- terra::resample(out, template, method = "near")
    return(out)
  }
  
  if (!fallback_extract) stop("Regional rasters are not aligned and fallback_extract=FALSE")
  
  if (verbose) message("Compiling using extract-based fallback (unaligned grids).")
  
  coords <- terra::crds(template)
  out_vals <- terra::values(template, mat = FALSE)
  out_vals[] <- NA_real_
  
  # For extract-based update, apply a max-like rule for overlaps.
  for (i in seq_along(region_rasters)) {
    if (verbose) message("  Region ", i, "/", length(region_rasters))
    v_add <- terra::extract(region_rasters[[i]], coords)[, 1]
    ok <- !is.na(v_add)
    if (!any(ok)) next
    
    # If out is NA, replace; else apply max (or sum if requested)
    if (fun == "max") {
      replace <- ok & (is.na(out_vals) | v_add > out_vals)
      out_vals[replace] <- v_add[replace]
    } else if (fun == "sum") {
      out_vals[ok] <- ifelse(is.na(out_vals[ok]), v_add[ok], out_vals[ok] + v_add[ok])
    } else {
      stop("Only fun='max' or fun='sum' supported in fallback mode.")
    }
  }
  
  out <- terra::rast(template)
  out[] <- out_vals
  out
}