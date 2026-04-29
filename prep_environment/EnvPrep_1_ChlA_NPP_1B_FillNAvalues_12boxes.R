####################################################################
## Fill NPP values based on chl-npp relationships in 12 lon-lat regions around Antarctica
## using monthly summer average files
####################################################################

## specify user and setup directory to look up data from
usr <- "VM"
source("0_SourceFile.R")

## set input and output folders
env.dir <- paste0(usr.main.dir,"data_environmental/derived/NPP")

##########################
library(terra)

monthly.cafe <- rast(file.path(env.dir,"NPP_monthly_OctMar_cafe.tif"))
monthly.cbpm <- rast(file.path(env.dir,"NPP_monthly_OctMar_cbpm.tif"))
monthly.eppl <- rast(file.path(env.dir,"NPP_monthly_OctMar_eppley.tif"))
monthly.vpmg <- rast(file.path(env.dir,"NPP_monthly_OctMar_vpmg.tif"))
chl.s <-        rast(file.path(env.dir,"NPP_monthly_OctMar_chla.tif"))

# ---- Standardize missing values ----
# Interpret 0 and negatives as missing for both NPP and chlorophyll.
monthly.cafe[monthly.cafe <= 0] <- NA
monthly.cbpm[monthly.cbpm <= 0] <- NA
monthly.eppl[monthly.eppl <= 0] <- NA
monthly.vpmg[monthly.vpmg <= 0] <- NA
chl.s[chl.s <= 0] <- NA

# ---- NPP infill function (per stack) ----
# Fills NA NPP using mean NPP within geographic boxes and similar chl bins.
# For each layer (e.g., month), compute mean NPP by (box_id, chl_bin) from observed cells,
# then impute NA cells by exact bin match or progressively widened bins (±1, ±2, ±5).
fill_npp_from_chl <- function(npp, chl,
                              bin_width = 0.01,
                              widen_bins = c(0L, 1L, 2L, 5L),
                              combine_box = list(`11` = 7),
                              verbose = TRUE) {
  stopifnot(inherits(npp, "SpatRaster"), inherits(chl, "SpatRaster"))
  stopifnot(nlyr(npp) >= 1L, nlyr(chl) >= 1L)
  
  # Ensure geometry is consistent (same extent, resolution, rows/cols).
  # CRS should already match; resample chl if grid differs.
  if (!compareGeom(npp, chl, stopOnError = FALSE)) {
    chl <- resample(chl, npp, method = "bilinear")
  }
  
  # Build fixed geographic boxes
  lon_r <- init(npp, "x")
  lat_r <- init(npp, "y")
  
  # (12 boxes) from longitudes/latitudes on the NPP grid.
  # Longitude codes: (-Inf,-90] -> 1; (-90,0] -> 2; (0,90] -> 3; (90,Inf] -> 4
  lon_rcl <- matrix(c(-Inf, -90, 1,
                      -90,   0, 2,
                      0,  90, 3,
                      90, Inf, 4), byrow = TRUE, ncol = 3)
  # Build fixed geographic boxes (12 boxes) from longitudes/latitudes on the NPP grid.
  # Longitude codes: (-Inf,-120] -> 1; (-120,60] -> 2;  (-60,0] -> 3;(0,60] -> 4; (60,120] -> 5; (120,Inf] -> 6
  # lon_rcl <- matrix(c(-Inf, -120, 1,
  #                     -120,  -60, 2,
  #                     - 60,    0, 3,
  #                        0,   60, 4,
  #                       60,  120, 5,
  #                      120,  Inf, 6), byrow = TRUE, ncol = 3)
  
  # Latitude  codes: (-Inf,-70] -> 3; (-70,-60] -> 2; (-60,Inf] -> 1
  lat_rcl <- matrix(c(-Inf, -70, 3,
                      -70,  -60, 2,
                      -60,  Inf, 1), byrow = TRUE, ncol = 3)
  
  lon_code <- classify(lon_r, lon_rcl, include.lowest = TRUE, right = TRUE)
  lat_code <- classify(lat_r, lat_rcl, include.lowest = TRUE, right = TRUE)
  
  # box_id = (lat_code - 1) * 4 + lon_code, so box_id ∈ {1..12}
  box_id_r <- app(c(lon_code, lat_code), fun = function(v) (v[2] - 1L) * 4L + v[1])
  box_id   <- values(box_id_r)
  
  # Output scaffold
  npp_out <- npp
  L_npp   <- nlyr(npp)
  L_chl   <- nlyr(chl)
  
  # Helper: aggregate mean NPP by (box_id, chl_bin) for one layer
  compute_mean_by_box_bin <- function(npp_vec, chl_bin, box_id) {
    ok <- !is.na(npp_vec) & !is.na(chl_bin) & !is.na(box_id)
    if (!any(ok)) {
      return(list(keys = integer(), means = numeric(), by_box = vector("list", 12L),
                  min_bin = 0L, width = 0L, stride = 10000000L))
    }
    stride <- 10000000L  # large enough to avoid key collisions
    keys <- box_id[ok] * stride + chl_bin[ok]
    
    df  <- data.frame(key = keys, npp = npp_vec[ok])
    agg <- aggregate(npp ~ key, data = df, FUN = function(x) mean(x, na.rm = TRUE))
    
    by_box <- vector("list", 12L)
    min_bin <- suppressWarnings(min(chl_bin, na.rm = TRUE)); if (!is.finite(min_bin)) min_bin <- 0L
    max_bin <- suppressWarnings(max(chl_bin, na.rm = TRUE)); if (!is.finite(max_bin)) max_bin <- min_bin
    width <- max(0L, max_bin - min_bin + 1L)
    for (b in 1:12) by_box[[b]] <- rep(NA_real_, width)
    
    k_box <- floor(agg$key / stride)
    k_bin <- agg$key - k_box * stride
    for (i in seq_len(nrow(agg))) {
      b  <- k_box[i]
      bb <- k_bin[i] - min_bin + 1L
      if (bb >= 1L && bb <= width) {
        by_box[[b]][bb] <- agg$npp[i]
      }
    }
    list(keys = agg$key, means = agg$npp, by_box = by_box,
         min_bin = min_bin, width = width, stride = stride)
  }
  
  # Helper: widened (box, bin) lookup for target cells
  widened_lookup <- function(target_box, target_bin, box_means, min_bin, width, widen_bins, combine_box) {
    idx_in_box <- which(!is.na(target_box) & !is.na(target_bin))
    res <- rep(NA_real_, length(target_box))
    if (!length(idx_in_box)) return(res)
    
    local_bin <- target_bin[idx_in_box] - min_bin + 1L
    ok_local <- local_bin >= 1L & local_bin <= width
    idx_in_box <- idx_in_box[ok_local]
    local_bin  <- local_bin[ok_local]
    
    for (w in widen_bins) {
      cand_idx1 <- pmax(1L, local_bin - w)
      cand_idx2 <- pmin(width, local_bin + w)
      
      bvec <- target_box[idx_in_box]
      # Extract means from per-box vectors
      r1 <- mapply(function(b, lb) box_means[[b]][lb], bvec, cand_idx1)
      r2 <- mapply(function(b, lb) box_means[[b]][lb], bvec, cand_idx2)
      pick <- ifelse(!is.na(r1), r1, r2)
      
      if (anyNA(pick) && length(combine_box)) {
        na_pos <- which(is.na(pick))
        for (p in na_pos) {
          b <- as.character(bvec[p])
          if (!is.null(combine_box[[b]])) {
            b2 <- combine_box[[b]]
            v1 <- box_means[[b2]][cand_idx1[p]]
            v2 <- box_means[[b2]][cand_idx2[p]]
            pick[p] <- ifelse(!is.na(v1), v1, v2)
          }
        }
      }
      
      to_fill <- which(is.na(res[idx_in_box]) & !is.na(pick))
      if (length(to_fill)) {
        res[idx_in_box[to_fill]] <- pick[to_fill]
      }
      if (all(!is.na(res[idx_in_box]))) break
    }
    res
  }
  
  # Layer-by-layer processing
  for (j in 1:L_npp) {
    if (verbose) message(sprintf("Layer %d/%d", j, L_npp))
    
    # One-to-one alignment: expect chl[[j]] corresponds to npp[[j]]
    chl_j <- chl[[j]]
    chl_vals <- values(chl_j)
    chl_bin  <- as.integer(floor(chl_vals / bin_width))
    
    npp_vec <- values(npp[[j]])
    
    fill_mask <- is.na(npp_vec) & !is.na(chl_bin) & !is.na(box_id)
    if (!any(fill_mask)) {
      npp_out[[j]] <- setValues(npp_out[[j]], npp_vec)
      next
    }
    
    agg_obj <- compute_mean_by_box_bin(npp_vec, chl_bin, box_id)
    
    keys_all <- box_id * agg_obj$stride + chl_bin
    filled0 <- rep(NA_real_, length(npp_vec))
    if (length(agg_obj$keys)) {
      m <- match(keys_all, agg_obj$keys)
      filled0 <- ifelse(!is.na(m), agg_obj$means[m], NA_real_)
    }
    
    target_idx <- which(fill_mask)
    widened <- widened_lookup(
      target_box = box_id[target_idx],
      target_bin = chl_bin[target_idx],
      box_means  = agg_obj$by_box,
      min_bin    = agg_obj$min_bin,
      width      = agg_obj$width,
      widen_bins = widen_bins,
      combine_box = combine_box
    )
    
    v_final <- filled0[target_idx]
    need_widen <- is.na(v_final)
    v_final[need_widen] <- widened[need_widen]
    
    npp_vec[target_idx] <- v_final
    npp_out[[j]] <- setValues(npp_out[[j]], npp_vec)
  }
  npp_out
}

# ---- Utility: align NPP stack to chlorophyll by layer names ----
# Restricts both stacks to the common set of layer names and orders chl to match NPP.
align_to_chl_names <- function(npp_stack, chl_stack) {
  npp_nm <- names(npp_stack)
  chl_nm <- names(chl_stack)
  common <- intersect(npp_nm, chl_nm)   # preserves NPP order while restricting to names in chl
  if (!length(common)) stop("No common layer names between NPP and chlorophyll.")
  npp_aligned <- npp_stack[[common]]
  chl_aligned <- chl_stack[[common]]
  list(npp = npp_aligned, chl = chl_aligned)
}

# ---- Per-model processing ----
# Each model is processed independently with layer-wise alignment to chlorophyll.
# Output uses LZW compression and 32-bit float unless double precision is required.

gdal_opts <- c("COMPRESS=LZW", "PREDICTOR=2")
dtype    <- "FLT4S"
na_flag  <- -9999

# 1) CAFE
aln <- align_to_chl_names(monthly.cafe, chl.s)
cafe.filled <- fill_npp_from_chl(
  npp = aln$npp,
  chl = aln$chl,
  bin_width  = 0.01,
  widen_bins = c(0L, 1L, 2L, 5L),
  combine_box = list(`11` = 7),
  verbose = TRUE
)
writeRaster(cafe.filled,
            filename = file.path(env.dir, "NPP_filled_monthly_OctMar_12boxes_cafe.tif"),
            overwrite = TRUE, gdal = gdal_opts, datatype = dtype, NAflag = na_flag)
writeRaster(aln$npp,
            filename = file.path(env.dir, "NPP_unfilled_monthly_OctMar_cafe.tif"),
            overwrite = TRUE, gdal = gdal_opts, datatype = dtype, NAflag = na_flag)

# 2) CBPM (has fewer layers; alignment trims to common set)
aln <- align_to_chl_names(monthly.cbpm, chl.s)
cbpm.filled <- fill_npp_from_chl(
  npp = aln$npp,
  chl = aln$chl,
  bin_width  = 0.01,
  widen_bins = c(0L, 1L, 2L, 5L),
  combine_box = list(`11` = 7),
  verbose = TRUE
)
writeRaster(cbpm.filled,
            filename = file.path(env.dir, "NPP_filled_monthly_OctMar_12boxes_cbpm.tif"),
            overwrite = TRUE, gdal = gdal_opts, datatype = dtype, NAflag = na_flag)
writeRaster(aln$npp,
            filename = file.path(env.dir, "NPP_unfilled_monthly_OctMar_12boxes_cbpm.tif"),
            overwrite = TRUE, gdal = gdal_opts, datatype = dtype, NAflag = na_flag)

# 3) EPPLEY
aln <- align_to_chl_names(monthly.eppl, chl.s)
eppl.filled <- fill_npp_from_chl(
  npp = aln$npp,
  chl = aln$chl,
  bin_width  = 0.01,
  widen_bins = c(0L, 1L, 2L, 5L),
  combine_box = list(`11` = 7),
  verbose = TRUE
)
writeRaster(eppl.filled,
            filename = file.path(env.dir, "NPP_filled_monthly_OctMar_12boxes_eppley.tif"),
            overwrite = TRUE, gdal = gdal_opts, datatype = dtype, NAflag = na_flag)
writeRaster(aln$npp,
            filename = file.path(env.dir, "NPP_unfilled_monthly_OctMar_eppley.tif"),
            overwrite = TRUE, gdal = gdal_opts, datatype = dtype, NAflag = na_flag)

# 4) VPMG
aln <- align_to_chl_names(monthly.vpmg, chl.s)
vpmg.filled <- fill_npp_from_chl(
  npp = aln$npp,
  chl = aln$chl,
  bin_width  = 0.01,
  widen_bins = c(0L, 1L, 2L, 5L),
  combine_box = list(`11` = 7),
  verbose = TRUE
)
writeRaster(vpmg.filled,
            filename = file.path(env.dir, "NPP_filled_monthly_OctMar_12boxes_vpmg.tif"),
            overwrite = TRUE, gdal = gdal_opts, datatype = dtype, NAflag = na_flag)
writeRaster(aln$npp,
            filename = file.path(env.dir, "NPP_unfilled_monthly_OctMar_vpmg.tif"),
            overwrite = TRUE, gdal = gdal_opts, datatype = dtype, NAflag = na_flag)

# ---- Optional: re-concatenate all filled models into a single stack ----
# npp.filled.all <- c(cafe.filled, cbpm.filled, eppl.filled, vpmg.filled)
# writeRaster(npp.filled.all,
#             filename = file.path(env.dir, "NPP_filled_monthly_12boxes_SO_ALL.tif"),
#             overwrite = TRUE, gdal = gdal_opts, datatype = dtype, NAflag = na_flag)


##########
## plot relationships between the filled NPP values and the underlying chl-a values
























































































npp.s <- c(monthly.cafe, monthly.cbpm, monthly.eppl, monthly.vpmg)
## Clean negatives once (vectorized)
npp.s[npp.s < 0] <- NA
chl.s[chl.s < 0] <- NA

#### crop raster to SO-extent
so_ext <- ext(-180, 180, -90, -50)
npp.crop <- crop(npp.s, so_ext)
chl.crop   <- crop(chl.s, so_ext)

fill_npp_from_chl <- function(npp, chl, lon_breaks = c(-Inf, -90, 0, 90, Inf),
                              lat_breaks = c(-Inf, -70, -60, Inf),
                              bin_width = 0.01,
                              widen_bins = c(0L, 1L, 2L, 5L),
                              combine_box = list(`11` = 7),
                              n_models = 4L,
                              layer_order = c("models_fast", "months_fast"),
                              verbose = TRUE) {
  # ---- Function goal ----
  # Fill NA NPP values using the relationship between NPP and chlorophyll-a (chl) in
  # nearby geographic boxes and similar chl values. For each (box, chl_bin) we compute
  # the mean NPP from observed data, then impute NA NPPs by matching their (box, chl_bin).
  # If no exact bin match exists, we widen the bin tolerance (±1, ±2, ±5 bins).
  #
  # Arguments:
  #   npp        : SpatRaster (≥1 layer) with NPP values (may contain NAs)
  #   chl        : SpatRaster (1+ layers) with chlorophyll-a values (same CRS as npp).
  #                Supports climatology (1 layer) or multi-layer (e.g., monthly).
  #   lon_breaks : (unused by current code) intended custom lon breaks for boxes
  #   lat_breaks : (unused by current code) intended custom lat breaks for boxes
  #   bin_width  : chl bin size; e.g. 0.01 means bins of width 0.01 chl units
  #   widen_bins : integer bin radii to widen search (0 => exact, 1 => ±0.01, 2 => ±0.02, etc.)
  #   combine_box: list mapping "sparse" boxes to a fallback box (e.g., box 11 -> 7)
  #   n_models   : number of NPP models concatenated in `npp` (used if chl has multiple layers)
  #   layer_order: how `npp` is stacked relative to months/models; choose "models_fast" (default)
  #                if your order is [model1_m1..m12, model2_m1..m12, ...]; choose "months_fast"
  #                if your order is [m1_model1, m1_model2, m1_model3, m1_model4, m2_model1, ...]
  #   verbose    : print progress messages
  #
  # Notes:
  # - This function vectorizes everything (no per-cell scanning) for speed.
  # - Memory: values() pulls full rasters into memory; OK for ~2M cells, but keep in mind.
  
  layer_order <- match.arg(layer_order)
  stopifnot(nlyr(npp) >= 1, nlyr(chl) >= 1)
  
  # ---- Geometry alignment ----
  # Ensure chl matches the npp grid (same extent/resolution/nrows/ncols). We assume CRS matches.
  # If geometry differs, resample chl to npp grid (bilinear is reasonable for continuous chl).
  # If CRS differs, project() should be used prior to this function (not handled here).
  if (!compareGeom(npp, chl, stopOnError = FALSE)) {
    chl <- resample(chl, npp, method = "bilinear")
  }
  
  # ---- Build geographic box IDs (1..12) from lon/lat ----
  # We create a box raster where each cell is assigned to one of 12 boxes:
  #   lon_code ∈ {1,2,3,4} (<= -90, (-90,0], (0,90], >90)
  #   lat_code ∈ {1,2,3}   (> -60, (-70,-60], <= -70)   [ordered north to south]
  #   box_id   = (lat_code - 1)*4 + lon_code  ∈ {1..12}
  #
  # NOTE: Although lon_breaks / lat_breaks are parameters, the code below uses fixed matrices.
  # They could be parameterized (keeping include.lowest/right semantics consistent).
  lon_r <- init(npp, "x")  # raster of longitudes for the npp grid
  lat_r <- init(npp, "y")  # raster of latitudes for the npp grid
  
  # classify longitude: intervals are (left, right] because right=TRUE and include.lowest=TRUE
  # So mapping is:
  #   (-Inf,-90] -> 1
  #   (-90,   0] -> 2
  #   (  0,  90] -> 3
  #   ( 90, Inf] -> 4
  lon_rcl <- matrix(c(-Inf, -90, 1,
                      -90,   0, 2,
                      0,  90, 3,
                      90, Inf, 4), byrow = TRUE, ncol = 3)
  
  # classify latitude (north to south):
  #   (-Inf,-70] -> code 3  (south of -70; "Lat70")
  #   (-70, -60] -> code 2  ( -70..-60; "Lat60")
  #   (-60, Inf] -> code 1  (> -60;    "Lat50")
  lat_rcl <- matrix(c(-Inf, -70, 3,  # south of -70 => "Lat70" (code 3)
                      -70,  -60, 2,  # -70..-60     => "Lat60" (code 2)
                      -60,  Inf, 1), byrow = TRUE, ncol = 3)  # > -60 => "Lat50" (code 1)
  
  lon_code <- classify(lon_r, lon_rcl, include.lowest = TRUE, right = TRUE)
  lat_code <- classify(lat_r, lat_rcl, include.lowest = TRUE, right = TRUE)
  
  # Combine lon_code and lat_code into a single box_id ∈ {1..12}
  # Formula: box_id = (lat_code - 1)*4 + lon_code
  box_id_r <- app(c(lon_code, lat_code), fun = function(v) (v[2] - 1L) * 4L + v[1])
  box_id   <- values(box_id_r)           # integer vector (1..12), constant across layers
  
  # ---- Extract vectors for fast, vectorized operations ----
  # (chl bins will be computed PER-LAYER below, because chl may be monthly.)
  # Note: chl_vals==NA => chl_bin==NA (propagated NA bins are intended), handled per-layer.
  
  # Prepare output raster (we will fill individual layers into this)
  npp_out <- npp
  L_npp   <- nlyr(npp_out)
  L_chl   <- nlyr(chl)
  
  # ---- Layer mapping sanity (only relevant if chl has multiple layers) ----
  if (L_chl > 1L) {
    if (layer_order == "models_fast") {
      # Expect NPP organized as: model1_m1..m12, model2_m1..m12, model3..., model4...
      # Optional: check that L_npp is a multiple of L_chl
      if ((L_npp %% L_chl) != 0L && verbose) {
        warning("nlyr(npp) is not a multiple of nlyr(chl); check layer_order or stacks.")
      }
    } else {
      # months_fast: NPP organized as: m1_model1, m1_model2, m1_model3, m1_model4, m2_model1, ...
      # Must know n_models (default 4)
      stopifnot(L_npp %% n_models == 0L)
      if ((L_npp / n_models) != L_chl && verbose) {
        warning("With 'months_fast', expected nlyr(chl) == nlyr(npp)/n_models; check inputs.")
      }
    }
  }
  
  # ---- Helper: aggregate mean NPP by (box_id, chl_bin) using observed (non-NA) NPP ----
  compute_mean_by_box_bin <- function(npp_vec, chl_bin, box_id) {
    # Mask: only use cells where NPP, chl_bin, and box_id are all available
    ok <- !is.na(npp_vec) & !is.na(chl_bin) & !is.na(box_id)
    if (!any(ok)) {
      # No observed data for this layer; return minimal structure
      return(list(keys = integer(), means = numeric(), by_box = vector("list", 12L),
                  min_bin = 0L, width = 0L, stride = 10^7L))
    }
    # Pack (box_id, chl_bin) into a single integer key for fast match/aggregate:
    #   key = box_id * stride + chl_bin
    # Choose a stride larger than the range of chl_bin to avoid collisions.
    # Here 1e7 is very conservative for typical chl ranges/bin_width.
    stride <- 10^7L
    keys <- box_id[ok] * stride + chl_bin[ok]
    
    # Aggregate mean NPP for each unique key
    # (aggregate is fine here; data.table could be faster, but this is portable)
    df <- data.frame(key = keys, npp = npp_vec[ok])
    agg <- aggregate(npp ~ key, data = df, FUN = function(x) mean(x, na.rm = TRUE))
    
    # Build per-box mean vectors indexed by local chl_bin (shifted so index starts at 1).
    # This enables fast widening searches without repeated subsetting.
    by_box <- vector("list", 12L)
    # Determine global chl_bin range to size the per-box vectors
    min_bin <- suppressWarnings(min(chl_bin, na.rm = TRUE)); if (!is.finite(min_bin)) min_bin <- 0L
    max_bin <- suppressWarnings(max(chl_bin, na.rm = TRUE)); if (!is.finite(max_bin)) max_bin <- min_bin
    width <- max(0L, max_bin - min_bin + 1L)  # number of distinct bin slots across entire raster
    
    # Initialize each box's vector to NA (to be filled where observations exist)
    for (b in 1:12) by_box[[b]] <- rep(NA_real_, width)
    
    # Unpack keys back to (box_id, chl_bin) to place means into by_box
    k_box <- floor(agg$key / stride)          # integer box_id
    k_bin <- agg$key - k_box * stride         # integer chl_bin
    for (i in seq_len(nrow(agg))) {
      b  <- k_box[i]
      bb <- k_bin[i] - min_bin + 1L          # local index into by_box[[b]]
      if (bb >= 1L && bb <= width) {
        by_box[[b]][bb] <- agg$npp[i]
      }
    }
    # Return everything needed for direct-exact and widened lookups
    list(keys = agg$key, means = agg$npp, by_box = by_box, min_bin = min_bin, width = width, stride = stride)
  }
  
  # ---- Helper: widened lookup for a set of target cells (vectorized) ----
  # For each target cell, try to find a mean NPP at (box, bin±w) for w in widen_bins.
  # If still NA, optionally "combine" box with a fallback box (e.g., 11 -> 7) at that same widened bin.
  widened_lookup <- function(target_box, target_bin, box_means, min_bin, width, widen_bins, combine_box) {
    # Work only on positions where we have both box and bin
    idx_in_box <- which(!is.na(target_box) & !is.na(target_bin))
    res <- rep(NA_real_, length(target_box))  # result vector for the subset
    if (!length(idx_in_box)) return(res)
    
    # Convert global bin to local [1..width] index (relative to min_bin)
    local_bin <- target_bin[idx_in_box] - min_bin + 1L
    ok_local <- local_bin >= 1L & local_bin <= width
    idx_in_box <- idx_in_box[ok_local]
    local_bin  <- local_bin[ok_local]
    
    # Try increasing radii defined in widen_bins (e.g., 0, 1, 2, 5)
    for (w in widen_bins) {
      # Candidate indices for left/right positions within bounds
      cand_idx1 <- pmax(1L, local_bin - w)     # left (or exact if w == 0)
      cand_idx2 <- pmin(width, local_bin + w)  # right (or exact if w == 0)
      
      # Fetch per-box means at candidate indices.
      # Prefer left/exact (r1) when available; else use right (r2).
      bvec <- target_box[idx_in_box]
      r1 <- mapply(function(b, lb) box_means[[b]][lb], bvec, cand_idx1)
      r2 <- mapply(function(b, lb) box_means[[b]][lb], bvec, cand_idx2)
      
      pick <- ifelse(!is.na(r1), r1, r2)
      
      # Box-combine rule: if still NA, try a fallback box (e.g., box 11 -> box 7)
      if (anyNA(pick) && length(combine_box)) {
        na_pos <- which(is.na(pick))
        for (p in na_pos) {
          b <- as.character(bvec[p])     # combine_box keys are character names
          if (!is.null(combine_box[[b]])) {
            b2 <- combine_box[[b]]
            v1 <- box_means[[b2]][cand_idx1[p]]
            v2 <- box_means[[b2]][cand_idx2[p]]
            pick[p] <- ifelse(!is.na(v1), v1, v2)
          }
        }
      }
      
      # Fill only positions still NA in res for this subset
      to_fill <- which(is.na(res[idx_in_box]) & !is.na(pick))
      if (length(to_fill)) {
        res[idx_in_box[to_fill]] <- pick[to_fill]
      }
      # Early exit if everything in this subset has been filled
      if (all(!is.na(res[idx_in_box]))) break
    }
    res
  }
  
  # ---- Main loop over NPP layers ----
  for (j in 1:nlyr(npp)) {
    if (verbose) message(sprintf("Layer %d/%d", j, nlyr(npp)))
    
    # -- Pick the matching chl layer for this NPP layer (supports 1-layer or multi-layer chl) --
    # If chl is a single layer (climatology), reuse it for all j.
    # If chl has multiple layers (e.g., months), map j -> k depending on `layer_order`.
    if (L_chl == 1L) {
      chl_j <- chl
    } else {
      if (layer_order == "models_fast") {
        # e.g., npp order: [cafe_1..12, cbpm_1..12, eppley_1..12, vpmg_1..12]
        # month index cycles every L_chl layers
        k <- ((j - 1L) %% L_chl) + 1L
      } else {
        # "months_fast": npp order e.g., [m1_cafe, m1_cbpm, m1_eppley, m1_vpmg, m2_cafe, ...]
        # month index advances every n_models layers
        k <- ((j - 1L) %/% n_models) %% L_chl + 1L
      }
      chl_j <- chl[[k]]
    }
    
    # -- Compute chl bins for this specific time-slice --
    chl_vals <- values(chl_j)                               # numeric, length = ncell
    chl_bin  <- as.integer(floor(chl_vals / bin_width))     # discretize chl into integer bins
    
    # Pull the j-th NPP layer as a vector for fast ops
    npp_vec <- values(npp[[j]])
    
    # Cells to fill are those where NPP is NA but both chl_bin and box_id are known
    fill_mask <- is.na(npp_vec) & !is.na(chl_bin) & !is.na(box_id)
    if (!any(fill_mask)) {
      npp_out[[j]] <- setValues(npp_out[[j]], npp_vec)
      next  # nothing to fill in this layer
    }
    
    # Precompute mean NPP by (box, bin) using observed cells in *this* layer (with its chl bins)
    agg_obj <- compute_mean_by_box_bin(npp_vec, chl_bin, box_id)
    
    # ---- Exact (box, bin) lookup using packed keys ----
    keys_all  <- box_id * agg_obj$stride + chl_bin
    filled0 <- rep(NA_real_, length(npp_vec))
    if (length(agg_obj$keys)) {
      m <- match(keys_all, agg_obj$keys)  # vectorized match
      filled0 <- ifelse(!is.na(m), agg_obj$means[m], NA_real_)
    }
    
    # ---- Widened lookup only for the cells that actually need filling ----
    target_idx <- which(fill_mask)
    widened <- widened_lookup(
      target_box = box_id[target_idx],
      target_bin = chl_bin[target_idx],
      box_means  = agg_obj$by_box,
      min_bin    = agg_obj$min_bin,
      width      = agg_obj$width,
      widen_bins = widen_bins,
      combine_box = combine_box
    )
    
    # Prefer exact match; fall back to widened value where exact is NA
    v_final <- filled0[target_idx]
    need_widen <- is.na(v_final)
    v_final[need_widen] <- widened[need_widen]
    
    # Write back filled values to the raster layer
    npp_vec[target_idx] <- v_final
    npp_out[[j]] <- setValues(npp_out[[j]], npp_vec)
  }
  # Return the filled NPP raster with all layers processed
  npp_out
}

# Fill (bin_width=0.01 reproduces your ±0.01, 0.02, 0.05 widening)
npp.filled <- fill_npp_from_chl(
  npp = npp.crop,
  chl = chl.crop,
  bin_width = 0.01,
  widen_bins = c(0L, 1L, 2L, 5L),
  # keep your special box 11 -> 7 combine rule
  combine_box = list(`11` = 7),
  n_models = 4L,
  layer_order = "models_fast",
  verbose = TRUE
)

# Save if needed
writeRaster(npp.filled, filename = paste0(env.dir,"NPP_filled_monthly_12boxes_SO.tif"), overwrite = TRUE)
writeRaster(npp.crop, filename = paste0(env.dir,"NPP_unfilled_monthly_SO.tif"), overwrite = TRUE)
# save(npp.filled,file=paste0(env.derived,file.npp,"_filled_",year,".Rdata"))
# save(npp.s.crop,file=paste0(env.derived,file.npp,"_original_",year,".Rdata"))














####################################################################
## Fill NPP values based on chl-npp relationships in 12 lon-lat regions around Antarctica
## using long-term climatology average files (1 chl-a, and 1 npp file per model)
####################################################################


fill_npp_from_chl <- function(npp, chl, lon_breaks = c(-Inf, -90, 0, 90, Inf),
                              lat_breaks = c(-Inf, -70, -60, Inf),
                              bin_width = 0.01,
                              widen_bins = c(0L, 1L, 2L, 5L),
                              combine_box = list(`11` = 7)) {
  # ---- Function goal ----
  # Fill NA NPP values using the relationship between NPP and chlorophyll-a (chl) in
  # nearby geographic boxes and similar chl values. For each (box, chl_bin) we compute
  # the mean NPP from observed data, then impute NA NPPs by matching their (box, chl_bin).
  # If no exact bin match exists, we widen the bin tolerance (±1, ±2, ±5 bins).
  #
  # Arguments:
  #   npp        : SpatRaster (≥1 layer) with NPP values (may contain NAs)
  #   chl        : SpatRaster (1 layer) with chlorophyll-a values (same CRS as npp)
  #   lon_breaks : (unused by current code) intended custom lon breaks for boxes
  #   lat_breaks : (unused by current code) intended custom lat breaks for boxes
  #   bin_width  : chl bin size; e.g. 0.01 means bins of width 0.01 chl units
  #   widen_bins : integer bin radii to widen search (0 => exact, 1 => ±0.01, 2 => ±0.02, etc.)
  #   combine_box: list mapping "sparse" boxes to a fallback box (e.g., box 11 -> 7)
  #
  # Notes:
  # - This function vectorizes everything (no per-cell scanning) for speed.
  # - Memory: values() pulls full rasters into memory; OK for your ~2M cells, but keep in mind.
  
  stopifnot(nlyr(chl) == 1, nlyr(npp) >= 1)
  
  # ---- Geometry alignment ----
  # Ensure chl matches the npp grid (same extent/resolution/nrows/ncols). We assume CRS matches.
  # If geometry differs, resample chl to npp grid (bilinear is reasonable for continuous chl).
  # If CRS differs, project() should be used prior to this function (not handled here).
  if (!compareGeom(npp, chl, stopOnError = FALSE)) {
    chl <- resample(chl, npp, method = "bilinear")
  }
  
  # ---- Build geographic box IDs (1..12) from lon/lat ----
  # We create a box raster where each cell is assigned to one of 12 boxes:
  #   lon_code ∈ {1,2,3,4} (<= -90, (-90,0], (0,90], >90)
  #   lat_code ∈ {1,2,3}   (> -60, (-70,-60], <= -70)   [ordered north to south]
  #   box_id   = (lat_code - 1)*4 + lon_code  ∈ {1..12}
  #
  # NOTE: Although lon_breaks / lat_breaks are parameters, the code below uses fixed matrices.
  # You could parameterize with these if desired (keeping include.lowest/right semantics consistent).
  lon_r <- init(npp, "x")  # raster of longitudes for the npp grid
  lat_r <- init(npp, "y")  # raster of latitudes for the npp grid
  
  # classify longitude: intervals are (left, right] because right=TRUE and include.lowest=TRUE
  # So mapping is:
  #   (-Inf,-90] -> 1
  #   (-90,   0] -> 2
  #   (  0,  90] -> 3
  #   ( 90, Inf] -> 4
  lon_rcl <- matrix(c(-Inf, -90, 1,
                      -90,   0, 2,
                      0,  90, 3,
                      90, Inf, 4), byrow = TRUE, ncol = 3)
  
  # classify latitude (north to south):
  #   (-Inf,-70] -> code 3  (south of -70; "Lat70")
  #   (-70, -60] -> code 2  ( -70..-60; "Lat60")
  #   (-60, Inf] -> code 1  (> -60;    "Lat50")
  lat_rcl <- matrix(c(-Inf, -70, 3,  # south of -70 => "Lat70" (code 3)
                      -70,  -60, 2,  # -70..-60     => "Lat60" (code 2)
                      -60,  Inf, 1), byrow = TRUE, ncol = 3)  # > -60 => "Lat50" (code 1)
  
  lon_code <- classify(lon_r, lon_rcl, include.lowest = TRUE, right = TRUE)
  lat_code <- classify(lat_r, lat_rcl, include.lowest = TRUE, right = TRUE)
  
  # Combine lon_code and lat_code into a single box_id ∈ {1..12}
  # Formula: box_id = (lat_code - 1)*4 + lon_code
  box_id_r <- app(c(lon_code, lat_code), fun = function(v) (v[2] - 1L) * 4L + v[1])
  
  # ---- Extract vectors for fast, vectorized operations ----
  chl_vals <- values(chl)                # numeric vector (length = number of cells)
  box_id   <- values(box_id_r)           # integer vector (1..12)
  chl_bin  <- as.integer(floor(chl_vals / bin_width))  # discretize chl into integer bins
  # Note: chl_vals==NA => chl_bin==NA (propagated NA bins are intended)
  
  # Prepare output raster (we will fill individual layers into this)
  npp_out <- npp
  
  # ---- Helper: aggregate mean NPP by (box_id, chl_bin) using observed (non-NA) NPP ----
  compute_mean_by_box_bin <- function(npp_vec) {
    # Mask: only use cells where NPP, chl_bin, and box_id are all available
    ok <- !is.na(npp_vec) & !is.na(chl_bin) & !is.na(box_id)
    if (!any(ok)) {
      # No observed data for this layer; return minimal structure
      return(list(keys = integer(), means = numeric(), by_box = vector("list", 12L)))
    }
    # Pack (box_id, chl_bin) into a single integer key for fast match/aggregate:
    #   key = box_id * stride + chl_bin
    # Choose a stride larger than the range of chl_bin to avoid collisions.
    # Here 1e7 is very conservative for typical chl ranges/bin_width.
    stride <- 10^7L
    keys <- box_id[ok] * stride + chl_bin[ok]
    
    # Aggregate mean NPP for each unique key
    # (aggregate is fine here; data.table could be faster, but this is portable)
    df <- data.frame(key = keys, npp = npp_vec[ok])
    agg <- aggregate(npp ~ key, data = df, FUN = function(x) mean(x, na.rm = TRUE))
    
    # Build per-box mean vectors indexed by local chl_bin (shifted so index starts at 1).
    # This enables fast widening searches without repeated subsetting.
    by_box <- vector("list", 12L)
    # Determine global chl_bin range to size the per-box vectors
    min_bin <- suppressWarnings(min(chl_bin, na.rm = TRUE)); if (!is.finite(min_bin)) min_bin <- 0L
    max_bin <- suppressWarnings(max(chl_bin, na.rm = TRUE)); if (!is.finite(max_bin)) max_bin <- min_bin
    width <- max(0L, max_bin - min_bin + 1L)  # number of distinct bin slots across entire raster
    
    # Initialize each box's vector to NA (to be filled where observations exist)
    for (b in 1:12) by_box[[b]] <- rep(NA_real_, width)
    
    # Unpack keys back to (box_id, chl_bin) to place means into by_box
    k_box <- floor(agg$key / stride)          # integer box_id
    k_bin <- agg$key - k_box * stride         # integer chl_bin
    for (i in seq_len(nrow(agg))) {
      b  <- k_box[i]
      bb <- k_bin[i] - min_bin + 1L          # local index into by_box[[b]]
      if (bb >= 1L && bb <= width) {
        by_box[[b]][bb] <- agg$npp[i]
      }
    }
    # Return everything needed for direct-exact and widened lookups
    list(keys = agg$key, means = agg$npp, by_box = by_box, min_bin = min_bin, width = width, stride = stride)
  }
  
  # ---- Helper: widened lookup for a set of target cells (vectorized) ----
  # For each target cell, try to find a mean NPP at (box, bin±w) for w in widen_bins.
  # If still NA, optionally "combine" box with a fallback box (e.g., 11 -> 7) at that same widened bin.
  widened_lookup <- function(target_box, target_bin, box_means, min_bin, width, widen_bins, combine_box) {
    # Work only on positions where we have both box and bin
    idx_in_box <- which(!is.na(target_box) & !is.na(target_bin))
    res <- rep(NA_real_, length(target_box))  # result vector for the subset
    if (!length(idx_in_box)) return(res)
    
    # Convert global bin to local [1..width] index (relative to min_bin)
    local_bin <- target_bin[idx_in_box] - min_bin + 1L
    ok_local <- local_bin >= 1L & local_bin <= width
    idx_in_box <- idx_in_box[ok_local]
    local_bin  <- local_bin[ok_local]
    
    # Try increasing radii defined in widen_bins (e.g., 0, 1, 2, 5)
    for (w in widen_bins) {
      # Candidate indices for left/right positions within bounds
      cand_idx1 <- pmax(1L, local_bin - w)     # left (or exact if w == 0)
      cand_idx2 <- pmin(width, local_bin + w)  # right (or exact if w == 0)
      
      # Fetch per-box means at candidate indices.
      # Prefer left/exact (r1) when available; else use right (r2).
      bvec <- target_box[idx_in_box]
      r1 <- mapply(function(b, lb) box_means[[b]][lb], bvec, cand_idx1)
      r2 <- mapply(function(b, lb) box_means[[b]][lb], bvec, cand_idx2)
      
      pick <- ifelse(!is.na(r1), r1, r2)
      
      # Box-combine rule: if still NA, try a fallback box (e.g., box 11 -> box 7)
      if (anyNA(pick) && length(combine_box)) {
        na_pos <- which(is.na(pick))
        for (p in na_pos) {
          b <- as.character(bvec[p])     # combine_box keys are character names
          if (!is.null(combine_box[[b]])) {
            b2 <- combine_box[[b]]
            v1 <- box_means[[b2]][cand_idx1[p]]
            v2 <- box_means[[b2]][cand_idx2[p]]
            pick[p] <- ifelse(!is.na(v1), v1, v2)
          }
        }
      }
      
      # Fill only positions still NA in res for this subset
      to_fill <- which(is.na(res[idx_in_box]) & !is.na(pick))
      if (length(to_fill)) {
        res[idx_in_box[to_fill]] <- pick[to_fill]
      }
      # Early exit if everything in this subset has been filled
      if (all(!is.na(res[idx_in_box]))) break
    }
    res
  }
  
  # ---- Main loop over NPP layers ----
  for (j in 1:nlyr(npp)) {
    message(sprintf("Layer %d/%d", j, nlyr(npp)))
    
    # Pull the j-th NPP layer as a vector for fast ops
    npp_vec <- values(npp[[j]])
    
    # Cells to fill are those where NPP is NA but both chl_bin and box_id are known
    fill_mask <- is.na(npp_vec) & !is.na(chl_bin) & !is.na(box_id)
    if (!any(fill_mask)) next  # nothing to fill in this layer
    
    # Precompute mean NPP by (box, bin) using observed cells in this layer
    agg_obj <- compute_mean_by_box_bin(npp_vec)
    
    # ---- Exact (box, bin) lookup using packed keys ----
    keys_all  <- box_id * agg_obj$stride + chl_bin
    filled0 <- rep(NA_real_, length(npp_vec))
    if (length(agg_obj$keys)) {
      m <- match(keys_all, agg_obj$keys)  # vectorized match
      filled0 <- ifelse(!is.na(m), agg_obj$means[m], NA_real_)
    }
    
    # ---- Widened lookup only for the cells that actually need filling ----
    target_idx <- which(fill_mask)
    widened <- widened_lookup(
      target_box = box_id[target_idx],
      target_bin = chl_bin[target_idx],
      box_means  = agg_obj$by_box,
      min_bin    = agg_obj$min_bin,
      width      = agg_obj$width,
      widen_bins = widen_bins,
      combine_box = combine_box
    )
    
    # Prefer exact match; fall back to widened value where exact is NA
    v_final <- filled0[target_idx]
    need_widen <- is.na(v_final)
    v_final[need_widen] <- widened[need_widen]
    
    # Write back filled values to the raster layer
    npp_vec[target_idx] <- v_final
    npp_out[[j]] <- setValues(npp_out[[j]], npp_vec)
  }
  
  # Return the filled NPP raster with all layers processed
  npp_out
}

# Clean negatives once (vectorized)
npp.s[npp.s < 0] <- NA
sum.chla[sum.chla < 0] <- NA

# Align chl to npp (resample is sufficient; both are EPSG:4326)
chl.s <- resample(sum.chla, npp.s, method = "bilinear")

# Southern Ocean crop once
so_ext <- ext(-180, 180, -90, -50)
npp.s.crop <- crop(npp.s, so_ext)
chl.crop   <- crop(chl.s, so_ext)

# Fill (bin_width=0.01 reproduces your ±0.01, 0.02, 0.05 widening)
npp.filled <- fill_npp_from_chl(
  npp = npp.s.crop,
  chl = chl.crop,
  bin_width = 0.01,
  widen_bins = c(0L, 1L, 2L, 5L),
  # keep your special box 11 -> 7 combine rule
  combine_box = list(`11` = 7)
)

# Save if needed
# writeRaster(npp.filled, filename = "NPP_filled_SO.tif", overwrite = TRUE)










