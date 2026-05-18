############################################################
# Purpose
# -------
# Starting from a "wide" cover-point annotation tables
# (one row per cell/image, one column per label), this script:
#
# 1) Loads environmental rasters (unscaled variables + polynomials/etc)
# 2) Extracts environmental values at sampled cells (2 km)
# 3) Merges environmental predictors with cover-point response data
# 4) Scales predictors (z-score) and saves scaling parameters
# 5) Creates a scaled environmental raster stack for prediction
# 6) Creates a dataframe of scaled values for all shelf cells (for prediction)
#
# Inputs (NEW format)
# -------------------
# - ASAID_cover_points_counts_by_cell_wide.csv
# - ASAID_cover_points_counts_by_image_wide.csv (optional, loaded but not required)
# - Circumpolar_EnvData_2km_shelf_mask_unscaled_variables.tif
# - Circumpolar_EnvData_2km_shelf_mask_unscaled_polynomials_etc.tif
#
# Outputs
# -------
# - cover_cell_metadata_env.rds
# - cover_cells_env.rds
# - cover_cell_metadata_env_scaled.rds
# - scaling_params_cover_2km.rds
# - Circumpolar_EnvData_2km_shelf_mask_scaled.tif
# - Circumpolar_EnvData_2km_shelf_mask_scaled_dataframe.rds
#
############################################################

############################
# 0) USER SETTINGS (EDIT ME)
############################

## specify user and setup directory to look up data from
#usr <- "VM"
#usr <- "SJ"
usr <- "JJ"
source("0_SourceFile.R")

## set folders
bio.dir      <- paste0(usr.dropbox.dir, "data_biological/")
env.derived  <- paste0(usr.dropbox.dir, "data_environmental/derived")
## Output folder for merged and scaled datasets (make sure this exists or is created before running)
output_dir <- paste0(usr.dropbox.dir,"data_products/modelling_files/circum_antarctic")
## Folder containing the downloaded ASAID annotation CSVs
annotation_dir <- paste0(bio.dir,"ASAID_image_annotations")
## Folder where annotation outputs are written
annotation_output_dir <- paste0(bio.dir,"ASAID_image_annotations/derived_outputs")

# Resolution label used in environmental filenames
res <- "2km"

# Environmental raster filenames (unscaled)
env_file_vars <- paste0("Circumpolar_EnvData_", res, "_shelf_mask_unscaled_variables.tif")
env_file_poly <- paste0("Circumpolar_EnvData_", res, "_shelf_mask_unscaled_polynomials_etc.tif")

# NEW annotation wide tables (cover points)
cover_cell_counts_file  <- file.path(annotation_output_dir, "ASAID_points_counts_by_cell_wide.csv")
cover_image_counts_file <- file.path(annotation_output_dir, "ASAID_points_counts_by_image_wide.csv")  # optional


# Label column that should be excluded from scorable totals (kept as a column)
unscorable_label_name <- "Unscorable"  # the *original* label text

# Optional: correlation threshold for auto-dropping predictors (set NULL to skip auto-drop)
cor_threshold <- NULL   # e.g., 0.7; NULL means "do not auto-drop"

# manual list of environmental variables to remove (to mimic your earlier workflow)
env_remove <- c("tpi5", "tpi11", "arag_mean", "no3_mean", "no3_sd", "po4_mean", "po4_sd")

# Predictors that are categorical codes and must NOT be z-scaled
categorical_vars <- c("geomorphology")




############################
# 1) PACKAGES
############################

needed_pkgs <- c("data.table", "dplyr", "terra")
# install.packages(setdiff(needed_pkgs, rownames(installed.packages())))

invisible(lapply(needed_pkgs, require, character.only = TRUE))


############################
# 2) LOAD COVER COUNTS (CELL + IMAGE)
############################
cover_cells_wide <- data.table::fread(cover_cell_counts_file)
message("Loaded cover cell-wide table: ", nrow(cover_cells_wide), " rows; ", ncol(cover_cells_wide), " columns.")

# Image-wide table is optional for this script (not required for cell-level extraction)
cover_images_wide <- NULL
if (file.exists(cover_image_counts_file)) {
  cover_images_wide <- data.table::fread(cover_image_counts_file)
  message("Loaded cover image-wide table: ", nrow(cover_images_wide), " rows; ", ncol(cover_images_wide), " columns.")
} else {
  message("Note: cover image-wide file not found (this is OK for cell-level env extraction).")
}

############################
# 3) LOAD ENVIRONMENTAL STACK
############################

env_path_vars <- file.path(env.derived, env_file_vars)
env_path_poly <- file.path(env.derived, env_file_poly)

env_stack <- terra::rast(c(env_path_vars, env_path_poly))
message("Loaded env_stack with ", terra::nlyr(env_stack), " layers.")


############################
# 4) CREATE CELL METADATA (COORDS) + EXTRACT ENV BY CELL ID
############################
# Use the FIRST layer as a reference raster for geometry (extent/res/CRS)
r_ref <- env_stack[[1]]

# Get xy coordinates (projected) for each sampled cell
cell_ids <- cover_cells_wide$cell_id

xy <- terra::xyFromCell(r_ref, cell_ids)  # returns matrix with columns x,y
cell_meta <- data.frame(
  cell_id = cell_ids,
  proj_coord_x = xy[,1],
  proj_coord_y = xy[,2]
)

# Convert projected coords to lon/lat for convenience
pts_xy <- terra::vect(cell_meta, geom = c("proj_coord_x", "proj_coord_y"), crs = terra::crs(r_ref))
pts_ll <- terra::project(pts_xy, "EPSG:4326")
ll <- terra::crds(pts_ll)

cell_meta$lon <- ll[,1]
cell_meta$lat <- ll[,2]

# Extract environmental values by cell index (fast & robust)
env_vals <- terra::extract(env_stack, cell_ids)

cell_metadata_env <- cbind(cell_meta, env_vals)

message("Extracted environmental values for ", nrow(cell_metadata_env), " sampled cells.")


############################
# 5) MERGE ENV + COVER RESPONSE (CELL-LEVEL)
############################

# Merge environmental predictors onto cover cell-wide table (sampled subset)
cover_cells_env <- cover_cells_wide |>
  dplyr::left_join(cell_metadata_env, by = "cell_id")

# Save unscaled outputs
saveRDS(cell_metadata_env, file.path(output_dir, paste0("cover_cell_metadata_env_", res, ".rds")))
saveRDS(cover_cells_env,    file.path(output_dir, paste0("cover_cells_env_", res, ".rds")))

message("Saved unscaled merged objects to: ", output_dir)


############################
# 6) SELECT PREDICTORS + SCALE (Z-SCORE), BUT DO NOT SCALE CATEGORICALS
############################

# Identify which columns are environmental predictors:
metadata_cols <- c("cell_id", "proj_coord_x", "proj_coord_y", "lon", "lat")
env_cols <- setdiff(names(cell_metadata_env), metadata_cols)

# Remove any variables in env_remove (if present)
env_cols_kept <- setdiff(env_cols, intersect(env_cols, env_remove))

# Ensure categorical vars are kept (unless you explicitly removed them)
categorical_vars <- intersect(categorical_vars, env_cols_kept)

# Optional: auto-drop highly correlated predictors (numeric-only, exclude categoricals)
if (!is.null(cor_threshold)) {
  message("Auto-dropping predictors with |cor| > ", cor_threshold, " (pairwise), numeric-only.")
  numeric_candidates <- setdiff(env_cols_kept, categorical_vars)
  x <- cell_metadata_env[, numeric_candidates, drop = FALSE]
  # keep only numeric columns (extra safety)
  is_num <- sapply(x, is.numeric)
  x <- x[, is_num, drop = FALSE]
  cmat <- stats::cor(x, use = "pairwise.complete.obs")
  # greedy drop: remove later variable in each high-correlation pair
  high <- which(abs(cmat) > cor_threshold & upper.tri(cmat), arr.ind = TRUE)
  drop_vars <- character(0)
  if (nrow(high) > 0) drop_vars <- unique(colnames(cmat)[high[, 2]])
  env_cols_kept <- setdiff(env_cols_kept, drop_vars)
  message("Dropped ", length(drop_vars), " predictors by correlation threshold.")
}

# Scale (z-score): (x - mean)/sd for numeric continuous predictors only
cell_metadata_env_scaled <- cell_metadata_env

scale_means <- rep(NA_real_, length(env_cols_kept))
scale_sds   <- rep(NA_real_, length(env_cols_kept))
names(scale_means) <- env_cols_kept
names(scale_sds)   <- env_cols_kept

for (v in env_cols_kept) {
  
  # Do NOT scale categorical predictors
  if (v %in% categorical_vars) {
    scale_means[v] <- NA_real_
    scale_sds[v]   <- NA_real_
    next
  }
  
  mu <- mean(cell_metadata_env_scaled[[v]], na.rm = TRUE)
  sdv <- sd(cell_metadata_env_scaled[[v]], na.rm = TRUE)
  scale_means[v] <- mu
  scale_sds[v]   <- sdv
  
  if (is.finite(sdv) && sdv > 0) {
    cell_metadata_env_scaled[[v]] <- (cell_metadata_env_scaled[[v]] - mu) / sdv
  } else {
    cell_metadata_env_scaled[[v]] <- 0
  }
}

# Keep a record of predictors included for modelling:
sel_predictors <- env_cols_kept

# Save scaled metadata + scaling parameters
saveRDS(cell_metadata_env_scaled, file.path(output_dir, paste0("cover_cell_metadata_env_scaled_", res, ".rds")))
saveRDS(
  list(
    predictors = sel_predictors,
    categorical_vars = categorical_vars,
    means = scale_means,
    sds = scale_sds,
    env_remove = env_remove,
    cor_threshold = cor_threshold
  ),
  file.path(output_dir, paste0("scaling_params_cover_", res, ".rds"))
)

message("Saved scaled cell metadata and scaling parameters (categoricals left unscaled).")

############################
# 7) SCALE THE ENVIRONMENTAL RASTER STACK FOR PREDICTIONS (DISK-BASED)
#    - does NOT scale categorical predictors (e.g., geomorphology)
############################

# Subset env_stack to the selected predictors (only those present in the raster)
sel_ra <- which(names(env_stack) %in% sel_predictors)
pred_stack_nc <- terra::subset(env_stack, sel_ra)

# terra temp settings
terra_tmp <- file.path(output_dir, "terra_tmp")
if (!dir.exists(terra_tmp)) dir.create(terra_tmp, recursive = TRUE)
terra::terraOptions(tempdir = terra_tmp, memfrac = 0.3)

# Process in groups (4 layers at a time)
idx_groups <- split(1:terra::nlyr(pred_stack_nc),
                    ceiling(seq_along(1:terra::nlyr(pred_stack_nc)) / 4))
tmp_files <- character(0)
for (g in seq_along(idx_groups)) {
  ## select the layers to process
  idx <- idx_groups[[g]]
  grp <- pred_stack_nc[[idx]]
  grp_scaled <- grp
  ## loop across each layer in the group and scale (except categoricals)
  for (j in seq_along(idx)) {
    lyr_name <- names(grp)[j]
    # Do NOT scale categorical layers
    if (lyr_name %in% categorical_vars) {
      grp_scaled[[j]] <- grp[[j]]
      next
    }
    ## load mean and sd values
    mu  <- scale_means[lyr_name]
    sdv <- scale_sds[lyr_name]
    ## scale
    if (!is.na(mu) && !is.na(sdv) && is.finite(sdv) && sdv > 0) {
      grp_scaled[[j]] <- (grp[[j]] - mu) / sdv
    } else {
      grp_scaled[[j]] <- grp[[j]] * 0
    }
  }
  ## write out the layers temporarily
  f_out <- file.path(env.derived, paste0("tmp_scaled_", res, "_group_", g, ".tif"))
  terra::writeRaster(grp_scaled, f_out, overwrite = TRUE,
                     wopt = list(gdal = c("COMPRESS=LZW")))
  tmp_files <- c(tmp_files, f_out)
  
  message("Wrote temp scaled file: ", f_out)
}

# Combine temp files into final scaled stack (file-backed)
scaled_raster_out <- file.path(env.derived, paste0("Circumpolar_EnvData_", res, "_shelf_mask_scaled.tif"))
scaled_all <- terra::rast(tmp_files)
terra::writeRaster(scaled_all, filename = scaled_raster_out, overwrite = TRUE,
                   wopt = list(gdal = c("COMPRESS=LZW")))

message("Wrote scaled prediction raster stack:\n  ", scaled_raster_out)

# Optional cleanup
# file.remove(tmp_files)


############################
# 8) CREATE A SCALED PREDICTION DATAFRAME (ONE ROW PER CELL)
############################

# Use the scaled stack we just wrote (or pred_stack_scaled in memory)
env_stack_scaled <- terra::rast(scaled_raster_out)

# Define "valid shelf cells" as those with non-NA depth (common convention)
# Adjust if your mask uses a different layer.
if (!("depth" %in% names(env_stack_scaled))) {
  stop("Expected a 'depth' layer to define valid shelf cells, but it is not present in the scaled stack.")
}

valid_cells <- which(!is.na(env_stack_scaled[["depth"]][]))

# Extract all scaled predictor values for valid cells
# (This can be memory-heavy; for 2km it is usually manageable.)
pred_vals <- terra::extract(env_stack_scaled, valid_cells)

pred_xy <- terra::xyFromCell(env_stack_scaled[[1]], valid_cells)
pred_df <- data.frame(
  cell_id = valid_cells,
  proj_coord_x = pred_xy[,1],
  proj_coord_y = pred_xy[,2],
  pred_vals
)

# Save dataframe for prediction workflows
pred_df_out <- file.path(env.derived, paste0("Circumpolar_EnvData_", res, "_shelf_mask_scaled_dataframe.rds"))
saveRDS(pred_df, pred_df_out)

message("Saved scaled prediction dataframe:\n  ", pred_df_out)


############################
# 9) DONE
############################

message("All steps completed successfully.")
message("Env outputs directory:\n  ", output_dir)