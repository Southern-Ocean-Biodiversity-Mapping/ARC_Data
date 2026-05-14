############################################################
# ASAID point-annotation processing: counts per 2 km raster cell
#
# What this script does
# ---------------------
# 1) Reads all "Points_*.csv" files exported from Squidle/ASAID.
# 2) Keeps only the fields needed to locate each annotation (lat/lon)
#    and identify its label (label.lineage_names).
# 3) Projects annotation coordinates to the raster CRS and assigns each
#    annotation to a 2x2 km raster cell (terra cell index).
# 4) Produces:
#    - A wide table: one row per raster cell, one column per label,
#      values = counts, plus totals (excluding Unscorable).
#    - A label frequency table: total occurrences per label across all cells.
#    - A lookup table mapping original label names -> safe column names.
#
# Notes
# -----
# - "Unscorable" is kept as its own column but is excluded from totals.
# - No taxonomic aggregation is performed: each unique label.lineage_names
#   becomes its own column.
# - Environmental variables are NOT extracted at this stage.
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
bio.dir      <- paste0(usr.main.dir, "data_biological/")
env.derived  <- paste0(usr.main.dir, "data_environmental/derived")

# Folder containing the downloaded ASAID annotation CSVs
annotation_dir <- paste0(bio.dir,"ASAID_image_annotations")

# File name of the raster used ONLY for cell indexing (no covariate extraction here)
raster_file <- "Circumpolar_EnvData_2km_shelf_mask_unscaled_variables.tif"

# Output folder (default: create a subfolder next to your annotation_dir)
output_dir <- file.path(annotation_dir, "derived_outputs")

# Pattern used to find the POINT annotation files only
points_pattern <- "^SquidleASAIDAnnotations_\\d{6}_Points_.*\\.csv$"

# Label that should be excluded from totals (but kept as a column)
unscorable_label <- "Unscorable"


############################
# 1) PACKAGES
############################

# Minimal set of packages for speed + clarity.
# If you don't have these installed, uncomment the install.packages() line.
needed_pkgs <- c("data.table", "dplyr", "tidyr", "stringr", "terra")

# install.packages(setdiff(needed_pkgs, rownames(installed.packages())))

invisible(lapply(needed_pkgs, require, character.only = TRUE))

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

############################
# 2) LOCATE INPUT FILES
############################
points_files <- list.files(annotation_dir, pattern = points_pattern, full.names = TRUE)

############################
# 3) READ AND COMBINE CSVs
############################

# Columns we need
keep_cols <- c(
  "point.pose.lat",
  "point.pose.lon",
  "label.lineage_names",
  "point.media.key",
  "point.media.deployment.key",          
  "point.media.deployment.campaign.key", 
  "point.media.data.equipment",          
  "timestamp"                             
)

read_one_points_file <- function(f) {
  
  # Read only the columns we need, and force consistent types across files
  dt <- data.table::fread(
    f,
    select = keep_cols,
    colClasses = list(
      numeric   = c("point.pose.lat", "point.pose.lon"),
      character = c(
        "label.lineage_names",
        "point.media.key",
        "point.media.deployment.key",
        "point.media.deployment.campaign.key",
        "point.media.data.equipment",
        "timestamp"
      )
    ),
    showProgress = FALSE
  )
  
  # Basic cleaning:
  # - Remove rows missing coordinates or label
  # - Trim whitespace in labels to avoid accidental duplicates
  dt <- dt[!is.na(point.pose.lat) &
             !is.na(point.pose.lon) &
             !is.na(label.lineage_names)]
  
  dt[, label.lineage_names := stringr::str_trim(label.lineage_names)]
  
  # Ensure timestamp is definitely character (extra belt-and-braces)
  dt[, timestamp := as.character(timestamp)]
  
  # Extract year
  dt[, year := substr(timestamp, 1, 4)]
  dt[, year := as.factor(year)]
  
  dt
}

## Reading and combining CSVs
all_pts <- data.table::rbindlist(lapply(points_files, read_one_points_file), use.names = TRUE)

message("Total annotations after basic cleaning: ", nrow(all_pts))


############################
# 4) LOAD RASTER AND ASSIGN CELL INDEX
############################

raster_path <- file.path(env.derived, raster_file)
r <- terra::rast(raster_path)
message("Raster loaded. CRS: ", terra::crs(r))

pts_vec_ll <- terra::vect(
  all_pts,
  geom = c("point.pose.lon", "point.pose.lat"),
  crs = "EPSG:4326"
)

pts_vec_xy <- terra::project(pts_vec_ll, terra::crs(r))

cell_id <- terra::cellFromXY(r, terra::crds(pts_vec_xy))
all_pts[, cell_id := cell_id]

n_outside <- sum(is.na(all_pts$cell_id))
if (n_outside > 0) {
  message("Dropping ", n_outside, " annotations outside raster extent.")
  all_pts <- all_pts[!is.na(cell_id), ]
}

message("Annotations with valid raster cell_id: ", nrow(all_pts))

############################
# 4b) CELL-LEVEL METADATA
############################

cell_metadata <- all_pts |>
  dplyr::as_tibble() |>
  dplyr::select(
    cell_id,
    surveyID   = point.media.deployment.campaign.key,
    transectID = point.media.deployment.key,
    gear       = point.media.data.equipment,
    year
  ) |>
  dplyr::group_by(cell_id) |>
  dplyr::summarise(
    surveyID   = dplyr::first(surveyID),
    transectID = dplyr::first(transectID),
    gear       = dplyr::first(gear),
    year       = dplyr::first(year),
    .groups = "drop"
  )

############################
# 5) COUNT ANNOTATIONS PER CELL PER LABEL (LONG)
############################

counts_long <- all_pts |>
  dplyr::as_tibble() |>
  dplyr::count(cell_id, label.lineage_names, name = "n")

############################
# 6) MAKE SAFE COLUMN NAMES
############################

labels <- sort(unique(counts_long$label.lineage_names))
safe_names <- make.names(labels, unique = TRUE)

label_lookup <- dplyr::tibble(
  label_lineage_names = labels,
  column_name = safe_names
)

counts_long2 <- counts_long |>
  dplyr::left_join(
    label_lookup,
    by = c("label.lineage_names" = "label_lineage_names")
  )

############################
# 7) WIDE TABLE: ONE ROW PER CELL
############################

counts_wide <- counts_long2 |>
  dplyr::select(cell_id, column_name, n) |>
  tidyr::pivot_wider(
    names_from  = column_name,
    values_from = n,
    values_fill = 0
  ) |>
  dplyr::arrange(cell_id) |>
  dplyr::left_join(cell_metadata, by = "cell_id")

############################
# 7a) DEFINE TRUE LABEL COLUMNS (FIXES ALL ERRORS)
############################

# Label columns are NUMERIC columns excluding cell_id
label_cols <- names(counts_wide)[sapply(counts_wide, is.numeric)]
label_cols <- setdiff(label_cols, "cell_id")

# Identify unscorable column safely
unscorable_col <- intersect(
  label_cols,
  label_lookup$column_name[
    label_lookup$label_lineage_names == unscorable_label
  ]
)

label_cols_excl_unscorable <- setdiff(label_cols, unscorable_col)

############################
# 7b) TOTALS PER CELL (SAFE WITH |>)
############################
# Compute totals explicitly
tot_excl <- if (length(label_cols_excl_unscorable) > 0) {
  rowSums(
    as.matrix(counts_wide[, label_cols_excl_unscorable, drop = FALSE]),
    na.rm = TRUE
  )
} else {
  rep(0, nrow(counts_wide))
}

tot_incl <- if (length(label_cols) > 0) {
  rowSums(
    as.matrix(counts_wide[, label_cols, drop = FALSE]),
    na.rm = TRUE
  )
} else {
  rep(0, nrow(counts_wide))
}

# Attach totals
counts_wide$cover_points_N <- tot_incl
counts_wide$cover_points_scorable <- tot_excl

# Reorder for readability
counts_wide <- counts_wide |>
  dplyr::relocate(
    cover_points_N,
    cover_points_scorable,
    .after = cell_id
  )

counts_wide <- counts_wide |>
  dplyr::relocate(surveyID, transectID, gear, year, .after = cell_id)

############################
# 7c) WIDE TABLE: ONE ROW PER IMAGE
############################

counts_long_image <- all_pts |>
  dplyr::as_tibble() |>
  dplyr::count(point.media.key, label.lineage_names, name = "n")

counts_wide_image <- counts_long_image |>
  dplyr::left_join(
    label_lookup,
    by = c("label.lineage_names" = "label_lineage_names")
  ) |>
  dplyr::select(point.media.key, column_name, n) |>
  tidyr::pivot_wider(
    names_from  = column_name,
    values_from = n,
    values_fill = 0
  ) |>
  dplyr::arrange(point.media.key)

############################
# 8) LABEL TOTALS ACROSS ALL CELLS
############################
# Sum each label column across all cells.
label_totals <- counts_wide |>
  dplyr::select(dplyr::all_of(label_cols)) |>
  dplyr::summarise(dplyr::across(dplyr::everything(), sum)) |>
  tidyr::pivot_longer(cols = dplyr::everything(),
                      names_to = "column_name",
                      values_to = "n_total") |>
  dplyr::left_join(label_lookup, by = "column_name") |>
  dplyr::arrange(dplyr::desc(n_total))

############################
# 9) WRITE OUTPUT FILES
############################

out_counts <- file.path(output_dir, "ASAID_points_counts_by_cell_wide.csv")
out_totals <- file.path(output_dir, "ASAID_label_totals_across_cells.csv")
out_lookup <- file.path(output_dir, "ASAID_label_column_lookup.csv")
out_runlog <- file.path(output_dir, "ASAID_run_summary.txt")
out_counts_image <- file.path(output_dir, "ASAID_points_counts_by_image_wide.csv")

data.table::fwrite(counts_wide, out_counts)
data.table::fwrite(label_totals, out_totals)
data.table::fwrite(label_lookup, out_lookup)
data.table::fwrite(counts_wide_image, out_counts_image)

# Write a small run summary (handy for reproducibility)
summary_lines <- c(
  "ASAID point-annotation processing summary",
  "========================================",
  paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  paste0("Annotation directory: ", normalizePath(annotation_dir, winslash = "/")),
  paste0("Raster used for cell indexing: ", normalizePath(raster_path, winslash = "/")),
  paste0("Number of point CSVs read: ", length(points_files)),
  paste0("Total annotations after cleaning: ", nrow(all_pts)),
  paste0("Annotations dropped outside raster: ", n_outside),
  paste0("Unique cells with annotations: ", dplyr::n_distinct(counts_long$cell_id)),
  paste0("Unique labels: ", dplyr::n_distinct(counts_long$label.lineage_names)),
  paste0("Unscorable label: ", unscorable_label),
  paste0("Unscorable column name (if present): ", ifelse(length(unscorable_col) == 1, unscorable_col, "NOT PRESENT")),
  "",
  "Outputs written:",
  paste0(" - ", out_counts),
  paste0(" - ", out_totals),
  paste0(" - ", out_lookup)
)

writeLines(summary_lines, out_runlog)

message("Done.")
message("Wrote: ", out_counts)
message("Wrote: ", out_totals)
message("Wrote: ", out_lookup)
message("Wrote: ", out_runlog)
message("Wrote: ", out_counts_image)