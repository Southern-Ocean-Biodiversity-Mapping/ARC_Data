############################################################
# ASAID mobiles annotation processing: counts per 2 km raster cell
#
# What this script does
# ---------------------
# 1) Reads all "Mobiles_*.csv" files exported from Squidle/ASAID.
# 2) Keeps only the fields needed to locate each annotation (lat/lon)
#    and identify its label (label.lineage_names).
# 3) Projects annotation coordinates to the raster CRS and assigns each
#    annotation to a 2x2 km raster cell (terra cell index).
# 4) Produces:
#    - A wide table: one row per raster cell, one column per label,
#      values = counts, plus total number of annotations per cell.
#    - A wide table: one row per image (point.media.key), one column per label,
#      values = counts, plus total number of annotations per image.
#    - A label frequency table: total occurrences per label across all cells.
#    - A lookup table mapping original label names -> safe column names.
#
# Notes
# -----
# - Mobiles data: no Unscorable handling is applied (totals are simple row sums).
# - No taxonomic aggregation is performed: each unique label.lineage_names
#   becomes its own column.
# - Environmental variables are NOT extracted at this stage (cell ID only).
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
bio.dir     <- paste0(usr.main.dir, "data_biological/")
env.derived <- paste0(usr.main.dir, "data_environmental/derived")

# Folder containing the downloaded ASAID annotation CSVs
annotation_dir <- paste0(bio.dir, "ASAID_image_annotations")

# File name of the raster used ONLY for cell indexing (no covariate extraction here)
raster_file <- "Circumpolar_EnvData_2km_shelf_mask_unscaled_variables.tif"

# Output folder (default: create a subfolder next to your annotation_dir)
output_dir <- file.path(annotation_dir, "derived_outputs")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Pattern used to find the MOBILES annotation files only
mobiles_pattern <- "^SquidleASAIDAnnotations_\\d{6}_Mobiles_.*\\.csv$"


############################
# 1) PACKAGES
############################

needed_pkgs <- c("data.table", "dplyr", "tidyr", "stringr", "terra")
# install.packages(setdiff(needed_pkgs, rownames(installed.packages())))
invisible(lapply(needed_pkgs, require, character.only = TRUE))


############################
# 2) LOCATE INPUT FILES
############################

mobiles_files <- list.files(annotation_dir, pattern = mobiles_pattern, full.names = TRUE)

if (length(mobiles_files) == 0) {
  stop("No mobiles-annotation files found. Check 'annotation_dir' and 'mobiles_pattern'.")
}

message("Found ", length(mobiles_files), " mobiles-annotation files.")


############################
# 3) READ AND COMBINE CSVs
############################

# Columns we need (based on the Mobiles export structure)
# (Mobiles files include extra polygon/bbox fields; we ignore them for counting.)
keep_cols <- c(
  "point.pose.lat",
  "point.pose.lon",
  "label.lineage_names",
  "point.media.key",
  "point.media.deployment.key",
  "timestamp"
)

read_one_mobiles_file <- function(f) {
  
  # Read only the columns we need, and force consistent types across files
  dt <- data.table::fread(
    f,
    select = keep_cols,
    colClasses = list(
      numeric   = c("point.pose.lat", "point.pose.lon"),
      character = c("label.lineage_names",
                    "point.media.key",
                    "point.media.deployment.key",
                    "timestamp")
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
  dt[, timestamp := as.character(timestamp)]
  
  dt
}

message("Reading and combining mobiles CSVs...")
all_mob <- data.table::rbindlist(lapply(mobiles_files, read_one_mobiles_file), use.names = TRUE)

message("Total mobiles annotations after cleaning: ", nrow(all_mob))


############################
# 4) LOAD RASTER AND ASSIGN CELL INDEX
############################

raster_path <- file.path(env.derived, raster_file)
if (!file.exists(raster_path)) {
  stop("Raster file not found:\n  ", raster_path)
}

r <- terra::rast(raster_path)
message("Raster loaded. CRS: ", terra::crs(r))

# Create a SpatVector from lon/lat in WGS84 (EPSG:4326)
pts_vec_ll <- terra::vect(
  all_mob,
  geom = c("point.pose.lon", "point.pose.lat"),
  crs = "EPSG:4326"
)

# Project points to raster CRS
pts_vec_xy <- terra::project(pts_vec_ll, terra::crs(r))

# Compute terra cell index for each annotation
cell_id <- terra::cellFromXY(r, terra::crds(pts_vec_xy))

# Attach cell_id back onto the table
all_mob[, cell_id := cell_id]

# Drop annotations outside raster extent
n_outside <- sum(is.na(all_mob$cell_id))
if (n_outside > 0) {
  message("Dropping ", n_outside, " mobiles annotations outside raster extent.")
  all_mob <- all_mob[!is.na(cell_id), ]
}

message("Mobiles annotations with valid cell_id: ", nrow(all_mob))


############################
# 5) COUNT ANNOTATIONS PER CELL PER LABEL (LONG)
############################

counts_long <- all_mob |>
  dplyr::as_tibble() |>
  dplyr::count(cell_id, label.lineage_names, name = "n")

message("Unique cells with mobiles annotations: ", dplyr::n_distinct(counts_long$cell_id))
message("Unique mobiles labels: ", dplyr::n_distinct(counts_long$label.lineage_names))


############################
# 6) MAKE SAFE COLUMN NAMES FOR WIDE OUTPUT
############################

labels <- sort(unique(counts_long$label.lineage_names))
safe_names <- make.names(labels, unique = TRUE)

label_lookup <- dplyr::tibble(
  label_lineage_names = labels,
  column_name = safe_names
)

counts_long2 <- counts_long |>
  dplyr::left_join(label_lookup, by = c("label.lineage_names" = "label_lineage_names"))


############################
# 7) WIDE TABLE: ONE ROW PER CELL, ONE COLUMN PER LABEL
############################

counts_wide_cell <- counts_long2 |>
  dplyr::select(cell_id, column_name, n) |>
  tidyr::pivot_wider(
    names_from = column_name,
    values_from = n,
    values_fill = 0
  ) |>
  dplyr::arrange(cell_id)

# Total mobiles annotations per cell (simple sum across label columns)
cell_label_cols <- setdiff(names(counts_wide_cell), "cell_id")

counts_wide_cell <- counts_wide_cell |>
  dplyr::mutate(
    n_total = if (length(cell_label_cols) > 0) {
      rowSums(dplyr::across(dplyr::all_of(cell_label_cols)))
    } else {
      0
    }
  ) |>
  dplyr::relocate(n_total, .after = cell_id)


############################
# 7b) WIDE TABLE: ONE ROW PER IMAGE, ONE COLUMN PER LABEL
############################

counts_long_image <- all_mob |>
  dplyr::as_tibble() |>
  dplyr::count(point.media.key, label.lineage_names, name = "n")

counts_long_image2 <- counts_long_image |>
  dplyr::left_join(label_lookup, by = c("label.lineage_names" = "label_lineage_names"))

counts_wide_image <- counts_long_image2 |>
  dplyr::select(point.media.key, column_name, n) |>
  tidyr::pivot_wider(
    names_from = column_name,
    values_from = n,
    values_fill = 0
  ) |>
  dplyr::arrange(point.media.key)

image_label_cols <- setdiff(names(counts_wide_image), "point.media.key")

counts_wide_image <- counts_wide_image |>
  dplyr::mutate(
    n_total = if (length(image_label_cols) > 0) {
      rowSums(dplyr::across(dplyr::all_of(image_label_cols)))
    } else {
      0
    }
  ) |>
  dplyr::relocate(n_total, .after = point.media.key)


############################
# 8) LABEL TOTALS ACROSS ALL CELLS
############################

label_totals <- counts_wide_cell |>
  dplyr::select(dplyr::all_of(cell_label_cols)) |>
  dplyr::summarise(dplyr::across(dplyr::everything(), sum)) |>
  tidyr::pivot_longer(cols = dplyr::everything(),
                      names_to = "column_name",
                      values_to = "n_total") |>
  dplyr::left_join(label_lookup, by = "column_name") |>
  dplyr::arrange(dplyr::desc(n_total))


############################
# 9) WRITE OUTPUT FILES
############################

out_cell   <- file.path(output_dir, "ASAID_mobiles_counts_by_cell_wide.csv")
out_image  <- file.path(output_dir, "ASAID_mobiles_counts_by_image_wide.csv")
out_totals <- file.path(output_dir, "ASAID_mobiles_label_totals_across_cells.csv")
out_lookup <- file.path(output_dir, "ASAID_mobiles_label_column_lookup.csv")
out_runlog <- file.path(output_dir, "ASAID_mobiles_run_summary.txt")

data.table::fwrite(counts_wide_cell, out_cell)
data.table::fwrite(counts_wide_image, out_image)
data.table::fwrite(label_totals, out_totals)
data.table::fwrite(label_lookup, out_lookup)

summary_lines <- c(
  "ASAID mobiles-annotation processing summary",
  "==========================================",
  paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  paste0("Annotation directory: ", normalizePath(annotation_dir, winslash = "/")),
  paste0("Raster used for cell indexing: ", normalizePath(raster_path, winslash = "/")),
  paste0("Number of mobiles CSVs read: ", length(mobiles_files)),
  paste0("Total mobiles annotations after cleaning: ", nrow(all_mob)),
  paste0("Mobiles annotations dropped outside raster: ", n_outside),
  paste0("Unique cells with mobiles annotations: ", dplyr::n_distinct(counts_long$cell_id)),
  paste0("Unique mobiles labels: ", dplyr::n_distinct(counts_long$label.lineage_names)),
  "",
  "Outputs written:",
  paste0(" - ", out_cell),
  paste0(" - ", out_image),
  paste0(" - ", out_totals),
  paste0(" - ", out_lookup),
  paste0(" - ", out_runlog)
)

writeLines(summary_lines, out_runlog)

message("Done.")
message("Wrote: ", out_cell)
message("Wrote: ", out_image)
message("Wrote: ", out_totals)
message("Wrote: ", out_lookup)
message("Wrote: ", out_runlog)