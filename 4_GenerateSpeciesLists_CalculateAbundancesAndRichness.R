############################################################
# Purpose
# -------
# Prepare biological response variables from ASAID cover-point
# annotation data for modelling.
#
# This script:
# 1) Loads cell-level cover data merged with environmental predictors
# 2) Summarises label prevalence and exports for manual curation
# 3) Applies user-defined label exclusions and merges
# 4) Generates response matrices:
#    - counts per label
#    - presence/absence per label
#    - proportional cover (optional)
# 5) Combines responses with scaled environmental predictors
# 6) Saves a structured object for downstream modelling workflows
#
# Inputs
# ------
# - cover_cells_env_2km.rds
# - cover_cell_metadata_env_scaled_2km.rds
# - scaling_params_cover_2km.rds
# - Curation Excel file (user-generated after reviewing prevalence)
#
# Outputs
# -------
# - cover_prevalence_summary_2km.xlsx
# - cover_modelling_inputs_2km.rds
#
############################################################


############################
# 0) USER SETTINGS
############################

usr <- "JJ"
source("0_SourceFile.R")

bio.dir      <- paste0(usr.main.dir, "data_biological/")
output_dir   <- paste0(usr.main.dir, "data_products/modelling_files/circum_antarctic")

res <- "2km"

cover_cells_env_file <- file.path(output_dir, paste0("cover_cells_env_", res, ".rds"))
cell_meta_scaled_file <- file.path(output_dir, paste0("cover_cell_metadata_env_scaled_", res, ".rds"))
cell_meta_file <- file.path(output_dir, paste0("cover_cell_metadata_env_", res, ".rds"))
scaling_params_file   <- file.path(output_dir, paste0("scaling_params_cover_", res, ".rds"))

prevalence_xlsx_out <- file.path(output_dir, paste0("cover_prevalence_summary_", res, ".xlsx"))
modelling_output_rds <- file.path(output_dir, paste0("cover_modelling_inputs_", res, ".rds"))

curation_xlsx <- file.path(output_dir, paste0("cover_label_curation_", res, ".xlsx"))

# Column names in curation sheet
cur_label_col   <- "Label"
cur_merge_col   <- "Merge_With_2pc_2km"
cur_exclude_col <- "Exclude_2km"

exclude_flag <- "x"

categorical_vars <- c("geomorphology")


############################
# 1) PACKAGES
############################

library(dplyr)
library(tidyr)
library(data.table)
library(readxl)
library(writexl)
library(purrr)


############################
# 2) LOAD DATA
############################

cover_cells_env <- readRDS(cover_cells_env_file) |> as_tibble()
cell_meta <- readRDS(cell_meta_file) |> as_tibble()
cell_meta_scaled <- readRDS(cell_meta_scaled_file) |> as_tibble()
scaling_params <- readRDS(scaling_params_file)

predictor_cols <- scaling_params$predictors

# Ensure only sampled cells
cover_cells_env <- cover_cells_env |> filter(cover_points_N > 0)


############################
# 3) IDENTIFY RESPONSE COLUMNS
############################
meta_cols <- c("cell_id", "surveyID", "transectID", "gear", "year",
               "proj_coord_x", "proj_coord_y", "lon", "lat")
total_cols <- c(
  "cover_points_N",
  "cover_points_scorable"
)

# ALSO explicitly remove predictors
predictor_cols <- scaling_params$predictors

non_label_cols <- unique(c(meta_cols, total_cols, predictor_cols))

label_cols <- setdiff(names(cover_cells_env), non_label_cols)

# keep only numeric columns
label_cols <- label_cols[sapply(cover_cells_env[label_cols], is.numeric)]

# SAFETY: remove obvious environmental leftovers
label_cols <- label_cols[!grepl(
  "mean|sd|tpi|depth|temp|sal|arag|po4|no3",
  label_cols
)]

if (length(label_cols) == 0) {
  stop("No label columns detected.")
}

############################
# 3b) CELL-LEVEL BIODIVERSITY METRICS
############################
# Select ONLY biological labels
label_cols_biota <- label_cols[grepl("^X1\\.1\\.Biota", label_cols)]

# Remove unscorable if present
label_cols_biota <- setdiff(
  label_cols_biota,
  grep("Unscorable", label_cols_biota, value = TRUE)
)

# Calculate species richness (presence/absence)
# and total abundance (sum of counts)
cover_cells_env <- cover_cells_env |>
  mutate(
    richness_raw = rowSums(
      across(all_of(label_cols_biota), ~ . > 0),
      na.rm = TRUE
    ),
    
    total_abundance = rowSums(
      across(all_of(label_cols_biota)),
      na.rm = TRUE
    )
  )

message("Calculated richness and total abundance (biota only, pre-curation).")

############################
# 4) PREVALENCE SUMMARY
############################

n_cells <- nrow(cover_cells_env)

summary_df <- data.frame(
  Label = label_cols,
  count_2km = sapply(label_cols, function(v)
    sum(cover_cells_env[[v]] > 0, na.rm = TRUE)
  ),
  prev_2km = sapply(label_cols, function(v)
    round(sum(cover_cells_env[[v]] > 0, na.rm = TRUE) / n_cells, 3)
  )
)

# Create full template structure for aggregating rare morphospecies labels
template_df <- summary_df %>%
  mutate(
    count_img = NA,
    prev_img = NA,
    Exclude_2km = "",
    Merge_With_2pc_2km = "",
    Notes = "",
    Final_labels_2km = ""
  ) %>%
  select(
    Label,
    count_img,
    prev_img,
    count_2km,
    prev_2km,
    Exclude_2km,
    Merge_With_2pc_2km,
    Notes,
    Final_labels_2km
  ) %>%
  arrange(Label)

# Write Excel file
writexl::write_xlsx(
  list(COVER_naming = template_df),
  prevalence_xlsx_out
)

message("Template curation file written to: ", prevalence_xlsx_out)


############################
# 5) LOAD AND APPLY CURATION
############################
if (!file.exists(curation_xlsx)) {
  stop(
    "Curation file not found.\n",
    "Review prevalence summary and create file at:\n", curation_xlsx
  )
}

# Load and standardise curation table
curation <- read_xlsx(curation_xlsx) |>
  as_tibble() |>
  mutate(across(everything(), ~ as.character(.))) |>
  filter(!is.na(Label)) |>
  filter(Label != "X")   # remove accidental Excel artifact

# Keep only biological labels
label_cols <- label_cols[grepl("^X1\\.1\\.Biota", label_cols)]
label_cols <- setdiff(label_cols, grep("Unscorable", label_cols, value = TRUE))

## LONG FORMAT
long_df <- cover_cells_env |>
  select(cell_id, all_of(label_cols)) |>
  pivot_longer(
    cols = -cell_id,
    names_to = "Label",
    values_to = "count"
  ) |>
  left_join(curation, by = "Label")

## APPLY EXCLUSIONS + MERGES
has_exclude <- cur_exclude_col %in% names(long_df)
has_merge   <- cur_merge_col %in% names(long_df)

long_df <- long_df |>
  mutate(
    exclude = if (has_exclude) {
      !is.na(.data[[cur_exclude_col]]) &
        tolower(trimws(.data[[cur_exclude_col]])) == exclude_flag
    } else {
      FALSE
    },
    
    new_label = case_when(
      exclude ~ NA_character_,
      has_merge &
        !is.na(.data[[cur_merge_col]]) &
        .data[[cur_merge_col]] != "" ~ .data[[cur_merge_col]],
      TRUE ~ Label
    )
  ) |>
  filter(!is.na(new_label))

## BACK TO WIDE FORMAT
cover_curated <- long_df |>
  group_by(cell_id, new_label) |>
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(
    id_cols = cell_id,
    names_from = new_label,
    values_from = count,
    values_fill = 0
  )

## simplify labels (first remove start, then remove "...")
names(cover_curated)[-1] <- gsub("X1\\.1\\.Biota\\...", "", names(cover_curated)[-1])
names(cover_curated)[-1] <- gsub("\\.\\.\\.", "_", names(cover_curated)[-1])
## replace "..3D" with "3D", and "..2D." with "2D"
names(cover_curated)[-1] <- gsub("\\.\\.3D", "3D", names(cover_curated)[-1])
names(cover_curated)[-1] <- gsub("\\.\\.2D.", "2D", names(cover_curated)[-1])

## REATTACH METADATA
meta_keep <- cover_cells_env |>
  select(
    cell_id,
    any_of(c(
      "surveyID",
      "transectID",
      "gear",
      "year",
      "proj_coord_x",
      "proj_coord_y",
      "lon",
      "lat",
      "cover_points_N",
      "cover_points_scorable",
      "richness_raw",
      "total_abundance"
    ))
  ) |>
  distinct()

cover_curated <- left_join(meta_keep, cover_curated, by = "cell_id")

############################
# 6) RESPONSE MATRICES
############################

response_cols <- setdiff(names(cover_curated), names(meta_keep))
response_cols <- response_cols[sapply(cover_curated[response_cols], is.numeric)]

# Counts
resp_counts <- cover_curated |> select(cell_id, all_of(response_cols))

# Presence/absence
resp_pa <- resp_counts
resp_pa[response_cols] <- lapply(resp_pa[response_cols], function(x) as.integer(x > 0))

# Proportions
resp_prop <- resp_counts
resp_prop[response_cols] <- lapply(response_cols, function(v) {
  resp_counts[[v]] / cover_curated$cover_points_scorable
})
names(resp_prop)[-1] <- response_cols


############################
# 7) FINAL MODELLING DATASET
############################

predictors <- cell_meta |>
  select(cell_id, all_of(predictor_cols)) |>
  mutate(across(any_of(categorical_vars), as.factor))

predictors_scaled <- cell_meta_scaled |>
  select(cell_id, all_of(predictor_cols)) |>
  mutate(across(any_of(categorical_vars), as.factor))

modelling_df <- resp_counts |>
  left_join(meta_keep, by = "cell_id") |>
  left_join(predictors_scaled, by = "cell_id")


############################
# 8) SAVE OUTPUT
############################
output <- list(
  metadata = list(
    resolution = res,
    predictors = predictor_cols,
    response_labels = response_cols,
    n_cells = nrow(resp_counts)
  ),
  responses = list(
    counts = resp_counts,
    presence_absence = resp_pa,
    proportions = resp_prop
  ),
  cell_metrics = meta_keep,   # Cell-level biodiversity summaries
  predictors = predictors,  # Environmental predictors
  predictors_scaled = predictors_scaled,  # Environmental predictors
  modelling_dataframe = modelling_df,
  prevalence_table = summary_df
)

saveRDS(output, modelling_output_rds)

message("Script completed successfully.")
message("Outputs saved to: ", modelling_output_rds)
message("Prevalence summary saved to: ", prevalence_xlsx_out)