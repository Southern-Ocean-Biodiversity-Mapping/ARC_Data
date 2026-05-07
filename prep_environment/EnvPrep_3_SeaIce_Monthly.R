# ============================================================
# Robust monthly Antarctic sea-ice concentration download + GeoTIFF export
# Option B: Use ERDDAP /files listings and download source NetCDFs directly
#
# Output GeoTIFFs:
#   SeaIce_monthly_YYYY_MM.tif
#
# Data sources (Southern Hemisphere / Antarctic):
#   1) Final CDR V4 monthly SH: nsidcG02202v4shmday [1](https://www.ecmwf.int/en/forecasts/dataset/ecmwf-reanalysis-v5)
#   2) NRT CDR V2 monthly SH : nsidcG10016v2shmday [2](https://data.nasa.gov/dataset/amsr-e-amsr2-unified-l3-daily-12-5-km-brightness-temperatures-sea-ice-concentration-motion-64a89)[3](https://nsidc.org/data/data-access-tool/AU_SI12/versions/1)
#
# Important note:
#   ERDDAP warns that variable names/metadata in source files may differ from griddap outputs. [1](https://www.ecmwf.int/en/forecasts/dataset/ecmwf-reanalysis-v5)
# ============================================================

suppressPackageStartupMessages({
  library(terra)
  library(ncdf4)
  library(lubridate)
})

# ----------------------------
# USER SETTINGS
# ----------------------------
seaice.dir <- "C:/Users/jjansen/UTAS Research Dropbox/Jan Jansen/Data/data_environmental/raw/SeaIce_monthly/"
dir.create(seaice.dir, recursive = TRUE, showWarnings = FALSE)

# Where to store downloaded NetCDFs and exported GeoTIFFs
nc_dir  <- file.path(seaice.dir, "source_nc")
tif_dir <- file.path(seaice.dir, "monthly_tif")
dir.create(nc_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(tif_dir, recursive = TRUE, showWarnings = FALSE)

# Start month requested
start_ym <- "2015-12"

# Output units: set TRUE to write % (0–100), FALSE to keep as fraction if provided that way
write_percent <- FALSE

# Download behaviour
overwrite_nc  <- FALSE
overwrite_tif <- TRUE

# ----------------------------
# ERDDAP "files" directory URLs (Southern Hemisphere / Antarctic)
# These pages list .nc files in the directory. [1](https://www.ecmwf.int/en/forecasts/dataset/ecmwf-reanalysis-v5)
# ----------------------------
files_v4_dir  <- "https://polarwatch.noaa.gov/erddap/files/nsidcG02202v4shmday/?C=N;O=D"
files_nrt_dir <- "https://polarwatch.noaa.gov/erddap/files/nsidcG10016v2shmday/?C=N;O=D"

# Base URLs for direct file download (IMPORTANT: must end with "/")
files_v4_base  <- "https://polarwatch.noaa.gov/erddap/files/nsidcG02202v4shmday/"
files_nrt_base <- "https://polarwatch.noaa.gov/erddap/files/nsidcG10016v2shmday/"

# ----------------------------
# Helper: robustly extract .nc filenames from the ERDDAP /files HTML listing
# We use rvest to parse HTML and then keep only clean filenames ending in .nc
# ----------------------------
get_nc_filenames <- function(listing_url) {
  if (!requireNamespace("xml2", quietly = TRUE)) install.packages("xml2")
  if (!requireNamespace("rvest", quietly = TRUE)) install.packages("rvest")
  
  doc <- xml2::read_html(listing_url)
  hrefs <- rvest::html_elements(doc, "a") |> rvest::html_attr("href")
  hrefs <- hrefs[!is.na(hrefs)]
  
  # Keep only plain .nc filenames (drop directories and any query strings/fragments)
  nc <- hrefs[grepl("\\.nc$", hrefs, ignore.case = TRUE)]
  nc <- sub("^.*/", "", nc)      # ensure just filename
  nc <- sub("\\?.*$", "", nc)    # drop query string if any
  nc <- sub("#.*$", "", nc)      # drop fragment if any
  
  unique(nc)
}

# ----------------------------
# Helper: parse YYYYMM from known filename patterns
# ----------------------------
parse_yyyymm <- function(fname) {
  if (grepl("^seaice_conc_monthly_sh_\\d{6}_", fname)) {
    return(sub("^seaice_conc_monthly_sh_(\\d{6}).*$", "\\1", fname))
  }
  if (grepl("^seaice_conc_monthly_icdr_sh_\\d{6}_", fname)) {
    return(sub("^seaice_conc_monthly_icdr_sh_(\\d{6}).*$", "\\1", fname))
  }
  NA_character_
}

# ----------------------------
# Helper: validate the downloaded file is NetCDF (not HTML)
# ----------------------------
is_valid_netcdf <- function(path) {
  if (!file.exists(path) || file.info(path)$size < 10000) return(FALSE)
  
  # Quick HTML sniff
  first <- readLines(path, n = 2, warn = FALSE)
  if (length(first) > 0 && grepl("<!DOCTYPE html>|<html", first[1], ignore.case = TRUE)) return(FALSE)
  
  # NetCDF open test
  ok <- try(ncdf4::nc_open(path), silent = TRUE)
  if (inherits(ok, "try-error")) return(FALSE)
  ncdf4::nc_close(ok)
  TRUE
}

# ----------------------------
# Helper: download using curl (handles redirects well)
# ----------------------------
download_file <- function(url, dest, overwrite = FALSE) {
  if (!overwrite && file.exists(dest) && file.info(dest)$size > 0) return(dest)
  
  if (!requireNamespace("curl", quietly = TRUE)) install.packages("curl")
  
  # download
  ok <- try(curl::curl_download(url, destfile = dest, quiet = TRUE, mode = "wb"), silent = TRUE)
  if (inherits(ok, "try-error")) {
    if (file.exists(dest)) file.remove(dest)
    return(NA_character_)
  }
  
  # validate
  if (!is_valid_netcdf(dest)) {
    if (file.exists(dest)) file.remove(dest)
    return(NA_character_)
  }
  
  dest
}

# ----------------------------
# Helper: pick SIC variable inside NetCDF robustly
# (Source-file variable names may differ from ERDDAP griddap names.) [1](https://www.ecmwf.int/en/forecasts/dataset/ecmwf-reanalysis-v5)
# ----------------------------
pick_sic_var <- function(ncfile) {
  nc <- ncdf4::nc_open(ncfile)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  
  v <- names(nc$var)
  if (length(v) == 0) stop("No variables found in: ", ncfile)
  
  # prefer variables containing seaice/ice + conc
  cand <- v[grepl("ice|seaice", v, ignore.case = TRUE) & grepl("conc", v, ignore.case = TRUE)]
  if (length(cand) == 0) cand <- v
  
  # among candidates, prefer monthly if present
  cand2 <- cand[grepl("month", cand, ignore.case = TRUE)]
  if (length(cand2) > 0) cand <- cand2
  
  # prefer cdr/icdr if present
  cand2 <- cand[grepl("cdr|icdr", cand, ignore.case = TRUE)]
  if (length(cand2) > 0) cand <- cand2
  
  cand[[1]]
}

# ============================================================
# 1) Build catalogue from listings
# ============================================================

cat("Reading ERDDAP listings...\n")

v4_names  <- get_nc_filenames(files_v4_dir)
nrt_names <- get_nc_filenames(files_nrt_dir)

# Build data frames with URLs and dates
mk_df <- function(fnames, base_url, source_tag) {
  yyyymm <- vapply(fnames, parse_yyyymm, character(1))
  keep <- !is.na(yyyymm)
  fnames <- fnames[keep]
  yyyymm <- yyyymm[keep]
  
  dates <- as.Date(paste0(yyyymm, "01"), format = "%Y%m%d")
  
  data.frame(
    source = source_tag,
    fname  = fnames,
    url    = paste0(base_url, fnames),  # <-- DIRECT FILE URL
    yyyymm = yyyymm,
    date   = dates,
    stringsAsFactors = FALSE
  )
}

df_v4  <- mk_df(v4_names,  files_v4_base,  "CDR_V4_SH")
df_nrt <- mk_df(nrt_names, files_nrt_base, "NRT_V2_SH")

df_all <- rbind(df_v4, df_nrt)

# Filter to requested date range
start_date <- as.Date(paste0(start_ym, "-01"))
end_date   <- floor_date(Sys.Date(), "month") - days(1)

df_all <- df_all[df_all$date >= start_date & df_all$date <= end_date, ]
if (nrow(df_all) == 0) stop("No monthly NetCDF files found in requested range.")

# Prefer final CDR V4 up to 2024-12, else prefer NRT
cutoff <- as.Date("2024-12-31")

df_pref <- do.call(rbind, lapply(split(df_all, df_all$date), function(d) {
  if (d$date[1] <= cutoff) {
    d2 <- d[d$source == "CDR_V4_SH", , drop = FALSE]
    if (nrow(d2) > 0) return(d2[1, , drop = FALSE])
  }
  d2 <- d[d$source == "NRT_V2_SH", , drop = FALSE]
  if (nrow(d2) > 0) return(d2[1, , drop = FALSE])
  d[1, , drop = FALSE]
}))

df_pref <- df_pref[order(df_pref$date), ]
keep_months <- c("12","01","02","03")
df_pref <- df_pref[format(df_pref$date, "%m") %in% keep_months, ]
cat("Months to download after DJFM filter: ", nrow(df_pref), "\n")

# Sanity: show first few URLs to confirm they include the filename
cat("\nFirst 3 download URLs:\n")
print(head(df_pref$url, 3))

# ============================================================
# 2) Download + convert to GeoTIFF
# ============================================================

terraOptions(progress = 1, memfrac = 0.6)

for (i in seq_len(nrow(df_pref))) {
  this <- df_pref[i, ]
  yy <- format(this$date, "%Y")
  mm <- format(this$date, "%m")
  
  nc_dest  <- file.path(nc_dir,  this$fname)
  tif_dest <- file.path(tif_dir, sprintf("SeaIce_monthly_%s_%s.tif", yy, mm))
  
  cat(sprintf("\n[%d/%d] %s (%s)\n", i, nrow(df_pref), format(this$date, "%Y-%m"), this$source))
  
  # Download and validate
  nc_local <- download_file(this$url, nc_dest, overwrite = overwrite_nc)
  if (is.na(nc_local)) {
    warning("Failed (or invalid NetCDF) for: ", this$url)
    next
  }
  
  # Choose variable and load raster
  sic_var <- pick_sic_var(nc_local)
  r <- try(terra::rast(nc_local, subds = sic_var), silent = TRUE)
  if (inherits(r, "try-error")) {
    r <- terra::rast(nc_local)
    if (terra::nlyr(r) > 1) r <- r[[1]]
  }
  
  # Optional conversion to percent
  if (write_percent) {
    rng <- terra::global(r, "range", na.rm = TRUE)
    if (!is.na(rng[1, 2]) && rng[1, 2] <= 1.5) r <- r * 100
  }
  
  # Write GeoTIFF
  terra::writeRaster(
    r, tif_dest,
    overwrite = overwrite_tif,
    datatype = "FLT4S",
    NAflag = -9999,
    gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2")
  )
  
  cat("Wrote: ", tif_dest, "\n")
}

cat("\nDone.\nGeoTIFFs in: ", tif_dir, "\nNetCDFs in: ", nc_dir, "\n")


################################################################################
################################################################################


# ============================================================
# Average sea-ice concentration maps (East Antarctica)
# Uses your GeoTIFFs: SeaIce_monthly_YYYY_MM.tif
#
# Strategy (robust):
#   1) Read monthly GeoTIFFs (projected grid)
#   2) Reproject each raster to lon/lat (EPSG:4326)
#   3) Crop to East Antarctica lon/lat box
#   4) Compute mean maps (overall and/or per calendar month)
# ============================================================

library(terra)

# ----------------------------
# USER SETTINGS
# ----------------------------

tif_dir <- "C:/Users/jjansen/UTAS Research Dropbox/Jan Jansen/Data/data_environmental/raw/SeaIce_monthly/monthly_tif"
out_dir <- file.path(tif_dir, "east_antarctica_means")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# East Antarctica definition (simple, reproducible lon/lat box)
# Adjust if you want a different sector
lon_min <- 100
lon_max <- 155
lat_min <- -69
lat_max <- -63   # try -60 for more polar-only

# Months of interest (Dec–Mar). Use 1:12 for all months.
target_months <- c(12, 1, 2, 3)

# Output unit preference: TRUE => percent 0–100, FALSE => keep native (often 0–1)
write_percent <- TRUE

# Reprojection resolution (degrees). 0.25° is plenty for 25-km source data.
# If you want smoother maps, use 0.1; if you want faster, use 0.5.
target_res_deg <- 0.25

# ----------------------------
# 1) LIST FILES + PARSE DATES
# ----------------------------

f <- list.files(tif_dir, pattern = "^SeaIce_monthly_\\d{4}_\\d{2}\\.tif$", full.names = TRUE)
stopifnot(length(f) > 0)

ym <- sub("^SeaIce_monthly_(\\d{4})_(\\d{2})\\.tif$", "\\1-\\2", basename(f))
dates <- as.Date(paste0(ym, "-01"))
mos <- as.integer(format(dates, "%m"))

# Sort chronologically just to be safe
ord <- order(dates)
f <- f[ord]; dates <- dates[ord]; mos <- mos[ord]

# ----------------------------
# 2) SET UP LON/LAT TEMPLATE (target grid)
# ----------------------------

ea_ext <- ext(lon_min, lon_max, lat_min, lat_max)

# Create a blank lon/lat raster template at the chosen resolution
tmpl <- rast(
  extent = ea_ext,
  resolution = target_res_deg,
  crs = "EPSG:4326"
)

# ----------------------------
# 3) HELPER: READ → OPTIONAL UNIT FIX → PROJECT → CROP
# ----------------------------

read_project_crop <- function(path) {
  r <- rast(path)
  
  # Some products are fraction (0–1), others percent (0–100).
  # We detect by max value and optionally convert to percent for consistency.
  rng <- global(r, "range", na.rm = TRUE)
  vmax <- rng[1, 2]
  
  if (write_percent) {
    # If vmax looks like fraction, convert to percent
    if (!is.na(vmax) && vmax <= 1.5) r <- r * 100
  }
  
  # Reproject to lon/lat template (bilinear is fine for concentration fields)
  r_ll <- project(r, tmpl, method = "bilinear")
  
  # Crop to East Antarctica box
  crop(r_ll, ea_ext)
}

# ----------------------------
# 4) BUILD STACK (cropped lon/lat rasters)
# ----------------------------

cat("Reprojecting + cropping", length(f), "rasters...\n")
r_list <- lapply(f, read_project_crop)

# Combine into one SpatRaster time series
r_ea <- rast(r_list)
time(r_ea) <- dates

# Sanity check: do we have valid cells?
valid_cells <- global(r_ea[[1]], "notNA")
if (valid_cells[1,1] == 0) {
  stop("East Antarctica crop produced 0 valid cells. Try adjusting lon/lat bounds or target_res_deg.")
}

# ----------------------------
# 5) OUTPUT A) OVERALL MEAN (all months/years)
# ----------------------------

mean_all <- app(r_ea, mean, na.rm = TRUE)
out_all <- file.path(out_dir, "SeaIce_mean_ALL_months_EastAnt.tif")
writeRaster(mean_all, out_all, overwrite = TRUE, datatype = "FLT4S",
            gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2"))
cat("Wrote:", out_all, "\n")

# ----------------------------
# 6) OUTPUT B) MEAN PER CALENDAR MONTH (Dec/Jan/Feb/Mar)
# ----------------------------

for (m in target_months) {
  idx <- which(mos == m)
  if (length(idx) == 0) {
    cat("No layers found for month", m, "\n")
    next
  }
  
  mean_m <- app(r_ea[[idx]], mean, na.rm = TRUE)
  
  out_m <- file.path(out_dir, sprintf("SeaIce_mean_month_%02d_EastAnt.tif", m))
  writeRaster(mean_m, out_m, overwrite = TRUE, datatype = "FLT4S",
              gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2"))
  cat("Wrote:", out_m, "(n =", length(idx), "months)\n")
}

cat("\nDone. Outputs in:\n", out_dir, "\n")

############################################################################################
#### plotting
tif.dir <- "C:/Users/jjansen/UTAS Research Dropbox/Jan Jansen/Data/data_environmental/raw/SeaIce_monthly/monthly_tif/"
## mean values
ice12 <- rast(paste0(tif.dir,"east_antarctica_means/SeaIce_mean_month_12_EastAnt.tif"))
ice01 <- rast(paste0(tif.dir,"east_antarctica_means/SeaIce_mean_month_01_EastAnt.tif"))
ice02 <- rast(paste0(tif.dir,"east_antarctica_means/SeaIce_mean_month_02_EastAnt.tif"))
ice03 <- rast(paste0(tif.dir,"east_antarctica_means/SeaIce_mean_month_03_EastAnt.tif"))

par(mfrow=c(4,1))
plot(ice12, main="Mean SIC - December")
plot(ice01, main="Mean SIC - January")
plot(ice02, main="Mean SIC - February")
plot(ice03, main="Mean SIC - March")

## monthly data from the last 5 years
xlim=c(100,155)
ylim=c(-69,-63)
tifs <- list.files(tif.dir, full.names = TRUE)
ice_12 <- rast(tifs[grep("_12.tif", tifs)])
ice_01 <- rast(tifs[grep("_01.tif", tifs)])
ice_02 <- rast(tifs[grep("_02.tif", tifs)])
ice_03 <- rast(tifs[grep("_03.tif", tifs)])

## reproject
ice_12_ll <- project(ice_12, "EPSG:4326")
ice_01_ll <- project(ice_01, "EPSG:4326")
ice_02_ll <- project(ice_02, "EPSG:4326")
ice_03_ll <- project(ice_03, "EPSG:4326")


par(mfrow=c(5,1))
plot(ice_01_ll[[6]], main="SIC - January 2021", xlim=xlim, ylim=ylim)
plot(ice_01_ll[[7]], main="SIC - January 2022", xlim=xlim, ylim=ylim)
plot(ice_01_ll[[8]], main="SIC - January 2023", xlim=xlim, ylim=ylim)
plot(ice_01_ll[[9]], main="SIC - January 2024", xlim=xlim, ylim=ylim)
plot(ice_01_ll[[10]],main="SIC - January 2025", xlim=xlim, ylim=ylim)

par(mfrow=c(5,1))
plot(ice_02_ll[[6]], main="SIC - February 2021", xlim=xlim, ylim=ylim)
plot(ice_02_ll[[7]], main="SIC - February 2022", xlim=xlim, ylim=ylim)
plot(ice_02_ll[[8]], main="SIC - February 2023", xlim=xlim, ylim=ylim)
plot(ice_02_ll[[9]], main="SIC - February 2024", xlim=xlim, ylim=ylim)
plot(ice_02_ll[[10]],main="SIC - February 2025", xlim=xlim, ylim=ylim)

par(mfrow=c(5,1))
plot(ice_03_ll[[6]], main="SIC - March 2021", xlim=xlim, ylim=ylim)
plot(ice_03_ll[[7]], main="SIC - March 2022", xlim=xlim, ylim=ylim)
plot(ice_03_ll[[8]], main="SIC - March 2023", xlim=xlim, ylim=ylim)
plot(ice_03_ll[[9]], main="SIC - March 2024", xlim=xlim, ylim=ylim)
plot(ice_03_ll[[10]],main="SIC - March 2025", xlim=xlim, ylim=ylim)





par(mfrow=c(5,1))
plot(ice_12_ll[[6]], main="SIC - December 2020", xlim=xlim, ylim=ylim)
plot(ice_12_ll[[7]], main="SIC - December 2021", xlim=xlim, ylim=ylim)
plot(ice_12_ll[[8]], main="SIC - December 2022", xlim=xlim, ylim=ylim)
plot(ice_12_ll[[9]], main="SIC - December 2023", xlim=xlim, ylim=ylim)
plot(ice_12_ll[[10]],main="SIC - December 2024", xlim=xlim, ylim=ylim)

























# ============================================================
# Robust: Mean monthly SIC (East Antarctica) + min/max 5% edge lines
# Input: SeaIce_monthly_YYYY_MM.tif (GeoTIFFs you downloaded & created)
# Output: PNG per month in /plots_east_antarctica
# ============================================================

library(terra)

# ----------------------------
# USER SETTINGS
# ----------------------------
tif_dir <- "C:/Users/jjansen/UTAS Research Dropbox/Jan Jansen/Data/data_environmental/raw/SeaIce_monthly/monthly_tif"
plot_dir <- file.path(tif_dir, "plots_east_antarctica")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# Months to process. Set 1:12 for all months.
target_months <- c(12, 1, 2, 3)

# East Antarctica sector (simple lon/lat box)
lon_min <- 0
lon_max <- 180
lat_min <- -90
lat_max <- -60   # use -60 if you want more polar only

# ----------------------------
# 1) LOAD FILES + PARSE DATES
# ----------------------------
f <- list.files(tif_dir, pattern = "^SeaIce_monthly_\\d{4}_\\d{2}\\.tif$", full.names = TRUE)
stopifnot(length(f) > 0)

ym <- sub("^SeaIce_monthly_(\\d{4})_(\\d{2})\\.tif$", "\\1-\\2", basename(f))
dates <- as.Date(paste0(ym, "-01"))
yrs <- as.integer(format(dates, "%Y"))
mos <- as.integer(format(dates, "%m"))

r_all <- rast(f)
time(r_all) <- dates

# ----------------------------
# 2) MASK TO EAST ANTARCTICA
# ----------------------------
east_poly_ll <- vect(
  matrix(c(
    lon_min, lat_min,
    lon_max, lat_min,
    lon_max, lat_max,
    lon_min, lat_max,
    lon_min, lat_min
  ), ncol = 2, byrow = TRUE),
  type = "polygons",
  crs = "EPSG:4326"
)

east_poly <- project(east_poly_ll, crs(r_all))
r_east <- mask(r_all, east_poly)

# ----------------------------
# 3) AUTO-DETECT SIC UNITS + THRESHOLD
# ----------------------------
rng <- global(r_east[[1]], "range", na.rm = TRUE)
vmax <- rng[1, 2]

thr <- if (!is.na(vmax) && vmax <= 1.5) 0.05 else 5
message("Detected max value ~", round(vmax, 3), " => using 5% threshold = ", thr)

# ----------------------------
# 4) HELPERS WITH SAFETY CHECKS
# ----------------------------

# Count non-NA cells in a layer
n_valid <- function(layer) {
  # global(..., "notNA") returns count of non-NA cells
  as.numeric(global(layer, "notNA", na.rm = FALSE))
}

# Compute extent area (km^2) where SIC >= threshold
extent_area_km2 <- function(layer, threshold) {
  if (n_valid(layer) == 0) return(NA_real_)
  m <- layer >= threshold
  a <- cellSize(layer, unit = "m")  # m^2
  as.numeric(global(ifel(m, a, NA), "sum", na.rm = TRUE)) / 1e6
}

# Safe contour extraction: returns NULL if it cannot make a contour
safe_edge_line <- function(layer, threshold) {
  if (n_valid(layer) == 0) return(NULL)
  
  # If threshold is above max or below min (within valid cells), contour will be empty
  rr <- global(layer, "range", na.rm = TRUE)
  if (is.na(rr[1,1]) || is.na(rr[1,2])) return(NULL)
  if (threshold < rr[1,1] || threshold > rr[1,2]) return(NULL)
  
  out <- try(as.contour(layer, levels = threshold), silent = TRUE)
  if (inherits(out, "try-error")) return(NULL)
  if (!inherits(out, "SpatVector")) return(NULL)
  if (nrow(out) == 0) return(NULL)
  
  out
}

# ----------------------------
# 5) LOOP MONTHS
# ----------------------------
terraOptions(progress = 1, memfrac = 0.6)

for (m in target_months) {
  
  idx <- which(mos == m)
  if (length(idx) < 2) {
    message("Skipping month ", m, " (need >=2 years; found ", length(idx), ")")
    next
  }
  
  r_m_all <- r_east[[idx]]
  d_m_all <- dates[idx]
  
  # ---- Drop layers that are all-NA after masking ----
  valid_counts <- vapply(1:nlyr(r_m_all), function(i) n_valid(r_m_all[[i]]), numeric(1))
  keep <- which(valid_counts > 0)
  
  if (length(keep) < 2) {
    message("Skipping month ", sprintf("%02d", m),
            ": after masking, <2 layers have data (keep=", length(keep), ").")
    next
  }
  
  r_m <- r_m_all[[keep]]
  d_m <- d_m_all[keep]
  
  message("\nMonth ", sprintf("%02d", m),
          ": layers=", nlyr(r_m),
          " (dropped ", length(idx) - length(keep), " empty layer(s))")
  
  # ---- Mean SIC for this calendar month ----
  mean_m <- app(r_m, fun = mean, na.rm = TRUE)
  
  # ---- Area at 5% threshold for each layer ----
  areas <- vapply(1:nlyr(r_m), function(i) extent_area_km2(r_m[[i]], thr), numeric(1))
  
  # If all NA (shouldn’t happen after keep filter), skip
  if (all(is.na(areas))) {
    message("Skipping month ", sprintf("%02d", m), ": all extent areas are NA after masking.")
    next
  }
  
  # Remove NA areas for min/max selection
  ok <- which(!is.na(areas))
  r_ok <- r_m[[ok]]
  d_ok <- d_m[ok]
  areas_ok <- areas[ok]
  
  i_min <- which.min(areas_ok)
  i_max <- which.max(areas_ok)
  
  layer_min <- r_ok[[i_min]]
  layer_max <- r_ok[[i_max]]
  d_min <- d_ok[i_min]
  d_max <- d_ok[i_max]
  
  message("  Min extent: ", format(d_min, "%Y-%m"),
          " area=", round(areas_ok[i_min]), " km^2")
  message("  Max extent: ", format(d_max, "%Y-%m"),
          " area=", round(areas_ok[i_max]), " km^2")
  
  # ---- Safe contour lines (may return NULL if no line) ----
  line_min <- safe_edge_line(layer_min, thr)
  line_max <- safe_edge_line(layer_max, thr)
  
  if (is.null(line_min)) message("  Note: no MIN contour line produced (likely no 5% crossing or empty).")
  if (is.null(line_max)) message("  Note: no MAX contour line produced (likely no 5% crossing or empty).")
  
  # ---- Plot ----
  out_png <- file.path(plot_dir, sprintf("EastAnt_meanSIC_minmaxEdge_%02d.png", m))
  png(out_png, width = 1800, height = 1600, res = 200)
  
  plot(
    mean_m,
    main = sprintf(
      "East Antarctica mean SIC (month=%02d)\nMin: %s | Max: %s | Threshold=%s",
      m, format(d_min, "%Y-%m"), format(d_max, "%Y-%m"), thr
    )
  )
  
  # overlay lines only if they exist
  if (!is.null(line_min)) plot(line_min, add = TRUE, col = "deepskyblue3", lwd = 2)
  if (!is.null(line_max)) plot(line_max, add = TRUE, col = "red3",         lwd = 2)
  
  legend(
    "bottomleft",
    legend = c(
      sprintf("Min 5%% extent (%s): %.0f km^2", format(d_min, "%Y-%m"), areas_ok[i_min]),
      sprintf("Max 5%% extent (%s): %.0f km^2", format(d_max, "%Y-%m"), areas_ok[i_max])
    ),
    col = c("deepskyblue3", "red3"),
    lwd = 2,
    bty = "n"
  )
  
  dev.off()
  message("Wrote: ", out_png)
}

message("\nAll done. Plots in: ", plot_dir)