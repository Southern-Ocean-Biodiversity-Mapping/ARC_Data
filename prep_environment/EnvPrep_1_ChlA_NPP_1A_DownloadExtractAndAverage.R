###########################################################
#### Download, Extract and Average NPP Files from the web
###########################################################
#### This script needs to only be run once at the start of the project (about 2h runtime)
###########################################################

## specify user and setup directory to look up data from
usr <- "VM"
source("prep_environment/EnvPrep_0_SourceFile.R")
## set input and output folders
npp.dir <- paste0(usr.main.dir,"data_environmental/raw/NPP_Files")
out.dir <- paste0(usr.main.dir,"data_environmental/derived/NPP")

#########################################
#### Download NPP files from the web
## setup folders
dir.create(file.path(npp.dir, "chla"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(npp.dir, "cafe"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(npp.dir, "vpmg"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(npp.dir, "cbpm"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(npp.dir, "eppley"), recursive = TRUE, showWarnings = FALSE)

### Chl-a
save.path  <- file.path(npp.dir, "chla")
links.file <- file.path(npp.dir, "links_chla.txt")
urls <- readLines(links.file)
## download all files (fastest using "wget"), takes less than 2s per monthly file, 20s per yearly file, about 15min total
system2(
  "wget",
  c("-q", "--show-progress", "--progress=bar:force:noscroll",   # progress bar only
    "--no-check-certificate",                                   # temporary TLS bypass
    "-P", save.path,                                            # output directory
    "-i", links.file                                            # read URLs from file
  )
)
all(file.exists(file.path(save.path, basename(urls))))
setdiff(basename(urls), list.files(save.path))

### cafe
save.path  <- file.path(npp.dir, "cafe")
links.file <- file.path(npp.dir, "links_cafe.txt")
urls <- readLines(links.file)
## download all files
system2(
  "wget",
  c("-q", "--show-progress", "--progress=bar:force:noscroll",   # progress bar only
    "--no-check-certificate",                                   # temporary TLS bypass
    "-P", save.path,                                            # output directory
    "-i", links.file                                            # read URLs from file
  )
)
all(file.exists(file.path(save.path, basename(urls))))
setdiff(basename(urls), list.files(save.path))

### vpmg
save.path  <- file.path(npp.dir, "vpmg")
links.file <- file.path(npp.dir, "links_vpmg.txt")
urls <- readLines(links.file)
## download all files
system2(
  "wget",
  c("-q", "--show-progress", "--progress=bar:force:noscroll",   # progress bar only
    "--no-check-certificate",                                   # temporary TLS bypass
    "-P", save.path,                                            # output directory
    "-i", links.file                                            # read URLs from file
  )
)
all(file.exists(file.path(save.path, basename(urls))))
setdiff(basename(urls), list.files(save.path))

### cbpm
save.path  <- file.path(npp.dir, "cbpm")
links.file <- file.path(npp.dir, "links_cbpm.txt")
urls <- readLines(links.file)
## download all files
system2(
  "wget",
  c("-q", "--show-progress", "--progress=bar:force:noscroll",   # progress bar only
    "--no-check-certificate",                                   # temporary TLS bypass
    "-P", save.path,                                            # output directory
    "-i", links.file                                            # read URLs from file
  )
)
all(file.exists(file.path(save.path, basename(urls))))
setdiff(basename(urls), list.files(save.path))

### eppley
save.path  <- file.path(npp.dir, "eppley")
links.file <- file.path(npp.dir, "links_eppley.txt")
urls <- readLines(links.file)
## download all files
system2(
  "wget",
  c("-q", "--show-progress", "--progress=bar:force:noscroll",   # progress bar only
    "--no-check-certificate",                                   # temporary TLS bypass
    "-P", save.path,                                            # output directory
    "-i", links.file                                            # read URLs from file
  )
)
all(file.exists(file.path(save.path, basename(urls))))
setdiff(basename(urls), list.files(save.path))

#########################################
#### extract all files
gz <- list.files(npp.dir, pattern = "\\.gz$", full.names = TRUE, recursive = TRUE)
n  <- length(gz); if (!n) stop("No .gz files found under: ", npp.dir)
done <- 0L
cat(sprintf("Unzipping %d files...\n", n))
for (f in gz) {
  if (system2("gunzip", c("-q", "-f", f)) == 0L) done <- done + 1L  # removes .gz on success
  cat(sprintf("\r[%d/%d] %s", done, n, basename(f))); flush.console()
}
cat("\nDone.\n")

#########################################
#### create single Southern Ocean tif file for each model containing all the data across time
library(terra)
models   <- c("cafe", "vpmg", "chla", "cbpm", "eppley")
# Southern Ocean crop window
so_ext   <- ext(-180, 180, -90, -50)

### Read HDFs -> set extent/CRS -> crop to Southern Ocean -> sanitize -> write stack
# Use a dedicated temp dir for terra intermediates (prevents filename collisions)
td <- file.path(npp.dir, ".terra_tmp")
dir.create(td, showWarnings = FALSE, recursive = TRUE)
terraOptions(tempdir = td, progress = 1)

# DOYs to include (Oct–Mar)
keep_doy <- c("274","275","305","306","335","336","001","032","060","061")

# Map DOY -> month label for band names
doy_to_month <- function(d) {
  if (d %in% c("274","275")) "October"
  else if (d %in% c("305","306")) "November"
  else if (d %in% c("335","336")) "December"
  else if (d == "001")             "January"
  else if (d == "032")             "February"
  else if (d %in% c("060","061"))  "March"
  else NA_character_
}

# Read one HDF, set georef, crop, and sanitize missing/invalid values
read_crop <- function(f) {
  r <- suppressWarnings(rast(f))  # HDF lacks embedded georef
  ext(r) <- c(-180, 180, -90, 90) # global 2160×4320 (~1/12°) grid
  crs(r) <- "EPSG:4326"
  r <- crop(r, so_ext, snap = "out")
  # Sanitize in one pass: set ≤−9999 or <0 to NA
  app(r, fun = function(x) { x[x <= -9999 | x < 0] <- NA_real_; x })
}

process_model <- function(model) {
  model_dir <- file.path(npp.dir, model)
  files <- list.files(model_dir, pattern = "\\.(\\d{7})\\.hdf$", full.names = TRUE)
  if (!length(files)) {
    message("No monthly *.YYYYDDD.hdf files in: ", model_dir)
    return(invisible(NULL))
  }
  
  yyyydoy <- sub(".*\\.(\\d{7})\\.hdf$", "\\1", files)
  doys    <- substr(yyyydoy, 5, 7)
  keep    <- doys %in% keep_doy
  files   <- files[keep]; yyyydoy <- yyyydoy[keep]; doys <- doys[keep]
  o <- order(as.integer(yyyydoy))
  files <- files[o]; yyyydoy <- yyyydoy[o]; doys <- doys[o]
  
  # Read, stack, and name: "YYYYDDD_Month"
  rs  <- lapply(files, read_crop)
  stk <- rast(rs)
  years <- substr(yyyydoy, 1, 4)
  names(stk) <- paste0(years, doys, "_", vapply(doys, doy_to_month, ""))
  
  # Write the monthly time‑series stack (Oct–Mar months across all years)
  out_stack <- file.path(out.dir, paste0("NPP_monthly_OctMar_", model, ".tif"))
  writeRaster(
    stk, out_stack, overwrite = TRUE,
    wopt = list(
      gdal = c("TILED=YES", "COMPRESS=LZW", "PREDICTOR=3", "BIGTIFF=IF_SAFER"), # GDAL/GTiff  [4](https://earthdatascience.org/courses/earth-analytics/multispectral-remote-sensing-modis/modis-data-in-R/)
      datatype = "FLT4S"
    )
  )
  invisible(out_stack)
}

# Run Step 1
paths <- lapply(models, process_model)
message("Step 1 complete. Wrote stacks to: ", out.dir)

### Read stack -> select Oct-2002..Mar-2020 -> pooled summer statistic -> write GeoTIFF
# Aggregation options for the single-band summer product:
stat        <- "mean"      # "mean" or "sum"
mean_policy <- "strict"    # "available" (mean over available) or "strict" (divide by expected)

# Date window (inclusive)
start_date <- as.Date("2002-10-01")
end_date   <- as.Date("2020-03-31")

# Convert a "YYYYDDD" string to Date
date_from_yyyydoy <- function(yyyydoy) {
  yy  <- as.integer(substr(yyyydoy, 1, 4))
  doy <- as.integer(substr(yyyydoy, 5, 7))
  as.Date(doy - 1, origin = paste0(yy, "-01-01"))
}

# Build pooled summer product from a *subset* of the stack
build_summer <- function(stk) {
  ssum <- app(stk, sum, na.rm = TRUE)
  if (stat == "sum") return(ssum)
  
  if (mean_policy == "available") {
    return(app(stk, mean, na.rm = TRUE))
  } else {
    den <- nlyr(stk)
    out <- ssum / den
    avail <- app(!is.na(stk), sum, na.rm = TRUE)  # NA where all layers missing
    mask(out, avail, maskvalues = 0)
  }
}

write_summer <- function(model) {
  in_stack <- file.path(out.dir, paste0("NPP_monthly_OctMar_", model, ".tif"))
  if (!file.exists(in_stack)) {
    message("Missing stack for model: ", model)
    return(invisible(NULL))
  }
  
  stk <- rast(in_stack)
  
  # Select bands whose "YYYYDDD" fall within [2002-10-01, 2020-03-31]
  nm       <- names(stk)
  yyyydoy  <- sub("^([0-9]{7}).*$", "\\1", nm)
  dates    <- date_from_yyyydoy(yyyydoy)
  sel      <- which(dates >= start_date & dates <= end_date)
  
  if (!length(sel)) {
    message("No layers in requested window for model: ", model)
    return(invisible(NULL))
  }
  
  stk_sub <- stk[[sel]]
  summer  <- build_summer(stk_sub)
  names(summer) <- paste0("OctMar_2002To2020_", stat, "_", mean_policy)
  
  out_file <- file.path(out.dir, paste0("NPP_climatology_OctMar_2002To2020_", model, "_", stat, "_", mean_policy, ".tif"))
  writeRaster(
    summer, out_file, overwrite = TRUE,
    wopt = list(
      gdal = c("TILED=YES", "COMPRESS=LZW", "PREDICTOR=3", "BIGTIFF=IF_SAFER"),  # GDAL/GTiff  [4](https://earthdatascience.org/courses/earth-analytics/multispectral-remote-sensing-modis/modis-data-in-R/)
      datatype = "FLT4S"
    )
  )
  invisible(out_file)
}

# Run Step 2
outs <- lapply(models, write_summer)
message("Step 2 complete. Climatologies written to: ", out.dir)






