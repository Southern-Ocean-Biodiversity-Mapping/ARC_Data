# run_full_pipeline_3D_to_2D.R
# Full pipeline: 3D sinking/advection -> bottom NPP -> 2D seafloor transport/sedimentation -> outputs
#
# 3D indexing upgrade:
# - index_mode_3d = "xy+z" uses 2D kNN + vertical-layer selection (fast, low memory)
# - na_action_3d  = "stop" stops particles that encounter non-finite u/v/w or bathymetry (prevents crashes)

suppressWarnings({
  library(terra)
  library(spatstat.geom)
  library(spatstat.random)
  library(ppmData)
})

# -------------------- USER SETTINGS --------------------

ptrackr_dir <- "/pvol3TB/Scripts/ptrackR_package_update2026/R/"
env_dir3 <- "/pvol3TB/data_environmental/ROMS_2k_files/"
out_dir  <- "/pvol3TB/FAM_outputs/"
cache_dir <- "/pvol3TB/tmp/ptrack_cache"
presence_tmp_dir <- "/pvol3TB/tmp"

npp.model.list <- c("cafe_12boxfilled",
                    "eppl_12boxfilled",
                    "vpmg_12boxfilled",
                    "cbpm_12boxfilled")

for(z in 1:4){
  npp.model <- npp.model.list[z]
# # NPP model label used in filenames
# npp.model <- "cafe_12boxfilled"
# #npp.model <- "eppl_12boxfilled"
# #npp.model <- "vpmg_12boxfilled"
# #npp.model <- "cbpm_12boxfilled"
# #npp.model <- "chla"

# Shared parameters
speed <- 200
time_steps_in_s <- 1800
regions <- 1:10

# 3D parameters
roms_slices_3d <- 84
looping_time_3d <- 0.25
runtime_3d <- 1
ROMSreduced_3d <- TRUE

# 3D indexing/NA behaviour
index_mode_3d <- "xy+z"   # "xy+z" (recommended) or "knn3d" (legacy)
na_action_3d  <- "stop"   # "stop" (recommended), "zero_velocity", or "error"

# 2D parameters
roms_slices_2d <- 113
looping_time_2d <- 0.25
runtime_2d <- 1
sedimentation_2d <- TRUE
particle_radius <- 0.0002
# only the decimal part for the character string
r_str <- sub("^0\\.", "", format(particle_radius, scientific = FALSE, trim = TRUE))

# Flux mode: "entry" (default), "total", "presence" (disk-backed exact P1)
flux_mode <- "entry"

# Presence-only parameters (only used for flux_mode="presence")
presence_buckets <- 512L
presence_flush_n <- 2e6
presence_cleanup <- TRUE

# PPP scaling divisor (as in your original script)
ppp_divisor <- 10000000 #5000000

# Reduce ROMS depth layers before 2D (as in your script)
reduce_depth_layers <- TRUE
layers_to_drop <- 4:31

# -------------------- SOURCING FUNCTIONS --------------------

source(file.path(ptrackr_dir, "setup_knn.R"))
source(file.path(ptrackr_dir, "buildparams.R"))
source(file.path(ptrackr_dir, "create_points_pattern.R"))

source(file.path(ptrackr_dir, "presence_flux_disk.R"))
source(file.path(ptrackr_dir, "trackit_2D.R"))

# 3D indexing utilities + upgraded 3D tracker
source(file.path(ptrackr_dir, "indexing_3D_xy_z.R"))
source(file.path(ptrackr_dir, "trackit_3D.R"))

# Wrapper (must include pass-through for index_mode_3d/na_action_3d; see file below)
source(file.path(ptrackr_dir, "loopit_2D3D_traj_ROMSreduced.R"))

# Analysis helpers
source(file.path(ptrackr_dir, "analysis_helpers_flux_sed.R"))
#source(file.path(ptrackr_dir, "analysis_helpers_compile.R"))
source(file.path(ptrackr_dir, "analysis_helpers_3D_backtransform.R"))
source(file.path(ptrackr_dir, "analysis_helpers_io_cache.R"))

# Ensure terra2im exists (qrbp/ppmData lineage) [1](https://rdrr.io/github/skiptoniam/qrbp/man/terra2im.html)[2](https://github.com/skiptoniam/ppmData)
if (!exists("terra2im", mode = "function")) {
  stop("terra2im() not found. Load the package that provides terra2im() before running (e.g., qrbp).")
}

# -------------------- STEP 1: 3D TRACKING PER REGION --------------------

# for (i in regions) {
#   message("===== 3D tracking: Region ", sprintf("%02d", i), " =====")
# 
#   # Load ROMS (expects Rdat6h)
#   load(fname_roms_region(env_dir3, i))
# 
#   # Load NPP particle start points (pts9: x,y,z,npp)
#   npp_pts_file <- paste0(env_dir3,
#                          "ocean_his_TrackingSetup_", npp.model,
#                          "_NPP9_Region", sprintf("%02d", i), ".Rdata")
#   load(npp_pts_file)
# 
#   start.time <- Sys.time()
#   track.3D <- loopit_2D3D(
#     pts_seeded = as.matrix(pts9[, 1:3]),
#     romsobject = Rdat6h,
#     projected = TRUE,
#     speed = speed,
#     domain = "3D",
#     runtime = runtime_3d,
#     roms_slices = roms_slices_3d,
#     looping_time = looping_time_3d,
#     ROMSreduced = ROMSreduced_3d,
#     time_steps_in_s = time_steps_in_s,
#     trajectories = FALSE,
#     index_mode_3d = index_mode_3d,
#     na_action_3d = na_action_3d
#   )
#   message("3D tracking took: ", format(Sys.time() - start.time))
# 
#   out_3d_file <- paste0(out_dir,
#                         "/tracking3D_Region", sprintf("%02d", i),
#                         "_", npp.model, "_NPP9_",
#                         speed, "mday_21days_red.Rdata")
#   save(track.3D, file = out_3d_file)
# 
#   rm(Rdat6h, pts9, track.3D)
#   gc(FALSE)
# }
# 
# # -------------------- STEP 2: 3D BACKTRANSFORM TO BOTTOM NPP RASTERS --------------------
# 
# for (i in regions) {
#   message("===== 3D backtransform: Region ", sprintf("%02d", i), " =====")
# 
#   # Reload pts9 to get particle NPP values
#   npp_pts_file <- paste0(env_dir3,
#                          "ocean_his_TrackingSetup_", npp.model,
#                          "_NPP9_Region", sprintf("%02d", i), ".Rdata")
#   load(npp_pts_file)
# 
#   # Load template NPP raster for region
#   npp_crop <- terra::rast(fname_npp_raster_region(env_dir3, i))
# 
#   # Load 3D tracking output
#   in_3d_file <- paste0(out_dir,
#                        "/tracking3D_Region", sprintf("%02d", i),
#                        "_", npp.model, "_NPP9_",
#                        speed, "mday_21days_red.Rdata")
#   load(in_3d_file)
# 
#   xyz_end <- track.3D$xyz_end
#   colnames(xyz_end) <- c("x", "y", "z")
# 
#   bt <- backtransform_3D_to_bottom_raster(
#     xyz_end = xyz_end,
#     npp_values = pts9[, 4],
#     npp_template = npp_crop
#   )
# 
#   xyz_end_npp <- bt$xyz_end_npp
#   cell_npp    <- bt$cell_npp
# 
#   save(xyz_end_npp, cell_npp,
#        file = paste0(out_dir,
#                      "/tracking3D_Region", sprintf("%02d", i),
#                      "_", npp.model, "_NPP9_",
#                      speed, "mday_21days_xyzbottomnpp.Rdata"))
# 
#   bottom_out <- paste0(out_dir,
#                        "/tracking3D_Region", sprintf("%02d", i),
#                        "_", npp.model,
#                        "_NPP9_", speed, "mday_21days_bottomnpp.tif")
#   terra::writeRaster(bt$bottom_npp_raster, bottom_out, overwrite = TRUE)
# 
#   rm(bt, track.3D, pts9, xyz_end, npp_crop)
#   gc(FALSE)
# }

# -------------------- STEP 3: 2D TRACKING FROM BOTTOM NPP PPP --------------------

for (i in regions) {
  message("===== 2D tracking: Region ", sprintf("%02d", i), " =====")

  ra9 <- terra::rast(paste0(out_dir,
                            "/tracking3D_Region", sprintf("%02d", i),
                            "_", npp.model,
                            "_NPP9_", speed, "mday_21days_bottomnpp.tif"))
  # PPP generation using terra2im + rpoispp (as before) [1](https://rdrr.io/github/skiptoniam/qrbp/man/terra2im.html)
  Z <- terra2im(ra9 / ppp_divisor)
  pp <- spatstat.random::rpoispp(Z)
  npp.pts <- cbind(pp$x, pp$y, 0)

  load(fname_roms_region(env_dir3, i))  # Rdat6h

  # Assign z values as depth of cells (as in your script)
  roms.coords <- cbind(as.vector(Rdat6h$x), as.vector(Rdat6h$y))
  x.range <- c(min(roms.coords[, 1]) - 1000, max(roms.coords[, 1]) + 1000)
  y.range <- c(min(roms.coords[, 2]) - 1000, max(roms.coords[, 2]) + 1000)

  empty.roms.ra <- terra::rast(ext = terra::ext(x.range[1], x.range[2], y.range[1], y.range[2]),
                               resolution = 2000)
  h.ra <- terra::setValues(empty.roms.ra, as.vector(Rdat6h$h))
  npp.pts.h <- terra::extract(h.ra, npp.pts[, 1:2, drop = FALSE])
  npp.pts[, 3] <- npp.pts.h[, 1]

  if (reduce_depth_layers) {
    Rdat6h$hh <- Rdat6h$hh[, , -layers_to_drop, drop = FALSE]
    Rdat6h$i_u <- Rdat6h$i_u[, , -layers_to_drop, , drop = FALSE]
    Rdat6h$i_v <- Rdat6h$i_v[, , -layers_to_drop, , drop = FALSE]
    Rdat6h$i_w <- Rdat6h$i_w[, , -layers_to_drop, , drop = FALSE]
  }

  start.time <- Sys.time()
  track.2D <- loopit_2D3D(
    pts_seeded = npp.pts,
    romsobject = Rdat6h,
    projected = TRUE,
    sedimentation = sedimentation_2d,
    speed = speed,
    particle_radius = particle_radius,
    domain = "2D",
    runtime = runtime_2d,
    looping_time = looping_time_2d,
    roms_slices = roms_slices_2d,
    trajectories = FALSE,
    flux_mode = flux_mode,
    presence_tmp_dir = presence_tmp_dir,
    presence_buckets = presence_buckets,
    presence_flush_n = presence_flush_n,
    presence_cleanup = presence_cleanup,
    time_steps_in_s = time_steps_in_s
  )
  message("2D tracking took: ", format(Sys.time() - start.time))

  out_2d_file <- paste0(out_dir,
                        "/tracking2D_Region", sprintf("%02d", i),
                        "_", npp.model,
                        "_NPP9_", speed, "mday_21days_",
                        "r", r_str, "_28days.Rdata")
  save(track.2D, npp.pts, file = out_2d_file)

  rm(Rdat6h, h.ra, empty.roms.ra, pp, Z)
  gc(FALSE)
}

# -------------------- STEP 4: 2D ANALYSIS OUTPUTS (SED + FLUX) --------------------

for (i in regions) {
  message("===== 2D analysis outputs: Region ", sprintf("%02d", i), " =====")

  in_2d_file <- paste0(out_dir,
                       "/tracking2D_Region", sprintf("%02d", i),
                       "_", npp.model,
                       "_NPP9_", speed, "mday_21days_",
                       "r", r_str, "_28days.Rdata")
  load(in_2d_file)

  npp_crop <- terra::rast(fname_npp_raster_region(env_dir3, i))
  load(fname_roms_region(env_dir3, i))

  map_knn_to_cell <- load_or_build_knn_map(
    cache_dir = cache_dir,
    region = i,
    speed = speed,
    npp_raster = npp_crop,
    roms_x = Rdat6h$x,
    roms_y = Rdat6h$y
  )

  sed_ra <- sedimentation_raster_from_pend(npp_raster = npp_crop, pend = track.2D$pend_sed)
  sed_out <- paste0(out_dir,
                    "/tracking2D_Region", sprintf("%02d", i),
                    "_", npp.model,
                    "_NPP9_", speed, "mday_21days_r", r_str, "_28days_sed.tif")
  terra::writeRaster(sed_ra, sed_out, overwrite = TRUE)

  flux_out <- paste0(out_dir,
                     "/tracking2D_Region", sprintf("%02d", i),
                     "_", npp.model,
                     "_NPP9_", speed, "mday_21days_r", r_str, "_28days_flux.tif")

  if (!is.null(track.2D$flux_counts)) {
    flux_ra <- flux_raster_from_knn_counts(npp_raster = npp_crop,
                                           flux_counts = track.2D$flux_counts,
                                           map_knn_to_cell = map_knn_to_cell)
    terra::writeRaster(flux_ra, flux_out, overwrite = TRUE)
  } else if (!is.null(track.2D$idx_list_2D)) {
    # Legacy fallback
    n_knn <- length(Rdat6h$x)
    flux_counts <- legacy_flux_counts_from_idx_list(track.2D$idx_list_2D, n_knn = n_knn, verbose = TRUE)
    flux_ra <- flux_raster_from_knn_counts(npp_raster = npp_crop,
                                           flux_counts = flux_counts,
                                           map_knn_to_cell = map_knn_to_cell)
    terra::writeRaster(flux_ra, flux_out, overwrite = TRUE)
  } else {
    message("No flux information found for Region ", sprintf("%02d", i))
  }

  rm(Rdat6h, track.2D, npp_crop, map_knn_to_cell)
  gc(FALSE)
}
# -------------------- STEP 5: CIRCUMPOLAR COMPILATION --------------------

uv.max <- terra::rast(paste0(env_dir3, "ocean_his_bottom_uv_max.tif"))

flux_list <- vector("list", length(regions))
sed_list  <- vector("list", length(regions))

for (k in seq_along(regions)) {
  i <- regions[k]
  flux_file <- paste0(out_dir,
                      "/tracking2D_Region", sprintf("%02d", i),
                      "_", npp.model,
                      "_NPP9_", speed, "mday_21days_r", r_str, "_28days_flux.tif")
  sed_file <- paste0(out_dir,
                     "/tracking2D_Region", sprintf("%02d", i),
                     "_", npp.model, 
                     "_NPP9_", speed, "mday_21days_r", r_str, "_28days_sed.tif")
  
  flux_list[[k]] <- if (file.exists(flux_file)) terra::rast(flux_file) else NULL
  sed_list[[k]]  <- if (file.exists(sed_file))  terra::rast(sed_file)  else NULL
}

flux_aligned <- lapply(flux_list, function(r) {
  resample(r, uv.max, method = "near")  # or "bilinear" if truly continuous
})
flux_circ <- do.call(mosaic, c(flux_aligned, list(fun = "max")))

sed_aligned <- lapply(sed_list, function(r) {
  resample(r, uv.max, method = "near")  # or "bilinear" if truly continuous
})
sed_circ <- do.call(mosaic, c(sed_aligned, list(fun = "max")))

terra::writeRaster(flux_circ,
                   file.path(out_dir, paste0("tracking2D_",npp.model,"_NPP9_", speed, "mday_21days_r", r_str, "_28days_flux_circumpolar.tif")),
                   overwrite = TRUE)
terra::writeRaster(sed_circ,
                   file.path(out_dir, paste0("tracking2D_",npp.model,"_NPP9_", speed, "mday_21days_r", r_str, "_28days_sed_circumpolar.tif")),
                   overwrite = TRUE)

message("Pipeline complete.")
}