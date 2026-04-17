# loopit_2D3D_traj_ROMSreduced.R (complete, with 3D indexing integration)

loopit_2D3D <- function(pts_seeded,
                        romsobject,
                        roms_slices = 1,
                        start_slice = 1,
                        domain = "2D",
                        trajectories = FALSE,
                        speed,
                        runtime = 10,
                        looping_time = 0.25,
                        sedimentation = FALSE,
                        sedimentation_mode = "legacy",
                        particle_radius = 0.00016,
                        time_steps_in_s = 1800,
                        uphill_restricted = NULL,
                        sed_at_max_speed = FALSE,
                        mean_move = FALSE,
                        projected = TRUE,
                        detailed_trajectories = FALSE,
                        ROMSreduced = FALSE,
                        cell_counting = FALSE,
                        flux_mode = c("entry", "total", "presence"),
                        presence_tmp_dir = "/pvol3TB/tmp",
                        presence_buckets = 256L,
                        presence_flush_n = 2e6,
                        presence_cleanup = TRUE,
                        index_mode_3d = c("xy+z", "knn3d"),
                        na_action_3d = c("stop", "zero_velocity", "error")) {
  
  flux_mode <- match.arg(flux_mode)
  index_mode_3d <- match.arg(index_mode_3d)
  na_action_3d <- match.arg(na_action_3d)
  
  pts <- pts_seeded
  id_vec <- seq_len(nrow(pts_seeded))
  
  # ---- kNN setup ----------------------------------------------------------
  kd.time <- Sys.time()
  
  # Backward compatible call to setup_knn (in case you didn't replace setup_knn.R)
  build_3d_tree <- (domain == "3D" && index_mode_3d == "knn3d")
  if ("build_3d" %in% names(formals(setup_knn))) {
    sknn <- with(romsobject, setup_knn(x, y, hh, build_3d = build_3d_tree))
  } else {
    sknn <- with(romsobject, setup_knn(x, y, hh))
    if (!build_3d_tree) {
      # If legacy setup_knn always builds kdtree, we just won't use it in xy+z mode.
    }
  }
  
  romsobject$kdtree <- sknn$kdtree
  romsobject$kdxy <- sknn$kdxy
  message("... kNN setup successfully")
  message("took: ", format(Sys.time() - kd.time))
  
  # Precompute hh_mat for xy+z (reuse across slices)
  hh_mat <- NULL
  if (domain == "3D" && index_mode_3d == "xy+z") {
    hh_mat <- hh_to_matrix(romsobject$hh)
  }
  
  # ---- ROMS params container ---------------------------------------------
  romsparams <- list(
    h = romsobject$h,
    i_u = romsobject$i_u,
    i_v = romsobject$i_v,
    i_w = romsobject$i_w,
    roms_ext = c(min(romsobject$x), max(romsobject$x), min(romsobject$y), max(romsobject$y))
  )
  
  # Optional max speed arrays (kept as-is; only used in 2D sedimentation)
  if (sed_at_max_speed && roms_slices == 4) {
    nrow_ <- dim(romsobject$i_u)[1]
    ncol_ <- dim(romsobject$i_u)[2]
    nlayer <- dim(romsobject$i_u)[3]
    i_u_max <- array(NA_real_, dim = c(nrow_, ncol_, nlayer))
    i_v_max <- array(NA_real_, dim = c(nrow_, ncol_, nlayer))
    for (irow in seq_len(nrow_)) {
      for (icol in seq_len(ncol_)) {
        for (ilayer in seq_len(nlayer)) {
          i_u_max[irow, icol, ilayer] <- max(abs(romsobject$i_u[irow, icol, ilayer, ]))
          i_v_max[irow, icol, ilayer] <- max(abs(romsobject$i_v[irow, icol, ilayer, ]))
        }
      }
    }
    romsparams$i_u_max <- i_u_max
    romsparams$i_v_max <- i_v_max
  }
  
  # Sedimentation params (2D)
  sedimentationparams <- if (domain == "2D") {
    buildparams(-speed / (60 * 60 * 24), time_step_in_s = time_steps_in_s, r = particle_radius)
  } else NULL
  
  # ---- Slice schedule -----------------------------------------------------
  curr_vector <- rep(seq_len(roms_slices), runtime)
  sliced_vector <- curr_vector[c(start_slice:length(curr_vector), seq_len(start_slice - 1))]
  total_loops <- length(sliced_vector)
  
  # ---- Output containers --------------------------------------------------
  xyz_end <- matrix(NA_real_, nrow = nrow(pts_seeded), ncol = 3L)
  colnames(xyz_end) <- c("x","y","z")
  stop_reason_global <- integer(nrow(pts_seeded))  # 0 = never stopped within runtime; else 1..4
  stop_reason_counts <- integer(4)                
  
  x_list <- list(); y_list <- list()
  x_sed_list <- list()
  y_sed_list <- list()
  depth_list <- if (domain == "3D") list() else NULL
  
  # Flux state (2D)
  flux_counts <- NULL
  prev_cell_idx <- NULL
  presence_writer <- NULL
  
  if (domain == "2D" && flux_mode == "presence") {
    max_cell_id <- prod(dim(romsobject$x))
    presence_writer <- presence_writer_init(
      tmp_dir = presence_tmp_dir,
      nbuckets = presence_buckets,
      max_cell_id = max_cell_id,
      flush_n = presence_flush_n,
      prefix = "presence_pairs"
    )
  }
  
  # Optional heavy objects
  if (trajectories) {
    idx_list_2D <- list()
    id_list <- list()
  } else {
    idx_list_2D <- NULL
    id_list <- NULL
  }
  
  # ---- Main loop ----------------------------------------------------------
  for (irun in seq_len(total_loops)) {
    
    if (irun == 1) message("starting # of particles: ", nrow(pts))
    message(irun, ".loop")
    s.time <- Sys.time()
    
    # Slice velocities
    if (roms_slices > 1) {
      romsparams$i_u <- romsobject$i_u[,,,sliced_vector[irun]]
      romsparams$i_v <- romsobject$i_v[,,,sliced_vector[irun]]
      romsparams$i_w <- romsobject$i_w[,,,sliced_vector[irun]]
    }
    
    if (trajectories) id_list[[irun]] <- id_vec
    
    if (domain == "3D") {
      obj <- trackit_3D(
        pts = pts,
        romsobject = romsobject,
        w_sink = speed,
        time = looping_time,
        romsparams = romsparams,
        loop_trackit = TRUE,
        time_steps_in_s = time_steps_in_s,
        projected = projected,
        index_mode_3d = index_mode_3d,
        na_action = na_action_3d,
        hh_mat = hh_mat,
        store_trajectory = FALSE,
        store_indices = FALSE
      )
      hh_mat <- obj$hh_mat
    } else {
      obj <- trackit_2D(
        pts = pts,
        romsobject = romsobject,
        w_sink = speed,
        time = looping_time,
        romsparams = romsparams,
        loop_trackit = TRUE,
        time_steps_in_s = time_steps_in_s,
        projected = projected,
        sedimentationparams = sedimentationparams,
        sedimentation = sedimentation,
        sedimentation_mode = sedimentation_mode,
        particle_radius = particle_radius,
        uphill_restricted = uphill_restricted,
        sed_at_max_speed = sed_at_max_speed,
        mean_move = mean_move,
        store_trajectory = (trajectories && detailed_trajectories),
        store_indices = trajectories,
        flux_mode = flux_mode,
        flux_counts = flux_counts,
        prev_cell_idx = prev_cell_idx,
        ids = id_vec,
        presence_writer = presence_writer
      )
      flux_counts <- obj$flux_counts
      prev_cell_idx <- obj$prev_cell_idx
    }
    
    if (domain == "2D" && !is.null(obj$stop_reason)) {
      sed_local <- which(obj$stop_reason == 2L)
      if (length(sed_local) > 0L) {
        stop_pos_sed <- obj$stop_pos[sed_local, , drop = FALSE]
        x_sed_list[[irun]] <- stop_pos_sed[, 1]
        y_sed_list[[irun]] <- stop_pos_sed[, 2]
      } else {
        x_sed_list[[irun]] <- numeric(0)
        y_sed_list[[irun]] <- numeric(0)
      }
    }
    
    # Update final positions in original indexing
    xyz_end[id_vec, ] <- obj$p_end
    
    # Store stop positions
    stopped_local <- which(obj$stopindex != 0L)
    if (length(stopped_local) > 0L) {
      stop_pos <- obj$stop_pos[stopped_local, , drop = FALSE]
      x_list[[irun]] <- stop_pos[,1]
      y_list[[irun]] <- stop_pos[,2]
      if (domain == "3D") depth_list[[irun]] <- stop_pos[,3]
    } else {
      x_list[[irun]] <- numeric(0)
      y_list[[irun]] <- numeric(0)
      if (domain == "3D") depth_list[[irun]] <- numeric(0)
    }
    
    if (trajectories && domain == "2D") {
      idx_list_2D[[irun]] <- obj$indices_2D
    }
    
    if (domain == "3D" && length(stopped_local) > 0L && !is.null(obj$stop_reason)) {
      # Map reasons from local particle indices to original global particle indices
      stopped_global_ids <- id_vec[stopped_local]
      reason_codes <- obj$stop_reason[stopped_local]
      
      stop_reason_global[stopped_global_ids] <- reason_codes
      stop_reason_counts <- stop_reason_counts + tabulate(reason_codes, nbins = 4)
    }
    
    # Reduce to still-active particles
    still_active <- which(obj$stopindex == 0L)
    id_vec <- id_vec[still_active]
    if (domain == "2D" && flux_mode %in% c("entry", "presence") && !is.null(prev_cell_idx)) {
      prev_cell_idx <- prev_cell_idx[still_active]
    }
    
    pts <- if (length(still_active) > 0L) obj$p_end[still_active, , drop = FALSE] else matrix(numeric(0), ncol = 3)
    
    message(format(Sys.time() - s.time))
    message(nrow(pts), " particles floating")
    
    if (nrow(pts) == 0L) break
  }
  
  # Finalise presence-only flux
  if (domain == "2D" && flux_mode == "presence") {
    flux_counts <- presence_writer_finalize(presence_writer, cleanup = presence_cleanup, verbose = TRUE)
  }
  
  pend <- if (domain == "3D") {
    cbind(unlist(x_list), unlist(y_list), unlist(depth_list))
  } else {
    cbind(unlist(x_list), unlist(y_list))
  }
  
  out <- list(
    pts = pts,
    pend = pend,
    xyz_end = xyz_end,
    flux_counts = flux_counts
  )
  
  if (trajectories) {
    out$idx_list_2D <- idx_list_2D
    out$id_list <- id_list
  }
  
  if (domain == "3D") {
    out$stop_reason_global <- stop_reason_global
    out$stop_reason_counts <- stop_reason_counts
    out$stop_reason_labels <- if (!is.null(obj$stop_reason_labels)) {
      obj$stop_reason_labels
    } else {
      c("bottom", "bad_velocity", "bad_bathymetry", "bad_position")
    }
  }
  if (domain == "2D") {
    out$pend_sed <- cbind(unlist(x_sed_list), unlist(y_sed_list))
  }

  out
}



#' # loopit_2D3D_traj_ROMSreduced.R
#' 
#' #' Loop particle tracking over consecutive ROMS time slices (2D or 3D).
#' #'
#' #' This wrapper runs trackit_2D or trackit_3D repeatedly over ROMS slices.
#' #' It supports:
#' #' - minimal-memory production runs (default): no trajectories, no index lists
#' #' - optional diagnostic outputs: trajectories and index lists
#' #' - streaming flux counting modes for 2D:
#' #'     "entry" (default), "total", "presence" (exact P1; disk-backed)
#' #'
#' #' @param pts_seeded Matrix of particles with 3 columns (x, y, z).
#' #' @param romsobject ROMS object with x, y, h, hh and velocity arrays.
#' #' @param roms_slices Number of ROMS time frames available in romsobject (4th dimension).
#' #' @param start_slice ROMS slice index to start with.
#' #' @param domain "2D" or "3D".
#' #' @param trajectories If TRUE, store index lists (heavy). Default FALSE.
#' #' @param speed Sinking speed (m/day).
#' #' @param runtime Number of times to loop through the ROMS slice sequence.
#' #' @param looping_time Duration (days) per wrapper loop step. Default 0.25 (6h).
#' #' @param sedimentation If TRUE (2D only), apply density-dependent settling.
#' #' @param particle_radius Particle radius for sedimentation model.
#' #' @param time_steps_in_s Time step in seconds used in tracking functions.
#' #' @param uphill_restricted Restrict uphill movement by depth difference (m).
#' #' @param sed_at_max_speed Use max velocity fields for sedimentation (requires 4 slices precomputed).
#' #' @param mean_move Use midpoint correction when uphill restriction triggers (slower).
#' #' @param projected TRUE if x/y are metres, FALSE if lon/lat.
#' #' @param detailed_trajectories If TRUE and trajectories=TRUE, store detailed x/y trajectories (very heavy).
#' #' @param ROMSreduced If TRUE, allow reducing depth layers during 3D runs (unchanged in this refactor).
#' #'
#' #' @param flux_mode For 2D flux counts: "entry" (default), "total", or "presence".
#' #' @param presence_tmp_dir Scratch directory used for disk-backed presence-only (P1) counting.
#' #' @param presence_buckets Number of bucket files for presence-only counting.
#' #' @param presence_flush_n Buffer size for presence-only pair writing.
#' #' @param presence_cleanup If TRUE, remove temporary presence bucket files after aggregation.
#' #'
#' #' @return A list with run outputs. Minimal outputs by default; heavy outputs only when requested.
#' loopit_2D3D <- function(pts_seeded,
#'                         romsobject,
#'                         roms_slices = 1,
#'                         start_slice = 1,
#'                         domain = "2D",
#'                         trajectories = FALSE,
#'                         speed,
#'                         runtime = 10,
#'                         looping_time = 0.25,
#'                         sedimentation = FALSE,
#'                         particle_radius = 0.00016,
#'                         time_steps_in_s = 1800,
#'                         uphill_restricted = NULL,
#'                         sed_at_max_speed = FALSE,
#'                         mean_move = FALSE,
#'                         projected = TRUE,
#'                         detailed_trajectories = FALSE,
#'                         ROMSreduced = FALSE,
#'                         cell_counting = FALSE,
#'                         flux_mode = c("entry", "total", "presence"),
#'                         presence_tmp_dir = "/pvol3TB/tmp",
#'                         presence_buckets = 256L,
#'                         presence_flush_n = 2e6,
#'                         presence_cleanup = TRUE) {
#'   
#'   flux_mode <- match.arg(flux_mode)
#'   
#'   pts <- pts_seeded
#' 
#'   index_mode_3d <- if (exists("index_mode_3d")) index_mode_3d else "xy+z"
#'   na_action_3d  <- if (exists("na_action_3d"))  na_action_3d  else "stop"
#'   
#'   # --- Build kdtrees once --------------------------------------------------
#'   # kd.time <- Sys.time()
#'   # sknn <- with(romsobject, setup_knn(x, y, hh))
#'   # romsobject$kdtree <- sknn$kdtree
#'   # romsobject$kdxy <- sknn$kdxy
#'   # message("... kdtree setup successfully")
#'   # message("took: ", format(Sys.time() - kd.time))
#' 
#'   kd.time <- Sys.time()
#'   build_3d_tree <- (domain == "3D" && index_mode_3d == "knn3d")
#'   sknn <- with(romsobject, setup_knn(x, y, hh, build_3d = build_3d_tree))
#'   romsobject$kdtree <- sknn$kdtree
#'   romsobject$kdxy <- sknn$kdxy
#'   message("... kNN setup successfully")
#'   print(paste0("took this long: ", Sys.time() - kd.time))
#'   
#'   # Precompute hh_mat for xy+z mode once (reuse across slices)
#'   hh_mat <- NULL
#'   if (domain == "3D" && index_mode_3d == "xy+z") {
#'     hh_mat <- hh_to_matrix(romsobject$hh)
#'   }
#'   
#'   # --- Stable global particle IDs -----------------------------------------
#'   id_vec <- seq_len(nrow(pts_seeded))
#'   
#'   # --- Minimal end-position store (matrix avoids data.frame copying) -------
#'   xyz_end <- matrix(NA_real_, nrow = nrow(pts_seeded), ncol = 3L)
#'   colnames(xyz_end) <- c("x", "y", "z")
#'   
#'   # --- Optional heavy stores ----------------------------------------------
#'   if (trajectories) {
#'     idx_list <- list()
#'     idx_list_2D <- list()
#'     id_list <- list()
#'   } else {
#'     idx_list <- NULL
#'     idx_list_2D <- NULL
#'     id_list <- NULL
#'   }
#'   
#'   if (detailed_trajectories) {
#'     ptrack_x_list <- list()
#'     ptrack_y_list <- list()
#'   } else {
#'     ptrack_x_list <- NULL
#'     ptrack_y_list <- NULL
#'   }
#'   
#'   # --- Sedimentation parameters (2D) --------------------------------------
#'   if (domain == "2D") {
#'     sedimentationparams <- buildparams(-speed / (60 * 60 * 24), time_step_in_s = time_steps_in_s, r = particle_radius)
#'   } else {
#'     sedimentationparams <- NULL
#'   }
#'   
#'   # --- ROMS parameter list -------------------------------------------------
#'   romsparams <- list()
#'   romsparams$h <- romsobject$h
#'   romsparams$i_u <- romsobject$i_u
#'   romsparams$i_v <- romsobject$i_v
#'   romsparams$i_w <- romsobject$i_w
#'   
#'   # Preserve existing behaviour: wrapper supplies explicit extent used for stop conditions.
#'   romsparams$roms_ext <- c(min(romsobject$x), max(romsobject$x), min(romsobject$y), max(romsobject$y))
#'   
#'   # Optional: max-speed sedimentation fields (unchanged; kept as in your version)
#'   if (sed_at_max_speed && roms_slices == 4) {
#'     nrow_ <- dim(romsobject$i_u)[1]
#'     ncol_ <- dim(romsobject$i_u)[2]
#'     nlayer <- dim(romsobject$i_u)[3]
#'     i_u_max <- array(NA_real_, dim = c(nrow_, ncol_, nlayer))
#'     i_v_max <- array(NA_real_, dim = c(nrow_, ncol_, nlayer))
#'     for (irow in seq_len(nrow_)) {
#'       for (icol in seq_len(ncol_)) {
#'         for (ilayer in seq_len(nlayer)) {
#'           i_u_max[irow, icol, ilayer] <- max(abs(romsobject$i_u[irow, icol, ilayer, ]))
#'           i_v_max[irow, icol, ilayer] <- max(abs(romsobject$i_v[irow, icol, ilayer, ]))
#'         }
#'       }
#'     }
#'     romsparams$i_u_max <- i_u_max
#'     romsparams$i_v_max <- i_v_max
#'   }
#'   
#'   # --- Slice schedule ------------------------------------------------------
#'   curr_vector <- rep(seq_len(roms_slices), runtime)
#'   sliced_vector <- curr_vector[c(start_slice:length(curr_vector), seq_len(start_slice - 1))]
#'   total_loops <- length(sliced_vector)
#'   
#'   # --- Flux accumulators ---------------------------------------------------
#'   flux_counts <- NULL
#'   prev_cell_idx <- NULL
#'   
#'   # Presence-only writer (P1) uses disk-backed buckets
#'   presence_writer <- NULL
#'   if (domain == "2D" && flux_mode == "presence") {
#'     # Maximum 2D cell index is the number of horizontal grid points used by kdxy.
#'     max_cell_id <- prod(dim(romsobject$x))
#'     presence_writer <- presence_writer_init(
#'       tmp_dir = presence_tmp_dir,
#'       nbuckets = presence_buckets,
#'       max_cell_id = max_cell_id,
#'       flush_n = presence_flush_n,
#'       prefix = "presence_pairs"
#'     )
#'   }
#'   
#'   # --- Stop position collectors -------------------------------------------
#'   x_list <- list()
#'   y_list <- list()
#'   depth_list <- if (domain == "3D") list() else NULL
#'   
#'   # --- Loop over ROMS slices ----------------------------------------------
#'   for (irun in seq_len(total_loops)) {
#'     
#'     if (irun == 1) message("starting # of particles: ", nrow(pts))
#'     message(irun, ".loop")
#'     s.time <- Sys.time()
#'     
#'     # Update ROMS slice velocities if multiple slices
#'     if (roms_slices > 1) {
#'       romsparams$i_u <- romsobject$i_u[,,,sliced_vector[irun]]
#'       romsparams$i_v <- romsobject$i_v[,,,sliced_vector[irun]]
#'       romsparams$i_w <- romsobject$i_w[,,,sliced_vector[irun]]
#'     }
#'     
#'     # Track IDs present in this wrapper loop if requested
#'     if (trajectories) id_list[[irun]] <- id_vec
#'     
#'     # Run tracking for currently active particles
#'     if (domain == "3D") {
#'       # obj <- trackit_3D(
#'       #   pts = pts,
#'       #   romsobject = romsobject,
#'       #   w_sink = speed,
#'       #   time = looping_time,
#'       #   romsparams = romsparams,
#'       #   loop_trackit = TRUE,
#'       #   time_steps_in_s = time_steps_in_s,
#'       #   projected = projected
#'       # )
#'       obj <- trackit_3D(
#'         pts = pts,
#'         romsobject = romsobject,
#'         w_sink = speed,
#'         time = looping_time,
#'         romsparams = romsparams,
#'         loop_trackit = TRUE,
#'         time_steps_in_s = time_steps_in_s,
#'         projected = projected,
#'         index_mode_3d = index_mode_3d,
#'         na_action = na_action_3d,
#'         hh_mat = hh_mat
#'       )
#'       hh_mat <- obj$hh_mat
#'     } else {
#'       obj <- trackit_2D(
#'         pts = pts,
#'         romsobject = romsobject,
#'         w_sink = speed,
#'         time = looping_time,
#'         romsparams = romsparams,
#'         loop_trackit = TRUE,
#'         time_steps_in_s = time_steps_in_s,
#'         projected = projected,
#'         sedimentationparams = sedimentationparams,
#'         sedimentation = sedimentation,
#'         particle_radius = particle_radius,
#'         uphill_restricted = uphill_restricted,
#'         sed_at_max_speed = sed_at_max_speed,
#'         mean_move = mean_move,
#'         store_trajectory = (trajectories && detailed_trajectories),
#'         store_indices = trajectories,
#'         flux_mode = flux_mode,
#'         flux_counts = flux_counts,
#'         prev_cell_idx = prev_cell_idx,
#'         ids = id_vec,
#'         presence_writer = presence_writer
#'       )
#'       
#'       # Carry streaming flux and entry-state across wrapper loops
#'       flux_counts <- obj$flux_counts
#'       prev_cell_idx <- obj$prev_cell_idx
#'     }
#'     
#'     # Final positions for all particles processed in this loop
#'     xyz_end[id_vec, ] <- obj$p_end
#'     
#'     # Record stop positions for particles that stopped in this loop step
#'     stopped_local <- which(obj$stopindex != 0L)
#'     if (length(stopped_local) > 0L) {
#'       stop_pos <- obj$stop_pos[stopped_local, , drop = FALSE]
#'       x_list[[irun]] <- stop_pos[, 1]
#'       y_list[[irun]] <- stop_pos[, 2]
#'       if (domain == "3D") depth_list[[irun]] <- stop_pos[, 3]
#'     } else {
#'       x_list[[irun]] <- numeric(0)
#'       y_list[[irun]] <- numeric(0)
#'       if (domain == "3D") depth_list[[irun]] <- numeric(0)
#'     }
#'     
#'     # Store indices if requested (2D)
#'     if (trajectories && domain == "2D") {
#'       idx_list_2D[[irun]] <- obj$indices_2D
#'     }
#'     
#'     # Store detailed trajectories if requested (2D)
#'     if (trajectories && detailed_trajectories && domain == "2D") {
#'       # ptrack is n_particles_active x 3 x ntime; store x/y series
#'       ptrack_x_list[[irun]] <- obj$ptrack[, 1, , drop = FALSE]
#'       ptrack_y_list[[irun]] <- obj$ptrack[, 2, , drop = FALSE]
#'     }
#'     
#'     # Remove stopped particles for next wrapper loop
#'     still_active <- which(obj$stopindex == 0L)
#'     id_vec <- id_vec[still_active]
#'     
#'     if (domain == "2D" && flux_mode == "entry" && !is.null(prev_cell_idx)) {
#'       prev_cell_idx <- prev_cell_idx[still_active]
#'     }
#'     
#'     # Prepare pts for next iteration
#'     if (length(still_active) > 0L) {
#'       pts <- obj$p_end[still_active, , drop = FALSE]
#'     } else {
#'       pts <- matrix(numeric(0), ncol = 3)
#'     }
#'     
#'     message(format(Sys.time() - s.time))
#'     message(nrow(pts), " particles floating")
#'     
#'     if (nrow(pts) == 0L) break
#'   }
#'   
#'   # Finalise presence-only (P1) if selected
#'   if (domain == "2D" && flux_mode == "presence") {
#'     flux_counts <- presence_writer_finalize(
#'       presence_writer,
#'       cleanup = presence_cleanup,
#'       verbose = TRUE
#'     )
#'   }
#'   
#'   # Combine stop positions across wrapper loops
#'   if (domain == "3D") {
#'     pend <- cbind(unlist(x_list), unlist(y_list), unlist(depth_list))
#'   } else {
#'     pend <- cbind(unlist(x_list), unlist(y_list))
#'   }
#'   
#'   message((nrow(pts_seeded) - nrow(pend)), " particle(s) still floating")
#'   
#'   # Build output
#'   out <- list(
#'     pts = pts,
#'     pend = pend,
#'     xyz_end = xyz_end,
#'     flux_counts = flux_counts
#'   )
#'   
#'   if (trajectories) {
#'     out$idx_list_2D <- idx_list_2D
#'     out$id_list <- id_list
#'   }
#'   if (trajectories && detailed_trajectories) {
#'     out$ptrack_x_list <- ptrack_x_list
#'     out$ptrack_y_list <- ptrack_y_list
#'   }
#'   
#'   out
#' }