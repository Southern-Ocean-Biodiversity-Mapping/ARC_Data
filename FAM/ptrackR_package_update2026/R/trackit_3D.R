# trackit_3D.R (updated with xy+z indexing + NA policy)

trackit_3D <- function(pts,
                       romsobject,
                       w_sink = 100,
                       time = 50,
                       romsparams = NULL,
                       loop_trackit = FALSE,
                       time_steps_in_s = 1800,
                       projected = TRUE,
                       store_trajectory = FALSE,
                       store_indices = FALSE,
                       index_mode_3d = c("xy+z", "knn3d"),
                       na_action = c("stop", "zero_velocity", "error"),
                       hh_mat = NULL) {
  
  index_mode_3d <- match.arg(index_mode_3d)
  na_action <- match.arg(na_action)
  
  # KD setup
  if (loop_trackit) {
    kdxy <- romsobject$kdxy
    kdtree <- romsobject$kdtree
  } else {
    # Build only 2D tree for xy+z unless legacy knn3d requested
    build_3d <- identical(index_mode_3d, "knn3d")
    sknn <- with(romsobject, setup_knn(x, y, hh, build_3d = build_3d))
    kdxy <- sknn$kdxy
    kdtree <- sknn$kdtree
  }
  
  # ROMS slice vars
  if (!is.null(romsparams)) {
    i_u <- romsparams$i_u
    i_v <- romsparams$i_v
    i_w <- romsparams$i_w
    h <- romsparams$h
    hh <- romsobject$hh
  } else {
    i_u <- romsobject$i_u[,,,1]
    i_v <- romsobject$i_v[,,,1]
    i_w <- romsobject$i_w[,,,1]
    h <- romsobject$h
    hh <- romsobject$hh
  }
  
  # Prepare hh matrix for fast (idx2,k) lookups
  if (index_mode_3d == "xy+z") {
    if (is.null(hh_mat)) {
      hh_mat <- hh_to_matrix(hh)
    }
    nxy <- nrow(hh_mat)
  }
  
  ntime <- as.integer(time * 24 * 2)
  w_sink_ms <- -w_sink/(60*60*24)
  
  n_particles <- nrow(pts)
  plast <- pts
  pnow <- pts
  
  if (store_trajectory) {
    ptrack <- array(NA_real_, dim = c(n_particles, 3L, ntime))
  } else {
    ptrack <- NULL
  }
  
  if (store_indices) {
    indices_2D <- vector("list", ntime)
    indices_3D <- vector("list", ntime)
  } else {
    indices_2D <- NULL
    indices_3D <- NULL
  }
  
  stopindex <- integer(n_particles)
  stop_pos <- matrix(NA_real_, nrow = n_particles, ncol = 3L)
  active <- rep(TRUE, n_particles)
  stop_reason <- integer(n_particles)   # 0 = not stopped in this call, otherwise 1..4
  reason_labels <- c("bottom", "bad_velocity", "bad_bathymetry", "bad_position")
  
  for (itime in seq_len(ntime)) {
    
    if (!any(active)) break
    
    a_idx <- which(active)
    plast_a <- plast[a_idx, , drop = FALSE]
    
    # 2D horizontal cell
    idx2 <- kdxy$query(plast_a[,1:2, drop = FALSE], k = 1, eps = 0, radius = 0)$nn.idx
    
    # Determine 3D cell index
    if (index_mode_3d == "knn3d") {
      if (is.null(kdtree)) stop("kdtree is NULL but index_mode_3d='knn3d'")
      idx3 <- kdtree$query(plast_a, k = 1, eps = 0, radius = 0)$nn.idx
    } else {
      # vertical layer selection
      k <- nearest_layer_index(plast_a[,3], idx2, hh_mat)
      idx3 <- idx2k_to_idx3(idx2, k, nxy)
    }
    
    if (store_indices) {
      indices_2D[[itime]] <- idx2
      indices_3D[[itime]] <- idx3
    }
    
    thisu <- i_u[idx3]
    thisv <- i_v[idx3]
    thisw <- i_w[idx3]
    thish <- h[idx2]
    
    # Handle non-finite velocities/bathymetry
    bad_vel <- !is.finite(thisu) | !is.finite(thisv) | !is.finite(thisw)
    bad_h   <- !is.finite(thish)
    
    if (na_action == "error" && (any(bad_vel) || any(bad_h))) {
      stop("Non-finite velocity or bathymetry encountered in 3D tracking.")
    } else if (na_action == "zero_velocity") {
      thisu[bad_vel] <- 0
      thisv[bad_vel] <- 0
      thisw[bad_vel] <- 0
    }
    
    # Update
    pnow_a <- plast_a
    
    if (!projected) {
      pnow_a[,1] <- plast_a[,1] + (thisu * time_steps_in_s) / (1.852 * 60 * 1000 * cos(plast_a[,2] * pi/180))
      pnow_a[,2] <- plast_a[,2] + (thisv * time_steps_in_s) / (1.852 * 60 * 1000)
    } else {
      pnow_a[,1] <- plast_a[,1] + (thisu * time_steps_in_s)
      pnow_a[,2] <- plast_a[,2] + (thisv * time_steps_in_s)
    }
    
    pnow_a[,3] <- pmin(0, plast_a[,3]) + ((thisw + w_sink_ms) * time_steps_in_s)
    
    # Stop condition: bottom hit OR invalid bathymetry/velocity propagation
    stopped_now <- (pnow_a[,3] <= -thish)
    
    # If na_action == "stop", stop particles that hit non-finite inputs or produce non-finite outputs
    if (na_action == "stop") {
      stopped_now <- stopped_now | bad_vel | bad_h | !is.finite(pnow_a[,3]) | !is.finite(pnow_a[,1]) | !is.finite(pnow_a[,2])
    }
    
    # Commit positions
    pnow[a_idx,] <- pnow_a
    plast[a_idx,] <- pnow_a
    
    if (store_trajectory) {
      ptrack[a_idx,,itime] <- pnow_a
    }
    
    # # Record stops
    # if (any(stopped_now, na.rm = TRUE)) {
    #   stop_global <- a_idx[which(stopped_now)]
    #   stopindex[stop_global] <- itime
    #   stop_pos[stop_global,] <- pnow[stop_global,]
    #   active[stop_global] <- FALSE
    # }

    # Identify non-finite positions (diagnostic)
    bad_pos <- !is.finite(pnow_a[,1]) | !is.finite(pnow_a[,2]) | !is.finite(pnow_a[,3])
    
    # Stop-set from bottom-hit
    stopped_bottom <- (pnow_a[,3] <= -thish)
    
    # Full stop condition (unchanged behaviour)
    stopped_now <- stopped_bottom
    if (na_action == "stop") {
      stopped_now <- stopped_now | bad_vel | bad_h | bad_pos
    }
    
    # Record stops + reasons (diagnostic only; does not change motion)
    if (any(stopped_now, na.rm = TRUE)) {
      stop_local <- which(stopped_now)
      stop_global <- a_idx[stop_local]
      
      stopindex[stop_global] <- itime
      stop_pos[stop_global, ] <- pnow[stop_global, ]
      
      # Reason codes with priority:
      # 3 = bad bathymetry, 2 = bad velocity, 4 = bad position, 1 = bottom
      reason_local <- integer(length(stop_local))
      
      # Apply priority conditions
      bad_h_local   <- bad_h[stop_local]
      bad_vel_local <- bad_vel[stop_local]
      bad_pos_local <- bad_pos[stop_local]
      bottom_local  <- stopped_bottom[stop_local]
      
      reason_local[bottom_local] <- 1L
      reason_local[bad_pos_local] <- 4L
      reason_local[bad_vel_local] <- 2L
      reason_local[bad_h_local] <- 3L
      
      stop_reason[stop_global] <- reason_local
      active[stop_global] <- FALSE
    }
    
  }
  
  if (store_trajectory && exists("itime", inherits = FALSE)) {
    ptrack <- ptrack[, , seq_len(itime), drop = FALSE]
  }
  
  list(
    p_end = plast,
    stopindex = stopindex,
    stop_pos = stop_pos,
    stop_reason = stop_reason,
    stop_reason_labels = reason_labels,
    indices = indices_3D,
    indices_2D = indices_2D,
    ptrack = ptrack,
    hh_mat = hh_mat
  )
}



#' # trackit_3D.R
#' 
#' #' Track particles through a 3D ROMS field (minimal-memory by default).
#' #'
#' #' Particles are advected by ROMS u/v/w and additionally sink at constant w_sink.
#' #' The default behaviour is optimised for large particle counts:
#' #' - no full trajectory storage unless requested
#' #' - no per-timestep index storage unless requested
#' #'
#' #' @param pts Matrix with columns (x, y, z).
#' #' @param romsobject ROMS object containing x, y, hh (3D), h (2D), and velocity arrays.
#' #' @param w_sink Sinking rate (m/day).
#' #' @param time Duration in days for this call.
#' #' @param romsparams Optional list providing sliced i_u/i_v/i_w and h/extent.
#' #' @param loop_trackit If TRUE, re-use kdtrees stored on romsobject.
#' #' @param time_steps_in_s Timestep in seconds.
#' #' @param projected If FALSE, interpret x/y as lon/lat and convert displacement; otherwise x/y are metres.
#' #' @param store_trajectory If TRUE, store full (x,y,z) trajectory array for this call.
#' #' @param store_indices If TRUE, store per-timestep 3D and 2D cell indices (heavy).
#' #'
#' #' @return List with:
#' #' - p_end: final positions for all input particles (n x 3)
#' #' - stopindex: timestep index when particle hit bottom (within this call), or 0
#' #' - stop_pos: bottom-hit positions (n x 3; NA for not stopped)
#' #' - indices (optional): list of 3D knn indices per timestep
#' #' - indices_2D (optional): list of 2D knn indices per timestep
#' #' - ptrack (optional): trajectory array (n x 3 x ntime)
#' trackit_3D <- function(pts,
#'                        romsobject,
#'                        w_sink = 100,
#'                        time = 50,
#'                        romsparams = NULL,
#'                        loop_trackit = FALSE,
#'                        time_steps_in_s = 1800,
#'                        projected = TRUE,
#'                        store_trajectory = FALSE,
#'                        store_indices = FALSE) {
#'   
#'   # --- KD-trees ------------------------------------------------------------
#'   if (loop_trackit) {
#'     kdtree <- romsobject$kdtree
#'     kdxy <- romsobject$kdxy
#'   } else {
#'     sknn <- with(romsobject, setup_knn(x, y, hh))
#'     kdtree <- sknn$kdtree
#'     kdxy <- sknn$kdxy
#'   }
#'   
#'   # --- ROMS slice parameters ----------------------------------------------
#'   if (!is.null(romsparams)) {
#'     i_u <- romsparams$i_u
#'     i_v <- romsparams$i_v
#'     i_w <- romsparams$i_w
#'     h <- romsparams$h
#'   } else {
#'     # Keep legacy behaviour for standalone calls:
#'     i_u <- romsobject$i_u[,,,1]
#'     i_v <- romsobject$i_v[,,,1]
#'     i_w <- romsobject$i_w[,,,1]
#'     h <- romsobject$h
#'   }
#'   
#'   # --- Time setup ----------------------------------------------------------
#'   ntime <- as.integer(time * 24 * 2)  # half-hour steps
#'   w_sink_ms <- -w_sink / (60 * 60 * 24)  # convert m/day to m/s
#'   
#'   n_particles <- nrow(pts)
#'   pnow <- pts
#'   plast <- pts
#'   
#'   # Optional storage (avoid unless requested)
#'   if (store_trajectory) {
#'     ptrack <- array(NA_real_, dim = c(n_particles, 3L, ntime))
#'   } else {
#'     ptrack <- NULL
#'   }
#'   
#'   if (store_indices) {
#'     indices <- vector("list", ntime)
#'     indices_2D <- vector("list", ntime)
#'   } else {
#'     indices <- NULL
#'     indices_2D <- NULL
#'   }
#'   
#'   stopindex <- integer(n_particles)                 # 0 means not stopped
#'   stop_pos <- matrix(NA_real_, nrow = n_particles, ncol = 3L)
#'   
#'   # Track only active particles to reduce work as particles stop
#'   active <- rep(TRUE, n_particles)
#'   
#'   for (itime in seq_len(ntime)) {
#'     
#'     if (!any(active)) break
#'     
#'     a_idx <- which(active)
#'     plast_a <- plast[a_idx, , drop = FALSE]
#'     
#'     # 3D nearest neighbour (current ROMS cell)
#'     dmap <- kdtree$query(plast_a, k = 1, eps = 0, radius = 0)
#'     idx3 <- dmap$nn.idx
#'     
#'     # 2D nearest neighbour (for bathymetry stop condition)
#'     idx2 <- kdxy$query(plast_a[, 1:2, drop = FALSE], k = 1, eps = 0, radius = 0)$nn.idx
#'     
#'     if (store_indices) {
#'       indices[[itime]] <- idx3
#'       indices_2D[[itime]] <- idx2
#'     }
#'     
#'     thisu <- i_u[idx3]
#'     thisv <- i_v[idx3]
#'     thisw <- i_w[idx3]
#'     
#'     pnow_a <- plast_a
#'     
#'     if (!projected) {
#'       pnow_a[, 1] <- plast_a[, 1] + (thisu * time_steps_in_s) / (1.852 * 60 * 1000 * cos(plast_a[, 2] * pi/180))
#'       pnow_a[, 2] <- plast_a[, 2] + (thisv * time_steps_in_s) / (1.852 * 60 * 1000)
#'     } else {
#'       pnow_a[, 1] <- plast_a[, 1] + (thisu * time_steps_in_s)
#'       pnow_a[, 2] <- plast_a[, 2] + (thisv * time_steps_in_s)
#'     }
#'     
#'     # Vertical update (sink + ROMS w)
#'     pnow_a[, 3] <- pmin(0, plast_a[, 3]) + ((thisw + w_sink_ms) * time_steps_in_s)
#'     
#'     # Stop when reaching seabed depth at the 2D position
#'     stopped_now <- pnow_a[, 3] <= -h[idx2]
#'     
#'     # Commit positions
#'     pnow[a_idx, ] <- pnow_a
#'     plast[a_idx, ] <- pnow_a
#'     
#'     if (store_trajectory) {
#'       ptrack[a_idx, , itime] <- pnow_a
#'     }
#'     
#'     # Record stop info for newly stopped particles
#'     if (any(stopped_now)) {
#'       stop_global <- a_idx[stopped_now]
#'       stopindex[stop_global] <- itime
#'       stop_pos[stop_global, ] <- pnow[stop_global, ]
#'       active[stop_global] <- FALSE
#'     }
#'   }
#'   
#'   if (store_trajectory && exists("itime", inherits = FALSE)) {
#'     ptrack <- ptrack[, , seq_len(itime), drop = FALSE]
#'   }
#'   
#'   list(
#'     p_end = plast,
#'     stopindex = stopindex,
#'     stop_pos = stop_pos,
#'     indices = indices,
#'     indices_2D = indices_2D,
#'     ptrack = ptrack
#'   )
#' }