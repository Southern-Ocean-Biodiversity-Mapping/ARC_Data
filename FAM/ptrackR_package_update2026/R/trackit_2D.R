# trackit_2D.R

#' Track particles through a ROMS field in 2D-space.
#'
#' This function advances particles horizontally using ROMS u/v components and
#' applies optional sedimentation rules. It can optionally accumulate flux counts
#' (cell visits) in a streaming manner without storing per-particle trajectories.
#'
#' Flux counting modes:
#' - "entry"   (default): count only when a particle enters a new cell.
#' - "total"             : count every timestep a particle occupies a cell.
#' - "presence"          : exact presence-only (P1): each particle counts at most once per cell
#'                         over the entire run; uses a disk-backed writer and is slower.
#'
#' @param pts Matrix with columns (x, y, z). z is not used for advection in 2D.
#' @param romsobject ROMS object with at least x, y, h and velocity fields.
#' @param w_sink Sinking rate (m/day). Used only for sedimentation parameterisation.
#' @param time Total duration in days for this call.
#' @param sedimentation TRUE/FALSE for density-dependent settling.
#' @param particle_radius Particle radius affecting sedimentation.
#' @param force_final_settling If TRUE, assign stopindex for any remaining particles at end.
#' @param romsparams ROMS parameters passed in from wrapper (recommended for performance).
#' @param sedimentationparams Precomputed sedimentation parameters (optional).
#' @param loop_trackit TRUE when called by loopit_2D3D to re-use kdtrees.
#' @param time_steps_in_s Timestep in seconds (default 1800 = 0.5h).
#' @param uphill_restricted If not NULL, restrict uphill movement beyond this depth difference.
#' @param sed_at_max_speed If TRUE, use precomputed max velocities when estimating sedimentation.
#' @param mean_move If TRUE, attempt midpoint correction for uphill restriction (slower).
#' @param projected If FALSE, interpret x/y as lon/lat and convert displacement; otherwise x/y are metres.
#'
#' @param store_trajectory If TRUE, store full (x,y,z) trajectory array. Expensive for many particles.
#' @param store_indices If TRUE, store list of 2D cell indices per timestep. Expensive for many particles.
#'
#' @param flux_mode One of "entry", "total", "presence". Default is "entry".
#' @param flux_counts Optional integer vector of existing counts to update (length = n_cells_2D).
#' @param prev_cell_idx Optional integer vector of previous cell indices for "entry" mode continuity.
#' @param ids Optional stable particle IDs (required for "presence").
#' @param presence_writer Disk-backed writer created by presence_writer_init() (required for "presence").
#'
#' @return A list containing:
#' - p_end: final positions for all input particles (n x 3)
#' - stopindex: timestep index (within this call) when particle stopped, or 0 if not stopped
#' - stop_pos: stop positions for particles that stopped (n x 3; NA for not stopped)
#' - flux_counts: updated counts (if flux_mode is "entry" or "total")
#' - prev_cell_idx: last cell index for remaining particles (for wrapper continuity)
#' - indices_2D: list of indices per timestep (only if store_indices=TRUE)
#' - ptrack: trajectory array (only if store_trajectory=TRUE)
trackit_2D <- function(pts,
                       romsobject,
                       w_sink = 100,
                       time = 50,
                       sedimentation = FALSE,
                       particle_radius = 0.00016,
                       sedimentation_mode = c("legacy", "fraction"),
                       force_final_settling = FALSE,
                       romsparams = NULL,
                       sedimentationparams = NULL,
                       loop_trackit = FALSE,
                       time_steps_in_s = 1800,
                       uphill_restricted = NULL,
                       sed_at_max_speed = FALSE,
                       mean_move = FALSE,
                       projected = TRUE,
                       store_trajectory = FALSE,
                       store_indices = FALSE,
                       flux_mode = c("entry", "total", "presence"),
                       flux_counts = NULL,
                       prev_cell_idx = NULL,
                       ids = NULL,
                       presence_writer = NULL) {
  
  flux_mode <- match.arg(flux_mode)
  sedimentation_mode <- match.arg(sedimentation_mode)
  
  # --- KD-trees ------------------------------------------------------------
  if (loop_trackit) {
    kdtree <- romsobject$kdtree
    kdxy <- romsobject$kdxy
  } else {
    sknn <- with(romsobject, setup_knn(x, y, hh))
    kdtree <- sknn$kdtree
    kdxy <- sknn$kdxy
  }
  
  # --- ROMS parameters -----------------------------------------------------
  if (!is.null(romsparams)) {
    i_u <- romsparams$i_u
    i_v <- romsparams$i_v
    h <- romsparams$h
    roms_ext <- romsparams$roms_ext
  } else {
    i_u <- romsobject$i_u
    i_v <- romsobject$i_v
    h <- romsobject$h
    
    # Preserve existing behaviour; wrapper normally supplies roms_ext explicitly.
    roms_ext <- c(min(romsobject$x), max(romsobject$x), min(romsobject$y), max(romsobject$y))
  }
  
  # --- Sedimentation parameters -------------------------------------------
  # w_sink is m/day; buildparams expects speed in m/s in your usage.
  w_sink_ms <- -w_sink / (60 * 60 * 24)
  
  if (!is.null(sedimentationparams)) {
    params <- sedimentationparams
  } else {
    params <- buildparams(w_sink_ms, time_step_in_s = time_steps_in_s, r = particle_radius)
  }
  
  # --- Time steps ----------------------------------------------------------
  ntime <- as.integer(time * 24 * 2)  # half-hour steps per day
  
  n_particles <- nrow(pts)
  pnow <- pts
  plast <- pts
  
  # Stable particle IDs (required for presence-only mode)
  if (is.null(ids)) ids <- seq_len(n_particles)
  if (flux_mode == "presence" && is.null(presence_writer)) {
    stop("presence_writer must be provided when flux_mode='presence'")
  }
  
  # Trajectories and indices are optional; avoid allocation unless requested.
  if (store_trajectory) {
    ptrack <- array(NA_real_, dim = c(n_particles, 3L, ntime))
  } else {
    ptrack <- NULL
  }
  
  if (store_indices) {
    indices_2D <- vector("list", ntime)
  } else {
    indices_2D <- NULL
  }
  
  # Stop bookkeeping
  stopindex <- integer(n_particles)                    # 0 means active
  stop_pos <- matrix(NA_real_, nrow = n_particles, ncol = 3L)
  stop_reason <- integer(n_particles)  # 0=not stopped in this call; 1..4=reason
  
  # Flux bookkeeping
  # Runtime note:
  #   "entry" and "total" are streaming operations and are typically fast.
  #   "presence" writes to disk and deduplicates later, so it is substantially slower.
  if (flux_mode %in% c("entry", "total")) {
    if (is.null(flux_counts)) {
      # flux_counts length must equal the number of 2D cells; caller should supply this.
      # If not supplied, it can be created after the first kdxy query using max cell index.
      flux_counts <- NULL
    }
  }
  
  if (flux_mode %in% c("entry", "presence") && is.null(prev_cell_idx)) {
    # Initialise previous cell positions from starting locations.
    prev_cell_idx <- kdxy$query(plast[, 1:2, drop = FALSE], k = 1, eps = 0, radius = 0)$nn.idx
  }

  # Active particle mask; shrink computations as particles stop within this call.
  active <- rep(TRUE, n_particles)
  
  for (itime in seq_len(ntime)) {
    
    if (!any(active)) break
    
    a_idx <- which(active)
    plast_a <- plast[a_idx, , drop = FALSE]
    
    # Map current positions to 2D cells (used for velocities)
    two_dim_pos <- kdxy$query(plast_a[, 1:2, drop = FALSE], k = 1, eps = 0, radius = 0)
    idx_for_roms <- two_dim_pos$nn.idx
    
    # Velocities at current cells
    thisu <- i_u[idx_for_roms]
    thisv <- i_v[idx_for_roms]
    thish <- h[idx_for_roms]
    
    # Advance positions
    pnow_a <- plast_a
    if (!projected) {
      pnow_a[, 1] <- plast_a[, 1] + (thisu * time_steps_in_s) / (1.852 * 60 * 1000 * cos(plast_a[, 2] * pi/180))
      pnow_a[, 2] <- plast_a[, 2] + (thisv * time_steps_in_s) / (1.852 * 60 * 1000)
    } else {
      pnow_a[, 1] <- plast_a[, 1] + (thisu * time_steps_in_s)
      pnow_a[, 2] <- plast_a[, 2] + (thisv * time_steps_in_s)
    }
    
    # Determine destination cell for the proposed new positions
    tdp_idx <- kdxy$query(pnow_a[, 1:2, drop = FALSE], k = 1, eps = 0, radius = 0)$nn.idx
    
    # Uphill restriction (optional)
    if (!is.null(uphill_restricted)) {
      uphill <- h[tdp_idx] < thish - uphill_restricted
      if (any(uphill)) {
        if (!mean_move) {
          pnow_a[uphill, ] <- plast_a[uphill, ]
          tdp_idx[uphill] <- idx_for_roms[uphill]
        } else {
          # Midpoint correction; repeated once as in original logic.
          pnow_a[uphill, ] <- (pnow_a[uphill, ] + plast_a[uphill, ]) / 2
          tdp_test <- kdxy$query(pnow_a[, 1:2, drop = FALSE], k = 1, eps = 0, radius = 0)$nn.idx
          still_uphill <- h[tdp_test] < thish - uphill_restricted
          
          if (any(still_uphill)) {
            pnow_a[still_uphill, ] <- (pnow_a[still_uphill, ] + plast_a[still_uphill, ]) / 2
            tdp_test2 <- kdxy$query(pnow_a[, 1:2, drop = FALSE], k = 1, eps = 0, radius = 0)$nn.idx
            stillstill_uphill <- h[tdp_test2] < thish - uphill_restricted
            if (any(stillstill_uphill)) {
              pnow_a[stillstill_uphill, ] <- plast_a[stillstill_uphill, ]
              tdp_test2[stillstill_uphill] <- idx_for_roms[stillstill_uphill]
            }
            tdp_idx <- tdp_test2
          } else {
            tdp_idx <- tdp_test
          }
        }
      }
    }
    
    # Stopping condition 1: outside domain extent
    stopped_out <- (pnow_a[, 1] < roms_ext[1] | pnow_a[, 1] > roms_ext[2] |
                      pnow_a[, 2] < roms_ext[3] | pnow_a[, 2] > roms_ext[4])
    
    # Stopping condition 2: sedimentation (optional)
    stopped_sed <- rep(FALSE, length(a_idx))
    if (sedimentation) {
      # Particle density per 2D cell among active particles at this step
      dens_all <- tabulate(tdp_idx)
      
      # Cells that contain >= 1 active particle
      cells_used <- which(dens_all > 0L)
      dens_used <- dens_all[cells_used]
      
      # Current speed squared per used cell
      if (!sed_at_max_speed) {
        vel_sq <- i_u[cells_used]^2 + i_v[cells_used]^2
      } else {
        # romsparams must supply i_u_max/i_v_max if sed_at_max_speed=TRUE
        vel_sq <- romsparams$i_u_max[cells_used]^2 + romsparams$i_v_max[cells_used]^2
      }
      
      # U_div term (McCave & Swift style)
      U_div <- 1 - (vel_sq / params$Ucsq)
      U_div[U_div < 0] <- 0
      
      # Number of particles to settle per cell at this step (legacy sign conventions retained)
      n_drop <- params$SedFunct(U_div, dens_used)
      
      dens_for_pt <- dens_all[tdp_idx]
      drop_for_pt <- n_drop[match(tdp_idx, cells_used)]
      
      if (sedimentation_mode == "legacy") {
        # Legacy approach from original code:
        # p = -(n_drop / dens)
        # This can exceed 1, which forces settling; no clamping is applied to match legacy behaviour.
        p_settle <- -(drop_for_pt / dens_for_pt)
      } else {
        # Alternative (non-legacy) fraction approach (kept available):
        # p = max(0, min(1, n_drop/dens)) after sign correction
        p_settle <- drop_for_pt / dens_for_pt
        p_settle <- pmax(0, pmin(1, p_settle))
      }
      
      p_settle[!is.finite(p_settle)] <- 0
      stopped_sed <- (runif(length(p_settle)) <= p_settle)
    }
    
    # Combine stopping and guard against non-finite velocity/positions.
    # If u/v are non-finite, the update can produce NaNs which then propagate into logical tests.
    bad_vel <- !is.finite(thisu) | !is.finite(thisv)
    bad_pos <- !is.finite(pnow_a[, 1]) | !is.finite(pnow_a[, 2])
    
    stopped_now <- stopped_out | stopped_sed | bad_vel | bad_pos
    stopped_now[is.na(stopped_now)] <- TRUE
    
    # Flux accounting (streaming)
    if (flux_mode == "total") {
      # Ensure flux_counts exists
      if (is.null(flux_counts)) flux_counts <- integer(max(tdp_idx))
      if (length(flux_counts) < max(tdp_idx)) length(flux_counts) <- max(tdp_idx)
      flux_counts <- flux_counts + tabulate(tdp_idx, nbins = length(flux_counts))
    } else if (flux_mode == "entry") {
      # Compare destination cell against previous cell indices for active particles.
      prev_a <- prev_cell_idx[a_idx]
      changed <- (tdp_idx != prev_a)
      if (any(changed)) {
        idx_changed <- tdp_idx[changed]
        if (is.null(flux_counts)) flux_counts <- integer(max(idx_changed))
        if (length(flux_counts) < max(idx_changed)) length(flux_counts) <- max(idx_changed)
        flux_counts <- flux_counts + tabulate(idx_changed, nbins = length(flux_counts))
      }
      # Update previous cell indices for active particles (even if unchanged)
      prev_cell_idx[a_idx] <- tdp_idx
    } else if (flux_mode == "presence") {
      # Write only entry events to reduce I/O volume; still exact for presence-only.
      prev_a <- prev_cell_idx[a_idx]
      changed <- (tdp_idx != prev_a)
      if (any(changed)) {
        presence_writer_add(presence_writer,
                            particle_id = ids[a_idx][changed],
                            cell_id = tdp_idx[changed])
      }
      prev_cell_idx[a_idx] <- tdp_idx
    }
    
    # Store indices if requested
    if (store_indices) {
      indices_2D[[itime]] <- tdp_idx
    }
    
    # Commit new positions to global arrays for active particles
    pnow[a_idx, ] <- pnow_a
    plast[a_idx, ] <- pnow_a
    
    # Store trajectory if requested
    if (store_trajectory) {
      ptrack[a_idx, , itime] <- pnow_a
    }
    
    # Record stop index and stop positions for particles stopping now (within this call)
    stop_local <- which(stopped_now)
    if (length(stop_local) > 0L) {
      stop_global <- a_idx[stop_local]
      
      stopindex[stop_global] <- itime
      stop_pos[stop_global, ] <- pnow[stop_global, ]
      active[stop_global] <- FALSE
      
      # Stop reason codes (diagnostic only):
      # 1 = out of domain
      # 2 = sedimentation
      # 3 = bad velocity (u/v not finite)
      # 4 = bad position (x/y not finite)
      reason_local <- integer(length(stop_local))
      
      out_local <- stopped_out[stop_local]
      sed_local <- stopped_sed[stop_local]
      vel_local <- bad_vel[stop_local]
      pos_local <- bad_pos[stop_local]
      
      reason_local[out_local] <- 1L
      reason_local[sed_local] <- 2L
      reason_local[vel_local] <- 3L
      reason_local[pos_local] <- 4L
      
      stop_reason[stop_global] <- reason_local
    }
    
  }
  
  # Force stopindex for any remaining particles (optional)
  if (force_final_settling) {
    still_active <- which(stopindex == 0L)
    if (length(still_active) > 0L) {
      stopindex[still_active] <- max(1L, itime)
      stop_pos[still_active, ] <- plast[still_active, ]
    }
  }
  
  # Trim trajectory array if stored (to actual itime reached)
  if (store_trajectory && exists("itime", inherits = FALSE)) {
    ptrack <- ptrack[, , seq_len(itime), drop = FALSE]
  }
  
  list(
    p_end = plast,
    stopindex = stopindex,
    stop_pos = stop_pos,
    stop_reason = stop_reason,
    stop_reason_labels = c("out_of_domain", "sedimentation", "bad_velocity", "bad_position"),
    flux_counts = flux_counts,
    prev_cell_idx = prev_cell_idx,
    indices_2D = indices_2D,
    ptrack = ptrack
  )
}