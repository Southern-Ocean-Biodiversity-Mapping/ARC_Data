# indexing_3D_xy_z.R
# Utilities for 3D indexing via 2D kNN + vertical layer selection ("xy+z").

#' Flatten hh (nx, ny, nz) into a matrix (nxy, nz) consistent with R's column-major order.
#' This lets us index column depths quickly by hh_mat[idx2, k].
#'
#' @param hh 3D numeric array (nx, ny, nz)
#' @return numeric matrix (nxy, nz)
hh_to_matrix <- function(hh) {
  d <- dim(hh)
  if (length(d) != 3L) stop("hh must be a 3D array (nx, ny, nz)")
  nxy <- d[1] * d[2]
  nz  <- d[3]
  matrix(hh, nrow = nxy, ncol = nz)
}

#' Select nearest vertical layer index for each particle based on depth.
#'
#' Uses an O(nz) scan per particle but vectorised over particles; nz is small (e.g. 31).
#'
#' @param z numeric vector of particle depths (negative downward, usually <= 0)
#' @param idx2 integer vector of horizontal cell indices (1..nxy)
#' @param hh_mat numeric matrix (nxy, nz) from hh_to_matrix()
#'
#' @return integer vector k (1..nz) nearest layer index per particle
nearest_layer_index <- function(z, idx2, hh_mat) {
  nz <- ncol(hh_mat)
  n  <- length(z)
  best_k <- integer(n)
  best_d <- rep(Inf, n)
  
  for (k in seq_len(nz)) {
    d <- abs(z - hh_mat[idx2, k])
    upd <- d < best_d
    if (any(upd)) {
      best_d[upd] <- d[upd]
      best_k[upd] <- k
    }
  }
  best_k
}

#' Convert (idx2, k) into linear 3D array indices for arrays shaped (nx, ny, nz).
#'
#' @param idx2 integer vector (1..nxy)
#' @param k integer vector (1..nz)
#' @param nxy integer scalar = nx*ny
#' @return integer vector idx3
idx2k_to_idx3 <- function(idx2, k, nxy) {
  idx2 + (k - 1L) * nxy
}