# presence_flux_disk.R

#' Initialise a disk-backed writer for exact presence-only (P1) flux counting.
#'
#' This writer stores (particle_id, cell_id) pairs in bucketed binary files.
#' Bucketing ensures each file is small enough to deduplicate in memory later.
#'
#' The on-disk format for each bucket is a sequence of records:
#'   [int32 n][int32 particle_id x n][int32 cell_id x n]
#'
#' @param tmp_dir Directory used to write bucket files.
#' @param nbuckets Number of bucket files to partition pairs across.
#' @param max_cell_id Maximum possible cell index (used for key packing during aggregation).
#' @param flush_n Number of pairs to buffer before flushing to disk.
#' @param prefix Filename prefix used for bucket files.
#'
#' @return An environment containing writer state.
presence_writer_init <- function(tmp_dir,
                                 nbuckets = 256L,
                                 max_cell_id,
                                 flush_n = 2e6,
                                 prefix = "presence_pairs") {
  if (!dir.exists(tmp_dir)) dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  nbuckets <- as.integer(nbuckets)
  max_cell_id <- as.integer(max_cell_id)
  
  if (nbuckets < 1L) stop("nbuckets must be >= 1")
  if (max_cell_id < 1L) stop("max_cell_id must be >= 1")
  
  # MULT is used later to pack (pid, cid) into a numeric key: key = pid * MULT + cid
  MULT <- as.numeric(max_cell_id) + 1
  
  files <- file.path(tmp_dir, sprintf("%s_bucket%04d.bin", prefix, seq_len(nbuckets)))
  
  st <- new.env(parent = emptyenv())
  st$tmp_dir <- tmp_dir
  st$nbuckets <- nbuckets
  st$max_cell_id <- max_cell_id
  st$MULT <- MULT
  st$flush_n <- as.integer(flush_n)
  st$prefix <- prefix
  st$files <- files
  
  # Buffers
  st$buf_pid <- integer(0)
  st$buf_cid <- integer(0)
  
  st
}

#' Add (particle_id, cell_id) pairs to the presence writer buffer.
#'
#' @param writer Presence writer environment from presence_writer_init().
#' @param particle_id Integer vector of particle IDs.
#' @param cell_id Integer vector of cell IDs (same length as particle_id).
#'
#' @return Invisibly returns writer.
presence_writer_add <- function(writer, particle_id, cell_id) {
  if (length(particle_id) != length(cell_id)) stop("particle_id and cell_id must have same length")
  if (!is.integer(particle_id)) particle_id <- as.integer(particle_id)
  if (!is.integer(cell_id)) cell_id <- as.integer(cell_id)
  
  writer$buf_pid <- c(writer$buf_pid, particle_id)
  writer$buf_cid <- c(writer$buf_cid, cell_id)
  
  if (length(writer$buf_pid) >= writer$flush_n) {
    presence_writer_flush(writer)
  }
  invisible(writer)
}

#' Flush buffered pairs to bucket files on disk.
#'
#' Uses a fast integer hash to bucket pairs. Data are written in compact binary
#' records to reduce overhead.
#'
#' @param writer Presence writer environment.
#' @return Invisibly returns writer.
presence_writer_flush <- function(writer) {
  n <- length(writer$buf_pid)
  if (n == 0L) return(invisible(writer))
  
  pid <- writer$buf_pid
  cid <- writer$buf_cid
  
  # Fast deterministic hash into buckets (uses double arithmetic safely for these magnitudes).
  # The constants are chosen for mixing; exact values are not important as long as stable.
  b <- as.integer((abs(pid * 1315423911 + cid * 2654435761) %% writer$nbuckets) + 1L)
  
  # Group by bucket using ordering to avoid creating huge split() lists.
  o <- order(b)
  b_sorted <- b[o]
  pid_sorted <- pid[o]
  cid_sorted <- cid[o]
  
  # Identify run boundaries per bucket
  runs <- which(c(TRUE, b_sorted[-1L] != b_sorted[-length(b_sorted)]))
  runs_end <- c(runs[-1L] - 1L, length(b_sorted))
  
  for (k in seq_along(runs)) {
    i1 <- runs[k]
    i2 <- runs_end[k]
    bucket <- b_sorted[i1]
    
    pid_k <- pid_sorted[i1:i2]
    cid_k <- cid_sorted[i1:i2]
    nk <- length(pid_k)
    
    con <- file(writer$files[bucket], open = "ab")
    on.exit(close(con), add = TRUE)
    
    writeBin(as.integer(nk), con, size = 4L, endian = "little")
    writeBin(as.integer(pid_k), con, size = 4L, endian = "little")
    writeBin(as.integer(cid_k), con, size = 4L, endian = "little")
    
    close(con)
  }
  
  # Clear buffers
  writer$buf_pid <- integer(0)
  writer$buf_cid <- integer(0)
  
  invisible(writer)
}

#' Aggregate exact presence-only (P1) counts per cell from bucket files.
#'
#' For each bucket file:
#'   1) read all records
#'   2) pack pairs into numeric keys: key = pid * MULT + cid
#'   3) unique(keys) to deduplicate within bucket
#'   4) tabulate cell IDs and add to global counts
#'
#' @param writer Presence writer environment.
#' @param cleanup If TRUE, remove bucket files after aggregation.
#' @param verbose If TRUE, print progress messages.
#'
#' @return Integer vector of length max_cell_id giving presence-only counts per cell.
presence_writer_finalize <- function(writer, cleanup = TRUE, verbose = TRUE) {
  # Ensure any buffered pairs are written
  presence_writer_flush(writer)
  
  max_cell <- writer$max_cell_id
  MULT <- writer$MULT
  counts <- integer(max_cell)
  
  for (f in writer$files) {
    if (!file.exists(f) || file.info(f)$size == 0) next
    
    if (verbose) message("Aggregating presence bucket: ", basename(f))
    
    con <- file(f, open = "rb")
    on.exit(close(con), add = TRUE)
    
    keys_chunks <- list()
    nch <- 0L
    
    repeat {
      n <- readBin(con, integer(), n = 1L, size = 4L, endian = "little")
      if (length(n) == 0L) break
      if (n <= 0L) next
      
      pid <- readBin(con, integer(), n = n, size = 4L, endian = "little")
      cid <- readBin(con, integer(), n = n, size = 4L, endian = "little")
      
      # Pack to numeric keys; safe under typical magnitudes (<= 2^53 exact integer in double).
      key <- as.numeric(pid) * MULT + as.numeric(cid)
      nch <- nch + 1L
      keys_chunks[[nch]] <- key
    }
    
    close(con)
    
    # Deduplicate within this bucket
    if (nch > 0L) {
      keys <- unique(unlist(keys_chunks, use.names = FALSE))
      cell_id <- as.integer(keys %% MULT)
      
      # tabulate requires positive integers; cell_id should be in 1..max_cell
      counts <- counts + tabulate(cell_id, nbins = max_cell)
    }
    
    rm(keys_chunks)
    gc(FALSE)
  }
  
  if (cleanup) {
    suppressWarnings(file.remove(writer$files))
  }
  
  counts
}