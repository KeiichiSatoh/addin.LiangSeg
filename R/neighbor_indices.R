#' Identify Daily-Life Neighbor Indices in a 3D Residential Array
#'
#' For each agent in a 3D residential array, returns the linear indices of
#' all residents who share the same 2D ground-level neighborhood, regardless
#' of which floor (z) they live on. Returns a flat 2D matrix that is easy to
#' use in downstream analysis.
#'
#' @details
#' This function models the assumption that residents in a multi-story
#' setting (e.g., a tower apartment building, or a stack of communities)
#' interact through ground-level activities such as shopping or walking,
#' rather than within their own floor. As a result, an agent at
#' \code{(r, c, z)} considers as potential neighbors:
#' \itemize{
#'   \item the eight surrounding cells on their own floor, and
#'   \item the corresponding \code{3 x 3} block on every other floor.
#' }
#'
#' Each agent's own cell is excluded from their neighbor list. Cells outside
#' the 2D boundary are represented as \code{NA}.
#'
#' Linear indices follow R's standard column-major order over the original
#' 3D array, so that \code{ar[result[i, ]]} retrieves the values of agent
#' \code{i}'s neighbors directly.
#'
#' This function is self-contained and does not depend on any external
#' helper functions.
#'
#' @param ar A 3-dimensional array where \code{ar[r, c, z]} represents the
#'   resident at row \code{r}, column \code{c}, floor \code{z}. The first
#'   two dimensions must each be at least 2.
#'
#' @return An integer matrix of dimensions
#'   \code{c(prod(dim(ar)), 8 + 9 * (dim(ar)[3] - 1))}.
#'   Row \code{i} contains the linear indices of all neighbors of the
#'   agent at linear position \code{i} in \code{ar}. Cells outside the
#'   2D boundary appear as \code{NA}.
#'
#'   Column names indicate the floor and relative position of each
#'   neighbor, e.g. \code{"floor1_top"}, \code{"floor2_self"}.
#'
#' @examples
#' # 3 x 3 footprint, 2 floors
#' ar <- array(1:18, dim = c(3, 3, 2))
#'
#' nbrs <- neighbor_indices(ar)
#' dim(nbrs)  # 18 rows (agents) x 17 columns (neighbors per agent)
#'
#' # Linear indices of the neighbors of the agent at ar[2, 2, 1]
#' # (linear index = 5):
#' nbrs[5, ]
#'
#' # Their actual values:
#' ar[nbrs[5, ]]
#'
#' @export
neighbor_indices <- function(ar){
  stopifnot(
    "'ar' must be a 3-dimensional array." =
      is.array(ar) && length(dim(ar)) == 3,
    "The first two dimensions of 'ar' must be at least 2." =
      dim(ar)[1] >= 2 && dim(ar)[2] >= 2
  )
  
  ar_dim <- dim(ar)
  nr <- ar_dim[1]
  nc <- ar_dim[2]
  nz <- ar_dim[3]
  N  <- nr * nc * nz
  
  # Linear-index array: idx_ar[r, c, z] = the linear index of ar[r, c, z]
  idx_ar <- array(seq_len(N), dim = ar_dim)
  
  # Position labels matching the 3rd-dimension order built below
  pos_labels <- c("self", "U", "UR", "R", "BR",
                  "B", "BL", "L", "UL")
  
  # --- Inline helper: 8-neighbor stack for one 2D slice ---------------
  # For a matrix `mat`, returns an [nr, nc, 9] array whose 3rd dimension
  # holds (in order): self, top, topright, right, bottomright,
  # bottom, bottomleft, left, topleft.
  # Out-of-boundary positions are filled with NA.
  build_neighbor_stack <- function(mat){
    right       <- cbind(mat[, 2:nc, drop = FALSE], NA)
    left        <- cbind(NA, mat[, 1:(nc - 1), drop = FALSE])
    top         <- rbind(NA, mat[1:(nr - 1), , drop = FALSE])
    bottom      <- rbind(mat[2:nr, , drop = FALSE], NA)
    topright    <- rbind(NA, right[1:(nr - 1), , drop = FALSE])
    topleft     <- rbind(NA, left[1:(nr - 1), , drop = FALSE])
    bottomright <- rbind(right[2:nr, , drop = FALSE], NA)
    bottomleft  <- rbind(left[2:nr, , drop = FALSE], NA)
    
    stack <- array(NA_integer_, dim = c(nr, nc, 9))
    stack[, , 1] <- mat
    stack[, , 2] <- top
    stack[, , 3] <- topright
    stack[, , 4] <- right
    stack[, , 5] <- bottomright
    stack[, , 6] <- bottom
    stack[, , 7] <- bottomleft
    stack[, , 8] <- left
    stack[, , 9] <- topleft
    stack
  }
  # --------------------------------------------------------------------
  
  # For each floor z2, build the index-neighbor stack once. This gives,
  # for every (r, c), the linear indices of the 3x3 block on floor z2.
  per_floor <- lapply(seq_len(nz), function(z2){
    build_neighbor_stack(idx_ar[, , z2])
  })
  
  # Build column blocks. For observer floor z, the column block for
  # floor z2 contains:
  #   - all 9 positions if z2 != z (other floors: keep self)
  #   - 8 positions    if z2 == z (own floor: drop self)
  # Each block is an [nr*nc, k] matrix (rows in column-major (r, c) order).
  build_block <- function(z, z2){
    block_3d <- per_floor[[z2]]
    keep <- if (z == z2) 2:9 else 1:9
    mat <- matrix(block_3d[, , keep],
                  nrow = nr * nc,
                  ncol = length(keep))
    colnames(mat) <- paste0("F", z2, "_", pos_labels[keep])
    mat
  }
  
  # For each observer floor z, concatenate column blocks across all z2.
  # Then stack the z-blocks vertically to form the final [N, ncol] matrix.
  per_observer_floor <- lapply(seq_len(nz), function(z){
    do.call(cbind, lapply(seq_len(nz), function(z2) build_block(z, z2)))
  })
  
  result <- do.call(rbind, per_observer_floor)
  
  # Sanity check on dimensions: should be [N, 8 + 9*(nz - 1)]
  stopifnot(
    nrow(result) == N,
    ncol(result) == 8 + 9 * (nz - 1)
  )
  
  result
}