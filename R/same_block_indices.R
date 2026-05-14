#' Identify Same-Block Resident Indices in a 3D Residential Array
#'
#' For each resident in a 3D residential array, returns the linear indices
#' of all other residents who belong to the same block. Block membership is
#' defined by the value of the cell on the first floor (\code{city[, , 1]}):
#' residents at \code{(r, c, z)} share a block iff they share the same
#' \code{city[r, c, 1]} value, regardless of floor.
#'
#' @details
#' This function models a community structure where block membership is a
#' ground-level property (e.g., neighborhood district, administrative
#' zone). All residents at the same horizontal position \code{(r, c)} on
#' any floor inherit the block ID stamped on the first floor.
#'
#' The function:
#' \enumerate{
#'   \item Replicates the first-floor block map across all floors.
#'   \item Groups all linear indices by block ID in a single pass
#'     (\code{O(N)} via \code{\link[base]{split}}).
#'   \item For each resident, returns the indices of all other members of
#'     their block.
#' }
#'
#' Cells with \code{NA} in \code{city[, , 1]} represent uninhabitable
#' positions; the corresponding rows in the output are filled entirely
#' with \code{NA}.
#'
#' Each resident's own index is excluded from their row.
#'
#' @param city A 3-dimensional array representing the city. Each element
#'   is a block ID (any atomic type that \code{\link[base]{split}}
#'   accepts as a factor). \code{NA} values mark uninhabitable cells.
#'   Block IDs are read from the first floor (\code{city[, , 1]}) and
#'   propagated to all floors.
#'
#' @return An integer matrix of dimensions \code{c(prod(dim(city)), K - 1)},
#'   where \code{K} is the size of the largest block. Row \code{i}
#'   contains the linear indices (in column-major order over
#'   \code{city}) of all other residents in the same block as resident
#'   \code{i}. Rows are padded with \code{NA} when a block has fewer
#'   than \code{K} members, or when the resident's cell is \code{NA}.
#'
#' @examples
#' # 3 x 3 footprint, 2 floors, 3 blocks labeled A/B/C on the 1st floor
#' city <- array(NA, dim = c(3, 3, 2))
#' city[, , 1] <- matrix(c("A", "A", "B",
#'                         "A", "B", "B",
#'                         "C", "C", "C"), nrow = 3, byrow = TRUE)
#' city[, , 2] <- city[, , 1]  # block IDs propagate (function does this internally too)
#'
#' nbrs <- same_block_indices(city)
#' dim(nbrs)
#'
#' # Indices of all same-block residents of city[1, 1, 1]:
#' nbrs[1, ]
#'
#' @seealso \code{\link{neighbor_indices}} for spatial (8-neighbor) lookups.
#'
#' @export
same_block_indices <- function(city){
  stopifnot(
    "'city' must be a 3-dimensional array." =
      is.array(city) && length(dim(city)) == 3
  )
  
  city_dim <- dim(city)
  N <- prod(city_dim)
  
  # Propagate first-floor block IDs to all floors.
  # Using rep() explicitly (rather than relying on array() recycling) makes
  # the intent unambiguous.
  city_1F <- city[, , 1]
  city_block <- array(rep(city_1F, times = city_dim[3]), dim = city_dim)
  
  # Group all linear indices by block ID in one pass.
  # split() automatically drops NA keys, so NA cells are excluded here.
  all_idx <- seq_len(N)
  block_members <- split(all_idx, city_block)
  
  # Largest block size determines K; output has K - 1 columns (self excluded).
  K <- max(lengths(block_members), 0L)
  ncol_out <- max(K - 1L, 0L)
  
  # Pre-fill output with NA. Rows corresponding to NA cells will stay NA.
  out <- matrix(NA_integer_, nrow = N, ncol = ncol_out)
  
  # For each block, fill in the rows of its members with the indices of
  # the other members (self excluded), padded with NA on the right.
  for (members in block_members) {
    if (length(members) <= 1L) next  # singletons have no same-block neighbors
    for (i in members) {
      others <- members[members != i]
      out[i, seq_along(others)] <- others
    }
  }
  
  out
}
