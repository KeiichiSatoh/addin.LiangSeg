#' Create a 3D city array
#'
#' Generates a 3D array representing a city composed of zones and lots,
#' with randomly assigned building heights.
#'
#' @param zone_dim Integer vector of length 2: number of zones (rows, cols).
#' @param lot_dim  Integer vector of length 2: lot dimensions within each zone.
#' @param NA_zone_ind Integer vector of zone indices to set as NA. Default NULL.
#' @param NA_lot_ind  Integer vector of lot indices (within a zone) to set as NA. Default NULL.
#' @param max_height Positive integer: maximum building height.
#' @param height_prop Numeric vector of length `max_height` giving the
#'   proportion of buildings at each height. Will be normalized to sum to 1.
#'
#' @return A 3D array of dimension `c(zone_dim[1]*lot_dim[1], zone_dim[2]*lot_dim[2], max_height)`.
#'   The first two axes are spatial coordinates, the third axis is height level.
#'   Values are zone numbers; cells above a building's height or in NA zones/lots are NA.
#'
#' @note Requires the external function `make_group_labels(n, prop, random_adjustment)`.
create_city <- function(zone_dim = c(2, 2),
                        lot_dim = c(3, 3),
                        NA_zone_ind = NULL,
                        NA_lot_ind = NULL,
                        max_height = 2,
                        height_prop = rep(1, max_height) / max_height) {
  
  #------------------------
  # Argument validation
  #------------------------
  stopifnot(
    length(zone_dim) == 2, all(zone_dim >= 1), all(zone_dim == as.integer(zone_dim)),
    length(lot_dim)  == 2, all(lot_dim  >= 1), all(lot_dim  == as.integer(lot_dim)),
    length(max_height) == 1, max_height >= 1, max_height == as.integer(max_height),
    length(height_prop) == max_height
  )
  
  if (any(is.na(height_prop)) || any(height_prop < 0)) {
    stop("'height_prop' must be non-negative and contain no NA")
  }
  if (sum(height_prop) <= 0) {
    stop("'height_prop' must have a positive sum")
  }
  height_prop <- height_prop / sum(height_prop)
  
  if (!is.null(NA_lot_ind) &&
      !all(NA_lot_ind %in% seq_len(lot_dim[1] * lot_dim[2]))) {
    stop("'NA_lot_ind' contains invalid indices")
  }
  if (!is.null(NA_zone_ind) && !all(NA_zone_ind > 0)) {
    stop("'NA_zone_ind' must contain positive values")
  }
  
  #------------------------
  # Build the zone-id matrix (vectorized via kronecker)
  #------------------------
  zone_ids <- matrix(seq_len(prod(zone_dim)),
                     nrow = zone_dim[1], ncol = zone_dim[2])
  city_0 <- kronecker(zone_ids,
                      matrix(1, nrow = lot_dim[1], ncol = lot_dim[2]))
  
  # Apply lot-level NA mask
  if (!is.null(NA_lot_ind)) {
    lot_mask <- matrix(1, nrow = lot_dim[1], ncol = lot_dim[2])
    lot_mask[NA_lot_ind] <- NA
    city_0 <- city_0 * kronecker(matrix(1, zone_dim[1], zone_dim[2]), lot_mask)
  }
  
  # Apply zone-level NA mask
  if (!is.null(NA_zone_ind)) {
    city_0[city_0 %in% NA_zone_ind] <- NA
  }
  
  #------------------------
  # Expand to 3D (broadcast city_0 across the height axis)
  #------------------------
  city <- array(city_0, dim = c(dim(city_0), max_height))
  
  #------------------------
  # Assign random heights
  #------------------------
  n_cell <- sum(!is.na(city_0))
  height <- rABM::make_group_labels(n = n_cell,
                                    prop = height_prop,
                                    random_adjustment = TRUE)
  height <- sample(height)  # randomize the height assignment
  
  height_mat <- matrix(NA_integer_, nrow = nrow(city_0), ncol = ncol(city_0))
  height_mat[!is.na(city_0)] <- height
  
  # Mask cells above each building's height (vectorized)
  height_array <- array(rep(height_mat, times = max_height), dim = dim(city))
  level_array  <- array(rep(seq_len(max_height),
                            each = prod(dim(city_0))),
                        dim = dim(city))
  city[level_array > height_array] <- NA
  
  city
}
