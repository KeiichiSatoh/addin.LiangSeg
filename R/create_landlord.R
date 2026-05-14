#' Create landlords and assign property ownership in a city
#'
#' Generates a set of landlords with specified aversion parameters and
#' assigns ownership of each cell in a city array. Ownership can be
#' assigned either uniformly at random or based on geographic distance
#' from each landlord's origin point, with support for "dispersed"
#' landlords who own property uniformly across the city.
#'
#' @param city A 3D array representing a city, typically produced by
#'   [create_city()]. The first two dimensions are spatial coordinates and
#'   the third is the height axis. NA cells are treated as non-housing and
#'   receive no owner.
#' @param n_landlord A positive integer giving the number of landlords to
#'   create. Default is `5`.
#' @param ethnic_aversion A non-negative numeric vector of length
#'   `n_landlord` representing each landlord's aversion to ethnic
#'   diversity (or similar ethnicity-based preference). The interpretation
#'   is application-specific; higher values typically indicate stronger
#'   aversion. Default is `c(1, 1, 1, 0, 0)`.
#' @param SES_aversion A non-negative numeric vector of length
#'   `n_landlord` representing each landlord's aversion based on
#'   socioeconomic status (SES). The interpretation is application-
#'   specific; higher values typically indicate stronger aversion.
#'   Default is `c(0, 0, 1, 1, 1)`.
#' @param dispersed_landlord An integer vector of length `n_landlord`
#'   containing 0 or 1. A value of `1` marks the landlord as "dispersed",
#'   meaning they own property uniformly across the city (rather than
#'   clustered near an origin point). Only relevant when
#'   `geo_difference = TRUE`.
#' @param geo_difference Logical. If `TRUE`, ownership is assigned based
#'   on geographic distance from each landlord's origin, with closer cells
#'   being more likely to be owned by that landlord. If `FALSE`, ownership
#'   is assigned uniformly at random. Default is `FALSE`.
#' @param geo_lambda A positive numeric scalar controlling the rate of
#'   distance decay in the ownership propensity. Higher values produce
#'   more spatially clustered ownership. Default is `1`. Only used when
#'   `geo_difference = TRUE`.
#' @param landlord_origin_vec An optional integer vector of length
#'   `n_landlord` giving the cell index (in column-major order over the
#'   first two dimensions of `city`) of each landlord's origin point.
#'   Must contain unique values within `1:prod(dim(city)[1:2])`. If `NULL`
#'   (the default), origins are sampled at random without replacement.
#'   Only used when `geo_difference = TRUE`.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{`landlord`}{A data frame with columns `ID`, `ethnic_aversion`,
#'       `SES_aversion`, and `dispersed`, one row per landlord.}
#'     \item{`ownership`}{A 3D integer array with the same dimensions as
#'       `city`. Each non-NA cell contains the ID of the landlord that
#'       owns it; cells that are NA in `city` are also NA here. Ownership
#'       is constant along the height axis (a building has a single owner
#'       across all of its levels).}
#'   }
#'
#' @details
#' When `geo_difference = TRUE`, ownership propensity for each landlord is
#' computed as `exp(-geo_lambda * d)`, where `d` is the Euclidean distance
#' from the landlord's origin to the target cell. For each cell, the owner
#' is sampled with probability proportional to these propensities across
#' landlords. Cells that coincide with a landlord's origin are assigned
#' deterministically to that landlord.
#'
#' Dispersed landlords (those with `dispersed_landlord == 1`) have their
#' propensity row replaced by its row mean, yielding a uniform propensity
#' across all cells while preserving the same scale as non-dispersed
#' landlords. This makes them compete on equal footing across the entire
#' city rather than clustering near an origin.
#'
#' @examples
#' \dontrun{
#' city <- create_city(zone_dim = c(2, 2), lot_dim = c(3, 3), max_height = 2)
#'
#' # Uniform random ownership
#' result1 <- create_landlord(city, n_landlord = 3,
#'                            ethnic_aversion = c(1, 0, 0),
#'                            SES_aversion    = c(0, 1, 0),
#'                            dispersed_landlord = c(0, 0, 0))
#'
#' # Geographically clustered ownership with two dispersed landlords
#' result2 <- create_landlord(city, n_landlord = 5,
#'                            ethnic_aversion = c(1, 1, 1, 0, 0),
#'                            SES_aversion    = c(0, 0, 1, 1, 1),
#'                            dispersed_landlord = c(0, 0, 0, 1, 1),
#'                            geo_difference = TRUE,
#'                            geo_lambda = 0.5)
#' }
#'
#' @seealso [create_city()], [create_resident()], [get_n_resident()]
#'
#' @export
create_landlord <- function(
    city,
    n_landlord = 5,
    ethnic_aversion = c(1, 1, 1, 0, 0),
    SES_aversion    = c(0, 0, 1, 1, 1),
    dispersed_landlord = c(1, 1, 1, 0, 0),
    geo_difference = FALSE,
    geo_lambda = 1,
    landlord_origin_vec = NULL) {
  
  #======================================================
  # Validation
  #======================================================
  # city
  if (!is.array(city)) {
    stop("'city' must be an array.")
  }
  
  # ethnic_aversion
  if (length(ethnic_aversion) != n_landlord) {
    stop("The length of 'ethnic_aversion' must be 'n_landlord'.")
  }
  if (any(ethnic_aversion < 0)) {
    stop("The elements of 'ethnic_aversion' must be non-negative.")
  }
  
  # SES_aversion
  if (length(SES_aversion) != n_landlord) {
    stop("The length of 'SES_aversion' must be 'n_landlord'.")
  }
  if (any(SES_aversion < 0)) {
    stop("The elements of 'SES_aversion' must be non-negative.")
  }
  
  # dispersed_landlord
  if (length(dispersed_landlord) != n_landlord) {
    stop("The length of 'dispersed_landlord' must be 'n_landlord'.")
  }
  if (!all(dispersed_landlord %in% c(0, 1))) {
    stop("The elements of 'dispersed_landlord' must be 0 or 1.")
  }
  
  #=========================
  # create df of landlord
  #=========================
  landlord <- data.frame(
    ID = seq_len(n_landlord),
    ethnic_aversion = ethnic_aversion,
    SES_aversion    = SES_aversion,
    dispersed = dispersed_landlord
  )
  
  #=========================
  # geo_difference
  #=========================
  city_dim <- dim(city)[1:2]
  
  if (geo_difference) {
    # create a distance matrix
    grid <- expand.grid(x = seq_len(city_dim[1]), y = seq_len(city_dim[2]))
    index_dist <- as.matrix(dist(grid))
    
    # validate landlord_origin_vec or create it
    if (!is.null(landlord_origin_vec)) {
      if (length(landlord_origin_vec) != n_landlord) {
        stop("The length of 'landlord_origin_vec' must be 'n_landlord'.")
      }
      if (any(landlord_origin_vec < 1 | landlord_origin_vec > prod(city_dim))) {
        stop("Elements of 'landlord_origin_vec' must be within valid city indices.")
      }
      if (anyDuplicated(landlord_origin_vec)) {
        stop("'landlord_origin_vec' must not contain duplicates.")
      }
    } else {
      landlord_origin_vec <- sample(seq_len(prod(city_dim)), size = n_landlord)
    }
    
    # calculate the ownership propensity (distance decay from each origin)
    owning_prop <- lapply(landlord_origin_vec, function(p) {
      exp(-geo_lambda * index_dist[p, ])
    })
    owning_prop <- do.call(rbind, owning_prop)
    rownames(owning_prop) <- seq_len(n_landlord)
    
    # For dispersed landlords, replace each row with its own mean so that
    # the propensity is uniform across cells but on the same scale as
    # non-dispersed landlords' distance-decayed propensities.
    if (any(dispersed_landlord == 1)) {
      disp_idx <- which(dispersed_landlord == 1)
      row_means <- rowMeans(owning_prop[disp_idx, , drop = FALSE])
      owning_prop[disp_idx, ] <- matrix(
        row_means,
        nrow = length(disp_idx),
        ncol = ncol(owning_prop)
      )
    }
    
    # decide the owners
    ownership <- vapply(seq_len(ncol(owning_prop)), function(i) {
      origin_match <- which(landlord_origin_vec == i)
      if (length(origin_match) == 1L) {
        origin_match
      } else {
        sample.int(n_landlord, size = 1, prob = owning_prop[, i])
      }
    }, integer(1))
    ownership <- matrix(ownership, nrow = city_dim[1], ncol = city_dim[2])
    
  } else {
    ownership <- matrix(
      sample(seq_len(n_landlord), size = prod(city_dim), replace = TRUE),
      nrow = city_dim[1],
      ncol = city_dim[2]
    )
  }
  
  # make 3D version of ownership (broadcast across the height axis)
  ownership_3D <- array(ownership, dim = dim(city))
  ownership_3D[is.na(city)] <- NA
  
  # output
  list(
    landlord = landlord,
    ownership = ownership_3D
  )
}
