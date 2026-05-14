#' Randomly Allocate Initial Residents to Houses in a City
#'
#' Randomly assigns a specified number of residents to available houses within
#' a city array. Houses marked as \code{NA} in the input array are excluded from
#' allocation, representing unavailable locations.
#'
#' @param city An array (created by \code{\link[base]{array}} or
#'   \code{\link[base]{matrix}}) representing the city. Each element corresponds
#'   to a house (cell). Elements with \code{NA} values are treated as unavailable
#'   and excluded from allocation. The array can have any dimensions.
#' @param n_resident A positive integer of length 1 specifying the number of
#'   residents to allocate. Must not exceed the number of available
#'   (non-\code{NA}) houses in \code{city}. Defaults to \code{10}.
#'
#' @return An integer vector of length \code{n_resident} containing the linear
#'   indices (positions in the flattened array, column-major order) of the
#'   houses to which residents are allocated. Indices are sampled without
#'   replacement from the available houses.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Validates that \code{n_resident} is a positive integer of length 1.
#'   \item Validates that \code{city} is an array (matrices are accepted, as
#'     they inherit from \code{array}).
#'   \item Constructs an index array of the same shape as \code{city}, masking
#'     positions where \code{city} is \code{NA}.
#'   \item Checks that the number of available houses is at least
#'     \code{n_resident}; otherwise, an error is thrown.
#'   \item Samples \code{n_resident} indices without replacement from the
#'     available houses using \code{\link[base]{sample.int}} to avoid the
#'     ambiguous behavior of \code{\link[base]{sample}} when given a
#'     length-one vector.
#' }
#'
#' The returned indices are linear (column-major) positions corresponding to
#' the flattened \code{city} array, which can be used directly to subset
#' \code{city} or related arrays of the same dimensions.
#'
#' @examples
#' # 2D city (5x5 grid) with some unavailable houses
#' city <- array(1, dim = c(5, 5))
#' city[1, 1] <- NA
#' city[3, 3] <- NA
#' initial_random_allocation(city, n_resident = 5)
#'
#' # 3D city
#' city_3d <- array(1, dim = c(4, 4, 2))
#' initial_random_allocation(city_3d, n_resident = 8)
#'
#' # Matrices are also accepted
#' city_mat <- matrix(1, nrow = 4, ncol = 4)
#' initial_random_allocation(city_mat, n_resident = 3)
#'
#' \dontrun{
#' # Error: too many residents for available houses
#' small_city <- array(1, dim = c(2, 2))
#' initial_random_allocation(small_city, n_resident = 10)
#'
#' # Error: non-integer n_resident
#' initial_random_allocation(array(1, dim = c(5, 5)), n_resident = 3.5)
#' }
#'
#' @export
initial_random_allocation <- function(city, n_resident = 10){
  # validation
  stopifnot(
    "'n_resident' must be a positive integer of length 1." =
      length(n_resident) == 1 && is.numeric(n_resident) &&
      n_resident > 0 && n_resident == as.integer(n_resident),
    "'city' must be an array." = is.array(city)
  )
  
  # preparation
  city_dim <- dim(city)
  house_index <- array(seq_len(prod(city_dim)), dim = city_dim)
  house_index[is.na(city)] <- NA
  house_index_vect <- na.omit(as.vector(house_index))
  
  # validate the number of houses
  if(length(house_index_vect) < n_resident){
    stop("Number of houses are less than 'n_resident'.")
  }
  
  # allocation (use sample.int to avoid sample()'s length-1 ambiguity)
  house_index_vect[sample.int(length(house_index_vect), size = n_resident)]
}