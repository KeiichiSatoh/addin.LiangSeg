#' Calculate the number of residents in a city
#'
#' Computes the number of occupied housing units in a city, given a vacancy
#' proportion. The total number of houses is determined by counting non-NA
#' cells in the city array.
#'
#' @param city An array representing a city, typically produced by
#'   [create_city()]. Non-NA cells are treated as housing units; NA cells
#'   are excluded from the count.
#' @param vacant_prop A single numeric value in `[0, 1]` giving the proportion
#'   of houses that are vacant. The remaining proportion `(1 - vacant_prop)`
#'   is occupied by residents.
#'
#' @return An integer scalar: the number of residents (i.e. occupied houses),
#'   computed as `round(n_houses * (1 - vacant_prop))`.
#'
#' @details
#' The function counts every non-NA cell in `city` as a distinct house.
#' For a 3D city array (with a height axis), this means each occupied
#' level of a building counts as one house.
#'
#' @examples
#' \dontrun{
#' city <- create_city(zone_dim = c(2, 2), lot_dim = c(3, 3), max_height = 2)
#' get_n_resident(city, vacant_prop = 0.1)
#' }
#'
#' @seealso [create_city()]
#'
#' @export
get_n_resident <- function(city, vacant_prop) {
  # validate the arguments
  if (!is.array(city)) {
    stop("'city' must be an array.")
  }
  if (!is.numeric(vacant_prop) || length(vacant_prop) != 1 ||
      vacant_prop < 0 || vacant_prop > 1) {
    stop("'vacant_prop' must be a single numeric value between 0 and 1.")
  }
  
  # number of total houses (non-NA cells)
  n_houses <- sum(!is.na(city))
  
  # number of residents
  n_resident <- round(n_houses * (1 - vacant_prop))
  
  n_resident
}
