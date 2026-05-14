#' Create a data frame of residents
#'
#' Constructs a data frame describing a population of residents, including
#' each resident's minority status, socioeconomic status (SES), and
#' preference parameters for ethnic and SES similarity in their neighbors.
#' Minority/majority labels are sampled internally based on a specified
#' minority proportion.
#'
#' @param n_resident A positive integer giving the number of residents to
#'   create. Default is `10`.
#' @param minority_prop A single numeric value in `(0, 1)` giving the
#'   target proportion of minority residents. Minority status is assigned
#'   by sampling from `rABM::make_group_labels()` using this proportion.
#'   Default is `0.4`.
#' @param preference_ethnicity A numeric vector of length `n_resident` with
#'   values in `[0, 1]` giving each resident's preference weight for ethnic
#'   similarity in neighbors. Default is `0.5` for all residents.
#' @param SES A numeric vector of length `n_resident` giving each
#'   resident's socioeconomic status. Must be non-negative.
#' @param preference_SES A numeric vector of length `n_resident` with
#'   values in `[0, 1]` giving each resident's preference weight for SES
#'   similarity in neighbors. Default is `0.5` for all residents.
#'
#' @return A data frame with `n_resident` rows and the following columns:
#'   \describe{
#'     \item{`ID`}{Integer resident identifier, from 1 to `n_resident`.}
#'     \item{`minority`}{Integer indicator: `0` for majority residents,
#'       `1` for minority residents.}
#'     \item{`SES`}{The resident's socioeconomic status.}
#'     \item{`preference_ethnicity`}{The resident's preference weight for
#'       ethnic similarity.}
#'     \item{`preference_SES`}{The resident's preference weight for SES
#'       similarity.}
#'   }
#'
#' @details
#' Minority status is assigned by calling
#' `rABM::make_group_labels(n = n_resident, prop = c(1 - minority_prop, minority_prop))`
#' and recoding the result so that majority residents are `0` and
#' minority residents are `1`.
#'
#' The interpretation of `preference_ethnicity` and `preference_SES` is
#' application-specific, but they are typically used as weights in a
#' utility function that scores potential residential locations based on
#' neighbor characteristics. Values close to `1` indicate strong preference
#' for similar neighbors; values close to `0` indicate indifference.
#'
#' @examples
#' \dontrun{
#' # Default population: 10 residents, 40% minority
#' residents <- create_resident()
#'
#' # Custom population: 6 residents, 50% minority, with varied preferences
#' residents <- create_resident(
#'   n_resident = 6,
#'   minority_prop = 0.5,
#'   preference_ethnicity = c(0.8, 0.8, 0.5, 0.5, 0.2, 0.2),
#'   SES = c(1, 2, 3, 1, 2, 3),
#'   preference_SES = rep(0.3, 6)
#' )
#' }
#'
#' @seealso [create_city()], [create_landlord()]
#'
#' @export
create_resident <- function(
    n_resident = 10,
    minority_prop = 0.4,
    preference_ethnicity = rep(0.5, n_resident),
    SES = c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5),
    preference_SES = rep(0.5, n_resident)
) {
  # validation
  if (!is.numeric(n_resident) || n_resident <= 0 || n_resident != as.integer(n_resident)) {
    stop("'n_resident' must be a positive integer.")
  }
  
  stopifnot(
    "The length of 'minority_prop' must be 1 and its range must be from 0 to 1." =
      length(minority_prop) == 1 && minority_prop > 0 && minority_prop < 1
  )
  
  if (length(preference_ethnicity) != n_resident) {
    stop("The length of 'preference_ethnicity' must be 'n_resident'.")
  }
  if (any(preference_ethnicity < 0 | preference_ethnicity > 1)) {
    stop("'preference_ethnicity' must be in [0, 1].")
  }
  
  if (length(SES) != n_resident) {
    stop("The length of 'SES' must be 'n_resident'.")
  }
  if (any(SES < 0)) {
    stop("'SES' must be non-negative.")
  }
  
  if (length(preference_SES) != n_resident) {
    stop("The length of 'preference_SES' must be 'n_resident'.")
  }
  if (any(preference_SES < 0 | preference_SES > 1)) {
    stop("'preference_SES' must be in [0, 1].")
  }
  
  # ethnicity
  minority <- rABM::make_group_labels(
    n = n_resident,
    prop = c(1 - minority_prop, minority_prop)
  )
  minority <- minority - 1
  
  # create a data.frame
  data.frame(
    ID = seq_len(n_resident),
    minority = minority,
    SES = SES,
    preference_ethnicity = preference_ethnicity,
    preference_SES = preference_SES
  )
}
