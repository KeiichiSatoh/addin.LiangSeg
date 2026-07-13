#' Randomly re-draw landlords' discrimination aversion scores by geographic type
#'
#' Re-samples the ethnic (\code{minority_aversion}) and/or socioeconomic
#' (\code{SES_aversion}) aversion scores of landlords in an existing
#' \code{Game} object, separately for geographically dispersed and
#' concentrated landlords (\code{G$landlord$dispersed}). Whether binary or
#' numeric sampling is used is inferred automatically from the current
#' values already stored in \code{G$landlord}, so the function can be
#' applied to a \code{Game} built with either
#' \code{landlord_hostile_input_type = "binary"} or \code{"numeric"} (see
#' \code{\link{set_segGame}}) without needing to specify the type explicitly.
#'
#' @param G A \code{Game} object (as returned by \code{\link{set_segGame}})
#'   containing a \code{landlord} field with columns \code{minority_aversion},
#'   \code{SES_aversion}, and \code{dispersed} (\code{1} = geographically
#'   dispersed, \code{0} = concentrated).
#' @param dispersed_aversion_minority Numeric vector of length 2, or
#'   \code{NULL}. Controls re-sampling of \code{minority_aversion} for
#'   dispersed landlords (\code{dispersed == 1}). Interpretation depends on
#'   the aversion type already present in \code{G$landlord$minority_aversion}
#'   (binary if it currently has 2 or fewer unique values, numeric
#'   otherwise): if binary, only the first element is used, as the
#'   probability that \code{minority_aversion = 1} (passed to
#'   \code{sample2(x = c(0,1), prob = c(1 - p, p))}); if numeric, the two
#'   elements are used as \code{c(mean, sd)} for \code{rnorm()}. Set to
#'   \code{NULL} to leave dispersed landlords' \code{minority_aversion}
#'   unchanged. Default: \code{c(0.5, 0.1)}.
#' @param dispersed_aversion_SES Numeric vector of length 2, or \code{NULL}.
#'   Same as \code{dispersed_aversion_minority}, but for \code{SES_aversion}
#'   of dispersed landlords (type inferred from
#'   \code{G$landlord$SES_aversion}). Default: \code{c(0.5, 0.1)}.
#' @param concentrated_aversion_minority Numeric vector of length 2, or
#'   \code{NULL}. Same as \code{dispersed_aversion_minority}, but for
#'   concentrated landlords (\code{dispersed == 0}). Default:
#'   \code{c(0.5, 0.1)}.
#' @param concentrated_aversion_SES Numeric vector of length 2, or
#'   \code{NULL}. Same as \code{dispersed_aversion_SES}, but for
#'   concentrated landlords (\code{dispersed == 0}). Default:
#'   \code{c(0.5, 0.1)}.
#'
#' @return A new \code{Game} object (a deep copy of \code{G} via
#'   \code{rABM::copy_obj()}), with \code{landlord$minority_aversion} and/or
#'   \code{landlord$SES_aversion} replaced for the requested
#'   dispersed/concentrated subsets. \code{G} itself is left unmodified.
#'   Groups for which the corresponding argument was \code{NULL} retain
#'   their original aversion values.
#'
#' @details
#' The function re-samples up to four independent slices of
#' \code{G$landlord}: (dispersed x minority), (dispersed x SES),
#' (concentrated x minority), (concentrated x SES). Each slice is controlled
#' by its own argument and can be skipped individually by passing
#' \code{NULL}, making it possible to, for example, only shift the
#' ethnic-aversion profile of dispersed landlords while leaving all other
#' aversion values untouched.
#'
#' Because the binary-vs-numeric decision is inferred from the *current*
#' number of unique values in \code{G$landlord}, this function is only
#' intended for updating aversion scores in place; it does not switch a
#' \code{Game} between \code{landlord_hostile_input_type = "binary"} and
#' \code{"numeric"} representations.
#'
#' @examples
#' \dontrun{
#' library(rABM)
#'
#' G <- set_segGame(landlord_hostile_input_type = "binary")
#'
#' # Make dispersed landlords much more ethnically averse, leave SES alone,
#' # and leave concentrated landlords entirely untouched
#' G2 <- modify_landlord_aversion(
#'   G,
#'   dispersed_aversion_minority     = c(0.9, NA),
#'   dispersed_aversion_SES          = NULL,
#'   concentrated_aversion_minority  = NULL,
#'   concentrated_aversion_SES       = NULL
#' )
#'
#' # Numeric aversion scores: shift means/sds for both groups
#' G_num <- set_segGame(landlord_hostile_input_type = "numeric")
#' G_num2 <- modify_landlord_aversion(
#'   G_num,
#'   dispersed_aversion_minority    = c(0.8, 0.05),
#'   dispersed_aversion_SES         = c(0.2, 0.05),
#'   concentrated_aversion_minority = c(0.3, 0.05),
#'   concentrated_aversion_SES      = c(0.3, 0.05)
#' )
#' }
#'
#' @seealso \code{\link{set_segGame}}, \code{\link[rABM]{copy_obj}}
#'
#' @export

modify_landlord_aversion <- function(G, 
                                     dispersed_aversion_minority = c(0.5, 0.1),
                                     dispersed_aversion_SES = c(0.5, 0.1),
                                     concentrated_aversion_minority = c(0.5, 0.1),
                                     concentrated_aversion_SES = c(0.5, 0.1)){
  # landlord
  landlord <- G$landlord
  
  
  # judge if aversion score is binrary or numeric
  aversion_type <- list(minority = "binary", SES = "binary")
  if(length(unique(landlord$minority_aversion)) > 2){
    aversion_type$minority <- "numeric"
  }
  if(length(unique(landlord$SES_aversion)) > 2){
    aversion_type$SES <- "numeric"
  }
  
  # dispersed ---------
  ### minority
  if(!is.null(dispersed_aversion_minority)){
    if(aversion_type$minority == "binary"){
      landlord$minority_aversion[landlord$dispersed==1] <- sample2(x = c(0, 1), 
                                                                   size = length(landlord$minority_aversion[landlord$dispersed==1]),
                                                                   replace = TRUE, 
                                                                   prob = c(1 - dispersed_aversion_minority[1], dispersed_aversion_minority[1]) 
      )
    }else{
      landlord$minority_aversion[landlord$dispersed==1] <- rnorm(n = length(landlord$minority_aversion[landlord$dispersed==1]),
                                                                 mean = dispersed_aversion_minority[1], 
                                                                 sd = dispersed_aversion_minority[2]) 
    }
  }
  
  ### SES
  if(!is.null(dispersed_aversion_SES)){
    if(aversion_type$SES == "binary"){
      landlord$SES_aversion[landlord$dispersed==1] <- sample2(x = c(0, 1), 
                                                              size = length(landlord$SES_aversion[landlord$dispersed==1]),
                                                              replace = TRUE, 
                                                              prob = c(1 - dispersed_aversion_SES[1], dispersed_aversion_SES[1]) 
      )
    }else{
      landlord$SES_aversion[landlord$dispersed==1] <- rnorm(n = length(landlord$SES_aversion[landlord$dispersed==1]),
                                                            mean = dispersed_aversion_SES[1], 
                                                            sd = dispersed_aversion_SES[2]) 
    }
  }
  
  # concentrated ---------
  ### minority
  if(!is.null(concentrated_aversion_minority)){
    if(aversion_type$minority == "binary"){
      landlord$minority_aversion[landlord$dispersed==0] <- sample2(x = c(0, 1), 
                                                                   size = length(landlord$minority_aversion[landlord$dispersed==0]),
                                                                   replace = TRUE, 
                                                                   prob = c(1 - concentrated_aversion_minority[1], concentrated_aversion_minority[1]) 
      )
    }else{
      landlord$minority_aversion[landlord$dispersed==0] <- rnorm(n = length(landlord$minority_aversion[landlord$dispersed==0]),
                                                                 mean = concentrated_aversion_minority[1], 
                                                                 sd = concentrated_aversion_minority[2])
    }
  }
  
  ### SES
  if(!is.null(concentrated_aversion_SES)){
    if(aversion_type$SES == "binary"){
      landlord$SES_aversion[landlord$dispersed==0] <- sample2(x = c(0, 1), 
                                                              size = length(landlord$SES_aversion[landlord$dispersed==0]),
                                                              replace = TRUE, 
                                                              prob = c(1 - concentrated_aversion_SES[1], concentrated_aversion_SES[1]) 
      )
    }else{
      landlord$SES_aversion[landlord$dispersed==0] <- rnorm(n = length(landlord$SES_aversion[landlord$dispersed==0]),
                                                            mean = concentrated_aversion_SES[1], 
                                                            sd = concentrated_aversion_SES[2]) 
    }
  }
  
  # attach the changes and return
  copied_G <- rABM::copy_obj(G)
  copied_G$landlord <- landlord
  copied_G
}
