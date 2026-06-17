#' Set up a Segregation ABM Game object
#'
#' Constructs a fully configured \code{Game} object for simulating residential
#' segregation dynamics. The model includes residents, landlords, and a
#' gridded city structure, and incorporates ethnic and socioeconomic
#' preferences of residents and discriminatory behaviour of landlords.
#'
#' @param city_zone_dim Integer vector of length 2. Number of zones along
#'   each dimension of the city grid (rows, columns). Default: \code{c(3, 3)}.
#' @param city_lot_dim Integer vector of length 2. Number of lots per zone
#'   along each dimension (rows, columns). Default: \code{c(5, 5)}.
#' @param city_max_height Positive integer. Maximum number of floors per lot.
#'   Default: \code{3}.
#' @param city_height_prop Numeric vector of length \code{city_max_height}.
#'   Relative proportions used to sample building heights. Values need not sum
#'   to 1 (they are normalised internally). Default: random draws from
#'   \code{rlnorm(city_max_height, 1, 1)}, evaluated at call time.
#' @param city_vacancy_rate Numeric scalar in \code{(0, 1)}. Proportion of
#'   housing units left vacant at initialisation. Default: \code{0.15}.
#' @param landlord_ratio_against_resident Positive numeric scalar. Number of
#'   landlords as a proportion of the number of residents (before rounding up).
#'   The result is rounded up and floored at 3. Default: \code{0.01}.
#' @param resident_minority_prop Numeric scalar in \code{(0, 1)}. Target
#'   proportion of minority residents. Default: \code{0.3}.
#' @param resident_SES_lambda_majo Positive numeric scalar. Poisson rate
#'   (\eqn{\lambda}) for the SES distribution of majority residents.
#'   Default: \code{2.5}.
#' @param resident_SES_lambda_mino Positive numeric scalar. Poisson rate
#'   (\eqn{\lambda}) for the SES distribution of minority residents.
#'   Default: \code{1.5}.
#' @param resident_preference_ethnicity_majo Numeric vector of length 2
#'   \code{c(mean, sd)}. Parameters of the normal distribution from which
#'   ethnic co-group preference weights are drawn for majority residents.
#'   Default: \code{c(5, 0)}.
#' @param resident_preference_ethnicity_mino Numeric vector of length 2
#'   \code{c(mean, sd)}. Same as above for minority residents.
#'   Default: \code{c(3, 1)}.
#' @param resident_preference_SES_majo Numeric vector of length 2
#'   \code{c(mean, sd)}. Parameters for the SES preference distribution of
#'   majority residents. Default: \code{c(5, 0)}.
#' @param resident_preference_SES_mino Numeric vector of length 2
#'   \code{c(mean, sd)}. Same as above for minority residents.
#'   Default: \code{c(5, 0)}.
#' @param landlord_hostile_input_type Character scalar, one of
#'   \code{"binary"} or \code{"numeric"}. Whether landlord aversion scores
#'   are binary (0/1 via bivariate Bernoulli) or continuous (via truncated
#'   multivariate normal). Default: \code{"binary"}.
#' @param landlord_ethnic_hostile_prop Numeric scalar in \code{[0, 1]}.
#'   When \code{landlord_hostile_input_type = "binary"}, the probability that
#'   a landlord is ethnically hostile. Default: \code{0.5}.
#' @param landlord_SES_hostile_prop Numeric scalar in \code{[0, 1]}.
#'   When \code{landlord_hostile_input_type = "binary"}, the probability that
#'   a landlord is SES-hostile. Default: \code{0.5}.
#' @param landlord_hostile_phi Numeric scalar. Correlation parameter
#'   (\eqn{\phi}) between ethnic and SES hostility in the bivariate Bernoulli
#'   model. Default: \code{0.1}.
#' @param landlord_ethnic_aversion_rnorm Numeric vector of length 2
#'   \code{c(mean, sd)}. When \code{landlord_hostile_input_type = "numeric"},
#'   parameters of the marginal normal distribution for ethnic aversion scores.
#'   Default: \code{c(0.5, 0.1)}.
#' @param landlord_SES_aversion_rnorm Numeric vector of length 2
#'   \code{c(mean, sd)}. When \code{landlord_hostile_input_type = "numeric"},
#'   parameters of the marginal normal distribution for SES aversion scores.
#'   Default: \code{c(0.5, 0.1)}.
#' @param landlord_aversion_cor Numeric scalar in \code{[-1, 1]}. Correlation
#'   between ethnic and SES aversion in the numeric model. Default: \code{0.1}.
#' @param landlord_dispersed_prop Numeric scalar in \code{[0, 1]}. Proportion
#'   of landlords designated as geographically dispersed (i.e., owning
#'   properties spread across zones). Default: \code{0.3}.
#' @param landlord_geo_difference Logical. If \code{TRUE}, landlords differ in
#'   their geographic ownership patterns (concentrated vs. dispersed).
#'   Default: \code{TRUE}.
#' @param landlord_geo_theta Positive numeric scalar. Concentration parameter
#'   controlling how strongly concentrated landlords cluster their holdings.
#'   Default: \code{3}.
#' @param block_eth_eff Numeric scalar. Effect size of block-level minority
#'   proportion on landlord decline probability when
#'   \code{block_eth_context = TRUE} is passed to \code{landlord_decline()} or
#'   \code{resident_move_seq()}. Default: \code{0.1}.
#' @param block_SES_eff Numeric scalar. Effect size of block-level SES
#'   composition (reserved for future use). Default: \code{0.1}.
#' @param default_select_resident_prop Numeric scalar in \code{(0, 1]}.
#'   Default proportion of residents selected per step by
#'   \code{select_resident()}. Stored in \code{G$act_defaults} and used as
#'   the default argument value of \code{select_resident(prop)}. Default:
#'   \code{0.25}.
#' @param default_criteria_resident Character vector. Default criteria used by
#'   \code{resident_choose_house()} and \code{resident_move_seq()} to
#'   determine which components enter the resident utility score. Any subset
#'   of \code{c("intercept", "ethnicity", "SES")}. Default: all three.
#' @param default_resident_softmax_beta Positive numeric scalar. Default
#'   inverse temperature (\eqn{\beta}) for the softmax transformation of
#'   resident utility scores. Stored in \code{G$act_defaults} and used as the
#'   default argument value of \code{beta} in \code{resident_choose_house()}
#'   and \code{resident_move_seq()}. Default: \code{1}.
#' @param default_criteraia_landlord Character vector. Default criteria used
#'   by \code{landlord_decline()} and \code{resident_move_seq()} to determine
#'   which components enter the landlord decline score. Any subset of
#'   \code{c("intercept", "ethnicity", "SES")}. Note: the argument name
#'   contains a typo (\code{criteraia}) for backward compatibility. Default:
#'   all three.
#' @param default_block_eth_context Logical. Default value of
#'   \code{block_eth_context} in \code{landlord_decline()} and
#'   \code{resident_move_seq()}. If \code{TRUE}, the block-level minority
#'   proportion modifies the landlord's decline probability. Default:
#'   \code{FALSE}.
#' @param default_include_landlord_seq Logical. Default value of
#'   \code{include_landlord} in \code{resident_move_seq()}. If \code{FALSE},
#'   landlord screening is skipped and all selected residents are always
#'   accepted. Useful for isolating the effect of resident preferences alone.
#'   Default: \code{TRUE}.
#'
#' @return A \code{Game} object (R6 class) configured with the following
#'   fields:
#'   \describe{
#'     \item{States}{
#'       \code{landlord}, \code{resident}, \code{resident_mat_ethnicity},
#'       \code{resident_mat_SES}, \code{house}, \code{house_neighbor_ind},
#'       \code{settings}, \code{record}, \code{act_defaults}
#'     }
#'     \item{Active states}{
#'       \code{house_neib_ethnicity}, \code{house_neib_SES},
#'       \code{house_intercept}, \code{block_minority_prop}
#'     }
#'     \item{Act functions (batch update)}{
#'       \code{select_resident}, \code{resident_choose_house},
#'       \code{landlord_decline}, \code{resident_move}
#'     }
#'     \item{Act functions (sequential update)}{
#'       \code{select_resident}, \code{resident_move_seq}
#'     }
#'     \item{Plot functions}{
#'       \code{plot_city}, \code{plot_ownership}
#'     }
#'     \item{Report functions}{
#'       \code{report_segregation}, \code{report_landlord_stat}
#'     }
#'   }
#'
#'   The \code{settings} state contains the following named elements:
#'   \code{resident_n}, \code{landlord_n}, \code{house_n} (number of
#'   habitable units, i.e. cells with a landlord), \code{cell_n},
#'   \code{city_dim}, \code{city_zone_dim}, \code{city_lot_dim},
#'   \code{city_max_height}, \code{block_eth_eff}, \code{block_SES_eff}.
#'
#'   The \code{act_defaults} state stores the six default parameter values
#'   listed above (\code{select_resident_prop}, \code{criteria_resident},
#'   \code{resident_softmax_beta}, \code{criteraia_landlord},
#'   \code{block_eth_context}, \code{include_landlord_seq}), and can be
#'   updated at run time to change the behaviour of Act functions without
#'   re-initialising the \code{Game} object.
#'
#' @details
#' \subsection{Model overview}{
#'   The model simulates residential sorting driven by two mechanisms:
#'   (1) residents' preferences for co-ethnic and same-SES neighbours, and
#'   (2) landlords' discriminatory screening of applicants by ethnicity and/or
#'   SES. Two update modes are provided: a batch mode and a sequential mode.
#'
#'   \strong{Batch mode} (\code{plan = c("select_resident",
#'   "resident_choose_house", "landlord_decline", "resident_move")}):
#'   All selected residents choose their candidate houses simultaneously,
#'   landlords then screen all applications at once, and accepted residents
#'   move. When two residents choose the same unit, the later one (after
#'   random reordering) stays put.
#'
#'   \strong{Sequential mode} (\code{plan = c("select_resident",
#'   "resident_move_seq")}):
#'   Selected residents are processed one at a time in random order. For each
#'   resident, the full cycle of (1) house choice, (2) landlord screening, and
#'   (3) movement is completed before the next resident is processed.
#'   Neighbourhood composition scores (Active states) are recomputed at each
#'   step, so later residents see the updated positions of earlier movers.
#'   Setting \code{include_landlord = FALSE} (or via
#'   \code{default_include_landlord_seq}) disables landlord screening
#'   entirely, so every resident always moves to their chosen unit. This is
#'   useful for isolating the effect of residential preferences alone.
#'   Both modes are theoretically equivalent when \code{include_landlord =
#'   TRUE} and can be used to verify convergence properties.
#' }
#'
#' \subsection{City structure}{
#'   The city is a 3-D array of dimension
#'   \code{(city_zone_dim[1] * city_lot_dim[1],}
#'   \code{city_zone_dim[2] * city_lot_dim[2],}
#'   \code{city_max_height)}.
#'   Each cell represents one housing unit. Zones group lots for block-level
#'   analysis.
#' }
#'
#' \subsection{Resident score computation and softmax choice}{
#'   The raw utility score for resident \eqn{i} considering house \eqn{j} is:
#'   \deqn{s_{ij} = w_{\text{intercept},i} \cdot v_{\text{intercept},j}
#'                + w_{\text{eth},i} \cdot v_{\text{eth},j}
#'                + w_{\text{SES},i} \cdot v_{\text{SES},j}}
#'   where \eqn{w} are preference weights and \eqn{v} are neighbourhood
#'   composition values from the active-state functions. Negative preference
#'   values (representing out-group preference) are supported.
#'
#'   Raw scores are converted to selection probabilities via the softmax
#'   transformation (implemented by \code{rABM::prob_softmax()}):
#'   \deqn{P(j \mid i) = \frac{\exp(\beta \cdot s_{ij})}
#'                             {\sum_{j'} \exp(\beta \cdot s_{ij'})}}
#'   The inverse temperature parameter \eqn{\beta} (argument \code{beta},
#'   default value stored in \code{act_defaults$resident_softmax_beta})
#'   controls the degree of rationality:
#'   \itemize{
#'     \item \eqn{\beta \to 0}: near-uniform (bounded rationality / random
#'       choice)
#'     \item \eqn{\beta = 1}: standard probabilistic choice
#'     \item \eqn{\beta \to \infty}: near-deterministic optimal choice
#'       (perfect rationality)
#'   }
#'   This formulation is consistent with McFadden's random utility model, in
#'   which Gumbel-distributed unobserved heterogeneity yields softmax choice
#'   probabilities. After transformation, occupied units and units without a
#'   landlord are masked to zero; \code{rABM::sample_weighted()} re-normalises
#'   internally before sampling.
#' }
#'
#' \subsection{Landlord decline model}{
#'   The log-odds of declining resident \eqn{i} for landlord \eqn{l} is:
#'   \deqn{\text{logit}(p_{\text{decline}}) =
#'         \sum_k \beta_k \cdot a_{lk} \cdot x_{ik} - 6}
#'   where \eqn{\beta_k} are fixed scaling coefficients
#'   (\code{intercept}: 1, \code{ethnicity}: 12, \code{SES}: 12/5),
#'   \eqn{a_{lk}} is landlord \eqn{l}'s aversion on dimension \eqn{k}, and
#'   \eqn{x_{ik}} is resident \eqn{i}'s profile on dimension \eqn{k}.
#'   For SES, the resident's score is reversed (\eqn{5 - \text{SES}}) so that
#'   higher SES residents face lower decline probability.
#'
#'   The baseline offset \code{-6} sets the decline probability to almost 0\%
#'   when no intercept is included. When \code{intercept} is included as a
#'   criterion (\eqn{\beta_{\text{intercept}} = 1},
#'   \eqn{a_{l,\text{intercept}} = 1}), a residual decline probability of
#'   approximately 0.7\% remains, representing random ``capricious'' refusals
#'   by landlords irrespective of resident attributes.
#'
#'   Landlord screening can be disabled entirely in the sequential mode by
#'   setting \code{include_landlord = FALSE} in \code{resident_move_seq()},
#'   or by changing \code{G$act_defaults$include_landlord_seq} at run time.
#' }
#'
#' @examples
#' \dontrun{
#' library(rABM)
#'
#' # Build the Game object with default settings
#' G <- set_segGame()
#'
#' # --- Batch mode ---
#' G_batch <- run_Game(
#'   G,
#'   plan = c("select_resident", "resident_choose_house",
#'            "landlord_decline", "resident_move"),
#'   times = 10,
#'   fields_to_save = G$notes$fields_to_save
#' )
#'
#' # --- Sequential mode with landlord screening (default) ---
#' G_seq <- run_Game(
#'   G,
#'   plan = c("select_resident", "resident_move_seq"),
#'   times = 10,
#'   fields_to_save = G$notes$fields_to_save
#' )
#'
#' # --- Sequential mode without landlord screening ---
#' # Option A: via argument
#' G$resident_move_seq(include_landlord = FALSE)
#'
#' # Option B: change act_defaults globally before running
#' G$act_defaults$include_landlord_seq <- FALSE
#' G_seq_noland <- run_Game(
#'   G,
#'   plan = c("select_resident", "resident_move_seq"),
#'   times = 10,
#'   fields_to_save = G$notes$fields_to_save
#' )
#'
#' # Visualise the city at the final time step
#' G_batch$plot_city(time = 11)
#'
#' # Report dissimilarity index over time
#' G_batch$report_segregation(level = "city", show = "ethnicity",
#'                            index = "D", log = 1:10)
#'
#' # Landlord-level descriptive statistics at current time
#' G_batch$report_landlord_stat()
#'
#' # Plot landlord ownership and aversion map
#' G_batch$plot_ownership(show = "ethnicity")
#' }
#'
#' @seealso
#'   \code{\link[rABM]{Game}}, \code{\link[rABM]{run_Game}},
#'   \code{\link[rABM]{value_of}}, \code{\link[rABM]{prob_softmax}},
#'   \code{\link[rABM]{sample_weighted}}, \code{\link{create_resident}},
#'   \code{\link{create_landlord}}, \code{\link{create_city}}
#'
#' @importFrom matrixStats rowMeans2 rowCumsums rowSums2 rowCounts
#' @importFrom ggplot2 ggplot aes geom_tile geom_line geom_vline geom_hline
#'   geom_text scale_fill_gradient2 scale_x_continuous scale_y_continuous
#'   coord_fixed labs theme element_blank element_rect xlab ylab ylim
#'   theme_bw
#' @importFrom segregation dissimilarity mutual_total
#'
#' @export
set_segGame <- function(
    city_zone_dim = c(3,3),
    city_lot_dim  = c(5,5),
    city_max_height = 3,
    city_height_prop = rlnorm(city_max_height, 1, 1), # User can change the distribution
    city_vacancy_rate = 0.15,
    landlord_ratio_against_resident = 0.01,
    resident_minority_prop = 0.3,
    resident_SES_lambda_majo = 2.5,
    resident_SES_lambda_mino = 1.5,
    resident_preference_ethnicity_majo = c(5,0),      # c(mean, sd) input to rnorm
    resident_preference_ethnicity_mino = c(3,1),      
    resident_preference_SES_majo = c(5,0),            
    resident_preference_SES_mino = c(5,0),            
    landlord_hostile_input_type = c("binary", "numeric"),
    landlord_ethnic_hostile_prop = 0.5,     
    landlord_SES_hostile_prop = 0.5,
    landlord_hostile_phi = 0.1,
    landlord_ethnic_aversion_rnorm = c(0.5, 0.1),
    landlord_SES_aversion_rnorm = c(0.5, 0.1),
    landlord_aversion_cor = 0.1,
    landlord_dispersed_prop = 0.3,
    landlord_geo_difference = TRUE,
    landlord_geo_theta = 3,
    block_eth_eff = 0.1,
    block_SES_eff = 0.1,
    default_select_resident_prop = 0.25,
    default_criteria_resident = c("intercept", "ethnicity", "SES"),
    default_resident_softmax_beta = 1,
    default_criteraia_landlord = c("intercept", "ethnicity", "SES"),
    default_block_eth_context = FALSE,
    default_include_landlord_seq = TRUE
    ){
  
  #==========================================
  # create a Game object
  #==========================================
  G <- Game()
  
  # note
  G$notes$fields_to_save <- c("landlord","resident","house", "record")
  
  #==========================================
  # State
  #==========================================
  
  
  # city
  city <- create_city(zone_dim = city_zone_dim, 
                      lot_dim = city_lot_dim, 
                      max_height = city_max_height, 
                      height_prop = city_height_prop)
  
  # calculate n -------------
  resident_n <- get_n_resident(city, vacant_prop = city_vacancy_rate)
  landlord_n <- ceiling(resident_n * landlord_ratio_against_resident)
  if(landlord_n < 3){
    landlord_n <- 3
  }
  
  # resident SES
  gr_tb <- table(make_group_labels(n = resident_n, prop = c(1 - resident_minority_prop,
                                                            resident_minority_prop)))
  resident_SES <- c(
    rpois(n = gr_tb[1], lambda = resident_SES_lambda_majo),
    rpois(n = gr_tb[2], lambda = resident_SES_lambda_mino)
  )
  
  # resident coeth pref
  resident_preference_ethnicity <- c(
    rnorm(n = gr_tb[1], 
          mean = resident_preference_ethnicity_majo[1], 
          sd = resident_preference_ethnicity_majo[2]),
    rnorm(n = gr_tb[2], 
          mean = resident_preference_ethnicity_mino[1], 
          sd = resident_preference_ethnicity_mino[2]))
  
  # resident SES pref
  resident_preference_SES <- c(
    rnorm(n = gr_tb[1], 
          mean = resident_preference_SES_majo[1], 
          sd = resident_preference_SES_majo[2]),
    rnorm(n = gr_tb[2], 
          mean = resident_preference_SES_mino[1], 
          sd = resident_preference_SES_mino[2]))
  
  
  # resident --------------
  resident <- create_resident(n_resident = resident_n, 
                              minority_prop = resident_minority_prop, 
                              preference_ethnicity = resident_preference_ethnicity, 
                              SES = resident_SES, 
                              preference_SES = resident_preference_SES, 
                              preference_intercept = rep(1, resident_n))
  
  
  # landlord ----------------
  # calculate landlord attribute
  landlord_hostile_input_type <- match.arg(landlord_hostile_input_type)
  if(landlord_hostile_input_type == "binary"){
    #---- binary -------------------
    landlord_aversion_out <- rbivbinom(n = landlord_n, 
                                       p = c(landlord_ethnic_hostile_prop, landlord_SES_hostile_prop),
                                       Sigma = matrix(c(1, landlord_hostile_phi, 
                                                        landlord_hostile_phi, 1), nrow = 2))
  }else{
    #---- numeric ------------------
    landlord_aversion_out <- rmvnorm(n = landlord_n, 
                                     mu = c(landlord_ethnic_aversion_rnorm[1],
                                            landlord_SES_aversion_rnorm[1]), 
                                     sigma = c(landlord_ethnic_aversion_rnorm[2],
                                               landlord_SES_aversion_rnorm[2]), 
                                     R = matrix(c(1, landlord_aversion_cor,
                                                  landlord_aversion_cor, 1), nrow = 2))
  }
  
  # landlord disparsed
  landlord_dispersed <- make_group_labels(
    landlord_n, 
    c(1 - landlord_dispersed_prop, landlord_dispersed_prop)) - 1
  landlord_dispersed <- sample(landlord_dispersed)
  
  # create landlord  
  landlord_out <- create_landlord(city = city, 
                                  n_landlord = landlord_n, 
                                  minority_aversion = landlord_aversion_out[,1],
                                  SES_aversion = landlord_aversion_out[,2],
                                  intercept_aversion = rep(1, landlord_n),
                                  dispersed_landlord = landlord_dispersed,
                                  geo_difference = landlord_geo_difference,
                                  geo_lambda = landlord_geo_theta)
  landlord <- landlord_out$landlord
  
  add_field(G, State(landlord))
  
  
  # house
  house_allocated <- initial_random_allocation(city = city, n_resident = resident_n)
  resident <- data.frame(resident, house = house_allocated)
  
  # ethnicity and SES matrix (for update calculation)
  resident_mat_ethnicity <- matrix(0, resident_n, 2, dimnames = list(seq_len(resident_n), c("majority", "minority")))
  resident_mat_ethnicity[resident$minority==0, "majority"] <- 1
  resident_mat_ethnicity[resident$minority==1, "minority"] <- 1
  
  resident_mat_SES <- matrix(0, resident_n, 6, 
                             dimnames = list(seq_len(resident_n), 0:5))
  for(j in 0:4){
    resident_mat_SES[resident$SES==j, c(as.character(j), as.character(j + 1))] <- 1
  }
  resident_mat_SES[resident$SES==5, "5"] <- 1
  
  add_field(G, 
            State(resident), 
            State(resident_mat_ethnicity),
            State(resident_mat_SES))
  
  # house -------------------
  # cache the neighborhood indices
  house <- data.frame(ID = 1:prod(dim(city)),
                      block = as.vector(city),
                      landlord = as.vector(landlord_out$ownership))
  house_neighbor_ind <- neighbor_indices(city)
  add_field(G, State(house), State(house_neighbor_ind))
  
  # settings ---------------
  # dimension of the city
  city_dim <- c(city_zone_dim[1] * city_lot_dim[1],
                city_zone_dim[2] * city_lot_dim[2])
  
  settings <- list(resident_n = resident_n, 
                   landlord_n = landlord_n,
                   house_n = nrow(na.exclude(house)),
                   cell_n = prod(dim(city)),
                   city_dim = city_dim,
                   city_zone_dim = city_zone_dim,
                   city_lot_dim  = city_lot_dim,
                   city_max_height = city_max_height,
                   block_eth_eff = block_eth_eff,
                   block_SES_eff = block_SES_eff)
  add_field(G, State(settings))
  
  #--- record ---------------
  record <- list(selected_house = c(),
                 selected_resident = seq_len(resident_n),
                 landlord_decision = c())
  add_field(G, State(record))
  
  #--- act_variables --------
  act_defaults = list(
    select_resident_prop = default_select_resident_prop,
    criteria_resident = default_criteria_resident,
    resident_softmax_beta = default_resident_softmax_beta,
    criteraia_landlord = default_criteraia_landlord,
    block_eth_context = default_block_eth_context,
    include_landlord_seq = default_include_landlord_seq
    )
  add_field(G, State(act_defaults))
  
  
  #===================================
  # Active state
  #===================================
  
  #------ same_ethnicity_prop in each neighborhood --------
  house_neib_ethnicity <- function(){
    # resident_ID
    resident_ID <- self$record$selected_resident
    
    # create a vector of house that indicates residents' ethnicity
    house_resident_ethnicity <- rep(NA, length = self$settings$cell_n)
    house_resident_ethnicity[self$resident$house] <- self$resident$minority
    
    # calculate the proportion of each ethnicity for each house
    neighbor_ethnicity <- t(apply(self$house_neighbor_ind, MARGIN = 1, 
                                  FUN = function(X){house_resident_ethnicity[X]}))
    neighbor_minority_prop <- matrixStats::rowMeans2(neighbor_ethnicity, na.rm = TRUE)
    neighbor_majority_prop <- 1 - neighbor_minority_prop
    neighbor_ethnicity_prop <- matrix(c(neighbor_majority_prop, neighbor_minority_prop), self$settings$cell_n, 2)
    neighbor_ethnicity_prop[is.nan(neighbor_ethnicity_prop)] <- 0
    
    # same ethnic proportion for each resident 
    same_eth_prop <- t(neighbor_ethnicity_prop %*% t(self$resident_mat_ethnicity[resident_ID, ,drop = FALSE]))
    same_eth_prop
  }
  
  add_field(G, Active(house_neib_ethnicity))
  
  #------ same SES_prop in each neghborhood --------
  house_neib_SES <- function(){
    # resident_ID
    resident_ID <- self$record$selected_resident
    
    # create a vector of house that indicates residents' SES
    house_resident_SES <- rep(NA, length = self$settings$cell_n)
    house_resident_SES[self$resident$house] <- self$resident$SES
    
    # calculate the proportion of each ethnicity for each house
    neighbor_SES <- t(apply(self$house_neighbor_ind, MARGIN = 1, 
                            FUN = function(X){house_resident_SES[X]}))
    SES_count <- lapply(0:5, function(i){
      matrixStats::rowCounts(neighbor_SES, value = i, na.rm = TRUE)})
    SES_count <- do.call(cbind, SES_count)  # Note: column labels are 0:5
    SES_prop <- SES_count/matrixStats::rowSums2(SES_count)
    SES_prop[is.nan(SES_prop)] <- 0
    
    # same SES(including + 1 level) proportion for each resident 
    same_SES_prop <- t(SES_prop %*% t(self$resident_mat_SES[resident_ID, ,drop = FALSE]))
    same_SES_prop
  }
  add_field(G, Active(house_neib_SES))
  
  
  #------- random ----------------------------
  house_intercept <- function(){
    # resident_ID
    resident_ID <- self$record$selected_resident
    
    mat <- matrix(1, length(resident_ID), self$settings$cell_n)
    mat/self$settings$cell_n
  }
  add_field(G, Active(house_intercept))
  
  #------- ethnic context by block ---------------------
  
  block_minority_prop <- function(){
    house <- data.frame(
      block = self$house$block,
      minority = 0
    )
    house$minority[self$resident$house] <- self$resident$minority
    
    # tabulate
    minority_prop <- prop.table(table(house$block, house$minority), 1)[,"1"]
    
    # return
    minority_prop
  }
  add_field(G, Active(block_minority_prop))
  
  
  #===================================
  # Plot
  #===================================
  
  # city with residents' attribute --------------
  plot_city <- function(time = NULL, show = "ethnicity"){
    if(is.null(time)){
      house_attr <- rep(NA, self$settings$cell_n)
      if(show == "ethnicity"){
        house_attr[self$resident$house] <- self$resident$minority
        house_attr[house_attr==0] <- -1
      }else if(show == "SES"){
        house_attr[self$resident$house] <- self$resident$SES    
      }else{
        stop("'show' must be either 'ethnicity' or 'SES'.")
      }
      city <- array(house_attr, dim = c(self$settings$city_dim, self$settings$city_max_height))
    }else{
      message("plotting the city of time ", time)
      house_attr <- rep(NA, self$settings$cell_n)
      if(show == "ethnicity"){
        house_attr[self$log[[time]]$resident$house] <- self$log[[time]]$resident$minority
        house_attr[house_attr==0] <- -1
      }else if(show == "SES"){
        house_attr[self$log[[time]]$resident$house] <- self$log[[time]]$resident$SES    
      }else{
        stop("'show' must be either 'ethnicity' or 'SES'.")
      }
      city <- array(house_attr, dim = c(self$settings$city_dim, self$settings$city_max_height))
    }
    
    # create an 2D map by calculate the average 
    city_2D <- apply(city, c(1,2), mean, na.rm = TRUE)
    city_2D[is.nan(city_2D)] <- 0
    city_2D_df <- reshape2::melt(city_2D)
    
    # plot
    # Var1 = Row index -> y-axis (scale_y_reverse: from upper to bottom）
    # Var2 = column index -> x-axis（from left to right）
    p <- ggplot(data = city_2D_df) + 
      geom_tile(aes(x = Var2, y = Var1, fill = value)) + 
      scale_fill_gradient2(low  = "#00A0E9",
                           mid  = "white",
                           high = "#E60012",
                           midpoint = 0) + 
      geom_vline(xintercept = seq(1.5, self$settings$city_dim[2], 1),
                 linewidth = .1, linetype = 'dashed', color = 'gray60') + 
      geom_hline(yintercept = seq(1.5, self$settings$city_dim[1], 1),
                 linewidth = .1, linetype = 'dashed', color = 'gray60') +
      geom_vline(xintercept = seq(self$settings$city_lot_dim[2] + 0.5, 
                                  self$settings$city_dim[2], 
                                  self$settings$city_lot_dim[2]),
                 linetype = 'solid', color = 'black') +
      geom_hline(yintercept = seq(self$settings$city_lot_dim[1] + 0.5, 
                                  self$settings$city_dim[1], 
                                  self$settings$city_lot_dim[1]),
                 linetype = 'solid', color = 'black') + 
      labs(x = "", y = "", fill = show) + 
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_reverse(expand = c(0, 0)) + 
      coord_fixed() + 
      theme(
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1, fill = NA)
      )
    
    print(p)
  }
  
  add_field(G, Plot(plot_city))
  
  # Ownership map
  plot_ownership <- function(show = "ethnicity"){
    # reconstruct the city
    city <- array(self$house$landlord, dim = c(self$settings$city_dim, self$settings$city_max_height))
    ownership <- city[,,1]
    ownership_df <- reshape2::melt(ownership)
    colnames(ownership_df)[3] <- "landlord"
    
    # integrate landlords' attribute
    ownership_df <- merge(ownership_df, self$landlord, by.x = "landlord", by.y = "ID", all.x = TRUE)
    if(show == "ethnicity"){
      ownership_df$fill <- ownership_df$minority_aversion
    }else if(show == "SES"){
      ownership_df$fill <- ownership_df$SES_aversion
    }else{
      stop("'show' must be either 'ethnicity' or 'SES'.")
    }
    
    # message
    message("Showing the landlords' aversion on ", show,".")
    
    # plot
     p <- ggplot(data = ownership_df, aes(x = Var2, y = Var1)) + 
      geom_tile(aes(fill = fill)) + 
      geom_text(aes(label = landlord), color = "black") + 
      scale_fill_gradient2(low  = "white",
                           high = "#E60012") + 
      geom_vline(xintercept = seq(1.5, self$settings$city_dim[2], 1),
                 linewidth = .1, linetype = 'dashed', color = 'gray60') + 
      geom_hline(yintercept = seq(1.5, self$settings$city_dim[1], 1),
                 linewidth = .1, linetype = 'dashed', color = 'gray60') +
      geom_vline(xintercept = seq(self$settings$city_lot_dim[2] + 0.5, 
                                  self$settings$city_dim[2], 
                                  self$settings$city_lot_dim[2]),
                 linetype = 'solid', color = 'black') +
      geom_hline(yintercept = seq(self$settings$city_lot_dim[1] + 0.5, 
                                  self$settings$city_dim[1], 
                                  self$settings$city_lot_dim[1]),
                 linetype = 'solid', color = 'black') + 
      labs(x = "", y = "", fill = "aversion") + 
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_reverse(expand = c(0, 0)) +
      coord_fixed() + 
      theme(
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1, fill = NA)
      )
    
    print(p)
  }
  add_field(G, Plot(plot_ownership))
  
  #======================================================
  # Act
  #======================================================
  #---- choose resident randomly ---------------
  select_resident <- function(prop = self$act_defaults$select_resident_prop, 
                              n = NULL){
    if(is.null(n)){
      n <- floor(self$settings$resident_n * prop)
    }
    
    # select agent ID
    id <- rABM::sample2(seq_len(self$settings$resident_n), size = n)
    
    # attach to record
    self$record$selected_resident <- id
  }
  add_field(G, Act(select_resident))
  
  
  #---- resident choose house ------------------
  resident_choose_house <- function(criteria = self$act_defaults$criteria_resident,
                                    beta = self$act_defaults$resident_softmax_beta){
    # validate criteria
    criteria <- match.arg(criteria, 
                          c("intercept", "ethnicity", "SES"), 
                          several.ok = TRUE)
    # selected resident
    resid_ID <- self$record$selected_resident
    
    # specification table for each criterion
    spec <- data.frame(
      criterion  = c("intercept", "ethnicity", "SES"),
      preference = c("preference_intercept", "preference_ethnicity", "preference_SES"),
      house      = c("house_intercept", "house_neib_ethnicity", "house_neib_SES"),
      stringsAsFactors = FALSE
    )
    spec <- spec[spec$criterion %in% criteria, , drop = FALSE]
    
    # retrieve the relevant preference
    resident_preference <- self$resident[resid_ID, spec$preference, drop = FALSE]
    
    # calculate score
    score_list <- lapply(seq_len(nrow(spec)), function(i){
      self[[spec$house[i]]] * resident_preference[, i]
    })
    score <- Reduce(`+`, score_list)
    
    # softmax transformation
    # beta: inverse temperature parameter.
    #   beta -> 0 : near-random selection (bounded rationality)
    #   beta = 1  : standard probabilistic choice
    #   beta -> large : near-deterministic optimal choice (perfect rationality)
    score2 <- rABM::prob_softmax(score, beta = beta)
    
    # mask the houses that are occupied by unmoving residents
    unmoving_resident <- self$resident$ID[!self$resident$ID %in% resid_ID]
    unmoving_house    <- self$resident$house[unmoving_resident]
    score2[, unmoving_house] <- 0
    
    # mask houses with no landlord
    score2[, is.na(self$house$landlord)] <- 0
    
    # weighted sample (sample_weighted re-normalises internally)
    selected_house <- rABM::sample_weighted(score2, size = 1)
    
    # attach to record
    self$record$selected_house <- data.frame(
      resident      = resid_ID,
      current_house = self$resident$house[resid_ID],
      candid_house  = selected_house
    )
  }
  add_field(G, Act(resident_choose_house))
  
  
  # owner decide whether to decline -------
  landlord_decline <- function(criteria = self$act_defaults$criteraia_landlord,
                               block_eth_context = self$act_defaults$block_eth_context){
    # validate criteria
    criteria <- match.arg(criteria, 
                          c("intercept", "ethnicity", "SES"), 
                          several.ok = TRUE)
    
    # specification table for each criterion
    spec <- data.frame(
      criterion = c("intercept", "ethnicity", "SES"),
      landlord_aversion  = c("intercept_aversion", "minority_aversion", "SES_aversion"),
      resident_profile   = c("intercept", "minority", "SES"),
      beta      = c(1, 12, 12/5),
      stringsAsFactors = FALSE
    )
    spec <- spec[spec$criterion %in% criteria, , drop = FALSE]
    
    # resident profile (aligned with selected_house row order via match)
    resident_ID <- self$record$selected_house$resident
    resident_profile <- self$resident[resident_ID, spec$resident_profile, drop = FALSE]
    
    if("SES" %in% criteria){
      resident_profile$SES <- 5 - resident_profile$SES
    }
    
    # landlords' aversion (aligned with selected_house row order)
    relevant_landlord   <- self$house$landlord[self$record$selected_house$candid_house]
    landlord_preference <- self$landlord[relevant_landlord, spec$landlord_aversion, drop = FALSE]
    
    # beta matrix (one beta per criterion, broadcast across rows)
    beta_mat <- matrix(spec$beta, 
                       nrow = nrow(resident_profile), 
                       ncol = length(spec$beta), 
                       byrow = TRUE)
    
    # calculate score
    # NOTE: -6 is the baseline offset to set the baseline declination probability to
    # almost 0% when no intercept is included.
    # When intercept is included (beta = 1, intercept_aversion = 1), a small residual
    # declination probability (~0.7%) remains, representing random "capricious" refusals
    # by landlords regardless of resident attributes.
    house_score <- matrixStats::rowSums2(
      as.matrix(resident_profile * landlord_preference * beta_mat)
    ) - 6
    
    # block context
    if(block_eth_context){
      # cache block_minority_prop
      block_minority_prop <- self$block_minority_prop
      
      # retrieve block id
      selected_block <- self$house$block[self$record$selected_house$candid_house]
      
      # minority prop
      minority_prop <- block_minority_prop[selected_block]
      
      # calculate score
      # Note: if both minority_prop and block_eth_eff is 1, the calculated score 
      # completely cancel out the decline probability. To make score adjust to logit, 
      # setting -12
      block_eth_score <- minority_prop * self$settings$block_eth_eff * (-12)
      
      # renewed score
      house_score <- house_score + block_eth_score
    }
    
    # probability of accepting (1) vs declining (0), per row
    prob     <- 1 / (1 + exp(-house_score))
    prob_mat <- cbind(1 - prob, prob)
    colnames(prob_mat) <- NULL
    decision <- rABM::sample_weighted(prob_mat = prob_mat, size = 1) - 1
    self$record$landlord_decision <- decision
    
    # keep only declined proposals (decision == 0)
    self$record$selected_house <- self$record$selected_house[decision == 0, , drop = FALSE]
    
    # return(invisible)
    invisible(self)
  }
  add_field(G, Act(landlord_decline))
  
  #---- resident move -------------------------
  resident_move <- function(){
    selected_house <- self$record$selected_house
    
    # do nothing if there is no selected house remained
    if(nrow(selected_house)==0){return(NULL)}
    
    # randamize the order
    selected_house <- selected_house[sample(seq_len(nrow(selected_house))), ]
    dup <- duplicated(selected_house$candid_house)
    
    # set duplicated candid house to be the current house
    selected_house$candid_house[dup] <- selected_house$current_house[dup]
    
    # move
    self$resident[selected_house$resident, "house"] <- selected_house$candid_house
  }
  add_field(G, Act(resident_move))
  
  
  resident_move_seq <- function(criteria_resident = self$act_defaults$criteria_resident,
                                criteria_landlord = self$act_defaults$criteraia_landlord,
                                block_eth_context = self$act_defaults$block_eth_context,
                                beta = self$act_defaults$resident_softmax_beta,
                                include_landlord = self$act_defaults$include_landlord_seq) {
    
    # validate criteria
    criteria_resident <- match.arg(criteria_resident,
                                   c("intercept", "ethnicity", "SES"),
                                   several.ok = TRUE)
    criteria_landlord <- match.arg(criteria_landlord,
                                   c("intercept", "ethnicity", "SES"),
                                   several.ok = TRUE)
    
    # specification table for resident house choice
    spec_res <- data.frame(
      criterion  = c("intercept", "ethnicity", "SES"),
      preference = c("preference_intercept", "preference_ethnicity", "preference_SES"),
      house      = c("house_intercept", "house_neib_ethnicity", "house_neib_SES"),
      stringsAsFactors = FALSE
    )
    spec_res <- spec_res[spec_res$criterion %in% criteria_resident, , drop = FALSE]
    
    # specification table for landlord decline decision
    spec_land <- data.frame(
      criterion         = c("intercept", "ethnicity", "SES"),
      landlord_aversion = c("intercept_aversion", "minority_aversion", "SES_aversion"),
      resident_profile  = c("intercept", "minority", "SES"),
      beta_land         = c(1, 12, 12/5),
      stringsAsFactors  = FALSE
    )
    spec_land <- spec_land[spec_land$criterion %in% criteria_landlord, , drop = FALSE]
    
    # get selected residents and randomize order
    resid_IDs <- self$record$selected_resident
    resid_IDs <- rABM::sample2(resid_IDs, size = length(resid_IDs))
    
    # initialize record vectors for this step
    record_resident <- integer(length(resid_IDs))
    record_current  <- integer(length(resid_IDs))
    record_candid   <- integer(length(resid_IDs))
    record_decision <- integer(length(resid_IDs))
    
    # sequential loop over each resident
    for(i in seq_along(resid_IDs)){
      resid_ID <- resid_IDs[i]
      
      #----------------------------------------------
      # Step 1: resident chooses a candidate house
      #----------------------------------------------
      
      # temporarily set selected_resident to this single agent
      # Active states (house_neib_ethnicity / house_neib_SES) use this internally
      # and are recomputed reflecting the current positions of all residents
      self$record$selected_resident <- resid_ID
      
      # retrieve preference for this resident (1-row data.frame)
      resident_preference <- self$resident[resid_ID, spec_res$preference, drop = FALSE]
      
      # calculate score (Active states recomputed here, reflecting current positions)
      score_list <- lapply(seq_len(nrow(spec_res)), function(k){
        self[[spec_res$house[k]]] * resident_preference[, k]
      })
      score <- Reduce(`+`, score_list)
      # score is a 1 x cell_n matrix here
      
      # softmax transformation
      # beta: inverse temperature parameter.
      #   beta -> 0 : near-random selection (bounded rationality)
      #   beta = 1  : standard probabilistic choice
      #   beta -> large : near-deterministic optimal choice (perfect rationality)
      score2 <- rABM::prob_softmax(score, beta = beta)
      
      # mask occupied houses (all residents except the current agent)
      other_residents <- self$resident$ID[self$resident$ID != resid_ID]
      occupied_house  <- self$resident$house[other_residents]
      score2[, occupied_house] <- 0
      
      # mask houses with no landlord
      score2[, is.na(self$house$landlord)] <- 0
      
      # weighted sample (sample_weighted re-normalises internally)
      candid_house  <- rABM::sample_weighted(score2, size = 1)
      current_house <- self$resident$house[resid_ID]
      
      # go NEXT if current_house is selected
      if(candid_house == current_house) next
      
      #----------------------------------------------
      # Step 2: landlord decides whether to decline
      # (skipped if include_landlord = FALSE: resident always moves)
      #----------------------------------------------
      
      if(isTRUE(include_landlord)){
        
        resident_profile <- self$resident[resid_ID, spec_land$resident_profile, drop = FALSE]
        
        if("SES" %in% criteria_landlord){
          resident_profile$SES <- 5 - resident_profile$SES
        }
        
        relevant_landlord   <- self$house$landlord[candid_house]
        landlord_preference <- self$landlord[relevant_landlord, spec_land$landlord_aversion, drop = FALSE]
        
        # calculate logit score
        # NOTE: -6 is the baseline offset to set the baseline declination probability to
        # almost 0% when no intercept is included.
        # When intercept is included (beta_land = 1, intercept_aversion = 1), a small
        # residual declination probability (~0.7%) remains, representing random
        # "capricious" refusals by landlords regardless of resident attributes.
        house_score <- sum(
          as.numeric(resident_profile) * as.numeric(landlord_preference) * spec_land$beta_land
        ) - 6
        
        # block ethnic context
        if(isTRUE(block_eth_context)){
          block_minority_prop <- self$block_minority_prop
          selected_block      <- self$house$block[candid_house]
          minority_prop       <- block_minority_prop[selected_block]
          block_eth_score     <- minority_prop * self$settings$block_eth_eff * (-12)
          house_score         <- house_score + block_eth_score
        }
        
        # probability of declining
        prob     <- 1 / (1 + exp(-house_score))
        prob_mat <- matrix(c(1 - prob, prob), nrow = 1)
        decision <- rABM::sample_weighted(prob_mat, size = 1) - 1
        # decision: 1 = declined, 0 = accepted
        
      }else{
        # landlord screening skipped: resident always accepted
        decision <- 0L
      }
      
      #----------------------------------------------
      # Step 3: move if accepted
      #----------------------------------------------
      
      if(decision == 0){
        self$resident[resid_ID, "house"] <- candid_house
      }
      
      # record
      record_resident[i] <- resid_ID
      record_current[i]  <- current_house
      record_candid[i]   <- candid_house
      record_decision[i] <- decision
    }
    
    # restore selected_resident to the full (randomized) set
    self$record$selected_resident <- resid_IDs
    
    # store results in record
    self$record$selected_house <- data.frame(
      resident      = record_resident,
      current_house = record_current,
      candid_house  = record_candid
    )
    self$record$landlord_decision <- record_decision
    
    invisible(self)
  }
  add_field(G, Act(resident_move_seq))
  
  #================================================
  # Report
  #================================================
  
  # report segregation ---------------------------
  report_segregation <- function(level = c("city", "landlord"), 
                                 show = c("ethnicity", "SES"),
                                 index = c("D","M","H"),
                                 log = NULL,
                                 plot = TRUE){
    # match arg
    level <- match.arg(level)
    show  <- match.arg(show)
    index <- match.arg(index)
    
    if(show == "ethnicity"){
      variable <- "minority"
    }else{
      variable <- "SES"
    }
    
    if(level == "city"){
      obs <- "block"
    }else{
      obs <- "landlord"
    }
    
    # From the current observation
    if(is.null(log)){
      house_variable <- rep(NA, length = self$settings$cell_n)
      house_variable[self$resident$house] <- self$resident[[variable]]
      
      tb_df <- as.data.frame(table(self$house[[obs]], house_variable))
      ind <- switch(index,
                    "D" = {segregation::dissimilarity(tb_df, unit = "Var1", group = "house_variable", weight = "Freq")$est},
                    "M" = {segregation::mutual_total(tb_df, unit = "Var1", group = "house_variable", weight = "Freq")$est[1]},
                    "H" = {segregation::mutual_total(tb_df, unit = "Var1", group = "house_variable", weight = "Freq")$est[2]}
      )
      
      cat(paste0("Reporting ", index, "-index at '", level, "' level on '", show, "'.", "\n"))
      return(ind)
    }else{
      # From log --------------
      house_list <- rABM::value_of(G = self, field_name = "house", log = log)
      resident_list <- rABM::value_of(G = self, field_name = "resident", log = log)
      time_list <- unlist(rABM::value_of(G = self, field_name = "time", log = log))
      
      ind <- numeric(length = length(house_list))
      for(t in seq_len(length(house_list))){
        house_variable <- rep(NA, length=self$settings$cell_n)
        house_variable[resident_list[[t]]$house] <- resident_list[[t]][[variable]]
        obs_variable <- house_list[[t]][[obs]]
        tb_df <- as.data.frame(table(obs_variable, house_variable))
        
        ind_t <- switch(index,
                        "D" = {segregation::dissimilarity(tb_df, unit = "obs_variable", group = "house_variable", weight = "Freq")$est},
                        "M" = {segregation::mutual_total(tb_df, unit = "obs_variable", group = "house_variable", weight = "Freq")$est[1]},
                        "H" = {segregation::mutual_total(tb_df, unit = "obs_variable", group = "house_variable", weight = "Freq")$est[2]}
        )
        ind[t] <- ind_t
      }
      
      # plot
      if(isTRUE(plot)){
        plot_df <- data.frame(time = time_list, ind = ind)
        p <- ggplot(data = plot_df) + geom_line(aes(x = time, y = ind), linewidth = 1) + 
          theme_bw() + xlab("time") + ylab(paste0(index, "-index")) + ylim(c(0,1))
        print(p)
      }
      
      # return
      cat(paste0("Reporting ", index, "-index at '", level, "' level on '", show, "'.", "\n"))
      names(ind) <- as.numeric(time_list)
      return(ind)
    }
  }
  
  add_field(G, Report(report_segregation))
  
  
  # report landlord's descriptive statistics
  report_landlord_stat <- function(log = NULL){
    if(is.null(log)){
      # house
      house <- self$house
      house <- data.frame(house, minority = NA, SES = NA)
      house[self$resident$house, c("minority","SES")] <- self$resident[ ,c("minority","SES")]
      
      # landlord
      landlord <- self$landlord
      out <- data.frame(landlord,
                        vacancy = 1- tapply(!is.na(house$minority), house$landlord, mean),
                        minority_prop = tapply(house$minority, house$landlord, mean, na.rm = TRUE),
                        SES_mean = tapply(house$SES, house$landlord, mean, na.rm = TRUE),
                        SES_sd   = tapply(house$SES, house$landlord, sd, na.rm = TRUE))
    }else{
      house_list <- rABM::value_of(G = self, field_name = "house", log = log)
      resident_list <- rABM::value_of(G = self, field_name = "resident", log = log)
      landlord_list <- rABM::value_of(G = self, field_name = "landlord", log = log)
      time_list <- unlist(rABM::value_of(G = self, field_name = "time", log = log))
      out_list <- vector("list", length = length(time_list))
      for(t in seq_len(length(out_list))){
        house <- house_list[[t]]
        house <- data.frame(house, minority = NA, SES = NA)
        house[resident_list[[t]]$house, c("minority", "SES")] <- resident_list[[t]][ ,c("minority", "SES")]
        
        # landlord
        landlord <- landlord_list[[t]]
        out_list[[t]] <- data.frame(time = time_list[[t]],
                                    landlord, 
                                    vacancy = 1 - tapply(!is.na(house$minority), house$landlord, mean),
                                    minority_prop = tapply(house$minority, house$landlord, mean, na.rm = TRUE),
                                    SES_mean = tapply(house$SES, house$landlord, mean, na.rm = TRUE),
                                    SES_sd   = tapply(house$SES, house$landlord, sd, na.rm = TRUE))
      }
      out <- do.call(rbind, out_list)
    }
    
    out
  }
  
  add_field(G, Report(report_landlord_stat))
  
  
  #================================================
  # output G
  #================================================
  
  G
}

