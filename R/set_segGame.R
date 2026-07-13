#' Set up a Segregation ABM Game object
#'
#' Constructs a fully configured \code{Game} object for simulating residential
#' segregation dynamics. The model includes residents, landlords, and a
#' gridded city structure, and incorporates ethnic and socioeconomic
#' preferences of residents, discriminatory behaviour of landlords, and
#' landlord ownership turnover over time.
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
#' @param resident_ownership_prop Numeric scalar in \code{[0, 1]}. Target
#'   proportion of residents who own their home
#'   (\code{resident$house_owning = 1}) versus rent (\code{= 0}), assigned
#'   independently of ethnicity and SES via \code{make_group_labels()}. Used
#'   by \code{select_resident()} when \code{separate_home_owners = TRUE} to
#'   stratify resident selection by ownership status. Default: \code{0.5}.
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
#' @param block_eth_threshold Numeric scalar. Minimum block-level minority
#'   proportion required for the block ethnic context effect to be applied.
#'   When \code{block_eth_context = TRUE} in \code{resident_move_seq()}, the
#'   minority proportion of the block containing the candidate house
#'   (\code{self$block_minority_prop[selected_block]}) is compared against
#'   this threshold; the penalty
#'   (\code{minority_prop * settings$block_eth_eff * (-12)}) is added to the
#'   landlord's decline score only when the block's minority proportion is at
#'   or above \code{block_eth_threshold}. Below the threshold, block-level
#'   context has no effect and only individual-level criteria
#'   (\code{criteria_landlord}) enter the decline score. Default \code{0}
#'   effectively applies the block effect unconditionally (matching prior
#'   behaviour). Note: this threshold is currently only applied in
#'   \code{resident_move_seq()} (sequential mode); \code{landlord_decline()}
#'   (batch mode) applies the block effect without a threshold gate whenever
#'   \code{block_eth_context = TRUE}. Stored in
#'   \code{G$act_defaults$block_eth_threshold} and can be modified after
#'   construction, e.g. \code{G$act_defaults$block_eth_threshold <- 0.3}.
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
#' @param default_criteria_landlord Character vector. Default criteria used
#'   by \code{landlord_decline()} and \code{resident_move_seq()} to determine
#'   which components enter the landlord decline score. Any subset of
#'   \code{c("intercept", "ethnicity", "SES")}. Default: all three.
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
#' @param sort_by_SES Logical. Default value of \code{sort_by_SES} in
#'   \code{select_resident()}, stored in \code{G$act_defaults}. If
#'   \code{TRUE}, the residents drawn for a given step are additionally
#'   re-ordered in descending order of SES (highest SES first) before being
#'   returned. This only affects processing order and is most relevant to
#'   \code{resident_move_seq()}, where residents are handled one at a time and
#'   later residents observe the updated positions of earlier ones. Default:
#'   \code{FALSE}.
#' @param convergence_tol Positive numeric scalar. Convergence tolerance used
#'   by the built-in \code{Stop} function \code{converged()}. Convergence is
#'   declared once the range (maximum \eqn{-} minimum) of the most recent
#'   \code{convergence_times} values of the city-level ethnic dissimilarity
#'   index (D) falls below this tolerance. Default: \code{0.05}.
#' @param convergence_times Positive integer. Number of most recent
#'   dissimilarity-index (D) values retained and checked by
#'   \code{converged()} when evaluating the convergence criterion. Default:
#'   \code{10}.
#' @param convergence_maxit Positive integer. Maximum number of time steps.
#'   \code{converged()} returns \code{TRUE} once \code{self$time} reaches this
#'   value, regardless of whether the tolerance-based criterion has been met,
#'   ensuring the simulation always terminates. Default: \code{500}.
#' @param convergence_burnin Non-negative integer. Number of initial time
#'   steps during which \code{converged()} always returns \code{FALSE},
#'   allowing the system to move away from its (typically random) initial
#'   configuration before convergence is evaluated. Default: \code{100}.
#' @param separate_home_owners Logical. Default value of
#'   \code{separate_home_owners} in \code{select_resident()}, stored in
#'   \code{G$act_defaults}. If \code{TRUE}, the \code{n} residents selected
#'   each step are drawn separately from home owners
#'   (\code{resident$house_owning == 1}) and non-owners
#'   (\code{house_owning == 0}) in proportions controlled by
#'   \code{home_owner_selection_prop}, instead of being drawn as a single
#'   uniform sample from all residents. Default: \code{FALSE}.
#' @param home_owner_selection_prop Numeric scalar in \code{[0, 1]}. When
#'   \code{separate_home_owners = TRUE} in \code{select_resident()}, the
#'   proportion of the \code{n} selected residents per step drawn from home
#'   owners; the remainder (\code{n - round(home_owner_selection_prop * n)})
#'   is drawn from non-owners. Has no effect when
#'   \code{separate_home_owners = FALSE}. Default: \code{0.25}.
#' @param landlord_change_para Positive numeric scalar. Governs the landlord
#'   ownership-turnover hazard used by \code{landlord_change()}. For each
#'   house, the cumulative probability that its landlord changes is modelled
#'   as a Rayleigh distribution (a Weibull distribution with \code{shape = 2})
#'   over \code{house$landlord_duration}, implemented as
#'   \code{pweibull(landlord_duration, shape = 2, scale =
#'   landlord_change_para * sqrt(2))}, with \code{landlord_change_para} acting
#'   as the Rayleigh scale parameter \eqn{\sigma}. Larger values make
#'   ownership more stable (the turnover probability rises more slowly with
#'   duration); as a rule of thumb, the cumulative change probability reaches
#'   \eqn{\approx} 39\% once \code{landlord_duration} reaches
#'   \code{landlord_change_para}, and \eqn{\approx} 86\% once it reaches
#'   \code{2 * landlord_change_para}. Default: \code{30} (time steps).
#'
#' @return A \code{Game} object (R6 class) configured with the following
#'   fields:
#'   \describe{
#'     \item{States}{
#'       \code{landlord}, \code{resident}, \code{resident_mat_ethnicity},
#'       \code{resident_mat_SES}, \code{house}, \code{house_neighbor_ind},
#'       \code{house_neighbor_ind_0idx} (0-indexed, NA-as-\code{-1} version of
#'       \code{house_neighbor_ind} used internally by the Rcpp-backed
#'       neighbourhood-proportion helpers), \code{settings}, \code{record},
#'       \code{act_defaults}, \code{convergence}
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
#'     \item{Act functions (landlord turnover)}{
#'       \code{landlord_change}
#'     }
#'     \item{Plot functions}{
#'       \code{plot_city}, \code{plot_ownership}
#'     }
#'     \item{Report functions}{
#'       \code{report_segregation}, \code{report_landlord_stat},
#'       \code{report_result_df}
#'     }
#'     \item{Stop functions}{
#'       \code{converged}
#'     }
#'   }
#'
#'   The \code{resident} state additionally includes a \code{house_owning}
#'   column (\code{1} = owner, \code{0} = renter; see
#'   \code{resident_ownership_prop}), and the \code{house} state includes a
#'   \code{landlord_duration} column (number of consecutive time steps the
#'   current landlord has held that house; see \code{landlord_change_para}).
#'
#'   The \code{settings} state contains the following named elements:
#'   \code{resident_n}, \code{landlord_n}, \code{house_n} (number of
#'   habitable units, i.e. cells with a landlord), \code{cell_n},
#'   \code{city_dim}, \code{city_zone_dim}, \code{city_lot_dim},
#'   \code{city_max_height}, \code{block_eth_eff}, \code{block_SES_eff},
#'   \code{convergence_tol}, \code{convergence_times},
#'   \code{convergence_maxit}, \code{convergence_D_prev} (initialised to
#'   \code{0}), and \code{convergence_burnin}.
#'
#'   The \code{act_defaults} state stores the eleven default parameter values
#'   listed above (\code{select_resident_prop}, \code{criteria_resident},
#'   \code{resident_softmax_beta}, \code{criteria_landlord},
#'   \code{block_eth_context}, \code{include_landlord_seq},
#'   \code{sort_by_SES}, \code{block_eth_threshold},
#'   \code{separate_home_owners}, \code{home_owner_selection_prop},
#'   \code{landlord_change_para}), and can be updated at run time to change
#'   the behaviour of Act functions without re-initialising the \code{Game}
#'   object.
#'
#'   The \code{convergence} state stores \code{prev_D}, the trailing history
#'   (up to \code{convergence_times} values) of the city-level dissimilarity
#'   index (D) used internally by \code{converged()}.
#'
#' @details
#' \subsection{Model overview}{
#'   The model simulates residential sorting driven by two mechanisms:
#'   (1) residents' preferences for co-ethnic and same-SES neighbours, and
#'   (2) landlords' discriminatory screening of applicants by ethnicity and/or
#'   SES. Two update modes are provided: a batch mode and a sequential mode.
#'   Independently of these, landlord ownership itself can turn over over time
#'   via \code{landlord_change()} (see below).
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
#'   Selected residents are processed one at a time in random order (unless
#'   \code{sort_by_SES = TRUE}, see below). For each resident, the full cycle
#'   of (1) house choice, (2) landlord screening, and (3) movement is
#'   completed before the next resident is processed. Neighbourhood
#'   composition scores (Active states) are recomputed at each step, so later
#'   residents see the updated positions of earlier movers. Setting
#'   \code{include_landlord = FALSE} (or via
#'   \code{default_include_landlord_seq}) disables landlord screening
#'   entirely, so every resident always moves to their chosen unit. This is
#'   useful for isolating the effect of residential preferences alone.
#'   Both modes are theoretically equivalent when \code{include_landlord =
#'   TRUE} and \code{block_eth_threshold = 0}, and can be used to verify
#'   convergence properties; see \code{block_eth_threshold} above for a case
#'   where the two modes currently diverge.
#' }
#'
#' \subsection{Selecting residents and processing order}{
#'   \code{select_resident()} draws \code{n} (or \code{prop} of) residents at
#'   random for the current step and stores their IDs in
#'   \code{self$record$selected_resident}.
#'
#'   When \code{separate_home_owners = TRUE} (default stored in
#'   \code{act_defaults$separate_home_owners}), the draw is stratified by
#'   \code{resident$house_owning}: a proportion \code{home_owner_selection_prop}
#'   of the \code{n} selected residents is drawn from home owners and the
#'   remainder from non-owners, then the combined set is shuffled. This is
#'   useful for scenarios where owners and renters should be represented in a
#'   fixed ratio each step regardless of their relative population sizes
#'   (e.g. to ensure owners — who may move less often in reality — are still
#'   adequately sampled). When \code{FALSE} (default), residents are drawn as
#'   a single uniform sample, ignoring ownership status.
#'
#'   When \code{sort_by_SES = TRUE} (default stored in
#'   \code{act_defaults$sort_by_SES}), the drawn residents are additionally
#'   re-ordered in descending order of SES before being stored. This does not
#'   change \emph{which} residents are selected, only the order in which they
#'   are subsequently processed. The ordering has no effect under batch mode
#'   (all selected residents are evaluated simultaneously), but it does affect
#'   \code{resident_move_seq()}, where residents earlier in the queue move
#'   (and thereby reshape neighbourhood composition) before residents later in
#'   the queue make their choice.
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
#'
#'   Internally, \code{house_neib_ethnicity} and \code{house_neib_SES} compute
#'   neighbourhood composition proportions via the Rcpp helpers
#'   \code{neib_ethnicity_prop()} / \code{neib_SES_prop()} over
#'   \code{house_neighbor_ind_0idx}, then combine them with resident-level
#'   ethnicity/SES matrices via ordinary (BLAS-backed) matrix multiplication.
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
#'   When \code{block_eth_context = TRUE}, a further penalty based on the
#'   block-level minority proportion is added to the score; see
#'   \code{block_eth_threshold} above for how this differs between the batch
#'   and sequential update modes.
#'
#'   Landlord screening can be disabled entirely in the sequential mode by
#'   setting \code{include_landlord = FALSE} in \code{resident_move_seq()},
#'   or by changing \code{G$act_defaults$include_landlord_seq} at run time.
#' }
#'
#' \subsection{Landlord ownership turnover}{
#'   Independently of resident movement, \code{landlord_change()} models
#'   landlords selling/losing their properties to a new, randomly drawn
#'   landlord (\code{rABM::sample2()} over \code{seq_len(landlord_n)}, with
#'   replacement, so the incoming landlord may coincide with the outgoing
#'   one). For each house, the probability that a change occurs this step is
#'   \code{pweibull(house$landlord_duration, shape = 2, scale =
#'   landlord_change_para * sqrt(2))} — a Rayleigh hazard that starts at 0 for
#'   a freshly-acquired house and increases with tenure (see
#'   \code{landlord_change_para} above). Houses that change owners have
#'   \code{landlord_duration} reset to \code{1}; unchanged houses have it
#'   incremented by \code{1}. This Act function is independent of resident
#'   preferences and of \code{landlord_decline()}/\code{resident_move_seq()}'s
#'   discrimination model — it does not affect who lives where, only who owns
#'   which house and (via \code{landlord$minority_aversion}/\code{SES_aversion})
#'   with what discriminatory disposition. It is not included in either
#'   default \code{plan}; add \code{"landlord_change"} to the \code{plan}
#'   argument of \code{run_Game()} to enable turnover.
#' }
#'
#' \subsection{Convergence criterion}{
#'   The \code{Game} object is equipped with a \code{Stop} function,
#'   \code{converged()}, that can be used to terminate a simulation
#'   automatically (e.g. when \code{times} is left unspecified or set to a
#'   large upper bound in \code{run_Game()}). At each time step, it computes
#'   the city-level ethnic dissimilarity index (D) between minority and
#'   majority residents across blocks (\code{self$house$block}) and appends
#'   it to a trailing history capped at \code{convergence_times} values
#'   (stored in \code{self$convergence$prev_D}).
#'
#'   \code{converged()} returns \code{FALSE} while \code{self$time} is less
#'   than or equal to \code{convergence_burnin}, allowing the system to move
#'   away from its initial (typically near-random) configuration. After the
#'   burn-in period, it returns \code{TRUE} as soon as either of the following
#'   holds:
#'   \enumerate{
#'     \item The D-index history contains at least \code{convergence_times}
#'       values and their range (maximum \eqn{-} minimum) is smaller than
#'       \code{convergence_tol}, indicating the segregation level has
#'       stabilised; or
#'     \item \code{self$time} has reached \code{convergence_maxit}, which
#'       guarantees the simulation terminates even if the tolerance-based
#'       criterion is never satisfied.
#'   }
#'   The relevant defaults (\code{convergence_tol}, \code{convergence_times},
#'   \code{convergence_maxit}, \code{convergence_burnin}) are stored in
#'   \code{G$settings} and can be modified after construction, e.g.
#'   \code{G$settings$convergence_tol <- 0.02}.
#' }
#'
#' \subsection{Plot functions in detail}{
#'   Both plot functions read directly from \code{self} (i.e. the live
#'   \code{Game} object's current state, or a single saved time point from
#'   \code{self$log}), and are meant for quick visual diagnostics during
#'   interactive exploration rather than for extracting numbers. Neither
#'   function returns a usable value (both simply \code{print()} a
#'   \pkg{ggplot2} object and return it invisibly via \code{print}); to get
#'   the underlying \code{ggplot} object for further customisation (e.g.
#'   adding \code{ggtitle()}, changing the colour scale, or faceting), copy
#'   the function body out of \code{set_segGame()} and adapt it, or simply
#'   reconstruct the input data frame (\code{city_2D_df} / \code{ownership_df})
#'   yourself and call \code{ggplot()} on it directly.
#'
#'   \describe{
#'     \item{\code{plot_city(time = NULL, show = "ethnicity")}}{
#'       Draws a top-down map of the city grid, one tile per lot.
#'       \itemize{
#'         \item \strong{Data source (\code{time})}: if \code{time = NULL}
#'           (default), the plot reflects the \emph{current} state of
#'           \code{self$resident} (i.e. wherever the simulation currently
#'           stands). If \code{time} is set to an integer, the plot instead
#'           reads \code{self$log[[time]]$resident}, i.e. a snapshot saved by
#'           \code{run_Game()}. This only works if \code{"resident"} was
#'           included in \code{fields_to_save} when the simulation was run
#'           (the default \code{G$notes$fields_to_save} already includes it),
#'           and if \code{times} in \code{run_Game()} was left short enough,
#'           or \code{log} indices requested, that the desired time point was
#'           actually retained. This lets a researcher compare snapshots at
#'           different points of the same run without re-simulating, e.g.
#'           \code{G$plot_city(time = 1)} vs. \code{G$plot_city(time = 50)}.
#'         \item \strong{Transformation (\code{show})}: with
#'           \code{show = "ethnicity"} (default), each occupied lot is coded
#'           \code{-1} (majority) or \code{+1} (minority); vacant/unowned
#'           cells are \code{NA}. With \code{show = "SES"}, each occupied lot
#'           is coded by the resident's raw SES value (\code{0}-\code{5}).
#'           Because a lot can contain multiple floors
#'           (\code{city_max_height} > 1), the function first builds a 3-D
#'           array (rows \eqn{\times} columns \eqn{\times} floors) and then
#'           collapses it to 2-D by taking the \strong{mean across floors}
#'           (\code{NA}s removed), so what is plotted is, technically, the
#'           average ethnicity code or average SES of the floors within each
#'           lot, not a single resident's value. \code{NaN} results (i.e. a
#'           lot entirely vacant/unowned) are recoded to \code{0} before
#'           plotting.
#'         \item \strong{Colour scale}: a diverging red-white-blue gradient
#'           centred at \code{0} (\code{scale_fill_gradient2()}). For
#'           \code{show = "ethnicity"} this cleanly separates
#'           majority-dominated lots (blue, negative) from
#'           minority-dominated lots (red, positive), with mixed/vacant lots
#'           near white. For \code{show = "SES"} the midpoint of \code{0}
#'           does not correspond to a meaningful "neutral" SES value (SES
#'           ranges \code{0}-\code{5}), so the colour scale is less
#'           informative for that case and researchers who want a proper SES
#'           gradient may prefer to rebuild the plot with
#'           \code{scale_fill_gradient()} or \code{scale_fill_viridis_c()}
#'           instead.
#'         \item \strong{Overlaid grid lines}: dashed grey lines mark
#'           individual lots; solid black lines mark zone boundaries, derived
#'           from \code{self$settings$city_lot_dim} and
#'           \code{self$settings$city_dim}. These will automatically adapt if
#'           \code{city_zone_dim} / \code{city_lot_dim} are changed at
#'           construction time.
#'       }
#'     }
#'     \item{\code{plot_ownership(show = "ethnicity")}}{
#'       Draws a map of landlord ownership footprints, with each tile
#'       labelled by its owning landlord's ID.
#'       \itemize{
#'         \item Unlike \code{plot_city()}, this always uses only the
#'           \strong{ground floor} (\code{city[,,1]}) as a proxy for
#'           ownership across the lot, on the assumption that ownership does
#'           not vary by floor; there is no \code{time} argument, so it
#'           always reflects the \emph{current} \code{self$house$landlord}
#'           assignment (there is no way to inspect ownership at a past log
#'           point through this function as written, and — since
#'           \code{landlord_change()} mutates ownership over time if included
#'           in the \code{plan} — this plot reflects turnover as well as the
#'           original geographic assignment).
#'         \item \strong{Transformation (\code{show})}: with
#'           \code{show = "ethnicity"} (default) each tile is coloured by
#'           that landlord's \code{minority_aversion} value; with
#'           \code{show = "SES"} it is coloured by \code{SES_aversion}. These
#'           are the same aversion values used inside \code{landlord_decline()}
#'           / \code{resident_move_seq()}, so this plot is a direct way to
#'           see \emph{where} the more discriminatory landlords' portfolios
#'           are concentrated in space (relevant to
#'           \code{landlord_dispersed_prop} and \code{landlord_geo_theta}).
#'         \item \strong{Colour scale}: a sequential white-to-red gradient
#'           (\code{scale_fill_gradient2(low = "white", high = "#E60012")}),
#'           i.e. darker red = more hostile landlord on the selected
#'           dimension. If \code{landlord_hostile_input_type = "binary"} was
#'           used at construction, this collapses to essentially two shades
#'           (aversion \code{0} or \code{1}); with
#'           \code{landlord_hostile_input_type = "numeric"} it shows a
#'           continuous gradient of aversion.
#'       }
#'     }
#'   }
#' }
#'
#' \subsection{Report functions in detail}{
#'   The three report functions differ from the plot functions in that they
#'   return ordinary R objects (numeric scalars/vectors or data frames),
#'   making them the natural entry point for exporting simulation output to
#'   external statistical software or for building custom summaries/plots not
#'   covered by \code{plot_city()} / \code{plot_ownership()}.
#'
#'   \describe{
#'     \item{\code{report_segregation(level, show, index, log, plot)}}{
#'       Computes an aggregate segregation index from the current
#'       (\code{log = NULL}) or logged (\code{log = <vector of time indices>})
#'       state of the city.
#'       \itemize{
#'         \item \strong{\code{level}}: \code{"city"} treats each
#'           \code{self$house$block} as the unit across which segregation is
#'           measured (i.e. how unevenly residents are sorted across
#'           neighbourhoods/blocks); \code{"landlord"} instead treats each
#'           landlord's portfolio as the unit (i.e. how unevenly residents
#'           are sorted across landlords' holdings, which is closer to a
#'           measure of \emph{landlord-driven} sorting).
#'         \item \strong{\code{show}}: \code{"ethnicity"} uses
#'           \code{resident$minority} as the group variable;
#'           \code{"SES"} uses \code{resident$SES} (all 6 SES categories, not
#'           collapsed).
#'         \item \strong{\code{index}}: which segregation statistic to
#'           compute via package \pkg{segregation}. \code{"D"} is the
#'           dissimilarity index (\eqn{0}-\eqn{1}; share of one group that
#'           would need to relocate to achieve even distribution — this is
#'           the same statistic used internally by \code{converged()}, but
#'           computed here on demand and configurable by \code{level}/
#'           \code{show}). \code{"M"} is Theil's information-theoretic mutual
#'           information index (\code{mutual_total()$est[1]}, unbounded
#'           above, sensitive to the number of groups/units). \code{"H"} is
#'           the associated normalised entropy index
#'           (\code{mutual_total()$est[2]}, bounded in \eqn{[0, 1]}, and
#'           generally the more comparable choice to \code{"D"} across
#'           differently-sized runs).
#'         \item \strong{\code{log}}: if \code{NULL} (default), a single
#'           numeric value for the \emph{current} state is returned. If set
#'           to a vector of saved time indices (e.g. \code{log = 1:50}), the
#'           function instead loops over \code{self$log} at those indices via
#'           \code{rABM::value_of()}, and returns a named numeric vector (one
#'           value per requested time point, named by the actual logged
#'           \code{time} value) — this requires that \code{"house"} and
#'           \code{"resident"} were both included in \code{fields_to_save}
#'           when the run was executed. This is the natural way to trace how
#'           segregation evolves over the course of a run, e.g. to visually
#'           confirm convergence or to compare trajectories across parameter
#'           settings.
#'         \item \strong{\code{plot}}: only takes effect when \code{log} is
#'           supplied; if \code{TRUE} (default) a \pkg{ggplot2} line plot of
#'           the index against time is drawn (y-axis fixed to \code{[0, 1]}
#'           via \code{ylim()} — note this clipping is only strictly correct
#'           for \code{index = "D"} or \code{"H"}; if you request
#'           \code{index = "M"} and its values exceed 1, the plot will
#'           silently clip them, so it is worth checking the returned numeric
#'           vector directly in that case rather than relying on the plot
#'           alone). Set \code{plot = FALSE} to suppress the figure and just
#'           get the numbers back, e.g. for stitching multiple runs together
#'           into a single custom \pkg{ggplot2}/\pkg{data.table} pipeline.
#'       }
#'     }
#'     \item{\code{report_landlord_stat(log = NULL)}}{
#'       Returns landlord-level descriptive statistics as a data frame with
#'       one row per landlord (or, when \code{log} is supplied, one row per
#'       landlord per requested time point, with an added \code{time}
#'       column, stacked via \code{rbind}).
#'       \itemize{
#'         \item Columns beyond the landlord's own attributes
#'           (\code{minority_aversion}, \code{SES_aversion}, etc.) are:
#'           \code{vacancy} (proportion of that landlord's units currently
#'           unoccupied), \code{minority_prop} (share of that landlord's
#'           tenants who are minority residents, \code{NA} if the landlord
#'           currently has no tenants), \code{SES_mean} and \code{SES_sd}
#'           (mean/sd of tenant SES for that landlord).
#'         \item \strong{\code{log}}: as with \code{report_segregation()},
#'           \code{NULL} gives the current snapshot; a vector of saved time
#'           indices gives a stacked panel data frame across time, suitable
#'           directly for a fixed-effects or random-effects panel regression
#'           (e.g. \code{plm}) of tenant composition on landlord aversion
#'           traits or on \code{landlord_dispersed}/geographic variables,
#'           analogous to the kind of panel set-up already used elsewhere in
#'           this workflow. Requires \code{"house"}, \code{"resident"}, and
#'           \code{"landlord"} to have been included in
#'           \code{fields_to_save}.
#'       }
#'     }
#'     \item{\code{report_result_df(time = NULL)}}{
#'       Returns a single "flat" resident-level data frame — one row per
#'       resident — with the resident's own attributes (including
#'       \code{resident_owning}, i.e. \code{house_owning}) plus the attributes
#'       of their currently-assigned house/block (\code{house_ID},
#'       \code{house_block}, \code{landlord_duration}) and landlord
#'       (\code{landlord_ID}, \code{landlord_minority_aversion},
#'       \code{landlord_SES_aversion}, \code{landlord_dispersed}) already
#'       merged in.
#'       \itemize{
#'         \item \strong{\code{time}}: if \code{NULL} (default), uses the
#'           live/current \code{self}. If set to an integer, temporarily
#'           substitutes \code{self <- self$log[[time]]} before assembling
#'           the data frame, i.e. it reconstructs the identical
#'           multi-level (resident-in-house-in-landlord) data frame for any
#'           single saved time point. Passing a \code{time} larger than
#'           \code{length(self$log)} raises an explicit error rather than
#'           returning \code{NULL}/an empty frame.
#'         \item This is the function best suited for exporting a single
#'           cross-section (or, by looping over several values of
#'           \code{time} and \code{rbind()}-ing the results, a full panel) to
#'           multilevel/hierarchical models (e.g. residents nested in
#'           landlords, via \code{lme4::lmer()} or \code{lavaan}), since it
#'           already carries the nesting structure (resident \eqn{\to} house
#'           \eqn{\to} landlord) as merged columns rather than requiring the
#'           researcher to join \code{self$resident}, \code{self$house}, and
#'           \code{self$landlord} by hand.
#'       }
#'     }
#'   }
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
#' # --- Sequential mode, processing residents in descending SES order ---
#' G$act_defaults$sort_by_SES <- TRUE
#'
#' # --- Sequential mode with landlord ownership turnover each step ---
#' G_turnover <- set_segGame(landlord_change_para = 20)
#' G_turnover_out <- run_Game(
#'   G_turnover,
#'   plan = c("select_resident", "resident_move_seq", "landlord_change"),
#'   times = 50,
#'   fields_to_save = G_turnover$notes$fields_to_save
#' )
#'
#' # --- Stratified resident selection by home-ownership status ---
#' G_owners <- set_segGame(separate_home_owners = TRUE,
#'                         home_owner_selection_prop = 0.5)
#'
#' # --- Run until convergence (or convergence_maxit), using the built-in
#' #     Stop function instead of a fixed number of steps ---
#' G_conv <- set_segGame(convergence_tol = 0.02, convergence_burnin = 50)
#' G_conv_out <- run_Game(
#'   G_conv,
#'   plan = c("select_resident", "resident_move_seq"),
#'   stop = "converged",
#'   fields_to_save = G_conv$notes$fields_to_save
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
#' # Flat resident/house/landlord data frame for the current time
#' G_batch$report_result_df()
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
    resident_ownership_prop = 0.5,
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
    block_eth_threshold = 0,
    default_select_resident_prop = 0.25,
    default_criteria_resident = c("intercept", "ethnicity", "SES"),
    default_resident_softmax_beta = 1,
    default_criteria_landlord = c("intercept", "ethnicity", "SES"),
    default_block_eth_context = FALSE,
    default_include_landlord_seq = TRUE,
    sort_by_SES = FALSE,
    convergence_tol = 0.05,
    convergence_times = 10,
    convergence_maxit = 500,
    convergence_burnin = 100,
    separate_home_owners = FALSE,
    home_owner_selection_prop = 0.25,
    landlord_change_para = 30
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
  
  # house ownership
  house_owning <- make_group_labels(nrow(resident), c(1 - resident_ownership_prop,  resident_ownership_prop)) - 1
  resident <- data.frame(resident, 
                         house_owning = sample(house_owning))
  
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
                      landlord = as.vector(landlord_out$ownership),
                      landlord_duration = 1)
  
  # attach the building ID
  city_dim <- dim(city)
  house_map <- array(house$ID, dim = city_dim)
  building_map <- array(NA, dim = city_dim)
  k <- 0
  for(i in 1:city_dim[1]){for(j in 1:city_dim[2]){
    k <- k + 1
    building_map[i,j, ] <- k
  }}
  house <- data.frame(house,
                      building = as.vector(building_map)) 
  house[is.na(house$block), ] <- NA

  # house neighbor indices
  house_neighbor_ind <- neighbor_indices(city)
  
  # Rcpp用に0-indexed化し、NAは-1で表現（初期化時に1回だけ変換）
  house_neighbor_ind_0idx <- house_neighbor_ind - 1L
  house_neighbor_ind_0idx[is.na(house_neighbor_ind_0idx)] <- -1L
  storage.mode(house_neighbor_ind_0idx) <- "integer"
  
  add_field(G, State(house), 
            State(house_neighbor_ind),
            State(house_neighbor_ind_0idx))
  
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
    criteria_landlord = default_criteria_landlord,
    block_eth_context = default_block_eth_context,
    include_landlord_seq = default_include_landlord_seq,
    sort_by_SES = sort_by_SES,
    block_eth_threshold = block_eth_threshold,
    separate_home_owners = separate_home_owners,
    home_owner_selection_prop = home_owner_selection_prop,
    landlord_change_para = landlord_change_para
    )
  add_field(G, State(act_defaults))
  
  
  #===================================
  # Active state
  #===================================
  
  #------ same_ethnicity_prop in each neighborhood --------
  house_neib_ethnicity <- function(){
    resident_ID <- self$record$selected_resident
    
    house_resident_ethnicity <- rep(NA_real_, length = self$settings$cell_n)
    house_resident_ethnicity[self$resident$house] <- self$resident$minority
    
    # C++で近傍割合を計算（apply()ループを置き換え）
    neighbor_ethnicity_prop <- neib_ethnicity_prop(self$house_neighbor_ind_0idx,
                                                   house_resident_ethnicity)
    
    # 以降は元のまま（BLAS最適化済みの行列積なのでRのままでよい）
    same_eth_prop <- t(neighbor_ethnicity_prop %*% t(self$resident_mat_ethnicity[resident_ID, , drop = FALSE]))
    same_eth_prop
  }
  
  
  add_field(G, Active(house_neib_ethnicity))
  
  #------ same SES_prop in each neghborhood --------
  house_neib_SES <- function(){
    resident_ID <- self$record$selected_resident
    
    house_resident_SES <- rep(NA_real_, length = self$settings$cell_n)
    house_resident_SES[self$resident$house] <- self$resident$SES
    
    SES_prop <- neib_SES_prop(self$house_neighbor_ind_0idx, house_resident_SES)
    
    same_SES_prop <- t(SES_prop %*% t(self$resident_mat_SES[resident_ID, , drop = FALSE]))
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
                              n = NULL, 
                              sort_by_SES = self$act_defaults$sort_by_SES,
                              separate_home_owners = self$act_defaults$separate_home_owners,
                              home_owner_selection_prop = self$act_defaults$home_owner_selection_prop){
    if(is.null(n)){
      n <- floor(self$settings$resident_n * prop)
    }
    
    
    # select agent ID (with randamization)
    if(isTRUE(separate_home_owners)){
      home_owners_ID <- self$resident$ID[self$resident$house_owning == 1]
      non_owners_ID  <- self$resident$ID[self$resident$house_owning == 0]
      home_owner_n <- round(home_owner_selection_prop * n)
      non_owner_n <- n - home_owner_n
      id <- sample(c(rABM::sample2(home_owners_ID, size = home_owner_n), 
              rABM::sample2(non_owners_ID, size = non_owner_n))) 
    }else{
      id <- rABM::sample2(seq_len(self$settings$resident_n), size = n)
    }
    
    # sort by SES
    if(isTRUE(sort_by_SES)){
      id_SES <- 6 - self$resident$SES[id]
      id <- sort_by(id, id_SES)
    }
    
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
  landlord_decline <- function(criteria = self$act_defaults$criteria_landlord,
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
                                criteria_landlord = self$act_defaults$criteria_landlord,
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
          
          ## Block context come into effect only if it is over the threshold
          if(minority_prop >= self$act_defaults$block_eth_threshold){
            block_eth_score     <- minority_prop * self$settings$block_eth_eff * (-12)
            house_score         <- house_score + block_eth_score
          }
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
  
  
  # landlord_change ------------------------------
  landlord_change <- function(){
    # decide whether the landlord changes
    change_prob <- pweibull(self$house$landlord_duration,
                            shape = 2,
                            scale = self$act_defaults$landlord_change_para * sqrt(2))
    
    owner_change <- runif(length(change_prob), min = 0, max = 1) < change_prob
    
    # change the owners (permitting unchange)
    self$house$landlord[which(owner_change)] <- rABM::sample2(
      seq_len(self$settings$landlord_n),
      sum(owner_change, na.rm = TRUE),
      replace = TRUE
    )
    
    # update the count
    self$house$landlord_duration[which(owner_change)]  <- 1
    self$house$landlord_duration[which(!owner_change)] <- self$house$landlord_duration[which(!owner_change)] + 1
  }
  add_field(G, Act(landlord_change))
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
  
  
  # report result df ------------------------
  report_result_df <- function(time = NULL){
    if(!is.null(time)){
      if(length(self$log) < time){stop("'time' is outside of the length of logged time.")}
      self <- self$log[[time]]
    }
    
    # merge resident and house
    df <- self$resident
    df <- data.frame(resident_ID = df$ID, 
                     resident_minority = df$minority,
                     resident_SES = df$SES,
                     resident_preference_ethnicity = df$preference_ethnicity,
                     resident_preference_SES = df$preference_SES, 
                     house = df$house,
                     resident_owning = df$house_owning)
    df2 <- merge(df, self$house, by.x = "house", by.y = "ID", all.x = TRUE, sort = FALSE)
    df3 <- cbind(df2[ ,2:7], house_ID = df2$house, house_block = df2$block, 
                 house_building = df2$building, 
                 landlord_ID = df2$landlord,
                 landlord_duration = df2$landlord_duration)
    # further merge landlord
    df4 <- merge(df3, self$landlord, by.x = "landlord_ID", by.y = "ID", all.x = TRUE, sort = FALSE)
    df5 <- data.frame(df4[ ,2:10],
                      landlord_ID = df4$landlord_ID, 
                      landlord_minority_aversion = df4$minority_aversion,
                      landlord_SES_aversion = df4$SES_aversion,
                      landlord_dispersed = df4$dispersed,
                      landlord_duration = df4$landlord_duration)
    # return
    df5
  }
  add_field(G, Report(report_result_df))
  
  #================================================
  # stop FUN
  #================================================
  
  # add convergence settings to the field settings
  G$settings$convergence_tol   <- convergence_tol
  G$settings$convergence_times <- convergence_times
  G$settings$convergence_maxit <- convergence_maxit 
  G$settings$convergence_D_prev <- 0
  G$settings$convergence_burnin <- convergence_burnin
  
  convergence <- list()
  convergence$prev_D <- c()
  add_field(G, State(convergence))
  
  # convergence evaluation
  converged <- function(){
    # skip if within the time of burnin
    if(self$time <= self$settings$convergence_burnin){return(FALSE)}
    
    # merge the dataset
    df <- merge(self$resident[ ,c("house", "minority")], self$house[ ,c("ID", "block")], 
                by.x = "house", by.y = "ID", all.x = TRUE)
    df$majority <- 1 - df$minority
    
    # calculate D index
    overall_sum_mino <- sum(df$minority)
    overall_sum_majo <- sum(df$majority)
    block_sum_mino <- tapply(df$minority, df$block, sum)
    block_sum_majo <- tapply(df$majority, df$block, sum)
    D <- sum(abs(block_sum_mino/overall_sum_mino - block_sum_majo/overall_sum_majo))/2
    
    # update D history (keep only the last convergence_times values)
    D_history <- tail(c(self$convergence$prev_D, D), self$settings$convergence_times)
    self$convergence$prev_D <- D_history
    
    # evaluation of the stop condition
    if(length(D_history) >= self$settings$convergence_times &&
       diff(range(D_history)) < self$settings$convergence_tol){
      return(TRUE)
    }
    
    self$time >= self$settings$convergence_maxit
  }
  
  
  add_field(G, Stop(converged))
  
  
  #================================================
  # output G
  #================================================
  
  G
}
