# Internal helper: adjust a single pairwise Pearson target correlation to the
# corresponding latent Gaussian copula correlation via numerical inversion.
# Uses pbivnorm::pbivnorm() for the bivariate normal CDF.
#
# @param rho_target  Target Pearson (phi) correlation for the binary pair.
# @param p1, p2      Marginal success probabilities for the two variables.
# @param tol         Convergence tolerance passed to stats::uniroot().
# @return Scalar copula correlation rho such that the expected binary Pearson
#         correlation equals rho_target.
.adjust_copula_cor <- function(rho_target, p1, p2, tol = 1e-6) {
  q1 <- qnorm(p1)
  q2 <- qnorm(p2)
  denom <- sqrt(p1 * (1 - p1) * p2 * (1 - p2))

  # Expected binary Pearson r as a function of latent copula correlation rho
  expected_r <- function(rho) {
    (pbivnorm::pbivnorm(q1, q2, rho) - p1 * p2) / denom
  }

  # Feasible range: rho in (-1, 1), but expected_r may not span all of
  # (-1, 1) for the binary case; bracket tightly and catch uniroot failures.
  lo <- -1 + 1e-6
  hi <-  1 - 1e-6

  r_lo <- expected_r(lo)
  r_hi <- expected_r(hi)

  if (rho_target < r_lo || rho_target > r_hi) {
    stop(sprintf(
      "rho_target = %.4f is outside the feasible Pearson range [%.4f, %.4f] ",
      "for p1 = %.4f, p2 = %.4f.",
      rho_target, r_lo, r_hi, p1, p2
    ))
  }

  uniroot(function(rho) expected_r(rho) - rho_target,
          interval = c(lo, hi), tol = tol)$root
}


#' Generate Random Samples from a Multivariate Binary Distribution
#'
#' Generates \code{n} rows of correlated binary random variables using either
#' a Gaussian copula (default, any number of variables \eqn{k}) or the direct
#' Pearson/\eqn{\phi} method (bivariate case only).
#'
#' For the copula method, \code{Sigma} is interpreted as the target
#' \strong{Pearson} correlation matrix of the binary outcomes.  Internally,
#' each pairwise entry is converted to the corresponding latent Gaussian copula
#' correlation via numerical inversion, so the simulated binary variables
#' match the requested correlations in expectation.
#'
#' @param n Integer. Number of observations to generate.
#' @param p Numeric vector of length \eqn{k}. Marginal success probabilities
#'   \eqn{P(X_j = 1)} for each variable \eqn{j = 1, \ldots, k}.
#'   All values must lie strictly in \eqn{(0, 1)}.
#' @param Sigma Numeric \eqn{k \times k} \strong{target Pearson} correlation
#'   matrix of the binary variables.  Must be symmetric with unit diagonal.
#'   \describe{
#'     \item{\code{method = "copula"}}{Each off-diagonal entry is the desired
#'       Pearson (\eqn{\phi}) correlation of the corresponding binary pair.
#'       The matrix is converted internally to latent copula correlations via
#'       \code{.adjust_copula_cor()}, so the resulting positive-definiteness
#'       of the latent matrix is checked after conversion.}
#'     \item{\code{method = "pearson"}}{Only \code{Sigma[1,2]} is used as the
#'       target \eqn{\phi}; \code{k} must equal 2.}
#'   }
#' @param method Character string, either \code{"copula"} (default) or
#'   \code{"pearson"}.  Partial matching is supported.  \code{"pearson"} is
#'   only valid when \code{k = 2}; an error is thrown otherwise.
#' @param print Logical. If \code{TRUE}, prints the sample Pearson
#'   correlation(s) of the generated binary matrix to the console and returns
#'   a named list; if \code{FALSE} (default), returns the binary matrix
#'   directly.
#'   \describe{
#'     \item{\code{method = "pearson"}}{Prints a single sample \eqn{\phi}
#'       (= Pearson \eqn{r}) value.}
#'     \item{\code{method = "copula"}}{Prints the full \eqn{k \times k}
#'       sample Pearson correlation matrix.}
#'   }
#'
#' @details
#' \strong{Gaussian copula method} (\code{method = "copula"}):
#'
#' Because thresholding a continuous latent variable attenuates correlations,
#' naively passing the target Pearson matrix as the copula matrix would yield
#' binary correlations that are systematically lower than requested.  To
#' correct for this, each off-diagonal entry \eqn{r_{ij}} of \code{Sigma} is
#' replaced by the latent copula correlation \eqn{\rho_{ij}} satisfying
#' \deqn{r_{ij} = \frac{P(Z_i \le \Phi^{-1}(p_i),\, Z_j \le \Phi^{-1}(p_j);\,
#'   \rho_{ij}) - p_i p_j}{\sqrt{p_i(1-p_i)\,p_j(1-p_j)}}}
#' where \eqn{\Phi^{-1}} is the standard normal quantile function and
#' \eqn{P(\cdot;\rho)} denotes the standard bivariate normal CDF with
#' correlation \eqn{\rho}.  The inversion is performed numerically via
#' \code{stats::uniroot()} for each pair \eqn{(i, j)}.
#'
#' After all pairwise corrections the adjusted matrix is checked for positive
#' definiteness; if it is not (which can occur when \code{Sigma} has entries
#' near the Fréchet bounds), an error is thrown.
#'
#' A latent vector \eqn{\mathbf{Z} \sim \mathcal{N}_k(\mathbf{0},
#' \Sigma_{\text{adj}})} is then drawn via \code{MASS::mvrnorm()} and
#' thresholded:
#' \deqn{X_j = \mathbf{1}\{Z_j \le \Phi^{-1}(p_j)\}}
#'
#' \strong{Pearson / \eqn{\phi} method} (\code{method = "pearson"},
#' \eqn{k = 2} only):
#'
#' For binary variables, Pearson's \eqn{r} and the \eqn{\phi} coefficient
#' coincide:
#' \deqn{\phi = r = \frac{P(X_1=1,X_2=1) - p_1 p_2}
#'       {\sqrt{p_1(1-p_1)\,p_2(1-p_2)}}}
#' The four joint cell probabilities are derived directly from
#' \eqn{\phi = } \code{Sigma[1,2]}:
#' \deqn{P(X_1=1, X_2=1) = \phi\sqrt{p_1(1-p_1)p_2(1-p_2)} + p_1 p_2}
#' \deqn{P(X_1=1, X_2=0) = p_1 - P(X_1=1, X_2=1)}
#' \deqn{P(X_1=0, X_2=1) = p_2 - P(X_1=1, X_2=1)}
#' \deqn{P(X_1=0, X_2=0) = 1 - \text{(above three)}}
#' An error is raised if any cell probability falls outside \eqn{[0, 1]}.
#'
#' @return
#' \describe{
#'   \item{\code{print = FALSE}}{An integer matrix of dimensions
#'     \code{n} \eqn{\times} \code{k} with column names
#'     \code{x1, x2, \ldots, xk}.  All entries are 0 or 1.}
#'   \item{\code{print = TRUE}}{A named list with components:
#'     \describe{
#'       \item{\code{X}}{The integer matrix described above.}
#'       \item{\code{cor}}{For \code{method = "pearson"}: a scalar sample
#'         \eqn{\phi}.  For \code{method = "copula"}: the \eqn{k \times k}
#'         sample Pearson correlation matrix of \code{X}.}
#'     }
#'   }
#' }
#'
#' @examples
#' set.seed(42)
#'
#' # --- Gaussian copula: trivariate case ---
#' # Sigma is the TARGET Pearson correlation matrix of the binary outcomes.
#' Sigma3 <- matrix(c(1.0, 0.6, 0.3,
#'                    0.6, 1.0, 0.2,
#'                    0.3, 0.2, 1.0), nrow = 3)
#' res3 <- rbivbinom(n = 2000, p = c(0.3, 0.5, 0.7), Sigma = Sigma3,
#'                   print = TRUE)
#' # res3$cor should be close to Sigma3
#'
#' # --- Gaussian copula: bivariate case ---
#' Sigma2 <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
#' res2 <- rbivbinom(n = 1000, p = c(0.4, 0.6), Sigma = Sigma2, print = TRUE)
#'
#' # --- Pearson method: bivariate only ---
#' Phi <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
#' res_p <- rbivbinom(n = 1000, p = c(0.4, 0.6), Sigma = Phi,
#'                    method = "pearson", print = TRUE)
#' res_p$cor   # should be close to 0.5
#'
#' @seealso
#' \code{\link[stats]{rbinom}} for independent binary samples;
#' \code{\link[MASS]{mvrnorm}} for multivariate normal sampling;
#' \code{\link[pbivnorm]{pbivnorm}} for the bivariate normal CDF.
#'
#' @importFrom MASS mvrnorm
#' @importFrom pbivnorm pbivnorm
#' @importFrom stats cor qnorm uniroot
#' @export
rbivbinom <- function(n, p, Sigma, method = c("copula", "pearson"),
                      print = FALSE) {
  method <- match.arg(method)
  k <- length(p)

  # ---- input checks ----------------------------------------------------------
  if (!is.matrix(Sigma) || nrow(Sigma) != k || ncol(Sigma) != k) {
    stop("'Sigma' must be a ", k, " x ", k, " matrix matching length(p).")
  }
  if (any(p <= 0 | p >= 1)) {
    stop("All elements of 'p' must be strictly in (0, 1).")
  }
  if (method == "pearson" && k != 2) {
    stop("method = \"pearson\" is only supported for k = 2 variables.")
  }

  # ---- Pearson / phi method (k = 2 only) ------------------------------------
  if (method == "pearson") {
    phi <- Sigma[1, 2]
    p1 <- p[1]; p2 <- p[2]

    p11 <- phi * sqrt(p1 * (1 - p1) * p2 * (1 - p2)) + p1 * p2
    p10 <- p1 - p11
    p01 <- p2 - p11
    p00 <- 1 - p11 - p10 - p01

    probs <- c(p00, p01, p10, p11)
    if (any(probs < 0 | probs > 1)) {
      stop("Under the specified Sigma (phi) and p, the cell probabilities are ",
           "outside [0, 1]. Choose a feasible correlation value.")
    }

    # cell index: 1=(0,0), 2=(0,1), 3=(1,0), 4=(1,1)
    idx <- sample(1:4, n, replace = TRUE, prob = probs)
    X <- cbind(x1 = as.integer(idx %in% c(3, 4)),
               x2 = as.integer(idx %in% c(2, 4)))

    if (!print) return(X)

    sample_phi <- cor(X[, 1], X[, 2])
    message("Sample phi (= Pearson r) between x1 and x2: ",
            round(sample_phi, 4))
    return(invisible(list(X = X, cor = sample_phi)))
  }

  # ---- Gaussian copula method (general k) ------------------------------------

  # Step 1: convert target Pearson correlations to latent copula correlations
  Sigma_adj <- diag(k)
  for (i in seq_len(k - 1)) {
    for (j in seq(i + 1, k)) {
      rho_adj <- .adjust_copula_cor(Sigma[i, j], p[i], p[j])
      Sigma_adj[i, j] <- rho_adj
      Sigma_adj[j, i] <- rho_adj
    }
  }

  # Step 2: check positive definiteness of the adjusted matrix
  ev <- eigen(Sigma_adj, symmetric = TRUE, only.values = TRUE)$values
  if (any(ev <= 0)) {
    stop("The adjusted latent copula correlation matrix is not positive ",
         "definite. Consider using a target Sigma closer to the identity.")
  }

  # Step 3: draw latent MVN and threshold
  thresholds <- qnorm(p)
  Z <- MASS::mvrnorm(n = n, mu = rep(0, k), Sigma = Sigma_adj)
  X <- sweep(Z, 2, thresholds, FUN = `<=`) * 1L
  colnames(X) <- paste0("x", seq_len(k))

  if (!print) return(X)

  sample_cor <- cor(X)
  message("Sample Pearson correlation matrix:")
  print(round(sample_cor, 4))
  return(invisible(list(X = X, cor = sample_cor)))
}
