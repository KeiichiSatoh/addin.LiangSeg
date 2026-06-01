#' Generate Random Samples from a Multivariate Normal Distribution
#'
#' A convenience wrapper around \code{\link[MASS]{mvrnorm}} that accepts a
#' mean vector, a vector of standard deviations, and a correlation matrix
#' instead of a covariance matrix.  The covariance matrix is constructed
#' internally as \eqn{\Sigma = \mathrm{diag}(\sigma) \, R \, \mathrm{diag}(\sigma)}.
#'
#' @param n     Integer. Number of observations to generate.
#' @param mu    Numeric vector of length \eqn{k}. Mean of each variable.
#' @param sigma Numeric vector of length \eqn{k} of strictly positive standard
#'   deviations.
#' @param R     Numeric \eqn{k \times k} correlation matrix.  Must be symmetric
#'   positive definite with unit diagonal.
#' @param print Logical. If \code{TRUE}, prints the sample means, standard
#'   deviations, and Pearson correlation matrix of the generated data to the
#'   console and returns a named list; if \code{FALSE} (default), returns the
#'   matrix directly.
#'
#' @details
#' The covariance matrix passed to \code{MASS::mvrnorm} is
#' \deqn{\Sigma = \mathrm{diag}(\sigma) \, R \, \mathrm{diag}(\sigma)}
#' so that \eqn{\mathrm{Var}(X_j) = \sigma_j^2} and
#' \eqn{\mathrm{Corr}(X_i, X_j) = R_{ij}}.
#'
#' @return
#' \describe{
#'   \item{\code{print = FALSE}}{A numeric matrix of dimensions
#'     \code{n} \eqn{\times} \code{k} with column names
#'     \code{x1, x2, \ldots, xk}.}
#'   \item{\code{print = TRUE}}{A named list with components:
#'     \describe{
#'       \item{\code{X}}{The numeric matrix described above.}
#'       \item{\code{mean}}{Sample mean vector of length \eqn{k}.}
#'       \item{\code{sd}}{Sample standard deviation vector of length \eqn{k}.}
#'       \item{\code{cor}}{Sample Pearson correlation matrix
#'         (\eqn{k \times k}).}
#'     }
#'   }
#' }
#'
#' @examples
#' set.seed(42)
#'
#' # Bivariate normal with different means, SDs, and a moderate correlation
#' R2 <- matrix(c(1.0, 0.6,
#'                0.6, 1.0), nrow = 2)
#' mat <- rmvnorm(n = 1000, mu = c(10, 20), sigma = c(2, 5), R = R2)
#' colMeans(mat)     # should be near c(10, 20)
#' apply(mat, 2, sd) # should be near c(2, 5)
#' cor(mat)          # should be near R2
#'
#' # With print = TRUE: prints summary and returns list
#' res <- rmvnorm(n = 1000, mu = c(10, 20), sigma = c(2, 5), R = R2,
#'                print = TRUE)
#' res$mean
#' res$sd
#' res$cor
#'
#' # Trivariate case
#' R3 <- matrix(c(1.0, 0.5, 0.3,
#'                0.5, 1.0, 0.2,
#'                0.3, 0.2, 1.0), nrow = 3)
#' res3 <- rmvnorm(n = 1000, mu = c(0, 5, 100), sigma = c(1, 2, 10), R = R3,
#'                 print = TRUE)
#'
#' @seealso
#' \code{\link[MASS]{mvrnorm}} for the underlying sampler.
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats cor sd
#' @export
rmvnorm <- function(n, mu, sigma, R, print = FALSE) {
  k <- length(mu)

  # ---- input checks ----------------------------------------------------------
  if (length(sigma) != k) {
    stop("'sigma' must have the same length as 'mu' (", k, ").")
  }
  if (any(sigma <= 0)) {
    stop("All elements of 'sigma' must be strictly positive.")
  }
  if (!is.matrix(R) || nrow(R) != k || ncol(R) != k) {
    stop("'R' must be a ", k, " x ", k, " matrix matching length(mu).")
  }
  if (!isTRUE(all.equal(diag(R), rep(1, k)))) {
    stop("'R' must be a correlation matrix with unit diagonal.")
  }
  if (!isTRUE(all.equal(R, t(R)))) {
    stop("'R' must be symmetric.")
  }
  ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  if (any(ev <= 0)) {
    stop("'R' must be positive definite.")
  }

  # ---- build covariance matrix and sample ------------------------------------
  D <- diag(sigma, nrow = k)
  Sigma <- D %*% R %*% D

  X <- MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma)
  colnames(X) <- paste0("x", seq_len(k))

  if (!print) return(X)

  # ---- print summary ---------------------------------------------------------
  sample_mean <- colMeans(X)
  sample_sd   <- apply(X, 2, sd)
  sample_cor  <- cor(X)

  message("Sample means:")
  print(round(sample_mean, 4))
  message("Sample standard deviations:")
  print(round(sample_sd, 4))
  message("Sample Pearson correlation matrix:")
  print(round(sample_cor, 4))

  return(invisible(list(X = X, mean = sample_mean, sd = sample_sd,
                        cor = sample_cor)))
}
