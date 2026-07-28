#' Factor-analytic environment covariance
#'
#' @description
#' A factor-analytic (FA) structure is the standard parsimonious model for the
#' genetic variance-covariance between environments in multi-environment trials
#' (Smith, Cullis & Thompson; Burgueno et al.; used by Mothukuri et al. 2025 and
#' Montesinos-Lopez et al. 2023). It models
#' \deqn{\Sigma_E = \Lambda \Lambda^\top + \Psi,}
#' where \eqn{\Lambda} is an \eqn{E \times k} matrix of loadings on `k` latent
#' factors and \eqn{\Psi} is a diagonal matrix of environment-specific
#' variances. A low-rank `k` captures the dominant correlation structure with
#' few parameters, which stabilises the coupled information matrix
#' ([met_information()]) and the simulation ([simulate_met()]) when the number
#' of environments is large.
#'
#' `fa_sigma_e()` either builds \eqn{\Sigma_E} directly from supplied loadings
#' and specific variances, or approximates a supplied full covariance by its
#' best rank-`k` FA form (the leading eigen-components plus a diagonal residual
#' that preserves the original variances).
#'
#' @param Sigma_E Optional full \eqn{E \times E} covariance to approximate. If
#'   supplied, a rank-`n_factors` FA approximation is returned.
#' @param n_factors Number of latent factors `k` (rank of \eqn{\Lambda\Lambda^\top}).
#' @param loadings Optional \eqn{E \times k} loadings matrix \eqn{\Lambda}
#'   (used when `Sigma_E` is not supplied).
#' @param psi Optional length-`E` (or scalar) vector of specific variances
#'   \eqn{\Psi}; defaults to a small positive value.
#'
#' @return An \eqn{E \times E} positive-definite covariance matrix of FA form.
#'
#' @references
#' Burgueno, J., Crossa, J., Cornelius, P. L., & Yang, R.-C. (2011).
#' Using factor analytic models for joining environments and genotypes without
#' crossover genotype x environment interaction. *Crop Science*, 51(4),
#' 1584-1596.
#'
#' @seealso [met_information()], [simulate_met()].
#' @export
fa_sigma_e <- function(Sigma_E = NULL, n_factors = 1L,
                       loadings = NULL, psi = NULL) {
  k <- as.integer(n_factors)
  if (k < 1L) stop("`n_factors` must be a positive integer.")

  if (!is.null(Sigma_E)) {
    S <- as.matrix(Sigma_E)
    if (nrow(S) != ncol(S)) stop("`Sigma_E` must be square.")
    E <- nrow(S)
    if (k >= E) return(S)                      # nothing to reduce
    eg  <- eigen((S + t(S)) / 2, symmetric = TRUE)
    Lam <- eg$vectors[, seq_len(k), drop = FALSE] %*%
           diag(sqrt(pmax(eg$values[seq_len(k)], 0)), k)
    Psi <- pmax(diag(S) - rowSums(Lam^2), 1e-8)  # residual variances, non-negative
    out <- Lam %*% t(Lam) + diag(Psi, E)
    dimnames(out) <- dimnames(S)
    return((out + t(out)) / 2)
  }

  if (is.null(loadings))
    stop("Supply either `Sigma_E` (to approximate) or `loadings` (to build).")
  Lam <- as.matrix(loadings)
  E <- nrow(Lam)
  if (ncol(Lam) != k)
    stop("`loadings` must have `n_factors` columns.")
  if (is.null(psi)) psi <- 1e-2
  if (length(psi) == 1L) psi <- rep(psi, E)
  if (length(psi) != E) stop("`psi` must have length E or be a scalar.")
  out <- Lam %*% t(Lam) + diag(psi, E)
  (out + t(out)) / 2
}
