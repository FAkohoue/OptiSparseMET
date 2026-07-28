#' Variance-component sensitivity of a MET design (P3.5)
#'
#' @description
#' Design optimality is conditional on the assumed variance components. A design
#' recommended under one variance ratio may not be best under another.
#' `sensitivity_varcomp()` sweeps the residual-to-genetic variance ratio
#' \eqn{\lambda = \sigma_e^2/\sigma_g^2} for a fixed allocation and reports the
#' across-TPE mean PEV and CDmean at each value (via [met_information()]), so a
#' programme can judge whether its design choice is robust to variance
#' misspecification.
#'
#' @inheritParams met_information
#' @param lambda_grid Numeric vector of \eqn{\sigma_e^2/\sigma_g^2} ratios to
#'   evaluate (with `sigma_g2` fixed at 1).
#'
#' @return A data frame with columns `lambda`, `mean_PEV`, and `CDmean`.
#'
#' @seealso [met_information()], [simulate_met()].
#' @export
sensitivity_varcomp <- function(allocation_matrix, G, Sigma_E = NULL,
                                lambda_grid = c(0.25, 0.5, 1, 2, 4, 8),
                                reps = NULL, env_efficiency = NULL,
                                max_dim = 6000L) {
  lambda_grid <- sort(unique(lambda_grid))
  rows <- lapply(lambda_grid, function(lam) {
    info <- met_information(allocation_matrix, G = G, Sigma_E = Sigma_E,
                            sigma_g2 = 1, sigma_e2 = lam,
                            reps = reps, env_efficiency = env_efficiency,
                            target = "across_tpe", max_dim = max_dim)
    data.frame(lambda = lam, mean_PEV = info$mean_PEV, CDmean = info$CDmean)
  })
  do.call(rbind, rows)
}
