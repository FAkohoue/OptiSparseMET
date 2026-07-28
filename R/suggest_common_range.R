#' Suggest a range of common-treatment counts to evaluate
#'
#' @description
#' The number of common treatments has three governing bounds: a hard cap
#' (`C < k`, the per-environment capacity), a connectivity-driven upper bound
#' (enough genotypes shared between environment pairs to estimate their genetic
#' correlations, net of what the sparse overlap already provides), and a coverage
#' reference (where unique lines can no longer be covered). `suggest_common_range()`
#' turns these into a grid to pass as `common_values` to
#' [tune_common_treatments()].
#'
#' @details
#' The connectivity target follows the sampling error of a correlation: to
#' estimate a genetic correlation `rho` with standard error `target_se`,
#' \deqn{\text{target\_shared} \approx \left(\frac{1-\rho^2}{\text{target\_se}}\right)^2}
#' connecting genotypes are needed per environment pair. Random sparse allocation
#' already supplies about \eqn{k^2/J} shared genotypes per pair; disjoint
#' (non-overlapping) designs supply none, so common treatments must provide the
#' whole amount. Weak correlations (`rho` near 0) need more, which is why the
#' optimum rises as environments become weakly correlated.
#'
#' The coverage reference \eqn{(kE - J)/(E-1)} is reported but is *not* used as a
#' hard cap: under weak correlation the optimum can lie beyond it, trading
#' coverage for connectivity.
#'
#' @param n_test_entries_per_environment Plots per environment, `k`.
#' @param n_treatments Total number of candidate lines, `J`.
#' @param n_environments Number of environments, `E`.
#' @param sparse_allocation `"random"` (sparse overlap ~ \eqn{k^2/J}) or
#'   `"disjoint"` (no sparse overlap; common treatments carry all connectivity).
#' @param rho Typical genetic correlation between environments (use a low value
#'   for a broad/heterogeneous TPE). Default 0.5.
#' @param target_se Target standard error for the estimated correlation. Default
#'   0.15. Lower values demand more connectivity.
#' @param n_points Approximate number of grid points. Default 8.
#'
#' @return A list with `common_values` (integer vector for
#'   [tune_common_treatments()]) and `rationale` (the target connectivity, the
#'   expected sparse overlap, the connectivity-driven and coverage bounds, and
#'   the hard cap).
#'
#' @seealso [tune_common_treatments()], [suggest_safe_k()].
#' @export
suggest_common_range <- function(n_test_entries_per_environment,
                                 n_treatments, n_environments,
                                 sparse_allocation = c("random", "disjoint"),
                                 rho = 0.5, target_se = 0.15, n_points = 8L) {
  sparse_allocation <- match.arg(sparse_allocation)
  k <- as.integer(n_test_entries_per_environment)
  J <- as.integer(n_treatments)
  E <- as.integer(n_environments)
  if (k < 2L || J < 2L || E < 2L) stop("Need k, J >= 2 and E >= 2.")

  target_shared <- ((1 - rho^2) / target_se)^2
  overlap0 <- if (sparse_allocation == "random") (k^2) / J else 0
  C_need   <- max(0, ceiling(target_shared - overlap0))
  C_cover  <- floor((k * E - J) / (E - 1))                 # coverage reference (may be <= 0)
  C_hard   <- min(k - 1L, J - 1L)

  C_upper <- min(C_hard, max(C_need, 6L))                  # cover the connectivity rise
  pts     <- max(2L, min(as.integer(n_points), C_upper + 1L))
  vals    <- sort(unique(round(seq(0, C_upper, length.out = pts))))

  list(
    common_values = as.integer(vals),
    rationale = list(
      target_shared          = target_shared,
      expected_sparse_overlap = overlap0,
      connectivity_upper     = C_need,
      coverage_reference     = C_cover,
      hard_cap               = C_hard,
      swept_to               = C_upper
    )
  )
}
