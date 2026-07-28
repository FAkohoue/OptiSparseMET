#' Recommend a replication level from accuracy, cost, and seed constraints (decision 4)
#'
#' @description
#' `recommend_replication()` sweeps candidate replication levels for a given
#' allocation and reports the realised accuracy, genetic gain, and design-based
#' PEV at each level (via [simulate_met()]), then recommends the level at the
#' point of diminishing returns. It is the resource-aware form of decision 4 in
#' Colmant et al. (2026): it evaluates replication at the programme's own spatial
#' noise (`sigma_e2`), relationship matrix, and -- optionally -- seed budget,
#' and it varies replication with the number of unique entries held fixed, so
#' the two are not confounded (as they were in the paper's field-size grid).
#'
#' @details
#' A replication level `r` scales the effective plots per tested cell
#' (`reps = allocation_matrix * r`), so fractional (p-rep) levels are allowed.
#' If `seed_available` (a per-line seed budget) and `seed_required_per_plot` are
#' supplied, a level is flagged infeasible when total required seed exceeds the
#' budget.
#'
#' @inheritParams simulate_met
#' @param replication_levels Numeric vector of replication levels to evaluate.
#' @param min_gain Minimum accuracy improvement (per unit replication) required
#'   to justify the next level; the recommendation is the largest level whose
#'   marginal accuracy gain still meets `min_gain`, else the smallest level.
#' @param seed_available Optional total number of seeds available per line
#'   (scalar) for a feasibility flag.
#' @param seed_required_per_plot Optional seeds needed per plot (scalar).
#'
#' @return A list with `table` (a data frame: replication, accuracy_mean,
#'   gain_mean, mean_PEV, feasible) and `recommended` (the chosen level).
#'
#' @seealso [simulate_met()], [assign_replication_by_seed()].
#' @export
recommend_replication <- function(allocation_matrix, G, Sigma_E = NULL,
                                  sigma_g2 = 1, sigma_e2 = 1,
                                  replication_levels = c(1, 1.5, 2, 3),
                                  n_sim = 30L, min_gain = 0.02,
                                  seed_available = NULL,
                                  seed_required_per_plot = NULL,
                                  seed = NULL, max_dim = 6000L) {
  replication_levels <- sort(unique(replication_levels))
  J <- nrow(allocation_matrix)

  rows <- lapply(replication_levels, function(r) {
    reps <- allocation_matrix * r
    sm <- simulate_met(allocation_matrix, G = G, Sigma_E = Sigma_E,
                       sigma_g2 = sigma_g2, sigma_e2 = sigma_e2,
                       reps = reps, n_sim = n_sim, seed = seed, max_dim = max_dim)
    feasible <- TRUE
    if (!is.null(seed_available) && !is.null(seed_required_per_plot)) {
      plots_per_line <- rowSums(reps)
      feasible <- all(plots_per_line * seed_required_per_plot <= seed_available)
    }
    data.frame(replication = r,
               accuracy_mean = sm$accuracy_mean,
               gain_mean = sm$gain_mean,
               mean_PEV = sm$mean_PEV,
               feasible = feasible)
  })
  tab <- do.call(rbind, rows)

  # Recommend the largest feasible level whose marginal accuracy gain per unit
  # replication still meets min_gain; fall back to the smallest feasible level.
  feas <- tab[tab$feasible, , drop = FALSE]
  if (nrow(feas) == 0L) {
    recommended <- NA_real_
  } else {
    recommended <- feas$replication[1]
    if (nrow(feas) > 1L) {
      for (i in 2:nrow(feas)) {
        dacc <- feas$accuracy_mean[i] - feas$accuracy_mean[i - 1]
        dr   <- feas$replication[i]  - feas$replication[i - 1]
        if (dacc / dr >= min_gain) recommended <- feas$replication[i] else break
      }
    }
  }

  list(table = tab, recommended = recommended)
}
