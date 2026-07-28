#' Explore the sparse-testing trade-off (individuals x per-environment x environments)
#'
#' @description
#' The single largest driver of accuracy in Colmant et al. (2026) was *how* the
#' target populations are sampled: how many individuals are available, how many
#' are tested per environment, and how many environments are used, under a fixed
#' plot budget. `sparsity_grid()` sweeps that three-way grid, builds a random
#' sparse allocation for each cell, and scores realised accuracy and genetic gain
#' by simulation -- turning "how sparse should we go?" into a comparison you can
#' read off a table.
#'
#' @param G Genomic relationship matrix (row/column names) with at least
#'   `max(ia_values)` genotypes.
#' @param Sigma_E Environment covariance with at least `max(noe_values)`
#'   environments (a broad-TPE structure makes the sweep informative).
#' @param ia_values Numbers of individuals *available* to sample from the TPG.
#' @param ipf_values Numbers of individuals *tested per environment*.
#' @param noe_values Numbers of environments (each >= 2).
#' @param plot_budget Optional cap on total plots (`ipf * noe`); cells exceeding
#'   it are skipped.
#' @inheritParams simulate_met
#'
#' @return A data frame (sorted by `accuracy_mean`, best first) with one row per
#'   feasible grid cell: `ia`, `ipf`, `noe`, `total_plots`, `sparsity_pct`
#'   (relative to a fully balanced design), `accuracy_mean`, `gain_mean`,
#'   `common_selected_mean`, `mean_PEV`, and `CDmean`.
#'
#' @references
#' Colmant, A., Pita, F., & Covarrubias-Pazaran, G. (2026). Some objective
#' functions and ideas to optimize experimental designs in artificial selection
#' programs. *Crop Science*, 66, e70337.
#'
#' @seealso [simulate_met()], [recommend_replication()].
#' @export
sparsity_grid <- function(G, Sigma_E,
                          ia_values, ipf_values, noe_values,
                          plot_budget = NULL,
                          sigma_g2 = 1, sigma_e2 = 1,
                          n_sim = 30L, select_fraction = 0.1,
                          bv_target = "across_tpe", target_envs = NULL,
                          seed = NULL, max_dim = 6000L) {
  G <- as.matrix(G)
  if (is.null(rownames(G))) rownames(G) <- colnames(G) <- paste0("L", seq_len(nrow(G)))
  Sigma_E <- as.matrix(Sigma_E)
  if (max(ia_values) > nrow(G))
    stop("`G` has fewer genotypes than max(ia_values).")
  if (max(noe_values) > nrow(Sigma_E))
    stop("`Sigma_E` has fewer environments than max(noe_values).")
  if (any(noe_values < 2)) stop("`noe_values` must all be >= 2.")

  base_seed <- if (is.null(seed)) 1L else seed
  rows <- list()

  for (ia in ia_values) for (ipf in ipf_values) for (noe in noe_values) {
    if (ipf > ia) next
    total_plots <- ipf * noe
    if (!is.null(plot_budget) && total_plots > plot_budget) next

    lines <- rownames(G)[seq_len(ia)]
    envs  <- paste0("E", seq_len(noe))
    Gsub    <- G[lines, lines, drop = FALSE]
    SigEsub <- Sigma_E[seq_len(noe), seq_len(noe), drop = FALSE]

    # Random sparse allocation: each environment gets ipf distinct lines; a line
    # may be tested in 0, 1 or more environments (untested lines are predicted).
    set.seed(base_seed + ia + 7L * ipf + 13L * noe)
    M <- matrix(0L, ia, noe, dimnames = list(lines, envs))
    for (e in seq_len(noe)) M[sample.int(ia, ipf), e] <- 1L
    if (all(rowSums(M) == 0) || all(colSums(M) == 0)) next

    sm <- tryCatch(
      simulate_met(M, G = Gsub, Sigma_E = SigEsub,
                   sigma_g2 = sigma_g2, sigma_e2 = sigma_e2,
                   n_sim = n_sim, select_fraction = select_fraction,
                   bv_target = bv_target, target_envs = target_envs,
                   seed = base_seed, max_dim = max_dim),
      error = function(e) NULL)
    if (is.null(sm)) next

    rows[[length(rows) + 1L]] <- data.frame(
      ia = ia, ipf = ipf, noe = noe,
      total_plots = total_plots,
      sparsity_pct = 100 * (1 - ipf / ia),
      accuracy_mean = sm$accuracy_mean,
      gain_mean = sm$gain_mean,
      common_selected_mean = sm$common_selected_mean,
      mean_PEV = sm$mean_PEV,
      CDmean = sm$CDmean,
      stringsAsFactors = FALSE)
  }

  if (!length(rows)) stop("No feasible grid cells (check plot_budget and ipf <= ia).")
  tab <- do.call(rbind, rows)
  tab[order(-tab$accuracy_mean), , drop = FALSE]
}
