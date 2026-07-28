#' Validate design-stage criteria against realised selection outcomes
#'
#' @description
#' Cheap design-stage criteria (mean PEV, CDmean, ...) are only useful if they
#' actually track the outcomes a breeder cares about -- selection accuracy and
#' genetic gain. Mothukuri et al. (2025) show that this link, while real, is
#' often weak and programme-specific, so it should be measured rather than assumed.
#'
#' `validate_design_criteria()` takes a set of candidate allocations, computes
#' for each both the design-based criteria (via [met_information()]) and the
#' realised outcomes (via [simulate_met()]), and reports their correlation across
#' designs. A strong negative correlation between `mean_PEV` and accuracy (or
#' positive between `CDmean` and accuracy) means the cheap criterion is a
#' trustworthy proxy for your programme; a weak one is a warning that you should
#' rely on simulation, not the criterion, when choosing a design.
#'
#' @param allocations A named list of \eqn{J \times E} allocation matrices
#'   (for example the `allocation_matrix` outputs of [allocate_sparse_met()]).
#'   Use a diverse set -- balanced, sparse, and optimised -- so the correlation
#'   is informative.
#' @inheritParams simulate_met
#' @param select_fraction Top fraction selected for the outcome metrics.
#'
#' @return A list with:
#' \describe{
#'   \item{`table`}{One row per design: `design`, the design-based `mean_PEV`
#'     and `CDmean`, and the realised `accuracy`, `gain`, `common_selected`, and
#'     `avg_rank`.}
#'   \item{`correlations`}{A data frame of Pearson correlations between each
#'     design-stage criterion and each realised outcome across the designs (with
#'     the number of designs, `n`).}
#' }
#'
#' @references
#' Mothukuri, S. R., Beyene, Y., Gultas, M., Burgueno, J., & Griebel, S. (2025).
#' Optimization of sparse phenotyping strategy in multi-environmental trials in
#' maize. *Theoretical and Applied Genetics*, 138, 62.
#'
#' @seealso [met_information()], [simulate_met()], [sensitivity_varcomp()].
#' @export
validate_design_criteria <- function(allocations, G, Sigma_E = NULL,
                                     sigma_g2 = 1, sigma_e2 = 1,
                                     env_efficiency = NULL,
                                     n_sim = 30L, select_fraction = 0.1,
                                     seed = NULL, max_dim = 6000L) {
  if (!is.list(allocations) || length(allocations) < 2L)
    stop("`allocations` must be a named list of at least two allocation matrices.")
  nms <- names(allocations)
  if (is.null(nms)) nms <- paste0("design_", seq_along(allocations))

  rows <- lapply(seq_along(allocations), function(i) {
    A <- allocations[[i]]
    info <- met_information(A, G = G, Sigma_E = Sigma_E,
                            sigma_g2 = sigma_g2, sigma_e2 = sigma_e2,
                            env_efficiency = env_efficiency,
                            target = "across_tpe", max_dim = max_dim)
    sm <- simulate_met(A, G = G, Sigma_E = Sigma_E,
                       sigma_g2 = sigma_g2, sigma_e2 = sigma_e2,
                       env_efficiency = env_efficiency,
                       n_sim = n_sim, select_fraction = select_fraction,
                       seed = seed, max_dim = max_dim)
    data.frame(
      design          = nms[i],
      mean_PEV        = info$mean_PEV,
      CDmean          = info$CDmean,
      accuracy        = sm$accuracy_mean,
      gain            = sm$gain_mean,
      common_selected = sm$common_selected_mean,
      avg_rank        = sm$avg_rank_mean,
      stringsAsFactors = FALSE
    )
  })
  tab <- do.call(rbind, rows)

  crit_cols <- c("mean_PEV", "CDmean")
  out_cols  <- c("accuracy", "gain", "common_selected", "avg_rank")
  cor_rows <- list()
  for (cc in crit_cols) for (oc in out_cols) {
    r <- suppressWarnings(stats::cor(tab[[cc]], tab[[oc]]))
    cor_rows[[length(cor_rows) + 1L]] <- data.frame(
      criterion = cc, outcome = oc, correlation = r,
      n = nrow(tab), stringsAsFactors = FALSE)
  }
  list(table = tab, correlations = do.call(rbind, cor_rows))
}
