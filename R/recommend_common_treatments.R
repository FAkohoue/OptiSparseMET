#' Recommend the number of common treatments from environmental kinship
#'
#' @description
#' Common (check-like) treatments tested in every environment provide the
#' cross-environment connectivity that lets a sparse MET be analysed jointly.
#' How many are needed depends on how correlated the environments are: strongly
#' correlated sites borrow information freely and need few, weakly correlated
#' sites need many. `recommend_common_treatments()` reads that correlation
#' straight from the breeder's **environmental kinship** (the environment
#' covariance `Sigma_E`, or a relationship matrix from
#' [build_environment_relationship()]), so the breeder does not have to guess a
#' connectivity level. It reports the typical correlation `rho`, the analytic
#' [suggest_common_range()] band, and -- if `G` and `Sigma_E` are supplied --
#' the simulated optimum from [tune_common_treatments()].
#'
#' @param env_relationship Environment covariance `Sigma_E` (default) or an
#'   environment correlation/relationship matrix (`is_covariance = FALSE`).
#' @param n_test_entries_per_environment,n_treatments Passed to
#'   [suggest_common_range()]; the per-environment capacity and total treatment
#'   pool.
#' @param n_environments Number of environments; defaults to `nrow(env_relationship)`.
#' @param sparse_allocation `"random"` or `"disjoint"` (see
#'   [suggest_common_range()]).
#' @param target_se,n_points Passed to [suggest_common_range()].
#' @param is_covariance If `TRUE` (default) `env_relationship` is a covariance and
#'   is converted to a correlation before averaging; set `FALSE` if it is already
#'   a correlation/similarity in \[-1, 1\].
#' @param evaluate If `TRUE` and `G`/`Sigma_E` are given, also run
#'   [tune_common_treatments()] over the suggested range and return the simulated
#'   optimum.
#' @param G,Sigma_E Genomic relationship and environment covariance for the
#'   optional simulation (`evaluate = TRUE`). Treatments/environments are taken
#'   from their dimnames.
#' @param ... Further arguments passed to [tune_common_treatments()].
#' @return A list with `rho` (typical environment correlation), `common_range`
#'   (suggested counts), `suggest` (the full [suggest_common_range()] object),
#'   and, when evaluated, `tune` (the [tune_common_treatments()] object) and
#'   `recommended` (its `best_estimated` optimum).
#' @seealso [suggest_common_range()], [tune_common_treatments()],
#'   [build_environment_relationship()].
#' @examples
#' SigE <- matrix(c(1, .3, .1, .3, 1, .2, .1, .2, 1), 3)
#' recommend_common_treatments(SigE, n_test_entries_per_environment = 60,
#'                             n_treatments = 300, n_environments = 3)$rho
#' @export
recommend_common_treatments <- function(env_relationship,
                                        n_test_entries_per_environment,
                                        n_treatments, n_environments = NULL,
                                        sparse_allocation = c("random", "disjoint"),
                                        target_se = 0.15, n_points = 8L,
                                        is_covariance = TRUE,
                                        evaluate = FALSE, G = NULL,
                                        Sigma_E = NULL, ...) {
  sparse_allocation <- match.arg(sparse_allocation)
  D <- as.matrix(env_relationship)
  if (is.null(n_environments)) n_environments <- nrow(D)

  R <- if (is_covariance) stats::cov2cor(D) else D
  rho <- mean(R[upper.tri(R)])
  if (!is.finite(rho)) rho <- 0
  rho <- max(0, min(0.999, rho))                 # connectivity uses a low-side rho

  sr <- suggest_common_range(n_test_entries_per_environment, n_treatments,
                             n_environments, sparse_allocation = sparse_allocation,
                             rho = rho, target_se = target_se, n_points = n_points)

  out <- list(rho = rho, common_range = sr$common_values, suggest = sr)

  if (evaluate && !is.null(G) && !is.null(Sigma_E)) {
    G <- as.matrix(G); Sigma_E <- as.matrix(Sigma_E)
    trts <- rownames(G); if (is.null(trts)) trts <- paste0("L", seq_len(nrow(G)))
    envs <- colnames(Sigma_E)
    if (is.null(envs)) envs <- paste0("E", seq_len(nrow(Sigma_E)))
    tc <- tune_common_treatments(
      treatments = trts, environments = envs,
      n_test_entries_per_environment = n_test_entries_per_environment,
      G = G, Sigma_E = Sigma_E, common_values = sr$common_values,
      sparse_allocation = sparse_allocation, ...)
    out$tune <- tc
    out$recommended <- tc$optima$best_estimated
  }
  out
}
