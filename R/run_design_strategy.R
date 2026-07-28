#' Run the sparse-MET decision-point framework end to end
#'
#' @description
#' `run_design_strategy()` chains the Colmant et al. (2026) decision points into
#' a single call so a strategy can be specified and evaluated in one place:
#' (1) choose environments, (2) choose the tested individuals, (3) allocate
#' entries across environments (optionally refined for G x E prediction), (4)
#' recommend replication, and finally evaluate the result. Each stage is
#' optional -- supply only the pieces relevant to the decision you are testing --
#' and each underlying tool remains usable on its own. For the plantable field
#' book, pass the returned allocation to [plan_sparse_met_design()].
#'
#' @param treatments Character vector of all candidate line IDs.
#' @param environments Character vector of environment names. If `NULL`, they are
#'   chosen from `env_relationship` (decision 1).
#' @param n_test_entries_per_environment Entries per environment (decision 3).
#' @param G,Sigma_E Genomic relationship and environment covariance, required for
#'   individual selection, G x E refinement, and evaluation.
#' @param env_relationship,n_environments,env_select_method Decision 1: choose
#'   `n_environments` from the `env_relationship` matrix
#'   (see [build_environment_relationship()] and [select_environments()]).
#' @param n_individuals,individual_criterion Decision 2: choose `n_individuals`
#'   tested lines from `treatments` by a training-set criterion
#'   (see [select_individuals()]). If `NULL`, all `treatments` are used.
#' @param allocation_method,target_replications,common_treatments,balance
#'   Decision 3: passed to [allocate_sparse_met()].
#' @param refine_gxe Logical. If `TRUE` (and `G`/`Sigma_E` supplied), refine the
#'   allocation with [optimize_allocation_gxe()].
#' @param recommend_reps Logical. If `TRUE`, run [recommend_replication()]
#'   (decision 4) on the allocation.
#' @param replication_levels,seed_available,seed_required_per_plot Passed to
#'   [recommend_replication()].
#' @param evaluate Logical. If `TRUE` (and `G`/`Sigma_E` supplied), evaluate the
#'   final allocation with [met_information()] and [simulate_met()].
#' @param n_sim,sigma_g2,sigma_e2,bv_target,target_envs Evaluation settings
#'   passed to [simulate_met()].
#' @param seed,max_dim Reproducibility and size controls.
#'
#' @return A list with `decisions` (chosen `environments`, `individuals`,
#'   `allocation_method`), `allocation_matrix` (final, after any refinement),
#'   `replication` (from [recommend_replication()] or `NULL`), and `evaluation`
#'   (`met_information` summary and `simulate_met` outcomes, or `NULL`).
#'
#' @seealso [select_environments()], [select_individuals()],
#'   [allocate_sparse_met()], [optimize_allocation_gxe()],
#'   [recommend_replication()], [simulate_met()], [plan_sparse_met_design()].
#' @export
run_design_strategy <- function(treatments,
                                environments = NULL,
                                n_test_entries_per_environment,
                                G = NULL, Sigma_E = NULL,
                                env_relationship = NULL, n_environments = NULL,
                                env_select_method = "representative",
                                n_individuals = NULL,
                                individual_criterion = "cdmean",
                                allocation_method = "random_balanced",
                                target_replications = NULL,
                                common_treatments = NULL,
                                balance = "none",
                                refine_gxe = FALSE,
                                recommend_reps = FALSE,
                                replication_levels = c(1, 1.5, 2),
                                seed_available = NULL,
                                seed_required_per_plot = NULL,
                                evaluate = TRUE, n_sim = 30L,
                                sigma_g2 = 1, sigma_e2 = 1,
                                bv_target = "across_tpe", target_envs = NULL,
                                seed = NULL, max_dim = 6000L) {

  treatments <- unique(as.character(treatments))

  ## Decision 1 -- environments -------------------------------------------------
  if (is.null(environments)) {
    if (is.null(env_relationship) || is.null(n_environments))
      stop("Supply `environments`, or `env_relationship` + `n_environments`.")
    environments <- select_environments(env_relationship, n_environments,
                                        method = env_select_method,
                                        seed = seed)$selected
  }

  ## Decision 2 -- individuals --------------------------------------------------
  individuals <- treatments
  if (!is.null(n_individuals)) {
    if (is.null(G)) stop("`G` is required to select individuals.")
    Gsel <- as.matrix(G)[intersect(treatments, rownames(G)),
                         intersect(treatments, rownames(G)), drop = FALSE]
    individuals <- select_individuals(Gsel, n_train = n_individuals,
                                      criterion = individual_criterion,
                                      seed = seed)$selected
  }

  ## Decision 3 -- allocation ---------------------------------------------------
  alloc <- allocate_sparse_met(
    treatments = individuals, environments = environments,
    allocation_method = allocation_method,
    n_test_entries_per_environment = n_test_entries_per_environment,
    target_replications = target_replications,
    common_treatments = common_treatments,
    balance = balance, seed = seed)
  M <- alloc$allocation_matrix

  ## Align Sigma_E with the environments actually retained in the allocation. --
  ## Decision 1 may select a subset of environments, so the supplied covariance
  ## (defined over the candidate environments) must be subset to match `M`.
  Sigma_E_M <- Sigma_E
  if (!is.null(Sigma_E)) {
    Sigma_E <- as.matrix(Sigma_E)
    env_M <- colnames(M)
    if (!is.null(colnames(Sigma_E)) && all(env_M %in% colnames(Sigma_E))) {
      Sigma_E_M <- Sigma_E[env_M, env_M, drop = FALSE]
    } else if (nrow(Sigma_E) != length(env_M)) {
      stop("`Sigma_E` must be named over the environments or match the ",
           "number of selected environments.")
    }
  }

  refined <- FALSE
  if (refine_gxe) {
    if (is.null(G) || is.null(Sigma_E))
      stop("`G` and `Sigma_E` are required for `refine_gxe = TRUE`.")
    Gm <- as.matrix(G)[rownames(M), rownames(M), drop = FALSE]
    M  <- optimize_allocation_gxe(M, G = Gm, Sigma_E = Sigma_E_M,
                                  sigma_g2 = sigma_g2, sigma_e2 = sigma_e2,
                                  seed = seed, max_dim = max_dim)$allocation_matrix
    refined <- TRUE
  }

  ## Decision 4 -- replication --------------------------------------------------
  replication <- NULL
  if (recommend_reps) {
    if (is.null(G) || is.null(Sigma_E))
      stop("`G` and `Sigma_E` are required for `recommend_reps = TRUE`.")
    Gm <- as.matrix(G)[rownames(M), rownames(M), drop = FALSE]
    replication <- recommend_replication(
      M, G = Gm, Sigma_E = Sigma_E_M, sigma_g2 = sigma_g2, sigma_e2 = sigma_e2,
      replication_levels = replication_levels, n_sim = n_sim,
      seed_available = seed_available,
      seed_required_per_plot = seed_required_per_plot,
      seed = seed, max_dim = max_dim)
  }

  ## Evaluation -----------------------------------------------------------------
  evaluation <- NULL
  if (evaluate && !is.null(G) && !is.null(Sigma_E)) {
    Gm  <- as.matrix(G)[rownames(M), rownames(M), drop = FALSE]
    info <- met_information(M, G = Gm, Sigma_E = Sigma_E_M,
                            sigma_g2 = sigma_g2, sigma_e2 = sigma_e2,
                            max_dim = max_dim)
    sim  <- simulate_met(M, G = Gm, Sigma_E = Sigma_E_M,
                         sigma_g2 = sigma_g2, sigma_e2 = sigma_e2,
                         n_sim = n_sim, bv_target = bv_target,
                         target_envs = target_envs, seed = seed, max_dim = max_dim)
    evaluation <- list(mean_PEV = info$mean_PEV, CDmean = info$CDmean,
                       accuracy_mean = sim$accuracy_mean, gain_mean = sim$gain_mean,
                       common_selected_mean = sim$common_selected_mean,
                       avg_rank_mean = sim$avg_rank_mean)
  }

  list(
    decisions = list(environments = environments,
                     individuals = individuals,
                     allocation_method = alloc$summary$allocation_method,
                     gxe_refined = refined),
    allocation_matrix = M,
    replication = replication,
    evaluation = evaluation
  )
}
