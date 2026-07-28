#' Robustly optimise a common-treatment set
#'
#' @description
#' Selects the size and composition of the global connectivity set in a sparse
#' multi-environment trial (MET). A common treatment is required to occur in
#' every environment. This function decides presence only; local plot
#' replication is a separate, downstream seed-feasible precision decision.
#'
#' The set is a design-stage decision, not a quantity that changes during field
#' deployment. When `n_common = NULL`, the function searches feasible sizes and
#' identities. Supplying `n_common` fixes only the number to be selected; the
#' function still selects the identities. Once the returned `selected` vector
#' is adopted as `common_treatments` in allocation or planning, those identities
#' are fixed as present in every environment of that released design.
#'
#' Candidate common-set sizes are compared by a multi-criterion, risk-averse
#' objective combining:
#' \itemize{
#'   \item prediction reliability of the target-population-of-environments mean;
#'   \item precision of pairwise environment connectivity;
#'   \item genetic effective sample size of the common set; and
#'   \item testing breadth retained for non-common treatments.
#' }
#' Reliability is integrated over a prior of variance-component and
#' environment-covariance scenarios. Pairwise connectivity is the number of
#' distinct shared treatments. Genetic effective sample size is retained as a
#' separate diversity criterion, not substituted for the observed
#' shared-treatment count. Correlation-dependent targets are
#' calculated separately for every environment pair. By default, the global
#' common set is chosen with maximin protection of the least-connected pair;
#' lower-tail CVaR and mean pairwise aggregation are available explicitly.
#' Family and genomic-group coverage, entry capacities, and one network-wide
#' seed inventory are enforced.
#'
#' Replication never increases the pairwise connectivity count: one genotype
#' shared by two environments remains one shared genotype whether it has one
#' plot or several plots at either environment. Use
#' [assign_replication_by_seed()] or [plan_sparse_met_design()] after allocation
#' to assign feasible local replication from the remaining seed inventory.
#'
#' @details
#' For a proposed common set with relationship submatrix \eqn{G_C}, its genetic
#' effective sample size is
#' \deqn{n_{\mathrm{eff}} =
#'   \{\operatorname{tr}(G_C)\}^2 / \operatorname{tr}(G_C^2).}
#' Thus duplicated or strongly related anchors reduce the diversity score
#' without altering the number of distinct shared identities. For environments
#' \eqn{a,b}, let \eqn{S_{ab}} be the set of distinct genotypes present in both
#' environments. Connectivity is \eqn{|S_{ab}|}; repeating a member does not
#' change it. The shared-treatment count is compared with the conservative
#' correlation-precision
#' requirement
#' \eqn{\max_s\{(1-\rho_{ab,s}^2)/\mathrm{target\_se}\}^2} across prior
#' scenarios \eqn{s}.
#'
#' Within a scenario, the planning-stage posterior covariance for common
#' treatment \eqn{j}, using one presence record per environment, is computed
#' from
#' \deqn{\left[
#'   \{\sigma_{g,s}^2 G_{jj}\Sigma_{E,s}\}^{-1} +
#'   \operatorname{diag}(\epsilon/\sigma_{e,s}^2)
#'   \right]^{-1}.}
#' Its TPE-mean PEV is the mean of that matrix. Scenario reliabilities are then
#' aggregated by their prior mean, worst case, or lower-tail CVaR. These compact
#' calculations make the common-set search scalable; the delivered whole-trial
#' design should still be evaluated with [met_information()] and
#' [simulate_met()].
#'
#' @param G Named genomic or pedigree relationship matrix.
#' @param Sigma_E Named environment covariance matrix.
#' @param treatment_info Data frame with columns `Treatment` and `Family`.
#' @param seed_info Data frame with columns `Treatment` and `SeedAvailable`.
#' @param seed_required_per_plot Positive environment-specific mandatory seed
#'   debit for one selected common presence, as a named vector or a data frame
#'   with columns `Environment` and `SeedRequiredPerPlot`. Supply the per-plot
#'   requirement when one plot is mandatory; multiply it by the mandatory local
#'   plot count when the local design requires more than one plot.
#' @param entry_capacities Named positive integer vector giving the number of
#'   distinct candidate treatments accommodated in each environment.
#' @param n_common Optional single common-set size. If supplied, the size is
#'   fixed but the function still selects the member identities. If `NULL`, a
#'   nested grid of feasible sizes and their candidate identities is optimised.
#' @param minimum_seed_buffer Non-negative seed reserve retained for every
#'   treatment after placing one plot of a selected common treatment in every
#'   environment. Default 0.
#' @param scenarios Prior from [robust_scenarios()]. Defaults to residual
#'   variances `c(0.5, 1, 2)` and environment-covariance shrinkage
#'   `c(0, 0.25, 0.5)`.
#' @param aggregate Risk aggregation: `"cvar"` (default), `"mean"`, or `"min"`.
#' @param cvar_alpha Lower-tail probability used by CVaR.
#' @param env_efficiency Named vector of local-design efficiencies in `(0, 1]`.
#' @param target_se Target standard error for pairwise environment correlations.
#' @param pair_aggregate How to protect connectivity across environment pairs.
#'   `"maximin"` (default) maximises the weakest pair, `"cvar"` maximises the
#'   lower-tail mean, and `"mean"` maximises average pairwise attainment.
#' @param pair_cvar_alpha Lower-tail fraction used only when
#'   `pair_aggregate = "cvar"`.
#' @param min_per_family Minimum number of common treatments per family.
#' @param n_groups Number of genomic groups to cover; defaults to the number of
#'   families, bounded by the feasible candidate pool.
#' @param objective_weights Named non-negative weights for `reliability`,
#'   `connectivity`, `genetic_diversity`, and `testing_breadth`.
#' @param require_full_coverage If `TRUE`, reject common-set sizes that make it
#'   arithmetically impossible to test every treatment at least once.
#' @param max_count_evaluations Maximum number of automatically generated
#'   common-set sizes. The grid spans the entire feasible interval.
#' @param seed Optional integer controlling genomic clustering.
#'
#' @return A list containing `selected`, `n_common`, `common_presence` (a binary
#'   common-treatment by environment matrix), `comparison`, `pareto`,
#'   `candidate_sets`, `selection_diagnostics`, `pairwise_connectivity`,
#'   `seed_ledger`, and `rationale`.
#'
#' @seealso [suggest_common_treatments()], [plan_sparse_met_design()],
#'   [robust_scenarios()]
#' @export
optimize_common_treatments <- function(
    G,
    Sigma_E,
    treatment_info,
    seed_info,
    seed_required_per_plot,
    entry_capacities,
    n_common = NULL,
    minimum_seed_buffer = 0,
    scenarios = NULL,
    aggregate = c("cvar", "mean", "min"),
    cvar_alpha = 0.25,
    env_efficiency = 1,
    target_se = 0.15,
    pair_aggregate = c("maximin", "cvar", "mean"),
    pair_cvar_alpha = 0.25,
    min_per_family = 1L,
    n_groups = NULL,
    objective_weights = c(
      reliability = 0.40,
      connectivity = 0.30,
      genetic_diversity = 0.15,
      testing_breadth = 0.15
    ),
    require_full_coverage = TRUE,
    max_count_evaluations = 15L,
    seed = NULL
) {
  aggregate <- match.arg(aggregate)
  pair_aggregate <- match.arg(pair_aggregate)
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
        abs(seed - round(seed)) > 1e-8)
      stop("`seed` must be one finite integer or NULL.")
    set.seed(as.integer(round(seed)))
  }

  ti <- as.data.frame(treatment_info, stringsAsFactors = FALSE)
  si <- as.data.frame(seed_info, stringsAsFactors = FALSE)
  if (!all(c("Treatment", "Family") %in% names(ti)))
    stop("`treatment_info` needs columns `Treatment` and `Family`.")
  if (!all(c("Treatment", "SeedAvailable") %in% names(si)))
    stop("`seed_info` needs columns `Treatment` and `SeedAvailable`.")
  trts <- as.character(ti$Treatment)
  if (!length(trts) || anyNA(trts) || any(!nzchar(trts)) ||
      anyDuplicated(trts))
    stop("`treatment_info$Treatment` must contain unique, non-missing IDs.")
  fam <- as.character(ti$Family)
  fam[is.na(fam) | !nzchar(fam)] <- "UNGROUPED"
  names(fam) <- trts
  seed_ids <- as.character(si$Treatment)
  if (anyNA(seed_ids) || any(!nzchar(seed_ids)) || anyDuplicated(seed_ids))
    stop("`seed_info$Treatment` must contain unique, non-missing IDs.")
  seed_available <- as.numeric(si$SeedAvailable)
  if (any(!is.finite(seed_available)) || any(seed_available < 0))
    stop("`SeedAvailable` must contain finite, non-negative values.")
  names(seed_available) <- seed_ids
  if (!all(trts %in% seed_ids))
    stop("`seed_info` must cover every treatment in `treatment_info`.")
  seed_available <- seed_available[trts]

  G <- as.matrix(G)
  if (!is.numeric(G) || nrow(G) != ncol(G) || any(!is.finite(G)) ||
      is.null(rownames(G)) || is.null(colnames(G)) ||
      anyDuplicated(rownames(G)) || anyDuplicated(colnames(G)) ||
      !all(trts %in% rownames(G)) || !all(trts %in% colnames(G)))
    stop("`G` must be a finite, named square matrix covering all treatments.")
  G <- G[trts, trts, drop = FALSE]
  if (!isTRUE(all.equal(G, t(G), tolerance = 1e-8)))
    stop("`G` must be symmetric.")
  if (any(diag(G) <= 0))
    stop("`G` must have a strictly positive diagonal.")
  eig_G <- eigen((G + t(G)) / 2, symmetric = TRUE,
                 only.values = TRUE)$values
  if (min(eig_G) < -1e-8 * max(1, max(abs(eig_G))))
    stop("`G` must be positive semidefinite.")

  Sigma_E <- as.matrix(Sigma_E)
  if (!is.numeric(Sigma_E) || nrow(Sigma_E) != ncol(Sigma_E) ||
      any(!is.finite(Sigma_E)) || is.null(rownames(Sigma_E)) ||
      is.null(colnames(Sigma_E)) || anyDuplicated(rownames(Sigma_E)) ||
      anyDuplicated(colnames(Sigma_E)))
    stop("`Sigma_E` must be a finite, named square matrix.")
  if (!isTRUE(all.equal(Sigma_E, t(Sigma_E), tolerance = 1e-8)))
    stop("`Sigma_E` must be symmetric.")
  if (any(diag(Sigma_E) <= 0))
    stop("`Sigma_E` must have a strictly positive diagonal.")
  eig_Sigma <- eigen((Sigma_E + t(Sigma_E)) / 2, symmetric = TRUE,
                     only.values = TRUE)$values
  if (min(eig_Sigma) < -1e-8 * max(1, max(abs(eig_Sigma))))
    stop("`Sigma_E` must be positive semidefinite.")
  envs <- rownames(Sigma_E)
  if (!setequal(envs, colnames(Sigma_E)))
    stop("`Sigma_E` row and column names must contain the same environments.")
  Sigma_E <- Sigma_E[envs, envs, drop = FALSE]
  E <- length(envs)
  J <- length(trts)
  if (E < 2L)
    stop("At least two environments are required to optimize connectivity.")

  .env_vector <- function(x, nm, integer = FALSE, lower = -Inf,
                          strictly = FALSE) {
    if (is.data.frame(x)) {
      value_col <- if (identical(nm, "seed_required_per_plot"))
        "SeedRequiredPerPlot" else NULL
      if (is.null(value_col) ||
          !all(c("Environment", value_col) %in% names(x)))
        stop("`", nm, "` has an unsupported data-frame format.")
      if (anyDuplicated(as.character(x$Environment)))
        stop("`", nm, "` contains duplicate environments.")
      vals <- as.numeric(x[[value_col]])
      names(vals) <- as.character(x$Environment)
      x <- vals
    }
    nx <- names(x)
    x <- as.numeric(x)
    if (length(x) == 1L) {
      x <- rep(x, E)
      names(x) <- envs
    } else {
      if (is.null(nx) || !all(envs %in% nx))
        stop("`", nm, "` must be scalar or named for every environment.")
      names(x) <- nx
      x <- x[envs]
    }
    bad <- !is.finite(x) | if (strictly) x <= lower else x < lower
    if (any(bad))
      stop("`", nm, "` contains invalid values.")
    if (integer && any(abs(x - round(x)) > 1e-8))
      stop("`", nm, "` must contain integer values.")
    if (integer) x <- as.integer(round(x))
    stats::setNames(x, envs)
  }

  seed_req <- .env_vector(seed_required_per_plot,
                          "seed_required_per_plot", lower = 0,
                          strictly = TRUE)
  entry_cap <- .env_vector(entry_capacities, "entry_capacities",
                           integer = TRUE, lower = 0, strictly = TRUE)
  efficiency <- .env_vector(env_efficiency, "env_efficiency",
                             lower = 0, strictly = TRUE)
  if (any(efficiency > 1))
    stop("`env_efficiency` values must not exceed 1.")
  if (!is.numeric(minimum_seed_buffer) ||
      length(minimum_seed_buffer) != 1L ||
      !is.finite(minimum_seed_buffer) || minimum_seed_buffer < 0)
    stop("`minimum_seed_buffer` must be one finite non-negative value.")
  if (!is.numeric(target_se) || length(target_se) != 1L ||
      !is.finite(target_se) || target_se <= 0)
    stop("`target_se` must be one finite positive value.")
  if (!is.numeric(pair_cvar_alpha) || length(pair_cvar_alpha) != 1L ||
      !is.finite(pair_cvar_alpha) || pair_cvar_alpha <= 0 ||
      pair_cvar_alpha > 1)
    stop("`pair_cvar_alpha` must be one value in (0, 1].")
  if (!is.numeric(cvar_alpha) || length(cvar_alpha) != 1L ||
      !is.finite(cvar_alpha) || cvar_alpha <= 0 || cvar_alpha > 1)
    stop("`cvar_alpha` must be in (0, 1].")
  if (!is.numeric(min_per_family) || length(min_per_family) != 1L ||
      !is.finite(min_per_family) ||
      abs(min_per_family - round(min_per_family)) > 1e-8 ||
      min_per_family < 0)
    stop("`min_per_family` must be one non-negative integer.")
  min_per_family <- as.integer(round(min_per_family))
  if (!is.logical(require_full_coverage) ||
      length(require_full_coverage) != 1L || is.na(require_full_coverage))
    stop("`require_full_coverage` must be `TRUE` or `FALSE`.")

  if (is.null(scenarios)) {
    scenarios <- robust_scenarios(
      sigma_g2 = 1,
      sigma_e2 = c(0.5, 1, 2),
      sigmaE_shrink = c(0, 0.25, 0.5)
    )
  }
  valid_scenario <- vapply(scenarios, function(sc)
    is.list(sc) &&
      all(c("sigma_g2", "sigma_e2", "sigmaE_shrink", "prob") %in% names(sc)) &&
      all(vapply(sc[c("sigma_g2", "sigma_e2", "sigmaE_shrink", "prob")],
                 function(z) is.numeric(z) && length(z) == 1L &&
                   is.finite(z), logical(1))) &&
      sc$sigma_g2 > 0 && sc$sigma_e2 > 0 &&
      sc$sigmaE_shrink >= 0 && sc$sigmaE_shrink <= 1 && sc$prob >= 0,
    logical(1))
  if (!length(scenarios) || !all(valid_scenario))
    stop("`scenarios` must be a valid non-empty robust-scenario prior.")
  scenario_sigma <- function(sc) {
    base <- if (!is.null(sc$Sigma_E)) as.matrix(sc$Sigma_E) else Sigma_E
    if (nrow(base) != E || ncol(base) != E ||
        is.null(rownames(base)) || is.null(colnames(base)) ||
        !setequal(rownames(base), envs) || !setequal(colnames(base), envs))
      stop("Every scenario-specific `Sigma_E` must cover the same environments.")
    base <- base[envs, envs, drop = FALSE]
    if (any(!is.finite(base)) ||
        !isTRUE(all.equal(base, t(base), tolerance = 1e-8)))
      stop("Every scenario-specific `Sigma_E` must be finite and symmetric.")
    (1 - sc$sigmaE_shrink) * base +
      sc$sigmaE_shrink * diag(diag(base))
  }
  scenario_prob <- vapply(scenarios, `[[`, numeric(1), "prob")
  if (sum(scenario_prob) <= 0)
    stop("Scenario probabilities must not all be zero.")
  scenario_prob <- scenario_prob / sum(scenario_prob)

  required_weights <- c("reliability", "connectivity",
                        "genetic_diversity", "testing_breadth")
  if (is.null(names(objective_weights)) ||
      !all(required_weights %in% names(objective_weights)))
    stop("`objective_weights` must name: ",
         paste(required_weights, collapse = ", "), ".")
  objective_weights <- as.numeric(objective_weights[required_weights])
  names(objective_weights) <- required_weights
  if (any(!is.finite(objective_weights)) ||
      any(objective_weights < 0) || sum(objective_weights) <= 0)
    stop("`objective_weights` must be finite, non-negative, and not all zero.")
  objective_weights <- objective_weights / sum(objective_weights)

  # Common-set selection funds exactly one presence plot per environment.
  # Additional local plots are assigned later and cannot increase connectivity.
  mandatory_seed <- sum(seed_req)
  feasible <- trts[
    seed_available >= mandatory_seed + minimum_seed_buffer - 1e-10
  ]
  if (!length(feasible))
    stop("No treatment can fund common presence across all ",
         "environments while retaining the requested seed buffer.")

  family_levels <- unique(fam[trts])
  family_feasible_n <- table(factor(fam[feasible], levels = family_levels))
  missing_family <- family_levels[family_feasible_n < min_per_family]
  if (length(missing_family))
    stop("Seed feasibility cannot provide `min_per_family` common treatments ",
         "for: ", paste(missing_family, collapse = ", "), ".")

  if (is.null(n_groups))
    n_groups <- length(family_levels)
  if (!is.numeric(n_groups) || length(n_groups) != 1L ||
      !is.finite(n_groups) || abs(n_groups - round(n_groups)) > 1e-8 ||
      n_groups < 1)
    stop("`n_groups` must be one positive integer.")
  n_groups <- as.integer(round(n_groups))
  n_groups <- min(n_groups, length(feasible))
  groups <- .cluster_G_groups(G[feasible, feasible, drop = FALSE], n_groups)

  centrality <- rowMeans(G[feasible, trts, drop = FALSE])
  centrality <- .common_rescale01(centrality)
  seed_slack <- (seed_available[feasible] - mandatory_seed -
                   minimum_seed_buffer) /
    pmax(seed_available[feasible], 1e-12)
  seed_slack <- .common_rescale01(seed_slack)
  base_rank <- 0.75 * centrality + 0.25 * seed_slack

  chosen <- character(0)
  for (ff in family_levels) {
    candidates <- feasible[fam[feasible] == ff]
    ord <- candidates[order(base_rank[candidates], decreasing = TRUE)]
    chosen <- c(chosen, utils::head(ord, min_per_family))
  }
  chosen <- unique(chosen)
  for (gg in sort(unique(groups))) {
    if (!any(groups[intersect(chosen, names(groups))] == gg)) {
      candidates <- setdiff(names(groups)[groups == gg], chosen)
      if (length(candidates))
        chosen <- c(chosen,
                    candidates[which.max(base_rank[candidates])])
    }
  }
  chosen <- unique(chosen)

  genetic_distance <- outer(diag(G), diag(G), `+`) - 2 * G
  genetic_distance[genetic_distance < 0 & genetic_distance > -1e-8] <- 0
  genetic_distance <- sqrt(pmax(genetic_distance, 0))
  while (length(chosen) < length(feasible)) {
    remaining <- setdiff(feasible, chosen)
    diversity <- if (!length(chosen)) {
      rep(1, length(remaining))
    } else {
      apply(genetic_distance[remaining, chosen, drop = FALSE], 1L, min)
    }
    diversity <- .common_rescale01(diversity)
    names(diversity) <- remaining
    fill_score <- 0.55 * diversity + 0.30 * centrality[remaining] +
      0.15 * seed_slack[remaining]
    chosen <- c(chosen, remaining[which.max(fill_score)])
  }
  nested_order <- chosen
  # The coverage construction above is the exact lower bound for this nested
  # order; using its prefix length avoids relying on an arbitrary GxE floor.
  family_group_prefix <- 0L
  for (ii in seq_along(nested_order)) {
    prefix <- nested_order[seq_len(ii)]
    family_ok <- all(table(factor(fam[prefix],
                                  levels = family_levels)) >= min_per_family)
    group_ok <- setequal(unique(groups[prefix]), unique(groups))
    if (family_ok && group_ok) {
      family_group_prefix <- ii
      break
    }
  }
  minimum_count <- max(1L, family_group_prefix)

  max_common <- min(length(feasible), min(entry_cap))
  if (isTRUE(require_full_coverage) && E > 1L) {
    coverage_limit <- floor((sum(entry_cap) - J) / (E - 1L))
    max_common <- min(max_common, coverage_limit)
  }
  if (max_common < minimum_count)
    stop("Entry capacities cannot accommodate family/group connectivity and ",
         "full treatment coverage. Feasible common count is at most ",
         max_common, " but coverage requires at least ", minimum_count, ".")

  if (!is.null(n_common)) {
    if (!is.numeric(n_common) || length(n_common) != 1L ||
        !is.finite(n_common) || abs(n_common - round(n_common)) > 1e-8)
      stop("`n_common` must be one finite integer or NULL.")
    counts <- as.integer(round(n_common))
  } else {
    if (!is.numeric(max_count_evaluations) ||
        length(max_count_evaluations) != 1L ||
        !is.finite(max_count_evaluations) ||
        abs(max_count_evaluations - round(max_count_evaluations)) > 1e-8 ||
        max_count_evaluations < 2)
      stop("`max_count_evaluations` must be an integer of at least 2.")
    max_count_evaluations <- as.integer(round(max_count_evaluations))
    n_values <- min(max_count_evaluations,
                    max_common - minimum_count + 1L)
    counts <- unique(as.integer(round(seq(minimum_count, max_common,
                                          length.out = n_values))))
  }
  if (!length(counts) || anyNA(counts) || any(counts < minimum_count) ||
      any(counts > max_common))
    stop("Every evaluated common count must lie in [", minimum_count, ", ",
         max_common, "].")
  counts <- sort(unique(counts))

  pair_index <- which(upper.tri(Sigma_E), arr.ind = TRUE)
  pair_names <- paste(envs[pair_index[, 1L]],
                      envs[pair_index[, 2L]], sep = "__")
  target_by_scenario <- vapply(seq_along(scenarios), function(ss) {
    sc <- scenarios[[ss]]
    Sig <- scenario_sigma(sc)
    RR <- stats::cov2cor(Sig)
    rho <- RR[pair_index]
    pmax(2, ((1 - rho^2) / target_se)^2)
  }, numeric(nrow(pair_index)))
  if (is.null(dim(target_by_scenario)))
    target_by_scenario <- matrix(target_by_scenario,
                                 nrow = nrow(pair_index))
  pair_target <- apply(target_by_scenario, 1L, max)
  names(pair_target) <- pair_names
  aggregate_pair_connectivity <- function(pair_fraction) {
    pair_fraction <- pmin(1, pmax(0, as.numeric(pair_fraction)))
    if (pair_aggregate == "maximin")
      return(min(pair_fraction))
    if (pair_aggregate == "mean")
      return(mean(pair_fraction))
    .common_cvar(
      pair_fraction,
      rep(1 / length(pair_fraction), length(pair_fraction)),
      pair_cvar_alpha
    )
  }
  tpe_variance_ok <- vapply(scenarios, function(sc) {
    Sig <- scenario_sigma(sc)
    is.finite(mean(Sig)) && mean(Sig) > 0
  }, logical(1))
  if (!all(tpe_variance_ok))
    stop("Every robust scenario must imply positive variance for the ",
         "across-environment TPE mean.")

  aggregate_values <- function(values) {
    if (aggregate == "mean")
      return(sum(values * scenario_prob))
    if (aggregate == "min")
      return(min(values))
    .common_cvar(values, scenario_prob, cvar_alpha)
  }

  reliability_cache <- new.env(parent = emptyenv())
  individual_reliability <- function(gdiag, presence_j) {
    cache_key <- paste(
      formatC(gdiag, digits = 14L, format = "fg", flag = "#"),
      paste(as.integer(presence_j), collapse = ","),
      sep = "|"
    )
    if (exists(cache_key, envir = reliability_cache, inherits = FALSE))
      return(get(cache_key, envir = reliability_cache, inherits = FALSE))
    values <- vapply(seq_along(scenarios), function(ss) {
      sc <- scenarios[[ss]]
      Sig <- scenario_sigma(sc)
      K <- sc$sigma_g2 * gdiag * Sig
      K <- .bend_pd((K + t(K)) / 2, eps = 1e-8)
      precision <- solve(K) +
        diag(as.numeric(presence_j * efficiency / sc$sigma_e2), E)
      post <- tryCatch(solve(precision),
                       error = function(e) .pinv_sym_dense(precision))
      prior_var <- mean(K)
      post_var <- mean(post)
      pmin(1, pmax(0, 1 - post_var / prior_var))
    }, numeric(1))
    out <- aggregate_values(values)
    assign(cache_key, out, envir = reliability_cache)
    out
  }

  evaluate_count <- function(count) {
    ids <- nested_order[seq_len(count)]
    Gc <- G[ids, ids, drop = FALSE]
    eig <- pmax(eigen((Gc + t(Gc)) / 2, symmetric = TRUE,
                      only.values = TRUE)$values, 0)
    genetic_ess <- if (sum(eig^2) > 0)
      sum(eig)^2 / sum(eig^2) else 1
    genetic_fraction <- min(1, genetic_ess / count)

    # Binary common presence: every selected identity occurs once in the
    # planning incidence matrix at every environment. Plot replication is not
    # selected here.
    R <- matrix(1L, nrow = count, ncol = E,
                dimnames = list(ids, envs))
    seed_used <- rep(sum(seed_req), count)
    names(seed_used) <- ids
    rel_j <- vapply(seq_len(count), function(jj)
      individual_reliability(Gc[jj, jj], R[jj, ]), numeric(1))

    # Every selected common genotype is present in every environment. Each
    # identity contributes exactly one unit to each pair, irrespective of
    # relatedness or the later number of local plots.
    pair_info <- rep(count, nrow(pair_index))
    names(pair_info) <- pair_names
    connectivity <- aggregate_pair_connectivity(pair_info / pair_target)

    robust_reliability <- mean(rel_j)
    testing_breadth <- (sum(entry_cap) - (E - 1L) * count) /
      sum(entry_cap)
    testing_breadth <- pmin(1, pmax(0, testing_breadth))
    score <- objective_weights["reliability"] * robust_reliability +
      objective_weights["connectivity"] * connectivity +
      objective_weights["genetic_diversity"] * genetic_fraction +
      objective_weights["testing_breadth"] * testing_breadth
    list(
      ids = ids,
      presence = R,
      seed_used = seed_used,
      score = unname(score),
      reliability = robust_reliability,
      connectivity = connectivity,
      min_pair_connectivity = min(pmin(1, pair_info / pair_target)),
      genetic_ess = genetic_ess,
      genetic_fraction = genetic_fraction,
      testing_breadth = testing_breadth,
      pair_info = pair_info
    )
  }

  fits <- lapply(counts, evaluate_count)
  comparison <- data.frame(
    n_common = counts,
    objective = vapply(fits, `[[`, numeric(1), "score"),
    robust_reliability = vapply(fits, `[[`, numeric(1), "reliability"),
    pair_connectivity = vapply(fits, `[[`, numeric(1), "connectivity"),
    minimum_pair_connectivity =
      vapply(fits, `[[`, numeric(1), "min_pair_connectivity"),
    genetic_effective_n = vapply(fits, `[[`, numeric(1), "genetic_ess"),
    genetic_effective_fraction =
      vapply(fits, `[[`, numeric(1), "genetic_fraction"),
    testing_breadth = vapply(fits, `[[`, numeric(1), "testing_breadth"),
    stringsAsFactors = FALSE
  )
  dominated <- vapply(seq_len(nrow(comparison)), function(ii) {
    if (nrow(comparison) == 1L) return(FALSE)
    other <- seq_len(nrow(comparison)) != ii
    all_metrics <- c("robust_reliability", "pair_connectivity",
                     "genetic_effective_fraction", "testing_breadth")
    metric_delta <- sweep(
      as.matrix(comparison[other, all_metrics, drop = FALSE]),
      2L,
      as.numeric(unlist(
        comparison[ii, all_metrics, drop = FALSE],
        use.names = FALSE
      )),
      `-`
    )
    any(vapply(seq_len(nrow(metric_delta)), function(jj) {
      all(metric_delta[jj, ] >= -1e-10) &&
        any(metric_delta[jj, ] > 1e-10)
    }, logical(1)))
  }, logical(1))
  comparison$pareto <- !dominated
  best_idx <- which.max(comparison$objective)
  best <- fits[[best_idx]]

  all_seed_used <- stats::setNames(rep(0, J), trts)
  all_seed_used[best$ids] <- best$seed_used
  seed_ledger <- data.frame(
    Treatment = trts,
    IsCommon = trts %in% best$ids,
    SeedAvailable = as.numeric(seed_available[trts]),
    SeedAllocatedToCommonPresence = as.numeric(all_seed_used[trts]),
    SeedRemainingAfterCommonPresence =
      as.numeric(seed_available[trts] - all_seed_used[trts]),
    Feasible = seed_available[trts] - all_seed_used[trts] >=
      minimum_seed_buffer - 1e-8,
    stringsAsFactors = FALSE
  )
  pairwise_connectivity <- data.frame(
    Environment1 = envs[pair_index[, 1L]],
    Environment2 = envs[pair_index[, 2L]],
    TargetDistinctSharedTreatments = as.numeric(pair_target),
    AchievedDistinctSharedTreatments = as.numeric(best$pair_info),
    TargetFraction = pmin(1, as.numeric(best$pair_info / pair_target)),
    stringsAsFactors = FALSE
  )
  selection_diagnostics <- data.frame(
    Treatment = feasible,
    Family = unname(fam[feasible]),
    GeneticGroup = unname(groups[feasible]),
    GeneticCentralityScore = unname(centrality[feasible]),
    SeedHeadroomScore = unname(seed_slack[feasible]),
    NestedRank = match(feasible, nested_order),
    Selected = feasible %in% best$ids,
    stringsAsFactors = FALSE
  )
  selection_diagnostics <- selection_diagnostics[
    order(selection_diagnostics$NestedRank), , drop = FALSE
  ]

  list(
    selected = best$ids,
    n_common = length(best$ids),
    common_presence = best$presence,
    comparison = comparison,
    pareto = comparison[comparison$pareto, , drop = FALSE],
    candidate_sets = stats::setNames(lapply(fits, `[[`, "ids"), counts),
    selection_diagnostics = selection_diagnostics,
    pairwise_connectivity = pairwise_connectivity,
    seed_ledger = seed_ledger,
    rationale = list(
      method = "robust_nested_common_set_presence_only",
      aggregate = aggregate,
      cvar_alpha = cvar_alpha,
      pair_aggregate = pair_aggregate,
      pair_cvar_alpha = pair_cvar_alpha,
      scenario_count = length(scenarios),
      feasible_candidate_pool = length(feasible),
      minimum_count_from_family_and_group_coverage = minimum_count,
      maximum_count_from_capacity_and_coverage = max_common,
      objective_weights = objective_weights,
      identity_selection =
        "constrained_greedy_genetic_diversity_centrality_seed_headroom",
      replication_decision =
        "downstream_seed_feasible_local_precision_not_connectivity",
      minimum_seed_buffer = minimum_seed_buffer,
      selected_count_was_user_fixed = !is.null(n_common)
    )
  )
}


.common_rescale01 <- function(x) {
  nx <- names(x)
  x <- as.numeric(x)
  if (!length(x)) return(x)
  rr <- range(x, finite = TRUE)
  out <- if (!all(is.finite(rr)) || diff(rr) <= .Machine$double.eps)
    rep(0.5, length(x)) else (x - rr[1L]) / diff(rr)
  names(out) <- nx
  out
}


.common_cvar <- function(values, probs, alpha) {
  ord <- order(values)
  values <- values[ord]
  probs <- probs[ord] / sum(probs)
  take <- pmin(probs, pmax(0, alpha - c(0, head(cumsum(probs), -1L))))
  sum(values * take) / alpha
}
