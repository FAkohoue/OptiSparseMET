#' Construct reference designs for sparse MET benchmarking
#'
#' `make_met_benchmark_designs()` constructs a full-testing upper bound and
#' equal-budget sparse comparators. The sparse comparators comprise naive
#' within-environment random sampling, the M3-inspired coverage-first allocator,
#' and strict M4 equireplication when its degree sequence is feasible. An
#' infeasible M4 is reported, not silently relabelled as an exact M4 design.
#'
#' @param treatments,environments Unique treatment and environment identifiers.
#' @param n_test_entries_per_environment Integer scalar or environment-named
#'   vector giving the number of test entries at each site.
#' @param target_replications Target number of environments per non-common
#'   treatment for M3 and M4.
#' @param common_treatments Treatments fixed as present at every site. Their
#'   local plot replication remains design-specific.
#' @param m4_fallback `"omit"` (default) omits an infeasible strict M4;
#'   `"approximate"` constructs an explicitly labelled approximate comparator.
#' @param seed Integer seed. The same seed is passed to reproducible allocators.
#' @param ... Further arguments passed to [allocate_sparse_met()], including
#'   genetic grouping and seed constraints.
#'
#' @return A list with `designs`, `diagnostics`, and the requested capacities.
#' @export
make_met_benchmark_designs <- function(
    treatments, environments, n_test_entries_per_environment,
    target_replications, common_treatments = NULL,
    m4_fallback = c("omit", "approximate"), seed = 1L, ...) {
  m4_fallback <- match.arg(m4_fallback)
  treatments <- as.character(treatments)
  environments <- as.character(environments)
  if (!length(treatments) || anyNA(treatments) ||
      any(!nzchar(treatments)) || anyDuplicated(treatments))
    stop("`treatments` must contain unique, non-missing, non-empty IDs.")
  if (length(environments) < 2L || anyNA(environments) ||
      any(!nzchar(environments)) || anyDuplicated(environments))
    stop("`environments` must contain at least two unique IDs.")
  common_treatments <- intersect(
    unique(as.character(common_treatments)), treatments
  )
  if (!is.numeric(target_replications) ||
      length(target_replications) != 1L ||
      !is.finite(target_replications) ||
      target_replications < 1 ||
      target_replications != as.integer(target_replications))
    stop("`target_replications` must be one positive integer.")
  target_replications <- as.integer(target_replications)
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != as.integer(seed))
    stop("`seed` must be one finite integer.")
  seed <- as.integer(seed)

  k <- .benchmark_named_numeric(
    n_test_entries_per_environment, environments,
    "n_test_entries_per_environment", integer = TRUE, lower = 1
  )
  if (any(k > length(treatments)))
    stop("An environment capacity cannot exceed the number of treatments.")
  if (any(k < length(common_treatments)))
    stop("Every environment must accommodate all common treatments.")

  full <- matrix(
    1L, length(treatments), length(environments),
    dimnames = list(treatments, environments)
  )

  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)

  random_sparse <- matrix(
    0L, length(treatments), length(environments),
    dimnames = list(treatments, environments)
  )
  if (length(common_treatments))
    random_sparse[common_treatments, ] <- 1L
  sparse_pool <- setdiff(treatments, common_treatments)
  for (e in environments) {
    n_add <- k[e] - length(common_treatments)
    if (n_add > 0L)
      random_sparse[sample(sparse_pool, n_add), e] <- 1L
  }

  m3 <- allocate_sparse_met(
    treatments = treatments,
    environments = environments,
    allocation_method = "random_balanced",
    n_test_entries_per_environment = k,
    target_replications = target_replications,
    common_treatments = common_treatments,
    seed = seed, ...
  )

  m4_check <- check_equireplicate_feasibility(
    n_treatments_total = length(treatments),
    n_environments = length(environments),
    n_test_entries_per_environment = unname(k),
    target_replications = target_replications,
    n_common_treatments = length(common_treatments)
  )
  designs <- list(
    full = full,
    random_sparse = random_sparse,
    M3_coverage_first = m3$allocation_matrix
  )
  m4_status <- if (isTRUE(m4_check$feasible)) "strict_exact" else "omitted"
  if (isTRUE(m4_check$feasible)) {
    m4 <- allocate_sparse_met(
      treatments = treatments,
      environments = environments,
      allocation_method = "equireplicate",
      n_test_entries_per_environment = k,
      target_replications = target_replications,
      common_treatments = common_treatments,
      allow_approximate = FALSE, seed = seed, ...
    )
    designs$M4_equireplicate <- m4$allocation_matrix
  } else if (m4_fallback == "approximate") {
    m4 <- allocate_sparse_met(
      treatments = treatments,
      environments = environments,
      allocation_method = "equireplicate",
      n_test_entries_per_environment = k,
      target_replications = target_replications,
      common_treatments = common_treatments,
      allow_approximate = TRUE, seed = seed, ...
    )
    designs$M4_approximate <- m4$allocation_matrix
    m4_status <- "approximate_explicitly_labelled"
  }

  coverage <- vapply(
    designs, function(A) mean(rowSums(A > 0) > 0), numeric(1)
  )
  diagnostics <- data.frame(
    design = names(designs),
    total_cells = vapply(designs, sum, numeric(1)),
    treatment_coverage = unname(coverage),
    min_environments_per_treatment = vapply(
      designs, function(A) min(rowSums(A > 0)), numeric(1)
    ),
    max_environments_per_treatment = vapply(
      designs, function(A) max(rowSums(A > 0)), numeric(1)
    ),
    stringsAsFactors = FALSE
  )

  structure(
    list(
      designs = designs,
      diagnostics = diagnostics,
      capacities = k,
      m4_feasibility = m4_check,
      m4_status = m4_status,
      seed = seed
    ),
    class = "met_benchmark_designs"
  )
}


#' Benchmark MET designs under statistical and operational uncertainty
#'
#' `benchmark_met_designs()` compares allocations using common genetic and
#' residual Monte Carlo draws. Hence, design contrasts are paired: each design
#' is challenged with the same realised genotype-by-environment effects and
#' cell-level standard-normal errors. The function evaluates prediction error
#' variance (PEV), coefficient of determination (CDmean), selection accuracy,
#' realised genetic gain, selection coincidence, rank quality, cost, plot
#' capacity, seed use, and treatment coverage.
#'
#' Environmental covariance candidates, variance-component values, and complete
#' site-loss patterns are crossed into a stress-test grid. No probability is
#' assigned to unsupported environmental kernels; scenario summaries are
#' uniform sensitivity summaries unless `scenario_probabilities` is supplied.
#'
#' @param designs A named list of allocation matrices, a result from
#'   [make_met_benchmark_designs()], or a named list whose elements contain
#'   `allocation_matrix` and optional `reps`.
#' @param G Named additive or hybrid relationship matrix.
#' @param Sigma_E Central named environment covariance matrix.
#' @param Sigma_E_candidates Optional named list of covariance candidates.
#' @param sigma_g2,sigma_e2 Positive variance-component values crossed into the
#'   scenario grid.
#' @param environment_dropout Named list of character vectors. Each vector
#'   identifies sites assumed to fail completely. The default retains all sites.
#' @param scenario_probabilities Optional non-negative probabilities for the
#'   complete crossed scenario grid. They are normalised to sum to one.
#' @param reps Optional named list of replication-count matrices, one per
#'   design. Replication stored inside a design element takes precedence.
#' @param env_efficiency Scalar or environment-named vector in `(0, 1]`.
#' @param cost_per_plot Non-negative scalar or environment-named vector.
#' @param fixed_plot_overhead Non-negative scalar or environment-named vector
#'   for checks, borders, or other physical plots not represented in `reps`.
#' @param environment_plot_capacity Optional scalar or environment-named vector
#'   of physical plot limits.
#' @param seed_available Optional treatment-named inventory, or a data frame
#'   with `Treatment` and `SeedAvailable`.
#' @param seed_required_per_plot Positive scalar, environment-named vector, or
#'   treatment-by-environment matrix.
#' @param minimum_seed_buffer Non-negative scalar or treatment-named vector.
#' @param n_sim Number of paired Monte Carlo replicates.
#' @param select_fraction Selected fraction.
#' @param seed Integer random seed.
#' @param max_dim Dense mixed-model dimension guard.
#'
#' @return A list with `design_summary`, `scenario_results`,
#'   `feasibility`, `scenarios`, and the aligned designs and replication counts.
#' @export
benchmark_met_designs <- function(
    designs, G, Sigma_E = NULL, Sigma_E_candidates = NULL,
    sigma_g2 = 1, sigma_e2 = c(0.5, 1, 2),
    environment_dropout = list(none = character()),
    scenario_probabilities = NULL, reps = NULL, env_efficiency = 1,
    cost_per_plot = 1, fixed_plot_overhead = 0,
    environment_plot_capacity = NULL,
    seed_available = NULL, seed_required_per_plot = 1,
    minimum_seed_buffer = 0, n_sim = 100L, select_fraction = 0.1,
    seed = 1L, max_dim = 6000L) {
  if (inherits(designs, "met_benchmark_designs"))
    designs <- designs$designs
  if (!is.list(designs) || length(designs) < 2L)
    stop("`designs` must be a named list of at least two designs.")
  if (is.null(names(designs)) || any(!nzchar(names(designs))) ||
      anyDuplicated(names(designs)))
    stop("`designs` must have unique, non-empty names.")
  if (!is.null(reps) &&
      (!is.list(reps) || is.null(names(reps)) ||
       !all(names(designs) %in% names(reps))))
    stop("`reps` must be NULL or a named list covering every design.")
  if (!is.numeric(n_sim) || length(n_sim) != 1L || !is.finite(n_sim) ||
      n_sim < 2L || n_sim != as.integer(n_sim))
    stop("`n_sim` must be an integer >= 2.")
  if (!is.numeric(select_fraction) || length(select_fraction) != 1L ||
      !is.finite(select_fraction) || select_fraction <= 0 ||
      select_fraction > 1)
    stop("`select_fraction` must be in (0, 1].")
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != as.integer(seed))
    stop("`seed` must be one finite integer.")

  extracted <- lapply(seq_along(designs), function(i) {
    z <- designs[[i]]
    A <- if (is.list(z) && !is.null(z$allocation_matrix))
      z$allocation_matrix else z
    R <- if (is.list(z) && !is.null(z$reps)) z$reps else
      if (!is.null(reps)) reps[[names(designs)[i]]] else A
    list(allocation = as.matrix(A), reps = as.matrix(R))
  })
  names(extracted) <- names(designs)
  ref_A <- extracted[[1L]]$allocation
  if (is.null(rownames(ref_A)) || is.null(colnames(ref_A)) ||
      anyDuplicated(rownames(ref_A)) || anyDuplicated(colnames(ref_A)))
    stop("Every allocation needs unique treatment and environment dimnames.")
  treatments <- rownames(ref_A)
  environments <- colnames(ref_A)
  for (nm in names(extracted)) {
    z <- extracted[[nm]]
    if (is.null(rownames(z$allocation)) || is.null(colnames(z$allocation)) ||
        !setequal(rownames(z$allocation), treatments) ||
        !setequal(colnames(z$allocation), environments))
      stop("All designs must contain the same named treatments and environments.")
    z$allocation <- z$allocation[treatments, environments, drop = FALSE]
    if (is.null(rownames(z$reps)) || is.null(colnames(z$reps)) ||
        !all(treatments %in% rownames(z$reps)) ||
        !all(environments %in% colnames(z$reps)))
      stop("Every replication matrix must cover all treatments and environments.")
    z$reps <- z$reps[treatments, environments, drop = FALSE]
    if (!is.numeric(z$allocation) || any(!is.finite(z$allocation)) ||
        any(z$allocation < 0))
      stop("Allocation matrices must be finite and non-negative.")
    z$allocation <- 1L * (z$allocation > 0)
    if (!is.numeric(z$reps) || any(!is.finite(z$reps)) ||
        any(z$reps < 0) ||
        any(z$reps[z$allocation > 0] <= 0) ||
        any(z$reps[z$allocation == 0] != 0))
      stop("Replication counts must be positive exactly in allocated cells.")
    extracted[[nm]] <- z
  }

  J <- length(treatments)
  E <- length(environments)
  if (!is.numeric(max_dim) || length(max_dim) != 1L ||
      !is.finite(max_dim) || max_dim < 1L ||
      max_dim != as.integer(max_dim))
    stop("`max_dim` must be one positive integer.")
  if (J * E > as.integer(max_dim))
    stop(
      "The benchmark has ", J * E,
      " genotype-by-environment effects, exceeding `max_dim = ",
      as.integer(max_dim), "`. Increase the guard only when dense-memory ",
      "requirements are acceptable."
    )
  G <- as.matrix(G)
  if (is.null(rownames(G)) || is.null(colnames(G)) ||
      !all(treatments %in% rownames(G)) ||
      !all(treatments %in% colnames(G)))
    stop("`G` must be a named matrix covering every treatment.")
  G <- G[treatments, treatments, drop = FALSE]
  if (!isTRUE(all.equal(G, t(G), tolerance = 1e-8)) ||
      any(!is.finite(G)))
    stop("`G` must be finite and symmetric after alignment.")

  env_efficiency <- .benchmark_named_numeric(
    env_efficiency, environments, "env_efficiency",
    lower = .Machine$double.eps, upper = 1
  )
  cost_per_plot <- .benchmark_named_numeric(
    cost_per_plot, environments, "cost_per_plot", lower = 0
  )
  fixed_plot_overhead <- .benchmark_named_numeric(
    fixed_plot_overhead, environments, "fixed_plot_overhead", lower = 0
  )
  plot_capacity <- if (is.null(environment_plot_capacity)) NULL else
    .benchmark_named_numeric(
      environment_plot_capacity, environments,
      "environment_plot_capacity", lower = 0
    )

  covariances <- .benchmark_covariance_candidates(
    Sigma_E, Sigma_E_candidates, environments
  )
  if (!is.numeric(sigma_g2) || !length(sigma_g2) ||
      any(!is.finite(sigma_g2)) || any(sigma_g2 <= 0) ||
      !is.numeric(sigma_e2) || !length(sigma_e2) ||
      any(!is.finite(sigma_e2)) || any(sigma_e2 <= 0))
    stop("`sigma_g2` and `sigma_e2` must contain finite positive values.")
  if (!is.list(environment_dropout) || !length(environment_dropout))
    stop("`environment_dropout` must be a non-empty list.")
  if (is.null(names(environment_dropout)))
    names(environment_dropout) <- paste0("dropout_", seq_along(environment_dropout))
  if (any(!nzchar(names(environment_dropout))) ||
      anyDuplicated(names(environment_dropout)))
    stop("`environment_dropout` must have unique, non-empty names.")
  environment_dropout <- lapply(environment_dropout, function(x) {
    x <- unique(as.character(x))
    if (any(!x %in% environments))
      stop("Every dropout environment must occur in the designs.")
    if (length(x) == E)
      stop("A dropout scenario cannot remove every environment.")
    x
  })

  scenario_grid <- expand.grid(
    covariance = names(covariances),
    sigma_g2 = as.numeric(sigma_g2),
    sigma_e2 = as.numeric(sigma_e2),
    dropout = names(environment_dropout),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  if (is.null(scenario_probabilities)) {
    scenario_grid$probability <- rep(1 / nrow(scenario_grid),
                                     nrow(scenario_grid))
  } else {
    if (!is.numeric(scenario_probabilities) ||
        length(scenario_probabilities) != nrow(scenario_grid) ||
        any(!is.finite(scenario_probabilities)) ||
        any(scenario_probabilities < 0) ||
        sum(scenario_probabilities) <= 0)
      stop("`scenario_probabilities` must be non-negative, finite, and have ",
           "one value per crossed scenario.")
    scenario_grid$probability <- scenario_probabilities /
      sum(scenario_probabilities)
  }
  scenario_grid$scenario <- paste0("S", seq_len(nrow(scenario_grid)))

  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))
  Zu <- matrix(stats::rnorm(J * E * n_sim), J * E, n_sim)
  Ze <- matrix(stats::rnorm(J * E * n_sim), J * E, n_sim)

  scenario_rows <- vector("list", nrow(scenario_grid) * length(extracted))
  rr <- 0L
  for (s in seq_len(nrow(scenario_grid))) {
    sc <- scenario_grid[s, , drop = FALSE]
    Sig <- covariances[[sc$covariance]]
    covU <- sc$sigma_g2 * kronecker(Sig, G)
    U <- .benchmark_psd_root(covU) %*% Zu
    dropped <- environment_dropout[[sc$dropout]]
    for (nm in names(extracted)) {
      rr <- rr + 1L
      z <- extracted[[nm]]
      A <- z$allocation
      R <- z$reps
      if (length(dropped)) {
        A[, dropped] <- 0L
        R[, dropped] <- 0
      }
      info <- met_information(
        A, G = G, Sigma_E = Sig, sigma_g2 = sc$sigma_g2,
        sigma_e2 = sc$sigma_e2, reps = R,
        env_efficiency = env_efficiency, max_dim = max_dim
      )
      outcome <- .benchmark_paired_outcomes(
        A, R, info, U, Ze, Sig, G, sc$sigma_g2, sc$sigma_e2,
        env_efficiency, select_fraction
      )
      scenario_rows[[rr]] <- data.frame(
        scenario = sc$scenario,
        design = nm,
        covariance = sc$covariance,
        sigma_g2 = sc$sigma_g2,
        sigma_e2 = sc$sigma_e2,
        dropout = sc$dropout,
        probability = sc$probability,
        mean_PEV = info$mean_PEV,
        CDmean = info$CDmean,
        accuracy = outcome$accuracy_mean,
        accuracy_se = outcome$accuracy_se,
        accuracy_ci95_lower = outcome$accuracy_ci95[1L],
        accuracy_ci95_upper = outcome$accuracy_ci95[2L],
        gain = outcome$gain_mean,
        gain_se = outcome$gain_se,
        gain_ci95_lower = outcome$gain_ci95[1L],
        gain_ci95_upper = outcome$gain_ci95[2L],
        common_selected = outcome$common_selected_mean,
        avg_rank = outcome$avg_rank_mean,
        stringsAsFactors = FALSE
      )
    }
  }
  scenario_results <- do.call(rbind, scenario_rows)

  seed_inputs <- .benchmark_seed_inputs(
    seed_available, seed_required_per_plot, minimum_seed_buffer,
    treatments, environments
  )
  feasibility_rows <- lapply(names(extracted), function(nm) {
    z <- extracted[[nm]]
    test_plots_by_environment <- colSums(z$reps)
    plots_by_environment <- test_plots_by_environment + fixed_plot_overhead
    capacity_ok <- is.null(plot_capacity) ||
      all(plots_by_environment <= plot_capacity + 1e-8)
    used_seed <- if (is.null(seed_inputs)) rep(NA_real_, J) else
      rowSums(z$reps * seed_inputs$required)
    seed_ok <- is.null(seed_inputs) ||
      all(used_seed <= seed_inputs$spendable + 1e-8)
    coverage <- rowSums(z$allocation) > 0
    data.frame(
      design = nm,
      test_plots = sum(z$reps),
      fixed_plot_overhead = sum(fixed_plot_overhead),
      plots = sum(plots_by_environment),
      cost = sum(plots_by_environment * cost_per_plot),
      treatment_coverage = mean(coverage),
      uncovered_treatments = sum(!coverage),
      min_environments_per_treatment = min(rowSums(z$allocation)),
      max_site_plots = max(plots_by_environment),
      plot_capacity_feasible = capacity_ok,
      seed_feasible = seed_ok,
      feasible = all(coverage) && capacity_ok && seed_ok,
      seed_used = if (is.null(seed_inputs)) NA_real_ else sum(used_seed),
      minimum_seed_remaining = if (is.null(seed_inputs)) NA_real_ else
        min(seed_inputs$available - used_seed),
      stringsAsFactors = FALSE
    )
  })
  feasibility <- do.call(rbind, feasibility_rows)

  summary_rows <- lapply(names(extracted), function(nm) {
    x <- scenario_results[scenario_results$design == nm, , drop = FALSE]
    p <- x$probability / sum(x$probability)
    data.frame(
      design = nm,
      mean_PEV = sum(p * x$mean_PEV),
      worst_mean_PEV = max(x$mean_PEV),
      CDmean = sum(p * x$CDmean),
      worst_CDmean = min(x$CDmean),
      accuracy = sum(p * x$accuracy),
      worst_accuracy = min(x$accuracy),
      gain = sum(p * x$gain),
      worst_gain = min(x$gain),
      common_selected = sum(p * x$common_selected),
      worst_common_selected = min(x$common_selected),
      avg_rank = sum(p * x$avg_rank),
      worst_avg_rank = max(x$avg_rank),
      stringsAsFactors = FALSE
    )
  })
  design_summary <- merge(
    do.call(rbind, summary_rows), feasibility, by = "design", sort = FALSE
  )
  design_summary$pareto_accuracy <- .benchmark_pareto_flag(
    design_summary$accuracy, design_summary$cost, design_summary$feasible
  )
  design_summary$pareto_gain <- .benchmark_pareto_flag(
    design_summary$gain, design_summary$cost, design_summary$feasible
  )
  cost_range <- range(design_summary$cost)
  design_summary$cost_minimisation_percent <- if (diff(cost_range) > 0)
    100 * (cost_range[2L] - design_summary$cost) / diff(cost_range) else 0
  full_idx <- match("full", design_summary$design)
  if (!is.na(full_idx)) {
    full_cost <- design_summary$cost[full_idx]
    full_accuracy <- design_summary$accuracy[full_idx]
    design_summary$cost_saving_vs_full_percent <-
      100 * (full_cost - design_summary$cost) / full_cost
    design_summary$accuracy_vs_full_percent <- if (
      is.finite(full_accuracy) && abs(full_accuracy) > 1e-12
    ) 100 * design_summary$accuracy / full_accuracy else NA_real_
  } else {
    design_summary$cost_saving_vs_full_percent <- NA_real_
    design_summary$accuracy_vs_full_percent <- NA_real_
  }
  design_summary$robust_accuracy_rank <- rank(
    -design_summary$worst_accuracy, ties.method = "min"
  )

  list(
    design_summary = design_summary,
    scenario_results = scenario_results,
    feasibility = feasibility,
    scenarios = scenario_grid,
    designs = lapply(extracted, `[[`, "allocation"),
    reps = lapply(extracted, `[[`, "reps"),
    paired_simulation = TRUE,
    n_sim = as.integer(n_sim),
    seed = as.integer(seed)
  )
}


#' Benchmark environmental covariance models against response evidence
#'
#' This function compares independence, each separate environmental kernel, an
#' explicitly labelled equal-weight single-kernel comparator, and the
#' historically calibrated covariance when response evidence is supplied. The
#' equal-weight matrix is included only as a benchmark; it is never adopted as
#' the package's no-history covariance.
#'
#' @param kernels Named environmental kernels.
#' @param target Optional empirical covariance/correlation matrix, or named list
#'   of matrices.
#' @param historical Optional historical long-format MET data.
#' @param truth Optional known covariance used only to score recovery in a
#'   simulation study.
#' @param genotype_col,environment_col,value_col,year_col Historical column
#'   names passed to [calibrate_environment_covariance()].
#' @param interaction_policy,n_boot,seed Passed to
#'   [calibrate_environment_covariance()].
#' @param ... Further calibration arguments.
#'
#' @return A list containing model scores, target-specific scores, the
#'   calibration result, candidates, and blocked cross-validation diagnostics.
#' @export
benchmark_environment_models <- function(
    kernels, target = NULL, historical = NULL, truth = NULL,
    genotype_col = "genotype", environment_col = "environment",
    value_col = "value", year_col = NULL,
    interaction_policy = c("evidence", "exclude", "include"),
    n_boot = 200L, seed = 1L, ...) {
  interaction_policy <- match.arg(interaction_policy)
  K <- .validate_environment_kernels(kernels)
  envs <- rownames(K[[1L]])
  calibration <- calibrate_environment_covariance(
    kernels = K, target = target, historical = historical,
    genotype_col = genotype_col, environment_col = environment_col,
    value_col = value_col, year_col = year_col,
    interaction_policy = interaction_policy,
    n_boot = n_boot, seed = seed, ...
  )
  identity <- diag(length(envs))
  dimnames(identity) <- list(envs, envs)
  equal_weight <- combine_environment_kernels(
    K, weights = stats::setNames(rep(1, length(K)), names(K))
  )
  candidates <- c(
    list(independent = identity, single_equal_weight = equal_weight),
    stats::setNames(K, paste0("separate_", names(K)))
  )
  if (calibration$status == "historically_calibrated")
    candidates$historically_calibrated <- calibration$Sigma_E

  references <- NULL
  reference_basis <- "none"
  if (!is.null(truth)) {
    references <- list(known_truth = .benchmark_align_covariance(truth, envs))
    reference_basis <- "known_truth"
  } else if (!is.null(calibration$target)) {
    if (is.list(target) && !is.matrix(target)) {
      references <- lapply(target, .benchmark_align_covariance, envs = envs)
      if (is.null(names(references)))
        names(references) <- paste0("target_", seq_along(references))
    } else {
      references <- list(empirical_target = calibration$target)
    }
    reference_basis <- "empirical_response"
  }

  if (is.null(references)) {
    score_by_target <- do.call(rbind, lapply(names(candidates), function(nm)
      data.frame(
        model = nm, reference = NA_character_, n_pairs = NA_integer_,
        rmse = NA_real_, mae = NA_real_, correlation = NA_real_,
        relative_frobenius_error = NA_real_, stringsAsFactors = FALSE
      )))
  } else {
    score_by_target <- do.call(rbind, lapply(names(candidates), function(nm)
      do.call(rbind, lapply(names(references), function(rn)
        .benchmark_covariance_score(candidates[[nm]], references[[rn]],
                                    nm, rn)))))
  }
  scores <- do.call(rbind, lapply(split(score_by_target,
                                        score_by_target$model), function(x)
    data.frame(
      model = x$model[1L],
      n_references = sum(!is.na(x$reference)),
      rmse = mean(x$rmse, na.rm = TRUE),
      mae = mean(x$mae, na.rm = TRUE),
      correlation = mean(x$correlation, na.rm = TRUE),
      relative_frobenius_error =
        mean(x$relative_frobenius_error, na.rm = TRUE),
      stringsAsFactors = FALSE
    )))
  for (nm in c("rmse", "mae", "correlation", "relative_frobenius_error"))
    scores[[nm]][is.nan(scores[[nm]])] <- NA_real_
  scores$rmse_rank <- rank(scores$rmse, na.last = "keep", ties.method = "min")

  list(
    scores = scores,
    score_by_target = score_by_target,
    reference_basis = reference_basis,
    calibration = calibration,
    candidates = candidates,
    cross_validation = calibration$cross_validation,
    interaction_evidence = calibration$interaction_evidence
  )
}


#' Simulate recovery of a known environmental covariance
#'
#' `simulate_environment_model_benchmark()` generates multi-year
#' genotype-by-environment responses under a known mixture of separate kernels.
#' It then fits the standard calibration workflow and reports covariance and
#' component-weight recovery. This is a demonstration of identifiability and
#' predictive validation, not evidence about a specific breeding programme.
#'
#' @param kernels Named environmental kernels.
#' @param true_weights Named non-negative weights for any kernels and
#'   `identity`. Omitted components receive zero weight.
#' @param G Optional genotype relationship matrix. If omitted, unrelated
#'   genotypes are simulated.
#' @param n_genotypes,n_years Positive integers.
#' @param sigma_g2,sigma_e2 Positive genetic and residual variances.
#' @param reps_per_cell Positive number of plot replicates represented by each
#'   simulated genotype-environment mean.
#' @param missing_rate Fraction of response cells set to missing completely at
#'   random.
#' @param interaction_policy,n_boot,seed Passed to the calibration benchmark.
#' @param ... Further arguments passed to [benchmark_environment_models()].
#'
#' @return A list containing the known covariance, simulated historical data,
#'   true and estimated weights, weight-recovery diagnostics, and the complete
#'   environmental-model benchmark.
#' @export
simulate_environment_model_benchmark <- function(
    kernels, true_weights, G = NULL, n_genotypes = 200L, n_years = 8L,
    sigma_g2 = 1, sigma_e2 = 1, reps_per_cell = 2,
    missing_rate = 0, interaction_policy = c("evidence", "exclude", "include"),
    n_boot = 200L, seed = 1L, ...) {
  interaction_policy <- match.arg(interaction_policy)
  K <- .validate_environment_kernels(kernels)
  envs <- rownames(K[[1L]])
  component_names <- c(names(K), "identity")
  if (!is.numeric(true_weights) || is.null(names(true_weights)) ||
      any(!names(true_weights) %in% component_names) ||
      any(!is.finite(true_weights)) || any(true_weights < 0) ||
      sum(true_weights) <= 0)
    stop("`true_weights` must be a named, non-negative vector over kernel ",
         "names and optionally `identity`.")
  w_true <- stats::setNames(rep(0, length(component_names)), component_names)
  w_true[names(true_weights)] <- true_weights
  w_true <- w_true / sum(w_true)
  identity <- diag(length(envs))
  dimnames(identity) <- list(envs, envs)
  all_components <- c(K, list(identity = identity))
  Sigma_true <- .normalise_environment_kernel(.bend_pd(
    Reduce(`+`, Map(`*`, all_components, as.list(w_true)))
  ))

  whole <- function(x, lower)
    is.numeric(x) && length(x) == 1L && is.finite(x) &&
      x >= lower && x == as.integer(x)
  if (!whole(n_genotypes, 3L) || !whole(n_years, 2L))
    stop("`n_genotypes` must be >= 3 and `n_years` must be >= 2.")
  if (!is.numeric(sigma_g2) || length(sigma_g2) != 1L ||
      !is.finite(sigma_g2) || sigma_g2 <= 0 ||
      !is.numeric(sigma_e2) || length(sigma_e2) != 1L ||
      !is.finite(sigma_e2) || sigma_e2 <= 0 ||
      !is.numeric(reps_per_cell) || length(reps_per_cell) != 1L ||
      !is.finite(reps_per_cell) || reps_per_cell <= 0)
    stop("Variances and `reps_per_cell` must be finite and positive.")
  if (!is.numeric(missing_rate) || length(missing_rate) != 1L ||
      !is.finite(missing_rate) || missing_rate < 0 || missing_rate >= 0.8)
    stop("`missing_rate` must be one value in [0, 0.8).")

  if (is.null(G)) {
    G <- diag(as.integer(n_genotypes))
    genotypes <- paste0("G", seq_len(n_genotypes))
    dimnames(G) <- list(genotypes, genotypes)
  } else {
    G <- as.matrix(G)
    if (nrow(G) != ncol(G) || is.null(rownames(G)) ||
        is.null(colnames(G)) || !setequal(rownames(G), colnames(G)) ||
        !isTRUE(all.equal(G, t(G), tolerance = 1e-8)) ||
        any(!is.finite(G)))
      stop("`G` must be a finite, symmetric, named square matrix.")
    genotypes <- rownames(G)
    G <- G[genotypes, genotypes, drop = FALSE]
    n_genotypes <- length(genotypes)
  }

  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))
  root <- .benchmark_psd_root(sigma_g2 * kronecker(Sigma_true, G))
  template <- expand.grid(
    genotype = genotypes, environment = envs,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  years <- lapply(seq_len(n_years), function(y) {
    u <- as.numeric(root %*% stats::rnorm(n_genotypes * length(envs)))
    value <- u + stats::rnorm(length(u), sd = sqrt(sigma_e2 / reps_per_cell))
    if (missing_rate > 0) {
      n_miss <- floor(missing_rate * length(value))
      if (n_miss) value[sample.int(length(value), n_miss)] <- NA_real_
    }
    data.frame(
      genotype = template$genotype,
      environment = template$environment,
      year = paste0("Y", y),
      value = value,
      stringsAsFactors = FALSE
    )
  })
  historical <- do.call(rbind, years)
  benchmark <- benchmark_environment_models(
    kernels = K, historical = historical, truth = Sigma_true,
    genotype_col = "genotype", environment_col = "environment",
    value_col = "value", year_col = "year",
    interaction_policy = interaction_policy,
    n_boot = n_boot, seed = seed, ...
  )
  estimated <- benchmark$calibration$weights
  weight_recovery <- data.frame(
    component = component_names,
    true_weight = unname(w_true),
    estimated_weight = if (is.null(estimated)) NA_real_ else
      unname(estimated[component_names]),
    stringsAsFactors = FALSE
  )
  weight_recovery$absolute_error <- abs(
    weight_recovery$estimated_weight - weight_recovery$true_weight
  )

  list(
    Sigma_E_true = Sigma_true,
    historical = historical,
    true_weights = w_true,
    estimated_weights = estimated,
    weight_recovery = weight_recovery,
    benchmark = benchmark,
    seed = as.integer(seed)
  )
}


#' Stress-test environmental kernels under missing covariates
#'
#' `benchmark_environment_missingness()` repeatedly masks observed weather,
#' soil, management, and geographic cells, rebuilds the kernels through the
#' package quality-control and imputation workflow, and scores the resulting
#' covariance candidates against response evidence or a known truth.
#'
#' @param weather,soil,management,geography Environmental blocks accepted by
#'   [build_environment_kernels()].
#' @param missing_rates Fractions of eligible covariate cells masked completely
#'   at random. Include zero for the unperturbed reference.
#' @param n_repeats Number of masks at each missingness level.
#' @param target,historical,truth Evidence inputs passed to
#'   [benchmark_environment_models()].
#' @param build_args Named list of additional arguments passed to
#'   [build_environment_kernels()].
#' @param benchmark_args Named list of additional arguments passed to
#'   [benchmark_environment_models()]. The default disables bootstrapping inside
#'   each perturbation.
#' @param seed Integer random seed.
#'
#' @return A list with run-level model scores, missingness-level summaries,
#'   masking diagnostics, a kernel-rebuild warning ledger, and the rebuilt
#'   kernel objects. Expected quality-control warnings caused by deliberate
#'   masking are captured in the ledger rather than emitted to the caller.
#' @export
benchmark_environment_missingness <- function(
    weather = NULL, soil = NULL, management = NULL, geography = NULL,
    missing_rates = c(0, 0.05, 0.10, 0.20), n_repeats = 20L,
    target = NULL, historical = NULL, truth = NULL,
    build_args = list(missing_action = "impute", impute = "median"),
    benchmark_args = list(n_boot = 0L), seed = 1L) {
  blocks <- Filter(Negate(is.null), list(
    weather = weather, soil = soil, management = management,
    geography = geography
  ))
  if (!length(blocks))
    stop("Supply at least one environmental data block.")
  if (!is.numeric(missing_rates) || !length(missing_rates) ||
      any(!is.finite(missing_rates)) ||
      any(missing_rates < 0 | missing_rates >= 1))
    stop("`missing_rates` must contain values in [0, 1).")
  missing_rates <- sort(unique(as.numeric(missing_rates)))
  if (!is.numeric(n_repeats) || length(n_repeats) != 1L ||
      !is.finite(n_repeats) || n_repeats < 1L ||
      n_repeats != as.integer(n_repeats))
    stop("`n_repeats` must be one positive integer.")
  if (!is.list(build_args) ||
      (length(build_args) &&
       (is.null(names(build_args)) || any(!nzchar(names(build_args))))) ||
      !is.list(benchmark_args) ||
      (length(benchmark_args) &&
       (is.null(names(benchmark_args)) || any(!nzchar(names(benchmark_args))))))
    stop("`build_args` and `benchmark_args` must be named lists.")
  reserved_build <- intersect(
    names(build_args), c("weather", "soil", "management", "geography")
  )
  if (length(reserved_build))
    stop("Environmental blocks must be supplied directly, not in `build_args`.")
  reserved_benchmark <- intersect(
    names(benchmark_args), c("kernels", "target", "historical", "truth")
  )
  if (length(reserved_benchmark))
    stop("Evidence inputs must be supplied directly, not in `benchmark_args`.")
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != as.integer(seed))
    stop("`seed` must be one finite integer.")

  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))

  score_rows <- list()
  mask_rows <- list()
  warning_rows <- list()
  rebuilt <- list()
  row_i <- 0L
  for (rate in missing_rates) {
    repeats <- if (rate == 0) 1L else as.integer(n_repeats)
    for (r in seq_len(repeats)) {
      row_i <- row_i + 1L
      run_id <- paste0("missing_", format(rate, trim = TRUE),
                       "_rep_", r)
      masked <- lapply(blocks, .benchmark_mask_environment_block, rate = rate)
      n_eligible <- sum(vapply(masked, attr, numeric(1), "n_eligible"))
      n_masked <- sum(vapply(masked, attr, numeric(1), "n_masked"))
      masked <- lapply(masked, function(x) {
        attr(x, "n_eligible") <- NULL
        attr(x, "n_masked") <- NULL
        x
      })
      build_call <- c(
        list(
          weather = masked$weather, soil = masked$soil,
          management = masked$management, geography = masked$geography
        ),
        build_args
      )
      build_warnings <- character()
      kb <- withCallingHandlers(
        do.call(build_environment_kernels, build_call),
        warning = function(w) {
          build_warnings <<- c(build_warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
      build_warnings <- unique(build_warnings)
      bench_call <- c(
        list(
          kernels = kb$kernels, target = target,
          historical = historical, truth = truth
        ),
        benchmark_args
      )
      bm <- do.call(benchmark_environment_models, bench_call)
      score <- bm$scores
      score$missing_rate <- rate
      score$replicate <- r
      score$run <- run_id
      score_rows[[row_i]] <- score
      mask_rows[[row_i]] <- data.frame(
        run = run_id, missing_rate = rate, replicate = r,
        eligible_cells = n_eligible, masked_cells = n_masked,
        realised_missing_rate = if (n_eligible)
          n_masked / n_eligible else 0,
        qc_warning_count = length(build_warnings),
        qc_warnings = paste(build_warnings, collapse = " | "),
        stringsAsFactors = FALSE
      )
      warning_rows[[row_i]] <- if (length(build_warnings)) {
        data.frame(
          run = run_id,
          warning = build_warnings,
          stringsAsFactors = FALSE
        )
      } else {
        data.frame(
          run = character(),
          warning = character(),
          stringsAsFactors = FALSE
        )
      }
      rebuilt[[run_id]] <- list(kernels = kb, benchmark = bm)
    }
  }
  run_scores <- do.call(rbind, score_rows)
  masking <- do.call(rbind, mask_rows)
  warning_ledger <- do.call(rbind, warning_rows)
  split_scores <- split(
    run_scores, interaction(run_scores$model, run_scores$missing_rate,
                            drop = TRUE)
  )
  summary <- do.call(rbind, lapply(split_scores, function(x)
    data.frame(
      model = x$model[1L],
      missing_rate = x$missing_rate[1L],
      n_repeats = nrow(x),
      rmse_mean = mean(x$rmse, na.rm = TRUE),
      rmse_sd = stats::sd(x$rmse, na.rm = TRUE),
      correlation_mean = mean(x$correlation, na.rm = TRUE),
      correlation_sd = stats::sd(x$correlation, na.rm = TRUE),
      relative_frobenius_error_mean =
        mean(x$relative_frobenius_error, na.rm = TRUE),
      stringsAsFactors = FALSE
    )))
  numeric_summary <- setdiff(names(summary),
                             c("model", "missing_rate", "n_repeats"))
  for (nm in numeric_summary)
    summary[[nm]][is.nan(summary[[nm]])] <- NA_real_

  list(
    run_scores = run_scores,
    summary = summary,
    masking = masking,
    warnings = warning_ledger,
    rebuilt = rebuilt,
    seed = as.integer(seed)
  )
}


#' Summarise stability of repeated MET design decisions
#'
#' @param common_sets Named list of common-treatment vectors from repeated
#'   perturbations, bootstrap samples, or variance scenarios.
#' @param site_capacities Optional named list or matrix of site-capacity
#'   recommendations. Rows are runs and columns are environments.
#' @param allocations Optional named list of aligned allocation matrices.
#'
#' @return Selection frequencies and pairwise Jaccard stability for common sets,
#'   plus capacity and allocation stability when supplied.
#' @export
summarize_design_stability <- function(
    common_sets, site_capacities = NULL, allocations = NULL) {
  if (!is.list(common_sets) || length(common_sets) < 2L)
    stop("`common_sets` must contain at least two repeated selections.")
  if (is.null(names(common_sets)))
    names(common_sets) <- paste0("run_", seq_along(common_sets))
  common_sets <- lapply(common_sets, function(x)
    unique(as.character(x[!is.na(x) & nzchar(as.character(x))])))
  universe <- sort(unique(unlist(common_sets, use.names = FALSE)))
  frequency <- data.frame(
    Treatment = universe,
    selection_frequency = vapply(
      universe, function(g) mean(vapply(
        common_sets, function(s) g %in% s, logical(1)
      )), numeric(1)
    ),
    n_selected = vapply(
      universe, function(g) sum(vapply(
        common_sets, function(s) g %in% s, logical(1)
      )), integer(1)
    ),
    stringsAsFactors = FALSE
  )
  jaccard <- .benchmark_pairwise_jaccard(common_sets)

  capacity_summary <- NULL
  capacity_runs <- NULL
  if (!is.null(site_capacities)) {
    capacity_runs <- if (is.list(site_capacities)) {
      envs <- unique(unlist(lapply(site_capacities, names),
                            use.names = FALSE))
      if (!length(envs))
        stop("List elements in `site_capacities` must be named.")
      out <- matrix(NA_real_, length(site_capacities), length(envs),
                    dimnames = list(names(site_capacities), envs))
      for (i in seq_along(site_capacities))
        out[i, names(site_capacities[[i]])] <- site_capacities[[i]]
      out
    } else {
      as.matrix(site_capacities)
    }
    if (!is.numeric(capacity_runs) || any(!is.finite(capacity_runs)) ||
        is.null(colnames(capacity_runs)))
      stop("`site_capacities` must be finite numeric values with site names.")
    capacity_summary <- data.frame(
      Environment = colnames(capacity_runs),
      median = apply(capacity_runs, 2L, stats::median),
      IQR = apply(capacity_runs, 2L, stats::IQR),
      minimum = apply(capacity_runs, 2L, min),
      maximum = apply(capacity_runs, 2L, max),
      coefficient_of_variation = apply(capacity_runs, 2L, function(x) {
        m <- mean(x)
        if (m == 0) NA_real_ else stats::sd(x) / m
      }),
      stringsAsFactors = FALSE
    )
  }

  allocation_frequency <- NULL
  allocation_jaccard <- NULL
  if (!is.null(allocations)) {
    if (!is.list(allocations) || length(allocations) < 2L)
      stop("`allocations` must contain at least two matrices.")
    ref <- as.matrix(allocations[[1L]])
    if (is.null(rownames(ref)) || is.null(colnames(ref)))
      stop("Allocation matrices require dimnames.")
    cells <- lapply(allocations, function(A) {
      A <- as.matrix(A)
      if (is.null(rownames(A)) || is.null(colnames(A)) ||
          !setequal(rownames(A), rownames(ref)) ||
          !setequal(colnames(A), colnames(ref)))
        stop("All allocation matrices must have aligned dimnames.")
      A <- A[rownames(ref), colnames(ref), drop = FALSE]
      which(A > 0)
    })
    allocation_jaccard <- .benchmark_pairwise_jaccard(cells)
    arr <- simplify2array(lapply(allocations, function(A) {
      A <- as.matrix(A)[rownames(ref), colnames(ref), drop = FALSE]
      1 * (A > 0)
    }))
    allocation_frequency <- apply(arr, c(1L, 2L), mean)
  }

  list(
    common_frequency = frequency,
    common_jaccard = jaccard,
    mean_common_jaccard = mean(jaccard[upper.tri(jaccard)]),
    minimum_common_jaccard = min(jaccard[upper.tri(jaccard)]),
    capacity_summary = capacity_summary,
    capacity_runs = capacity_runs,
    allocation_frequency = allocation_frequency,
    allocation_jaccard = allocation_jaccard
  )
}


.benchmark_named_numeric <- function(
    x, ids, nm, integer = FALSE, lower = -Inf, upper = Inf) {
  if (!is.numeric(x) || !length(x) || any(!is.finite(x)))
    stop("`", nm, "` must contain finite numeric values.")
  if (length(x) == 1L) {
    x <- rep(x, length(ids))
    names(x) <- ids
  } else if (!is.null(names(x))) {
    if (anyNA(names(x)) || any(!nzchar(names(x))) || anyDuplicated(names(x)))
      stop("Named `", nm, "` must have unique, non-empty names.")
    if (!all(ids %in% names(x)))
      stop("Named `", nm, "` must cover all required identifiers.")
    x <- x[ids]
  } else if (length(x) == length(ids)) {
    names(x) <- ids
  } else {
    stop("`", nm, "` must be scalar, aligned, or named.")
  }
  if (any(x < lower | x > upper))
    stop("`", nm, "` contains values outside the allowed range.")
  if (integer && any(abs(x - round(x)) > 1e-8))
    stop("`", nm, "` must contain integers.")
  if (integer) x <- as.integer(round(x))
  stats::setNames(x, ids)
}


.benchmark_align_covariance <- function(M, envs) {
  M <- as.matrix(M)
  if (!is.numeric(M) || nrow(M) != ncol(M) ||
      is.null(rownames(M)) || is.null(colnames(M)) ||
      !all(envs %in% rownames(M)) || !all(envs %in% colnames(M)))
    stop("Every environmental covariance must be a named square matrix ",
         "covering all environments.")
  M <- M[envs, envs, drop = FALSE]
  if (any(!is.finite(M)) ||
      !isTRUE(all.equal(M, t(M), tolerance = 1e-8)) ||
      any(diag(M) <= 0))
    stop("Environmental covariance matrices must be finite, symmetric, and ",
         "have positive diagonal elements.")
  ee <- eigen((M + t(M)) / 2, symmetric = TRUE,
              only.values = TRUE)$values
  if (min(ee) < -1e-8 * max(1, max(abs(ee))))
    stop("Environmental covariance matrices must be positive semidefinite.")
  M
}


.benchmark_covariance_candidates <- function(Sigma_E, candidates, envs) {
  if (is.null(candidates)) {
    if (is.null(Sigma_E)) {
      M <- diag(length(envs))
      dimnames(M) <- list(envs, envs)
      return(list(independent = M))
    }
    return(list(central = .benchmark_align_covariance(Sigma_E, envs)))
  }
  if (!is.list(candidates) || !length(candidates))
    stop("`Sigma_E_candidates` must be a non-empty named list.")
  if (is.null(names(candidates)) || any(!nzchar(names(candidates))) ||
      anyDuplicated(names(candidates)))
    stop("`Sigma_E_candidates` needs unique, non-empty names.")
  out <- lapply(candidates, .benchmark_align_covariance, envs = envs)
  if (!is.null(Sigma_E) && !"central" %in% names(out))
    out <- c(
      list(central = .benchmark_align_covariance(Sigma_E, envs)),
      out
    )
  out
}


.benchmark_psd_root <- function(M) {
  M <- (as.matrix(M) + t(as.matrix(M))) / 2
  z <- tryCatch(t(chol(M)), error = function(e) NULL)
  if (!is.null(z)) return(z)
  ee <- eigen(M, symmetric = TRUE)
  scale <- max(1, max(abs(ee$values)))
  if (min(ee$values) < -1e-8 * scale)
    stop("A benchmark covariance is not positive semidefinite.")
  ee$vectors %*% (t(ee$vectors) * sqrt(pmax(ee$values, 0)))
}


.benchmark_summary_ci <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x))
    return(list(mean = NA_real_, sd = NA_real_, se = NA_real_,
                ci = c(NA_real_, NA_real_)))
  s <- if (length(x) > 1L) stats::sd(x) else 0
  se <- s / sqrt(length(x))
  mult <- if (length(x) > 1L) stats::qt(0.975, length(x) - 1L) else 0
  list(mean = mean(x), sd = s, se = se,
       ci = mean(x) + c(-1, 1) * mult * se)
}


.benchmark_paired_outcomes <- function(
    A, R, info, U, Ze, Sigma_E, G, sigma_g2, sigma_e2,
    env_efficiency, select_fraction) {
  J <- nrow(A)
  E <- ncol(A)
  n_sim <- ncol(U)
  d <- as.numeric(sweep(R, 2L, env_efficiency / sigma_e2, `*`))
  present <- as.numeric(A) > 0
  pres_idx <- which(present)
  env_of <- rep(seq_len(E), each = J)
  Cinv <- tryCatch(solve(info$C_uu), error = function(e)
    .pinv_sym_dense(info$C_uu))
  w <- rep(1 / E, E)
  n_sel <- max(1L, floor(select_fraction * J))
  acc <- gain <- common <- avg_rank <- numeric(n_sim)
  for (s in seq_len(n_sim)) {
    u <- U[, s]
    true_bv <- as.numeric(matrix(u, J, E) %*% w)
    ybar <- numeric(J * E)
    ybar[pres_idx] <- u[pres_idx] +
      Ze[pres_idx, s] / sqrt(d[pres_idx])
    rhs <- numeric(J * E)
    for (e in seq_len(E)) {
      idx <- which(env_of == e & present)
      if (!length(idx)) next
      mu <- sum(d[idx] * ybar[idx]) / sum(d[idx])
      rhs[idx] <- d[idx] * (ybar[idx] - mu)
    }
    pred <- as.numeric(matrix(Cinv %*% rhs, J, E) %*% w)
    acc[s] <- suppressWarnings(stats::cor(pred, true_bv))
    selected <- order(pred, decreasing = TRUE)[seq_len(n_sel)]
    best <- order(true_bv, decreasing = TRUE)[seq_len(n_sel)]
    gain[s] <- mean(true_bv[selected]) - mean(true_bv)
    common[s] <- length(intersect(selected, best)) / n_sel
    avg_rank[s] <- mean(rank(-true_bv, ties.method = "average")[selected])
  }
  a <- .benchmark_summary_ci(acc)
  g <- .benchmark_summary_ci(gain)
  list(
    accuracy_mean = a$mean, accuracy_sd = a$sd, accuracy_se = a$se,
    accuracy_ci95 = pmin(1, pmax(-1, a$ci)),
    gain_mean = g$mean, gain_sd = g$sd, gain_se = g$se,
    gain_ci95 = g$ci,
    common_selected_mean = mean(common),
    avg_rank_mean = mean(avg_rank),
    n_selected = n_sel
  )
}


.benchmark_seed_inputs <- function(
    available, required, buffer, treatments, environments) {
  if (is.null(available)) return(NULL)
  if (is.data.frame(available)) {
    .validate_cols(available, c("Treatment", "SeedAvailable"),
                   "seed_available")
    available <- stats::setNames(
      as.numeric(available$SeedAvailable),
      as.character(available$Treatment)
    )
  }
  available <- .benchmark_named_numeric(
    available, treatments, "seed_available", lower = 0
  )
  buffer <- .benchmark_named_numeric(
    buffer, treatments, "minimum_seed_buffer", lower = 0
  )
  if (any(buffer > available))
    stop("`minimum_seed_buffer` cannot exceed available seed.")
  if (is.matrix(required)) {
    if (is.null(rownames(required)) || is.null(colnames(required)) ||
        !all(treatments %in% rownames(required)) ||
        !all(environments %in% colnames(required)))
      stop("Matrix `seed_required_per_plot` must cover treatments and sites.")
    required <- required[treatments, environments, drop = FALSE]
    if (!is.numeric(required) || any(!is.finite(required)) ||
        any(required <= 0))
      stop("`seed_required_per_plot` must be finite and positive.")
  } else {
    req_env <- .benchmark_named_numeric(
      required, environments, "seed_required_per_plot",
      lower = .Machine$double.eps
    )
    required <- matrix(
      rep(req_env, each = length(treatments)),
      length(treatments), length(environments),
      dimnames = list(treatments, environments)
    )
  }
  list(
    available = available,
    spendable = available - buffer,
    buffer = buffer,
    required = required
  )
}


.benchmark_covariance_score <- function(candidate, reference,
                                        model, reference_name) {
  candidate <- .benchmark_align_covariance(candidate, rownames(reference))
  reference <- .benchmark_align_covariance(reference, rownames(reference))
  idx <- upper.tri(reference) & is.finite(reference) & is.finite(candidate)
  a <- candidate[idx]
  b <- reference[idx]
  cor_ab <- if (length(a) > 2L && stats::sd(a) > 0 && stats::sd(b) > 0)
    suppressWarnings(stats::cor(a, b)) else NA_real_
  denom <- sqrt(sum(b^2))
  data.frame(
    model = model,
    reference = reference_name,
    n_pairs = length(a),
    rmse = sqrt(mean((a - b)^2)),
    mae = mean(abs(a - b)),
    correlation = cor_ab,
    relative_frobenius_error = if (denom > 0)
      sqrt(sum((a - b)^2)) / denom else NA_real_,
    stringsAsFactors = FALSE
  )
}


.benchmark_pairwise_jaccard <- function(sets) {
  n <- length(sets)
  nm <- names(sets)
  if (is.null(nm)) nm <- paste0("run_", seq_len(n))
  out <- matrix(1, n, n, dimnames = list(nm, nm))
  if (n > 1L) for (i in seq_len(n - 1L)) for (j in (i + 1L):n) {
    union_n <- length(union(sets[[i]], sets[[j]]))
    value <- if (union_n == 0L) 1 else
      length(intersect(sets[[i]], sets[[j]])) / union_n
    out[i, j] <- out[j, i] <- value
  }
  out
}


.benchmark_pareto_flag <- function(benefit, cost, feasible = NULL) {
  if (is.null(feasible)) feasible <- rep(TRUE, length(benefit))
  if (length(feasible) != length(benefit) || anyNA(feasible))
    stop("Pareto feasibility indicators must be complete and aligned.")
  valid <- is.finite(benefit) & is.finite(cost) & as.logical(feasible)
  out <- rep(FALSE, length(benefit))
  if (any(valid))
    out[valid] <- .pareto_flag(benefit[valid], cost[valid])
  out
}


.benchmark_mask_environment_block <- function(x, rate) {
  d <- as.data.frame(x, check.names = FALSE, stringsAsFactors = FALSE)
  original_rownames <- rownames(x)
  eligible_columns <- setdiff(names(d), "environment")
  n_eligible <- nrow(d) * length(eligible_columns)
  n_mask <- floor(rate * n_eligible)
  if (n_mask > 0L) {
    chosen <- sample.int(n_eligible, n_mask)
    row_index <- ((chosen - 1L) %% nrow(d)) + 1L
    col_index <- ((chosen - 1L) %/% nrow(d)) + 1L
    for (i in seq_along(chosen))
      d[row_index[i], eligible_columns[col_index[i]]] <- NA
  }
  if (!is.null(original_rownames))
    rownames(d) <- original_rownames
  attr(d, "n_eligible") <- n_eligible
  attr(d, "n_masked") <- n_mask
  d
}
