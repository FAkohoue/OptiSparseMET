#' Criterion-driven optimisation of a sparse MET allocation
#'
#' @description
#' Searches for the genotype-by-environment allocation that maximises a
#' [design_objective()] -- any weighted combination of statistical quality
#' (reliability / CDmean), expected genetic gain (single-trait or multi-trait
#' index), and resource cost -- by exchange, simulated annealing, evolutionary
#' mutation-selection, or guarded exact binary search. Unlike the greedy
#' constructors in [allocate_sparse_met()], the
#' allocation itself is optimised against the criterion on the coupled
#' [met_information()] matrix. The scoring can be made robust to uncertain
#' variance components by passing [robust_scenarios()].
#'
#' @details
#' The move set is chosen by `preserve`, which controls what is held fixed:
#' \describe{
#'   \item{`"margins"`}{Margin-preserving swaps only: both the replication of
#'     every genotype and the size of every environment stay fixed, so a fixed
#'     resource budget and (if started from an M4 design) equal replication are
#'     retained while quality is improved.}
#'   \item{`"replication"`}{Relocations: each genotype keeps its number of
#'     environments (equal replication preserved, M4-style) but environment
#'     sizes may change -- criterion-driven allocation with site-specific
#'     capacities.}
#'   \item{`"none"`}{Relocate/add/remove: plots may change (bounded by
#'     `objective$budget`), so cost is optimised jointly with quality and gain --
#'     the fully simultaneous mode.}
#' }
#' All evaluations are deterministic; components are normalised by the starting
#' design so the objective weights are on a comparable scale.
#'
#' @param allocation_matrix Starting genotype-by-environment 0/1 matrix
#'   (dimnames required), e.g. from [allocate_sparse_met()].
#' @param G,Sigma_E Relationship inputs for [met_information()].
#' @param objective A list of objective settings: `weights` (list with `gain`,
#'   `reliability`, `cost`), `prop`, `sigma_g`, `cost_per_plot`, `trait_weights`,
#'   `trait_gencov`, `budget`.
#' @param sigma_g2,sigma_e2 Nominal variance components (used when `robust` is
#'   `NULL`, and to normalise the objective reference).
#' @param tpe_weights Target-population weights passed to [met_information()].
#' @param env_efficiency Optional fixed site-efficiency vector. For genuine
#'   allocation-replication-layout optimisation, use `design_evaluator`.
#' @param design_evaluator Optional function taking one candidate allocation
#'   matrix and returning a list with any of `reps`, `env_efficiency`,
#'   `local_information`, `cost_per_plot`, `fixed_plot_overhead`, `feasible`,
#'   `fieldbooks`, and `metadata`. It is called once per unique candidate (results are cached),
#'   allowing the annealing search to score actual plantable local layouts.
#' @param preserve `"margins"`, `"replication"`, or `"none"` (see Details).
#' @param allocation_criterion `"weighted"`, `"mean_pev"`, `"cdmean"`, or
#'   `"expected_gain"`; passed to [design_objective()].
#' @param search_method Search engine. `"annealing"` is the scalable default;
#'   `"exchange"` accepts improving moves only; `"genetic"` uses an elitist
#'   mutation-selection population; and `"exact"` exhaustively evaluates all
#'   feasible binary designs for very small problems. `"mip"` is an alias for
#'   the dependency-free exact binary search and is intentionally guarded by
#'   `exact_max_cells`.
#' @param common_set Whether globally common treatments are `"fixed"`, forbidden
#'   (`"none"`), or allowed to change in size and identity
#'   (`"optimize_jointly"`).
#' @param common_treatments Optional fixed common-treatment IDs. With
#'   `common_set = "fixed"`, `NULL` infers them from all-environment rows of the
#'   starting design.
#' @param common_count_range Allowed minimum and maximum number of global common
#'   treatments during joint optimisation. Defaults to the feasible range.
#' @param common_weight Non-negative weight for a connectivity/diversity utility
#'   during joint common-set search. It is constant or zero outside
#'   `common_set = "optimize_jointly"`.
#' @param minimum_treatment_environments,maximum_treatment_environments Scalar or
#'   treatment-specific bounds on the number of environments assigned to each
#'   treatment.
#' @param exact_max_cells Maximum `nrow * ncol` accepted by exact/`"mip"`
#'   enumeration.
#' @param robust Optional list from [robust_scenarios()]; if given, designs are
#'   scored by [robust_design_score()].
#' @param robust_aggregate,cvar_alpha Aggregation for robust scoring.
#' @param seed_available Optional network-wide seed inventory, as a named
#'   numeric vector or a data frame with columns `Treatment` and
#'   `SeedAvailable`. Optimisation moves that would overspend a treatment are
#'   rejected.
#' @param seed_required_per_environment Seed consumed by one allocation in each
#'   environment: a positive scalar, a named numeric vector, or a data frame
#'   with columns `Environment` and `SeedRequiredPerPlot`.
#' @param minimum_seed_buffer Non-negative seed reserve retained for every
#'   treatment. Default 0.
#' @param environment_capacities Optional scalar or per-environment upper bounds
#'   on the number of allocated entries. This prevents unconstrained add/relocate
#'   moves from exceeding contracted site capacity.
#' @param minimum_environment_entries Scalar or per-environment lower bounds.
#'   Default 1 keeps every environment connected.
#' @param n_starts Number of random restarts. Default 5.
#' @param iters Annealing iterations per restart. Default 200.
#' @param cooling Geometric cooling factor per iteration. Default 0.97.
#' @param seed Optional RNG seed.
#' @param max_dim Guard passed to [met_information()].
#' @param verbose Print per-restart progress.
#' @return A list with the best `allocation_matrix`, `score`, raw objective
#'   `components`, starting score, restart `trace`, full `trajectory`, move and
#'   evaluator `diagnostics`, preservation and robustness settings, optional
#'   `seed_summary`, the best `design_evaluation`, and a validated
#'   `sparse_met_design` in `design`.
#' @seealso [design_objective()], [robust_design_score()], [pareto_designs()],
#'   [allocate_sparse_met()].
#' @examples
#' set.seed(1)
#' G <- crossprod(matrix(rnorm(24 * 60), 60, 24)) / 60 + diag(24) * 0.3
#' dimnames(G) <- list(paste0("L", 1:24), paste0("L", 1:24))
#' SigE <- diag(3) * 0.5 + 0.5
#' dimnames(SigE) <- list(paste0("E", 1:3), paste0("E", 1:3))
#' M <- matrix(0L, 24, 3, dimnames = list(rownames(G), colnames(SigE)))
#' for (e in 1:3) M[sample(24, 12), e] <- 1L
#' opt <- optimize_design(M, G, SigE, preserve = "margins",
#'                        n_starts = 2, iters = 40, seed = 1)
#' opt$score >= opt$score_start
#' @export
optimize_design <- function(allocation_matrix, G, Sigma_E = NULL,
                            objective = list(),
                            sigma_g2 = 1, sigma_e2 = 1,
                            tpe_weights = NULL, env_efficiency = NULL,
                            design_evaluator = NULL,
                            preserve = c("margins", "replication", "none"),
                            allocation_criterion = c("weighted", "mean_pev",
                                                     "cdmean", "expected_gain"),
                            search_method = c("annealing", "exchange", "genetic",
                                              "exact", "mip"),
                            common_set = c("fixed", "optimize_jointly", "none"),
                            common_treatments = NULL,
                            common_count_range = NULL,
                            common_weight = 0.10,
                            minimum_treatment_environments = 0L,
                            maximum_treatment_environments = NULL,
                            robust = NULL,
                            robust_aggregate = c("mean", "cvar", "min"),
                            cvar_alpha = 0.25,
                            seed_available = NULL,
                            seed_required_per_environment = NULL,
                            minimum_seed_buffer = 0,
                            environment_capacities = NULL,
                            minimum_environment_entries = 1L,
                            n_starts = 5L, iters = 200L, cooling = 0.97,
                            exact_max_cells = 16L,
                            seed = NULL, max_dim = 6000L, verbose = FALSE) {
  preserve <- match.arg(preserve)
  allocation_criterion <- match.arg(allocation_criterion)
  search_method <- match.arg(search_method)
  if (search_method == "mip") search_method <- "exact"
  common_set <- match.arg(common_set)
  robust_aggregate <- match.arg(robust_aggregate)
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
        abs(seed - round(seed)) > 1e-8)
      stop("`seed` must be one finite integer or NULL.")
    set.seed(as.integer(round(seed)))
  }
  if (length(minimum_seed_buffer) != 1L ||
      !is.numeric(minimum_seed_buffer) ||
      !is.finite(minimum_seed_buffer) || minimum_seed_buffer < 0)
    stop("`minimum_seed_buffer` must be one finite non-negative number.")
  if (is.null(seed_available) && minimum_seed_buffer > 0)
    stop("`seed_available` is required when `minimum_seed_buffer` is positive.")
  if (!is.null(design_evaluator) && !is.function(design_evaluator))
    stop("`design_evaluator` must be NULL or a function.")
  if (!is.numeric(common_weight) || length(common_weight) != 1L ||
      !is.finite(common_weight) || common_weight < 0)
    stop("`common_weight` must be one finite non-negative value.")
  if (!is.numeric(exact_max_cells) || length(exact_max_cells) != 1L ||
      !is.finite(exact_max_cells) || exact_max_cells < 1L ||
      exact_max_cells != as.integer(exact_max_cells))
    stop("`exact_max_cells` must be a positive integer.")
  exact_max_cells <- as.integer(exact_max_cells)

  M0 <- allocation_matrix
  storage.mode(M0) <- "integer"
  if (is.null(rownames(M0)) || is.null(colnames(M0)))
    stop("`allocation_matrix` needs row and column names.")
  if (anyDuplicated(rownames(M0)) || anyDuplicated(colnames(M0)) ||
      anyNA(M0) || !all(M0 %in% c(0L, 1L)))
    stop("`allocation_matrix` must be a named, unique, finite 0/1 matrix.")
  if (any(colSums(M0) == 0L))
    stop("Every environment in `allocation_matrix` must contain an entry.")

  normalise_treatment_bound <- function(x, arg, default) {
    if (is.null(x)) x <- default
    if (!is.numeric(x) || !length(x) || any(!is.finite(x)) ||
        any(abs(x - round(x)) > 1e-8))
      stop(sprintf("`%s` must contain finite integer values.", arg))
    if (!is.null(names(x)) && any(names(x) != "")) {
      if (anyDuplicated(names(x)) || !all(rownames(M0) %in% names(x)))
        stop(sprintf("Named `%s` must cover every treatment.", arg))
      x <- x[rownames(M0)]
    } else if (length(x) == 1L) x <- rep(x, nrow(M0))
    else if (length(x) != nrow(M0))
      stop(sprintf("`%s` must be scalar or have one value per treatment.", arg))
    stats::setNames(as.integer(round(x)), rownames(M0))
  }
  minimum_treatment_environments <- normalise_treatment_bound(
    minimum_treatment_environments, "minimum_treatment_environments", 0L)
  maximum_treatment_environments <- normalise_treatment_bound(
    maximum_treatment_environments, "maximum_treatment_environments", ncol(M0))
  if (any(minimum_treatment_environments < 0L) ||
      any(maximum_treatment_environments > ncol(M0)) ||
      any(minimum_treatment_environments > maximum_treatment_environments))
    stop("Treatment bounds must satisfy 0 <= minimum <= maximum <= ncol(allocation_matrix).")

  inferred_common <- rownames(M0)[rowSums(M0) == ncol(M0)]
  if (common_set == "fixed") {
    if (is.null(common_treatments)) common_treatments <- inferred_common
    common_treatments <- unique(as.character(common_treatments))
    if (!all(common_treatments %in% rownames(M0)))
      stop("`common_treatments` contains IDs absent from the allocation.")
    if (length(common_treatments) &&
        any(rowSums(M0[common_treatments, , drop = FALSE]) != ncol(M0)))
      stop("Every fixed common treatment must occur in every starting environment.")
  } else {
    common_treatments <- character(0)
  }
  if (is.null(common_count_range)) {
    common_count_range <- if (common_set == "fixed")
      c(length(common_treatments), min(colSums(M0))) else if (common_set == "none")
        c(0L, 0L) else c(0L, min(colSums(M0)))
  }
  if (!is.numeric(common_count_range) || length(common_count_range) != 2L ||
      any(!is.finite(common_count_range)) ||
      any(abs(common_count_range - round(common_count_range)) > 1e-8))
    stop("`common_count_range` must contain two finite integers.")
  common_count_range <- as.integer(round(common_count_range))
  if (common_count_range[1L] < 0L ||
      common_count_range[2L] < common_count_range[1L] ||
      common_count_range[2L] > min(colSums(M0)))
    stop("Invalid `common_count_range` for the starting environment sizes.")
  if (length(n_starts) != 1L || !is.finite(n_starts) ||
      n_starts < 1 || abs(n_starts - round(n_starts)) > 1e-8 ||
      length(iters) != 1L || !is.finite(iters) ||
      iters < 1 || abs(iters - round(iters)) > 1e-8)
    stop("`n_starts` and `iters` must be positive integers.")
  n_starts <- as.integer(round(n_starts))
  iters <- as.integer(round(iters))
  if (length(cooling) != 1L || !is.finite(cooling) ||
      cooling <= 0 || cooling > 1)
    stop("`cooling` must be in (0, 1].")

  normalise_environment_bound <- function(x, arg, default) {
    if (is.null(x)) x <- default
    if (!is.numeric(x) || !length(x) || any(!is.finite(x)) ||
        any(abs(x - round(x)) > 1e-8))
      stop(sprintf("`%s` must contain finite integer values.", arg))
    if (!is.null(names(x)) && any(names(x) != "")) {
      if (anyDuplicated(names(x)) || !all(colnames(M0) %in% names(x)))
        stop(sprintf("Named `%s` must cover every environment.", arg))
      x <- x[colnames(M0)]
    } else if (length(x) == 1L) {
      x <- rep(x, ncol(M0))
    } else if (length(x) != ncol(M0)) {
      stop(sprintf("`%s` must be scalar or have one value per environment.", arg))
    }
    stats::setNames(as.integer(round(x)), colnames(M0))
  }
  environment_capacities <- normalise_environment_bound(
    environment_capacities, "environment_capacities",
    rep(nrow(M0), ncol(M0)))
  minimum_environment_entries <- normalise_environment_bound(
    minimum_environment_entries, "minimum_environment_entries",
    rep(1L, ncol(M0)))
  if (any(minimum_environment_entries < 1L) ||
      any(environment_capacities < minimum_environment_entries) ||
      any(environment_capacities > nrow(M0)))
    stop("Environment bounds must satisfy 1 <= minimum <= capacity <= nrow(allocation_matrix).")

  seed_budget <- seed_cost <- NULL
  if (!is.null(seed_available)) {
    if (is.data.frame(seed_available)) {
      needed <- c("Treatment", "SeedAvailable")
      if (!all(needed %in% names(seed_available)))
        stop("Data-frame `seed_available` needs Treatment and SeedAvailable columns.")
      if (anyDuplicated(seed_available$Treatment))
        stop("`seed_available$Treatment` must be unique.")
      seed_vec <- stats::setNames(as.numeric(seed_available$SeedAvailable),
                                  as.character(seed_available$Treatment))
    } else {
      seed_vec <- as.numeric(seed_available)
      names(seed_vec) <- names(seed_available)
    }
    if (is.null(names(seed_vec)) || any(names(seed_vec) == "") ||
        anyDuplicated(names(seed_vec)) ||
        !all(rownames(M0) %in% names(seed_vec)) ||
        any(!is.finite(seed_vec)) || any(seed_vec < 0))
      stop("`seed_available` must be a uniquely named, finite non-negative inventory covering all treatments.")
    seed_budget <- seed_vec[rownames(M0)] - minimum_seed_buffer
    if (any(seed_budget < 0))
      stop("`minimum_seed_buffer` exceeds available seed for at least one treatment.")

    if (is.null(seed_required_per_environment))
      stop("`seed_required_per_environment` is required with `seed_available`.")
    if (is.data.frame(seed_required_per_environment)) {
      needed <- c("Environment", "SeedRequiredPerPlot")
      if (!all(needed %in% names(seed_required_per_environment)))
        stop("Seed-requirement data frame needs Environment and SeedRequiredPerPlot columns.")
      if (anyDuplicated(seed_required_per_environment$Environment))
        stop("Seed-requirement environments must be unique.")
      req <- stats::setNames(
        as.numeric(seed_required_per_environment$SeedRequiredPerPlot),
        as.character(seed_required_per_environment$Environment))
    } else {
      req <- as.numeric(seed_required_per_environment)
      names(req) <- names(seed_required_per_environment)
    }
    if (!is.null(names(req)) && any(names(req) != "")) {
      if (anyDuplicated(names(req)) ||
          !all(colnames(M0) %in% names(req)))
        stop("Named seed requirements must cover every environment.")
      req <- req[colnames(M0)]
    } else if (length(req) == 1L) {
      req <- rep(req, ncol(M0))
    } else if (length(req) != ncol(M0)) {
      stop("Seed requirements must be scalar or have one value per environment.")
    }
    if (any(!is.finite(req)) || any(req <= 0))
      stop("Seed requirements must be finite and strictly positive.")
    seed_cost <- stats::setNames(as.numeric(req), colnames(M0))
  } else if (!is.null(seed_required_per_environment)) {
    stop("`seed_available` is required when seed requirements are supplied.")
  }

  seed_used <- function(M)
    rowSums(sweep(M, 2L, seed_cost, `*`))
  design_feasible <- function(M) {
    loads <- colSums(M)
    treatment_loads <- rowSums(M)
    load_ok <- all(loads >= minimum_environment_entries &
                     loads <= environment_capacities)
    treatment_ok <- all(treatment_loads >= minimum_treatment_environments &
                          treatment_loads <= maximum_treatment_environments)
    common_now <- sum(treatment_loads == ncol(M))
    common_ok <- common_now >= common_count_range[1L] &&
      common_now <= common_count_range[2L]
    if (common_set == "fixed" && length(common_treatments))
      common_ok <- common_ok &&
        all(treatment_loads[common_treatments] == ncol(M))
    if (common_set == "none") common_ok <- common_now == 0L
    seed_ok <- is.null(seed_budget) ||
      all(seed_used(M) <= seed_budget + 1e-10)
    load_ok && treatment_ok && common_ok && seed_ok
  }
  if (!design_feasible(M0))
    stop("The starting allocation violates an environment or seed constraint.")

  obj <- utils::modifyList(
    list(weights = list(gain = 1, reliability = 0, cost = 0),
         prop = 0.1, sigma_g = 1, cost_per_plot = 1,
         trait_weights = NULL, trait_gencov = NULL,
         R_T = NULL, multitrait = "exact", budget = NULL),
    objective)

  evaluator_cache <- new.env(parent = emptyenv(), hash = TRUE)
  evaluator_calls <- 0L
  evaluator_cache_hits <- 0L
  evaluator_failures <- 0L
  evaluate_candidate <- function(M) {
    key <- paste0(as.integer(M), collapse = "")
    if (exists(key, envir = evaluator_cache, inherits = FALSE)) {
      evaluator_cache_hits <<- evaluator_cache_hits + 1L
      return(get(key, envir = evaluator_cache, inherits = FALSE))
    }
    evaluator_calls <<- evaluator_calls + 1L
    ans <- if (is.null(design_evaluator)) {
      list(reps = NULL, env_efficiency = env_efficiency,
           local_information = NULL, feasible = TRUE,
           cost_per_plot = obj$cost_per_plot, fixed_plot_overhead = 0)
    } else tryCatch(design_evaluator(M), error = function(e) {
      evaluator_failures <<- evaluator_failures + 1L
      structure(list(feasible = FALSE, error = conditionMessage(e)),
                class = "design_evaluator_error")
    })
    if (!is.list(ans))
      stop("`design_evaluator` must return a list.")
    if (is.null(ans$feasible)) ans$feasible <- TRUE
    if (!is.logical(ans$feasible) || length(ans$feasible) != 1L ||
        is.na(ans$feasible))
      stop("`design_evaluator()$feasible` must be TRUE or FALSE.")
    if (is.null(ans$cost_per_plot)) ans$cost_per_plot <- obj$cost_per_plot
    if (is.null(ans$fixed_plot_overhead)) ans$fixed_plot_overhead <- 0
    assign(key, ans, envir = evaluator_cache)
    ans
  }

  # Reference for normalisation: components of the starting design at nominal
  # variance parameters.
  base_eval <- evaluate_candidate(M0)
  if (!isTRUE(base_eval$feasible))
    stop("The starting allocation cannot be converted to a feasible local design",
         if (!is.null(base_eval$error)) paste0(": ", base_eval$error) else ".")
  base <- design_objective(M0, G = G, Sigma_E = Sigma_E,
                           sigma_g2 = sigma_g2, sigma_e2 = sigma_e2,
                           reps = base_eval$reps,
                           env_efficiency = base_eval$env_efficiency,
                           tpe_weights = tpe_weights,
                           local_information = base_eval$local_information,
                           prop = obj$prop, sigma_g = obj$sigma_g,
                           trait_weights = obj$trait_weights,
                           trait_gencov = obj$trait_gencov,
                           R_T = obj$R_T, multitrait = obj$multitrait,
                           cost_per_plot = base_eval$cost_per_plot,
                           fixed_plot_overhead = base_eval$fixed_plot_overhead,
                           criterion = allocation_criterion,
                           weights = obj$weights, ref = NULL,
                           budget = NULL, max_dim = max_dim)
  ref <- list(gain = base$gain, reliability = base$reliability,
              mean_PEV = base$mean_PEV,
              cost = if (base$cost > 0) base$cost else 1)

  score_fun <- function(M) {
    if (!design_feasible(M)) return(-Inf)
    ev <- evaluate_candidate(M)
    if (!isTRUE(ev$feasible)) return(-Inf)
    ans <- if (is.null(robust)) {
      design_objective(M, G = G, Sigma_E = Sigma_E,
                       sigma_g2 = sigma_g2, sigma_e2 = sigma_e2,
                       reps = ev$reps, env_efficiency = ev$env_efficiency,
                       tpe_weights = tpe_weights,
                       local_information = ev$local_information,
                       prop = obj$prop, sigma_g = obj$sigma_g,
                       trait_weights = obj$trait_weights,
                       trait_gencov = obj$trait_gencov,
                       R_T = obj$R_T, multitrait = obj$multitrait,
                       cost_per_plot = ev$cost_per_plot,
                       fixed_plot_overhead = ev$fixed_plot_overhead,
                       criterion = allocation_criterion,
                       weights = obj$weights, ref = ref,
                       budget = obj$budget, max_dim = max_dim)$score
    } else {
      robust_design_score(M, G = G, Sigma_E = Sigma_E, scenarios = robust,
                          aggregate = robust_aggregate, cvar_alpha = cvar_alpha,
                          prop = obj$prop, sigma_g = obj$sigma_g,
                          trait_weights = obj$trait_weights,
                          trait_gencov = obj$trait_gencov,
                          R_T = obj$R_T, multitrait = obj$multitrait,
                          reps = ev$reps,
                          env_efficiency = ev$env_efficiency,
                          tpe_weights = tpe_weights,
                          local_information = ev$local_information,
                          cost_per_plot = ev$cost_per_plot,
                          fixed_plot_overhead = ev$fixed_plot_overhead,
                          criterion = allocation_criterion,
                          weights = obj$weights, ref = ref,
                          budget = obj$budget, max_dim = max_dim)$score
    }
    if (common_set == "optimize_jointly" && common_weight > 0)
      ans <- ans + common_weight * .joint_common_utility(M, G, Sigma_E)
    ans
  }

  score_start <- score_fun(M0)
  if (!is.finite(score_start))
    stop("The starting design is infeasible under the objective budget or ",
         "produces a non-finite score.")
  best_M <- M0; best_score <- score_start; trace <- numeric(0)
  proposals <- accepted <- improved <- infeasible_proposals <- 0L
  trajectory <- list()

  if (search_method == "exact") {
    nc <- length(M0)
    if (nc > exact_max_cells)
      stop("Exact/`mip` search requires nrow * ncol <= `exact_max_cells` (",
           exact_max_cells, "); use annealing, exchange, or genetic search.")
    row_target <- rowSums(M0); col_target <- colSums(M0)
    n_candidates <- 2^nc
    for (code in 0:(n_candidates - 1)) {
      proposals <- proposals + 1L
      bits <- as.integer(intToBits(code))[seq_len(nc)]
      candidate <- matrix(bits, nrow(M0), ncol(M0),
                          dimnames = dimnames(M0))
      if (preserve == "margins" &&
          (!identical(as.integer(rowSums(candidate)), as.integer(row_target)) ||
           !identical(as.integer(colSums(candidate)), as.integer(col_target)))) next
      if (preserve == "replication" &&
          !identical(as.integer(rowSums(candidate)), as.integer(row_target))) next
      if (!design_feasible(candidate)) next
      new_s <- score_fun(candidate)
      accepted <- accepted + 1L
      if (is.finite(new_s) && new_s > best_score) {
        improved <- improved + 1L
        best_score <- new_s; best_M <- candidate
      }
    }
    trace <- best_score
    trajectory[[1L]] <- data.frame(restart = 1L, final_score = best_score,
                                    best_score = best_score)
  } else if (search_method == "genetic") {
    population_size <- max(4L, n_starts)
    population <- vector("list", population_size)
    population[[1L]] <- M0
    for (i in 2:population_size) {
      candidate <- .perturb_design(M0, preserve, obj$budget, 5L + i,
                                   common_set)
      population[[i]] <- if (design_feasible(candidate)) candidate else M0
    }
    scores <- vapply(population, score_fun, numeric(1))
    for (generation in seq_len(iters)) {
      elite_n <- max(2L, ceiling(population_size / 2))
      elite <- order(scores, decreasing = TRUE)[seq_len(elite_n)]
      next_population <- population[elite]
      while (length(next_population) < population_size) {
        parent <- population[[sample(elite, 1L)]]
        child <- .perturb_design(parent, preserve, obj$budget,
                                 sample.int(4L, 1L), common_set)
        proposals <- proposals + 1L
        if (!design_feasible(child)) {
          infeasible_proposals <- infeasible_proposals + 1L
          child <- parent
        } else accepted <- accepted + 1L
        next_population[[length(next_population) + 1L]] <- child
      }
      population <- next_population
      scores <- vapply(population, score_fun, numeric(1))
      gi <- which.max(scores)
      if (scores[gi] > best_score) {
        improved <- improved + 1L
        best_score <- scores[gi]; best_M <- population[[gi]]
      }
      trace <- c(trace, best_score)
    }
    trajectory <- lapply(seq_along(population), function(i)
      data.frame(restart = i, final_score = scores[i], best_score = best_score))
  } else {
    trajectory <- vector("list", n_starts)
    for (st in seq_len(n_starts)) {
      # First restart begins at the supplied design; later ones are perturbed.
      cur_M <- if (st == 1L) M0 else
        .perturb_design(M0, preserve, obj$budget, 10L, common_set)
      if (!design_feasible(cur_M)) cur_M <- M0
      cur_s <- score_fun(cur_M)
      restart_best <- cur_s
      temp <- max(abs(cur_s), 1e-3) * 0.5

      for (it in seq_len(iters)) {
        proposals <- proposals + 1L
        prop_M <- .propose_move(cur_M, preserve, obj$budget, common_set)
        if (is.null(prop_M)) { temp <- temp * cooling; next }
        if (!design_feasible(prop_M)) {
          infeasible_proposals <- infeasible_proposals + 1L
          temp <- temp * cooling; next
        }
        new_s <- score_fun(prop_M)
        d <- new_s - cur_s
        accept <- is.finite(new_s) && (d > 0 ||
          (search_method == "annealing" && temp > 0 &&
             stats::runif(1) < exp(d / temp)))
        if (accept) {
          accepted <- accepted + 1L
          cur_M <- prop_M; cur_s <- new_s
          if (cur_s > restart_best) restart_best <- cur_s
          if (cur_s > best_score) {
            improved <- improved + 1L
            best_score <- cur_s; best_M <- cur_M
          }
        }
        temp <- temp * cooling
      }
      trace <- c(trace, restart_best)
      trajectory[[st]] <- data.frame(restart = st, final_score = cur_s,
                                      best_score = restart_best)
      if (verbose) message(sprintf("restart %d: score = %.5f (best = %.5f)",
                                   st, cur_s, best_score))
    }
  }

  best_eval <- evaluate_candidate(best_M)
  comp <- design_objective(best_M, G = G, Sigma_E = Sigma_E,
                           sigma_g2 = sigma_g2, sigma_e2 = sigma_e2,
                           reps = best_eval$reps,
                           env_efficiency = best_eval$env_efficiency,
                           tpe_weights = tpe_weights,
                           local_information = best_eval$local_information,
                           prop = obj$prop, sigma_g = obj$sigma_g,
                           trait_weights = obj$trait_weights,
                           trait_gencov = obj$trait_gencov,
                           R_T = obj$R_T, multitrait = obj$multitrait,
                           cost_per_plot = best_eval$cost_per_plot,
                           fixed_plot_overhead = best_eval$fixed_plot_overhead,
                           criterion = allocation_criterion,
                           weights = obj$weights, ref = ref,
                           budget = obj$budget, max_dim = max_dim)

  seed_summary <- if (is.null(seed_budget)) NULL else {
    used <- seed_used(best_M)
    data.frame(Treatment = rownames(best_M),
               SeedAvailable = as.numeric(seed_budget + minimum_seed_buffer),
               MinimumBuffer = rep(minimum_seed_buffer, nrow(best_M)),
               SeedAllocated = as.numeric(used),
               SeedRemaining = as.numeric(seed_budget +
                                            minimum_seed_buffer - used),
               SeedSpendableRemaining = as.numeric(seed_budget - used),
               stringsAsFactors = FALSE)
  }

  optimizer_diagnostics <- list(
    search_method = search_method,
    exact_certified = identical(search_method, "exact"),
    proposals = proposals, accepted = accepted,
    acceptance_rate = accepted / max(1L, proposals),
    improvements = improved,
    infeasible_proposals = infeasible_proposals,
    unique_design_evaluations = evaluator_calls,
    evaluator_cache_hits = evaluator_cache_hits,
    evaluator_failures = evaluator_failures,
    common_set = common_set,
    common_treatments = rownames(best_M)[rowSums(best_M) == ncol(best_M)],
    n_common = sum(rowSums(best_M) == ncol(best_M)))
  design <- sparse_met_design(
    best_M,
    reps = if (is.null(best_eval$reps)) best_M else best_eval$reps,
    fieldbooks = best_eval$fieldbooks,
    G = G, Sigma_E = Sigma_E, tpe_weights = tpe_weights,
    sigma_g2 = sigma_g2, sigma_e2 = sigma_e2,
    diagnostics = optimizer_diagnostics,
    provenance = list(engine = "optimize_design",
                      preserve = preserve,
                      allocation_criterion = allocation_criterion,
                      search_method = search_method,
                      common_set = common_set,
                      robust = !is.null(robust),
                      seed = seed))

  list(allocation_matrix = best_M, score = best_score,
       components = comp[c("reliability", "mean_PEV", "gain", "plots", "cost")],
       score_start = score_start, trace = trace,
       preserve = preserve, allocation_criterion = allocation_criterion,
       search_method = search_method, common_set = common_set,
       common_treatments = optimizer_diagnostics$common_treatments,
       robust = !is.null(robust),
       seed_summary = seed_summary,
       design_evaluation = best_eval,
       design = design,
       trajectory = do.call(rbind, trajectory),
       diagnostics = optimizer_diagnostics)
}


# ---- move generators --------------------------------------------------------

# One random move respecting `preserve`; returns a modified matrix or NULL.
.propose_move <- function(M, preserve, budget, common_set = "fixed") {
  if (common_set == "optimize_jointly" && stats::runif(1) < 0.30) {
    common_move <- sample(c("identity", "promote", "demote"), 1L)
    ans <- switch(common_move,
                  identity = .move_common_identity(M),
                  promote = .move_common_promote(M),
                  demote = .move_common_demote(M))
    if (!is.null(ans)) return(ans)
  }
  if (preserve == "margins") return(.move_swap(M))
  if (preserve == "replication") return(.move_relocate(M))
  # preserve == "none": mix relocate / add / remove
  u <- stats::runif(1)
  if (u < 0.5) .move_relocate(M)
  else if (u < 0.75) .move_add(M, budget)
  else .move_remove(M)
}

# Margin-preserving swap: a in e1 (not e2) and b in e2 (not e1) exchange envs.
.move_swap <- function(M) {
  E <- ncol(M); if (E < 2L) return(NULL)
  es <- sample.int(E, 2L); e1 <- es[1L]; e2 <- es[2L]
  A <- which(M[, e1] == 1L & M[, e2] == 0L)
  B <- which(M[, e2] == 1L & M[, e1] == 0L)
  if (!length(A) || !length(B)) return(NULL)
  a <- A[sample.int(length(A), 1L)]; b <- B[sample.int(length(B), 1L)]
  M[a, e1] <- 0L; M[a, e2] <- 1L; M[b, e2] <- 0L; M[b, e1] <- 1L
  M
}

# Relocate: move one genotype from an environment it is in to one it is not
# (its replication is unchanged; environment sizes change).
.move_relocate <- function(M) {
  E <- ncol(M); if (E < 2L) return(NULL)
  g <- sample.int(nrow(M), 1L)
  inn <- which(M[g, ] == 1L); out <- which(M[g, ] == 0L)
  if (!length(inn) || !length(out)) return(NULL)
  e1 <- inn[sample.int(length(inn), 1L)]; e2 <- out[sample.int(length(out), 1L)]
  if (sum(M[, e1]) <= 1L) return(NULL)          # keep environment non-empty
  M[g, e1] <- 0L; M[g, e2] <- 1L
  M
}

# Add one plot (bounded by budget).
.move_add <- function(M, budget) {
  if (!is.null(budget) && sum(M != 0) >= budget) return(NULL)
  zeros <- which(M == 0L)
  if (!length(zeros)) return(NULL)
  M[zeros[sample.int(length(zeros), 1L)]] <- 1L
  M
}

# Remove one plot, keeping every environment non-empty.
.move_remove <- function(M) {
  ones <- which(M == 1L)
  if (length(ones) <= ncol(M)) return(NULL)
  idx <- ones[sample.int(length(ones), 1L)]
  e <- ((idx - 1L) %/% nrow(M)) + 1L
  if (sum(M[, e]) <= 1L) return(NULL)
  M[idx] <- 0L
  M
}

# Apply several random moves to diversify a restart.
.perturb_design <- function(M, preserve, budget, n, common_set = "fixed") {
  for (i in seq_len(n)) {
    m <- .propose_move(M, preserve, budget, common_set)
    if (!is.null(m)) M <- m
  }
  M
}


# Swap the identity of one global common treatment while preserving every
# environment margin exactly.
.move_common_identity <- function(M) {
  common <- which(rowSums(M) == ncol(M))
  other <- which(rowSums(M) > 0L & rowSums(M) < ncol(M))
  if (!length(common) || !length(other)) return(NULL)
  a <- sample(common, 1L); b <- sample(other, 1L)
  tmp <- M[a, ]; M[a, ] <- M[b, ]; M[b, ] <- tmp
  M
}


# Promote a partially observed treatment to the global common set. Donor cells
# are taken within the missing environments, preserving all environment sizes.
.move_common_promote <- function(M) {
  rs <- rowSums(M)
  candidates <- which(rs > 0L & rs < ncol(M))
  if (!length(candidates)) return(NULL)
  b <- sample(candidates, 1L)
  missing <- which(M[b, ] == 0L)
  out <- M
  for (e in missing) {
    donors <- which(out[, e] == 1L & rowSums(out) > 1L &
                      rowSums(out) < ncol(out))
    donors <- setdiff(donors, b)
    if (!length(donors)) return(NULL)
    d <- sample(donors, 1L)
    out[d, e] <- 0L; out[b, e] <- 1L
  }
  out
}


# Demote one common treatment in one environment and transfer that cell to a
# non-common recipient, again preserving the environment margin.
.move_common_demote <- function(M) {
  common <- which(rowSums(M) == ncol(M))
  if (!length(common)) return(NULL)
  a <- sample(common, 1L)
  env_order <- sample.int(ncol(M))
  for (e in env_order) {
    recipients <- which(M[, e] == 0L & rowSums(M) < ncol(M) - 1L)
    if (!length(recipients)) next
    b <- sample(recipients, 1L)
    M[a, e] <- 0L; M[b, e] <- 1L
    return(M)
  }
  NULL
}


# Unit-scale joint common-set utility. Pairwise overlap protects covariance
# estimation; genetic effective size discourages redundant common anchors.
.joint_common_utility <- function(M, G, Sigma_E, target_se = 0.15) {
  E <- ncol(M); J <- nrow(M)
  if (E < 2L || J < 1L) return(0)
  overlap <- crossprod(M)
  achieved <- overlap[upper.tri(overlap)]
  if (is.null(Sigma_E)) {
    targets <- rep(max(2, ceiling(sqrt(J))), length(achieved))
  } else {
    S <- as.matrix(Sigma_E)[colnames(M), colnames(M), drop = FALSE]
    rho <- stats::cov2cor(S)[upper.tri(S)]
    targets <- pmin(J, pmax(2, ceiling(((1 - rho^2) / target_se)^2)))
  }
  attainment <- pmin(1, achieved / pmax(1, targets))
  common <- which(rowSums(M) == E)
  diversity <- 0
  if (length(common)) {
    Gc <- as.matrix(G)[rownames(M)[common], rownames(M)[common], drop = FALSE]
    neff <- sum(diag(Gc))^2 / max(sum(Gc^2), .Machine$double.eps)
    diversity <- min(1, neff / length(common))
  }
  0.50 * min(attainment) + 0.25 * mean(attainment) + 0.25 * diversity
}
