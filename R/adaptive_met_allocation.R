#' Sequentially add the most informative MET observations
#'
#' @description
#' Builds or updates a sparse multi-environment allocation one
#' genotype-by-environment cell at a time. At each step it evaluates every
#' feasible unobserved cell and adds the one with the largest marginal gain in
#' the requested prediction criterion. This supports rolling multi-year plans
#' and within-season augmentation after observations, seed, or site capacity
#' change.
#'
#' @param observed_allocation Named genotype-by-environment 0/1 matrix describing
#'   observations already collected or irrevocably committed.
#' @param G,Sigma_E Relationship matrices used by [design_objective()].
#' @param target_environment_sizes Scalar or named vector giving the maximum
#'   final number of allocated treatments in each environment.
#' @param batch_size Number of new cells to recommend. `NULL` fills every
#'   environment to its target size.
#' @param allocation_criterion `"mean_pev"`, `"cdmean"`, or
#'   `"expected_gain"`.
#' @param robust Optional scenario list from [robust_scenarios()].
#' @param robust_aggregate,cvar_alpha Robust aggregation controls.
#' @param candidate_mask Optional named 0/1 matrix identifying cells that may be
#'   added. Already observed cells need not be marked.
#' @param required_common_treatments Treatments that must be present in every
#'   environment before discretionary additions are selected.
#' @param minimum_treatment_environments Minimum final environment coverage per
#'   treatment. Coverage deficits are filled before discretionary additions.
#' @param sigma_g2,sigma_e2,tpe_weights,env_efficiency Statistical assumptions.
#' @param cost_per_plot Optional scalar or environment-specific plot cost.
#' @param seed_available Optional named treatment-level seed inventory.
#' @param seed_required_per_environment Optional scalar or environment-specific
#'   seed cost for one new observation; required with `seed_available`.
#' @param minimum_seed_buffer Seed reserve that cannot be spent.
#' @param seed,max_dim Reproducibility and information-matrix size control.
#'
#' @return A list with the updated `allocation_matrix`, ordered
#'   `recommendations`, criterion scores before and after, an audit `trajectory`,
#'   seed use, and a validated `sparse_met_design` object.
#' @seealso [optimize_design()], [allocate_sparse_met()], [met_information()]
#' @export
adaptive_met_allocation <- function(
    observed_allocation,
    G,
    Sigma_E = NULL,
    target_environment_sizes,
    batch_size = NULL,
    allocation_criterion = c("mean_pev", "cdmean", "expected_gain"),
    robust = NULL,
    robust_aggregate = c("mean", "cvar", "min"),
    cvar_alpha = 0.25,
    candidate_mask = NULL,
    required_common_treatments = NULL,
    minimum_treatment_environments = 0L,
    sigma_g2 = 1,
    sigma_e2 = 1,
    tpe_weights = NULL,
    env_efficiency = NULL,
    cost_per_plot = 1,
    seed_available = NULL,
    seed_required_per_environment = NULL,
    minimum_seed_buffer = 0,
    seed = NULL,
    max_dim = 6000L) {
  allocation_criterion <- match.arg(allocation_criterion)
  robust_aggregate <- match.arg(robust_aggregate)
  if (!is.null(seed)) set.seed(seed)

  M <- as.matrix(observed_allocation)
  if (!is.numeric(M) || anyNA(M) || !all(M %in% c(0, 1)) ||
      is.null(rownames(M)) || is.null(colnames(M)) ||
      anyDuplicated(rownames(M)) || anyDuplicated(colnames(M)))
    stop("`observed_allocation` must be a uniquely named finite 0/1 matrix.")
  storage.mode(M) <- "integer"
  J <- nrow(M); E <- ncol(M); trts <- rownames(M); envs <- colnames(M)
  if (J < 1L || E < 2L) stop("At least one treatment and two environments are required.")

  normalise_env <- function(x, arg, integer = FALSE, positive = FALSE) {
    if (!is.numeric(x) || !length(x) || any(!is.finite(x)))
      stop("`", arg, "` must contain finite numeric values.")
    if (!is.null(names(x)) && any(names(x) != "")) {
      if (anyDuplicated(names(x)) || !all(envs %in% names(x)))
        stop("Named `", arg, "` must cover every environment.")
      x <- x[envs]
    } else if (length(x) == 1L) x <- rep(x, E)
    else if (length(x) != E)
      stop("`", arg, "` must be scalar or have one value per environment.")
    if (integer && any(abs(x - round(x)) > 1e-8))
      stop("`", arg, "` must contain integers.")
    if (positive && any(x <= 0)) stop("`", arg, "` must be strictly positive.")
    stats::setNames(if (integer) as.integer(round(x)) else as.numeric(x), envs)
  }
  capacities <- normalise_env(target_environment_sizes,
                              "target_environment_sizes", integer = TRUE)
  if (any(capacities < colSums(M)) || any(capacities > J))
    stop("Target environment sizes must lie between current sizes and nrow(allocation).")
  cost_per_plot <- normalise_env(cost_per_plot, "cost_per_plot")
  if (any(cost_per_plot < 0)) stop("`cost_per_plot` cannot be negative.")

  if (is.null(candidate_mask)) candidate_mask <- matrix(
    1L, J, E, dimnames = dimnames(M))
  candidate_mask <- as.matrix(candidate_mask)
  if (!identical(dim(candidate_mask), dim(M)) ||
      is.null(rownames(candidate_mask)) || is.null(colnames(candidate_mask)) ||
      !setequal(rownames(candidate_mask), trts) ||
      !setequal(colnames(candidate_mask), envs))
    stop("`candidate_mask` must be a named matrix aligned with the allocation.")
  candidate_mask <- candidate_mask[trts, envs, drop = FALSE]
  if (anyNA(candidate_mask) || !all(candidate_mask %in% c(0, 1)))
    stop("`candidate_mask` must contain only 0/1 values.")

  min_rep <- minimum_treatment_environments
  if (!is.numeric(min_rep) || !length(min_rep) || any(!is.finite(min_rep)))
    stop("`minimum_treatment_environments` must be numeric.")
  if (!is.null(names(min_rep)) && any(names(min_rep) != "")) {
    if (anyDuplicated(names(min_rep)) || !all(trts %in% names(min_rep)))
      stop("Named minimum treatment coverage must cover every treatment.")
    min_rep <- min_rep[trts]
  } else if (length(min_rep) == 1L) min_rep <- rep(min_rep, J)
  else if (length(min_rep) != J)
    stop("Minimum treatment coverage must be scalar or treatment-specific.")
  if (any(abs(min_rep - round(min_rep)) > 1e-8) ||
      any(min_rep < 0 | min_rep > E))
    stop("Minimum treatment coverage must contain integers in [0, n_environment].")
  min_rep <- stats::setNames(as.integer(round(min_rep)), trts)

  required_common_treatments <- unique(as.character(required_common_treatments))
  if (length(required_common_treatments) &&
      !all(required_common_treatments %in% trts))
    stop("`required_common_treatments` contains unknown IDs.")

  seed_cost <- seed_budget <- NULL
  if (!is.null(seed_available)) {
    if (is.data.frame(seed_available)) {
      if (!all(c("Treatment", "SeedAvailable") %in% names(seed_available)) ||
          anyDuplicated(seed_available$Treatment))
        stop("Seed data must contain unique Treatment and SeedAvailable columns.")
      seed_available <- stats::setNames(
        as.numeric(seed_available$SeedAvailable),
        as.character(seed_available$Treatment))
    }
    if (!is.numeric(seed_available) || is.null(names(seed_available)) ||
        anyDuplicated(names(seed_available)) ||
        !all(trts %in% names(seed_available)) ||
        any(!is.finite(seed_available)) || any(seed_available < 0))
      stop("`seed_available` must be a named finite non-negative vector covering treatments.")
    if (is.null(seed_required_per_environment))
      stop("`seed_required_per_environment` is required with seed inventory.")
    if (is.data.frame(seed_required_per_environment)) {
      if (!all(c("Environment", "SeedRequiredPerPlot") %in%
               names(seed_required_per_environment)) ||
          anyDuplicated(seed_required_per_environment$Environment))
        stop("Seed requirements must contain unique Environment and SeedRequiredPerPlot columns.")
      seed_required_per_environment <- stats::setNames(
        as.numeric(seed_required_per_environment$SeedRequiredPerPlot),
        as.character(seed_required_per_environment$Environment))
    }
    seed_cost <- normalise_env(seed_required_per_environment,
                               "seed_required_per_environment", positive = TRUE)
    seed_budget <- seed_available[trts] - minimum_seed_buffer
    if (any(seed_budget < 0)) stop("The minimum seed buffer exceeds available seed.")
  }
  seed_used <- function(X) if (is.null(seed_cost))
    stats::setNames(rep(0, J), trts) else
      stats::setNames(as.numeric(X %*% seed_cost), trts)

  total_open <- sum(capacities - colSums(M))
  if (is.null(batch_size)) batch_size <- total_open
  if (!is.numeric(batch_size) || length(batch_size) != 1L ||
      !is.finite(batch_size) || batch_size < 0 ||
      batch_size != as.integer(batch_size))
    stop("`batch_size` must be NULL or one non-negative integer.")
  batch_size <- min(as.integer(batch_size), total_open)

  base <- design_objective(
    M, G = G, Sigma_E = Sigma_E, sigma_g2 = sigma_g2, sigma_e2 = sigma_e2,
    tpe_weights = tpe_weights, env_efficiency = env_efficiency,
    cost_per_plot = cost_per_plot, criterion = allocation_criterion)
  ref <- list(gain = base$gain, reliability = base$reliability,
              mean_PEV = base$mean_PEV,
              cost = if (base$cost > 0) base$cost else 1)
  score_design <- function(X) {
    if (is.null(robust))
      design_objective(
        X, G = G, Sigma_E = Sigma_E, sigma_g2 = sigma_g2,
        sigma_e2 = sigma_e2, tpe_weights = tpe_weights,
        env_efficiency = env_efficiency, cost_per_plot = cost_per_plot,
        criterion = allocation_criterion, ref = ref)$score
    else robust_design_score(
      X, G = G, Sigma_E = Sigma_E, scenarios = robust,
      aggregate = robust_aggregate, cvar_alpha = cvar_alpha,
      tpe_weights = tpe_weights, env_efficiency = env_efficiency,
      cost_per_plot = cost_per_plot, criterion = allocation_criterion,
      ref = ref)$score
  }
  score_before <- score_design(M)
  trajectory <- list()
  recommendations <- list()
  additions <- 0L

  feasible_cells <- function(X, restrict_rows = NULL) {
    idx <- which(X == 0L & candidate_mask == 1L, arr.ind = TRUE)
    if (!nrow(idx)) return(idx)
    idx <- idx[colSums(X)[idx[, 2L]] < capacities[idx[, 2L]], , drop = FALSE]
    if (!is.null(restrict_rows) && nrow(idx))
      idx <- idx[idx[, 1L] %in% restrict_rows, , drop = FALSE]
    if (!is.null(seed_budget) && nrow(idx)) {
      used <- seed_used(X)
      ok <- used[idx[, 1L]] + seed_cost[idx[, 2L]] <=
        seed_budget[idx[, 1L]] + 1e-10
      idx <- idx[ok, , drop = FALSE]
    }
    idx
  }
  add_best <- function(reason, restrict_rows = NULL) {
    idx <- feasible_cells(M, restrict_rows)
    if (!nrow(idx)) return(FALSE)
    current <- score_design(M)
    scores <- numeric(nrow(idx))
    for (ii in seq_len(nrow(idx))) {
      X <- M; X[idx[ii, 1L], idx[ii, 2L]] <- 1L
      scores[ii] <- score_design(X)
    }
    labels <- paste(rownames(M)[idx[, 1L]], colnames(M)[idx[, 2L]], sep = "\r")
    best <- order(-scores, labels)[1L]
    r <- idx[best, 1L]; e <- idx[best, 2L]
    M[r, e] <<- 1L
    additions <<- additions + 1L
    recommendations[[additions]] <<- data.frame(
      Step = additions, Treatment = rownames(M)[r], Environment = colnames(M)[e],
      Reason = reason, ScoreBefore = current, ScoreAfter = scores[best],
      MarginalGain = scores[best] - current, Cost = cost_per_plot[e],
      stringsAsFactors = FALSE)
    trajectory[[additions]] <<- data.frame(
      Step = additions, Score = scores[best], RemainingCapacity =
        sum(capacities - colSums(M)), stringsAsFactors = FALSE)
    TRUE
  }

  mandatory_missing <- if (length(required_common_treatments))
    sum(M[required_common_treatments, , drop = FALSE] == 0L) else 0L
  M_mandatory <- M
  if (length(required_common_treatments))
    M_mandatory[required_common_treatments, ] <- 1L
  coverage_missing <- sum(pmax(0L, min_rep - rowSums(M_mandatory)))
  if (mandatory_missing + coverage_missing > batch_size)
    stop("`batch_size` is too small to satisfy required common and coverage constraints.")

  if (length(required_common_treatments)) {
    while (any(M[required_common_treatments, , drop = FALSE] == 0L)) {
      deficient <- which(rowSums(M) < E & trts %in% required_common_treatments)
      if (!add_best("required_common", deficient))
        stop("Required common treatments cannot be completed under capacity/seed constraints.")
    }
  }
  while (any(rowSums(M) < min_rep)) {
    deficient <- which(rowSums(M) < min_rep)
    if (!add_best("minimum_coverage", deficient))
      stop("Minimum treatment coverage is infeasible under capacity/seed constraints.")
  }
  while (additions < batch_size && sum(colSums(M) < capacities) > 0L) {
    if (!add_best("marginal_information")) break
  }

  rec <- if (length(recommendations)) do.call(rbind, recommendations) else
    data.frame(Step = integer(), Treatment = character(), Environment = character(),
               Reason = character(), ScoreBefore = numeric(), ScoreAfter = numeric(),
               MarginalGain = numeric(), Cost = numeric())
  traj <- if (length(trajectory)) do.call(rbind, trajectory) else
    data.frame(Step = integer(), Score = numeric(), RemainingCapacity = integer())
  score_after <- score_design(M)
  diagnostics <- list(
    engine = "adaptive_sequential", additions = additions,
    requested_batch_size = batch_size, allocation_criterion = allocation_criterion,
    robust = !is.null(robust), capacities = capacities)
  design <- sparse_met_design(
    M, reps = M, G = G, Sigma_E = Sigma_E, tpe_weights = tpe_weights,
    sigma_g2 = sigma_g2, sigma_e2 = sigma_e2, diagnostics = diagnostics,
    provenance = list(engine = "adaptive_met_allocation", seed = seed))
  list(allocation_matrix = M, recommendations = rec,
       score_before = score_before, score_after = score_after,
       trajectory = traj, seed_used = seed_used(M), diagnostics = diagnostics,
       design = design)
}
