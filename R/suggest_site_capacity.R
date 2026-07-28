#' Suggest the optimal number of plots for an unconstrained site
#'
#' @description
#' When a location is offered with no hard limit on field capacity ("use as much
#' space as you like"), more plots test more of the target population and raise
#' prediction accuracy and genetic gain -- but with diminishing returns, and
#' every extra plot costs resources. `suggest_site_capacity()` sweeps candidate
#' plot counts, simulates realised accuracy and gain, and recommends the plot
#' number at the point of diminishing returns, so the breeder maximises gain
#' while not wasting resources. Three scopes are supported:
#' \describe{
#'   \item{`scope = "single"`}{One generous site varies while the other sites
#'     are held at their real capacities. Returns the plot number for that site.}
#'   \item{`scope = "all"`}{*Every* site is unconstrained: the candidate capacity
#'     is applied to all `n_env` sites at once, and the recommendation is the
#'     common per-site capacity to use everywhere.}
#'   \item{`scope = "subset"`}{Only `focal_envs` vary; every other environment
#'     is held at its named value in `site_capacities`. This is the appropriate
#'     mode when several sites are unconstrained but partner sites have fixed
#'     plot commitments.}
#' }
#'
#' @details
#' For each candidate capacity `n`, a nested sparse design is built: the focal
#' environment(s) test `n` entries while all non-focal capacities stay fixed.
#' [simulate_met()] then returns realised selection accuracy and genetic gain
#' for the whole design under a common Monte Carlo stream. The marginal gain in
#' accuracy per additional *trial-wide physical plot* is tracked; the
#' recommendation is the largest consecutive candidate whose marginal gain
#' still meets `min_gain`. Beyond that point extra plots buy little accuracy
#' and may be better spent elsewhere. When `cost_per_plot` is supplied,
#' gain-per-cost is also reported.
#'
#' When `robust` is supplied, accuracy and gain are aggregated across its
#' variance/covariance scenarios before capacity selection. Maximin aggregation
#' does not use scenario probabilities and is therefore suitable for separate
#' main-effect and interaction structures whose relative weights are unknown.
#'
#' Each plot holds one entry unless `plots_per_entry > 1` (e.g. replicated
#' checks), which scales the plot and cost accounting but not the number of
#' distinct entries.
#'
#' @param G Genomic relationship matrix (row/column names) with at least
#'   `max(candidate_plots)` genotypes.
#' @param Sigma_E Environment covariance with at least `1 + other_n_env`
#'   environments.
#' @param candidate_plots Integer vector of candidate entry counts to evaluate
#'   per site (e.g. `seq(20, 200, by = 20)`).
#' @param scope `"single"` (default) optimises one generous site with the others
#'   fixed; `"all"` optimises every site together; `"subset"` varies only
#'   `focal_envs`.
#' @param focal_envs For `scope = "subset"`, environment names whose capacity is
#'   set to each candidate value.
#' @param site_capacities For `scope = "subset"`, a named vector covering every
#'   environment. Values for non-focal environments remain fixed.
#' @param n_env Total number of sites, used when `scope = "all"`. Defaults to
#'   `1 + other_n_env`.
#' @param other_n_env For `scope = "single"`: number of other sites in the trial
#'   (each at `other_capacity`). Default 2.
#' @param other_capacity For `scope = "single"`: entries tested at each of the
#'   other sites. Defaults to the median of `candidate_plots`; set it to the real
#'   capacities of the constrained sites for a realistic answer.
#' @param plots_per_entry Plots consumed per entry (for cost/plot accounting):
#'   a scalar or a named/vector value per environment. Default 1.
#' @param check_plots_per_site Check plots added at each site (e.g.
#'   `n_checks * n_blocks`), so `total_plots` reflects the real field size
#'   including checks, not just test entries. May be scalar or per-environment.
#'   Default 0. See
#'   [test_entry_capacity()] to convert an offered plot budget.
#' @param cost_per_plot Optional cost of one plot; enables a `gain_per_cost`
#'   column. Default `NULL`.
#' @param min_gain Minimum marginal accuracy gain per added plot required to keep
#'   growing the site (used by `select = "diminishing"`). Default `0.002`.
#' @param select How to pick the optimum **within the candidate interval**.
#'   `"max"` (default) chooses the capacity that **maximises accuracy** within the
#'   interval -- for an unrestricted site, the breeder just supplies the interval
#'   and the package returns the best size in it, with no further choice needed.
#'   `"diminishing"` returns the knee (where marginal gain drops below `min_gain`)
#'   if you would rather not over-invest, and `"target"` the smallest capacity
#'   reaching `target_accuracy`. `candidate_plots` *is* the interval.
#' @param target_accuracy Target for `select = "target"`. Default 0.7.
#' @param robust Optional list from [robust_scenarios()]. When supplied, every
#'   candidate capacity is simulated under every variance/covariance scenario.
#'   This allows unweighted environmental interaction kernels to affect the
#'   capacity recommendation even when no historical MET covariance is
#'   available.
#' @param robust_aggregate Aggregation across `robust` scenarios: `"min"`
#'   (default; maximin), `"mean"`, or `"cvar"`.
#' @param cvar_alpha Lower-tail fraction used when
#'   `robust_aggregate = "cvar"`.
#' @inheritParams simulate_met
#'
#' @return A list with `table` (one row per candidate: `plots` = per-site
#'   entries, `total_plots` across the trial, `accuracy_mean`, `gain_mean`,
#'   `marginal_accuracy_per_plot`, and optionally `gain_per_cost`) and
#'   `recommended_plots` (the diminishing-returns optimum; the per-site capacity
#'   to use everywhere when `scope = "all"`, at the focal site when
#'   `scope = "single"`, or at all focal sites when `scope = "subset"`), plus
#'   robust-analysis status, aggregation rule, and scenario count.
#'
#' @references
#' Colmant, A., Pita, F., & Covarrubias-Pazaran, G. (2026). Some objective
#' functions and ideas to optimize experimental designs in artificial selection
#' programs. *Crop Science*, 66, e70337.
#'
#' @seealso [allocate_sparse_met()], [sparsity_grid()], [recommend_replication()].
#' @examples
#' set.seed(1)
#' G <- crossprod(matrix(rnorm(80 * 150), 150, 80)) / 150 + diag(80) * 0.3
#' dimnames(G) <- list(paste0("L", 1:80), paste0("L", 1:80))
#' SigE <- diag(3) * 0.5 + 0.5
#' dimnames(SigE) <- list(paste0("E", 1:3), paste0("E", 1:3))
#' res <- suggest_site_capacity(G, SigE,
#'                              candidate_plots = c(20, 40, 60, 80),
#'                              other_n_env = 2, other_capacity = 30,
#'                              n_sim = 10, seed = 1)
#' res$recommended_plots
#' @export
suggest_site_capacity <- function(G, Sigma_E, candidate_plots,
                                  scope = c("single", "all", "subset"),
                                  n_env = NULL,
                                  other_n_env = 2L, other_capacity = NULL,
                                  focal_envs = NULL,
                                  site_capacities = NULL,
                                  plots_per_entry = 1,
                                  check_plots_per_site = 0,
                                  cost_per_plot = NULL,
                                  sigma_g2 = 1, sigma_e2 = 1,
                                  n_sim = 30L, select_fraction = 0.1,
                                  min_gain = 0.002,
                                  select = c("max", "diminishing", "target"),
                                  target_accuracy = 0.7,
                                  seed = NULL, max_dim = 6000L,
                                  robust = NULL,
                                  robust_aggregate = c("min", "mean", "cvar"),
                                  cvar_alpha = 0.25) {
  scope <- match.arg(scope)
  robust_aggregate <- match.arg(robust_aggregate)
  G <- as.matrix(G)
  Sigma_E <- as.matrix(Sigma_E)

  if (!is.numeric(G) || nrow(G) != ncol(G) || any(!is.finite(G)))
    stop("`G` must be a finite numeric square matrix.")
  if (is.null(rownames(G)) && is.null(colnames(G))) {
    rownames(G) <- colnames(G) <- paste0("L", seq_len(nrow(G)))
  } else if (is.null(rownames(G)) || is.null(colnames(G)) ||
             anyDuplicated(rownames(G)) || anyDuplicated(colnames(G)) ||
             !setequal(rownames(G), colnames(G))) {
    stop("`G` must have unique, matching row and column names.")
  } else {
    G <- G[rownames(G), rownames(G), drop = FALSE]
  }
  if (!is.numeric(Sigma_E) || nrow(Sigma_E) != ncol(Sigma_E) ||
      any(!is.finite(Sigma_E)) ||
      !isTRUE(all.equal(Sigma_E, t(Sigma_E), tolerance = 1e-8)))
    stop("`Sigma_E` must be a finite, symmetric numeric square matrix.")
  if (!is.numeric(candidate_plots) || !length(candidate_plots) ||
      any(!is.finite(candidate_plots)) ||
      any(abs(candidate_plots - round(candidate_plots)) > 1e-8))
    stop("`candidate_plots` must contain finite integer values.")
  candidate_plots <- sort(unique(as.integer(candidate_plots)))
  if (any(candidate_plots < 1L)) stop("`candidate_plots` must be positive.")
  if (max(candidate_plots) > nrow(G))
    stop("`G` has fewer genotypes than max(candidate_plots).")

  sigma_names <- rownames(Sigma_E)
  if (xor(is.null(rownames(Sigma_E)), is.null(colnames(Sigma_E))) ||
      (!is.null(sigma_names) &&
       (anyDuplicated(sigma_names) || anyDuplicated(colnames(Sigma_E)) ||
        !setequal(sigma_names, colnames(Sigma_E)))))
    stop("`Sigma_E` must have unique, matching row and column names.")
  if (!is.null(sigma_names))
    Sigma_E <- Sigma_E[sigma_names, sigma_names, drop = FALSE]

  if (scope == "all") {
    if (is.null(n_env)) n_env <- 1L + as.integer(other_n_env)
    if (length(n_env) != 1L || !is.finite(n_env) ||
        abs(n_env - round(n_env)) > 1e-8)
      stop("`n_env` must be one integer.")
    n_env <- as.integer(round(n_env))
    if (n_env < 2L) stop("`n_env` must be at least 2 for scope = 'all'.")
    if (n_env > nrow(Sigma_E))
      stop("`Sigma_E` needs at least `n_env` environments.")
    n_env_total <- n_env
    envs <- if (is.null(sigma_names)) paste0("E", seq_len(n_env_total)) else
      sigma_names[seq_len(n_env_total)]
    base_capacity <- rep(NA_integer_, n_env_total)
    focal_idx <- seq_len(n_env_total)
  } else if (scope == "single") {
    if (length(other_n_env) != 1L || !is.finite(other_n_env) ||
        abs(other_n_env - round(other_n_env)) > 1e-8 ||
        other_n_env < 0)
      stop("`other_n_env` must be one non-negative integer.")
    other_n_env <- as.integer(round(other_n_env))
    n_env_total <- 1L + other_n_env
    if (n_env_total > nrow(Sigma_E))
      stop("`Sigma_E` needs at least 1 + other_n_env environments.")
    if (is.null(other_capacity))
      other_capacity <- as.integer(round(stats::median(candidate_plots)))
    if (!is.numeric(other_capacity) || any(!is.finite(other_capacity)) ||
        any(abs(other_capacity - round(other_capacity)) > 1e-8) ||
        !(length(other_capacity) %in% c(1L, other_n_env)))
      stop("`other_capacity` must be one integer or one value per other environment.")
    other_capacity <- rep(as.integer(round(other_capacity)),
                          length.out = other_n_env)
    if (any(other_capacity < 1L | other_capacity > nrow(G)))
      stop("Every `other_capacity` must be between 1 and nrow(G).")
    envs <- if (is.null(sigma_names)) paste0("E", seq_len(n_env_total)) else
      sigma_names[seq_len(n_env_total)]
    base_capacity <- c(NA_integer_, other_capacity)
    focal_idx <- 1L
  } else {
    if (is.null(site_capacities) || !is.numeric(site_capacities) ||
        is.null(names(site_capacities)) || any(names(site_capacities) == "") ||
        anyDuplicated(names(site_capacities)) ||
        any(!is.finite(site_capacities)) ||
        any(abs(site_capacities - round(site_capacities)) > 1e-8) ||
        any(site_capacities < 1L | site_capacities > nrow(G)))
      stop("`site_capacities` must be a uniquely named integer vector between 1 and nrow(G).")
    if (is.null(focal_envs) || !length(focal_envs) ||
        anyNA(focal_envs) || anyDuplicated(focal_envs) ||
        !all(focal_envs %in% names(site_capacities)))
      stop("`focal_envs` must be unique names in `site_capacities`.")
    envs <- names(site_capacities)
    n_env_total <- length(envs)
    if (n_env_total > nrow(Sigma_E))
      stop("`Sigma_E` has fewer environments than `site_capacities`.")
    if (!is.null(sigma_names) && !all(envs %in% sigma_names))
      stop("Named `Sigma_E` does not contain every environment in `site_capacities`.")
    base_capacity <- as.integer(site_capacities)
    names(base_capacity) <- envs
    focal_idx <- match(focal_envs, envs)
  }

  normalise_by_environment <- function(x, arg, lower, strict = FALSE) {
    if (!is.numeric(x) || !length(x) || any(!is.finite(x)))
      stop(sprintf("`%s` must contain finite numeric values.", arg))
    if (!is.null(names(x)) && any(names(x) != "")) {
      if (anyDuplicated(names(x)) || !all(envs %in% names(x)))
        stop(sprintf("Named `%s` must cover every environment exactly once.", arg))
      x <- x[envs]
    } else if (length(x) == 1L) {
      x <- rep(x, n_env_total)
    } else if (length(x) != n_env_total) {
      stop(sprintf("`%s` must be scalar or have one value per environment.", arg))
    }
    bad <- if (strict) x <= lower else x < lower
    if (any(bad))
      stop(sprintf("`%s` values must be %s %s.", arg,
                   if (strict) "greater than" else "at least", lower))
    stats::setNames(as.numeric(x), envs)
  }

  plots_per_entry <- normalise_by_environment(
    plots_per_entry, "plots_per_entry", 0, strict = TRUE)
  check_plots_per_site <- normalise_by_environment(
    check_plots_per_site, "check_plots_per_site", 0)
  if (!is.null(cost_per_plot) &&
      (length(cost_per_plot) != 1L || !is.numeric(cost_per_plot) ||
       !is.finite(cost_per_plot) || cost_per_plot < 0))
    stop("`cost_per_plot` must be one finite non-negative number.")
  if (!is.numeric(min_gain) || length(min_gain) != 1L ||
      !is.finite(min_gain) || min_gain < 0)
    stop("`min_gain` must be one finite non-negative number.")
  if (!is.numeric(target_accuracy) || length(target_accuracy) != 1L ||
      !is.finite(target_accuracy) ||
      target_accuracy < -1 || target_accuracy > 1)
    stop("`target_accuracy` must be in [-1, 1].")
  if (!is.null(robust) && (!is.list(robust) || !length(robust)))
    stop("`robust` must be NULL or a non-empty list from `robust_scenarios()`.")
  if (!is.numeric(cvar_alpha) || length(cvar_alpha) != 1L ||
      !is.finite(cvar_alpha) || cvar_alpha < 0 || cvar_alpha > 1)
    stop("`cvar_alpha` must be one finite value in [0, 1].")

  if (!is.null(seed) &&
      (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
       abs(seed - round(seed)) > 1e-8))
    stop("`seed` must be one finite integer or NULL.")
  base_seed <- if (is.null(seed)) 1L else as.integer(round(seed))
  n_lines <- nrow(G)
  lines <- rownames(G)
  SigEsub <- if (is.null(sigma_names)) {
    Sigma_E[seq_len(n_env_total), seq_len(n_env_total), drop = FALSE]
  } else {
    Sigma_E[envs, envs, drop = FALSE]
  }
  dimnames(SigEsub) <- list(envs, envs)

  capacity_scenarios <- if (is.null(robust)) {
    list(list(
      sigma_g2 = sigma_g2, sigma_e2 = sigma_e2,
      sigmaE_shrink = 0, prob = 1, Sigma_E = SigEsub,
      Sigma_E_candidate = "supplied_central"
    ))
  } else robust
  valid_scenario <- vapply(capacity_scenarios, function(sc)
    is.list(sc) &&
      all(c("sigma_g2", "sigma_e2", "sigmaE_shrink", "prob") %in%
            names(sc)) &&
      all(vapply(sc[c("sigma_g2", "sigma_e2", "sigmaE_shrink", "prob")],
                 function(x) is.numeric(x) && length(x) == 1L &&
                   is.finite(x), logical(1))) &&
      sc$sigma_g2 > 0 && sc$sigma_e2 > 0 &&
      sc$sigmaE_shrink >= 0 && sc$sigmaE_shrink <= 1 && sc$prob >= 0,
    logical(1))
  if (!all(valid_scenario))
    stop("Every robust scenario needs finite positive variances, shrinkage ",
         "in [0, 1], and a finite non-negative probability.")
  scenario_probs <- vapply(capacity_scenarios, `[[`, numeric(1), "prob")
  if (sum(scenario_probs) <= 0)
    stop("Robust scenario probabilities must not all be zero.")
  scenario_probs <- scenario_probs / sum(scenario_probs)

  align_scenario_covariance <- function(sc) {
    S <- if (is.null(sc$Sigma_E)) Sigma_E else as.matrix(sc$Sigma_E)
    if (!is.numeric(S) || nrow(S) != ncol(S) || any(!is.finite(S)) ||
        !isTRUE(all.equal(S, t(S), tolerance = 1e-8)))
      stop("Every scenario `Sigma_E` must be a finite symmetric square matrix.")
    if (!is.null(rownames(S))) {
      if (is.null(colnames(S)) || anyDuplicated(rownames(S)) ||
          !setequal(rownames(S), colnames(S)) ||
          !all(envs %in% rownames(S)))
        stop("Every named scenario `Sigma_E` must contain all environments.")
      S <- S[envs, envs, drop = FALSE]
    } else {
      if (nrow(S) < n_env_total)
        stop("A scenario `Sigma_E` has too few environments.")
      S <- S[seq_len(n_env_total), seq_len(n_env_total), drop = FALSE]
      dimnames(S) <- list(envs, envs)
    }
    (1 - sc$sigmaE_shrink) * S +
      sc$sigmaE_shrink * diag(diag(S))
  }
  scenario_covariances <- lapply(capacity_scenarios,
                                 align_scenario_covariance)
  aggregate_capacity_metric <- function(values) {
    switch(
      robust_aggregate,
      min = min(values),
      mean = sum(values * scenario_probs),
      cvar = .cvar_lower(values, scenario_probs, cvar_alpha)
    )
  }

  # One fixed random order per environment makes the capacity sweep nested:
  # every larger candidate contains the entries evaluated at smaller candidates.
  set.seed(base_seed)
  entry_order <- lapply(seq_len(n_env_total),
                        function(e) sample.int(n_lines, n_lines))

  rows <- lapply(candidate_plots, function(n) {
    M <- matrix(0L, n_lines, n_env_total, dimnames = list(lines, envs))
    capacities <- base_capacity
    capacities[focal_idx] <- n
    for (e in seq_len(n_env_total))
      M[entry_order[[e]][seq_len(capacities[e])], e] <- 1L
    total_plots <- sum(capacities * plots_per_entry +
                         check_plots_per_site)

    scenario_metrics <- lapply(seq_along(capacity_scenarios), function(i) {
      scenario <- capacity_scenarios[[i]]
      sm <- simulate_met(
        M, G = G, Sigma_E = scenario_covariances[[i]],
        sigma_g2 = scenario$sigma_g2, sigma_e2 = scenario$sigma_e2,
        n_sim = n_sim, select_fraction = select_fraction,
        seed = base_seed, max_dim = max_dim
      )
      c(accuracy = sm$accuracy_mean, gain = sm$gain_mean)
    })
    scenario_metrics <- do.call(rbind, scenario_metrics)

    data.frame(plots = n,
               total_plots = total_plots,
               accuracy_mean =
                 aggregate_capacity_metric(scenario_metrics[, "accuracy"]),
               gain_mean =
                 aggregate_capacity_metric(scenario_metrics[, "gain"]),
               stringsAsFactors = FALSE)
  })
  tab <- do.call(rbind, rows)
  if (is.null(tab) || nrow(tab) < 1L)
    stop("No candidate capacities could be evaluated.")
  tab <- tab[order(tab$plots), , drop = FALSE]

  # Use the actual trial-wide plot increment: it correctly accounts for several
  # focal sites and environment-specific replication.
  d_acc  <- c(NA_real_, diff(tab$accuracy_mean))
  d_plot <- c(NA_real_, diff(tab$total_plots))
  tab$marginal_accuracy_per_plot <- d_acc / d_plot
  if (!is.null(cost_per_plot))
    tab$gain_per_cost <- tab$gain_mean / (tab$total_plots * cost_per_plot)

  # Recommendation within the user's candidate interval:
  #  "diminishing" (default) -- grow while each step buys >= min_gain accuracy
  #     per plot, then stop (the knee);
  #  "max" -- exploit free capacity: the candidate with the highest accuracy
  #     (i.e. the top of the interval you allowed);
  #  "target" -- the smallest capacity reaching `target_accuracy`.
  select <- match.arg(select)
  recommended <- switch(select,
    max = tab$plots[which.max(tab$accuracy_mean)],
    target = {
      ok <- tab$plots[is.finite(tab$accuracy_mean) &
                        tab$accuracy_mean >= target_accuracy]
      if (length(ok)) min(ok) else tab$plots[which.max(tab$accuracy_mean)]
    },
    diminishing = {
      rec <- tab$plots[1]
      if (nrow(tab) > 1L) for (i in 2:nrow(tab)) {
        if (is.finite(tab$marginal_accuracy_per_plot[i]) &&
            tab$marginal_accuracy_per_plot[i] >= min_gain) rec <- tab$plots[i]
        else break
      }
      rec
    })

  rownames(tab) <- NULL
  list(
    table = tab, recommended_plots = recommended, select = select,
    robust = !is.null(robust), robust_aggregate = robust_aggregate,
    n_scenarios = length(capacity_scenarios)
  )
}
