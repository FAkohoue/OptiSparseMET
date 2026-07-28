#' Criterion-driven optimisation of a sparse MET allocation
#'
#' @description
#' Searches for the genotype-by-environment allocation that maximises a
#' [design_objective()] -- any weighted combination of statistical quality
#' (reliability / CDmean), expected genetic gain (single-trait or multi-trait
#' index), and resource cost -- by simulated-annealing exchange with multiple
#' restarts. Unlike the greedy constructor in [allocate_sparse_met()], the
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
#' @param preserve `"margins"`, `"replication"`, or `"none"` (see Details).
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
#' @return A list with `allocation_matrix` (best found), `score`, `components`
#'   (raw reliability/gain/plots/cost of the best design), `score_start`,
#'   `trace` (best score per restart), `preserve`, `robust` (logical), and
#'   `seed_summary` when seed constraints are active.
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
                            preserve = c("margins", "replication", "none"),
                            robust = NULL,
                            robust_aggregate = c("mean", "cvar", "min"),
                            cvar_alpha = 0.25,
                            seed_available = NULL,
                            seed_required_per_environment = NULL,
                            minimum_seed_buffer = 0,
                            environment_capacities = NULL,
                            minimum_environment_entries = 1L,
                            n_starts = 5L, iters = 200L, cooling = 0.97,
                            seed = NULL, max_dim = 6000L, verbose = FALSE) {
  preserve <- match.arg(preserve)
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

  M0 <- allocation_matrix
  storage.mode(M0) <- "integer"
  if (is.null(rownames(M0)) || is.null(colnames(M0)))
    stop("`allocation_matrix` needs row and column names.")
  if (anyDuplicated(rownames(M0)) || anyDuplicated(colnames(M0)) ||
      anyNA(M0) || !all(M0 %in% c(0L, 1L)))
    stop("`allocation_matrix` must be a named, unique, finite 0/1 matrix.")
  if (any(colSums(M0) == 0L))
    stop("Every environment in `allocation_matrix` must contain an entry.")
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
    load_ok <- all(loads >= minimum_environment_entries &
                     loads <= environment_capacities)
    seed_ok <- is.null(seed_budget) ||
      all(seed_used(M) <= seed_budget + 1e-10)
    load_ok && seed_ok
  }
  if (!design_feasible(M0))
    stop("The starting allocation violates an environment or seed constraint.")

  obj <- utils::modifyList(
    list(weights = list(gain = 1, reliability = 0, cost = 0),
         prop = 0.1, sigma_g = 1, cost_per_plot = 1,
         trait_weights = NULL, trait_gencov = NULL,
         R_T = NULL, multitrait = "exact", budget = NULL),
    objective)

  # Reference for normalisation: components of the starting design at nominal
  # variance parameters.
  base <- design_objective(M0, G = G, Sigma_E = Sigma_E,
                           sigma_g2 = sigma_g2, sigma_e2 = sigma_e2,
                           prop = obj$prop, sigma_g = obj$sigma_g,
                           trait_weights = obj$trait_weights,
                           trait_gencov = obj$trait_gencov,
                           R_T = obj$R_T, multitrait = obj$multitrait,
                           cost_per_plot = obj$cost_per_plot,
                           weights = obj$weights, ref = NULL,
                           budget = NULL, max_dim = max_dim)
  ref <- list(gain = base$gain, reliability = base$reliability,
              cost = if (base$cost > 0) base$cost else 1)

  score_fun <- function(M) {
    if (!design_feasible(M)) return(-Inf)
    if (is.null(robust)) {
      design_objective(M, G = G, Sigma_E = Sigma_E,
                       sigma_g2 = sigma_g2, sigma_e2 = sigma_e2,
                       prop = obj$prop, sigma_g = obj$sigma_g,
                       trait_weights = obj$trait_weights,
                       trait_gencov = obj$trait_gencov,
                       R_T = obj$R_T, multitrait = obj$multitrait,
                       cost_per_plot = obj$cost_per_plot,
                       weights = obj$weights, ref = ref,
                       budget = obj$budget, max_dim = max_dim)$score
    } else {
      robust_design_score(M, G = G, Sigma_E = Sigma_E, scenarios = robust,
                          aggregate = robust_aggregate, cvar_alpha = cvar_alpha,
                          prop = obj$prop, sigma_g = obj$sigma_g,
                          trait_weights = obj$trait_weights,
                          trait_gencov = obj$trait_gencov,
                          R_T = obj$R_T, multitrait = obj$multitrait,
                          cost_per_plot = obj$cost_per_plot,
                          weights = obj$weights, ref = ref,
                          budget = obj$budget, max_dim = max_dim)$score
    }
  }

  score_start <- score_fun(M0)
  if (!is.finite(score_start))
    stop("The starting design is infeasible under the objective budget or ",
         "produces a non-finite score.")
  best_M <- M0; best_score <- score_start; trace <- numeric(0)

  for (st in seq_len(n_starts)) {
    # First restart begins at the supplied design; later ones are perturbed.
    cur_M <- if (st == 1L) M0 else .perturb_design(M0, preserve, obj$budget, 10L)
    if (!design_feasible(cur_M)) cur_M <- M0
    cur_s <- score_fun(cur_M)
    restart_best <- cur_s
    temp  <- max(abs(cur_s), 1e-3) * 0.5

    for (it in seq_len(iters)) {
      prop_M <- .propose_move(cur_M, preserve, obj$budget)
      if (is.null(prop_M)) { temp <- temp * cooling; next }
      if (!design_feasible(prop_M)) { temp <- temp * cooling; next }
      new_s <- score_fun(prop_M)
      d <- new_s - cur_s
      if (is.finite(new_s) && (d > 0 || stats::runif(1) < exp(d / temp))) {
        cur_M <- prop_M; cur_s <- new_s
        if (cur_s > restart_best) restart_best <- cur_s
        if (cur_s > best_score) { best_score <- cur_s; best_M <- cur_M }
      }
      temp <- temp * cooling
    }
    trace <- c(trace, restart_best)
    if (verbose) message(sprintf("restart %d: score = %.5f (best = %.5f)",
                                 st, cur_s, best_score))
  }

  comp <- design_objective(best_M, G = G, Sigma_E = Sigma_E,
                           sigma_g2 = sigma_g2, sigma_e2 = sigma_e2,
                           prop = obj$prop, sigma_g = obj$sigma_g,
                           trait_weights = obj$trait_weights,
                           trait_gencov = obj$trait_gencov,
                           R_T = obj$R_T, multitrait = obj$multitrait,
                           cost_per_plot = obj$cost_per_plot,
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

  list(allocation_matrix = best_M, score = best_score,
       components = comp[c("reliability", "mean_PEV", "gain", "plots", "cost")],
       score_start = score_start, trace = trace,
       preserve = preserve, robust = !is.null(robust),
       seed_summary = seed_summary)
}


# ---- move generators --------------------------------------------------------

# One random move respecting `preserve`; returns a modified matrix or NULL.
.propose_move <- function(M, preserve, budget) {
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
.perturb_design <- function(M, preserve, budget, n) {
  for (i in seq_len(n)) {
    m <- .propose_move(M, preserve, budget)
    if (!is.null(m)) M <- m
  }
  M
}
