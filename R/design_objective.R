#' Evaluate a sparse MET design under a combined objective
#'
#' @description
#' A single, deterministic scoring function for a candidate allocation that
#' combines the three things a breeder trades off: statistical quality
#' (reliability / CDmean, or mean PEV), expected genetic gain (single-trait or a
#' multi-trait economic index), and resource cost (plots). It is built on the
#' coupled two-level [met_information()] matrix, so no Monte-Carlo simulation is
#' needed inside an optimisation loop -- the objective is exact given the
#' variance parameters. This is the objective maximised by [optimize_design()].
#'
#' @details
#' The scalar score (higher is better) is a weighted combination
#' \deqn{s = w_{gain}\,\tilde{\Delta G} + w_{rel}\,\tilde{r^2} - w_{cost}\,\tilde{c},}
#' where each component is divided by a reference value (`ref`) so the terms are
#' on a comparable scale; [optimize_design()] supplies `ref` from the starting
#' design automatically. A design whose plot count exceeds `budget` is infeasible
#' (`score = -Inf`). For a single trait, gain and reliability are monotonically
#' related, so weighting both is redundant; the distinct trade-off is
#' benefit-vs-cost (and the multi-trait index vs cost).
#'
#' @param allocation_matrix Genotype-by-environment 0/1 (or replication-count)
#'   matrix with dimnames.
#' @param G,Sigma_E,sigma_g2,sigma_e2,reps,env_efficiency,max_dim Passed to
#'   [met_information()] (across-TPE target).
#' @param prop,sigma_g,trait_weights,trait_gencov Passed to
#'   [expected_genetic_gain()]. Supplying `trait_weights` (e.g. from
#'   [desired_gain_weights()]) with `trait_gencov` (the trait genetic covariance
#'   \eqn{\Sigma_T}) turns on the multi-trait index objective.
#' @param R_T Trait residual covariance for the exact multi-trait path (defaults
#'   to `sigma_e2 * I`).
#' @param multitrait `"exact"` (default) uses the full trait-covariance Kronecker
#'   MME via [met_information_mt()] for the index reliability; `"approx"` uses the
#'   single-trait reliability times the index SD.
#' @param cost_per_plot Cost of one physical plot.
#'   Default 1.
#' @param weights Named list of non-negative weights `gain`, `reliability`,
#'   `cost`. Default `list(gain = 1, reliability = 0, cost = 0)`.
#' @param ref Named list of reference values (`gain`, `reliability`, `cost`) used
#'   to normalise the components; `NULL` (default) uses 1 for each.
#' @param budget Optional maximum number of plots; designs above it are
#'   infeasible.
#' @return A list with `score` and the raw components `reliability`, `mean_PEV`,
#'   `gain`, `plots`, `cost`, and `feasible`.
#' @seealso [optimize_design()], [expected_genetic_gain()], [met_information()].
#' @examples
#' set.seed(1)
#' G <- crossprod(matrix(rnorm(30 * 60), 60, 30)) / 60 + diag(30) * 0.3
#' dimnames(G) <- list(paste0("L", 1:30), paste0("L", 1:30))
#' SigE <- diag(3) * 0.5 + 0.5
#' dimnames(SigE) <- list(paste0("E", 1:3), paste0("E", 1:3))
#' M <- matrix(0L, 30, 3, dimnames = list(rownames(G), colnames(SigE)))
#' for (e in 1:3) M[sample(30, 15), e] <- 1L
#' design_objective(M, G, SigE,
#'                  weights = list(gain = 1, reliability = 0, cost = 0.2))$score
#' @export
design_objective <- function(allocation_matrix, G, Sigma_E = NULL,
                             sigma_g2 = 1, sigma_e2 = 1,
                             reps = NULL, env_efficiency = NULL,
                             prop = 0.1, sigma_g = 1,
                             trait_weights = NULL, trait_gencov = NULL,
                             R_T = NULL, multitrait = c("exact", "approx"),
                             cost_per_plot = 1,
                             weights = list(gain = 1, reliability = 0, cost = 0),
                             ref = NULL, budget = NULL, max_dim = 6000L) {
  multitrait <- match.arg(multitrait)
  is_mt <- !is.null(trait_weights) && !is.null(trait_gencov)
  if (!is.numeric(cost_per_plot) || length(cost_per_plot) != 1L ||
      !is.finite(cost_per_plot) || cost_per_plot < 0)
    stop("`cost_per_plot` must be one finite non-negative number.")
  if (!is.null(budget) &&
      (!is.numeric(budget) || length(budget) != 1L ||
       !is.finite(budget) || budget < 0))
    stop("`budget` must be NULL or one finite non-negative number.")
  if (!is.list(weights) || (length(weights) && is.null(names(weights))))
    stop("`weights` must be a named list.")
  supplied_weights <- weights[intersect(names(weights),
                                        c("gain", "reliability", "cost"))]
  if (any(vapply(supplied_weights, function(x)
    !is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0,
    logical(1))))
    stop("Objective weights must be finite non-negative scalars.")

  if (is_mt && multitrait == "exact") {
    # Exact multi-trait reliability from the trait-covariance Kronecker MME.
    mt <- met_information_mt(allocation_matrix, G = G, Sigma_E = Sigma_E,
                             Sigma_T = trait_gencov, R_T = R_T,
                             index_weights = trait_weights, sigma_e2 = sigma_e2,
                             reps = reps, env_efficiency = env_efficiency,
                             max_dim = max_dim)
    reliability <- mt$CDmean_index
    mean_PEV <- mt$mean_PEV_index
    gain <- selection_intensity(prop) *
      sqrt(max(0, min(1, reliability))) * mt$sigma_index
  } else {
    info <- met_information(allocation_matrix, G = G, Sigma_E = Sigma_E,
                            sigma_g2 = sigma_g2, sigma_e2 = sigma_e2,
                            reps = reps, env_efficiency = env_efficiency,
                            target = "across_tpe", max_dim = max_dim)
    reliability <- info$CDmean
    mean_PEV <- info$mean_PEV
    gain <- expected_genetic_gain(reliability = reliability, sigma_g = sigma_g,
                                  prop = prop, trait_weights = trait_weights,
                                  trait_gencov = trait_gencov)$gain
  }
  plot_matrix <- if (is.null(reps)) allocation_matrix else as.matrix(reps)
  if (!all(dim(plot_matrix) == dim(allocation_matrix)) ||
      !is.numeric(plot_matrix) || any(!is.finite(plot_matrix)) ||
      any(plot_matrix < 0))
    stop("Plot counts from `reps`/`allocation_matrix` must be a finite, ",
         "non-negative matrix matching the allocation dimensions.")
  plots <- sum(plot_matrix)
  cost  <- plots * cost_per_plot

  out <- list(reliability = reliability, mean_PEV = mean_PEV,
              gain = gain, plots = plots, cost = cost, feasible = TRUE)

  if (!is.null(budget) && plots > budget) {
    out$feasible <- FALSE
    out$score <- -Inf
    return(out)
  }
  out$score <- .combine_design_score(out, weights, ref)
  out
}


# Weighted, reference-normalised scalarisation (higher = better).
.combine_design_score <- function(components, weights, ref = NULL) {
  wv <- function(nm) if (is.null(weights[[nm]])) 0 else weights[[nm]]
  rv <- function(nm, d) if (is.null(ref) || is.null(ref[[nm]]) ||
                            !is.finite(ref[[nm]]) || ref[[nm]] == 0) d else ref[[nm]]
  g <- components$gain        / rv("gain", 1)
  r <- components$reliability / rv("reliability", 1)
  c <- components$cost        / rv("cost", 1)
  wv("gain") * g + wv("reliability") * r - wv("cost") * c
}
