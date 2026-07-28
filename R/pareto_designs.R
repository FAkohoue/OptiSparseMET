#' Trade-off (Pareto) frontier of designs over benefit vs cost
#'
#' @description
#' There is rarely a single best design once resources matter: a bigger design
#' predicts better but costs more plots. `pareto_designs()` maps that trade-off
#' by re-running [optimize_design()] across a range of cost weights (an
#' \eqn{\varepsilon}-constraint-style sweep), then flags the **non-dominated**
#' (Pareto-optimal) designs -- those for which no other design is both cheaper
#' and at least as good. The breeder reads the frontier and picks the point that
#' matches their budget, rather than committing to one weighting up front. The
#' sweep can itself be robust by passing [robust_scenarios()].
#'
#' @param allocation_matrix Starting allocation (dimnames required).
#' @param G,Sigma_E Relationship inputs.
#' @param cost_weights Vector of cost weights to sweep (benefit weight fixed at
#'   1). Larger values buy cheaper, smaller designs. Default `seq(0, 1, 0.25)`.
#' @param benefit `"gain"` (default) or `"reliability"` -- the benefit objective
#'   traded against cost.
#' @param preserve Move set for [optimize_design()]; `"none"` (default) lets plot
#'   count vary so cost is actually optimised.
#' @param objective Base objective settings (see [optimize_design()]); the
#'   benefit and cost weights are set by the sweep.
#' @param ... Further arguments passed to [optimize_design()] (e.g. `sigma_g2`,
#'   `robust`, `n_starts`, `iters`, `seed`).
#' @return A list with `frontier` (a data frame: `cost_weight`, `reliability`,
#'   `gain`, `plots`, `cost`, `score`, `pareto`) and `designs` (the list of
#'   allocation matrices, one per cost weight).
#' @seealso [optimize_design()], [design_objective()].
#' @examples
#' set.seed(1)
#' G <- crossprod(matrix(rnorm(24 * 60), 60, 24)) / 60 + diag(24) * 0.3
#' dimnames(G) <- list(paste0("L", 1:24), paste0("L", 1:24))
#' SigE <- diag(3) * 0.5 + 0.5
#' dimnames(SigE) <- list(paste0("E", 1:3), paste0("E", 1:3))
#' M <- matrix(0L, 24, 3, dimnames = list(rownames(G), colnames(SigE)))
#' for (e in 1:3) M[sample(24, 12), e] <- 1L
#' pf <- pareto_designs(M, G, SigE, cost_weights = c(0, 0.5, 1),
#'                      n_starts = 1, iters = 30, seed = 1)
#' pf$frontier
#' @export
pareto_designs <- function(allocation_matrix, G, Sigma_E = NULL,
                           cost_weights = seq(0, 1, by = 0.25),
                           benefit = c("gain", "reliability"),
                           preserve = "none", objective = list(), ...) {
  benefit <- match.arg(benefit)

  runs <- lapply(cost_weights, function(cw) {
    w <- list(gain = if (benefit == "gain") 1 else 0,
              reliability = if (benefit == "reliability") 1 else 0,
              cost = cw)
    obj <- utils::modifyList(objective, list(weights = w))
    optimize_design(allocation_matrix, G = G, Sigma_E = Sigma_E,
                    objective = obj, preserve = preserve, ...)
  })

  fr <- data.frame(
    cost_weight = cost_weights,
    reliability = vapply(runs, function(r) r$components$reliability, numeric(1)),
    gain        = vapply(runs, function(r) r$components$gain, numeric(1)),
    plots       = vapply(runs, function(r) r$components$plots, numeric(1)),
    cost        = vapply(runs, function(r) r$components$cost, numeric(1)),
    score       = vapply(runs, function(r) r$score, numeric(1)),
    stringsAsFactors = FALSE)

  ben <- if (benefit == "gain") fr$gain else fr$reliability
  fr$pareto <- .pareto_flag(ben, fr$cost)

  list(frontier = fr, designs = lapply(runs, function(r) r$allocation_matrix))
}


# Flag non-dominated points: maximise benefit, minimise cost. Point i is
# dominated if some j is >= in benefit and <= in cost, with at least one strict.
.pareto_flag <- function(benefit, cost) {
  n <- length(benefit)
  keep <- rep(TRUE, n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (j == i) next
      if (benefit[j] >= benefit[i] && cost[j] <= cost[i] &&
          (benefit[j] > benefit[i] || cost[j] < cost[i])) {
        keep[i] <- FALSE; break
      }
    }
  }
  keep
}
