#' Build a prior over variance components for robust design
#'
#' @description
#' Optimal designs are usually *locally* optimal -- best only at the assumed
#' variance parameters. A robust (Bayesian) design instead performs well across a
#' plausible range of those parameters. `robust_scenarios()` builds that range as
#' a weighted set of scenarios varying the residual/genetic variance ratio (the
#' heritability regime) and, optionally, the trust placed in the assumed
#' environment covariance `Sigma_E` (shrinking it toward independence to reflect
#' uncertain G x E correlations).
#'
#' @param sigma_e2 Vector of residual variances to consider (with `sigma_g2`,
#'   these set the heritability regime). Default `c(0.5, 1, 2)`.
#' @param sigma_g2 Genetic variance(s). Default 1.
#' @param sigmaE_shrink Vector of shrinkage weights in \[0, 1\]; a scenario uses
#'   `(1 - s) * Sigma_E + s * diag(diag(Sigma_E))`, so `s = 0` trusts the
#'   covariance fully and `s = 1` assumes independent environments. Default 0.
#' @param probs Optional prior probabilities (recycled/normalised); defaults to
#'   uniform.
#' @param Sigma_E_candidates Optional named list of alternative environment
#'   covariance matrices (for example calibration bootstrap or modality
#'   sensitivity candidates). Each is crossed with the variance and shrinkage
#'   grid. If omitted, the `Sigma_E` supplied to the scoring function is used.
#' @return A list of scenarios, each a list with `sigma_g2`, `sigma_e2`,
#'   `sigmaE_shrink`, and `prob`.
#' @seealso [robust_design_score()], [optimize_design()].
#' @examples
#' robust_scenarios(sigma_e2 = c(0.5, 1, 2), sigmaE_shrink = c(0, 0.5))
#' @export
robust_scenarios <- function(sigma_e2 = c(0.5, 1, 2), sigma_g2 = 1,
                             sigmaE_shrink = 0, probs = NULL,
                             Sigma_E_candidates = NULL) {
  if (!is.numeric(sigma_g2) || !length(sigma_g2) ||
      any(!is.finite(sigma_g2)) || any(sigma_g2 <= 0))
    stop("`sigma_g2` must contain finite positive values.")
  if (!is.numeric(sigma_e2) || !length(sigma_e2) ||
      any(!is.finite(sigma_e2)) || any(sigma_e2 <= 0))
    stop("`sigma_e2` must contain finite positive values.")
  if (!is.numeric(sigmaE_shrink) || !length(sigmaE_shrink) ||
      any(!is.finite(sigmaE_shrink)))
    stop("`sigmaE_shrink` must contain finite numeric values.")
  if (any(sigmaE_shrink < 0 | sigmaE_shrink > 1))
    stop("`sigmaE_shrink` must be in [0, 1].")
  candidates <- if (is.null(Sigma_E_candidates)) list(default = NULL) else
    Sigma_E_candidates
  if (!is.list(candidates) || !length(candidates))
    stop("`Sigma_E_candidates` must be NULL or a non-empty list.")
  candidate_ok <- vapply(candidates, function(M) {
    if (is.null(M)) return(TRUE)
    M <- as.matrix(M)
    if (!is.numeric(M) || nrow(M) != ncol(M) || any(!is.finite(M)) ||
        !isTRUE(all.equal(M, t(M), tolerance = 1e-8)) ||
        any(diag(M) <= 0))
      return(FALSE)
    ee <- eigen((M + t(M)) / 2, symmetric = TRUE,
                only.values = TRUE)$values
    min(ee) >= -1e-8 * max(1, max(abs(ee)))
  }, logical(1))
  if (!all(candidate_ok))
    stop("Every `Sigma_E_candidates` element must be a finite, symmetric, ",
         "positive-semidefinite square matrix with positive diagonal.")
  if (is.null(names(candidates)))
    names(candidates) <- paste0("candidate_", seq_along(candidates))
  grid <- expand.grid(sigma_g2 = sigma_g2, sigma_e2 = sigma_e2,
                      sigmaE_shrink = sigmaE_shrink,
                      candidate = seq_along(candidates),
                      KEEP.OUT.ATTRS = FALSE)
  n <- nrow(grid)
  if (is.null(probs)) probs <- rep(1, n)
  if (!is.numeric(probs) || !length(probs) || any(!is.finite(probs)))
    stop("`probs` must contain finite numeric values.")
  probs <- rep(probs, length.out = n)
  if (any(probs < 0) || sum(probs) == 0)
    stop("`probs` must be non-negative and not all zero.")
  probs <- probs / sum(probs)
  lapply(seq_len(n), function(i)
    list(sigma_g2 = grid$sigma_g2[i], sigma_e2 = grid$sigma_e2[i],
         sigmaE_shrink = grid$sigmaE_shrink[i], prob = probs[i],
         Sigma_E = candidates[[grid$candidate[i]]],
         Sigma_E_candidate = names(candidates)[grid$candidate[i]]))
}


#' Robust (Bayesian / maximin / CVaR) score of a design over a prior
#'
#' @description
#' Evaluates [design_objective()] under every scenario from [robust_scenarios()]
#' and aggregates the scores into one robust value. `"mean"` gives the
#' Bayes-optimal expected objective over the prior; `"min"` gives the worst-case
#' (maximin) objective; `"cvar"` gives the conditional value-at-risk -- the mean
#' of the worst `cvar_alpha` fraction of scenarios -- a tunable risk-averse
#' middle ground.
#'
#' @param allocation_matrix,G,Sigma_E The design and relationship inputs.
#' @param scenarios A list from [robust_scenarios()].
#' @param aggregate `"mean"` (Bayes), `"min"` (maximin), or `"cvar"`.
#' @param cvar_alpha Tail fraction for `"cvar"` (e.g. 0.25 = worst quarter).
#' @param ... Further arguments passed to [design_objective()] (e.g. `weights`,
#'   `ref`, `prop`, `cost_per_plot`, `budget`, `trait_weights`).
#' @return A list with `score` (aggregated), `aggregate`, and per-scenario
#'   `scenario_scores`, `scenario_probs`, and covariance-candidate labels.
#' @seealso [robust_scenarios()], [optimize_design()].
#' @export
robust_design_score <- function(allocation_matrix, G, Sigma_E, scenarios,
                                aggregate = c("mean", "cvar", "min"),
                                cvar_alpha = 0.25, ...) {
  aggregate <- match.arg(aggregate)
  if (!is.list(scenarios) || !length(scenarios))
    stop("`scenarios` must be a non-empty list from `robust_scenarios()`.")
  if (!is.numeric(cvar_alpha) || length(cvar_alpha) != 1L ||
      !is.finite(cvar_alpha) || cvar_alpha < 0 || cvar_alpha > 1)
    stop("`cvar_alpha` must be one finite value in [0, 1].")
  valid_scenario <- vapply(scenarios, function(sc)
    is.list(sc) &&
      all(c("sigma_g2", "sigma_e2", "sigmaE_shrink", "prob") %in% names(sc)) &&
      all(vapply(sc[c("sigma_g2", "sigma_e2", "sigmaE_shrink", "prob")],
                 function(x) is.numeric(x) && length(x) == 1L &&
                   is.finite(x), logical(1))) &&
      sc$sigma_g2 > 0 && sc$sigma_e2 > 0 &&
      sc$sigmaE_shrink >= 0 && sc$sigmaE_shrink <= 1 && sc$prob >= 0,
    logical(1))
  if (!all(valid_scenario))
    stop("Every scenario needs finite positive variances, shrinkage in [0, 1], ",
         "and a finite non-negative probability.")
  Sigma_E <- if (is.null(Sigma_E)) NULL else as.matrix(Sigma_E)

  scores <- vapply(scenarios, function(sc) {
    base_Sigma <- if (!is.null(sc$Sigma_E)) as.matrix(sc$Sigma_E) else Sigma_E
    if (!is.null(base_Sigma) && !is.null(Sigma_E) &&
        (!identical(dim(base_Sigma), dim(Sigma_E)) ||
         (!is.null(rownames(Sigma_E)) &&
          !setequal(rownames(base_Sigma), rownames(Sigma_E)))))
      stop("A scenario `Sigma_E` candidate is not aligned with `Sigma_E`.")
    if (!is.null(base_Sigma) && !is.null(Sigma_E) &&
        !is.null(rownames(Sigma_E)))
      base_Sigma <- base_Sigma[rownames(Sigma_E), colnames(Sigma_E),
                               drop = FALSE]
    SigE_sc <- if (is.null(base_Sigma)) NULL else
      (1 - sc$sigmaE_shrink) * base_Sigma +
      sc$sigmaE_shrink * diag(diag(base_Sigma))
    design_objective(allocation_matrix, G = G, Sigma_E = SigE_sc,
                     sigma_g2 = sc$sigma_g2, sigma_e2 = sc$sigma_e2, ...)$score
  }, numeric(1))
  probs <- vapply(scenarios, function(s) s$prob, numeric(1))
  if (sum(probs) <= 0)
    stop("Scenario probabilities must not all be zero.")
  probs <- probs / sum(probs)

  agg <- switch(aggregate,
    mean = sum(scores * probs),
    min  = min(scores),
    cvar = .cvar_lower(scores, probs, cvar_alpha))

  list(score = agg, aggregate = aggregate,
       scenario_scores = scores, scenario_probs = probs,
       scenario_covariance = vapply(
         scenarios,
         function(s) if (is.null(s$Sigma_E_candidate))
           "supplied_central" else as.character(s$Sigma_E_candidate),
         character(1)
       ))
}


# Conditional value-at-risk of the lower (worst) tail: probability-weighted mean
# of the worst `alpha` fraction of scores.
.cvar_lower <- function(scores, probs, alpha) {
  if (alpha <= 0) return(min(scores))
  if (alpha >= 1) return(sum(scores * probs))
  ord <- order(scores)                       # ascending: worst first
  s <- scores[ord]; p <- probs[ord]
  cum <- cumsum(p)
  k <- which(cum >= alpha)[1]
  if (is.na(k)) k <- length(s)
  w <- p[seq_len(k)]
  w[k] <- w[k] - (cum[k] - alpha)            # trim the last slice to reach alpha
  sum(s[seq_len(k)] * w) / alpha
}
