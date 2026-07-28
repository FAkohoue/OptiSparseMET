#' Select a representative training subset of the TPG by CDmean (decision 2)
#'
#' @description
#' `select_individuals()` chooses `n_train` genotypes that best represent the
#' target population of genotypes (TPG) for genomic prediction, by maximising the
#' mean coefficient of determination (CDmean; Rincent et al. 2012) for predicting
#' a target set. It is the principled tool for decision 2 in Colmant et al.
#' (2026), replacing random sampling.
#'
#' @details
#' For a candidate training set \eqn{T} (one observation per training genotype),
#' with \eqn{\lambda = \sigma_e^2/\sigma_g^2}, centring matrix
#' \eqn{M = I - \mathbf{1}\mathbf{1}^\top/|T|}, and incidence \eqn{Z} mapping
#' training observations to all genotypes,
#' \deqn{\mathrm{CD}_i = 1 - \lambda\,[(Z^\top M Z + \lambda G^{-1})^{-1}]_{ii}
#'   / G_{ii}.}
#' The mean CD over the `target` genotypes is maximised by an exchange algorithm
#' (swap one training genotype out, one in, keep improving swaps).
#'
#' @param G Genomic (or pedigree) relationship matrix with row/column names.
#' @param n_train Number of genotypes to select for the training/calibration set.
#' @param lambda Variance ratio \eqn{\sigma_e^2/\sigma_g^2}. Default 1.
#' @param target Optional character vector of genotype names over which the
#'   criterion is evaluated. Defaults to all genotypes in `G`.
#' @param criterion Training-set relationship measurement to optimise (Mothukuri
#'   et al. 2025; STPGA, Akdemir 2017). One of `"cdmean"` (default), `"cdmax"`,
#'   `"pevmean"`, `"pevmax"`, `"aopt"` (trace of PEV), `"dopt"` (log-determinant
#'   of PEV), `"goptpev"` (largest PEV eigenvalue), or `"neg_dist"` (negative
#'   mean genetic distance to the training set). CD-type and `neg_dist` are
#'   maximised; the PEV/A/D/GOPT criteria are minimised. Note (per the paper)
#'   that these criteria predict realised accuracy only weakly, so validate the
#'   chosen criterion for your programme with [validate_design_criteria()].
#' @param method `"exchange"` (default; searches by swapping members) or
#'   `"random"` (control). `"cdmean_exchange"` is accepted as a deprecated alias.
#' @param iters Number of exchange iterations.
#' @param seed Optional RNG seed.
#'
#' @return A list with `selected` (genotype names), the optimised `criterion` and
#'   its `score`, `measures` (all eight relationship measurements for the chosen
#'   set), `mean_CD` (kept for backward compatibility), and `method`.
#'
#' @seealso [select_environments()], [met_information()].
#' @export
select_individuals <- function(G, n_train, lambda = 1,
                               target = NULL,
                               criterion = c("cdmean", "cdmax", "pevmean",
                                             "pevmax", "aopt", "dopt", "goptpev",
                                             "neg_dist"),
                               method = c("exchange", "random",
                                          "cdmean_exchange"),
                               iters = 500L, seed = NULL) {
  criterion <- match.arg(criterion)
  method    <- match.arg(method)
  if (method == "cdmean_exchange") method <- "exchange"  # back-compatible alias
  G <- as.matrix(G)
  p <- nrow(G)
  if (is.null(rownames(G))) rownames(G) <- colnames(G) <- paste0("G", seq_len(p))
  ids <- rownames(G)
  if (n_train < 1L || n_train >= p) stop("`n_train` must satisfy 1 <= n_train < nrow(G).")
  if (!is.null(seed)) set.seed(seed)

  Ginv  <- tryCatch(solve(G), error = function(e) .pinv_sym_dense(G))
  gdiag <- diag(G)
  # Genetic distance matrix for the distance-based criterion.
  Dgen  <- sqrt(pmax(outer(gdiag, gdiag, "+") - 2 * G, 0))
  tgt <- if (is.null(target)) seq_len(p) else match(target, ids)
  if (anyNA(tgt)) stop("Some `target` names are not in rownames(G).")

  # All eight training-set relationship measurements (Mothukuri et al. 2025;
  # STPGA, Akdemir 2017) for a candidate training set `sel`, evaluated on the
  # target set. PEV is the target block of (Z'MZ + lambda G^{-1})^{-1}.
  criteria <- function(sel) {
    nt  <- length(sel)
    inT <- numeric(p); inT[sel] <- 1
    ZMZ <- diag(inT) - outer(inT, inT) / nt
    Cinv <- tryCatch(solve(ZMZ + lambda * Ginv),
                     error = function(e) .pinv_sym_dense(ZMZ + lambda * Ginv))
    PEV  <- Cinv[tgt, tgt, drop = FALSE]
    dpev <- diag(PEV)
    CD   <- 1 - lambda * dpev / gdiag[tgt]
    # mean distance from each target individual to the training set
    meandist <- mean(rowMeans(Dgen[tgt, sel, drop = FALSE]))
    list(
      cdmean  = mean(CD),
      cdmax   = max(CD),
      pevmean = mean(dpev),
      pevmax  = max(dpev),
      aopt    = sum(dpev),                                   # trace of PEV
      dopt    = as.numeric(determinant(PEV, logarithm = TRUE)$modulus),
      goptpev = max(eigen(PEV, symmetric = TRUE, only.values = TRUE)$values),
      neg_dist = -meandist
    )
  }

  # Score to MAXIMISE: CD-type and neg_dist are maximised as-is; PEV/A/D/GOPT
  # type criteria are minimised, so negate them.
  maximise_as_is <- criterion %in% c("cdmean", "cdmax", "neg_dist")
  score_of <- function(cr) if (maximise_as_is) cr[[criterion]] else -cr[[criterion]]

  sel <- sample.int(p, n_train)
  cr  <- criteria(sel)

  if (method == "exchange") {
    cur <- score_of(cr)
    for (i in seq_len(iters)) {
      out_pos <- sel[sample.int(length(sel), 1L)]
      in_pool <- setdiff(seq_len(p), sel)
      if (!length(in_pool)) break
      in_pos  <- in_pool[sample.int(length(in_pool), 1L)]
      cand    <- c(setdiff(sel, out_pos), in_pos)
      cand_cr <- criteria(cand)
      new     <- score_of(cand_cr)
      if (new > cur) { sel <- cand; cur <- new; cr <- cand_cr }
    }
  }

  list(
    selected   = ids[sort(sel)],
    criterion  = criterion,
    score      = cr[[criterion]],
    measures   = cr,          # all eight relationship measurements for the chosen set
    mean_CD    = cr$cdmean,   # backward-compatible field
    method     = method
  )
}
