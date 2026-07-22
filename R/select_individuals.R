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
#' @param target Optional character vector of genotype names over which CDmean is
#'   averaged. Defaults to all genotypes in `G`.
#' @param method `"cdmean_exchange"` (default) or `"random"` (control).
#' @param iters Number of exchange iterations.
#' @param seed Optional RNG seed.
#'
#' @return A list with `selected` (genotype names), `mean_CD`, and `method`.
#'
#' @seealso [select_environments()], [met_information()].
#' @export
select_individuals <- function(G, n_train, lambda = 1,
                               target = NULL,
                               method = c("cdmean_exchange", "random"),
                               iters = 500L, seed = NULL) {
  method <- match.arg(method)
  G <- as.matrix(G)
  p <- nrow(G)
  if (is.null(rownames(G))) rownames(G) <- colnames(G) <- paste0("G", seq_len(p))
  ids <- rownames(G)
  if (n_train < 1L || n_train >= p) stop("`n_train` must satisfy 1 <= n_train < nrow(G).")
  if (!is.null(seed)) set.seed(seed)

  Ginv <- tryCatch(solve(G), error = function(e) .pinv_sym_dense(G))
  gdiag <- diag(G)
  tgt <- if (is.null(target)) seq_len(p) else match(target, ids)
  if (anyNA(tgt)) stop("Some `target` names are not in rownames(G).")

  mean_cd <- function(sel) {
    nt  <- length(sel)
    inT <- numeric(p); inT[sel] <- 1
    ZMZ <- diag(inT) - outer(inT, inT) / nt
    Cinv <- tryCatch(solve(ZMZ + lambda * Ginv),
                     error = function(e) .pinv_sym_dense(ZMZ + lambda * Ginv))
    CD <- 1 - lambda * diag(Cinv) / gdiag
    mean(CD[tgt])
  }

  sel <- sample.int(p, n_train)

  if (method == "cdmean_exchange") {
    cur <- mean_cd(sel)
    for (i in seq_len(iters)) {
      out_pos <- sel[sample.int(length(sel), 1L)]
      in_pool <- setdiff(seq_len(p), sel)
      if (!length(in_pool)) break
      in_pos  <- in_pool[sample.int(length(in_pool), 1L)]
      cand    <- c(setdiff(sel, out_pos), in_pos)
      new     <- mean_cd(cand)
      if (new > cur) { sel <- cand; cur <- new }
    }
    score <- cur
  } else {
    score <- mean_cd(sel)
  }

  list(selected = ids[sort(sel)], mean_CD = score, method = method)
}
