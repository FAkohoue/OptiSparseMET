#' Robust consensus of several relationship matrices
#'
#' @description
#' Given the same relationship matrix estimated over several years (or replicates,
#' data sources), a *consensus* summarises them into one. The naive choice -- the
#' arithmetic (Euclidean) mean -- is neither robust to an anomalous year nor
#' geometrically appropriate for covariance/relationship matrices, which live on
#' the manifold of symmetric positive-definite (SPD) matrices. `consensus_relationship()`
#' offers better options:
#' \describe{
#'   \item{`"rv_weighted"` (default)}{A STATIS-style consensus: each matrix is
#'     weighted by how much it agrees with the others (the leading eigenvector of
#'     the between-matrix RV-coefficient matrix), so an outlier year is
#'     automatically down-weighted. The **robust** choice.}
#'   \item{`"geometric"`}{The log-Euclidean mean \eqn{\exp(\overline{\log D_i})},
#'     the principled mean on the SPD manifold (Arsigny et al. 2007); avoids the
#'     variance "swelling" of Euclidean averaging.}
#'   \item{`"median"`}{Element-wise median across matrices (resistant to
#'     outliers), projected to the nearest positive-definite matrix.}
#'   \item{`"mean"`}{The plain arithmetic mean (Euclidean), provided for
#'     comparison.}
#' }
#'
#' @param matrices A list of relationship matrices over the same individuals.
#' @param method `"rv_weighted"`, `"geometric"`, `"median"`, or `"mean"`.
#' @param eps Eigenvalue floor for the SPD operations. Default 1e-8.
#' @param use_off_diagonal For STATIS weighting, compare only unique
#'   off-diagonal relationships. This prevents a common unit diagonal from
#'   artificially inflating agreement.
#' @param center Centre the vectorised relationships before computing the
#'   STATIS/RV agreement. Default `TRUE`.
#' @return The consensus matrix. For `"rv_weighted"` the per-matrix weights are
#'   attached as `attr(., "weights")`.
#' @references
#' Arsigny, V., Fillard, P., Pennec, X., Ayache, N. (2007). Geometric means in a
#' novel vector space structure on symmetric positive-definite matrices.
#' Lavit, C. et al. (1994). The ACT (STATIS method).
#' @seealso [assess_envirotype_stability()], [tune_relationship_matrix()].
#' @examples
#' set.seed(1)
#' mk <- function() { A <- matrix(rnorm(5 * 8), 5); tcrossprod(A) / 8 + diag(5) }
#' Ds <- list(mk(), mk(), mk())
#' Cw <- consensus_relationship(Ds, "rv_weighted")
#' attr(Cw, "weights")
#' @export
consensus_relationship <- function(matrices,
                                   method = c("rv_weighted", "geometric",
                                              "median", "mean"),
                                   eps = 1e-8,
                                   use_off_diagonal = TRUE,
                                   center = TRUE) {
  method <- match.arg(method)
  if (!length(matrices)) stop("`matrices` is empty.")
  if (!is.numeric(eps) || length(eps) != 1L ||
      !is.finite(eps) || eps <= 0)
    stop("`eps` must be one finite positive number.")
  mats <- lapply(matrices, function(m) {
    m <- as.matrix(m)
    if (!is.numeric(m) || nrow(m) != ncol(m) || any(!is.finite(m)))
      stop("Every relationship matrix must be finite, numeric, and square.")
    if (!isTRUE(all.equal(m, t(m), tolerance = 1e-8)))
      stop("Every relationship matrix must be symmetric.")
    ee <- eigen(m, symmetric = TRUE, only.values = TRUE)$values
    if (min(ee) < -1e-8 * max(abs(ee), 1))
      stop("Every relationship matrix must be positive semidefinite.")
    if (sum(m * m) <= 0)
      stop("Relationship matrices must contain non-zero information.")
    m
  })
  ids <- rownames(mats[[1]])
  if (!is.null(ids)) {
    if (is.null(colnames(mats[[1]])) || anyDuplicated(ids) ||
        anyDuplicated(colnames(mats[[1]])) ||
        !setequal(ids, colnames(mats[[1]])))
      stop("Named matrices need unique, matching row and column names.")
    mats <- lapply(mats, function(m) {
      if (is.null(rownames(m)) || is.null(colnames(m)) ||
          !setequal(ids, rownames(m)) || !setequal(ids, colnames(m)))
        stop("All matrices must cover exactly the same individuals.")
      m[ids, ids, drop = FALSE]
    })
  } else {
    if (any(vapply(mats, function(m)
      !is.null(rownames(m)) || !is.null(colnames(m)), logical(1))))
      stop("Either all consensus matrices must be named or none may be named.")
    dims <- vapply(mats, nrow, integer(1))
    if (any(dims != dims[1L]))
      stop("Unnamed consensus matrices must have identical dimensions.")
  }
  k <- length(mats)

  C <- switch(method,
    mean = Reduce(`+`, mats) / k,
    median = {
      arr <- simplify2array(mats)
      .bend_pd(apply(arr, c(1, 2), stats::median, na.rm = TRUE), eps)
    },
    geometric = {
      logs <- lapply(mats, function(m) .logm_spd(m, eps))
      .expm_sym(Reduce(`+`, logs) / k)
    },
    rv_weighted = {
      if (k == 1L) {
        out <- mats[[1L]]
        attr(out, "weights") <- stats::setNames(1, names(mats))
        out
      } else {
      vectorise <- function(A) {
        raw <- if (use_off_diagonal) A[upper.tri(A)] else as.numeric(A)
        z <- raw
        if (center && length(z) > 1L && stats::sd(z) > 0)
          z <- z - mean(z)
        # Degenerate all-zero off-diagonals contain no agreement information;
        # use the diagonal only as a documented fallback rather than failing.
        if (!length(z) || sum(z * z) <= .Machine$double.eps)
          z <- diag(A)
        z
      }
      vv <- lapply(mats, vectorise)
      rv <- function(A, B) sum(A * B) / sqrt(sum(A * A) * sum(B * B))
      S <- diag(1, k)
      for (i in seq_len(k)) for (j in seq_len(k)) if (i < j)
        S[i, j] <- S[j, i] <- rv(vv[[i]], vv[[j]])
      w <- abs(eigen(S, symmetric = TRUE)$vectors[, 1]); w <- w / sum(w)
      out <- Reduce(`+`, Map(function(m, wi) wi * m, mats, w))
      attr(out, "weights") <- stats::setNames(w, names(mats))
      attr(out, "rv_agreement") <- S
      out
      }
    })

  wts <- attr(C, "weights")
  rv_agreement <- attr(C, "rv_agreement")
  if (!is.null(ids)) dimnames(C) <- list(ids, ids)
  C <- (C + t(C)) / 2
  if (!is.null(wts)) attr(C, "weights") <- wts
  if (!is.null(rv_agreement)) attr(C, "rv_agreement") <- rv_agreement
  C
}


# ---- SPD helpers ------------------------------------------------------------

.logm_spd <- function(M, eps = 1e-8) {
  e <- eigen((M + t(M)) / 2, symmetric = TRUE)
  v <- pmax(e$values, eps)
  e$vectors %*% (log(v) * t(e$vectors))
}

.expm_sym <- function(M) {
  e <- eigen((M + t(M)) / 2, symmetric = TRUE)
  e$vectors %*% (exp(e$values) * t(e$vectors))
}

.bend_pd <- function(M, eps = 1e-8) {
  dn <- dimnames(M)
  e <- eigen((M + t(M)) / 2, symmetric = TRUE)
  v <- pmax(e$values, eps * max(abs(e$values), 1))
  out <- e$vectors %*% (v * t(e$vectors))
  dimnames(out) <- dn
  out
}
