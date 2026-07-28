#' Exact multi-trait MET information matrix (trait-covariance Kronecker term)
#'
#' @description
#' The exact multi-trait extension of [met_information()]. The additive genetic
#' values form a genotype-by-environment-by-trait array with the separable
#' covariance \eqn{\Sigma_T \otimes \Sigma_E \otimes G}, and the within-plot
#' residuals have trait covariance \eqn{R_T}. Rather than approximating a
#' multi-trait index by the single-trait reliability (as
#' [expected_genetic_gain()] does by default), `met_information_mt()` solves the
#' full multi-trait mixed-model equations and returns the exact prediction-error
#' variance and reliability of a desired-gain / economic index.
#'
#' @details
#' When every plot measures all traits (the usual case), the multi-trait model
#' with separable structure decouples exactly by the **canonical transformation**
#' (Hayes & Hill 1981; Thompson): simultaneously diagonalising \eqn{\Sigma_T} and
#' \eqn{R_T} turns the problem into `T` independent single-trait analyses, each
#' solved by [met_information()] at the corresponding canonical genetic variance
#' \eqn{\lambda_k} (with unit residual variance). The canonical
#' prediction-error variances are recombined into the original-scale
#' trait PEV matrix per genotype, and hence the reliability of any linear index
#' \eqn{H = w'g}. This is exact (verified against a direct multi-trait MME solve),
#' not an approximation, and reuses the fast single-trait engine `T` times.
#'
#' @param allocation_matrix Genotype-by-environment 0/1 (or replication) matrix,
#'   shared across traits (every plot measures all traits).
#' @param G Genomic relationship matrix (row/column names).
#' @param Sigma_E Environment covariance (defaults to identity).
#' @param Sigma_T Trait **genetic** (co)variance matrix (`T x T`, dimnames
#'   optional).
#' @param R_T Trait **residual** (co)variance matrix (`T x T`); defaults to
#'   `sigma_e2 * I`.
#' @param index_weights Optional index weight vector `w` (length `T`), typically
#'   from [desired_gain_weights()]. At least one weight must be non-zero. If
#'   `NULL`, per-trait reliabilities are returned but no index summary.
#' @param sigma_e2 Residual scale used only when `R_T` is `NULL`. Default 1.
#' @param reps,env_efficiency,max_dim Passed to [met_information()].
#' @return A list with `CDmean_index` (mean index reliability), `rel_index_per_line`,
#'   `PEV_index_per_line`, `mean_PEV_index`, `CD_per_trait` (a `T`-vector of mean
#'   per-trait reliabilities), `sigma_index` (\eqn{\sqrt{w'\Sigma_T w}}),
#'   `lambda` (canonical genetic variances), and `Q` (canonical transform).
#' @references
#' Hayes, J. F., & Hill, W. G. (1981). Modification of estimates of parameters in
#' the construction of genetic selection indices ('bending'). *Biometrics*, 37,
#' 483-493. Thompson, R. (1976). The estimation of maternal genetic variances.
#' *Biometrics*, 32, 903-917.
#' @seealso [met_information()], [desired_gain_weights()], [expected_genetic_gain()].
#' @examples
#' set.seed(1)
#' G <- crossprod(matrix(rnorm(8 * 20), 20, 8)) / 20 + diag(8) * 0.4
#' dimnames(G) <- list(paste0("L", 1:8), paste0("L", 1:8))
#' SigE <- matrix(c(1, .5, .5, 1), 2); dimnames(SigE) <- list(c("E1","E2"), c("E1","E2"))
#' SigT <- matrix(c(1, .3, .3, .6), 2)
#' M <- matrix(0L, 8, 2, dimnames = list(rownames(G), colnames(SigE)))
#' for (e in 1:2) M[sample(8, 5), e] <- 1L
#' met_information_mt(M, G, SigE, Sigma_T = SigT, index_weights = c(1, -0.5))$CDmean_index
#' @export
met_information_mt <- function(allocation_matrix, G, Sigma_E = NULL,
                               Sigma_T, R_T = NULL, index_weights = NULL,
                               sigma_e2 = 1, reps = NULL, env_efficiency = NULL,
                               max_dim = 6000L) {
  Sigma_T <- as.matrix(Sigma_T)
  Tn <- nrow(Sigma_T)
  if (!is.numeric(Sigma_T) || Tn < 1L || ncol(Sigma_T) != Tn ||
      any(!is.finite(Sigma_T)) ||
      !isTRUE(all.equal(Sigma_T, t(Sigma_T), tolerance = 1e-8)))
    stop("`Sigma_T` must be a finite symmetric square matrix.")
  trait_names <- rownames(Sigma_T)
  if (!is.null(trait_names) || !is.null(colnames(Sigma_T))) {
    if (is.null(trait_names) || is.null(colnames(Sigma_T)) ||
        anyDuplicated(trait_names) || anyDuplicated(colnames(Sigma_T)) ||
        !setequal(trait_names, colnames(Sigma_T)))
      stop("Named `Sigma_T` must have unique matching row and column names.")
    Sigma_T <- Sigma_T[trait_names, trait_names, drop = FALSE]
  }
  if (!is.numeric(sigma_e2) || length(sigma_e2) != 1L ||
      !is.finite(sigma_e2) || sigma_e2 <= 0)
    stop("`sigma_e2` must be one finite positive number.")
  if (is.null(R_T)) R_T <- diag(sigma_e2, Tn)
  R_T <- as.matrix(R_T)
  if (!is.numeric(R_T) || nrow(R_T) != Tn || ncol(R_T) != Tn ||
      any(!is.finite(R_T)) ||
      !isTRUE(all.equal(R_T, t(R_T), tolerance = 1e-8)))
    stop("`R_T` must be finite, symmetric, and the same dimension as `Sigma_T`.")
  if (!is.null(trait_names) &&
      (!is.null(rownames(R_T)) || !is.null(colnames(R_T)))) {
    if (is.null(rownames(R_T)) || is.null(colnames(R_T)) ||
        anyDuplicated(rownames(R_T)) || anyDuplicated(colnames(R_T)) ||
        !setequal(trait_names, rownames(R_T)) ||
        !setequal(trait_names, colnames(R_T)))
      stop("Named `R_T` must cover the traits in `Sigma_T`.")
    R_T <- R_T[trait_names, trait_names, drop = FALSE]
  }
  validate_psd <- function(M, nm) {
    ee <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
    if (min(ee) < -1e-8 * max(abs(ee), 1))
      stop("`", nm, "` must be positive semidefinite.")
    .bend_pd(M)
  }
  Sigma_T <- validate_psd(Sigma_T, "Sigma_T")
  R_T <- validate_psd(R_T, "R_T")

  lines <- rownames(allocation_matrix)
  Gm <- as.matrix(G)
  if (!is.null(lines) && !is.null(rownames(Gm))) Gm <- Gm[lines, lines, drop = FALSE]
  J <- nrow(Gm)
  if (is.null(Sigma_E)) Sigma_E <- diag(ncol(allocation_matrix))
  Sigma_E <- as.matrix(Sigma_E)
  s <- mean(Sigma_E)                       # = omega' Sigma_E omega, omega = 1/E
  if (!is.finite(s) || s <= 0)
    stop("`Sigma_E` must imply positive across-TPE variance.")

  if (!is.null(index_weights)) {
    if (!is.numeric(index_weights) || length(index_weights) != Tn ||
        any(!is.finite(index_weights)))
      stop("`index_weights` must contain one finite value per trait.")
    if (!is.null(names(index_weights)) && !is.null(trait_names)) {
      if (anyDuplicated(names(index_weights)) ||
          !setequal(trait_names, names(index_weights)))
        stop("Named `index_weights` must cover every trait.")
      index_weights <- index_weights[trait_names]
    }
    if (all(index_weights == 0))
      stop("`index_weights` must contain at least one non-zero weight.")
  }

  # ---- canonical transformation: Q R_T Q' = I, Q Sigma_T Q' = diag(lambda) ----
  # Guard a singular/non-PD residual trait covariance by bending it to nearest PD.
  L    <- tryCatch(t(chol(R_T)), error = function(e) t(chol(.bend_pd(R_T))))
  Linv <- solve(L)
  B    <- Linv %*% Sigma_T %*% t(Linv)
  eig  <- eigen((B + t(B)) / 2, symmetric = TRUE)
  lambda <- eig$values
  if (any(lambda <= 0))
    stop("`Sigma_T` and `R_T` yield non-positive canonical variances; check they are valid covariance matrices.")
  Q    <- t(eig$vectors) %*% Linv
  Qinv <- solve(Q)

  # ---- single-trait PEV per canonical trait (raw across-TPE PEV) ----
  pev <- matrix(NA_real_, Tn, J)           # canonical trait x genotype
  for (k in seq_len(Tn)) {
    info_k <- met_information(allocation_matrix, G = G, Sigma_E = Sigma_E,
                              sigma_g2 = lambda[k], sigma_e2 = 1,
                              reps = reps, env_efficiency = env_efficiency,
                              target = "across_tpe", max_dim = max_dim)
    pev[k, ] <- info_k$PEV_per_line
  }

  # ---- recombine to original-scale trait PEV and index reliability ----
  rel_trait <- matrix(NA_real_, Tn, J)
  rel_index <- rep(NA_real_, J)
  pev_index <- rep(NA_real_, J)
  w <- if (!is.null(index_weights)) as.numeric(index_weights) else NULL
  for (i in seq_len(J)) {
    PEV_orig <- Qinv %*% diag(pev[, i], Tn) %*% t(Qinv)   # T x T, scale s*Gii baked in
    denom <- s * Gm[i, i]
    for (t in seq_len(Tn))
      rel_trait[t, i] <- pmin(1, pmax(
        0, 1 - PEV_orig[t, t] / (denom * Sigma_T[t, t])))
    if (!is.null(w)) {
      PEV_H <- as.numeric(t(w) %*% PEV_orig %*% w)
      VarH  <- as.numeric(denom * (t(w) %*% Sigma_T %*% w))
      rel_index[i] <- pmin(1, pmax(0, 1 - PEV_H / VarH))
      pev_index[i] <- PEV_H
    }
  }
  names(rel_index) <- names(pev_index) <- lines

  list(
    CDmean_index      = if (is.null(w)) NA_real_ else mean(rel_index),
    rel_index_per_line = rel_index,
    PEV_index_per_line = pev_index,
    mean_PEV_index     = if (is.null(w)) NA_real_ else mean(pev_index),
    CD_per_trait       = stats::setNames(rowMeans(rel_trait), trait_names),
    sigma_index        = if (is.null(w)) NA_real_ else
                           sqrt(as.numeric(t(w) %*% Sigma_T %*% w)),
    lambda             = lambda,
    Q                  = Q
  )
}
