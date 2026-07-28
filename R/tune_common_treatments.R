#' Tune the number of common treatments across environments
#'
#' @description
#' Common treatments (lines tested in every environment) trade off against unique
#' coverage: more of them stabilise the estimated genotype-by-environment
#' covariance and connect an otherwise-disconnected design, but consume plots
#' that could test unique lines. `tune_common_treatments()` sweeps the number of
#' common treatments and reports four curves so the optimum can be read off
#' directly rather than guessed.
#'
#' @details
#' The crucial point is that a design criterion which assumes the environment
#' covariance \eqn{\Sigma_E} is *known* almost always prefers zero common
#' treatments -- their value is that they let you *estimate* \eqn{\Sigma_E} and
#' connect the design. The function therefore reports, for each candidate count:
#' \describe{
#'   \item{`accuracy_known`}{Simulated accuracy with \eqn{\Sigma_E} assumed known
#'     (via [simulate_met()]); typically decreasing in the count -- the coverage
#'     cost.}
#'   \item{`accuracy_estimated`}{Simulated accuracy with \eqn{\Sigma_E}
#'     re-estimated from the data each replicate; this is the objective that can
#'     show an interior optimum, and it rises when correlations are weak (weak
#'     correlations are noisier to estimate, so more common treatments help).}
#'   \item{`min_shared`, `shared_var`}{Connectivity: the minimum and variance of
#'     the number of genotypes shared between environment pairs (the G x E
#'     estimation objective; increasing, with a knee).}
#'   \item{`coverage`}{Number of unique lines actually tested (decreasing).}
#' }
#' The gap between `accuracy_known` and `accuracy_estimated` is the diagnostic:
#' it is widest when environments are weakly correlated, and a
#' `best_estimated` that runs to the top of `common_values` is a signal to
#' partition the environments into mega-environments rather than add more common
#' treatments.
#'
#' @param treatments Character vector of candidate line IDs.
#' @param environments Character vector of environment names.
#' @param n_test_entries_per_environment Plots (entries) per environment, `k`.
#' @param G,Sigma_E Genomic relationship and environment covariance.
#' @param common_values Integer vector of common-treatment counts to evaluate
#'   (each `< n_test_entries_per_environment`).
#' @param sparse_allocation How the non-common lines are spread: `"random"`
#'   (random per environment) or `"disjoint"` (non-overlapping blocks -- the
#'   worst case for connectivity, where common treatments matter most).
#' @param sigma_g2,sigma_e2 Variance components.
#' @param n_sim Simulation replicates per count.
#' @param select_fraction Top fraction for the gain statistic.
#' @param seed,max_dim Reproducibility and size controls.
#'
#' @return A list with `table` (one row per count: `n_common`, `accuracy_known`,
#'   `accuracy_estimated`, `gain_estimated`, `min_shared`, `shared_var`,
#'   `coverage`) and `optima` (the counts maximising `accuracy_known` and
#'   `accuracy_estimated`, and the connectivity `knee`).
#'
#' @seealso [simulate_met()], [allocate_sparse_met()], [select_environments()].
#' @export
tune_common_treatments <- function(treatments, environments,
                                   n_test_entries_per_environment,
                                   G, Sigma_E, common_values,
                                   sparse_allocation = c("random", "disjoint"),
                                   sigma_g2 = 1, sigma_e2 = 1,
                                   n_sim = 40L, select_fraction = 0.1,
                                   seed = NULL, max_dim = 6000L) {
  sparse_allocation <- match.arg(sparse_allocation)
  treatments <- unique(as.character(treatments))
  environments <- unique(as.character(environments))
  J <- length(treatments); E <- length(environments); k <- n_test_entries_per_environment
  G <- as.matrix(G)[treatments, treatments, drop = FALSE]
  Sigma_E <- as.matrix(Sigma_E)
  if (any(common_values >= k)) stop("all `common_values` must be < n_test_entries_per_environment.")
  base_seed <- if (is.null(seed)) 1L else seed

  build <- function(C) {
    set.seed(base_seed + 101L * C)
    M <- matrix(0L, J, E, dimnames = list(treatments, environments))
    if (C > 0) M[seq_len(C), ] <- 1L
    ksp  <- k - C
    pool <- if (C < J) (C + 1L):J else integer(0)
    if (ksp > 0 && length(pool)) {
      if (sparse_allocation == "random") {
        for (e in seq_len(E)) M[sample(pool, min(ksp, length(pool))), e] <- 1L
      } else {                              # disjoint blocks
        for (e in seq_len(E)) {
          blk <- pool[((e - 1L) * ksp + 1L):(e * ksp)]
          blk <- blk[!is.na(blk) & blk <= J]
          if (length(blk)) M[blk, e] <- 1L
        }
      }
    }
    M
  }

  # Realized accuracy with Sigma_E RE-ESTIMATED from the data each replicate.
  sim_est <- function(M) {
    d       <- as.numeric(sweep(M, 2, 1 / sigma_e2, `*`))   # column-major (e-1)J+j
    present <- as.numeric(M) > 0
    gd      <- diag(G); Ginv <- tryCatch(solve(G), error = function(e) .pinv_sym_dense(G))
    Lc      <- t(chol(sigma_g2 * kronecker(Sigma_E, G)))
    n_sel   <- max(1L, floor(select_fraction * J))
    set.seed(base_seed)
    accs <- numeric(n_sim); gains <- numeric(n_sim)
    for (s in seq_len(n_sim)) {
      u  <- as.numeric(Lc %*% stats::rnorm(J * E))
      U  <- matrix(u, J, E)
      Y  <- matrix(NA_real_, J, E)
      pres_cells <- which(present == 1)
      Y[pres_cells] <- U[pres_cells] + stats::rnorm(length(pres_cells)) * sqrt(sigma_e2)
      # method-of-moments Sigma_E from shared genotypes
      Se <- diag(E)
      for (a in seq_len(E)) for (b in seq_len(E)) {
        sh <- which(!is.na(Y[, a]) & !is.na(Y[, b]))
        if (length(sh) >= 2) Se[a, b] <- mean(Y[sh, a] * Y[sh, b] / gd[sh])
        else if (a != b) Se[a, b] <- 0
      }
      Se <- (Se + t(Se)) / 2
      eg <- eigen(Se, symmetric = TRUE); Se <- eg$vectors %*% diag(pmax(eg$values, 1e-3), E) %*% t(eg$vectors)
      SeInv <- solve(Se)
      Cuu <- diag(d) + kronecker(SeInv, Ginv)
      ybar <- numeric(J * E); ybar[pres_cells] <- Y[pres_cells]   # observed phenotypes
      rhs <- numeric(J * E)
      env_of <- rep(seq_len(E), each = J)
      for (e in seq_len(E)) {
        idx <- which(env_of == e & present == 1); if (!length(idx)) next
        muhat <- sum(d[idx] * ybar[idx]) / sum(d[idx])
        Cuu[idx, idx] <- Cuu[idx, idx] - outer(d[idx], d[idx]) / sum(d[idx])
        rhs[idx] <- d[idx] * (ybar[idx] - muhat)
      }
      uhat <- as.numeric(tryCatch(solve(Cuu, rhs),
        error = function(e) .pinv_sym_dense(Cuu) %*% rhs))
      trueBV <- rowMeans(U); predBV <- rowMeans(matrix(uhat, J, E))
      accs[s]  <- stats::cor(predBV, trueBV)
      sel <- order(predBV, decreasing = TRUE)[seq_len(n_sel)]
      gains[s] <- mean(trueBV[sel]) - mean(trueBV)
    }
    c(acc = mean(accs), gain = mean(gains))
  }

  rows <- lapply(common_values, function(C) {
    M <- build(C)
    sm_known <- simulate_met(M, G = G, Sigma_E = Sigma_E, sigma_g2 = sigma_g2,
                             sigma_e2 = sigma_e2, n_sim = n_sim,
                             select_fraction = select_fraction, seed = base_seed,
                             max_dim = max_dim)
    est <- sim_est(M)
    shared <- crossprod(M)                     # E x E genotypes shared per pair
    off <- shared[upper.tri(shared)]
    data.frame(n_common = C,
               accuracy_known     = sm_known$accuracy_mean,
               accuracy_estimated = est[["acc"]],
               gain_estimated     = est[["gain"]],
               min_shared = min(off), shared_var = stats::var(off),
               coverage = sum(rowSums(M) > 0), stringsAsFactors = FALSE)
  })
  tab <- do.call(rbind, rows)

  # connectivity knee: where the marginal gain in min_shared drops below 20% of
  # the first step (a simple, transparent elbow rule).
  ms <- tab$min_shared; knee <- tab$n_common[1]
  if (length(ms) > 2) {
    dm <- diff(ms); thr <- 0.2 * max(dm[1], 1e-9)
    below <- which(dm < thr)
    if (length(below)) knee <- tab$n_common[below[1]]
  }
  list(
    table = tab,
    optima = list(
      best_known     = tab$n_common[which.max(tab$accuracy_known)],
      best_estimated = tab$n_common[which.max(tab$accuracy_estimated)],
      connectivity_knee = knee
    )
  )
}
