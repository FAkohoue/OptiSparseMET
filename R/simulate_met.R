#' Simulate a sparse MET and report realized accuracy and genetic gain
#'
#' @description
#' `simulate_met()` closes the loop between design and outcome. It simulates
#' genotype-by-environment effects with a genuine across-environment genetic
#' correlation structure (\eqn{u \sim \mathrm{MVN}(0, \sigma_g^2\,\Sigma_E
#' \otimes G)}), generates plot means for the cells an allocation actually
#' tests, predicts breeding values by GBLUP through the coupled MET information
#' matrix (see [met_information()]), and reports the realized **selection
#' accuracy** (correlation between predicted and true across-TPE breeding
#' values) and **genetic gain** (mean true merit of the top fraction selected on
#' prediction). This is the quantity breeders care about, and it lets a design's
#' A / D / CDmean proxies be checked against the outcome they are meant to
#' improve.
#'
#' Unlike a no-G×E simulation, the correlation structure in `Sigma_E` means
#' cross-environment connectivity genuinely affects accuracy, so balanced,
#' sparse, and near-balanced allocations can be compared fairly.
#'
#' @inheritParams met_information
#' @param n_sim Number of independent simulation replicates.
#' @param select_fraction Top fraction selected (on predicted BV) for the
#'   genetic-gain statistic.
#' @param seed Optional RNG seed.
#'
#' @return A list with `accuracy_mean`, `accuracy_sd`, `gain_mean`, `gain_sd`,
#'   `n_sim`, and the design-based `mean_PEV` / `CDmean` from [met_information()].
#'
#' @seealso [met_information()], [allocate_sparse_met()].
#' @export
simulate_met <- function(allocation_matrix, G, Sigma_E = NULL,
                         sigma_g2 = 1, sigma_e2 = 1,
                         reps = NULL, env_efficiency = NULL,
                         n_sim = 50L, select_fraction = 0.1,
                         seed = NULL, max_dim = 6000L) {

  info <- met_information(allocation_matrix, G = G, Sigma_E = Sigma_E,
                          sigma_g2 = sigma_g2, sigma_e2 = sigma_e2,
                          reps = reps, env_efficiency = env_efficiency,
                          target = "across_tpe", max_dim = max_dim)

  lines <- rownames(allocation_matrix)
  envs  <- colnames(allocation_matrix)
  J <- length(lines); E <- length(envs)

  if (is.null(Sigma_E)) Sigma_E <- diag(E)
  Sigma_E <- as.matrix(Sigma_E)
  if (is.null(reps)) reps <- allocation_matrix
  reps <- as.matrix(reps)
  if (is.null(env_efficiency)) env_efficiency <- rep(1, E)
  Gsub <- as.matrix(G[lines, lines, drop = FALSE])

  d       <- as.numeric(sweep(reps, 2, env_efficiency / sigma_e2, `*`))
  present <- as.numeric(allocation_matrix) > 0
  env_of  <- rep(seq_len(E), each = J)

  Cuu_inv <- solve(info$C_uu)

  # Cholesky factor of the prior GxE covariance for simulating u.
  covU <- sigma_g2 * kronecker(Sigma_E, Gsub)
  Lc   <- t(chol(covU))                       # lower-triangular; Lc %*% z ~ N(0, covU)

  if (!is.null(seed)) set.seed(seed)
  n_sel <- max(1L, floor(select_fraction * J))
  accs  <- numeric(n_sim)
  gains <- numeric(n_sim)
  pres_idx <- which(present == 1)

  for (s in seq_len(n_sim)) {
    u      <- as.numeric(Lc %*% stats::rnorm(J * E))
    trueBV <- rowMeans(matrix(u, nrow = J, ncol = E))    # average over environments

    ybar <- numeric(J * E)
    ybar[pres_idx] <- u[pres_idx] +
      stats::rnorm(length(pres_idx)) * sqrt(sigma_e2 / d[pres_idx])

    rhs <- numeric(J * E)
    for (e in seq_len(E)) {
      idx <- which(env_of == e & present == 1)
      if (!length(idx)) next
      muhat <- sum(d[idx] * ybar[idx]) / sum(d[idx])     # absorbed environment mean
      rhs[idx] <- d[idx] * (ybar[idx] - muhat)
    }

    uhat   <- as.numeric(Cuu_inv %*% rhs)
    predBV <- rowMeans(matrix(uhat, nrow = J, ncol = E))

    accs[s]  <- stats::cor(predBV, trueBV)
    sel      <- order(predBV, decreasing = TRUE)[seq_len(n_sel)]
    gains[s] <- mean(trueBV[sel]) - mean(trueBV)
  }

  list(
    accuracy_mean = mean(accs),
    accuracy_sd   = stats::sd(accs),
    gain_mean     = mean(gains),
    gain_sd       = stats::sd(gains),
    n_sim         = n_sim,
    mean_PEV      = info$mean_PEV,
    CDmean        = info$CDmean
  )
}
