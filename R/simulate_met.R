#' Simulate a sparse MET and report realised accuracy and genetic gain
#'
#' @description
#' `simulate_met()` closes the loop between design and outcome. It simulates
#' genotype-by-environment effects with a genuine across-environment genetic
#' correlation structure (\eqn{u \sim \mathrm{MVN}(0, \sigma_g^2\,\Sigma_E
#' \otimes G)}), generates plot means for the cells an allocation actually
#' tests, predicts breeding values by GBLUP through the coupled MET information
#' matrix (see [met_information()]), and reports the realised **selection
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
#' @param bv_target Which breeding value the design is scored against (Colmant
#'   et al. 2026): `"across_tpe"` (default; average of all environments),
#'   `"environment_specific"` (a single environment), or `"mega_environment"`
#'   (average of a subset). This can change which design wins, so broad- and
#'   specific-adaptation strategies can be compared honestly.
#' @param target_envs Environment name(s) or index/indices defining the target
#'   for `"environment_specific"` (exactly one) or `"mega_environment"` (a
#'   subset). Ignored for `"across_tpe"`.
#' @param seed Optional RNG seed.
#'
#' @return A list with the realised accuracy (`accuracy_mean`, `_sd`, `_se`, and
#'   Monte Carlo `_ci95`), genetic gain / selection differential
#'   (`gain_mean`, `_sd`, `_se`, and `_ci95`), and the breeder-facing
#'   selection metrics of Mothukuri et al. (2025): `common_selected_mean`/`_sd`
#'   (fraction of the truly-best set that prediction also selects, higher is
#'   better) and `avg_rank_mean`/`_sd` (mean true rank of the predicted-selected
#'   lines; best possible is `(n_selected + 1) / 2`). Also `n_selected`, `n_sim`,
#'   and the design-based `mean_PEV` / `CDmean` from [met_information()].
#'
#' @seealso [met_information()], [allocate_sparse_met()].
#' @export
simulate_met <- function(allocation_matrix, G, Sigma_E = NULL,
                         sigma_g2 = 1, sigma_e2 = 1,
                         reps = NULL, env_efficiency = NULL,
                         n_sim = 50L, select_fraction = 0.1,
                         bv_target = c("across_tpe", "environment_specific",
                                       "mega_environment"),
                         target_envs = NULL,
                         seed = NULL, max_dim = 6000L) {

  bv_target <- match.arg(bv_target)
  if (!is.numeric(n_sim) || length(n_sim) != 1L || !is.finite(n_sim) ||
      n_sim < 2L || n_sim != as.integer(n_sim))
    stop("`n_sim` must be an integer >= 2.")
  n_sim <- as.integer(n_sim)
  if (!is.numeric(select_fraction) || length(select_fraction) != 1L ||
      !is.finite(select_fraction) || select_fraction <= 0 ||
      select_fraction > 1)
    stop("`select_fraction` must be in (0, 1].")
  if (!is.numeric(sigma_g2) || length(sigma_g2) != 1L ||
      !is.finite(sigma_g2) || sigma_g2 <= 0)
    stop("`sigma_g2` must be a finite positive scalar.")
  if (!is.numeric(sigma_e2) || length(sigma_e2) != 1L ||
      !is.finite(sigma_e2) || sigma_e2 <= 0)
    stop("`sigma_e2` must be a finite positive scalar.")

  info <- met_information(allocation_matrix, G = G, Sigma_E = Sigma_E,
                          sigma_g2 = sigma_g2, sigma_e2 = sigma_e2,
                          reps = reps, env_efficiency = env_efficiency,
                          target = "across_tpe", max_dim = max_dim)

  lines <- rownames(allocation_matrix)
  envs  <- colnames(allocation_matrix)
  J <- length(lines); E <- length(envs)
  if (J < 2L)
    stop("`simulate_met()` requires at least two genotypes.")

  # Breeding-value target as a weight vector over environments (Colmant et al.
  # 2026: broad vs specific adaptation). "across_tpe" averages all environments;
  # "environment_specific" targets one; "mega_environment" averages a subset.
  tgt_idx <- if (is.null(target_envs)) NULL else {
    ti <- if (is.character(target_envs)) match(target_envs, envs) else as.integer(target_envs)
    if (anyNA(ti)) stop("`target_envs` not found among the environment names.")
    ti
  }
  w <- switch(bv_target,
    across_tpe = rep(1 / E, E),
    environment_specific = {
      if (is.null(tgt_idx) || length(tgt_idx) != 1L)
        stop("`bv_target = \"environment_specific\"` needs exactly one `target_envs`.")
      wv <- numeric(E); wv[tgt_idx] <- 1; wv
    },
    mega_environment = {
      if (is.null(tgt_idx) || length(tgt_idx) < 1L)
        stop("`bv_target = \"mega_environment\"` needs `target_envs` (a subset).")
      wv <- numeric(E); wv[tgt_idx] <- 1 / length(tgt_idx); wv
    })

  if (is.null(Sigma_E)) {
    Sigma_E <- diag(E)
    dimnames(Sigma_E) <- list(envs, envs)
  }
  Sigma_E <- as.matrix(Sigma_E)
  if (!is.null(rownames(Sigma_E)) && !is.null(colnames(Sigma_E)))
    Sigma_E <- Sigma_E[envs, envs, drop = FALSE]
  if (is.null(reps)) {
    reps <- allocation_matrix
  } else if (!is.null(rownames(reps)) && !is.null(colnames(reps))) {
    reps <- reps[lines, envs, drop = FALSE]
  }
  reps <- as.matrix(reps)
  if (is.null(env_efficiency)) {
    env_efficiency <- stats::setNames(rep(1, E), envs)
  } else if (!is.null(names(env_efficiency))) {
    env_efficiency <- env_efficiency[envs]
  }
  Gsub <- as.matrix(G[lines, lines, drop = FALSE])

  d       <- as.numeric(sweep(reps, 2, env_efficiency / sigma_e2, `*`))
  present <- as.numeric(allocation_matrix) > 0
  env_of  <- rep(seq_len(E), each = J)

  # PD-safe inverse of the coefficient matrix: plain solve when well-conditioned,
  # else a jittered solve, else a symmetric pseudo-inverse (a singular Sigma_E
  # makes C_uu ill-conditioned).
  Cuu_inv <- tryCatch(solve(info$C_uu), error = function(e) {
    Cm <- info$C_uu
    tryCatch(solve(Cm + diag(1e-8 * mean(diag(Cm)), nrow(Cm))),
             error = function(e2) .pinv_sym_dense(Cm))
  })

  # Matrix square root of the prior GxE covariance for simulating u. Use the
  # Cholesky when possible; if the covariance is singular or non-PD (e.g. two
  # perfectly correlated environments, a rank-deficient Sigma_E), fall back to a
  # jittered Cholesky and then to a symmetric eigen-based root so simulation
  # never fails on redundant environments.
  covU <- sigma_g2 * kronecker(Sigma_E, Gsub)
  Lc <- tryCatch(t(chol(covU)), error = function(e) NULL)
  if (is.null(Lc)) {
    jit <- 1e-8 * mean(diag(covU))
    Lc <- tryCatch(t(chol(covU + diag(jit, nrow(covU)))), error = function(e) NULL)
  }
  if (is.null(Lc)) {
     e <- eigen((covU + t(covU)) / 2, symmetric = TRUE)
    vals <- pmax(e$values, 0)
    Lc <- e$vectors %*% (t(e$vectors) * sqrt(vals))   # symmetric PSD square root
  }

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
        abs(seed - round(seed)) > 1e-8)
      stop("`seed` must be one finite integer or NULL.")
    set.seed(as.integer(round(seed)))
  }
  n_sel <- max(1L, floor(select_fraction * J))
  accs      <- numeric(n_sim)
  gains     <- numeric(n_sim)
  common    <- numeric(n_sim)   # coincidence of selected sets (Mothukuri et al. 2025)
  avg_rank  <- numeric(n_sim)   # mean true rank of the predicted-selected lines
  pres_idx  <- which(present == 1)

  for (s in seq_len(n_sim)) {
    u      <- as.numeric(Lc %*% stats::rnorm(J * E))
    trueBV <- as.numeric(matrix(u, nrow = J, ncol = E) %*% w)   # target-weighted BV

    ybar <- numeric(J * E)
    ybar[pres_idx] <- u[pres_idx] +
      stats::rnorm(length(pres_idx)) * sqrt(1 / d[pres_idx])

    rhs <- numeric(J * E)
    for (e in seq_len(E)) {
      idx <- which(env_of == e & present == 1)
      if (!length(idx)) next
      muhat <- sum(d[idx] * ybar[idx]) / sum(d[idx])     # absorbed environment mean
      rhs[idx] <- d[idx] * (ybar[idx] - muhat)
    }

    uhat   <- as.numeric(Cuu_inv %*% rhs)
    predBV <- as.numeric(matrix(uhat, nrow = J, ncol = E) %*% w)

    accs[s]   <- stats::cor(predBV, trueBV)
    sel_pred  <- order(predBV, decreasing = TRUE)[seq_len(n_sel)]
    sel_true  <- order(trueBV, decreasing = TRUE)[seq_len(n_sel)]
    gains[s]  <- mean(trueBV[sel_pred]) - mean(trueBV)   # realized selection differential
    # Coincidence: fraction of the truly-best set that prediction also selects.
    common[s]   <- length(intersect(sel_pred, sel_true)) / n_sel
    # Average rank: mean true rank (1 = best) of the predicted-selected lines.
    rank_true   <- rank(-trueBV, ties.method = "average")
    avg_rank[s] <- mean(rank_true[sel_pred])
  }

  accuracy_se <- stats::sd(accs) / sqrt(n_sim)
  gain_se <- stats::sd(gains) / sqrt(n_sim)
  mc_multiplier <- stats::qt(0.975, df = n_sim - 1L)
  accuracy_ci95 <- mean(accs) + c(-1, 1) * mc_multiplier * accuracy_se
  accuracy_ci95 <- pmin(1, pmax(-1, accuracy_ci95))
  gain_ci95 <- mean(gains) + c(-1, 1) * mc_multiplier * gain_se

  list(
    accuracy_mean     = mean(accs),
    accuracy_sd       = stats::sd(accs),
    accuracy_se       = accuracy_se,
    accuracy_ci95     = accuracy_ci95,
    gain_mean         = mean(gains),
    gain_sd           = stats::sd(gains),
    gain_se           = gain_se,
    gain_ci95         = gain_ci95,
    # Breeder-facing outcome metrics (Mothukuri et al. 2025):
    common_selected_mean = mean(common),   # in [0, 1]; higher is better
    common_selected_sd   = stats::sd(common),
    avg_rank_mean        = mean(avg_rank),  # best possible = (n_sel + 1) / 2
    avg_rank_sd          = stats::sd(avg_rank),
    n_selected           = n_sel,
    n_sim         = n_sim,
    mean_PEV      = info$mean_PEV,
    CDmean        = info$CDmean
  )
}
