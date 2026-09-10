#' Fit a genetic environment covariance from historical MET responses
#'
#' @description
#' Fits an environment covariance directly from replicated or adjusted
#' genotype-by-environment responses by restricted maximum likelihood (REML).
#' Missing genotype-environment cells are handled in the likelihood, rather
#' than by pairwise correlations. Available covariance models are diagonal,
#' factor analytic, and unstructured.
#'
#' @param data Long-format historical MET data.
#' @param genotype_col,environment_col,response_col Column names.
#' @param environments Optional environment order; defaults to observed order.
#' @param model Covariance model: `"fa"`, `"unstructured"`, or `"diagonal"`.
#' @param rank Factor-analytic rank, used for `model = "fa"`.
#' @param residual_variance Known residual variance: a scalar or one value per
#'   environment. If `NULL`, it is estimated from within-cell replication;
#'   this requires at least one replicated genotype-environment cell.
#' @param control Named list passed to [stats::nlminb()].
#' @return A `historical_met_fit` object containing `Sigma_E`, its correlation
#'   form, residual variances, FA loadings/specific variances when relevant,
#'   REML log likelihood, information criteria, convergence diagnostics, and
#'   the aggregated data used by the fit.
#' @export
fit_historical_met <- function(
    data, genotype_col = "Genotype", environment_col = "Environment",
    response_col = "Value", environments = NULL,
    model = c("fa", "unstructured", "diagonal"), rank = 2L,
    residual_variance = NULL,
    control = list(iter.max = 500L, eval.max = 1000L,
                   rel.tol = 1e-8)) {
  model <- match.arg(model)
  if (!is.data.frame(data) ||
      !all(c(genotype_col, environment_col, response_col) %in% names(data)))
    stop("`data` does not contain the requested genotype, environment, and response columns.")
  d <- data.frame(
    genotype = as.character(data[[genotype_col]]),
    environment = as.character(data[[environment_col]]),
    response = as.numeric(data[[response_col]]), stringsAsFactors = FALSE)
  d <- d[!is.na(d$genotype) & nzchar(d$genotype) &
           !is.na(d$environment) & nzchar(d$environment) &
           is.finite(d$response), , drop = FALSE]
  if (!nrow(d)) stop("No finite historical MET observations remain.")
  if (is.null(environments)) environments <- unique(d$environment)
  environments <- unique(as.character(environments))
  if (length(environments) < 2L || anyNA(environments) ||
      any(!nzchar(environments)))
    stop("At least two valid environments are required.")
  d <- d[d$environment %in% environments, , drop = FALSE]
  E <- length(environments)
  if (!all(environments %in% d$environment))
    stop("Every requested environment must have observations.")
  if (!is.numeric(rank) || length(rank) != 1L || !is.finite(rank) ||
      rank != as.integer(rank) || rank < 1L || rank >= E)
    stop("`rank` must be an integer in 1, ..., E - 1.")
  rank <- as.integer(rank)

  key <- interaction(d$genotype, d$environment, drop = TRUE,
                     lex.order = TRUE, sep = "\r")
  split_idx <- split(seq_len(nrow(d)), key)
  split_y <- lapply(split_idx, function(i) d$response[i])
  meta <- d[vapply(split_idx, `[`, integer(1), 1L),
            c("genotype", "environment"), drop = FALSE]
  agg <- data.frame(
    genotype = meta$genotype, environment = meta$environment,
    value = vapply(split_y, mean, numeric(1)),
    n = vapply(split_y, length, integer(1)),
    within_var = vapply(split_y, function(x)
      if (length(x) > 1L) stats::var(x) else NA_real_, numeric(1)),
    stringsAsFactors = FALSE)
  agg$env_index <- match(agg$environment, environments)

  if (is.null(residual_variance)) {
    ok <- agg$n > 1L & is.finite(agg$within_var)
    df <- sum(agg$n[ok] - 1L)
    if (df <= 0L)
      stop("`residual_variance` is required when historical cells are unreplicated.")
    pooled <- max(sum((agg$n[ok] - 1L) * agg$within_var[ok]) / df,
                  .Machine$double.eps)
    residual_variance <- rep(pooled, E)
  }
  residual_variance <- .normalise_environment_variance(
    residual_variance, environments, "residual_variance")

  by_genotype <- split(agg, agg$genotype)
  by_genotype <- by_genotype[vapply(by_genotype, nrow, integer(1)) >= 1L]
  n_obs <- nrow(agg)
  if (n_obs <= E)
    stop("Historical data do not contain residual degrees of freedom for REML.")

  wide <- matrix(NA_real_, length(by_genotype), E,
                 dimnames = list(names(by_genotype), environments))
  for (g in seq_along(by_genotype))
    wide[g, by_genotype[[g]]$env_index] <- by_genotype[[g]]$value
  S0 <- stats::cov(wide, use = "pairwise.complete.obs")
  if (any(!is.finite(S0))) {
    vv <- vapply(seq_len(E), function(e) stats::var(wide[, e], na.rm = TRUE),
                 numeric(1))
    vv[!is.finite(vv) | vv <= 0] <- stats::var(agg$value)
    S0 <- diag(pmax(vv, 1e-4), E)
  }
  S0 <- .bend_pd((S0 + t(S0)) / 2, eps = 1e-5)

  if (model == "diagonal") {
    theta0 <- log(sqrt(diag(S0)))
    unpack <- function(theta) diag(exp(2 * theta), E)
    n_covpar <- E
  } else if (model == "unstructured") {
    L0 <- t(chol(S0)); pos <- which(lower.tri(L0, diag = TRUE), arr.ind = TRUE)
    theta0 <- vapply(seq_len(nrow(pos)), function(i) {
      v <- L0[pos[i, 1], pos[i, 2]]
      if (pos[i, 1] == pos[i, 2]) log(v) else v
    }, numeric(1))
    unpack <- function(theta) {
      L <- matrix(0, E, E)
      for (i in seq_len(nrow(pos)))
        L[pos[i, 1], pos[i, 2]] <- if (pos[i, 1] == pos[i, 2])
          exp(theta[i]) else theta[i]
      tcrossprod(L)
    }
    n_covpar <- length(theta0)
  } else {
    ee <- eigen(S0, symmetric = TRUE)
    Lambda0 <- ee$vectors[, seq_len(rank), drop = FALSE] %*%
      diag(sqrt(pmax(ee$values[seq_len(rank)] * 0.8, 1e-4)), rank)
    load_pos <- do.call(rbind, lapply(seq_len(rank), function(k)
      cbind(row = k:E, col = k)))
    load_theta <- vapply(seq_len(nrow(load_pos)), function(i) {
      v <- Lambda0[load_pos[i, "row"], load_pos[i, "col"]]
      if (load_pos[i, "row"] == load_pos[i, "col"])
        log(max(abs(v), 1e-3)) else v
    }, numeric(1))
    psi0 <- pmax(diag(S0 - tcrossprod(Lambda0)), diag(S0) * 0.05, 1e-5)
    theta0 <- c(load_theta, log(psi0))
    unpack_parts <- function(theta) {
      Lambda <- matrix(0, E, rank)
      for (i in seq_len(nrow(load_pos)))
        Lambda[load_pos[i, "row"], load_pos[i, "col"]] <-
          if (load_pos[i, "row"] == load_pos[i, "col"])
            exp(theta[i]) else theta[i]
      psi <- exp(theta[nrow(load_pos) + seq_len(E)])
      list(Sigma = tcrossprod(Lambda) + diag(psi),
           loadings = Lambda, specific = psi)
    }
    unpack <- function(theta) unpack_parts(theta)$Sigma
    n_covpar <- length(theta0)
  }

  reml_deviance <- function(theta) {
    Sigma <- unpack(theta)
    XtViX <- matrix(0, E, E); XtViy <- numeric(E)
    yViy <- 0; logdetV <- 0
    for (dg in by_genotype) {
      idx <- dg$env_index
      V <- Sigma[idx, idx, drop = FALSE] +
        diag(residual_variance[idx] / dg$n, length(idx))
      ch <- tryCatch(chol(V), error = function(e) NULL)
      if (is.null(ch)) return(1e100)
      Vi <- chol2inv(ch)
      y <- dg$value
      XtViX[idx, idx] <- XtViX[idx, idx, drop = FALSE] + Vi
      XtViy[idx] <- XtViy[idx] + as.numeric(Vi %*% y)
      yViy <- yViy + as.numeric(crossprod(y, Vi %*% y))
      logdetV <- logdetV + 2 * sum(log(diag(ch)))
    }
    chx <- tryCatch(chol(XtViX), error = function(e) NULL)
    if (is.null(chx)) return(1e100)
    beta <- chol2inv(chx) %*% XtViy
    quad <- yViy - as.numeric(crossprod(XtViy, beta))
    val <- logdetV + 2 * sum(log(diag(chx))) + quad +
      (n_obs - E) * log(2 * pi)
    if (!is.finite(val)) 1e100 else val
  }

  fit <- stats::nlminb(theta0, reml_deviance, control = control)
  Sigma_hat <- unpack(fit$par)
  dimnames(Sigma_hat) <- list(environments, environments)
  sd_hat <- sqrt(diag(Sigma_hat))
  Corr <- Sigma_hat / outer(sd_hat, sd_hat)
  parts <- if (model == "fa") unpack_parts(fit$par) else NULL
  if (!is.null(parts)) rownames(parts$loadings) <- environments
  structure(list(
    Sigma_E = Sigma_hat, correlation = Corr,
    residual_variance = residual_variance,
    model = model, rank = if (model == "fa") rank else NA_integer_,
    loadings = if (is.null(parts)) NULL else parts$loadings,
    specific_variances = if (is.null(parts)) NULL else
      stats::setNames(parts$specific, environments),
    logLik = -0.5 * fit$objective,
    AIC = fit$objective + 2 * n_covpar,
    BIC = fit$objective + log(n_obs - E) * n_covpar,
    convergence = fit$convergence, message = fit$message,
    iterations = fit$iterations, evaluations = fit$evaluations,
    n_observations = n_obs, n_genotypes = length(by_genotype),
    environments = environments, aggregated_data = agg,
    call = match.call()), class = "historical_met_fit")
}

#' @export
print.historical_met_fit <- function(x, ...) {
  cat("Historical MET covariance fit\n")
  cat("  Model:", x$model,
      if (identical(x$model, "fa")) paste0("(", x$rank, ")") else "", "\n")
  cat("  Environments:", length(x$environments),
      " Genotypes:", x$n_genotypes, " Cells:", x$n_observations, "\n")
  cat("  REML logLik:", format(x$logLik, digits = 6),
      " Convergence:", x$convergence, "\n")
  invisible(x)
}
