# Matrix-free PCG implementation used by met_information() for large METs.
.pcg_solve <- function(Afun, b, diagonal, tol, maxit) {
  x <- numeric(length(b))
  r <- b
  z <- r / pmax(diagonal, .Machine$double.eps)
  p <- z
  rz <- sum(r * z)
  target <- tol * max(sqrt(sum(b * b)), 1)
  if (sqrt(sum(r * r)) <= target)
    return(list(x = x, converged = TRUE, iterations = 0L,
                residual = sqrt(sum(r * r))))
  converged <- FALSE; it_done <- 0L
  for (it in seq_len(maxit)) {
    Ap <- Afun(p)
    denom <- sum(p * Ap)
    if (!is.finite(denom) || denom <= 0) break
    alpha <- rz / denom
    x <- x + alpha * p
    r <- r - alpha * Ap
    rn <- sqrt(sum(r * r)); it_done <- it
    if (rn <= target) { converged <- TRUE; break }
    z <- r / pmax(diagonal, .Machine$double.eps)
    rz_new <- sum(r * z)
    if (!is.finite(rz_new) || rz == 0) break
    p <- z + (rz_new / rz) * p
    rz <- rz_new
  }
  list(x = x, converged = converged, iterations = it_done,
       residual = sqrt(sum(r * r)))
}

.met_information_pcg <- function(
    allocation_matrix, Gsub, Sigma_E, Ginv, SigEinv,
    sigma_g2, sigma_e2, reps, env_efficiency, target, tpe_weights,
    local_information, trace_samples, solver_tol, solver_maxit) {
  lines <- rownames(allocation_matrix); envs <- colnames(allocation_matrix)
  J <- nrow(allocation_matrix); E <- ncol(allocation_matrix)
  present <- allocation_matrix > 0
  blocks <- vector("list", E); names(blocks) <- envs
  data_diag <- matrix(0, J, E)

  if (is.null(local_information)) {
    dmat <- sweep(reps, 2L, env_efficiency / sigma_e2, `*`) * present
    for (e in seq_len(E)) {
      de <- dmat[, e]; sd <- sum(de)
      blocks[[e]] <- list(d = de, sum = sd)
      if (sd > 0) data_diag[, e] <- de - de^2 / sd
    }
    data_apply <- function(x) {
      X <- matrix(x, J, E); Y <- matrix(0, J, E)
      for (e in seq_len(E)) {
        de <- blocks[[e]]$d; sd <- blocks[[e]]$sum
        if (sd > 0) Y[, e] <- de * X[, e] - de * sum(de * X[, e]) / sd
      }
      as.numeric(Y)
    }
  } else {
    if (!is.list(local_information) || is.null(names(local_information)) ||
        any(names(local_information) == "") || anyDuplicated(names(local_information)) ||
        !all(envs %in% names(local_information)))
      stop("`local_information` must be a uniquely named list covering every environment.")
    for (e in seq_len(E)) {
      L0 <- as.matrix(local_information[[envs[e]]])
      if (!is.numeric(L0) || nrow(L0) != ncol(L0) || any(!is.finite(L0)) ||
          is.null(rownames(L0)) || is.null(colnames(L0)) ||
          anyDuplicated(rownames(L0)) || anyDuplicated(colnames(L0)) ||
          !setequal(rownames(L0), colnames(L0)) ||
          any(!rownames(L0) %in% lines) ||
          !isTRUE(all.equal(L0, t(L0), tolerance = 1e-8)))
        stop("Every local-information matrix must be finite, symmetric, named, and square.")
      L0 <- L0[rownames(L0), rownames(L0), drop = FALSE]
      ee <- eigen((L0 + t(L0)) / 2, symmetric = TRUE,
                  only.values = TRUE)$values
      if (min(ee) < -1e-8 * max(abs(ee), 1))
        stop("A local-information matrix is not positive semidefinite.")
      absent <- intersect(rownames(L0), lines[!present[, e]])
      if (length(absent) && max(abs(L0[absent, , drop = FALSE])) > 1e-8)
        stop("A local-information matrix contains information for absent lines.")
      Le <- matrix(0, J, J, dimnames = list(lines, lines))
      ids <- rownames(L0); Le[ids, ids] <- (L0 + t(L0)) / 2
      blocks[[e]] <- Le; data_diag[, e] <- diag(Le)
    }
    data_apply <- function(x) {
      X <- matrix(x, J, E); Y <- matrix(0, J, E)
      for (e in seq_len(E)) Y[, e] <- blocks[[e]] %*% X[, e]
      as.numeric(Y)
    }
  }

  prior_apply <- function(x) {
    X <- matrix(x, J, E)
    as.numeric((Ginv %*% X %*% t(SigEinv)) / sigma_g2)
  }
  Afun <- function(x) data_apply(x) + prior_apply(x)
  prior_diag <- outer(diag(Ginv), diag(SigEinv)) / sigma_g2
  precond_diag <- as.numeric(data_diag + prior_diag)

  probes <- .with_local_seed(1L,
    matrix(sample(c(-1, 1), J * trace_samples, replace = TRUE),
           J, trace_samples))
  pev_diag <- numeric(J)
  converged <- logical(trace_samples)
  iterations <- integer(trace_samples)
  residuals <- numeric(trace_samples)
  for (k in seq_len(trace_samples)) {
    z <- probes[, k]
    b <- as.numeric(outer(z, tpe_weights))
    sol <- .pcg_solve(Afun, b, precond_diag, solver_tol, solver_maxit)
    hz <- as.numeric(matrix(sol$x, J, E) %*% tpe_weights)
    pev_diag <- pev_diag + z * hz
    converged[k] <- sol$converged
    iterations[k] <- sol$iterations
    residuals[k] <- sol$residual
  }
  PEV_line <- pmax(pev_diag / trace_samples, 0)
  names(PEV_line) <- lines
  target_SigE <- as.numeric(crossprod(tpe_weights,
                                      Sigma_E %*% tpe_weights))
  priorvar_line <- sigma_g2 * diag(Gsub) * target_SigE
  CD_line <- pmin(1, pmax(0, 1 - PEV_line / priorvar_line))

  out_mean_PEV <- mean(PEV_line)
  out_CD <- mean(CD_line)
  if (target == "environment_specific") {
    idx_present <- which(as.numeric(present) == 1L)
    cell_trace <- 0
    cell_probes <- .with_local_seed(2L,
      matrix(sample(c(-1, 1), length(idx_present) * trace_samples,
                    replace = TRUE), length(idx_present), trace_samples))
    for (k in seq_len(trace_samples)) {
      z <- cell_probes[, k]
      b <- numeric(J * E); b[idx_present] <- z
      sol <- .pcg_solve(Afun, b, precond_diag, solver_tol, solver_maxit)
      cell_trace <- cell_trace + sum(z * sol$x[idx_present])
    }
    out_mean_PEV <- cell_trace / (trace_samples * length(idx_present))
    out_CD <- NA_real_
  }

  list(
    mean_PEV = out_mean_PEV, CDmean = out_CD,
    PEV_per_line = PEV_line, CD_per_line = CD_line,
    C_uu = NULL, data_information = NULL,
    J = J, E = E, target = target, tpe_weights = tpe_weights,
    sigma_e2 = sigma_e2,
    local_information_used = !is.null(local_information),
    solver_used = "pcg", approximation = TRUE,
    solver_diagnostics = list(
      trace_samples = trace_samples,
      converged_fraction = mean(converged),
      mean_iterations = mean(iterations), max_iterations = max(iterations),
      max_residual = max(residuals), tolerance = solver_tol))
}
