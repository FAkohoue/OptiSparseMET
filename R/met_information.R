#' Across-environment (MET-level) information matrix and prediction reliability
#'
#' @description
#' `met_information()` couples the two design levels that `OptiSparseMET`
#' targets. It builds the joint mixed-model information matrix for a
#' genotype-by-environment (G x E) model over the whole multi-environment trial,
#' combining **which lines are tested in which environments** (the allocation
#' incidence, Level 1) with **how precisely each environment measures a plot**
#' (a within-environment design-efficiency factor, Level 2). From it, it returns
#' the prediction-error variance (PEV) and reliability (CDmean) of genotype
#' performance across the target population of environments (TPE).
#'
#' This is the quantity the allocation exists to optimise, so it is the natural
#' objective for [allocate_sparse_met()] (see `objective = "optcontrib_GxE"`) and
#' the metric that demonstrates the two levels are optimised consistently rather
#' than in isolation.
#'
#' @details
#' The G x E effects \eqn{u} (one per line-by-environment cell) have prior
#' covariance \eqn{\sigma_g^2\,(\Sigma_E \otimes G)}, so prior precision
#' \eqn{\sigma_g^{-2}(\Sigma_E^{-1} \otimes G^{-1})}. Each present cell
#' \eqn{(j,e)} contributes diagonal information
#' \eqn{d_{je} = \text{reps}_{je}\, \epsilon_e / \sigma_{\mathrm{res},e}^2}, where
#' \eqn{\epsilon_e \in (0,1]} is the within-environment efficiency factor (from a
#' local design's A / CDmean evaluation). For final evaluation, a full local
#' treatment-information matrix may replace this scalar approximation and retain
#' off-diagonal contrast information. Environment means are absorbed as fixed
#' effects, giving the effect information matrix
#' \deqn{C_{uu} = D - D X (X^\top D X)^{-1} X^\top D + \sigma_g^{-2}
#'   (\Sigma_E^{-1} \otimes G^{-1}).}
#' (The absorption form is verified to match the full mixed-model-equations
#' solution.) The across-TPE breeding value of line \eqn{j} is
#' \eqn{\bar u_j = \sum_e w_e u_{je}}, where non-negative target weights sum to
#' one. Its PEV is \eqn{a_j^\top C_{uu}^{-1} a_j}, with \eqn{a_j} placing the
#' weights on line \eqn{j}'s environment cells, and its reliability is
#' \eqn{\text{CD}_j = 1 - \text{PEV}_j / \operatorname{Var}(\bar u_j)}.
#'
#' @param allocation_matrix A \eqn{J \times E} 0/1 incidence matrix (lines as
#'   rows, environments as columns) with row and column names. Typically
#'   `alloc$allocation_matrix` from [allocate_sparse_met()].
#' @param G Genomic (or pedigree) relationship matrix with row/column names
#'   covering the lines in `allocation_matrix`.
#' @param Sigma_E Optional \eqn{E \times E} environment covariance (genetic
#'   correlations across environments). Defaults to the identity (independent
#'   environments, unit variance).
#' @param sigma_g2 Genetic variance scalar.
#' @param sigma_e2 Residual variance scalar or a vector with one value per
#'   environment. A named vector is aligned to `allocation_matrix` columns.
#' @param reps Optional \eqn{J \times E} matrix of plot counts per cell.
#'   Defaults to `allocation_matrix` (one plot per present cell).
#' @param env_efficiency Optional length-\eqn{E} vector of within-environment
#'   design efficiency factors \eqn{\epsilon_e \in (0,1]} (the coupling to the
#'   local design). Defaults to 1 for every environment.
#' @param target Either `"across_tpe"` (reliability of the TPE-average breeding
#'   value, the default) or `"environment_specific"` (per-cell PEV summary).
#' @param max_dim Integer threshold on \eqn{J \times E}. With
#'   `solver = "auto"`, larger problems use the iterative matrix-free path.
#' @param solver `"auto"`, `"dense"`, or `"pcg"`. The PCG path avoids the
#'   dense \eqn{JE \times JE} matrix and estimates target PEV diagonals with
#'   deterministic Hutchinson probes.
#' @param trace_samples Number of probe vectors for `solver = "pcg"`.
#' @param solver_tol,solver_maxit Relative tolerance and iteration limit for
#'   each preconditioned conjugate-gradient solve.
#' @param tpe_weights Optional non-negative weights defining the target
#'   population of environments. Values are normalised to sum to one. The
#'   default gives every environment equal weight.
#' @param local_information Optional named list of treatment-information
#'   matrices, one per environment. Each matrix may cover all lines or a named
#'   subset of the lines present in that environment and must already have
#'   nuisance fixed and random effects absorbed. When supplied, these matrices
#'   replace the scalar `reps * env_efficiency / sigma_e2` approximation and
#'   therefore preserve the local layout's full off-diagonal information.
#'
#' @return A list with `mean_PEV`, `CDmean`, per-line PEV and CD, aligned TPE
#'   weights and residual variances, dimensions, target, solver identity,
#'   approximation flag, and solver diagnostics. The exact dense path also
#'   returns the assembled `C_uu` and `data_information` matrices; these are
#'   `NULL` on the matrix-free PCG path.
#'
#' @seealso [allocate_sparse_met()], [simulate_met()].
#' @export
met_information <- function(allocation_matrix, G, Sigma_E = NULL,
                            sigma_g2 = 1, sigma_e2 = 1,
                            reps = NULL, env_efficiency = NULL,
                            target = c("across_tpe", "environment_specific"),
                            max_dim = 6000L, tpe_weights = NULL,
                            local_information = NULL,
                            solver = c("auto", "dense", "pcg"),
                            trace_samples = 40L, solver_tol = 1e-7,
                            solver_maxit = 1000L) {
  target <- match.arg(target)
  solver <- match.arg(solver)
  allocation_matrix <- as.matrix(allocation_matrix)
  if (!is.numeric(allocation_matrix) || any(!is.finite(allocation_matrix)) ||
      any(allocation_matrix < 0))
    stop("`allocation_matrix` must be a finite, non-negative numeric matrix.")
  if (is.null(rownames(allocation_matrix)) || is.null(colnames(allocation_matrix)))
    stop("`allocation_matrix` must have row (line) and column (environment) names.")
  if (anyDuplicated(rownames(allocation_matrix)) ||
      anyDuplicated(colnames(allocation_matrix)))
    stop("`allocation_matrix` row and column names must be unique.")
  if (!is.numeric(sigma_g2) || length(sigma_g2) != 1L ||
      !is.finite(sigma_g2) || sigma_g2 <= 0)
    stop("`sigma_g2` must be a finite positive scalar.")
  if (!is.numeric(max_dim) || length(max_dim) != 1L ||
      !is.finite(max_dim) || max_dim < 1)
    stop("`max_dim` must be a positive scalar.")

  lines <- rownames(allocation_matrix)
  envs  <- colnames(allocation_matrix)
  J <- length(lines); E <- length(envs)
  if (J < 1L || E < 1L) stop("allocation_matrix must be non-empty.")
  sigma_e2 <- .normalise_environment_variance(sigma_e2, envs)
  tpe_weights <- .normalise_tpe_weights(tpe_weights, envs)
  solver_used <- if (solver == "auto")
    if (J * E <= max_dim) "dense" else "pcg" else solver
  if (solver_used == "dense" && J * E > max_dim)
    stop("J*E = ", J * E, " exceeds max_dim (", max_dim, "). Increase max_dim ",
         "or use `solver = \"pcg\"`.")
  if (!is.numeric(trace_samples) || length(trace_samples) != 1L ||
      !is.finite(trace_samples) || trace_samples < 2L ||
      trace_samples != as.integer(trace_samples))
    stop("`trace_samples` must be an integer >= 2.")
  if (!is.numeric(solver_tol) || length(solver_tol) != 1L ||
      !is.finite(solver_tol) || solver_tol <= 0 ||
      !is.numeric(solver_maxit) || length(solver_maxit) != 1L ||
      !is.finite(solver_maxit) || solver_maxit < 1L ||
      solver_maxit != as.integer(solver_maxit))
    stop("Invalid `solver_tol` or `solver_maxit`.")

  G <- as.matrix(G)
  if (!is.numeric(G) || nrow(G) != ncol(G) || any(!is.finite(G)))
    stop("`G` must be a finite numeric square matrix.")
  if (is.null(rownames(G)) || is.null(colnames(G)) ||
      anyDuplicated(rownames(G)) || anyDuplicated(colnames(G)) ||
      !all(lines %in% rownames(G)) || !all(lines %in% colnames(G)))
    stop("`G` must have unique row/column names covering all allocation lines.")
  Gsub <- as.matrix(G[lines, lines, drop = FALSE])
  if (!isTRUE(all.equal(Gsub, t(Gsub), tolerance = 1e-8)))
    stop("`G` must be symmetric after alignment to allocation lines.")
  if (any(diag(Gsub) <= 0))
    stop("`G` must have strictly positive diagonal elements.")

  if (is.null(Sigma_E)) {
    Sigma_E <- diag(E)
    dimnames(Sigma_E) <- list(envs, envs)
  }
  Sigma_E <- as.matrix(Sigma_E)
  if (nrow(Sigma_E) != E || ncol(Sigma_E) != E)
    stop("`Sigma_E` must be E x E.")
  if (!is.numeric(Sigma_E) || any(!is.finite(Sigma_E)))
    stop("`Sigma_E` must be a finite numeric matrix.")
  if (!is.null(rownames(Sigma_E)) || !is.null(colnames(Sigma_E))) {
    if (is.null(rownames(Sigma_E)) || is.null(colnames(Sigma_E)) ||
        anyDuplicated(rownames(Sigma_E)) || anyDuplicated(colnames(Sigma_E)) ||
        !all(envs %in% rownames(Sigma_E)) ||
        !all(envs %in% colnames(Sigma_E)))
      stop("Named `Sigma_E` must uniquely cover all allocation environments.")
    Sigma_E <- Sigma_E[envs, envs, drop = FALSE]
  } else {
    dimnames(Sigma_E) <- list(envs, envs)
  }
  if (!isTRUE(all.equal(Sigma_E, t(Sigma_E), tolerance = 1e-8)))
    stop("`Sigma_E` must be symmetric.")

  if (is.null(reps)) {
    reps <- allocation_matrix
  } else {
    reps <- as.matrix(reps)
    if (!is.null(rownames(reps)) || !is.null(colnames(reps))) {
      if (is.null(rownames(reps)) || is.null(colnames(reps)) ||
          !all(lines %in% rownames(reps)) || !all(envs %in% colnames(reps)))
        stop("Named `reps` must cover all allocation lines and environments.")
      reps <- reps[lines, envs, drop = FALSE]
    }
  }
  reps <- as.matrix(reps)
  if (!all(dim(reps) == c(J, E))) stop("`reps` must be J x E.")
  if (!is.numeric(reps) || any(!is.finite(reps)) || any(reps < 0))
    stop("`reps` must contain finite, non-negative plot counts.")
  present_matrix <- allocation_matrix > 0
  if (any(reps[present_matrix] <= 0) || any(reps[!present_matrix] != 0))
    stop("`reps` must be positive exactly where `allocation_matrix` is present.")

  if (is.null(env_efficiency)) {
    env_efficiency <- stats::setNames(rep(1, E), envs)
  } else if (!is.null(names(env_efficiency))) {
    if (!all(envs %in% names(env_efficiency)))
      stop("Named `env_efficiency` must cover all environments.")
    env_efficiency <- env_efficiency[envs]
  }
  if (length(env_efficiency) != E) stop("`env_efficiency` must have length E.")
  if (!is.numeric(env_efficiency) || any(!is.finite(env_efficiency)) ||
      any(env_efficiency <= 0 | env_efficiency > 1))
    stop("`env_efficiency` values must be finite and in (0, 1].")

  # Covariance precision requires a positive-definite matrix. A Moore-Penrose
  # inverse is not appropriate here: a zero-variance direction would receive
  # zero rather than infinite precision. Reject materially indefinite inputs
  # and bend singular/near-singular PSD matrices by a negligible eigenvalue
  # floor before inversion.
  covariance_precision <- function(M, nm) {
    M <- (M + t(M)) / 2
    ee <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
    scale <- max(abs(ee), 1)
    if (min(ee) < -1e-8 * scale)
      stop("`", nm, "` must be positive semidefinite.")
    solve(.bend_pd(M, eps = 1e-8))
  }
  Ginv    <- covariance_precision(Gsub, "G")
  SigEinv <- covariance_precision(Sigma_E, "Sigma_E")

  if (solver_used == "pcg")
    return(.met_information_pcg(
      allocation_matrix = allocation_matrix, Gsub = Gsub,
      Sigma_E = Sigma_E, Ginv = Ginv, SigEinv = SigEinv,
      sigma_g2 = sigma_g2, sigma_e2 = sigma_e2, reps = reps,
      env_efficiency = env_efficiency, target = target,
      tpe_weights = tpe_weights, local_information = local_information,
      trace_samples = as.integer(trace_samples), solver_tol = solver_tol,
      solver_maxit = as.integer(solver_maxit)))

  # Prior precision, u indexed as (e - 1) * J + j  => kronecker(SigEinv, Ginv).
  prior_prec <- (1 / sigma_g2) * kronecker(SigEinv, Ginv)

  # Cell information d_{je} and env-mean absorption.
  d <- as.numeric(sweep(reps, 2, env_efficiency / sigma_e2, `*`))  # column-major = (e-1)*J+j
  present <- as.numeric(present_matrix)
  # X'D X is diagonal E x E: sum of d over present cells per environment.
  env_of <- rep(seq_len(E), each = J)
  d_present <- d * present
  XtDX_diag <- tapply(d_present, env_of, sum)
  XtDX_diag <- as.numeric(XtDX_diag[as.character(seq_len(E))])
  XtDX_diag[is.na(XtDX_diag)] <- 0

  # Absorption term A = D X (X'DX)^{-1} X'D. With X a present-cell env indicator
  # and D diagonal, A is block-diagonal by environment: within env e,
  # A[c, c'] = d_c d_c' / XtDX_diag[e]  (present cells only).
  data_information <- matrix(0, J * E, J * E)
  if (is.null(local_information)) {
    diag(data_information) <- d_present
    for (e in seq_len(E)) {
      idx <- which(env_of == e & present == 1)
      if (length(idx) == 0L || XtDX_diag[e] <= 0) next
      data_information[idx, idx] <- data_information[idx, idx] -
        outer(d[idx], d[idx]) / XtDX_diag[e]
    }
  } else {
    if (!is.list(local_information) || is.null(names(local_information)) ||
        any(names(local_information) == "") ||
        anyDuplicated(names(local_information)) ||
        !all(envs %in% names(local_information)))
      stop("`local_information` must be a uniquely named list covering every environment.")
    for (e in seq_len(E)) {
      Le0 <- as.matrix(local_information[[envs[e]]])
      if (!is.numeric(Le0) || nrow(Le0) != ncol(Le0) ||
          any(!is.finite(Le0)) || is.null(rownames(Le0)) ||
          is.null(colnames(Le0)) || anyDuplicated(rownames(Le0)) ||
          anyDuplicated(colnames(Le0)) ||
          !setequal(rownames(Le0), colnames(Le0)) ||
          any(!rownames(Le0) %in% lines))
        stop("Every local-information matrix must be a finite, named square matrix over allocation lines.")
      Le0 <- Le0[rownames(Le0), rownames(Le0), drop = FALSE]
      if (!isTRUE(all.equal(Le0, t(Le0), tolerance = 1e-8)))
        stop("Local-information matrix for `", envs[e], "` is not symmetric.")
      ee <- eigen((Le0 + t(Le0)) / 2, symmetric = TRUE,
                  only.values = TRUE)$values
      if (min(ee) < -1e-8 * max(abs(ee), 1))
        stop("Local-information matrix for `", envs[e],
             "` is not positive semidefinite.")
      absent <- intersect(rownames(Le0), lines[!present_matrix[, e]])
      if (length(absent) &&
          max(abs(Le0[absent, , drop = FALSE])) > 1e-8)
        stop("Local-information matrix for `", envs[e],
             "` contains information for absent lines.")
      Le <- matrix(0, J, J, dimnames = list(lines, lines))
      ids <- rownames(Le0)
      Le[ids, ids] <- (Le0 + t(Le0)) / 2
      idx <- (e - 1L) * J + seq_len(J)
      data_information[idx, idx] <- Le
    }
  }
  Cuu <- data_information + prior_prec

  Cuu_inv <- tryCatch(solve(Cuu), error = function(e) .pinv_sym_dense(Cuu))

  # Across-TPE: BV_j = mean_e u_{je}; positions for line j are j + (0:(E-1))*J.
  PEV_line <- numeric(J)
  priorvar_line <- numeric(J)
  target_SigE <- as.numeric(crossprod(
    tpe_weights, Sigma_E %*% tpe_weights
  ))
  if (!is.finite(target_SigE) || target_SigE <= 0)
    stop("The across-TPE variance implied by `Sigma_E` must be positive.")
  for (j in seq_len(J)) {
    idx <- j + (seq_len(E) - 1L) * J
    PEV_line[j] <- as.numeric(crossprod(
      tpe_weights, Cuu_inv[idx, idx, drop = FALSE] %*% tpe_weights
    ))
    priorvar_line[j] <- sigma_g2 * Gsub[j, j] * target_SigE
  }
  names(PEV_line) <- lines
  CD_line <- pmin(1, pmax(0, 1 - PEV_line / priorvar_line))

  if (target == "environment_specific") {
    cell_PEV <- diag(Cuu_inv)
    out_mean_PEV <- mean(cell_PEV[present == 1])
    out_CD       <- NA_real_
  } else {
    out_mean_PEV <- mean(PEV_line)
    out_CD       <- mean(CD_line)
  }

  list(
    mean_PEV     = out_mean_PEV,
    CDmean       = out_CD,
    PEV_per_line = PEV_line,
    CD_per_line  = CD_line,
    C_uu         = Cuu,
    data_information = data_information,
    J            = J,
    E            = E,
    target       = target,
    tpe_weights  = tpe_weights,
    sigma_e2     = sigma_e2,
    local_information_used = !is.null(local_information),
    solver_used = "dense",
    approximation = FALSE,
    solver_diagnostics = NULL
  )
}
