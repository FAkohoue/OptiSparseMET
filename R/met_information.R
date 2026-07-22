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
#' \eqn{d_{je} = \text{reps}_{je}\, \epsilon_e / \sigma_e^2}, where
#' \eqn{\epsilon_e \in (0,1]} is the within-environment efficiency factor (from a
#' local design's A / CDmean evaluation). Environment means are absorbed as fixed
#' effects, giving the effect information matrix
#' \deqn{C_{uu} = D - D X (X^\top D X)^{-1} X^\top D + \sigma_g^{-2}
#'   (\Sigma_E^{-1} \otimes G^{-1}).}
#' (The absorption form is verified to match the full mixed-model-equations
#' solution.) The across-TPE breeding value of line \eqn{j} is
#' \eqn{\bar u_j = E^{-1}\sum_e u_{je}}; its PEV is
#' \eqn{a_j^\top C_{uu}^{-1} a_j} with \eqn{a_j} averaging line \eqn{j}'s cells,
#' and its reliability is \eqn{\text{CD}_j = 1 - \text{PEV}_j / \operatorname{Var}(\bar u_j)}.
#'
#' @param allocation_matrix A \eqn{J \times E} 0/1 incidence matrix (lines as
#'   rows, environments as columns) with row and column names. Typically
#'   `alloc$allocation_matrix` from [allocate_sparse_met()].
#' @param G Genomic (or pedigree) relationship matrix with row/column names
#'   covering the lines in `allocation_matrix`.
#' @param Sigma_E Optional \eqn{E \times E} environment covariance (genetic
#'   correlations across environments). Defaults to the identity (independent
#'   environments, unit variance).
#' @param sigma_g2,sigma_e2 Genetic and residual variance scalars.
#' @param reps Optional \eqn{J \times E} matrix of plot counts per cell.
#'   Defaults to `allocation_matrix` (one plot per present cell).
#' @param env_efficiency Optional length-\eqn{E} vector of within-environment
#'   design efficiency factors \eqn{\epsilon_e \in (0,1]} (the coupling to the
#'   local design). Defaults to 1 for every environment.
#' @param target Either `"across_tpe"` (reliability of the TPE-average breeding
#'   value, the default) or `"environment_specific"` (per-cell PEV summary).
#' @param max_dim Integer guard on \eqn{J \times E}; above it the dense inverse
#'   is refused with an informative error.
#'
#' @return A list with `mean_PEV`, `CDmean`, `PEV_per_line` (across-TPE), the
#'   assembled `C_uu`, and metadata (`J`, `E`, `target`).
#'
#' @seealso [allocate_sparse_met()], [simulate_met()].
#' @export
met_information <- function(allocation_matrix, G, Sigma_E = NULL,
                            sigma_g2 = 1, sigma_e2 = 1,
                            reps = NULL, env_efficiency = NULL,
                            target = c("across_tpe", "environment_specific"),
                            max_dim = 6000L) {
  target <- match.arg(target)
  if (is.null(rownames(allocation_matrix)) || is.null(colnames(allocation_matrix)))
    stop("`allocation_matrix` must have row (line) and column (environment) names.")

  lines <- rownames(allocation_matrix)
  envs  <- colnames(allocation_matrix)
  J <- length(lines); E <- length(envs)
  if (J < 1L || E < 1L) stop("allocation_matrix must be non-empty.")
  if (J * E > max_dim)
    stop("J*E = ", J * E, " exceeds max_dim (", max_dim, "). Increase max_dim ",
         "only if you can afford a dense (JE x JE) inverse.")

  if (is.null(rownames(G)) || !all(lines %in% rownames(G)))
    stop("`G` must have row/column names covering all lines in allocation_matrix.")
  Gsub <- as.matrix(G[lines, lines, drop = FALSE])

  if (is.null(Sigma_E)) Sigma_E <- diag(E)
  Sigma_E <- as.matrix(Sigma_E)
  if (nrow(Sigma_E) != E || ncol(Sigma_E) != E)
    stop("`Sigma_E` must be E x E.")

  if (is.null(reps)) reps <- allocation_matrix
  reps <- as.matrix(reps)
  if (!all(dim(reps) == c(J, E))) stop("`reps` must be J x E.")

  if (is.null(env_efficiency)) env_efficiency <- rep(1, E)
  if (length(env_efficiency) != E) stop("`env_efficiency` must have length E.")

  # Inverses (Cholesky with pseudo-inverse fallback).
  safe_inv <- function(M) {
    out <- try(solve(M), silent = TRUE)
    if (inherits(out, "try-error")) .pinv_sym_dense(M) else out
  }
  Ginv    <- safe_inv(Gsub)
  SigEinv <- safe_inv(Sigma_E)

  # Prior precision, u indexed as (e - 1) * J + j  => kronecker(SigEinv, Ginv).
  prior_prec <- (1 / sigma_g2) * kronecker(SigEinv, Ginv)

  # Cell information d_{je} and env-mean absorption.
  d <- as.numeric(sweep(reps, 2, env_efficiency / sigma_e2, `*`))  # column-major = (e-1)*J+j
  present <- as.numeric(allocation_matrix) > 0
  # X'D X is diagonal E x E: sum of d over present cells per environment.
  env_of <- rep(seq_len(E), each = J)
  d_present <- d * present
  XtDX_diag <- tapply(d_present, env_of, sum)
  XtDX_diag <- as.numeric(XtDX_diag[as.character(seq_len(E))])
  XtDX_diag[is.na(XtDX_diag)] <- 0

  # Absorption term A = D X (X'DX)^{-1} X'D. With X a present-cell env indicator
  # and D diagonal, A is block-diagonal by environment: within env e,
  # A[c, c'] = d_c d_c' / XtDX_diag[e]  (present cells only).
  Cuu <- diag(d) + prior_prec
  for (e in seq_len(E)) {
    idx <- which(env_of == e & present == 1)
    if (length(idx) == 0L || XtDX_diag[e] <= 0) next
    Cuu[idx, idx] <- Cuu[idx, idx] - outer(d[idx], d[idx]) / XtDX_diag[e]
  }

  Cuu_inv <- safe_inv(Cuu)

  # Across-TPE: BV_j = mean_e u_{je}; positions for line j are j + (0:(E-1))*J.
  PEV_line <- numeric(J)
  priorvar_line <- numeric(J)
  mean_SigE <- mean(Sigma_E)
  for (j in seq_len(J)) {
    idx <- j + (seq_len(E) - 1L) * J
    PEV_line[j]      <- sum(Cuu_inv[idx, idx]) / (E^2)
    priorvar_line[j] <- sigma_g2 * Gsub[j, j] * mean_SigE
  }
  names(PEV_line) <- lines
  CD_line <- 1 - PEV_line / priorvar_line

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
    J            = J,
    E            = E,
    target       = target
  )
}
