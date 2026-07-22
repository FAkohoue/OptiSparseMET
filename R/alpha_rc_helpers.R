# ==============================================================================
# alpha_rc_helpers.R
# Internal helper functions shared across:
#   - alpha_rc_stream()
#   - evaluate_design_efficiency()
#   - optimize_alpha_rc()
#
# These functions are NOT exported and are NOT intended for direct user calls.
# They are prefixed with "." to signal internal status.
#
# Contents:
#   Matrix helpers
#     .make_sparse_incidence()   - build a sparse 0/1 incidence matrix
#     .pinv_sym_dense()          - Moore-Penrose pseudoinverse (symmetric)
#     .safe_logdet_psd_dense()   - log-determinant of a PSD matrix
#     .pairwise_diff_mean_var()  - mean pairwise contrast variance from C^{-1}
#     .solve_C()                 - unified solver for dense and sparse C
#     .ar1_precision_sparse()    - AR1 precision matrix (tridiagonal sparse)
#     .trace_subinv_est()        - Hutchinson stochastic trace estimator
#
#   RNG sandbox
#     .with_local_seed()         - evaluate an expression under an isolated seed
#
#   Spatial helpers
#     .build_neighbor_pairs()    - enumerate plot neighbour pairs by Chebyshev
#     .score_dispersion()        - total relatedness of neighbouring plots
# ==============================================================================


# ==============================================================================
# MATRIX HELPERS
# ==============================================================================

# ------------------------------------------------------------------------------
# .make_sparse_incidence
# ------------------------------------------------------------------------------
# Build a sparse 0/1 incidence matrix from a character (or factor) vector.
#
# Each unique non-NA level in `levels_vec` becomes one column. Row i has a 1
# in the column corresponding to levels_vec[i], and 0 elsewhere. NA entries
# produce all-zero rows.
#
# Used to build Z matrices (replicate, block, row, column, treatment incidence)
# in evaluate_design_efficiency().
#
# Arguments:
#   levels_vec  - character or factor vector of group labels (length n).
#                 NA values are silently treated as unassigned (zero row).
#
# Returns:
#   A named list:
#     $M       - sparse Matrix of dimension n x length(unique non-NA levels),
#                column names set to the unique level labels.
#     $levels  - character vector of unique non-NA levels (column order of M).
#
# Example:
#   .make_sparse_incidence(c("A", "B", "A", NA, "C"))
#   # 5 x 3 sparse matrix; row 4 is all-zero
# ------------------------------------------------------------------------------
.make_sparse_incidence <- function(levels_vec) {
  nn <- length(levels_vec)
  lv <- unique(levels_vec[!is.na(levels_vec)])
  if (length(lv) == 0) {
    return(list(
      M      = Matrix::Matrix(0, nrow = nn, ncol = 0, sparse = TRUE),
      levels = character(0)
    ))
  }
  j  <- match(levels_vec, lv)
  ok <- !is.na(j)
  M  <- Matrix::sparseMatrix(
    i    = which(ok),
    j    = j[ok],
    x    = 1,
    dims = c(nn, length(lv))
  )
  colnames(M) <- lv
  list(M = M, levels = lv)
}


# ------------------------------------------------------------------------------
# .pinv_sym_dense
# ------------------------------------------------------------------------------
# Compute the Moore-Penrose pseudoinverse of a real symmetric matrix.
#
# Uses the spectral decomposition A = V diag(lambda) V'. Only eigenvalues
# strictly greater than `tol` are inverted; the remaining (near-zero)
# eigenvalues contribute zero to the pseudoinverse. This is equivalent to
# truncated eigenvalue inversion.
#
# Used as a fallback in evaluate_design_efficiency() when the sparse Cholesky
# factorisation of K fails (near-singular genomic relationship matrices).
#
# Arguments:
#   A    - symmetric numeric matrix (dense or coercible to dense).
#   tol  - eigenvalue threshold below which eigenvalues are treated as zero.
#           Default 1e-10.
#
# Returns:
#   Numeric matrix of the same dimensions as A, the pseudoinverse A^+.
#   Returns the zero matrix of the same dimensions if all eigenvalues <= tol.
# ------------------------------------------------------------------------------
.pinv_sym_dense <- function(A, tol = 1e-10) {
  eg   <- eigen(as.matrix(A), symmetric = TRUE)
  vals <- eg$values
  vecs <- eg$vectors
  keep <- vals > tol
  if (!any(keep)) return(matrix(0, nrow(A), ncol(A)))
  vecs[, keep, drop = FALSE] %*%
    diag(1 / vals[keep], nrow = sum(keep)) %*%
    t(vecs[, keep, drop = FALSE])
}


# ------------------------------------------------------------------------------
# .safe_logdet_psd_dense
# ------------------------------------------------------------------------------
# Compute the log-determinant of a positive semi-definite (PSD) matrix.
#
# Symmetrises the matrix first to guard against floating-point asymmetry,
# then extracts eigenvalues and sums the logs of those that exceed `tol`.
# Eigenvalues at or below `tol` are treated as structural zeros and excluded,
# making this safe for rank-deficient (semidefinite) matrices.
#
# Used in evaluate_design_efficiency() to compute the D-criterion:
#   D_criterion = exp( log_det(H V H) / (p - 1) )
# where H is the centering matrix and V is the treatment variance-covariance
# submatrix extracted from C^{-1}.
#
# Arguments:
#   M    - numeric matrix, assumed symmetric PSD.
#   tol  - eigenvalue threshold for numerical rank determination. Default 1e-10.
#
# Returns:
#   Numeric scalar: sum of log(eigenvalues > tol).
#   Returns NA_real_ if no eigenvalue exceeds tol (rank-0 matrix).
# ------------------------------------------------------------------------------
.safe_logdet_psd_dense <- function(M, tol_rel = 1e-8, expected_rank = NULL) {
  # Log pseudo-determinant of a symmetric PSD matrix: sum of logs of the
  # positive eigenvalues. Uses a RELATIVE tolerance (tol_rel * max eigenvalue)
  # so the result is invariant to the overall scale of M -- an absolute
  # threshold could either keep a numerical-noise "zero" eigenvalue (collapsing
  # the determinant) or drop a genuinely small one (inflating it).
  #
  # When expected_rank is supplied (e.g. p - 1 for a centered contrast
  # covariance HVH), exactly that many top eigenvalues are used, guaranteeing the
  # eigenvalue count matches the divisor used by the caller. If fewer positive
  # eigenvalues survive the tolerance than expected_rank, the matrix is
  # rank-deficient beyond its structural null space and NA is returned.
  ev <- eigen((M + t(M)) / 2, symmetric = TRUE, only.values = TRUE)$values
  ev <- sort(ev, decreasing = TRUE)
  mx <- ev[1]
  if (!is.finite(mx) || mx <= 0) return(NA_real_)
  ev_pos <- ev[ev > tol_rel * mx]
  if (!is.null(expected_rank)) {
    if (length(ev_pos) < expected_rank) return(NA_real_)
    ev_pos <- ev_pos[seq_len(expected_rank)]
  }
  if (length(ev_pos) == 0) return(NA_real_)
  sum(log(ev_pos))
}


# ------------------------------------------------------------------------------
# .pairwise_diff_mean_var
# ------------------------------------------------------------------------------
# Compute the mean pairwise contrast variance from a variance-covariance
# submatrix V = C^{-1}[trt, trt].
#
# For a vector of treatment estimates tau_hat with variance-covariance matrix V,
# the variance of the contrast (tau_i - tau_j) is:
#
#   Var(tau_i - tau_j) = V_ii + V_jj - 2 * V_ij
#
# This function computes the mean of this quantity over all p*(p-1)/2 unique
# pairs (i, j) with i < j. This is the A-criterion in criterion form:
#
#   A_criterion = mean_{i < j} Var(tau_i - tau_j)
#               = (2 / (p*(p-1))) * sum_{i < j} (V_ii + V_jj - 2*V_ij)
#
# Lower values indicate more precise estimation of treatment contrasts.
#
# Used in evaluate_design_efficiency() for the FIXED_TREATMENT_BLUE_CONTRAST
# computation mode.
#
# Arguments:
#   V  - square numeric matrix (p x p), the treatment variance-covariance
#        submatrix extracted from C^{-1}.
#
# Returns:
#   Numeric scalar. NA_real_ if p < 2.
# ------------------------------------------------------------------------------
.pairwise_diff_mean_var <- function(V) {
  p <- nrow(V)
  if (p < 2) return(NA_real_)
  d <- outer(diag(V), diag(V), "+") - 2 * V
  mean(d[upper.tri(d)])
}


# ------------------------------------------------------------------------------
# .solve_C
# ------------------------------------------------------------------------------
# Unified linear system solver for both dense and sparse C matrices.
#
# The mixed-model coefficient matrix C can be either a base R dense matrix
# (from `as.matrix()`) or a sparse Matrix object (from the Matrix package),
# depending on the spatial_engine setting and design size. This wrapper
# dispatches to the appropriate solver so calling code does not need to branch.
#
# For sparse C, uses Matrix::solve() with sparse = FALSE to return a dense
# solution matrix (the right-hand side B is typically sparse but the solution
# is dense). For dense C, uses base R solve().
#
# Used in evaluate_design_efficiency() to extract submatrices of C^{-1}
# corresponding to treatment effects.
#
# Arguments:
#   Cmat  - the mixed model coefficient matrix; either a base R matrix or a
#           Matrix sparse matrix.
#   B     - right-hand side matrix (dense or sparse).
#
# Returns:
#   Dense numeric matrix: the solution X such that Cmat %*% X = B.
# ------------------------------------------------------------------------------
.solve_C <- function(Cmat, B) {
  if (inherits(Cmat, "Matrix")) {
    Matrix::solve(Cmat, B, sparse = FALSE)
  } else {
    solve(Cmat, as.matrix(B))
  }
}


# ------------------------------------------------------------------------------
# .ar1_precision_sparse
# ------------------------------------------------------------------------------
# Construct the precision (inverse covariance) matrix of a stationary AR1
# process as a sparse tridiagonal matrix.
#
# For a stationary AR1 process of length n with autocorrelation parameter rho,
# the marginal variance is sigma^2 / (1 - rho^2). The precision matrix
# Q = Sigma^{-1} / sigma^2 (scaled to unit marginal variance) is tridiagonal:
#
#   Let a = 1 / (1 - rho^2)
#
#   Diagonal entries:
#     Q[1,1]   = Q[n,n]   = a               (edge)
#     Q[i,i]              = (1 + rho^2) * a  (interior, i = 2, ..., n-1)
#
#   Off-diagonal entries (all):
#     Q[i, i+1] = Q[i+1, i] = -rho * a
#
# This matrix is used to build the residual precision matrix Q in
# evaluate_design_efficiency() for AR1 and AR1xAR1 residual structures:
#
#   AR1:      Q_residual = (1/sigma_e^2) * Q_AR1(rho_row) %x% I_cols
#   AR1xAR1:  Q_residual = (1/sigma_e^2) * Q_AR1(rho_col) %x% Q_AR1(rho_row)
#
# Arguments:
#   nn   - positive integer. Length of the AR1 process (number of rows or
#          columns in the field).
#   rho  - numeric in (-1, 1). AR1 autocorrelation parameter.
#
# Returns:
#   A symmetric sparse tridiagonal Matrix of dimension nn x nn.
#   For nn = 1, returns Matrix::Diagonal(1, 1).
#
# Notes:
#   The function does not validate rho. The caller (evaluate_design_efficiency)
#   is responsible for checking |rho| < 1.
# ------------------------------------------------------------------------------
.ar1_precision_sparse <- function(nn, rho) {
  if (nn <= 0) stop("n must be >= 1")
  if (nn == 1) return(Matrix::Diagonal(1, 1))
  a <- 1 / (1 - rho^2)
  d <- rep((1 + rho^2) * a, nn); d[1] <- a; d[nn] <- a
  o <- rep(-rho * a, nn - 1)
  Matrix::sparseMatrix(
    i    = c(seq_len(nn),     seq_len(nn - 1), 2:nn),
    j    = c(seq_len(nn),     2:nn,            seq_len(nn - 1)),
    x    = c(d, o, o),
    dims = c(nn, nn)
  )
}


# ------------------------------------------------------------------------------
# .build_residual_precision
# ------------------------------------------------------------------------------
# Construct the residual precision matrix Q = R^{-1} for a set of plots given
# their (Row, Column) coordinates on an n_rows x n_cols field grid.
#
# Structures:
#   "IID"      : Q = (1/sigma_e2) I
#   "AR1"      : row-only AR1. Covariance = kronecker(I_cols, Sigma_row);
#                precision = kronecker(I_cols, Qrow).
#   "AR1xAR1"  : separable. Covariance = kronecker(Sigma_col, Sigma_row);
#                precision = kronecker(Qcol, Qrow).
#
# CRITICAL indexing note: because the inner (fast-varying) Kronecker factor is
# Qrow (dimension n_rows), the linear index of a plot at (Row r, Column c) into
# the Kronecker grid is
#       grid_index = (c - 1) * n_rows + r.
# Using (r - 1) * n_cols + c (the previous implementation) transposes the row/
# column autocorrelation axes and scrambles adjacency on non-square fields.
# This helper is verified against a hand-built AR1xAR1 precision on a non-square
# grid in tests/testthat/test-ar1-precision.R.
#
# For a fully populated field (one plot per cell) the returned Q is the exact
# reordered full-grid precision. When some cells are empty, Q is the subset
# Qgrid[present, present], i.e. the precision of the observed plots conditional
# on the empty cells -- a deliberate, documented modelling approximation carried
# over from the original evaluators.
# ------------------------------------------------------------------------------
.ensure_full_column_rank <- function(X, Q) {
  # Guard against a singular fixed-effects design (e.g. intercept retained
  # alongside a full set of treatment dummies when check_as_fixed = FALSE and
  # there are no check rows). Rank is assessed on X'QX (p x p). The intercept is
  # dropped first (the usual redundant column); any remaining dependencies are
  # removed by QR pivoting. Returns X with full column rank.
  rk <- function(M) as.integer(Matrix::rankMatrix(
    Matrix::crossprod(M, Q %*% M), method = "qr"))
  if (rk(X) == ncol(X)) return(X)

  if ("(Intercept)" %in% colnames(X)) {
    X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
    message("Fixed-effects design was rank-deficient; dropped '(Intercept)'.")
    if (rk(X) == ncol(X)) return(X)
  }

  qrc     <- qr(as.matrix(Matrix::crossprod(X, Q %*% X)))
  keep    <- sort(qrc$pivot[seq_len(qrc$rank)])
  dropped <- setdiff(colnames(X), colnames(X)[keep])
  if (length(dropped))
    message("Fixed-effects design rank-deficient; dropped column(s): ",
            paste(dropped, collapse = ", "))
  X[, keep, drop = FALSE]
}


.ar1_cov <- function(nn, rho) {
  # AR1 correlation matrix, rho^|i-j| (dense; used by covariance-based structures).
  i <- seq_len(nn)
  rho^abs(outer(i, i, "-"))
}


.spatial_kernel <- function(d, kernel, range, matern_nu = 1.5) {
  # Isotropic correlation as a function of Euclidean plot distance.
  if (range <= 0) stop("`kernel_range` must be positive.")
  switch(kernel,
    exponential = exp(-d / range),
    gaussian    = exp(-(d / range)^2),
    matern      = {
      if (isTRUE(all.equal(matern_nu, 0.5))) {
        exp(-d / range)
      } else if (isTRUE(all.equal(matern_nu, 1.5))) {
        a <- sqrt(3) * d / range; (1 + a) * exp(-a)
      } else if (isTRUE(all.equal(matern_nu, 2.5))) {
        a <- sqrt(5) * d / range; (1 + a + (a^2) / 3) * exp(-a)
      } else {
        stop("`matern_nu` must be 0.5, 1.5 or 2.5.")
      }
    },
    stop("Unknown kernel: ", kernel)
  )
}


.build_residual_precision <- function(row_idx, col_idx, n_rows, n_cols,
                                      residual_structure, rho_row, rho_col,
                                      sigma_e2, nugget = 0,
                                      kernel_range = NULL, matern_nu = 1.5,
                                      max_dense_solve = 3000L) {
  nn <- length(row_idx)
  ri <- as.integer(row_idx)
  ci <- as.integer(col_idx)

  if (residual_structure == "IID") {
    return(Matrix::Diagonal(nn, 1 / sigma_e2))
  }

  grid_index <- (ci - 1L) * n_rows + ri

  # ---- Fast Kronecker-precision path (no nugget) -----------------------------
  if (residual_structure %in% c("AR1", "AR1xAR1")) {
    if (residual_structure == "AR1") {
      Qrow  <- .ar1_precision_sparse(n_rows, rho_row)
      Qgrid <- Matrix::kronecker(Matrix::Diagonal(n_cols, 1), Qrow) * (1 / sigma_e2)
    } else {
      Qrow  <- .ar1_precision_sparse(n_rows, rho_row)
      Qcol  <- .ar1_precision_sparse(n_cols, rho_col)
      Qgrid <- Matrix::kronecker(Qcol, Qrow) * (1 / sigma_e2)
    }
    return(Qgrid[grid_index, grid_index, drop = FALSE])
  }

  # ---- Covariance-based structures -------------------------------------------
  # These build the covariance of the OBSERVED plots and invert it (the marginal
  # form). For a fully populated field this coincides with subsetting the
  # full-grid precision; when cells are empty the marginal form used here is the
  # statistically correct one. Cost is O(n^3), hence the size guard.
  if (nn > max_dense_solve)
    stop("Structure '", residual_structure, "' needs a dense ", nn, " x ", nn,
         " solve, above max_dense_solve (", max_dense_solve, ").")

  if (nugget < 0 || nugget >= 1)
    stop("`nugget` must satisfy 0 <= nugget < 1.")

  if (residual_structure == "AR1xAR1_nugget") {
    Kfull <- kronecker(.ar1_cov(n_cols, rho_col), .ar1_cov(n_rows, rho_row))
    Ksub  <- Kfull[grid_index, grid_index, drop = FALSE]
  } else if (residual_structure %in% c("exponential", "gaussian", "matern")) {
    if (is.null(kernel_range))
      stop("`kernel_range` is required for kernel structure '", residual_structure, "'.")
    d <- as.matrix(stats::dist(cbind(as.numeric(ri), as.numeric(ci))))
    Ksub <- .spatial_kernel(d, residual_structure, kernel_range, matern_nu)
  } else {
    stop("Unknown residual_structure: ", residual_structure)
  }

  Sigma <- (1 - nugget) * Ksub + nugget * diag(nn)
  Qd <- tryCatch(solve(Sigma), error = function(e) .pinv_sym_dense(Sigma))
  Matrix::Matrix((Qd + t(Qd)) / 2 / sigma_e2, sparse = FALSE)
}


# ------------------------------------------------------------------------------
# .pspline_block
# ------------------------------------------------------------------------------
# Tensor-product P-spline (SpATS-style) spatial surface, returned in the
# mixed-model reparameterisation so it can be used as a proper random effect.
#
# A tensor B-spline basis is built over the row and column coordinates, with
# second-order difference penalties on each margin:
#     P = lambda_row * (I_kc %x% D_r) + lambda_col * (D_c %x% I_kr)
# matching the coefficient ordering of Z, whose row for plot (r, c) is
#     kronecker(Bc[c, ], Br[r, ]).
#
# That penalty is IMPROPER: its null space (constant, linear row, linear column
# and their product) is unpenalised, so using it directly leaves those
# directions nearly free and they compete with the genotype effects. We
# therefore eigen-decompose P and keep only the positive-eigenvalue (range)
# space, giving a proper random effect with diagonal precision. The null-space
# directions are intentionally dropped: the intercept and the Row/Column random
# effects already represent them. With this reparameterisation, letting the
# smoothing parameters grow shrinks the surface away and the design criteria
# converge to the no-spline model (verified numerically).
# ------------------------------------------------------------------------------
.pspline_block <- function(row_idx, col_idx, n_rows, n_cols,
                           knots_row = 8L, knots_col = 8L,
                           lambda_row = 1, lambda_col = 1,
                           degree = 3L, tol = 1e-8) {
  if (!requireNamespace("splines", quietly = TRUE))
    stop("Package 'splines' is required for the P-spline spatial surface.")

  kr <- max(as.integer(knots_row), degree + 1L)
  kc <- max(as.integer(knots_col), degree + 1L)
  kr <- min(kr, n_rows); kc <- min(kc, n_cols)
  if (kr < degree + 1L || kc < degree + 1L)
    stop("Field is too small for a P-spline surface; use a residual structure instead.")

  Br <- splines::bs(seq_len(n_rows), df = kr, degree = degree, intercept = TRUE)
  Bc <- splines::bs(seq_len(n_cols), df = kc, degree = degree, intercept = TRUE)
  Br <- matrix(as.numeric(Br), nrow = n_rows)
  Bc <- matrix(as.numeric(Bc), nrow = n_cols)

  dpen <- function(k, order = 2L) {
    D <- diag(k)
    for (i in seq_len(order)) D <- diff(D)
    crossprod(D)
  }
  P <- lambda_row * kronecker(diag(kc), dpen(kr)) +
       lambda_col * kronecker(dpen(kc), diag(kr))

  eg   <- eigen((P + t(P)) / 2, symmetric = TRUE)
  keep <- eg$values > tol * max(eg$values, 1)
  if (!any(keep)) stop("P-spline penalty has no positive eigenvalues.")
  U   <- eg$vectors[, keep, drop = FALSE]
  Lam <- eg$values[keep]

  ri <- as.integer(row_idx); ci <- as.integer(col_idx)
  Zfull <- t(vapply(seq_along(ri), function(i) kronecker(Bc[ci[i], ], Br[ri[i], ]),
                    numeric(kr * kc)))
  list(Z = Zfull %*% U, lambda = Lam)
}


# ------------------------------------------------------------------------------
# .trace_subinv_est
# ------------------------------------------------------------------------------
# Estimate the trace of a submatrix of C^{-1} using the Hutchinson stochastic
# trace estimator with Rademacher random vectors.
#
# For a large symmetric positive definite matrix C and an index set `idx`,
# computing the full submatrix C^{-1}[idx, idx] is expensive when length(idx)
# is large (>= eff_full_max in evaluate_design_efficiency). This estimator
# approximates trace(C^{-1}[idx, idx]) without forming C^{-1} explicitly.
#
# The Hutchinson estimator:
#   trace(A) ~= (1/m) * sum_{k=1}^{m}  z_k' A z_k
#
# where z_k are iid Rademacher vectors (entries in {-1, +1} with equal
# probability). Applied here to A = C^{-1}[idx, idx]:
#
#   For each k:
#     1. Draw z_k ~ Rademacher(length(idx))
#     2. Form rhs: a zero vector of length nrow(C) with z_k at positions idx
#     3. Solve: sol_k = C^{-1} rhs_k
#     4. Accumulate: z_k' sol_k[idx]  =  z_k' C^{-1}[idx,idx] z_k
#
# The estimator is unbiased: E[z' A z] = trace(A) for Rademacher z.
# Variance decreases as O(1/m); more samples reduce variance at linear cost.
#
# Used in evaluate_design_efficiency() for the _APPROX computation modes
# (FIXED_TREATMENT_BLUE_APPROX, RANDOM_TREATMENT_PEV_APPROX) when the number
# of treatments exceeds `eff_full_max`.
#
# Arguments:
#   Cmat        - the mixed model coefficient matrix (dense or sparse).
#   idx         - integer vector. Row/column indices into Cmat corresponding
#                 to the treatment effects whose trace is estimated.
#   m           - positive integer. Number of Rademacher probe vectors.
#                 Default 80. Higher values reduce estimator variance.
#   seed_local  - integer. Seed for the Rademacher draws, run in an isolated
#                 RNG scope via on.exit(). Does not affect the global seed.
#
# Returns:
#   Numeric scalar. Estimated trace(C^{-1}[idx, idx]).
#
# Reference:
#   Hutchinson, M.F. (1990). A stochastic estimator of the trace of the
#   influence matrix for Laplacian smoothing splines. Communications in
#   Statistics - Simulation and Computation, 19(2), 433-450.
# ------------------------------------------------------------------------------
.trace_subinv_est <- function(Cmat, idx, m = 80, seed_local = 1) {
  old <- NULL
  has <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (has) old <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  set.seed(seed_local)
  on.exit({
    if (has) assign(".Random.seed", old, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)

  p   <- length(idx)
  acc <- 0
  for (ii in seq_len(m)) {
    z   <- sample(c(-1, 1), p, replace = TRUE)
    rhs <- Matrix::Matrix(0, nrow = nrow(Cmat), ncol = 1, sparse = TRUE)
    rhs[idx, 1] <- z
    sol <- .solve_C(Cmat, rhs)
    acc <- acc + sum(z * sol[idx, 1])
  }
  acc / m
}


# ==============================================================================
# RNG SANDBOX
# ==============================================================================

# ------------------------------------------------------------------------------
# .with_local_seed
# ------------------------------------------------------------------------------
# Evaluate an expression under a specified seed without permanently altering
# the global RNG state.
#
# Saves the current .Random.seed, sets seed_local, evaluates expr, then
# restores the original .Random.seed via on.exit(). This ensures that calling
# code with a fixed seed (e.g. clustering, dispersion optimisation) does not
# consume random numbers from the main design seed stream.
#
# Used in alpha_rc_stream() for:
#   - k-means clustering initialisation (cluster_seed)
#   - genetic dispersion optimisation (dispersion_seed)
#
# Used in evaluate_design_efficiency() via .trace_subinv_est() for the
# Hutchinson estimator seed.
#
# Arguments:
#   seed_local  - integer. Seed to set before evaluating expr.
#   expr        - R expression to evaluate under seed_local. Evaluated with
#                 force() to prevent lazy evaluation issues.
#
# Returns:
#   The return value of expr.
#
# Notes:
#   If .Random.seed does not exist in .GlobalEnv before the call (i.e. no RNG
#   has been used yet in the session), the saved state is NULL and the seed is
#   simply removed on exit rather than restored.
# ------------------------------------------------------------------------------
.with_local_seed <- function(seed_local, expr) {
  old <- NULL
  has <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (has) old <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  set.seed(seed_local)
  on.exit({
    if (has) assign(".Random.seed", old, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  force(expr)
}


# ==============================================================================
# SPATIAL HELPERS
# ==============================================================================

# ------------------------------------------------------------------------------
# .build_neighbor_pairs
# ------------------------------------------------------------------------------
# Enumerate all pairs of plots that are within Chebyshev distance `radius` of
# each other on a 2D field grid.
#
# The Chebyshev distance between two plots at (r1, c1) and (r2, c2) is:
#   d_Cheb = max(|r1 - r2|, |c1 - c2|)
#
# Two plots are neighbours if d_Cheb <= radius. This includes diagonal
# neighbours: radius = 1 captures up to 8 immediate neighbours; radius = 2
# captures up to 24 neighbours.
#
# The function returns each pair (i, j) with i < j exactly once (upper
# triangle only), iterating over all plots i in order and collecting all
# j > i that satisfy the distance criterion.
#
# Used in alpha_rc_stream() by apply_genetic_dispersion() to score the total
# relatedness of neighbouring plots under the dispersion matrix.
#
# Arguments:
#   row     - integer vector of length n. Row coordinate of each plot.
#   col     - integer vector of length n. Column coordinate of each plot.
#   radius  - positive integer. Chebyshev radius defining the neighbourhood.
#             Default 1 (8-connected neighbourhood).
#
# Returns:
#   Integer matrix with two columns named "i" and "j", one row per
#   neighbouring pair (i < j). Returns a 0-row matrix with the same column
#   names if n < 2 or no pairs fall within radius.
# ------------------------------------------------------------------------------
.build_neighbor_pairs <- function(row, col, radius = 1) {
  n <- length(row)
  if (n < 2) return(matrix(integer(0), ncol = 2,
                            dimnames = list(NULL, c("i", "j"))))
  out <- vector("list", 0); kk <- 1L
  for (i in seq_len(n - 1)) {
    dr <- abs(row[i] - row[(i + 1):n])
    dc <- abs(col[i] - col[(i + 1):n])
    ok <- pmax(dr, dc) <= radius
    if (any(ok)) {
      jj        <- ((i + 1):n)[ok]
      out[[kk]] <- cbind(i = rep(i, length(jj)), j = jj)
      kk        <- kk + 1L
    }
  }
  if (length(out) == 0) return(matrix(integer(0), ncol = 2,
                                      dimnames = list(NULL, c("i", "j"))))
  ans <- do.call(rbind, out)
  colnames(ans) <- c("i", "j")
  ans
}


# ------------------------------------------------------------------------------
# .score_dispersion
# ------------------------------------------------------------------------------
# Compute the total pairwise relatedness of neighbouring non-check plots under
# a genomic or pedigree relationship matrix.
#
# For each pair (i, j) in `pairs`, if both plots i and j are non-check entries
# with known line IDs in Ksub, the relatedness Ksub[line_i, line_j] is added
# to the score. Pairs involving at least one check or one unrecognised line
# contribute zero.
#
# This score is the objective function minimised by the dispersion
# optimisation in alpha_rc_stream(). A lower score means genetically similar
# lines are placed further apart in the field, reducing spatial confounding
# between genetic and environmental effects.
#
# Used in alpha_rc_stream() by apply_genetic_dispersion() at each swap
# proposal to decide whether to accept the swap (accept if new score < current
# score).
#
# Arguments:
#   trt_vec   - character vector of length n. Treatment label for each plot
#               (NA for unused plots).
#   is_check  - logical vector of length n. TRUE if plot is a check.
#   Ksub      - numeric matrix. Submatrix of the relationship matrix restricted
#               to the non-check entry treatments in the current field layout.
#               rownames and colnames must match the treatment labels of
#               non-check entries.
#   pairs     - integer matrix with columns "i" and "j" as returned by
#               .build_neighbor_pairs(). Each row is a neighbouring plot pair.
#
# Returns:
#   Numeric scalar. Sum of Ksub[line_i, line_j] over all valid non-check
#   neighbour pairs. Returns 0 if pairs has 0 rows or no valid pairs exist.
# ------------------------------------------------------------------------------
.score_dispersion <- function(trt_vec, is_check, Ksub, pairs) {
  if (nrow(pairs) == 0) return(0)
  line_levels <- colnames(Ksub)
  line_idx    <- rep(NA_integer_, length(trt_vec))
  line_idx[!is_check] <- match(trt_vec[!is_check], line_levels)
  ii <- pairs[, "i"]; jj <- pairs[, "j"]
  li <- line_idx[ii]; lj <- line_idx[jj]
  ok <- !is.na(li) & !is.na(lj)
  if (!any(ok)) return(0)
  sum(Ksub[cbind(li[ok], lj[ok])])
}
