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
.safe_logdet_psd_dense <- function(M, tol = 1e-10) {
  ev <- eigen((M + t(M)) / 2, symmetric = TRUE, only.values = TRUE)$values
  ev <- ev[ev > tol]
  if (length(ev) == 0) return(NA_real_)
  sum(log(ev))
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
