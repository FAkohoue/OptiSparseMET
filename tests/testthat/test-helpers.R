library(Matrix)

# -- .make_sparse_incidence ----------------------------------------------------

test_that(".make_sparse_incidence returns correct structure", {
  v   <- c("A", "B", "A", NA, "C")
  res <- OptiSparseMET:::.make_sparse_incidence(v)
  
  expect_type(res, "list")
  expect_true("M" %in% names(res))
  expect_true("levels" %in% names(res))
  expect_equal(nrow(res$M), 5L)
  expect_equal(ncol(res$M), 3L)
  expect_equal(sort(res$levels), c("A", "B", "C"))
})

test_that(".make_sparse_incidence: NA rows are all-zero", {
  v   <- c("A", NA, "B")
  res <- OptiSparseMET:::.make_sparse_incidence(v)
  expect_true(all(res$M[2, ] == 0))
})

test_that(".make_sparse_incidence: empty levels returns 0-column matrix", {
  res <- OptiSparseMET:::.make_sparse_incidence(c(NA, NA))
  expect_equal(ncol(res$M), 0L)
  expect_equal(length(res$levels), 0L)
})

test_that(".make_sparse_incidence: each row sums to 0 or 1", {
  v   <- c("X", "Y", "X", NA, "Z")
  res <- OptiSparseMET:::.make_sparse_incidence(v)
  # as.matrix() required: rowSums() on a sparse Matrix with one column
  # treats it as a vector and fails with a dimension error.
  row_sums <- rowSums(as.matrix(res$M))
  expect_true(all(row_sums %in% c(0, 1)))
})

# -- .ar1_precision_sparse ----------------------------------------------------

test_that(".ar1_precision_sparse returns a square sparse matrix", {
  Q <- OptiSparseMET:::.ar1_precision_sparse(5L, 0.3)
  expect_true(inherits(Q, "Matrix"))
  expect_equal(dim(Q), c(5L, 5L))
})

test_that(".ar1_precision_sparse: n=1 returns 1x1 identity", {
  Q <- OptiSparseMET:::.ar1_precision_sparse(1L, 0.5)
  expect_equal(as.matrix(Q), matrix(1, 1, 1))
})

test_that(".ar1_precision_sparse is symmetric positive definite", {
  Q  <- OptiSparseMET:::.ar1_precision_sparse(6L, 0.4)
  Qd <- as.matrix(Q)
  expect_equal(Qd, t(Qd))
  ev <- eigen(Qd, only.values = TRUE)$values
  expect_true(all(ev > 0))
})

test_that(".ar1_precision_sparse rho=0 gives identity * 1/(1-0^2) = I", {
  Q  <- as.matrix(OptiSparseMET:::.ar1_precision_sparse(4L, 0))
  expect_equal(Q, diag(4L), tolerance = 1e-10)
})

# -- .pairwise_diff_mean_var ---------------------------------------------------

test_that(".pairwise_diff_mean_var returns NA for p < 2", {
  V <- matrix(1, 1, 1)
  expect_true(is.na(OptiSparseMET:::.pairwise_diff_mean_var(V)))
})

test_that(".pairwise_diff_mean_var: diagonal V gives non-negative result", {
  V   <- diag(c(1, 2, 3))
  val <- OptiSparseMET:::.pairwise_diff_mean_var(V)
  expect_true(is.numeric(val))
  expect_true(val >= 0)
})

test_that(".pairwise_diff_mean_var: equal variances give deterministic result", {
  # For V = s^2 * I, Var(ti - tj) = 2*s^2 for all pairs
  s2  <- 0.5
  V   <- diag(rep(s2, 4L))
  val <- OptiSparseMET:::.pairwise_diff_mean_var(V)
  expect_equal(val, 2 * s2, tolerance = 1e-10)
})

# -- .pinv_sym_dense -----------------------------------------------------------

test_that(".pinv_sym_dense inverts a full-rank symmetric matrix", {
  A    <- crossprod(matrix(rnorm(9), 3, 3)) + diag(3)
  Ainv <- OptiSparseMET:::.pinv_sym_dense(A)
  expect_equal(A %*% Ainv, diag(3), tolerance = 1e-8)
})

test_that(".pinv_sym_dense returns zero matrix for zero input", {
  A    <- matrix(0, 3, 3)
  Ainv <- OptiSparseMET:::.pinv_sym_dense(A)
  expect_equal(Ainv, matrix(0, 3, 3))
})

# -- .safe_logdet_psd_dense ---------------------------------------------------

test_that(".safe_logdet_psd_dense returns correct log-determinant", {
  A      <- diag(c(2, 3, 4))
  logdet <- OptiSparseMET:::.safe_logdet_psd_dense(A)
  expect_equal(logdet, log(2) + log(3) + log(4), tolerance = 1e-10)
})

test_that(".safe_logdet_psd_dense returns NA for zero matrix", {
  A <- matrix(0, 3, 3)
  expect_true(is.na(OptiSparseMET:::.safe_logdet_psd_dense(A)))
})

# -- .build_neighbor_pairs ----------------------------------------------------

test_that(".build_neighbor_pairs returns matrix with columns i and j", {
  row   <- c(1, 1, 2, 2)
  col   <- c(1, 2, 1, 2)
  pairs <- OptiSparseMET:::.build_neighbor_pairs(row, col, radius = 1)
  
  expect_true(is.matrix(pairs))
  expect_equal(colnames(pairs), c("i", "j"))
})

test_that(".build_neighbor_pairs: all pairs have i < j", {
  row   <- c(1, 1, 2, 2, 3)
  col   <- c(1, 2, 1, 2, 1)
  pairs <- OptiSparseMET:::.build_neighbor_pairs(row, col, radius = 1)
  if (nrow(pairs) > 0) expect_true(all(pairs[, "i"] < pairs[, "j"]))
})

test_that(".build_neighbor_pairs: n < 2 returns empty matrix", {
  pairs <- OptiSparseMET:::.build_neighbor_pairs(1L, 1L, radius = 1)
  expect_equal(nrow(pairs), 0L)
})

test_that(".build_neighbor_pairs: larger radius captures more pairs", {
  row <- 1:5; col <- 1:5
  p1  <- OptiSparseMET:::.build_neighbor_pairs(row, col, radius = 1)
  p2  <- OptiSparseMET:::.build_neighbor_pairs(row, col, radius = 2)
  expect_true(nrow(p2) >= nrow(p1))
})

# -- .with_local_seed ---------------------------------------------------------

test_that(".with_local_seed isolates RNG state", {
  set.seed(42L)
  r_before <- runif(1)
  
  set.seed(42L)
  val <- OptiSparseMET:::.with_local_seed(99L, runif(1))
  r_after <- runif(1)
  
  # The local seed should produce different values than seed 42
  set.seed(99L); r_local <- runif(1)
  expect_equal(val, r_local)
  
  # The outer RNG state should be restored: r_after equals what runif(1)
  # would give after set.seed(42) + runif(1)
  expect_equal(r_after, r_before, tolerance = 1e-15)
})