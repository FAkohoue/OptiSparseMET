# Tests for the extended spatial models (nugget, distance kernels, P-spline)
# and compare_spatial_models(). Invariants verified independently in numpy.

ar1_cov_ref <- function(n, rho) outer(seq_len(n), seq_len(n),
                                      function(i, j) rho^abs(i - j))

cells_fb <- function(n_rows, n_cols) {
  expand.grid(Row = seq_len(n_rows), Column = seq_len(n_cols))
}

# ---- S1: AR1xAR1 + nugget ----------------------------------------------------

test_that("nugget = 0 reproduces pure AR1xAR1 and nugget -> 1 approaches IID", {
  nr <- 3L; nc <- 4L; rr <- 0.6; rc <- 0.3; se2 <- 2
  cl <- cells_fb(nr, nc)

  q_pure <- as.matrix(.build_residual_precision(cl$Row, cl$Column, nr, nc,
                                                "AR1xAR1", rr, rc, se2))
  q_nug0 <- as.matrix(.build_residual_precision(cl$Row, cl$Column, nr, nc,
                                                "AR1xAR1_nugget", rr, rc, se2,
                                                nugget = 0))
  expect_equal(unname(q_nug0), unname(q_pure), tolerance = 1e-8)

  q_nug <- as.matrix(.build_residual_precision(cl$Row, cl$Column, nr, nc,
                                               "AR1xAR1_nugget", rr, rc, se2,
                                               nugget = 0.999))
  expect_equal(unname(q_nug), diag(1 / se2, nr * nc), tolerance = 1e-2)
})

test_that("nugget precision is symmetric and positive definite", {
  nr <- 3L; nc <- 4L
  cl <- cells_fb(nr, nc)
  Q <- as.matrix(.build_residual_precision(cl$Row, cl$Column, nr, nc,
                                           "AR1xAR1_nugget", 0.6, 0.3, 1,
                                           nugget = 0.25))
  expect_equal(Q, t(Q), tolerance = 1e-10)
  expect_true(all(eigen(Q, symmetric = TRUE, only.values = TRUE)$values > 0))
})

test_that("nugget outside [0,1) is rejected", {
  cl <- cells_fb(3L, 3L)
  expect_error(.build_residual_precision(cl$Row, cl$Column, 3L, 3L,
                                         "AR1xAR1_nugget", .5, .5, 1, nugget = 1))
})

# ---- S2: isotropic distance kernels -----------------------------------------

test_that("distance kernels give valid correlation matrices", {
  nr <- 3L; nc <- 4L
  cl <- cells_fb(nr, nc)
  for (k in c("exponential", "gaussian", "matern")) {
    Q <- as.matrix(.build_residual_precision(cl$Row, cl$Column, nr, nc, k,
                                             0, 0, 1, nugget = 0.1,
                                             kernel_range = 2, matern_nu = 1.5))
    expect_equal(Q, t(Q), tolerance = 1e-10, info = k)
    expect_true(all(eigen(Q, symmetric = TRUE, only.values = TRUE)$values > 0),
                info = k)
  }
})

test_that("exponential kernel with a tiny range approaches IID", {
  nr <- 3L; nc <- 3L; se2 <- 1.5
  cl <- cells_fb(nr, nc)
  Q <- as.matrix(.build_residual_precision(cl$Row, cl$Column, nr, nc,
                                           "exponential", 0, 0, se2,
                                           nugget = 0, kernel_range = 1e-4))
  expect_equal(unname(Q), diag(1 / se2, nr * nc), tolerance = 1e-6)
})

test_that("kernel structures require a range", {
  cl <- cells_fb(3L, 3L)
  expect_error(.build_residual_precision(cl$Row, cl$Column, 3L, 3L,
                                         "exponential", 0, 0, 1))
})

# ---- S3: P-spline block ------------------------------------------------------

test_that("P-spline block is reparameterised to a proper (positive) penalty", {
  nr <- 6L; nc <- 8L
  cl <- cells_fb(nr, nc)
  sp <- .pspline_block(cl$Row, cl$Column, nr, nc,
                       knots_row = 5, knots_col = 6,
                       lambda_row = 1, lambda_col = 1)
  expect_equal(nrow(sp$Z), nr * nc)
  expect_equal(ncol(sp$Z), length(sp$lambda))
  # all retained eigenvalues strictly positive => proper random effect
  expect_true(all(sp$lambda > 0))
  # a 2nd-order difference penalty on a tensor basis has a 4-dim null space
  # (constant, linear row, linear column, bilinear), which is dropped.
  expect_equal(ncol(sp$Z), 5 * 6 - 4)
})

test_that("larger smoothing parameters shrink the spline surface away", {
  nr <- 6L; nc <- 8L
  cl <- cells_fb(nr, nc)
  lo <- .pspline_block(cl$Row, cl$Column, nr, nc, 5, 6, 1, 1)$lambda
  hi <- .pspline_block(cl$Row, cl$Column, nr, nc, 5, 6, 1e4, 1e4)$lambda
  # penalty eigenvalues scale with lambda, so the prior precision grows
  expect_gt(min(hi), min(lo))
})

# ---- S4: compare_spatial_models() -------------------------------------------

test_that("compare_spatial_models returns one row per model", {
  fb <- data.frame(
    Plot = 1:12,
    Row = rep(1:3, times = 4), Column = rep(1:4, each = 3),
    Rep = rep(1:2, each = 6), IBlock = rep(1:4, each = 3),
    BlockInRep = rep(1:2, times = 6),
    Treatment = c(paste0("G", 1:6), paste0("G", 1:6)), Check = FALSE,
    stringsAsFactors = FALSE)
  vc <- list(sigma_e2 = 2, sigma_g2 = 1, sigma_rep2 = 1,
             sigma_ib2 = 1, sigma_r2 = 1, sigma_c2 = 1)

  res <- compare_spatial_models(
    field_book = fb, n_rows = 3, n_cols = 4,
    check_treatments = character(0),
    models = c("IID", "AR1xAR1", "AR1xAR1_nugget", "exponential"),
    treatment_effect = "random", prediction_type = "IID", varcomp = vc)

  expect_equal(nrow(res), 4L)
  expect_true(all(c("model", "A_criterion", "CDmean", "mean_PEV") %in% names(res)))
  expect_true(all(is.na(res$error)))
  expect_true(all(res$mean_PEV > 0))
})
