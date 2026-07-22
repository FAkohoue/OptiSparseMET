# Ground-truth tests for .build_residual_precision() (P0.1).
#
# The residual precision for a separable AR1xAR1 field is
#   Q = (1/sigma_e2) * solve(kronecker(Sigma_col, Sigma_row))
# where Sigma_row[i,j] = rho_row^|i-j| and Sigma_col[i,j] = rho_col^|i-j|.
# The linear index of a plot at (Row r, Column c) into that Kronecker grid is
# (c - 1) * n_rows + r. These tests build the oracle by hand on a NON-SQUARE
# grid with rho_row != rho_col, the case that exposed the original transposed
# index. They fail on the pre-fix code and pass after it.

ar1_cov <- function(n, rho) outer(seq_len(n), seq_len(n),
                                  function(i, j) rho^abs(i - j))

test_that("AR1xAR1 precision matches a hand-built oracle on a non-square grid", {
  n_rows <- 3L; n_cols <- 4L
  rho_row <- 0.5; rho_col <- 0.2; sigma_e2 <- 2

  Sr <- ar1_cov(n_rows, rho_row)
  Sc <- ar1_cov(n_cols, rho_col)
  oracle_full <- solve(kronecker(Sc, Sr)) / sigma_e2   # 12 x 12

  # A scrambled plot order over all 12 cells (r fast, then permuted).
  cells <- expand.grid(Row = seq_len(n_rows), Column = seq_len(n_cols))
  set.seed(42); ord <- sample(nrow(cells))
  cells <- cells[ord, ]
  gi <- (cells$Column - 1L) * n_rows + cells$Row       # oracle index

  Q <- .build_residual_precision(cells$Row, cells$Column, n_rows, n_cols,
                                 "AR1xAR1", rho_row, rho_col, sigma_e2)
  Q <- as.matrix(Q)
  expected <- oracle_full[gi, gi]

  expect_equal(unname(Q), unname(expected), tolerance = 1e-10)
})

test_that("row-only AR1 precision matches kronecker(I_cols, Qrow) oracle", {
  n_rows <- 3L; n_cols <- 4L
  rho_row <- 0.6; sigma_e2 <- 1.5

  Sr <- ar1_cov(n_rows, rho_row)
  oracle_full <- kronecker(diag(n_cols), solve(Sr)) / sigma_e2

  cells <- expand.grid(Row = seq_len(n_rows), Column = seq_len(n_cols))
  gi <- (cells$Column - 1L) * n_rows + cells$Row

  Q <- as.matrix(.build_residual_precision(cells$Row, cells$Column,
                                           n_rows, n_cols, "AR1", rho_row, 0, sigma_e2))
  expect_equal(unname(Q), unname(oracle_full[gi, gi]), tolerance = 1e-10)
})

test_that("row and column autocorrelation axes are distinct (not transposed)", {
  # On a non-square grid, swapping rho_row and rho_col must change Q. The old
  # buggy index made the two axes effectively interchangeable.
  n_rows <- 3L; n_cols <- 5L
  cells <- expand.grid(Row = seq_len(n_rows), Column = seq_len(n_cols))

  Q1 <- as.matrix(.build_residual_precision(cells$Row, cells$Column,
                                            n_rows, n_cols, "AR1xAR1", 0.7, 0.1, 1))
  Q2 <- as.matrix(.build_residual_precision(cells$Row, cells$Column,
                                            n_rows, n_cols, "AR1xAR1", 0.1, 0.7, 1))
  expect_false(isTRUE(all.equal(Q1, Q2)))
})

test_that("IID structure returns a scaled identity", {
  Q <- .build_residual_precision(c(1L, 2L, 3L), c(1L, 1L, 1L),
                                 3L, 1L, "IID", 0, 0, 4)
  expect_equal(as.matrix(Q), diag(1 / 4, 3), tolerance = 1e-12)
})
