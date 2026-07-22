# Tests for .safe_logdet_psd_dense() (P0.3): scale-robust pseudo-determinant
# with an exact expected rank.

# Fixed orthonormal basis so eigenvectors are known and reproducible.
Qmat <- qr.Q(qr(matrix(c(1, 0.3, -0.2,
                         0.3, 1, 0.1,
                        -0.2, 0.1, 1), 3, 3)))
psd_with_eigs <- function(lams) Qmat %*% diag(lams) %*% t(Qmat)

test_that("pseudo-determinant drops the structural zero and matches the oracle", {
  M <- psd_with_eigs(c(4, 1, 0))          # rank 2, one structural zero
  ld <- .safe_logdet_psd_dense(M, expected_rank = 2)
  expect_equal(ld, log(4) + log(1), tolerance = 1e-10)
  expect_equal(exp(ld / 2), 2, tolerance = 1e-10)   # geometric mean of {4,1}
})

test_that("result scales linearly with the matrix (scale invariance of D)", {
  M  <- psd_with_eigs(c(4, 1, 0))
  d1 <- exp(.safe_logdet_psd_dense(M, expected_rank = 2) / 2)
  d2 <- exp(.safe_logdet_psd_dense(1000 * M, expected_rank = 2) / 2)
  expect_equal(d2, 1000 * d1, tolerance = 1e-8)      # 2000 vs 2
})

test_that("a genuinely small eigenvalue is retained, not dropped as zero", {
  M  <- psd_with_eigs(c(4, 1e-6, 0))
  ld <- .safe_logdet_psd_dense(M, expected_rank = 2)
  # Value is sqrt(4 * 1e-6) = 2e-3; loosen tolerance because eigen() recovers the
  # tiny eigenvalue with ~1e-9 relative error, which the log amplifies.
  expect_equal(exp(ld / 2), sqrt(4 * 1e-6), tolerance = 1e-6)  # 2e-3
})

test_that("returns NA when fewer positive eigenvalues than expected_rank", {
  M <- psd_with_eigs(c(4, 0, 0))          # rank 1, expected 2
  expect_true(is.na(.safe_logdet_psd_dense(M, expected_rank = 2)))
})
