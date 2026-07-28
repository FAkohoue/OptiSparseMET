# Tests for robust consensus of relationship matrices.

mk_spd <- function(n, seed) {
  set.seed(seed); A <- matrix(rnorm(n * (n + 3)), n); tcrossprod(A) / (n + 3) + diag(n)
}

test_that("'mean' equals the arithmetic mean; all methods stay symmetric PD", {
  Ds <- list(mk_spd(5, 1), mk_spd(5, 2), mk_spd(5, 3))
  Cm <- consensus_relationship(Ds, "mean")
  expect_equal(Cm, Reduce(`+`, Ds) / 3, tolerance = 1e-10)
  for (m in c("mean", "geometric", "median", "rv_weighted")) {
    C <- consensus_relationship(Ds, m)
    expect_equal(C, t(C), tolerance = 1e-8)
    expect_true(all(eigen(C, symmetric = TRUE, only.values = TRUE)$values > -1e-8))
  }
})

test_that("STATIS (rv_weighted) down-weights an anomalous year", {
  # two similar years and one very different (outlier) year
  base <- mk_spd(6, 10)
  Ds <- list(y1 = base + 0.01 * mk_spd(6, 11),
             y2 = base + 0.01 * mk_spd(6, 12),
             y3 = mk_spd(6, 99) * 3)                 # anomalous
  C <- consensus_relationship(Ds, "rv_weighted")
  w <- attr(C, "weights")
  expect_length(w, 3L)
  expect_equal(sum(w), 1, tolerance = 1e-8)
  # the outlier year gets the smallest weight
  expect_true(which.min(w) == 3L)
})

test_that("geometric mean differs from the Euclidean mean (SPD geometry)", {
  Ds <- list(diag(c(1, 4)), diag(c(4, 1)))
  Cmean <- consensus_relationship(Ds, "mean")        # diag(2.5, 2.5)
  Cgeo  <- consensus_relationship(Ds, "geometric")   # diag(2, 2) = sqrt(1*4)
  expect_equal(diag(Cgeo), c(2, 2), tolerance = 1e-6)
  expect_false(isTRUE(all.equal(diag(Cmean), diag(Cgeo))))
})

test_that("empty / mismatched inputs error", {
  expect_error(consensus_relationship(list()), "empty")
  A <- matrix(1, 2, 2, dimnames = list(c("a","b"), c("a","b")))
  B <- matrix(1, 2, 2, dimnames = list(c("a","c"), c("a","c")))
  expect_error(consensus_relationship(list(A, B)), "same individuals")
})
