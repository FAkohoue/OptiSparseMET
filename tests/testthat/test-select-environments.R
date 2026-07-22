# Tests for select_environments() (P4.1 / decision 1).

mkD <- function() {
  E <- c("E1","E2","E3","E4","E5")
  D <- matrix(c(
    1.00, 0.98, 0.30, 0.20, 0.10,
    0.98, 1.00, 0.30, 0.20, 0.10,
    0.30, 0.30, 1.00, 0.40, 0.20,
    0.20, 0.20, 0.40, 1.00, 0.50,
    0.10, 0.10, 0.20, 0.50, 1.00), 5, 5, byrow = TRUE,
    dimnames = list(E, E))
  D
}

test_that("optcontrib (spread) avoids the near-duplicate environment pair", {
  D <- mkD()
  sel <- select_environments(D, n = 2, method = "optcontrib")$selected
  expect_false(all(c("E1", "E2") %in% sel))   # E1,E2 are near-identical
  expect_equal(length(unique(sel)), 2L)
})

test_that("all methods return n distinct environments", {
  D <- mkD()
  for (m in c("representative", "optcontrib", "kmeans", "hclust", "random")) {
    sel <- select_environments(D, n = 3, method = m, seed = 1)$selected
    expect_equal(length(unique(sel)), 3L, info = m)
    expect_true(all(sel %in% rownames(D)), info = m)
  }
})

test_that("representative coverage beats a redundant selection", {
  D <- mkD()
  rep_sel <- select_environments(D, n = 2, method = "representative")
  cov_redundant <- sum(apply(D[, c("E1","E2"), drop = FALSE], 1, max))
  expect_gte(rep_sel$score, cov_redundant)   # coverage-optimal >= redundant pair
})

test_that("random selection is reproducible with a seed", {
  D <- mkD()
  a <- select_environments(D, n = 3, method = "random", seed = 42)$selected
  b <- select_environments(D, n = 3, method = "random", seed = 42)$selected
  expect_identical(a, b)
})
