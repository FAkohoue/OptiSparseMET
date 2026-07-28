# Tests for suggest_common_range() -- the C-range calculator.

test_that("range respects the hard cap and starts at zero", {
  r <- suggest_common_range(n_test_entries_per_environment = 24,
                            n_treatments = 80, n_environments = 4)
  expect_equal(min(r$common_values), 0L)
  expect_true(max(r$common_values) <= 24 - 1)     # C < k
  expect_true(all(r$common_values == sort(r$common_values)))
})

test_that("disjoint designs need more common treatments than random ones", {
  rand <- suggest_common_range(24, 80, 4, sparse_allocation = "random")
  disj <- suggest_common_range(24, 80, 4, sparse_allocation = "disjoint")
  # random already has ~k^2/J shared per pair, so it needs fewer commons
  expect_gt(disj$rationale$connectivity_upper, rand$rationale$connectivity_upper)
  expect_equal(disj$rationale$expected_sparse_overlap, 0)
})

test_that("weak correlations demand a larger connectivity target", {
  strong <- suggest_common_range(24, 80, 4, sparse_allocation = "disjoint", rho = 0.6)
  weak   <- suggest_common_range(24, 80, 4, sparse_allocation = "disjoint", rho = 0.05)
  expect_gt(weak$rationale$target_shared, strong$rationale$target_shared)
  expect_gte(weak$rationale$connectivity_upper, strong$rationale$connectivity_upper)
})

test_that("a large program with ample sparse overlap barely needs common treatments", {
  r <- suggest_common_range(n_test_entries_per_environment = 200,
                            n_treatments = 1000, n_environments = 8,
                            sparse_allocation = "random", rho = 0.6)
  # sparse overlap (~k^2/J = 40) exceeds the target, so no commons are required
  expect_equal(r$rationale$connectivity_upper, 0)
  expect_true(max(r$common_values) <= 6)
})

test_that("the suggested range plugs into tune_common_treatments", {
  set.seed(1)
  J <- 60L; B <- matrix(rnorm(J * J), J, J); G <- tcrossprod(B) / J + diag(J) * 0.4
  dimnames(G) <- list(paste0("L", 1:J), paste0("L", 1:J))
  SigE <- diag(4); dimnames(SigE) <- list(paste0("E", 1:4), paste0("E", 1:4))
  rr <- suggest_common_range(20, J, 4, sparse_allocation = "disjoint", n_points = 4)
  res <- tune_common_treatments(
    treatments = rownames(G), environments = colnames(SigE),
    n_test_entries_per_environment = 20, G = G, Sigma_E = SigE,
    common_values = rr$common_values, sparse_allocation = "disjoint",
    n_sim = 10, seed = 1)
  expect_equal(nrow(res$table), length(rr$common_values))
})
