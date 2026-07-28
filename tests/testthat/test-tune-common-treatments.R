# Tests for tune_common_treatments(). Invariants verified independently in numpy:
#   - a disconnected (disjoint) design has an INTERIOR optimum number of common
#     treatments under estimated Sigma_E (C* ~ 6 in the reference run);
#   - weakly-correlated environments push that optimum HIGHER, not lower;
#   - a known-Sigma_E criterion does not reward common treatments (prefers few).

mk <- function(J = 60L, seed = 0) {
  set.seed(seed)
  B <- matrix(rnorm(J * J), J, J); G <- tcrossprod(B) / J + diag(J) * 0.4
  dimnames(G) <- list(paste0("L", seq_len(J)), paste0("L", seq_len(J)))
  strong <- matrix(c(1, .7, .3, 0,  .7, 1, .4, .1,
                     .3, .4, 1, .6,  0, .1, .6, 1), 4, 4)
  weak <- diag(4) * 0.9 + 0.1
  dimnames(strong) <- dimnames(weak) <- list(paste0("E", 1:4), paste0("E", 1:4))
  list(G = G, strong = strong, weak = weak)
}

test_that("tune_common_treatments returns the four diagnostic curves", {
  d <- mk()
  res <- tune_common_treatments(
    treatments = rownames(d$G), environments = colnames(d$strong),
    n_test_entries_per_environment = 20, G = d$G, Sigma_E = d$strong,
    common_values = c(0, 4, 8, 12), sparse_allocation = "disjoint",
    n_sim = 15, seed = 1)
  expect_true(all(c("n_common", "accuracy_known", "accuracy_estimated",
                    "min_shared", "shared_var", "coverage") %in% names(res$table)))
  expect_equal(nrow(res$table), 4L)
  # connectivity increases and coverage decreases with more common treatments
  expect_true(all(diff(res$table$min_shared) >= 0))
  expect_true(all(diff(res$table$coverage) <= 0))
  expect_true(res$optima$best_estimated %in% res$table$n_common)
})

test_that("common treatments help under a disconnected design (estimated Sigma_E)", {
  d <- mk()
  res <- tune_common_treatments(
    treatments = rownames(d$G), environments = colnames(d$strong),
    n_test_entries_per_environment = 24, G = d$G, Sigma_E = d$strong,
    common_values = c(0, 6, 12), sparse_allocation = "disjoint",
    n_sim = 30, seed = 1)
  # A disconnected design has an interior optimum under estimated Sigma_E: the
  # best count is greater than zero (numpy reference: C* ~ 6).
  expect_gt(res$optima$best_estimated, 0)
  # ...and the known-Sigma_E criterion does NOT reward common treatments as much:
  # its optimum is at or below the estimated one.
  expect_lte(res$optima$best_known, res$optima$best_estimated)
})

test_that("weakly-correlated environments push the optimum higher than strong ones", {
  d <- mk()
  cvals <- c(0, 6, 12, 16)
  common_args <- list(
    treatments = rownames(d$G), environments = colnames(d$strong),
    n_test_entries_per_environment = 24, G = d$G,
    common_values = cvals, sparse_allocation = "disjoint", n_sim = 20, seed = 1)
  rs <- do.call(tune_common_treatments, c(common_args, list(Sigma_E = d$strong)))
  rw <- do.call(tune_common_treatments, c(common_args, list(Sigma_E = d$weak)))
  expect_gte(rw$optima$best_estimated, rs$optima$best_estimated)
})
