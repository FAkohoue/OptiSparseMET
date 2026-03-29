# ==============================================================================
# test-check_balanced_incomplete_feasibility.R
# ==============================================================================

test_that("exact feasibility detected: difference == 0", {
  res <- check_balanced_incomplete_feasibility(
    n_treatments_total             = 120,
    n_environments                 = 4,
    n_test_entries_per_environment = 65,
    target_replications            = 2,
    n_common_treatments            = 10
  )
  expect_true(res$feasible)
  expect_equal(res$difference, 0L)
  expect_equal(res$n_sparse_treatments, 110L)
  expect_equal(res$total_sparse_slots, 4L * 55L)
  expect_equal(res$required_sparse_slots, 110L * 2L)
})

test_that("slot deficit detected: difference negative", {
  res <- check_balanced_incomplete_feasibility(
    n_treatments_total             = 120,
    n_environments                 = 4,
    n_test_entries_per_environment = 60,
    target_replications            = 2,
    n_common_treatments            = 10
  )
  expect_false(res$feasible)
  expect_equal(res$difference, -20L)
})

test_that("slot surplus detected: difference positive", {
  res <- check_balanced_incomplete_feasibility(
    n_treatments_total             = 120,
    n_environments                 = 4,
    n_test_entries_per_environment = 70,
    target_replications            = 2,
    n_common_treatments            = 10
  )
  expect_false(res$feasible)
  expect_equal(res$difference, 20L)
})

test_that("heterogeneous environment capacities handled", {
  res <- check_balanced_incomplete_feasibility(
    n_treatments_total             = 100,
    n_environments                 = 4,
    n_test_entries_per_environment = c(40, 45, 40, 45),
    target_replications            = 3,
    n_common_treatments            = 5
  )
  expect_equal(length(res$k_sparse), 4L)
  expect_equal(res$k_sparse, c(35L, 40L, 35L, 40L))
  expect_equal(res$total_sparse_slots, 150L)
})

test_that("no common treatments: identity reduces to standard BIBD", {
  res <- check_balanced_incomplete_feasibility(
    n_treatments_total             = 12,
    n_environments                 = 4,
    n_test_entries_per_environment = 3,
    target_replications            = 1,
    n_common_treatments            = 0
  )
  expect_true(res$feasible)
  expect_equal(res$n_sparse_treatments, 12L)
})

test_that("message field is a character string", {
  res <- check_balanced_incomplete_feasibility(
    n_treatments_total = 10, n_environments = 2,
    n_test_entries_per_environment = 5, target_replications = 1
  )
  expect_type(res$message, "character")
  expect_gt(nchar(res$message), 0)
})

test_that("stops when environment capacity less than n_common_treatments", {
  expect_error(
    check_balanced_incomplete_feasibility(
      n_treatments_total             = 20,
      n_environments                 = 3,
      n_test_entries_per_environment = 3,
      target_replications            = 2,
      n_common_treatments            = 5
    )
  )
})

test_that("stops when n_common_treatments exceeds n_treatments_total", {
  expect_error(
    check_balanced_incomplete_feasibility(
      n_treatments_total = 10, n_environments = 3,
      n_test_entries_per_environment = 5, target_replications = 1,
      n_common_treatments = 15
    )
  )
})
