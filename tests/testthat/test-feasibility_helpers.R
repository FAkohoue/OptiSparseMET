# Tests for feasibility_helpers.R
#
# Covers:
# - min_k_for_full_coverage()
# - suggest_safe_k()
# - warn_if_k_too_small()
# - .check_full_coverage_feasibility()

# ============================================================
# min_k_for_full_coverage()
# ============================================================

test_that("min_k_for_full_coverage returns correct minimum with no common treatments", {
  out <- min_k_for_full_coverage(
    n_treatments_total  = 80,
    n_environments      = 3,
    n_common_treatments = 0
  )

  expect_true(is.list(out))
  expect_equal(out$n_sparse_treatments, 80)
  expect_equal(out$min_sparse_slots_per_environment, 27)
  expect_equal(out$min_total_entries_per_environment, 27)
})

test_that("min_k_for_full_coverage returns correct minimum with common treatments", {
  out <- min_k_for_full_coverage(
    n_treatments_total  = 80,
    n_environments      = 3,
    n_common_treatments = 5
  )

  expect_true(is.list(out))
  expect_equal(out$n_sparse_treatments, 75)
  expect_equal(out$min_sparse_slots_per_environment, 25)
  expect_equal(out$min_total_entries_per_environment, 30)
})

test_that("min_k_for_full_coverage errors on invalid inputs", {
  expect_error(
    min_k_for_full_coverage(
      n_treatments_total  = 0,
      n_environments      = 3,
      n_common_treatments = 0
    )
  )

  expect_error(
    min_k_for_full_coverage(
      n_treatments_total  = 80,
      n_environments      = 0,
      n_common_treatments = 0
    )
  )

  expect_error(
    min_k_for_full_coverage(
      n_treatments_total  = 80,
      n_environments      = 3,
      n_common_treatments = 81
    )
  )
})

# ============================================================
# suggest_safe_k()
# ============================================================

test_that("suggest_safe_k returns minimum plus buffer", {
  trt <- paste0("L", sprintf("%03d", 1:80))
  env <- c("Env1", "Env2", "Env3")

  out <- suggest_safe_k(
    treatments        = trt,
    environments      = env,
    common_treatments = NULL,
    buffer            = 3
  )

  # minimum is ceil(80/3) = 27, plus buffer 3 -> 30
  expect_equal(out, 30)
})

test_that("suggest_safe_k accounts for common treatments", {
  trt    <- paste0("L", sprintf("%03d", 1:80))
  env    <- c("Env1", "Env2", "Env3")
  common <- trt[1:5]

  out <- suggest_safe_k(
    treatments        = trt,
    environments      = env,
    common_treatments = common,
    buffer            = 2
  )

  # sparse = 75, ceil(75/3) = 25, plus 5 common = 30, plus buffer 2 -> 32
  expect_equal(out, 32)
})

test_that("suggest_safe_k errors on invalid buffer", {
  trt <- paste0("L", sprintf("%03d", 1:80))
  env <- c("Env1", "Env2", "Env3")

  expect_error(
    suggest_safe_k(
      treatments   = trt,
      environments = env,
      buffer       = -1
    )
  )
})

# ============================================================
# warn_if_k_too_small()
# ============================================================

test_that("warn_if_k_too_small warns when capacity is insufficient", {
  trt <- paste0("L", sprintf("%03d", 1:80))
  env <- c("Env1", "Env2", "Env3")

  expect_warning(
    warn_if_k_too_small(
      treatments                     = trt,
      environments                   = env,
      n_test_entries_per_environment = 20
    )
  )
})

test_that("warn_if_k_too_small does not warn when capacity is sufficient", {
  trt <- paste0("L", sprintf("%03d", 1:80))
  env <- c("Env1", "Env2", "Env3")

  expect_no_warning(
    warn_if_k_too_small(
      treatments                     = trt,
      environments                   = env,
      n_test_entries_per_environment = 30
    )
  )
})

test_that("warn_if_k_too_small accounts for common treatments", {
  trt    <- paste0("L", sprintf("%03d", 1:80))
  env    <- c("Env1", "Env2", "Env3")
  common <- trt[1:5]

  # 29 total entries per env: 29 - 5 = 24 sparse slots per env
  # total sparse slots = 72 < 75 sparse treatments -> should warn
  expect_warning(
    warn_if_k_too_small(
      treatments                     = trt,
      environments                   = env,
      n_test_entries_per_environment = 29,
      common_treatments              = common
    )
  )
})

# ============================================================
# .check_full_coverage_feasibility()
# ============================================================

test_that(".check_full_coverage_feasibility passes when capacity is sufficient", {
  trt <- paste0("L", sprintf("%03d", 1:80))
  env <- c("Env1", "Env2", "Env3")

  expect_invisible(
    .check_full_coverage_feasibility(
      treatments                     = trt,
      environments                   = env,
      n_test_entries_per_environment = 30,
      common_treatments              = NULL
    )
  )
})

test_that(".check_full_coverage_feasibility errors when capacity is insufficient", {
  trt <- paste0("L", sprintf("%03d", 1:80))
  env <- c("Env1", "Env2", "Env3")

  expect_error(
    .check_full_coverage_feasibility(
      treatments                     = trt,
      environments                   = env,
      n_test_entries_per_environment = 20,
      common_treatments              = NULL
    )
  )
})

test_that(".check_full_coverage_feasibility accounts for common treatments correctly", {
  trt    <- paste0("L", sprintf("%03d", 1:80))
  env    <- c("Env1", "Env2", "Env3")
  common <- trt[1:5]

  # 29 total -> 24 sparse slots per env -> 72 total < 75 sparse -> error
  expect_error(
    .check_full_coverage_feasibility(
      treatments                     = trt,
      environments                   = env,
      n_test_entries_per_environment = 29,
      common_treatments              = common
    )
  )

  # 30 total -> 25 sparse slots per env -> 75 total >= 75 sparse -> OK
  expect_invisible(
    .check_full_coverage_feasibility(
      treatments                     = trt,
      environments                   = env,
      n_test_entries_per_environment = 30,
      common_treatments              = common
    )
  )
})

test_that(".check_full_coverage_feasibility accepts vector k values", {
  trt <- paste0("L", sprintf("%03d", 1:80))
  env <- c("Env1", "Env2", "Env3")

  # 25 + 25 + 30 = 80 sparse slots >= 80 treatments -> exactly feasible
  expect_invisible(
    .check_full_coverage_feasibility(
      treatments                     = trt,
      environments                   = env,
      n_test_entries_per_environment = c(25, 25, 30),
      common_treatments              = NULL
    )
  )
})

test_that(".check_full_coverage_feasibility errors when k vector length is invalid", {
  trt <- paste0("L", sprintf("%03d", 1:80))
  env <- c("Env1", "Env2", "Env3")

  expect_error(
    .check_full_coverage_feasibility(
      treatments                     = trt,
      environments                   = env,
      n_test_entries_per_environment = c(30, 30),   # length 2, need 3
      common_treatments              = NULL
    )
  )
})

test_that(".check_full_coverage_feasibility errors when an environment has fewer slots than common treatments", {
  trt    <- paste0("L", sprintf("%03d", 1:20))
  env    <- c("Env1", "Env2", "Env3")
  common <- trt[1:5]

  # Env2 has capacity 4 < 5 common treatments -> error
  expect_error(
    .check_full_coverage_feasibility(
      treatments                     = trt,
      environments                   = env,
      n_test_entries_per_environment = c(5, 4, 6),
      common_treatments              = common
    )
  )
})
