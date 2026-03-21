# Tests for allocate_sparse_met() and check_balanced_incomplete_feasibility()
#
# Each test block is self-contained: the helper data is reconstructed at the
# top of each test so that failures are isolated and reproducible without
# depending on shared state.
#
# Naming convention: "function — what is being tested — expected outcome"

# ============================================================
# allocate_sparse_met(): output structure
# ============================================================

test_that("allocate_sparse_met returns a correctly structured list for random_balanced", {
  x <- make_example_sparsemet_data()
  
  out <- allocate_sparse_met(
    treatments                     = x$treatments,
    environments                   = x$environments,
    allocation_method              = "random_balanced",
    n_test_entries_per_environment = 30,
    target_replications            = 1,
    seed                           = 123
  )
  
  expect_true(is.list(out))
  expect_true(is.matrix(out$allocation_matrix))
  expect_equal(nrow(out$allocation_matrix), length(x$treatments))
  expect_equal(ncol(out$allocation_matrix), length(x$environments))
  expect_true(all(out$allocation_matrix %in% c(0L, 1L)))
  
  expect_true(is.data.frame(out$allocation_long))
  expect_equal(
    nrow(out$allocation_long),
    length(x$treatments) * length(x$environments)
  )
  
  expect_equal(length(out$environment_sizes), length(x$environments))
  expect_true(all(out$environment_sizes == 30))
})

# ============================================================
# allocate_sparse_met(): M3 / M4 aliases
# ============================================================

test_that("allocate_sparse_met translates M3 alias to random_balanced internally", {
  x <- make_example_sparsemet_data()
  
  out <- allocate_sparse_met(
    treatments                     = x$treatments,
    environments                   = x$environments,
    allocation_method              = "M3",
    n_test_entries_per_environment = 30,
    target_replications            = 1,
    seed                           = 123
  )
  
  expect_true(is.list(out))
  expect_equal(out$summary$allocation_method, "random_balanced")
})

test_that("allocate_sparse_met translates M4 alias to balanced_incomplete internally", {
  x <- make_example_sparsemet_data()
  
  # 75 sparse treatments x 1 rep = 3 environments x 25 sparse slots:
  # exact feasibility holds, so allow_approximate = FALSE is valid here.
  out <- allocate_sparse_met(
    treatments                     = x$treatments,
    environments                   = x$environments,
    allocation_method              = "M4",
    n_test_entries_per_environment = 30,
    target_replications            = 1,
    common_treatments              = x$common_treatments,
    allow_approximate              = FALSE,
    seed                           = 123
  )
  
  expect_true(is.list(out))
  expect_equal(out$summary$allocation_method, "balanced_incomplete")
})

# ============================================================
# allocate_sparse_met(): common treatments
# ============================================================

test_that("allocate_sparse_met assigns common treatments to every environment", {
  x <- make_example_sparsemet_data()
  
  out <- allocate_sparse_met(
    treatments                     = x$treatments,
    environments                   = x$environments,
    allocation_method              = "random_balanced",
    n_test_entries_per_environment = 30,
    target_replications            = 1,
    common_treatments              = x$common_treatments,
    seed                           = 123
  )
  
  mat <- out$allocation_matrix
  
  common_rows <- mat[x$common_treatments, , drop = FALSE]
  expect_true(
    all(rowSums(common_rows) == length(x$environments)),
    info = "Each common treatment must be assigned to every environment"
  )
  
  common_flags <- out$allocation_long[
    out$allocation_long$Treatment %in% x$common_treatments &
      out$allocation_long$Assigned == 1,
    "IsCommonTreatment"
  ]
  expect_true(
    all(common_flags),
    info = "IsCommonTreatment flag must be TRUE for all common treatment assignments"
  )
})

# ============================================================
# allocate_sparse_met(): replication
# ============================================================

test_that("allocate_sparse_met respects minimum coverage for random_balanced", {
  x <- make_example_sparsemet_data()
  
  out <- allocate_sparse_met(
    treatments                     = x$treatments,
    environments                   = x$environments,
    allocation_method              = "random_balanced",
    n_test_entries_per_environment = 30,
    target_replications            = 1,
    seed                           = 123
  )
  
  reps <- out$line_replications[
    !names(out$line_replications) %in% x$common_treatments
  ]
  
  expect_true(
    all(reps >= 1),
    info = "Every non-common treatment must appear in at least one environment"
  )
})

# ============================================================
# allocate_sparse_met(): balanced_incomplete — exact feasibility error
# ============================================================

test_that("allocate_sparse_met stops when exact balanced_incomplete is infeasible and allow_approximate = FALSE", {
  x <- make_example_sparsemet_data()
  
  # 80 treatments, 3 environments, k = 17 per environment, r = 2:
  # required slots = 80 x 2 = 160; available slots = 3 x 17 = 51.
  # Also impossible to assign all treatments even once, since 51 < 80.
  expect_error(
    allocate_sparse_met(
      treatments                     = x$treatments,
      environments                   = x$environments,
      allocation_method              = "balanced_incomplete",
      n_test_entries_per_environment = 17,
      target_replications            = 2,
      allow_approximate              = FALSE,
      seed                           = 123
    ),
    info = "Should error when exact balance is impossible and total capacity cannot cover all treatments"
  )
})

# ============================================================
# allocate_sparse_met(): balanced_incomplete — approximate fallback
# ============================================================

test_that("allocate_sparse_met completes with approximate balanced_incomplete when all treatments can still be assigned", {
  x <- make_example_sparsemet_data()
  
  # 80 treatments, 3 environments, k = 30 per environment:
  # available slots = 90, so all 80 treatments can be assigned at least once.
  # But target_replications = 2 would require 160 slots, so exact balance is
  # impossible. Approximate allocation should therefore succeed.
  out <- allocate_sparse_met(
    treatments                     = x$treatments,
    environments                   = x$environments,
    allocation_method              = "balanced_incomplete",
    n_test_entries_per_environment = 30,
    target_replications            = 2,
    allow_approximate              = TRUE,
    seed                           = 123
  )
  
  expect_true(is.list(out))
  expect_true(is.matrix(out$allocation_matrix))
  expect_equal(nrow(out$allocation_matrix), length(x$treatments))
  expect_equal(ncol(out$allocation_matrix), length(x$environments))
  
  expect_true(all(out$environment_sizes <= 30))
  expect_true(all(out$line_replications >= 1))
})

# ============================================================
# check_balanced_incomplete_feasibility(): output structure
# ============================================================

test_that("check_balanced_incomplete_feasibility returns a correctly structured list", {
  chk <- check_balanced_incomplete_feasibility(
    n_treatments_total             = 80,
    n_environments                 = 4,
    n_test_entries_per_environment = 20,
    target_replications            = 1,
    n_common_treatments            = 0
  )
  
  expect_true(is.list(chk))
  expect_true(all(c("feasible", "message", "difference",
                    "total_sparse_slots", "required_sparse_slots",
                    "n_sparse_treatments", "k_sparse") %in% names(chk)))
  expect_true(is.logical(chk$feasible))
  expect_true(is.character(chk$message))
})

test_that("check_balanced_incomplete_feasibility correctly identifies an exact feasible case", {
  chk <- check_balanced_incomplete_feasibility(
    n_treatments_total             = 80,
    n_environments                 = 4,
    n_test_entries_per_environment = 20,
    target_replications            = 1,
    n_common_treatments            = 0
  )
  
  expect_true(chk$feasible)
  expect_equal(chk$difference, 0)
  expect_equal(chk$total_sparse_slots, 80)
  expect_equal(chk$required_sparse_slots, 80)
})

test_that("check_balanced_incomplete_feasibility correctly identifies a slot deficit", {
  chk <- check_balanced_incomplete_feasibility(
    n_treatments_total             = 80,
    n_environments                 = 4,
    n_test_entries_per_environment = 20,
    target_replications            = 2,
    n_common_treatments            = 0
  )
  
  expect_false(chk$feasible)
  expect_equal(chk$difference, -80)
})

test_that("check_balanced_incomplete_feasibility accounts correctly for common treatments", {
  chk <- check_balanced_incomplete_feasibility(
    n_treatments_total             = 80,
    n_environments                 = 4,
    n_test_entries_per_environment = 20,
    target_replications            = 1,
    n_common_treatments            = 5
  )
  
  expect_false(chk$feasible)
  expect_equal(chk$n_sparse_treatments, 75)
  expect_equal(chk$total_sparse_slots, 60)
  expect_equal(chk$difference, -15)
})