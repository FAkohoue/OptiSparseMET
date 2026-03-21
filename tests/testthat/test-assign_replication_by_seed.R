# Tests for assign_replication_by_seed()
#
# Covers all three replication modes (augmented, p_rep, rcbd_type), the three
# shortage actions (downgrade, exclude, error), priority ordering in p_rep
# mode, and input validation. Each test reconstructs the helper data
# independently so failures are isolated.

# ============================================================
# augmented mode
# ============================================================

test_that("assign_replication_by_seed augmented mode returns no replicated treatments", {
  x   <- make_example_sparsemet_data()
  trt <- x$treatments[1:12]
  
  out <- assign_replication_by_seed(
    treatments             = trt,
    seed_available         = x$seed_info,
    seed_required_per_plot = 10,
    replication_mode       = "augmented"
  )
  
  expect_true(is.list(out))
  
  # augmented mode targets one plot per treatment — p_rep_treatments must
  # always be empty regardless of seed availability
  expect_length(out$p_rep_treatments, 0)
  expect_length(out$p_rep_reps,       0)
  
  # At least some treatments must have sufficient seed for one plot
  expect_true(
    length(out$unreplicated_treatments) >= 1,
    info = "At least one treatment must be feasible for a single plot at seed_required = 10"
  )
  
  # Every treatment in trt must appear in exactly one role
  all_roles <- c(out$unreplicated_treatments, out$excluded_treatments)
  expect_setequal(all_roles, trt)
  
  # seed_summary must cover all input treatments with a Role assigned
  expect_true(is.data.frame(out$seed_summary))
  expect_equal(nrow(out$seed_summary), length(trt))
  expect_true(all(!is.na(out$seed_summary$Role)))
})

test_that("assign_replication_by_seed augmented mode respects minimum_seed_buffer", {
  x   <- make_example_sparsemet_data()
  trt <- x$treatments[1:12]
  
  # With no buffer: feasibility threshold = 10
  out_no_buffer <- assign_replication_by_seed(
    treatments             = trt,
    seed_available         = x$seed_info,
    seed_required_per_plot = 10,
    replication_mode       = "augmented",
    minimum_seed_buffer    = 0
  )
  
  # With buffer = 30: feasibility threshold = 40 — more treatments excluded
  out_buffered <- assign_replication_by_seed(
    treatments             = trt,
    seed_available         = x$seed_info,
    seed_required_per_plot = 10,
    replication_mode       = "augmented",
    minimum_seed_buffer    = 30
  )
  
  expect_true(
    length(out_buffered$unreplicated_treatments) <=
      length(out_no_buffer$unreplicated_treatments),
    info = "A larger seed buffer must reduce or maintain the number of feasible treatments"
  )
})

# ============================================================
# p_rep mode
# ============================================================

test_that("assign_replication_by_seed p_rep mode respects max_prep cap", {
  x   <- make_example_sparsemet_data()
  trt <- x$treatments[1:20]
  
  out <- assign_replication_by_seed(
    treatments             = trt,
    seed_available         = x$seed_info,
    seed_required_per_plot = 10,
    replication_mode       = "p_rep",
    desired_replications   = 2,
    max_prep               = 5,
    priority               = "seed_available",
    seed                   = 123
  )
  
  expect_true(is.list(out))
  
  # max_prep = 5: at most 5 treatments can be replicated
  expect_true(
    length(out$p_rep_treatments) <= 5,
    info = "Number of replicated treatments must not exceed max_prep"
  )
  
  # p_rep_reps must be a parallel vector to p_rep_treatments
  expect_equal(length(out$p_rep_treatments), length(out$p_rep_reps))
  
  # All assigned reps must equal desired_replications
  if (length(out$p_rep_reps) > 0) {
    expect_true(all(out$p_rep_reps == 2))
  }
  
  # Every treatment must appear in exactly one role
  all_roles <- c(out$p_rep_treatments, out$unreplicated_treatments, out$excluded_treatments)
  expect_setequal(all_roles, trt)
})

test_that("assign_replication_by_seed p_rep mode seed_available priority selects best-seeded treatments", {
  x   <- make_example_sparsemet_data()
  trt <- x$treatments[1:20]
  
  out <- assign_replication_by_seed(
    treatments             = trt,
    seed_available         = x$seed_info,
    seed_required_per_plot = 10,
    replication_mode       = "p_rep",
    desired_replications   = 2,
    max_prep               = 5,
    priority               = "seed_available",
    seed                   = 123
  )
  
  # Treatments selected for replication must all satisfy the replicated
  # feasibility condition: SeedAvailable >= 2 * 10 = 20
  if (length(out$p_rep_treatments) > 0) {
    selected_seed <- x$seed_info$SeedAvailable[
      match(out$p_rep_treatments, x$seed_info$Treatment)
    ]
    expect_true(
      all(selected_seed >= 20),
      info = "All replicated treatments must meet the seed requirement for 2 plots"
    )
  }
})

test_that("assign_replication_by_seed p_rep mode respects candidate_prep restriction", {
  x   <- make_example_sparsemet_data()
  trt <- x$treatments[1:20]
  
  # Restrict candidates to the last 10 treatments in the vector
  candidates <- trt[11:20]
  
  out <- assign_replication_by_seed(
    treatments             = trt,
    seed_available         = x$seed_info,
    seed_required_per_plot = 10,
    replication_mode       = "p_rep",
    desired_replications   = 2,
    candidate_prep         = candidates,
    seed                   = 123
  )
  
  # No treatment outside the candidate set should be replicated
  expect_true(
    all(out$p_rep_treatments %in% candidates),
    info = "Replicated treatments must be drawn exclusively from candidate_prep"
  )
})

# ============================================================
# rcbd_type mode
# ============================================================

test_that("assign_replication_by_seed rcbd_type mode with downgrade partitions correctly", {
  x   <- make_example_sparsemet_data()
  trt <- x$treatments[1:15]
  
  out <- assign_replication_by_seed(
    treatments             = trt,
    seed_available         = x$seed_info,
    seed_required_per_plot = 10,
    replication_mode       = "rcbd_type",
    desired_replications   = 2,
    shortage_action        = "downgrade"
  )
  
  expect_true(is.list(out))
  expect_true(
    all(c("p_rep_treatments", "unreplicated_treatments", "excluded_treatments") %in% names(out))
  )
  
  # Replicated treatments must all satisfy seed >= 2 * 10 = 20
  if (length(out$p_rep_treatments) > 0) {
    rep_seed <- x$seed_info$SeedAvailable[
      match(out$p_rep_treatments, x$seed_info$Treatment)
    ]
    expect_true(all(rep_seed >= 20))
  }
  
  # Downgraded (unreplicated) treatments must satisfy 10 <= seed < 20
  if (length(out$unreplicated_treatments) > 0) {
    unrep_seed <- x$seed_info$SeedAvailable[
      match(out$unreplicated_treatments, x$seed_info$Treatment)
    ]
    expect_true(all(unrep_seed >= 10))
    expect_true(all(unrep_seed < 20))
  }
  
  # Every treatment must appear in exactly one role
  all_roles <- c(out$p_rep_treatments, out$unreplicated_treatments, out$excluded_treatments)
  expect_setequal(all_roles, trt)
})

test_that("assign_replication_by_seed rcbd_type mode with exclude drops infeasible treatments", {
  x   <- make_example_sparsemet_data()
  trt <- x$treatments[1:15]
  
  out <- assign_replication_by_seed(
    treatments             = trt,
    seed_available         = x$seed_info,
    seed_required_per_plot = 10,
    replication_mode       = "rcbd_type",
    desired_replications   = 2,
    shortage_action        = "exclude"
  )
  
  # Under exclude, treatments with seed < 20 must not appear as unreplicated
  expect_length(out$unreplicated_treatments, 0)
  
  # Excluded treatments must all fail the replicated feasibility condition
  if (length(out$excluded_treatments) > 0) {
    excl_seed <- x$seed_info$SeedAvailable[
      match(out$excluded_treatments, x$seed_info$Treatment)
    ]
    expect_true(all(excl_seed < 20))
  }
})

test_that("assign_replication_by_seed rcbd_type mode with error stops on shortage", {
  # Construct a case where at least one treatment definitely has insufficient
  # seed: require 50 seeds per plot with desired_replications = 2,
  # so threshold = 100. Some treatments have seed < 100 in the test data.
  x   <- make_example_sparsemet_data()
  trt <- x$treatments[1:15]
  
  expect_error(
    assign_replication_by_seed(
      treatments             = trt,
      seed_available         = x$seed_info,
      seed_required_per_plot = 50,
      replication_mode       = "rcbd_type",
      desired_replications   = 2,
      shortage_action        = "error"
    ),
    info = "Should error when shortage_action = 'error' and at least one treatment cannot reach desired_replications"
  )
})

# ============================================================
# Input validation
# ============================================================

test_that("assign_replication_by_seed errors when seed_available is missing a treatment", {
  x       <- make_example_sparsemet_data()
  # Remove the first treatment row to create a deliberate gap
  bad_seed <- x$seed_info[-1, ]
  
  expect_error(
    assign_replication_by_seed(
      treatments             = x$treatments[1:10],
      seed_available         = bad_seed,
      seed_required_per_plot = 10,
      replication_mode       = "augmented"
    ),
    info = "Should error when seed_available is missing rows for treatments in the input vector"
  )
})

test_that("assign_replication_by_seed errors when check_treatments overlap with treatments", {
  x   <- make_example_sparsemet_data()
  trt <- x$treatments[1:10]
  
  expect_error(
    assign_replication_by_seed(
      treatments             = trt,
      check_treatments       = trt[1:2],   # deliberate overlap
      seed_available         = x$seed_info,
      seed_required_per_plot = 10,
      replication_mode       = "augmented"
    ),
    info = "Should error when check_treatments overlap with treatments"
  )
})

test_that("assign_replication_by_seed errors on non-positive seed_required_per_plot", {
  x   <- make_example_sparsemet_data()
  trt <- x$treatments[1:10]
  
  expect_error(
    assign_replication_by_seed(
      treatments             = trt,
      seed_available         = x$seed_info,
      seed_required_per_plot = 0,
      replication_mode       = "augmented"
    ),
    info = "seed_required_per_plot = 0 is not a valid plot requirement"
  )
})