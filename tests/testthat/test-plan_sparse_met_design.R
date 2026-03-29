# ============================================================
# Tests for plan_sparse_met_design() and combine_met_fieldbooks()
# Compatible with automatic replication detection + efficiency capture
# ============================================================

# ============================================================
# plan_sparse_met_design(): output structure
# ============================================================

test_that("plan_sparse_met_design returns a correctly structured list", {

  x <- make_example_sparsemet_data()

  out <- plan_sparse_met_design(
    treatments                     = x$treatments,
    environments                   = x$environments,
    allocation_method              = "random_balanced",
    n_test_entries_per_environment = 30,
    target_replications            = 1,
    common_treatments              = x$common_treatments,
    env_design_specs               = x$env_design_specs,
    treatment_info                 = x$treatment_info,
    seed_info                      = x$seed_info,
    seed_required_per_plot         = x$seed_required_per_plot,
    seed                           = 123
  )

  expect_true(is.list(out))

  expected_keys <- c(
    "sparse_allocation",
    "environment_designs",
    "combined_field_book",
    "environment_summary",
    "group_environment_summary",
    "efficiency_summary",
    "summary",
    "seed_used"
  )

  expect_setequal(names(out), expected_keys)

  expect_true(is.list(out$environment_designs))
  expect_true(is.data.frame(out$combined_field_book))
  expect_true(is.data.frame(out$environment_summary))
  expect_true(is.data.frame(out$efficiency_summary))
})


# ============================================================
# plan_sparse_met_design(): environment coverage
# ============================================================

test_that("plan_sparse_met_design produces a design for every environment", {

  x <- make_example_sparsemet_data()

  out <- plan_sparse_met_design(
    treatments                     = x$treatments,
    environments                   = x$environments,
    allocation_method              = "random_balanced",
    n_test_entries_per_environment = 30,
    target_replications            = 1,
    common_treatments              = x$common_treatments,
    env_design_specs               = x$env_design_specs,
    treatment_info                 = x$treatment_info,
    seed                           = 123
  )

  expect_setequal(
    names(out$environment_designs),
    x$environments
  )

  expect_setequal(
    out$environment_summary$Environment,
    x$environments
  )
})


# ============================================================
# plan_sparse_met_design(): combined_field_book columns
# ============================================================

test_that("combined_field_book contains required MET-level columns", {

  x <- make_example_sparsemet_data()

  out <- plan_sparse_met_design(
    treatments                     = x$treatments,
    environments                   = x$environments,
    allocation_method              = "random_balanced",
    n_test_entries_per_environment = 30,
    target_replications            = 1,
    common_treatments              = x$common_treatments,
    env_design_specs               = x$env_design_specs,
    treatment_info                 = x$treatment_info,
    seed_info                      = x$seed_info,
    seed_required_per_plot         = x$seed_required_per_plot,
    seed                           = 123
  )

  fb <- out$combined_field_book

  required_cols <- c(
    "Environment",
    "Treatment",
    "Family",
    "Gcluster",
    "Block",
    "Plot",
    "Row",
    "Column",
    "IsCommonTreatment",
    "LocalDesign",
    "ReplicationMode",
    "SparseMethod"
  )

  expect_true(all(required_cols %in% names(fb)))
})


# ============================================================
# combined_field_book covers all environments
# ============================================================

test_that("combined_field_book covers all environments", {

  x <- make_example_sparsemet_data()

  out <- plan_sparse_met_design(
    treatments                     = x$treatments,
    environments                   = x$environments,
    allocation_method              = "random_balanced",
    n_test_entries_per_environment = 30,
    target_replications            = 1,
    common_treatments              = x$common_treatments,
    env_design_specs               = x$env_design_specs,
    treatment_info                 = x$treatment_info,
    seed                           = 123
  )

  expect_setequal(
    unique(out$combined_field_book$Environment),
    x$environments
  )
})


# ============================================================
# common treatment flag correctness
# ============================================================

test_that("IsCommonTreatment flag is correct", {

  x <- make_example_sparsemet_data()

  out <- plan_sparse_met_design(
    treatments                     = x$treatments,
    environments                   = x$environments,
    allocation_method              = "random_balanced",
    n_test_entries_per_environment = 30,
    target_replications            = 1,
    common_treatments              = x$common_treatments,
    env_design_specs               = x$env_design_specs,
    treatment_info                 = x$treatment_info,
    seed                           = 123
  )

  fb <- out$combined_field_book

  expect_true(
    all(
      fb$IsCommonTreatment[
        fb$Treatment %in% x$common_treatments
      ]
    )
  )

  expect_true(
    all(
      !fb$IsCommonTreatment[
        !fb$Treatment %in% x$common_treatments
      ]
    )
  )
})


# ============================================================
# efficiency summary structure
# ============================================================

test_that("efficiency_summary has expected structure", {

  x <- make_example_sparsemet_data()

  out <- plan_sparse_met_design(
    treatments                     = x$treatments,
    environments                   = x$environments,
    allocation_method              = "random_balanced",
    n_test_entries_per_environment = 30,
    target_replications            = 1,
    common_treatments              = x$common_treatments,
    env_design_specs               = x$env_design_specs,
    treatment_info                 = x$treatment_info,
    seed                           = 123
  )

  eff <- out$efficiency_summary

  expect_true(is.data.frame(eff))

  required_cols <- c(
    "Environment",
    "LocalDesign",
    "Metric",
    "Value",
    "ValueType"
  )

  expect_true(all(required_cols %in% names(eff)))
})


# ============================================================
# efficiency columns in environment_summary
# ============================================================

test_that("environment_summary contains efficiency columns", {

  x <- make_example_sparsemet_data()

  out <- plan_sparse_met_design(
    treatments                     = x$treatments,
    environments                   = x$environments,
    allocation_method              = "random_balanced",
    n_test_entries_per_environment = 30,
    target_replications            = 1,
    common_treatments              = x$common_treatments,
    env_design_specs               = x$env_design_specs,
    treatment_info                 = x$treatment_info,
    seed                           = 123
  )

  env_sum <- out$environment_summary

  expect_true(all(
    c("has_efficiency", "eff_A", "eff_D", "eff_mean_PEV") %in% names(env_sum)
  ))
})


# ============================================================
# missing environment specification error
# ============================================================

test_that("error if env_design_specs missing environment", {

  x <- make_example_sparsemet_data()

  bad_specs <- x$env_design_specs[-1]

  expect_error(
    plan_sparse_met_design(
      treatments       = x$treatments,
      environments     = x$environments,
      env_design_specs = bad_specs
    )
  )
})


# ============================================================
# combine_met_fieldbooks: basic structure
# ============================================================

test_that("combine_met_fieldbooks returns expected structure", {

  fb_E1 <- data.frame(
    Treatment = c("L001", "L002", "CHK1"),
    Family    = c("F1", "F2", "CHECK"),
    Block     = 1L,
    Plot      = 1:3,
    Row       = 1L,
    Column    = 1:3
  )

  fb_E2 <- data.frame(
    Treatment = c("L003", "L004", "CHK1"),
    Family    = c("F3", "F4", "CHECK"),
    Block     = 1L,
    Plot      = 1:3,
    Row       = 1L,
    Column    = 1:3
  )

  out <- combine_met_fieldbooks(
    field_books       = list(E1 = fb_E1, E2 = fb_E2),
    local_designs     = c(E1 = "met_prep_famoptg", E2 = "met_prep_famoptg"),
    replication_modes = c(E1 = "augmented",        E2 = "augmented"),
    sparse_method     = "balanced_incomplete",
    common_treatments = "CHK1"
  )

  expect_true(is.data.frame(out))
  expect_equal(nrow(out), 6)
  expect_setequal(unique(out$Environment), c("E1", "E2"))
})


# ============================================================
# combine_met_fieldbooks: row names reset correctly
# ============================================================

test_that("combine_met_fieldbooks resets row names", {

  fb_E1 <- data.frame(
    Treatment = c("L001", "CHK1"),
    Family    = c("F1", "CHECK"),
    Block     = 1L,
    Plot      = 1:2,
    Row       = 1L,
    Column    = 1:2,
    row.names = c("a", "b")
  )

  fb_E2 <- data.frame(
    Treatment = c("L002", "CHK1"),
    Family    = c("F2", "CHECK"),
    Block     = 1L,
    Plot      = 1:2,
    Row       = 1L,
    Column    = 1:2,
    row.names = c("c", "d")
  )

  out <- combine_met_fieldbooks(
    field_books = list(E1 = fb_E1, E2 = fb_E2)
  )

  expect_equal(
    rownames(out),
    as.character(seq_len(nrow(out)))
  )
})


# ============================================================
# combine_met_fieldbooks: heterogeneous columns preserved
# ============================================================

test_that("combine_met_fieldbooks preserves heterogeneous columns", {

  fb_E1 <- data.frame(
    Treatment = "L001",
    Family    = "F1",
    Plot      = 1
  )

  fb_E2 <- data.frame(
    Treatment       = "L002",
    Family          = "F2",
    Plot            = 1,
    SpatialResidual = 0.2
  )

  out <- combine_met_fieldbooks(
    field_books = list(E1 = fb_E1, E2 = fb_E2)
  )

  expect_true("SpatialResidual" %in% names(out))
  expect_true(is.na(out$SpatialResidual[out$Environment == "E1"]))
})
