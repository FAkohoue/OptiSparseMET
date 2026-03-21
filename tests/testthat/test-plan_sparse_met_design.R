# Tests for plan_sparse_met_design() and combine_met_fieldbooks()

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
    "summary",
    "seed_used"
  )
  
  expect_equal(
    length(setdiff(expected_keys, names(out))),
    0
  )
  
  expect_true(is.list(out$environment_designs))
  expect_true(is.data.frame(out$combined_field_book))
  expect_true(is.data.frame(out$environment_summary))
})

# ============================================================
# environment coverage
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
  
  expect_equal(
    sort(names(out$environment_designs)),
    sort(x$environments)
  )
  
  expect_equal(
    sort(out$environment_summary$Environment),
    sort(x$environments)
  )
})

# ============================================================
# combined_field_book columns
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
    "SparseMethod"
  )
  
  expect_equal(
    length(setdiff(required_cols, names(fb))),
    0
  )
})

# ============================================================
# environment column coverage
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
  
  expect_true(
    setequal(
      unique(out$combined_field_book$Environment),
      x$environments
    )
  )
})

# ============================================================
# common treatment flag
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
# missing env spec error
# ============================================================

test_that("errors when env_design_specs missing environment", {
  
  x <- make_example_sparsemet_data()
  
  incomplete_specs <-
    x$env_design_specs[
      setdiff(x$environments, x$environments[1])
    ]
  
  expect_error(
    
    plan_sparse_met_design(
      treatments                     = x$treatments,
      environments                   = x$environments,
      allocation_method              = "random_balanced",
      n_test_entries_per_environment = 30,
      target_replications            = 1,
      env_design_specs               = incomplete_specs,
      seed                           = 123
    )
    
  )
})

# ============================================================
# combine_met_fieldbooks tests remain unchanged
# ============================================================