# Tests for OptiSparseMET_example_data
#
# Verifies that the bundled example dataset loads correctly and that all
# components have the structure required to pass them directly to package
# functions. Tests are ordered from coarse (does it load, does it have the
# right keys) to fine (are internal IDs consistent across components).

# ============================================================
# Loading and top-level structure
# ============================================================

test_that("OptiSparseMET_example_data loads and contains all expected components", {
  data("OptiSparseMET_example_data", package = "OptiSparseMET")

  expect_true(exists("OptiSparseMET_example_data"))
  expect_true(is.list(OptiSparseMET_example_data))

  expected_names <- c(
    "treatments",
    "environments",
    "treatment_info",
    "common_treatments",
    "seed_info",
    "seed_required_per_plot",
    "OptiSparseMET_GRM",
    "OptiSparseMET_A",
    "OptiSparseMET_K",
    "sparse_example_args_random_balanced",
    "sparse_example_args_balanced_incomplete",
    "env_design_specs"
  )

  missing <- setdiff(expected_names, names(OptiSparseMET_example_data))
  expect_equal(
    length(missing), 0,
    info = paste0("Missing components: ", paste(missing, collapse = ", "))
  )
})

# ============================================================
# Core vectors: treatments and environments
# ============================================================

test_that("OptiSparseMET_example_data treatments and environments are valid character vectors", {
  data("OptiSparseMET_example_data", package = "OptiSparseMET")
  x <- OptiSparseMET_example_data

  expect_true(is.character(x$treatments))
  expect_true(length(x$treatments) > 0)
  expect_true(all(!is.na(x$treatments)))
  expect_equal(length(x$treatments), length(unique(x$treatments)),
               info = "Treatment IDs must be unique")

  expect_true(is.character(x$environments))
  expect_true(length(x$environments) >= 2,
              info = "At least two environments are required for a MET")
  expect_true(all(!is.na(x$environments)))
  expect_equal(length(x$environments), length(unique(x$environments)),
               info = "Environment names must be unique")
})

test_that("OptiSparseMET_example_data common_treatments is a subset of treatments", {
  data("OptiSparseMET_example_data", package = "OptiSparseMET")
  x <- OptiSparseMET_example_data

  expect_true(is.character(x$common_treatments))
  expect_true(length(x$common_treatments) >= 1)
  expect_true(
    all(x$common_treatments %in% x$treatments),
    info = "Every common treatment must appear in the full treatments vector"
  )
})

# ============================================================
# treatment_info
# ============================================================

test_that("OptiSparseMET_example_data treatment_info has required columns and covers all treatments", {
  data("OptiSparseMET_example_data", package = "OptiSparseMET")
  x <- OptiSparseMET_example_data

  expect_true(is.data.frame(x$treatment_info))
  expect_true(all(c("Treatment", "Family") %in% names(x$treatment_info)),
              info = "treatment_info must have Treatment and Family columns")

  missing_trt <- setdiff(x$treatments, x$treatment_info$Treatment)
  expect_equal(
    length(missing_trt), 0,
    info = paste0("Treatments missing from treatment_info: ",
                  paste(utils::head(missing_trt, 5), collapse = ", "))
  )

  expect_true(all(!is.na(x$treatment_info$Family)))
  expect_true(all(nchar(x$treatment_info$Family) > 0))
})

# ============================================================
# seed_info and seed_required_per_plot
# ============================================================

test_that("OptiSparseMET_example_data seed_info covers all treatments with positive values", {
  data("OptiSparseMET_example_data", package = "OptiSparseMET")
  x <- OptiSparseMET_example_data

  expect_true(is.data.frame(x$seed_info))
  expect_true(all(c("Treatment", "SeedAvailable") %in% names(x$seed_info)),
              info = "seed_info must have Treatment and SeedAvailable columns")

  missing_seed <- setdiff(x$treatments, x$seed_info$Treatment)
  expect_equal(
    length(missing_seed), 0,
    info = paste0("Treatments missing from seed_info: ",
                  paste(utils::head(missing_seed, 5), collapse = ", "))
  )

  expect_true(
    all(x$seed_info$SeedAvailable > 0),
    info = "All SeedAvailable values must be positive"
  )
})

test_that("OptiSparseMET_example_data seed_required_per_plot covers all environments with positive values", {
  data("OptiSparseMET_example_data", package = "OptiSparseMET")
  x <- OptiSparseMET_example_data

  expect_true(is.data.frame(x$seed_required_per_plot))
  expect_true(
    all(c("Environment", "SeedRequiredPerPlot") %in% names(x$seed_required_per_plot)),
    info = "seed_required_per_plot must have Environment and SeedRequiredPerPlot columns"
  )

  missing_env <- setdiff(x$environments, x$seed_required_per_plot$Environment)
  expect_equal(
    length(missing_env), 0,
    info = paste0("Environments missing from seed_required_per_plot: ",
                  paste(missing_env, collapse = ", "))
  )

  expect_true(
    all(x$seed_required_per_plot$SeedRequiredPerPlot > 0),
    info = "All SeedRequiredPerPlot values must be positive"
  )
})

# ============================================================
# Relationship matrices
# ============================================================

test_that("OptiSparseMET_example_data relationship matrices are square with correct IDs", {
  data("OptiSparseMET_example_data", package = "OptiSparseMET")
  x <- OptiSparseMET_example_data

  for (mat_name in c("OptiSparseMET_GRM", "OptiSparseMET_A", "OptiSparseMET_K")) {
    mat <- x[[mat_name]]

    expect_true(is.matrix(mat),
                info = paste(mat_name, "must be a matrix"))
    expect_equal(nrow(mat), ncol(mat),
                 info = paste(mat_name, "must be square"))
    expect_equal(nrow(mat), length(x$treatments),
                 info = paste(mat_name, "must have one row/column per treatment"))
    expect_equal(rownames(mat), x$treatments,
                 info = paste(mat_name, "rownames must match treatments in order"))
    expect_equal(colnames(mat), x$treatments,
                 info = paste(mat_name, "colnames must match treatments in order"))
    expect_true(
      all(diag(mat) > 0),
      info = paste(mat_name, "diagonal must be strictly positive")
    )
  }
})

# ============================================================
# Pre-built argument lists
# ============================================================

test_that("OptiSparseMET_example_data sparse argument lists have required fields and consistent IDs", {
  data("OptiSparseMET_example_data", package = "OptiSparseMET")
  x <- OptiSparseMET_example_data

  required_fields <- c("treatments", "environments", "allocation_method",
                       "n_test_entries_per_environment", "seed")

  for (args_name in c("sparse_example_args_random_balanced",
                      "sparse_example_args_balanced_incomplete")) {
    args <- x[[args_name]]

    expect_true(is.list(args),
                info = paste(args_name, "must be a list"))
    expect_true(
      all(required_fields %in% names(args)),
      info = paste(args_name, "is missing required fields:",
                   paste(setdiff(required_fields, names(args)), collapse = ", "))
    )
    expect_equal(args$treatments,   x$treatments,
                 info = paste(args_name, "treatments must match x$treatments"))
    expect_equal(args$environments, x$environments,
                 info = paste(args_name, "environments must match x$environments"))
  }

  expect_equal(
    x$sparse_example_args_random_balanced$allocation_method,
    "random_balanced"
  )
  expect_equal(
    x$sparse_example_args_balanced_incomplete$allocation_method,
    "balanced_incomplete"
  )
})

# ============================================================
# env_design_specs
# ============================================================

test_that("OptiSparseMET_example_data env_design_specs covers all environments and has required fields", {
  data("OptiSparseMET_example_data", package = "OptiSparseMET")
  x <- OptiSparseMET_example_data

  expect_true(is.list(x$env_design_specs))

  missing_specs <- setdiff(x$environments, names(x$env_design_specs))
  expect_equal(
    length(missing_specs), 0,
    info = paste0("Environments missing from env_design_specs: ",
                  paste(missing_specs, collapse = ", "))
  )

  for (env in x$environments) {
    spec <- x$env_design_specs[[env]]

    expect_true("design" %in% names(spec),
                info = paste("env_design_specs[[", env, "]] must have a 'design' field"))
    expect_true(
      spec$design %in% c("met_prep_famoptg", "met_alpha_rc_stream"),
      info = paste("env_design_specs[[", env, "]]$design must be 'met_prep_famoptg' or 'met_alpha_rc_stream'")
    )
    expect_true("check_treatments" %in% names(spec),
                info = paste(env, "spec must include check_treatments"))
    expect_true("check_families" %in% names(spec),
                info = paste(env, "spec must include check_families"))
    expect_equal(
      length(spec$check_treatments),
      length(spec$check_families),
      info = paste(env, "check_treatments and check_families must have equal length")
    )
  }
})
