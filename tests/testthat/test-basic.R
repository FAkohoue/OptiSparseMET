test_that("package loads and all exported functions exist", {
  expect_true(isNamespaceLoaded("OptiSparseMET"))
  
  exported <- getNamespaceExports("OptiSparseMET")
  expect_true("met_prep_famoptg"                %in% exported)
  expect_true("met_evaluate_famoptg_efficiency" %in% exported)
  expect_true("met_optimize_famoptg"            %in% exported)
  expect_true("met_alpha_rc_stream"             %in% exported)
  expect_true("met_evaluate_alpha_efficiency"   %in% exported)
  expect_true("met_optimize_alpha_rc"           %in% exported)
})

test_that("internal helpers exist and are not exported", {
  # Helpers should be accessible via ::: but not via ::
  expect_true(is.function(OptiSparseMET:::.make_sparse_incidence))
  expect_true(is.function(OptiSparseMET:::.ar1_precision_sparse))
  expect_true(is.function(OptiSparseMET:::.solve_C))
  expect_true(is.function(OptiSparseMET:::.pinv_sym_dense))
  expect_true(is.function(OptiSparseMET:::.safe_logdet_psd_dense))
  expect_true(is.function(OptiSparseMET:::.pairwise_diff_mean_var))
  expect_true(is.function(OptiSparseMET:::.trace_subinv_est))
  expect_true(is.function(OptiSparseMET:::.with_local_seed))
  expect_true(is.function(OptiSparseMET:::.build_neighbor_pairs))
  expect_true(is.function(OptiSparseMET:::.score_dispersion))
  
  # Confirm they are not in the public namespace
  expect_false(".make_sparse_incidence" %in% getNamespaceExports("OptiSparseMET"))
  expect_false(".ar1_precision_sparse"  %in% getNamespaceExports("OptiSparseMET"))
})

test_that("example dataset loads and has expected structure", {
  data("OptiSparseMET_example_data", package = "OptiSparseMET")
  expect_true(exists("OptiSparseMET_example_data"))
  expect_true(is.list(OptiSparseMET_example_data))
  
  expected_names <- c(
    "treatments", "environments", "common_treatments",
    "treatment_info", "seed_info", "seed_required_per_plot",
    "env_design_specs"
  )
  expect_true(all(expected_names %in% names(OptiSparseMET_example_data)))
})