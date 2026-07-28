test_that("reference benchmark designs distinguish full, random, M3, and M4", {
  treatments <- paste0("L", 1:12)
  environments <- paste0("E", 1:4)
  out <- make_met_benchmark_designs(
    treatments, environments,
    n_test_entries_per_environment = 6,
    target_replications = 2,
    seed = 7
  )

  expect_s3_class(out, "met_benchmark_designs")
  expect_true(all(c("full", "random_sparse", "M3_coverage_first",
                    "M4_equireplicate") %in% names(out$designs)))
  expect_true(out$m4_feasibility$feasible)
  expect_equal(
    rowSums(out$designs$M4_equireplicate),
    setNames(rep(2, 12), treatments)
  )
  expect_equal(
    colSums(out$designs$M4_equireplicate),
    setNames(rep(6, 4), environments)
  )
  expect_equal(out$m4_status, "strict_exact")
})

test_that("infeasible strict M4 is reported rather than mislabelled", {
  out <- make_met_benchmark_designs(
    paste0("L", 1:11), paste0("E", 1:4),
    n_test_entries_per_environment = 6,
    target_replications = 2,
    m4_fallback = "omit", seed = 2
  )
  expect_false(out$m4_feasibility$feasible)
  expect_false("M4_equireplicate" %in% names(out$designs))
  expect_equal(out$m4_status, "omitted")
})

test_that("MET benchmark uses paired draws and reports robust feasibility", {
  treatments <- paste0("L", 1:12)
  environments <- paste0("E", 1:4)
  designs <- make_met_benchmark_designs(
    treatments, environments, 6, 2, seed = 3
  )
  X <- matrix(seq_len(12 * 18) / 100, 12, 18)
  G <- tcrossprod(X) + diag(12) * 0.3
  dimnames(G) <- list(treatments, treatments)
  Sigma <- matrix(0.35, 4, 4)
  diag(Sigma) <- 1
  dimnames(Sigma) <- list(environments, environments)

  out <- benchmark_met_designs(
    designs, G = G, Sigma_E = Sigma,
    sigma_e2 = 1,
    environment_dropout = list(none = character(), lose_E4 = "E4"),
    environment_plot_capacity = 6,
    cost_per_plot = c(E1 = 1, E2 = 1, E3 = 1.2, E4 = 1.2),
    seed_available = setNames(rep(20, 12), treatments),
    seed_required_per_plot = 1,
    n_sim = 8, seed = 11
  )

  expect_true(out$paired_simulation)
  expect_equal(nrow(out$scenario_results),
               length(out$designs) * 2L)
  expect_true(all(c("accuracy", "worst_accuracy", "gain",
                    "cost", "feasible", "pareto_accuracy",
                    "cost_minimisation_percent") %in%
                  names(out$design_summary)))
  expect_false(out$feasibility$plot_capacity_feasible[
    out$feasibility$design == "full"
  ])
  expect_true(all(
    !out$design_summary$pareto_accuracy | out$design_summary$feasible
  ))
  expect_true(all(out$design_summary$cost_minimisation_percent >= 0 &
                    out$design_summary$cost_minimisation_percent <= 100))
})

test_that("environment-model benchmark labels equal weighting as a comparator", {
  env <- paste0("E", 1:5)
  K1 <- outer(1:5, 1:5, function(i, j) exp(-abs(i - j)))
  K2 <- outer(c(1, 3, 2, 5, 4), c(1, 3, 2, 5, 4),
              function(i, j) exp(-abs(i - j)))
  dimnames(K1) <- dimnames(K2) <- list(env, env)
  truth <- 0.75 * K1 + 0.15 * K2 + 0.10 * diag(5)
  dimnames(truth) <- list(env, env)

  out <- benchmark_environment_models(
    list(weather = K1, soil = K2),
    target = truth, truth = truth, n_boot = 0, seed = 1
  )
  expect_true(all(c("independent", "single_equal_weight",
                    "separate_weather", "separate_soil",
                    "historically_calibrated") %in% out$scores$model))
  expect_equal(out$reference_basis, "known_truth")
  expect_equal(out$calibration$status, "historically_calibrated")
  expect_lt(
    out$scores$rmse[out$scores$model == "historically_calibrated"],
    out$scores$rmse[out$scores$model == "single_equal_weight"]
  )
})

test_that("known-truth simulation returns weight and covariance recovery", {
  env <- paste0("E", 1:5)
  K1 <- outer(1:5, 1:5, function(i, j) exp(-abs(i - j) / 2))
  K2 <- outer(c(1, 3, 2, 5, 4), c(1, 3, 2, 5, 4),
              function(i, j) exp(-abs(i - j) / 2))
  dimnames(K1) <- dimnames(K2) <- list(env, env)
  out <- simulate_environment_model_benchmark(
    list(weather = K1, soil = K2),
    true_weights = c(weather = 0.7, soil = 0.2, identity = 0.1),
    n_genotypes = 30, n_years = 4,
    n_boot = 0, seed = 5
  )

  expect_equal(nrow(out$historical), 30 * 5 * 4)
  expect_equal(sum(out$true_weights), 1)
  expect_true(all(c("component", "true_weight", "estimated_weight",
                    "absolute_error") %in% names(out$weight_recovery)))
  expect_equal(out$benchmark$reference_basis, "known_truth")
})

test_that("environmental missingness benchmark rebuilds every perturbation", {
  env <- paste0("E", 1:5)
  weather <- data.frame(
    environment = env,
    temperature = c(24, 25, 23, 27, 26),
    rainfall = c(500, 430, 610, 390, 470)
  )
  soil <- data.frame(
    environment = env,
    pH = c(5.5, 6.0, 5.8, 6.2, 5.9),
    clay = c(30, 25, 35, 20, 28)
  )
  base <- build_environment_kernels(weather = weather, soil = soil)
  truth <- 0.7 * base$kernels$weather + 0.3 * base$kernels$soil
  out <- benchmark_environment_missingness(
    weather = weather, soil = soil,
    missing_rates = c(0, 0.1), n_repeats = 2,
    truth = truth, seed = 4
  )

  expect_equal(nrow(out$masking), 3L)
  expect_equal(out$masking$masked_cells[out$masking$missing_rate == 0], 0)
  expect_true(all(c("qc_warning_count", "qc_warnings") %in%
                    names(out$masking)))
  expect_gt(sum(out$masking$qc_warning_count), 0)
  expect_gt(nrow(out$warnings), 0)
  expect_true(all(c("model", "missing_rate", "rmse_mean") %in%
                    names(out$summary)))
  expect_length(out$rebuilt, 3L)
})

test_that("design stability reports common, capacity, and allocation stability", {
  A1 <- matrix(c(1, 0, 1, 1, 0, 1), 3, 2,
               dimnames = list(paste0("L", 1:3), c("E1", "E2")))
  A2 <- matrix(c(1, 1, 0, 1, 0, 1), 3, 2,
               dimnames = dimnames(A1))
  out <- summarize_design_stability(
    common_sets = list(a = c("L1", "L2"), b = c("L1", "L3")),
    site_capacities = list(a = c(E1 = 20, E2 = 22),
                           b = c(E1 = 21, E2 = 22)),
    allocations = list(a = A1, b = A2)
  )

  expect_equal(out$common_frequency$selection_frequency[
    out$common_frequency$Treatment == "L1"
  ], 1)
  expect_equal(out$common_jaccard["a", "b"], 1 / 3)
  expect_equal(nrow(out$capacity_summary), 2)
  expect_equal(dim(out$allocation_frequency), c(3, 2))
})

test_that("Pareto plot maps maximum and minimum cost to 0 and 100 percent", {
  fr <- data.frame(
    design = c("full", "middle", "small"),
    accuracy = c(0.84, 0.80, 0.70),
    cost = c(10, 6, 1),
    pareto_accuracy = TRUE
  )
  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp)
  out <- plot_pareto_frontier(fr, metric = "accuracy", label_designs = TRUE)
  grDevices::dev.off()
  unlink(tmp)

  expect_equal(out$cost_minimisation_percent[out$cost == 10], 0)
  expect_equal(out$cost_minimisation_percent[out$cost == 1], 100)
  expect_equal(out$performance_maximisation_percent[out$accuracy == 0.84], 84)
})

test_that("infeasible designs cannot define the plotted Pareto frontier", {
  fr <- data.frame(
    design = c("feasible", "infeasible"),
    accuracy = c(0.75, 0.95),
    cost = c(6, 2),
    feasible = c(TRUE, FALSE),
    pareto_accuracy = c(TRUE, TRUE)
  )
  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp)
  out <- plot_pareto_frontier(fr, metric = "accuracy")
  grDevices::dev.off()
  unlink(tmp)

  expect_true(out$pareto_plotted[out$design == "feasible"])
  expect_false(out$pareto_plotted[out$design == "infeasible"])
  expect_false(out$frontier_eligible[out$design == "infeasible"])
})
