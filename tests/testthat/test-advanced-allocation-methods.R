advanced_fixture <- function(J = 8L, E = 3L) {
  trt <- paste0("L", seq_len(J)); env <- paste0("E", seq_len(E))
  G <- diag(J) * 0.8 + 0.2
  dimnames(G) <- list(trt, trt)
  Sigma_E <- diag(E) * 0.65 + 0.35
  dimnames(Sigma_E) <- list(env, env)
  list(trt = trt, env = env, G = G, Sigma_E = Sigma_E)
}

test_that("design_objective exposes explicit prediction criteria", {
  x <- advanced_fixture(6, 3)
  M <- matrix(c(1, 1, 0, 1, 0, 1,
                1, 0, 1, 0, 1, 1,
                0, 1, 1, 1, 1, 0), 6, 3,
              dimnames = list(x$trt, x$env))
  pev <- design_objective(M, x$G, x$Sigma_E, criterion = "mean_pev")
  cd <- design_objective(M, x$G, x$Sigma_E, criterion = "cdmean")
  gain <- design_objective(M, x$G, x$Sigma_E,
                           criterion = "expected_gain")
  expect_equal(pev$score, -pev$mean_PEV)
  expect_equal(cd$score, cd$reliability)
  expect_equal(gain$score, gain$gain)
})

test_that("optimizers expose exchange, genetic, and exact engines", {
  x <- advanced_fixture(6, 3)
  M <- matrix(0L, 6, 3, dimnames = list(x$trt, x$env))
  M[cbind(c(1, 2, 3, 4, 5, 6), c(1, 2, 3, 1, 2, 3))] <- 1L
  M[cbind(c(1, 2, 3, 4, 5, 6), c(2, 3, 1, 2, 3, 1))] <- 1L

  ex <- optimize_design(M, x$G, x$Sigma_E, search_method = "exchange",
                        allocation_criterion = "mean_pev",
                        n_starts = 1, iters = 8, seed = 1)
  expect_equal(rowSums(ex$allocation_matrix), rowSums(M))
  expect_equal(colSums(ex$allocation_matrix), colSums(M))
  expect_gte(ex$score, ex$score_start)

  ga <- optimize_design(M, x$G, x$Sigma_E, search_method = "genetic",
                        allocation_criterion = "cdmean",
                        n_starts = 4, iters = 4, seed = 2)
  expect_equal(ga$search_method, "genetic")
  expect_gte(ga$score, ga$score_start)

  y <- advanced_fixture(3, 2)
  M2 <- matrix(c(1, 0, 1, 0, 1, 1), 3, 2,
               dimnames = list(y$trt, y$env))
  exact <- optimize_design(M2, y$G, y$Sigma_E,
                           search_method = "mip", exact_max_cells = 6,
                           allocation_criterion = "mean_pev",
                           n_starts = 1, iters = 1)
  expect_true(exact$diagnostics$exact_certified)
  expect_equal(exact$diagnostics$search_method, "exact")
})

test_that("joint common-set optimization respects requested count bounds", {
  x <- advanced_fixture(8, 3)
  out <- allocate_sparse_met(
    x$trt, x$env, allocation_method = "prediction_optimal",
    n_test_entries_per_environment = 4,
    G = x$G, Sigma_E = x$Sigma_E,
    common_set = "optimize_jointly",
    common_control = list(count_range = c(1, 2), weight = 0.2),
    search_method = "exchange",
    optimizer_control = list(n_starts = 1, iters = 12), seed = 3)
  n_common <- sum(rowSums(out$allocation_matrix) == length(x$env))
  expect_gte(n_common, 1L)
  expect_lte(n_common, 2L)
  expect_equal(out$summary$allocation_method, "prediction_optimal")
  expect_equal(out$summary$constructor_method, "random_balanced")
})

test_that("robust scenarios support operational uncertainty and site loss", {
  x <- advanced_fixture(6, 3)
  M <- matrix(1L, 6, 3, dimnames = list(x$trt, x$env))
  scenarios <- robust_scenarios(
    sigma_e2 = 1, sigmaE_shrink = 0,
    operational_scenarios = list(
      baseline = list(),
      lose_E2 = list(site_loss = "E2"),
      expensive_E3 = list(cost_per_plot = c(E1 = 1, E2 = 1, E3 = 3)),
      shifted_tpe = list(tpe_weights = c(E1 = 0.7, E2 = 0.2, E3 = 0.1))))
  ans <- robust_design_score(
    M, x$G, x$Sigma_E, scenarios, aggregate = "cvar",
    criterion = "mean_pev")
  expect_length(ans$scenario_scores, 4L)
  expect_true(all(is.finite(ans$scenario_scores)))
  expect_setequal(ans$operational_scenario,
                  c("baseline", "lose_E2", "expensive_E3", "shifted_tpe"))
})

test_that("adaptive allocation returns an auditable feasible batch", {
  x <- advanced_fixture(7, 3)
  M <- matrix(0L, 7, 3, dimnames = list(x$trt, x$env))
  out <- adaptive_met_allocation(
    M, x$G, x$Sigma_E, target_environment_sizes = 4,
    batch_size = 12, required_common_treatments = "L1",
    minimum_treatment_environments = 1, seed = 5)
  expect_equal(unname(colSums(out$allocation_matrix)), rep(4L, 3))
  expect_true(all(out$allocation_matrix["L1", ] == 1L))
  expect_true(all(rowSums(out$allocation_matrix) >= 1L))
  expect_equal(nrow(out$recommendations), 12L)
  expect_gte(out$score_after, out$score_before)
  expect_s3_class(out$design, "sparse_met_design")
})

test_that("robust_prediction is available through the breeder-facing allocator", {
  x <- advanced_fixture(7, 3)
  out <- allocate_sparse_met(
    x$trt, x$env, allocation_method = "robust_prediction",
    n_test_entries_per_environment = 4,
    G = x$G, Sigma_E = x$Sigma_E,
    search_method = "exchange",
    optimizer_control = list(n_starts = 1, iters = 5,
                             robust_sigma_e2 = c(0.8, 1.2),
                             sigmaE_shrink = c(0, 0.5)), seed = 7)
  expect_true(out$advanced_optimization$robust)
  expect_equal(out$summary$allocation_method, "robust_prediction")
  expect_equal(unname(colSums(out$allocation_matrix)), rep(4L, 3))
})

test_that("run_design_strategy forwards advanced allocation controls", {
  x <- advanced_fixture(7, 3)
  out <- run_design_strategy(
    x$trt, x$env, n_test_entries_per_environment = 4,
    G = x$G, Sigma_E = x$Sigma_E,
    allocation_method = "prediction_optimal",
    allocation_criterion = "cdmean", search_method = "exchange",
    optimizer_control = list(n_starts = 1, iters = 5),
    evaluate = FALSE, seed = 9)
  expect_equal(out$decisions$allocation_method, "prediction_optimal")
  expect_equal(out$decisions$allocation_criterion, "cdmean")
  expect_equal(out$decisions$search_method, "exchange")
  expect_equal(unname(colSums(out$allocation_matrix)), rep(4L, 3))
})
