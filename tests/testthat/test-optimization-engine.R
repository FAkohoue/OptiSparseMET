# Tests for the criterion-driven optimal-design engine: selection intensity and
# expected genetic gain, the combined design objective, robust/Bayesian scoring,
# the exchange optimiser, and the Pareto explorer. All deterministic (numpy-
# verified reference values where relevant).

mk_opt <- function(J = 24L, E = 3L, seed = 1) {
  set.seed(seed)
  B <- matrix(rnorm(J * 60), 60, J)
  G <- crossprod(B) / 60 + diag(J) * 0.3
  dimnames(G) <- list(paste0("L", seq_len(J)), paste0("L", seq_len(J)))
  SigE <- diag(E) * 0.5 + 0.5
  dimnames(SigE) <- list(paste0("E", seq_len(E)), paste0("E", seq_len(E)))
  M <- matrix(0L, J, E, dimnames = list(rownames(G), colnames(SigE)))
  for (e in seq_len(E)) M[sample.int(J, 12L), e] <- 1L
  list(G = G, SigE = SigE, M = M)
}

# ---- Layer 1: selection intensity + expected genetic gain -------------------

test_that("selection_intensity matches standard values", {
  expect_equal(selection_intensity(0.10), 1.754983, tolerance = 1e-5)
  expect_equal(selection_intensity(0.05), 2.062713, tolerance = 1e-5)
  expect_equal(selection_intensity(0.01), 2.665214, tolerance = 1e-5)
  expect_error(selection_intensity(0), "in \\(0, 1\\)")
  expect_error(selection_intensity(1), "in \\(0, 1\\)")
})

test_that("expected_genetic_gain follows the breeder's equation", {
  g <- expected_genetic_gain(reliability = 0.5, sigma_g = 1, prop = 0.1)
  expect_equal(g$gain, 1.754983 * sqrt(0.5), tolerance = 1e-5)  # = 1.24096
  expect_equal(g$accuracy, sqrt(0.5), tolerance = 1e-8)
  # monotone increasing in reliability
  lo <- expected_genetic_gain(reliability = 0.2, prop = 0.1)$gain
  hi <- expected_genetic_gain(reliability = 0.8, prop = 0.1)$gain
  expect_gt(hi, lo)
})

test_that("multi-trait economic index gain uses sqrt(w' Sigma_T w)", {
  ST <- matrix(c(1, .3, .3, .5), 2)
  g <- expected_genetic_gain(reliability = 0.5, prop = 0.1,
                             trait_weights = c(2, -1), trait_gencov = ST)
  expect_equal(g$sigma, sqrt(sum(c(2, -1) * (ST %*% c(2, -1)))), tolerance = 1e-6)
  expect_equal(g$gain, 1.754983 * sqrt(0.5) * g$sigma, tolerance = 1e-5) # 2.25432
})

test_that("expected_genetic_gain can read reliability from met_information", {
  d <- mk_opt()
  info <- met_information(d$M, G = d$G, Sigma_E = d$SigE)
  g <- expected_genetic_gain(info = info, prop = 0.1)$gain
  expect_equal(g, selection_intensity(0.1) * sqrt(info$CDmean), tolerance = 1e-8)
})

# ---- Layer 2: design_objective ----------------------------------------------

test_that("design_objective returns components and honours a budget", {
  d <- mk_opt()
  o <- design_objective(d$M, d$G, d$SigE,
                        weights = list(gain = 1, reliability = 0, cost = 0))
  expect_true(all(c("reliability", "gain", "plots", "cost", "score") %in% names(o)))
  expect_equal(o$plots, sum(d$M != 0))
  expect_true(o$feasible)
  # a budget below the design's plot count makes it infeasible
  tight <- design_objective(d$M, d$G, d$SigE, budget = sum(d$M) - 1)
  expect_false(tight$feasible)
  expect_identical(tight$score, -Inf)
})

test_that("design_objective counts physical replication as plots", {
  d <- mk_opt()
  reps <- 2 * d$M
  o <- design_objective(d$M, d$G, d$SigE, reps = reps)
  expect_equal(o$plots, sum(reps))
  expect_false(design_objective(
    d$M, d$G, d$SigE, reps = reps,
    budget = sum(reps) - 1)$feasible)
})

test_that("a higher cost weight lowers the score of the same design", {
  d <- mk_opt()
  ref <- list(gain = 1, reliability = 1, cost = sum(d$M))
  s0 <- design_objective(d$M, d$G, d$SigE, weights = list(gain = 1, cost = 0),
                         ref = ref)$score
  s2 <- design_objective(d$M, d$G, d$SigE, weights = list(gain = 1, cost = 2),
                         ref = ref)$score
  expect_lt(s2, s0)                     # the cost term subtracts
})

# ---- Layer 3: robust scoring ------------------------------------------------

test_that("robust_scenarios builds a normalised prior", {
  sc <- robust_scenarios(sigma_e2 = c(0.5, 1, 2), sigmaE_shrink = c(0, 0.5))
  expect_equal(length(sc), 6L)
  expect_equal(sum(vapply(sc, function(s) s$prob, numeric(1))), 1, tolerance = 1e-8)
  expect_error(robust_scenarios(sigmaE_shrink = 1.5), "in \\[0, 1\\]")
  expect_error(robust_scenarios(sigma_g2 = 0), "positive")
  expect_error(robust_scenarios(sigma_e2 = NA_real_), "positive")
})

test_that("robust_design_score aggregates mean/min/cvar consistently", {
  d <- mk_opt()
  sc <- robust_scenarios(sigma_e2 = c(0.5, 1, 2), sigmaE_shrink = c(0, 0.5))
  mean_s <- robust_design_score(d$M, d$G, d$SigE, sc, aggregate = "mean")
  min_s  <- robust_design_score(d$M, d$G, d$SigE, sc, aggregate = "min")$score
  cvar_s <- robust_design_score(d$M, d$G, d$SigE, sc, aggregate = "cvar",
                                cvar_alpha = 0.5)$score
  expect_length(mean_s$scenario_scores, 6L)
  # worst-case <= cvar(worst half) <= mean
  expect_lte(min_s, cvar_s + 1e-9)
  expect_lte(cvar_s, mean_s$score + 1e-9)
  expect_error(
    robust_design_score(d$M, d$G, d$SigE, list(), aggregate = "mean"),
    "non-empty")
  expect_error(
    robust_design_score(d$M, d$G, d$SigE, sc, aggregate = "cvar",
                        cvar_alpha = 1.5),
    "in \\[0, 1\\]")
})

# ---- Layer 4: optimize_design -----------------------------------------------

test_that("optimize_design never worsens the starting score", {
  d <- mk_opt()
  opt <- optimize_design(d$M, d$G, d$SigE, preserve = "margins",
                         n_starts = 2, iters = 40, seed = 1)
  expect_gte(opt$score, opt$score_start)
})

test_that("preserve = 'margins' keeps both replication and environment sizes", {
  d <- mk_opt()
  opt <- optimize_design(d$M, d$G, d$SigE, preserve = "margins",
                         n_starts = 2, iters = 50, seed = 2)
  expect_equal(colSums(opt$allocation_matrix), colSums(d$M))
  expect_equal(rowSums(opt$allocation_matrix), rowSums(d$M))
})

test_that("preserve = 'replication' keeps replication but may resize sites", {
  d <- mk_opt()
  opt <- optimize_design(d$M, d$G, d$SigE, preserve = "replication",
                         n_starts = 2, iters = 50, seed = 3)
  expect_equal(rowSums(opt$allocation_matrix), rowSums(d$M))  # equal replication
})

test_that("optimize_design accepts a robust objective", {
  d <- mk_opt()
  sc <- robust_scenarios(sigma_e2 = c(0.5, 2), sigmaE_shrink = c(0, 0.5))
  opt <- optimize_design(d$M, d$G, d$SigE, preserve = "margins",
                         robust = sc, robust_aggregate = "cvar",
                         n_starts = 1, iters = 30, seed = 4)
  expect_true(is.finite(opt$score))
  expect_true(opt$robust)
})

test_that("optimize_design preserves network seed and site capacities", {
  d <- mk_opt(J = 18L, E = 3L)
  available <- stats::setNames(rep(8, nrow(d$M)), rownames(d$M))
  seed_cost <- c(E1 = 1, E2 = 2, E3 = 3)
  opt <- optimize_design(
    d$M, d$G, d$SigE, preserve = "none",
    seed_available = available,
    seed_required_per_environment = seed_cost,
    minimum_seed_buffer = 1,
    environment_capacities = colSums(d$M),
    n_starts = 1, iters = 30, seed = 8)
  used <- rowSums(sweep(opt$allocation_matrix, 2, seed_cost, `*`))
  expect_true(all(used <= available - 1 + 1e-8))
  expect_true(all(colSums(opt$allocation_matrix) <= colSums(d$M)))
  expect_true(is.data.frame(opt$seed_summary))
})

# ---- Layer 5: pareto_designs ------------------------------------------------

test_that(".pareto_flag identifies dominated points", {
  expect_equal(.pareto_flag(c(1, 2, 3), c(1, 2, 3)), c(TRUE, TRUE, TRUE))
  # (ben 2, cost 3) dominated by (ben 2, cost 2)
  expect_equal(.pareto_flag(c(1, 2, 2), c(1, 2, 3)), c(TRUE, TRUE, FALSE))
})

test_that("pareto_designs returns a frontier with a pareto flag", {
  d <- mk_opt()
  pf <- pareto_designs(d$M, d$G, d$SigE, cost_weights = c(0, 0.5, 1),
                       preserve = "none", n_starts = 1, iters = 25, seed = 1)
  expect_true(all(c("cost_weight", "gain", "cost", "pareto") %in% names(pf$frontier)))
  expect_equal(nrow(pf$frontier), 3L)
  expect_type(pf$frontier$pareto, "logical")
  expect_length(pf$designs, 3L)
})
