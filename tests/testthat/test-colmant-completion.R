# Tests for the remaining Colmant (2026) decision-point tools: BV targets,
# sparsity grid, environment-relationship builder, and the orchestrator.

mk <- function(J = 40L, E = 4L, seed = 11) {
  set.seed(seed)
  B <- matrix(rnorm(J * J), J, J); G <- tcrossprod(B) / J + diag(J) * 0.3
  dimnames(G) <- list(paste0("L", seq_len(J)), paste0("L", seq_len(J)))
  # two mega-environments: {E1,E2} correlated, {E3,E4} correlated
  SigE <- matrix(c(1, .8, .1, 0,  .8, 1, 0, .1,
                   .1, 0, 1, .8,   0, .1, .8, 1), 4, 4)
  dimnames(SigE) <- list(paste0("E", 1:4), paste0("E", 1:4))
  list(G = G, SigE = SigE, J = J, E = E)
}

# ---- BV targets: broad vs specific adaptation --------------------------------

test_that("mega-environment target rewards the design that tests those environments", {
  d <- mk(); ln <- rownames(d$G)
  A <- matrix(0L, d$J, 4, dimnames = list(ln, colnames(d$SigE))); A[, c(1, 2)] <- 1L
  B <- matrix(0L, d$J, 4, dimnames = list(ln, colnames(d$SigE))); B[, c(3, 4)] <- 1L

  # Under the mega-env {E1,E2} target, A (which tests it) beats B by a lot.
  accA <- simulate_met(A, G = d$G, Sigma_E = d$SigE, n_sim = 40,
                       bv_target = "mega_environment", target_envs = c("E1", "E2"),
                       seed = 1)$accuracy_mean
  accB <- simulate_met(B, G = d$G, Sigma_E = d$SigE, n_sim = 40,
                       bv_target = "mega_environment", target_envs = c("E1", "E2"),
                       seed = 1)$accuracy_mean
  expect_gt(accA, accB)
})

test_that("environment_specific needs exactly one target environment", {
  d <- mk(); ln <- rownames(d$G)
  M <- matrix(1L, d$J, 4, dimnames = list(ln, colnames(d$SigE)))
  expect_error(simulate_met(M, G = d$G, Sigma_E = d$SigE, n_sim = 5,
                            bv_target = "environment_specific"))
  expect_error(simulate_met(M, G = d$G, Sigma_E = d$SigE, n_sim = 5,
                            bv_target = "environment_specific",
                            target_envs = c("E1", "E2")))
  res <- simulate_met(M, G = d$G, Sigma_E = d$SigE, n_sim = 10,
                      bv_target = "environment_specific", target_envs = "E1", seed = 1)
  expect_true(is.finite(res$accuracy_mean))
})

# ---- build_environment_relationship -----------------------------------------

test_that("genetic-correlation E matrix is a correlation matrix over environments", {
  set.seed(2)
  GxE <- matrix(rnorm(50 * 4), 50, 4, dimnames = list(NULL, paste0("E", 1:4)))
  GxE[, 2] <- GxE[, 1] + rnorm(50, 0, 0.1)   # E1, E2 nearly identical
  D <- build_environment_relationship(GxE, source = "genetic_correlation")
  expect_equal(dim(D), c(4L, 4L))
  expect_equal(diag(D), rep(1, 4), tolerance = 1e-8, ignore_attr = TRUE)
  expect_gt(D["E1", "E2"], 0.9)              # captured the strong correlation
})

test_that("enviromic Gaussian kernel yields a valid similarity matrix", {
  set.seed(3)
  cov_mat <- matrix(rnorm(6 * 5), 6, 5, dimnames = list(paste0("E", 1:6), NULL))
  D <- build_environment_relationship(cov_mat, source = "enviromic", kernel = "gaussian")
  expect_equal(dim(D), c(6L, 6L))
  expect_equal(diag(D), rep(1, 6), tolerance = 1e-8, ignore_attr = TRUE)  # self-similarity = 1
  expect_true(all(D >= 0 & D <= 1))
  # feeds select_environments
  sel <- select_environments(D, n = 3, method = "representative")$selected
  expect_equal(length(sel), 3L)
})

# ---- sparsity_grid ----------------------------------------------------------

test_that("sparsity_grid sweeps the grid and returns ranked accuracy", {
  d <- mk(J = 60L)
  tab <- sparsity_grid(
    G = d$G, Sigma_E = d$SigE,
    ia_values = c(40, 60), ipf_values = c(20, 40), noe_values = c(2, 4),
    plot_budget = 200, n_sim = 15, seed = 1)
  expect_true(all(c("ia", "ipf", "noe", "total_plots", "sparsity_pct",
                    "accuracy_mean") %in% names(tab)))
  expect_true(all(tab$total_plots <= 200))
  # sorted best-first
  expect_equal(tab$accuracy_mean, sort(tab$accuracy_mean, decreasing = TRUE))
})

# ---- run_design_strategy orchestrator ---------------------------------------

test_that("run_design_strategy chains selection, allocation and evaluation", {
  d <- mk()
  set.seed(4)
  GxE <- matrix(rnorm(d$J * d$E), d$J, d$E, dimnames = list(NULL, colnames(d$SigE)))
  Denv <- build_environment_relationship(GxE, source = "genetic_correlation")

  out <- run_design_strategy(
    treatments = rownames(d$G),
    n_test_entries_per_environment = 20,
    G = d$G, Sigma_E = d$SigE,
    env_relationship = Denv, n_environments = 3,
    allocation_method = "random_balanced",
    evaluate = TRUE, n_sim = 15, seed = 1)

  expect_equal(length(out$decisions$environments), 3L)
  expect_true(is.matrix(out$allocation_matrix))
  expect_equal(ncol(out$allocation_matrix), 3L)
  expect_true(is.finite(out$evaluation$accuracy_mean))
  expect_true(is.finite(out$evaluation$mean_PEV))
})

test_that("run_design_strategy applies individual selection and G x E refinement", {
  d <- mk()
  out <- run_design_strategy(
    treatments = rownames(d$G),
    environments = colnames(d$SigE),
    n_test_entries_per_environment = 15,
    G = d$G, Sigma_E = d$SigE,
    n_individuals = 30, individual_criterion = "cdmean",
    refine_gxe = TRUE, evaluate = TRUE, n_sim = 10, seed = 1, max_dim = 2000)
  expect_equal(length(out$decisions$individuals), 30L)
  expect_true(out$decisions$gxe_refined)
  expect_equal(nrow(out$allocation_matrix), 30L)
})
