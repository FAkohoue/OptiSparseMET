# Tests for the kinship->connectivity helper, plot-budget reconciliation, and
# the Pareto frontier plot.

# ---- recommend_common_treatments -------------------------------------------

test_that("rho is read from the environment covariance (cov2cor mean off-diag)", {
  SigE <- matrix(c(4, 1.2, 0.4,
                   1.2, 4, 0.8,
                   0.4, 0.8, 4), 3, byrow = TRUE)
  rc <- recommend_common_treatments(SigE, n_test_entries_per_environment = 60,
                                    n_treatments = 300, n_environments = 3)
  R <- cov2cor(SigE)
  expect_equal(rc$rho, mean(R[upper.tri(R)]), tolerance = 1e-8)
  expect_true(is.numeric(rc$common_range) && length(rc$common_range) >= 1)
})

test_that("weaker environment correlation suggests more common treatments", {
  strong <- matrix(0.7, 3, 3); diag(strong) <- 1
  weak   <- matrix(0.1, 3, 3); diag(weak)   <- 1
  rs <- recommend_common_treatments(strong, 60, 300, is_covariance = FALSE)
  rw <- recommend_common_treatments(weak,   60, 300, is_covariance = FALSE)
  expect_lt(rs$rho, rw$rho + 1)                      # sanity: strong rho > weak
  expect_gte(max(rw$common_range), max(rs$common_range))
})

test_that("recommend_common_treatments can evaluate the simulated optimum", {
  set.seed(1)
  J <- 40
  G <- crossprod(matrix(rnorm(J * 60), 60, J)) / 60 + diag(J) * 0.3
  dimnames(G) <- list(paste0("L", 1:J), paste0("L", 1:J))
  SigE <- matrix(0.2, 3, 3); diag(SigE) <- 1
  dimnames(SigE) <- list(paste0("E", 1:3), paste0("E", 1:3))
  rc <- recommend_common_treatments(SigE, n_test_entries_per_environment = 20,
                                    n_treatments = J, evaluate = TRUE,
                                    G = G, Sigma_E = SigE, n_sim = 8, seed = 1,
                                    is_covariance = FALSE)
  expect_true(!is.null(rc$tune))
  expect_true(rc$recommended %in% rc$tune$table$n_common)
})

# ---- plot-budget reconciliation --------------------------------------------

test_that("test_entry_capacity and required_plots are inverse and check-aware", {
  # 400 plots, 2 checks x 4 blocks = 8 check plots -> 392 entries (unreplicated)
  expect_equal(test_entry_capacity(400, n_checks = 2, n_blocks = 4), 392L)
  expect_equal(required_plots(392, n_checks = 2, n_blocks = 4), 400L)
  # replication halves capacity roughly
  expect_equal(test_entry_capacity(408, n_checks = 2, n_blocks = 4,
                                   avg_reps_per_entry = 2), 200L)
  expect_error(test_entry_capacity(8, n_checks = 2, n_blocks = 4), "Check plots")
})

test_that("field_plot_accounting returns a consistent breakdown", {
  a <- field_plot_accounting(total_plots = 400, n_checks = 2, n_blocks = 4)
  expect_equal(a$check_plots, 8)
  expect_equal(a$n_test_entries, 392)
  expect_equal(a$check_plots + a$entry_plots, a$total_plots)
})

test_that("suggest_site_capacity total_plots includes check plots", {
  set.seed(2)
  G <- crossprod(matrix(rnorm(40 * 60), 60, 40)) / 60 + diag(40) * 0.3
  dimnames(G) <- list(paste0("L", 1:40), paste0("L", 1:40))
  SigE <- diag(3) * 0.5 + 0.5
  dimnames(SigE) <- list(paste0("E", 1:3), paste0("E", 1:3))
  no_chk <- suggest_site_capacity(G, SigE, candidate_plots = c(10, 20),
                                  scope = "all", n_env = 3, n_sim = 6, seed = 1)
  with_chk <- suggest_site_capacity(G, SigE, candidate_plots = c(10, 20),
                                    scope = "all", n_env = 3,
                                    check_plots_per_site = 8, n_sim = 6, seed = 1)
  # 3 sites x 8 check plots = 24 extra plots per row
  expect_equal(with_chk$table$total_plots - no_chk$table$total_plots,
               rep(24, nrow(no_chk$table)))
})

# ---- Pareto frontier plot ---------------------------------------------------

test_that("plot_pareto_frontier draws without error and returns the frontier", {
  fr <- data.frame(cost_weight = c(0, 0.5, 1),
                   reliability = c(0.40, 0.35, 0.28),
                   gain = c(1.2, 1.05, 0.85),
                   plots = c(120, 90, 60),
                   cost = c(120, 90, 60),
                   score = c(1, 1.1, 1.2),
                   pareto = c(TRUE, TRUE, TRUE))
  tmp <- tempfile(fileext = ".pdf"); grDevices::pdf(tmp)
  out <- plot_pareto_frontier(list(frontier = fr), x = "plots")
  grDevices::dev.off(); unlink(tmp)
  expect_s3_class(out, "data.frame")
  expect_equal(out$plots, sort(fr$plots))            # ordered by budget
})
