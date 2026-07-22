# End-to-end tests for the new spatial options THROUGH the evaluators and
# optimisers (the helper-level maths is covered in test-spatial-structures.R).

# 6 x 8 grid, 12 genotypes replicated 4 times, no checks.
alpha_grid_fb <- function(n_rows = 6L, n_cols = 8L, n_geno = 12L) {
  N <- n_rows * n_cols
  data.frame(
    Plot       = seq_len(N),
    Row        = rep(seq_len(n_rows), times = n_cols),
    Column     = rep(seq_len(n_cols), each = n_rows),
    Rep        = rep(1:4, each = N / 4),
    IBlock     = rep(1:8, each = N / 8),
    BlockInRep = rep(1:2, length.out = N),
    Treatment  = rep(paste0("G", seq_len(n_geno)), length.out = N),
    Check      = FALSE,
    stringsAsFactors = FALSE
  )
}
alpha_vc <- list(sigma_e2 = 2, sigma_g2 = 1, sigma_rep2 = 0.5,
                 sigma_ib2 = 0.3, sigma_r2 = 0.2, sigma_c2 = 0.2)

ev_alpha <- function(fb, ...) {
  met_evaluate_alpha_efficiency(
    field_book = fb, n_rows = 6, n_cols = 8,
    check_treatments = character(0),
    treatment_effect = "random", prediction_type = "IID",
    varcomp = alpha_vc, ...)
}

# ---- P-spline end-to-end: the lambda -> Inf invariant -----------------------

test_that("a heavily penalised P-spline reproduces the no-spline model", {
  fb <- alpha_grid_fb()
  none <- ev_alpha(fb, residual_structure = "IID", spatial_random = "none")
  stiff <- ev_alpha(fb, residual_structure = "IID", spatial_random = "pspline",
                    spline_knots_row = 5, spline_knots_col = 5,
                    spline_lambda_row = 1e8, spline_lambda_col = 1e8)
  # Verified independently in numpy: as the smoothing parameters grow the
  # surface shrinks away and the criteria converge to the no-spline model.
  expect_equal(stiff$mean_PEV, none$mean_PEV, tolerance = 1e-5)
  expect_equal(stiff$CDmean,   none$CDmean,   tolerance = 1e-5)
})

test_that("a flexible P-spline surface absorbs signal and raises PEV", {
  fb <- alpha_grid_fb()
  none  <- ev_alpha(fb, residual_structure = "IID", spatial_random = "none")
  loose <- ev_alpha(fb, residual_structure = "IID", spatial_random = "pspline",
                    spline_knots_row = 5, spline_knots_col = 5,
                    spline_lambda_row = 0.01, spline_lambda_col = 0.01)
  expect_gt(loose$mean_PEV, none$mean_PEV)
})

# ---- New residual structures end-to-end (alpha evaluator) -------------------

test_that("nugget = 0 through the evaluator matches pure AR1xAR1", {
  fb <- alpha_grid_fb()
  pure <- ev_alpha(fb, residual_structure = "AR1xAR1", rho_row = .5, rho_col = .3)
  nug0 <- ev_alpha(fb, residual_structure = "AR1xAR1_nugget",
                   rho_row = .5, rho_col = .3, nugget = 0)
  expect_equal(nug0$mean_PEV, pure$mean_PEV, tolerance = 1e-6)
})

test_that("every new residual structure yields finite, sane criteria", {
  fb <- alpha_grid_fb()
  cases <- list(
    list(residual_structure = "AR1xAR1_nugget", rho_row = .5, rho_col = .3, nugget = .2),
    list(residual_structure = "exponential", kernel_range = 3, nugget = .1),
    list(residual_structure = "gaussian",    kernel_range = 3, nugget = .1),
    list(residual_structure = "matern",      kernel_range = 3, nugget = .1, matern_nu = 1.5),
    list(residual_structure = "matern",      kernel_range = 3, nugget = .1, matern_nu = 2.5)
  )
  for (a in cases) {
    eff <- do.call(ev_alpha, c(list(fb = fb), a))
    lab <- paste(a$residual_structure, a$matern_nu)
    expect_true(is.finite(eff$mean_PEV), info = lab)
    expect_gt(eff$mean_PEV, 0)
    expect_lt(eff$CDmean, 1)
  }
})

test_that("an invalid Matern smoothness is rejected", {
  fb <- alpha_grid_fb()
  expect_error(ev_alpha(fb, residual_structure = "matern",
                        kernel_range = 3, matern_nu = 1.234))
})

# ---- famoptg evaluator with the new options ---------------------------------

test_that("famoptg evaluator accepts the new spatial structures", {
  n_rows <- 6L; n_cols <- 6L; N <- n_rows * n_cols
  fb <- data.frame(
    Plot      = seq_len(N),
    Row       = rep(seq_len(n_rows), times = n_cols),
    Column    = rep(seq_len(n_cols), each = n_rows),
    Block     = rep(1:6, each = N / 6),
    Treatment = rep(paste0("G", 1:9), length.out = N),
    Family    = "F1",
    stringsAsFactors = FALSE
  )
  vc <- list(sigma_e2 = 2, sigma_g2 = 1, sigma_b2 = 0.4,
             sigma_r2 = 0.2, sigma_c2 = 0.2)
  run <- function(...) met_evaluate_famoptg_efficiency(
    field_book = fb, n_rows = n_rows, n_cols = n_cols,
    check_treatments = character(0),
    treatment_effect = "random", prediction_type = "IID",
    varcomp = vc, ...)

  nug <- run(residual_structure = "AR1xAR1_nugget", rho_row = .4, rho_col = .4,
             nugget = .2)
  ker <- run(residual_structure = "exponential", kernel_range = 2, nugget = .1)
  spl <- run(residual_structure = "IID", spatial_random = "pspline",
             spline_knots_row = 5, spline_knots_col = 5,
             spline_lambda_row = 1e8, spline_lambda_col = 1e8)
  none <- run(residual_structure = "IID")

  for (e in list(nug, ker, spl)) expect_true(is.finite(e$mean_PEV))
  # stiff spline again collapses to the no-spline model
  expect_equal(spl$mean_PEV, none$mean_PEV, tolerance = 1e-5)
})

# ---- optimiser accepts the new spatial options ------------------------------

test_that("met_optimize_alpha_rc runs under the new spatial structures", {
  # Reuse the shared, known-valid alpha geometry from helper-test-data.R so the
  # test exercises the spatial options rather than design feasibility.
  fa <- make_alpha_args()
  fa$verbose <- NULL

  run <- function(...) do.call(met_optimize_alpha_rc, c(fa, list(
    treatment_effect = "fixed", method = "RS", n_restarts = 2L,
    criterion = "A", verbose_opt = FALSE, ...)))

  nug <- run(residual_structure = "AR1xAR1_nugget",
             rho_row = 0.4, rho_col = 0.3, nugget = 0.2)
  expect_true(is.finite(nug$optimization$best_score))
  expect_gt(nug$optimization$best_score, 0)

  ker <- run(residual_structure = "exponential", kernel_range = 3, nugget = 0.1)
  expect_true(is.finite(ker$optimization$best_score))
})
