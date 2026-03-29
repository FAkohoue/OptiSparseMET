# ==============================================================================
# test-met_evaluate_alpha_efficiency.R
# Tests for met_evaluate_alpha_efficiency() -- uses sigma_rep2 + sigma_ib2.
# Key distinction from met_evaluate_famoptg_efficiency(): field_book must have
# Rep/IBlock columns and varcomp must contain sigma_rep2 and sigma_ib2.
# ==============================================================================

make_design_alpha <- function(seed = 1L) {
  do.call(met_alpha_rc_stream, make_alpha_args(seed = seed))
}

eval_alpha <- function(design, ...) {
  met_evaluate_alpha_efficiency(
    field_book       = design$field_book,
    n_rows           = 6L, n_cols = 10L,
    check_treatments = c("CHK1", "CHK2", "CHK3"),
    ...
  )
}

# ── Return structure ──────────────────────────────────────────────────────────

test_that("returns a list with required fields", {
  eff <- eval_alpha(make_design_alpha(), treatment_effect = "fixed")
  for (f in c("A_criterion", "A_efficiency", "D_criterion", "D_efficiency",
              "mode", "treatment_effect", "residual_structure_used")) {
    expect_true(f %in% names(eff), info = f)
  }
})

test_that("model string contains Rep and IBlock (not just Block)", {
  eff <- eval_alpha(make_design_alpha(), treatment_effect = "fixed")
  expect_true(grepl("Rep", eff$model))
  expect_true(grepl("IBlock", eff$model))
})

test_that("mode is FIXED_TREATMENT_BLUE_CONTRAST for fixed effects", {
  eff <- eval_alpha(make_design_alpha(), treatment_effect = "fixed")
  expect_equal(eff$mode, "FIXED_TREATMENT_BLUE_CONTRAST")
})

# ── Fixed effects: A and D criteria ──────────────────────────────────────────

test_that("A_efficiency = 1 / A_criterion", {
  eff <- eval_alpha(make_design_alpha(), treatment_effect = "fixed")
  expect_equal(eff$A_efficiency, 1 / eff$A_criterion, tolerance = 1e-10)
})

test_that("D_criterion is positive and D_efficiency = 1 / D_criterion", {
  eff <- eval_alpha(make_design_alpha(), treatment_effect = "fixed")
  expect_gt(eff$D_criterion, 0)
  expect_equal(eff$D_efficiency, 1 / eff$D_criterion, tolerance = 1e-10)
})

test_that("A and D aliases match efficiency fields", {
  eff <- eval_alpha(make_design_alpha(), treatment_effect = "fixed")
  expect_equal(eff$A, eff$A_efficiency)
  expect_equal(eff$D, eff$D_efficiency)
})

# ── Random effects: CDmean ────────────────────────────────────────────────────

test_that("CDmean in [0, 1] for random effects IID", {
  eff <- eval_alpha(make_design_alpha(),
                    treatment_effect = "random", prediction_type = "IID")
  expect_gte(eff$CDmean, 0)
  expect_lte(eff$CDmean, 1)
})

test_that("CDmean = 1 - mean_PEV / sigma_g2", {
  vc  <- alpha_varcomp()
  eff <- eval_alpha(make_design_alpha(),
                    treatment_effect = "random", prediction_type = "IID",
                    varcomp = vc)
  expect_equal(eff$CDmean, 1 - eff$mean_PEV / vc$sigma_g2, tolerance = 1e-10)
})

test_that("D_criterion is NA for random effects", {
  eff <- eval_alpha(make_design_alpha(),
                    treatment_effect = "random", prediction_type = "IID")
  expect_true(is.na(eff$D_criterion))
})

# ── sigma_rep2 + sigma_ib2 requirement ───────────────────────────────────────

test_that("stops when sigma_rep2 missing from varcomp", {
  vc <- alpha_varcomp()
  vc$sigma_rep2 <- NULL
  expect_error(
    eval_alpha(make_design_alpha(), treatment_effect = "fixed", varcomp = vc)
  )
})

test_that("stops when sigma_ib2 missing from varcomp", {
  vc <- alpha_varcomp()
  vc$sigma_ib2 <- NULL
  expect_error(
    eval_alpha(make_design_alpha(), treatment_effect = "fixed", varcomp = vc)
  )
})

test_that("stops when sigma_b2 supplied instead of sigma_rep2 + sigma_ib2", {
  vc <- list(sigma_e2 = 1, sigma_g2 = 1, sigma_b2 = 1,
             sigma_r2 = 0.1, sigma_c2 = 0.1)
  expect_error(
    eval_alpha(make_design_alpha(), treatment_effect = "fixed", varcomp = vc)
  )
})

# ── Residual structures ───────────────────────────────────────────────────────

test_that("AR1xAR1 gives different A_criterion than IID", {
  eff_iid <- eval_alpha(make_design_alpha(),
                        treatment_effect = "fixed", residual_structure = "IID")
  eff_ar1 <- eval_alpha(make_design_alpha(),
                        treatment_effect = "fixed",
                        residual_structure = "AR1xAR1",
                        rho_row = 0.3, rho_col = 0.2)
  expect_false(isTRUE(all.equal(eff_iid$A_criterion, eff_ar1$A_criterion)))
})

# ── Validation errors ─────────────────────────────────────────────────────────

test_that("stops when required field_book columns missing", {
  d  <- make_design_alpha()
  fb <- d$field_book
  fb$Rep <- NULL
  expect_error(
    met_evaluate_alpha_efficiency(
      field_book = fb, n_rows = 6L, n_cols = 10L,
      check_treatments = c("CHK1", "CHK2", "CHK3"),
      treatment_effect = "fixed"
    )
  )
})

test_that("stops when treatment_effect = random and prediction_type = none", {
  expect_error(
    eval_alpha(make_design_alpha(),
               treatment_effect = "random", prediction_type = "none")
  )
})