# ==============================================================================
# test-met_evaluate_famoptg_efficiency.R
# Tests for met_evaluate_famoptg_efficiency() -- uses sigma_b2 varcomp.
# ==============================================================================

make_design_famoptg <- function(seed = 1L) {
  do.call(met_prep_famoptg, make_famoptg_args(seed = seed))
}

eval_famoptg <- function(design, ...) {
  met_evaluate_famoptg_efficiency(
    field_book       = design$field_book,
    n_rows           = 6L, n_cols = 6L,
    check_treatments = c("CHK1", "CHK2"),
    ...
  )
}

# ── Return structure ──────────────────────────────────────────────────────────

test_that("returns a list with required fields", {
  eff <- eval_famoptg(make_design_famoptg(), treatment_effect = "fixed")
  expect_true(is.list(eff))
  for (f in c("A_criterion", "A_efficiency", "mode",
              "treatment_effect", "residual_structure_used")) {
    expect_true(f %in% names(eff), info = f)
  }
})

test_that("mode is FIXED_TREATMENT_BLUE_CONTRAST for fixed effects", {
  eff <- eval_famoptg(make_design_famoptg(), treatment_effect = "fixed")
  expect_equal(eff$mode, "FIXED_TREATMENT_BLUE_CONTRAST")
})

# ── Fixed effects: A and D criteria ──────────────────────────────────────────

test_that("A_efficiency = 1 / A_criterion (mathematical identity)", {
  eff <- eval_famoptg(make_design_famoptg(), treatment_effect = "fixed")
  expect_equal(eff$A_efficiency, 1 / eff$A_criterion, tolerance = 1e-10)
})

test_that("D_criterion is positive and D_efficiency = 1 / D_criterion", {
  eff <- eval_famoptg(make_design_famoptg(), treatment_effect = "fixed")
  expect_gt(eff$D_criterion, 0)
  expect_equal(eff$D_efficiency, 1 / eff$D_criterion, tolerance = 1e-10)
})

test_that("A and D aliases match efficiency fields", {
  eff <- eval_famoptg(make_design_famoptg(), treatment_effect = "fixed")
  expect_equal(eff$A, eff$A_efficiency)
  expect_equal(eff$D, eff$D_efficiency)
})

# ── Random effects: CDmean ────────────────────────────────────────────────────

test_that("CDmean computed for random effects", {
  eff <- eval_famoptg(make_design_famoptg(),
                      treatment_effect = "random", prediction_type = "IID")
  expect_true(!is.null(eff$CDmean))
  expect_gte(eff$CDmean, 0)
  expect_lte(eff$CDmean, 1)
})

test_that("CDmean = 1 - mean_PEV / sigma_g2 (mathematical identity)", {
  vc  <- famoptg_varcomp()
  eff <- eval_famoptg(make_design_famoptg(),
                      treatment_effect = "random", prediction_type = "IID",
                      varcomp = vc)
  expect_equal(eff$CDmean, 1 - eff$mean_PEV / vc$sigma_g2, tolerance = 1e-10)
})

test_that("mode is RANDOM_TREATMENT_PEV for random effects", {
  eff <- eval_famoptg(make_design_famoptg(),
                      treatment_effect = "random", prediction_type = "IID")
  expect_equal(eff$mode, "RANDOM_TREATMENT_PEV")
})

test_that("D_criterion is NA for random effects", {
  eff <- eval_famoptg(make_design_famoptg(),
                      treatment_effect = "random", prediction_type = "IID")
  expect_true(is.na(eff$D_criterion))
})

# ── Residual structures ───────────────────────────────────────────────────────

test_that("AR1xAR1 gives different A_criterion than IID", {
  eff_iid <- eval_famoptg(make_design_famoptg(),
                          treatment_effect = "fixed", residual_structure = "IID")
  eff_ar1 <- eval_famoptg(make_design_famoptg(),
                          treatment_effect = "fixed",
                          residual_structure = "AR1xAR1",
                          rho_row = 0.3, rho_col = 0.2)
  expect_false(isTRUE(all.equal(eff_iid$A_criterion, eff_ar1$A_criterion)))
})

# ── sigma_b2 requirement ──────────────────────────────────────────────────────

test_that("stops when sigma_b2 missing from varcomp", {
  vc <- famoptg_varcomp()
  vc$sigma_b2 <- NULL
  expect_error(
    eval_famoptg(make_design_famoptg(), treatment_effect = "fixed", varcomp = vc)
  )
})

test_that("stops when sigma_rep2 supplied instead of sigma_b2", {
  vc <- famoptg_varcomp()
  vc$sigma_b2  <- NULL
  vc$sigma_rep2 <- 0.5
  expect_error(
    eval_famoptg(make_design_famoptg(), treatment_effect = "fixed", varcomp = vc)
  )
})

# ── Validation errors ─────────────────────────────────────────────────────────

test_that("stops when required field_book columns missing", {
  d  <- make_design_famoptg()
  fb <- d$field_book
  fb$Block <- NULL
  expect_error(
    met_evaluate_famoptg_efficiency(
      field_book = fb, n_rows = 6L, n_cols = 6L,
      check_treatments = c("CHK1", "CHK2"),
      treatment_effect = "fixed"
    )
  )
})

test_that("stops when treatment_effect = random and prediction_type = none", {
  expect_error(
    eval_famoptg(make_design_famoptg(),
                 treatment_effect = "random", prediction_type = "none")
  )
})

test_that("stops when rho_row outside (-1, 1)", {
  expect_error(
    eval_famoptg(make_design_famoptg(),
                 treatment_effect = "fixed",
                 residual_structure = "AR1xAR1",
                 rho_row = 1.5, rho_col = 0.3)
  )
})