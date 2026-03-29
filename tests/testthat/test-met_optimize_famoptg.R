# ==============================================================================
# test-met_optimize_famoptg.R
# Tests for met_optimize_famoptg() -- RS optimizer for met_prep_famoptg designs.
# ==============================================================================

base_opt_args <- function(...) {
  fa        <- make_famoptg_args()
  fa$verbose <- NULL           # met_optimize_famoptg has no verbose parameter
  overrides <- list(...)
  defaults  <- list(
    treatment_effect   = "fixed",
    residual_structure = "IID",
    criterion          = "A",
    n_restarts         = 3L,
    verbose_opt        = FALSE
  )
  c(fa, modifyList(defaults, overrides))
}

# ── Return structure ──────────────────────────────────────────────────────────

test_that("returns met_prep_famoptg-compatible list plus efficiency and optimization slots", {
  out <- do.call(met_optimize_famoptg, base_opt_args())
  expect_true(is.list(out))
  expect_true(is.matrix(out$layout_matrix))
  expect_s3_class(out$field_book, "data.frame")
  expect_true(is.list(out$efficiency))
  expect_true(is.list(out$optimization))
})

test_that("optimization slot has required fields", {
  out <- do.call(met_optimize_famoptg, base_opt_args())
  opt <- out$optimization
  for (f in c("method", "criterion", "best_score",
              "score_history", "master_seed", "n_restarts", "n_failed")) {
    expect_true(f %in% names(opt), info = f)
  }
})

test_that("method is always RS", {
  out <- do.call(met_optimize_famoptg, base_opt_args())
  expect_equal(out$optimization$method, "RS")
})

test_that("score_history length equals n_restarts", {
  out <- do.call(met_optimize_famoptg, base_opt_args(n_restarts = 5L))
  expect_length(out$optimization$score_history, 5L)
})

# ── Criterion directions ──────────────────────────────────────────────────────

test_that("A criterion: best_score <= all non-NA score_history values", {
  out  <- do.call(met_optimize_famoptg, base_opt_args(n_restarts = 5L))
  hist <- out$optimization$score_history
  expect_true(out$optimization$best_score <= min(hist, na.rm = TRUE) + 1e-10)
})

test_that("CDmean criterion: best_score is positive and <= 1", {
  out <- do.call(met_optimize_famoptg,
                 base_opt_args(
                   treatment_effect = "random",
                   prediction_type  = "IID",
                   criterion        = "CDmean",
                   varcomp          = famoptg_varcomp()
                 ))
  expect_gte(out$optimization$best_score, 0)
  expect_lte(out$optimization$best_score, 1)
})

test_that("CDmean score_history values are all positive when not NA", {
  out <- do.call(met_optimize_famoptg,
                 base_opt_args(
                   treatment_effect = "random",
                   prediction_type  = "IID",
                   criterion        = "CDmean",
                   varcomp          = famoptg_varcomp(),
                   n_restarts       = 4L
                 ))
  hist <- out$optimization$score_history[!is.na(out$optimization$score_history)]
  expect_true(all(hist >= 0))
})

# ── Structural guarantees preserved by optimizer ──────────────────────────────

test_that("checks appear in every block of the optimised design", {
  out <- do.call(met_optimize_famoptg, base_opt_args())
  fb  <- out$field_book
  for (chk in c("CHK1", "CHK2")) {
    expect_equal(sort(unique(fb$Block[fb$Treatment == chk])), 1:3)
  }
})

test_that("p-rep treatments never duplicated within a block", {
  out  <- do.call(met_optimize_famoptg, base_opt_args())
  fb   <- out$field_book
  args <- make_famoptg_args()
  for (trt in args$p_rep_treatments) {
    trt_blocks <- fb$Block[fb$Treatment == trt]
    expect_equal(length(trt_blocks), length(unique(trt_blocks)))
  }
})

# ── Pre-flight errors ─────────────────────────────────────────────────────────

test_that("stops when CDmean used with fixed effects", {
  expect_error(
    do.call(met_optimize_famoptg,
            base_opt_args(treatment_effect = "fixed", criterion = "CDmean"))
  )
})

test_that("stops when random effects used with prediction_type = none", {
  expect_error(
    do.call(met_optimize_famoptg,
            base_opt_args(treatment_effect = "random",
                          prediction_type  = "none",
                          criterion        = "A"))
  )
})

test_that("stops when n_restarts < 1", {
  expect_error(
    do.call(met_optimize_famoptg, base_opt_args(n_restarts = 0L))
  )
})