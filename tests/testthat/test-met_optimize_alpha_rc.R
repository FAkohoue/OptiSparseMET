# ==============================================================================
# test-met_optimize_alpha_rc.R
# Tests for met_optimize_alpha_rc() -- RS/SA/GA optimizer for
# met_alpha_rc_stream designs.
# ==============================================================================

base_opt_alpha <- function(method = "RS", criterion = "A", ...) {
  fa        <- make_alpha_args()
  fa$verbose <- NULL           # met_optimize_alpha_rc has no verbose parameter
  overrides <- list(...)
  defaults  <- list(
    treatment_effect   = "fixed",
    residual_structure = "IID",
    method             = method,
    criterion          = criterion,
    n_restarts         = 2L,
    verbose_opt        = FALSE
  )
  c(fa, modifyList(defaults, overrides))
}

# ── Return structure ──────────────────────────────────────────────────────────

test_that("RS: returns met_alpha_rc_stream list plus efficiency and optimization", {
  out <- do.call(met_optimize_alpha_rc, base_opt_alpha("RS"))
  expect_true(is.list(out))
  expect_true(is.matrix(out$layout_matrix))
  expect_s3_class(out$field_book, "data.frame")
  expect_true(is.list(out$efficiency))
  expect_true(is.list(out$optimization))
})

test_that("optimization slot has required fields for all methods", {
  for (m in c("RS", "SA", "GA")) {
    args <- base_opt_alpha(m)
    if (m == "SA") {
      args$sa_max_iter   <- 5L
      args$sa_temp_start <- 1.0
      args$sa_temp_end   <- 0.01
    }
    if (m == "GA") {
      args$ga_pop_size      <- 4L
      args$ga_n_generations <- 3L
    }
    out <- do.call(met_optimize_alpha_rc, args)
    opt <- out$optimization
    for (f in c("method", "criterion", "best_score",
                "score_history", "master_seed", "n_failed")) {
      expect_true(f %in% names(opt), info = paste(m, f))
    }
    expect_equal(opt$method, m)
  }
})

# ── Method-specific params slots ──────────────────────────────────────────────

test_that("RS: sa_params and ga_params are NULL", {
  out <- do.call(met_optimize_alpha_rc, base_opt_alpha("RS"))
  expect_null(out$optimization$sa_params)
  expect_null(out$optimization$ga_params)
  expect_equal(out$optimization$n_restarts, 2L)
})

test_that("SA: sa_params present, ga_params NULL", {
  args <- base_opt_alpha("SA")
  args$sa_max_iter <- 5L; args$sa_temp_start <- 1.0; args$sa_temp_end <- 0.01
  out <- do.call(met_optimize_alpha_rc, args)
  expect_true(is.list(out$optimization$sa_params))
  expect_null(out$optimization$ga_params)
})

test_that("GA: ga_params present, sa_params NULL", {
  args <- base_opt_alpha("GA")
  args$ga_pop_size <- 4L; args$ga_n_generations <- 3L
  out <- do.call(met_optimize_alpha_rc, args)
  expect_true(is.list(out$optimization$ga_params))
  expect_null(out$optimization$sa_params)
})

# ── SA score_history is a list of vectors ─────────────────────────────────────

test_that("SA: score_history is a list of numeric vectors", {
  args <- base_opt_alpha("SA")
  args$sa_max_iter <- 5L; args$sa_temp_start <- 1.0; args$sa_temp_end <- 0.01
  out <- do.call(met_optimize_alpha_rc, args)
  expect_true(is.list(out$optimization$score_history))
})

# ── Criterion directions ──────────────────────────────────────────────────────

test_that("A criterion RS: best_score <= all non-NA history values", {
  out  <- do.call(met_optimize_alpha_rc, base_opt_alpha("RS", "A", n_restarts = 4L))
  hist <- out$optimization$score_history
  expect_true(out$optimization$best_score <= min(hist, na.rm = TRUE) + 1e-10)
})

test_that("CDmean RS: best_score in [0, 1]", {
  args <- base_opt_alpha("RS", "CDmean")
  args$treatment_effect <- "random"
  args$prediction_type  <- "IID"
  args$varcomp          <- alpha_varcomp()
  out  <- do.call(met_optimize_alpha_rc, args)
  expect_gte(out$optimization$best_score, 0)
  expect_lte(out$optimization$best_score, 1)
})

test_that("CDmean score_history values are >= 0 when not NA", {
  args <- base_opt_alpha("RS", "CDmean", n_restarts = 3L)
  args$treatment_effect <- "random"
  args$prediction_type  <- "IID"
  args$varcomp          <- alpha_varcomp()
  out  <- do.call(met_optimize_alpha_rc, args)
  hist <- out$optimization$score_history
  expect_true(all(hist[!is.na(hist)] >= 0))
})

# ── Design constraints preserved ──────────────────────────────────────────────

test_that("entries appear exactly once per replicate in optimised design", {
  out     <- do.call(met_optimize_alpha_rc, base_opt_alpha("RS"))
  fb      <- out$field_book[!is.na(out$field_book$Treatment), ]
  entry_fb <- fb[!fb$Check, ]
  for (rep_id in unique(entry_fb$Rep)) {
    trts <- entry_fb$Treatment[entry_fb$Rep == rep_id]
    expect_equal(length(trts), length(unique(trts)))
  }
})

# ── Pre-flight errors ─────────────────────────────────────────────────────────

test_that("stops when CDmean used with fixed effects", {
  expect_error(
    do.call(met_optimize_alpha_rc,
            base_opt_alpha("RS", "CDmean", treatment_effect = "fixed"))
  )
})

test_that("stops when sa_temp_end >= sa_temp_start", {
  args <- base_opt_alpha("SA")
  args$sa_max_iter <- 5L; args$sa_temp_start <- 0.5; args$sa_temp_end <- 1.0
  expect_error(do.call(met_optimize_alpha_rc, args))
})

test_that("stops when ga_pop_size < 4", {
  args <- base_opt_alpha("GA")
  args$ga_pop_size <- 2L; args$ga_n_generations <- 3L
  expect_error(do.call(met_optimize_alpha_rc, args))
})

test_that("stops when ga_elitism >= ga_pop_size", {
  args <- base_opt_alpha("GA")
  args$ga_pop_size <- 4L; args$ga_n_generations <- 3L; args$ga_elitism <- 4L
  expect_error(do.call(met_optimize_alpha_rc, args))
})