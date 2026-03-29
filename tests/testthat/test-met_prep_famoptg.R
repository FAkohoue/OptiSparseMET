# ==============================================================================
# test-met_prep_famoptg.R
# Tests for met_prep_famoptg() -- the MET-context repeated-check block
# design constructor.
# ==============================================================================

# ── Return structure ──────────────────────────────────────────────────────────

test_that("returns list with layout_matrix, field_book, seed_used", {
  args <- make_famoptg_args()
  out  <- do.call(met_prep_famoptg, args)
  expect_true(is.list(out))
  expect_true(is.matrix(out$layout_matrix))
  expect_s3_class(out$field_book, "data.frame")
  expect_length(out$seed_used, 1L)
  expect_false("efficiency" %in% names(out))
})

test_that("layout_matrix has correct dimensions", {
  args <- make_famoptg_args()
  out  <- do.call(met_prep_famoptg, args)
  expect_equal(dim(out$layout_matrix), c(6L, 6L))
})

test_that("field_book has all required columns", {
  args <- make_famoptg_args()
  out  <- do.call(met_prep_famoptg, args)
  expect_true(all(c("Treatment", "Family", "Gcluster",
                    "Block", "Plot", "Row", "Column") %in% names(out$field_book)))
})

test_that("used plots in field_book match non-NA cells in layout_matrix", {
  args <- make_famoptg_args()
  out  <- do.call(met_prep_famoptg, args)
  expect_equal(nrow(out$field_book), sum(!is.na(out$layout_matrix)))
})

# ── Reproducibility ──────────────────────────────────────────────────────────

test_that("same seed produces identical field_book", {
  args <- make_famoptg_args(seed = 42L)
  out1 <- do.call(met_prep_famoptg, args)
  out2 <- do.call(met_prep_famoptg, args)
  expect_identical(out1$field_book, out2$field_book)
})

test_that("different seeds produce different layouts", {
  out1 <- do.call(met_prep_famoptg, make_famoptg_args(seed = 1L))
  out2 <- do.call(met_prep_famoptg, make_famoptg_args(seed = 2L))
  expect_false(identical(out1$field_book$Treatment, out2$field_book$Treatment))
})

# ── Block structure constraints ───────────────────────────────────────────────

test_that("check treatments appear in every block", {
  args <- make_famoptg_args()
  out  <- do.call(met_prep_famoptg, args)
  fb   <- out$field_book
  checks <- c("CHK1", "CHK2")
  for (chk in checks) {
    blocks_with_chk <- unique(fb$Block[fb$Treatment == chk])
    expect_equal(sort(blocks_with_chk), 1:3)
  }
})

test_that("p-rep treatments appear in exactly p_rep_reps blocks", {
  args <- make_famoptg_args()
  out  <- do.call(met_prep_famoptg, args)
  fb   <- out$field_book
  for (i in seq_along(args$p_rep_treatments)) {
    trt    <- args$p_rep_treatments[i]
    expect_equal(sum(fb$Treatment == trt), args$p_rep_reps[i])
  }
})

test_that("p-rep treatments never appear twice in the same block", {
  args <- make_famoptg_args()
  out  <- do.call(met_prep_famoptg, args)
  fb   <- out$field_book
  for (trt in args$p_rep_treatments) {
    trt_blocks <- fb$Block[fb$Treatment == trt]
    expect_equal(length(trt_blocks), length(unique(trt_blocks)))
  }
})

test_that("unreplicated treatments appear exactly once", {
  args <- make_famoptg_args()
  out  <- do.call(met_prep_famoptg, args)
  fb   <- out$field_book
  for (trt in args$unreplicated_treatments) {
    expect_equal(sum(fb$Treatment == trt), 1L)
  }
})

# ── Augmented design class ───────────────────────────────────────────────────

test_that("augmented design: all non-checks unreplicated", {
  args <- make_famoptg_args()
  args$p_rep_treatments <- character(0)
  args$p_rep_reps       <- integer(0)
  args$p_rep_families   <- character(0)
  args$unreplicated_treatments <- paste0("U", 1:24)
  args$unreplicated_families   <- rep("F1", 24)
  args$n_rows <- 6L; args$n_cols <- 6L  # 2*3 + 24 = 30 -> adjust
  args$n_rows <- 5L; args$n_cols <- 6L  # 30 exactly
  out <- do.call(met_prep_famoptg, args)
  fb  <- out$field_book
  non_check <- fb[!fb$Treatment %in% c("CHK1", "CHK2"), ]
  counts    <- table(non_check$Treatment)
  expect_true(all(counts == 1L))
})

# ── Gcluster column ──────────────────────────────────────────────────────────

test_that("Gcluster is NA for all plots when cluster_source = 'Family'", {
  args <- make_famoptg_args()
  args$cluster_source <- "Family"
  out  <- do.call(met_prep_famoptg, args)
  expect_true(all(is.na(out$field_book$Gcluster)))
})

# ── check_placement modes ─────────────────────────────────────────────────────

test_that("systematic check_placement runs without error", {
  args <- make_famoptg_args()
  args$check_placement <- "systematic"
  expect_no_error(do.call(met_prep_famoptg, args))
})

test_that("optimal check_placement runs without error", {
  args <- make_famoptg_args()
  args$check_placement    <- "optimal"
  args$check_opt_attempts <- 5L
  expect_no_error(do.call(met_prep_famoptg, args))
})

# ── Validation errors ─────────────────────────────────────────────────────────

test_that("stops when check_treatments overlap with p_rep_treatments", {
  args <- make_famoptg_args()
  args$p_rep_treatments[1] <- "CHK1"
  expect_error(do.call(met_prep_famoptg, args))
})

test_that("stops when p_rep_reps exceeds n_blocks", {
  args <- make_famoptg_args()
  args$p_rep_reps[1] <- 99L
  expect_error(do.call(met_prep_famoptg, args))
})

test_that("stops when check_families length mismatches check_treatments", {
  args <- make_famoptg_args()
  args$check_families <- c("CHECK")  # length 1, treatments length 2
  expect_error(do.call(met_prep_famoptg, args))
})

test_that("stops when duplicate check_treatments supplied", {
  args <- make_famoptg_args()
  args$check_treatments <- c("CHK1", "CHK1")
  args$check_families   <- c("CHECK", "CHECK")
  expect_error(do.call(met_prep_famoptg, args))
})
