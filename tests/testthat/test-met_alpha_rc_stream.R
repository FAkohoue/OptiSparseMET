# ==============================================================================
# test-met_alpha_rc_stream.R
# Tests for met_alpha_rc_stream() -- the MET-context alpha row-column
# stream constructor.
# ==============================================================================

# ── Return structure ──────────────────────────────────────────────────────────

test_that("returns list with layout_matrix, field_book, design_info, seed_used", {
  out <- do.call(met_alpha_rc_stream, make_alpha_args())
  expect_true(is.list(out))
  expect_true(is.matrix(out$layout_matrix))
  expect_s3_class(out$field_book, "data.frame")
  expect_true(is.list(out$design_info))
  expect_length(out$seed_used, 1L)
  expect_false("efficiency" %in% names(out))
})

test_that("layout_matrix has correct dimensions", {
  out <- do.call(met_alpha_rc_stream, make_alpha_args())
  expect_equal(dim(out$layout_matrix), c(6L, 10L))
})

test_that("field_book has all required columns", {
  out <- do.call(met_alpha_rc_stream, make_alpha_args())
  expect_true(all(c("Plot", "Row", "Column", "Rep", "IBlock",
                    "BlockInRep", "Treatment", "Family",
                    "Gcluster", "Check") %in% names(out$field_book)))
})

test_that("total_used_plots in design_info matches non-NA layout cells", {
  out <- do.call(met_alpha_rc_stream, make_alpha_args())
  expect_equal(
    sum(!is.na(out$layout_matrix)),
    out$design_info$total_used_plots
  )
})

# ── Reproducibility ──────────────────────────────────────────────────────────

test_that("same seed produces identical field_book", {
  out1 <- do.call(met_alpha_rc_stream, make_alpha_args(seed = 7L))
  out2 <- do.call(met_alpha_rc_stream, make_alpha_args(seed = 7L))
  expect_identical(out1$field_book, out2$field_book)
})

# ── Design structure constraints ──────────────────────────────────────────────

test_that("each entry appears exactly once per replicate", {
  out <- do.call(met_alpha_rc_stream, make_alpha_args())
  fb  <- out$field_book[!is.na(out$field_book$Treatment), ]
  entry_fb <- fb[!fb$Check, ]
  for (rep_id in unique(entry_fb$Rep)) {
    rep_trts <- entry_fb$Treatment[entry_fb$Rep == rep_id]
    expect_equal(length(rep_trts), length(unique(rep_trts)))
  }
})

test_that("check treatments appear in every incomplete block", {
  out    <- do.call(met_alpha_rc_stream, make_alpha_args())
  fb     <- out$field_book[!is.na(out$field_book$Treatment), ]
  checks <- c("CHK1", "CHK2", "CHK3")
  for (blk in unique(fb$IBlock[!is.na(fb$IBlock)])) {
    blk_checks <- fb$Treatment[fb$IBlock == blk & fb$Check]
    expect_true(all(checks %in% blk_checks))
  }
})

test_that("n_reps in design_info matches argument", {
  out <- do.call(met_alpha_rc_stream, make_alpha_args())
  expect_equal(out$design_info$n_reps, 2L)
})

test_that("block sizes respect min_block_size and max_block_size", {
  out <- do.call(met_alpha_rc_stream, make_alpha_args())
  fb  <- out$field_book[!is.na(out$field_book$Treatment), ]
  block_sizes <- as.integer(table(fb$IBlock))
  expect_true(all(block_sizes >= 8L & block_sizes <= 12L))
})

# ── Gcluster column ──────────────────────────────────────────────────────────

test_that("Gcluster is NA for all plots when cluster_source = 'Family'", {
  args <- make_alpha_args()
  args$cluster_source <- "Family"
  out  <- do.call(met_alpha_rc_stream, args)
  expect_true(all(is.na(out$field_book$Gcluster)))
})

test_that("check plots always have Gcluster = NA", {
  args <- make_alpha_args()
  out  <- do.call(met_alpha_rc_stream, args)
  fb   <- out$field_book[!is.na(out$field_book$Treatment), ]
  expect_true(all(is.na(fb$Gcluster[fb$Check])))
})

# ── Block-size parameterisation ───────────────────────────────────────────────

test_that("user-fixed n_blocks_per_rep respected", {
  args <- make_alpha_args()
  args$n_blocks_per_rep <- 2L
  args$min_block_size   <- NULL
  args$max_block_size   <- NULL
  out <- do.call(met_alpha_rc_stream, args)
  expect_equal(out$design_info$n_blocks_per_rep, 2L)
})

test_that("stops when n_blocks_per_rep infeasible", {
  args <- make_alpha_args()
  args$n_blocks_per_rep <- 99L
  expect_error(do.call(met_alpha_rc_stream, args))
})

# ── Validation errors ─────────────────────────────────────────────────────────

test_that("stops when check_treatments overlap with entry_treatments", {
  args <- make_alpha_args()
  args$entry_treatments[1] <- "CHK1"
  expect_error(do.call(met_alpha_rc_stream, args))
})

test_that("stops when min_block_size > max_block_size", {
  args <- make_alpha_args()
  args$min_block_size <- 15L
  args$max_block_size <- 10L
  expect_error(do.call(met_alpha_rc_stream, args))
})

test_that("stops when field too small for required plots", {
  args <- make_alpha_args()
  args$n_rows <- 2L
  args$n_cols <- 2L
  expect_error(do.call(met_alpha_rc_stream, args))
})