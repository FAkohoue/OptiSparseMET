# Test for the fixed-effects rank-deficiency guard (P0.4).
#
# With treatment_effect = "fixed", check_as_fixed = FALSE, and no check rows, the
# intercept equals the sum of the full Line dummy set, so X'QX is singular and
# the previous code errored inside .solve_C. The guard drops the intercept and
# the evaluation proceeds.

make_all_entry_fb <- function() {
  trt <- rep(paste0("E", 1:6), times = 2)   # 6 entries, each twice = 12 plots
  data.frame(
    Plot       = seq_len(12),
    Row        = rep(1:4, each = 3),
    Column     = rep(1:3, times = 4),
    Rep        = rep(1:2, each = 6),
    IBlock     = rep(1:4, each = 3),
    BlockInRep = rep(1:2, times = 6),
    Treatment  = trt,
    Check      = FALSE,
    stringsAsFactors = FALSE
  )
}

vc <- list(sigma_e2 = 1, sigma_g2 = 1, sigma_rep2 = 0.2,
           sigma_ib2 = 0.1, sigma_r2 = 0.05, sigma_c2 = 0.05)

test_that("singular fixed design is repaired instead of erroring", {
  expect_message(
    eff <- met_evaluate_alpha_efficiency(
      field_book = make_all_entry_fb(), n_rows = 4, n_cols = 3,
      check_treatments = character(0),
      treatment_effect = "fixed", check_as_fixed = FALSE,
      varcomp = vc, residual_structure = "IID"
    ),
    regexp = "rank-deficient"
  )
  expect_true(is.finite(eff$A_criterion))
  expect_gt(eff$A_criterion, 0)
})
