# Tests for the K_ii-scaled per-line reliability / CDmean (P0.2).
#
# Oracle relationship: CD_i = 1 - PEV_i / (sigma_g2 * K_ii)  =>
#   PEV_i = sigma_g2 * K_ii * (1 - CD_i), hence
#   mean_PEV = mean( sigma_g2 * K_ii * (1 - CD_i) ).
# With a constant diagonal K_ii = d this reduces to
#   mean_PEV = sigma_g2 * d * (1 - CDmean).
# The old (buggy) code used K_ii = 1, so this identity holds ONLY when the
# denominator carries K_ii. Testing with d != 1 therefore distinguishes the fix
# from the bug without needing an external PEV oracle.

# Minimal hand-built alpha field book: 4 rows x 3 cols = 12 plots.
# 6 entries (E1..E6, once each) + 2 checks (C1, C2, three times each).
make_fb <- function() {
  trt <- c("C1","E1","E2",
           "C2","E3","E4",
           "C1","E5","E6",
           "C2","C1","C2")
  fb <- data.frame(
    Plot       = seq_len(12),
    Row        = rep(1:4, each = 3),
    Column     = rep(1:3, times = 4),
    Rep        = rep(1:2, each = 6),
    IBlock     = rep(1:4, each = 3),
    BlockInRep = rep(1:2, times = 6),
    Treatment  = trt,
    Check      = trt %in% c("C1","C2"),
    stringsAsFactors = FALSE
  )
  fb
}

make_K <- function(d) {
  ent <- paste0("E", 1:6)
  K <- matrix(0.1, 6, 6, dimnames = list(ent, ent))
  diag(K) <- d
  K
}

vc <- list(sigma_e2 = 1, sigma_g2 = 0.5, sigma_rep2 = 0.2,
           sigma_ib2 = 0.1, sigma_r2 = 0.05, sigma_c2 = 0.05)

test_that("CDmean carries K_ii in the denominator (d = 1.5)", {
  d  <- 1.5
  eff <- met_evaluate_alpha_efficiency(
    field_book = make_fb(), n_rows = 4, n_cols = 3,
    check_treatments = c("C1","C2"),
    treatment_effect = "random", prediction_type = "GBLUP",
    K = make_K(d), varcomp = vc, residual_structure = "IID"
  )
  # Identity that only holds when K_ii is in the denominator:
  expect_equal(eff$mean_PEV, vc$sigma_g2 * d * (1 - eff$CDmean), tolerance = 1e-8)
  # Per-line consistency
  expect_equal(mean(eff$CD_per_line), eff$CDmean, tolerance = 1e-10)
  expect_true(all(eff$CD_per_line <= 1))
})

test_that("backward compatible when diag(K) = 1", {
  eff <- met_evaluate_alpha_efficiency(
    field_book = make_fb(), n_rows = 4, n_cols = 3,
    check_treatments = c("C1","C2"),
    treatment_effect = "random", prediction_type = "GBLUP",
    K = make_K(1), varcomp = vc, residual_structure = "IID"
  )
  expect_equal(eff$CDmean, 1 - eff$mean_PEV / vc$sigma_g2, tolerance = 1e-10)
})

test_that("IID prediction keeps the simple form (K_ii = 1)", {
  eff <- met_evaluate_alpha_efficiency(
    field_book = make_fb(), n_rows = 4, n_cols = 3,
    check_treatments = c("C1","C2"),
    treatment_effect = "random", prediction_type = "IID",
    varcomp = vc, residual_structure = "IID"
  )
  expect_equal(eff$CDmean, 1 - eff$mean_PEV / vc$sigma_g2, tolerance = 1e-10)
})
