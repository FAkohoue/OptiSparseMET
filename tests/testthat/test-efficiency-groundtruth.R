# Ground-truth efficiency tests (P2.1 / P2.2): compare the evaluator against
# INDEPENDENT oracles, not against its own formulas.
#
# Setup: a balanced "CRD" -- t treatments each replicated n times, all plots in
# a single field cell (Row = Column = Rep = IBlock = 1). The nuisance random
# effects then reduce to single-level random intercepts that cancel in treatment
# contrasts, giving closed-form values:
#   A_criterion (fixed) = 2 * sigma_e2 / n     (mean pairwise contrast variance)
#   D_criterion (fixed) =     sigma_e2 / n     (geometric mean of contrast eigenvalues)
# Both are independent of the nuisance variances. Verified by an explicit dense
# mixed-model-equations solve in numpy (see project notes).

single_cell_fb <- function(t, n) {
  trt <- rep(paste0("G", seq_len(t)), each = n)
  N   <- t * n
  data.frame(
    Plot = seq_len(N), Row = 1L, Column = 1L,
    Rep = 1L, IBlock = 1L, BlockInRep = 1L,
    Treatment = trt, Check = FALSE, stringsAsFactors = FALSE
  )
}

test_that("A- and D-criteria match the balanced-CRD closed form (fixed effects)", {
  t <- 4L; n <- 3L; sigma_e2 <- 2
  vc <- list(sigma_e2 = sigma_e2, sigma_g2 = 1,
             sigma_rep2 = 0.3, sigma_ib2 = 0.7, sigma_r2 = 0.5, sigma_c2 = 0.9)
  eff <- met_evaluate_alpha_efficiency(
    field_book = single_cell_fb(t, n), n_rows = 1, n_cols = 1,
    check_treatments = character(0),
    treatment_effect = "fixed", varcomp = vc, residual_structure = "IID"
  )
  expect_equal(eff$A_criterion, 2 * sigma_e2 / n, tolerance = 1e-8)  # 4/3
  expect_equal(eff$D_criterion,     sigma_e2 / n, tolerance = 1e-8)  # 2/3
})

test_that("relative A-efficiency is exactly 1 for the orthogonal balanced design", {
  # A_efficiency_rel compares the design against an orthogonal CRD using the same
  # total plots: (2 * sigma_e2 * p / nn) / A_criterion. A balanced CRD IS that
  # reference design, so the ratio must be exactly 1. It is also the inverse
  # criterion's honest counterpart: A_efficiency is 1 / A_criterion, not a
  # relative efficiency.
  t <- 4L; n <- 3L; sigma_e2 <- 2
  vc <- list(sigma_e2 = sigma_e2, sigma_g2 = 1,
             sigma_rep2 = 0.3, sigma_ib2 = 0.7, sigma_r2 = 0.5, sigma_c2 = 0.9)
  eff <- met_evaluate_alpha_efficiency(
    field_book = single_cell_fb(t, n), n_rows = 1, n_cols = 1,
    check_treatments = character(0),
    treatment_effect = "fixed", varcomp = vc, residual_structure = "IID"
  )
  expect_equal(eff$A_efficiency_rel, 1, tolerance = 1e-8)
  # and the inverse-criterion field remains what it claims to be
  expect_equal(eff$A_efficiency, 1 / eff$A_criterion, tolerance = 1e-10)
})

test_that("random-IID mean PEV / CDmean match the closed-form MME value", {
  # Single-cell balanced design: 4 genotypes x 3 reps, all nuisance factors
  # single-level. The genotype PEV is 0.55 and CDmean 0.45 (verified by an
  # independent dense MME solve; genotype PEV diagonal = c(0.55, 0.55, 0.55, 0.55)).
  # This also guards the fixed-effects offset in the PEV extraction: the genotype
  # block sits AFTER the fixed block in the coefficient matrix.
  vc <- list(sigma_e2 = 2, sigma_g2 = 1,
             sigma_rep2 = 1, sigma_ib2 = 1, sigma_r2 = 1, sigma_c2 = 1)
  eff <- met_evaluate_alpha_efficiency(
    field_book = single_cell_fb(4L, 3L), n_rows = 1, n_cols = 1,
    check_treatments = character(0),
    treatment_effect = "random", prediction_type = "IID",
    varcomp = vc, residual_structure = "IID")
  expect_equal(eff$mean_PEV, 0.55, tolerance = 1e-6)
  expect_equal(eff$CDmean,   0.45, tolerance = 1e-6)
  expect_equal(eff$CD_per_line, rep(0.45, 4), tolerance = 1e-6)
})

test_that("more replication strictly lowers random-IID mean PEV", {
  vc <- list(sigma_e2 = 2, sigma_g2 = 1,
             sigma_rep2 = 1, sigma_ib2 = 1, sigma_r2 = 1, sigma_c2 = 1)
  run <- function(n) met_evaluate_alpha_efficiency(
    field_book = single_cell_fb(4L, n), n_rows = 1, n_cols = 1,
    check_treatments = character(0),
    treatment_effect = "random", prediction_type = "IID",
    varcomp = vc, residual_structure = "IID")$mean_PEV
  expect_lt(run(4L), run(2L))
})
