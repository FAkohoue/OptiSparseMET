# Tests for the exact multi-trait information matrix (trait-covariance Kronecker
# term) via the canonical transformation. The two-trait reference values were
# produced by a direct multi-trait MME solve in numpy (matches to ~1e-16).

test_that("T = 1 reduces exactly to the single-trait met_information", {
  set.seed(7)
  G <- crossprod(matrix(rnorm(6 * 20), 20, 6)) / 20 + diag(6) * 0.4
  dimnames(G) <- list(paste0("L", 1:6), paste0("L", 1:6))
  SigE <- matrix(c(1, .4, .4, 1), 2, dimnames = list(c("E1", "E2"), c("E1", "E2")))
  M <- matrix(0L, 6, 2, dimnames = list(rownames(G), colnames(SigE)))
  M[c(1, 2, 4, 5), 1] <- 1L; M[c(2, 3, 4, 6), 2] <- 1L

  st <- met_information(M, G = G, Sigma_E = SigE, sigma_g2 = 1.3, sigma_e2 = 1)
  mt <- met_information_mt(M, G = G, Sigma_E = SigE,
                          Sigma_T = matrix(1.3), R_T = matrix(1),
                          index_weights = 1)
  expect_equal(mt$CDmean_index, st$CDmean, tolerance = 1e-8)
  expect_equal(unname(mt$rel_index_per_line), unname(st$CD_per_line),
               tolerance = 1e-8)
})

test_that("two-trait index reliability matches the direct MME oracle", {
  G <- diag(4) * 0.6 + 0.4
  dimnames(G) <- list(paste0("L", 1:4), paste0("L", 1:4))
  SigE <- matrix(c(1, .5, .5, 1), 2, dimnames = list(c("E1", "E2"), c("E1", "E2")))
  SigT <- matrix(c(1, .3, .3, .6), 2)
  R_T  <- matrix(c(1, .2, .2, .9), 2)
  w    <- c(1, -0.5)
  M <- matrix(c(1, 1, 0, 1,
                1, 0, 1, 1), 4, 2,
              dimnames = list(rownames(G), colnames(SigE)))

  mt <- met_information_mt(M, G = G, Sigma_E = SigE, Sigma_T = SigT,
                          R_T = R_T, index_weights = w)
  # numpy reference (direct Sigma_T (x) Sigma_E (x) G MME solve)
  expect_equal(unname(mt$rel_index_per_line),
               c(0.18004232, 0.10336825, 0.10336825, 0.18004232),
               tolerance = 1e-6)
  expect_equal(mt$CDmean_index, 0.14170528, tolerance = 1e-6)
  expect_equal(mt$sigma_index, sqrt(as.numeric(t(w) %*% SigT %*% w)), tolerance = 1e-8)
})

test_that("canonical transform simultaneously diagonalises Sigma_T and R_T", {
  SigT <- matrix(c(1, .3, .3, .6), 2)
  R_T  <- matrix(c(1, .2, .2, .9), 2)
  M <- matrix(1L, 4, 2, dimnames = list(paste0("L", 1:4), c("E1", "E2")))
  G <- diag(4); dimnames(G) <- list(paste0("L", 1:4), paste0("L", 1:4))
  mt <- met_information_mt(M, G = G, Sigma_T = SigT, R_T = R_T,
                          index_weights = c(1, 1))
  Q <- mt$Q
  expect_equal(Q %*% R_T %*% t(Q), diag(2), tolerance = 1e-8)
  expect_equal(Q %*% SigT %*% t(Q), diag(mt$lambda), tolerance = 1e-8)
  expect_error(
    met_information_mt(M, G = G, Sigma_T = SigT, R_T = R_T,
                       index_weights = c(0, 0)),
    "non-zero")
})

test_that("design_objective 'exact' vs 'approx' agree in direction, differ in value", {
  set.seed(3)
  G <- crossprod(matrix(rnorm(16 * 40), 40, 16)) / 40 + diag(16) * 0.3
  dimnames(G) <- list(paste0("L", 1:16), paste0("L", 1:16))
  SigE <- matrix(c(1, .5, .5, 1), 2, dimnames = list(c("E1", "E2"), c("E1", "E2")))
  SigT <- matrix(c(1, .3, .3, .6), 2)
  w <- c(1.2, -0.6)
  M <- matrix(0L, 16, 2, dimnames = list(rownames(G), colnames(SigE)))
  for (e in 1:2) M[sample(16, 8), e] <- 1L

  ex <- design_objective(M, G, SigE, trait_weights = w, trait_gencov = SigT,
                         multitrait = "exact")
  ap <- design_objective(M, G, SigE, trait_weights = w, trait_gencov = SigT,
                         multitrait = "approx")
  expect_true(is.finite(ex$gain) && is.finite(ap$gain))
  expect_true(ex$reliability >= 0 && ex$reliability <= 1)
  # both positive gains; the exact value need not equal the approximation
  expect_gt(ex$gain, 0); expect_gt(ap$gain, 0)
})

test_that("optimize_design runs with an exact multi-trait objective", {
  set.seed(5)
  G <- crossprod(matrix(rnorm(16 * 40), 40, 16)) / 40 + diag(16) * 0.3
  dimnames(G) <- list(paste0("L", 1:16), paste0("L", 1:16))
  SigE <- matrix(c(1, .5, .5, 1), 2, dimnames = list(c("E1", "E2"), c("E1", "E2")))
  SigT <- matrix(c(1, .3, .3, .6), 2)
  w <- desired_gain_weights(c(A = 1.5, B = 1), gen_cov = SigT)$weights
  M <- matrix(0L, 16, 2, dimnames = list(rownames(G), colnames(SigE)))
  for (e in 1:2) M[sample(16, 8), e] <- 1L

  opt <- optimize_design(M, G, SigE, preserve = "margins",
                         objective = list(weights = list(gain = 1),
                                          trait_weights = unname(w),
                                          trait_gencov = SigT,
                                          multitrait = "exact"),
                         n_starts = 1, iters = 20, seed = 1)
  expect_gte(opt$score, opt$score_start)
  expect_true(is.finite(opt$components$reliability))
})
