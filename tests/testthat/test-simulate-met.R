# Tests for simulate_met() (P3.3): realized accuracy / genetic gain under a
# G x E model. Numpy-verified direction: a fully balanced design yields higher
# across-TPE accuracy than a sparse one.

mk_sim <- function() {
  set.seed(11)
  J <- 40L; E <- 4L
  B <- matrix(rnorm(J * J), J, J); G <- tcrossprod(B) / J + diag(J) * 0.3
  dimnames(G) <- list(paste0("L", 1:J), paste0("L", 1:J))
  SigE <- matrix(c(1, .6, .2, -.1,
                   .6, 1, .3,  0,
                   .2, .3, 1, .5,
                  -.1,  0, .5, 1), 4, 4)
  list(G = G, SigE = SigE, J = J, E = E)
}

test_that("simulate_met survives a singular Sigma_E (redundant environments)", {
  set.seed(1)
  J <- 20L
  G <- crossprod(matrix(rnorm(J * 40), 40, J)) / 40 + diag(J) * 0.3
  dimnames(G) <- list(paste0("L", 1:J), paste0("L", 1:J))
  # E2 is identical to E1 -> Sigma_E is singular (rank 2), Cholesky fails
  SigE <- matrix(c(1, 1, .2,
                   1, 1, .2,
                   .2, .2, 1), 3, dimnames = list(paste0("E", 1:3), paste0("E", 1:3)))
  M <- matrix(0L, J, 3, dimnames = list(rownames(G), colnames(SigE)))
  for (e in 1:3) M[sample.int(J, 10), e] <- 1L
  sm <- simulate_met(M, G = G, Sigma_E = SigE, n_sim = 5, seed = 1)
  expect_true(is.finite(sm$accuracy_mean))       # eigen-root fallback, no crash
})

test_that("accuracy is a valid correlation and gain is positive on average", {
  d <- mk_sim()
  full <- matrix(1L, d$J, d$E,
                 dimnames = list(rownames(d$G), paste0("E", seq_len(d$E))))
  res <- simulate_met(full, G = d$G, Sigma_E = d$SigE,
                      n_sim = 30, seed = 1)
  expect_true(res$accuracy_mean > -1 && res$accuracy_mean < 1)
  expect_gt(res$gain_mean, 0)
})

test_that("balanced design beats a sparse design on realized accuracy", {
  d <- mk_sim()
  ln <- rownames(d$G); en <- paste0("E", seq_len(d$E))
  full <- matrix(1L, d$J, d$E, dimnames = list(ln, en))
  set.seed(3)
  sparse <- matrix(0L, d$J, d$E, dimnames = list(ln, en))
  for (j in seq_len(d$J)) sparse[j, sample(d$E, 2)] <- 1L   # each line in 2 of 4 envs

  acc_full   <- simulate_met(full,   G = d$G, Sigma_E = d$SigE, n_sim = 40, seed = 5)$accuracy_mean
  acc_sparse <- simulate_met(sparse, G = d$G, Sigma_E = d$SigE, n_sim = 40, seed = 5)$accuracy_mean
  expect_gt(acc_full, acc_sparse)
})

test_that("residual information depends on reps divided by residual variance", {
  d <- mk_sim()
  M <- matrix(1L, d$J, d$E,
              dimnames = list(rownames(d$G), paste0("E", seq_len(d$E))))
  r1 <- simulate_met(M, d$G, d$SigE, sigma_e2 = 1,
                     reps = M, n_sim = 12, seed = 19)
  r2 <- simulate_met(M, d$G, d$SigE, sigma_e2 = 2,
                     reps = 2 * M, n_sim = 12, seed = 19)
  expect_equal(r1$accuracy_mean, r2$accuracy_mean, tolerance = 1e-12)
  expect_equal(r1$gain_mean, r2$gain_mean, tolerance = 1e-12)
  expect_equal(r1$mean_PEV, r2$mean_PEV, tolerance = 1e-12)
})

test_that("simulation reports finite Monte Carlo uncertainty intervals", {
  d <- mk_sim()
  M <- matrix(1L, d$J, d$E,
              dimnames = list(rownames(d$G), paste0("E", seq_len(d$E))))
  out <- simulate_met(M, d$G, d$SigE, n_sim = 10, seed = 2)
  expect_true(all(is.finite(c(out$accuracy_se, out$gain_se,
                              out$accuracy_ci95, out$gain_ci95))))
  expect_lte(out$accuracy_ci95[1], out$accuracy_mean)
  expect_gte(out$accuracy_ci95[2], out$accuracy_mean)
  expect_lte(out$gain_ci95[1], out$gain_mean)
  expect_gte(out$gain_ci95[2], out$gain_mean)
  expect_error(simulate_met(M, d$G, d$SigE, n_sim = 1), "integer >= 2")
})
