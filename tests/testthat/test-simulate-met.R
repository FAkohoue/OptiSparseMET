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
