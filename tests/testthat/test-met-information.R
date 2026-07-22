# Tests for met_information() (P3.2): the across-environment coupled
# information matrix. Oracle values from an independent numpy MME solve.

mk <- function() {
  G <- matrix(c(1, .2, .2, .2, 1, .2, .2, .2, 1), 3, 3,
              dimnames = list(paste0("L", 1:3), paste0("L", 1:3)))
  SigE <- matrix(c(1, .5, .5, 1), 2, 2)
  alloc <- matrix(c(1, 1, 0,
                    1, 0, 1), nrow = 3,
                  dimnames = list(paste0("L", 1:3), paste0("E", 1:2)))
  list(G = G, SigE = SigE, alloc = alloc)
}

test_that("across-TPE mean_PEV and CDmean match the numpy MME oracle", {
  d <- mk()
  res <- met_information(d$alloc, G = d$G, Sigma_E = d$SigE,
                         sigma_g2 = 1, sigma_e2 = 1)
  expect_equal(res$mean_PEV, 0.6225, tolerance = 1e-4)
  expect_equal(res$CDmean,   0.17,   tolerance = 1e-4)
})

test_that("more connectivity (balanced) lowers across-TPE PEV", {
  d <- mk()
  full <- matrix(1L, 3, 2, dimnames = dimnames(d$alloc))
  sparse_pev <- met_information(d$alloc, G = d$G, Sigma_E = d$SigE)$mean_PEV
  full_pev   <- met_information(full,    G = d$G, Sigma_E = d$SigE)$mean_PEV
  expect_lt(full_pev, sparse_pev)
})

test_that("lower within-environment efficiency raises PEV (two-level coupling)", {
  d <- mk()
  full <- matrix(1L, 3, 2, dimnames = dimnames(d$alloc))
  hi <- met_information(full, G = d$G, Sigma_E = d$SigE,
                        env_efficiency = c(1, 1))$mean_PEV
  lo <- met_information(full, G = d$G, Sigma_E = d$SigE,
                        env_efficiency = c(0.5, 0.5))$mean_PEV
  expect_gt(lo, hi)
})
