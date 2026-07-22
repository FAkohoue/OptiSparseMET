# Tests for select_individuals() (P4.2 / decision 2). CDmean behaviour verified
# independently in numpy: CD in [0,1], rises with training size, exchange
# improves over random.

mkG <- function(p = 30) {
  set.seed(7)
  B <- matrix(rnorm(p * p), p, p)
  G <- tcrossprod(B) / p + diag(p) * 0.4
  dimnames(G) <- list(paste0("G", seq_len(p)), paste0("G", seq_len(p)))
  G
}

test_that("CDmean exchange selects n_train genotypes and beats random", {
  G <- mkG()
  ex  <- select_individuals(G, n_train = 10, method = "cdmean_exchange",
                            iters = 500, seed = 5)
  rnd <- select_individuals(G, n_train = 10, method = "random", seed = 5)
  expect_equal(length(ex$selected), 10L)
  expect_true(all(ex$selected %in% rownames(G)))
  expect_gte(ex$mean_CD, rnd$mean_CD)          # optimisation does not hurt
  expect_true(ex$mean_CD >= 0 && ex$mean_CD <= 1)
})

test_that("mean CDmean increases with training set size", {
  G <- mkG()
  cd_small <- select_individuals(G, n_train = 5,  method = "random", seed = 1)$mean_CD
  cd_large <- select_individuals(G, n_train = 20, method = "random", seed = 1)$mean_CD
  expect_gt(cd_large, cd_small)
})
