# Tests for optimize_allocation_gxe() (P4.3) and recommend_replication() (P4.4).

mk <- function(J = 24L, E = 4L) {
  set.seed(9)
  B <- matrix(rnorm(J * J), J, J); G <- tcrossprod(B) / J + diag(J) * 0.4
  dimnames(G) <- list(paste0("L", seq_len(J)), paste0("L", seq_len(J)))
  SigE <- matrix(c(1, .6, .2, -.1,
                   .6, 1, .3,  0,
                   .2, .3, 1, .5,
                  -.1,  0, .5, 1), 4, 4)
  ln <- rownames(G); en <- paste0("E", seq_len(E))
  # equireplicate sparse start: each line in 2 of 4 environments
  set.seed(3)
  M <- matrix(0L, J, E, dimnames = list(ln, en))
  for (j in seq_len(J)) M[j, sample(E, 2)] <- 1L
  list(G = G, SigE = SigE, M = M)
}

test_that("G x E allocation optimisation lowers mean PEV and preserves margins", {
  d <- mk()
  res <- optimize_allocation_gxe(d$M, G = d$G, Sigma_E = d$SigE,
                                 iter = 60, seed = 1, max_dim = 500)
  expect_lte(res$mean_PEV_after, res$mean_PEV_before)          # never worsens
  expect_equal(rowSums(res$allocation_matrix), rowSums(d$M))   # replication kept
  expect_equal(colSums(res$allocation_matrix), colSums(d$M))   # env sizes kept
})

test_that("recommend_replication sweeps levels and accuracy rises with replication", {
  d <- mk()
  rr <- recommend_replication(d$M, G = d$G, Sigma_E = d$SigE,
                              replication_levels = c(1, 2),
                              n_sim = 20, seed = 2, max_dim = 500)
  expect_true(all(c("replication", "accuracy_mean", "mean_PEV") %in% names(rr$table)))
  acc <- rr$table$accuracy_mean[order(rr$table$replication)]
  expect_gte(acc[2], acc[1])                                    # more reps -> >= accuracy
  expect_true(rr$recommended %in% rr$table$replication)
})

test_that("seed budget flags infeasible replication", {
  d <- mk()
  rr <- recommend_replication(d$M, G = d$G, Sigma_E = d$SigE,
                              replication_levels = c(1, 5),
                              n_sim = 10, seed = 2, max_dim = 500,
                              seed_available = 3, seed_required_per_plot = 1)
  # each line is in 2 envs; at r = 5 that is 10 plots/line * 1 seed > 3 -> infeasible
  expect_false(rr$table$feasible[rr$table$replication == 5])
})
