# Tests for P1: equireplicate allocation, near-balanced concurrence / env-pair
# balance, and the strict-BIBD feasibility conditions.

test_that("equireplicate is the canonical name and M4 is a silent alias", {
  trts <- paste0("L", 1:60); envs <- paste0("E", 1:4)   # J*=60, I=4, r=2, k=30
  a_new <- allocate_sparse_met(trts, envs, "equireplicate",
                               n_test_entries_per_environment = 30,
                               target_replications = 2, seed = 1)
  expect_equal(a_new$summary$allocation_method, "equireplicate")

  # "M4" is accepted as an alias and resolves to equireplicate without any
  # deprecation message.
  msgs <- capture_messages(
    a_m4 <- allocate_sparse_met(trts, envs, "M4",
                                n_test_entries_per_environment = 30,
                                target_replications = 2, seed = 1)
  )
  expect_equal(a_m4$summary$allocation_method, "equireplicate")
  expect_false(any(grepl("deprecated", msgs)))

  # The former "balanced_incomplete" name has been removed entirely.
  expect_error(
    allocate_sparse_met(trts, envs, "balanced_incomplete",
                        n_test_entries_per_environment = 30,
                        target_replications = 2, seed = 1)
  )
})

test_that("balance preserves equal replication and equal environment size", {
  trts <- paste0("L", 1:60); envs <- paste0("E", 1:4)
  a <- allocate_sparse_met(trts, envs, "equireplicate",
                           n_test_entries_per_environment = 30,
                           target_replications = 2,
                           balance = "env_pair", balance_iter = 3000,
                           balance_seed = 7, seed = 1)
  M <- a$allocation_matrix
  expect_true(all(rowSums(M) == 2))          # equal replication
  expect_true(all(colSums(M) == 30))         # equal environment size
})

test_that("env_pair balancing reduces (or holds) env-pair variance", {
  trts <- paste0("L", 1:100); envs <- paste0("E", 1:5)  # J*=100, I=5, r=2, k=40
  set.seed(1)
  a0 <- allocate_sparse_met(trts, envs, "equireplicate",
                            n_test_entries_per_environment = 40,
                            target_replications = 2, balance = "none", seed = 1)
  a1 <- allocate_sparse_met(trts, envs, "equireplicate",
                            n_test_entries_per_environment = 40,
                            target_replications = 2, balance = "env_pair",
                            balance_iter = 4000, balance_seed = 3, seed = 1)
  ep_var <- function(M) { X <- crossprod(M); stats::var(X[upper.tri(X)]) }
  expect_lte(ep_var(a1$allocation_matrix), ep_var(a0$allocation_matrix) + 1e-9)
  expect_false(is.null(a1$balance_report))
})

test_that("feasibility reports strict-BIBD conditions correctly", {
  # J*=110, I=4, r=2: slot identity holds but Fisher (I >= J*) fails.
  f1 <- check_equireplicate_feasibility(
    n_treatments_total = 110, n_environments = 4,
    n_test_entries_per_environment = 55, target_replications = 2)
  expect_true(f1$feasible)
  expect_false(f1$strict_bibd_possible)
  expect_false(f1$fisher_ok)

  # Fano-type symmetric BIBD: v=b=7, r=k=3, lambda=1 -> strict BIBD possible.
  f2 <- check_equireplicate_feasibility(
    n_treatments_total = 7, n_environments = 7,
    n_test_entries_per_environment = 3, target_replications = 3)
  expect_true(f2$feasible)
  expect_true(f2$fisher_ok)
  expect_equal(f2$bibd_lambda, 1)
  expect_true(f2$strict_bibd_possible)
})
