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

test_that("strict equireplication realizes unequal site loads with groups", {
  trts <- paste0("L", 1:12)
  envs <- paste0("E", 1:4)
  loads <- stats::setNames(c(8L, 6L, 4L, 6L), envs)
  ti <- data.frame(
    Treatment = trts,
    Family = rep(paste0("F", 1:3), each = 4),
    stringsAsFactors = FALSE)
  a <- allocate_sparse_met(
    treatments = trts, environments = envs,
    allocation_method = "equireplicate",
    n_test_entries_per_environment = loads,
    target_replications = 2L,
    allocation_group_source = "Family",
    treatment_info = ti,
    allow_approximate = FALSE,
    seed = 9)
  expect_true(all(rowSums(a$allocation_matrix) == 2L))
  expect_equal(colSums(a$allocation_matrix), loads)
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

test_that("adaptive connectivity uses distinct shared treatments, not plots", {
  trts <- paste0("L", 1:60)
  envs <- paste0("E", 1:4)
  latent <- rbind(
    E1 = c(1.0, 0.0),
    E2 = c(0.0, 1.0),
    E3 = c(0.9, 0.1),
    E4 = c(0.1, 0.9)
  )
  Sigma_E <- tcrossprod(latent) + diag(0.2, 4)
  dimnames(Sigma_E) <- list(envs, envs)

  out <- allocate_sparse_met(
    treatments = trts,
    environments = envs,
    allocation_method = "equireplicate",
    n_test_entries_per_environment = 30,
    target_replications = 2,
    Sigma_E = Sigma_E,
    pair_aggregate = "maximin",
    balance_iter = 5000,
    balance_seed = 8,
    seed = 3
  )

  expect_true(all(rowSums(out$allocation_matrix) == 2L))
  expect_true(all(colSums(out$allocation_matrix) == 30L))
  expect_equal(
    out$pairwise_connectivity$AchievedDistinctSharedTreatments,
    as.numeric(out$overlap_matrix[upper.tri(out$overlap_matrix)])
  )
  expect_gte(
    out$pair_refinement_report$after_score[["primary"]],
    out$pair_refinement_report$before_score[["primary"]]
  )
  rho <- stats::cov2cor(Sigma_E)[upper.tri(Sigma_E)]
  targets <- out$pairwise_connectivity$TargetDistinctSharedTreatments
  expect_gte(targets[which.min(abs(rho))],
             targets[which.max(abs(rho))])
})

test_that("balancing preserves group counts and treatment seed budgets", {
  trts <- paste0("L", 1:60)
  envs <- paste0("E", 1:4)
  ti <- data.frame(Treatment = trts,
                   Family = rep(paste0("F", 1:6), each = 10))
  available <- stats::setNames(rep(10, length(trts)), trts)
  costs <- stats::setNames(1:4, envs)
  args <- list(
    treatments = trts, environments = envs,
    allocation_method = "random_balanced",
    n_test_entries_per_environment = 30,
    target_replications = 2,
    allocation_group_source = "Family", treatment_info = ti,
    seed_available = available,
    seed_required_per_environment = costs,
    seed = 12)
  a0 <- do.call(allocate_sparse_met, c(args, list(balance = "none")))
  a1 <- do.call(allocate_sparse_met, c(
    args, list(balance = "env_pair", balance_iter = 2000,
               balance_seed = 13)))
  group_counts <- function(M)
    vapply(seq_len(ncol(M)), function(e)
      paste(rowsum(M[, e], ti$Family)[, 1], collapse = ","),
      character(1))
  expect_equal(group_counts(a1$allocation_matrix),
               group_counts(a0$allocation_matrix))
  used <- rowSums(sweep(a1$allocation_matrix, 2, costs, `*`))
  expect_true(all(used <= available + 1e-8))
  expect_true(isTRUE(a1$balance_report$constraints_preserved[["groups"]]))
  expect_true(isTRUE(a1$balance_report$constraints_preserved[["seed_budget"]]))
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
