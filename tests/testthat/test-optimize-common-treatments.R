make_common_optimizer_data <- function() {
  ids <- paste0("G", seq_len(12))
  envs <- c("E1", "E2", "E3")
  G <- diag(12)
  G[G == 0] <- 0.05
  diag(G) <- 1
  dimnames(G) <- list(ids, ids)
  Sigma_E <- matrix(c(
    1.0, 0.45, 0.20,
    0.45, 1.0, 0.30,
    0.20, 0.30, 1.0
  ), 3, 3, dimnames = list(envs, envs))
  treatment_info <- data.frame(
    Treatment = ids,
    Family = rep(c("F1", "F2", "F3"), each = 4),
    stringsAsFactors = FALSE
  )
  seed_info <- data.frame(
    Treatment = ids,
    SeedAvailable = rep(100, 12),
    stringsAsFactors = FALSE
  )
  list(
    ids = ids, envs = envs, G = G, Sigma_E = Sigma_E,
    treatment_info = treatment_info, seed_info = seed_info,
    q = stats::setNames(c(2, 3, 4), envs),
    capacity = stats::setNames(c(8L, 8L, 8L), envs)
  )
}

test_that("common-set optimisation returns binary presence only", {
  x <- make_common_optimizer_data()
  out <- optimize_common_treatments(
    G = x$G, Sigma_E = x$Sigma_E,
    treatment_info = x$treatment_info, seed_info = x$seed_info,
    seed_required_per_plot = x$q,
    entry_capacities = x$capacity,
    n_common = 3L, n_groups = 1L,
    seed = 4
  )

  expect_equal(dim(out$common_presence), c(3L, 3L))
  expect_true(all(out$common_presence == 1L))
  expect_equal(nrow(out$pairwise_connectivity), 3L)
  expect_equal(
    out$pairwise_connectivity$AchievedDistinctSharedTreatments,
    rep(3, 3L)
  )
  expect_true(all(out$pairwise_connectivity$TargetFraction >= 0 &
                    out$pairwise_connectivity$TargetFraction <= 1))
})

test_that("common seed is charged once at each environment where present", {
  x <- make_common_optimizer_data()
  out <- optimize_common_treatments(
    G = x$G, Sigma_E = x$Sigma_E,
    treatment_info = x$treatment_info, seed_info = x$seed_info,
    seed_required_per_plot = x$q,
    entry_capacities = x$capacity,
    n_common = 3L, n_groups = 1L,
    minimum_seed_buffer = 0,
    seed = 4
  )

  expected <- as.numeric(out$common_presence %*% x$q)
  ledger <- out$seed_ledger[
    match(rownames(out$common_presence), out$seed_ledger$Treatment), ]
  expect_equal(ledger$SeedAllocatedToCommonPresence, expected)
  expect_true(all(out$seed_ledger$Feasible))
  expect_true(all(
    out$seed_ledger$SeedAllocatedToCommonPresence[
      !out$seed_ledger$IsCommon
    ] == 0
  ))
})

test_that("automatic count search is capacity-feasible and Pareto-auditable", {
  x <- make_common_optimizer_data()
  out <- optimize_common_treatments(
    G = x$G, Sigma_E = x$Sigma_E,
    treatment_info = x$treatment_info, seed_info = x$seed_info,
    seed_required_per_plot = x$q,
    entry_capacities = x$capacity,
    n_groups = 1L,
    seed = 9
  )

  expect_true(out$n_common %in% out$comparison$n_common)
  expect_true(all(out$comparison$n_common <= 6L))
  expect_true(all(out$comparison$testing_breadth >= 0 &
                    out$comparison$testing_breadth <= 1))
  expect_true(any(out$comparison$pareto))
  expect_equal(out$rationale$pair_aggregate, "maximin")
  expect_equal(
    out$comparison$pair_connectivity,
    out$comparison$minimum_pair_connectivity
  )
  expect_equal(out$pareto,
               out$comparison[out$comparison$pareto, , drop = FALSE])
  expect_equal(sum(out$selection_diagnostics$Selected), out$n_common)
  expect_equal(out$selection_diagnostics$NestedRank,
               seq_len(nrow(out$selection_diagnostics)))
})

test_that("breeders can choose maximin, CVaR, or mean pair protection", {
  x <- make_common_optimizer_data()
  fit <- function(pair_aggregate) optimize_common_treatments(
    G = x$G, Sigma_E = x$Sigma_E,
    treatment_info = x$treatment_info, seed_info = x$seed_info,
    seed_required_per_plot = x$q,
    entry_capacities = x$capacity,
    n_common = 3L, n_groups = 1L,
    pair_aggregate = pair_aggregate,
    pair_cvar_alpha = 0.5,
    seed = 7
  )
  out <- lapply(c("maximin", "cvar", "mean"), fit)
  names(out) <- c("maximin", "cvar", "mean")

  expect_equal(out$maximin$rationale$pair_aggregate, "maximin")
  expect_equal(out$cvar$rationale$pair_aggregate, "cvar")
  expect_equal(out$mean$rationale$pair_aggregate, "mean")
  expect_equal(
    out$maximin$comparison$pair_connectivity,
    min(out$maximin$pairwise_connectivity$TargetFraction)
  )
  expect_equal(
    out$mean$comparison$pair_connectivity,
    mean(out$mean$pairwise_connectivity$TargetFraction)
  )
  expect_gte(
    out$cvar$comparison$pair_connectivity,
    min(out$cvar$pairwise_connectivity$TargetFraction)
  )
  expect_lte(
    out$cvar$comparison$pair_connectivity,
    mean(out$cvar$pairwise_connectivity$TargetFraction)
  )
})

test_that("family coverage failure is explicit rather than silently relaxed", {
  x <- make_common_optimizer_data()
  x$seed_info$SeedAvailable[x$treatment_info$Family == "F3"] <- 0
  expect_error(
    optimize_common_treatments(
      G = x$G, Sigma_E = x$Sigma_E,
      treatment_info = x$treatment_info, seed_info = x$seed_info,
      seed_required_per_plot = x$q,
      entry_capacities = x$capacity,
      min_per_family = 1L
    ),
    "cannot provide.*F3"
  )
})

test_that("p-rep environment specs enforce user-fixed treatment reps", {
  ids <- paste0("L", 1:6)
  envs <- c("E1", "E2")
  ti <- data.frame(
    Treatment = ids,
    Family = rep(c("F1", "F2"), each = 3),
    stringsAsFactors = FALSE
  )
  si <- data.frame(
    Treatment = ids,
    SeedAvailable = rep(20, 6),
    stringsAsFactors = FALSE
  )
  base_spec <- list(
    design = "met_prep_famoptg",
    replication_mode = "p_rep",
    desired_replications = 2L,
    max_prep = 1L,
    max_extra_replication_plots = 1L,
    shortage_action = "downgrade",
    check_treatments = "CHK1",
    check_families = "CHECK",
    n_blocks = 3L,
    n_rows = 3L,
    n_cols = 3L,
    order = "row",
    serpentine = FALSE,
    cluster_source = "Family",
    use_dispersion = FALSE
  )
  specs <- list(E1 = base_spec, E2 = base_spec)
  specs$E1$fixed_treatment_reps <- c(L1 = 1L, L2 = 2L)
  specs$E2$fixed_treatment_reps <- c(L1 = 2L, L2 = 1L)

  out <- plan_sparse_met_design(
    treatments = ids, environments = envs,
    allocation_method = "random_balanced",
    n_test_entries_per_environment = c(E1 = 4L, E2 = 4L),
    common_treatments = c("L1", "L2"),
    env_design_specs = specs,
    treatment_info = ti, seed_info = si,
    seed_required_per_plot = c(E1 = 1, E2 = 1),
    minimum_seed_buffer = 0,
    seed = 10
  )

  tab1 <- table(out$environment_designs$E1$field_book$Treatment)
  tab2 <- table(out$environment_designs$E2$field_book$Treatment)
  expect_equal(as.integer(tab1[c("L1", "L2")]), c(1, 2))
  expect_equal(as.integer(tab2[c("L1", "L2")]), c(2, 1))
})

test_that("three reps consume two physical extra-plot slots", {
  ids <- paste0("L", 1:4)
  ti <- data.frame(Treatment = ids, Family = "F1",
                   stringsAsFactors = FALSE)
  si <- data.frame(Treatment = ids, SeedAvailable = 20,
                   stringsAsFactors = FALSE)
  spec <- list(
    design = "met_prep_famoptg",
    replication_mode = "p_rep",
    desired_replications = 2L,
    max_prep = 2L,
    max_extra_replication_plots = 1L,
    fixed_treatment_reps = c(L1 = 3L),
    shortage_action = "downgrade",
    check_treatments = "CHK1",
    check_families = "CHECK",
    n_blocks = 3L, n_rows = 3L, n_cols = 3L,
    cluster_source = "Family", use_dispersion = FALSE
  )
  expect_error(
    plan_sparse_met_design(
      treatments = ids, environments = c("E1", "E2"),
      allocation_method = "random_balanced",
      n_test_entries_per_environment = c(E1 = 3L, E2 = 3L),
      common_treatments = "L1",
      env_design_specs = list(E1 = spec, E2 = spec),
      treatment_info = ti, seed_info = si,
      seed_required_per_plot = c(E1 = 1, E2 = 1),
      seed = 2
    ),
    "requires 2 extra plots"
  )
})
