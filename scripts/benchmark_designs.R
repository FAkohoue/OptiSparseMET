# Reproducible public benchmark for OptiSparseMET
#
# Run manually after installing the package:
#   Rscript scripts/benchmark_designs.R
#
# This script writes no files. It uses synthetic data with known truth and does
# not represent evidence for a particular crop or breeding programme.

library(OptiSparseMET)

benchmark_seed <- 2026L
set.seed(benchmark_seed)

# 1. Known-truth environmental covariance recovery ----------------------------

environments <- paste0("E", 1:6)
K_weather <- outer(
  1:6, 1:6, function(i, j) exp(-abs(i - j) / 1.5)
)
K_soil <- outer(
  c(1, 3, 2, 6, 4, 5),
  c(1, 3, 2, 6, 4, 5),
  function(i, j) exp(-abs(i - j) / 1.5)
)
dimnames(K_weather) <- dimnames(K_soil) <-
  list(environments, environments)

environment_benchmark <- simulate_environment_model_benchmark(
  kernels = list(weather = K_weather, soil = K_soil),
  true_weights = c(weather = 0.65, soil = 0.25, identity = 0.10),
  n_genotypes = 150,
  n_years = 8,
  reps_per_cell = 2,
  missing_rate = 0.05,
  n_boot = 200,
  seed = benchmark_seed
)

message("\nKnown and estimated environmental component weights")
print(environment_benchmark$weight_recovery)
message("\nEnvironmental covariance recovery")
print(environment_benchmark$benchmark$scores)
message("\nBlocked environmental validation")
print(environment_benchmark$benchmark$cross_validation)

# 2. Equal-budget allocation comparators --------------------------------------

treatments <- paste0("L", 1:76)
common <- treatments[1:4]
trial_environments <- paste0("S", 1:4)

# With four common treatments, 72 sparse treatments, and two environments per
# sparse treatment: 72 * 2 = 144 sparse slots = 4 * (40 - 4).
reference_designs <- make_met_benchmark_designs(
  treatments = treatments,
  environments = trial_environments,
  n_test_entries_per_environment = 40,
  target_replications = 2,
  common_treatments = common,
  m4_fallback = "omit",
  seed = benchmark_seed
)

message("\nReference-design construction")
print(reference_designs$diagnostics)
print(reference_designs$m4_feasibility)

# 3. Construct the OptiSparseMET candidate ------------------------------------

marker_surrogate <- matrix(rnorm(length(treatments) * 100),
                           length(treatments), 100)
G <- tcrossprod(marker_surrogate) / 100 + diag(length(treatments)) * 0.20
dimnames(G) <- list(treatments, treatments)

Sigma_central <- matrix(0.35, 4, 4)
diag(Sigma_central) <- 1
dimnames(Sigma_central) <- list(trial_environments, trial_environments)

Sigma_weak <- matrix(0.10, 4, 4)
diag(Sigma_weak) <- 1
dimnames(Sigma_weak) <- list(trial_environments, trial_environments)

# Twenty candidates can fund four sites and are eligible to become common.
# The others can fund three sites: enough for sparse testing but not complete
# testing. This creates a real network-wide seed constraint.
seed_inventory <- setNames(
  ifelse(seq_along(treatments) <= 20, 5, 3),
  treatments
)
treatment_info <- data.frame(
  Treatment = treatments,
  Family = rep(paste0("F", 1:4), length.out = length(treatments)),
  stringsAsFactors = FALSE
)
seed_info <- data.frame(
  Treatment = treatments,
  SeedAvailable = unname(seed_inventory),
  stringsAsFactors = FALSE
)

optimised_common <- optimize_common_treatments(
  G = G,
  Sigma_E = Sigma_central,
  treatment_info = treatment_info,
  seed_info = seed_info,
  seed_required_per_plot = setNames(rep(1, 4), trial_environments),
  entry_capacities = setNames(rep(40L, 4), trial_environments),
  minimum_seed_buffer = 0,
  n_groups = 4,
  pair_aggregate = "maximin",
  seed = benchmark_seed
)

optimised_allocation <- allocate_sparse_met(
  treatments = treatments,
  environments = trial_environments,
  allocation_method = "random_balanced",
  n_test_entries_per_environment = 40,
  target_replications = 2,
  common_treatments = optimised_common$selected,
  seed_available = seed_inventory,
  seed_required_per_environment =
    setNames(rep(1, 4), trial_environments),
  minimum_seed_buffer = 0,
  balance = "both",
  balance_iter = 5000,
  Sigma_E = Sigma_central,
  pair_aggregate = "maximin",
  allocation_group_source = "Family",
  treatment_info = treatment_info,
  min_env_per_group = 2,
  seed = benchmark_seed
)
optimised_reps <- optimised_allocation$allocation_matrix
reference_designs$designs$OptiSparseMET_optimised <- list(
  allocation_matrix = optimised_allocation$allocation_matrix,
  reps = optimised_reps
)

# 4. Paired quantitative-genetic and operational benchmark --------------------

design_benchmark <- benchmark_met_designs(
  designs = reference_designs,
  G = G,
  Sigma_E = Sigma_central,
  Sigma_E_candidates = list(
    central = Sigma_central,
    weakly_correlated = Sigma_weak
  ),
  sigma_g2 = 1,
  sigma_e2 = c(0.5, 1, 2),
  environment_dropout = list(
    none = character(),
    fourth_site_lost = "S4"
  ),
  cost_per_plot = c(S1 = 1.0, S2 = 1.0, S3 = 1.2, S4 = 1.3),
  environment_plot_capacity = c(S1 = 40, S2 = 40, S3 = 40, S4 = 40),
  seed_available = seed_inventory,
  seed_required_per_plot = 1,
  minimum_seed_buffer = 0,
  n_sim = 200,
  select_fraction = 0.10,
  seed = benchmark_seed
)

message("\nCross-scenario design summary")
print(design_benchmark$design_summary)
message("\nOptimised common-set decision")
print(optimised_common$comparison)
print(optimised_common$selected)
message("\nOperational feasibility")
print(design_benchmark$feasibility)
message("\nScenario-specific results")
print(design_benchmark$scenario_results)

# Uncomment in an interactive graphics session.
# plot_pareto_frontier(
#   design_benchmark, metric = "accuracy", label_designs = TRUE
# )
