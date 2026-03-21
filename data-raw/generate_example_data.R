# ============================================================
# OptiSparseMET_example_data generator
#
# Produces a single named list saved to data/ as an .RData file.
# All components share the same 120 treatment IDs and 4 environments,
# so any combination of components can be passed to package functions
# without modification.
#
# Run this script once from the package root to regenerate the data object.
# ============================================================

set.seed(123)

# ------------------------------------------------------------
# 1. Core treatment and environment structure
# ------------------------------------------------------------

n_treatments <- 120L
treatments   <- paste0("L", sprintf("%03d", seq_len(n_treatments)))
environments <- c("Env_A", "Env_B", "Env_C", "Env_D")

# 16 families assigned randomly so that family group sizes are unequal —
# this produces a realistic imbalance for demonstrating group-guided
# allocation and adjacency control.
treatment_info <- data.frame(
  Treatment = treatments,
  Family    = paste0("F", sprintf("%02d", sample(1:16, n_treatments, replace = TRUE))),
  stringsAsFactors = FALSE
)

# First 8 treatments designated as common — forced into all environments
# before sparse allocation. Chosen from the head of the treatment vector
# for reproducibility; in practice these would be elite connectors or
# benchmark lines.
common_treatments <- treatments[1:8]

# ------------------------------------------------------------
# 2. Seed availability
# ------------------------------------------------------------

# Seed quantities drawn from a range wide enough to produce a mix of
# replicated, unreplicated, and excluded treatments when passed to
# assign_replication_by_seed() at typical per-plot requirements.
# Even spacing (by = 2) avoids ties in seed-priority ranking.
seed_info <- data.frame(
  Treatment     = treatments,
  SeedAvailable = sample(seq(20, 140, by = 2), n_treatments, replace = TRUE),
  stringsAsFactors = FALSE
)

# Per-plot seed requirements differ across environments to reflect that
# plot size and sowing density vary by location. Env_C has the lowest
# requirement (6), which means more treatments can be replicated there.
seed_required_per_plot <- data.frame(
  Environment        = environments,
  SeedRequiredPerPlot = c(8, 10, 6, 9),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 3. Relationship matrices
# ------------------------------------------------------------

# Helper: generate a positive semi-definite matrix from random feature scores.
# tcrossprod(Xs) / n_features approximates a genomic or pedigree relationship
# matrix. The diagonal is inflated by diag_scale to ensure strict positive
# definiteness — required for Cholesky factorization in efficiency evaluation
# and for well-conditioned eigendecomposition in clustering.
make_psd_matrix <- function(ids, n_features = 20, diag_scale = 1) {
  X  <- matrix(rnorm(length(ids) * n_features), nrow = length(ids), ncol = n_features)
  Xs <- scale(X, center = TRUE, scale = TRUE)
  K  <- tcrossprod(Xs) / n_features
  diag(K) <- diag(K) + diag_scale
  rownames(K) <- ids
  colnames(K) <- ids
  K
}

# GRM: more features -> finer genomic differentiation between lines.
# A:   fewer features and larger diagonal -> broader pedigree clusters.
# K:   intermediate; used for prediction efficiency and dispersion scoring.
# All three use the same treatment IDs as row/column names, so no id_map
# is needed when passing them directly to package functions.
OptiSparseMET_GRM <- make_psd_matrix(treatments, n_features = 25, diag_scale = 1.0)
OptiSparseMET_A   <- make_psd_matrix(treatments, n_features = 15, diag_scale = 1.2)
OptiSparseMET_K   <- make_psd_matrix(treatments, n_features = 20, diag_scale = 1.1)

# ------------------------------------------------------------
# 4. Pre-built argument lists for allocate_sparse_met()
# ------------------------------------------------------------

# Random balanced: environment capacities are intentionally unequal
# (45, 50, 40, 48) to demonstrate that random_balanced handles
# heterogeneous sizes while balanced_incomplete requires equal capacity.
# target_replications = 1 keeps the example tractable for illustration.
sparse_example_args_random_balanced <- list(
  treatments                     = treatments,
  environments                   = environments,
  allocation_method              = "random_balanced",
  n_test_entries_per_environment = c(45L, 50L, 40L, 48L),
  target_replications            = 1L,
  common_treatments              = common_treatments,
  allow_approximate              = TRUE,
  seed                           = 123
)

# Balanced incomplete: equal capacity per environment (38) chosen so that
# the slot identity J* x r = I x k* is exactly satisfied after removing
# 8 common treatments: 112 x 1 = 4 x 28. allow_approximate = FALSE
# verifies that exact feasibility holds for this configuration.
sparse_example_args_balanced_incomplete <- list(
  treatments                     = treatments,
  environments                   = environments,
  allocation_method              = "balanced_incomplete",
  n_test_entries_per_environment = rep(38L, length(environments)),
  target_replications            = 1L,
  common_treatments              = common_treatments,
  allow_approximate              = FALSE,
  seed                           = 123
)

# ------------------------------------------------------------
# 5. Environment-specific local design specifications
# ------------------------------------------------------------

# Each environment uses a different design engine and replication mode,
# making env_design_specs suitable for demonstrating the full range of
# plan_sparse_met_design() behaviour in a single call.

env_design_specs <- list(
  
  # Env_A: p-rep design via prep_famoptg().
  # Up to 12 treatments replicated twice (max_prep = 12); remaining
  # allocated treatments receive one plot. Treatments with insufficient
  # seed for 2 plots are downgraded to 1 plot rather than excluded.
  # 10 x 12 grid accommodates 2 checks x 4 blocks + 12 x 2 + remaining
  # unreplicated entries.
  Env_A = list(
    design               = "prep_famoptg",
    replication_mode     = "p_rep",
    desired_replications = 2L,
    max_prep             = 12L,
    shortage_action      = "downgrade",
    check_treatments     = c("CHK1", "CHK2"),
    check_families       = c("CHECK", "CHECK"),
    n_blocks             = 4L,
    n_rows               = 10L,
    n_cols               = 12L,
    order                = "row",
    serpentine           = TRUE,
    cluster_source       = "Family",
    eval_efficiency      = FALSE,
    use_dispersion       = FALSE
  ),
  
  # Env_B: alpha row-column stream design via alpha_rc_stream().
  # n_reps = 2 with a 10 x 12 grid; min_entry_slots_per_block = 8 ensures
  # incomplete blocks are not too small after inserting 2 checks per block.
  # entry_treatments and entry_families are set internally by
  # plan_sparse_met_design() and must not appear here.
  Env_B = list(
    design                    = "alpha_rc_stream",
    check_treatments          = c("CHK1", "CHK2"),
    check_families            = c("CHECK", "CHECK"),
    n_reps                    = 2L,
    n_rows                    = 10L,
    n_cols                    = 12L,
    order                     = "row",
    serpentine                = TRUE,
    min_entry_slots_per_block = 8L,
    cluster_source            = "Family",
    eval_efficiency           = FALSE,
    use_dispersion            = FALSE,
    verbose                   = FALSE
  ),
  
  # Env_C: augmented repeated-check design via prep_famoptg().
  # All allocated entries are unreplicated (replication_mode = "augmented").
  # Column-major traversal without serpentine for contrast with Env_A.
  # A 10 x 10 grid is appropriate for a large unreplicated screening set.
  Env_C = list(
    design           = "prep_famoptg",
    replication_mode = "augmented",
    check_treatments = c("CHK1", "CHK2"),
    check_families   = c("CHECK", "CHECK"),
    n_blocks         = 4L,
    n_rows           = 10L,
    n_cols           = 10L,
    order            = "column",
    serpentine       = FALSE,
    cluster_source   = "Family",
    eval_efficiency  = FALSE,
    use_dispersion   = FALSE
  ),
  
  # Env_D: RCBD-type repeated-check design via prep_famoptg().
  # All allocated entries targeted for 2 replications (rcbd_type mode).
  # Treatments with insufficient seed are downgraded to 1 plot.
  # 4 blocks x (2 checks + entries/2) fills the 10 x 12 grid.
  Env_D = list(
    design               = "prep_famoptg",
    replication_mode     = "rcbd_type",
    desired_replications = 2L,
    shortage_action      = "downgrade",
    check_treatments     = c("CHK1", "CHK2"),
    check_families       = c("CHECK", "CHECK"),
    n_blocks             = 4L,
    n_rows               = 10L,
    n_cols               = 12L,
    order                = "row",
    serpentine           = FALSE,
    cluster_source       = "Family",
    eval_efficiency      = FALSE,
    use_dispersion       = FALSE
  )
)

# ------------------------------------------------------------
# 6. Assemble final example object
# ------------------------------------------------------------

OptiSparseMET_example_data <- list(
  treatments                          = treatments,
  environments                        = environments,
  treatment_info                      = treatment_info,
  common_treatments                   = common_treatments,
  seed_info                           = seed_info,
  seed_required_per_plot              = seed_required_per_plot,
  OptiSparseMET_GRM                   = OptiSparseMET_GRM,
  OptiSparseMET_A                     = OptiSparseMET_A,
  OptiSparseMET_K                     = OptiSparseMET_K,
  sparse_example_args_random_balanced     = sparse_example_args_random_balanced,
  sparse_example_args_balanced_incomplete = sparse_example_args_balanced_incomplete,
  env_design_specs                    = env_design_specs
)

# ------------------------------------------------------------
# 7. Save to data/
# ------------------------------------------------------------

# bzip2 compression and version = 2 ensure compatibility with R >= 3.5.
# Regenerate by sourcing this file from the package root directory.
save(
  OptiSparseMET_example_data,
  file     = "data/OptiSparseMET_example_data.RData",
  compress = "bzip2",
  version  = 2
)