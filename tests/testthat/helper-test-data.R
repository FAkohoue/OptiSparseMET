# Helper function for unit tests: generate a self-consistent minimal dataset
# that exercises the full OptiSparseMET pipeline without the runtime cost of
# the full 120-treatment example data. All components share the same treatment
# IDs, environments, and check set so they can be passed to any package
# function in combination without modification.
#
# Usage: call make_example_sparsemet_data() at the top of each test file
# that needs a realistic but lightweight dataset.

make_example_sparsemet_data <- function(seed = 123) {
  set.seed(seed)
  
  # 80 treatments across 3 environments keeps test runtime short while still
  # producing non-trivial allocation and field layouts.
  n_treatments <- 80L
  treatments   <- paste0("L", sprintf("%03d", seq_len(n_treatments)))
  environments <- c("Env1", "Env2", "Env3")
  
  # 10 families with random unequal group sizes — sufficient to test
  # family-guided allocation and within-block adjacency control without
  # requiring the larger 16-family structure of the production example data.
  treatment_info <- data.frame(
    Treatment = treatments,
    Family    = paste0("F", sprintf("%02d", sample(1:10, n_treatments, replace = TRUE))),
    stringsAsFactors = FALSE
  )
  
  # 5 common treatments — small enough to leave meaningful sparse allocation
  # capacity in each environment after they are removed from the pool.
  common_treatments <- treatments[1:5]
  
  # Seed range chosen to produce all three replication outcomes
  # (replicated, unreplicated, excluded) at the per-plot requirements below.
  # Even spacing avoids ties in seed-priority ranking under "seed_available"
  # priority mode.
  seed_info <- data.frame(
    Treatment     = treatments,
    SeedAvailable = sample(seq(20, 120, by = 2), n_treatments, replace = TRUE),
    stringsAsFactors = FALSE
  )
  
  # Env2 has the highest per-plot requirement (10), so it will produce the
  # most excluded treatments when seed_info is passed to
  # assign_replication_by_seed(). Env3 has the lowest (6), useful for
  # testing the downgrade path with a larger feasible pool.
  seed_required_per_plot <- data.frame(
    Environment         = environments,
    SeedRequiredPerPlot = c(8, 10, 6),
    stringsAsFactors = FALSE
  )
  
  env_design_specs <- list(
    
    # Env1: p-rep design via prep_famoptg().
    # Up to 8 treatments replicated twice; remainder unreplicated or
    # downgraded if seed is insufficient for 2 plots. Grid 6 x 6 = 36 cells
    # accommodates 2 checks x 4 blocks + 8 x 2 + remaining single-plot entries.
    Env1 = list(
      design               = "prep_famoptg",
      replication_mode     = "p_rep",
      desired_replications = 2L,
      max_prep             = 8L,
      shortage_action      = "downgrade",
      check_treatments     = c("CHK1", "CHK2"),
      check_families       = c("CHECK", "CHECK"),
      n_blocks             = 4L,
      n_rows               = 6L,
      n_cols               = 6L,
      order                = "row",
      serpentine           = TRUE,
      cluster_source       = "Family",
      eval_efficiency      = FALSE,
      use_dispersion       = FALSE
    ),
    
    # Env2: alpha row-column stream design via alpha_rc_stream().
    # n_reps = 2 on an 8 x 10 grid. min_entry_slots_per_block = 6 ensures
    # incomplete blocks are not too small after inserting 2 checks per block.
    # entry_treatments and entry_families are injected internally by
    # plan_sparse_met_design() and must not appear in this spec.
    Env2 = list(
      design                    = "alpha_rc_stream",
      check_treatments          = c("CHK1", "CHK2"),
      check_families            = c("CHECK", "CHECK"),
      n_reps                    = 2L,
      n_rows                    = 8L,
      n_cols                    = 10L,
      order                     = "row",
      serpentine                = TRUE,
      min_entry_slots_per_block = 6L,
      cluster_source            = "Family",
      eval_efficiency           = FALSE,
      use_dispersion            = FALSE,
      verbose                   = FALSE
    ),
    
    # Env3: augmented repeated-check design via prep_famoptg().
    # All allocated entries unreplicated — exercises the augmented path
    # and the column-major non-serpentine traversal branch.
    # 5 x 6 = 30 cells fits 2 checks x 4 blocks + remaining single-plot entries.
    Env3 = list(
      design           = "prep_famoptg",
      replication_mode = "augmented",
      check_treatments = c("CHK1", "CHK2"),
      check_families   = c("CHECK", "CHECK"),
      n_blocks         = 4L,
      n_rows           = 5L,
      n_cols           = 6L,
      order            = "column",
      serpentine       = FALSE,
      cluster_source   = "Family",
      eval_efficiency  = FALSE,
      use_dispersion   = FALSE
    )
  )
  
  list(
    treatments             = treatments,
    environments           = environments,
    treatment_info         = treatment_info,
    common_treatments      = common_treatments,
    seed_info              = seed_info,
    seed_required_per_plot = seed_required_per_plot,
    env_design_specs       = env_design_specs
  )
}