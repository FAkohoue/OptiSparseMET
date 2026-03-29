# ==============================================================================
# helper-test-data.R
# Shared fixtures for the OptiSparseMET test suite.
#
# All make_* functions return minimal valid inputs for their respective
# function. Tests modify copies rather than duplicating boilerplate.
# This file is auto-sourced by testthat before each test file.
# ==============================================================================

# -- Allocation fixture --------------------------------------------------------
make_alloc_args <- function(seed = 1L) {
  list(
    treatments   = paste0("L", sprintf("%03d", 1:30)),
    environments = c("E1", "E2", "E3"),
    allocation_method              = "random_balanced",
    n_test_entries_per_environment = 12L,
    target_replications            = 1L,
    seed                           = seed
  )
}

# -- Seed replication fixture --------------------------------------------------
make_seed_args <- function() {
  trt <- paste0("L", sprintf("%03d", 1:10))
  list(
    treatments = trt,
    seed_available = data.frame(
      Treatment     = trt,
      SeedAvailable = c(50, 45, 40, 35, 30, 20, 18, 15, 12, 10),
      stringsAsFactors = FALSE
    ),
    seed_required_per_plot = 10
  )
}

# -- met_prep_famoptg fixture --------------------------------------------------
# 2 checks x 3 blocks = 6 check plots
# 6 p-rep treatments x 2 reps   = 12 p-rep plots
# 18 unreplicated treatments x 1 = 18 unrep plots
# Total = 36 plots -> 6 x 6 field
make_famoptg_args <- function(seed = 1L) {
  list(
    check_treatments        = c("CHK1", "CHK2"),
    check_families          = c("CHECK", "CHECK"),
    p_rep_treatments        = paste0("P", 1:6),
    p_rep_reps              = rep(2L, 6L),
    p_rep_families          = rep(c("F1", "F2"), 3),
    unreplicated_treatments = paste0("U", 1:18),
    unreplicated_families   = rep(c("F1", "F2", "F3"), 6),
    n_blocks                = 3L,
    n_rows                  = 6L,
    n_cols                  = 6L,
    seed                    = seed,
    verbose                 = FALSE
  )
}

# -- met_alpha_rc_stream fixture -----------------------------------------------
# 3 checks, 20 entries, 2 reps
# min/max block size 10/12 -> block sizes 10-11 entries + 3 checks = 13-14 total
# used plots = 2 * (20 + b*3); auto b ~ ceil(20/9)=3 -> 2*(20+9)=58 plots in 6x10=60
make_alpha_args <- function(seed = 1L) {
  list(
    check_treatments = c("CHK1", "CHK2", "CHK3"),
    check_families   = c("CHECK", "CHECK", "CHECK"),
    entry_treatments = paste0("G", 1:20),
    entry_families   = rep(c("F1", "F2"), 10),
    n_reps           = 2L,
    n_rows           = 6L,
    n_cols           = 10L,
    min_block_size   = 8L,
    max_block_size   = 12L,
    seed             = seed,
    verbose          = FALSE
  )
}

# -- Efficiency varcomp helpers ------------------------------------------------
famoptg_varcomp <- function() {
  list(sigma_e2 = 1, sigma_g2 = 1, sigma_b2 = 1, sigma_r2 = 0.1, sigma_c2 = 0.1)
}

alpha_varcomp <- function() {
  list(sigma_e2 = 1, sigma_g2 = 1, sigma_rep2 = 0.5,
       sigma_ib2 = 0.3, sigma_r2 = 0.1, sigma_c2 = 0.1)
}

# -- Relationship matrix helper ------------------------------------------------
make_sym_pd <- function(n, seed = 99L) {
  set.seed(seed)
  raw <- matrix(rnorm(n * n), n, n)
  K   <- crossprod(raw) / n
  diag(K) <- diag(K) + 0.5
  K
}