# benchmark_designs.R  (P3.6)
# ------------------------------------------------------------------------------
# Reproducible benchmarks of OptiSparseMET against reference tools. Run in an
# environment where the optional packages are installed. This script is NOT part
# of R CMD check; it documents how to validate the criteria and optimizers
# against independent implementations.
#
#   * sommer  -- cross-check design-based PEV against a fitted GBLUP model
#   * odw / DiGGer / blocksdesign -- compare within-field A-optimality and
#     realized accuracy for alpha / p-rep layouts
#
# Usage: Rscript scripts/benchmark_designs.R
# ------------------------------------------------------------------------------

library(OptiSparseMET)

## ---- 1. Design-based PEV vs sommer (GBLUP) -----------------------------------
if (requireNamespace("sommer", quietly = TRUE)) {
  message("Benchmark 1: met_information PEV vs sommer GBLUP PEV")
  # 1. Build a small allocation with allocate_sparse_met().
  # 2. Compute across-TPE PEV with met_information().
  # 3. Simulate one MET dataset (simulate_met internals) and fit the SAME model
  #    in sommer with FIXED variance components (not estimated), then compare
  #    mean PEV. They should agree within numerical tolerance.
  # (Fill in with your program's G and Sigma_E.)
} else {
  message("Skipping sommer benchmark: package not installed.")
}

## ---- 2. Within-field A-optimality vs odw / DiGGer ----------------------------
if (requireNamespace("odw", quietly = TRUE) ||
    requireNamespace("DiGGer", quietly = TRUE)) {
  message("Benchmark 2: met_optimize_alpha_rc A-criterion vs reference engine")
  # 1. Construct an alpha row-column design with met_alpha_rc_stream().
  # 2. Evaluate A_criterion with met_evaluate_alpha_efficiency().
  # 3. Build a comparable design with odw / DiGGer and compare A-optimality and
  #    realized accuracy (via simulate_met on the resulting incidence).
} else {
  message("Skipping odw/DiGGer benchmark: packages not installed.")
}

## ---- 3. Allocation: balanced vs sparse vs near-balanced -----------------------
message("Benchmark 3: realized accuracy across allocation strategies")
# Use simulate_met() to compare accuracy/gain for:
#   * fully balanced (all lines in all environments),
#   * random sparse,
#   * equireplicate with balance = "env_pair",
#   * optimize_allocation_gxe() refinement,
# under a broad-TPE Sigma_E. Expect near-balanced / optimized >= random sparse.
