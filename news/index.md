# Changelog

## OptiSparseMET 0.1.0

### Initial release

`OptiSparseMET` introduces a unified framework for sparse
multi-environment trial (MET) design, jointly addressing treatment
allocation across environments and within-environment field design under
shared statistical, genetic, and logistical constraints.

The package links allocation structure, genetic connectivity, seed
availability, and spatial design assumptions in a single reproducible
workflow compatible with mixed-model inference.

------------------------------------------------------------------------

### Across-environment allocation

- Added
  [`allocate_sparse_met()`](https://FAkohoue.github.io/OptiSparseMET/reference/allocate_sparse_met.md)
  for constructing treatment-by-environment incidence matrices using
  sparse testing principles.

  - Supports `"random_balanced"` (M3) for flexible approximate balance
    with coverage-first guarantees.
  - Supports `"balanced_incomplete"` (M4) for BIBD-inspired uniform
    replication structure with enforced equal replication and equal
    environment sizes across the trial.
  - Accepts `"M3"` and `"M4"` as convenient aliases.
  - Supports unequal environment capacities under `random_balanced`.
  - Supports common treatments to ensure design-based connectivity
    across environments.
  - Returns `$allocation_matrix` (binary treatment-by-environment
    incidence matrix) from which pairwise co-occurrence can be computed
    post-hoc as `out$allocation_matrix %*% t(out$allocation_matrix)`.

- Added
  [`check_balanced_incomplete_feasibility()`](https://FAkohoue.github.io/OptiSparseMET/reference/check_balanced_incomplete_feasibility.md)
  to verify the slot identity (J\* × r = I × k\*) before attempting
  balanced incomplete allocation, confirming that equal replication is
  achievable for the chosen dimensions.

- Added
  [`derive_allocation_groups()`](https://FAkohoue.github.io/OptiSparseMET/reference/derive_allocation_groups.md)
  to construct grouping structures from:

  - family membership labels
  - genomic relationship matrix (GRM)
  - pedigree relationship matrix (A matrix)

  Group-guided allocation improves genetic connectedness and stability
  of cross-environment inference.

------------------------------------------------------------------------

### Capacity and feasibility helpers

- Added
  [`suggest_safe_k()`](https://FAkohoue.github.io/OptiSparseMET/reference/suggest_safe_k.md)
  to propose a safe uniform value for `n_test_entries_per_environment`
  given treatment count, environment count, common treatments, and a
  user-defined buffer.
- Added
  [`min_k_for_full_coverage()`](https://FAkohoue.github.io/OptiSparseMET/reference/min_k_for_full_coverage.md)
  to compute the minimum per-environment capacity required for every
  non-common treatment to be assigned at least once.
- Added
  [`warn_if_k_too_small()`](https://FAkohoue.github.io/OptiSparseMET/reference/warn_if_k_too_small.md)
  to provide a non-fatal diagnostic warning when the chosen capacity is
  insufficient for full treatment coverage.

These helpers prevent the most common failure mode: passing a capacity
too small to assign all treatments before
[`allocate_sparse_met()`](https://FAkohoue.github.io/OptiSparseMET/reference/allocate_sparse_met.md)
is called.

------------------------------------------------------------------------

### Seed-aware replication planning

- Added
  [`assign_replication_by_seed()`](https://FAkohoue.github.io/OptiSparseMET/reference/assign_replication_by_seed.md)
  to determine feasible replication levels based on available seed
  quantities and per-plot seed requirements.
  - Supports `"augmented"`, `"p_rep"`, and `"rcbd_type"` replication
    modes.
  - Supports `shortage_action` values `"error"`, `"downgrade"`, and
    `"exclude"` for handling treatments with insufficient seed.
  - Returns a role-partitioned list (`p_rep_treatments`,
    `unreplicated_treatments`, `excluded_treatments`) suitable for
    direct input to
    [`met_prep_famoptg()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_prep_famoptg.md).

------------------------------------------------------------------------

### Within-environment field design engines

#### Block-based repeated-check designs

- Added
  [`met_prep_famoptg()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_prep_famoptg.md)
  for constructing augmented, partially replicated (p-rep), and
  RCBD-type repeated-check block designs.

  Key structural guarantees:

  - Check treatments appear in every block.
  - Replicated (p-rep) treatments appear at most once per block.
  - Unreplicated treatments appear exactly once across the field.
  - Optional genetic dispersion optimisation using GRM or A matrix.
  - Optional within-environment efficiency evaluation (A, D, CDmean).

- Added
  [`met_evaluate_famoptg_efficiency()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_evaluate_famoptg_efficiency.md)
  for evaluating the statistical efficiency of
  [`met_prep_famoptg()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_prep_famoptg.md)
  designs under:

  - Fixed or random treatment effects.
  - IID, AR1, or AR1×AR1 residual covariance structures.
  - A-optimality, D-efficiency, CDmean, and mean PEV criteria.
  - Requires `sigma_b2` (block variance) in `varcomp`.

- Added
  [`met_optimize_famoptg()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_optimize_famoptg.md)
  for criterion-driven optimisation of
  [`met_prep_famoptg()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_prep_famoptg.md)
  designs via Random Restart.

#### Row-column alpha designs

- Added
  [`met_alpha_rc_stream()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_alpha_rc_stream.md)
  for generating alpha row-column stream designs suitable for large
  structured fields.

  Key features:

  - Repeated checks in every incomplete block.
  - Configurable block sizes via `min_block_size` / `max_block_size` or
    fixed `n_blocks_per_rep`.
  - Row-major, column-major, and serpentine traversal orders.
  - Optional genetic grouping from family, GRM, or A matrix.
  - Optional within-environment efficiency evaluation.

- Added
  [`met_evaluate_alpha_efficiency()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_evaluate_alpha_efficiency.md)
  for evaluating the statistical efficiency of
  [`met_alpha_rc_stream()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_alpha_rc_stream.md)
  designs under the same criteria as
  [`met_evaluate_famoptg_efficiency()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_evaluate_famoptg_efficiency.md).

  - Requires `sigma_rep2` (replicate variance) and `sigma_ib2`
    (incomplete block within replicate variance) in `varcomp`,
    distinguishing it from the block-based evaluator.

- Added
  [`met_optimize_alpha_rc()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_optimize_alpha_rc.md)
  for criterion-driven optimisation of
  [`met_alpha_rc_stream()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_alpha_rc_stream.md)
  designs via Random Restart, Simulated Annealing, or Genetic Algorithm.

------------------------------------------------------------------------

### Pipeline orchestration

- Added
  [`plan_sparse_met_design()`](https://FAkohoue.github.io/OptiSparseMET/reference/plan_sparse_met_design.md)
  providing an end-to-end MET design workflow that integrates:

  - Across-environment allocation via
    [`allocate_sparse_met()`](https://FAkohoue.github.io/OptiSparseMET/reference/allocate_sparse_met.md).
  - Per-environment seed feasibility via
    [`assign_replication_by_seed()`](https://FAkohoue.github.io/OptiSparseMET/reference/assign_replication_by_seed.md).
  - Environment-specific field design via
    [`met_prep_famoptg()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_prep_famoptg.md)
    or
    [`met_alpha_rc_stream()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_alpha_rc_stream.md),
    selected by `design` in `env_design_specs`.
  - Mixed design strategies across environments in a single call.

  Each environment’s design engine is specified via `env_design_specs`,
  a named list where `design = "met_prep_famoptg"` or
  `design = "met_alpha_rc_stream"` selects the constructor.

- Added
  [`combine_met_fieldbooks()`](https://FAkohoue.github.io/OptiSparseMET/reference/combine_met_fieldbooks.md)
  to stack environment-level field books produced by
  [`met_prep_famoptg()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_prep_famoptg.md)
  or
  [`met_alpha_rc_stream()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_alpha_rc_stream.md)
  into a unified MET-level field book with standard metadata columns
  (`Environment`, `LocalDesign`, `ReplicationMode`, `SparseMethod`,
  `IsCommonTreatment`). Handles heterogeneous column sets across
  environments by filling missing columns with `NA`.

------------------------------------------------------------------------

### Statistical framework

Implements the sparse testing identity:

``` R
N = J × r = I × k
```

making explicit the tradeoff between number of treatments (J), number of
environments (I), replication depth (r), and treatments per environment
(k).

Design outputs are compatible with mixed-model analysis:

``` R
y = Xβ + Zg + e
```

where the covariance of g may be proportional to a GRM or pedigree A
matrix, and the residual covariance structure may be IID, AR1, or
AR1×AR1.

Supports design strategies that improve: - cross-environment genetic
connectivity - G×E estimation stability - genomic prediction performance
(CDmean criterion; Rincent et al. 2012) - precision of BLUP estimates

------------------------------------------------------------------------

### Efficiency diagnostics

Within-environment designs support optional efficiency evaluation,
summarized via
[`plan_sparse_met_design()`](https://FAkohoue.github.io/OptiSparseMET/reference/plan_sparse_met_design.md)
across all environments:

- A-criterion and A-efficiency
- D-criterion and D-efficiency
- Mean prediction error variance (mean PEV)
- CDmean (coefficient of determination for genomic prediction)

Metrics are reported in `$efficiency_summary` (long format) and
`$environment_summary` (wide format with `has_efficiency`, `eff_A`,
`eff_D`, `eff_mean_PEV` columns).

------------------------------------------------------------------------

### Infrastructure

- Initial package structure with all exports documented via roxygen2.
- Bundled example dataset `OptiSparseMET_example_data` with 120
  treatments, 4 environments, GRM, pedigree A matrix, prediction matrix
  K, seed availability data, pre-built allocation argument lists, and
  environment-specific design specifications.
- Vignette describing the statistical framework, two-stage pipeline, and
  worked examples.
- Unit test suite covering all 13 exported functions plus internal
  helpers.
- pkgdown configuration for the documentation website.
- GitHub Actions workflows for R CMD check and pkgdown deployment.
