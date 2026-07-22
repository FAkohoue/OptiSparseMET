# OptiSparseMET 0.2.0 (development)

## Correctness fixes (P0)

* **Random-treatment PEV/CDmean fixed-effect offset (major).** In
  `met_evaluate_alpha_efficiency()` and `met_evaluate_famoptg_efficiency()`, the
  genotype index used to extract prediction-error variances from the coefficient
  matrix `C = [[X'QX, X'QZ], [Z'QX, Z'QZ + G^{-1}]]` was indexed within the
  random block `Z`, but the block is offset by the number of fixed-effect
  columns. As a result the reported `mean_PEV` / `CDmean` / `CD_per_line`
  averaged over a window shifted by `ncol(X)` -- including a nuisance random
  effect's PEV and dropping the last genotype(s). This affected every
  random-treatment (genomic-selection) evaluation and the `CDmean`-criterion
  optimisers. Fixed by offsetting the extraction index by `ncol(X)`; verified
  against an independent dense MME solve (single-cell example: corrected
  `mean_PEV = 0.55`, previously `0.6625`).

* **AR1 / AR1×AR1 residual precision**: fixed a transposed Kronecker index that
  swapped the row/column autocorrelation axes and scrambled spatial adjacency on
  non-square fields. Residual-precision construction is now centralised in the
  tested helper `.build_residual_precision()` and verified against a hand-built
  AR1×AR1 precision (`test-ar1-precision.R`). Affects every AR1-based A/D/CDmean
  value in `met_evaluate_alpha_efficiency()` and `met_evaluate_famoptg_efficiency()`.
* **CDmean scaling**: per-line reliability now divides by `sigma_g2 * K_ii`
  (diagonal of the relationship matrix) instead of `sigma_g2` alone, correcting
  the value for VanRaden GRMs and inbred/pedigree material where `diag(K) != 1`.
* **D-criterion**: `.safe_logdet_psd_dense()` now uses a relative eigenvalue
  tolerance and an exact expected rank, making the pseudo-determinant
  scale-robust.
* **Rank-deficiency guard**: a singular fixed-effects design (e.g. intercept plus
  a full dummy set with no check rows) is now repaired instead of erroring.

## Spatial models

* Residual structures extended beyond `IID` / `AR1` / `AR1xAR1`:
  * **`"AR1xAR1_nugget"`** adds an independent measurement-error (nugget)
    component. Pure AR1xAR1 often fits field data poorly, so this is usually the
    more realistic applied model.
  * **`"exponential"`, `"gaussian"`, `"matern"`** give isotropic correlation as a
    function of Euclidean plot distance, for irregular or non-rectangular
    layouts. These build the covariance of the observed plots and invert it (the
    marginal form), so they also handle designs with empty cells correctly.
* New `spatial_random = "pspline"` fits a **SpATS-style two-dimensional
  penalised-spline surface** as a random effect (Rodriguez-Alvarez et al. 2018)
  rather than as a residual covariance. The tensor B-spline penalty is improper
  (its null space is constant plus linear trends), so the block is supplied in
  the mixed-model reparameterisation: only the positive-eigenvalue space is kept,
  giving a proper random effect whose surface shrinks away as the smoothing
  parameters grow. Available in both evaluators and both optimisers.
* New `compare_spatial_models()` scores one design under several spatial models
  side by side -- the spatial analogue of `sensitivity_varcomp()`. A design that
  wins under every model is a safe choice; one that wins under a single model is
  a warning.

## Validation and coupling (P3)

* `met_evaluate_alpha_efficiency()` / `met_evaluate_famoptg_efficiency()` now
  also return `A_efficiency_rel`, the classical **relative** A-efficiency in
  (0, 1] (vs an orthogonal CRD on the same plots). The existing `A_efficiency`
  / `$A` fields are documented as the inverse criterion, not relative efficiency.
* New `sensitivity_varcomp()`: sweep the residual/genetic variance ratio and
  report how across-TPE PEV and CDmean respond, so a design choice can be judged
  for robustness to variance misspecification.
* The genuine contrast-based CDmean of Rincent et al. (2012) is available via
  `select_individuals()`; the evaluators' per-line quantity is now correctly
  labeled as reliability (see P0 CDmean fix) rather than mislabeled as the
  Rincent CDmean.

* New `met_information()`: the across-environment (MET-level) coupled
  information matrix. It combines the allocation incidence (Level 1) with a
  within-environment efficiency factor (Level 2) under a G x E model
  (`Sigma_E ⊗ G`) and returns the across-TPE mean PEV and CDmean. This backs the
  "coupled two-level" claim with code; the fixed-effect absorption is verified
  to match a full mixed-model-equations solve.
* New `simulate_met()`: simulates a MET with a genuine across-environment
  genetic correlation structure, predicts by GBLUP, and reports realized
  selection accuracy and genetic gain -- linking A/D/CDmean design proxies to
  the outcome they are meant to improve (replaces the previous no-G×E
  vignette demonstration).

## Decision-point framework (P4, in progress)

* New `select_environments()` (decision 1): choose a representative subset of
  environments from a TPE. Default `"representative"` maximises TPE coverage
  (k-medoid), correcting the spread-only `-q'Dq` criterion (available as
  `"optcontrib"`); `"kmeans"`, `"hclust"`, and `"random"` baselines included.
* New `select_individuals()` (decision 2): choose a representative training
  subset of the TPG by maximising CDmean (Rincent et al. 2012) via an exchange
  algorithm.
* New `optimize_allocation_gxe()` (decision 3): refine an equireplicate
  allocation to minimise across-TPE PEV under the coupled G x E information
  matrix -- the reproducible, CDmean-based replacement for the paper's
  unbounded `-q'Dq` / connectedness objectives, preserving replication and
  environment-size margins.
* New `recommend_replication()` (decision 4): sweep replication levels for an
  allocation, report realized accuracy / gain / PEV via `simulate_met()`, honour
  a seed budget, and recommend the diminishing-returns level -- de-confounding
  replication from the number of unique entries.

## Polish

* `met_prep_famoptg()` now emits the field-dimension auto-correction as an
  informational message under `verbose` (matching its documented "silently
  adjusts" behaviour) instead of a `warning()`.
* `plan_sparse_met_design()` accepts `allocation_method = "equireplicate"`
  (and the `"M4"` alias).

## Allocation (P1)

* Renamed `allocation_method = "balanced_incomplete"` to `"equireplicate"`
  (the old name was a misnomer and has been **removed**, not deprecated;
  `"M4"` remains as an alias): the constructor guarantees equal replication and
  equal environment size, not a strict BIBD. The bundled example data argument
  list is likewise renamed `sparse_example_args_equireplicate`.
* New `balance` argument (`"env_pair"`, `"line_pair"`, `"both"`) drives the
  design toward a near-balanced (regular-graph) structure while preserving the
  replication and environment-size margins; returns a `balance_report`.
* `check_equireplicate_feasibility()` (and the retained
  `check_balanced_incomplete_feasibility()`) now report the strict-BIBD
  existence conditions: `bibd_lambda`, `bibd_lambda_integer`, `fisher_ok`, and
  `strict_bibd_possible`.
* New `construct_exact_bibd()` builds a true BIBD via `crossdes` when one
  provably exists (small-J regime).

# OptiSparseMET 0.1.0

## Initial release

`OptiSparseMET` introduces a unified framework for sparse multi-environment
trial (MET) design, jointly addressing treatment allocation across environments
and within-environment field design under shared statistical, genetic, and
logistical constraints.

The package links allocation structure, genetic connectivity, seed
availability, and spatial design assumptions in a single reproducible workflow
compatible with mixed-model inference.

---

## Across-environment allocation

- Added `allocate_sparse_met()` for constructing treatment-by-environment
  incidence matrices using sparse testing principles.
  - Supports `"random_balanced"` (M3) for flexible approximate balance with
    coverage-first guarantees.
  - Supports `"balanced_incomplete"` (M4) for BIBD-inspired uniform
    replication structure with enforced equal replication and equal
    environment sizes across the trial.
  - Accepts `"M3"` and `"M4"` as convenient aliases.
  - Supports unequal environment capacities under `random_balanced`.
  - Supports common treatments to ensure design-based connectivity across
    environments.
  - Returns `$allocation_matrix` (binary treatment-by-environment incidence
    matrix) from which pairwise co-occurrence can be computed post-hoc as
    `out$allocation_matrix %*% t(out$allocation_matrix)`.

- Added `check_balanced_incomplete_feasibility()` to verify the slot identity
  (J* × r = I × k*) before attempting balanced incomplete allocation,
  confirming that equal replication is achievable for the chosen dimensions.

- Added `derive_allocation_groups()` to construct grouping structures from:
  - family membership labels
  - genomic relationship matrix (GRM)
  - pedigree relationship matrix (A matrix)

  Group-guided allocation improves genetic connectedness and stability of
  cross-environment inference.

---

## Capacity and feasibility helpers

- Added `suggest_safe_k()` to propose a safe uniform value for
  `n_test_entries_per_environment` given treatment count, environment count,
  common treatments, and a user-defined buffer.
- Added `min_k_for_full_coverage()` to compute the minimum per-environment
  capacity required for every non-common treatment to be assigned at least
  once.
- Added `warn_if_k_too_small()` to provide a non-fatal diagnostic warning when
  the chosen capacity is insufficient for full treatment coverage.

These helpers prevent the most common failure mode: passing a capacity too
small to assign all treatments before `allocate_sparse_met()` is called.

---

## Seed-aware replication planning

- Added `assign_replication_by_seed()` to determine feasible replication
  levels based on available seed quantities and per-plot seed requirements.
  - Supports `"augmented"`, `"p_rep"`, and `"rcbd_type"` replication modes.
  - Supports `shortage_action` values `"error"`, `"downgrade"`, and
    `"exclude"` for handling treatments with insufficient seed.
  - Returns a role-partitioned list (`p_rep_treatments`,
    `unreplicated_treatments`, `excluded_treatments`) suitable for direct
    input to `met_prep_famoptg()`.

---

## Within-environment field design engines

### Block-based repeated-check designs

- Added `met_prep_famoptg()` for constructing augmented, partially replicated
  (p-rep), and RCBD-type repeated-check block designs.

  Key structural guarantees:
  - Check treatments appear in every block.
  - Replicated (p-rep) treatments appear at most once per block.
  - Unreplicated treatments appear exactly once across the field.
  - Optional genetic dispersion optimisation using GRM or A matrix.
  - Optional within-environment efficiency evaluation (A, D, CDmean).

- Added `met_evaluate_famoptg_efficiency()` for evaluating the statistical
  efficiency of `met_prep_famoptg()` designs under:
  - Fixed or random treatment effects.
  - IID, AR1, or AR1×AR1 residual covariance structures.
  - A-optimality, D-efficiency, CDmean, and mean PEV criteria.
  - Requires `sigma_b2` (block variance) in `varcomp`.

- Added `met_optimize_famoptg()` for criterion-driven optimisation of
  `met_prep_famoptg()` designs via Random Restart.

### Row-column alpha designs

- Added `met_alpha_rc_stream()` for generating alpha row-column stream designs
  suitable for large structured fields.

  Key features:
  - Repeated checks in every incomplete block.
  - Configurable block sizes via `min_block_size` / `max_block_size` or fixed
    `n_blocks_per_rep`.
  - Row-major, column-major, and serpentine traversal orders.
  - Optional genetic grouping from family, GRM, or A matrix.
  - Optional within-environment efficiency evaluation.

- Added `met_evaluate_alpha_efficiency()` for evaluating the statistical
  efficiency of `met_alpha_rc_stream()` designs under the same criteria as
  `met_evaluate_famoptg_efficiency()`.
  - Requires `sigma_rep2` (replicate variance) and `sigma_ib2` (incomplete
    block within replicate variance) in `varcomp`, distinguishing it from the
    block-based evaluator.

- Added `met_optimize_alpha_rc()` for criterion-driven optimisation of
  `met_alpha_rc_stream()` designs via Random Restart, Simulated Annealing, or
  Genetic Algorithm.

---

## Pipeline orchestration

- Added `plan_sparse_met_design()` providing an end-to-end MET design workflow
  that integrates:
  - Across-environment allocation via `allocate_sparse_met()`.
  - Per-environment seed feasibility via `assign_replication_by_seed()`.
  - Environment-specific field design via `met_prep_famoptg()` or
    `met_alpha_rc_stream()`, selected by `design` in `env_design_specs`.
  - Mixed design strategies across environments in a single call.

  Each environment's design engine is specified via `env_design_specs`, a
  named list where `design = "met_prep_famoptg"` or
  `design = "met_alpha_rc_stream"` selects the constructor.

- Added `combine_met_fieldbooks()` to stack environment-level field books
  produced by `met_prep_famoptg()` or `met_alpha_rc_stream()` into a unified
  MET-level field book with standard metadata columns (`Environment`,
  `LocalDesign`, `ReplicationMode`, `SparseMethod`, `IsCommonTreatment`).
  Handles heterogeneous column sets across environments by filling missing
  columns with `NA`.

---

## Statistical framework

Implements the sparse testing identity:

    N = J × r = I × k

making explicit the tradeoff between number of treatments (J), number of
environments (I), replication depth (r), and treatments per environment (k).

Design outputs are compatible with mixed-model analysis:

    y = Xβ + Zg + e

where the covariance of g may be proportional to a GRM or pedigree A matrix,
and the residual covariance structure may be IID, AR1, or AR1×AR1.

Supports design strategies that improve:
- cross-environment genetic connectivity
- G×E estimation stability
- genomic prediction performance (CDmean criterion; Rincent et al. 2012)
- precision of BLUP estimates

---

## Efficiency diagnostics

Within-environment designs support optional efficiency evaluation, summarized
via `plan_sparse_met_design()` across all environments:

- A-criterion and A-efficiency
- D-criterion and D-efficiency
- Mean prediction error variance (mean PEV)
- CDmean (coefficient of determination for genomic prediction)

Metrics are reported in `$efficiency_summary` (long format) and
`$environment_summary` (wide format with `has_efficiency`, `eff_A`, `eff_D`,
`eff_mean_PEV` columns).

---

## Infrastructure

- Initial package structure with all exports documented via roxygen2.
- Bundled example dataset `OptiSparseMET_example_data` with 120 treatments,
  4 environments, GRM, pedigree A matrix, prediction matrix K, seed
  availability data, pre-built allocation argument lists, and
  environment-specific design specifications.
- Vignette describing the statistical framework, two-stage pipeline, and
  worked examples.
- Unit test suite covering all 13 exported functions plus internal helpers.
- pkgdown configuration for the documentation website.
- GitHub Actions workflows for R CMD check and pkgdown deployment.