# Example data for OptiSparseMET

Synthetic dataset for demonstrating and testing sparse multi-environment
trial allocation and within-environment field design construction across
the functions of `OptiSparseMET`. All components are internally
consistent: treatment IDs, family labels, seed inventories, relationship
matrices, and environment specifications refer to the same set of
synthetic lines and environments, so any combination of components can
be passed to the package functions without modification.

## Usage

``` r
OptiSparseMET_example_data
```

## Format

A named list with the following components:

- `treatments`:

  Character vector of test treatment IDs. These are the candidate lines
  to be allocated across environments. Does not include check
  treatments.

- `environments`:

  Character vector of environment names. Used as the `environments`
  argument in
  [`allocate_sparse_met()`](https://FAkohoue.github.io/OptiSparseMET/reference/allocate_sparse_met.md)
  and as keys in `env_design_specs`.

- `treatment_info`:

  Data frame with one row per treatment in `treatments`. Contains at
  minimum columns `Treatment` and `Family`. Used as the `treatment_info`
  argument in
  [`allocate_sparse_met()`](https://FAkohoue.github.io/OptiSparseMET/reference/allocate_sparse_met.md)
  when `allocation_group_source = "Family"`.

- `common_treatments`:

  Character vector of treatment IDs forced into all environments before
  sparse allocation. A subset of `treatments`. Used as the
  `common_treatments` argument in
  [`allocate_sparse_met()`](https://FAkohoue.github.io/OptiSparseMET/reference/allocate_sparse_met.md)
  and
  [`combine_met_fieldbooks()`](https://FAkohoue.github.io/OptiSparseMET/reference/combine_met_fieldbooks.md).

- `seed_info`:

  Data frame with columns `Treatment` and `SeedAvailable`. Contains one
  row per treatment in `treatments` and gives the seed quantity
  available for each line. Used as the `seed_available` argument in
  [`assign_replication_by_seed()`](https://FAkohoue.github.io/OptiSparseMET/reference/assign_replication_by_seed.md).

- `seed_required_per_plot`:

  Data frame with columns `Environment` and `SeedRequiredPerPlot`. Gives
  the per-plot seed requirement for each environment, reflecting that
  seed consumption may differ across locations. Used as the
  `seed_required_per_plot` argument in
  [`assign_replication_by_seed()`](https://FAkohoue.github.io/OptiSparseMET/reference/assign_replication_by_seed.md)
  after subsetting to the environment of interest.

- `OptiSparseMET_GRM`:

  Numeric matrix. Genomic relationship matrix with row and column names
  matching the IDs in `treatments`. Used as the `GRM` argument in
  [`allocate_sparse_met()`](https://FAkohoue.github.io/OptiSparseMET/reference/allocate_sparse_met.md)
  when `allocation_group_source = "GRM"`, and in
  [`alpha_rc_stream()`](https://FAkohoue.github.io/OptiSparseMET/reference/alpha_rc_stream.md)
  or
  [`prep_famoptg()`](https://FAkohoue.github.io/OptiSparseMET/reference/prep_famoptg.md)
  when `cluster_source = "GRM"` or `dispersion_source = "GRM"`.

- `OptiSparseMET_A`:

  Numeric matrix. Pedigree-based numerator relationship matrix with row
  and column names matching the IDs in `treatments`. Used as the `A`
  argument when `allocation_group_source = "A"`, `cluster_source = "A"`,
  or `dispersion_source = "A"`.

- `OptiSparseMET_K`:

  Numeric matrix. Relationship matrix intended for mixed-model
  efficiency evaluation and dispersion optimization. Used as the `K`
  argument in
  [`alpha_rc_stream()`](https://FAkohoue.github.io/OptiSparseMET/reference/alpha_rc_stream.md)
  and
  [`prep_famoptg()`](https://FAkohoue.github.io/OptiSparseMET/reference/prep_famoptg.md)
  when `prediction_type %in% c("GBLUP", "PBLUP")` or
  `dispersion_source = "K"`.

- `sparse_example_args_random_balanced`:

  Named list of arguments ready to pass to
  [`allocate_sparse_met()`](https://FAkohoue.github.io/OptiSparseMET/reference/allocate_sparse_met.md)
  with `allocation_method = "random_balanced"`. All list names
  correspond directly to function parameter names.

- `sparse_example_args_balanced_incomplete`:

  Named list of arguments ready to pass to
  [`allocate_sparse_met()`](https://FAkohoue.github.io/OptiSparseMET/reference/allocate_sparse_met.md)
  with `allocation_method = "balanced_incomplete"`. All list names
  correspond directly to function parameter names.

- `env_design_specs`:

  Named list with one element per environment in `environments`. Each
  element is itself a named list of arguments specifying the local field
  design for that environment, suitable for passing to
  [`prep_famoptg()`](https://FAkohoue.github.io/OptiSparseMET/reference/prep_famoptg.md)
  or
  [`alpha_rc_stream()`](https://FAkohoue.github.io/OptiSparseMET/reference/alpha_rc_stream.md).

## Source

Synthetically generated for package documentation and testing.

## Details

All data are synthetic and generated solely for use in package examples,
unit tests, and vignettes. No real breeding trial data are included.
Relationship matrices (`OptiSparseMET_GRM`, `OptiSparseMET_A`,
`OptiSparseMET_K`) are positive semi-definite by construction and have
row and column names consistent with `treatments`, so they can be used
directly without an `id_map`. The `seed_info` and
`seed_required_per_plot` components are calibrated so that a range of
replication outcomes (fully replicated, partially replicated, and
excluded treatments) arise naturally when passed to
[`assign_replication_by_seed()`](https://FAkohoue.github.io/OptiSparseMET/reference/assign_replication_by_seed.md),
making the data useful for demonstrating all three replication modes.
