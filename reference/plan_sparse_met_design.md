# Plan a sparse multi-environment trial design and assemble a combined field book

`plan_sparse_met_design()` is the end-to-end workflow wrapper for sparse
multi-environment trial (MET) design. It first allocates treatments
across environments using
[`allocate_sparse_met()`](https://FAkohoue.github.io/OptiSparseMET/reference/allocate_sparse_met.md),
then builds a local field design for each environment using the
environment-specific specification in `env_design_specs`, and finally
stacks all local field books into one combined MET field book using
[`combine_met_fieldbooks()`](https://FAkohoue.github.io/OptiSparseMET/reference/combine_met_fieldbooks.md).

## Usage

``` r
plan_sparse_met_design(
  treatments,
  environments,
  allocation_method = c("random_balanced", "balanced_incomplete", "M3", "M4"),
  n_test_entries_per_environment,
  target_replications = NULL,
  common_treatments = NULL,
  env_design_specs,
  treatment_info = NULL,
  seed_info = NULL,
  seed_required_per_plot = NULL,
  allocation_group_source = c("none", "Family", "GRM", "A"),
  GRM = NULL,
  A = NULL,
  id_map = NULL,
  group_method = c("kmeans", "hclust"),
  group_seed = 1,
  group_attempts = 25,
  n_pcs_use = Inf,
  min_groups_per_environment = NULL,
  min_env_per_group = NULL,
  balance_groups_across_env = TRUE,
  force_group_connectivity = TRUE,
  allow_approximate = FALSE,
  seed = NULL
)
```

## Arguments

- treatments:

  Character vector of treatment IDs.

- environments:

  Character vector of environment names.

- allocation_method:

  Character scalar passed to
  [`allocate_sparse_met()`](https://FAkohoue.github.io/OptiSparseMET/reference/allocate_sparse_met.md).

- n_test_entries_per_environment:

  Integer scalar or integer vector passed to
  [`allocate_sparse_met()`](https://FAkohoue.github.io/OptiSparseMET/reference/allocate_sparse_met.md).

- target_replications:

  Optional positive integer passed to
  [`allocate_sparse_met()`](https://FAkohoue.github.io/OptiSparseMET/reference/allocate_sparse_met.md).

- common_treatments:

  Optional character vector of common treatments.

- env_design_specs:

  Named list with one element per environment. Each element must contain
  at least a `design` field selecting the local design engine
  (`"prep_famoptg"` or `"alpha_rc_stream"`). All other fields are passed
  through to the corresponding local design function, except a small set
  of pipeline-consumed fields used internally for seed-aware
  replication.

- treatment_info:

  Optional data frame with at least `Treatment` and optionally `Family`.

- seed_info:

  Optional data frame used by
  [`assign_replication_by_seed()`](https://FAkohoue.github.io/OptiSparseMET/reference/assign_replication_by_seed.md)
  when seed-aware replication is requested for a local design.

- seed_required_per_plot:

  Optional scalar, named vector, or named list used by
  [`assign_replication_by_seed()`](https://FAkohoue.github.io/OptiSparseMET/reference/assign_replication_by_seed.md)
  when seed-aware replication is requested.

- allocation_group_source:

  Character scalar passed to
  [`allocate_sparse_met()`](https://FAkohoue.github.io/OptiSparseMET/reference/allocate_sparse_met.md).

- GRM:

  Optional genomic relationship matrix passed to
  [`allocate_sparse_met()`](https://FAkohoue.github.io/OptiSparseMET/reference/allocate_sparse_met.md).

- A:

  Optional pedigree relationship matrix passed to
  [`allocate_sparse_met()`](https://FAkohoue.github.io/OptiSparseMET/reference/allocate_sparse_met.md).

- id_map:

  Optional treatment-to-matrix ID map passed to
  [`allocate_sparse_met()`](https://FAkohoue.github.io/OptiSparseMET/reference/allocate_sparse_met.md).

- group_method:

  Character scalar passed to
  [`allocate_sparse_met()`](https://FAkohoue.github.io/OptiSparseMET/reference/allocate_sparse_met.md).

- group_seed:

  Integer seed passed to
  [`allocate_sparse_met()`](https://FAkohoue.github.io/OptiSparseMET/reference/allocate_sparse_met.md).

- group_attempts:

  Integer passed to
  [`allocate_sparse_met()`](https://FAkohoue.github.io/OptiSparseMET/reference/allocate_sparse_met.md).

- n_pcs_use:

  Integer or `Inf` passed to
  [`allocate_sparse_met()`](https://FAkohoue.github.io/OptiSparseMET/reference/allocate_sparse_met.md).

- min_groups_per_environment:

  Optional integer passed to
  [`allocate_sparse_met()`](https://FAkohoue.github.io/OptiSparseMET/reference/allocate_sparse_met.md).

- min_env_per_group:

  Optional integer passed to
  [`allocate_sparse_met()`](https://FAkohoue.github.io/OptiSparseMET/reference/allocate_sparse_met.md).

- balance_groups_across_env:

  Logical passed to
  [`allocate_sparse_met()`](https://FAkohoue.github.io/OptiSparseMET/reference/allocate_sparse_met.md).

- force_group_connectivity:

  Logical passed to
  [`allocate_sparse_met()`](https://FAkohoue.github.io/OptiSparseMET/reference/allocate_sparse_met.md).

- allow_approximate:

  Logical passed to
  [`allocate_sparse_met()`](https://FAkohoue.github.io/OptiSparseMET/reference/allocate_sparse_met.md).

- seed:

  Optional integer seed for reproducibility.

## Value

A named list with:

- `sparse_allocation`:

  Output from
  [`allocate_sparse_met()`](https://FAkohoue.github.io/OptiSparseMET/reference/allocate_sparse_met.md).

- `environment_designs`:

  Named list of local design outputs, one per environment.

- `combined_field_book`:

  Combined MET field book from
  [`combine_met_fieldbooks()`](https://FAkohoue.github.io/OptiSparseMET/reference/combine_met_fieldbooks.md).

- `environment_summary`:

  Environment-level summary table including design metadata, plot
  counts, treatment counts, and extracted efficiency fields such as
  `eff_model`, `eff_A`, `eff_D`, and `eff_mean_PEV` when available.

- `group_environment_summary`:

  Group-by-environment summary propagated from
  `sparse_allocation$group_by_environment`.

- `efficiency_summary`:

  Long-format table of extracted efficiency outputs with columns
  `Environment`, `LocalDesign`, `Metric`, `Value`, and `ValueType`.
  Empty when no efficiency outputs are detected.

- `summary`:

  High-level summary list for the overall MET design.

- `seed_used`:

  Seed used internally.

## Details

Compared with calling the lower-level functions manually,
`plan_sparse_met_design()` ensures that:

- every environment has a design specification

- allocation metadata are propagated into the combined field book

- local design metadata are harmonized across heterogeneous design
  engines

- seed-aware replication can be applied environment by environment

- local design efficiency outputs, when present as `design$efficiency`,
  are automatically captured into `environment_summary` and
  `efficiency_summary`

- environment specs are validated before design construction

- [`prep_famoptg()`](https://FAkohoue.github.io/OptiSparseMET/reference/prep_famoptg.md)
  replication mode can be inferred automatically when not supplied
  explicitly

- a registry-based compatibility layer allows additional future design
  engines to be added with minimal changes
