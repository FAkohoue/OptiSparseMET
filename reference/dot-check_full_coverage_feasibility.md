# Enforce full-coverage feasibility for sparse MET allocation

Internal function called by
[`allocate_sparse_met()`](https://FAkohoue.github.io/OptiSparseMET/reference/allocate_sparse_met.md)
to verify that the supplied `n_test_entries_per_environment` provides
enough sparse slots to assign every non-common treatment at least once.
Stops with an informative error when the condition is not met.

## Usage

``` r
.check_full_coverage_feasibility(
  treatments,
  environments,
  n_test_entries_per_environment,
  common_treatments = character(0)
)
```

## Arguments

- treatments:

  Character vector of all treatment IDs.

- environments:

  Character vector of environment names.

- n_test_entries_per_environment:

  Integer scalar or integer vector.

- common_treatments:

  Character vector of common treatment IDs.

## Value

Invisibly returns `TRUE` when the condition is satisfied.

## Details

Unlike
[`warn_if_k_too_small()`](https://FAkohoue.github.io/OptiSparseMET/reference/warn_if_k_too_small.md),
this function stops rather than warns and is not intended for direct
use.
