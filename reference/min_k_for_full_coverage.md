# Compute the minimum per-environment capacity for full treatment coverage

Let \\J\\ be the total number of treatments, \\C\\ the number of common
treatments assigned to every environment, and \\I\\ the number of
environments. The number of sparse (non-common) treatments is:

\$\$J^\* = J - C\$\$

Common treatments consume \\C\\ slots in every environment before sparse
allocation begins. The number of sparse-allocatable slots per
environment is:

\$\$k_e^\* = k_e - C\$\$

For every sparse treatment to appear in at least one environment, the
total number of sparse slots across all environments must satisfy:

\$\$\sum\_{e=1}^{I} k_e^\* \geq J^\*\$\$

Under equal per-environment capacity, the minimum number of sparse slots
per environment is:

\$\$k^\*\_{\min} = \left\lceil J^\* / I \right\rceil\$\$

and the corresponding minimum total entries per environment is:

\$\$k\_{\min} = k^\*\_{\min} + C\$\$

This function returns both quantities. The value
`min_total_entries_per_environment` is what should be passed to
`n_test_entries_per_environment` to guarantee full coverage at minimum
capacity.

## Usage

``` r
min_k_for_full_coverage(
  n_treatments_total,
  n_environments,
  n_common_treatments = 0
)
```

## Arguments

- n_treatments_total:

  Positive integer. Total number of treatments including common
  treatments. Corresponds to `length(treatments)` in
  [`allocate_sparse_met()`](https://FAkohoue.github.io/OptiSparseMET/reference/allocate_sparse_met.md).

- n_environments:

  Positive integer. Number of environments in the MET. Must be at least
  1.

- n_common_treatments:

  Non-negative integer, default `0`. Number of treatments assigned to
  every environment before sparse allocation. Must not exceed
  `n_treatments_total`.

## Value

A named list with three components:

- `n_sparse_treatments`:

  Integer. Number of non-common treatments \\J^\* = J - C\\. These are
  the treatments subject to sparse allocation.

- `min_sparse_slots_per_environment`:

  Integer. Minimum number of sparse-allocatable slots required per
  environment under equal capacities: \\\lceil J^\* / I \rceil\\. This
  is the minimum value of \\k_e^\* = k_e - C\\.

- `min_total_entries_per_environment`:

  Integer. Minimum total entries per environment including common
  treatments: \\\lceil J^\* / I \rceil + C\\. Pass this as
  `n_test_entries_per_environment` to guarantee that every non-common
  treatment can be assigned at least once.

## Details

Returns the minimum value of `n_test_entries_per_environment` under a
uniform capacity assumption that guarantees every non-common treatment
can be assigned to at least one environment. This is the floor below
which any sparse allocation must leave some treatments unassigned.

## Examples

``` r
## No common treatments: minimum k is ceil(80 / 3) = 27
min_k_for_full_coverage(
  n_treatments_total  = 80,
  n_environments      = 3,
  n_common_treatments = 0
)
#> $n_sparse_treatments
#> [1] 80
#> 
#> $min_sparse_slots_per_environment
#> [1] 27
#> 
#> $min_total_entries_per_environment
#> [1] 27
#> 

## With 5 common treatments: J* = 75, minimum sparse k = ceil(75 / 3) = 25,
## minimum total k = 25 + 5 = 30
min_k_for_full_coverage(
  n_treatments_total  = 80,
  n_environments      = 3,
  n_common_treatments = 5
)
#> $n_sparse_treatments
#> [1] 75
#> 
#> $min_sparse_slots_per_environment
#> [1] 25
#> 
#> $min_total_entries_per_environment
#> [1] 30
#> 
```
