# Warn when per-environment capacity is insufficient for full treatment coverage

The feasibility condition checked is:

\$\$\sum\_{e=1}^{I} (k_e - C) \geq J^\*\$\$

where \\k_e\\ is the capacity of environment \\e\\, \\C\\ is the number
of common treatments, and \\J^\*\\ is the number of non-common
treatments. When the condition fails, the warning message includes the
current total sparse slots, the number of sparse treatments, and the
minimum uniform capacity that would satisfy the condition.

Use this function interactively before calling
[`allocate_sparse_met()`](https://FAkohoue.github.io/OptiSparseMET/reference/allocate_sparse_met.md)
when you want a non-fatal check, or inside wrapper functions that should
guide the user toward a valid configuration without stopping.

## Usage

``` r
warn_if_k_too_small(
  treatments,
  environments,
  n_test_entries_per_environment,
  common_treatments = NULL
)
```

## Arguments

- treatments:

  Character vector of treatment IDs. Duplicates are silently removed.

- environments:

  Character vector of environment names. Duplicates are silently
  removed.

- n_test_entries_per_environment:

  Integer scalar or integer vector. Total number of entries per
  environment including common treatments. If a scalar, the same
  capacity is applied to all environments. If a vector, its length must
  equal the number of environments. All values must be at least
  `length(common_treatments)`.

- common_treatments:

  Optional character vector of common treatment IDs. Values not present
  in `treatments` are silently dropped. If `NULL`, no common treatments
  are assumed.

## Value

Invisibly returns `NULL`. Called for its side effect (warning) when the
capacity is insufficient.

## Details

Checks whether the supplied `n_test_entries_per_environment` provides
enough sparse slots to assign every non-common treatment at least once.
If not, emits a warning that states the deficit and the minimum capacity
needed to resolve it. Execution continues regardless — this function is
a diagnostic aid, not a hard stop.

## Examples

``` r
trt <- paste0("L", sprintf("%03d", 1:80))
env <- c("Env1", "Env2", "Env3")

## k = 20: total sparse slots = 3 x 20 = 60 < 80 — warns and suggests k = 27
warn_if_k_too_small(
  treatments                     = trt,
  environments                   = env,
  n_test_entries_per_environment = 20
)
#> Warning: Current design is infeasible for full treatment coverage: 60 sparse slots available for 80 non-common treatments. For a uniform design, use at least 27 entries per environment.

## k = 27: total sparse slots = 3 x 27 = 81 >= 80 — no warning
warn_if_k_too_small(
  treatments                     = trt,
  environments                   = env,
  n_test_entries_per_environment = 27
)

## Heterogeneous capacities: total = 20 + 25 + 20 = 65 < 80 — warns
warn_if_k_too_small(
  treatments                     = trt,
  environments                   = env,
  n_test_entries_per_environment = c(20, 25, 20)
)
#> Warning: Current design is infeasible for full treatment coverage: 65 sparse slots available for 80 non-common treatments. For a uniform design, use at least 27 entries per environment.
```
