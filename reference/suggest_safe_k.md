# Suggest a safe uniform per-environment capacity for sparse MET allocation

The suggested value is:

\$\$k\_{\text{safe}} = \left\lceil J^\* / I \right\rceil + C + b\$\$

where \\J^\* = J - C\\ is the number of non-common treatments, \\I\\ is
the number of environments, \\C\\ is the number of common treatments,
and \\b\\ is `buffer`. Setting `buffer = 0` returns the strict minimum
from
[`min_k_for_full_coverage()`](https://FAkohoue.github.io/OptiSparseMET/reference/min_k_for_full_coverage.md).
A buffer of 3 to 5 is typically sufficient to absorb rounding and ensure
[`allocate_sparse_met()`](https://FAkohoue.github.io/OptiSparseMET/reference/allocate_sparse_met.md)
runs without needing `allow_approximate = TRUE`.

## Usage

``` r
suggest_safe_k(treatments, environments, common_treatments = NULL, buffer = 3L)
```

## Arguments

- treatments:

  Character vector of treatment IDs. Duplicates are silently removed.

- environments:

  Character vector of environment names. Duplicates are silently
  removed.

- common_treatments:

  Optional character vector of common treatment IDs. Values not present
  in `treatments` are silently dropped. If `NULL`, no common treatments
  are assumed.

- buffer:

  Non-negative integer, default `3`. Extra slots added on top of the
  strict minimum capacity. A small buffer reduces the risk of
  infeasibility from rounding and makes the configuration more robust to
  minor changes in trial dimensions.

## Value

Integer scalar. A safe uniform value for
`n_test_entries_per_environment`.

## Details

Returns a single integer suitable for passing to
`n_test_entries_per_environment` in
[`allocate_sparse_met()`](https://FAkohoue.github.io/OptiSparseMET/reference/allocate_sparse_met.md).
The value is the minimum capacity that guarantees full treatment
coverage plus a user-defined buffer, making it appropriate for examples,
tests, and demonstration workflows where a comfortably feasible setting
is preferred over the exact minimum.

## Examples

``` r
trt <- paste0("L", sprintf("%03d", 1:80))
env <- c("Env1", "Env2", "Env3")

## No common treatments, default buffer of 3:
## ceil(80 / 3) + 0 + 3 = 27 + 3 = 30
suggest_safe_k(
  treatments        = trt,
  environments      = env,
  common_treatments = NULL,
  buffer            = 3
)
#> [1] 30

## With 5 common treatments, buffer of 2:
## ceil(75 / 3) + 5 + 2 = 25 + 5 + 2 = 32
suggest_safe_k(
  treatments        = trt,
  environments      = env,
  common_treatments = trt[1:5],
  buffer            = 2
)
#> [1] 32
```
