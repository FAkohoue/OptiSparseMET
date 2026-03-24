# Feasibility helpers for sparse MET allocation

Utility functions for assessing whether a proposed sparse MET
configuration has sufficient per-environment capacity to assign every
non-common treatment at least once, and for deriving safe values for
`n_test_entries_per_environment` prior to calling
[`allocate_sparse_met()`](https://FAkohoue.github.io/OptiSparseMET/reference/allocate_sparse_met.md).

## Details

These helpers address a specific failure mode: passing
`n_test_entries_per_environment` that is too small to cover all
treatments, which causes
[`allocate_sparse_met()`](https://FAkohoue.github.io/OptiSparseMET/reference/allocate_sparse_met.md)
to either error, produce incomplete allocation, or require
`allow_approximate = TRUE` in ways that leave some treatments
unassigned. Calling these helpers before allocation lets the user verify
feasibility and choose a safe capacity without trial and error.

The three exported helpers serve distinct purposes.
[`min_k_for_full_coverage()`](https://FAkohoue.github.io/OptiSparseMET/reference/min_k_for_full_coverage.md)
computes the arithmetic minimum.
[`suggest_safe_k()`](https://FAkohoue.github.io/OptiSparseMET/reference/suggest_safe_k.md)
adds a user-defined buffer on top of that minimum and returns a single
ready-to-use integer.
[`warn_if_k_too_small()`](https://FAkohoue.github.io/OptiSparseMET/reference/warn_if_k_too_small.md)
checks an existing configuration and emits a warning with a corrective
suggestion if the capacity is insufficient. The internal function
[`.check_full_coverage_feasibility()`](https://FAkohoue.github.io/OptiSparseMET/reference/dot-check_full_coverage_feasibility.md)
is called by
[`allocate_sparse_met()`](https://FAkohoue.github.io/OptiSparseMET/reference/allocate_sparse_met.md)
to enforce the same condition with a hard stop.
