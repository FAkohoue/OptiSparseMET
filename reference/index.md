# Package index

## Across-environment allocation

Functions for distributing treatments across environments and for
verifying the feasibility and balance of the resulting incidence
structure.

- [`allocate_sparse_met()`](https://FAkohoue.github.io/OptiSparseMET/reference/allocate_sparse_met.md)
  : Allocate test treatments across environments for sparse MET designs
- [`check_balanced_incomplete_feasibility()`](https://FAkohoue.github.io/OptiSparseMET/reference/check_balanced_incomplete_feasibility.md)
  : Evaluate feasibility of an exact balanced incomplete sparse MET
  allocation
- [`derive_allocation_groups()`](https://FAkohoue.github.io/OptiSparseMET/reference/derive_allocation_groups.md)
  : Derive allocation group labels for sparse MET treatment assignment

## Feasibility and capacity helpers

Pre-flight diagnostics to verify that the chosen per-environment
capacity is sufficient to assign every non-common treatment at least
once before allocation begins.

- [`feasibility_helpers`](https://FAkohoue.github.io/OptiSparseMET/reference/feasibility_helpers.md)
  : Feasibility helpers for sparse MET allocation
- [`min_k_for_full_coverage()`](https://FAkohoue.github.io/OptiSparseMET/reference/min_k_for_full_coverage.md)
  : Compute the minimum per-environment capacity for full treatment
  coverage
- [`suggest_safe_k()`](https://FAkohoue.github.io/OptiSparseMET/reference/suggest_safe_k.md)
  : Suggest a safe uniform per-environment capacity for sparse MET
  allocation
- [`warn_if_k_too_small()`](https://FAkohoue.github.io/OptiSparseMET/reference/warn_if_k_too_small.md)
  : Warn when per-environment capacity is insufficient for full
  treatment coverage

## Seed-aware replication planning

Partition treatments into replicated, unreplicated, and excluded roles
based on available seed quantities and per-plot seed requirements.

- [`assign_replication_by_seed()`](https://FAkohoue.github.io/OptiSparseMET/reference/assign_replication_by_seed.md)
  : Classify treatments into replication roles based on seed
  availability

## Within-environment field design

Construct field layouts within each environment. Two engines are
available: block-based repeated-check designs and alpha row-column
stream designs.

- [`met_prep_famoptg()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_prep_famoptg.md)
  : Construct a repeated-check block design with flexible replication
- [`met_alpha_rc_stream()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_alpha_rc_stream.md)
  : Construct a stream-based repeated-check alpha row-column design

## Efficiency evaluation

Evaluate the statistical efficiency of within-environment designs under
A-optimality, D-optimality, and CDmean criteria before field deployment.

- [`met_evaluate_famoptg_efficiency()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_evaluate_famoptg_efficiency.md)
  : Evaluate the statistical efficiency of a repeated-check block design
- [`met_evaluate_alpha_efficiency()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_evaluate_alpha_efficiency.md)
  : Evaluate the statistical efficiency of an alpha-lattice design

## Design optimisation

Search for higher-efficiency field arrangements using criterion-driven
optimisation. Supports Random Restart, Simulated Annealing, and Genetic
Algorithm methods.

- [`met_optimize_famoptg()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_optimize_famoptg.md)
  : Search for a criterion-optimal repeated-check block design
- [`met_optimize_alpha_rc()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_optimize_alpha_rc.md)
  : Search for a criterion-optimal alpha-lattice design

## Pipeline and assembly

End-to-end pipeline orchestration and assembly of environment-level
field books into a single MET-level field book.

- [`plan_sparse_met_design()`](https://FAkohoue.github.io/OptiSparseMET/reference/plan_sparse_met_design.md)
  : Plan a sparse multi-environment trial design and assemble a combined
  field book
- [`combine_met_fieldbooks()`](https://FAkohoue.github.io/OptiSparseMET/reference/combine_met_fieldbooks.md)
  : Combine environment-level field books into a single MET field book

## Example data

Bundled example dataset for illustrating and testing the full
OptiSparseMET workflow.

- [`OptiSparseMET_example_data`](https://FAkohoue.github.io/OptiSparseMET/reference/OptiSparseMET_example_data.md)
  : Example data for OptiSparseMET
