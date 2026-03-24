# Evaluate feasibility of an exact balanced incomplete sparse MET allocation

For a balanced incomplete allocation to be exact, the slot identity:

\$\$J^\* \times r = \sum\_{e=1}^{I} k_e^\*\$\$

must hold, where \\J^\*\\ is the number of non-common treatments, \\r\\
is the target number of environments per treatment, \\I\\ is the number
of environments, and \\k_e^\*\\ is the number of sparse-allocatable
slots in environment \\e\\ after common treatments have been subtracted
from its total capacity. When environments have equal capacity \\k^\*\\,
this reduces to the standard BIBD identity \\J^\* \times r = I \times
k^\*\\.

The function evaluates this condition and returns both a logical
indicator and the signed difference \\\sum k_e^\* - J^\* \times r\\,
which quantifies the degree of imbalance when exact feasibility fails. A
positive difference means there are more available slots than required —
some treatments will receive an extra replication in an approximate
allocation. A negative difference means slots are insufficient — some
treatments will be under-replicated. Both cases can be handled by
[`allocate_sparse_met()`](https://FAkohoue.github.io/OptiSparseMET/reference/allocate_sparse_met.md)
with `allow_approximate = TRUE`, but the magnitude of the difference
informs whether the resulting imbalance is practically acceptable.

## Usage

``` r
check_balanced_incomplete_feasibility(
  n_treatments_total,
  n_environments,
  n_test_entries_per_environment,
  target_replications,
  n_common_treatments = 0L
)
```

## Arguments

- n_treatments_total:

  Positive integer. Total number of test treatments, including any
  common treatments. This is the full candidate pool before any
  subdivision into common and sparse subsets.

- n_environments:

  Integer greater than or equal to 2. Number of environments in the
  trial.

- n_test_entries_per_environment:

  Integer scalar or integer vector. Total number of test treatments
  assigned to each environment, including common treatments. If a
  scalar, the same capacity is applied to all environments. If a vector,
  its length must equal `n_environments`. All values must be positive
  and must not be less than `n_common_treatments`.

- target_replications:

  Positive integer. Desired number of environments in which each
  non-common treatment should appear. This is the \\r\\ term in the slot
  identity. For exact feasibility, the product \\J^\* \times r\\ must
  equal the total number of sparse slots \\\sum k_e^\*\\.

- n_common_treatments:

  Non-negative integer, default `0`. Number of treatments forced into
  all environments before sparse allocation. These treatments consume
  capacity in every environment and are excluded from the sparse
  allocation pool. Setting this to `0` corresponds to a design with no
  forced common entries.

## Value

A named list with the following components:

- `feasible`:

  Logical. `TRUE` if and only if the slot identity holds exactly, i.e.
  `difference == 0`.

- `n_sparse_treatments`:

  Integer. Number of non-common treatments \\J^\* = J - C\\. These are
  the treatments subject to sparse allocation.

- `k_sparse`:

  Integer vector of length `n_environments`. Sparse- allocatable slots
  per environment \\k_e^\* = k_e - C\\.

- `total_sparse_slots`:

  Integer. Total available slots for sparse allocation:
  \\\sum\_{e=1}^{I} k_e^\*\\.

- `required_sparse_slots`:

  Integer. Slots required for exact balance: \\J^\* \times r\\.

- `difference`:

  Integer. Signed difference \\\sum k_e^\* - J^\* \times r\\. Zero when
  allocation is exactly feasible. Positive when slots exceed
  requirement; negative when slots are insufficient.

- `message`:

  Character scalar. Human-readable summary of the feasibility check,
  including the values of the key quantities and the direction of any
  imbalance.

## Details

`check_balanced_incomplete_feasibility()` determines whether the
parameters of a proposed sparse MET design admit an exact balanced
incomplete allocation — that is, whether the total number of
sparse-allocatable treatment slots across environments exactly equals
the number of non-common treatments multiplied by the target
replication. The function is a diagnostic companion to
[`allocate_sparse_met()`](https://FAkohoue.github.io/OptiSparseMET/reference/allocate_sparse_met.md)
with `allocation_method = "balanced_incomplete"`, and should be called
before allocation when the user wants to verify feasibility or
understand the magnitude of any imbalance before deciding whether to set
`allow_approximate = TRUE`.

Common treatments are assigned to every environment before sparse
allocation begins. Their count is subtracted from each environment's
total capacity to obtain the sparse-allocatable slots:

\$\$k_e^\* = k_e - C\$\$

where \\k_e\\ is the total number of test entries assigned to
environment \\e\\ (`n_test_entries_per_environment`) and \\C\\ is the
number of common treatments (`n_common_treatments`). The number of
non-common treatments is:

\$\$J^\* = J - C\$\$

where \\J\\ is `n_treatments_total`. If any environment has \\k_e \<
C\\, that environment cannot accommodate all common treatments and the
function stops with an error before evaluating the balance condition.

This function performs no allocation. It only evaluates the arithmetic
condition and returns the components needed to diagnose feasibility. For
the actual construction of the incidence matrix, see
[`allocate_sparse_met()`](https://FAkohoue.github.io/OptiSparseMET/reference/allocate_sparse_met.md).

## References

Montesinos-López, O. A., Mosqueda-González, B. A., Salinas-Ruiz, J.,
Montesinos-López, A., & Crossa, J. (2023). Sparse multi-trait genomic
prediction under balanced incomplete block design. *The Plant Genome*,
16, e20305.

## Examples

``` r
## Example 1: exact feasibility — 110 sparse treatments x 2 reps = 220,
## 4 environments x 55 sparse slots each = 220.
check_balanced_incomplete_feasibility(
  n_treatments_total             = 120,
  n_environments                 = 4,
  n_test_entries_per_environment = 65,
  target_replications            = 2,
  n_common_treatments            = 10
)
#> $feasible
#> [1] TRUE
#> 
#> $n_sparse_treatments
#> [1] 110
#> 
#> $k_sparse
#> [1] 55 55 55 55
#> 
#> $total_sparse_slots
#> [1] 220
#> 
#> $required_sparse_slots
#> [1] 220
#> 
#> $difference
#> [1] 0
#> 
#> $message
#> [1] "Exact balanced incomplete allocation is feasible: 110 sparse treatments x 2 replications = 220 sparse slots, matching the available total."
#> 
# feasible = TRUE, difference = 0

## Example 2: slot deficit — 110 sparse treatments x 2 reps = 220,
## 4 environments x 50 sparse slots each = 200. The deficit of 20 means
## approximately 20 treatments will receive only 1 replication under an
## approximate allocation.
check_balanced_incomplete_feasibility(
  n_treatments_total             = 120,
  n_environments                 = 4,
  n_test_entries_per_environment = 60,
  target_replications            = 2,
  n_common_treatments            = 10
)
#> $feasible
#> [1] FALSE
#> 
#> $n_sparse_treatments
#> [1] 110
#> 
#> $k_sparse
#> [1] 50 50 50 50
#> 
#> $total_sparse_slots
#> [1] 200
#> 
#> $required_sparse_slots
#> [1] 220
#> 
#> $difference
#> [1] -20
#> 
#> $message
#> [1] "Exact balanced incomplete allocation is not feasible: available sparse slots = 200, required sparse slots = 220, difference = -20."
#> 
# feasible = FALSE, difference = -20

## Example 3: slot surplus — 110 sparse treatments x 2 reps = 220,
## 4 environments x 60 sparse slots each = 240. The surplus of 20 means
## approximately 20 treatments will receive a third replication under an
## approximate allocation.
check_balanced_incomplete_feasibility(
  n_treatments_total             = 120,
  n_environments                 = 4,
  n_test_entries_per_environment = 70,
  target_replications            = 2,
  n_common_treatments            = 10
)
#> $feasible
#> [1] FALSE
#> 
#> $n_sparse_treatments
#> [1] 110
#> 
#> $k_sparse
#> [1] 60 60 60 60
#> 
#> $total_sparse_slots
#> [1] 240
#> 
#> $required_sparse_slots
#> [1] 220
#> 
#> $difference
#> [1] 20
#> 
#> $message
#> [1] "Exact balanced incomplete allocation is not feasible: available sparse slots = 240, required sparse slots = 220, difference = 20."
#> 
# feasible = FALSE, difference = 20

## Example 4: heterogeneous environment capacities — useful when
## environments differ in the number of plots available.
check_balanced_incomplete_feasibility(
  n_treatments_total             = 100,
  n_environments                 = 4,
  n_test_entries_per_environment = c(40, 45, 40, 45),
  target_replications            = 3,
  n_common_treatments            = 5
)
#> $feasible
#> [1] FALSE
#> 
#> $n_sparse_treatments
#> [1] 95
#> 
#> $k_sparse
#> [1] 35 40 35 40
#> 
#> $total_sparse_slots
#> [1] 150
#> 
#> $required_sparse_slots
#> [1] 285
#> 
#> $difference
#> [1] -135
#> 
#> $message
#> [1] "Exact balanced incomplete allocation is not feasible: available sparse slots = 150, required sparse slots = 285, difference = -135."
#> 
```
