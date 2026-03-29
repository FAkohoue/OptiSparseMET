# Allocate test treatments across environments for sparse MET designs

Two allocation strategies are available. `"random_balanced"` implements
an M3-inspired stochastic allocation that approximates balance without
requiring exact BIBD parameters to be satisfiable – appropriate when
environment capacities differ or when the trial dimensions do not admit
an exact balanced solution. Unlike the original M3 of Montesinos-Lopez
et al. (2023), which allocates location by location and may silently
leave some lines with fewer than \\r\\ replications, this implementation
uses a two-phase coverage-first strategy that additionally guarantees
every non-common treatment appears in at least one environment before
replication filling begins.

The `"balanced_incomplete"` method implements the M4 allocation of
Montesinos-Lopez et al. (2023): every non-common treatment appears in
exactly \\r\\ environments (equal replication) and every environment
receives exactly \\k^\*\\ sparse treatments (equal environment sizes),
so the resource identity \\J^\* \times r = I \times k^\*\\ holds
exactly. This equal-replication, equal-environment-size guarantee is
what distinguishes M4 from M3 in the paper, and is always the goal in
plant breeding programs where thousands of lines are tested across a few
environments.

`allow_approximate = FALSE` (the default) is the standard M4 path: the
slot identity must hold exactly, and the function stops with an
informative error if it cannot be met, so the caller always knows
whether equal replication was achieved. Construction first tries
[`crossdes::find.BIB()`](https://rdrr.io/pkg/crossdes/man/find.BIB.html)
(if crossdes is installed), then falls back to a greedy load-balanced
constructor. `allow_approximate = TRUE` relaxes the slot identity and
allows minor replication imbalances across lines; it is a fallback for
exploratory use, not the intended primary path.

Before allocation begins, the function calls
[`.check_full_coverage_feasibility()`](https://FAkohoue.github.io/OptiSparseMET/reference/dot-check_full_coverage_feasibility.md)
to verify that the total number of sparse slots across all environments
is sufficient to assign every non-common treatment at least once. If it
is not, the function stops with an informative error naming the deficit
and the minimum capacity that would resolve it. Use
[`suggest_safe_k()`](https://FAkohoue.github.io/OptiSparseMET/reference/suggest_safe_k.md)
or
[`min_k_for_full_coverage()`](https://FAkohoue.github.io/OptiSparseMET/reference/min_k_for_full_coverage.md)
to choose a feasible `n_test_entries_per_environment` before calling
this function.

## Usage

``` r
allocate_sparse_met(
  treatments,
  environments,
  allocation_method = c("random_balanced", "balanced_incomplete", "M3", "M4"),
  n_test_entries_per_environment,
  target_replications = NULL,
  common_treatments = NULL,
  allocation_group_source = c("none", "Family", "GRM", "A"),
  treatment_info = NULL,
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

  Character vector of test treatment IDs to allocate across
  environments. Check treatments should not be included here; they are
  handled separately within the within-environment design functions.
  Duplicate values are silently removed.

- environments:

  Character vector of environment names. Must contain at least two
  elements. Duplicate values are silently removed.

- allocation_method:

  Character scalar. Sparse allocation strategy. Accepted values are
  `"random_balanced"` (M3-inspired stochastic allocation) and
  `"balanced_incomplete"` (M4-type BIBD allocation). The aliases `"M3"`
  and `"M4"` are also accepted and translated internally to their
  canonical names before any further processing.

- n_test_entries_per_environment:

  Integer scalar or integer vector. Total number of test treatments per
  environment including common treatments. If a scalar, it applies
  uniformly to all environments. If a vector, its length must equal the
  number of environments. All values must be positive and must be at
  least `length(common_treatments)`. Use
  [`suggest_safe_k()`](https://FAkohoue.github.io/OptiSparseMET/reference/suggest_safe_k.md)
  to choose a value that guarantees full treatment coverage.

- target_replications:

  Optional positive integer. Target number of environments in which each
  non-common treatment should appear. For `"random_balanced"`, this is a
  soft target that guides the phase-two filling. For
  `"balanced_incomplete"`, this is the strict replication level required
  for an exact balanced solution. If `NULL`, the function infers a value
  as `floor(total_sparse_slots / n_sparse_treatments)`, with a minimum
  of 1.

- common_treatments:

  Optional character vector of treatment IDs to assign to all
  environments. Placed before sparse allocation begins and excluded from
  the sparse pool. Values not present in `treatments` are silently
  dropped.

- allocation_group_source:

  Character scalar. Controls whether and how genetic group structure
  guides allocation. `"none"` disables group-guided allocation.
  `"Family"` uses `treatment_info$Family`. `"GRM"` derives clusters from
  `GRM`. `"A"` derives clusters from `A`. Active in both allocation
  phases when set to anything other than `"none"`.

- treatment_info:

  Optional data frame. Required when
  `allocation_group_source = "Family"`. Must contain columns `Treatment`
  and `Family`. When `allocation_group_source %in% c("GRM", "A")`, this
  argument is optional but used to anchor the number of clusters when
  supplied.

- GRM:

  Optional numeric matrix. Genomic relationship matrix. Required when
  `allocation_group_source = "GRM"`. Row and column names must match
  treatment IDs or be reachable through `id_map`.

- A:

  Optional numeric matrix. Pedigree-based numerator relationship matrix.
  Required when `allocation_group_source = "A"`. Same naming
  requirements as `GRM`.

- id_map:

  Optional data frame with columns `Treatment` and `LineID`. Required
  only when `allocation_group_source %in% c("GRM", "A")` and treatment
  IDs do not match the row/column names of the relationship matrix.

- group_method:

  Character scalar. Clustering algorithm applied to the principal
  components of `GRM` or `A`. `"kmeans"` or `"hclust"`. Ignored when
  `allocation_group_source %in% c("none", "Family")`.

- group_seed:

  Integer. Seed for k-means initialization. Active only when
  `allocation_group_source %in% c("GRM", "A")` and
  `group_method = "kmeans"`.

- group_attempts:

  Integer. Number of random restarts for k-means. Active only when
  `allocation_group_source %in% c("GRM", "A")` and
  `group_method = "kmeans"`.

- n_pcs_use:

  Integer or `Inf`. Number of leading principal components retained for
  matrix-based clustering. `Inf` retains all components corresponding to
  positive eigenvalues. Ignored when
  `allocation_group_source %in% c("none", "Family")`.

- min_groups_per_environment:

  Optional positive integer. Minimum number of allocation groups that
  should be represented in each environment, where feasible given
  available treatments. Active in phase two when
  `allocation_group_source` is not `"none"`.

- min_env_per_group:

  Optional positive integer. Minimum number of environments in which
  each allocation group should appear, where feasible. Used in phase two
  when `force_group_connectivity = TRUE`.

- balance_groups_across_env:

  Logical. If `TRUE`, phase two preferentially assigns treatments from
  groups whose current allocation count is below their proportional
  target. Active when `allocation_group_source` is not `"none"`.

- force_group_connectivity:

  Logical. If `TRUE`, phase two preferentially assigns treatments from
  groups that have not yet appeared in `min_env_per_group` environments.
  Active when `allocation_group_source` is not `"none"` and
  `min_env_per_group` is not `NULL`.

- allow_approximate:

  Logical, default `FALSE`. When `FALSE` and
  `allocation_method = "balanced_incomplete"`, the slot identity \\J^\*
  \times r = I \times k^\*\\ must hold exactly; the function stops with
  an error if it does not. This is the standard M4 path and guarantees
  equal replication for every non-common treatment. When `TRUE`, the
  slot identity is not enforced and minor replication imbalances are
  accepted; this is a relaxed fallback for exploratory use, not the
  primary mode. For `"random_balanced"`, this argument has no effect.

- seed:

  Optional integer. Random seed for reproducibility. Controls the random
  order in which sparse treatments are processed in phase one and the
  stochastic tie-breaking in phase two under `"random_balanced"`, as
  well as tie-breaking in the strict exact constructor. If `NULL`, no
  seed is set and results may differ across runs; the seed used
  internally is returned as `seed_used`.

## Value

A named list with the following components:

- `allocation_matrix`:

  Binary integer matrix of dimension `n_treatments x n_environments`
  with `dimnames` set to `treatments` and `environments`. Entry `[i, e]`
  is `1L` if treatment `i` is assigned to environment `e`, and `0L`
  otherwise. Every non-common treatment is guaranteed to have row sum at
  least 1.

- `allocation_long`:

  Long-format data frame with one row per treatment-by-environment
  combination. Columns: `Treatment`, `Environment`, `Assigned` (integer
  0/1), `IsCommonTreatment` (logical), and `AllocationGroup` (character,
  present when `allocation_group_source` is not `"none"`).

- `overlap_matrix`:

  Square integer matrix of dimension `n_environments x n_environments`
  giving the number of treatments shared between each pair of
  environments. Diagonal entries give the total treatment count per
  environment.

- `line_replications`:

  Named integer vector of length `n_treatments`. Row sums of
  `allocation_matrix`: the number of environments each treatment enters.

- `environment_sizes`:

  Named integer vector of length `n_environments`. Column sums of
  `allocation_matrix`: the number of treatments assigned to each
  environment.

- `group_assignment`:

  Data frame with columns `Treatment` and `AllocationGroup`, one row per
  treatment. `NULL` when `allocation_group_source = "none"`.

- `group_by_environment`:

  Data frame summarizing the count of assigned treatments from each
  allocation group in each environment. Columns: `Environment`,
  `AllocationGroup`, `n_treatments`. `NULL` when
  `allocation_group_source = "none"`.

- `group_overlap_matrix`:

  Square integer matrix of dimension `n_environments x n_environments`
  giving the number of shared allocation groups between each pair of
  environments. `NULL` when `allocation_group_source = "none"`.

- `summary`:

  Named list with allocation metadata: `allocation_method`,
  `allocation_group_source`, `target_replications`,
  `n_treatments_total`, `n_sparse_treatments`, `n_common_treatments`,
  `total_sparse_slots`, `environment_sizes`, `min_replication`,
  `max_replication`, `mean_replication`, `min_sparse_replication`,
  `max_sparse_replication`, `mean_sparse_replication`,
  `min_common_replication`, `max_common_replication`, and
  `mean_common_replication`.

- `seed_used`:

  The integer seed passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) internally, or
  `NULL` if no seed was supplied.

## Details

`allocate_sparse_met()` constructs the treatment-by-environment
incidence structure for a sparse multi-environment trial. Given a set of
candidate lines and a set of environments, the function determines which
lines enter which environments subject to replication targets, capacity
constraints, optional common treatment requirements, and optional
genetic structure constraints. The function guarantees that every
non-common treatment is assigned to at least one environment before
attempting to fill additional replication slots.

### Relation to sparse testing theory

The allocation strategies implemented here correspond to M3 and M4 in
Montesinos-Lopez et al. (2023). The underlying resource identity is:

\$\$N = J \times r = I \times k\$\$

where \\J\\ is the number of lines, \\I\\ is the number of environments,
\\k\\ is the number of lines per environment, \\r\\ is the number of
environments each line enters, and \\N\\ is the total number of
line-by-environment observations. This identity makes the tradeoff
between coverage (\\k\\) and replication depth (\\r\\) explicit given
fixed \\N\\.

### Common treatments

Treatments in `common_treatments` are assigned to every environment
before sparse allocation begins and are excluded from the sparse
allocation pool. The per-environment sparse capacity is:

\$\$k_e^\* = k_e - C\$\$

where \\C\\ is the number of common treatments and \\k_e\\ is the total
capacity of environment \\e\\. Common treatments establish direct
cross-environment connectivity that does not depend on the relationship
structure in the covariance model, and provide stable references for
estimating environment-level effects. They are most important when
environments are weakly correlated.

### Two-phase allocation

For `"random_balanced"` and approximate `"balanced_incomplete"`, phase
one iterates over sparse treatments in random order and assigns each to
one environment, preferring environments where the treatment's genetic
group is not yet represented when `allocation_group_source` is active.
Phase two iterates over environments in decreasing order of remaining
capacity and fills each to its target, guided by line-level replication
deficit scores and group-balance penalties.

### Within-environment layout

Once allocation is complete, the set of treatments assigned to each
environment is passed to the within-environment design function.
[`met_prep_famoptg()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_prep_famoptg.md)
constructs repeated-check augmented, partially replicated (p-rep), and
RCBD-type block designs.
[`met_alpha_rc_stream()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_alpha_rc_stream.md)
constructs stream-based alpha row-column designs for fixed-grid field
geometry. Both accept output from
[`assign_replication_by_seed()`](https://FAkohoue.github.io/OptiSparseMET/reference/assign_replication_by_seed.md)
directly. The end-to-end pipeline that coordinates allocation and local
design construction is
[`plan_sparse_met_design()`](https://FAkohoue.github.io/OptiSparseMET/reference/plan_sparse_met_design.md).

### Allocation groups

When `allocation_group_source` is not `"none"`, both phases use genetic
group membership to guide assignments. Family-based grouping reads
labels from `treatment_info$Family`. GRM-based and A-based grouping
derive clusters from the leading principal components of the respective
relationship matrix. In phase two, the allocator penalizes assignments
that concentrate a group in too few environments
(`balance_groups_across_env`) and gives preference to groups that have
not yet reached `min_env_per_group` (`force_group_connectivity`).
Line-level replication targets are preserved within these group-level
constraints.

### M4 allocation and the role of allow_approximate

The M4 method (Montesinos-Lopez et al., 2023) enforces two conditions:

1.  Equal replication: every non-common treatment appears in exactly
    \\r\\ environments.

2.  Equal environment sizes: every environment receives exactly \\k^\*\\
    sparse treatments, so \\J^\* \times r = I \times k^\*\\.

`allow_approximate = FALSE` (the default) enforces both conditions
strictly. If the slot identity \\J^\* \times r = I \times k^\*\\ does
not hold for the chosen `n_test_entries_per_environment` and
`target_replications`, the function stops with an informative error. Use
[`check_balanced_incomplete_feasibility()`](https://FAkohoue.github.io/OptiSparseMET/reference/check_balanced_incomplete_feasibility.md)
to verify the slot identity before calling, or adjust \\k\\ and \\r\\ so
that \\J^\* \times r = I \times k^\*\\. Construction first tries
[`crossdes::find.BIB()`](https://rdrr.io/pkg/crossdes/man/find.BIB.html)
(if crossdes is installed), then falls back to a greedy load-balanced
constructor that guarantees equal replication even when crossdes is
absent.

`allow_approximate = TRUE` relaxes the slot identity: the function
constructs the most balanced allocation it can without stopping on
infeasibility, accepting that some lines may receive more or fewer
replications than \\r\\. This mode is useful for exploratory analysis
but does not provide the equal-replication guarantee that is the
defining property of M4.

## References

Montesinos-Lopez, O. A., Mosqueda-Gonzalez, B. A., Salinas-Ruiz, J.,
Montesinos-Lopez, A., & Crossa, J. (2023). Sparse multi-trait genomic
prediction under balanced incomplete block design. *The Plant Genome*,
16, e20305. [doi:10.1002/tpg2.20305](https://doi.org/10.1002/tpg2.20305)

## See also

[`suggest_safe_k()`](https://FAkohoue.github.io/OptiSparseMET/reference/suggest_safe_k.md)
and
[`min_k_for_full_coverage()`](https://FAkohoue.github.io/OptiSparseMET/reference/min_k_for_full_coverage.md)
for choosing a feasible `n_test_entries_per_environment` before calling
this function.
[`check_balanced_incomplete_feasibility()`](https://FAkohoue.github.io/OptiSparseMET/reference/check_balanced_incomplete_feasibility.md)
for verifying the slot identity before attempting a
`"balanced_incomplete"` allocation.
[`derive_allocation_groups()`](https://FAkohoue.github.io/OptiSparseMET/reference/derive_allocation_groups.md)
for inspecting the group structure that guides allocation when
`allocation_group_source` is not `"none"`.
[`met_prep_famoptg()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_prep_famoptg.md)
and
[`met_alpha_rc_stream()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_alpha_rc_stream.md)
for the within-environment design functions that consume the allocation
output.
[`plan_sparse_met_design()`](https://FAkohoue.github.io/OptiSparseMET/reference/plan_sparse_met_design.md)
for the end-to-end two-stage MET pipeline.

## Examples

``` r
treatments <- paste0("L", sprintf("%03d", 1:120))
envs       <- c("Env1", "Env2", "Env3", "Env4")
fam        <- rep(paste0("F", 1:6), each = 20)

treatment_info <- data.frame(
  Treatment = treatments,
  Family    = fam,
  stringsAsFactors = FALSE
)

## Check a safe per-environment capacity before running allocation
suggest_safe_k(treatments, envs, buffer = 3)  # 33
#> [1] 33

## Example 1: random balanced allocation with family-guided group spreading
out1 <- allocate_sparse_met(
  treatments                     = treatments,
  environments                   = envs,
  allocation_method              = "random_balanced",
  n_test_entries_per_environment = 33,
  target_replications            = 1,
  allocation_group_source        = "Family",
  treatment_info                 = treatment_info,
  min_groups_per_environment     = 4,
  min_env_per_group              = 2,
  seed                           = 123
)

out1$summary
#> $allocation_method
#> [1] "random_balanced"
#> 
#> $allocation_group_source
#> [1] "Family"
#> 
#> $target_replications
#> [1] 1
#> 
#> $n_treatments_total
#> [1] 120
#> 
#> $n_sparse_treatments
#> [1] 120
#> 
#> $n_common_treatments
#> [1] 0
#> 
#> $total_sparse_slots
#> [1] 132
#> 
#> $environment_sizes
#> Env1 Env2 Env3 Env4 
#>   33   33   33   33 
#> 
#> $min_replication
#> [1] 1
#> 
#> $max_replication
#> [1] 2
#> 
#> $mean_replication
#> [1] 1.1
#> 
#> $min_sparse_replication
#> [1] 1
#> 
#> $max_sparse_replication
#> [1] 2
#> 
#> $mean_sparse_replication
#> [1] 1.1
#> 
#> $min_common_replication
#> [1] NA
#> 
#> $max_common_replication
#> [1] NA
#> 
#> $mean_common_replication
#> [1] NA
#> 
#> $n_groups
#> [1] 6
#> 
head(out1$allocation_long)
#>   Treatment Environment Assigned IsCommonTreatment AllocationGroup
#> 1      L001        Env1        1             FALSE              F1
#> 2      L002        Env1        0             FALSE              F1
#> 3      L003        Env1        0             FALSE              F1
#> 4      L004        Env1        1             FALSE              F1
#> 5      L005        Env1        0             FALSE              F1
#> 6      L006        Env1        1             FALSE              F1
head(out1$group_by_environment)
#>   Environment AllocationGroup n_treatments
#> 1        Env1              F1            7
#> 2        Env2              F1            4
#> 3        Env3              F1            7
#> 4        Env4              F1            4
#> 5        Env1              F2            6
#> 6        Env2              F2            6

## Example 2: M4 balanced incomplete allocation (paper method).
## Equal replication (r=2) and equal environment sizes (k*=55 per env).
## Slot identity: J* x r = I x k* => 110 x 2 = 4 x 55 = 220. Valid.
## allow_approximate = FALSE enforces the slot identity strictly, stopping
## with an error if it is not met (the default safe behaviour).
out2 <- allocate_sparse_met(
  treatments                     = treatments,
  environments                   = envs,
  allocation_method              = "balanced_incomplete",
  n_test_entries_per_environment = 65,
  target_replications            = 2,
  common_treatments              = treatments[1:10],
  allow_approximate              = FALSE,
  seed                           = 123
)

out2$summary
#> $allocation_method
#> [1] "balanced_incomplete"
#> 
#> $allocation_group_source
#> [1] "none"
#> 
#> $target_replications
#> [1] 2
#> 
#> $n_treatments_total
#> [1] 120
#> 
#> $n_sparse_treatments
#> [1] 110
#> 
#> $n_common_treatments
#> [1] 10
#> 
#> $total_sparse_slots
#> [1] 220
#> 
#> $environment_sizes
#> Env1 Env2 Env3 Env4 
#>   65   65   65   65 
#> 
#> $min_replication
#> [1] 2
#> 
#> $max_replication
#> [1] 4
#> 
#> $mean_replication
#> [1] 2.166667
#> 
#> $min_sparse_replication
#> [1] 2
#> 
#> $max_sparse_replication
#> [1] 2
#> 
#> $mean_sparse_replication
#> [1] 2
#> 
#> $min_common_replication
#> [1] 4
#> 
#> $max_common_replication
#> [1] 4
#> 
#> $mean_common_replication
#> [1] 4
#> 
#> $n_groups
#> [1] 0
#> 
# Every sparse line appears in exactly 2 environments
range(out2$line_replications[!(names(out2$line_replications) %in%
                                treatments[1:10])])
#> [1] 2 2
```
