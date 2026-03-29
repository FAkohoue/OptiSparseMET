# OptiSparseMET

![OptiSparseMET logo](reference/figures/logo.png)

------------------------------------------------------------------------

## Overview

`OptiSparseMET` is an R framework for constructing sparse
multi-environment trials (METs) that jointly addresses treatment
allocation across environments and field design within environments –
under shared statistical, genetic, and logistical constraints.

The package targets a structural challenge common to modern breeding
programs: the number of candidate lines routinely exceeds what any
single environment can accommodate, environments are heterogeneous
rather than interchangeable, and seed availability imposes limits that
purely theoretical designs ignore. Standard MET tools typically address
one of these problems at a time. `OptiSparseMET` integrates all of them
within a single workflow.

------------------------------------------------------------------------

## Conceptual Framework

Sparse MET design is a two-level problem and should be treated as such.

**Level 1 — Across-environment allocation** determines which treatments
appear in which environments, how many times each treatment is
replicated across the trial, and whether the resulting incidence
structure preserves sufficient genetic connectivity for valid
cross-environment inference.

**Level 2 — Within-environment design** determines blocking structure,
spatial layout, replication within each environment, and local control
of field heterogeneity.

These two levels are not merely sequential steps — they are
statistically coupled, and optimizing them independently produces
inferior designs. The linkage operates through four mechanisms:

**1. The incidence matrix couples both levels inside the information
matrix.** In the linear mixed model $y = X\beta + Zg + e$, the precision
of all genetic value estimates is governed by the coefficient matrix
$C = Z^{\top}V^{- 1}Z - Z^{\top}V^{- 1}X\left( X^{\top}V^{- 1}X \right)^{- 1}X^{\top}V^{- 1}Z$,
where $V = ZKZ^{\top}\sigma_{g}^{2} + R\sigma_{e}^{2}$. The allocation
decision determines the sparsity pattern of $Z$ (which lines appear
where); the within-environment blocking structure determines $R$ (the
residual covariance). Both enter $V$ and therefore $C^{- 1}$. Neither
can be optimized in isolation because they interact inside the inversion
of $V$.

**2. Allocation fixes which within-environment designs are feasible.**
Once allocation assigns $k_{e}$ lines to environment $e$, the
within-environment design must arrange exactly those $k_{e}$ treatments
across the available $n_{\text{rows}} \times n_{\text{cols}}$ field. If
$k_{e}$ is incompatible with the blocking structure — for example, not a
multiple of the target block size, or exceeding field capacity — the
design is infeasible regardless of how statistically ideal the
allocation was. Allocation and field geometry must be co-designed.

**3. Block efficiency propagates into cross-environment inference.** The
precision of a genetic value estimate for line $j$ in environment $e$ is
proportional to $e_{j}\, r_{j}^{(e)}$, where $r_{j}^{(e)}$ is the number
of plots and $e_{j} \in (0,1\rbrack$ is the efficiency factor of the
within-environment design relative to a completely randomized layout. A
poor block design reduces $e_{j}$, inflating the variance of each BLUP.
These inflated variances propagate into cross-environment covariance
estimates, degrading G×E inference and genetic correlation estimation
even when the allocation incidence structure is perfectly balanced.

**4. CDmean — the genomic prediction criterion — depends on both
levels.** The CDmean criterion,
$\text{CDmean} = 1 - \overline{\text{PEV}}/\sigma_{g}^{2}$, where PEV
depends on both $Z$ (allocation) and $R^{- 1}$ (blocking), cannot be
maximized by fixing either level independently. Spreading genetically
diverse lines across environments improves the genomic connectivity
captured in $Z^{\top}R^{- 1}Z$; efficient blocking sharpens $R^{- 1}$.
Both contributions are necessary.

`OptiSparseMET` formalizes the link between the two levels: the
allocation output specifies exactly which lines enter each environment,
and the within-environment design engine receives precisely that set,
ensuring that the incidence structure and the blocking structure are
optimized consistently within the same statistical framework.

------------------------------------------------------------------------

## Statistical Foundations

### Sparse testing identity

The theoretical basis follows Montesinos-Lopez et al. (2023), who
formalize the resource identity underlying balanced sparse designs:

$$N = J \times r = I \times k$$

| Symbol | Meaning                              |
|--------|--------------------------------------|
| $J$    | Total number of treatments           |
| $I$    | Number of environments               |
| $k$    | Number of treatments per environment |
| $r$    | Number of environments per treatment |

Given fixed total resources $N$, this identity makes the tradeoff
between coverage breadth ($k$) and replication depth ($r$) explicit.

### Allocation strategies

Two strategies are available, corresponding to M3 and M4 in
Montesinos-Lopez et al. (2023):

| Strategy    | Argument              | Properties                                                                                                                                                                 |
|-------------|-----------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| M3-inspired | `random_balanced`     | Coverage-first stochastic allocation; guarantees every treatment appears at least once; tolerates unequal environment sizes                                                |
| M4 BIBD     | `balanced_incomplete` | Enforces equal replication and equal environment sizes (slot identity $J^{*} \times r = I \times k^{*}$); pairwise co-occurrence is computed and returned but not enforced |

#### M4 – equal replication and equal environment sizes

The M4 method (Montesinos-Lopez et al. 2023) enforces two structural
guarantees:

1.  **Equal replication** – every non-common treatment appears in
    exactly $r$ environments.
2.  **Equal environment sizes** – every environment receives exactly
    $k^{*}$ sparse treatments, so the resource identity
    $J^{*} \times r = I \times k^{*}$ holds exactly ($J^{*} = J - C$,
    $k^{*} = k - C$, where $C$ is the number of common treatments).

These are the guarantees that distinguish M4 from M3 in the paper. In
plant breeding programs where thousands of lines are tested across a few
environments, equal replication means every candidate is evaluated the
same number of times – a fundamental fairness and precision requirement.

`allow_approximate = FALSE` (the default) enforces both conditions
strictly. If the slot identity $J^{*} \times r = I \times k^{*}$ does
not hold for the chosen dimensions, the function stops with a clear
error before any allocation is attempted. Use
[`check_balanced_incomplete_feasibility()`](https://FAkohoue.github.io/OptiSparseMET/reference/check_balanced_incomplete_feasibility.md)
to verify the slot identity first, or adjust $k$ and $r$ so that the
identity holds. Construction proceeds via
[`crossdes::find.BIB()`](https://rdrr.io/pkg/crossdes/man/find.BIB.html)
when the [crossdes](https://CRAN.R-project.org/package=crossdes) package
is installed, falling back to a greedy load-balanced constructor
otherwise.

`allow_approximate = TRUE` relaxes the slot identity and allows minor
replication imbalances. This is a fallback for exploratory use, not the
primary mode.

### Genetic connectedness

`OptiSparseMET` addresses genetic disconnectedness through three
mechanisms: **common treatments** forced into every environment
establish model-free cross-environment connectivity; **family-based
allocation** distributes each family group across environments; and
**GRM/A-based allocation** uses genomic or pedigree relationships to
prevent genetic clustering.

### Seed constraints

[`assign_replication_by_seed()`](https://FAkohoue.github.io/OptiSparseMET/reference/assign_replication_by_seed.md)
takes a data frame of available seed quantities and a per-plot seed
requirement and returns a replication plan that respects those
constraints – making designs deployable rather than merely theoretically
optimal.

------------------------------------------------------------------------

## Main Functions

### Across-environment allocation

| Function                                                                                                                                 | Description                                                                                   |
|------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------|
| [`allocate_sparse_met()`](https://FAkohoue.github.io/OptiSparseMET/reference/allocate_sparse_met.md)                                     | Distribute treatments across environments (M3 or M4) with enforced equal replication under M4 |
| [`check_balanced_incomplete_feasibility()`](https://FAkohoue.github.io/OptiSparseMET/reference/check_balanced_incomplete_feasibility.md) | Verify the slot identity $J^{*} \times r = I \times k^{*}$ before attempting M4 allocation    |
| [`derive_allocation_groups()`](https://FAkohoue.github.io/OptiSparseMET/reference/derive_allocation_groups.md)                           | Derive grouping structure from family labels, GRM, or pedigree matrix                         |

### Feasibility and capacity helpers

| Function                                                                                                     | Description                                              |
|--------------------------------------------------------------------------------------------------------------|----------------------------------------------------------|
| [`suggest_safe_k()`](https://FAkohoue.github.io/OptiSparseMET/reference/suggest_safe_k.md)                   | Propose a safe `n_test_entries_per_environment` value    |
| [`min_k_for_full_coverage()`](https://FAkohoue.github.io/OptiSparseMET/reference/min_k_for_full_coverage.md) | Compute the minimum capacity for full treatment coverage |
| [`warn_if_k_too_small()`](https://FAkohoue.github.io/OptiSparseMET/reference/warn_if_k_too_small.md)         | Non-fatal pre-flight capacity check                      |

Call
[`suggest_safe_k()`](https://FAkohoue.github.io/OptiSparseMET/reference/suggest_safe_k.md)
or
[`min_k_for_full_coverage()`](https://FAkohoue.github.io/OptiSparseMET/reference/min_k_for_full_coverage.md)
before
[`allocate_sparse_met()`](https://FAkohoue.github.io/OptiSparseMET/reference/allocate_sparse_met.md)
whenever trial dimensions change.

### Seed-aware replication

| Function                                                                                                           | Description                                                                                       |
|--------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------|
| [`assign_replication_by_seed()`](https://FAkohoue.github.io/OptiSparseMET/reference/assign_replication_by_seed.md) | Partition treatments into replicated, unreplicated, and excluded roles based on seed availability |

### Within-environment field design

| Function                                                                                                                     | Description                                                                                                                                       |
|------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------|
| [`met_prep_famoptg()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_prep_famoptg.md)                               | Augmented, p-rep, and RCBD-type repeated-check block designs                                                                                      |
| [`met_alpha_rc_stream()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_alpha_rc_stream.md)                         | Alpha row-column stream designs for fixed-grid field deployment                                                                                   |
| [`met_evaluate_famoptg_efficiency()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_evaluate_famoptg_efficiency.md) | A/D/CDmean efficiency evaluation for [`met_prep_famoptg()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_prep_famoptg.md) designs       |
| [`met_evaluate_alpha_efficiency()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_evaluate_alpha_efficiency.md)     | A/D/CDmean efficiency evaluation for [`met_alpha_rc_stream()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_alpha_rc_stream.md) designs |
| [`met_optimize_famoptg()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_optimize_famoptg.md)                       | Random Restart optimisation for [`met_prep_famoptg()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_prep_famoptg.md) designs            |
| [`met_optimize_alpha_rc()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_optimize_alpha_rc.md)                     | RS/SA/GA optimisation for [`met_alpha_rc_stream()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_alpha_rc_stream.md) designs            |

### Pipeline and assembly

| Function                                                                                                   | Description                                                 |
|------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------|
| [`plan_sparse_met_design()`](https://FAkohoue.github.io/OptiSparseMET/reference/plan_sparse_met_design.md) | End-to-end two-stage pipeline in a single call              |
| [`combine_met_fieldbooks()`](https://FAkohoue.github.io/OptiSparseMET/reference/combine_met_fieldbooks.md) | Stack environment-level field books into one MET field book |

------------------------------------------------------------------------

## Workflow

    STEP 0  Verify capacity
            suggest_safe_k()  OR  min_k_for_full_coverage()
            |
    STEP 1  Allocate treatments across environments
            allocate_sparse_met()
            |
    STEP 2  Define replication based on seed availability
            assign_replication_by_seed()
            |
    STEP 3  Build within-environment field designs
            met_prep_famoptg()  OR  met_alpha_rc_stream()
            |
    STEP 4  Assemble the combined MET field book
            combine_met_fieldbooks()

Or run the entire pipeline in one call:
[`plan_sparse_met_design()`](https://FAkohoue.github.io/OptiSparseMET/reference/plan_sparse_met_design.md).

------------------------------------------------------------------------

## 5.5 Pipeline inputs: required and optional

Before running any pipeline function, it helps to know exactly what each
function needs. The tables below list every input for the four main
pipeline functions, distinguishing what is strictly required from what
is optional.

### `allocate_sparse_met()` — across-environment allocation

**Required inputs**

| Argument                         | Type             | Description                                                   |
|----------------------------------|------------------|---------------------------------------------------------------|
| `treatments`                     | character vector | All candidate line IDs (J total)                              |
| `environments`                   | character vector | Environment names (≥ 2)                                       |
| `allocation_method`              | character        | `"random_balanced"` (M3) or `"balanced_incomplete"` (M4)      |
| `n_test_entries_per_environment` | integer          | Total entries per environment including common treatments (k) |

**Optional inputs**

| Argument                     | Default  | Description                                                                                           |
|------------------------------|----------|-------------------------------------------------------------------------------------------------------|
| `target_replications`        | inferred | Target environments per sparse line (r); computed from slot identity if NULL                          |
| `common_treatments`          | none     | Lines forced into every environment before sparse allocation                                          |
| `allow_approximate`          | `FALSE`  | `FALSE` = strict equal replication; `TRUE` = relaxed fallback                                         |
| `allocation_group_source`    | `"none"` | Genetic grouping: `"Family"`, `"GRM"`, or `"A"`                                                       |
| `treatment_info`             | NULL     | Data frame with `Treatment` and `Family` columns (required when `allocation_group_source = "Family"`) |
| `GRM`                        | NULL     | Genomic relationship matrix (required when `allocation_group_source = "GRM"`)                         |
| `A`                          | NULL     | Pedigree relationship matrix (required when `allocation_group_source = "A"`)                          |
| `min_groups_per_environment` | NULL     | Minimum genetic groups per environment                                                                |
| `min_env_per_group`          | NULL     | Minimum environments per genetic group                                                                |
| `seed`                       | NULL     | Integer seed for reproducibility                                                                      |

### `assign_replication_by_seed()` — seed-aware replication

**Required inputs**

| Argument                 | Type             | Description                                                    |
|--------------------------|------------------|----------------------------------------------------------------|
| `treatments`             | character vector | All candidate line IDs                                         |
| `seed_available`         | data frame       | Must contain `Treatment` and `SeedAvailable` columns           |
| `seed_required_per_plot` | integer          | Seeds needed per plot (scalar or named vector per environment) |
| `replication_mode`       | character        | `"augmented"`, `"p_rep"`, or `"rcbd_type"`                     |

**Optional inputs**

| Argument               | Default            | Description                                                                    |
|------------------------|--------------------|--------------------------------------------------------------------------------|
| `desired_replications` | 2                  | Target number of plots per replicated line                                     |
| `shortage_action`      | `"downgrade"`      | What to do when seed is insufficient: `"downgrade"`, `"exclude"`, or `"error"` |
| `max_prep`             | NULL               | Maximum number of p-rep treatments (p-rep mode only)                           |
| `priority`             | `"seed_available"` | Selection criterion for p-rep candidates                                       |
| `minimum_seed_buffer`  | 0                  | Extra seeds reserved per line beyond the plot requirement                      |
| `seed`                 | NULL               | Integer seed for reproducibility                                               |

### `met_prep_famoptg()` — block-based field design

**Required inputs**

| Argument           | Type             | Description                                                 |
|--------------------|------------------|-------------------------------------------------------------|
| `check_treatments` | character vector | Check (control) treatment IDs; appear in every block        |
| `check_families`   | character vector | Family labels for checks; same length as `check_treatments` |
| `n_blocks`         | integer          | Number of incomplete blocks                                 |
| `n_rows`           | integer          | Number of field rows                                        |
| `n_cols`           | integer          | Number of field columns                                     |

At least one of `p_rep_treatments` or `unreplicated_treatments` must be
supplied.

**Optional inputs**

| Argument                  | Default   | Description                                                             |
|---------------------------|-----------|-------------------------------------------------------------------------|
| `p_rep_treatments`        | NULL      | Treatments to replicate; typically `rep_plan$p_rep_treatments`          |
| `p_rep_reps`              | NULL      | Replication count per p-rep line; typically `rep_plan$p_rep_reps`       |
| `p_rep_families`          | NULL      | Family labels for p-rep treatments                                      |
| `unreplicated_treatments` | NULL      | Treatments to appear once; typically `rep_plan$unreplicated_treatments` |
| `unreplicated_families`   | NULL      | Family labels for unreplicated treatments                               |
| `replication_mode`        | `"p_rep"` | `"p_rep"`, `"augmented"`, or `"rcbd_type"`                              |
| `cluster_source`          | `"none"`  | Genetic dispersion grouping: `"none"`, `"Family"`, `"GRM"`, `"A"`       |
| `eval_efficiency`         | `FALSE`   | Compute A, D, CDmean efficiency metrics                                 |
| `order`                   | `"row"`   | Plot traversal order: `"row"`, `"col"`, or `"serpentine"`               |
| `seed`                    | NULL      | Integer seed for reproducibility                                        |

### `met_alpha_rc_stream()` — row-column alpha design

**Required inputs**

| Argument           | Type             | Description                                           |
|--------------------|------------------|-------------------------------------------------------|
| `check_treatments` | character vector | Check treatment IDs; appear in every incomplete block |
| `check_families`   | character vector | Family labels for checks                              |
| `entry_treatments` | character vector | Entry (non-check) treatment IDs                       |
| `entry_families`   | character vector | Family labels for entries                             |
| `n_reps`           | integer          | Number of field replicates                            |
| `n_rows`           | integer          | Number of field rows                                  |
| `n_cols`           | integer          | Number of field columns                               |

**Optional inputs**

| Argument          | Default  | Description                                                       |
|-------------------|----------|-------------------------------------------------------------------|
| `min_block_size`  | 6        | Minimum entries (excluding checks) per incomplete block           |
| `max_block_size`  | NULL     | Maximum entries per incomplete block                              |
| `cluster_source`  | `"none"` | Genetic dispersion grouping: `"none"`, `"Family"`, `"GRM"`, `"A"` |
| `eval_efficiency` | `FALSE`  | Compute A, D, CDmean efficiency metrics                           |
| `order`           | `"row"`  | Plot traversal order: `"row"`, `"col"`, or `"serpentine"`         |
| `serpentine`      | `FALSE`  | Reverse alternating rows/columns for physical continuity          |
| `seed`            | NULL     | Integer seed for reproducibility                                  |

### Minimum working example

The absolute minimum to run the full pipeline from allocation to field
book:

``` r
library(OptiSparseMET)

## Minimum inputs: just lines, environments, and field dimensions
treatments <- paste0("L", sprintf("%03d", 1:120))
envs       <- c("E1", "E2", "E3", "E4")

## Stage 0: verify k
k <- suggest_safe_k(treatments, envs, buffer = 3)  # 33

## Stage 1: M3 allocation (no common treatments, no grouping)
alloc <- allocate_sparse_met(
  treatments                     = treatments,
  environments                   = envs,
  allocation_method              = "random_balanced",
  n_test_entries_per_environment = k,
  seed                           = 1
)

## Stage 2: seed plan (uniform seed, no shortage)
seed_df <- data.frame(
  Treatment     = treatments,
  SeedAvailable = 100L
)
rep_plan <- assign_replication_by_seed(
  treatments             = treatments,
  seed_available         = seed_df,
  seed_required_per_plot = 10L,
  replication_mode       = "augmented"
)

## Stage 3: within-environment design (checks + unreplicated entries)
design <- met_prep_famoptg(
  check_treatments        = c("CHK1", "CHK2"),
  check_families          = c("CHECK", "CHECK"),
  unreplicated_treatments = rep_plan$unreplicated_treatments,
  unreplicated_families   = rep("F1", length(rep_plan$unreplicated_treatments)),
  n_blocks = 4L, n_rows = 10L, n_cols = 12L,
  seed     = 1
)

## Stage 4: combine
met_book <- combine_met_fieldbooks(
  field_books = list(E1 = design$field_book)
)
```

------------------------------------------------------------------------

## Quick Start

``` r
library(OptiSparseMET)

treatments <- paste0("L", sprintf("%03d", 1:120))
envs       <- c("E1", "E2", "E3", "E4")
common     <- treatments[1:10]    # 10 common treatments -> J* = 110

# Step 0: verify minimum capacity before allocating
# J=120, C=10, I=4: J*=110, min k* = ceil(110/4) = 28, min k = 38
# suggest_safe_k adds a buffer of 3 on top of the minimum
suggest_safe_k(treatments, envs,
               common_treatments = common,
               buffer = 3)           # returns 41

# Before M4 allocation, verify the slot identity holds.
# With k=65, r=2, C=10: k*=55, J*=110. Identity: 110*2 = 4*55 = 220. Holds.
check_balanced_incomplete_feasibility(
  n_treatments_total             = 120,
  n_environments                 = 4,
  n_test_entries_per_environment = 65,
  target_replications            = 2,
  n_common_treatments            = 10
)
# feasible = TRUE: 110*2 = 4*55 = 220 (difference = 0)

# Step 1a: M3-inspired random balanced allocation
alloc_m3 <- allocate_sparse_met(
  treatments                     = treatments,
  environments                   = envs,
  allocation_method              = "random_balanced",
  n_test_entries_per_environment = 41,
  target_replications            = 1,
  common_treatments              = common,
  seed                           = 123
)

alloc_m3$summary$min_sparse_replication  # every treatment in >= 1 environment
alloc_m3$summary$mean_sparse_replication

# Step 1b: M4 BIBD allocation -- slot identity must hold exactly.
# J*=110, I=4, r=2: k* = 110*2/4 = 55, k = 55+10 = 65.
# Slot identity: 110*2 = 4*55 = 220. Satisfied.
alloc_m4 <- allocate_sparse_met(
  treatments                     = treatments,
  environments                   = envs,
  allocation_method              = "balanced_incomplete",
  n_test_entries_per_environment = 65,
  target_replications            = 2,
  common_treatments              = common,
  allow_approximate              = FALSE,
  seed                           = 123
)

alloc_m4$summary$min_sparse_replication  # 2: every sparse treatment in exactly 2 envs
alloc_m4$summary$max_sparse_replication  # 2: equal replication confirmed

# Step 2: seed-aware replication plan
seed_df <- data.frame(
  Treatment     = treatments,
  SeedAvailable = sample(10:100, length(treatments), replace = TRUE)
)

rep_plan <- assign_replication_by_seed(
  treatments             = treatments,
  seed_available         = seed_df,
  seed_required_per_plot = 10,
  replication_mode       = "p_rep",
  desired_replications   = 2,
  max_prep               = 15,
  shortage_action        = "downgrade"
)

rep_plan$p_rep_treatments        # treatments receiving 2 plots
rep_plan$unreplicated_treatments # treatments receiving 1 plot
rep_plan$excluded_treatments     # treatments with insufficient seed

# Step 3: within-environment field design
# For augmented / p-rep / RCBD-type block designs:
design_E1 <- met_prep_famoptg(
  check_treatments        = c("CHK1", "CHK2"),
  check_families          = c("CHECK", "CHECK"),
  p_rep_treatments        = rep_plan$p_rep_treatments,
  p_rep_reps              = rep_plan$p_rep_reps,
  p_rep_families          = rep("F1", length(rep_plan$p_rep_treatments)),
  unreplicated_treatments = rep_plan$unreplicated_treatments,
  unreplicated_families   = rep("F1", length(rep_plan$unreplicated_treatments)),
  n_blocks = 4L, n_rows = 10L, n_cols = 12L,
  seed     = 123
)

# For alpha row-column designs:
design_E2 <- met_alpha_rc_stream(
  check_treatments = c("CHK1", "CHK2"),
  check_families   = c("CHECK", "CHECK"),
  entry_treatments = rep_plan$p_rep_treatments,
  entry_families   = rep("F1", length(rep_plan$p_rep_treatments)),
  n_reps = 2L, n_rows = 8L, n_cols = 10L,
  seed   = 123
)

# Step 4: combine into single MET field book
met_book <- combine_met_fieldbooks(
  field_books       = list(E1 = design_E1$field_book,
                           E2 = design_E2$field_book),
  local_designs     = c(E1 = "met_prep_famoptg",
                        E2 = "met_alpha_rc_stream"),
  replication_modes = c(E1 = "p_rep", E2 = "met_alpha_rc_stream"),
  sparse_method     = "random_balanced",
  common_treatments = common
)

head(met_book[, 1:8])
```

### End-to-end pipeline

For standard workflows,
[`plan_sparse_met_design()`](https://FAkohoue.github.io/OptiSparseMET/reference/plan_sparse_met_design.md)
handles steps 1–4 in a single call using `env_design_specs`, a named
list where each environment is mapped to its local design arguments. Set
`design = "met_prep_famoptg"` or `design = "met_alpha_rc_stream"` in
each spec:

``` r
env_specs <- list(
  E1 = list(
    design               = "met_prep_famoptg",
    replication_mode     = "p_rep",
    desired_replications = 2L,
    max_prep             = 15L,
    shortage_action      = "downgrade",
    check_treatments     = c("CHK1", "CHK2"),
    check_families       = c("CHECK", "CHECK"),
    n_blocks = 4L, n_rows = 10L, n_cols = 12L
  ),
  E2 = list(
    design           = "met_alpha_rc_stream",
    check_treatments = c("CHK1", "CHK2"),
    check_families   = c("CHECK", "CHECK"),
    n_reps = 2L, n_rows = 8L, n_cols = 10L
  )
)

out <- plan_sparse_met_design(
  treatments                     = treatments,
  environments                   = envs[1:2],
  allocation_method              = "random_balanced",
  n_test_entries_per_environment = 41,
  target_replications            = 1,
  common_treatments              = common,
  env_design_specs               = env_specs,
  seed_info                      = seed_df,
  seed_required_per_plot         = data.frame(
    Environment         = envs[1:2],
    SeedRequiredPerPlot = c(10, 10)
  ),
  seed = 123
)

out$combined_field_book  # full MET field book
out$environment_summary  # per-environment design summary
out$efficiency_summary   # efficiency metrics (when eval_efficiency = TRUE)
```

------------------------------------------------------------------------

## 6.4 Slot identity feasibility by J\*, I, and r

The slot identity $J^{*} \times r = I \times k^{*}$ requires that
$J^{*} \times r$ be exactly divisible by $I$. Whether this is achievable
for a given combination of sparse treatments ($J^{*}$), environments
($I$), and replication ($r$) depends on the shared factors of these
three numbers.

### The divisibility rule

For the slot identity to hold, $I$ must divide $J^{*} \times r$ exactly.
Every prime factor of $I$ that is absent from $J^{*}$ must be supplied
by $r$. This has a practical consequence for the most common case in
plant breeding:

- **$I = 4$ environments, $J^{*}$ odd**: $J^{*} \times 1 = \text{odd}$
  (not divisible by 4); $J^{*} \times 2 = 2 \times \text{odd}$
  (divisible by 2 but not by 4 for odd $J^{*}$); only $r = 4$ guarantees
  divisibility. But $r = 4$ gives $k^{*} = J^{*}$ — full replication —
  which defeats the purpose of sparse testing. **Practical fix**: adjust
  $C$ by 1 so that $J^{*} = J - C$ becomes even.

- **$I = 4$ environments, $J^{*}$ even but not divisible by 4**: $r = 2$
  always works (e.g. $J^{*} = 110$: $110 \times 2/4 = 55$).

- **$I = 3$ environments**: feasibility depends on divisibility by 3. If
  $J^{*}$ is divisible by 3, any $r$ works. Otherwise $r$ must be a
  multiple of 3.

- **$I = 6$ environments**: requires divisibility by $2 \times 3 = 6$.
  Odd $J^{*}$ not divisible by 3 requires $r$ divisible by 6.

### Feasibility table: $r = 2$

The table shows $k^{*}$ when the slot identity holds, and `--` when it
does not for that combination. Add $C$ (common treatments) to $k^{*}$ to
get the `n_test_entries_per_environment` argument.

| $J^{*}$ | $I = 3$ | $I = 4$ | $I = 5$ | $I = 6$ | $I = 7$ | $I = 8$ | $I = 9$ | $I = 10$ |
|--------:|--------:|--------:|--------:|--------:|--------:|--------:|--------:|---------:|
|      60 |      40 |      30 |      24 |      20 |       – |      15 |       – |       12 |
|      70 |       – |      35 |      28 |       – |      20 |       – |       – |       14 |
|      75 |      50 |       – |      30 |      25 |       – |       – |       – |       15 |
|      80 |       – |      40 |      32 |       – |       – |      20 |       – |       16 |
|      90 |      60 |      45 |      36 |      30 |       – |       – |      20 |       18 |
|     100 |       – |      50 |      40 |       – |       – |      25 |       – |       20 |
|     110 |       – |      55 |      44 |       – |       – |       – |       – |       22 |
|     112 |       – |      56 |       – |       – |      32 |      28 |       – |        – |
|     120 |      80 |      60 |      48 |      40 |       – |      30 |       – |       24 |
|     150 |     100 |      75 |      60 |      50 |       – |       – |       – |       30 |
|     200 |       – |     100 |      80 |       – |       – |      50 |       – |       40 |

### Feasibility table: $r = 3$

| $J^{*}$ | $I = 3$ | $I = 4$ | $I = 5$ | $I = 6$ | $I = 7$ | $I = 8$ | $I = 9$ | $I = 10$ |
|--------:|--------:|--------:|--------:|--------:|--------:|--------:|--------:|---------:|
|      60 |      60 |      45 |      36 |      30 |       – |       – |      20 |       18 |
|      70 |      70 |       – |      42 |      35 |      30 |       – |       – |       21 |
|      75 |      75 |       – |      45 |       – |       – |       – |      25 |        – |
|      80 |      80 |      60 |      48 |      40 |       – |      30 |       – |       24 |
|      90 |      90 |       – |      54 |      45 |       – |       – |      30 |       27 |
|     100 |     100 |      75 |      60 |      50 |       – |       – |       – |       30 |
|     110 |     110 |       – |      66 |      55 |       – |       – |       – |       33 |
|     112 |     112 |      84 |       – |      56 |      48 |      42 |       – |        – |
|     120 |     120 |      90 |      72 |      60 |       – |      45 |      40 |       36 |
|     150 |     150 |       – |      90 |      75 |       – |       – |      50 |       45 |
|     200 |     200 |     150 |     120 |     100 |       – |      75 |       – |       60 |

### What to do when your combination gives `--`

Use
[`check_balanced_incomplete_feasibility()`](https://FAkohoue.github.io/OptiSparseMET/reference/check_balanced_incomplete_feasibility.md)
to diagnose the problem and try one of these adjustments:

1.  **Adjust $C$ by 1**: adding or removing one common treatment changes
    $J^{*}$ by 1, which may make it divisible by $I$ for the chosen $r$.
2.  **Try $r = 2$ instead of $r = 1$**, or $r = 3$ instead of $r = 2$ —
    the extra factor may resolve the divisibility.
3.  **Use `random_balanced` (M3)** if exact equal replication is not
    essential. M3 does not require the slot identity and tolerates odd
    $J^{*}$ freely.
4.  **Use `allow_approximate = TRUE`** as a fallback — the allocation
    proceeds with the closest possible balance, accepting minor
    replication differences.

``` r
## Quick check: is your combination feasible?
## J* = 75 (odd), I = 4, r = 2 -- should give --
check_balanced_incomplete_feasibility(
  n_treatments_total             = 83,   # J = J* + C = 75 + 8
  n_environments                 = 4,
  n_test_entries_per_environment = 30,   # k* guess: 30 - 8 = 22, 4*22=88 != 75*2=150
  target_replications            = 2,
  n_common_treatments            = 8
)
## feasible = FALSE -> adjust

## Fix: change C from 8 to 9 -> J* = 74 (even), r=2: 74*2/4 = 37
check_balanced_incomplete_feasibility(
  n_treatments_total             = 83,
  n_environments                 = 4,
  n_test_entries_per_environment = 46,   # k* = 37, k = 37+9 = 46
  target_replications            = 2,
  n_common_treatments            = 9
)
## feasible = TRUE
```

------------------------------------------------------------------------

## Design Strategy Notes

**Use `random_balanced`** when environment capacities differ
substantially, when exact BIBD parameters are not achievable, or when
some stochasticity in allocation is acceptable. Unlike the original M3
of Montesinos-Lopez et al. (2023), this implementation guarantees every
treatment appears in at least one environment before replication filling
begins.

**Use `balanced_incomplete` with `allow_approximate = FALSE`** (the
default) when equal replication is a hard requirement – every sparse
treatment must appear in exactly $r$ environments. Verify the slot
identity $J^{*} \times r = I \times k^{*}$ first with
[`check_balanced_incomplete_feasibility()`](https://FAkohoue.github.io/OptiSparseMET/reference/check_balanced_incomplete_feasibility.md).
The function stops with a clear error if the identity does not hold, so
you always know whether the equal-replication guarantee was met. Install
the `crossdes` package for guaranteed construction when available.

**Use `balanced_incomplete` with `allow_approximate = TRUE`** as a
fallback when the slot identity cannot be satisfied for the chosen
dimensions but you still want to attempt a balanced allocation. Some
lines will receive more or fewer replications than $r$. This is an
exploratory mode, not the primary path.

**Use
[`met_prep_famoptg()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_prep_famoptg.md)**
for designs requiring repeated checks in every block, partial
replication (p-rep), or RCBD-type structure. The optimizer
[`met_optimize_famoptg()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_optimize_famoptg.md)
accepts A, D, and CDmean criteria.

**Use
[`met_alpha_rc_stream()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_alpha_rc_stream.md)**
for fixed-grid field deployments where spatial row-column control is a
priority. The optimizer
[`met_optimize_alpha_rc()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_optimize_alpha_rc.md)
supports Random Restart, Simulated Annealing, and Genetic Algorithm
methods.

**Use GRM/A-based grouping** when genomic prediction is a primary
objective, when family labels are too coarse to capture relevant genetic
structure, or when related lines risk clustering into the same
environments.

**Include common treatments** whenever environments are weakly
correlated, when genetic connectivity cannot otherwise be guaranteed, or
when a stable reference set is needed for cross-environment
benchmarking.

------------------------------------------------------------------------

## Installation

Install from GitHub with vignettes (recommended):

``` r
install.packages("remotes")
remotes::install_github("FAkohoue/OptiSparseMET",
  build_vignettes = TRUE,
  dependencies    = TRUE
)
```

Install without vignettes for a faster install:

``` r
remotes::install_github("FAkohoue/OptiSparseMET",
  build_vignettes = FALSE,
  dependencies    = TRUE
)
```

For guaranteed exact BIBD construction via `crossdes`:

``` r
install.packages("crossdes")
```

------------------------------------------------------------------------

## Documentation

Full documentation, function reference, and tutorials are available at:

<https://FAkohoue.github.io/OptiSparseMET/>

To read the vignette after installation:

``` r
vignette("OptiSparseMET-introduction", package = "OptiSparseMET")
```

------------------------------------------------------------------------

## Citation

If you use `OptiSparseMET` in published research, please cite:

    Akohoue, F. (2026).
    OptiSparseMET: Sparse Multi-Environment Trial Design with Flexible Local
    Field Layout. R package version 0.1.0.
    https://github.com/FAkohoue/OptiSparseMET

------------------------------------------------------------------------

## Reference

Montesinos-Lopez O.A., Mosqueda-Gonzalez B.A., Salinas-Ruiz J.,
Montesinos-Lopez A., Crossa J. (2023). Sparse multi-trait genomic
prediction under balanced incomplete block design. *The Plant Genome*,
16, e20305. <https://doi.org/10.1002/tpg2.20305>

------------------------------------------------------------------------

## Contributing

Issues, bug reports, and feature suggestions are welcome:
<https://github.com/FAkohoue/OptiSparseMET/issues>

------------------------------------------------------------------------

## License

MIT License (c) Felicien Akohoue
