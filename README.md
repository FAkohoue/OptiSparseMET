# OptiSparseMET

<p align="center">
  <img src="man/figures/logo.svg" alt="OptiSparseMET logo" width="100%">
</p>

<!-- badges: start -->
[![R-CMD-check](https://github.com/FAkohoue/OptiSparseMET/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/FAkohoue/OptiSparseMET/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/FAkohoue/OptiSparseMET/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/FAkohoue/OptiSparseMET/actions/workflows/pkgdown.yaml)
<!-- badges: end -->

---

## Overview

`OptiSparseMET` is an R framework for constructing sparse multi-environment
trials (METs) that jointly addresses treatment allocation across environments
and field design within environments — under shared statistical, genetic, and
logistical constraints.

The package targets a structural challenge common to modern breeding programs:
the number of candidate lines routinely exceeds what any single environment
can accommodate, environments are heterogeneous rather than interchangeable,
and seed availability imposes limits that purely theoretical designs ignore.
Standard MET tools typically address one of these problems at a time.
`OptiSparseMET` integrates all of them within a single workflow.

---

## Conceptual Framework

Sparse MET design is a two-level problem and should be treated as such.

**Level 1 — Across-environment allocation** determines which treatments appear
in which environments, how many times each treatment is replicated across the
trial, and whether the resulting incidence structure preserves sufficient
genetic connectivity for valid cross-environment inference.

**Level 2 — Within-environment design** determines blocking structure, spatial
layout, replication within each environment, and local control of field
heterogeneity.

These two levels are statistically linked. Allocation decisions constrain
which field layouts are feasible within environments; field layout choices
affect the precision with which allocation-level effects can be estimated.
Optimizing each level in isolation is not equivalent to joint optimization.
`OptiSparseMET` formalizes the link between them.

---

## Statistical Foundations

### Sparse testing identity

The theoretical basis follows Montesinos-López et al. (2023), who formalize
the resource identity underlying balanced sparse designs:

$$N = J \times r = I \times k$$

| Symbol | Meaning |
|--------|---------|
| $J$ | Total number of treatments |
| $I$ | Number of environments |
| $k$ | Number of treatments per environment |
| $r$ | Number of environments per treatment |

Given fixed total resources $N$, this identity makes the tradeoff between
coverage breadth ($k$) and replication depth ($r$) explicit.

### Allocation strategies

| Strategy | Argument | Properties |
|----------|----------|------------|
| M3 | `random_balanced` | Stochastic allocation with approximate balance; tolerates unequal environment sizes |
| M4 | `balanced_incomplete` | BIBD-inspired allocation with enforced equal replication and uniform co-occurrence |

### Genetic connectedness

`OptiSparseMET` addresses genetic disconnectedness through three mechanisms:
**common treatments** forced into every environment; **family-based
allocation** that distributes each family group across environments; and
**GRM/A-based allocation** that uses genomic or pedigree relationships to
prevent genetic clustering across environments.

### Seed constraints

`assign_replication_by_seed()` takes a data frame of available seed quantities
and a per-plot seed requirement and returns a replication plan that respects
those constraints — making designs deployable rather than merely theoretically
optimal.

---

## Main Functions

### Across-environment allocation

| Function | Description |
|----------|-------------|
| `allocate_sparse_met()` | Distribute treatments across environments (M3 or M4) |
| `check_balanced_incomplete_feasibility()` | Verify BIBD feasibility given trial dimensions |
| `derive_allocation_groups()` | Derive grouping structure from family, GRM, or A matrix |

### Feasibility and capacity helpers

| Function | Description |
|----------|-------------|
| `suggest_safe_k()` | Propose a safe `n_test_entries_per_environment` value |
| `min_k_for_full_coverage()` | Compute minimum capacity for full treatment coverage |
| `warn_if_k_too_small()` | Non-fatal pre-flight capacity check |

Call `suggest_safe_k()` or `min_k_for_full_coverage()` before
`allocate_sparse_met()` whenever trial dimensions change.

### Within-environment field design

| Function | Description |
|----------|-------------|
| `prep_famoptg()` | Augmented, p-rep, or repeated-check block designs |
| `alpha_rc_stream()` | Alpha row-column designs for fixed-grid field deployment |

### Pipeline and assembly

| Function | Description |
|----------|-------------|
| `plan_sparse_met_design()` | End-to-end pipeline in a single call |
| `combine_met_fieldbooks()` | Stack environment-level field books into one MET field book |

---

## Workflow

```
STEP 0  Verify capacity
        suggest_safe_k()  OR  min_k_for_full_coverage()
        ↓
STEP 1  Allocate treatments across environments
        allocate_sparse_met()
        ↓
STEP 2  Define replication based on seed availability
        assign_replication_by_seed()
        ↓
STEP 3  Build within-environment field designs
        prep_famoptg()  OR  alpha_rc_stream()
        ↓
STEP 4  Assemble the combined MET field book
        combine_met_fieldbooks()
```

---

## Quick Start

```r
library(OptiSparseMET)

treatments <- paste0("L", sprintf("%03d", 1:120))
envs       <- c("E1", "E2", "E3", "E4")

# Step 0: verify capacity
suggest_safe_k(treatments, envs, buffer = 3)  # returns 33

# Step 1: allocate
alloc <- allocate_sparse_met(
  treatments                     = treatments,
  environments                   = envs,
  allocation_method              = "balanced_incomplete",
  n_test_entries_per_environment = 50,
  target_replications            = 2,
  common_treatments              = treatments[1:10],
  allow_approximate              = TRUE,
  seed                           = 123
)

# Step 2: seed-aware replication
seed_df <- data.frame(
  Treatment     = treatments,
  SeedAvailable = sample(10:100, length(treatments), replace = TRUE)
)

rep_plan <- assign_replication_by_seed(
  treatments             = treatments,
  seed_available         = seed_df,
  seed_required_per_plot = 10
)

# Step 3: within-environment design
design_E1 <- prep_famoptg(...)

# Step 4: combine into MET field book
met <- combine_met_fieldbooks(...)
```

---

## Design Strategy Notes

**Use `random_balanced`** when environment capacities differ substantially,
when exact BIBD parameters are not achievable, or when some stochasticity in
allocation is acceptable.

**Use `balanced_incomplete`** when environments are comparable in size, when
equal replication is a hard requirement, and when the downstream analysis
relies on the precision guarantees that uniform co-occurrence provides.

**Use GRM/A-based grouping** when genomic prediction is a primary objective,
when family labels are too coarse to capture relevant genetic structure, or
when related lines risk clustering into the same environments.

**Include common treatments** whenever environments are weakly correlated,
when genetic connectivity cannot otherwise be guaranteed, or when a stable
reference set is needed for cross-environment benchmarking.

---

## Installation

Install from GitHub with vignettes (recommended):

```r
install.packages("remotes")
remotes::install_github("FAkohoue/OptiSparseMET",
  build_vignettes = TRUE,
  dependencies    = TRUE
)
```

Install without vignettes for a faster install:

```r
remotes::install_github("FAkohoue/OptiSparseMET",
  build_vignettes = FALSE,
  dependencies    = TRUE
)
```

---

## Documentation

Full documentation, function reference, and tutorials are available at:

<https://FAkohoue.github.io/OptiSparseMET/>

To read the vignette after installation:

```r
vignette("OptiSparseMET-introduction", package = "OptiSparseMET")
```

---

## Citation

If you use `OptiSparseMET` in published research, please cite:

```
Akohoue, F. (2026).
OptiSparseMET: Sparse Multi-Environment Trial Design with Flexible Local
Field Layout Construction. R package version 0.1.0.
https://github.com/FAkohoue/OptiSparseMET
```

---

## Reference

Montesinos-López O.A., Mosqueda-González B.A., Salinas-Ruiz J.,
Montesinos-López A., Crossa J. (2023). Sparse multi-trait genomic prediction
under balanced incomplete block design. *The Plant Genome*, 16, e20305.

---

## Contributing

Issues, bug reports, and feature suggestions are welcome:
<https://github.com/FAkohoue/OptiSparseMET/issues>

---

## License

MIT License © Félicien Akohoue