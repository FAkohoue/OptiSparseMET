# OptiSparseMET: A Unified Framework for Sparse Multi-Environment Trial Design

## 1. Introduction

The fundamental constraint in multi-environment trials is that the
number of candidate lines $J$ routinely exceeds the number of plots
available across the $I$ environments in the trial. Complete evaluation
of all lines in all environments is therefore infeasible, and some form
of incomplete allocation is required.

Sparse MET design formalizes this allocation problem. Rather than
assigning lines to environments arbitrarily or by convenience, a
principled sparse design controls replication, co-occurrence, and
genetic connectivity simultaneously — so that the resulting data support
the mixed-model inference the trial was built for. `OptiSparseMET`
addresses this at two levels: across-environment allocation and
within-environment field design, treated jointly rather than
sequentially.

This package is built for breeding workflows where design quality must
satisfy both statistical and operational constraints:

- limited field capacity per environment
- unequal environment sizes
- strong family or genomic structure among candidates
- seed availability constraints on replication
- need for mixed-model compatible layouts
- need for a unified MET-level field book

------------------------------------------------------------------------

## 2. Statistical Framework

The underlying model is a standard linear mixed model:

$$y = X\beta + Zg + e$$

where $g \sim N\left( 0,\, K\sigma_{g}^{2} \right)$ and
$e \sim N\left( 0,\, R\sigma_{e}^{2} \right)$. The matrix $K$ is either
a genomic relationship matrix (GRM) or a pedigree-based numerator
relationship matrix $A$, depending on the available data. The covariance
structure of $R$ can be identity (IID), or structured as AR1 or
AR1$\times$AR1 to account for spatial gradients within environments.

Design choices affect this model at multiple points. The incidence
structure of $Z$ — determined by allocation — governs the estimability
and precision of $\widehat{g}$. The blocking structure within
environments governs $R$ and the degree to which local field
heterogeneity is absorbed rather than confounded with genetic effects.
Both levels need to be well-constructed for the model to return reliable
estimates.

------------------------------------------------------------------------

## 3. Two-Level Design Structure

Sparse MET design in `OptiSparseMET` is intentionally divided into two
linked levels.

### 3.1 Across-Environment Allocation

The resource identity underlying balanced sparse allocation is:

$$J \times r = I \times k$$

where $J$ is the total number of treatments, $r$ is the number of
environments each treatment enters, $I$ is the number of environments,
and $k$ is the number of treatments tested per environment. Given fixed
total resources $N = J \times r = I \times k$, this identity makes the
tradeoff between coverage breadth ($k$) and replication depth ($r$)
explicit: increasing one necessarily decreases the other.

Allocation strategies in `OptiSparseMET` operate by constructing
treatment-by-environment incidence matrices that satisfy or approximate
this constraint while respecting genetic structure constraints and
common treatment requirements.

### 3.2 Within-Environment Design

Each environment receives an independent field design constructed by one
of two engines:

- [`prep_famoptg()`](https://FAkohoue.github.io/OptiSparseMET/reference/prep_famoptg.md)
  for block-based repeated-check designs, including augmented and p-rep
  layouts
- [`alpha_rc_stream()`](https://FAkohoue.github.io/OptiSparseMET/reference/alpha_rc_stream.md)
  for row-column alpha designs under large-scale spatial structure

The choice between them depends on field geometry, check structure,
replication goals, and the spatial model anticipated at analysis.

------------------------------------------------------------------------

## 4. Pipeline Overview

The package workflow proceeds in four stages:

- [`allocate_sparse_met()`](https://FAkohoue.github.io/OptiSparseMET/reference/allocate_sparse_met.md)
  distributes treatments across environments
- [`assign_replication_by_seed()`](https://FAkohoue.github.io/OptiSparseMET/reference/assign_replication_by_seed.md)
  verifies replication feasibility when seed constraints are relevant
- local design engines construct within-environment layouts
- [`combine_met_fieldbooks()`](https://FAkohoue.github.io/OptiSparseMET/reference/combine_met_fieldbooks.md)
  merges all layouts into a single MET-level field book

``` r
# install if needed: install.packages("DiagrammeR")
library(DiagrammeR)

grViz("
digraph OptiSparseMET {
  graph [layout = dot, rankdir = LR]
  node [shape = rectangle, style = rounded, fontname = Helvetica]

  A [label = 'Input data\nTreatments\nEnvironments\nFamilies / GRM / A\nSeed inventory']
  B [label = 'Across-environment allocation\nallocate_sparse_met()\n(M3 or M4)']
  C [label = 'Allocation matrix']
  D [label = 'Seed-aware replication\nassign_replication_by_seed()']
  E [label = 'Within-environment design\nprep_famoptg()\nor\nalpha_rc_stream()']
  F [label = 'Local field books']
  G [label = 'combine_met_fieldbooks()']
  H [label = 'Outputs\ncombined_field_book\nenvironment_summary\nefficiency_summary']

  A -> B -> C -> D -> E -> F -> G -> H
}
")
```

------------------------------------------------------------------------

## 5. Allocation Strategies

### 5.1 Random Balanced (M3)

The `random_balanced` strategy uses stochastic allocation to achieve
approximate balance without requiring exact BIBD parameters to be
satisfiable. It is appropriate when environment capacities are
heterogeneous, when the trial dimensions do not admit an exact balanced
solution, or when some flexibility in the incidence structure is
acceptable. Balance is approximate rather than exact, but the
construction is reproducible through the `seed` argument.

### 5.2 Balanced Incomplete (M4)

The `balanced_incomplete` strategy enforces equal replication across
environments and approximately uniform pairwise co-occurrence of
treatments, following BIBD principles applied at the MET level. It is
appropriate when environments are comparable in capacity, equal
replication is a hard requirement, and uniform co-occurrence is needed
for precision of treatment comparisons. When exact slot equality holds
and `allow_approximate = FALSE`, strict replication equality among
non-common treatments is enforced. When exact balance cannot be
satisfied, `allow_approximate = TRUE` relaxes the requirement.

------------------------------------------------------------------------

## 6. M3 versus M4

| Feature                 | M3 `random_balanced`       | M4 `balanced_incomplete` |
|-------------------------|----------------------------|--------------------------|
| Replication equality    | Approximate                | Exact or near-exact      |
| Pairwise co-occurrence  | Stochastic                 | More uniform             |
| Feasibility requirement | Flexible                   | Stricter                 |
| Typical use case        | Heterogeneous environments | Uniform MET networks     |
| Balance guarantee       | Heuristic                  | Design-based             |

In practice, M3 is often easier to use in real breeding pipelines. M4 is
preferable when exact balance is central to the design objective.

``` r
library(OptiSparseMET)

data("OptiSparseMET_example_data")
x <- OptiSparseMET_example_data

alloc_M3 <- allocate_sparse_met(
  treatments                     = x$treatments,
  environments                   = x$environments,
  allocation_method              = "random_balanced",
  n_test_entries_per_environment = 40,
  target_replications            = 2,
  common_treatments              = x$common_treatments,
  seed                           = 123
)

alloc_M4 <- allocate_sparse_met(
  treatments                     = x$treatments,
  environments                   = x$environments,
  allocation_method              = "balanced_incomplete",
  n_test_entries_per_environment = 40,
  target_replications            = 2,
  common_treatments              = x$common_treatments,
  allow_approximate              = TRUE,
  seed                           = 123
)

alloc_M3$summary
alloc_M4$summary
```

------------------------------------------------------------------------

## 7. Genetic Structure-Aware Allocation

Standard allocation treats lines as exchangeable. This is inadequate
when genomic prediction is a downstream objective, family structure is
strong, or environments risk receiving genetically clustered subsets of
entries. `OptiSparseMET` incorporates genetic structure through three
mechanisms.

### 7.1 Family-Based Allocation

Family-based allocation distributes related materials across
environments rather than concentrating them in a few sites, preventing
systematic differences in the genetic composition of environments that
would otherwise inflate GxE estimates.

### 7.2 GRM-Based Allocation

Genomic relationship matrix allocation uses marker-derived relatedness
to spread genetically similar lines across environments, improving the
conditions under which GBLUP models are estimated and reducing the risk
of genetic clustering.

### 7.3 Pedigree-Based Allocation

When genomic data are unavailable, the pedigree-based numerator
relationship matrix $A$ plays the same role as the GRM in guiding
allocation away from genetic clustering.

All three mechanisms improve cross-environment connectedness,
comparability of genetic effects across sites, stability of GxE
estimation, and prediction accuracy under GBLUP and PBLUP models.

------------------------------------------------------------------------

## 8. Common Treatments

Let $C$ denote the number of common treatments — treatments assigned to
every environment before sparse allocation begins. The effective
per-environment sparse capacity is:

$$k_{e}^{*} = k_{e} - C$$

Common treatments play two key roles. They establish direct
cross-environment connectivity that does not depend on model assumptions
about $K$, and they stabilize estimation of environment-level effects.
Without them, cross-environment linkage relies entirely on the
covariance structure of $K$, which is model-based rather than
design-based. Common treatments are most important when environments are
weakly correlated.

------------------------------------------------------------------------

## 9. Feasibility Helpers

Before allocation, it is important to verify that the chosen
per-environment capacity is sufficient to assign every non-common
treatment at least once. Passing a capacity that is too small causes
[`allocate_sparse_met()`](https://FAkohoue.github.io/OptiSparseMET/reference/allocate_sparse_met.md)
to stop with an informative error. The package provides three helpers
for this pre-flight check.

| Function                                                                                                     | Purpose                                                                        |
|--------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------|
| [`min_k_for_full_coverage()`](https://FAkohoue.github.io/OptiSparseMET/reference/min_k_for_full_coverage.md) | Compute the minimum entries per environment required to place all treatments   |
| [`suggest_safe_k()`](https://FAkohoue.github.io/OptiSparseMET/reference/suggest_safe_k.md)                   | Propose a safe uniform value of `n_test_entries_per_environment` with a buffer |
| [`warn_if_k_too_small()`](https://FAkohoue.github.io/OptiSparseMET/reference/warn_if_k_too_small.md)         | Emit a non-fatal warning when the chosen capacity is insufficient              |

These helpers are especially useful when many common treatments reduce
effective sparse capacity.

``` r
library(OptiSparseMET)

data("OptiSparseMET_example_data")
x <- OptiSparseMET_example_data

# Suggest a safe k with a small buffer above the strict minimum
suggest_safe_k(
  treatments        = x$treatments,
  environments      = x$environments,
  common_treatments = x$common_treatments,
  buffer            = 3
)

# Compute the strict minimum directly
min_k_for_full_coverage(
  n_treatments_total  = length(x$treatments),
  n_environments      = length(x$environments),
  n_common_treatments = length(x$common_treatments)
)

# Non-fatal check on a proposed capacity
warn_if_k_too_small(
  treatments                     = x$treatments,
  environments                   = x$environments,
  n_test_entries_per_environment = 40,
  common_treatments              = x$common_treatments
)
```

------------------------------------------------------------------------

## 10. Seed-Aware Replication

A design that requests replication without accounting for seed
availability is not deployable. The feasibility condition for assigning
replication $r_{i}$ to treatment $i$ is:

$$s_{i} \geq r_{i} \times q_{e}$$

where $s_{i}$ is the seeds available for treatment $i$, $r_{i}$ is the
assigned replication level, and $q_{e}$ is the seeds required per plot.
[`assign_replication_by_seed()`](https://FAkohoue.github.io/OptiSparseMET/reference/assign_replication_by_seed.md)
evaluates this condition and assigns each treatment one of three roles —
replicated, unreplicated (downgraded from the target), or excluded —
depending on the design mode and available seed. This step is especially
relevant before calling
[`prep_famoptg()`](https://FAkohoue.github.io/OptiSparseMET/reference/prep_famoptg.md)
in p-rep or RCBD-type repeated-check designs.

------------------------------------------------------------------------

## 11. Within-Environment Design Engines

### 11.1 `prep_famoptg()`

This function builds repeated-check block designs supporting three main
design classes:

- **Augmented repeated-check design**: checks repeated in every block,
  all test entries unreplicated
- **Partially replicated (p-rep) design**: a subset of entries
  replicated, the remainder unreplicated
- **RCBD-type repeated-check design**: all non-check entries replicated
  at a common level

Core rules enforced across all classes: checks appear in every block;
replicated non-check entries appear at most once per block; unreplicated
entries appear exactly once in the full design. Optionally, the function
minimizes same-group adjacency within blocks, applies dispersion
optimization using a relationship matrix, and evaluates mixed-model
efficiency metrics.

### 11.2 `alpha_rc_stream()`

This function builds stream-based alpha row-column designs over a fixed
grid. It is best suited when fields are large, row and column trends
matter, and spatial models such as AR1 or AR1$\times$AR1 are planned at
the analysis stage. Each entry appears exactly once per replicate and
checks appear once in every incomplete block. Replicate boundaries are
determined by position in the traversal stream rather than by geometric
subdivision of the grid.

------------------------------------------------------------------------

## 12. Design Engine Selection Guide

``` r
library(DiagrammeR)

grViz("
digraph design_choice {
  graph [layout = dot, rankdir = TB]
  node [shape = diamond, fontname = Helvetica]

  A [label = 'Is field geometry strongly row-column structured?']
  B [label = 'Is spatial modeling planned?']
  C [label = 'Use alpha_rc_stream()']
  D [label = 'Need augmented or p-rep repeated-check design?']
  E [label = 'Use prep_famoptg()']

  A -> B [label = 'yes']
  B -> C [label = 'yes']
  B -> D [label = 'no']
  A -> D [label = 'no']
  D -> E [label = 'yes']
  D -> C [label = 'no']
}
")
```

Use
[`prep_famoptg()`](https://FAkohoue.github.io/OptiSparseMET/reference/prep_famoptg.md)
when repeated checks are central, when a p-rep or augmented design is
needed, or when block structure is the primary local control mechanism.

Use
[`alpha_rc_stream()`](https://FAkohoue.github.io/OptiSparseMET/reference/alpha_rc_stream.md)
when the field is large and rectangular, when row and column variation
is substantial, or when a spatial residual model is anticipated at
analysis.

------------------------------------------------------------------------

## 13. Example: Allocation

``` r
library(OptiSparseMET)

data("OptiSparseMET_example_data")
x <- OptiSparseMET_example_data

alloc <- allocate_sparse_met(
  treatments                     = x$treatments,
  environments                   = x$environments,
  allocation_method              = "balanced_incomplete",
  n_test_entries_per_environment = 40,
  target_replications            = 2,
  common_treatments              = x$common_treatments,
  allow_approximate              = TRUE,
  seed                           = 123
)

alloc$summary
head(alloc$allocation_long)
alloc$overlap_matrix
```

------------------------------------------------------------------------

## 14. Example: Full Pipeline

``` r
out <- plan_sparse_met_design(
  treatments                     = x$treatments,
  environments                   = x$environments,
  allocation_method              = "random_balanced",
  n_test_entries_per_environment = 40,
  common_treatments              = x$common_treatments,
  env_design_specs               = x$env_design_specs,
  treatment_info                 = x$treatment_info,
  seed_info                      = x$seed_info,
  seed_required_per_plot         = x$seed_required_per_plot,
  seed                           = 123
)

names(out)
head(out$combined_field_book)
head(out$environment_summary)
```

------------------------------------------------------------------------

## 15. Environment-Specific Design Choice

One of the main strengths of the package is that the local design engine
can differ across environments. This is controlled through
`env_design_specs`, a named list with one element per environment. Each
element contains the design arguments for that environment plus a
`design` field specifying which engine to use.

``` r
env_design_specs <- list(

  # p-rep repeated-check design
  Env1 = list(
    design               = "prep_famoptg",
    replication_mode     = "p_rep",
    desired_replications = 2,
    check_treatments     = c("CHK1", "CHK2"),
    check_families       = c("CHECK", "CHECK"),
    n_blocks             = 4,
    n_rows               = 10,
    n_cols               = 10,
    cluster_source       = "Family"
  ),

  # alpha row-column stream design
  Env2 = list(
    design           = "alpha_rc_stream",
    check_treatments = c("CHK1", "CHK2"),
    check_families   = c("CHECK", "CHECK"),
    n_reps           = 2,
    n_rows           = 10,
    n_cols           = 10,
    cluster_source   = "Family"
  ),

  # augmented repeated-check design
  Env3 = list(
    design               = "prep_famoptg",
    replication_mode     = "augmented",
    check_treatments     = c("CHK1", "CHK2"),
    check_families       = c("CHECK", "CHECK"),
    n_blocks             = 4,
    n_rows               = 10,
    n_cols               = 10,
    cluster_source       = "Family"
  )
)
```

This allows the MET-level design to adapt to environment-specific
operational and logistical realities, including different field
geometries, different check sets, and different replication
requirements.

------------------------------------------------------------------------

## 16. Efficiency Evaluation

When local designs are built with efficiency evaluation enabled,
[`plan_sparse_met_design()`](https://FAkohoue.github.io/OptiSparseMET/reference/plan_sparse_met_design.md)
extracts the resulting metrics and stores them in `environment_summary`.
This makes it possible to compare design quality across environments
within the same MET before committing resources to field deployment.

Typical metrics include A-optimality-based summaries, D-efficiency, mean
prediction error variance (PEV), and the residual structure used in the
evaluation model.

``` r
out$environment_summary
```

------------------------------------------------------------------------

## 17. Connectivity and Genomic Prediction: Simulation Illustration

The benefit of better cross-environment connectivity can be explored
through simulation. The following code generates sparse allocations
under M3 and M4, simulates phenotypic data under a genomic model, and
fits a GBLUP model to each dataset. In many settings, the more regular
co-occurrence structure of M4 leads to lower average PEV, more stable
BLUPs, and stronger cross-environment connectivity. The exact magnitude
of the advantage depends on the covariance structure and degree of
imbalance.

``` r
library(sommer)
library(MASS)

set.seed(1)

n_geno <- 200
n_env  <- 4

G0 <- matrix(rnorm(n_geno^2), n_geno)
G0 <- tcrossprod(G0) / n_geno

simulate_MET <- function(allocation_matrix) {
  geno_ids <- rownames(allocation_matrix)
  env_ids  <- colnames(allocation_matrix)

  g_eff <- MASS::mvrnorm(1, mu = rep(0, length(geno_ids)), Sigma = G0)
  e_eff <- rnorm(length(env_ids), 0, 1)

  dat <- expand.grid(Geno = geno_ids, Env = env_ids, stringsAsFactors = FALSE)
  dat$present <- as.vector(allocation_matrix)
  dat <- dat[dat$present == 1, ]

  dat$y <- g_eff[match(dat$Geno, geno_ids)] +
            e_eff[match(dat$Env,  env_ids)]  +
            rnorm(nrow(dat), 0, 1)
  dat
}

alloc_M3 <- allocate_sparse_met(
  treatments                     = paste0("L", 1:n_geno),
  environments                   = paste0("E", 1:n_env),
  allocation_method              = "random_balanced",
  n_test_entries_per_environment = 80,
  target_replications            = 2,
  seed                           = 1
)

alloc_M4 <- allocate_sparse_met(
  treatments                     = paste0("L", 1:n_geno),
  environments                   = paste0("E", 1:n_env),
  allocation_method              = "balanced_incomplete",
  n_test_entries_per_environment = 80,
  target_replications            = 2,
  allow_approximate              = TRUE,
  seed                           = 1
)

dat_M3 <- simulate_MET(alloc_M3$allocation_matrix)
dat_M4 <- simulate_MET(alloc_M4$allocation_matrix)

fit_M3 <- sommer::mmer(
  y ~ Env,
  random = ~ vsr(Geno, Gu = G0),
  data   = dat_M3
)

fit_M4 <- sommer::mmer(
  y ~ Env,
  random = ~ vsr(Geno, Gu = G0),
  data   = dat_M4
)
```

------------------------------------------------------------------------

## 18. Design Component Summary

| Component                     | Governing principle                                                     |
|-------------------------------|-------------------------------------------------------------------------|
| Across-environment allocation | Balance replication depth against coverage breadth                      |
| Genetic connectivity          | Common treatments provide model-free cross-environment linkage          |
| Genetic structure             | GRM / $A$ / family grouping prevents clustering and supports prediction |
| Within-environment design     | Block and spatial structure control local heterogeneity                 |
| Replication feasibility       | Seed inventory determines what is deployable                            |
| Efficiency evaluation         | Mixed-model diagnostics quantify design quality prior to deployment     |

------------------------------------------------------------------------

## 19. Conclusions

`OptiSparseMET` treats MET design as a two-level problem and provides
machinery to solve both levels in a consistent statistical framework.
Allocation decisions are tied to within-environment design through
shared constraints on genetic structure and connectivity. Seed
feasibility is evaluated against replication targets before a design is
finalized. Efficiency can be evaluated locally and summarized globally
across environments.

The result is a pipeline that produces designs which are statistically
grounded and operationally deployable, rather than optimizing one at the
expense of the other.

------------------------------------------------------------------------

## 20. References

Montesinos-López, O. A., Mosqueda-González, B. A., Salinas-Ruiz, J.,
Montesinos-López, A., & Crossa, J. (2023). Sparse multi-trait genomic
prediction under balanced incomplete block design. *The Plant Genome*,
16, e20305.
