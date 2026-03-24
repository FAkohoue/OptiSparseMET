# OptiSparseMET

`OptiSparseMET` is an R framework for constructing sparse multi-environment trials (METs) that jointly addresses treatment allocation across environments and field design within environments — under shared statistical, genetic, and logistical constraints.

The package targets a structural challenge common to modern breeding programs: the number of candidate lines routinely exceeds what any single environment can accommodate, environments are heterogeneous rather than interchangeable, and seed availability imposes limits that purely theoretical designs ignore. Standard MET tools typically address one of these problems at a time. `OptiSparseMET` integrates all of them within a single workflow.

---

## Conceptual Framework

Sparse MET design is a two-level problem and should be treated as such.

**Level 1 — Across-environment allocation** determines which treatments appear in which environments, how many times each treatment is replicated across the trial, and whether the resulting incidence structure preserves sufficient genetic connectivity for valid cross-environment inference.

**Level 2 — Within-environment design** determines blocking structure, spatial layout, replication within each environment, and local control of field heterogeneity.

These two levels are statistically linked. Allocation decisions constrain which field layouts are feasible within environments; field layout choices affect the precision with which allocation-level effects can be estimated. Optimizing each level in isolation and patching them together is not equivalent to joint optimization. `OptiSparseMET` formalizes the link between them.

---

## Statistical Foundations

### Sparse Testing Identity

The theoretical basis for sparse allocation follows Montesinos-López et al. (2023), who formalize the resource identity underlying balanced sparse designs:

$$N = J \times r = I \times k$$

| Symbol | Meaning |
|---|---|
| $J$ | total number of treatments |
| $I$ | number of environments |
| $k$ | number of treatments per environment |
| $r$ | number of environments per treatment |

Given fixed total resources $N$, this identity makes the tradeoff between coverage breadth ($k$) and replication depth ($r$) explicit: increasing one necessarily decreases the other.

### Allocation Strategies

Two allocation strategies are implemented.

| Strategy | Argument | Properties |
|---|---|---|
| M3 | `random_balanced` | Stochastic allocation with approximate balance; tolerates unequal environment sizes |
| M4 | `balanced_incomplete` | BIBD-inspired allocation with enforced equal replication and uniform co-occurrence |

The `balanced_incomplete` strategy constructs a balanced incomplete incidence structure at the MET level. When the design dimensions admit an exact solution, equal replication is guaranteed. When exact balance is not achievable, approximate balance is available via `allow_approximate = TRUE`.

### Feasibility Helpers

Before running allocation, the package provides utilities to verify that the chosen per-environment capacity is sufficient to assign every non-common treatment at least once. Passing a capacity that is too small causes allocation to either fail or leave some treatments unassigned.

| Function | Purpose |
|---|---|
| `min_k_for_full_coverage()` | Compute the minimum entries per environment needed to place all treatments |
| `suggest_safe_k()` | Propose a safe uniform value of `n_test_entries_per_environment` given the trial dimensions |
| `warn_if_k_too_small()` | Emit a non-fatal warning when the chosen capacity is insufficient |

Call `suggest_safe_k()` or `min_k_for_full_coverage()` before `allocate_sparse_met()` when the trial dimensions are new or when the number of treatments and environments has changed.

### Genetic Connectedness

Genetic disconnectedness between environments — where different genetic material is tested in different locations with no shared lines — partially confounds environment effects and genetic effects. `OptiSparseMET` addresses this through three mechanisms.

**Common treatments** are forced into every environment before sparse allocation begins, establishing a baseline of direct cross-environment connectivity that does not depend on model assumptions.

**Family-based allocation** distributes each family group across environments rather than concentrating it in a subset of them, preventing systematic differences in the genetic composition of environments that would inflate GxE estimates.

**GRM/A-based allocation** uses a genomic relationship matrix (GRM) or pedigree-based numerator relationship matrix (A) to guide how related lines are spread across environments, avoiding genetic clustering that degrades prediction accuracy under GBLUP and PBLUP models.

### Within-Environment Field Design

`prep_famoptg()` builds augmented, p-rep, or repeated-check RCBD-type designs. Checks appear in every block, test entries appear at most once per block, and adjacency between related entries is minimized. When a relationship matrix is supplied, spatial dispersion of related lines within the field layout is additionally controlled.

`alpha_rc_stream()` constructs alpha-type row-column designs for large-scale field deployment with explicit spatial structure. These are appropriate when field dimensions impose simultaneous row and column blocking and when AR1 or AR1×AR1 spatial covariance models are anticipated at the analysis stage.

### Mixed-Model Efficiency

Designs can be evaluated under mixed-model frameworks before field deployment. Supported models include fixed-effects contrast precision, IID random effects, GBLUP and PBLUP using GRM or A, and spatial models with AR1 or AR1×AR1 error structure. The intent is to let the user verify that a proposed design will support the planned inference before committing resources.

### Seed Constraints

`assign_replication_by_seed()` takes a data frame of available seed quantities and a per-plot seed requirement, and returns a replication plan that respects those constraints. Lines with insufficient seed for the target replication are flagged or excluded. This is what makes designs produced by the package deployable rather than merely theoretically optimal.

---

## Workflow

```
STEP 0: Verify capacity
        ↓
        suggest_safe_k()  OR  min_k_for_full_coverage()

STEP 1: Allocate treatments across environments
        ↓
        allocate_sparse_met()

STEP 2: Define replication based on seed availability
        ↓
        assign_replication_by_seed()

STEP 3: Build within-environment field designs
        ↓
        prep_famoptg()  OR  alpha_rc_stream()

STEP 4: Assemble the combined MET field book
        ↓
        combine_met_fieldbooks()
```

---

## Minimal Reproducible Example

```r
library(OptiSparseMET)

treatments <- paste0("L", sprintf("%03d", 1:120))
envs       <- c("E1", "E2", "E3", "E4")

# Step 0: verify that the chosen k is sufficient for full coverage
suggest_safe_k(treatments, envs, buffer = 3)  # returns 33

# Step 1: allocation
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

# Step 3: within-environment design (see prep_famoptg() documentation)
design_E1 <- prep_famoptg(...)

# Step 4: combine into MET field book
met <- combine_met_fieldbooks(...)
```

---

## Design Strategy Notes

**Use `random_balanced`** when environment capacities differ substantially, when exact BIBD parameters are not achievable given the trial dimensions, or when a degree of stochasticity in allocation is acceptable.

**Use `balanced_incomplete`** when environments are comparable in size, when equal replication across environments is a hard requirement, and when the downstream analysis will rely on the precision guarantees that uniform co-occurrence provides.

**Use GRM/A-based grouping** when genomic prediction is a primary objective, when family labels are too coarse to capture relevant genetic structure, or when there is a risk that related lines will cluster into the same environments under simpler allocation rules.

**Include common treatments** whenever environments are weakly correlated with each other, when genetic connectivity cannot otherwise be guaranteed, or when a stable reference set is needed for cross-environment benchmarking.

---

## Installation

**Build Vignettes**

```r
# install.packages("remotes")
remotes::install_github("FAkohoue/OptiSparseMET", build_vignettes = TRUE,
  dependencies = TRUE
)
```

**Without Vignettes** 

```r
# install.packages("remotes")
remotes::install_github("FAkohoue/OptiSparseMET", build_vignettes = FALSE,
  dependencies = TRUE
)
```
---

## Documentation

Full documentation and tutorials are available at: https://FAkohoue.github.io/OptiSparseMET/

---

## Citation

Akohoue, F. (2026). OptiSparseMET: Sparse Multi-Environment Trial Design with Flexible Local Field Layout Construction. https://github.com/FAkohoue/OptiSparseMET

## Contributing

Issues and suggestions are welcome: https://github.com/FAkohoue/OptiSparseMET/issues

---

## Reference

Montesinos-López, O. A., Mosqueda-González, B. A., Salinas-Ruiz, J., Montesinos-López, A., & Crossa, J. (2023). Sparse multi-trait genomic prediction under balanced incomplete block design. *The Plant Genome*, 16, e20305.