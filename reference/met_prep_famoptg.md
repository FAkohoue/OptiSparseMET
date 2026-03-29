# Construct a repeated-check block design with flexible replication

`met_prep_famoptg()` is the OptiSparseMET version of `prep_famoptg()`
from the OptiDesign package. The `met_` prefix avoids namespace
conflicts when both packages are loaded simultaneously. All arguments,
return values, and internal logic are identical to `prep_famoptg()` in
OptiDesign.

`met_prep_famoptg()` builds a repeated-check block design for plant
breeding and agronomic experiments. Check treatments are included in
every block, while non-check treatments are allocated according to
user-specified replication levels. Depending on the treatment structure
supplied, the same function generates three design classes:

- **Augmented repeated-check design**: checks repeated in every block,
  all test entries unreplicated (`p_rep_treatments = NULL`).

- **Partially replicated (p-rep) repeated-check design**: a mixture of
  replicated and unreplicated non-check entries.

- **RCBD-type repeated-check design**: all non-check entries given the
  same replication across blocks.

**Design evaluation has been separated from construction.** This
function returns the field book, layout matrix, and seed only. To
compute A, D, and CDmean optimality criteria call
[`met_evaluate_famoptg_efficiency()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_evaluate_famoptg_efficiency.md)
on the returned `field_book`. To search for an optimised design across
many randomisations call
[`met_optimize_famoptg()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_optimize_famoptg.md).

## Usage

``` r
met_prep_famoptg(
  check_treatments,
  check_families,
  p_rep_treatments,
  p_rep_reps,
  p_rep_families,
  unreplicated_treatments,
  unreplicated_families,
  n_blocks,
  n_rows,
  n_cols,
  order = "column",
  serpentine = FALSE,
  seed = NULL,
  attempts = 1000,
  warn_and_correct = TRUE,
  fix_rows = TRUE,
  cluster_source = c("Family", "GRM", "A"),
  GRM = NULL,
  A = NULL,
  id_map = NULL,
  cluster_method = c("kmeans", "hclust"),
  cluster_seed = 1,
  cluster_attempts = 25,
  n_pcs_use = Inf,
  check_placement = c("random", "systematic", "optimal"),
  check_opt_attempts = 200,
  use_dispersion = FALSE,
  dispersion_source = c("K", "A", "GRM"),
  dispersion_radius = 1,
  dispersion_iters = 2000,
  dispersion_seed = 1,
  K = NULL,
  line_id_map = NULL,
  verbose = TRUE
)
```

## Arguments

- check_treatments:

  Character vector of check treatment identifiers. Checks appear in
  every block. Must not overlap with `p_rep_treatments` or
  `unreplicated_treatments` and must not contain duplicates.

- check_families:

  Character vector of the same length as `check_treatments`. Family
  labels for checks, used for adjacency scoring. Always stored in the
  `Family` column; `Gcluster` is always `NA` for check plots.

- p_rep_treatments:

  Character vector of replicated non-check treatment identifiers. Each
  entry appears in exactly `p_rep_reps[i]` distinct blocks. May be
  `NULL` or `character(0)` for augmented designs.

- p_rep_reps:

  Integer vector of the same length as `p_rep_treatments`. Number of
  blocks each p-rep treatment appears in. Must satisfy
  `1 <= p_rep_reps[i] <= n_blocks` for all \\i\\.

- p_rep_families:

  Character vector of the same length as `p_rep_treatments`. Family
  labels for p-rep entries.

- unreplicated_treatments:

  Character vector of treatment identifiers that appear exactly once in
  the full design. May be `NULL` or `character(0)`.

- unreplicated_families:

  Character vector of the same length as `unreplicated_treatments`.
  Family labels for unreplicated entries.

- n_blocks:

  Positive integer. Number of experimental blocks.

- n_rows:

  Positive integer. Number of field rows.

- n_cols:

  Positive integer. Number of field columns.

- order:

  Character. Field traversal direction: `"column"` fills column by
  column; `"row"` fills row by row.

- serpentine:

  Logical. If `TRUE`, alternate columns (when `order = "column"`) or
  alternate rows (when `order = "row"`) traverse in reverse direction,
  producing a boustrophedon planting path.

- seed:

  Optional integer. If `NULL`, a seed is drawn from
  `1:.Machine$integer.max` and returned as `seed_used`. Controls all
  randomised steps: p-rep block assignment, unreplicated entry
  distribution, within-block shuffling, check placement, and dispersion
  optimisation.

- attempts:

  Positive integer. Number of within-block shuffle attempts to minimise
  adjacent same-group pairs. Default 1000.

- warn_and_correct:

  Logical. If `TRUE`, silently adjust field dimensions when
  `n_rows * n_cols` does not match the required total plots. If `FALSE`,
  stop with an error.

- fix_rows:

  Logical. Used only when `warn_and_correct = TRUE`. If `TRUE`, hold
  `n_rows` fixed and expand `n_cols`; otherwise hold `n_cols` fixed and
  expand `n_rows`.

- cluster_source:

  Character. Grouping source for adjacency scoring: `"Family"`, `"GRM"`,
  or `"A"`.

- GRM:

  Optional square numeric matrix with rownames and colnames equal to
  line IDs. Required when `cluster_source = "GRM"` or
  `dispersion_source = "GRM"`.

- A:

  Optional square numeric matrix. Required when `cluster_source = "A"`
  or `dispersion_source = "A"`.

- id_map:

  Optional data frame with columns `Treatment` and `LineID`. Required
  when treatment labels do not match relationship matrix rownames.

- cluster_method:

  Character. Clustering algorithm: `"kmeans"` or `"hclust"` (Ward's D2
  linkage).

- cluster_seed:

  Integer. Seed for k-means initialisation, run in an isolated RNG
  scope.

- cluster_attempts:

  Positive integer. Number of k-means random restarts (`nstart`).

- n_pcs_use:

  Positive number or `Inf`. Leading PCs retained for matrix-based
  clustering. Actual number used is
  `min(n_pcs_use, n_positive_eigenvalues - 1)`.

- check_placement:

  Character. Check position rule within blocks:

  `"random"`

  :   Check positions drawn uniformly.

  `"systematic"`

  :   Check positions spaced evenly through the block.

  `"optimal"`

  :   Keeps the layout that maximises mean nearest-neighbour distance
      between check plots in the field across `check_opt_attempts`
      candidates.

- check_opt_attempts:

  Positive integer. Number of candidate layouts evaluated when
  `check_placement = "optimal"`. Default 200.

- use_dispersion:

  Logical. If `TRUE`, apply a post-layout swap-based dispersion
  optimisation to reduce pairwise relatedness among neighbouring
  non-check plots.

- dispersion_source:

  Character. Matrix used for dispersion scoring: `"K"`, `"A"`, or
  `"GRM"`. The corresponding argument must be non-`NULL`.

- dispersion_radius:

  Positive integer. Chebyshev distance radius defining the neighbourhood
  for dispersion scoring. Default 1.

- dispersion_iters:

  Non-negative integer. Number of swap proposals for dispersion
  optimisation. Default 2000.

- dispersion_seed:

  Integer. Seed for the dispersion step, run in an isolated RNG scope.

- K:

  Optional square numeric matrix. Used for dispersion scoring when
  `dispersion_source = "K"`.

- line_id_map:

  Optional data frame with columns `Treatment` and `LineID`. Required
  when treatment labels do not match `rownames(K)` or the dispersion
  matrix.

- verbose:

  Logical. If `TRUE`, prints a one-line summary of field dimensions,
  block count, and treatment counts.

## Value

A named list with three components:

- `layout_matrix`:

  Character matrix of dimension `n_rows x n_cols`. Used cells contain
  treatment IDs; unused cells (when field is larger than required) are
  `NA`.

- `field_book`:

  Data frame with one row per assigned plot. Columns: `Treatment`,
  `Family` (original user-supplied label), `Gcluster` (`NA` when
  `cluster_source = "Family"` or for check plots; otherwise the
  matrix-derived cluster label prefixed `"G"` or `"A"`), `Block`
  (integer 1 to `n_blocks`), `Plot` (sequential integer), `Row`,
  `Column`.

- `seed_used`:

  Integer. The random seed used.

## Details

### Design structure and allocation rules

Let:

- \\c\\ = number of check treatments

- \\b\\ = `n_blocks`

- \\v_p\\ = number of p-rep treatments with replication counts \\r_1,
  r_2, \ldots, r\_{v_p}\\

- \\v_u\\ = number of unreplicated treatments

Total plots required: \$\$\text{total} = b \times c + \sum\_{i=1}^{v_p}
r_i + v_u\$\$

Allocation rules enforced by construction:

- Checks:

  Appear in every block. Each check contributes \\b\\ plots total.

- P-rep treatments:

  Each treatment \\i\\ appears in exactly \\r_i\\ blocks, always in
  distinct blocks - never twice in the same block. This is the core
  p-rep constraint.

- Unreplicated treatments:

  Appear exactly once in the full design, distributed as evenly as
  possible across blocks.

### Three design classes

**Augmented repeated-check design** - obtained when `p_rep_treatments`
is empty and all non-check entries are supplied via
`unreplicated_treatments`. Checks are the only repeated treatments.

**P-rep repeated-check design** - obtained when some non-check entries
in `p_rep_treatments` have `p_rep_reps > 1` and others are unreplicated.
The most general use case: a mixture of replicated candidate entries and
single-plot entries alongside repeated checks.

**RCBD-type repeated-check design** - obtained when all non-check
treatments are supplied via `p_rep_treatments` with a common replication
count greater than 1. If that count equals `n_blocks`, every non-check
treatment appears once in every block - the closest repeated-check
analogue of a classical RCBD.

### Field dimension handling

If the supplied `n_rows x n_cols` field does not exactly match the
required total plots, the function either stops (when
`warn_and_correct = FALSE`) or silently adjusts one dimension
(`warn_and_correct = TRUE`). When adjusting, `fix_rows = TRUE` holds
`n_rows` fixed and expands `n_cols`; `fix_rows = FALSE` does the
reverse. Surplus cells remain `NA` in the layout matrix.

### Within-block arrangement

Treatments are shuffled within each block over `attempts` random
permutations to minimise the count of adjacent same-group treatment
pairs in the 1D block ordering. The group label used for adjacency
scoring is determined by `cluster_source` (see below). Checks are
inserted at positions controlled by `check_placement`:

- `"random"`:

  Check positions drawn uniformly.

- `"systematic"`:

  Check positions spaced evenly through the block.

- `"optimal"`:

  Runs `check_opt_attempts` candidate layouts and keeps the one
  maximising the mean nearest-neighbour distance between check plots in
  the field - encouraging spatial spread of checks.

### Grouping and adjacency scoring

`cluster_source` determines the group label used to penalise adjacent
same-group treatments within blocks:

- `"Family"`:

  Uses `p_rep_families` and `unreplicated_families` directly. `Gcluster`
  column in the field book is `NA` for all plots.

- `"GRM"`:

  Derives clusters from PCA of the genomic relationship matrix `GRM`.
  Number of clusters = `length(unique(p_rep_families))`. Cluster labels
  are prefixed `"G"`.

- `"A"`:

  Derives clusters from PCA of the pedigree relationship matrix `A`.
  Cluster labels are prefixed `"A"`.

Check treatments always use `check_families` for grouping and always
have `Gcluster = NA`.

### Dispersion optimisation

When `use_dispersion = TRUE`, a local swap search runs after layout
construction. At each of `dispersion_iters` iterations two non-check
plots are selected at random; the swap is accepted if it reduces the
total neighbourhood relatedness score - the sum of relationship matrix
values between non-check neighbour pairs within Chebyshev distance
`dispersion_radius`. The matrix used is selected by `dispersion_source`.

## See also

[`met_evaluate_famoptg_efficiency()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_evaluate_famoptg_efficiency.md)
to compute A, D, and CDmean optimality criteria on the returned
`field_book`.
[`met_optimize_famoptg()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_optimize_famoptg.md)
to search for a criterion-optimal design using Random Restart.

## Examples

``` r
## Augmented repeated-check design: checks repeated, entries unreplicated
design_aug <- met_prep_famoptg(
  check_treatments        = c("CHK1", "CHK2"),
  check_families          = c("CHECK", "CHECK"),
  p_rep_treatments        = character(0),
  p_rep_reps              = integer(0),
  p_rep_families          = character(0),
  unreplicated_treatments = paste0("E", 1:80),
  unreplicated_families   = rep(paste0("F", 1:4), 20),
  n_blocks = 5, n_rows = 10, n_cols = 20
)
#> Warning: Field size (10 x 20 = 200) does not match required plots (90). Adjusting dimensions.
#> Field = 10 x 9 | total plots = 90 | checks/block = 2 | p-rep entries = 0 | unreplicated = 80 | n_blocks = 5
dim(design_aug$layout_matrix)
#> [1] 10  9
head(design_aug$field_book)
#>   Treatment Family Gcluster Block Plot Row Column
#> 1       E24     F4     <NA>     1    1   1      1
#> 2       E37     F1     <NA>     1    2   2      1
#> 3       E40     F4     <NA>     1    3   3      1
#> 4      CHK1  CHECK     <NA>     1    4   4      1
#> 5       E15     F3     <NA>     1    5   5      1
#> 6       E25     F1     <NA>     1    6   6      1

## P-rep design: some entries replicated twice, others once
design_prep <- met_prep_famoptg(
  check_treatments        = c("CHK1", "CHK2"),
  check_families          = c("CHECK", "CHECK"),
  p_rep_treatments        = paste0("P", 1:20),
  p_rep_reps              = rep(2L, 20),
  p_rep_families          = rep(paste0("F", 1:4), 5),
  unreplicated_treatments = paste0("U", 1:60),
  unreplicated_families   = rep(paste0("F", 1:4), 15),
  n_blocks = 5, n_rows = 15, n_cols = 20,
  check_placement = "systematic"
)
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#> Field = 15 x 8 | total plots = 120 | checks/block = 2 | p-rep entries = 20 | unreplicated = 60 | n_blocks = 5
table(design_prep$field_book$Block)
#> 
#>  1  2  3  4  5 
#> 22 23 20 24 21 

## Evaluate separately
eff <- met_evaluate_famoptg_efficiency(
  field_book         = design_prep$field_book,
  n_rows = 15, n_cols = 20,
  check_treatments   = c("CHK1", "CHK2"),
  treatment_effect   = "fixed",
  residual_structure = "AR1xAR1",
  rho_row = 0.10, rho_col = 0.10
)
eff$A_criterion
#> [1] 3.067164
```
