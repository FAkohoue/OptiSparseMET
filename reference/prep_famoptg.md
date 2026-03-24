# Construct a repeated-check block design with flexible replication structure

The allocation rules that govern all three design classes are:

- Check treatments appear exactly once in every block.

- Replicated non-check treatments appear at most once per block; their
  replicates are distributed across distinct blocks.

- Unreplicated treatments appear exactly once in the full design.

A treatment with replication \\r\\ and \\b \geq r\\ blocks therefore
occupies \\r\\ distinct blocks. When \\r = b\\, the treatment appears in
every block, which is the RCBD analogue within this framework.

Within each block, the 1D ordering of treatments is shuffled to reduce
adjacency between entries from the same genetic group, using up to
`attempts` shuffle proposals per block. This adjacency control operates
on group labels derived from `cluster_source`. Optional dispersion
optimization applies a swap-based local search after the full layout is
placed, using pairwise relatedness from a relationship matrix to
penalize spatial proximity of similar entries in the 2D field grid.

## Usage

``` r
prep_famoptg(
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
  eval_efficiency = FALSE,
  treatment_effect = c("random", "fixed"),
  prediction_type = c("none", "IID", "GBLUP", "PBLUP"),
  K = NULL,
  line_id_map = NULL,
  varcomp = list(sigma_e2 = 1, sigma_g2 = 1, sigma_b2 = 1, sigma_r2 = 1, sigma_c2 = 1),
  check_as_fixed = TRUE,
  residual_structure = c("IID", "AR1", "AR1xAR1"),
  rho_row = 0,
  rho_col = 0,
  spatial_engine = c("auto", "sparse", "dense"),
  dense_max_n = 5000,
  eff_trace_samples = 80,
  eff_full_max = 400,
  check_placement = c("random", "systematic", "optimal"),
  check_opt_attempts = 200,
  use_dispersion = FALSE,
  dispersion_source = c("K", "A", "GRM"),
  dispersion_radius = 1,
  dispersion_iters = 2000,
  dispersion_seed = 1
)
```

## Arguments

- check_treatments:

  Character vector of check treatment IDs. Checks appear exactly once in
  every block and define the repeated reference structure of the design.
  Required for all design classes. Must not overlap with
  `p_rep_treatments` or `unreplicated_treatments`.

- check_families:

  Character vector of the same length as `check_treatments`. Family or
  group labels for check treatments. Required regardless of
  `cluster_source` because checks always carry a group label in the
  output field book. Must align element-wise with `check_treatments`.

- p_rep_treatments:

  Character vector of replicated non-check treatment IDs. Each treatment
  in this vector appears `p_rep_reps[i]` times in the full design,
  distributed across that many distinct blocks. Supply `character(0)` or
  `NULL` for an augmented design with no replicated non-check entries.
  Must align element-wise with `p_rep_reps` and `p_rep_families`.

- p_rep_reps:

  Integer vector of replication counts, one per element of
  `p_rep_treatments`. Each value must satisfy
  `p_rep_reps[i] <= n_blocks`. A treatment with replication \\r\\ is
  placed in \\r\\ distinct blocks. Supply `integer(0)` when
  `p_rep_treatments` is empty.

- p_rep_families:

  Character vector of the same length as `p_rep_treatments`. Family or
  group labels for replicated non-check treatments. Used directly when
  `cluster_source = "Family"`. When `cluster_source %in% c("GRM", "A")`,
  still used to anchor the number of clusters to the number of distinct
  families among non-check treatments. Supply `character(0)` when
  `p_rep_treatments` is empty.

- unreplicated_treatments:

  Character vector of non-check treatments assigned exactly one plot
  each. Supply `character(0)` or `NULL` when all non-check treatments
  are replicated. In combination with empty `p_rep_treatments` and
  repeated checks, these entries define an augmented repeated-check
  design. Must align element-wise with `unreplicated_families`.

- unreplicated_families:

  Character vector of the same length as `unreplicated_treatments`.
  Family or group labels for unreplicated treatments. Required whenever
  `unreplicated_treatments` is non-empty. Supply `character(0)` when
  `unreplicated_treatments` is empty.

- n_blocks:

  Positive integer. Number of experimental blocks. Determines how many
  times checks are repeated and sets the upper bound on replication for
  any non-check treatment: `p_rep_reps[i] <= n_blocks`. When a
  replicated treatment should appear in every block, its replication
  must equal `n_blocks`.

- n_rows:

  Positive integer. Number of rows in the field grid. Together with
  `n_cols`, determines total field capacity and the 2D coordinates in
  `layout_matrix` and `field_book`. When `warn_and_correct = TRUE` and
  `fix_rows = TRUE`, `n_rows` is held fixed and `n_cols` is adjusted.

- n_cols:

  Positive integer. Number of columns in the field grid. When
  `warn_and_correct = TRUE` and `fix_rows = FALSE`, `n_cols` is held
  fixed and `n_rows` is adjusted.

- order:

  Character scalar. Direction in which the field grid is traversed when
  mapping the ordered treatment sequence to spatial coordinates. `"row"`
  fills row by row; `"column"` fills column by column. Should reflect
  the physical planting or harvesting direction. Interacts with
  `serpentine`.

- serpentine:

  Logical. If `TRUE`, alternate rows (when `order = "row"`) or alternate
  columns (when `order = "column"`) reverse traversal direction,
  producing a boustrophedon mapping from treatment sequence to field
  coordinates. Does not affect block composition; affects only the
  spatial positions assigned to each treatment.

- seed:

  Optional integer. Random seed controlling block assignment of
  replicated treatments, within-block shuffling, check placement
  candidates, and dispersion optimization when `dispersion_seed` is
  `NULL`. If `NULL`, a seed is generated internally and returned as
  `seed_used`.

- attempts:

  Positive integer. Maximum number of shuffle proposals per block used
  to reduce same-group adjacency in the 1D block ordering. Larger values
  increase search effort without changing the algorithm. Most relevant
  when group structure is highly imbalanced within blocks or when many
  similar entries occur together.

- warn_and_correct:

  Logical. If `FALSE`, the function stops when `n_rows * n_cols` does
  not equal the total required plot count. If `TRUE`, one dimension is
  expanded upward to accommodate all plots, with surplus cells left as
  `NA` in `layout_matrix`. Which dimension is adjusted is controlled by
  `fix_rows`.

- fix_rows:

  Logical. Active only when `warn_and_correct = TRUE`. If `TRUE`,
  `n_rows` is held fixed and `n_cols` is adjusted. If `FALSE`, `n_cols`
  is held fixed and `n_rows` is adjusted. Use `TRUE` when the number of
  field rows is physically constrained.

- cluster_source:

  Character scalar. Source of group labels for within-block adjacency
  control. `"Family"` uses the supplied family vectors. `"GRM"` derives
  cluster labels from the eigenstructure of `GRM`. `"A"` derives cluster
  labels from the eigenstructure of `A`. Determines which of `GRM`, `A`,
  `id_map`, `cluster_method`, `cluster_seed`, `cluster_attempts`, and
  `n_pcs_use` become active.

- GRM:

  Optional numeric matrix. Genomic relationship matrix. Required when
  `cluster_source = "GRM"` or when `use_dispersion = TRUE` and
  `dispersion_source = "GRM"`. Must be square with row and column names
  matching treatment IDs or reachable through `id_map`.

- A:

  Optional numeric matrix. Pedigree-based numerator relationship matrix.
  Required when `cluster_source = "A"` or when `use_dispersion = TRUE`
  and `dispersion_source = "A"`. Same naming requirements as `GRM`.

- id_map:

  Optional data frame with columns `Treatment` and `LineID`. Required
  only when `cluster_source %in% c("GRM", "A")` and treatment IDs in the
  design do not match the row names of the selected relationship matrix.

- cluster_method:

  Character scalar. Clustering algorithm applied to PCA scores in
  matrix-based grouping. `"kmeans"` uses k-means with `cluster_attempts`
  random restarts seeded by `cluster_seed`. `"hclust"` uses Ward's
  criterion hierarchical clustering. Ignored when
  `cluster_source = "Family"`.

- cluster_seed:

  Integer. Seed for k-means initialization. Active only when
  `cluster_source %in% c("GRM", "A")` and `cluster_method = "kmeans"`.

- cluster_attempts:

  Integer. Number of random restarts for k-means. Active only when
  `cluster_source %in% c("GRM", "A")` and `cluster_method = "kmeans"`.
  Larger values reduce the risk of poor local optima.

- n_pcs_use:

  Integer or `Inf`. Number of leading principal components retained for
  PCA-based clustering. `Inf` retains all components corresponding to
  positive eigenvalues. Smaller values preserve only broad genetic
  structure. Must be at least 2. Ignored when
  `cluster_source = "Family"`.

- eval_efficiency:

  Logical. If `TRUE`, mixed-model efficiency metrics are computed on the
  final design. If `FALSE`, all efficiency arguments are ignored and
  `efficiency` is returned as `NULL`.

- treatment_effect:

  Character scalar. Specifies whether non-check treatments enter the
  efficiency model as fixed or random effects. `"fixed"` computes
  precision of pairwise treatment contrasts. `"random"` computes average
  prediction error variance. Active only when `eval_efficiency = TRUE`.

- prediction_type:

  Character scalar. Random-effect model for efficiency evaluation.
  `"none"` skips. `"IID"` assumes independent entry effects. `"GBLUP"`
  and `"PBLUP"` use the relationship matrix `K`. Active only when
  `eval_efficiency = TRUE` and `treatment_effect = "random"`.

- K:

  Optional numeric matrix. Relationship matrix serving two purposes: (1)
  random-effect efficiency when
  `prediction_type %in% c("GBLUP", "PBLUP")`; (2) dispersion
  optimization when `use_dispersion = TRUE` and
  `dispersion_source = "K"`. Ignored otherwise. Row and column names
  must match treatment IDs or be reachable through `line_id_map`.

- line_id_map:

  Optional data frame with columns `Treatment` and `LineID`. Required
  only when `K` is active and treatment labels differ from
  `rownames(K)`. Serves the same role for `K` that `id_map` serves for
  `GRM` and `A`.

- varcomp:

  Named list of variance components for efficiency evaluation. Must
  contain `sigma_e2`, `sigma_g2`, `sigma_b2`, `sigma_r2`, and
  `sigma_c2`. All values except `sigma_g2` must be strictly positive.
  `sigma_g2` must be positive when `treatment_effect = "random"`. Active
  only when `eval_efficiency = TRUE`. The ratio of `sigma_b2` to
  `sigma_e2` has the largest effect on efficiency metrics for block
  designs.

- check_as_fixed:

  Logical. If `TRUE`, check treatments are included as fixed indicator
  columns in the mixed-model coefficient matrix during efficiency
  evaluation. Active only when `eval_efficiency = TRUE`.

- residual_structure:

  Character scalar. Residual covariance structure for efficiency
  evaluation. `"IID"` assumes independent residuals. `"AR1"` introduces
  row-wise autocorrelation controlled by `rho_row`. `"AR1xAR1"`
  introduces both row-wise and column-wise autocorrelation controlled by
  `rho_row` and `rho_col`. Active only when `eval_efficiency = TRUE`.

- rho_row:

  Numeric. AR1 autocorrelation parameter in the row direction. Must
  satisfy \\\|\rho\_{\text{row}}\| \< 1\\. Active when
  `eval_efficiency = TRUE` and
  `residual_structure %in% c("AR1", "AR1xAR1")`.

- rho_col:

  Numeric. AR1 autocorrelation parameter in the column direction. Must
  satisfy \\\|\rho\_{\text{col}}\| \< 1\\. Active when
  `eval_efficiency = TRUE` and `residual_structure = "AR1xAR1"`.

- spatial_engine:

  Character scalar. Computational strategy for the residual precision
  matrix when `residual_structure != "IID"`. `"auto"` selects `"dense"`
  when the number of observed plots is at most `dense_max_n` and
  `"sparse"` otherwise. `"dense"` and `"sparse"` force the respective
  implementation. Active only when `eval_efficiency = TRUE` and
  `residual_structure != "IID"`.

- dense_max_n:

  Integer. Plot count threshold used by `spatial_engine = "auto"` to
  switch between dense and sparse computation. Active only when
  `eval_efficiency = TRUE`, `residual_structure != "IID"`, and
  `spatial_engine = "auto"`.

- eff_trace_samples:

  Integer. Number of Hutchinson trace samples used when the target
  treatment dimension exceeds `eff_full_max` and exact efficiency
  extraction is replaced by stochastic approximation. Larger values
  reduce approximation variance at the cost of runtime. Active only when
  `eval_efficiency = TRUE`.

- eff_full_max:

  Integer. Maximum treatment dimension for which exact efficiency
  extraction via linear solve is attempted. Above this threshold the
  function switches to stochastic trace estimation. Active only when
  `eval_efficiency = TRUE`.

- check_placement:

  Character scalar. Strategy for positioning checks within blocks.
  `"random"` draws check positions freely within each block.
  `"systematic"` spaces checks approximately evenly in the 1D block
  ordering. `"optimal"` evaluates `check_opt_attempts` candidate layouts
  and selects the one maximizing the mean minimum Euclidean distance
  between check positions in the 2D grid.

- check_opt_attempts:

  Positive integer. Number of candidate layouts evaluated when
  `check_placement = "optimal"`. Larger values increase the chance of
  finding a spatially well-dispersed check arrangement at higher runtime
  cost. Ignored unless `check_placement = "optimal"`.

- use_dispersion:

  Logical. If `TRUE`, a swap-based local search is applied after layout
  construction to reduce spatial concentration of genetically similar
  non-check entries. The pairwise relatedness matrix is selected by
  `dispersion_source`. If `FALSE`, all dispersion arguments are ignored.

- dispersion_source:

  Character scalar. Relationship matrix used for pairwise relatedness
  scoring in dispersion optimization. `"K"` uses `K` (with optional
  `line_id_map`), `"A"` uses `A`, `"GRM"` uses `GRM`. Active only when
  `use_dispersion = TRUE`.

- dispersion_radius:

  Positive integer. Chebyshev neighborhood radius for dispersion
  scoring. Two plots are neighbors when \\\max(\|\Delta r\|, \|\Delta
  c\|) \leq \texttt{dispersion\\radius}\\. `1` considers immediate
  neighbors only; larger values extend the local window over which
  relatedness is penalized. Active only when `use_dispersion = TRUE`.

- dispersion_iters:

  Non-negative integer. Number of swap proposals in the dispersion
  search. Larger values allow more extensive search at higher runtime
  cost. Active only when `use_dispersion = TRUE`.

- dispersion_seed:

  Optional integer. Seed for the dispersion swap search, applied
  independently of `seed`. When `NULL`, `seed_used` is applied. Use when
  the initial layout and the dispersion refinement should be
  reproducible independently.

## Value

A named list with the following components:

- `layout_matrix`:

  Character matrix of dimension `n_rows x n_cols`. Each cell contains a
  treatment ID or `NA` for surplus unused positions.

- `field_book`:

  Data frame with one row per assigned plot. Columns: `Treatment`,
  `Family`, `Gcluster`, `Block`, `Plot`, `Row`, `Column`. `Gcluster` is
  `NA` for check treatments and for all treatments when
  `cluster_source = "Family"`.

- `efficiency`:

  `NULL` if `eval_efficiency = FALSE`; otherwise a named list of
  efficiency metrics including `mode`, `A`, `D`, `mean_PEV` or
  `mean_VarDiff`, `n_trt` or `n_lines`, `residual_structure_used`,
  `spatial_engine_used`, and `notes`.

- `seed_used`:

  Integer. The random seed used internally, returned for reproducibility
  records.

## Details

`prep_famoptg()` builds a block field design in which check treatments
appear in every block and non-check treatments are allocated across
blocks according to specified replication levels. The same function
produces three design classes depending on the supplied treatment
structure: an augmented repeated-check design when all non-check
treatments are unreplicated, a partially replicated (p-rep) design when
some non-check treatments are replicated and others are not, and an
RCBD-type repeated-check design when all non-check treatments receive a
common replication greater than 1. Optionally, the function applies
post-layout dispersion optimization to reduce local relatedness among
non-check entries, and evaluates the resulting design under a
mixed-model framework.

Serpentine traversal uses `mod()` from the **pracma** package. When the
adjusted field size exceeds the plot count, surplus cells appear as `NA`
in `layout_matrix`; the `field_book` contains only assigned plots.

The function separates construction logic, grouping logic, dispersion
logic, and efficiency logic. This means that, for example, family labels
can drive adjacency control while a genomic matrix drives dispersion
optimization, and efficiency evaluation can use a third, dedicated
relationship matrix `K`. None of these layers need to use the same
source.

## Design classes

**Augmented repeated-check design**

Obtained when `p_rep_treatments` is empty or `NULL` and all non-check
entries are supplied through `unreplicated_treatments`. Checks are
repeated in every block; all test entries appear once. Appropriate for
early-stage screening when the number of candidates exceeds what can be
replicated.

**Partially replicated (p-rep) repeated-check design**

Obtained when some non-check treatments are in `p_rep_treatments` with
replication greater than 1 and others remain in
`unreplicated_treatments`. Checks are repeated in every block; selected
candidates are replicated while others receive a single plot.

**RCBD-type repeated-check design**

Obtained when all non-check treatments are in `p_rep_treatments` with a
common replication greater than 1 and `unreplicated_treatments` is
empty. If the common replication equals `n_blocks`, every non-check
treatment appears once in every block alongside the checks, producing
the closest analogue to a classical RCBD within this repeated-check
framework. If the common replication is less than `n_blocks`, balance is
maintained but the design is no longer a strict classical RCBD.

## Construction sequence

1.  Inputs are validated and field dimensions are reconciled with the
    total plot count required.

2.  Group labels are derived from family labels or from PCA-based
    clustering of `GRM` or `A`.

3.  Replicated non-check entries are assigned to blocks, subject to the
    constraint that no treatment appears twice in the same block.

4.  Unreplicated entries are distributed across blocks.

5.  Checks are inserted into each block using the chosen placement
    strategy.

6.  Within each block, treatment order is shuffled to reduce same-group
    adjacency in the 1D block sequence.

7.  The ordered treatment sequence is mapped to the field grid according
    to `order` and `serpentine`.

8.  Optionally, a swap-based dispersion search reduces local relatedness
    in the 2D layout.

9.  Optionally, mixed-model efficiency metrics are computed on the final
    design.

## Grouping and why it matters

Group labels serve two purposes. During block construction, they define
which entries are considered similar for the purpose of adjacency
control in the 1D block ordering. During optional dispersion
optimization, a relationship matrix is used to score pairwise
relatedness between neighboring plots in the 2D grid and a swap search
reduces that score.

Three grouping modes are available. `"Family"` reads labels directly
from the supplied family vectors. `"GRM"` derives cluster labels from
the eigenstructure of the genomic relationship matrix via PCA followed
by k-means or hierarchical clustering. `"A"` applies the same procedure
to the pedigree relationship matrix. Matrix-based grouping is preferable
when family labels are too coarse to capture the relevant genetic
structure, or when dispersing genomically similar materials more
precisely is a priority.

## Argument dependencies

**Grouping mode**

When `cluster_source = "Family"`, labels come from `check_families`,
`p_rep_families`, and `unreplicated_families`; `GRM`, `A`, `id_map`,
`cluster_method`, `cluster_seed`, `cluster_attempts`, and `n_pcs_use`
are ignored. When `cluster_source = "GRM"`, `GRM` is required; `id_map`
is needed only if treatment IDs differ from `rownames(GRM)`. When
`cluster_source = "A"`, the same applies with `A`. In both matrix-based
modes, `cluster_method`, `cluster_seed`, `cluster_attempts`, and
`n_pcs_use` are active.

**Efficiency evaluation**

When `eval_efficiency = FALSE`, all efficiency arguments are ignored.
When `eval_efficiency = TRUE` and `treatment_effect = "fixed"`,
precision of fixed-effect treatment contrasts is computed and
`prediction_type`, `K`, and `line_id_map` are ignored. When
`treatment_effect = "random"`, `prediction_type` becomes active:
`"none"` skips random-effect efficiency; `"IID"` treats entry effects as
independent; `"GBLUP"` and `"PBLUP"` require `K`, with `line_id_map`
needed when treatment IDs differ from `rownames(K)`.

**Residual structure**

`"IID"` ignores `rho_row` and `rho_col`. `"AR1"` uses `rho_row` only.
`"AR1xAR1"` uses both. All residual arguments are inactive when
`eval_efficiency = FALSE`.

**Check placement**

`"optimal"` activates `check_opt_attempts` and evaluates multiple
candidate layouts, selecting the one that maximizes the mean minimum
distance between check positions in the 2D grid. `"random"` and
`"systematic"` ignore `check_opt_attempts`.

**Dispersion optimization**

When `use_dispersion = FALSE`, all dispersion arguments are ignored.
When `use_dispersion = TRUE`, `dispersion_source` selects the matrix
used for pairwise relatedness scoring: `"K"` requires `K` (with optional
`line_id_map`), `"A"` requires `A`, `"GRM"` requires `GRM`.

**Field dimension correction**

When `warn_and_correct = FALSE`, the function stops if `n_rows * n_cols`
does not match the total required plot count. When
`warn_and_correct = TRUE`, one dimension is adjusted upward to
accommodate all plots; `fix_rows` determines which dimension is held
fixed.

## Examples

``` r
data("OptiSparseMET_example_data", package = "OptiSparseMET")
x <- OptiSparseMET_example_data

env_trt <- x$treatments[1:24]
env_fam <- x$treatment_info$Family[match(env_trt, x$treatment_info$Treatment)]

## Example 1: p-rep repeated-check design — first 8 entries replicated
## twice, remaining 16 unreplicated. Checks in every block.
out_prep <- prep_famoptg(
  check_treatments         = c("CHK1", "CHK2"),
  check_families           = c("CHECK", "CHECK"),
  p_rep_treatments         = env_trt[1:8],
  p_rep_reps               = rep(2L, 8),
  p_rep_families           = env_fam[1:8],
  unreplicated_treatments  = env_trt[9:24],
  unreplicated_families    = env_fam[9:24],
  n_blocks                 = 4,
  n_rows                   = 6,
  n_cols                   = 8,
  order                    = "row",
  serpentine               = TRUE,
  seed                     = 123,
  cluster_source           = "Family",
  eval_efficiency          = FALSE,
  use_dispersion           = FALSE
)
#> Warning: Field size (6 x 8 = 48) does not match required (40). Adjusting dimensions.

dim(out_prep$layout_matrix)   # 6 x 8
#> [1] 6 7
head(out_prep$field_book)
#>        Treatment Family Gcluster Block Plot Row Column
#> plot_1      L008    F06     <NA>     1    1   1      1
#> plot_2      L012    F14     <NA>     1    2   1      2
#> plot_3      CHK2  CHECK     <NA>     1    3   1      3
#> plot_4      L015    F10     <NA>     1    4   1      4
#> plot_5      L002    F15     <NA>     1    5   1      5
#> plot_6      L007    F02     <NA>     1    6   1      6
out_prep$seed_used
#> [1] 123

if (FALSE) { # \dontrun{
## Example 2: augmented repeated-check design — no replicated entries,
## all 24 test entries unreplicated.
out_aug <- prep_famoptg(
  check_treatments         = c("CHK1", "CHK2"),
  check_families           = c("CHECK", "CHECK"),
  p_rep_treatments         = character(0),
  p_rep_reps               = integer(0),
  p_rep_families           = character(0),
  unreplicated_treatments  = env_trt,
  unreplicated_families    = env_fam,
  n_blocks                 = 4,
  n_rows                   = 5,
  n_cols                   = 6,
  order                    = "row",
  serpentine               = FALSE,
  seed                     = 123,
  cluster_source           = "Family",
  eval_efficiency          = FALSE,
  use_dispersion           = FALSE
)

dim(out_aug$layout_matrix)
head(out_aug$field_book)

## Example 3: RCBD-type repeated-check design — all 24 entries replicated
## 4 times (once per block). Equivalent to a classical RCBD with repeated
## checks.
out_rcbd <- prep_famoptg(
  check_treatments         = c("CHK1", "CHK2"),
  check_families           = c("CHECK", "CHECK"),
  p_rep_treatments         = env_trt,
  p_rep_reps               = rep(4L, length(env_trt)),
  p_rep_families           = env_fam,
  unreplicated_treatments  = character(0),
  unreplicated_families    = character(0),
  n_blocks                 = 4,
  n_rows                   = 14,
  n_cols                   = 8,
  order                    = "row",
  serpentine               = TRUE,
  seed                     = 123,
  cluster_source           = "Family",
  eval_efficiency          = FALSE,
  use_dispersion           = FALSE
)

dim(out_rcbd$layout_matrix)
head(out_rcbd$field_book)
} # }
```
