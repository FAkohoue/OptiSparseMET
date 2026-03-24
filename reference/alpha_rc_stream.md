# Construct a fixed-grid alpha row-column design using stream-based incomplete blocking

The function is designed for field situations where the physical grid is
fixed before design construction, planting or harvesting follows a
continuous field-book order rather than rectangular replicate
boundaries, and repeated checks in every incomplete block are required.
Replicate boundaries are determined by position in the traversal stream,
not by geometric subdivision of the grid. This means replicate 2 begins
exactly where replicate 1 ends in stream order, and all surplus cells
accumulate at the tail of the stream as unused positions.

Given those constraints, the function proceeds as follows:

1.  Generates the global field stream from `n_rows`, `n_cols`, `order`,
    and `serpentine`.

2.  Determines the number of incomplete blocks supportable per replicate
    given entry count, check count, and `min_entry_slots_per_block`.

3.  Partitions each replicate segment into incomplete blocks.

4.  Assigns checks and entries to blocks subject to the capacity plan.

5.  Improves entry arrangement heuristically to reduce local
    concentration of genetically or structurally similar material.

6.  Places surplus cells as trailing `NA` at the end of the stream.

7.  Optionally refines placement via dispersion optimization.

8.  Optionally evaluates the resulting design under a mixed-model
    framework.

## Usage

``` r
alpha_rc_stream(
  check_treatments,
  check_families,
  entry_treatments,
  entry_families,
  n_reps,
  n_rows,
  n_cols,
  order = "column",
  serpentine = FALSE,
  seed = NULL,
  attempts = 5000,
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
  min_entry_slots_per_block = 8,
  max_blocks_per_rep = NULL,
  eval_efficiency = FALSE,
  treatment_effect = c("random", "fixed"),
  prediction_type = c("none", "IID", "GBLUP", "PBLUP"),
  K = NULL,
  line_id_map = NULL,
  varcomp = list(sigma_e2 = 1, sigma_g2 = 1, sigma_rep2 = 1, sigma_ib2 = 1, sigma_r2 = 1,
    sigma_c2 = 1),
  check_as_fixed = TRUE,
  residual_structure = c("IID", "AR1", "AR1xAR1"),
  rho_row = 0,
  rho_col = 0,
  spatial_engine = c("auto", "sparse", "dense"),
  dense_max_n = 5000,
  eff_trace_samples = 80,
  eff_full_max = 400,
  check_placement = c("systematic", "random"),
  check_position_pattern = c("spread", "corners_first"),
  use_dispersion = FALSE,
  dispersion_source = c("K", "A", "GRM"),
  dispersion_radius = 1,
  dispersion_iters = 2000,
  dispersion_seed = 1,
  verbose = TRUE
)
```

## Arguments

- check_treatments:

  Character vector of check IDs. Every check is placed exactly once in
  every incomplete block. Checks define the within-block benchmark
  structure across both blocks and replicates, and their count directly
  affects the remaining entry capacity per block.

- check_families:

  Character vector of the same length as `check_treatments`. Family or
  group labels for checks. Used directly when
  `cluster_source = "Family"` and stored in the output field book under
  all settings. Must align exactly with `check_treatments`.

- entry_treatments:

  Character vector of test entry IDs. Each entry appears exactly once
  per replicate. Must not overlap with `check_treatments`.

- entry_families:

  Character vector of the same length as `entry_treatments`. Family or
  group labels for entries. Used directly when
  `cluster_source = "Family"`. When `cluster_source %in% c("GRM", "A")`,
  these labels inform the target number of clusters but do not define
  them.

- n_reps:

  Integer. Number of replicates. Each replicate is a contiguous segment
  of the global stream of length `E + b * C`. Increasing `n_reps`
  increases the total number of used plots and reduces the stream length
  available per replicate, which may reduce the number of supportable
  incomplete blocks.

- n_rows:

  Integer. Number of field rows. Fixed before design construction and
  never altered by this function. Together with `n_cols`, determines
  total field capacity and the global traversal stream.

- n_cols:

  Integer. Number of field columns. Fixed before design construction and
  never altered by this function.

- order:

  Character. Global traversal direction. `"row"` generates a row-major
  stream; `"column"` generates a column-major stream. This determines
  where replicate and incomplete block boundaries fall and where
  trailing `NA` cells are placed. Should reflect the operational
  movement direction in the field.

- serpentine:

  Logical. If `TRUE`, alternate rows (when `order = "row"`) or alternate
  columns (when `order = "column"`) reverse traversal direction,
  producing a boustrophedon stream. Use when field-book order follows
  serpentine movement. Changes stream order and consequently replicate
  and incomplete block boundaries.

- seed:

  Optional integer. Random seed for reproducibility. Controls entry
  allocation, within-block arrangement, random check placement, and
  dispersion optimization unless a separate `dispersion_seed` is
  supplied. If `NULL`, a seed is generated internally and returned in
  `seed_used`.

- attempts:

  Integer. Number of swap proposals used when improving entry allocation
  across blocks within each replicate. Larger values increase search
  effort without changing the underlying algorithm. Relevant mainly when
  family structure is strongly imbalanced or when the number of entries
  per replicate is large.

- warn_and_correct:

  Logical. Retained for interface continuity. Field dimensions are fixed
  in this function and are not altered regardless of this setting.

- fix_rows:

  Logical. Retained for interface continuity. Field dimensions are fixed
  in this function.

- cluster_source:

  Character. Grouping source for entry arrangement. `"Family"` uses
  `entry_families` and `check_families` directly. `"GRM"` derives
  clusters from the leading principal components of the genomic
  relationship matrix. `"A"` derives clusters from the pedigree
  relationship matrix. Determines which of `GRM`, `A`, `id_map`,
  `cluster_method`, `cluster_seed`, `cluster_attempts`, and `n_pcs_use`
  become active.

- GRM:

  Optional numeric matrix. Genomic relationship matrix. Required when
  `cluster_source = "GRM"` or when `use_dispersion = TRUE` and
  `dispersion_source = "GRM"`. Row and column names must match entry IDs
  or be reachable through `id_map`.

- A:

  Optional numeric matrix. Pedigree-based numerator relationship matrix.
  Required when `cluster_source = "A"` or when `use_dispersion = TRUE`
  and `dispersion_source = "A"`.

- id_map:

  Optional data frame with columns `Treatment` and `LineID`. Required
  only when `cluster_source %in% c("GRM", "A")` and treatment IDs in the
  field book do not match the row/column names of the relationship
  matrix.

- cluster_method:

  Character. Clustering method applied to the principal components of
  `GRM` or `A`. `"kmeans"` or `"hclust"`. Ignored when
  `cluster_source = "Family"`.

- cluster_seed:

  Integer. Seed for k-means initialization. Active only when
  `cluster_source %in% c("GRM", "A")` and `cluster_method = "kmeans"`.

- cluster_attempts:

  Integer. Number of random restarts for k-means. Active only when
  `cluster_source %in% c("GRM", "A")` and `cluster_method = "kmeans"`.

- n_pcs_use:

  Integer or `Inf`. Number of leading principal components retained for
  matrix-based clustering. `Inf` retains all components. Ignored when
  `cluster_source = "Family"`. Smaller values retain only broad
  structure; larger values preserve finer genetic differentiation.

- min_entry_slots_per_block:

  Integer. Minimum number of entry slots permitted in any incomplete
  block after checks are placed. Constrains the maximum number of blocks
  per replicate: if adding another block would reduce entry slots below
  this threshold, the block count is capped. Particularly important when
  check count is large relative to block capacity.

- max_blocks_per_rep:

  Optional integer. Upper bound on the number of incomplete blocks per
  replicate. If `NULL`, block count is derived solely from capacity and
  `min_entry_slots_per_block`. Use when operational or analytical
  considerations favor fewer, larger blocks.

- eval_efficiency:

  Logical. If `TRUE`, mixed-model efficiency diagnostics are computed on
  the non-`NA` plots of the resulting design. If `FALSE`, all
  efficiency-related arguments are ignored and `efficiency` is returned
  as `NULL`.

- treatment_effect:

  Character. Specifies how entries enter the efficiency model. `"fixed"`
  computes precision of pairwise treatment contrasts. `"random"`
  computes average prediction error variance of BLUP-type predictors.
  Active only when `eval_efficiency = TRUE`.

- prediction_type:

  Character. Random-effect model for efficiency evaluation. `"none"`
  skips random-effect evaluation. `"IID"` assumes independent entry
  effects. `"GBLUP"` and `"PBLUP"` use the relationship matrix `K`.
  Active only when `eval_efficiency = TRUE` and
  `treatment_effect = "random"`.

- K:

  Optional numeric matrix. Relationship matrix used for random-effect
  efficiency (`prediction_type %in% c("GBLUP", "PBLUP")`) or for
  dispersion optimization (`use_dispersion = TRUE` and
  `dispersion_source = "K"`). Ignored otherwise.

- line_id_map:

  Optional data frame with columns `Treatment` and `LineID`. Required
  only when `K` is active and treatment labels in the field book differ
  from `rownames(K)`. Serves the same purpose for `K` that `id_map`
  serves for `GRM` and `A`.

- varcomp:

  Named list of variance components used in efficiency evaluation. Must
  contain `sigma_e2`, `sigma_g2`, `sigma_rep2`, `sigma_ib2`, `sigma_r2`,
  and `sigma_c2`. Active only when `eval_efficiency = TRUE`. Use values
  consistent with the anticipated trial model; the relative magnitude of
  `sigma_ib2` to `sigma_e2` has the largest effect on efficiency metrics
  for incomplete block designs.

- check_as_fixed:

  Logical. If `TRUE`, checks are included as fixed indicators in the
  efficiency model. Active only when `eval_efficiency = TRUE`.

- residual_structure:

  Character. Residual covariance structure assumed in the efficiency
  model. `"IID"` ignores spatial correlation. `"AR1"` introduces
  row-wise autocorrelation controlled by `rho_row`. `"AR1xAR1"`
  introduces both row-wise and column-wise autocorrelation controlled by
  `rho_row` and `rho_col`. Active only when `eval_efficiency = TRUE`.

- rho_row:

  Numeric. AR1 autocorrelation parameter in the row direction. Active
  when `eval_efficiency = TRUE` and
  `residual_structure %in% c("AR1", "AR1xAR1")`.

- rho_col:

  Numeric. AR1 autocorrelation parameter in the column direction. Active
  when `eval_efficiency = TRUE` and `residual_structure = "AR1xAR1"`.

- spatial_engine:

  Character. Retained for interface compatibility. Efficiency
  computations use sparse Matrix operations internally regardless of
  this setting. Accepted values are `"auto"`, `"sparse"`, and `"dense"`.

- dense_max_n:

  Integer. Retained for interface compatibility.

- eff_trace_samples:

  Integer. Number of Hutchinson trace samples used when the number of
  target treatments exceeds `eff_full_max` and exact efficiency
  extraction is replaced by an approximation. Larger values reduce
  approximation variance at the cost of runtime.

- eff_full_max:

  Integer. Maximum treatment dimension for which exact efficiency
  extraction is attempted. Above this threshold the function switches to
  stochastic trace estimation.

- check_placement:

  Character. Method for positioning checks within blocks. `"systematic"`
  spaces checks approximately evenly across the block stream. `"random"`
  draws check positions freely. Does not interact with
  `check_position_pattern`.

- check_position_pattern:

  Retained for interface compatibility. Not used by the stream-based
  placement logic.

- use_dispersion:

  Logical. If `TRUE`, a post-hoc swap search is applied to reduce local
  concentration of genetically similar entries within replicates, using
  pairwise relatedness from the matrix selected by `dispersion_source`.
  If `FALSE`, all dispersion arguments are ignored.

- dispersion_source:

  Character. Relationship matrix used for dispersion scoring. `"K"` uses
  `K` (with optional `line_id_map`), `"A"` uses `A`, `"GRM"` uses `GRM`.
  Active only when `use_dispersion = TRUE`.

- dispersion_radius:

  Integer. Neighborhood radius for dispersion scoring. Two plots are
  neighbors when \\\max(\|\Delta r\|, \|\Delta c\|) \leq
  \texttt{dispersion\\radius}\\. Use `1` for immediate neighbors only;
  larger values extend the local window over which relatedness is
  penalized.

- dispersion_iters:

  Integer. Number of swap proposals in the dispersion search. Larger
  values increase optimization effort and runtime.

- dispersion_seed:

  Optional integer. Seed for the dispersion swap search, applied
  independently of `seed`. Useful when the initial design and the
  dispersion refinement should be reproducible separately.

- verbose:

  Logical. If `TRUE`, prints the derived replicate segment lengths and
  incomplete block structure to the console. Use for diagnostics when
  the automatic block plan needs to be verified.

## Value

A named list with the following components:

- `layout_matrix`:

  Character matrix of dimension `n_rows x n_cols`. Each cell contains a
  treatment ID or `NA` for unused positions.

- `field_book`:

  Data frame with one row per field cell, including unused cells.
  Columns include row, column, replicate, incomplete block, treatment
  ID, entry or check status, and family label.

- `efficiency`:

  `NULL` if `eval_efficiency = FALSE`; otherwise a named list of
  efficiency metrics computed from non-`NA` plots under the specified
  mixed model.

- `seed_used`:

  Integer. The random seed used internally, returned for reproducibility
  records.

- `design_info`:

  Named list summarizing replicate segment lengths, incomplete block
  sizes, number of used and unused cells, and the key settings that
  determined the block plan.

## Details

`alpha_rc_stream()` builds a row-column field design over a fixed
`n_rows x n_cols` grid. The full grid is traversed in a single ordered
stream, partitioned into contiguous replicate segments, and each segment
is further subdivided into incomplete blocks of possibly unequal size.
Within every replicate, each entry appears exactly once. Within every
incomplete block, every check appears exactly once. Field cells that
cannot be assigned to any replicate are left as trailing `NA` positions
at the end of the stream.

## Stream-based versus rectangle-based layout

The distinction from conventional row-column design functions is that
replicates are defined by stream position, not by rectangular spatial
partitions. This matters when operational field movement is continuous
(row-wise or column-wise), when replicate boundaries are logistical
rather than geometric, and when all unused cells should fall at the end
of the sequence rather than within replicate areas. The incomplete block
structure is consequently a function of stream segment length and
capacity constraints, not of field shape.

## Plot count identity

Let `E` be the number of entries, `C` the number of checks, `R` the
number of replicates, and `b` the number of incomplete blocks per
replicate. The number of used plots per replicate is:

\$\$E + b \times C\$\$

and the total number of used plots across the trial is:

\$\$R \times (E + b \times C)\$\$

This total must not exceed the fixed field capacity `n_rows x n_cols`.
The function chooses `b` to be as large as possible subject to
`min_entry_slots_per_block`, `max_blocks_per_rep` if supplied, and the
capacity constraint. Efficiency evaluation and dispersion optimization
operate only on non-`NA` plots.

## Argument dependencies

**Grouping source**

When `cluster_source = "Family"`, grouping comes directly from
`check_families` and `entry_families`; `GRM`, `A`, `id_map`,
`cluster_method`, `cluster_seed`, `cluster_attempts`, and `n_pcs_use`
are all ignored. When `cluster_source = "GRM"`, `GRM` is required and
`id_map` is needed only if entry names differ from `rownames(GRM)`;
grouping is derived by PCA followed by clustering. The same applies for
`cluster_source = "A"` with matrix `A`.

**Efficiency evaluation**

When `eval_efficiency = FALSE`, all efficiency arguments are ignored.
When `eval_efficiency = TRUE` and `treatment_effect = "fixed"`,
fixed-treatment contrast precision is computed and `prediction_type`,
`K`, and `line_id_map` are ignored. When `treatment_effect = "random"`,
`prediction_type` becomes active. `"none"` skips random-effect
efficiency; `"IID"` treats entry effects as independent and ignores `K`
and `line_id_map`; `"GBLUP"` and `"PBLUP"` require `K`, with
`line_id_map` needed if treatment IDs differ from `rownames(K)`.

**Residual structure**

When `residual_structure = "IID"`, `rho_row` and `rho_col` are ignored.
`"AR1"` uses only `rho_row`. `"AR1xAR1"` uses both `rho_row` and
`rho_col`. All residual structure arguments are inactive when
`eval_efficiency = FALSE`.

**Dispersion optimization**

When `use_dispersion = FALSE`, all dispersion arguments are ignored.
When `use_dispersion = TRUE`, `dispersion_source` selects the
relationship matrix used for pairwise similarity scoring: `"K"` requires
`K` (with optional `line_id_map`), `"A"` requires `A`, `"GRM"` requires
`GRM`.

**Interface-continuity arguments**

`warn_and_correct`, `fix_rows`, `spatial_engine`, `dense_max_n`, and
`check_position_pattern` are retained for interface consistency with
related functions. They do not alter the stream-based design logic.

## References

Montesinos-López, O. A., Mosqueda-González, B. A., Salinas-Ruiz, J.,
Montesinos-López, A., & Crossa, J. (2023). Sparse multi-trait genomic
prediction under balanced incomplete block design. *The Plant Genome*,
16, e20305.

## Examples

``` r
data("OptiSparseMET_example_data", package = "OptiSparseMET")
x <- OptiSparseMET_example_data

env_trt <- x$treatments[1:24]
env_fam <- x$treatment_info$Family[match(env_trt, x$treatment_info$Treatment)]

out <- alpha_rc_stream(
  check_treatments          = c("CHK1", "CHK2"),
  check_families            = c("CHECK", "CHECK"),
  entry_treatments          = env_trt,
  entry_families            = env_fam,
  n_reps                    = 2,
  n_rows                    = 8,
  n_cols                    = 8,
  order                     = "row",
  serpentine                = TRUE,
  seed                      = 123,
  cluster_source            = "Family",
  min_entry_slots_per_block = 6,
  eval_efficiency           = FALSE,
  check_placement           = "random",
  use_dispersion            = FALSE,
  verbose                   = FALSE
)

dim(out$layout_matrix)
#> [1] 8 8
head(out$field_book)
#>   Treatment Family Gcluster Check PlotStream Rep IBlock BlockInRep Row Column
#> 1      CHK1  CHECK     <NA>  TRUE          1   1      1          1   1      1
#> 2      L014    F09     <NA> FALSE          2   1      1          1   1      2
#> 3      L005    F03     <NA> FALSE          3   1      1          1   1      3
#> 4      L012    F14     <NA> FALSE          4   1      1          1   1      4
#> 5      L011    F04     <NA> FALSE          5   1      1          1   1      5
#> 6      CHK2  CHECK     <NA>  TRUE          6   1      1          1   1      6
out$design_info$n_blocks_per_rep
#> [1] 4
```
