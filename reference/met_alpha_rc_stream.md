# Construct a stream-based repeated-check alpha row-column design

`met_alpha_rc_stream()` is the OptiSparseMET version of
`alpha_rc_stream()` from the OptiDesign package. The `met_` prefix
avoids namespace conflicts when both packages are loaded simultaneously.
All arguments, return values, and internal logic are identical to
`alpha_rc_stream()` in OptiDesign.

`met_alpha_rc_stream()` builds a randomised alpha-lattice incomplete
block design for breeding and agronomic experiments on a fixed field
grid of `n_rows x n_cols` plots. The field is converted to a
one-dimensional planting stream according to `order` and `serpentine`,
partitioned into `n_reps` contiguous replicate segments, and each
replicate is further divided into incomplete blocks. Every incomplete
block contains all check treatments plus a subset of entry treatments.
Entry treatments appear exactly once per replicate. Trailing unused
field cells remain `NA` in the layout matrix and field book.

**Design evaluation has been separated from construction.** This
function returns the field book, layout matrix, and design metadata
only. To compute A-, D-, and CDmean optimality criteria, call
[`met_evaluate_alpha_efficiency()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_evaluate_alpha_efficiency.md)
on the returned `field_book`. To search for an optimised design using
Random Restart, Simulated Annealing, or a Genetic Algorithm, call
[`met_optimize_alpha_rc()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_optimize_alpha_rc.md).

## Usage

``` r
met_alpha_rc_stream(
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
  n_blocks_per_rep = NULL,
  min_block_size = NULL,
  max_block_size = NULL,
  check_placement = c("systematic", "random"),
  check_position_pattern = c("spread", "corners_first"),
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
  every incomplete block of every replicate. Must not overlap with
  `entry_treatments` and must not contain duplicates.

- check_families:

  Character vector of the same length as `check_treatments`. Family
  labels for checks, used for adjacency scoring. Always stored in the
  `Family` column; `Gcluster` is always `NA` for check plots regardless
  of `cluster_source`.

- entry_treatments:

  Character vector of non-check treatment identifiers. Each entry
  appears exactly once per replicate. Must not contain duplicates.

- entry_families:

  Character vector of the same length as `entry_treatments`. Family
  labels for entries. Used for adjacency scoring when
  `cluster_source = "Family"` and to set the number of clusters when
  `cluster_source %in% c("GRM", "A")`.

- n_reps:

  Positive integer. Number of contiguous replicate segments.

- n_rows:

  Positive integer. Number of field rows.

- n_cols:

  Positive integer. Number of field columns.

- order:

  Character. Stream traversal direction: `"row"` fills row by row;
  `"column"` fills column by column.

- serpentine:

  Logical. If `TRUE`, alternate rows (when `order = "row"`) or alternate
  columns (when `order = "column"`) traverse in reverse direction.

- seed:

  Optional integer. If `NULL`, a seed is drawn from
  `1:.Machine$integer.max` and returned as `seed_used`. Controls all
  randomised steps including entry randomisation, within-block
  shuffling, check placement, and dispersion optimisation.

- attempts:

  Positive integer. Number of within-block shuffling attempts to
  minimise adjacent same-group pairs. Exits early when score reaches
  zero.

- warn_and_correct:

  Logical. Retained for interface compatibility. Field capacity
  violations always raise an error.

- fix_rows:

  Logical. Retained for interface compatibility. Field geometry is
  always fixed.

- cluster_source:

  Character. Grouping source for adjacency scoring. `"Family"` uses
  `entry_families` directly. `"GRM"` and `"A"` derive clusters from PCA
  of the respective relationship matrix.

- GRM:

  Optional square numeric matrix with rownames and colnames equal to
  line IDs. Required when `cluster_source = "GRM"` or when
  `use_dispersion = TRUE` and `dispersion_source = "GRM"`.

- A:

  Optional square numeric matrix. Required when `cluster_source = "A"`
  or `dispersion_source = "A"`.

- id_map:

  Optional data frame with columns `Treatment` and `LineID`. Required
  when treatment labels do not match relationship matrix rownames.

- cluster_method:

  Character. Clustering algorithm: `"kmeans"` uses
  [`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html) with
  `cluster_attempts` random starts; `"hclust"` uses Ward's D2 linkage.

- cluster_seed:

  Integer. Seed for k-means initialisation, run in an isolated RNG
  scope.

- cluster_attempts:

  Positive integer. Number of k-means random restarts. Ignored when
  `cluster_method = "hclust"`.

- n_pcs_use:

  Positive number or `Inf`. Number of leading PCs for matrix-based
  clustering. Actual number used is
  `min(n_pcs_use, n_positive_eigenvalues - 1)`.

- n_blocks_per_rep:

  Optional positive integer. If supplied, the exact number of incomplete
  blocks per replicate. Validated against the feasible range derived
  from `min_block_size`, `max_block_size`, and field capacity. If
  `NULL`, the block count is set automatically to \\b\_\text{max}\\.

- min_block_size:

  Optional positive integer \\\geq c\\. Minimum **total** block size
  (checks + entries). Sets an upper bound on \\b\\.

- max_block_size:

  Optional positive integer \\\geq c\\. Maximum **total** block size
  (checks + entries). Sets a lower bound on \\b\\.

- check_placement:

  Character. Check position rule within blocks. `"systematic"`
  distributes checks evenly; `"random"` draws positions uniformly
  without replacement.

- check_position_pattern:

  Character. Retained for interface compatibility.

- use_dispersion:

  Logical. If `TRUE`, apply a post-layout swap-based local dispersion
  optimisation to reduce pairwise relatedness among neighbouring
  non-check plots. Requires a relationship matrix specified by
  `dispersion_source`. Default `FALSE`.

- dispersion_source:

  Character. Matrix used for dispersion scoring when
  `use_dispersion = TRUE`: `"K"`, `"A"`, or `"GRM"`. The corresponding
  argument must be non-`NULL`.

- dispersion_radius:

  Positive integer. Chebyshev distance radius defining the neighbourhood
  for dispersion scoring.

- dispersion_iters:

  Non-negative integer. Number of random swap proposals for dispersion
  optimisation.

- dispersion_seed:

  Integer. Seed for the dispersion step, run in an isolated RNG scope.
  Defaults to `seed_used` when `NULL`.

- K:

  Optional square numeric matrix. Used for dispersion scoring when
  `dispersion_source = "K"`.

- line_id_map:

  Optional data frame with columns `Treatment` and `LineID`. Required
  when treatment labels do not match `rownames(K)` or the dispersion
  matrix.

- verbose:

  Logical. If `TRUE`, prints a one-line summary of field usage, trailing
  `NA` plots, replicate sizes, block count and origin, block sizes in
  replicate 1, and active block-size bounds.

## Value

A named list with four components:

- `layout_matrix`:

  Character matrix of dimension `n_rows x n_cols`. Used cells contain
  treatment IDs; unused trailing cells are `NA`.

- `field_book`:

  Data frame with one row per field plot. Columns: `Plot` (stream
  position), `Row`, `Column`, `Rep` (`NA` for unused), `IBlock` (global
  block index), `BlockInRep` (block index within rep), `Treatment` (`NA`
  for unused), `Family`, `Gcluster` (`NA` when
  `cluster_source = "Family"` or for check plots), `Check` (logical).

- `design_info`:

  Named list: `n_rows`, `n_cols`, `total_plots`, `n_reps`, `n_checks`,
  `n_entries`, `n_blocks_per_rep` (resolved), `n_blocks_per_rep_user`,
  `min_block_size`, `max_block_size`, `min_entry_slots_per_block`,
  `max_entry_slots_per_block`, `rep_sizes`, `total_used_plots`,
  `trailing_na_plots`, `block_plan`, `block_meta`, `order`,
  `serpentine`.

- `seed_used`:

  Integer. The random seed used.

## Details

### Design structure

Let:

- \\r\\ = number of replicates (`n_reps`)

- \\c\\ = number of checks (`length(check_treatments)`)

- \\v\\ = number of entries (`length(entry_treatments)`)

- \\b\\ = number of incomplete blocks per replicate

Each replicate contains all \\v\\ entries exactly once and all \\c\\
checks repeated once in each of the \\b\\ incomplete blocks. Entries are
split across blocks as evenly as possible - block entry counts differ by
at most one. The number of used plots per replicate is:

\$\$\text{plots per rep} = v + b \times c\$\$

The total used plots in the design is:

\$\$\text{total used} = r \times (v + b \times c)\$\$

If the fixed field contains more cells than this total, unused trailing
cells remain `NA`. The function stops if the required plots exceed field
capacity.

### Block-size parameterisation

Block-size constraints are expressed in terms of **total block size**
(checks + entries) because in a repeated-check design the whole block is
the operational unit. The mapping is:

\$\$\text{min\\entry\\slots} = \text{min\\block\\size} - c\$\$
\$\$\text{max\\entry\\slots} = \text{max\\block\\size} - c\$\$

For example, with \\c = 3\\ checks, `min_block_size = 19` and
`max_block_size = 20` means each block holds between 16 and 17 entries
plus 3 check positions. Both bounds must be \\\geq c\\.

### Block-count determination

After translating total-block-size limits into entry-slot limits, the
number of incomplete blocks per replicate \\b\\ must satisfy a feasible
range \\\[b\_\text{min},\\ b\_\text{max}\]\\:

- Lower bound from `max_block_size`:

  \\b \geq \lceil v / (\text{max\\block\\size} - c)\rceil\\

- Upper bound from `min_block_size`:

  \\b \leq \lfloor v / (\text{min\\block\\size} - c)\rfloor\\

- Upper bound from field capacity:

  \\b \leq \lfloor (\text{total\\plots}/r - v) / c \rfloor\\

The function stops with an informative error if \\b\_\text{min} \>
b\_\text{max}\\. At least one of `n_blocks_per_rep`, `min_block_size`,
or `max_block_size` must be supplied.

**User-fixed mode** (`n_blocks_per_rep` supplied): the value is
validated against the feasible range and all block-size constraints.

**Automatic mode** (`n_blocks_per_rep = NULL`): \\b\\ is set to
\\b\_\text{max}\\ - the largest feasible number of blocks, producing the
smallest feasible blocks, which generally gives the best local control.

### Field stream construction

The field is mapped to a one-dimensional stream by iterating over rows
(`order = "row"`) or columns (`order = "column"`). When
`serpentine = TRUE`, alternate rows or columns reverse direction,
producing a boustrophedon planting path. Replicate segments occupy
contiguous portions of the stream and incomplete blocks occupy
contiguous sub-segments within each replicate.

### Within-block arrangement

Entries are randomised independently within each replicate, then split
across blocks. Checks are inserted at positions determined by
`check_placement`:

- `"systematic"`: positions computed as
  `round(seq(1, block_size, length.out = n_checks + 2))[2:(n_checks + 1)]`,
  distributing checks as evenly as possible.

- `"random"`: positions drawn uniformly without replacement.

Within-block order is optimised over `attempts` random shuffles to
minimise the count of adjacent same-group treatment pairs. The search
exits early when the adjacency score reaches zero.

### Grouping and adjacency scoring

The grouping label used for adjacency scoring is determined by
`cluster_source`:

- `"Family"`: uses `entry_families` and `check_families` directly.
  `Gcluster` is `NA` for all plots.

- `"GRM"` or `"A"`: derives clusters from the leading principal
  components of the relationship matrix. The number of clusters equals
  `length(unique(entry_families))`. Cluster labels are prefixed `"G"`
  (GRM) or `"A"` (A matrix). Check treatments always use
  `check_families` and always have `Gcluster = NA`.

### Dispersion optimisation

When `use_dispersion = TRUE`, a local swap search is applied after
layout construction. At each of `dispersion_iters` iterations, two
non-check plots are selected at random; the swap is accepted if it
reduces the total neighbourhood relatedness score - the sum of
relationship matrix values between all non-check neighbour pairs within
Chebyshev distance `dispersion_radius`. The matrix used is selected by
`dispersion_source`.

## See also

[`met_evaluate_alpha_efficiency()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_evaluate_alpha_efficiency.md)
to compute A, D, and CDmean optimality criteria on the returned
`field_book`.
[`met_optimize_alpha_rc()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_optimize_alpha_rc.md)
to search for a criterion-optimal design using Random Restart, Simulated
Annealing, or a Genetic Algorithm.

## Examples

``` r
## Basic construction: 167 entries, 3 checks, 3 reps
design <- met_alpha_rc_stream(
  check_treatments = c("CHK1", "CHK2", "CHK3"),
  check_families   = c("CHECK", "CHECK", "CHECK"),
  entry_treatments = paste0("G", 1:167),
  entry_families   = rep(paste0("F", 1:7), length.out = 167),
  n_reps           = 3,
  n_rows           = 30,
  n_cols           = 20,
  min_block_size   = 19,
  max_block_size   = 20
)
#> Fixed field = 30 x 20; used plots = 591; trailing NA plots = 9; replicate used sizes = {197, 197, 197}; blocks/rep = 10 (derived); block sizes in rep 1 = {20, 20, 20, 20, 20, 20, 20, 19, 19, 19}; min block size = 19; max block size = 20

dim(design$layout_matrix)             # 30 x 20
#> [1] 30 20
head(design$field_book)
#>   Plot Row Column Rep IBlock BlockInRep Treatment Family Gcluster Check
#> 1    1   1      1   1      1          1      G159     F5     <NA> FALSE
#> 2    2   2      1   1      1          1       G65     F2     <NA> FALSE
#> 3    3   3      1   1      1          1       G26     F5     <NA> FALSE
#> 4    4   4      1   1      1          1       G14     F7     <NA> FALSE
#> 5    5   5      1   1      1          1       G29     F1     <NA> FALSE
#> 6    6   6      1   1      1          1      CHK1  CHECK     <NA>  TRUE
design$design_info$n_blocks_per_rep
#> [1] 10
design$design_info$trailing_na_plots
#> [1] 9

## Evaluate separately
eff <- met_evaluate_alpha_efficiency(
  field_book         = design$field_book,
  n_rows             = 30,
  n_cols             = 20,
  check_treatments   = c("CHK1", "CHK2", "CHK3"),
  treatment_effect   = "fixed",
  residual_structure = "AR1xAR1",
  rho_row = 0.10, rho_col = 0.10
)
eff$A_criterion
#> [1] 0.7456078
eff$D_criterion
#> [1] 0.3528107
```
