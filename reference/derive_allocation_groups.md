# Derive allocation group labels for sparse MET treatment assignment

Four grouping modes are supported. `"none"` assigns all treatments to a
single group labelled `"ALL"`, which disables group-guided allocation
without requiring any change to the allocation function call. `"Family"`
reads group labels directly from `treatment_info$Family`, one label per
treatment. `"GRM"` and `"A"` derive cluster labels from the
eigenstructure of the genomic or pedigree relationship matrix
respectively, using PCA followed by k-means or hierarchical clustering.

When the number of clusters is not determined by the user directly, it
is anchored to the number of distinct family labels among the supplied
treatments if `treatment_info` is available, or otherwise approximated
as \\\max(2,\\ \lfloor\sqrt{n}\rfloor)\\ where \\n\\ is the number of
treatments.

## Usage

``` r
derive_allocation_groups(
  treatments,
  allocation_group_source = c("none", "Family", "GRM", "A"),
  treatment_info = NULL,
  GRM = NULL,
  A = NULL,
  id_map = NULL,
  group_method = c("kmeans", "hclust"),
  group_seed = 1,
  group_attempts = 25,
  n_pcs_use = Inf
)
```

## Arguments

- treatments:

  Character vector of treatment IDs. Duplicate values are silently
  deduplicated. Must contain at least one element.

- allocation_group_source:

  Character scalar. Grouping mode. One of:

  - `"none"`: all treatments assigned to a single group `"ALL"`.

  - `"Family"`: group labels read from `treatment_info$Family`.

  - `"GRM"`: cluster labels derived from `GRM` via PCA and clustering.

  - `"A"`: cluster labels derived from `A` via PCA and clustering.

- treatment_info:

  Optional data frame. Required when
  `allocation_group_source = "Family"`. Must contain columns `Treatment`
  and `Family`. When `allocation_group_source %in% c("GRM", "A")`, this
  argument is optional but, if supplied with a `Family` column, is used
  to anchor the number of clusters to the number of distinct families
  among the supplied treatments. The function stops if any treatment in
  `treatments` is absent from `treatment_info$Treatment` when
  `allocation_group_source = "Family"`.

- GRM:

  Optional numeric matrix. Genomic relationship matrix. Required when
  `allocation_group_source = "GRM"`. Must be square with row and column
  names. Row names must match treatment IDs in `treatments` or be
  reachable through `id_map`.

- A:

  Optional numeric matrix. Pedigree-based numerator relationship matrix.
  Required when `allocation_group_source = "A"`. Same structural
  requirements as `GRM`.

- id_map:

  Optional data frame with columns `Treatment` and `LineID`. Required
  only when treatment IDs in `treatments` do not match the row names of
  `GRM` or `A`. The function uses `LineID` to look up the corresponding
  matrix rows. Ignored when
  `allocation_group_source %in% c("none", "Family")`.

- group_method:

  Character scalar. Clustering algorithm applied to the PCA scores.
  `"kmeans"` uses k-means with `group_attempts` random restarts.
  `"hclust"` uses Ward's criterion hierarchical clustering. Ignored when
  `allocation_group_source %in% c("none", "Family")`.

- group_seed:

  Integer. Random seed passed to k-means initialization. Active only
  when `allocation_group_source %in% c("GRM", "A")` and
  `group_method = "kmeans"`. Has no effect on hierarchical clustering.

- group_attempts:

  Integer. Number of random restarts for k-means. Larger values reduce
  the risk of converging to a poor local optimum. Active only when
  `allocation_group_source %in% c("GRM", "A")` and
  `group_method = "kmeans"`.

- n_pcs_use:

  Integer or `Inf`. Number of leading principal components retained for
  clustering. `Inf` retains all components corresponding to positive
  eigenvalues, up to \\n - 1\\. Smaller integer values retain only the
  leading components, preserving broad structure and discarding finer
  differentiation. Must be at least 2. Ignored when
  `allocation_group_source %in% c("none", "Family")`.

## Value

A data frame with one row per element of `treatments` (after
deduplication) and the following columns:

- `Treatment`:

  Character. Treatment ID, in the order they appear in `treatments`
  after deduplication.

- `AllocationGroup`:

  Character. Derived group label. `"ALL"` under `"none"`; the family
  label string under `"Family"`; a string of the form `"GRP_G{k}"` under
  `"GRM"` or `"GRP_A{k}"` under `"A"`, where `{k}` is the integer
  cluster index.

## Details

`derive_allocation_groups()` assigns a group label to each treatment
prior to sparse allocation across environments. These labels are then
used by
[`allocate_sparse_met()`](https://FAkohoue.github.io/OptiSparseMET/reference/allocate_sparse_met.md)
to guide the incidence structure so that genetic groups — defined by
family membership or by clusters derived from a relationship matrix —
are distributed across environments rather than concentrated in a subset
of them. The function is called internally by
[`allocate_sparse_met()`](https://FAkohoue.github.io/OptiSparseMET/reference/allocate_sparse_met.md)
when `allocation_group_source` is not `"none"`, but can also be called
directly to inspect or audit the group structure before running
allocation.

### Matrix-based grouping

When `allocation_group_source %in% c("GRM", "A")`, the function extracts
the treatment-level submatrix, performs eigendecomposition, and retains
eigenvectors corresponding to positive eigenvalues (threshold \\\>
10^{-10}\\) as principal components. Component scores are scaled by the
square root of the corresponding eigenvalues before clustering, which
weights components proportionally to their contribution to variance in
the relationship matrix.

The number of components retained is controlled by `n_pcs_use`. Setting
`n_pcs_use = Inf` retains all positive-eigenvalue components up to
\\\min(n\_{\text{pos}},\\ n - 1)\\. Smaller values retain only the
leading components, preserving broad genetic structure at the cost of
finer differentiation. At least 2 informative components must be
available; the function stops with an error if this condition is not
met.

K-means clustering uses `group_attempts` random restarts seeded by
`group_seed`. Hierarchical clustering uses Ward's minimum variance
criterion (`method = "ward.D2"`) and is not affected by `group_seed` or
`group_attempts`. Resulting cluster labels are prefixed `"GRP_G"` for
GRM clusters and `"GRP_A"` for pedigree clusters.

### ID matching

By default, treatment IDs in `treatments` are matched directly to row
names of the relationship matrix. When field-book treatment labels
differ from matrix row names, supply `id_map` with columns `Treatment`
and `LineID`; the function uses this map to resolve the correspondence
before extracting the submatrix.

## Examples

``` r
treatments <- paste0("L", sprintf("%03d", 1:12))

treatment_info <- data.frame(
  Treatment = treatments,
  Family    = rep(c("F1", "F2", "F3"), each = 4),
  stringsAsFactors = FALSE
)

## Example 1: family-based groups — labels read directly from treatment_info.
grp_fam <- derive_allocation_groups(
  treatments              = treatments,
  allocation_group_source = "Family",
  treatment_info          = treatment_info
)
grp_fam
#>    Treatment AllocationGroup
#> 1       L001              F1
#> 2       L002              F1
#> 3       L003              F1
#> 4       L004              F1
#> 5       L005              F2
#> 6       L006              F2
#> 7       L007              F2
#> 8       L008              F2
#> 9       L009              F3
#> 10      L010              F3
#> 11      L011              F3
#> 12      L012              F3
# 12 rows; AllocationGroup is "F1", "F2", or "F3" according to family

## Example 2: no grouping — all treatments assigned to "ALL".
## Use this to disable group-guided allocation without changing any other
## argument in allocate_sparse_met().
grp_none <- derive_allocation_groups(
  treatments              = treatments,
  allocation_group_source = "none"
)
unique(grp_none$AllocationGroup)  # "ALL"
#> [1] "ALL"

## Example 3: GRM-based clustering — groups derived from eigenstructure
## of a genomic relationship matrix. Number of clusters anchored to the
## number of families in treatment_info when supplied.
set.seed(1)
n <- length(treatments)
raw <- matrix(rnorm(n * n), n, n)
GRM <- crossprod(raw) / n
diag(GRM) <- diag(GRM) + 0.1
rownames(GRM) <- colnames(GRM) <- treatments

grp_grm <- derive_allocation_groups(
  treatments              = treatments,
  allocation_group_source = "GRM",
  GRM                     = GRM,
  treatment_info          = treatment_info,
  group_method            = "kmeans",
  group_seed              = 42,
  group_attempts          = 25,
  n_pcs_use               = Inf
)
grp_grm
#>    Treatment AllocationGroup
#> 1       L001          GRP_G2
#> 2       L002          GRP_G3
#> 3       L003          GRP_G3
#> 4       L004          GRP_G2
#> 5       L005          GRP_G2
#> 6       L006          GRP_G1
#> 7       L007          GRP_G3
#> 8       L008          GRP_G1
#> 9       L009          GRP_G3
#> 10      L010          GRP_G2
#> 11      L011          GRP_G2
#> 12      L012          GRP_G3
# AllocationGroup values are "GRP_G1", "GRP_G2", "GRP_G3"
# (3 clusters, anchored to 3 families in treatment_info)
```
