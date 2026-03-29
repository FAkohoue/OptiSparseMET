# Combine environment-level field books into a single MET field book

Environment-level field books produced by
[`met_prep_famoptg()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_prep_famoptg.md)
or
[`met_alpha_rc_stream()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_alpha_rc_stream.md)
may differ in their column sets depending on the design options used in
each environment. A direct
[`rbind()`](https://rdrr.io/r/base/cbind.html) over such a list fails
unless all data frames share identical columns.
`combine_met_fieldbooks()` resolves this by computing the union of all
column names across environments, padding missing columns with `NA` in
each data frame before stacking, and placing a standard set of columns
at the front of the result regardless of the order in which they appear
in the individual field books.

The following columns are guaranteed to be present in the output even
when absent in the input, filled with `NA` where the source field book
did not contain them: `Treatment`, `Family`, `Gcluster`, `Block`,
`Plot`, `Row`, `Column`. All other columns from any input field book are
preserved and appear after the standard column set.

MET-level metadata is added to every row before stacking:

- `Environment`: taken from the name of the list element.

- `LocalDesign`: the design type used in that environment, from
  `local_designs` if supplied, otherwise `NA`.

- `ReplicationMode`: the replication mode used in that environment, from
  `replication_modes` if supplied, otherwise `NA`.

- `SparseMethod`: the across-environment allocation strategy, from
  `sparse_method` if supplied, otherwise `NA`.

- `IsCommonTreatment`: logical flag indicating whether the treatment in
  each row appears in `common_treatments`.

## Usage

``` r
combine_met_fieldbooks(
  field_books,
  local_designs = NULL,
  replication_modes = NULL,
  sparse_method = NULL,
  common_treatments = NULL
)
```

## Arguments

- field_books:

  Named list of data frames, one per environment. List names are used as
  the values of the `Environment` column in the output and must
  therefore be unique and non-empty. Each element must be a data frame;
  the function stops with an error if any element is not.

- local_designs:

  Optional named character vector or named list. Gives the local design
  type used in each environment, typically `"met_prep_famoptg"` or
  `"met_alpha_rc_stream"` for built-in engines. Names must match the
  names of `field_books`. Environments not present in `local_designs`
  receive `NA` in the `LocalDesign` column. If `NULL`, all environments
  receive `NA`.

- replication_modes:

  Optional named character vector or named list. Gives the replication
  mode used in each environment, corresponding to the `replication_mode`
  argument passed to
  [`assign_replication_by_seed()`](https://FAkohoue.github.io/OptiSparseMET/reference/assign_replication_by_seed.md)
  for that environment. Names must match the names of `field_books`.
  Environments not present in `replication_modes` receive `NA` in the
  `ReplicationMode` column. If `NULL`, all environments receive `NA`.

- sparse_method:

  Optional character scalar. The across-environment allocation strategy
  used in
  [`allocate_sparse_met()`](https://FAkohoue.github.io/OptiSparseMET/reference/allocate_sparse_met.md),
  e.g. `"balanced_incomplete"` or `"random_balanced"`. Applied uniformly
  to all rows of the combined field book. If `NULL`, the `SparseMethod`
  column is filled with `NA`.

- common_treatments:

  Optional character vector of treatment IDs forced into all
  environments. Used to populate the `IsCommonTreatment` logical column
  by matching each row's `Treatment` value against this vector. If
  `NULL`, `IsCommonTreatment` is `FALSE` for all rows.

## Value

A data frame with one row per plot across all environments. Columns
appear in the following order: `Environment`, `LocalDesign`,
`ReplicationMode`, `SparseMethod`, `IsCommonTreatment`, `Treatment`,
`Family`, `Gcluster`, `Block`, `Plot`, `Row`, `Column`, followed by any
additional columns present in the input field books in the order they
are encountered. Row names are reset to `NULL`. Columns absent in a
given environment's field book are filled with `NA` for all rows
belonging to that environment.

## Details

`combine_met_fieldbooks()` stacks a named list of environment-level
field books into one combined multi-environment trial (MET) field book.
Each environment's field book is augmented with MET-level metadata
columns before stacking, and columns that are present in some
environments but absent in others are filled with `NA` rather than
causing an error. The result is a single flat data frame suitable for
export as a field book or as input to downstream mixed-model analysis
pipelines.

## See also

[`met_prep_famoptg()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_prep_famoptg.md)
and
[`met_alpha_rc_stream()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_alpha_rc_stream.md)
which produce the environment-level field books that serve as input to
this function.
[`plan_sparse_met_design()`](https://FAkohoue.github.io/OptiSparseMET/reference/plan_sparse_met_design.md)
which calls this function internally as the final step of the two-stage
MET pipeline.

## Examples

``` r
## Minimal example: two environments with the same columns.
fb_E1 <- data.frame(
  Treatment = c("L001", "L002", "CHK1"),
  Family    = c("F1", "F2", "CHECK"),
  Block     = c(1L, 1L, 1L),
  Plot      = 1:3,
  Row       = c(1L, 1L, 1L),
  Column    = 1:3,
  stringsAsFactors = FALSE
)

fb_E2 <- data.frame(
  Treatment = c("L003", "L004", "CHK1"),
  Family    = c("F3", "F4", "CHECK"),
  Block     = c(1L, 1L, 1L),
  Plot      = 1:3,
  Row       = c(1L, 1L, 1L),
  Column    = 1:3,
  stringsAsFactors = FALSE
)

met <- combine_met_fieldbooks(
  field_books        = list(E1 = fb_E1, E2 = fb_E2),
  local_designs      = c(E1 = "met_prep_famoptg", E2 = "met_prep_famoptg"),
  replication_modes  = c(E1 = "augmented", E2 = "augmented"),
  sparse_method      = "balanced_incomplete",
  common_treatments  = "CHK1"
)

nrow(met)                        # 6 -- three plots per environment
#> [1] 6
unique(met$Environment)          # "E1" "E2"
#> [1] "E1" "E2"
met$IsCommonTreatment            # TRUE only for CHK1 rows
#> [1] FALSE FALSE  TRUE FALSE FALSE  TRUE
head(met[, 1:8])
#>   Environment      LocalDesign ReplicationMode        SparseMethod
#> 1          E1 met_prep_famoptg       augmented balanced_incomplete
#> 2          E1 met_prep_famoptg       augmented balanced_incomplete
#> 3          E1 met_prep_famoptg       augmented balanced_incomplete
#> 4          E2 met_prep_famoptg       augmented balanced_incomplete
#> 5          E2 met_prep_famoptg       augmented balanced_incomplete
#> 6          E2 met_prep_famoptg       augmented balanced_incomplete
#>   IsCommonTreatment Treatment Family Gcluster
#> 1             FALSE      L001     F1       NA
#> 2             FALSE      L002     F2       NA
#> 3              TRUE      CHK1  CHECK       NA
#> 4             FALSE      L003     F3       NA
#> 5             FALSE      L004     F4       NA
#> 6              TRUE      CHK1  CHECK       NA

## Heterogeneous columns: E2 has an extra column absent in E1.
## combine_met_fieldbooks() fills the missing column with NA for E1 rows.
fb_E2$SpatialResidual <- rnorm(3)

met2 <- combine_met_fieldbooks(
  field_books = list(E1 = fb_E1, E2 = fb_E2)
)

"SpatialResidual" %in% names(met2)                         # TRUE
#> [1] TRUE
is.na(met2$SpatialResidual[met2$Environment == "E1"])      # all TRUE
#> [1] TRUE TRUE TRUE
```
