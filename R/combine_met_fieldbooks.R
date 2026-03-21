#' Combine environment-level field books into a single MET field book
#'
#' `combine_met_fieldbooks()` stacks a named list of environment-level field
#' books into one combined multi-environment trial (MET) field book. Each
#' environment's field book is augmented with MET-level metadata columns before
#' stacking, and columns that are present in some environments but absent in
#' others are filled with `NA` rather than causing an error. The result is a
#' single flat data frame suitable for export as a field book or as input to
#' downstream mixed-model analysis pipelines.
#'
#' @description
#' Environment-level field books produced by `prep_famoptg()` or
#' `alpha_rc_stream()` may differ in their column sets depending on the design
#' options used in each environment. A direct `rbind()` over such a list fails
#' unless all data frames share identical columns. `combine_met_fieldbooks()`
#' resolves this by computing the union of all column names across environments,
#' padding missing columns with `NA` in each data frame before stacking, and
#' placing a standard set of columns at the front of the result regardless of
#' the order in which they appear in the individual field books.
#'
#' The following columns are guaranteed to be present in the output even when
#' absent in the input, filled with `NA` where the source field book did not
#' contain them: `Treatment`, `Family`, `Gcluster`, `Block`, `Plot`, `Row`,
#' `Column`. All other columns from any input field book are preserved and
#' appear after the standard column set.
#'
#' MET-level metadata is added to every row before stacking:
#'
#' - `Environment`: taken from the name of the list element.
#' - `LocalDesign`: the design type used in that environment, from
#'   `local_designs` if supplied, otherwise `NA`.
#' - `ReplicationMode`: the replication mode used in that environment, from
#'   `replication_modes` if supplied, otherwise `NA`.
#' - `SparseMethod`: the across-environment allocation strategy, from
#'   `sparse_method` if supplied, otherwise `NA`.
#' - `IsCommonTreatment`: logical flag indicating whether the treatment in each
#'   row appears in `common_treatments`.
#'
#' @param field_books Named list of data frames, one per environment. List
#'   names are used as the values of the `Environment` column in the output
#'   and must therefore be unique and non-empty. Each element must be a data
#'   frame; the function stops with an error if any element is not.
#'
#' @param local_designs Optional named character vector or named list. Gives
#'   the local design type used in each environment, e.g. `"augmented"`,
#'   `"p_rep"`, or `"alpha_rc"`. Names must match the names of `field_books`.
#'   Environments not present in `local_designs` receive `NA` in the
#'   `LocalDesign` column. If `NULL`, all environments receive `NA`.
#'
#' @param replication_modes Optional named character vector or named list.
#'   Gives the replication mode used in each environment, corresponding to the
#'   `replication_mode` argument passed to `assign_replication_by_seed()` for
#'   that environment. Names must match the names of `field_books`.
#'   Environments not present in `replication_modes` receive `NA` in the
#'   `ReplicationMode` column. If `NULL`, all environments receive `NA`.
#'
#' @param sparse_method Optional character scalar. The across-environment
#'   allocation strategy used in `allocate_sparse_met()`, e.g.
#'   `"balanced_incomplete"` or `"random_balanced"`. Applied uniformly to all
#'   rows of the combined field book. If `NULL`, the `SparseMethod` column is
#'   filled with `NA`.
#'
#' @param common_treatments Optional character vector of treatment IDs forced
#'   into all environments. Used to populate the `IsCommonTreatment` logical
#'   column by matching each row's `Treatment` value against this vector. If
#'   `NULL`, `IsCommonTreatment` is `FALSE` for all rows.
#'
#' @return A data frame with one row per plot across all environments. Columns
#'   appear in the following order: `Environment`, `LocalDesign`,
#'   `ReplicationMode`, `SparseMethod`, `IsCommonTreatment`, `Treatment`,
#'   `Family`, `Gcluster`, `Block`, `Plot`, `Row`, `Column`, followed by any
#'   additional columns present in the input field books in the order they are
#'   encountered. Row names are reset to `NULL`. Columns absent in a given
#'   environment's field book are filled with `NA` for all rows belonging to
#'   that environment.
#'
#' @examples
#' ## Minimal example: two environments with the same columns.
#' fb_E1 <- data.frame(
#'   Treatment = c("L001", "L002", "CHK1"),
#'   Family    = c("F1", "F2", "CHECK"),
#'   Block     = c(1L, 1L, 1L),
#'   Plot      = 1:3,
#'   Row       = c(1L, 1L, 1L),
#'   Column    = 1:3,
#'   stringsAsFactors = FALSE
#' )
#'
#' fb_E2 <- data.frame(
#'   Treatment = c("L003", "L004", "CHK1"),
#'   Family    = c("F3", "F4", "CHECK"),
#'   Block     = c(1L, 1L, 1L),
#'   Plot      = 1:3,
#'   Row       = c(1L, 1L, 1L),
#'   Column    = 1:3,
#'   stringsAsFactors = FALSE
#' )
#'
#' met <- combine_met_fieldbooks(
#'   field_books        = list(E1 = fb_E1, E2 = fb_E2),
#'   local_designs      = c(E1 = "augmented", E2 = "augmented"),
#'   replication_modes  = c(E1 = "augmented", E2 = "augmented"),
#'   sparse_method      = "balanced_incomplete",
#'   common_treatments  = "CHK1"
#' )
#'
#' nrow(met)                        # 6 — three plots per environment
#' unique(met$Environment)          # "E1" "E2"
#' met$IsCommonTreatment            # TRUE only for CHK1 rows
#' head(met[, 1:8])
#'
#' ## Heterogeneous columns: E2 has an extra column absent in E1.
#' ## combine_met_fieldbooks() fills the missing column with NA for E1 rows.
#' fb_E2$SpatialResidual <- rnorm(3)
#'
#' met2 <- combine_met_fieldbooks(
#'   field_books = list(E1 = fb_E1, E2 = fb_E2)
#' )
#'
#' "SpatialResidual" %in% names(met2)         # TRUE
#' is.na(met2$SpatialResidual[met2$Environment == "E1"])  # all TRUE
#'
#' @export
#' 
#' 
combine_met_fieldbooks <- function(
    field_books,
    local_designs = NULL,
    replication_modes = NULL,
    sparse_method = NULL,
    common_treatments = NULL
) {
  
  if (!is.list(field_books) || length(field_books) < 1 || is.null(names(field_books))) {
    stop("`field_books` must be a named list of data frames.")
  }
  
  envs <- names(field_books)
  out_list <- vector("list", length(field_books))
  names(out_list) <- envs
  
  for (env in envs) {
    fb <- field_books[[env]]
    
    if (!is.data.frame(fb)) {
      stop("Each element of `field_books` must be a data frame.")
    }
    
    # Ensure core columns exist
    core_cols <- c("Treatment", "Family", "Gcluster", "Block", "Plot", "Row", "Column")
    for (cc in core_cols) {
      if (!cc %in% names(fb)) {
        fb[[cc]] <- NA
      }
    }
    
    # Add MET-level metadata
    fb$Environment <- env
    
    fb$LocalDesign <- if (!is.null(local_designs)) {
      unname(local_designs[[env]])
    } else {
      NA_character_
    }
    
    fb$ReplicationMode <- if (!is.null(replication_modes)) {
      unname(replication_modes[[env]])
    } else {
      NA_character_
    }
    
    fb$SparseMethod <- if (is.null(sparse_method)) NA_character_ else sparse_method
    fb$IsCommonTreatment <- if (is.null(common_treatments)) FALSE else fb$Treatment %in% common_treatments
    
    # Put standard columns first
    front_cols <- c(
      "Environment", "LocalDesign", "ReplicationMode", "SparseMethod",
      "IsCommonTreatment", "Treatment", "Family", "Gcluster",
      "Block", "Plot", "Row", "Column"
    )
    
    remain_cols <- setdiff(names(fb), front_cols)
    fb <- fb[, c(front_cols, remain_cols), drop = FALSE]
    
    out_list[[env]] <- fb
  }
  
  # Union of all column names across environments
  all_cols <- unique(unlist(lapply(out_list, names)))
  
  # Add missing columns to each field book
  out_list <- lapply(out_list, function(df) {
    missing_cols <- setdiff(all_cols, names(df))
    if (length(missing_cols) > 0) {
      for (mc in missing_cols) {
        df[[mc]] <- NA
      }
    }
    df[, all_cols, drop = FALSE]
  })
  
  out <- do.call(rbind, out_list)
  rownames(out) <- NULL
  out
}