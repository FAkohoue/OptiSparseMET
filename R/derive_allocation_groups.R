#' Derive allocation group labels for sparse MET treatment assignment
#'
#' `derive_allocation_groups()` assigns a group label to each treatment prior
#' to sparse allocation across environments. These labels are then used by
#' [allocate_sparse_met()] to guide the incidence structure so that genetic
#' groups -- defined by family membership or by clusters derived from a
#' relationship matrix -- are distributed across environments rather than
#' concentrated in a subset of them. The function is called internally by
#' [allocate_sparse_met()] when `allocation_group_source` is not `"none"`,
#' but can also be called directly to inspect or audit the group structure
#' before running allocation.
#'
#' @description
#' Four grouping modes are supported. `"none"` assigns all treatments to a
#' single group labelled `"ALL"`, which disables group-guided allocation
#' without requiring any change to the allocation function call. `"Family"`
#' reads group labels directly from `treatment_info$Family`, one label per
#' treatment. `"GRM"` and `"A"` derive cluster labels from the eigenstructure
#' of the genomic or pedigree relationship matrix respectively, using PCA
#' followed by k-means or hierarchical clustering.
#'
#' When the number of clusters is not determined by the user directly, it is
#' anchored to the number of distinct family labels among the supplied
#' treatments if `treatment_info` is available, or otherwise approximated as
#' \eqn{\max(2,\, \lfloor\sqrt{n}\rfloor)} where \eqn{n} is the number of
#' treatments.
#'
#' @details
#' ## Matrix-based grouping
#'
#' When `allocation_group_source %in% c("GRM", "A")`, the function extracts
#' the treatment-level submatrix, performs eigendecomposition, and retains
#' eigenvectors corresponding to positive eigenvalues (threshold
#' \eqn{> 10^{-10}}) as principal components. Component scores are scaled by
#' the square root of the corresponding eigenvalues before clustering, which
#' weights components proportionally to their contribution to variance in the
#' relationship matrix.
#'
#' The number of components retained is controlled by `n_pcs_use`. Setting
#' `n_pcs_use = Inf` retains all positive-eigenvalue components up to
#' \eqn{\min(n_{\text{pos}},\, n - 1)}. Smaller values retain only the leading
#' components, preserving broad genetic structure at the cost of finer
#' differentiation. At least 2 informative components must be available;
#' the function stops with an error if this condition is not met.
#'
#' K-means clustering uses `group_attempts` random restarts seeded by
#' `group_seed`. Hierarchical clustering uses Ward's minimum variance criterion
#' (`method = "ward.D2"`) and is not affected by `group_seed` or
#' `group_attempts`. Resulting cluster labels are prefixed `"GRP_G"` for GRM
#' clusters and `"GRP_A"` for pedigree clusters.
#'
#' ## ID matching
#'
#' By default, treatment IDs in `treatments` are matched directly to row names
#' of the relationship matrix. When field-book treatment labels differ from
#' matrix row names, supply `id_map` with columns `Treatment` and `LineID`;
#' the function uses this map to resolve the correspondence before extracting
#' the submatrix.
#'
#' @param treatments Character vector of treatment IDs. Duplicate values are
#'   silently deduplicated. Must contain at least one element.
#'
#' @param allocation_group_source Character scalar. Grouping mode. One of:
#' \describe{
#'   \item{`"none"`}{All treatments assigned to a single group `"ALL"`.}
#'   \item{`"Family"`}{Group labels read from `treatment_info$Family`.}
#'   \item{`"GRM"`}{Cluster labels derived from `GRM` via PCA and clustering.}
#'   \item{`"A"`}{Cluster labels derived from `A` via PCA and clustering.}
#' }
#'
#' @param treatment_info Optional data frame. Required when
#'   `allocation_group_source = "Family"`. Must contain columns `Treatment`
#'   and `Family`. When `allocation_group_source %in% c("GRM", "A")`, this
#'   argument is optional but, if supplied with a `Family` column, is used to
#'   anchor the number of clusters to the number of distinct families among
#'   the supplied treatments. The function stops if any treatment in
#'   `treatments` is absent from `treatment_info$Treatment` when
#'   `allocation_group_source = "Family"`.
#'
#' @param GRM Optional numeric matrix. Genomic relationship matrix. Required
#'   when `allocation_group_source = "GRM"`. Must be square with row and
#'   column names. Row names must match treatment IDs in `treatments` or be
#'   reachable through `id_map`.
#'
#' @param A Optional numeric matrix. Pedigree-based numerator relationship
#'   matrix. Required when `allocation_group_source = "A"`. Same structural
#'   requirements as `GRM`.
#'
#' @param id_map Optional data frame with columns `Treatment` and `LineID`.
#'   Required only when treatment IDs in `treatments` do not match the row
#'   names of `GRM` or `A`. The function uses `LineID` to look up the
#'   corresponding matrix rows. Ignored when
#'   `allocation_group_source %in% c("none", "Family")`.
#'
#' @param group_method Character scalar. Clustering algorithm applied to the
#'   PCA scores. `"kmeans"` uses k-means with `group_attempts` random
#'   restarts. `"hclust"` uses Ward's criterion hierarchical clustering.
#'   Ignored when `allocation_group_source %in% c("none", "Family")`.
#'
#' @param group_seed Integer. Random seed passed to k-means initialization.
#'   Active only when `allocation_group_source %in% c("GRM", "A")` and
#'   `group_method = "kmeans"`. Has no effect on hierarchical clustering.
#'
#' @param group_attempts Integer. Number of random restarts for k-means.
#'   Larger values reduce the risk of converging to a poor local optimum.
#'   Active only when `allocation_group_source %in% c("GRM", "A")` and
#'   `group_method = "kmeans"`.
#'
#' @param n_pcs_use Integer or `Inf`. Number of leading principal components
#'   retained for clustering. `Inf` retains all components corresponding to
#'   positive eigenvalues, up to \eqn{n - 1}. Smaller integer values retain
#'   only the leading components. Must be at least 2. Ignored when
#'   `allocation_group_source %in% c("none", "Family")`.
#'
#' @return A data frame with one row per element of `treatments` (after
#'   deduplication) and the following columns:
#' \describe{
#'   \item{`Treatment`}{Character. Treatment ID, in the order they appear in
#'     `treatments` after deduplication.}
#'   \item{`AllocationGroup`}{Character. Derived group label. `"ALL"` under
#'     `"none"`; the family label string under `"Family"`; a string of the
#'     form `"GRP_G{k}"` under `"GRM"` or `"GRP_A{k}"` under `"A"`, where
#'     `{k}` is the integer cluster index.}
#' }
#'
#' @seealso [allocate_sparse_met()] which calls this function internally when
#'   `allocation_group_source` is not `"none"`. Call `derive_allocation_groups()`
#'   directly to inspect or audit the group structure before running allocation.
#'
#' @examples
#' treatments <- paste0("L", sprintf("%03d", 1:12))
#'
#' treatment_info <- data.frame(
#'   Treatment = treatments,
#'   Family    = rep(c("F1", "F2", "F3"), each = 4),
#'   stringsAsFactors = FALSE
#' )
#'
#' ## Example 1: family-based groups
#' grp_fam <- derive_allocation_groups(
#'   treatments              = treatments,
#'   allocation_group_source = "Family",
#'   treatment_info          = treatment_info
#' )
#' grp_fam
#' # AllocationGroup is "F1", "F2", or "F3"
#'
#' ## Example 2: no grouping
#' grp_none <- derive_allocation_groups(
#'   treatments              = treatments,
#'   allocation_group_source = "none"
#' )
#' unique(grp_none$AllocationGroup)  # "ALL"
#'
#' ## Example 3: GRM-based clustering
#' set.seed(1)
#' n   <- length(treatments)
#' raw <- matrix(rnorm(n * n), n, n)
#' GRM <- crossprod(raw) / n
#' diag(GRM) <- diag(GRM) + 0.1
#' rownames(GRM) <- colnames(GRM) <- treatments
#'
#' grp_grm <- derive_allocation_groups(
#'   treatments              = treatments,
#'   allocation_group_source = "GRM",
#'   GRM                     = GRM,
#'   treatment_info          = treatment_info,
#'   group_method            = "kmeans",
#'   group_seed              = 42,
#'   group_attempts          = 25,
#'   n_pcs_use               = Inf
#' )
#' grp_grm
#' # AllocationGroup values are "GRP_G1", "GRP_G2", "GRP_G3"
#'
#' @importFrom stats kmeans hclust cutree dist setNames
#' @export
derive_allocation_groups <- function(
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
) {

  allocation_group_source <- match.arg(allocation_group_source)
  group_method            <- match.arg(group_method)

  if (missing(treatments) || length(treatments) < 1)
    stop("`treatments` must contain at least one treatment ID.")

  treatments <- unique(as.character(treatments))

  if (!(is.numeric(n_pcs_use) && length(n_pcs_use) == 1 &&
        (is.finite(n_pcs_use) || is.infinite(n_pcs_use)) && n_pcs_use > 0))
    stop("`n_pcs_use` must be a single positive number or Inf.")

  # ---- 1. No grouping ----
  if (allocation_group_source == "none") {
    return(data.frame(
      Treatment       = treatments,
      AllocationGroup = rep("ALL", length(treatments)),
      stringsAsFactors = FALSE
    ))
  }

  # ---- 2. Family-based grouping ----
  if (allocation_group_source == "Family") {
    if (is.null(treatment_info) ||
        !is.data.frame(treatment_info) ||
        !all(c("Treatment", "Family") %in% names(treatment_info)))
      stop(
        "When `allocation_group_source = 'Family'`, `treatment_info` must be a data.frame ",
        "with columns `Treatment` and `Family`."
      )

    treatment_info$Treatment <- as.character(treatment_info$Treatment)
    treatment_info$Family    <- as.character(treatment_info$Family)

    fam <- treatment_info$Family[match(treatments, treatment_info$Treatment)]

    if (any(is.na(fam))) {
      miss <- treatments[is.na(fam)]
      stop("Some treatments have no family in `treatment_info`: ",
           paste(utils::head(miss, 10), collapse = ", "))
    }

    return(data.frame(
      Treatment       = treatments,
      AllocationGroup = fam,
      stringsAsFactors = FALSE
    ))
  }

  # ---- 3. Matrix-based grouping (GRM or A) ----
  Kc <- if (allocation_group_source == "GRM") GRM else A

  if (is.null(Kc))
    stop("Relationship matrix required for `allocation_group_source = '",
         allocation_group_source, "'.")
  if (is.null(rownames(Kc)) || is.null(colnames(Kc)))
    stop("Relationship matrix must have row names and column names.")

  if (is.null(id_map)) {
    line_ids <- stats::setNames(treatments, treatments)
  } else {
    if (!is.data.frame(id_map) || !all(c("Treatment", "LineID") %in% names(id_map)))
      stop("`id_map` must be a data.frame with columns `Treatment` and `LineID`.")
    line_ids <- stats::setNames(as.character(id_map$LineID), as.character(id_map$Treatment))
  }

  treatment_line_ids <- unname(line_ids[treatments])
  missing_ids        <- setdiff(treatment_line_ids, rownames(Kc))
  if (length(missing_ids) > 0)
    stop("Some treatment IDs are missing from the selected relationship matrix.")

  Ksub <- Kc[treatment_line_ids, treatment_line_ids, drop = FALSE]
  eg   <- eigen(Ksub, symmetric = TRUE)
  pos  <- which(eg$values > 1e-10)

  if (length(pos) < 2)
    stop("The selected relationship matrix has too few positive eigenvalues for clustering.")

  # Determine number of clusters
  if (!is.null(treatment_info) && is.data.frame(treatment_info) &&
      all(c("Treatment", "Family") %in% names(treatment_info))) {

    treatment_info$Treatment <- as.character(treatment_info$Treatment)
    treatment_info$Family    <- as.character(treatment_info$Family)
    fam <- treatment_info$Family[match(treatments, treatment_info$Treatment)]
    fam <- fam[!is.na(fam)]
    k_groups <- if (length(fam) > 0) max(2L, length(unique(fam))) else
      max(2L, round(sqrt(length(treatments))))
  } else {
    k_groups <- max(2L, round(sqrt(length(treatments))))
  }

  max_possible_pcs <- min(length(pos), nrow(Ksub) - 1L)
  n_pcs <- if (is.infinite(n_pcs_use)) max_possible_pcs else
    min(as.integer(n_pcs_use), max_possible_pcs)

  if (n_pcs < 2)
    stop("Fewer than 2 informative PCs are available for clustering.")

  pcs <- eg$vectors[, pos[seq_len(n_pcs)], drop = FALSE]
  pcs <- sweep(pcs, 2, sqrt(eg$values[pos[seq_len(n_pcs)]]), `*`)

  if (group_method == "kmeans") {
    set.seed(group_seed)
    cl <- stats::kmeans(pcs, centers = k_groups, nstart = group_attempts)$cluster
  } else {
    d  <- stats::dist(pcs)
    hc <- stats::hclust(d, method = "ward.D2")
    cl <- stats::cutree(hc, k = k_groups)
  }

  prefix <- if (allocation_group_source == "GRM") "GRP_G" else "GRP_A"

  data.frame(
    Treatment       = treatments,
    AllocationGroup = paste0(prefix, cl),
    stringsAsFactors = FALSE
  )
}
