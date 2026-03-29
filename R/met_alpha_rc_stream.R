# ==============================================================================
# met_alpha_rc_stream.R
#
# OptiSparseMET version of alpha_rc_stream() from OptiDesign.
# The met_ prefix prevents namespace conflicts when both packages are loaded.
# Function body is identical to alpha_rc_stream(). Only the name changes.
# ==============================================================================

#' Construct a stream-based repeated-check alpha row-column design
#'
#' @description
#' `met_alpha_rc_stream()` is the OptiSparseMET version of `alpha_rc_stream()`
#' from the OptiDesign package. The `met_` prefix avoids namespace conflicts
#' when both packages are loaded simultaneously. All arguments, return values,
#' and internal logic are identical to `alpha_rc_stream()` in OptiDesign.
#'
#' `met_alpha_rc_stream()` builds a randomised alpha-lattice incomplete block
#' design for breeding and agronomic experiments on a fixed field grid of
#' `n_rows x n_cols` plots. The field is converted to a one-dimensional
#' planting stream according to `order` and `serpentine`, partitioned into
#' `n_reps` contiguous replicate segments, and each replicate is further
#' divided into incomplete blocks. Every incomplete block contains all check
#' treatments plus a subset of entry treatments. Entry treatments appear
#' exactly once per replicate. Trailing unused field cells remain `NA` in the
#' layout matrix and field book.
#'
#' **Design evaluation has been separated from construction.** This function
#' returns the field book, layout matrix, and design metadata only. To compute
#' A-, D-, and CDmean optimality criteria, call
#' [met_evaluate_alpha_efficiency()] on the returned `field_book`. To search
#' for an optimised design using Random Restart, Simulated Annealing, or a
#' Genetic Algorithm, call [met_optimize_alpha_rc()].
#'
#' @details
#' ## Design structure
#'
#' Let:
#' - \eqn{r} = number of replicates (`n_reps`)
#' - \eqn{c} = number of checks (`length(check_treatments)`)
#' - \eqn{v} = number of entries (`length(entry_treatments)`)
#' - \eqn{b} = number of incomplete blocks per replicate
#'
#' Each replicate contains all \eqn{v} entries exactly once and all \eqn{c}
#' checks repeated once in each of the \eqn{b} incomplete blocks. Entries are
#' split across blocks as evenly as possible - block entry counts differ by at
#' most one. The number of used plots per replicate is:
#'
#' \deqn{\text{plots per rep} = v + b \times c}
#'
#' The total used plots in the design is:
#'
#' \deqn{\text{total used} = r \times (v + b \times c)}
#'
#' If the fixed field contains more cells than this total, unused trailing
#' cells remain `NA`. The function stops if the required plots exceed field
#' capacity.
#'
#' ## Block-size parameterisation
#'
#' Block-size constraints are expressed in terms of **total block size**
#' (checks + entries) because in a repeated-check design the whole block is
#' the operational unit. The mapping is:
#'
#' \deqn{\text{min\_entry\_slots} = \text{min\_block\_size} - c}
#' \deqn{\text{max\_entry\_slots} = \text{max\_block\_size} - c}
#'
#' For example, with \eqn{c = 3} checks, `min_block_size = 19` and
#' `max_block_size = 20` means each block holds between 16 and 17 entries plus
#' 3 check positions. Both bounds must be \eqn{\geq c}.
#'
#' ## Block-count determination
#'
#' After translating total-block-size limits into entry-slot limits, the number
#' of incomplete blocks per replicate \eqn{b} must satisfy a feasible range
#' \eqn{[b_\text{min},\, b_\text{max}]}:
#'
#' \describe{
#'   \item{Lower bound from `max_block_size`}{
#'     \eqn{b \geq \lceil v / (\text{max\_block\_size} - c)\rceil}
#'   }
#'   \item{Upper bound from `min_block_size`}{
#'     \eqn{b \leq \lfloor v / (\text{min\_block\_size} - c)\rfloor}
#'   }
#'   \item{Upper bound from field capacity}{
#'     \eqn{b \leq \lfloor (\text{total\_plots}/r - v) / c \rfloor}
#'   }
#' }
#'
#' The function stops with an informative error if
#' \eqn{b_\text{min} > b_\text{max}}. At least one of `n_blocks_per_rep`,
#' `min_block_size`, or `max_block_size` must be supplied.
#'
#' **User-fixed mode** (`n_blocks_per_rep` supplied): the value is validated
#' against the feasible range and all block-size constraints.
#'
#' **Automatic mode** (`n_blocks_per_rep = NULL`): \eqn{b} is set to
#' \eqn{b_\text{max}} - the largest feasible number of blocks, producing the
#' smallest feasible blocks, which generally gives the best local control.
#'
#' ## Field stream construction
#'
#' The field is mapped to a one-dimensional stream by iterating over rows
#' (`order = "row"`) or columns (`order = "column"`). When `serpentine = TRUE`,
#' alternate rows or columns reverse direction, producing a boustrophedon
#' planting path. Replicate segments occupy contiguous portions of the stream
#' and incomplete blocks occupy contiguous sub-segments within each replicate.
#'
#' ## Within-block arrangement
#'
#' Entries are randomised independently within each replicate, then split
#' across blocks. Checks are inserted at positions determined by
#' `check_placement`:
#' - `"systematic"`: positions computed as
#'   `round(seq(1, block_size, length.out = n_checks + 2))[2:(n_checks + 1)]`,
#'   distributing checks as evenly as possible.
#' - `"random"`: positions drawn uniformly without replacement.
#'
#' Within-block order is optimised over `attempts` random shuffles to minimise
#' the count of adjacent same-group treatment pairs. The search exits early
#' when the adjacency score reaches zero.
#'
#' ## Grouping and adjacency scoring
#'
#' The grouping label used for adjacency scoring is determined by
#' `cluster_source`:
#'
#' - `"Family"`: uses `entry_families` and `check_families` directly.
#'   `Gcluster` is `NA` for all plots.
#' - `"GRM"` or `"A"`: derives clusters from the leading principal components
#'   of the relationship matrix. The number of clusters equals
#'   `length(unique(entry_families))`. Cluster labels are prefixed `"G"`
#'   (GRM) or `"A"` (A matrix). Check treatments always use `check_families`
#'   and always have `Gcluster = NA`.
#'
#' ## Dispersion optimisation
#'
#' When `use_dispersion = TRUE`, a local swap search is applied after layout
#' construction. At each of `dispersion_iters` iterations, two non-check plots
#' are selected at random; the swap is accepted if it reduces the total
#' neighbourhood relatedness score - the sum of relationship matrix values
#' between all non-check neighbour pairs within Chebyshev distance
#' `dispersion_radius`. The matrix used is selected by `dispersion_source`.
#'
#' @param check_treatments Character vector of check treatment identifiers.
#'   Checks appear in every incomplete block of every replicate. Must not
#'   overlap with `entry_treatments` and must not contain duplicates.
#'
#' @param check_families Character vector of the same length as
#'   `check_treatments`. Family labels for checks, used for adjacency scoring.
#'   Always stored in the `Family` column; `Gcluster` is always `NA` for
#'   check plots regardless of `cluster_source`.
#'
#' @param entry_treatments Character vector of non-check treatment identifiers.
#'   Each entry appears exactly once per replicate. Must not contain duplicates.
#'
#' @param entry_families Character vector of the same length as
#'   `entry_treatments`. Family labels for entries. Used for adjacency scoring
#'   when `cluster_source = "Family"` and to set the number of clusters when
#'   `cluster_source %in% c("GRM", "A")`.
#'
#' @param n_reps Positive integer. Number of contiguous replicate segments.
#'
#' @param n_rows Positive integer. Number of field rows.
#'
#' @param n_cols Positive integer. Number of field columns.
#'
#' @param order Character. Stream traversal direction: `"row"` fills row by
#'   row; `"column"` fills column by column.
#'
#' @param serpentine Logical. If `TRUE`, alternate rows (when `order = "row"`)
#'   or alternate columns (when `order = "column"`) traverse in reverse
#'   direction.
#'
#' @param seed Optional integer. If `NULL`, a seed is drawn from
#'   `1:.Machine$integer.max` and returned as `seed_used`. Controls all
#'   randomised steps including entry randomisation, within-block shuffling,
#'   check placement, and dispersion optimisation.
#'
#' @param attempts Positive integer. Number of within-block shuffling attempts
#'   to minimise adjacent same-group pairs. Exits early when score reaches
#'   zero.
#'
#' @param warn_and_correct Logical. Retained for interface compatibility.
#'   Field capacity violations always raise an error.
#'
#' @param fix_rows Logical. Retained for interface compatibility. Field
#'   geometry is always fixed.
#'
#' @param cluster_source Character. Grouping source for adjacency scoring.
#'   `"Family"` uses `entry_families` directly. `"GRM"` and `"A"` derive
#'   clusters from PCA of the respective relationship matrix.
#'
#' @param GRM Optional square numeric matrix with rownames and colnames equal
#'   to line IDs. Required when `cluster_source = "GRM"` or when
#'   `use_dispersion = TRUE` and `dispersion_source = "GRM"`.
#'
#' @param A Optional square numeric matrix. Required when
#'   `cluster_source = "A"` or `dispersion_source = "A"`.
#'
#' @param id_map Optional data frame with columns `Treatment` and `LineID`.
#'   Required when treatment labels do not match relationship matrix rownames.
#'
#' @param cluster_method Character. Clustering algorithm: `"kmeans"` uses
#'   `stats::kmeans()` with `cluster_attempts` random starts; `"hclust"` uses
#'   Ward's D2 linkage.
#'
#' @param cluster_seed Integer. Seed for k-means initialisation, run in an
#'   isolated RNG scope.
#'
#' @param cluster_attempts Positive integer. Number of k-means random restarts.
#'   Ignored when `cluster_method = "hclust"`.
#'
#' @param n_pcs_use Positive number or `Inf`. Number of leading PCs for
#'   matrix-based clustering. Actual number used is
#'   `min(n_pcs_use, n_positive_eigenvalues - 1)`.
#'
#' @param n_blocks_per_rep Optional positive integer. If supplied, the exact
#'   number of incomplete blocks per replicate. Validated against the feasible
#'   range derived from `min_block_size`, `max_block_size`, and field capacity.
#'   If `NULL`, the block count is set automatically to \eqn{b_\text{max}}.
#'
#' @param min_block_size Optional positive integer \eqn{\geq c}. Minimum
#'   **total** block size (checks + entries). Sets an upper bound on \eqn{b}.
#'
#' @param max_block_size Optional positive integer \eqn{\geq c}. Maximum
#'   **total** block size (checks + entries). Sets a lower bound on \eqn{b}.
#'
#' @param check_placement Character. Check position rule within blocks.
#'   `"systematic"` distributes checks evenly; `"random"` draws positions
#'   uniformly without replacement.
#'
#' @param check_position_pattern Character. Retained for interface
#'   compatibility.
#'
#' @param use_dispersion Logical. If `TRUE`, apply a post-layout swap-based
#'   local dispersion optimisation to reduce pairwise relatedness among
#'   neighbouring non-check plots. Requires a relationship matrix specified
#'   by `dispersion_source`. Default `FALSE`.
#'
#' @param dispersion_source Character. Matrix used for dispersion scoring when
#'   `use_dispersion = TRUE`: `"K"`, `"A"`, or `"GRM"`. The corresponding
#'   argument must be non-`NULL`.
#'
#' @param dispersion_radius Positive integer. Chebyshev distance radius
#'   defining the neighbourhood for dispersion scoring.
#'
#' @param dispersion_iters Non-negative integer. Number of random swap
#'   proposals for dispersion optimisation.
#'
#' @param dispersion_seed Integer. Seed for the dispersion step, run in an
#'   isolated RNG scope. Defaults to `seed_used` when `NULL`.
#'
#' @param K Optional square numeric matrix. Used for dispersion scoring when
#'   `dispersion_source = "K"`.
#'
#' @param line_id_map Optional data frame with columns `Treatment` and
#'   `LineID`. Required when treatment labels do not match `rownames(K)` or
#'   the dispersion matrix.
#'
#' @param verbose Logical. If `TRUE`, prints a one-line summary of field
#'   usage, trailing `NA` plots, replicate sizes, block count and origin,
#'   block sizes in replicate 1, and active block-size bounds.
#'
#' @return A named list with four components:
#' \describe{
#'   \item{`layout_matrix`}{Character matrix of dimension `n_rows x n_cols`.
#'     Used cells contain treatment IDs; unused trailing cells are `NA`.}
#'   \item{`field_book`}{Data frame with one row per field plot. Columns:
#'     `Plot` (stream position), `Row`, `Column`, `Rep` (`NA` for unused),
#'     `IBlock` (global block index), `BlockInRep` (block index within rep),
#'     `Treatment` (`NA` for unused), `Family`, `Gcluster` (`NA` when
#'     `cluster_source = "Family"` or for check plots), `Check` (logical).}
#'   \item{`design_info`}{Named list: `n_rows`, `n_cols`, `total_plots`,
#'     `n_reps`, `n_checks`, `n_entries`, `n_blocks_per_rep` (resolved),
#'     `n_blocks_per_rep_user`, `min_block_size`, `max_block_size`,
#'     `min_entry_slots_per_block`, `max_entry_slots_per_block`, `rep_sizes`,
#'     `total_used_plots`, `trailing_na_plots`, `block_plan`, `block_meta`,
#'     `order`, `serpentine`.}
#'   \item{`seed_used`}{Integer. The random seed used.}
#' }
#'
#' @seealso [met_evaluate_alpha_efficiency()] to compute A, D, and CDmean
#'   optimality criteria on the returned `field_book`.
#'   [met_optimize_alpha_rc()] to search for a criterion-optimal design using
#'   Random Restart, Simulated Annealing, or a Genetic Algorithm.
#'
#' @examples
#' ## Basic construction: 167 entries, 3 checks, 3 reps
#' design <- met_alpha_rc_stream(
#'   check_treatments = c("CHK1", "CHK2", "CHK3"),
#'   check_families   = c("CHECK", "CHECK", "CHECK"),
#'   entry_treatments = paste0("G", 1:167),
#'   entry_families   = rep(paste0("F", 1:7), length.out = 167),
#'   n_reps           = 3,
#'   n_rows           = 30,
#'   n_cols           = 20,
#'   min_block_size   = 19,
#'   max_block_size   = 20
#' )
#'
#' dim(design$layout_matrix)             # 30 x 20
#' head(design$field_book)
#' design$design_info$n_blocks_per_rep
#' design$design_info$trailing_na_plots
#'
#' ## Evaluate separately
#' eff <- met_evaluate_alpha_efficiency(
#'   field_book         = design$field_book,
#'   n_rows             = 30,
#'   n_cols             = 20,
#'   check_treatments   = c("CHK1", "CHK2", "CHK3"),
#'   treatment_effect   = "fixed",
#'   residual_structure = "AR1xAR1",
#'   rho_row = 0.10, rho_col = 0.10
#' )
#' eff$A_criterion
#' eff$D_criterion
#'
#' @importFrom stats dist cutree hclust kmeans
#' @export
met_alpha_rc_stream <- function(
    check_treatments,
    check_families,
    entry_treatments,
    entry_families,
    n_reps,
    n_rows,
    n_cols,
    order                  = "column",
    serpentine             = FALSE,
    seed                   = NULL,
    attempts               = 5000,
    warn_and_correct       = TRUE,
    fix_rows               = TRUE,
    cluster_source         = c("Family", "GRM", "A"),
    GRM                    = NULL,
    A                      = NULL,
    id_map                 = NULL,
    cluster_method         = c("kmeans", "hclust"),
    cluster_seed           = 1,
    cluster_attempts       = 25,
    n_pcs_use              = Inf,
    n_blocks_per_rep       = NULL,
    min_block_size         = NULL,
    max_block_size         = NULL,
    check_placement        = c("systematic", "random"),
    check_position_pattern = c("spread", "corners_first"),
    use_dispersion         = FALSE,
    dispersion_source      = c("K", "A", "GRM"),
    dispersion_radius      = 1,
    dispersion_iters       = 2000,
    dispersion_seed        = 1,
    K                      = NULL,       # kept for dispersion_source = "K"
    line_id_map            = NULL,
    verbose                = TRUE
) {

  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Package 'Matrix' is required.")

  # ============================================================
  # 0. RNG
  # ============================================================
  seed_used <- if (is.null(seed)) sample.int(.Machine$integer.max, 1) else seed
  set.seed(seed_used)

  # ============================================================
  # 1. VALIDATION
  # ============================================================
  cluster_source         <- match.arg(cluster_source)
  cluster_method         <- match.arg(cluster_method)
  check_placement        <- match.arg(check_placement)
  check_position_pattern <- match.arg(check_position_pattern)
  dispersion_source      <- match.arg(dispersion_source)

  if (!order %in% c("row", "column"))
    stop("Invalid 'order'. Use 'row' or 'column'.")
  if (!is.logical(serpentine) || length(serpentine) != 1)
    stop("serpentine must be TRUE or FALSE.")
  if (length(n_reps) != 1 || is.na(n_reps) || n_reps < 1)
    stop("n_reps must be >= 1.")
  if (length(n_rows) != 1 || is.na(n_rows) || n_rows < 1)
    stop("n_rows must be >= 1.")
  if (length(n_cols) != 1 || is.na(n_cols) || n_cols < 1)
    stop("n_cols must be >= 1.")
  if (length(check_treatments) != length(check_families))
    stop("Length of check_families must match length of check_treatments.")
  if (length(entry_treatments) != length(entry_families))
    stop("Length of entry_families must match length of entry_treatments.")
  if (anyDuplicated(check_treatments)) stop("Duplicate check_treatments found.")
  if (anyDuplicated(entry_treatments)) stop("Duplicate entry_treatments found.")
  if (length(intersect(check_treatments, entry_treatments)) > 0)
    stop("A treatment cannot be both a check and an entry.")

  if (!is.null(n_blocks_per_rep)) {
    if (length(n_blocks_per_rep) != 1 || is.na(n_blocks_per_rep) || n_blocks_per_rep < 1)
      stop("n_blocks_per_rep must be a single integer >= 1 or NULL.")
    n_blocks_per_rep <- as.integer(n_blocks_per_rep)
  }
  if (!is.null(min_block_size)) {
    if (length(min_block_size) != 1 || is.na(min_block_size) || min_block_size < 1)
      stop("min_block_size must be a single integer >= 1 or NULL.")
    min_block_size <- as.integer(min_block_size)
  }
  if (!is.null(max_block_size)) {
    if (length(max_block_size) != 1 || is.na(max_block_size) || max_block_size < 1)
      stop("max_block_size must be a single integer >= 1 or NULL.")
    max_block_size <- as.integer(max_block_size)
  }
  if (!is.null(min_block_size) && !is.null(max_block_size) && min_block_size > max_block_size)
    stop("min_block_size cannot be greater than max_block_size.")

  if (dispersion_radius < 1)  stop("dispersion_radius must be >= 1.")
  if (dispersion_iters  < 0)  stop("dispersion_iters must be >= 0.")
  if (!(is.numeric(n_pcs_use) && length(n_pcs_use) == 1 &&
        (is.finite(n_pcs_use) || is.infinite(n_pcs_use)) && n_pcs_use > 0))
    stop("n_pcs_use must be a single positive number or Inf.")

  total_plots <- as.integer(n_rows * n_cols)
  n_checks    <- length(check_treatments)
  n_entries   <- length(entry_treatments)

  if (!is.null(min_block_size) && min_block_size < n_checks)
    stop(paste0("min_block_size (", min_block_size,
                ") cannot be smaller than the number of checks per block (", n_checks, ")."))
  if (!is.null(max_block_size) && max_block_size < n_checks)
    stop(paste0("max_block_size (", max_block_size,
                ") cannot be smaller than the number of checks per block (", n_checks, ")."))

  min_entry_slots_per_block <- if (is.null(min_block_size)) NULL else (min_block_size - n_checks)
  max_entry_slots_per_block <- if (is.null(max_block_size)) NULL else (max_block_size - n_checks)

  if (!is.null(min_entry_slots_per_block) && min_entry_slots_per_block < 0)
    stop("Derived min_entry_slots_per_block is negative. Check min_block_size.")
  if (!is.null(max_entry_slots_per_block) && max_entry_slots_per_block < 0)
    stop("Derived max_entry_slots_per_block is negative. Check max_block_size.")
  if (!is.null(min_entry_slots_per_block) && min_entry_slots_per_block < 1)
    min_entry_slots_per_block <- 1L
  if (!is.null(max_entry_slots_per_block) && max_entry_slots_per_block < 1)
    stop(paste0("max_block_size = ", max_block_size,
                " leaves fewer than 1 entry slot per block after allocating ", n_checks, " checks."))

  # ============================================================
  # 2. GLOBAL LOOKUPS
  # ============================================================
  family_lookup <- setNames(
    c(check_families, entry_families),
    c(check_treatments, entry_treatments)
  )

  gcluster_lookup <- setNames(
    rep(NA_character_, length(c(check_treatments, entry_treatments))),
    c(check_treatments, entry_treatments)
  )

  # ============================================================
  # 3. CLUSTERING FROM GRM / A (OPTIONAL)
  # ============================================================
  if (cluster_source %in% c("GRM", "A")) {
    Kc <- if (cluster_source == "GRM") GRM else A

    if (is.null(Kc))
      stop(paste0("cluster_source='", cluster_source, "' selected but matrix is NULL."))
    if (is.null(rownames(Kc)) || is.null(colnames(Kc)))
      stop("GRM/A must have rownames and colnames.")

    k_clusters <- length(unique(entry_families))
    if (k_clusters < 2) stop("Need at least 2 unique families for matrix-derived clustering.")

    if (is.null(id_map)) {
      line_ids <- setNames(entry_treatments, entry_treatments)
    } else {
      if (!is.data.frame(id_map) || !all(c("Treatment", "LineID") %in% names(id_map)))
        stop("id_map must contain columns Treatment and LineID.")
      line_ids <- setNames(id_map$LineID, id_map$Treatment)
    }

    ids         <- unname(line_ids[entry_treatments])
    missing_ids <- setdiff(ids, rownames(Kc))
    if (length(missing_ids) > 0) stop("Some entry LineIDs are missing in GRM/A.")

    Ksub <- Kc[ids, ids, drop = FALSE]
    eg   <- eigen(Ksub, symmetric = TRUE)
    pos  <- which(eg$values > 1e-10)
    if (length(pos) < 2) stop("GRM/A has too few positive eigenvalues for clustering.")

    max_possible_pcs <- min(length(pos), nrow(Ksub) - 1)
    n_pcs <- if (is.infinite(n_pcs_use)) max_possible_pcs else
      min(as.integer(n_pcs_use), max_possible_pcs)
    if (n_pcs < 2) stop("Fewer than 2 PCs available for clustering.")

    pcs <- eg$vectors[, pos[seq_len(n_pcs)], drop = FALSE]
    pcs <- sweep(pcs, 2, sqrt(eg$values[pos[seq_len(n_pcs)]]), `*`)

    if (cluster_method == "kmeans") {
      clust <- .with_local_seed(
        cluster_seed,
        stats::kmeans(pcs, centers = k_clusters, nstart = cluster_attempts)$cluster
      )
    } else {
      hc    <- stats::hclust(stats::dist(pcs), method = "ward.D2")
      clust <- stats::cutree(hc, k = k_clusters)
    }

    prefix <- if (cluster_source == "GRM") "G" else "A"
    gcluster_lookup[entry_treatments] <- paste0(prefix, clust)
  }

  get_adj_group <- function(trt) {
    if (cluster_source == "Family") {
      family_lookup[trt]
    } else if (trt %in% check_treatments) {
      family_lookup[trt]
    } else {
      gcluster_lookup[trt]
    }
  }

  # ============================================================
  # 4. LOCAL HELPERS
  # ============================================================
  build_positions <- function(n_rows, n_cols,
                              order = c("column", "row"),
                              serpentine = FALSE) {
    order <- match.arg(order)
    out <- vector("list", n_rows * n_cols); kk <- 1L

    if (order == "row") {
      for (r in seq_len(n_rows)) {
        cols <- seq_len(n_cols)
        if (isTRUE(serpentine) && (r %% 2L == 0L)) cols <- rev(cols)
        for (c in cols) { out[[kk]] <- c(Row = r, Column = c); kk <- kk + 1L }
      }
    } else {
      for (c in seq_len(n_cols)) {
        rows <- seq_len(n_rows)
        if (isTRUE(serpentine) && (c %% 2L == 0L)) rows <- rev(rows)
        for (r in rows) { out[[kk]] <- c(Row = r, Column = c); kk <- kk + 1L }
      }
    }
    ans        <- as.data.frame(do.call(rbind, out))
    ans$Row    <- as.integer(ans$Row)
    ans$Column <- as.integer(ans$Column)
    ans$Plot   <- seq_len(nrow(ans))
    ans
  }

  split_integer <- function(n, k) {
    if (k < 1) stop("k must be >= 1.")
    base <- n %/% k; rem <- n %% k
    ans  <- rep(base, k)
    if (rem > 0) ans[seq_len(rem)] <- ans[seq_len(rem)] + 1L
    as.integer(ans)
  }

  resolve_n_blocks_per_rep <- function(total_plots, n_reps, n_entries, n_checks,
                                       min_entry_slots_per_block = NULL,
                                       max_entry_slots_per_block = NULL,
                                       n_blocks_per_rep = NULL) {
    if (is.null(min_entry_slots_per_block) &&
        is.null(max_entry_slots_per_block) &&
        is.null(n_blocks_per_rep))
      stop("Provide n_blocks_per_rep or at least one of min_block_size / max_block_size.")

    lower_bound <- 1L
    if (!is.null(max_entry_slots_per_block))
      lower_bound <- max(lower_bound, ceiling(n_entries / max_entry_slots_per_block))

    upper_bound <- max(1L, n_entries)
    if (!is.null(min_entry_slots_per_block)) {
      if (min_entry_slots_per_block > n_entries)
        stop(paste0("The minimum block size leaves ", min_entry_slots_per_block,
                    " entry slots per block, which exceeds number of entries (", n_entries, ")."))
      upper_bound <- floor(n_entries / min_entry_slots_per_block)
    }

    if (n_checks > 0L) {
      max_b_by_field <- floor((total_plots / n_reps - n_entries) / n_checks)
      if (is.na(max_b_by_field) || max_b_by_field < 1L)
        stop(paste0("Field is too small to allocate ", n_entries, " entries across ", n_reps,
                    " replicates with repeated checks in every block."))
      upper_bound <- min(upper_bound, max_b_by_field)
    }

    if (lower_bound > upper_bound)
      stop(paste0("No feasible n_blocks_per_rep satisfies the total block-size constraints. ",
                  "After removing ", n_checks, " repeated checks per block, the entry-slot bounds imply at least ",
                  lower_bound, " block(s) per replicate and at most ", upper_bound, "."))

    if (!is.null(n_blocks_per_rep)) {
      if (n_blocks_per_rep < lower_bound || n_blocks_per_rep > upper_bound)
        stop(paste0("Requested n_blocks_per_rep = ", n_blocks_per_rep,
                    " is not feasible. Feasible range is [", lower_bound, ", ", upper_bound, "]."))
      entry_counts <- split_integer(n_entries, n_blocks_per_rep)
      if (!is.null(min_entry_slots_per_block) && any(entry_counts < min_entry_slots_per_block))
        stop("Requested n_blocks_per_rep violates min_block_size after accounting for checks.")
      if (!is.null(max_entry_slots_per_block) && any(entry_counts > max_entry_slots_per_block))
        stop("Requested n_blocks_per_rep violates max_block_size after accounting for checks.")
      used_per_rep <- n_entries + n_blocks_per_rep * n_checks
      total_used   <- n_reps * used_per_rep
      if (total_used > total_plots)
        stop(paste0("Requested n_blocks_per_rep = ", n_blocks_per_rep,
                    " does not fit the fixed field capacity. Required plots = ",
                    total_used, ", available = ", total_plots, "."))
      return(as.integer(n_blocks_per_rep))
    }

    resolved     <- as.integer(upper_bound)
    entry_counts <- split_integer(n_entries, resolved)
    if (!is.null(min_entry_slots_per_block) && any(entry_counts < min_entry_slots_per_block))
      stop("Internal error: derived block count violates min_block_size.")
    if (!is.null(max_entry_slots_per_block) && any(entry_counts > max_entry_slots_per_block))
      stop("Internal error: derived block count violates max_block_size.")
    resolved
  }

  make_block_plan <- function(n_reps, n_entries, n_checks, n_blocks_per_rep) {
    entry_counts <- split_integer(n_entries, n_blocks_per_rep)
    block_sizes  <- entry_counts + n_checks
    do.call(rbind, lapply(seq_len(n_reps), function(rr) {
      data.frame(Rep = rr, BlockInRep = seq_len(n_blocks_per_rep),
                 EntryCapacity = entry_counts, BlockSize = block_sizes,
                 stringsAsFactors = FALSE)
    }))
  }

  arrange_block_entries <- function(entries, checks, attempts,
                                    check_placement = c("systematic", "random")) {
    check_placement <- match.arg(check_placement)
    total_len <- length(entries) + length(checks)
    if (total_len == 0) return(character(0))
    if (length(checks) == 0)
      return(if (length(entries) > 1) sample(entries) else entries)

    insert_checks <- function(entries_ord) {
      out <- rep(NA_character_, total_len)
      if (check_placement == "systematic") {
        pos <- round(seq(1, total_len, length.out = length(checks) + 2))[2:(length(checks) + 1)]
      } else {
        pos <- sort(sample(seq_len(total_len), length(checks)))
      }
      out[pos]       <- checks
      out[is.na(out)] <- entries_ord
      out
    }

    score_adj <- function(vec) {
      if (length(vec) <= 1) return(0)
      grp <- vapply(vec, get_adj_group, character(1))
      sum(grp[-1] == grp[-length(grp)])
    }

    best <- c(checks, entries); best_score <- Inf
    for (aa in seq_len(max(1L, attempts))) {
      ent  <- if (length(entries) > 1) sample(entries) else entries
      cand <- insert_checks(ent)
      sc   <- score_adj(cand)
      if (sc < best_score) { best_score <- sc; best <- cand; if (best_score == 0) break }
    }
    best
  }

  apply_genetic_dispersion <- function(fb, check_treatments, Kdisp,
                                       line_id_map, radius, iters, seed_local) {
    .with_local_seed(seed_local, {
      trt      <- as.character(fb$Treatment)
      is_check <- trt %in% check_treatments
      movable  <- which(!is_check & !is.na(trt))
      if (length(movable) < 2 || iters <= 0) return(fb)

      non_trt <- unique(trt[!is_check & !is.na(trt)])
      if (is.null(line_id_map)) {
        line_ids <- setNames(non_trt, non_trt)
      } else {
        if (!is.data.frame(line_id_map) ||
            !all(c("Treatment", "LineID") %in% names(line_id_map)))
          stop("line_id_map must be a data.frame with columns: Treatment, LineID")
        line_ids <- setNames(line_id_map$LineID, line_id_map$Treatment)
      }
      ids  <- unname(line_ids[non_trt])
      miss <- setdiff(ids, rownames(Kdisp))
      if (length(miss) > 0)
        stop("Some non-check LineIDs are missing in dispersion matrix rownames/colnames.")

      Ksub             <- Kdisp[ids, ids, drop = FALSE]
      colnames(Ksub)   <- non_trt
      rownames(Ksub)   <- non_trt

      pairs     <- .build_neighbor_pairs(fb$Row, fb$Column, radius = radius)
      best_trt  <- trt
      best_sc   <- .score_dispersion(best_trt, is_check, Ksub, pairs)

      for (it in seq_len(iters)) {
        ij   <- sample(movable, 2, replace = FALSE)
        cand <- best_trt; cand[ij] <- cand[rev(ij)]
        sc   <- .score_dispersion(cand, is_check, Ksub, pairs)
        if (sc < best_sc) { best_sc <- sc; best_trt <- cand }
      }
      fb$Treatment <- best_trt
      fb
    })
  }

  # ============================================================
  # 5. GLOBAL STREAM -> REPLICATES -> BLOCKS
  # ============================================================
  pos_all <- build_positions(n_rows, n_cols, order = order, serpentine = serpentine)

  resolved_blocks_per_rep <- resolve_n_blocks_per_rep(
    total_plots = total_plots, n_reps = n_reps,
    n_entries = n_entries, n_checks = n_checks,
    min_entry_slots_per_block = min_entry_slots_per_block,
    max_entry_slots_per_block = max_entry_slots_per_block,
    n_blocks_per_rep = n_blocks_per_rep
  )

  block_plan <- make_block_plan(n_reps, n_entries, n_checks, resolved_blocks_per_rep)

  rep_used_size    <- sum(block_plan$BlockSize[block_plan$Rep == 1])
  rep_sizes        <- rep(rep_used_size, n_reps)
  total_used_plots <- sum(rep_sizes)

  if (total_used_plots > total_plots)
    stop(paste0("The fixed field (", n_rows, " x ", n_cols, " = ", total_plots,
                " plots) is too small. Required used plots = ", total_used_plots, "."))

  rep_starts <- cumsum(c(1L, head(rep_sizes, -1)))
  rep_ends   <- cumsum(rep_sizes)

  if (verbose) {
    msg <- paste0(
      "Fixed field = ", n_rows, " x ", n_cols,
      "; used plots = ", total_used_plots,
      "; trailing NA plots = ", total_plots - total_used_plots,
      "; replicate used sizes = {", paste(rep_sizes, collapse = ", "), "}",
      "; blocks/rep = ", resolved_blocks_per_rep,
      if (!is.null(n_blocks_per_rep)) " (user-fixed)" else " (derived)",
      "; block sizes in rep 1 = {",
      paste(block_plan$BlockSize[block_plan$Rep == 1], collapse = ", "), "}"
    )
    if (!is.null(min_block_size)) msg <- paste0(msg, "; min block size = ", min_block_size)
    if (!is.null(max_block_size)) msg <- paste0(msg, "; max block size = ", max_block_size)
    message(msg)
  }

  # ============================================================
  # 6. BUILD FIELD BOOK TEMPLATE
  # ============================================================
  field_book             <- pos_all
  field_book$Rep         <- NA_integer_
  field_book$IBlock      <- NA_integer_
  field_book$BlockInRep  <- NA_integer_
  field_book$Treatment   <- NA_character_
  field_book$Family      <- NA_character_
  field_book$Gcluster    <- NA_character_
  field_book$Check       <- FALSE

  for (rr in seq_len(n_reps)) {
    idx <- rep_starts[rr]:rep_ends[rr]
    field_book$Rep[idx] <- rr
  }

  global_block_counter <- 1L
  block_meta_list      <- list()

  for (rr in seq_len(n_reps)) {
    rep_idx         <- rep_starts[rr]:rep_ends[rr]
    rep_block_sizes <- block_plan$BlockSize[block_plan$Rep == rr]
    block_starts    <- cumsum(c(1L, head(rep_block_sizes, -1)))
    block_ends      <- cumsum(rep_block_sizes)

    for (bb in seq_along(rep_block_sizes)) {
      local_idx  <- block_starts[bb]:block_ends[bb]
      global_idx <- rep_idx[local_idx]
      field_book$IBlock[global_idx]     <- global_block_counter
      field_book$BlockInRep[global_idx] <- bb

      block_meta_list[[length(block_meta_list) + 1L]] <- data.frame(
        Rep           = rr,
        IBlock        = global_block_counter,
        BlockInRep    = bb,
        BlockSize     = rep_block_sizes[bb],
        EntryCapacity = block_plan$EntryCapacity[block_plan$Rep == rr][bb],
        StartStream   = min(global_idx),
        EndStream     = max(global_idx),
        stringsAsFactors = FALSE
      )
      global_block_counter <- global_block_counter + 1L
    }
  }

  block_meta <- do.call(rbind, block_meta_list)

  # ============================================================
  # 7. ALLOCATE TREATMENTS WITHIN EACH REPLICATE
  # ============================================================
  for (rr in seq_len(n_reps)) {
    rep_entries     <- sample(entry_treatments)
    rep_block_sizes <- block_plan$BlockSize[block_plan$Rep == rr]
    rep_entry_caps  <- block_plan$EntryCapacity[block_plan$Rep == rr]
    entry_split     <- split(rep_entries, rep(seq_along(rep_entry_caps), times = rep_entry_caps))
    rep_idx         <- rep_starts[rr]:rep_ends[rr]
    block_starts    <- cumsum(c(1L, head(rep_block_sizes, -1)))
    block_ends      <- cumsum(rep_block_sizes)

    for (bb in seq_along(rep_block_sizes)) {
      local_idx  <- block_starts[bb]:block_ends[bb]
      global_idx <- rep_idx[local_idx]
      ent_here   <- if (bb <= length(entry_split)) entry_split[[bb]] else character(0)
      block_trt  <- arrange_block_entries(ent_here, check_treatments, attempts, check_placement)

      field_book$Treatment[global_idx] <- block_trt
      field_book$Family[global_idx]    <- family_lookup[block_trt]
      field_book$Gcluster[global_idx]  <- ifelse(
        cluster_source == "Family" | block_trt %in% check_treatments,
        NA_character_, gcluster_lookup[block_trt]
      )
      field_book$Check[global_idx] <- block_trt %in% check_treatments
    }
  }

  # ============================================================
  # 8. OPTIONAL DISPERSION OPTIMIZATION
  # ============================================================
  if (isTRUE(use_dispersion)) {
    Kdisp <- switch(dispersion_source, "K" = K, "A" = A, "GRM" = GRM)
    if (is.null(Kdisp))
      stop("use_dispersion=TRUE but selected dispersion matrix is NULL.")
    if (is.null(rownames(Kdisp)) || is.null(colnames(Kdisp)))
      stop("Dispersion matrix must have rownames and colnames.")

    seed_disp_used <- if (!is.null(dispersion_seed)) dispersion_seed else seed_used
    used_idx <- which(!is.na(field_book$Treatment))
    fb_used  <- field_book[used_idx, , drop = FALSE]

    fb_used <- apply_genetic_dispersion(
      fb = fb_used, check_treatments = check_treatments,
      Kdisp = Kdisp, line_id_map = line_id_map,
      radius = dispersion_radius, iters = dispersion_iters,
      seed_local = seed_disp_used
    )

    field_book$Treatment[used_idx] <- fb_used$Treatment
    field_book$Family[used_idx]    <- family_lookup[fb_used$Treatment]
    field_book$Gcluster[used_idx]  <- ifelse(
      cluster_source == "Family" | fb_used$Treatment %in% check_treatments,
      NA_character_, gcluster_lookup[fb_used$Treatment]
    )
    field_book$Check[used_idx] <- fb_used$Treatment %in% check_treatments
  }

  # ============================================================
  # 9. LAYOUT MATRIX
  # ============================================================
  layout_matrix <- matrix(NA_character_, nrow = n_rows, ncol = n_cols)
  for (ii in seq_len(nrow(field_book)))
    layout_matrix[field_book$Row[ii], field_book$Column[ii]] <- field_book$Treatment[ii]

  # ============================================================
  # 10. DESIGN INFO
  # ============================================================
  design_info <- list(
    n_rows                    = n_rows,
    n_cols                    = n_cols,
    total_plots               = total_plots,
    n_reps                    = n_reps,
    n_checks                  = n_checks,
    n_entries                 = n_entries,
    n_blocks_per_rep          = resolved_blocks_per_rep,
    n_blocks_per_rep_user     = n_blocks_per_rep,
    min_block_size            = min_block_size,
    max_block_size            = max_block_size,
    min_entry_slots_per_block = min_entry_slots_per_block,
    max_entry_slots_per_block = max_entry_slots_per_block,
    rep_sizes                 = rep_sizes,
    total_used_plots          = total_used_plots,
    trailing_na_plots         = total_plots - total_used_plots,
    block_plan                = block_plan,
    block_meta                = block_meta,
    order                     = order,
    serpentine                = serpentine
  )

  field_book <- field_book[, c("Plot", "Row", "Column", "Rep", "IBlock",
                                "BlockInRep", "Treatment", "Family",
                                "Gcluster", "Check")]

  list(
    layout_matrix = layout_matrix,
    field_book    = field_book,
    design_info   = design_info,
    seed_used     = seed_used
  )
}
