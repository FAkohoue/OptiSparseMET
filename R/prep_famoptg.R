#' Construct a repeated-check block design with flexible replication structure
#'
#' `prep_famoptg()` builds a block field design in which check treatments
#' appear in every block and non-check treatments are allocated across blocks
#' according to specified replication levels. The same function produces three
#' design classes depending on the supplied treatment structure: an augmented
#' repeated-check design when all non-check treatments are unreplicated, a
#' partially replicated (p-rep) design when some non-check treatments are
#' replicated and others are not, and an RCBD-type repeated-check design when
#' all non-check treatments receive a common replication greater than 1.
#' Optionally, the function applies post-layout dispersion optimization to
#' reduce local relatedness among non-check entries, and evaluates the
#' resulting design under a mixed-model framework.
#'
#' @description
#' The allocation rules that govern all three design classes are:
#'
#' - Check treatments appear exactly once in every block.
#' - Replicated non-check treatments appear at most once per block; their
#'   replicates are distributed across distinct blocks.
#' - Unreplicated treatments appear exactly once in the full design.
#'
#' A treatment with replication \eqn{r} and \eqn{b \geq r} blocks therefore
#' occupies \eqn{r} distinct blocks. When \eqn{r = b}, the treatment appears in
#' every block, which is the RCBD analogue within this framework.
#'
#' Within each block, the 1D ordering of treatments is shuffled to reduce
#' adjacency between entries from the same genetic group, using up to `attempts`
#' shuffle proposals per block. This adjacency control operates on group labels
#' derived from `cluster_source`. Optional dispersion optimization applies a
#' swap-based local search after the full layout is placed, using pairwise
#' relatedness from a relationship matrix to penalize spatial proximity of
#' similar entries in the 2D field grid.
#'
#' @section Design classes:
#'
#' **Augmented repeated-check design**
#'
#' Obtained when `p_rep_treatments` is empty or `NULL` and all non-check
#' entries are supplied through `unreplicated_treatments`. Checks are repeated
#' in every block; all test entries appear once. Appropriate for early-stage
#' screening when the number of candidates exceeds what can be replicated.
#'
#' **Partially replicated (p-rep) repeated-check design**
#'
#' Obtained when some non-check treatments are in `p_rep_treatments` with
#' replication greater than 1 and others remain in `unreplicated_treatments`.
#' Checks are repeated in every block; selected candidates are replicated while
#' others receive a single plot.
#'
#' **RCBD-type repeated-check design**
#'
#' Obtained when all non-check treatments are in `p_rep_treatments` with a
#' common replication greater than 1 and `unreplicated_treatments` is empty. If
#' the common replication equals `n_blocks`, every non-check treatment appears
#' once in every block alongside the checks, producing the closest analogue to
#' a classical RCBD within this repeated-check framework. If the common
#' replication is less than `n_blocks`, balance is maintained but the design is
#' no longer a strict classical RCBD.
#'
#' @section Construction sequence:
#'
#' 1. Inputs are validated and field dimensions are reconciled with the total
#'    plot count required.
#' 2. Group labels are derived from family labels or from PCA-based clustering
#'    of `GRM` or `A`.
#' 3. Replicated non-check entries are assigned to blocks, subject to the
#'    constraint that no treatment appears twice in the same block.
#' 4. Unreplicated entries are distributed across blocks.
#' 5. Checks are inserted into each block using the chosen placement strategy.
#' 6. Within each block, treatment order is shuffled to reduce same-group
#'    adjacency in the 1D block sequence.
#' 7. The ordered treatment sequence is mapped to the field grid according to
#'    `order` and `serpentine`.
#' 8. Optionally, a swap-based dispersion search reduces local relatedness in
#'    the 2D layout.
#' 9. Optionally, mixed-model efficiency metrics are computed on the final
#'    design.
#'
#' @section Grouping and why it matters:
#'
#' Group labels serve two purposes. During block construction, they define which
#' entries are considered similar for the purpose of adjacency control in the 1D
#' block ordering. During optional dispersion optimization, a relationship matrix
#' is used to score pairwise relatedness between neighboring plots in the 2D
#' grid and a swap search reduces that score.
#'
#' Three grouping modes are available. `"Family"` reads labels directly from the
#' supplied family vectors. `"GRM"` derives cluster labels from the eigenstructure
#' of the genomic relationship matrix via PCA followed by k-means or hierarchical
#' clustering. `"A"` applies the same procedure to the pedigree relationship
#' matrix. Matrix-based grouping is preferable when family labels are too coarse
#' to capture the relevant genetic structure, or when dispersing genomically
#' similar materials more precisely is a priority.
#'
#' @section Argument dependencies:
#'
#' **Grouping mode**
#'
#' When `cluster_source = "Family"`, labels come from `check_families`,
#' `p_rep_families`, and `unreplicated_families`; `GRM`, `A`, `id_map`,
#' `cluster_method`, `cluster_seed`, `cluster_attempts`, and `n_pcs_use` are
#' ignored. When `cluster_source = "GRM"`, `GRM` is required; `id_map` is
#' needed only if treatment IDs differ from `rownames(GRM)`. When
#' `cluster_source = "A"`, the same applies with `A`. In both matrix-based
#' modes, `cluster_method`, `cluster_seed`, `cluster_attempts`, and `n_pcs_use`
#' are active.
#'
#' **Efficiency evaluation**
#'
#' When `eval_efficiency = FALSE`, all efficiency arguments are ignored. When
#' `eval_efficiency = TRUE` and `treatment_effect = "fixed"`, precision of
#' fixed-effect treatment contrasts is computed and `prediction_type`, `K`, and
#' `line_id_map` are ignored. When `treatment_effect = "random"`,
#' `prediction_type` becomes active: `"none"` skips random-effect efficiency;
#' `"IID"` treats entry effects as independent; `"GBLUP"` and `"PBLUP"` require
#' `K`, with `line_id_map` needed when treatment IDs differ from `rownames(K)`.
#'
#' **Residual structure**
#'
#' `"IID"` ignores `rho_row` and `rho_col`. `"AR1"` uses `rho_row` only.
#' `"AR1xAR1"` uses both. All residual arguments are inactive when
#' `eval_efficiency = FALSE`.
#'
#' **Check placement**
#'
#' `"optimal"` activates `check_opt_attempts` and evaluates multiple candidate
#' layouts, selecting the one that maximizes the mean minimum distance between
#' check positions in the 2D grid. `"random"` and `"systematic"` ignore
#' `check_opt_attempts`.
#'
#' **Dispersion optimization**
#'
#' When `use_dispersion = FALSE`, all dispersion arguments are ignored. When
#' `use_dispersion = TRUE`, `dispersion_source` selects the matrix used for
#' pairwise relatedness scoring: `"K"` requires `K` (with optional
#' `line_id_map`), `"A"` requires `A`, `"GRM"` requires `GRM`.
#'
#' **Field dimension correction**
#'
#' When `warn_and_correct = FALSE`, the function stops if `n_rows * n_cols`
#' does not match the total required plot count. When `warn_and_correct = TRUE`,
#' one dimension is adjusted upward to accommodate all plots; `fix_rows`
#' determines which dimension is held fixed.
#'
#' @details
#' Serpentine traversal uses `mod()` from the **pracma** package. When the
#' adjusted field size exceeds the plot count, surplus cells appear as `NA` in
#' `layout_matrix`; the `field_book` contains only assigned plots.
#'
#' The function separates construction logic, grouping logic, dispersion logic,
#' and efficiency logic. This means that, for example, family labels can drive
#' adjacency control while a genomic matrix drives dispersion optimization, and
#' efficiency evaluation can use a third, dedicated relationship matrix `K`.
#' None of these layers need to use the same source.
#'
#' @param check_treatments Character vector of check treatment IDs. Checks
#'   appear exactly once in every block and define the repeated reference
#'   structure of the design. Required for all design classes. Must not overlap
#'   with `p_rep_treatments` or `unreplicated_treatments`.
#'
#' @param check_families Character vector of the same length as
#'   `check_treatments`. Family or group labels for check treatments. Required
#'   regardless of `cluster_source` because checks always carry a group label
#'   in the output field book. Must align element-wise with `check_treatments`.
#'
#' @param p_rep_treatments Character vector of replicated non-check treatment
#'   IDs. Each treatment in this vector appears `p_rep_reps[i]` times in the
#'   full design, distributed across that many distinct blocks. Supply
#'   `character(0)` or `NULL` for an augmented design with no replicated
#'   non-check entries. Must align element-wise with `p_rep_reps` and
#'   `p_rep_families`.
#'
#' @param p_rep_reps Integer vector of replication counts, one per element of
#'   `p_rep_treatments`. Each value must satisfy `p_rep_reps[i] <= n_blocks`.
#'   A treatment with replication \eqn{r} is placed in \eqn{r} distinct blocks.
#'   Supply `integer(0)` when `p_rep_treatments` is empty.
#'
#' @param p_rep_families Character vector of the same length as
#'   `p_rep_treatments`. Family or group labels for replicated non-check
#'   treatments. Used directly when `cluster_source = "Family"`. When
#'   `cluster_source %in% c("GRM", "A")`, still used to anchor the number of
#'   clusters to the number of distinct families among non-check treatments.
#'   Supply `character(0)` when `p_rep_treatments` is empty.
#'
#' @param unreplicated_treatments Character vector of non-check treatments
#'   assigned exactly one plot each. Supply `character(0)` or `NULL` when all
#'   non-check treatments are replicated. In combination with empty
#'   `p_rep_treatments` and repeated checks, these entries define an augmented
#'   repeated-check design. Must align element-wise with
#'   `unreplicated_families`.
#'
#' @param unreplicated_families Character vector of the same length as
#'   `unreplicated_treatments`. Family or group labels for unreplicated
#'   treatments. Required whenever `unreplicated_treatments` is non-empty.
#'   Supply `character(0)` when `unreplicated_treatments` is empty.
#'
#' @param n_blocks Positive integer. Number of experimental blocks. Determines
#'   how many times checks are repeated and sets the upper bound on replication
#'   for any non-check treatment: `p_rep_reps[i] <= n_blocks`. When a
#'   replicated treatment should appear in every block, its replication must
#'   equal `n_blocks`.
#'
#' @param n_rows Positive integer. Number of rows in the field grid. Together
#'   with `n_cols`, determines total field capacity and the 2D coordinates in
#'   `layout_matrix` and `field_book`. When `warn_and_correct = TRUE` and
#'   `fix_rows = TRUE`, `n_rows` is held fixed and `n_cols` is adjusted.
#'
#' @param n_cols Positive integer. Number of columns in the field grid. When
#'   `warn_and_correct = TRUE` and `fix_rows = FALSE`, `n_cols` is held fixed
#'   and `n_rows` is adjusted.
#'
#' @param order Character scalar. Direction in which the field grid is
#'   traversed when mapping the ordered treatment sequence to spatial
#'   coordinates. `"row"` fills row by row; `"column"` fills column by column.
#'   Should reflect the physical planting or harvesting direction. Interacts
#'   with `serpentine`.
#'
#' @param serpentine Logical. If `TRUE`, alternate rows (when `order = "row"`)
#'   or alternate columns (when `order = "column"`) reverse traversal direction,
#'   producing a boustrophedon mapping from treatment sequence to field
#'   coordinates. Does not affect block composition; affects only the spatial
#'   positions assigned to each treatment.
#'
#' @param seed Optional integer. Random seed controlling block assignment of
#'   replicated treatments, within-block shuffling, check placement candidates,
#'   and dispersion optimization when `dispersion_seed` is `NULL`. If `NULL`,
#'   a seed is generated internally and returned as `seed_used`.
#'
#' @param attempts Positive integer. Maximum number of shuffle proposals per
#'   block used to reduce same-group adjacency in the 1D block ordering.
#'   Larger values increase search effort without changing the algorithm. Most
#'   relevant when group structure is highly imbalanced within blocks or when
#'   many similar entries occur together.
#'
#' @param warn_and_correct Logical. If `FALSE`, the function stops when
#'   `n_rows * n_cols` does not equal the total required plot count. If `TRUE`,
#'   one dimension is expanded upward to accommodate all plots, with surplus
#'   cells left as `NA` in `layout_matrix`. Which dimension is adjusted is
#'   controlled by `fix_rows`.
#'
#' @param fix_rows Logical. Active only when `warn_and_correct = TRUE`. If
#'   `TRUE`, `n_rows` is held fixed and `n_cols` is adjusted. If `FALSE`,
#'   `n_cols` is held fixed and `n_rows` is adjusted. Use `TRUE` when the
#'   number of field rows is physically constrained.
#'
#' @param cluster_source Character scalar. Source of group labels for
#'   within-block adjacency control. `"Family"` uses the supplied family
#'   vectors. `"GRM"` derives cluster labels from the eigenstructure of `GRM`.
#'   `"A"` derives cluster labels from the eigenstructure of `A`. Determines
#'   which of `GRM`, `A`, `id_map`, `cluster_method`, `cluster_seed`,
#'   `cluster_attempts`, and `n_pcs_use` become active.
#'
#' @param GRM Optional numeric matrix. Genomic relationship matrix. Required
#'   when `cluster_source = "GRM"` or when `use_dispersion = TRUE` and
#'   `dispersion_source = "GRM"`. Must be square with row and column names
#'   matching treatment IDs or reachable through `id_map`.
#'
#' @param A Optional numeric matrix. Pedigree-based numerator relationship
#'   matrix. Required when `cluster_source = "A"` or when
#'   `use_dispersion = TRUE` and `dispersion_source = "A"`. Same naming
#'   requirements as `GRM`.
#'
#' @param id_map Optional data frame with columns `Treatment` and `LineID`.
#'   Required only when `cluster_source %in% c("GRM", "A")` and treatment IDs
#'   in the design do not match the row names of the selected relationship
#'   matrix.
#'
#' @param cluster_method Character scalar. Clustering algorithm applied to PCA
#'   scores in matrix-based grouping. `"kmeans"` uses k-means with
#'   `cluster_attempts` random restarts seeded by `cluster_seed`. `"hclust"`
#'   uses Ward's criterion hierarchical clustering. Ignored when
#'   `cluster_source = "Family"`.
#'
#' @param cluster_seed Integer. Seed for k-means initialization. Active only
#'   when `cluster_source %in% c("GRM", "A")` and `cluster_method = "kmeans"`.
#'
#' @param cluster_attempts Integer. Number of random restarts for k-means.
#'   Active only when `cluster_source %in% c("GRM", "A")` and
#'   `cluster_method = "kmeans"`. Larger values reduce the risk of poor local
#'   optima.
#'
#' @param n_pcs_use Integer or `Inf`. Number of leading principal components
#'   retained for PCA-based clustering. `Inf` retains all components
#'   corresponding to positive eigenvalues. Smaller values preserve only broad
#'   genetic structure. Must be at least 2. Ignored when
#'   `cluster_source = "Family"`.
#'
#' @param eval_efficiency Logical. If `TRUE`, mixed-model efficiency metrics
#'   are computed on the final design. If `FALSE`, all efficiency arguments are
#'   ignored and `efficiency` is returned as `NULL`.
#'
#' @param treatment_effect Character scalar. Specifies whether non-check
#'   treatments enter the efficiency model as fixed or random effects. `"fixed"`
#'   computes precision of pairwise treatment contrasts. `"random"` computes
#'   average prediction error variance. Active only when
#'   `eval_efficiency = TRUE`.
#'
#' @param prediction_type Character scalar. Random-effect model for efficiency
#'   evaluation. `"none"` skips. `"IID"` assumes independent entry effects.
#'   `"GBLUP"` and `"PBLUP"` use the relationship matrix `K`. Active only when
#'   `eval_efficiency = TRUE` and `treatment_effect = "random"`.
#'
#' @param K Optional numeric matrix. Relationship matrix serving two purposes:
#'   (1) random-effect efficiency when `prediction_type %in% c("GBLUP",
#'   "PBLUP")`; (2) dispersion optimization when `use_dispersion = TRUE` and
#'   `dispersion_source = "K"`. Ignored otherwise. Row and column names must
#'   match treatment IDs or be reachable through `line_id_map`.
#'
#' @param line_id_map Optional data frame with columns `Treatment` and
#'   `LineID`. Required only when `K` is active and treatment labels differ
#'   from `rownames(K)`. Serves the same role for `K` that `id_map` serves for
#'   `GRM` and `A`.
#'
#' @param varcomp Named list of variance components for efficiency evaluation.
#'   Must contain `sigma_e2`, `sigma_g2`, `sigma_b2`, `sigma_r2`, and
#'   `sigma_c2`. All values except `sigma_g2` must be strictly positive.
#'   `sigma_g2` must be positive when `treatment_effect = "random"`. Active
#'   only when `eval_efficiency = TRUE`. The ratio of `sigma_b2` to `sigma_e2`
#'   has the largest effect on efficiency metrics for block designs.
#'
#' @param check_as_fixed Logical. If `TRUE`, check treatments are included as
#'   fixed indicator columns in the mixed-model coefficient matrix during
#'   efficiency evaluation. Active only when `eval_efficiency = TRUE`.
#'
#' @param residual_structure Character scalar. Residual covariance structure
#'   for efficiency evaluation. `"IID"` assumes independent residuals. `"AR1"`
#'   introduces row-wise autocorrelation controlled by `rho_row`. `"AR1xAR1"`
#'   introduces both row-wise and column-wise autocorrelation controlled by
#'   `rho_row` and `rho_col`. Active only when `eval_efficiency = TRUE`.
#'
#' @param rho_row Numeric. AR1 autocorrelation parameter in the row direction.
#'   Must satisfy \eqn{|\rho_{\text{row}}| < 1}. Active when
#'   `eval_efficiency = TRUE` and `residual_structure %in% c("AR1",
#'   "AR1xAR1")`.
#'
#' @param rho_col Numeric. AR1 autocorrelation parameter in the column
#'   direction. Must satisfy \eqn{|\rho_{\text{col}}| < 1}. Active when
#'   `eval_efficiency = TRUE` and `residual_structure = "AR1xAR1"`.
#'
#' @param spatial_engine Character scalar. Computational strategy for the
#'   residual precision matrix when `residual_structure != "IID"`. `"auto"`
#'   selects `"dense"` when the number of observed plots is at most
#'   `dense_max_n` and `"sparse"` otherwise. `"dense"` and `"sparse"` force
#'   the respective implementation. Active only when `eval_efficiency = TRUE`
#'   and `residual_structure != "IID"`.
#'
#' @param dense_max_n Integer. Plot count threshold used by `spatial_engine =
#'   "auto"` to switch between dense and sparse computation. Active only when
#'   `eval_efficiency = TRUE`, `residual_structure != "IID"`, and
#'   `spatial_engine = "auto"`.
#'
#' @param eff_trace_samples Integer. Number of Hutchinson trace samples used
#'   when the target treatment dimension exceeds `eff_full_max` and exact
#'   efficiency extraction is replaced by stochastic approximation. Larger
#'   values reduce approximation variance at the cost of runtime. Active only
#'   when `eval_efficiency = TRUE`.
#'
#' @param eff_full_max Integer. Maximum treatment dimension for which exact
#'   efficiency extraction via linear solve is attempted. Above this threshold
#'   the function switches to stochastic trace estimation. Active only when
#'   `eval_efficiency = TRUE`.
#'
#' @param check_placement Character scalar. Strategy for positioning checks
#'   within blocks. `"random"` draws check positions freely within each block.
#'   `"systematic"` spaces checks approximately evenly in the 1D block
#'   ordering. `"optimal"` evaluates `check_opt_attempts` candidate layouts
#'   and selects the one maximizing the mean minimum Euclidean distance between
#'   check positions in the 2D grid.
#'
#' @param check_opt_attempts Positive integer. Number of candidate layouts
#'   evaluated when `check_placement = "optimal"`. Larger values increase the
#'   chance of finding a spatially well-dispersed check arrangement at higher
#'   runtime cost. Ignored unless `check_placement = "optimal"`.
#'
#' @param use_dispersion Logical. If `TRUE`, a swap-based local search is
#'   applied after layout construction to reduce spatial concentration of
#'   genetically similar non-check entries. The pairwise relatedness matrix
#'   is selected by `dispersion_source`. If `FALSE`, all dispersion arguments
#'   are ignored.
#'
#' @param dispersion_source Character scalar. Relationship matrix used for
#'   pairwise relatedness scoring in dispersion optimization. `"K"` uses `K`
#'   (with optional `line_id_map`), `"A"` uses `A`, `"GRM"` uses `GRM`.
#'   Active only when `use_dispersion = TRUE`.
#'
#' @param dispersion_radius Positive integer. Chebyshev neighborhood radius
#'   for dispersion scoring. Two plots are neighbors when
#'   \eqn{\max(|\Delta r|, |\Delta c|) \leq \texttt{dispersion\_radius}}.
#'   `1` considers immediate neighbors only; larger values extend the local
#'   window over which relatedness is penalized. Active only when
#'   `use_dispersion = TRUE`.
#'
#' @param dispersion_iters Non-negative integer. Number of swap proposals in
#'   the dispersion search. Larger values allow more extensive search at higher
#'   runtime cost. Active only when `use_dispersion = TRUE`.
#'
#' @param dispersion_seed Optional integer. Seed for the dispersion swap
#'   search, applied independently of `seed`. When `NULL`, `seed_used` is
#'   applied. Use when the initial layout and the dispersion refinement should
#'   be reproducible independently.
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{`layout_matrix`}{Character matrix of dimension `n_rows x n_cols`.
#'     Each cell contains a treatment ID or `NA` for surplus unused positions.}
#'   \item{`field_book`}{Data frame with one row per assigned plot. Columns:
#'     `Treatment`, `Family`, `Gcluster`, `Block`, `Plot`, `Row`, `Column`.
#'     `Gcluster` is `NA` for check treatments and for all treatments when
#'     `cluster_source = "Family"`.}
#'   \item{`efficiency`}{`NULL` if `eval_efficiency = FALSE`; otherwise a
#'     named list of efficiency metrics including `mode`, `A`, `D`,
#'     `mean_PEV` or `mean_VarDiff`, `n_trt` or `n_lines`, `residual_structure_used`,
#'     `spatial_engine_used`, and `notes`.}
#'   \item{`seed_used`}{Integer. The random seed used internally, returned for
#'     reproducibility records.}
#' }
#'
#' @examples
#' data("OptiSparseMET_example_data", package = "OptiSparseMET")
#' x <- OptiSparseMET_example_data
#'
#' env_trt <- x$treatments[1:24]
#' env_fam <- x$treatment_info$Family[match(env_trt, x$treatment_info$Treatment)]
#'
#' ## Example 1: p-rep repeated-check design — first 8 entries replicated
#' ## twice, remaining 16 unreplicated. Checks in every block.
#' out_prep <- prep_famoptg(
#'   check_treatments         = c("CHK1", "CHK2"),
#'   check_families           = c("CHECK", "CHECK"),
#'   p_rep_treatments         = env_trt[1:8],
#'   p_rep_reps               = rep(2L, 8),
#'   p_rep_families           = env_fam[1:8],
#'   unreplicated_treatments  = env_trt[9:24],
#'   unreplicated_families    = env_fam[9:24],
#'   n_blocks                 = 4,
#'   n_rows                   = 6,
#'   n_cols                   = 8,
#'   order                    = "row",
#'   serpentine               = TRUE,
#'   seed                     = 123,
#'   cluster_source           = "Family",
#'   eval_efficiency          = FALSE,
#'   use_dispersion           = FALSE
#' )
#'
#' dim(out_prep$layout_matrix)   # 6 x 8
#' head(out_prep$field_book)
#' out_prep$seed_used
#'
#' \dontrun{
#' ## Example 2: augmented repeated-check design — no replicated entries,
#' ## all 24 test entries unreplicated.
#' out_aug <- prep_famoptg(
#'   check_treatments         = c("CHK1", "CHK2"),
#'   check_families           = c("CHECK", "CHECK"),
#'   p_rep_treatments         = character(0),
#'   p_rep_reps               = integer(0),
#'   p_rep_families           = character(0),
#'   unreplicated_treatments  = env_trt,
#'   unreplicated_families    = env_fam,
#'   n_blocks                 = 4,
#'   n_rows                   = 5,
#'   n_cols                   = 6,
#'   order                    = "row",
#'   serpentine               = FALSE,
#'   seed                     = 123,
#'   cluster_source           = "Family",
#'   eval_efficiency          = FALSE,
#'   use_dispersion           = FALSE
#' )
#'
#' dim(out_aug$layout_matrix)
#' head(out_aug$field_book)
#'
#' ## Example 3: RCBD-type repeated-check design — all 24 entries replicated
#' ## 4 times (once per block). Equivalent to a classical RCBD with repeated
#' ## checks.
#' out_rcbd <- prep_famoptg(
#'   check_treatments         = c("CHK1", "CHK2"),
#'   check_families           = c("CHECK", "CHECK"),
#'   p_rep_treatments         = env_trt,
#'   p_rep_reps               = rep(4L, length(env_trt)),
#'   p_rep_families           = env_fam,
#'   unreplicated_treatments  = character(0),
#'   unreplicated_families    = character(0),
#'   n_blocks                 = 4,
#'   n_rows                   = 14,
#'   n_cols                   = 8,
#'   order                    = "row",
#'   serpentine               = TRUE,
#'   seed                     = 123,
#'   cluster_source           = "Family",
#'   eval_efficiency          = FALSE,
#'   use_dispersion           = FALSE
#' )
#'
#' dim(out_rcbd$layout_matrix)
#' head(out_rcbd$field_book)
#' }
#'
#' @importFrom stats runif setNames
#' @importFrom utils head
#' @importFrom pracma mod
#' @export
#' 
prep_famoptg <- function(
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
    varcomp = list(
      sigma_e2 = 1,
      sigma_g2 = 1,
      sigma_b2 = 1,
      sigma_r2 = 1,
      sigma_c2 = 1
    ),
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
) {
  
  # ============================================================
  # 0. RNG POLICY
  # ============================================================
  # Reproducible if seed provided; otherwise random but record the seed used.
  seed_used <- seed
  if (is.null(seed_used)) {
    seed_used <- sample.int(.Machine$integer.max, 1)
  }
  set.seed(seed_used)
  
  # Utility: set a seed locally without hijacking the global RNG stream
  .with_local_seed <- function(seed_local, expr) {
    old_seed <- NULL
    has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (has_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    
    set.seed(seed_local)
    on.exit({
      if (has_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      }
    }, add = TRUE)
    
    force(expr)
  }
  
  # ============================================================
  # 1. INPUT VALIDATION AND ARGUMENT NORMALIZATION
  # ============================================================
  cluster_source <- match.arg(cluster_source)
  cluster_method <- match.arg(cluster_method)
  treatment_effect <- match.arg(treatment_effect)
  prediction_type <- match.arg(prediction_type)
  residual_structure <- match.arg(residual_structure)
  spatial_engine <- match.arg(spatial_engine)
  check_placement <- match.arg(check_placement)
  dispersion_source <- match.arg(dispersion_source)
  
  if (length(check_treatments) != length(check_families)) {
    stop("Length of check_families must match length of check_treatments.")
  }
  if (length(p_rep_treatments) != length(p_rep_reps) ||
      length(p_rep_treatments) != length(p_rep_families)) {
    stop("Lengths of p_rep_treatments, p_rep_reps, and p_rep_families must all match.")
  }
  if (length(unreplicated_treatments) != length(unreplicated_families)) {
    stop("Length of unreplicated_families must match length of unreplicated_treatments.")
  }
  if (any(p_rep_reps > n_blocks)) {
    stop("Each p-rep treatment replication count must not exceed n_blocks.")
  }
  if (n_blocks < 1) {
    stop("n_blocks must be at least 1.")
  }
  if (!order %in% c("column", "row")) {
    stop("Invalid 'order'. Use 'row' or 'column'.")
  }
  if (!(
    is.numeric(n_pcs_use) && length(n_pcs_use) == 1 &&
    (is.finite(n_pcs_use) || is.infinite(n_pcs_use)) && n_pcs_use > 0
  )) {
    stop("n_pcs_use must be a single positive number or Inf.")
  }
  if (!is.list(varcomp) ||
      !all(c("sigma_e2", "sigma_g2", "sigma_b2", "sigma_r2", "sigma_c2") %in% names(varcomp))) {
    stop("varcomp must be a named list with sigma_e2, sigma_g2, sigma_b2, sigma_r2, sigma_c2.")
  }
  if (residual_structure != "IID") {
    if (!is.numeric(rho_row) || length(rho_row) != 1 || abs(rho_row) >= 1) {
      stop("rho_row must satisfy |rho_row| < 1.")
    }
    if (!is.numeric(rho_col) || length(rho_col) != 1 || abs(rho_col) >= 1) {
      stop("rho_col must satisfy |rho_col| < 1.")
    }
  }
  if (dispersion_radius < 1) stop("dispersion_radius must be >= 1.")
  if (dispersion_iters < 0) stop("dispersion_iters must be >= 0.")
  if (check_opt_attempts < 1) stop("check_opt_attempts must be >= 1.")
  
  # ============================================================
  # 2. FIELD DIMENSION ACCOUNTING AND CORRECTION
  # ============================================================
  total_checks <- n_blocks * length(check_treatments)
  total_prep <- sum(p_rep_reps)
  total_unrep <- length(unreplicated_treatments)
  total_required <- total_checks + total_prep + total_unrep
  
  field_size <- n_rows * n_cols
  if (field_size != total_required) {
    if (warn_and_correct) {
      warning(
        paste0(
          "Field size (", n_rows, " x ", n_cols, " = ", field_size,
          ") does not match required (", total_required, "). Adjusting dimensions."
        )
      )
      if (fix_rows) {
        n_cols <- ceiling(total_required / n_rows)
      } else {
        n_rows <- ceiling(total_required / n_cols)
      }
      field_size <- n_rows * n_cols
    } else {
      stop(
        paste0(
          "Provided field size (", field_size,
          ") does not match required (", total_required,
          "). Adjust n_rows/n_cols or enable warn_and_correct."
        )
      )
    }
  }
  
  # ============================================================
  # 3. CLUSTER PREPARATION (FAMILY OR GRM/A-BASED)
  # ============================================================
  family_lookup <- setNames(
    c(check_families, p_rep_families, unreplicated_families),
    c(check_treatments, p_rep_treatments, unreplicated_treatments)
  )
  
  checks_trt <- check_treatments
  noncheck_trt <- c(p_rep_treatments, unreplicated_treatments)
  
  gcluster_lookup <- setNames(
    rep(NA_character_, length(c(checks_trt, noncheck_trt))),
    c(checks_trt, noncheck_trt)
  )
  
  if (cluster_source %in% c("GRM", "A")) {
    Kc <- if (cluster_source == "GRM") GRM else A
    if (is.null(Kc)) stop(paste0("cluster_source='", cluster_source, "' selected but matrix is NULL."))
    if (is.null(rownames(Kc)) || is.null(colnames(Kc))) stop("GRM/A must have rownames and colnames.")
    
    noncheck_fams <- c(p_rep_families, unreplicated_families)
    k_clusters <- length(unique(noncheck_fams))
    if (k_clusters < 2) stop("Non-check treatments have <2 unique families; clustering is not meaningful.")
    
    if (is.null(id_map)) {
      line_ids <- setNames(noncheck_trt, noncheck_trt)
    } else {
      if (!is.data.frame(id_map) || !all(c("Treatment", "LineID") %in% names(id_map))) {
        stop("id_map must be a data.frame with columns: Treatment, LineID")
      }
      line_ids <- setNames(id_map$LineID, id_map$Treatment)
    }
    
    noncheck_line_ids <- unname(line_ids[noncheck_trt])
    missing <- setdiff(noncheck_line_ids, rownames(Kc))
    if (length(missing) > 0) stop("Some non-check LineIDs are not found in GRM/A rownames.")
    
    Ksub <- Kc[noncheck_line_ids, noncheck_line_ids, drop = FALSE]
    eg <- eigen(Ksub, symmetric = TRUE)
    pos <- which(eg$values > 1e-10)
    if (length(pos) < 2) stop("GRM/A has too few positive eigenvalues for PCA clustering.")
    
    max_possible_pcs <- min(length(pos), nrow(Ksub) - 1)
    if (is.infinite(n_pcs_use)) {
      n_pcs <- max_possible_pcs
    } else {
      n_pcs_req <- as.integer(n_pcs_use)
      n_pcs <- min(n_pcs_req, max_possible_pcs)
      if (n_pcs_req > max_possible_pcs) {
        warning(
          paste0(
            "Requested n_pcs_use=", n_pcs_use,
            " but only ", max_possible_pcs, " PCs are available; using ", n_pcs, "."
          )
        )
      }
    }
    if (n_pcs < 2) stop("After bounding, fewer than 2 PCs are available for clustering.")
    
    pcs <- eg$vectors[, pos[seq_len(n_pcs)], drop = FALSE]
    pcs <- sweep(pcs, 2, sqrt(eg$values[pos[seq_len(n_pcs)]]), `*`)
    
    if (cluster_method == "kmeans") {
      clust <- .with_local_seed(cluster_seed, {
        stats::kmeans(pcs, centers = k_clusters, nstart = cluster_attempts)$cluster
      })
    } else {
      d <- stats::dist(pcs)
      hc <- stats::hclust(d, method = "ward.D2")
      clust <- stats::cutree(hc, k = k_clusters)
    }
    
    prefix <- if (cluster_source == "GRM") "G" else "A"
    gcluster_lookup[noncheck_trt] <- paste0(prefix, clust)
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
  
  place_checks_in_block <- function(checks, others, mode = c("random", "systematic")) {
    mode <- match.arg(mode)
    total_len <- length(checks) + length(others)
    
    if (mode == "random") {
      return(sample(c(checks, others)))
    }
    
    pos <- round(seq(1, total_len, length.out = length(checks) + 2))[2:(length(checks) + 1)]
    pos <- pmin(pmax(pos, 1), total_len)
    
    out <- rep(NA_character_, total_len)
    out[pos] <- checks
    out[is.na(out)] <- sample(others)
    out
  }
  
  # ============================================================
  # 4. BLOCK CONSTRUCTION
  # ============================================================
  p_rep_assignments <- vector("list", length(p_rep_treatments))
  names(p_rep_assignments) <- p_rep_treatments
  
  available_blocks <- vector("list", length(p_rep_treatments))
  for (i in seq_along(p_rep_treatments)) {
    available_blocks[[i]] <- sample(seq_len(n_blocks))
  }
  
  for (i in seq_along(p_rep_treatments)) {
    trt_i <- p_rep_treatments[i]
    reps_i <- p_rep_reps[i]
    if (reps_i > length(available_blocks[[i]])) {
      stop(paste0("Not enough unique blocks available for p-rep treatment '", trt_i, "'."))
    }
    assigned <- sample(available_blocks[[i]], reps_i)
    available_blocks[[i]] <- setdiff(available_blocks[[i]], assigned)
    p_rep_assignments[[trt_i]] <- assigned
  }
  
  unrep_treatments_shuffled <- sample(unreplicated_treatments)
  block_unrep_list <- split(
    unrep_treatments_shuffled,
    rep(seq_len(n_blocks), length.out = length(unrep_treatments_shuffled))
  )
  
  build_blocks_once <- function(seed_offset = 0) {
    if (!is.null(seed_used)) {
      set.seed(seed_used + seed_offset)
    }
    
    blocks <- vector("list", n_blocks)
    
    for (b in seq_len(n_blocks)) {
      others <- character(0)
      
      for (trt in names(p_rep_assignments)) {
        if (is.element(b, p_rep_assignments[[trt]])) {
          others <- c(others, trt)
        }
      }
      if (b <= length(block_unrep_list)) {
        others <- c(others, block_unrep_list[[b]])
      }
      
      if (check_placement %in% c("random", "systematic")) {
        block_treatments <- place_checks_in_block(check_treatments, others, mode = check_placement)
      } else {
        block_treatments <- place_checks_in_block(check_treatments, others, mode = "random")
      }
      
      block_family <- vapply(block_treatments, function(trt) family_lookup[trt], character(1))
      block_gcluster <- vapply(block_treatments, function(trt) gcluster_lookup[trt], character(1))
      block_adj <- vapply(block_treatments, function(trt) get_adj_group(trt), character(1))
      
      valid_order <- FALSE
      attempt <- 1
      trt_sh <- block_treatments
      fam_sh <- block_family
      gcl_sh <- block_gcluster
      
      while (!valid_order && attempt <= attempts) {
        ord <- sample(seq_along(block_treatments))
        trt_sh <- block_treatments[ord]
        fam_sh <- block_family[ord]
        gcl_sh <- block_gcluster[ord]
        adj_sh <- block_adj[ord]
        
        if (!any(adj_sh[-1] == adj_sh[-length(adj_sh)])) {
          valid_order <- TRUE
        }
        attempt <- attempt + 1
      }
      
      if (!valid_order) {
        warning(
          paste0(
            "Could not fully avoid adjacent groups in block ", b,
            " after ", attempts, " attempts."
          )
        )
      }
      
      blocks[[b]] <- data.frame(
        Treatment = trt_sh,
        Family = fam_sh,
        Gcluster = gcl_sh,
        Block = b,
        stringsAsFactors = FALSE
      )
    }
    
    blocks
  }
  
  build_positions <- function(n_rows, n_cols, order, serpentine) {
    positions <- vector("list", n_rows * n_cols)
    kpos <- 1
    
    if (order == "row") {
      for (r in seq_len(n_rows)) {
        cols <- seq_len(n_cols)
        if (serpentine && mod(r, 2) == 0) cols <- rev(cols)
        for (c in cols) {
          positions[[kpos]] <- c(Row = r, Column = c)
          kpos <- kpos + 1
        }
      }
    } else {
      for (c in seq_len(n_cols)) {
        rows <- seq_len(n_rows)
        if (serpentine && mod(c, 2) == 0) rows <- rev(rows)
        for (r in rows) {
          positions[[kpos]] <- c(Row = r, Column = c)
          kpos <- kpos + 1
        }
      }
    }
    
    pos_mat <- as.data.frame(do.call(rbind, positions))
    pos_mat$Row <- as.integer(pos_mat$Row)
    pos_mat$Column <- as.integer(pos_mat$Column)
    pos_mat
  }
  
  # ============================================================
  # 5. CHECK PLACEMENT STRATEGY (INCLUDING OPTIMIZATION)
  # ============================================================
  if (check_placement != "optimal") {
    blocks <- build_blocks_once(seed_offset = 0)
  } else {
    best_score <- -Inf
    best_blocks <- NULL
    
    for (k in seq_len(check_opt_attempts)) {
      cand_blocks <- build_blocks_once(seed_offset = k)
      fb_cand <- do.call(rbind, cand_blocks)
      rownames(fb_cand) <- paste0("plot_", seq_len(nrow(fb_cand)))
      
      n_assigned_c <- nrow(fb_cand)
      pos_mat <- build_positions(n_rows, n_cols, order, serpentine)
      fb_cand$Row <- pos_mat$Row[seq_len(n_assigned_c)]
      fb_cand$Column <- pos_mat$Column[seq_len(n_assigned_c)]
      
      is_chk <- fb_cand$Treatment %in% check_treatments
      chk_rc <- unique(cbind(fb_cand$Row[is_chk], fb_cand$Column[is_chk]))
      
      if (nrow(chk_rc) <= 1) {
        score <- 0
      } else {
        d <- as.matrix(stats::dist(chk_rc))
        diag(d) <- Inf
        score <- mean(apply(d, 1, min))
      }
      
      if (score > best_score) {
        best_score <- score
        best_blocks <- cand_blocks
      }
    }
    blocks <- best_blocks
  }
  
  final_data <- do.call(rbind, blocks)
  n_assigned <- nrow(final_data)
  rownames(final_data) <- paste0("plot_", seq_len(n_assigned))
  
  pos_mat <- build_positions(n_rows, n_cols, order, serpentine)
  final_data$Plot <- seq_len(n_assigned)
  final_data$Row <- pos_mat$Row[seq_len(n_assigned)]
  final_data$Column <- pos_mat$Column[seq_len(n_assigned)]
  
  layout_matrix <- matrix(NA, nrow = n_rows, ncol = n_cols)
  for (i in seq_len(n_assigned)) {
    layout_matrix[final_data$Row[i], final_data$Column[i]] <- final_data$Treatment[i]
  }
  
  # ============================================================
  # 6. OPTIONAL GENETIC DISPERSION LOCAL SEARCH
  # ============================================================
  build_neighbor_pairs <- function(row, col, radius = 1) {
    n <- length(row)
    if (n < 2) return(matrix(integer(0), ncol = 2, dimnames = list(NULL, c("i", "j"))))
    
    pairs <- vector("list", 0)
    k <- 1
    for (i in seq_len(n - 1)) {
      dr <- abs(row[i] - row[(i + 1):n])
      dc <- abs(col[i] - col[(i + 1):n])
      ok <- pmax(dr, dc) <= radius
      if (any(ok)) {
        jj <- (i + 1):n
        jj <- jj[ok]
        pairs[[k]] <- cbind(i = rep(i, length(jj)), j = jj)
        k <- k + 1
      }
    }
    if (length(pairs) == 0) return(matrix(integer(0), ncol = 2, dimnames = list(NULL, c("i", "j"))))
    out <- do.call(rbind, pairs)
    colnames(out) <- c("i", "j")
    out
  }
  
  score_dispersion <- function(trt_vec, is_check, Ksub, pairs) {
    if (nrow(pairs) == 0) return(0)
    
    line_levels <- colnames(Ksub)
    line_idx <- rep(NA_integer_, length(trt_vec))
    line_idx[!is_check] <- match(trt_vec[!is_check], line_levels)
    
    ii <- pairs[, "i"]; jj <- pairs[, "j"]
    li <- line_idx[ii]; lj <- line_idx[jj]
    ok <- !is.na(li) & !is.na(lj)
    if (!any(ok)) return(0)
    
    sum(Ksub[cbind(li[ok], lj[ok])])
  }
  
  apply_genetic_dispersion <- function(
    fb,
    check_treatments,
    Kdisp,
    line_id_map,
    radius,
    iters,
    seed_local
  ) {
    .with_local_seed(seed_local, {
      trt <- as.character(fb$Treatment)
      is_check <- trt %in% check_treatments
      movable <- which(!is_check)
      
      if (length(movable) < 2 || iters <= 0) return(fb)
      
      non_trt <- unique(trt[!is_check])
      
      if (is.null(line_id_map)) {
        line_ids <- setNames(non_trt, non_trt)
      } else {
        if (!is.data.frame(line_id_map) || !all(c("Treatment", "LineID") %in% names(line_id_map))) {
          stop("line_id_map must be a data.frame with columns: Treatment, LineID")
        }
        line_ids <- setNames(line_id_map$LineID, line_id_map$Treatment)
      }
      
      ids <- unname(line_ids[non_trt])
      miss <- setdiff(ids, rownames(Kdisp))
      if (length(miss) > 0) stop("Some non-check LineIDs are missing in dispersion matrix rownames/colnames.")
      
      Ksub <- Kdisp[ids, ids, drop = FALSE]
      colnames(Ksub) <- non_trt
      rownames(Ksub) <- non_trt
      
      pairs <- build_neighbor_pairs(fb$Row, fb$Column, radius = radius)
      best_trt <- trt
      best_score <- score_dispersion(best_trt, is_check, Ksub, pairs)
      
      for (it in seq_len(iters)) {
        ij <- sample(movable, 2, replace = FALSE)
        cand <- best_trt
        cand[ij] <- cand[rev(ij)]
        sc <- score_dispersion(cand, is_check, Ksub, pairs)
        if (sc < best_score) {
          best_score <- sc
          best_trt <- cand
        }
      }
      
      fb$Treatment <- best_trt
      fb
    })
  }
  
  if (isTRUE(use_dispersion)) {
    Kdisp <- switch(
      dispersion_source,
      "K" = K,
      "A" = A,
      "GRM" = GRM
    )
    
    if (is.null(Kdisp)) stop("use_dispersion=TRUE but selected dispersion matrix is NULL.")
    if (is.null(rownames(Kdisp)) || is.null(colnames(Kdisp))) stop("Dispersion matrix must have rownames and colnames.")
    
    seed_disp_used <- dispersion_seed
    if (is.null(seed_disp_used)) seed_disp_used <- seed_used
    
    final_data <- apply_genetic_dispersion(
      fb = final_data,
      check_treatments = check_treatments,
      Kdisp = Kdisp,
      line_id_map = line_id_map,
      radius = dispersion_radius,
      iters = dispersion_iters,
      seed_local = seed_disp_used
    )
    
    layout_matrix[,] <- NA
    for (i in seq_len(n_assigned)) {
      layout_matrix[final_data$Row[i], final_data$Column[i]] <- final_data$Treatment[i]
    }
  }
  
  # ============================================================
  # 7. OPTIONAL EFFICIENCY EVALUATION
  # ============================================================
  efficiency <- NULL
  
  if (isTRUE(eval_efficiency) && (treatment_effect == "fixed" || prediction_type != "none")) {
    
    if (!requireNamespace("Matrix", quietly = TRUE)) {
      stop("Package 'Matrix' is required for efficiency evaluation.")
    }
    
    ar1_precision_sparse <- function(nn, rho) {
      if (nn <= 0) stop("n must be >= 1")
      if (nn == 1) return(Matrix::Diagonal(1, 1))
      if (abs(rho) >= 1) stop("AR1 requires |rho| < 1")
      
      a <- 1 / (1 - rho^2)
      d <- rep((1 + rho^2) * a, nn)
      d[1] <- a
      d[nn] <- a
      o <- rep(-rho * a, nn - 1)
      
      Matrix::sparseMatrix(
        i = c(seq_len(nn), seq_len(nn - 1), 2:nn),
        j = c(seq_len(nn), 2:nn, seq_len(nn - 1)),
        x = c(d, o, o),
        dims = c(nn, nn)
      )
    }
    
    make_sparse_incidence <- function(levels_vec) {
      nn <- length(levels_vec)
      lv <- unique(levels_vec[!is.na(levels_vec)])
      if (length(lv) == 0) {
        return(list(M = Matrix::Matrix(0, nrow = nn, ncol = 0, sparse = TRUE), levels = character(0)))
      }
      j <- match(levels_vec, lv)
      ok <- !is.na(j)
      M <- Matrix::sparseMatrix(i = which(ok), j = j[ok], x = 1, dims = c(nn, length(lv)))
      colnames(M) <- lv
      list(M = M, levels = lv)
    }
    
    pinv_sym_dense <- function(A, tol = 1e-10) {
      eg <- eigen(A, symmetric = TRUE)
      vals <- eg$values
      vecs <- eg$vectors
      keep <- vals > tol
      if (!any(keep)) return(matrix(0, nrow(A), ncol(A)))
      vecs[, keep, drop = FALSE] %*% diag(1 / vals[keep]) %*% t(vecs[, keep, drop = FALSE])
    }
    
    solve_C <- function(Cmat, B) {
      out <- try({
        fac <- Matrix::Cholesky(Cmat, LDL = TRUE, Imult = 0)
        Matrix::solve(fac, B)
      }, silent = TRUE)
      if (inherits(out, "try-error")) out <- Matrix::solve(Cmat, B)
      out
    }
    
    trace_subinv_est <- function(Cmat, idx, m = 80, seed_local = 1) {
      .with_local_seed(seed_local, {
        p <- length(idx)
        if (p == 0) return(NA_real_)
        nn <- nrow(Cmat)
        acc <- 0
        for (k in seq_len(m)) {
          z <- sample(c(-1, 1), p, replace = TRUE)
          u <- Matrix::sparseVector(i = idx, x = z, length = nn)
          x <- solve_C(Cmat, u)
          acc <- acc + as.numeric(Matrix::crossprod(u, x))
        }
        acc / m
      })
    }
    
    safe_logdet_psd_dense <- function(A, tol = 1e-10) {
      eg <- eigen(A, symmetric = TRUE)
      lam <- eg$values
      lam_pos <- lam[lam > tol]
      if (length(lam_pos) == 0) return(-Inf)
      sum(log(lam_pos))
    }
    
    ar1_cov_dense <- function(nn, rho) {
      idx <- seq_len(nn)
      outer(idx, idx, function(i, j) rho^abs(i - j))
    }
    
    pairwise_diff_mean_var <- function(V) {
      p <- nrow(V)
      if (p < 2) return(NA_real_)
      one <- rep(1, p)
      num <- p * sum(diag(V)) - as.numeric(t(one) %*% V %*% one)
      (2 * num) / (p * (p - 1))
    }
    
    fb <- final_data
    nn <- nrow(fb)
    
    trt <- as.character(fb$Treatment)
    blk <- as.character(fb$Block)
    rw  <- as.character(fb$Row)
    cl  <- as.character(fb$Column)
    
    is_check <- trt %in% check_treatments
    
    sigma_e2 <- varcomp$sigma_e2
    sigma_g2 <- varcomp$sigma_g2
    sigma_b2 <- varcomp$sigma_b2
    sigma_r2 <- varcomp$sigma_r2
    sigma_c2 <- varcomp$sigma_c2
    
    if (any(c(sigma_e2, sigma_b2, sigma_r2, sigma_c2) <= 0)) {
      stop("Variance components sigma_e2, sigma_b2, sigma_r2, sigma_c2 must be > 0.")
    }
    if (treatment_effect == "random" && sigma_g2 <= 0) {
      stop("sigma_g2 must be > 0 when treatment_effect='random'.")
    }
    
    has_holes <- n_assigned < (n_rows * n_cols)
    
    spatial_engine_use <- spatial_engine
    if (spatial_engine_use == "auto") {
      spatial_engine_use <- if (nn <= dense_max_n) "dense" else "sparse"
    }
    
    eff_notes <- character(0)
    residual_use <- residual_structure
    if (residual_structure != "IID" && has_holes) {
      eff_notes <- c(eff_notes, "Spatial residual requested and layout has holes; AR1 submatrix on observed plots was used.")
    }
    
    if (residual_use == "IID") {
      Q <- (1 / sigma_e2) * Matrix::Diagonal(nn, 1)
      
    } else if (spatial_engine_use == "dense") {
      
      if (residual_use == "AR1") {
        Rrow <- ar1_cov_dense(n_rows, rho_row)
        Rcol <- diag(n_cols)
      } else if (residual_use == "AR1xAR1") {
        Rrow <- ar1_cov_dense(n_rows, rho_row)
        Rcol <- ar1_cov_dense(n_cols, rho_col)
      } else stop("Unsupported residual_structure.")
      
      Rgrid <- sigma_e2 * kronecker(Rcol, Rrow)
      grid_index <- (as.integer(rw) - 1) * n_cols + as.integer(cl)
      Robs <- Rgrid[grid_index, grid_index, drop = FALSE]
      Q <- Matrix::Matrix(solve(Robs), sparse = FALSE)
      
    } else {
      
      Qrow <- ar1_precision_sparse(n_rows, rho_row)
      if (residual_use == "AR1") {
        Qcol <- Matrix::Diagonal(n_cols, 1)
      } else if (residual_use == "AR1xAR1") {
        Qcol <- ar1_precision_sparse(n_cols, rho_col)
      } else stop("Unsupported residual_structure.")
      
      Qgrid <- Matrix::kronecker(Qcol, Qrow) * (1 / sigma_e2)
      grid_index <- (as.integer(rw) - 1) * n_cols + as.integer(cl)
      Q <- Qgrid[grid_index, grid_index, drop = FALSE]
    }
    
    X_int <- Matrix::Matrix(rep(1, nn), ncol = 1, sparse = TRUE)
    colnames(X_int) <- "(Intercept)"
    X <- X_int
    
    if (isTRUE(check_as_fixed)) {
      chk_vec <- ifelse(is_check, trt, NA)
      chk_inc <- make_sparse_incidence(chk_vec)
      if (ncol(chk_inc$M) > 0) {
        keep <- intersect(check_treatments, colnames(chk_inc$M))
        chkM <- chk_inc$M[, keep, drop = FALSE]
        chkM <- chkM[, check_treatments[check_treatments %in% keep], drop = FALSE]
        colnames(chkM) <- paste0("Check_", colnames(chkM))
        X <- cbind(X, chkM)
      }
    }
    
    if (treatment_effect == "fixed") {
      trt_fix_vec <- ifelse(!is_check, trt, NA)
      trt_inc <- make_sparse_incidence(trt_fix_vec)
      if (ncol(trt_inc$M) < 2) stop("Not enough fixed non-check treatments to compute efficiency.")
      colnames(trt_inc$M) <- paste0("Line_", colnames(trt_inc$M))
      X <- cbind(X, trt_inc$M)
    }
    
    Zb <- make_sparse_incidence(blk)$M
    Zr <- make_sparse_incidence(rw)$M
    Zc <- make_sparse_incidence(cl)$M
    colnames(Zb) <- paste0("Block_", colnames(Zb))
    colnames(Zr) <- paste0("Row_",   colnames(Zr))
    colnames(Zc) <- paste0("Col_",   colnames(Zc))
    
    Z_list <- list(Zb = Zb, Zr = Zr, Zc = Zc)
    
    Zg <- NULL
    if (treatment_effect == "random") {
      g_vec <- ifelse(!is_check, trt, NA)
      g_inc <- make_sparse_incidence(g_vec)
      Zg <- g_inc$M
      if (ncol(Zg) < 2) stop("Not enough random non-check treatments to compute efficiency.")
      colnames(Zg) <- paste0("Line_", colnames(Zg))
      Z_list$Zg <- Zg
    }
    
    Z <- do.call(cbind, Z_list)
    
    XtQX <- Matrix::crossprod(X, Q %*% X)
    XtQZ <- Matrix::crossprod(X, Q %*% Z)
    ZtQX <- t(XtQZ)
    ZtQZ <- Matrix::crossprod(Z, Q %*% Z)
    
    Ginv <- Matrix::Diagonal(ncol(Z), 0)
    
    idx0 <- 0
    pB <- ncol(Zb); pR <- ncol(Zr); pC <- ncol(Zc)
    
    if (pB > 0) {
      ii <- (idx0 + 1):(idx0 + pB)
      Ginv[ii, ii] <- Matrix::Diagonal(pB, 1 / sigma_b2)
      idx0 <- idx0 + pB
    }
    if (pR > 0) {
      ii <- (idx0 + 1):(idx0 + pR)
      Ginv[ii, ii] <- Matrix::Diagonal(pR, 1 / sigma_r2)
      idx0 <- idx0 + pR
    }
    if (pC > 0) {
      ii <- (idx0 + 1):(idx0 + pC)
      Ginv[ii, ii] <- Matrix::Diagonal(pC, 1 / sigma_c2)
      idx0 <- idx0 + pC
    }
    
    trt_idx <- integer(0)
    
    if (treatment_effect == "random") {
      pG <- ncol(Zg)
      trt_idx <- (idx0 + 1):(idx0 + pG)
      
      if (prediction_type == "IID") {
        Ginv[trt_idx, trt_idx] <- Matrix::Diagonal(pG, 1 / sigma_g2)
      } else {
        if (is.null(K)) stop("K must be provided for prediction_type='GBLUP'/'PBLUP'.")
        if (is.null(rownames(K)) || is.null(colnames(K))) stop("K must have rownames and colnames.")
        
        line_trt_order <- sub("^Line_", "", colnames(Zg))
        
        if (is.null(line_id_map)) {
          line_ids <- setNames(line_trt_order, line_trt_order)
        } else {
          if (!is.data.frame(line_id_map) || !all(c("Treatment", "LineID") %in% names(line_id_map))) {
            stop("line_id_map must be a data.frame with columns: Treatment, LineID")
          }
          line_ids <- setNames(line_id_map$LineID, line_id_map$Treatment)
        }
        
        ids <- unname(line_ids[line_trt_order])
        miss <- setdiff(ids, rownames(K))
        if (length(miss) > 0) stop("Some line IDs are missing in K rownames/colnames.")
        
        Ksub2 <- K[ids, ids, drop = FALSE]
        
        Kinv_try <- try({
          Km <- Matrix::Matrix(Ksub2, sparse = FALSE)
          facK <- Matrix::Cholesky(Km, LDL = FALSE, Imult = 0)
          Matrix::solve(facK, Matrix::Diagonal(nrow(Km), 1))
        }, silent = TRUE)
        
        if (!inherits(Kinv_try, "try-error")) {
          Kinv <- Kinv_try
        } else {
          if (nrow(Ksub2) > 2500) warning("K is large and not Cholesky-factorable; dense generalized inverse may be slow.")
          Kinv <- Matrix::Matrix(pinv_sym_dense(Ksub2), sparse = FALSE)
        }
        
        Ginv[trt_idx, trt_idx] <- (1 / sigma_g2) * Kinv
      }
    }
    
    Cmat <- rbind(cbind(XtQX, XtQZ), cbind(ZtQX, ZtQZ + Ginv))
    if (spatial_engine_use == "dense") {
      Cmat <- as.matrix(Cmat)
    } else {
      Cmat <- Matrix::Matrix(Cmat, sparse = TRUE)
    }
    
    eff <- list(
      model = "mu + Treatment + Block + Row + Column + e",
      treatment_effect = treatment_effect,
      prediction_type = if (treatment_effect == "random") prediction_type else NA_character_,
      residual_structure_requested = residual_structure,
      residual_structure_used = residual_use,
      spatial_engine_used = spatial_engine_use,
      notes = eff_notes
    )
    
    if (treatment_effect == "fixed") {
      
      xnames <- colnames(X)
      trt_cols <- grep("^Line_", xnames)
      if (length(trt_cols) < 2) stop("Not enough fixed non-check treatments to compute efficiency.")
      p <- length(trt_cols)
      
      if (p <= eff_full_max) {
        B <- Matrix::sparseMatrix(i = trt_cols, j = seq_along(trt_cols), x = 1, dims = c(nrow(Cmat), p))
        Xsol <- solve_C(Cmat, B)
        Vsub <- as.matrix(Xsol[trt_cols, , drop = FALSE])
        
        mean_var_diff <- pairwise_diff_mean_var(Vsub)
        A_opt <- 1 / mean_var_diff
        
        H <- diag(p) - matrix(1 / p, p, p)
        Vctr <- H %*% Vsub %*% H
        logdet <- safe_logdet_psd_dense(Vctr)
        D_opt <- if (is.finite(logdet)) exp(-logdet / (p - 1)) else NA_real_
        
        eff$mode <- "FIXED_TREATMENT_BLUE_CONTRAST"
        eff$A <- A_opt
        eff$D <- D_opt
        eff$mean_VarDiff <- mean_var_diff
        eff$n_trt <- p
        
      } else {
        tr_est <- trace_subinv_est(Cmat, trt_cols, m = eff_trace_samples, seed_local = 1)
        mean_var_coef <- tr_est / p
        A_opt <- 1 / mean_var_coef
        
        eff$mode <- "FIXED_TREATMENT_BLUE_APPROX"
        eff$A <- A_opt
        eff$D <- NA_real_
        eff$mean_VarCoef <- mean_var_coef
        eff$mean_Var <- mean_var_coef
        eff$n_trt <- p
        eff$notes <- c(
          eff$notes,
          paste0(
            "Fixed-treatment target dimension (", p,
            ") > eff_full_max (", eff_full_max,
            "); A via stochastic trace estimator; D = NA."
          )
        )
      }
      
    } else {
      
      p <- length(trt_idx)
      if (p < 2) stop("Not enough random non-check treatments to compute efficiency.")
      
      if (p <= eff_full_max) {
        B <- Matrix::sparseMatrix(i = trt_idx, j = seq_along(trt_idx), x = 1, dims = c(nrow(Cmat), p))
        Xsol <- solve_C(Cmat, B)
        PEVsub <- as.matrix(Xsol[trt_idx, , drop = FALSE])
        
        A_pred <- 1 / mean(diag(PEVsub))
        logdet <- safe_logdet_psd_dense(PEVsub)
        D_pred <- if (is.finite(logdet)) exp(-logdet / p) else NA_real_
        
        eff$mode <- paste0("RANDOM_TREATMENT_BLUP_", prediction_type)
        eff$A <- A_pred
        eff$D <- D_pred
        eff$mean_PEV <- mean(diag(PEVsub))
        eff$n_lines <- p
        
      } else {
        tr_est <- trace_subinv_est(Cmat, trt_idx, m = eff_trace_samples, seed_local = 1)
        mean_pev <- tr_est / p
        A_pred <- 1 / mean_pev
        
        eff$mode <- paste0("RANDOM_TREATMENT_BLUP_", prediction_type, "_APPROX")
        eff$A <- A_pred
        eff$D <- NA_real_
        eff$mean_PEV <- mean_pev
        eff$n_lines <- p
        eff$notes <- c(
          eff$notes,
          paste0(
            "Random-treatment target dimension (", p,
            ") > eff_full_max (", eff_full_max,
            "); A via stochastic trace estimator; D = NA."
          )
        )
      }
    }
    
    efficiency <- eff
  }
  
  # ============================================================
  # 8. RETURN OBJECT
  # ============================================================
  list(
    layout_matrix = layout_matrix,
    field_book = final_data,
    efficiency = efficiency,
    seed_used = seed_used
  )
}
