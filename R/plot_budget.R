#' Test-entry capacity of an offered plot budget
#'
#' @description
#' Field partners offer space in **total plots**, but the allocation capacity
#' (`n_test_entries_per_environment`) and [suggest_site_capacity()] count **test
#' entries**, which exclude the check plots and replication a within-environment
#' design adds. Because checks appear in every block, a design uses
#' \eqn{b \times c} check plots (blocks times checks) plus the replicated and
#' unreplicated entry plots. `test_entry_capacity()` returns the largest number
#' of distinct test entries that fit in `total_plots`.
#'
#' @param total_plots Total plots offered at the site.
#' @param n_checks Number of check treatments (each appears in every block).
#' @param n_blocks Number of blocks (checks contribute `n_checks * n_blocks`
#'   plots).
#' @param avg_reps_per_entry Average plots per test entry: 1 for an
#'   augmented/unreplicated design, higher for p-rep or replicated entries
#'   (e.g. 1.3 if 30 percent of entries are duplicated). Default 1.
#' @return Integer number of distinct test entries.
#' @seealso [required_plots()], [field_plot_accounting()],
#'   [suggest_site_capacity()], [allocate_sparse_met()].
#' @examples
#' test_entry_capacity(400, n_checks = 2, n_blocks = 4)   # 2 checks x 4 blocks -> 392
#' @export
test_entry_capacity <- function(total_plots, n_checks = 0, n_blocks = 1,
                                avg_reps_per_entry = 1) {
  if (avg_reps_per_entry <= 0) stop("`avg_reps_per_entry` must be positive.")
  check_plots <- n_checks * n_blocks
  entry_plots <- total_plots - check_plots
  if (entry_plots < 1)
    stop("Check plots (n_checks * n_blocks = ", check_plots,
         ") already meet or exceed `total_plots`.")
  as.integer(floor(entry_plots / avg_reps_per_entry))
}


#' Total plots required for a number of test entries
#'
#' @description
#' The inverse of [test_entry_capacity()]: the total plots a design needs given
#' its test entries, check plots (`n_checks * n_blocks`), and replication.
#'
#' @param n_test_entries Number of distinct test entries (including common
#'   treatments).
#' @param n_checks,n_blocks,avg_reps_per_entry See [test_entry_capacity()].
#' @return Integer total plots required.
#' @seealso [test_entry_capacity()], [field_plot_accounting()].
#' @examples
#' required_plots(392, n_checks = 2, n_blocks = 4)   # -> 400
#' @export
required_plots <- function(n_test_entries, n_checks = 0, n_blocks = 1,
                           avg_reps_per_entry = 1) {
  as.integer(n_checks * n_blocks + ceiling(n_test_entries * avg_reps_per_entry))
}


#' Field plot accounting (checks, entries, replication)
#'
#' @description
#' Returns the full plot breakdown of a site's field, reconciling total plots
#' with test entries, check plots, and replication. Supply either `total_plots`
#' or `n_test_entries`.
#'
#' @param total_plots Total plots offered (optional if `n_test_entries` given).
#' @param n_test_entries Number of distinct test entries (optional if
#'   `total_plots` given).
#' @param n_checks,n_blocks,avg_reps_per_entry See [test_entry_capacity()].
#' @return A one-row data frame with `total_plots`, `check_plots`, `entry_plots`,
#'   `n_test_entries`, and `avg_reps_per_entry`.
#' @seealso [test_entry_capacity()], [required_plots()].
#' @examples
#' field_plot_accounting(total_plots = 400, n_checks = 2, n_blocks = 4)
#' @export
field_plot_accounting <- function(total_plots = NULL, n_test_entries = NULL,
                                  n_checks = 0, n_blocks = 1,
                                  avg_reps_per_entry = 1) {
  if (is.null(total_plots) && is.null(n_test_entries))
    stop("Supply `total_plots` or `n_test_entries`.")
  check_plots <- n_checks * n_blocks
  if (is.null(n_test_entries))
    n_test_entries <- test_entry_capacity(total_plots, n_checks, n_blocks,
                                          avg_reps_per_entry)
  entry_plots <- ceiling(n_test_entries * avg_reps_per_entry)
  if (is.null(total_plots)) total_plots <- check_plots + entry_plots
  data.frame(total_plots = total_plots, check_plots = check_plots,
             entry_plots = entry_plots, n_test_entries = n_test_entries,
             avg_reps_per_entry = avg_reps_per_entry,
             stringsAsFactors = FALSE)
}
