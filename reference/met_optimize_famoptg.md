# Search for a criterion-optimal repeated-check block design

`met_optimize_famoptg()` is the OptiSparseMET version of
`optimize_famoptg()` from the OptiDesign package. The `met_` prefix
avoids namespace conflicts when both packages are loaded simultaneously.
All arguments, return values, and internal logic are identical to
`optimize_famoptg()` in OptiDesign.

`met_optimize_famoptg()` wraps
[`met_prep_famoptg()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_prep_famoptg.md)
and
[`met_evaluate_famoptg_efficiency()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_evaluate_famoptg_efficiency.md)
in a Random Restart (RS) optimisation loop that searches for the design
with the best optimality criterion value among `n_restarts` independent
randomisations.

**Why Random Restart only?**
[`met_prep_famoptg()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_prep_famoptg.md)
enforces the p-rep constraint - that no replicated treatment appears
twice in the same block

- by construction at every call. Permutation-based methods such as
  Simulated Annealing or a Genetic Algorithm would require block-aware
  swap logic to preserve this constraint, making them substantially more
  complex for modest criterion improvement. RS generates fully valid
  designs at every restart with no risk of constraint violation.

**Constraint preservation**: every candidate is produced by
[`met_prep_famoptg()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_prep_famoptg.md),
so the following structural guarantees hold by construction for all
candidates, regardless of `n_restarts`:

- Check treatments appear in every block.

- P-rep treatments appear in exactly `p_rep_reps[i]` blocks each, always
  in distinct blocks - never twice in the same block.

- Unreplicated treatments appear exactly once.

- Family/cluster adjacency is minimised within blocks.

- Optional genetic dispersion is applied.

**Integrity guarantee**: every candidate additionally passes
`.check_famoptg_integrity()` before being scored or stored as the best.
The running best is updated only when both the criterion score improves
and an integrity re-check passes. A final integrity check is performed
before returning. If no valid design is found an emergency fallback
returns a single freshly constructed valid design.

## Usage

``` r
met_optimize_famoptg(
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
  check_placement = c("random", "systematic", "optimal"),
  check_opt_attempts = 200,
  use_dispersion = FALSE,
  dispersion_source = c("K", "A", "GRM"),
  dispersion_radius = 1,
  dispersion_iters = 2000,
  dispersion_seed = 1,
  K = NULL,
  line_id_map = NULL,
  treatment_effect = c("random", "fixed"),
  prediction_type = c("IID", "GBLUP", "PBLUP", "none"),
  check_as_fixed = TRUE,
  residual_structure = c("IID", "AR1", "AR1xAR1"),
  rho_row = 0,
  rho_col = 0,
  varcomp = list(sigma_e2 = 1, sigma_g2 = 1, sigma_b2 = 1, sigma_r2 = 1, sigma_c2 = 1),
  spatial_engine = c("auto", "sparse", "dense"),
  dense_max_n = 5000,
  eff_trace_samples = 80,
  eff_full_max = 400,
  criterion = c("A", "D", "both", "CDmean"),
  n_restarts = 50L,
  max_failure_rate = 0.5,
  verbose_opt = TRUE
)
```

## Arguments

- check_treatments, check_families, p_rep_treatments, p_rep_reps,
  p_rep_families, unreplicated_treatments, unreplicated_families,
  n_blocks, n_rows, n_cols, order, serpentine, seed, attempts,
  warn_and_correct, fix_rows, cluster_source, GRM, A, id_map,
  cluster_method, cluster_seed, cluster_attempts, n_pcs_use,
  check_placement, check_opt_attempts, use_dispersion,
  dispersion_source, dispersion_radius, dispersion_iters,
  dispersion_seed, K, line_id_map:

  Passed directly to
  [`met_prep_famoptg()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_prep_famoptg.md).
  See that function for full documentation.

- treatment_effect, prediction_type, check_as_fixed, residual_structure,
  rho_row, rho_col, varcomp, spatial_engine, dense_max_n,
  eff_trace_samples, eff_full_max:

  Passed directly to
  [`met_evaluate_famoptg_efficiency()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_evaluate_famoptg_efficiency.md).
  See that function for full documentation. Note that `varcomp` here
  uses `sigma_b2` (block variance) instead of `sigma_rep2` and
  `sigma_ib2`.

- criterion:

  Character. Optimality criterion to drive the search:

  `"A"`

  :   Minimise `A_criterion` (mean pairwise contrast variance for fixed
      effects; mean PEV for random effects). Valid for both
      `treatment_effect = "fixed"` and `"random"`. Recommended default.

  `"D"`

  :   Minimise `D_criterion` (geometric mean of contrast covariance
      eigenvalues). Fixed effects only. Falls back to `A_criterion` with
      a warning for random effects.

  `"both"`

  :   Minimise mean of `A_criterion` and `D_criterion`. Fixed effects
      only.

  `"CDmean"`

  :   Maximise CDmean (mean coefficient of determination for GEBV
      prediction; Rincent et al. 2012). Requires
      `treatment_effect = "random"` and
      `prediction_type %in% c("IID", "GBLUP", "PBLUP")`. Best suited for
      optimising designs for genomic selection.

- n_restarts:

  Positive integer. Number of independent random restarts. Each restart
  calls
  [`met_prep_famoptg()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_prep_famoptg.md)
  with a new seed and evaluates the resulting design. Default 50.

- max_failure_rate:

  Numeric in \\\[0, 1\]\\. Maximum tolerated fraction of restarts that
  may fail construction or integrity checking before the function stops
  with a diagnostic error. Default `0.5`. Below the threshold a warning
  is issued. Increase for very constrained field geometries where some
  failure is expected.

- verbose_opt:

  Logical. If `TRUE`, prints per-restart progress messages showing
  current score, running best, and failure status. Default `TRUE`.

## Value

The return value of
[`met_prep_famoptg()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_prep_famoptg.md)
for the best valid design found, augmented with two additional
components:

- `efficiency`:

  The full output of
  [`met_evaluate_famoptg_efficiency()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_evaluate_famoptg_efficiency.md)
  for the best design, including all applicable criterion values:
  `A_criterion`, `D_criterion`, `A_efficiency`, `D_efficiency`,
  `CDmean`, and `CD_per_line`.

- `optimization`:

  Named list of optimisation metadata:

  `method`

  :   Character. Always `"RS"`.

  `criterion`

  :   Character. As supplied.

  `best_score`

  :   Numeric. Best criterion value found, in natural direction: lower
      is better for `"A"`, `"D"`, `"both"`; higher is better for
      `"CDmean"` (positive value, negation is internal only).

  `score_history`

  :   Numeric vector of length `n_restarts`. Criterion value for each
      restart (`NA` for failed restarts). For `criterion = "CDmean"`,
      values are positive CDmean (higher = better). For all other
      criteria, lower = better.

  `master_seed`

  :   Integer. The master random seed used.

  `n_restarts`

  :   Integer. As supplied.

  `n_failed`

  :   Integer. Number of restarts that failed construction or integrity
      checking.

If no valid optimised design is found after all restarts, a warning is
issued and the function returns a single freshly constructed valid
design via emergency fallback. If even that fails after 10 attempts, the
function stops with a diagnostic error pointing to the likely cause.

## Details

### Optimisation target

All restarts minimise an internal score (lower = better internally):

|             |                                         |                            |
|-------------|-----------------------------------------|----------------------------|
| `criterion` | Internal score                          | Direction reported to user |
| `"A"`       | `A_criterion`                           | Lower is better            |
| `"D"`       | `D_criterion`                           | Lower is better            |
| `"both"`    | Mean of `A_criterion` and `D_criterion` | Lower is better            |
| `"CDmean"`  | Negated CDmean                          | Higher CDmean is better    |

For `"CDmean"` the internal negation is transparent to the user:
`best_score` and `score_history` are always reported as positive CDmean
values (higher = better).

### Integrity checking

Five structural constraints are verified for every candidate:

1.  No non-check entry appears more than once in any single block.

2.  Each p-rep treatment appears in exactly `p_rep_reps[i]` blocks.

3.  Each unreplicated treatment appears exactly once.

4.  All check treatments appear in every block.

5.  No p-rep treatment occupies the same block twice (core p-rep
    constraint).

Checks 2 and 5 together constitute the p-rep guarantee. Failure of any
check discards the candidate and counts toward `n_failed`.

### Failure handling

Failed restarts (construction error or integrity failure) are counted.
If the failure rate exceeds `max_failure_rate`, the function stops with
a diagnostic message. Below the threshold a warning is issued. If no
valid design is found after all restarts, up to 10 emergency fallback
attempts are made. If all fail, the function stops with an informative
error.

## References

Rincent, R., Laloe, D., Nicolas, S., Altmann, T., Brunel, D., Revilla,
P., ... & Moreau, L. (2012). Maximizing the reliability of genomic
selection by optimizing the calibration set of reference individuals.
*Genetics*, 192(2), 715-728.

Jones, B., Allen-Moyer, K., & Goos, P. (2021). A-optimal versus
D-optimal design of screening experiments. *Journal of Quality
Technology*, 53(4), 369-382.

## See also

[`met_prep_famoptg()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_prep_famoptg.md)
for single-design construction without optimisation.
[`met_evaluate_famoptg_efficiency()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_evaluate_famoptg_efficiency.md)
for standalone criterion evaluation.
[`met_optimize_alpha_rc()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_optimize_alpha_rc.md)
for the equivalent optimizer for
[`met_alpha_rc_stream()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_alpha_rc_stream.md)
designs (also supports SA and GA methods).

## Examples

``` r
## Optimise a p-rep design by A-criterion (fixed effects)
result <- met_optimize_famoptg(
  check_treatments        = c("CHK1", "CHK2"),
  check_families          = c("CHECK", "CHECK"),
  p_rep_treatments        = paste0("P", 1:20),
  p_rep_reps              = rep(2L, 20),
  p_rep_families          = rep(paste0("F", 1:4), 5),
  unreplicated_treatments = paste0("U", 1:60),
  unreplicated_families   = rep(paste0("F", 1:4), 15),
  n_blocks           = 5,
  n_rows             = 15,
  n_cols             = 20,
  treatment_effect   = "fixed",
  residual_structure = "AR1xAR1",
  rho_row            = 0.10,
  rho_col            = 0.10,
  criterion          = "A",
  n_restarts         = 50
)
#> 
#> [met_optimize_famoptg] criterion = 'A' | n_restarts = 50 | master seed = 787830449
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 1/50 | score = 3.152510 | best = 3.152510
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 2/50 | score = 3.170390 | best = 3.152510
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 3/50 | score = 3.087090 | best = 3.087090
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 4/50 | score = 3.040015 | best = 3.040015
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 5/50 | score = 3.150978 | best = 3.040015
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 6/50 | score = 3.127811 | best = 3.040015
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 7/50 | score = 3.111434 | best = 3.040015
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 8/50 | score = 3.073177 | best = 3.040015
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 9/50 | score = 3.169440 | best = 3.040015
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 10/50 | score = 3.122313 | best = 3.040015
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 11/50 | score = 3.119799 | best = 3.040015
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 12/50 | score = 3.086633 | best = 3.040015
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 13/50 | score = 3.119774 | best = 3.040015
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 14/50 | score = 3.107902 | best = 3.040015
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 15/50 | score = 3.115780 | best = 3.040015
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 16/50 | score = 3.132496 | best = 3.040015
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 17/50 | score = 3.124873 | best = 3.040015
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 18/50 | score = 3.121083 | best = 3.040015
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 19/50 | score = 3.079913 | best = 3.040015
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 20/50 | score = 3.111494 | best = 3.040015
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 21/50 | score = 3.090669 | best = 3.040015
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 22/50 | score = 3.058088 | best = 3.040015
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 23/50 | score = 3.140534 | best = 3.040015
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 24/50 | score = 3.212416 | best = 3.040015
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 25/50 | score = 3.062542 | best = 3.040015
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 26/50 | score = 3.202487 | best = 3.040015
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 27/50 | score = 3.103123 | best = 3.040015
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 28/50 | score = 3.027256 | best = 3.027256
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 29/50 | score = 3.148211 | best = 3.027256
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 30/50 | score = 3.241331 | best = 3.027256
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 31/50 | score = 3.149520 | best = 3.027256
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 32/50 | score = 3.144043 | best = 3.027256
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 33/50 | score = 3.098636 | best = 3.027256
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 34/50 | score = 3.299052 | best = 3.027256
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 35/50 | score = 3.080314 | best = 3.027256
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 36/50 | score = 3.075660 | best = 3.027256
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 37/50 | score = 3.068523 | best = 3.027256
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 38/50 | score = 3.222485 | best = 3.027256
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 39/50 | score = 2.997918 | best = 2.997918
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 40/50 | score = 3.019537 | best = 2.997918
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 41/50 | score = 3.169334 | best = 2.997918
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 42/50 | score = 3.203584 | best = 2.997918
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 43/50 | score = 3.060075 | best = 2.997918
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 44/50 | score = 3.154353 | best = 2.997918
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 45/50 | score = 3.302120 | best = 2.997918
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 46/50 | score = 3.049646 | best = 2.997918
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 47/50 | score = 3.109698 | best = 2.997918
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 48/50 | score = 3.083572 | best = 2.997918
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 49/50 | score = 3.145700 | best = 2.997918
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#>   [RS] 50/50 | score = 3.207284 | best = 2.997918
#> 
#> [met_optimize_famoptg] Done. Best A = 2.997918 (lower is better) | valid restarts = 50/50

result$efficiency$A_criterion      # best A-criterion found
#> [1] 2.997918
result$optimization$best_score     # same value, lower is better
#> [1] 2.997918
result$optimization$score_history  # per-restart scores
#>  [1] 3.152510 3.170390 3.087090 3.040015 3.150978 3.127811 3.111434 3.073177
#>  [9] 3.169440 3.122313 3.119799 3.086633 3.119774 3.107902 3.115780 3.132496
#> [17] 3.124873 3.121083 3.079913 3.111494 3.090669 3.058088 3.140534 3.212416
#> [25] 3.062542 3.202487 3.103123 3.027256 3.148211 3.241331 3.149520 3.144043
#> [33] 3.098636 3.299052 3.080314 3.075660 3.068523 3.222485 2.997918 3.019537
#> [41] 3.169334 3.203584 3.060075 3.154353 3.302120 3.049646 3.109698 3.083572
#> [49] 3.145700 3.207284
result$optimization$n_failed       # failed restarts
#> [1] 0

if (FALSE) { # \dontrun{
## Optimise an augmented design by CDmean (GBLUP)
result_cdmean <- met_optimize_famoptg(
  check_treatments        = c("CHK1", "CHK2"),
  check_families          = c("CHECK", "CHECK"),
  p_rep_treatments        = character(0),
  p_rep_reps              = integer(0),
  p_rep_families          = character(0),
  unreplicated_treatments = paste0("E", 1:120),
  unreplicated_families   = rep(paste0("F", 1:6), 20),
  n_blocks           = 6,
  n_rows             = 13,
  n_cols             = 10,
  treatment_effect   = "random",
  prediction_type    = "GBLUP",
  K                  = my_kinship_matrix,
  varcomp            = list(
    sigma_g2 = 0.4, sigma_e2 = 0.6,
    sigma_b2 = 0.1, sigma_r2 = 0.02, sigma_c2 = 0.02
  ),
  criterion          = "CDmean",
  n_restarts         = 50
)

result_cdmean$efficiency$CDmean        # mean GEBV reliability [0,1]
result_cdmean$optimization$best_score  # positive CDmean, higher = better
result_cdmean$optimization$n_failed    # failed restarts
} # }
```
