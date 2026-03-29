# Evaluate the statistical efficiency of a repeated-check block design

`met_evaluate_famoptg_efficiency()` is the OptiSparseMET version of
`evaluate_famoptg_efficiency()` from the OptiDesign package. The `met_`
prefix avoids namespace conflicts when both packages are loaded
simultaneously. All arguments, return values, and internal logic are
identical to `evaluate_famoptg_efficiency()` in OptiDesign.

`met_evaluate_famoptg_efficiency()` takes a `field_book` produced by
[`met_prep_famoptg()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_prep_famoptg.md)
and computes optimality criteria under a user-specified mixed model. It
is fully decoupled from design construction so that a single design can
be evaluated multiple times under different model assumptions without
rebuilding the layout.

This function is the sibling of
[`met_evaluate_alpha_efficiency()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_evaluate_alpha_efficiency.md)
for
[`met_alpha_rc_stream()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_alpha_rc_stream.md)
designs. The key structural difference is that
[`met_prep_famoptg()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_prep_famoptg.md)
designs have a single flat blocking structure (one set of `n_blocks`
blocks with no replicate or incomplete-block nesting), so the random
effect model here contains **Block + Row + Column** rather than the
**Rep + IBlock(Rep) + Row + Column** structure of alpha-lattice designs.

## Usage

``` r
met_evaluate_famoptg_efficiency(
  field_book,
  n_rows,
  n_cols,
  check_treatments,
  treatment_effect = c("random", "fixed"),
  prediction_type = c("IID", "GBLUP", "PBLUP", "none"),
  check_as_fixed = TRUE,
  residual_structure = c("IID", "AR1", "AR1xAR1"),
  rho_row = 0,
  rho_col = 0,
  varcomp = list(sigma_e2 = 1, sigma_g2 = 1, sigma_b2 = 1, sigma_r2 = 1, sigma_c2 = 1),
  K = NULL,
  line_id_map = NULL,
  spatial_engine = c("auto", "sparse", "dense"),
  dense_max_n = 5000,
  eff_trace_samples = 80,
  eff_full_max = 400
)
```

## Arguments

- field_book:

  Data frame produced by
  [`met_prep_famoptg()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_prep_famoptg.md).
  Must contain columns: `Treatment`, `Family`, `Block`, `Plot`, `Row`,
  `Column`.

- n_rows:

  Positive integer. Number of field rows (must match the original
  design).

- n_cols:

  Positive integer. Number of field columns (must match the original
  design).

- check_treatments:

  Character vector of check treatment identifiers (must match those used
  in
  [`met_prep_famoptg()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_prep_famoptg.md)).

- treatment_effect:

  Character. Whether entry treatments are modelled as `"fixed"`
  (BLUE-based A and D criteria) or `"random"` (PEV-based criterion and
  CDmean). Note: `"fixed"` does not produce CDmean; `"random"` does not
  produce D-criterion.

- prediction_type:

  Character. Random-effect prediction model:

  `"IID"`

  :   \\G^{-1}\_\text{entry} = \sigma_g^{-2} I\\. No relationship matrix
      required.

  `"GBLUP"` / `"PBLUP"`

  :   \\G^{-1}\_\text{entry} = \sigma_g^{-2} K^{-1}\\. Requires `K`.

  `"none"`

  :   Invalid when `treatment_effect = "random"`.

- check_as_fixed:

  Logical. If `TRUE`, checks are included as fixed effects in the model.
  Default `TRUE`.

- residual_structure:

  Character. Residual covariance structure: `"IID"`, `"AR1"` (row AR1
  only), or `"AR1xAR1"` (separable row x column AR1).

- rho_row:

  Numeric in \\(-1, 1)\\. AR1 autocorrelation along rows. Used when
  `residual_structure %in% c("AR1", "AR1xAR1")`.

- rho_col:

  Numeric in \\(-1, 1)\\. AR1 autocorrelation along columns. Used only
  when `residual_structure = "AR1xAR1"`.

- varcomp:

  Named list of variance components. Required components:

  `sigma_e2`

  :   Residual variance.

  `sigma_g2`

  :   Entry genetic variance. Used as denominator of CDmean when
      `treatment_effect = "random"`.

  `sigma_b2`

  :   Block variance. Corresponds to the single flat blocking level in
      [`met_prep_famoptg()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_prep_famoptg.md)
      designs (no replicate or incomplete-block nesting). This differs
      from
      [`met_evaluate_alpha_efficiency()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_evaluate_alpha_efficiency.md)
      which uses `sigma_rep2` and `sigma_ib2`.

  `sigma_r2`

  :   Row variance.

  `sigma_c2`

  :   Column variance.

- K:

  Optional square numeric matrix with rownames and colnames. Required
  when `prediction_type %in% c("GBLUP", "PBLUP")`.

- line_id_map:

  Optional data frame with columns `Treatment` and `LineID`. Required
  when treatment labels do not match `rownames(K)`.

- spatial_engine:

  Character. Computational engine: `"auto"` uses `"dense"` when used
  plots \\\leq\\ `dense_max_n`, otherwise `"sparse"`.

- dense_max_n:

  Integer. Threshold for `spatial_engine = "auto"`.

- eff_trace_samples:

  Positive integer. Number of Rademacher vectors for the Hutchinson
  stochastic trace estimator (large designs only).

- eff_full_max:

  Positive integer. Maximum number of treatments for exact \\C^{-1}\\
  submatrix extraction. Designs with more treatments use the stochastic
  estimator (`mode` suffix `_APPROX`).

## Value

A named list containing model metadata and criterion values:

- `model`:

  Character. Model string.

- `treatment_effect`:

  Character. As supplied.

- `prediction_type`:

  Character or `NA`.

- `residual_structure_requested`:

  Character. As supplied.

- `residual_structure_used`:

  Character. As resolved.

- `spatial_engine_used`:

  Character. `"dense"` or `"sparse"`.

- `mode`:

  Character. One of `"FIXED_TREATMENT_BLUE_CONTRAST"`,
  `"FIXED_TREATMENT_BLUE_APPROX"`, `"RANDOM_TREATMENT_PEV"`,
  `"RANDOM_TREATMENT_PEV_APPROX"`.

- `n_trt`:

  Integer. Number of treatment columns evaluated.

- `n_contrasts`:

  Integer. `n_trt - 1`. Fixed mode only.

- `A_criterion`:

  Numeric. Mean pairwise contrast variance (fixed) or mean PEV (random).
  Lower is better.

- `D_criterion`:

  Numeric or `NA`. Fixed full mode only. Lower is better.

- `A_efficiency`:

  Numeric. `1 / A_criterion`. Higher is better.

- `D_efficiency`:

  Numeric or `NA`. `1 / D_criterion`. Higher is better.

- `A`:

  Numeric. Alias for `A_efficiency` (backward compatibility).

- `D`:

  Numeric or `NA`. Alias for `D_efficiency` (backward compatibility).

- `mean_VarDiff`:

  Numeric. Mean pairwise contrast variance. Fixed full mode only.

- `PEV_criterion`:

  Numeric. Mean PEV. Random mode only.

- `mean_PEV`:

  Numeric. Alias for `PEV_criterion`. Random mode only.

- `CDmean`:

  Numeric. Mean coefficient of determination: \\1 - \text{mean\\PEV} /
  \sigma_g^2\\. Range \\\[0, 1\]\\. Higher is better. Random mode only.

- `CD_per_line`:

  Numeric vector or `NA`. Per-line coefficient of determination.
  Available in full mode only. Random mode only.

## Details

### Mixed model

\$\$y = X\beta + Zu + e\$\$

**Fixed part** \\X\\: intercept, optional check fixed effects
(`check_as_fixed = TRUE`), entry fixed effects when
`treatment_effect = "fixed"`.

**Random part** \\Z\\: block, row, column, and - when
`treatment_effect = "random"` - entry treatment effects. There is no
replicate or incomplete-block term because
[`met_prep_famoptg()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_prep_famoptg.md)
has no such nesting.

**Residual structure** \\R = \sigma_e^2 \Sigma\\:

- `"IID"`: \\\Sigma = I\\

- `"AR1"`: row-only AR1, \\Q\_\text{AR1}(\rho\_\text{row}) \otimes
  I\_\text{cols}\\

- `"AR1xAR1"`: separable row x column AR1,
  \\Q\_\text{AR1}(\rho\_\text{col}) \otimes
  Q\_\text{AR1}(\rho\_\text{row})\\

### Mixed model coefficient matrix

\$\$C = \begin{pmatrix} X^\top Q X & X^\top Q Z \\ Z^\top Q X & Z^\top Q
Z + G^{-1} \end{pmatrix}\$\$

where \\Q = R^{-1}\\ and \\G^{-1}\\ is block-diagonal with:

- Block random effects: \\\sigma_b^{-2} I\\

- Row random effects: \\\sigma_r^{-2} I\\

- Column random effects: \\\sigma_c^{-2} I\\

- Entry random effects (when `treatment_effect = "random"`):
  \\\sigma_g^{-2} I\\ (IID) or \\\sigma_g^{-2} K^{-1}\\ (GBLUP/PBLUP)

### Optimality criteria

#### Fixed treatment effects

**A-criterion** (lower = better): \$\$A\_\text{criterion} =
\bar{v}\_\text{diff} = \frac{2}{p(p-1)} \sum\_{i\<j}
\text{Var}(\hat{\tau}\_i - \hat{\tau}\_j)\$\$

**D-criterion** (lower = better): \$\$D\_\text{criterion} =
\exp\\\left(\frac{\log\det(HVH)}{p-1}\right)\$\$ where \\H = I_p -
p^{-1}J_p\\ and \\V\\ is the treatment variance-covariance submatrix of
\\C^{-1}\\.

**Efficiency forms** (higher = better): \$\$A\_\text{efficiency} = 1 /
A\_\text{criterion}, \quad D\_\text{efficiency} = 1 /
D\_\text{criterion}\$\$

#### Random treatment effects (genomic prediction)

**PEV-criterion** (lower = better): mean prediction error variance.
\$\$\text{PEV}\_\text{criterion} = \frac{1}{p} \sum\_{i=1}^{p}
\text{Var}(\hat{u}\_i - u_i)\$\$

**CDmean** (higher = better): mean coefficient of determination for GEBV
prediction (Rincent et al. 2012, *Genetics* 192:715-728):
\$\$\text{CDmean} = 1 -
\frac{\text{PEV}\_\text{criterion}}{\sigma_g^2}\$\$

**CD per line**: per-line coefficient of determination: \$\$CD_i = 1 -
\text{PEV}\_i / \sigma_g^2\$\$ Available in full mode only; `NA` in
`_APPROX` mode.

### Large-design approximation

When the number of treatments exceeds `eff_full_max`, a Hutchinson
stochastic trace estimator with `eff_trace_samples` Rademacher vectors
replaces exact submatrix extraction. The `mode` field carries the suffix
`_APPROX` and `D_criterion`, `D_efficiency`, and `CD_per_line` are `NA`.

## References

Rincent, R., Laloe, D., Nicolas, S., Altmann, T., Brunel, D., Revilla,
P., ... & Moreau, L. (2012). Maximizing the reliability of genomic
selection by optimizing the calibration set of reference individuals.
*Genetics*, 192(2), 715-728.

## See also

[`met_prep_famoptg()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_prep_famoptg.md)
to construct the design.
[`met_optimize_famoptg()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_optimize_famoptg.md)
to search for a criterion-optimal design.
[`met_evaluate_alpha_efficiency()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_evaluate_alpha_efficiency.md)
for the equivalent function for
[`met_alpha_rc_stream()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_alpha_rc_stream.md)
designs (different random effect structure: uses `sigma_rep2` and
`sigma_ib2` instead of `sigma_b2`).

## Examples

``` r
design <- met_prep_famoptg(
  check_treatments        = c("CHK1", "CHK2"),
  check_families          = c("CHECK", "CHECK"),
  p_rep_treatments        = paste0("P", 1:20),
  p_rep_reps              = rep(2L, 20),
  p_rep_families          = rep(paste0("F", 1:4), 5),
  unreplicated_treatments = paste0("U", 1:60),
  unreplicated_families   = rep(paste0("F", 1:4), 15),
  n_blocks = 5, n_rows = 15, n_cols = 20
)
#> Warning: Field size (15 x 20 = 300) does not match required plots (110). Adjusting dimensions.
#> Field = 15 x 8 | total plots = 120 | checks/block = 2 | p-rep entries = 20 | unreplicated = 60 | n_blocks = 5

## Fixed-effect evaluation (BLUE contrasts, IID residuals)
eff_fixed <- met_evaluate_famoptg_efficiency(
  field_book       = design$field_book,
  n_rows = 15, n_cols = 20,
  check_treatments = c("CHK1", "CHK2"),
  treatment_effect = "fixed"
)
eff_fixed$A_criterion   # lower is better
#> [1] 3.294417
eff_fixed$D_criterion   # lower is better
#> [1] 1.173881
eff_fixed$A_efficiency  # higher is better
#> [1] 0.3035438

## AR1xAR1 spatial model
eff_spatial <- met_evaluate_famoptg_efficiency(
  field_book         = design$field_book,
  n_rows = 15, n_cols = 20,
  check_treatments   = c("CHK1", "CHK2"),
  treatment_effect   = "fixed",
  residual_structure = "AR1xAR1",
  rho_row = 0.10, rho_col = 0.10
)
eff_spatial$A_criterion
#> [1] 3.212715

if (FALSE) { # \dontrun{
## Random-effect evaluation with GBLUP and CDmean
eff_gblup <- met_evaluate_famoptg_efficiency(
  field_book       = design$field_book,
  n_rows = 15, n_cols = 20,
  check_treatments = c("CHK1", "CHK2"),
  treatment_effect = "random",
  prediction_type  = "GBLUP",
  K                = my_kinship_matrix,
  varcomp          = list(
    sigma_g2 = 0.4, sigma_e2 = 0.6,
    sigma_b2 = 0.1, sigma_r2 = 0.02, sigma_c2 = 0.02
  )
)
eff_gblup$CDmean       # mean GEBV prediction reliability [0, 1]
eff_gblup$CD_per_line  # per-line reliability vector
eff_gblup$mean_PEV     # mean prediction error variance
} # }
```
