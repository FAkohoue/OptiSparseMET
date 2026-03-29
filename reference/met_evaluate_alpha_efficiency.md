# Evaluate the statistical efficiency of an alpha-lattice design

`met_evaluate_alpha_efficiency()` is the OptiSparseMET version of
`evaluate_alpha_efficiency()` from the OptiDesign package. The `met_`
prefix avoids namespace conflicts when both packages are loaded
simultaneously. All arguments, return values, and internal logic are
identical to `evaluate_alpha_efficiency()` in OptiDesign.

`met_evaluate_alpha_efficiency()` takes a `field_book` produced by
[`met_alpha_rc_stream()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_alpha_rc_stream.md)
and computes optimality criteria for the design under a user-specified
mixed model. It is fully decoupled from design construction, so a single
design can be evaluated multiple times under different model assumptions
without rebuilding the layout.

## Usage

``` r
met_evaluate_alpha_efficiency(
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
  varcomp = list(sigma_e2 = 1, sigma_g2 = 1, sigma_rep2 = 1, sigma_ib2 = 1, sigma_r2 = 1,
    sigma_c2 = 1),
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
  [`met_alpha_rc_stream()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_alpha_rc_stream.md).
  Must contain columns: `Plot`, `Row`, `Column`, `Rep`, `IBlock`,
  `BlockInRep`, `Treatment`, `Check`.

- n_rows:

  Positive integer. Number of field rows (must match the original
  design).

- n_cols:

  Positive integer. Number of field columns (must match the original
  design).

- check_treatments:

  Character vector of check treatment identifiers (must match those used
  in
  [`met_alpha_rc_stream()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_alpha_rc_stream.md)).

- treatment_effect:

  Character. Whether entry treatments are modelled as `"fixed"`
  (BLUE-based A and D criteria) or `"random"` (PEV-based criterion and
  CDmean). Note: `"fixed"` does not produce CDmean; `"random"` does not
  produce D-criterion.

- prediction_type:

  Character. Random-effect prediction model for entries. `"IID"` assumes
  independent entries (\\G^{-1} = \sigma_g^{-2} I\\). `"GBLUP"` and
  `"PBLUP"` use the supplied `K` matrix (\\G^{-1} = \sigma_g^{-2}
  K^{-1}\\). `"none"` is invalid when `treatment_effect = "random"` and
  raises an error.

- check_as_fixed:

  Logical. If `TRUE`, checks are included as fixed effects in the model.

- residual_structure:

  Character. Residual covariance structure. `"IID"` (independence),
  `"AR1"` (row AR1 only), or `"AR1xAR1"` (separable row x column AR1).

- rho_row:

  Numeric in \\(-1, 1)\\. AR1 autocorrelation parameter along rows. Used
  when `residual_structure %in% c("AR1", "AR1xAR1")`.

- rho_col:

  Numeric in \\(-1, 1)\\. AR1 autocorrelation parameter along columns.
  Used only when `residual_structure = "AR1xAR1"`.

- varcomp:

  Named list of variance components. Required components: `sigma_e2`
  (residual), `sigma_g2` (entry genetic), `sigma_rep2` (replicate),
  `sigma_ib2` (incomplete block), `sigma_r2` (row), `sigma_c2` (column).
  All components must be present even if the corresponding effect is
  absent from the model. `sigma_g2` is used as the denominator of CDmean
  when `treatment_effect = "random"`. Note: this function uses
  `sigma_rep2` and `sigma_ib2` for the Rep and IBlock(Rep) random
  effects, which differs from
  [`met_evaluate_famoptg_efficiency()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_evaluate_famoptg_efficiency.md)
  which uses `sigma_b2`.

- K:

  Optional square numeric matrix with rownames and colnames. Required
  when `prediction_type %in% c("GBLUP", "PBLUP")`.

- line_id_map:

  Optional data frame with columns `Treatment` and `LineID`. Required
  when treatment labels do not match `rownames(K)`.

- spatial_engine:

  Character. Computational engine for matrix operations. `"auto"` uses
  `"dense"` when used plots \\\leq\\ `dense_max_n`, otherwise
  `"sparse"`.

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

  Character. Computation mode: `"FIXED_TREATMENT_BLUE_CONTRAST"`,
  `"FIXED_TREATMENT_BLUE_APPROX"`, `"RANDOM_TREATMENT_PEV"`, or
  `"RANDOM_TREATMENT_PEV_APPROX"`.

- `n_trt`:

  Integer. Number of treatment columns evaluated.

- `n_contrasts`:

  Integer. `n_trt - 1`. Fixed mode only.

- `A_criterion`:

  Numeric. Mean pairwise contrast variance (fixed) or mean PEV (random).
  Lower is better.

- `D_criterion`:

  Numeric or `NA`. Geometric mean of contrast covariance eigenvalues.
  Fixed full mode only. Lower is better.

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

  Numeric. Mean coefficient of determination for GEBV prediction: \\1 -
  \text{mean\\PEV} / \sigma_g^2\\. Range \\\[0,1\]\\. Higher is better.
  Random mode only.

- `CD_per_line`:

  Numeric vector or `NA`. Per-line coefficient of determination: \\1 -
  \text{PEV}\_i / \sigma_g^2\\. Available in full mode only (`NA` in
  `_APPROX` mode). Random mode only.

## Details

### Mixed model

The model is: \$\$y = X\beta + Zu + e\$\$

**Fixed part** \\X\\: intercept, optional check fixed effects
(`check_as_fixed = TRUE`), entry treatment fixed effects when
`treatment_effect = "fixed"`.

**Random part** \\Z\\: replicate, incomplete block within replicate,
row, column, and - when `treatment_effect = "random"` - entry treatment
effects.

**Residual structure** \\R = \sigma_e^2 \Sigma\\:

- `"IID"`: \\\Sigma = I\\

- `"AR1"`: row-only AR1, \\Q\_\text{AR1}(\rho\_\text{row}) \otimes
  I\_\text{cols}\\

- `"AR1xAR1"`: separable row x column AR1,
  \\Q\_\text{AR1}(\rho\_\text{col}) \otimes
  Q\_\text{AR1}(\rho\_\text{row})\\

AR1 precision matrices are tridiagonal sparse matrices. For an
\\n\\-dimensional AR1 process with parameter \\\rho\\, let \\a =
1/(1-\rho^2)\\:

- Interior diagonal: \\(1+\rho^2) \times a\\

- Edge diagonal (positions 1 and \\n\\): \\a\\

- Off-diagonal: \\-\rho \times a\\

### Mixed model coefficient matrix

\$\$C = \begin{pmatrix} X^\top Q X & X^\top Q Z \\ Z^\top Q X & Z^\top Q
Z + G^{-1} \end{pmatrix}\$\$

where \\Q = R^{-1}\\ is the residual precision matrix.

### Random-effect structure G

- `"IID"`: \\G^{-1}\_\text{entry} = \sigma_g^{-2} I\\

- `"GBLUP"` / `"PBLUP"`: \\G^{-1}\_\text{entry} = \sigma_g^{-2}
  K^{-1}\\, computed via sparse Cholesky; falls back to Moore-Penrose
  pseudoinverse for near-singular `K`.

### Optimality criteria

#### Fixed treatment effects

**A-criterion** (lower = better): \$\$A\_\text{criterion} =
\bar{v}\_\text{diff} = \frac{2}{p(p-1)} \sum\_{i\<j}
\text{Var}(\hat{\tau}\_i - \hat{\tau}\_j)\$\$

**D-criterion** (lower = better): \$\$D\_\text{criterion} =
\exp\\\left(\frac{\log\det(HVH)}{p-1}\right)\$\$ where \\H = I_p -
p^{-1}J_p\\ is the centering matrix and \\V\\ is the treatment
variance-covariance submatrix.

**Efficiency forms** (higher = better): \$\$A\_\text{efficiency} = 1 /
A\_\text{criterion}, \quad D\_\text{efficiency} = 1 /
D\_\text{criterion}\$\$

#### Random treatment effects (genomic prediction)

**PEV-criterion** (lower = better): mean prediction error variance
across all entry lines. \$\$\text{PEV}\_\text{criterion} =
\frac{1}{p}\sum\_{i=1}^{p} \text{Var}(\hat{u}\_i - u_i)\$\$

**CDmean** (higher = better): mean coefficient of determination for
genomic breeding value prediction (Rincent et al. 2012, *Genetics*
192:715-728). \$\$\text{CDmean} = 1 -
\frac{\text{PEV}\_\text{criterion}}{\sigma_g^2}\$\$

CDmean measures the proportion of genetic variance explained by GEBV
predictions on average across lines. Values close to 1 indicate highly
reliable predictions; values close to 0 indicate near-uninformative
predictions. CDmean is only defined when `treatment_effect = "random"`.

**CD per line**: when the number of treatments does not exceed
`eff_full_max`, a per-line coefficient of determination vector is also
returned: \$\$CD_i = 1 - \text{PEV}\_i / \sigma_g^2\$\$

### Large-design approximation

When the number of treatments exceeds `eff_full_max`, a Hutchinson
stochastic trace estimator with `eff_trace_samples` Rademacher vectors
replaces exact submatrix extraction. The `mode` field carries the suffix
`_APPROX` and `D_criterion`, `D_efficiency`, and `CD_per_line` are `NA`.

## References

Rincent, R., Laloe, D., Nicolas, S., Altmann, T., Brunel, D., Revilla,
P., ... & Moreau, L. (2012). Maximizing the reliability of genomic
selection by optimizing the calibration set of reference individuals:
comparison of methods in two diverse groups of maize inbreds (*Zea mays*
L.). *Genetics*, 192(2), 715-728.

## See also

[`met_alpha_rc_stream()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_alpha_rc_stream.md)
to construct the design.
[`met_optimize_alpha_rc()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_optimize_alpha_rc.md)
to search for a criterion-optimal design.
[`met_evaluate_famoptg_efficiency()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_evaluate_famoptg_efficiency.md)
for the equivalent function for
[`met_prep_famoptg()`](https://FAkohoue.github.io/OptiSparseMET/reference/met_prep_famoptg.md)
designs (uses `sigma_b2` instead of `sigma_rep2`

- `sigma_ib2`).

## Examples

``` r
design <- met_alpha_rc_stream(
  check_treatments = c("CHK1", "CHK2", "CHK3"),
  check_families   = c("CHECK", "CHECK", "CHECK"),
  entry_treatments = paste0("G", 1:167),
  entry_families   = rep(paste0("F", 1:7), length.out = 167),
  n_reps = 3, n_rows = 30, n_cols = 20,
  min_block_size = 19, max_block_size = 20
)
#> Fixed field = 30 x 20; used plots = 591; trailing NA plots = 9; replicate used sizes = {197, 197, 197}; blocks/rep = 10 (derived); block sizes in rep 1 = {20, 20, 20, 20, 20, 20, 20, 19, 19, 19}; min block size = 19; max block size = 20

## Fixed-effect evaluation (BLUE contrasts, AR1xAR1 spatial)
eff_fixed <- met_evaluate_alpha_efficiency(
  field_book         = design$field_book,
  n_rows = 30, n_cols = 20,
  check_treatments   = c("CHK1", "CHK2", "CHK3"),
  treatment_effect   = "fixed",
  residual_structure = "AR1xAR1",
  rho_row = 0.10, rho_col = 0.10
)
eff_fixed$A_criterion   # lower is better
#> [1] 0.7436041
eff_fixed$D_criterion   # lower is better
#> [1] 0.3520789
eff_fixed$A_efficiency  # higher is better
#> [1] 1.344802

## Random-effect evaluation (GBLUP with CDmean)
if (FALSE) { # \dontrun{
eff_gblup <- met_evaluate_alpha_efficiency(
  field_book         = design$field_book,
  n_rows = 30, n_cols = 20,
  check_treatments   = c("CHK1", "CHK2", "CHK3"),
  treatment_effect   = "random",
  prediction_type    = "GBLUP",
  K                  = my_kinship_matrix,
  varcomp            = list(sigma_g2 = 0.4, sigma_e2 = 0.6,
                            sigma_rep2 = 0.1, sigma_ib2 = 0.05,
                            sigma_r2 = 0.02, sigma_c2 = 0.02),
  residual_structure = "AR1xAR1",
  rho_row = 0.10, rho_col = 0.10
)
eff_gblup$CDmean       # mean prediction reliability [0, 1]
eff_gblup$CD_per_line  # per-line reliability vector
eff_gblup$mean_PEV     # mean prediction error variance
} # }
```
