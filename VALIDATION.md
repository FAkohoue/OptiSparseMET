# OptiSparseMET validation summary (v0.2.0)

This document records how the scientifically critical parts of `OptiSparseMET`
are verified. Every numerically important quantity is checked against an
**independent oracle** — a closed-form value, a hand-built matrix, or a separate
implementation — rather than against the package's own formula. The oracles were
derived and checked in an independent numerical environment (NumPy) and are
encoded as `testthat` tests.

Run `devtools::document()` then `devtools::test()` to reproduce.

## Correctness fixes and their ground-truth tests

| Area | What is verified | Oracle | Test file |
|------|------------------|--------|-----------|
| AR1 / AR1×AR1 residual precision | The residual precision on a non-square grid matches the separable AR1 model; row/column axes are not transposed | Hand-built `solve(kronecker(Sigma_col, Sigma_row))`; confirmed `ar1_precision == inv(ar1_cov)` | `test-ar1-precision.R` |
| A-criterion (fixed) | Mean pairwise contrast variance | Balanced-CRD closed form `2·sigma_e2/n` | `test-efficiency-groundtruth.R` |
| D-criterion (fixed) | Geometric mean of contrast eigenvalues | Balanced-CRD closed form `sigma_e2/n` | `test-efficiency-groundtruth.R` |
| PEV / CDmean (random) | Mean PEV and reliability | Independent dense MME solve (`mean_PEV = 0.55`, `CDmean = 0.45`) | `test-efficiency-groundtruth.R` |
| CDmean scaling | `CD_i = 1 − PEV_i/(sigma_g2·K_ii)` uses the relationship diagonal | Identity that holds only when K_ii is in the denominator (tested with `diag(K)=1.5`) | `test-cdmean-kdiag.R` |
| D pseudo-determinant | Relative tolerance, exact rank, scale invariance | Known eigenvalue spectra | `test-logdet.R` |
| Rank-deficiency guard | Singular fixed design is repaired, not fatal | Contrast-invariance | `test-rank-guard.R` |

## Allocation

| Property | Verification | Test file |
|----------|-------------|-----------|
| `equireplicate` rename; `balanced_incomplete` deprecated | Message emitted; identical structure | `test-allocation-balance.R` |
| Balancing preserves equal replication and equal environment size | Row/column sums unchanged after swaps (proven in the swap prototype) | `test-allocation-balance.R` |
| `env_pair` balancing reduces env-pair variance (to 0 where achievable) | Prototype reached variance 0; test checks non-increase | `test-allocation-balance.R` |
| Feasibility reports strict-BIBD conditions | `strict_bibd_possible` FALSE for J≫I (Fisher fails); TRUE for the v=b=7 Fano design | `test-allocation-balance.R` |

## Coupled evaluation, simulation, and decision points

| Function | Verification | Test file |
|----------|-------------|-----------|
| `met_information()` | Fixed-effect absorption matches a full MME solve; balanced < sparse PEV; lower efficiency raises PEV; oracle values `mean_PEV = 0.6225`, `CDmean = 0.17` | `test-met-information.R` |
| `simulate_met()` | Balanced design beats sparse on realized accuracy (0.85 vs 0.73 in NumPy); accuracy in [-1,1]; positive gain | `test-simulate-met.R` |
| `select_environments()` | Spread method avoids near-duplicate environments; all methods return distinct sets | `test-select-environments.R` |
| `select_individuals()` | CDmean in [0,1], rises with training size, exchange beats random | `test-select-individuals.R` |
| `optimize_allocation_gxe()` | Mean PEV never worsens; margins preserved | `test-gxe-allocation-replication.R` |
| `recommend_replication()` | Accuracy rises with replication; seed budget flags infeasible levels | `test-gxe-allocation-replication.R` |
| `sensitivity_varcomp()` | PEV monotone increasing, CDmean monotone decreasing, in the variance ratio | `test-sensitivity.R` |

## External benchmarks

`scripts/benchmark_designs.R` documents optional cross-checks against `sommer`
(GBLUP PEV) and `odw` / `DiGGer` (within-field A-optimality), to run in an
environment where those packages are installed. These are not part of
`R CMD check`.

## Notes

* The `man/*.Rd` files and `NAMESPACE` are regenerated from roxygen with
  `devtools::document()`; the exported set is defined by the `@export` tags.
* A strict balanced incomplete block design (constant pairwise concurrence λ)
  generally cannot exist for sparse testing (J ≫ I); the package targets the
  achievable near-balanced (regular-graph) design and reports the strict-BIBD
  existence conditions rather than claiming a BIBD.
