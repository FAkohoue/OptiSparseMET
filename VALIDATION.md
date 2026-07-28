# OptiSparseMET validation protocol (version 0.2.0)

## Purpose

This document defines the evidence required before an `OptiSparseMET` design is
released for planting. The package addresses a coupled problem: allocation
across environments, replication within environments, field layout, genetic
connectedness, environmental uncertainty, and finite seed stocks all affect the
final information matrix. Validation must therefore test numerical correctness,
structural feasibility, and operational consistency.

Correctness tests use an independent oracle wherever possible. An oracle is a
closed-form result, a hand-constructed matrix, or an implementation that does
not call the package formula under test. Self-consistency checks remain useful
as smoke tests, but they are not treated as proof of correctness.

No command result is inferred from the source tree. The package maintainer runs
the verification sequence manually and retains the complete output.

## 1. Numerical correctness

| Area | Required evidence | Independent oracle |
|---|---|---|
| First-order autoregressive residuals | The precision matrix for a non-square field matches the specified row and column covariance | Direct inversion of the separable row-by-column covariance |
| Fixed-effect A-criterion | Mean pairwise contrast variance is correct | Balanced completely randomised design: \(2\sigma_e^2/n\) |
| Fixed-effect D-criterion | Contrast pseudo-determinant is correct and scale-stable | Balanced completely randomised design: \(\sigma_e^2/n\) |
| Prediction-error variance and CDmean | Mixed-model equations and reliability scaling are correct | Independent dense mixed-model-equation solution |
| Relationship diagonal | \(CD_i=1-\operatorname{PEV}_i/(\sigma_g^2K_{ii})\) uses each treatment's prior variance | A relationship matrix with a non-unit diagonal |
| Rank deficiency | Equivalent fixed-effect parameterisations yield the same contrasts | Direct comparison after removing the redundant intercept |
| Positive-semidefinite matrices | Invalid dimensions, names, asymmetry, and indefiniteness are detected | Hand-constructed valid and invalid matrices |

CDmean denotes the mean coefficient of determination, and PEV denotes
prediction-error variance.

## 2. Across-environment allocation

The allocation tests must establish the following properties.

1. `random_balanced` covers every treatment and respects the requested site
   capacities.
2. `equireplicate` gives every sparse treatment the requested environment
   replication and reproduces each named site load exactly when the degree
   sequence is realisable.
3. Line-pair and environment-pair balancing preserve row and column margins
   while reducing their stated concurrence criterion.
4. Strict balanced incomplete block design (BIBD) conditions are reported
   separately from the achievable equireplicate design.
5. A common treatment is present at every site, but its number of plots may
   differ among sites.
6. The common-set search compares nested candidate counts and audits genetic,
   family, environmental, seed, and capacity consequences.
7. Genotype-by-environment refinement preserves treatment replication and site
   size while seeking a lower across-environment mean PEV.

## 3. Seed feasibility

Seed is one network-wide inventory for each treatment. Validation must
reconstruct consumption from the delivered fieldbooks and verify that:

1. a treatment-site cell consumes seed only when that treatment is present at
   that site;
2. an absent treatment-site cell consumes exactly 0 g;
3. a common treatment is charged at every site because it is present at every
   site;
4. all local replicated plots are included in the charge;
5. environment-specific grams per plot are respected;
6. cumulative consumption never exceeds inventory less the declared reserve;
   and
7. the default reserve is 0 g, while any positive reserve is an explicit
   programme policy.

The delivered seed ledger, rather than an intermediate allocation object, is
the release record.

## 4. Environmental evidence

Environmental validation has four distinct targets.

1. **Data integrity.** Environment-date keys are unique; physical ranges,
   coverage, provenance, imputation, and source precedence are auditable.
2. **Feature construction.** Weather features respect crop stages and time;
   soil features respect depth; management variables preserve treatment
   meaning; and missing values are never silently replaced by zero.
3. **Kernel construction.** Weather, soil, management, and geographic kernels
   are aligned and positive semidefinite. Weather-by-soil,
   weather-by-management, and soil-by-management kernels represent distinct
   interaction hypotheses rather than automatic consequences of main-effect
   similarity.
4. **Evidence branch.** Historical multi-environment trial (MET) responses may
   calibrate trait-specific covariance and defensible kernel weights. Without
   such responses, the central genetic covariance is the identity matrix;
   environmental kernels enter as separate, unweighted maximin scenarios.

Dedicated variable-level interactions require four additional checks:

1. Every requested variable survived environmental quality control and matches
   one exact covariate column.
2. Tensor expansion, residual information, effective rank, bandwidths,
   positive semidefiniteness, and strong-heredity parent mappings are audited.
3. Historical support requires positive blocked-validation improvement, a
   multiplicity-adjusted paired sign test, and a positive fitted covariance
   contribution.
4. Screened or rejected interactions remain separate maximin structures; they
   are not removed from robustness analysis.

Without historical MET responses, no dedicated interaction is labelled
supported or unsupported and no environmental weight is defined.

Environmental grouping must infer both the number of groups and the algorithm.
The candidate solution is assessed by separation, minimum group size,
agreement across algorithms, and reproducibility across modality- or
year-specific relationships. A stable solution may be used as a hard
environmental stratum. A provisional solution is reported with caution. An
unstable solution must fall back to one unpartitioned group rather than force a
fragile mega-environment structure.

## 5. Coupled evaluation and simulation

`met_information()` must reproduce an independent mixed-model-equation solution
for the joint genetic-by-environment information matrix. Actual plot
replication, local design efficiency, residual variance, the genetic
relationship matrix, and the environmental covariance must all affect PEV in
the expected direction.

`simulate_met()` must use the same residual-information convention and report
finite Monte Carlo uncertainty intervals. Comparative tests should establish
directional sanity without treating simulation as a mathematical proof. For
example, a well-connected design should not systematically perform worse than
a deliberately disconnected control under its stated data-generating model.

Robust design comparisons must preserve common random numbers, nested capacity
sweeps, fixed partner-site capacities, and a common uncertainty scenario set.
Physical cost is the number of delivered plots, not the number of non-zero
treatment-site cells.

## 6. Benchmarking and release evidence

`scripts/benchmark_designs.R` is a reproducible public demonstration of the
package benchmark. It uses synthetic data and must not be interpreted as crop-
or programme-specific evidence.

A release benchmark must include:

1. complete testing as an information upper bound;
2. naive random sparse allocation as a weak control;
3. M3 coverage-first allocation;
4. strict M4 equireplication when its degree sequence is feasible;
5. the proposed optimised design;
6. paired Monte Carlo comparison, so every design faces the same realised
   genetic effects and residual standard-normal draws;
7. central and alternative environmental covariance structures;
8. plausible genetic and residual variance combinations;
9. complete loss of each operationally vulnerable site;
10. treatment coverage, physical plots, cost, and network-wide seed
    feasibility; and
11. average and worst-scenario PEV, CDmean, accuracy, gain, selection
    coincidence, and rank quality.

`benchmark_environment_models()` compares independence, each separate
environmental kernel, an explicitly labelled equal-weight single-kernel
comparator, and the historically calibrated covariance. Equal weighting is a
benchmark only. Without historical responses or known truth, environmental
models are not ranked because no defensible response target exists.

`simulate_environment_model_benchmark()` tests recovery of a known covariance
and known component weights. `benchmark_environment_missingness()` masks
covariate cells, rebuilds the quality-control and imputation workflow, and
quantifies degradation. `summarize_design_stability()` reports common-set,
site-capacity, and allocation reproducibility.

Independent comparison with other mixed-model or field-design software remains
valuable when the same covariance, fixed effects, residual structure, and
variance components can be imposed. External agreement supports implementation
validity, but it does not replace the analytic oracles above.

## 7. Manual verification sequence

Run these commands from a clean R session at the package root:

```r
devtools::document()
testthat::test_local()
devtools::check(document = FALSE, vignettes = FALSE)
devtools::install()
```

For an end-to-end validation dataset, inspect the site sizes, allocation
margins, environment-group status, seed ledger, common-treatment replication
matrix, minimum reserve, scenario-specific benchmark, design stability,
Pareto frontier, and simulation intervals. The Pareto graph must use one
vertical endpoint, 0--100% cost minimisation below, and decreasing real cost
above. Only designs passing coverage, physical-capacity, and seed checks may
define that frontier. Repeat hybrid validation with dominance disabled to
verify the additive-only path. Script completion alone is insufficient.

The design is ready for a plant breeder only when all numerical checks pass and
the delivered fieldbooks agree exactly with the capacity and seed ledgers.
