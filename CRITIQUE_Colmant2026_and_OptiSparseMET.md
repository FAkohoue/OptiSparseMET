# Critical appraisal of Colmant et al. (2026) and the OptiSparseMET response

## Purpose and scope

Colmant et al. (2026) provide a useful decision framework for sparse
multi-environment trials (METs). Their central contribution is to treat
environment choice, genotype choice, allocation, replication, and field layout
as related design decisions. This appraisal identifies the parts of that
framework that `OptiSparseMET` adopts, the assumptions it does not adopt, and
the methodological safeguards implemented in response.

The target population of environments (TPE) is the set of environments for
which the breeding programme seeks inference. The target population of
genotypes (TPG) is the corresponding set of candidate genotypes. Prediction
error variance (PEV) and the mean coefficient of determination (CDmean) measure
the precision of genetic predictions.

This is a methodological appraisal, not a claim that one software package is
universally superior. Empirical superiority requires independent comparison
with historical trials, realistic variance estimates, real marker or pedigree
relationships, and current breeding designs.

## 1. Contributions retained from the paper

Four ideas are retained.

1. Sparse MET design should be expressed as a sequence of explicit decisions.
2. The number of environments, number of genotypes, replication, and field
   size must be compared under a common resource budget.
3. Genetic and environmental relationships can inform design.
4. Random allocation is a necessary benchmark rather than an adequate default.

The paper's \(-q'Dq\) optimal-contribution criterion is also retained where its
interpretation is appropriate. With \(D\) defined as similarity and \(q\)
indicating a selected subset, minimising \(q'Dq\) favours mutual dissimilarity.
It is therefore a spread criterion.

## 2. Limits of the spread objective

Spread and representativeness are not synonymous. A set of mutually dissimilar
environments can consist of atypical extremes and fail to represent the centre
of the TPE. Similarly, a genetically diverse subset can be valuable for
diversity conservation without being the most informative training set for the
remaining TPG.

`OptiSparseMET` therefore uses:

- facility-location coverage as the default environment-representativeness
  objective;
- CDmean for genotype training-set representativeness;
- \(-q'Dq\) as an explicit spread baseline; and
- MET-level PEV for allocation refinement when the necessary covariance model
  is defensible.

This separation avoids interpreting a diversity proxy as direct prediction
accuracy.

## 3. Genetic-by-environment targets

A main-effect-only breeding-value target can conceal the value of allocation
when genotype-by-environment interaction is material. Conversely, a complex
interaction model is not automatically justified by small or weak historical
data.

The package distinguishes three prediction targets:

1. an across-TPE mean;
2. environment-specific genetic values; and
3. group-specific genetic values when environmental grouping is supported.

`met_information()` combines the treatment-by-environment incidence matrix,
actual replication, local design efficiency, residual variance, genetic
relationship, and environmental covariance. `simulate_met()` evaluates the
corresponding targets under declared data-generating models and reports Monte
Carlo uncertainty.

## 4. Environmental covariance and interactions

Environmental similarity is not automatically genetic covariance. Weather,
soil, management, and geography describe exposure; historical MET responses
describe trait-specific genetic response.

`build_environment_kernels()` therefore preserves separate main-effect kernels.
Cross-modality kernels represent weather-by-soil,
weather-by-management, soil-by-management, or declared higher-order
interactions. These interaction kernels do not assume that environments similar
under main effects must be similar under interactions, or that environments
dissimilar under main effects must be dissimilar under interactions.

Two evidence branches follow.

### Historical MET responses are available

Historical responses may calibrate trait-specific environmental covariance and
support empirical kernel contributions. The calibration must be cross-validated
and regularised because environment networks are usually small.

### Historical MET responses are unavailable

The identity matrix is the central genetic covariance. Environmental main and
interaction kernels remain separate, unweighted maximin scenarios. The user is
not asked to choose arbitrary weights, and an equal-weight average is not
silently substituted.

This no-history branch does not claim that environments are biologically
independent. It states that trait-specific genetic covariance has not been
identified and evaluates whether the design remains acceptable across
plausible environmental structures.

## 5. Environmental grouping

A fixed number of mega-environments is rarely defensible before the data are
examined. `infer_mega_environments()` estimates the number of groups across an
admissible range and compares hierarchical clustering with finite-descent
k-means. Selection uses separation, minimum group size, algorithm agreement,
and reproducibility across modality- or year-specific relationship matrices.

The output is stable, provisional, or unstable. An unstable result returns one
unpartitioned group. A provisional environmental cluster is not labelled a
genetic mega-environment without historical response evidence.

This fallback is preferable to forcing a biologically attractive but
statistically unstable classification.

## 6. Allocation and connectedness

In breeding networks, the number of treatments commonly exceeds the number of
environments. A strict balanced incomplete block design (BIBD), which requires
constant pairwise concurrence, is then often impossible. Equal treatment
replication and exact site sizes are achievable under a valid degree sequence,
but they do not by themselves prove a BIBD.

`OptiSparseMET` accordingly distinguishes:

- coverage-first `random_balanced` allocation;
- `equireplicate` allocation with exact treatment and site margins;
- line-pair and environment-pair concurrence improvement;
- strict BIBD feasibility reporting and optional exact construction when
  mathematically possible;
- common treatments for direct cross-environment links; and
- PEV-based genetic-by-environment refinement.

Common treatments are not assigned a fixed count or a common replication
scalar. Their number, identity, and site-specific replication can be optimised
jointly. Presence at every site is mandatory; equal replication is not.

## 7. Seed and physical feasibility

An allocation is not deployable if the same seed inventory is implicitly spent
again at every environment. The package treats each treatment's seed as one
network-wide stock.

A treatment consumes seed only where it is present. An absent
treatment-by-site cell consumes 0 g. A common treatment consumes seed at every
site. Local replication creates additional physical plots and therefore
additional seed demand. The default reserve is 0 g; a positive reserve must be
declared as a programme policy.

The final seed ledger and plot counts are reconstructed from delivered
fieldbooks. This guards against discrepancies between an incidence matrix and
the design actually planted.

## 8. Within-environment design

Across-environment allocation and local field layout cannot be optimised in
isolation. Local blocking and spatial correlation determine each site's
information contribution. `OptiSparseMET` supports augmented, partially
replicated, randomised complete block, and alpha row-column designs under
independent and autoregressive residual structures.

The package evaluates fixed-effect contrasts and random-effect prediction using
A-, D-, PEV-, and CDmean-based criteria. Numerical guards address residual-axis
orientation, relationship-matrix scaling, pseudo-determinant rank, and
rank-deficient fixed-effect parameterisations.

## 9. Remaining evidence requirements

The implemented framework still requires empirical scrutiny. In particular:

1. historical rice METs should provide realistic genetic and residual variance
   estimates;
2. parental marker data should validate additive and optional dominance
   testcross relationships;
3. environmental feature windows should be checked against crop phenology;
4. no-history maximin decisions should be compared with subsequent realised
   response;
5. competing breeding designs should be evaluated at equal physical cost; and
6. all numerical claims should be reproduced by the manual protocol in
   `VALIDATION.md`.

The appropriate conclusion is therefore conditional. `OptiSparseMET` provides
a stronger and more auditable design framework than isolated heuristic
allocation, but its value for plant breeders must be demonstrated with the
programme's own environments, germplasm, seed stocks, and historical outcomes.
