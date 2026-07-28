# OptiSparseMET 0.2.0 (development)

## Statistical, environmental, and operational hardening

* Added `build_variable_interaction_kernels()` for pre-specified within- and
  cross-modality variable hypotheses. It matches exact post-quality-control
  columns, constructs row-wise tensor features, residualises them against their
  parent main effects, builds Gaussian bandwidth-ensemble kernels, enforces
  strong-heredity parent retention, and returns variable, hierarchy, dimension,
  residual-information, bandwidth, and rank audits.
* `build_environment_kernels()` now accepts `variable_interactions` and the
  compact `variable_interaction_control` list. Omnibus modality interactions
  and dedicated variable interactions remain distinct structures.
* Added `assess_variable_interactions()` for component-wise leave-year-out or
  leave-environment-out validation, Holm multiplicity control, and bootstrap
  inclusion stability. Without historical MET evidence, no support decision or
  weight is created.
* `calibrate_environment_covariance()` now accepts
  `eligible_interactions`. Component-wise screened kernels can therefore be
  excluded from the central covariance while remaining available as separate
  structural uncertainty candidates.
* Added a dedicated environmental-interaction vignette and expanded the
  README, introductory vignette, benchmarking vignette, breeder guide, and
  regression documentation.
* Added a public benchmarking and validation framework. It constructs complete,
  naive-random, M3, and feasible strict-M4 comparators; evaluates designs with
  paired genetic and residual Monte Carlo draws; crosses covariance,
  variance-component, and complete site-loss scenarios; and reports accuracy,
  gain, prediction-error variance (PEV), mean coefficient of determination
  (CDmean), plots, cost, seed, and coverage in separate ledgers.
* Added known-truth environmental covariance recovery, response-based comparison
  of separate kernels with an explicitly labelled equal-weight single-kernel
  comparator, repeated environmental-covariate masking and imputation stress
  tests, and common-set, site-capacity, and allocation stability summaries.
  Expected low-coverage warnings created by deliberate masking are captured in
  an auditable per-run ledger instead of being emitted repeatedly.
  Equal weighting remains a benchmark only and is never adopted as the
  no-history covariance.
* Pairwise common-set connectivity now exposes `"maximin"`, `"cvar"`, and
  `"mean"` aggregation. Maximin is the default, so the global common-set size
  and correlation-adaptive allocation of additional distinct shared treatments
  protect the weakest environment pair unless the breeder explicitly chooses
  another rule.
* Corrected the definition of cross-environment connectivity: it is now based
  on distinct shared genotypes, never repeated plots of the same genotype.
  Optional `Sigma_E`-guided allocation refinement preserves treatment and site
  margins while allowing weakly correlated pairs to share more distinct
  treatments than strongly correlated pairs. Replication remains a separate
  seed-feasible precision decision.
* Simplified `optimize_common_treatments()` by removing `candidate_counts`,
  `min_reps_per_environment`, `max_reps_per_environment`,
  `common_extra_plot_capacities`, and `replication_penalty`. Supply `n_common`
  to fix the global core size, or leave it `NULL` for the automatic feasible
  search. The former `common_reps` output is replaced by the unambiguous binary
  `common_presence` matrix.
* Reworked `plot_pareto_frontier()` to display one statistical endpoint at a
  time. Cost minimisation runs from 0% to 100% on the lower axis, corresponding
  actual cost decreases on the upper axis, and only observed feasible
  non-dominated designs are connected. Infeasible comparators remain visible
  but cannot define the frontier.
* Added `infer_mega_environments()` to infer, rather than require, the number of
  groups in a multi-environment trial (MET). The procedure compares hierarchical and
  kernel-principal-coordinate k-means partitions across a feasible range,
  screens candidates by silhouette and minimum group size, and requires
  adjusted Rand index stability across algorithms, modalities, years, and
  block-balanced relationship bootstraps. Unsupported partitions return one
  explicitly unpartitioned network; single-block results are labelled
  provisional rather than interpreted as genetic mega-environments.
* Mega-environment k-means now uses a finite-descent Hartigan relocation
  algorithm with guaranteed non-empty initial groups. Every accepted move
  strictly reduces within-cluster sum of squares, so every start terminates
  without an iteration-limit failure. Multiple starts, deterministic seeds,
  supported alternatives, and the unpartitioned fallback are retained.
* Environmental interaction kernels now include weather-by-soil,
  weather-by-management, and soil-by-management. Their default construction is
  the centred functional analysis of variance (functional-ANOVA) product, which
  reduces constant/main-effect
  leakage while preserving positive semidefiniteness. Explicit higher-order
  terms are supported by `add_environment_kernel_interactions()`; the legacy
  raw Hadamard product remains available with `mode = "product"`.
* Interaction construction no longer implies interaction activation.
  `calibrate_environment_covariance()` now compares main-only and
  main-plus-interaction models by paired blocked validation and fixes
  unsupported interaction weights at zero in the central covariance. Every
  interaction kernel remains an unweighted structural uncertainty scenario
  when evidence is absent or insufficient: lack of response evidence is not
  treated as absence of physical interaction. The design workflow follows the
  central-model evidence gate and applies maximin design across main and
  interaction structures.
  Calibration now defaults to `ridge = 0`, so prior-weight shrinkage is opt-in.
* Corrected MET simulation residual sampling so that
  `reps * efficiency / sigma_e2` is applied exactly once. Simulation now
  reports Monte Carlo standard errors and 95% intervals for accuracy and gain.
* Physical plot cost now uses supplied replication counts, rather than merely
  counting non-zero treatment-by-environment cells.
* Added network-wide seed constraints to allocation, local partially
  replicated (p-rep) planning, and criterion-driven optimisation. A genotype
  consumes seed only at sites where it is present, and consumption accumulates
  across all delivered plots. Final MET plans reconstruct the ledger from the
  field books; the optional per-treatment reserve defaults to 0 g.
* Seed-constrained construction assigns treatments with the fewest affordable
  environments first, then fills the most constrained environment first.
  Post-allocation balance remains active with genetic groups and seed budgets
  by accepting only within-group, seed-feasible swaps.
* Strengthened matrix, naming, dimension, symmetry, positive-semidefinite,
  replication, dosage, and covariance validation throughout the MET and
  relationship-matrix core.
* `check_equireplicate_feasibility()` now reports slot equality and binary
  degree-sequence realizability separately, rather than treating matching slot
  totals alone as sufficient.
* Additive relationships remain the default. Dominance is explicit and
  optional for hybrid or specific combining ability targets; hybrid additive
  relationships can be derived from the parental genomic relationship matrix
  (GRM) with stable automatic hybrid identifiers.
* `suggest_site_capacity()` now supports a focal subset of unconstrained sites
  while partner capacities remain fixed. Candidate designs are nested, use a
  common Monte Carlo stream, and count environment-specific replication and
  check overhead. It can now evaluate every capacity under a
  `robust_scenarios()` set using maximin, mean, or lower-tail conditional
  value-at-risk (CVaR) aggregation.
* `optimize_design()` can enforce network seed budgets and per-environment
  lower/upper capacity bounds.
* `pracma` remains an imported dependency and the package continues to use
  `pracma::mod()` for the field-layout arithmetic that requires it.

## Environmental covariance, weather, soil, and provenance

* Added `build_environment_kernels()` and
  `combine_environment_kernels()`. Weather, soil, management, geography, and
  available cross-modality interaction products are now normalised and retained
  as separate modalities, preventing high-dimensional weather from receiving
  accidental weight by column count. The interaction kernels are structural
  hypotheses whether or not historical responses are available. Gaussian
  blocks use a bandwidth ensemble.
  Multi-kernel combinations now require explicit weights rather than silently
  assigning equal weights.
* Added `calibrate_environment_covariance()`: non-negative simplex weights are
  fitted against historical MET genetic correlations with
  leave-environment-out and optional leave-year-out validation, bootstrap
  weights, and covariance sensitivity candidates. Without historical genetic
  information, no weights are defined: the central covariance is identity and
  the modality kernels remain separate sensitivity scenarios.
* Added `consensus_environment_kernels()`, a median-of-within-kernel-ranks
  integration with no elicited modality weights, for descriptive environmental
  clustering only.
* Added `qc_environmental_data()` with duplicate-key checks, physical-range
  handling, coverage thresholds, explicit median/linear imputation,
  missingness indicators, provenance, and row-level imputation ledgers.
  Supplied and fetched source blocks now merge cell by cell.
* `fetch_weather_series()` now preserves, validates, and orders actual POWER
  response dates; computes day-after-start from dates; and attaches request,
  coverage, package-version, cache, and MD5 provenance.
* Added `calibrate_weather_series()` for station-based additive,
  precipitation-ratio, or linear bias correction. Historical envirotyping now
  retains Q10/Q50/Q90, trends, observed-year counts, station diagnostics, and
  bootstrap confidence intervals for across-year kernel stability.
* Added `rice_weather_features()` for phenology- or thermal-time-aligned
  vapour-pressure deficit (VPD), diurnal range,
  precipitation-minus-evapotranspiration (P-ET) water balance, hot nights,
  heat/cold/root-water stress
  severity, and consecutive dry/wet/hot spells.
* `fetch_soilgrids()` now supports multi-depth and multi-quantile retrieval plus
  neighbourhood mean/SD extraction. New `soil_profile_features()` derives
  thickness-weighted root-zone properties, Q0.95-Q0.05 uncertainty widths,
  spatial heterogeneity, available water capacity, and texture log ratios.
* `consensus_relationship()` now uses centred off-diagonal relationships for
  STATIS weights by default, so a common diagonal cannot inflate apparent
  agreement.
* `robust_scenarios()` can cross variance/shrinkage assumptions with alternative
  `Sigma_E` matrices. These candidates feed robust common-treatment and
  final-design evaluation.

## Common treatments, capacity, robustness

* New `optimize_common_treatments()` replaces mandatory user prespecification
  of common-set size with robust selection of size and identities. It returns a
  binary common-presence matrix and integrates variance-component
  and `Sigma_E` uncertainty with Bayes, maximin, or lower-tail CVaR criteria;
  evaluates genetic effective sample size as a distinct diversity criterion;
  and enforces family/group coverage, entry capacity, full population coverage,
  and one network-wide seed inventory. Pairwise connectivity is the literal
  number of distinct shared treatments. Local replication is assigned
  downstream from remaining seed and never increases connectivity.
* `plan_sparse_met_design()` now accepts per-environment
  `fixed_treatment_reps` inside p-rep design specifications. Optimised common
  replication is enforced first, and unused p-rep capacity remains available
  to non-common treatments.
* New `suggest_common_treatments()` suggests both the number and identities of
  common treatments when the breeder does not fix the number. The suggestion
  uses seed feasibility at every site, family and genetic-group
  representativeness, and a connectivity target for estimating G×E and
  between-site correlations. A supplied number remains authoritative.
* `suggest_site_capacity()` now defaults to `select = "max"` -- for an
  unrestricted site the breeder supplies only the capacity interval and the
  package returns the size that maximises accuracy within it (no rule to choose).
* Robustness pass (numerical-core audit): guarded further factorizations against
  singular inputs -- the Pesek-Baker desired-gain solve (singular trait
  covariance), the `tune_common_treatments()` coefficient solve, and the
  `met_information_mt()` residual-covariance Cholesky (now bends to nearest PD)
  all fall back gracefully, matching the earlier `simulate_met()` hardening. The
  core index patterns (`seq_len` throughout) and the `met_evaluate_*` /
  `met_information()` inverses were reviewed and are already guarded.

## Genomic-matrix tooling (additive + dominance; AGHmatrix / ASRgenomics or native)

* `build_relationship_matrix()` gains `method = "AGHmatrix"`
  (`AGHmatrix::Gmatrix`), a `ploidy` argument, a `tuneup` option, and
  `relationship = "additive" | "dominance"` -- the dominance matrix (Vitezica
  2013; `Gmatrix(method = "Vitezica")` or a built-in estimator) is important for
  hybrids/SCA.
* New `build_hybrid_relationship()`: additive genomic relationship among
  hybrids/testcrosses derived from the parental G and the cross design
  (`0.25(G[a,c]+G[a,d]+G[b,c]+G[b,d])`), so non-genotyped hybrids and the integer
  0/1/2 requirement of AGHmatrix are both handled.
* New `combine_relationship_matrices()`: variance-weighted combination of
  component matrices (e.g. additive + dominance) into one relationship for total
  genetic merit.
* New `tune_relationship_matrix()`: makes a matrix well-conditioned and
  invertible. `method = "asrgenomics"` uses `ASRgenomics::G.tuneup()`; the
  **native** `"bend"` (nearest positive-definite via eigenvalue flooring) and
  `"blend"` need no external package, so ASRgenomics is optional.
* New `kinship_pca()`: population-structure PCs via `ASRgenomics::kinship.pca()`
  (eigen fallback), for structure-aware allocation grouping and diagnostics.
* `simulate_met()` hardened to a PD-safe covariance root and pseudo-inverse, so
  a singular `Sigma_E` (e.g. two perfectly correlated environments) degrades
  gracefully instead of erroring.


## Planning-stage envirotyping (design before the season)

* `fetch_weather_series()` / weather fetch now **retry** (`max_tries`) with a
  longer timeout, so intermittent NASA POWER timeouts no longer drop sites;
  `historical_envirotype()` exposes `max_tries`.
* New `consensus_relationship()`: a **robust** consensus of several relationship
  matrices -- STATIS/RV-weighted (down-weights anomalous years, the default),
  log-Euclidean geometric mean (the principled SPD-manifold mean), element-wise
  median, or the plain arithmetic mean -- rather than only averaging.
* New `assess_envirotype_stability()`: builds the environment relationship matrix
  **per year**, reports how stable it is across years (Mantel-type correlation),
  and returns a robust **consensus** `D` (via `consensus_relationship()`,
  STATIS by default) to design on when stability is low.
* `historical_envirotype()` now also returns `combined` = the typical profile
  **plus** the interannual variability (columns suffixed `_iav`), so both
  within-season (`sd` features) and between-year variation can *enter* the
  environment relationship, not merely be reported.
* `historical_envirotype()` supports **site-specific growing windows** (`window`
  as a per-site data frame `environment`/`start_md`/`end_md`, or columns on
  `sites`), including windows that cross the calendar year (e.g. Peru Nov->Feb);
  and `station_daily` to **merge on-station manually-collected daily weather**
  with the downloaded series (station values take precedence, downloaded fills
  the gaps).
* `fetch_weather_series()` gains `cache_dir` (cache responses for instant
  re-runs) and `workers` (parallel fetch via \pkg{future.apply}); it warns if
  more than 20 daily parameters are requested. `historical_envirotype()` passes
  both through.
* SoilGrids documentation corrected: it is a single static model, so
  time-resolved soil (organic carbon, N, pH, moisture) needs the NASA POWER
  `GWET*` variables, satellite soil moisture, or user-supplied dated
  measurements (`soil =`), which may be one set per year.

* New `historical_envirotype()`: at design time the upcoming season's weather
  does not exist, so it characterises each candidate site from **several past
  seasons** over the intended calendar window -- returning the **typical**
  environmental profile (across-year mean, for descriptive site comparison)
  and the **interannual variability** (across-year SD/CV, the site's
  environmental risk, used to construct separate stress-test scenarios).
  Enviromic covariates alone are no longer described as predicting genetic
  `Sigma_E`; historical MET genetic information is required to calibrate that
  covariance.

## Temporal envirotyping (keep within-season dynamics)

* New `envirotype_features()`: turns a daily weather series into rich,
  fixed-length per-environment features instead of a single season mean --
  by time window / crop stage, with variability (`sd`, `cv`), stress-day counts,
  growing degree days, and an optional B-spline **functional** representation of
  the trajectory (on normalised phenological time). This lets the environmental
  kinship reflect how sites differ in their day-to-day variation and its timing,
  not just their averages. Stages can be defined by expert crop physiology
  (named day ranges), by equal intervals (`windows = K`), or by day-after-start
  cut points. (Dynamic time warping was considered and deliberately not used: it
  warps away the timing signal that envirotyping needs.)
* New `fetch_weather_series()`: retrieve the per-site **daily** NASA POWER series
  (long format) to feed `envirotype_features()`.

## Configurable weather fetch and supplied+fetched merging

* `build_enviromic_covariates()` weather fetch is now fully configurable:
  `weather_pars` requests any NASA POWER parameter (wind, dew point, soil
  moisture GWET*, evapotranspiration, PAR, ...), `weather_stats` sets how each is
  aggregated over the window (mean/sum/min/max/sd/median), and
  `weather_temporal` selects `"daily"`, `"monthly"`, or `"climatology"`.
* Supplied and fetched data are now **merged**: your own `weather`/`soil`
  columns take precedence and only the variables you are missing are fetched, so
  a partly-collected dataset is completed automatically into one matrix.
* `enviromic_variable_catalog()` extended to document the additional NASA POWER
  weather parameters (codes, meanings, units) alongside the SoilGrids properties.

## Management covariates

* `build_enviromic_covariates()` gains `scale_dummies` (default `TRUE`): set
  `FALSE` to leave 0/1 dummy columns (one-hot management, indicators) as-is while
  still standardising the continuous weather/soil/dose covariates.
* New `add_management_interactions()`: build management interaction covariates
  (e.g. fertiliser dose within each fertiliser type, planting mode by irrigation)
  to add to the management data before assembling the environmental covariates.
* Interaction selection: `list_candidate_interactions()` enumerates every
  possible interaction among covariate columns (labelled by type and column
  cost) so terms can be chosen deliberately; `screen_interactions()` ranks
  environmental covariates/interactions by variance explained in observed G x E
  (double-centre + SVD of the genotype-by-environment matrix; `"anova"` R^2 by
  default, optional random-forest importance via \pkg{ranger}); and
  `agronomic_interaction_hypotheses()` returns a curated, hypothesis-driven
  shortlist of known environmental/management interactions. Non-linear covariate
  effects are otherwise captured implicitly by the Gaussian environmental kernel
  in `build_environment_relationship()`.
* New `nest_dose_within_type()`: enforce the correct encoding of management
  **doses**, which must never be used ungrouped -- 100 kg NPK is not 100 kg urea,
  and 10 mm drip is not 10 mm manual spray. It replaces each dose with
  type-nested dose columns (dose carried only within its own fertiliser /
  pesticide / irrigation type, zero otherwise) and drops the raw dose, so a dose
  is always interpreted relative to the product/method that delivered it.

## Connectivity, plot-budget, and trade-off visualisation

* `recommend_common_treatments()`: reads the typical environment correlation
  `rho` directly from the environmental kinship (`Sigma_E` or an environment
  relationship matrix) and returns the suggested number of common treatments --
  optionally the simulated optimum via `tune_common_treatments()` -- so
  connectivity need not be guessed.
* `test_entry_capacity()` / `required_plots()` / `field_plot_accounting()`:
  convert between a partner's offered **total plots** and the number of **test
  entries** it holds, accounting for check plots (`n_checks * n_blocks`) and
  replication. `suggest_site_capacity()` gains `check_plots_per_site` so its
  `total_plots` reflects the real field size, not just test entries.
* `plot_pareto_frontier()`: plots a `pareto_designs()` result -- budget (plots or
  cost) against expected genetic gain and reliability, with the Pareto-optimal
  designs highlighted.

## Criterion-driven optimal-design engine

A deterministic engine (built on the coupled `met_information()` matrix, so no
Monte-Carlo inside the optimisation loop) that optimises the allocation itself,
rather than only refining a greedy construction:

* `selection_intensity()` and `expected_genetic_gain()`: the breeder's-equation
  gain `dG = i * r * sigma`, single-trait or as a multi-trait selection index
  (`trait_weights` + `trait_gencov`). (For one trait, gain is a monotone
  transform of CDmean; the gain scale becomes distinct for multi-trait indices.)
* `met_information_mt()`: the **exact** multi-trait extension of
  `met_information()` -- the full trait-covariance Kronecker MME
  (`Sigma_T (x) Sigma_E (x) G` genetic, `R_T` residual). It returns the exact
  prediction-error variance and reliability of a desired-gain / economic index
  via the canonical transformation (Hayes-Hill 1981), which reduces the
  multi-trait system to independent single-trait analyses and recombines them
  (verified against a direct multi-trait MME solve to machine precision).
  `design_objective()` / `optimize_design()` use it automatically when a trait
  index is supplied (`multitrait = "exact"`, the default), so the whole
  optimisation engine is now exact for multi-trait targets rather than the
  earlier first-order approximation.
* `desired_gain_weights()`: multi-trait index weights from breeder-defined
  *desired gains* instead of subjective economic weights -- the classical
  Pesek-Baker index (`b = G^-1 d`, built-in and deterministic) or extracted from
  the \pkg{DGQGSI} package (desired-gain and quadratic genomic indices). Feed its
  `weights`/`gen_cov` into the optimisation engine as `trait_weights`/
  `trait_gencov` to optimise a design for a desired-gain multi-trait target.
* `design_objective()`: one deterministic score combining statistical quality
  (CDmean / mean PEV), expected genetic gain, and resource cost with tunable
  weights and an optional plot `budget`.
* `optimize_design()`: simulated-annealing exchange with multiple restarts that
  maximises `design_objective()`. Move sets (`preserve = "margins"` /
  `"replication"` / `"none"`) keep a fixed budget or M4 equal replication, or let
  plots vary so cost is optimised jointly with quality and gain.
* `robust_scenarios()` + `robust_design_score()`: robust/Bayesian scoring over a
  prior on the variance components (heritability regime and trust in `Sigma_E`),
  aggregated as the Bayes expectation (`"mean"`), worst case (`"min"`), or CVaR
  (`"cvar"`). `optimize_design()` accepts `robust=` to optimise the robust score
  directly -- criterion, gain, cost and robustness simultaneously.
* `pareto_designs()`: the benefit-vs-cost trade-off frontier, flagging the
  non-dominated designs so the breeder can pick the point matching their budget.


## From raw inputs to a design (markers, pedigree, GPS, management)

* New `build_relationship_matrix()`: VanRaden genomic `G` from markers, pedigree
  `A` (Henderson recursion, with inbreeding), and single-step `H` blending both
  (Legarra et al. 2009) so markers and pedigree are used jointly and
  non-genotyped lines are included.
* `build_environment_relationship()` gains a `geo` option using great-circle
  (haversine) distance, correct for site GPS coordinates at continental scale.
* New `build_enviromic_covariates()`: assemble weather, soil, and management
  (planting system, fertilisation, ...) into the environment-by-covariate matrix,
  with optional fetching of weather (NASA POWER) and soil (ISRIC SoilGrids) from
  GPS coordinates. Fetching is gated behind optional packages and never runs in
  tests or examples.
* New `cluster_environments()`: partition sites into mega-environments (with the
  within-cluster covariance) to drive `bv_target = "mega_environment"` evaluation
  when environments are weakly to moderately correlated.
* `build_environment_relationship()` gains `variables` and `weights` for the
  enviromic source, so a breeder can restrict (and weight) the environmental
  kinship to the specific weather/soil covariates they judge to be the *limiting
  factors*, rather than using every assembled covariate.
* `allocate_sparse_met()` now accepts a **named** per-environment capacity vector
  (matched to environment names, any order) in `n_test_entries_per_environment`,
  so sites need not be equal-sized: a breeder can occupy exactly the space each
  partner offers per location. Both allocation methods now support these
  capacities when their respective feasibility conditions hold.
* New `enviromic_variable_catalog()`: lists every fetchable weather/soil variable
  with its exact column name, description and units, so a breeder knows which
  names to pass to the `variables` argument (SoilGrids codes like `phh2o`/`bdod`
  are otherwise opaque). The `variables` "not found" error now also lists the
  available column names.
* `allocate_sparse_met()` `"equireplicate"` (M4) now supports **unequal
  environment sizes** while preserving equal replication. Equal replication --
  not equal environment size -- is M4's defining strength, so a breeder can keep
  the stronger M4 method and still occupy site-specific capacities, provided the
  equal-replication degree sequence is realisable (a Gale-Ryser condition checked
  internally; the strict constructor fills sites by largest remaining capacity to
  realize it).
* New `suggest_site_capacity()`: for a location offered with *no* field-capacity
  limit, sweeps candidate plot counts, simulates realised accuracy and genetic
  gain, and recommends the plot number at the point of diminishing returns --
  maximising gain while optimising resource use. `scope = "all"` optimises every
  site together when all sites are unconstrained, returning the common per-site
  capacity to use everywhere.
* New `fetch_soilgrids()`: production-grade SoilGrids retrieval with selectable
  backend -- `"wcs"` (OGC Web Coverage Service, default), `"webdav"` (direct VRT
  access over GDAL `/vsicurl`), `"rest"` (legacy beta endpoint), or `"local"`
  (local raster files). Reprojects site coordinates onto the SoilGrids
  Interrupted Goode Homolosine grid (EPSG:152160) via \pkg{terra} and reads only
  the needed cells. Applies the SoilGrids per-property conversion factors (e.g.
  bdod/nitrogen /100, clay/sand/silt/soc/etc. /10) so values are returned in
  conventional units. Supports the six standard depth intervals and the mean /
  Q0.05 / Q0.5 / Q0.95 statistics. `build_enviromic_covariates()` gains
  `soil_backend`, `soil_properties`, `soil_depth`, `soil_quantile`, and
  `soil_local_paths` to route through it; the default is now the robust WCS
  backend rather than the fragile REST endpoint.

## Design-strategy testing (after Mothukuri et al. 2025)

* `simulate_met()` now also returns breeder-facing selection outcomes:
  `common_selected_mean` (fraction of the truly-best set that prediction also
  selects) and `avg_rank_mean` (mean true rank of the predicted-selected lines),
  alongside accuracy and genetic gain.
* `select_individuals()` gains a `criterion` argument exposing the eight
  training-set relationship measurements (CDMEAN, CDMAX, PEVMEAN, PEVMAX, AOPT,
  DOPT, GOPTPEV, negative distance) and returns all eight for the chosen set. The
  previous `method = "cdmean_exchange"` remains as an alias.
* New `validate_design_criteria()`: scores a set of candidate allocations on both
  the design-stage criteria (`met_information()`) and realised outcomes
  (`simulate_met()`), and reports their correlation across designs -- so a
  programme can check whether optimising a cheap criterion actually improves
  selection, rather than assuming it does.
* New `fa_sigma_e()`: build or low-rank-approximate a factor-analytic environment
  covariance for `met_information()` / `simulate_met()`.

## Correctness fixes (P0)

* **Random-treatment PEV/CDmean fixed-effect offset (major).** In
  `met_evaluate_alpha_efficiency()` and `met_evaluate_famoptg_efficiency()`, the
  genotype index used to extract prediction-error variances from the coefficient
  matrix `C = [[X'QX, X'QZ], [Z'QX, Z'QZ + G^{-1}]]` was indexed within the
  random block `Z`, but the block is offset by the number of fixed-effect
  columns. As a result the reported `mean_PEV` / `CDmean` / `CD_per_line`
  averaged over a window shifted by `ncol(X)` -- including a nuisance random
  effect's PEV and dropping the last genotype(s). This affected every
  random-treatment (genomic-selection) evaluation and the `CDmean`-criterion
  optimisers. Fixed by offsetting the extraction index by `ncol(X)`; verified
  against an independent dense MME solve (single-cell example: corrected
  `mean_PEV = 0.55`, previously `0.6625`).

* **AR1 / AR1×AR1 residual precision**: fixed a transposed Kronecker index that
  swapped the row/column autocorrelation axes and scrambled spatial adjacency on
  non-square fields. Residual-precision construction is now centralised in the
  tested helper `.build_residual_precision()` and verified against a hand-built
  AR1×AR1 precision (`test-ar1-precision.R`). Affects every AR1-based A/D/CDmean
  value in `met_evaluate_alpha_efficiency()` and `met_evaluate_famoptg_efficiency()`.
* **CDmean scaling**: per-line reliability now divides by `sigma_g2 * K_ii`
  (diagonal of the relationship matrix) instead of `sigma_g2` alone, correcting
  the value for VanRaden GRMs and inbred/pedigree material where `diag(K) != 1`.
* **D-criterion**: `.safe_logdet_psd_dense()` now uses a relative eigenvalue
  tolerance and an exact expected rank, making the pseudo-determinant
  scale-robust.
* **Rank-deficiency guard**: a singular fixed-effects design (e.g. intercept plus
  a full dummy set with no check rows) is now repaired instead of erroring.

## Spatial models

* Residual structures extended beyond `IID` / `AR1` / `AR1xAR1`:
  * **`"AR1xAR1_nugget"`** adds an independent measurement-error (nugget)
    component. Pure AR1xAR1 often fits field data poorly, so this is usually the
    more realistic applied model.
  * **`"exponential"`, `"gaussian"`, `"matern"`** give isotropic correlation as a
    function of Euclidean plot distance, for irregular or non-rectangular
    layouts. These build the covariance of the observed plots and invert it (the
    marginal form), so they also handle designs with empty cells correctly.
* New `spatial_random = "pspline"` fits a **SpATS-style two-dimensional
  penalised-spline surface** as a random effect (Rodriguez-Alvarez et al. 2018)
  rather than as a residual covariance. The tensor B-spline penalty is improper
  (its null space is constant plus linear trends), so the block is supplied in
  the mixed-model reparameterisation: only the positive-eigenvalue space is kept,
  giving a proper random effect whose surface shrinks away as the smoothing
  parameters grow. Available in both evaluators and both optimisers.
* New `compare_spatial_models()` scores one design under several spatial models
  side by side -- the spatial analogue of `sensitivity_varcomp()`. A design that
  wins under every model is a safe choice; one that wins under a single model is
  a warning.

## Validation and coupling (P3)

* `met_evaluate_alpha_efficiency()` / `met_evaluate_famoptg_efficiency()` now
  also return `A_efficiency_rel`, the classical **relative** A-efficiency in
  (0, 1] (vs an orthogonal CRD on the same plots). The existing `A_efficiency`
  / `$A` fields are documented as the inverse criterion, not relative efficiency.
* New `sensitivity_varcomp()`: sweep the residual/genetic variance ratio and
  report how across-TPE PEV and CDmean respond, so a design choice can be judged
  for robustness to variance misspecification.
* The genuine contrast-based CDmean of Rincent et al. (2012) is available via
  `select_individuals()`; the evaluators' per-line quantity is now correctly
  labelled as reliability (see P0 CDmean fix) rather than mislabelled as the
  Rincent CDmean.

* New `met_information()`: the across-environment (MET-level) coupled
  information matrix. It combines the allocation incidence (Level 1) with a
  within-environment efficiency factor (Level 2) under a G x E model
  (`Sigma_E ⊗ G`) and returns the across-TPE mean PEV and CDmean. This backs the
  "coupled two-level" claim with code; the fixed-effect absorption is verified
  to match a full mixed-model-equations solve.
* New `simulate_met()`: simulates a MET with a genuine across-environment
  genetic correlation structure, predicts by GBLUP, and reports realised
  selection accuracy and genetic gain -- linking A/D/CDmean design proxies to
  the outcome they are meant to improve (replaces the previous no-G×E
  vignette demonstration).

## Decision-point framework — completion (Colmant et al. 2026)

* `simulate_met()` gains a `bv_target` argument (`"across_tpe"`,
  `"environment_specific"`, `"mega_environment"`) so designs can be scored for
  broad *or* specific adaptation. This can change which strategy wins.
* New `sparsity_grid()`: sweep individuals-available x tested-per-environment x
  number-of-environments under a plot budget, scoring realised accuracy and gain
  by simulation — the paper's largest finding, made explorable.
* New `build_environment_relationship()`: construct the environment similarity
  matrix for `select_environments()` from historical genetic correlations or
  enviromic covariates (including unseen environments).
* New `suggest_common_range()`: derive the grid of common-treatment counts for
  `tune_common_treatments()` from the connectivity target (sampling error of a
  genetic correlation), the sparse overlap already present, and the coverage
  reference -- so the sweep range is calculated, not guessed.
* New `tune_common_treatments()`: sweep the number of common treatments and
  report four curves -- accuracy with the environment covariance assumed known,
  accuracy with it re-estimated each replicate, cross-environment connectivity,
  and unique coverage -- to locate the optimum. Because a known-covariance
  criterion almost always prefers zero common treatments, the estimated-covariance
  curve is what reveals the true (often interior) optimum; the two diverge most,
  and the optimum rises, when environments are weakly correlated.
* New `run_design_strategy()`: an end-to-end orchestrator chaining environment
  selection, individual selection, allocation (with optional G x E refinement),
  replication, and evaluation into one call, while each tool remains usable
  alone.

## Decision-point framework (P4, in progress)

* New `select_environments()` (decision 1): choose a representative subset of
  environments from a TPE. Default `"representative"` maximises TPE coverage
  (k-medoid), correcting the spread-only `-q'Dq` criterion (available as
  `"optcontrib"`); `"kmeans"`, `"hclust"`, and `"random"` baselines included.
* New `select_individuals()` (decision 2): choose a representative training
  subset of the TPG by maximising CDmean (Rincent et al. 2012) via an exchange
  algorithm.
* New `optimize_allocation_gxe()` (decision 3): refine an equireplicate
  allocation to minimise across-TPE PEV under the coupled G x E information
  matrix -- the reproducible, CDmean-based replacement for the paper's
  unbounded `-q'Dq` / connectedness objectives, preserving replication and
  environment-size margins.
* New `recommend_replication()` (decision 4): sweep replication levels for an
  allocation, report realised accuracy / gain / PEV via `simulate_met()`, honour
  a seed budget, and recommend the diminishing-returns level -- de-confounding
  replication from the number of unique entries.

## Polish

* `met_prep_famoptg()` now emits the field-dimension auto-correction as an
  informational message under `verbose` (matching its documented "silently
  adjusts" behaviour) instead of a `warning()`.
* `plan_sparse_met_design()` accepts `allocation_method = "equireplicate"`
  (and the `"M4"` alias).

## Allocation (P1)

* Renamed `allocation_method = "balanced_incomplete"` to `"equireplicate"`
  (the old name was a misnomer and has been **removed**, not deprecated;
  `"M4"` remains as an alias): the constructor guarantees equal replication and
  exact requested environment-size margins, not a strict BIBD. The bundled example data argument
  list is likewise renamed `sparse_example_args_equireplicate`.
* New `balance` argument (`"env_pair"`, `"line_pair"`, `"both"`) drives the
  design toward a near-balanced (regular-graph) structure while preserving the
  replication and environment-size margins; returns a `balance_report`.
* `check_equireplicate_feasibility()` (and the retained
  `check_balanced_incomplete_feasibility()`) now report the strict-BIBD
  existence conditions: `bibd_lambda`, `bibd_lambda_integer`, `fisher_ok`, and
  `strict_bibd_possible`.
* New `construct_exact_bibd()` builds a true BIBD via `crossdes` when one
  provably exists (small-J regime).

# OptiSparseMET 0.1.0

## Initial release

`OptiSparseMET` introduces a unified framework for sparse multi-environment
trial (MET) design, jointly addressing treatment allocation across environments
and within-environment field design under shared statistical, genetic, and
logistical constraints.

The package links allocation structure, genetic connectivity, seed
availability, and spatial design assumptions in a single reproducible workflow
compatible with mixed-model inference.

---

## Across-environment allocation

- Added `allocate_sparse_met()` for constructing treatment-by-environment
  incidence matrices using sparse testing principles.
  - Supports `"random_balanced"` (M3) for flexible approximate balance with
    coverage-first guarantees.
  - Supports `"balanced_incomplete"` (M4) for BIBD-inspired uniform
    replication structure with enforced equal replication and equal
    environment sizes across the trial.
  - Accepts `"M3"` and `"M4"` as convenient aliases.
  - Supports unequal environment capacities under `random_balanced`.
  - Supports common treatments to ensure design-based connectivity across
    environments.
  - Returns `$allocation_matrix` (binary treatment-by-environment incidence
    matrix) from which pairwise co-occurrence can be computed post-hoc as
    `out$allocation_matrix %*% t(out$allocation_matrix)`.

- Added `check_balanced_incomplete_feasibility()` to verify the slot identity
  (J* × r = I × k*) before attempting balanced incomplete allocation,
  confirming that equal replication is achievable for the chosen dimensions.

- Added `derive_allocation_groups()` to construct grouping structures from:
  - family membership labels
  - genomic relationship matrix (GRM)
  - pedigree relationship matrix (A matrix)

  Group-guided allocation improves genetic connectedness and stability of
  cross-environment inference.

---

## Capacity and feasibility helpers

- Added `suggest_safe_k()` to propose a safe uniform value for
  `n_test_entries_per_environment` given treatment count, environment count,
  common treatments, and a user-defined buffer.
- Added `min_k_for_full_coverage()` to compute the minimum per-environment
  capacity required for every non-common treatment to be assigned at least
  once.
- Added `warn_if_k_too_small()` to provide a non-fatal diagnostic warning when
  the chosen capacity is insufficient for full treatment coverage.

These helpers prevent the most common failure mode: passing a capacity too
small to assign all treatments before `allocate_sparse_met()` is called.

---

## Seed-aware replication planning

- Added `assign_replication_by_seed()` to determine feasible replication
  levels based on available seed quantities and per-plot seed requirements.
  - Supports `"augmented"`, `"p_rep"`, and `"rcbd_type"` replication modes.
  - Supports `shortage_action` values `"error"`, `"downgrade"`, and
    `"exclude"` for handling treatments with insufficient seed.
  - Returns a role-partitioned list (`p_rep_treatments`,
    `unreplicated_treatments`, `excluded_treatments`) suitable for direct
    input to `met_prep_famoptg()`.

---

## Within-environment field design engines

### Block-based repeated-check designs

- Added `met_prep_famoptg()` for constructing augmented, partially replicated
  (p-rep), and RCBD-type repeated-check block designs.

  Key structural guarantees:
  - Check treatments appear in every block.
  - Replicated (p-rep) treatments appear at most once per block.
  - Unreplicated treatments appear exactly once across the field.
  - Optional genetic dispersion optimisation using GRM or A matrix.
  - Optional within-environment efficiency evaluation (A, D, CDmean).

- Added `met_evaluate_famoptg_efficiency()` for evaluating the statistical
  efficiency of `met_prep_famoptg()` designs under:
  - Fixed or random treatment effects.
  - IID, AR1, or AR1×AR1 residual covariance structures.
  - A-optimality, D-efficiency, CDmean, and mean PEV criteria.
  - Requires `sigma_b2` (block variance) in `varcomp`.

- Added `met_optimize_famoptg()` for criterion-driven optimisation of
  `met_prep_famoptg()` designs via Random Restart.

### Row-column alpha designs

- Added `met_alpha_rc_stream()` for generating alpha row-column stream designs
  suitable for large structured fields.

  Key features:
  - Repeated checks in every incomplete block.
  - Configurable block sizes via `min_block_size` / `max_block_size` or fixed
    `n_blocks_per_rep`.
  - Row-major, column-major, and serpentine traversal orders.
  - Optional genetic grouping from family, GRM, or A matrix.
  - Optional within-environment efficiency evaluation.

- Added `met_evaluate_alpha_efficiency()` for evaluating the statistical
  efficiency of `met_alpha_rc_stream()` designs under the same criteria as
  `met_evaluate_famoptg_efficiency()`.
  - Requires `sigma_rep2` (replicate variance) and `sigma_ib2` (incomplete
    block within replicate variance) in `varcomp`, distinguishing it from the
    block-based evaluator.

- Added `met_optimize_alpha_rc()` for criterion-driven optimisation of
  `met_alpha_rc_stream()` designs via Random Restart, Simulated Annealing, or
  Genetic Algorithm.

---

## Pipeline orchestration

- Added `plan_sparse_met_design()` providing an end-to-end MET design workflow
  that integrates:
  - Across-environment allocation via `allocate_sparse_met()`.
  - Per-environment seed feasibility via `assign_replication_by_seed()`.
  - Environment-specific field design via `met_prep_famoptg()` or
    `met_alpha_rc_stream()`, selected by `design` in `env_design_specs`.
  - Mixed design strategies across environments in a single call.

  Each environment's design engine is specified via `env_design_specs`, a
  named list where `design = "met_prep_famoptg"` or
  `design = "met_alpha_rc_stream"` selects the constructor.

- Added `combine_met_fieldbooks()` to stack environment-level field books
  produced by `met_prep_famoptg()` or `met_alpha_rc_stream()` into a unified
  MET-level field book with standard metadata columns (`Environment`,
  `LocalDesign`, `ReplicationMode`, `SparseMethod`, `IsCommonTreatment`).
  Handles heterogeneous column sets across environments by filling missing
  columns with `NA`.

---

## Statistical framework

Implements the sparse testing identity:

    N = J × r = I × k

making explicit the tradeoff between number of treatments (J), number of
environments (I), replication depth (r), and treatments per environment (k).

Design outputs are compatible with mixed-model analysis:

    y = Xβ + Zg + e

where the covariance of g may be proportional to a GRM or pedigree A matrix,
and the residual covariance structure may be IID, AR1, or AR1×AR1.

Supports design strategies that improve:
- cross-environment genetic connectivity
- G×E estimation stability
- genomic prediction performance (CDmean criterion; Rincent et al. 2012)
- precision of BLUP estimates

---

## Efficiency diagnostics

Within-environment designs support optional efficiency evaluation, summarised
via `plan_sparse_met_design()` across all environments:

- A-criterion and A-efficiency
- D-criterion and D-efficiency
- Mean prediction error variance (mean PEV)
- CDmean (coefficient of determination for genomic prediction)

Metrics are reported in `$efficiency_summary` (long format) and
`$environment_summary` (wide format with `has_efficiency`, `eff_A`, `eff_D`,
`eff_mean_PEV` columns).

---

## Infrastructure

- Initial package structure with all exports documented via roxygen2.
- Bundled example dataset `OptiSparseMET_example_data` with 120 treatments,
  4 environments, GRM, pedigree A matrix, prediction matrix K, seed
  availability data, pre-built allocation argument lists, and
  environment-specific design specifications.
- Vignette describing the statistical framework, two-stage pipeline, and
  worked examples.
- Unit test suite covering all 13 exported functions plus internal helpers.
- pkgdown configuration for the documentation website.
- GitHub Actions workflows for R CMD check and pkgdown deployment.
