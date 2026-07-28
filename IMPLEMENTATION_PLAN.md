# OptiSparseMET implementation and verification ledger

## Scope

This ledger records the methodological capabilities implemented for version
0.2.0 and the evidence still required before release. It replaces the earlier
prospective implementation plan. An item described as implemented is present in
the source tree; it is not claimed to have passed the current manual test cycle
until the package maintainer records that result in `VALIDATION.md`.

The implementation is organised around six linked requirements:

1. correct information criteria;
2. feasible and connected across-environment allocation;
3. explicit genetic and environmental covariance;
4. network-wide seed and plot accounting;
5. robust optimisation under structural uncertainty; and
6. reproducible breeder-facing outputs.

## 1. Information criteria and mixed-model calculations

Implemented:

- separable first-order autoregressive residual precision with explicit row and
  column orientation;
- fixed-effect A- and D-criteria based on treatment contrasts;
- relative pseudo-determinant tolerance and rank checks;
- PEV and CDmean from mixed-model equations;
- reliability scaling by each treatment's relationship-matrix diagonal;
- rank-deficiency protection for equivalent fixed-effect parameterisations;
- MET-level information based on the allocation, actual replication, local
  efficiency, genetic relationship, and environmental covariance; and
- Monte Carlo simulation with uncertainty intervals.

Required verification:

- analytic balanced-design oracles;
- independent dense mixed-model-equation solutions;
- non-unit relationship diagonal cases;
- scale and rank stress tests; and
- directional comparisons under deliberately weak and strong connectivity.

## 2. Across-environment allocation

Implemented:

- `random_balanced` coverage-first allocation;
- `equireplicate` allocation with exact named site loads when the degree
  sequence is realisable;
- separate reporting of strict balanced incomplete block design feasibility;
- optional exact BIBD construction in the small-design regime;
- line-pair and environment-pair concurrence improvement;
- genetically and environmentally stratified allocation;
- a global common core based on binary presence, with local replication
  assigned separately from remaining seed; and
- PEV-based genetic-by-environment swap refinement that preserves allocation
  margins.

The common-treatment count is not a package constant.
`suggest_common_treatments()` constructs a feasible starting set.
`optimize_common_treatments()` compares nested counts and selects the number
and identities of common treatments. Each distinct shared genotype counts once;
local repeated plots do not increase connectivity.

Required verification:

- exact row and column margins;
- full treatment coverage;
- concurrence improvement without margin drift;
- common presence at every site;
- local replication that is feasible with remaining seed and physical plots;
- explicit infeasibility when capacities cannot satisfy the constraints; and
- Pareto comparisons based on actual plots.

## 3. Genetic relationships

Implemented:

- marker-based and pedigree-based additive relationship matrices;
- relationship alignment, scaling, tuning, and consensus;
- parental construction of additive hybrid relationship;
- optional marker-based dominance relationship; and
- named additive and dominance kernel combinations.

Additive relationship is the default for all genotype classes. Dominance is
never added automatically and is appropriate only when dominance variation is
part of the prediction target, principally in hybrids and testcrosses.

Required verification:

- marker and treatment name alignment;
- symmetry and positive-semidefinite checks;
- parental-to-hybrid identities;
- additive-only equivalence when dominance is disabled; and
- sensitivity to defensible variance-component combinations.

## 4. Environmental evidence

Implemented:

- environmental data quality control with physical ranges, coverage,
  provenance, and imputation ledgers;
- weather acquisition and crop-stage feature construction;
- depth-aware soil features;
- management and geographic covariates;
- separate weather, soil, management, and geographic kernels;
- functional-analysis-of-variance cross-modality interaction kernels;
- higher-order declared interaction hypotheses;
- bandwidth ensembles and positive-semidefinite normalisation;
- historical MET calibration of trait-specific covariance;
- a no-history branch with identity central covariance and separate unweighted
  maximin scenarios;
- consensus integration without arbitrary user weights;
- automated inference of the number of environmental groups; and
- stable, provisional, and unstable grouping outcomes.

Main-effect and interaction kernels are deliberately separate. Similarity
under weather, soil, or management main effects does not determine similarity
under their interactions. Historical responses may estimate trait-specific
contributions. Without those responses, no synthetic weighted average is
presented as genetic covariance.

Required verification:

- cell-wise precedence when observed and fetched data are merged;
- crop-stage and year alignment;
- no implicit zero imputation;
- invariance to feature-block column count after normalisation;
- reproducible group inference across algorithms and starts;
- finite k-means termination;
- modality-balanced stability assessment; and
- one-group fallback when grouping is unsupported.

## 5. Seed, replication, and fieldbooks

Implemented:

- one network-wide seed inventory per treatment;
- environment-specific grams per plot;
- treatment-specific reserve, with 0 g as the default;
- zero consumption for absent treatment-site cells;
- consumption at every occupied site, including all sites for common
  treatments;
- seed-feasible replication assignment;
- local augmented, partially replicated, randomised complete block, and alpha
  row-column field designs;
- reconstruction of seed and plot ledgers from the delivered fieldbooks; and
- exact reporting of physical plots rather than binary incidence.

Required verification:

- equality between fieldbook-derived and ledger-derived consumption;
- no overspent treatment;
- correct handling of unequal grams per plot;
- preservation of fixed common replication;
- use of unclaimed partial-replication positions by non-common treatments; and
- explicit failure when seed cannot satisfy mandatory presence.

## 6. Robust optimisation and breeder decisions

Implemented:

- variance-component sensitivity analysis;
- common scenario sets across competing designs;
- no-history maximin comparison across environmental kernels;
- nested site-capacity sweeps with fixed partner capacities;
- feasible Pareto frontiers for precision, gain, plots, seed, and
  representativeness;
- environment and individual subset selection;
- replication recommendation; and
- end-to-end assembly of heterogeneous local MET fieldbooks.

The Colmant \(-q'Dq\) objective is retained for environment or individual
spread. It is not used as a substitute for PEV when a defensible coupled
information matrix is available. Environmental grouping is not forced when
the evidence is unstable.

Required verification:

- identical scenario sets and random numbers for design comparisons;
- monotone capacity construction where nesting is requested;
- fixed non-focal site capacities during site-specific sweeps;
- reported optimiser convergence diagnostics;
- reproducible random seeds;
- Pareto dominance checks; and
- additive-only and additive-plus-dominance hybrid comparisons.

## 7. Documentation and release control

The package documentation must define each acronym at first use, distinguish
implemented facts from recommendations, and state the historical and
no-history branches explicitly. Function documentation is maintained in
roxygen source and synchronised manually with `man/*.Rd` until the maintainer
runs `devtools::document()`.

The release sequence is:

1. run the commands in `VALIDATION.md`;
2. resolve every failure and unexpected warning;
3. inspect the end-to-end design ledgers and fieldbooks;
4. compare the additive-only and optional-dominance paths;
5. archive the complete console output;
6. review the final breeder guide and vignette; and
7. release only the source state associated with that evidence.

The implementation is intended to support a breeder's decision, not conceal it.
Every released recommendation must remain traceable to its covariance
assumptions, environmental evidence, plot capacity, and seed inventory.
