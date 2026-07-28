# Five decisions in robust sparse multi-environment trial design

## Purpose

Colmant et al. (2026) organise multi-environment trial design around five
decisions: which environments to use, which individuals to test, how to
allocate them across environments, how much to replicate them, and how to
arrange them within fields. `OptiSparseMET` retains this useful structure but
replaces proxy objectives where a direct information criterion is available.

The target population of environments (TPE) is the set of environments for
which the breeding programme seeks prediction and selection decisions. The
target population of genotypes (TPG) is the set of candidate genotypes for
which those decisions are required. A genomic relationship matrix (GRM), or a
pedigree numerator relationship matrix when markers are unavailable, describes
genetic relatedness.

Two principles apply to all five decisions.

1. Prediction-error variance (PEV) and the mean coefficient of determination
   (CDmean) are preferred when the required covariance model is defensible.
2. Structural uncertainty is propagated explicitly. Historical MET responses
   may calibrate trait-specific environmental covariance. In their absence,
   the package does not ask the user to invent kernel weights.

## Decision 1: which environments should represent the TPE?

`select_environments()` accepts an aligned environmental similarity matrix and
provides five methods.

- `representative`, the default, maximises coverage of the complete TPE through
  a facility-location criterion.
- `optcontrib` implements the Colmant spread criterion \(-q'Dq\). It selects
  mutually dissimilar environments and is retained as a diversity baseline.
- `kmeans` and `hclust` select medoids from clustering solutions.
- `random` provides a control.

The distinction between spread and representativeness is substantive.
Minimising \(q'Dq\) can select atypical extremes that are mutually dissimilar
but collectively poor representatives of the central TPE. The package
therefore does not label the spread solution as automatically optimal.

The number of environments should be selected from a capacity or
representativeness curve, not assumed in advance. `suggest_site_capacity()`,
`sparsity_grid()`, and `pareto_designs()` support this comparison under fixed
partner-site capacities and actual plot costs.

### Environmental relationships

`build_environment_kernels()` constructs separate weather, soil, management,
and geographic kernels. Cross-modality kernels may represent
weather-by-soil, weather-by-management, soil-by-management, or a declared
higher-order interaction. These are distinct hypotheses: similarity under
main effects neither proves nor excludes similarity under interactions.

When valid historical MET responses are available,
`calibrate_environment_covariance()` may estimate trait-specific covariance and
the associated calibration workflow may support empirical kernel
contributions. Without historical responses, the
identity matrix is the central genetic covariance and the environmental
kernels remain separate, unweighted maximin sensitivity scenarios.

### Environmental groups

`infer_mega_environments()` compares hierarchical clustering and a finite-descent
k-means algorithm across an admissible range of group counts. It evaluates
silhouette separation, minimum group size, algorithm agreement, and
reproducibility across year- or modality-specific relationship matrices.

The result is labelled `stable`, `provisional`, or `unstable`. An unstable
solution returns one unpartitioned group. Moreover, an environmental cluster
becomes a genetic mega-environment only when historical response data support
the corresponding genotype-by-environment pattern.

## Decision 2: which individuals should represent the TPG?

`select_individuals()` uses the relationship matrix to choose a training or
testing subset. Its CDmean method targets prediction reliability for the
unselected TPG. Its spread method retains the \(-q'Dq\) logic as a diversity
baseline, and random selection provides a control.

The number of individuals is evaluated jointly with environment count,
replication, and field capacity through `sparsity_grid()` and simulation. This
joint comparison prevents a nominal gain from larger fields being confused
with a gain from more unique genotypes.

Additive relationship is the default. Dominance relationship is included only
when the biological target requires it, principally for hybrids or testcross
specific combining ability.

## Decision 3: how should individuals be allocated across environments?

`allocate_sparse_met()` first constructs a feasible incidence matrix.
`random_balanced` prioritises coverage under unequal site sizes.
`equireplicate` gives each sparse treatment the same environment replication
and reproduces the named site loads when the degree sequence is realisable.
Line-pair and environment-pair balancing improve the achievable concurrence
structure without claiming a strict balanced incomplete block design where one
cannot exist.

Common treatments provide direct connectivity. Their number is not hard-coded.
`suggest_common_treatments()` supplies a feasible initial set, whereas
`optimize_common_treatments()` compares nested counts and selects:

1. the number of common treatments;
2. their identities; and
3. their binary presence in every environment.

Pairwise connectivity is the number of distinct shared treatments. Repeated
plots of one genotype do not increase it.

`optimize_allocation_gxe()` then refines the feasible allocation by treatment
swaps that preserve row and column margins. The objective is the
across-environment mean PEV from `met_information()`. This direct criterion is
preferred to applying \(-q'Dq\) to a genetic-by-environment product matrix,
because the direct information matrix incorporates actual replication, local
precision, and the covariance target.

## Decision 4: how much replication should each treatment receive?

`recommend_replication()` compares replication levels, partial replication
fractions, and field sizes under the programme's variance and seed assumptions.
`plan_sparse_met_design()` assigns local replication after binary allocation,
using only the remaining seed inventory and physical plots. Local replication
may differ by genotype and site, but it does not alter connectivity.

Replication is charged in physical plots and grams of seed. Each treatment has
one network-wide seed inventory. It consumes seed only at sites where it is
present; an absent treatment-site cell consumes 0 g. A common treatment is
charged at every site. The default reserve is 0 g, and a positive reserve is
used only when the programme declares one.

## Decision 5: how should treatments be arranged within each field?

The within-environment engine supports augmented, partially replicated,
randomised complete block, and alpha row-column structures. Local optimisation
may use A-, D-, or CDmean-based criteria under independent, first-order
autoregressive, or separable row-by-column residual structures.

The correct local design depends on field geometry, residual correlation,
family structure, treatment relationship, and replication. The delivered
fieldbook therefore determines both physical plot accounting and the local
information contribution used by the MET-level evaluation.

## Integrated decision sequence

`plan_sparse_met_design()` coordinates the operational sequence:

1. validate treatment, environment, relationship, capacity, and seed inputs;
2. interpret historical or no-history environmental evidence;
3. infer environmental strata and apply the stability fallback;
4. determine feasible common and sparse treatment structures;
5. allocate treatments without spending seed inventory more than once;
6. construct local fieldbooks;
7. reconstruct plot and seed ledgers from those fieldbooks; and
8. compare feasible designs under common covariance and uncertainty scenarios.

This sequence converts the five conceptual decisions into one auditable design
for planting. The final choice remains a breeder's decision, but its statistical
precision, environmental assumptions, seed use, and field cost are explicit.
