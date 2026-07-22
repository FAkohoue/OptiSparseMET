# Answering the five decision points — a corrected, package-native design

Prepared for F. Akohoue — 2026-07-21. Companion to `IMPLEMENTATION_PLAN.md` and
`CRITIQUE_Colmant2026_and_OptiSparseMET.md`.

Goal: make `OptiSparseMET` give a **genuine, criterion-driven, and validated**
answer to each decision point in Colmant et al. (2026), Figure 1 — while
*correcting* the weaknesses in their objective functions and interpretation.
Every method below is presented as a **modified version** of the paper's idea,
with its own package name and an explicit statement of what we changed and why.
Nothing claims to reproduce the paper intact.

The five decision points:
1. Which & how many **locations** best represent the TPE?
2. Which & how many **individuals** best represent the TPG?
3. How to **allocate** a representative set of entries **across environments**?
4. How many **replicates** per entry should be planted?
5. How to **spread entries in plots within** an environment?

---

## Two systemic corrections that make all five answers sound

These fix interpretation errors in the paper that otherwise contaminate every
decision, so we apply them package-wide.

**S1 — Do not force a single "TPE-average" breeding-value target.** The paper
defines accuracy against the *average* true BV across environments, estimated
with a **main-effect-only** model, even though data are simulated with strong
G×E. That structurally favors broad adaptation and makes allocation look
unimportant. Our evaluation harness (`simulate_met()`, `met_information()`) will
support three explicit targets — **across-TPE mean**, **environment-specific**,
and **mega-environment/cluster** BVs — and will fit a **G×E-aware** model
(factor-analytic / unstructured `Sigma_E`), not main-effect-only. Every "which
is better" answer is reported per target so specific adaptation is not silently
penalized.

**S2 — Replace unbounded "spread" proxies with the quantity that determines
accuracy.** The paper's workhorse `-q'Dq` (and the underspecified connectedness
Eqs 2–3) maximize mutual *dissimilarity*, which is not the same as
representativeness and is not the same as prediction reliability. Wherever the
paper uses `-q'Dq`, we (a) keep it as a clearly labeled **spread baseline**, and
(b) add a **CDmean/PEV-based objective** computed from a real information matrix
as the recommended default. This is the single most important scientific upgrade
and is reused across decisions 1, 2, and 3.

---

## Decision 1 — Which & how many locations best represent the TPE?

**Paper's approach.** Select `n` environments with the modified
optimal-contribution function `f(q) = -q'Dq` (D = environment relationship),
compared to k-means, hierarchical clustering, and random; "how many" chosen from
a plot-budget grid.

**Weaknesses (from our critique).** `-q'Dq` maximizes spread, so it can pick
mutually uncorrelated *outlier* environments that miss the TPE's typical
conditions; the method comparison was uncontrolled (different inputs, unequal
compute); `E` was assumed known; and "how many" was never tied to a criterion.

**Our modified answer — `select_environments()` (modified optimal-contribution,
representativeness-corrected).**
- **Objective (default): representativeness-CDmean.** Choose the environment
  subset that maximizes the reliability of predicting the *unselected* remainder
  of the TPE — an environment-level analogue of Rincent-style training-set
  optimization. This directly encodes "represents the TPE," unlike `-q'Dq`.
- **Objective (baseline): `-q'Dq` (spread) + optional mean-representativeness
  term** `w1·coverage − w2·q'Dq`, documented as diversity-maximizing.
- **Baselines retained:** k-means, hclust, random — for honest benchmarking under
  *identical inputs and compute budget* (fixing the paper's uncontrolled
  comparison).
- **"How many": a criterion, not a budget.** Return a representativeness (or
  predicted-TPE-accuracy) curve vs `n`, with an elbow/marginal-gain rule and a
  target-reliability stopping option.
- **`E` construction — `build_environment_relationship()`:** from historical MET
  genetic correlations *or* enviromic covariates (weather/soil/remote sensing),
  with the trait-specificity caveat documented (an enviromic `E` is a general,
  not trait-specific, relationship).

**Naming:** "Modified optimal-contribution environment selection (Colmant et al.
2026, decision 1), representativeness-corrected." **Validation:** on a `D` with a
planted redundant cluster, the corrected objective avoids picking two
near-identical environments that `random` would; "how many" curve is monotone in
predicted-TPE accuracy under the broad-TPE simulation.

---

## Decision 2 — Which & how many individuals best represent the TPG?

**Paper's approach.** Framed as the sparsity dimension (`ia × ipf × noe` grid);
individuals mostly sampled at random, with `-q'Dq` available for selection.

**Weaknesses.** Random sampling is a weak answer to "which individuals
represent"; the representativeness objective for the TPG is undeveloped; "how
many" is only a simulation grid.

**Our modified answer — `select_individuals()` (CDmean training-set
optimization) + `sparsity_grid()`.**
- **Which: CDmean (Rincent et al. 2012) over the GRM.** Select the subset of the
  TPG that maximizes CDmean for predicting the rest of the TPG — the established,
  correct tool for "most representative/informative individuals." This is a
  principled upgrade over random and over raw `-q'Dq`.
- **Diversity alternative:** `-q'Dq(G)` exposed as a founder-contribution/spread
  baseline with the spread-vs-representativeness caveat.
- **How many: `sparsity_grid()`** explores `ia × ipf × noe` under a total-plot
  budget and returns accuracy/Δg from the S1 harness for each cell, so a program
  locates its own sweet spot rather than reading a paper-specific number.

**Naming:** "CDmean-based TPG representative-set selection (modified from Colmant
decision 2; objective from Rincent et al. 2012)." **Validation:** selected set
yields higher predicted accuracy on the held-out TPG than random of equal size;
grid respects the plot budget.

---

## Decision 3 — How to allocate a representative set of entries across environments?

**Paper's approach.** `-q'Dq` with `D = G ⊗ E` optimized by a GA (`evola`); plus
the connectedness objective functions (Eqs 2–3).

**Weaknesses.** Eqs 2–3 are underspecified/non-reproducible ("details on
request"); `-q'Dq(G⊗E)` is again a spread proxy; and the biased target +
main-effect estimation made allocation look nearly irrelevant. Capacity, equal
replication, common treatments, and seed limits are not jointly handled.

**Our modified answer — `allocate_sparse_met(objective = "optcontrib_GxE")`,
optimizing across-TPE CDmean.**
- **Objective (default): across-TPE CDmean/PEV from `met_information()`** — the
  MET-level information matrix combining the allocation incidence `Z` with
  `G ⊗ Sigma_E`. We optimize the quantity that actually governs prediction
  accuracy, **replacing the underspecified Eqs 2–3** with a well-defined,
  reproducible criterion.
- **Structural constraints/objectives:** near-balanced equireplicate concurrence
  and **environment-pair connectivity balance** (from the allocation phase),
  which is the connectedness property the paper gestured at but never balanced.
- **Baseline:** `-q'Dq(G⊗E)` spread objective retained and labeled.
- **Constraints honored jointly:** per-environment capacity, common treatments,
  and **seed availability** (`assign_replication_by_seed()`), making the
  allocation deployable — something the paper omits.
- **Fair evaluation (S1):** assessed under a G×E model and per target, so the
  "does allocation matter?" question is answered without the paper's
  model-misspecification.

**Naming:** "Modified connectedness allocation (Colmant decision 3): CDmean-based
G×E objective replacing Eqs 2–3." **Validation:** the objective improves
across-TPE PEV vs random allocation on a broad-TPE case; on a narrow TPE the gain
shrinks (matching the paper's qualitative finding, now on a sound basis).

---

## Decision 4 — How many replicates per entry should be planted?

**Paper's approach.** Grid of field size (100–900 plots) × replication
(1.0–2.0); conclusion: replication matters little once a GRM is used; ~400 plots
needed under high spatial noise.

**Weaknesses.** Field size is confounded with the number of unique entries;
conclusions are specific to the extreme spatial-noise setting; no seed
constraint; replication reported as a scalar average.

**Our modified answer — `recommend_replication()` (resource- and noise-aware),
built on `assign_replication_by_seed()`.**
- **Criterion-driven, not a static grid.** Given the program's **own** spatial-
  variance level, GRM, and **seed availability**, sweep replication level and
  p-rep fraction and return the accuracy/Δg-per-unit-cost curve; recommend the
  level at the diminishing-returns point.
- **De-confound size vs replication (fixing the paper).** Hold the number of
  unique entries fixed while varying replication, and separately vary field size,
  so the two effects are reported independently.
- **Seed realism.** Recommendations are feasible under seed limits — the
  package's existing strength, now wired into the replication decision.

**Naming:** "Resource-aware replication optimization (modified from Colmant
decision 4; confounding of field size and replication corrected; seed constraints
added)." **Validation:** on a low-noise scenario the recommended replication
differs from the high-noise scenario (showing noise-dependence the single-setting
paper could not); recommendations never exceed seed availability.

---

## Decision 5 — How to spread entries in plots within an environment?

**Paper's approach.** Compare random, A-optimality (`odw`), and a
genetic-by-spatial connectedness design (GbySpat, Eqs 2–3); report accuracy and
compute time.

**Weaknesses.** GbySpat (Eqs 2–3) is not reproducible; A-optimality via `odw` is
expensive. This is the decision the package already covers best — but with the
P0 correctness bugs (AR1 indexing, CDmean scaling) it is not yet trustworthy.

**Our modified answer — the existing within-environment engine, corrected and
extended.**
- **Constructors:** `met_prep_famoptg()` (augmented / p-rep / RCBD) and
  `met_alpha_rc_stream()` (alpha row-column) — already present.
- **Evaluation (after P0 fixes):** A / D / **corrected** CDmean under IID / AR1 /
  AR1×AR1 residuals — a correct, reproducible replacement for the paper's
  A-optimality vs GbySpat comparison.
- **A defined `GbySpat` objective (modified, reproducible):** minimize a
  spatially-weighted genomic-relationship clustering penalty — i.e. spread
  genetically similar lines apart in the field — or, equivalently, **directly
  maximize CDmean under an AR1×AR1 residual**. This gives the paper's
  "genetic-by-spatial" idea an explicit, documented formula instead of
  "details on request."
- **Search:** RS / SA / GA (`met_optimize_*`), reporting compute time like the
  paper so the accuracy-vs-cost trade-off is transparent.

**Naming:** "Genetic-by-spatial within-field optimization (modified from Colmant
decision 5): explicit spatially-weighted relationship objective / AR1×AR1
CDmean." **Validation:** the closed-form A/D oracles (P2.1) and the AR1 test
(P0.1) guarantee the criteria are correct; optimized designs beat random and the
family-clustered adverse control on realized accuracy (matching the paper's
direction, now bug-free).

---

## How the five answers chain into one pipeline

`plan_sparse_met_design()` will orchestrate the corrected decisions end to end,
each step feeding the next and each reporting its own criterion value:

```
select_environments()      # D1: which/how many locations  (repr-CDmean)
        v
select_individuals()       # D2: which/how many individuals (CDmean)  + sparsity_grid()
        v
allocate_sparse_met(objective="optcontrib_GxE")  # D3: across-env allocation (across-TPE CDmean)
        v
recommend_replication()    # D4: replication x field size x seed
        v
met_prep_famoptg() / met_alpha_rc_stream() + met_optimize_*()   # D5: within-field spatial
        v
combine_met_fieldbooks()   # unified MET field book
```

Because each step optimizes a **defined, validated** criterion (and the whole
chain can be scored by `simulate_met()` for realized accuracy/Δg under a G×E
model), the package answers each of the paper's questions *genuinely* — and,
where the paper's objective or interpretation was weak, answers it **better**,
under an honest "modified from Colmant et al. (2026)" label.

---

## Traceability to the implementation plan

| Decision | New/《changed》 functions | Plan items |
|---|---|---|
| D1 locations | `select_environments()`, `build_environment_relationship()` | P4.1 (+ P3.2 metric) |
| D2 individuals | `select_individuals()`, `sparsity_grid()` | P4.2 (+ P3.3) |
| D3 allocation | `allocate_sparse_met(objective="optcontrib_GxE")` | P4.3 (+ P1, P3.2) |
| D4 replication | `recommend_replication()` 《`assign_replication_by_seed()`》 | new item P4.4 (+ P3.3) |
| D5 within-field | 《`met_prep_famoptg`/`met_alpha_rc_stream`/`met_optimize_*`》, `GbySpat` objective | P0.1–P0.4, P2, P3.4 |
| Systemic S1/S2 | `simulate_met()`, `met_information()` | P3.2, P3.3 |

(Add **P4.4 `recommend_replication()`** to `IMPLEMENTATION_PLAN.md` Phase P4.)
