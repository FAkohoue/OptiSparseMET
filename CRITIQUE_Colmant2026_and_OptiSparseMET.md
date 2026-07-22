# Critical review: Colmant et al. (2026) and the OptiSparseMET package

Prepared for F. Akohoue — 2026-07-21

Scope: (A) a critical read of *Colmant, Pita & Covarrubias-Pazaran (2026), "Some
objective functions and ideas to optimize experimental designs in artificial
selection programs," Crop Science 66:e70337*; and (B) an assessment of what is
wrong, weak, or missing in `OptiSparseMET` if the goal is a package that
quantitative-genetics experts in plant breeding would trust and adopt.

---

## Part A — Critique of the paper

The paper is a *framework and simulation* paper, not a software paper. It
decomposes MET design into five sequential decision points and proposes an
objective function for each, validated by stochastic simulation. Its strengths
are the clean decision-point framing, the modified optimal-contribution idea
(−q′Dq) for environment/entry selection, and honest self-reported limitations.
The weaknesses that matter for anyone building on it:

**1. Sequential decomposition is greedy by construction.** The authors
optimize each decision independently and concede the result is not globally
optimal. There is no feedback between stages: environments chosen at step 1
constrain entry allocation at step 3, which constrains within-field design at
step 5, but information never flows backward. The paper offers no bound on how
far the greedy solution sits from the joint optimum, so "sufficiently close to
guide practical decisions" is an assertion, not a demonstrated property.

**2. A single genetic architecture, entirely simulated.** All conclusions rest
on one TPE/TPG: rG ~ MVN(0.3, 0.3), one complex trait, additive + dominance
infinitesimal model, no epistasis, 50 QTL/chromosome, and no real or semi-real
data. The authors acknowledge scenario-specificity but still generalize
("we expect our objective function to work well under both scenarios"). Without
varying the genetic-correlation structure, trait complexity, or heritability,
the robustness of every headline number is untested.

**3. The noise regime is extreme and drives the results.** Plot-level
reliability ≈ 0.1 with spatial variance set to 50·σ_g² and residual 5·σ_g². The
"≈400 plots needed" and "spatial optimization matters" conclusions are partly
artifacts of this deliberately harsh spatial setting; under moderate
heterogeneity the thresholds would differ substantially. This is not wrong, but
it is under-explored and easy to over-read.

**4. −q′Dq is a spread/redundancy criterion, not a representativeness proof.**
Minimizing q′Dq minimizes the sum of pairwise covariances among selected
environments — it maximizes mutual dissimilarity. That is not the same as
maximizing representativeness of the TPE: a set of maximally uncorrelated
environments can be a set of mutually extreme outliers that misses the TPE
centroid. The function ignores each environment's own variance (diagonal) and
its similarity to the *unselected* remainder. The link to "representativeness"
is intuitive but not formally derived.

**5. Uncontrolled method comparisons.** k-means and hierarchical clustering are
fed a table of genetic values; −q′Dq is fed a covariance matrix D; the GA runs
50 generations of optimization while clustering is one-shot. The methods differ
in inputs *and* in computational budget, so the reported ~20% advantage
conflates the objective function with the search effort and the input
representation.

**6. The connectedness objective functions (Eqs 2–3) are not reproducible.**
The genetic-by-spatial functions — the paper's most novel contribution — depend
on a "4D mapping function" for the inter-plot relationship matrix R_f whose
"details ... can be provided as complementary documentation on request." As
written, Eqs 2–3 are linear forms 1′(·)1 with no derivation connecting them to
prediction-error variance or to any inferential quantity. Compared with
A-optimality (which has an unambiguous PEV meaning), their theoretical grounding
is asserted rather than shown, and they cannot be independently implemented from
the text.

**7. The accuracy target biases the conclusions.** "True breeding value across
the TPE" is defined as the *average* of environment-specific true BVs under an
unstructured G×E model, while accuracy is estimated from a **main-effect-only**
model. Two consequences: (i) averaging BVs across environments with genetic
correlations that include negative values penalizes genotypes with strong but
opposite-sign G×E, structurally favoring broad main-effect strategies; and
(ii) fitting a main-effect model to data simulated with strong G×E is a
model misspecification. This plausibly explains the finding that entry
allocation had only "moderate" effect and that "random allocation was often
near-optimal" — both may be partly artifacts of the estimation model, not
properties of the designs.

**8. Weak inferential reporting.** Differences are reported as means ± SE over
100 random-sampling replicates, with no confidence intervals or significance
tests, and treatments in Figure 2 differ in several factors at once (ia, ipf,
noe, sparsity are correlated). Headline phrases like "up to a 20% gain"
report best-case cells rather than an average effect with uncertainty.

**9. Cost and time claims are implementation-bound.** Runtime comparisons
(seconds vs. ~28 min vs. ~56 min) reflect specific packages (`evola`, `odw`) on
one machine, not algorithmic complexity; and the cost model ("cost scales with
trial size") ignores per-location fixed costs, which the authors note only in
passing.

**Bottom line on the paper:** the decision-point framing and the
optimal-contribution idea are genuinely useful, and the limitations section is
commendably candid. But the quantitative claims are single-scenario,
model-misspecified in the estimation step, statistically under-analyzed, and —
for the flagship connectedness functions — not reproducible from the paper
alone. Treat its numbers as illustrative, not as design constants.

---

## Part B — What is wrong, weak, or missing in OptiSparseMET

**Framing first.** `OptiSparseMET` is *not* an implementation of Colmant et al.
(2026). It is built on Montesinos-Lopez et al. (2023): M3/M4 (BIBD) allocation,
within-field designs (augmented, p-rep, RCBD, alpha row-column), A/D/CDmean
evaluation, and RS/SA/GA optimization. That is a legitimate and different scope.
But two things follow: (1) if the intent is to be the reference implementation
of the paper's framework, the core contributions are absent; and (2)
independent of the paper, several correctness and rigor issues currently block
"expert-grade" status.

### B1. Correctness bugs (fix before anything else)

**(a) AR1 / AR1×AR1 residual grid indexing is inconsistent — likely a real
bug.** In both `met_evaluate_alpha_efficiency.R` and
`met_evaluate_famoptg_efficiency.R`:

```r
Qgrid      <- Matrix::kronecker(Qcol, Qrow) * (1 / sigma_e2)
grid_index <- (as.integer(Row) - 1L) * n_cols + as.integer(Column)
Q          <- Qgrid[grid_index, grid_index]
```

`kronecker(Qcol, Qrow)` has inner (fast) dimension `n_rows`, so the linear
index of plot (row r, col c) must be `(c - 1) * n_rows + r`. The code uses
`(r - 1) * n_cols + c`, which is the indexing for `kronecker(Qrow, Qcol)`. The
two agree only when `n_rows == n_cols` and only on the diagonal. For any
non-square field, or for AR1×AR1 with `rho_row != rho_col`, this swaps the row
and column autocorrelation axes and scrambles spatial adjacency, so the residual
precision applied to each plot is wrong. Because every A/D/CDmean value under an
AR1 model flows through this Q, and the SA/GA optimizers can be driven by an AR1
objective, the error propagates into both evaluation and optimization. Add a
unit test on a small non-square grid (e.g. 3×4) comparing `Q` against a manually
constructed AR1×AR1 precision, then fix the index to `(Column-1)*n_rows + Row`
(or build `kronecker(Qrow, Qcol)`).

**(b) "CDmean" is not the Rincent et al. (2012) CDmean.** The code computes
`CDmean = 1 - mean_PEV / sigma_g2` and cites Rincent 2012. The Rincent CDmean is
a *contrast-based* coefficient of determination that removes fixed effects
(via M = I − X(XᵀX)⁻¹Xᵀ) and normalizes prediction-error variance by the
relationship matrix, i.e. CD_i involves K_ii and the variance ratio, not a flat
division by σ_g². For inbred lines under VanRaden scaling, diag(K) is routinely
≠ 1 (can approach 2), so the current quantity is mis-scaled and can fall outside
[0, 1]. Either implement the genuine Rincent contrast form, or rename the output
to "mean reliability (K_ii = 1 approximation)" and drop the Rincent citation.
Expert users will check this, and a mislabeled criterion undermines credibility.

**(c) The test suite is largely tautological.** The efficiency tests assert
`expect_equal(A_efficiency, 1/A_criterion)` and
`expect_equal(CDmean, 1 - mean_PEV/sigma_g2)` — they re-derive the code's own
formulas and therefore validate nothing about statistical correctness. There is
no test comparing A/D/PEV against a closed-form small design, against `sommer`
/ `asreml` / `odw` output, or against a hand-computed mixed-model information
matrix. The current suite would pass even with bug (a) present. Add
ground-truth tests: a tiny RCBD/CRD where A- and D-optimality are known
analytically, and a PEV cross-check against a mixed-model fit.

### B2. Weaknesses in scientific rigor

**(d) Optimized proxies are never linked to realized accuracy or gain.** The
package optimizes A/D/CDmean — design-based surrogates — but ships no harness
showing these translate into higher selection accuracy or genetic gain, which
is the quantity breeders actually care about (and the entire point of the
paper). The only simulation, in the vignette, simulates data with **no G×E**
(one shared genetic value `g_eff` across environments plus an environment main
effect). Under that model, cross-environment connectivity cannot matter the way
the package's connectivity machinery claims, so the vignette both fails to
validate the tool and is quietly misleading. A credible package needs a
simulation module (AlphaSimR-style or GRM-covariance-based *with* a genetic
correlation structure across environments) that reports accuracy and Δg for the
designs it produces.

**(e) Variance components are treated as known and unquestioned.** All
optimality output is conditional on user-supplied σ² and ρ values, with no
estimation path, no default-elicitation guidance, and no sensitivity analysis.
Design optimality is famously sensitive to the assumed variance ratios; an
expert-grade tool should at minimum sweep a plausible range and report whether
the chosen design is robust to misspecification.

**(f) Metaheuristics lack diagnostics and benchmarking.** RS/SA/GA are provided
without convergence traces, restart-to-restart variance reporting, seed
sensitivity, or any benchmark against exact optima on small problems or against
established engines (`odw`, `DiGGer`, `blocksdesign`, `od`). There is currently
no evidence the optimizers find near-optimal designs rather than merely better-
than-random ones.

**(g) The headline "coupled two-level" claim is not implemented.** The README
argues at length that across-environment allocation (Level 1) and
within-environment design (Level 2) are statistically coupled inside a single
information matrix C⁻¹. But the efficiency code evaluates **one field_book at a
time** — it never assembles an across-environment information matrix that
combines the allocation incidence Z with a G×E covariance. So the central
conceptual selling point is asserted in prose but absent from the code; the two
levels are optimized *sequentially*, exactly the thing the README says is
inferior. Delivering on this claim requires a genuine MET-level C⁻¹ (block
structure over environments with a G⊗E or factor-analytic covariance).

### B3. Missing capabilities (relative to the paper and to expert practice)

**(h) Environment selection from a TPE (paper decisions 1–2) is entirely
absent.** No −q′Dq environment selection, no k-means/hclust baselines, no
construction of the environment-relationship matrix E, no enviromic integration.
The package takes the set of environments as given.

**(i) The TPG sparsity dimension (how many individuals to sample) is absent.**
The ia/ipf/noe trade-off that the paper identifies as the single largest driver
of accuracy has no counterpart in the package.

**(j) There is no objective-function-driven, G×E-aware across-environment
allocation.** Allocation is combinatorial BIBD; genetic connectivity is handled
only heuristically (common treatments, family/GRM grouping) and never optimized
against an information or connectedness criterion. The paper's −q′Dq with
D = G⊗E is not available.

**(k) The genetic-by-spatial connectedness functions (paper Eqs 2–3) are not
implemented** — understandable, since the paper does not fully specify them, but
worth noting if alignment with the paper is a goal.

**(l) No benchmarking against balanced designs or against reference software,
no benchmark datasets, and no validation report.** For adoption by experts the
package needs: reproducible comparisons to `odw`/`DiGGer` on shared examples,
at least one real or semi-real dataset, and a documented validation showing the
criteria and optimizers behave correctly.

**(m) Software-maturity gaps.** Version 0.1.0 / lifecycle "experimental"; git
history is a flat series of "update" commits; no NEWS-linked validation; no
multi-trait or unbalanced-evaluation support. None of these is fatal, but
collectively they signal pre-adoption maturity.

### What the package already does well (keep and build on)

Efficient sparse-matrix machinery and a Hutchinson trace estimator for
large designs; clean decoupling of design construction from evaluation;
feasibility/capacity helpers; and — notably — **seed-constrained replication
planning**, which the paper ignores entirely and which is a genuine practical
advantage for real programs. Documentation coverage (README, vignette, man
pages) is strong.

---

## Part C — Addenda from code review (2026-07-21, Q&A follow-ups)

### C1. A / D "efficiency" naming (from the A-optimality direction question)

There is **no directional error**: `A_criterion` is mean pairwise contrast
variance (lower is better, matching Colmant et al.), and the optimizer correctly
minimizes it. But two labeling hazards remain:

- **`$A` / `$D` aliases point to the reciprocal** (`A_efficiency = 1/A_criterion`,
  higher is better), so a user reading `$A` expecting "the A-optimality value
  from the paper" gets the opposite ordering. Deprecate the bare `$A`/`$D`
  aliases; force explicit `A_criterion` vs `A_efficiency`.
- **"Efficiency = 1/criterion" is a misnomer.** Classical A-efficiency is a
  *relative* quantity in (0, 1] (variance under an ideal/orthogonal reference
  ÷ variance under this design). `1/A_criterion` is an unbounded inverse
  variance. Rename (e.g. `A_value`/`A_inv`) or compute a true relative
  efficiency against a reference design.

### C2. A / D / CDmean formula correctness (from the formula question)

- **A-criterion: correct.** `Vsub` = treatment block of C⁻¹ = Var(BLUEs);
  `.pairwise_diff_mean_var` = mean of V_ii+V_jj−2V_ij. Caveat: it averages over
  checks *and* entries together (documentable choice); only as correct as R⁻¹
  (AR1 bug).
- **D-criterion: correct in formulation, numerically fragile.** Pseudo-
  determinant of the centered contrast covariance HVH via
  `sum(log(ev > tol))`, then geometric mean over q = p−1. Risk: the divisor
  `q = p−1` is hardcoded but the number of eigenvalues summed depends on an
  **absolute** `tol = 1e-10`. If the structural zero eigenvalue drifts above
  tol, D collapses toward 0; if a genuine small eigenvalue falls below tol, D is
  inflated. Use a **relative** tolerance (`tol · max(ev)`) or explicitly drop
  exactly one eigenvalue.
- **CDmean: only correct when diag(K) = 1.** Code computes
  `1 − mean_PEV/σ_g²`, i.e. it assumes K_ii = 1 for every line. Per-line
  reliability is `CD_i = 1 − PEV_i/(σ_g²·K_ii)`; for a VanRaden GRM or inbred /
  pedigree material with F > 0, diag(K) ≈ 1+F (up to ~2), so the value is
  mis-scaled and can leave [0,1]. Also it is the mean of individual
  reliabilities, not Rincent's contrast-based generalized CD. Fix: divide by
  `σ_g²·diag(K)_i`; implement the true contrast form or relabel honestly.
- **Rank-deficiency guard (robustness):** in the fixed branch with
  `check_as_fixed = FALSE`, the intercept is retained alongside a full set of
  `Line_` dummies → X is rank-deficient → C singular → `.solve_C` errors.
  Add a guard (drop intercept or reference level).

### C3. Fidelity to Montesinos-López et al. (2023) M3/M4 (from the M3/M4 question)

**M3 (`random_balanced`) — "inspired by", with three deviations from the
paper's M3:**

1. Forces full coverage (every sparse line in ≥ 1 environment). Paper's M3 does
   *not* guarantee coverage (Fig. 1c shows lines in 3/2/1 or possibly 0 envs).
   Documented enhancement.
2. Tolerates unequal environment sizes (`k_vec` per environment). Paper's M3
   keeps every location the same size k = ⌈J·r/I⌉.
3. Lets lines exceed r replications to fill capacity. Paper's M3 caps at r
   (lines already in r locations are removed from the candidate pool).

None is wrong, but the cumulative drift means results are not comparable to the
paper's M3; the label should stay "M3-inspired," as the README mostly does.

**M4 (`balanced_incomplete`) — NOT a BIBD; this is the headline finding.**

- The paper's M4 is a genuine **balanced** incomplete block design: every *pair*
  of lines co-occurs in the same location an equal number of times (constant λ),
  built with `crossdes::find.BIB()`.
- The package's M4 (constructor 7A) is a **greedy least-loaded** assignment that
  guarantees only **equal replication (r)** and **equal block size (k*)**. It
  computes pairwise concurrence (`overlap_matrix = tALLOC · ALLOC`) but **never
  balances or enforces it**. So it is an equireplicate, equal-block-size
  incomplete block design (a regular-graph design at best) — *not* a BIBD. The
  defining BIBD property (constant λ) is neither enforced nor optimized.
- Why this matters scientifically: genomic-prediction connectivity under sparse
  testing is driven primarily by *which pairs of lines co-occur across
  environments* — exactly the concurrence structure a BIBD balances. The package
  delivers the easy half (equal replication) and leaves the half that actually
  drives prediction accuracy (balanced concurrence) unoptimized, while marketing
  it as "BIBD (M4)." This is the most substantive M-method gap.
- Regime note (in the package's favor): a true BIBD generally **cannot exist**
  for breeding sizes — Fisher's inequality needs #blocks ≥ #treatments (I ≥ J),
  and λ = r(k−1)/(J−1) must be a positive integer; with J ≫ I (e.g. 1000 lines,
  4 environments) neither holds. So the relaxation is pragmatic and defensible —
  but then it should be **renamed** (e.g. "equireplicate incomplete" /
  "M4-approx"), not called BIBD. Conversely, in the *small*-J regime where a
  BIBD does exist (e.g. v=7, b=7, r=3, k=3, λ=1), the greedy constructor will
  not reliably find it, whereas the paper's `find.BIB()` would — so the package
  under-delivers there.

**`check_balanced_incomplete_feasibility()` — checks the wrong condition for its
name.** It validates only the slot-count identity J*·r = I·k* (`feasible =
TRUE` iff `difference == 0`). It does **not** check the actual BIBD existence
conditions: (a) λ = r(k*−1)/(J*−1) an integer, and (b) Fisher's I ≥ J*. It will
therefore report "feasible = TRUE" for many parameter sets where no BIBD exists.
Either add the λ-integrality and Fisher checks, or rename to
`check_equireplicate_feasibility()`.

**Evaluation mismatch.** The M-methods (and the paper) exist to serve
*multi-environment, multi-trait* genomic prediction with a G⊗Σ_E interaction
covariance. The package's efficiency evaluation is single-environment,
single-trait, so it never evaluates the across-environment prediction the
allocation is designed for — reinforcing Part B's B2(d)/B2(g) points.

**Credit where due:** the package correctly extends the resource identity to
include common treatments (J* = J − C), which the paper does not cover, and it
verifies slot feasibility before constructing — both good practices.

---

## Priority roadmap to "expert-grade"

1. **Fix bug B1(a)** (AR1 indexing) and add an analytic AR1×AR1 unit test on a
   non-square grid. Highest impact, lowest effort.
2. **Correct or rename CDmean** (B1(b)) to match Rincent et al. (2012), or drop
   the citation and relabel as reliability.
3. **Replace tautological tests with ground-truth tests** (B1(c)): closed-form
   A/D optima and a PEV cross-check against `sommer`/`asreml`.
4. **Add a G×E simulation-validation module** (B2(d)) reporting accuracy and Δg,
   and fix the no-G×E vignette so the demonstrated benefit is real.
5. **Implement the MET-level information matrix** (B2(g)) so the "coupled
   two-level" claim is backed by code — or soften the README to match the
   current per-environment scope.
6. **Add variance-component sensitivity analysis** (B2(e)) and optimizer
   convergence/benchmark diagnostics (B2(f)).
7. **If alignment with Colmant et al. is intended:** add environment selection
   (−q′Dq, with k-means/hclust baselines), the TPG sparsity trade-off, and
   objective-function-driven G×E allocation (B3 h–j).
8. **Benchmark against `odw`/`DiGGer` on shared datasets** and publish a short
   validation report (B3 l).
9. **Rename or fix M4 (C3):** either construct a true BIBD via
   `crossdes::find.BIB()` when feasible (and fall back gracefully), or add a
   concurrence-balancing objective; and rename `balanced_incomplete` /
   `check_balanced_incomplete_feasibility()` to reflect that only equal
   replication + equal block size are guaranteed.
10. **Harden the criteria (C1–C2):** deprecate `$A`/`$D` aliases, relative
    tolerance in the D pseudo-determinant, `K_ii` in CDmean, and a
    rank-deficiency guard in the fixed branch.

The single most important structural point: the package currently ships
*design-optimality proxies* but never demonstrates they yield *better selection*,
and its most-advertised idea (coupled across/within optimization) is not the
thing the code actually computes. Closing those two gaps — validated accuracy
and a true MET information matrix — is what would move it from a competent
design constructor to a tool experts would rely on.
